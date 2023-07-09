#!/bin/bash
{% include 'slurm_header.sh' %}
{{ environment_commands }}

# To load custom E3SM Diags environment, comment out line above using {# ... #}
# and uncomment lines below

#module load anaconda3/2019.03
#source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh
#conda activate e3sm_diags_env_dev

# Turn on debug output if needed
debug={{ debug }}
if [[ "${debug,,}" == "true" ]]; then
  set -x
fi

# Make sure UVCDAT doesn't prompt us about anonymous logging
export UVCDAT_ANONYMOUS_LOG=False

# Script dir
cd {{ scriptDir }}

# Get jobid
id=${SLURM_JOBID}

# Update status file
STARTTIME=$(date +%s)
echo "RUNNING ${id}" > {{ prefix }}.status

# Basic definitions
case="{{ case }}"
short="{{ short_name }}"
www="{{ www }}"
y1={{ year1 }}
y2={{ year2 }}
Y1="{{ '%04d' % (year1) }}"
Y2="{{ '%04d' % (year2) }}"
run_type="{{ run_type }}"
tag="{{ tag }}"

results_dir=${tag}_${Y1}-${Y2}

# Create temporary workdir
workdir=`mktemp -d tmp.${id}.XXXX`
cd ${workdir}

# Create local links to input climo files
tsDir={{ output }}/post/atm/{{ grid }}/clim/{{ '%dyr' % (ypf) }}
mkdir -p climo
#cd climo
ln -s ${tsDir}/*.nc ./climo
#cd ..
# Create symbolic links to input files
#input={{ input }}/{{ input_subdir }}
#for (( year=${y1}; year<=${y2}; year++ ))
#do
#  YYYY=`printf "%04d" ${year}`
#  for file in ${input}/${case}.{{ input_files }}.${YYYY}-*.nc
#  do
#    ln -s ${file} .
#  done
#  for file in ${input}/${case}.{{ input_files2 }}.${YYYY}-*.nc
#  do
#    ln -s ${file} .
#  done
#done

{%- if frequency != 'monthly' %}
# For non-monthly input files, need to add the last file of the previous year
year={{ year1 - 1 }}
YYYY=`printf "%04d" ${year}`
mapfile -t files < <( ls ${input}/{{ case }}.{{ input_files }}.${YYYY}-*.nc 2> /dev/null )
{% raw -%}
if [ ${#files[@]} -ne 0 ]
then
  ln -s ${files[-1]} .
fi
{%- endraw %}
# as well as first file of next year to ensure that first and last years are complete
year={{ year2 + 1 }}
YYYY=`printf "%04d" ${year}`
mapfile -t files < <( ls ${input}/{{ case }}.{{ input_files }}.${YYYY}-*.nc 2> /dev/null )
{% raw -%}
if [ ${#files[@]} -ne 0 ]
then
  ln -s ${files[0]} .
fi
{%- endraw %}
{%- endif %}

#cd ..

# Run E3SM chem Diags
echo
echo ===== RUN E3SM CHEM DIAGS  =====
echo

# Prepare configuration file
cat > e3sm_chem_diags.py << EOF
#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from calendar import monthrange
import pandas as pd

path = './climo/'
pathout = './'

short_name = '${short}'
startyear = '${Y1}'
endyear = '${Y2}'

filename = short_name+'*_ANN_*.nc'
seasons = ['ANN','DJF','MAM','JJA','SON']

varname = ["O3","OH","HO2","H2O2","CH2O","CH3O2","CH3OOH","NO","NO2","NO3","N2O5",
           "HNO3","HO2NO2","PAN","CO","C2H6","C3H8","C2H4","ROHO2","CH3COCH3","C2H5O2",
           "C2H5OOH","CH3CHO","CH3CO3","ISOP","ISOPO2","MVKMACR",
           "MVKO2","E90","N2OLNZ","NOYLNZ","CH4LNZ", "H2OLNZ","DMS","SO2","H2SO4","SOAG",
           "NH3","HCL"]

WD_list = ["C2H5OOH","CH2O","CH3CHO","CH3OOH","H2O2","H2SO4","HNO3","HO2NO2","SO2"]
trop_list =['O3','N2OLNZ','CH4LNZ']
layer = ['','_L1','_L2','_L3','_L4','_trop']

refer_in = xr.open_mfdataset(path+filename)
variablelist = list(refer_in.keys())

rearth = 6.37122e6 # Earth radius: m
unit_covet = 1.e-9*365*24*3600 # kg/sec -> Tg/year

area_rad = refer_in['area']         # radian (ncol)
area = area_rad * rearth * rearth  # m2

for ss in range(len(seasons)):
    linehead = '<h> E3SM chem tendency check (TD? units: Tg/year)</h>'
    linehead = linehead + '<pre> (WD:wet deposition; TD:tendency; TR: tendency after the stratosphere reset) </pre>'
    linehead = linehead + '<pre> (O:processes outside of chemistry; I:implicit solver; R: reaction rate reset; E:explicit solver; A:aero_model_gasaerexch; L:Linoz; N:reset negative values to zero; U:setting upper boundary values; B:setting lower boundary values; S:surface emission; D:dry deposition) </pre>'
    linehead = linehead + '<pre> (L1: top_of_model to 100 hPa; L2: 100 to 267 hPa; L3: 267 to 856 hPa; L4: 856 hPa to surface) </pre>'
    linehead = linehead + '<pre>'+short_name+'</pre>'
    linehead = linehead + '<pre>Simulation period: '+ startyear +' - '+ endyear + '</pre>'

    line_ann = linehead + '<p> Season: '+seasons[ss]+' </p>'
    line_ann = line_ann + '<pre> Chemistry           TDO            TDI            TRI            TDE            TRE            TDA            TDL            TDN            TDU            TDB            TDS             TDD            WD          total_TD   </pre>'
    line_txt = 'Chemistry           TDO            TDI            TRI            TDE            TRE            TDA            TDL            TDN            TDU            TDB            TDS             TDD            WD          total_TD   \n'

    fileout_ann = open(pathout+'chem_clim_'+seasons[ss]+'.html',"w")
    fileout_txt = open(pathout+'chem_clim_'+seasons[ss]+'.txt',"w")

    h0_in = xr.open_mfdataset(path+short_name+'*'+seasons[ss]+'*.nc')

    for var in range(len(varname)-1):
        if varname[var] in trop_list:
            total_layer = len(layer)
        else:
            total_layer = len(layer)-1

        if varname[var]+'_2DTDI' in variablelist:
            print(varname[var])
        else:
            print(varname[var]+' not exist')
            continue

        for ll in range(total_layer):

            TDB = h0_in[varname[var]+'_2DTDB'+layer[ll]] #kg/m2/sec
            TDD = h0_in[varname[var]+'_2DTDD'+layer[ll]]
            TDE = h0_in[varname[var]+'_2DTDE'+layer[ll]]
            TDI = h0_in[varname[var]+'_2DTDI'+layer[ll]]
            TRI = h0_in[varname[var]+'_2DTRI'+layer[ll]]
            TRE = h0_in[varname[var]+'_2DTRE'+layer[ll]]
            TDA = h0_in[varname[var]+'_2DTDA'+layer[ll]]
            TDL = h0_in[varname[var]+'_2DTDL'+layer[ll]]
            TDN = h0_in[varname[var]+'_2DTDN'+layer[ll]]
            TDO = h0_in[varname[var]+'_2DTDO'+layer[ll]]
            TDS = h0_in[varname[var]+'_2DTDS'+layer[ll]]
            TDU = h0_in[varname[var]+'_2DTDU'+layer[ll]]

            if varname[var] in WD_list:
                WD = h0_in['WD'+layer[ll]+'_'+varname[var]]#kg/m2/sec
                total_td = (WD+TDO+TDE+TDI+TRI+TRE+TDA+TDL+TDN+TDU+TDB+TDS+TDD)
                WD_total = ((WD*area).sum(axis=1)).mean()
            else:
                total_td = (TDO+TDE+TDI+TRI+TRE+TDA+TDL+TDN+TDU+TDB+TDS+TDD)

    # annunal mean
            TDO_total = ((TDO*area).sum(axis=1)).mean() #kg
            TDE_total = ((TDE*area).sum(axis=1)).mean()
            TDI_total = ((TDI*area).sum(axis=1)).mean()
            TRI_total = ((TRI*area).sum(axis=1)).mean()
            TRE_total = ((TRE*area).sum(axis=1)).mean()
            TDA_total = ((TDA*area).sum(axis=1)).mean()
            TDL_total = ((TDL*area).sum(axis=1)).mean()
            TDN_total = ((TDN*area).sum(axis=1)).mean()
            TDU_total = ((TDU*area).sum(axis=1)).mean()
            TDB_total = ((TDB*area).sum(axis=1)).mean()
            TDS_total = ((TDS*area).sum(axis=1)).mean()
            TDD_total = ((TDD*area).sum(axis=1)).mean()
            TD_total = ((total_td*area).sum(axis=1)).mean()

    # write out annual chem tendency
            line_ann = line_ann + '<pre> '+ format(varname[var]+layer[ll],'12s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDO_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDI_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TRI_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDE_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TRE_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDA_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDL_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDN_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDU_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDB_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDS_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDD_total)*unit_covet)
            if varname[var] in WD_list:
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(WD_total)*unit_covet)
            else:
                line_ann = line_ann + '     ----------'
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TD_total)*unit_covet)+'</pre>'

            line_txt = line_txt + '     '+format(varname[var]+layer[ll],'12s')
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDO_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDI_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TRI_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDE_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TRE_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDA_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDL_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDN_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDU_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDB_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDS_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDD_total)*unit_covet)
            if varname[var] in WD_list:
                line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(WD_total)*unit_covet)
            else:
                line_txt = line_txt + '     ----------'
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TD_total)*unit_covet) +'\n'

    fileout_ann.write(line_ann)
    fileout_txt.write(line_txt)
    fileout_ann.close()
    fileout_txt.close()
    h0_in.close()

for ss in range(len(seasons)):
    linehead = '<h> E3SM chem tendency check (units: Tg/year)</h>'
    linehead = linehead + '<pre> (TDI/TDE:tendency due to implicit/explicit solver; TEP/TIP:chemistry production rate from explicit/implicit solver before reset; TEL/TIL:chemistry loss rate from explicit/implicit solver before reset; </pre>'
    linehead = linehead + '<pre> TRI/TRE:tendency due to vmr and/or reaction rate reset; CEP/CIP:chemistry production rate from explicit/implicit solver after reset; CEL/CIL:chemistry loss rate from explicit/implicit solver after reset; </pre>'
    linehead = linehead + '<pre> MPP/MPL:Micheal Prather calculation for ozone prodcution and loss; L2 DIFF:L2-norm relative difference) </pre>'
    linehead = linehead + '<pre> Note: Implicit solver for troposphere only </pre>'
    linehead = linehead + '<pre> (L1: top_of_model to 100 hPa; L2: 100 to 267 hPa; L3: 267 to 856 hPa; L4: 856 hPa to surface) </pre>'
    linehead = linehead + '<pre>'+short_name+'</pre>'
    linehead = linehead + '<pre>Simulation period: '+ startyear +' - '+ endyear + '</pre>'
    line_ann = linehead + '<p> Season: '+seasons[ss]+' </p>'

    line_ann = line_ann + '<pre> Chemistry           TDI            TIP            TIL            NET1         L2 DIFF          TRI            CIP            CIL          NET2          L2 DIFF          MPP            MPL            NET3          L2 DIFF        </pre>'
    line_ann = line_ann + '<pre>                                                                (TIP+TIL)     (TDI;NET1)                                                  (CIP+CIL)     (TDI+TRI;NET2)                                (MPP+MPL)     (TDI+TRI;NET3)   </pre>'
    fileout_ann = open(pathout+'chem_prodloss_'+seasons[ss]+'.html',"w")
    h0_in = xr.open_mfdataset(path+short_name+'*'+seasons[ss]+'*.nc')

    for var in range(len(varname)):
        total_layer = len(layer)-1

        for ll in range(total_layer):

            if varname[var] == 'O3':
                TDI = h0_in[varname[var]+'_2DTDI'+layer[ll]] #kg/m2/sec
                TRI = h0_in[varname[var]+'_2DTRI'+layer[ll]] #kg/m2/sec
                TIP = h0_in[varname[var]+'_2DTIP'+layer[ll]]
                TIL = h0_in[varname[var]+'_2DTIL'+layer[ll]]
                CIP = h0_in[varname[var]+'_2DCIP'+layer[ll]]
                CIL = h0_in[varname[var]+'_2DCIL'+layer[ll]]
                MPP = h0_in[varname[var]+'_2DMPP'+layer[ll]]
                MPL = h0_in[varname[var]+'_2DMPL'+layer[ll]]
                total_net1 = TIP-TIL
                total_net2 = CIP-CIL
                total_net3 = MPP-MPL

                TDI_total = ((TDI*area).sum(axis=1)).mean() #kg
                TRI_total = ((TRI*area).sum(axis=1)).mean() #kg
                TIP_total = ((TIP*area).sum(axis=1)).mean()
                TIL_total = ((-TIL*area).sum(axis=1)).mean()
                NET1      = (((TIP-TIL)*area).sum(axis=1)).mean()
                CIP_total = ((CIP*area).sum(axis=1)).mean()
                CIL_total = ((-CIL*area).sum(axis=1)).mean()
                NET2      = (((CIP-CIL)*area).sum(axis=1)).mean()
                MPP_total = ((MPP*area).sum(axis=1)).mean()
                MPL_total = ((-MPL*area).sum(axis=1)).mean()
                NET3      = (((MPP-MPL)*area).sum(axis=1)).mean()
    #closure check 
                td_temp = total_net1.sum(axis=0)
                temp = TDI.sum(axis=0)
                a = np.square((td_temp-temp)*area)
                b = np.square(temp*area)
                SQDF1 = np.sqrt(a.sum()/b.sum())
                td_temp = total_net2.sum(axis=0)
                temp = (TDI+TRI).sum(axis=0)
                a = np.square((td_temp-temp)*area)
                b = np.square(temp*area)
                SQDF2 = np.sqrt(a.sum()/b.sum())
                td_temp = total_net3.sum(axis=0)
                a = np.square((td_temp-temp)*area)
                SQDF3 = np.sqrt(a.sum()/b.sum())
   # write out annual chem tendency
                line_ann = line_ann + '<pre> '+ format(varname[var]+layer[ll],'12s')
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDI_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TIP_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TIL_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(NET1)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(SQDF1))
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TRI_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(CIP_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(CIL_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(NET2)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(SQDF2))
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(MPP_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(MPL_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(NET3)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(SQDF3))+'</pre>'

            elif varname[var] == 'CO':
                TDE = h0_in[varname[var]+'_2DTDE'+layer[ll]] #kg/m2/sec
                TRE = h0_in[varname[var]+'_2DTRE'+layer[ll]] #kg/m2/sec
                TEP = h0_in[varname[var]+'_2DTEP'+layer[ll]]
                TEL = h0_in[varname[var]+'_2DTEL'+layer[ll]]
                CEP = h0_in[varname[var]+'_2DCEP'+layer[ll]]
                CEL = h0_in[varname[var]+'_2DCEL'+layer[ll]]
                total_net4 = TEP+TEL
                total_net5 = CEP+CEL
    # annunal mean
                TDE_total = ((TDE*area).sum(axis=1)).mean() #kg
                TRE_total = ((TRE*area).sum(axis=1)).mean() #kg
                TEP_total = ((TEP*area).sum(axis=1)).mean()
                TEL_total = ((TEL*area).sum(axis=1)).mean()
                NET4      = (((TEP+TEL)*area).sum(axis=1)).mean()
                CEP_total = ((CEP*area).sum(axis=1)).mean()
                CEL_total = ((CEL*area).sum(axis=1)).mean()
                NET5      = (((CEP+CEL)*area).sum(axis=1)).mean()

                td_temp = total_net4.sum(axis=0)
                temp = TDE.sum(axis=0)
                a = np.square((td_temp-temp)*area)
                b = np.square(temp*area)
                SQDF4 = np.sqrt(a.sum()/b.sum())
                td_temp = total_net5.sum(axis=0)
                temp = (TDE+TRE).sum(axis=0)
                a = np.square((td_temp-temp)*area)
                b = np.square(temp*area)
                SQDF5 = np.sqrt(a.sum()/b.sum())

                if ll == 0:
                    line_ann = line_ann + '<pre>   </pre>'
                    line_ann = line_ann + '<pre> Chemistry           TDE            TEP            TEL            NET1         L2 DIFF          TRE            CEP            CEL           NET2          L2 DIFF        </pre>'
                    line_ann = line_ann + '<pre>                                                                (TEP+TEL)     (TDE;NET1)                                                 (CEP+CEL)    (TDE+TRE;NET2)     </pre>'

                line_ann = line_ann + '<pre> '+ format(varname[var]+layer[ll],'12s')
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDE_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TEP_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TEL_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(NET4)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(SQDF4))
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TRE_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(CEP_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(CEL_total)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(NET5)*unit_covet)
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(SQDF5))+'</pre>'

            elif varname[var] == 'r_lch4':
                if ll == 0:
                    rch4 = h0_in[varname[var]+'_2D']
                else:
                    rch4 = h0_in[varname[var]+layer[ll]] #kg/m2/sec

    # annunal mean
                rch4_total = ((rch4*area).sum(axis=1)).mean() #kg

                if ll == 0:
                    line_ann = line_ann + '<pre>   </pre>'
                    line_ann = line_ann + '<pre> Chemistry       reaction rate (loss)      </pre>'

                line_ann = line_ann + '<pre> '+ format('CH4'+layer[ll],'12s')
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(rch4_total)*unit_covet) +'</pre>'

    fileout_ann.write(line_ann)
    fileout_ann.close()
    h0_in.close()

EOF

# Run diagnostics
command="python -u e3sm_chem_diags.py"
time ${command}
if [ $? != 0 ]; then
  cd ..
  echo 'ERROR (1)' > {{ prefix }}.status
  exit 1
fi

# Copy output to web server
echo
echo ===== COPY FILES TO WEB SERVER =====
echo

# Create top-level directory
f=${www}/${case}/e3sm_chem_diags/plots/
mkdir -p ${f}
if [ $? != 0 ]; then
  cd ..
  echo 'ERROR (2)' > {{ prefix }}.status
  exit 1
fi

# Copy files
cp *.html ${www}/${case}/e3sm_chem_diags/plots/
cp *.txt ${www}/${case}/e3sm_chem_diags/plots/

# Change file permissions
chmod -R go+rX,go-w ${www}/${case}/e3sm_chem_diags/plots/

if [ $? != 0 ]; then
  cd ..
  echo 'ERROR (3)' > {{ prefix }}.status
  exit 1
fi
cd ..
if [[ "${debug,,}" != "true" ]]; then
  rm -rf ${workdir}
fi

# Update status file and exit
{% raw %}
ENDTIME=$(date +%s)
ELAPSEDTIME=$(($ENDTIME - $STARTTIME))
{% endraw %}
echo ==============================================
echo "Elapsed time: $ELAPSEDTIME seconds"
echo ==============================================
echo 'OK' > {{ prefix }}.status
exit 0

