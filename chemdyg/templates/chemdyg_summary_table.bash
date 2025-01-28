#!/bin/bash
{% include 'inclusions/slurm_header.bash' %}
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
run_type="{{ run_type }}"
tag="{{ tag }}"

# Create temporary workdir
workdir=`mktemp -d tmp.${id}.XXXX`
cd ${workdir}

# Create local links to input climo files
#tsDir={{ output }}/post/atm/{{ grid }}/ts/monthly/{{ '%dyr' % (ypf) }}
mkdir -p ts
#cd ts
#ln -s ${tsDir}/*.nc .
#cd ..
# Create symbolic links to input files
input={{ input }}/{{ input_subdir }}
eamfile={{ input_files }}
for (( year=${y1}; year<=${y2}; year++ ))
do
  YYYY=`printf "%04d" ${year}`
  for file in ${input}/${case}.eam.h0.${YYYY}-*.nc
  do
    ln -s ${file} ./ts
  done
  for file in ${input}/${case}.${eamfile}.${YYYY}-*.nc
  do
    ln -s ${file} ./ts
  done
done

#{%- if frequency != 'monthly' %}
## For non-monthly input files, need to add the last file of the previous year
#year={{ year1 - 1 }}
#YYYY=`printf "%04d" ${year}`
#mapfile -t files < <( ls ${input}/{{ case }}.{{ input_files }}.${YYYY}-*.nc 2> /dev/null )
#{% raw -%}
#if [ ${#files[@]} -ne 0 ]
#then
#  ln -s ${files[-1]} .
#fi
#{%- endraw %}
# as well as first file of next year to ensure that first and last years are complete
#year={{ year2 + 1 }}
#YYYY=`printf "%04d" ${year}`
#mapfile -t files < <( ls ${input}/{{ case }}.{{ input_files }}.${YYYY}-*.nc 2> /dev/null )
#{% raw -%}
#if [ ${#files[@]} -ne 0 ]
#then
#  ln -s ${files[0]} .
#fi
#{%- endraw %}
#{%- endif %}

#cd ..

# Run E3SM chem Diags
echo
echo ===== RUN E3SM CHEM DIAGS  =====
echo

# Prepare configuration file
cat > chem_summary_table.py << EOF
#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from calendar import monthrange
import pandas as pd

path = './ts/'
pathout = './'

short_name = '${short}'
startyear = '${y1}'
endyear = '${y2}'

filename = short_name+'.eam.h0.*.nc'
filenameh1 = short_name+'.${eamfile}.*.nc'

varname = ["O3","CO","CH4LNZ","NO"]
layer = ['']

h0_in = xr.open_mfdataset(path+filename)
h1_in = xr.open_mfdataset(path+filenameh1)

variablelist = list(h0_in.keys())

timeperiod = len(h0_in['time'])
startdate = str(np.array(h0_in['time'].dt.year[0]))+'-01-01'

time_range_month = pd.date_range(startdate,  periods=timeperiod, freq='ME')
h0_in['time'] = time_range_month
h1_in['time'] = time_range_month

rearth = 6.37122e6 # Earth radius: m
unit_covet = 1.e-9*12 # kg/month -> Tg/year

area_rad = h0_in['area'][0]         # radian (ncol)
area = area_rad * rearth * rearth  # m2
lat = h0_in['lat'][0]
NH = area.where(lat >= 0)
SH = area.where(lat < 0)

lev = h0_in['lev']
time = h0_in['time']

year = np.array(time.dt.year)
month = np.array(time.dt.month)
linehead = '<h> E3SM chemistry high-level summary</h>'
linehead = linehead + '<pre>The * sign indicates that the calculation involves the stratosphere, while calculations without the sign primarily pertain to the troposhere. </pre>'
linehead = linehead + '<pre>'+short_name+'</pre>'
linehead = linehead + '<pre>Simulation period: '+ startyear +' - '+ endyear + '</pre>'
line_ann = linehead + '<p> Season: ANN </p>'

line_ann = line_ann + '<pre> Chemistry                                   Global        N. Hemisphere   S. Hemisphere </pre>'

fileout_ann = open(pathout+'chem_summary_table.html',"w")

dt = np.zeros(timeperiod)
for i in range(len(time)):
    dt[i] = monthrange(2001,month[i])[1]*3600*24

dt_array = xr.DataArray(dt, coords=[h0_in['time']], dims=["time"])

STE_time = np.zeros(len(time))
STE_time_NH = np.zeros(len(time))
STE_time_SH = np.zeros(len(time))
for var in range(len(varname)):
    total_layer = len(layer)

    for ll in range(total_layer):

        if varname[var] == 'O3':
            MSD = h1_in[varname[var]+'_2DMSD'] #kg/m2

            SCO = h0_in['SCO'] *2.1415e-14
            TCO = h0_in['TCO'] *2.1415e-14 #DU to Tg
            TDD = h0_in[varname[var]+'_2DTDD'+layer[ll]]
            CIP = h0_in[varname[var]+'_2DCIP'+layer[ll]] #kg/m2/sec
            CIL = h0_in[varname[var]+'_2DCIL'+layer[ll]] #kg/m2/sec
            total_net = CIP-CIL
            TOZ = SCO+TCO

            MSD_total = ((MSD*area).sum(axis=1)).mean() #kg
            TDD_total = (dt*(TDD*area).sum(axis=1)).mean() #kg
            CIP_total = (dt*(CIP*area).sum(axis=1)).mean()
            CIL_total = (dt*(-CIL*area).sum(axis=1)).mean()
            NET       = (dt*((CIP-CIL)*area).sum(axis=1)).mean()
            SCO_total = (SCO*area).sum(axis=1).mean() #kg
            TCO_total = (TCO*area).sum(axis=1).mean() #kg
            TOZ_total = (TOZ*area).sum(axis=1).mean() #kg
            # NH
            MSD_NH = ((MSD*NH).sum(axis=1)).mean() #kg
            TDD_NH = (dt*(TDD*NH).sum(axis=1)).mean() #kg
            CIP_NH = (dt*(CIP*NH).sum(axis=1)).mean()
            CIL_NH = (dt*(-CIL*NH).sum(axis=1)).mean()
            NET_NH = (dt*((CIP-CIL)*NH).sum(axis=1)).mean()
            SCO_NH = (SCO*NH).sum(axis=1).mean() #kg
            TCO_NH = (TCO*NH).sum(axis=1).mean() #kg
            TOZ_NH = (TOZ*NH).sum(axis=1).mean() #kg
            # SH
            MSD_SH = ((MSD*SH).sum(axis=1)).mean() #kg
            TDD_SH = (dt*(TDD*SH).sum(axis=1)).mean() #kg
            CIP_SH = (dt*(CIP*SH).sum(axis=1)).mean()
            CIL_SH = (dt*(-CIL*SH).sum(axis=1)).mean()
            NET_SH = (dt*((CIP-CIL)*SH).sum(axis=1)).mean()
            SCO_SH = (SCO*SH).sum(axis=1).mean() #kg
            TCO_SH = (TCO*SH).sum(axis=1).mean() #kg
            TOZ_SH = (TOZ*SH).sum(axis=1).mean() #kg
            # calculate STE
            for i in range(len(time)):
                MSDt = h1_in[varname[var]+'_2DMSD_trop'][i,:] #kg/m2

                TDBt = h0_in[varname[var]+'_2DTDB_trop'][i,:]
                TDDt = h0_in[varname[var]+'_2DTDD_trop'][i,:]
                TDEt = h0_in[varname[var]+'_2DTDE_trop'][i,:]
                TDIt = h0_in[varname[var]+'_2DTDI_trop'][i,:]
                TDAt = h0_in[varname[var]+'_2DTDA_trop'][i,:]
                TDLt = h0_in[varname[var]+'_2DTDL_trop'][i,:]
                TDNt = h0_in[varname[var]+'_2DTDN_trop'][i,:]
                TDOt = h0_in[varname[var]+'_2DTDO_trop'][i,:]
                TDSt = h0_in[varname[var]+'_2DTDS_trop'][i,:]
                TDUt = h0_in[varname[var]+'_2DTDU_trop'][i,:]

                total_td = (TDOt+TDEt+TDIt+TDAt+TDLt+TDNt+TDUt+TDBt+TDSt+TDDt)

                MSD_total = (MSDt*area).sum()
                td_temp = total_td*dt[i]
                TTD_total = (td_temp*area).sum()

                if i == 0:
                    STE = 'nan'
                    STE_NH = 'nan'
                    STE_SH = 'nan'
                else:
                    temp = MSD_old+td_temp
                    STE = ((MSDt-temp)*area).sum()
                    STE_NH = ((MSDt-temp)*NH).sum()
                    STE_SH = ((MSDt-temp)*SH).sum()

                STE_time[i] = STE
                STE_time_NH[i] = STE_NH
                STE_time_SH[i] = STE_SH

                MSD_old = MSDt

            STE_mean = STE_time[1::].mean()
            STE_mean_NH = STE_time_NH[1::].mean()
            STE_mean_SH = STE_time_SH[1::].mean()

   # write out annual chem tendency
            line_ann = line_ann + '<pre> '+ format('O3 burden (Tg)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TOZ_total))
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TOZ_NH))
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TOZ_SH)) +'</pre>'
            line_ann = line_ann + '<pre> '+ format('O3 surface deposition (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDD_total)*unit_covet) 
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDD_NH)*unit_covet) 
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDD_SH)*unit_covet) +'</pre>'
            line_ann = line_ann + '<pre> '+ format('O3 production (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(CIP_total)*unit_covet) 
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(CIP_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(CIP_SH)*unit_covet) +'</pre>'
            line_ann = line_ann + '<pre> '+ format('O3 loss (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(CIL_total)*unit_covet) 
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(CIL_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(CIL_SH)*unit_covet) +'</pre>'
            line_ann = line_ann + '<pre> '+ format('O3 net change (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(NET)*unit_covet) 
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(NET_NH)*unit_covet) 
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(NET_SH)*unit_covet) +'</pre>'
            line_ann = line_ann + '<pre> '+ format('O3 tropospheric burden (Tg)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TCO_total)) 
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TCO_NH))
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TCO_SH)) +'</pre>'
            line_ann = line_ann + '<pre> '+ format('O3 stratospheric burden (Tg)*','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(SCO_total)) 
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(SCO_NH))
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(SCO_SH)) +'</pre>'
            line_ann = line_ann + '<pre> '+ format('O3 STE (Tg/yr)*','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(STE_mean)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(STE_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(STE_SH)*unit_covet)+'</pre>'

        elif varname[var] == 'CO':
            MSD = h1_in[varname[var]+'_2DMSD'] #kg/m2
            TDS = h0_in[varname[var]+'_2DTDS'+layer[ll]] #kg/m2/sec
            TDD = h0_in[varname[var]+'_2DTDD'+layer[ll]]
            CEP = h0_in[varname[var]+'_2DCEP'+layer[ll]]
            CEL = h0_in[varname[var]+'_2DCEL'+layer[ll]]
            total_net2 = CEP+CEL
    # annunal mean
            MSD_total = ((MSD*area).sum(axis=1)).mean() #kg
            TDS_total = (dt*(TDS*area).sum(axis=1)).mean() #kg
            TDD_total = (dt*(TDD*area).sum(axis=1)).mean() #kg
            CEP_total = (dt*(CEP*area).sum(axis=1)).mean()
            CEL_total = (dt*(CEL*area).sum(axis=1)).mean()
            NET2      = (dt*((CEP+CEL)*area).sum(axis=1)).mean()
           #NH
            MSD_NH    = ((MSD*NH).sum(axis=1)).mean() #kg
            TDS_NH    = (dt*(TDS*NH).sum(axis=1)).mean() #kg
            TDD_NH    = (dt*(TDD*NH).sum(axis=1)).mean() #kg
            CEP_NH    = (dt*(CEP*NH).sum(axis=1)).mean()
            CEL_NH    = (dt*(CEL*NH).sum(axis=1)).mean()
            NET2_NH   = (dt*((CEP+CEL)*NH).sum(axis=1)).mean()
           #SH
            MSD_SH    = ((MSD*SH).sum(axis=1)).mean() #kg
            TDS_SH    = (dt*(TDS*SH).sum(axis=1)).mean() #kg
            TDD_SH    = (dt*(TDD*SH).sum(axis=1)).mean() #kg
            CEP_SH    = (dt*(CEP*SH).sum(axis=1)).mean()
            CEL_SH    = (dt*(CEL*SH).sum(axis=1)).mean()
            NET2_SH   = (dt*((CEP+CEL)*SH).sum(axis=1)).mean()

            line_ann = line_ann + '<pre> '+ format('CO burden (Tg)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(MSD_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(MSD_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(MSD_SH)*unit_covet)+'</pre>'
            line_ann = line_ann + '<pre> '+ format('CO emission (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDS_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDS_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDS_SH)*unit_covet)+'</pre>'
            line_ann = line_ann + '<pre> '+ format('CO surface deposition (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDD_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDD_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDD_SH)*unit_covet)+'</pre>'
            line_ann = line_ann + '<pre> '+ format('CO production (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(CEP_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(CEP_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(CEP_SH)*unit_covet)+'</pre>'
            line_ann = line_ann + '<pre> '+ format('CO loss (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(CEL_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(CEL_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(CEL_SH)*unit_covet)+'</pre>'
            line_ann = line_ann + '<pre> '+ format('CO net change (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(NET2)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(NET2_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(NET2_SH)*unit_covet)+'</pre>'
        elif varname[var] == 'CH4LNZ':
            if varname[var]+'_2DTDS' in variablelist:
                print(varname[var])
            else:
                print(varname[var]+' not exist')
                continue

            MSD = h1_in[varname[var]+'_2DMSD'] #kg/m2
            TDS = h0_in[varname[var]+'_2DTDS'+layer[ll]] #kg/m2/sec
            TDD = h0_in[varname[var]+'_2DTDD'+layer[ll]]
            CEL = h0_in['r_lch4_2D']
    # annunal mean
            MSD_total = ((MSD*area).sum(axis=1)).mean() #kg
            TDS_total = (dt*(TDS*area).sum(axis=1)).mean() #kg
            TDD_total = (dt*(TDD*area).sum(axis=1)).mean() #kg
            CEL_total = (dt*(CEL*area).sum(axis=1)).mean()
           #NH
            MSD_NH = ((MSD*NH).sum(axis=1)).mean() #kg
            TDS_NH = (dt*(TDS*NH).sum(axis=1)).mean() #kg
            TDD_NH = (dt*(TDD*NH).sum(axis=1)).mean() #kg
            CEL_NH = (dt*(CEL*NH).sum(axis=1)).mean()
           #SH
            MSD_SH = ((MSD*SH).sum(axis=1)).mean() #kg
            TDS_SH = (dt*(TDS*SH).sum(axis=1)).mean() #kg
            TDD_SH = (dt*(TDD*SH).sum(axis=1)).mean() #kg
            CEL_SH = (dt*(CEL*SH).sum(axis=1)).mean()

            line_ann = line_ann + '<pre> '+ format('CH4 burden (Tg)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(MSD_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(MSD_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(MSD_SH)*unit_covet)+'</pre>'
            line_ann = line_ann + '<pre> '+ format('CH4 emission (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDS_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDS_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDS_SH)*unit_covet)+'</pre>'
            line_ann = line_ann + '<pre> '+ format('CH4 surface deposition (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDD_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDD_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDD_SH)*unit_covet)+'</pre>'
            line_ann = line_ann + '<pre> '+ format('CH4 production (Tg/yr)','37s')
            line_ann = line_ann + '      --------     '+'</pre>'
            line_ann = line_ann + '<pre> '+ format('CH4 loss (Tg/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(CEL_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(CEL_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(CEL_SH)*unit_covet)+'</pre>'
            line_ann = line_ann + '<pre> '+ format('CH4 net change (Tg/yr)','37s')
            line_ann = line_ann + '      --------     '+'</pre>'
        elif varname[var] == 'NO':
            MSD = h1_in[varname[var]+'_2DMSD'] + h1_in[varname[var]+'2_2DMSD'] #kg/m2
            TDS = h0_in[varname[var]+'_2DTDS'+layer[ll]] + h0_in[varname[var]+'2_2DTDS'+layer[ll]]  #kg/m2/sec
            TDD = h0_in[varname[var]+'_2DTDD'+layer[ll]] + h0_in[varname[var]+'2_2DTDD'+layer[ll]]
            if 'NO_TDLgt' in variablelist:
                LGT = h0_in['NO_TDLgt'] # kg N/m2/sec
                ACF = h0_in['NO2_TDAcf'] 
    # annunal mean
            MSD_total = ((MSD*area).sum(axis=1)).mean() #kg
            TDS_total = (dt*(TDS*area).sum(axis=1)).mean() #kg
            TDD_total = (dt*(TDD*area).sum(axis=1)).mean() #kg
            MSD_NH = ((MSD*NH).sum(axis=1)).mean() #kg
            TDS_NH = (dt*(TDS*NH).sum(axis=1)).mean() #kg
            TDD_NH = (dt*(TDD*NH).sum(axis=1)).mean() #kg
            MSD_SH = ((MSD*SH).sum(axis=1)).mean() #kg
            TDS_SH = (dt*(TDS*SH).sum(axis=1)).mean() #kg
            TDD_SH = (dt*(TDD*SH).sum(axis=1)).mean() #kg
            if 'NO_TDLgt' in variablelist:
                LGT_3d = LGT.copy()
                ACF_3d = ACF.copy()
                LGT_NH = LGT.copy()
                ACF_NH = ACF.copy()
                LGT_SH = LGT.copy()
                ACF_SH = ACF.copy()
                for k in range(len(lev)):
                    LGT_3d[:,k,:] = LGT[:,k,:]*area
                    ACF_3d[:,k,:] = ACF[:,k,:]*area
                    LGT_NH[:,k,:] = LGT[:,k,:]*NH
                    ACF_NH[:,k,:] = ACF[:,k,:]*NH
                    LGT_SH[:,k,:] = LGT[:,k,:]*SH
                    ACF_SH[:,k,:] = ACF[:,k,:]*SH
             
                LGT_total = (dt*LGT_3d.sum(axis=1).sum(axis=1)).mean() 
                ACF_total = (dt*ACF_3d.sum(axis=1).sum(axis=1)).mean() 
                LGT_Ntotal = (dt*LGT_NH.sum(axis=1).sum(axis=1)).mean() 
                ACF_Ntotal = (dt*ACF_NH.sum(axis=1).sum(axis=1)).mean() 
                LGT_Stotal = (dt*LGT_SH.sum(axis=1).sum(axis=1)).mean() 
                ACF_Stotal = (dt*ACF_SH.sum(axis=1).sum(axis=1)).mean() 

            line_ann = line_ann + '<pre> '+ format('NOx burden (Tg N)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(MSD_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(MSD_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(MSD_SH)*unit_covet)+'</pre>'
            line_ann = line_ann + '<pre> '+ format('NOx emission (Tg N/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDS_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDS_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDS_SH)*unit_covet)+'</pre>'
            if 'NO_TDLgt' in variablelist:
                line_ann = line_ann + '<pre> '+ format('NOx lightning emis (Tg N/yr)','37s')
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(LGT_total)*unit_covet)
                line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(LGT_Ntotal)*unit_covet)
                line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(LGT_Stotal)*unit_covet)+'</pre>'
                line_ann = line_ann + '<pre> '+ format('NOx Aircraft emis (Tg N/yr)','37s')
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(ACF_total)*unit_covet)
                line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(ACF_Ntotal)*unit_covet)
                line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(ACF_Stotal)*unit_covet)+'</pre>'
            line_ann = line_ann + '<pre> '+ format('NOx surface deposition (Tg N/yr)','37s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDD_total)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDD_NH)*unit_covet)
            line_ann = line_ann + '       '+"{0:+.3e}".format(np.array(TDD_SH)*unit_covet)+'</pre>'

fileout_ann.write(line_ann)
fileout_ann.close()

EOF

# Run diagnostics
command="python -u chem_summary_table.py"
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
tmp_Y1=`printf "%04d" ${y1}`
tmp_Y2=`printf "%04d" ${y2}`
f=${www}/${case}/e3sm_chem_diags_${tmp_Y1}_${tmp_Y2}/plots/
mkdir -p ${f}
if [ $? != 0 ]; then
  cd ..
  echo 'ERROR (2)' > {{ prefix }}.status
  exit 1
fi

# Copy files
mv *.html ${f}

# Change file permissions
chmod -R go+rX,go-w ${f}

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

