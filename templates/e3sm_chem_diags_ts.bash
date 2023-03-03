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
run_type="{{ run_type }}"
tag="{{ tag }}"
years={{ ypf }}

# Create temporary workdir
workdir=`mktemp -d tmp.${id}.XXXX`
cd ${workdir}

# Create local links to input climo files
#tsDir={{ output }}/post/atm/{{ grid }}/ts/monthly/{{ '%dyr' % (ypf) }}
#mkdir -p ts
#cd ts
#ln -s ${tsDir}/*.nc .
#cd ..
# Create symbolic links to input files
input={{ input }}/{{ input_subdir }}
for (( year=${y1}; year<=${y2}; year++ ))
do
  YYYY=`printf "%04d" ${year}`
  for file in ${input}/${case}.eam.h0.${YYYY}-*.nc
  do
    ln -s ${file} .
  done
  for file in ${input}/${case}.eam.h1.${YYYY}-*.nc
  do
    ln -s ${file} .
  done
done

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
cat > e3sm_chem_diags_budget.py << EOF
#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from calendar import monthrange
import pandas as pd

path = './'
pathout = './'

short_name = '${short}'
startyear = '${y1}'
endyear = '${y2}'

filename = short_name+'.eam.h0.*.nc'
filenameh1 = short_name+'.eam.h1.*.nc'

varname = ["O3","OH","HO2","H2O2","CH2O","CH3O2","CH3OOH","NO","NO2","NO3","N2O5",
           "HNO3","HO2NO2","PAN","CO","C2H6","C3H8","C2H4","ROHO2","CH3COCH3","C2H5O2",
           "C2H5OOH","CH3CHO","CH3CO3","ISOP","ISOPO2","MVKMACR",
           "MVKO2","E90","N2OLNZ","NOYLNZ","CH4LNZ", "H2OLNZ","DMS","SO2","H2SO4","SOAG",
           "NH3","HCL"]

mol_mass = [
   47.9982000000000,        17.0068000000000,        33.0062000000000,
   34.0136000000000,        30.0252000000000,        47.0320000000000,
   48.0394000000000,        30.0061400000000,        46.0055400000000,
   62.0049400000000,        108.010480000000,        63.0123400000000,
   79.0117400000000,        121.047940000000,        28.0104000000000,
   30.0664000000000,        44.0922000000000,        28.0516000000000,
   77.0572000000000,        58.0768000000000,        61.0578000000000,
   62.0652000000000,        44.0510000000000,        75.0424000000000,
   68.1142000000000,        117.119800000000,        70.0878000000000,
   119.093400000000,        47.9982000000000,        28.0134800000000,
   14.0067400000000,        16.0406000000000,        18.0142000000000,
   62.1324000000000,        64.0648000000000,        98.0784000000000,
   12.0110000000000,        17.0289400000000,        36.4601000000000 ]

WD_list = ["C2H5OOH","CH2O","CH3CHO","CH3OOH","H2O2","H2SO4","HNO3","HO2NO2","SO2"]
trop_list =['O3','N2OLNZ','CH4LNZ']
layer = ['','_L1','_L2','_L3','_L4','_trop']

h0_in = xr.open_mfdataset(path+filename)
h1_in = xr.open_mfdataset(path+filenameh1)

variablelist = list(h1_in.keys())

timeperiod = len(h0_in['time'])
startdate = str(np.array(h0_in['time'].dt.year[0]))+'-01-01'

time_range_month = pd.date_range(startdate,  periods=timeperiod, freq='M')
h0_in['time'] = time_range_month
h1_in['time'] = time_range_month

rearth = 6.37122e6 # Earth radius: m

area_rad = h0_in['area']         # radian (ncol)
area = area_rad * rearth * rearth  # m2
time = h0_in['time']
year = np.array(time.dt.year)
month = np.array(time.dt.month)

line_budget = '<h> E3SM chem closure check and burden</h>'
line_budget = line_budget+'<pre> (MSD:averaged concentration after dry deposition; unit:Tg) </pre>'
line_budget = line_budget+'<pre> (VMR:averaged volume mixing ratio; unit: mol/mol) </pre>'
line_budget = line_budget+'<pre> L1: top_of_model to 100 hPa; L2: 100 to 267 hPa; L3: 267 to 856 hPa; L4; 856 hPa to surface</pre>'
line_budget = line_budget + '<pre>'+short_name+'</pre>'
line_budget = line_budget + '<pre>Simulation period: '+ startyear +' - '+ endyear + '</pre>'
line_budget = line_budget+'<pre> Chemistry     L2-norm rel. diff.     MSD          VMR  </pre>'
line_txt = 'Chemistry     L2-norm rel. diff.     MSD          VMR  \n'

fileout_budget = open(pathout+'chem_budget.html',"w")
fileout_txt = open(pathout+'chem_budget.txt',"w")

dt = np.zeros(timeperiod)
for i in range(len(time)):
    dt[i] = monthrange(2001,month[i])[1]*3600*24

dt_array = xr.DataArray(dt, coords=[h0_in['time']], dims=["time"])
mass = h0_in['MASS']

for var in range(len(varname)):
    if varname[var] in trop_list:
        total_layer = len(layer)
    else:
        total_layer = len(layer)-1

    if varname[var]+'_2DMSD' in variablelist:
        print(varname[var])
    else:
        continue

    for ll in range(total_layer):

        MSD = h1_in[varname[var]+'_2DMSD'+layer[ll]] #kg/m2

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
        else:
            total_td = (TDO+TDE+TDI+TRI+TRE+TDA+TDL+TDN+TDU+TDB+TDS+TDD)

        MSD_total = ((MSD*area).sum(axis=1)).mean()
        if ll == 0:
            MSD_vmr = (MSD_total/mol_mass[var]) / (mass[:,:,:].sum(axis=1).sum(axis=1).mean()/28.96)
        elif ll == 1:
            MSD_vmr = (MSD_total/mol_mass[var]) / (mass[:,0:26,:].sum(axis=1).sum(axis=1).mean()/28.96)
        elif ll == 2:
            MSD_vmr = (MSD_total/mol_mass[var]) / (mass[:,26:38,:].sum(axis=1).sum(axis=1).mean()/28.96)
        elif ll == 3:
            MSD_vmr = (MSD_total/mol_mass[var]) / (mass[:,38:58,:].sum(axis=1).sum(axis=1).mean()/28.96)
        elif ll == 4:
            MSD_vmr = (MSD_total/mol_mass[var]) / (mass[:,58::,:].sum(axis=1).sum(axis=1).mean()/28.96)

    # closure check (skip the first month)
        td_temp = (total_td[1::,:]*dt_array[1::]).sum(axis=0)
        temp = MSD[0,:]+td_temp
        a = np.square((MSD[-1,:]-temp)*area)
        b = np.square(MSD[-1,:]*area)
        SQDF = np.sqrt(a.sum()/b.sum())

    # write out closure check
        line_budget = line_budget + '<pre> '+ format(varname[var]+layer[ll],'12s')
        line_budget = line_budget + '     '+"{0:+.3e}".format(SQDF.values)
        line_budget = line_budget + '     '+"{0:+.3e}".format(MSD_total.values*1.e-9)
        if ll == 5:
            line_budget = line_budget + '     ----------' +'</pre>'
        else:
            line_budget = line_budget + '     '+"{0:+.3e}".format(np.array(MSD_vmr)) + '</pre>'

        line_txt = line_txt + format(varname[var]+layer[ll],'12s')
        line_txt = line_txt + '     '+"{0:+.3e}".format(SQDF.values)
        line_txt = line_txt + '     '+"{0:+.3e}".format(MSD_total.values*1.e-9)
        if ll == 5:
            line_txt = line_txt + '     ----------' +'\n'
        else:
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(MSD_vmr)) + '\n'

fileout_budget.write(line_budget)
fileout_txt.write(line_txt)

fileout_budget.close()
fileout_txt.close()

EOF

# Run diagnostics
command="python -u e3sm_chem_diags_budget.py"
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

