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
ncfile_save="{{ ncfile_save }}"
# diagnostics_base_path is set by zppy using the mache package
cmipDir="{{ diagnostics_base_path }}/observations/Atm/ChemDyg_inputs"
if [[ "${ncfile_save}" == "true" ]]; then
   results_dir={{ output }}/post/atm/ncfiles
   if [[ -d ${results_dir} ]]; then
      echo "directory exists."
   else
      mkdir -p ${results_dir}
   fi
fi

# Create temporary workdir
workdir=`mktemp -d tmp.${id}.XXXX`
cd ${workdir}

# Create local links to input climo files
tsDir={{ output }}/post/atm/{{ grid }}/ts/monthly/{{ '%dyr' % (ypf) }}
mkdir -p ts
#cd ts
ln -s ${tsDir}/TCO*.nc ./ts
ln -s ${cmipDir}/TCO_CMIP6/*TCO.nc ./ts
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
cat > CMIP6_TCO_comparison.py << EOF
#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from calendar import monthrange
import pandas as pd
import pylab

path = './ts/'
pathout = './'

short_name = '${short}'
startyear = '${Y1}'
endyear = '${Y2}'

filename = 'TCO_'+startyear+'01_'+endyear+'12.nc'

model_in_tco = xr.open_dataset(path+filename)
rearth = 6.37122e6
model_tco = model_in_tco['TCO']
model_area = model_in_tco['area']*rearth*rearth

time = len(model_in_tco['time'])
lat = len(model_in_tco['lat'])
lon = len(model_in_tco['lon'])
o3_time = []

for i in range(time):
    temp1 = 2.1415e-14 * (model_tco[i,:,:]*model_area[:,:]).sum()#DU to Tg
    o3_time.append(temp1)

y = xr.concat(o3_time,dim="time")
E3SM_TCO = xr.DataArray(y, coords=[model_in_tco['time'][0:len(o3_time)]], dims=["time"], name='TCO')

#--------- plot time series TCO ----------

CESM_in = xr.open_mfdataset(path+'CESM_TCO.nc')
GFDL_in = xr.open_mfdataset(path+'GFDL_TCO.nc')
GISS_in = xr.open_mfdataset(path+'GISS_TCO.nc')
MRI_in = xr.open_mfdataset(path+'MRI_TCO.nc')
UKESM_in = xr.open_mfdataset(path+'UKESM_TCO.nc')

timeperiod = len(CESM_in['time'])
timeperiod_year = round(timeperiod/12)

time_range_month = pd.date_range('1850-01-01',  periods=timeperiod, freq='M')
CESM_in['time'] = time_range_month
GFDL_in['time'] = time_range_month
GISS_in['time'] = time_range_month
MRI_in['time'] = time_range_month
UKESM_in['time'] = time_range_month

CESM_TCO = CESM_in['TCO']
GFDL_TCO = GFDL_in['TCO']
GISS_TCO = GISS_in['TCO']
MRI_TCO = MRI_in['TCO']
UKESM_TCO = UKESM_in['TCO']
year_start = (${y1}-1850)*12
t_end = len(E3SM_TCO)+year_start
t_end_year = round(t_end/12)
E3SM_long = np.zeros(t_end)
E3SM_long[0:year_start] = 'NAN'
E3SM_long[year_start:t_end] = E3SM_TCO[:].values
E3SM_xr = xr.DataArray(E3SM_long, coords=[CESM_in['time'][0:t_end]], dims=["time"], name='TCO')
E3SM_xr.to_netcdf(pathout+'E3SM_cmip_'+startyear+'-'+endyear+'.nc')

time_range_year = pd.date_range('1850-01-01',  periods=timeperiod_year, freq='Y')
CESM_ANN = CESM_TCO.groupby('time.year').mean('time')
CESM_std = CESM_TCO.groupby('time.year').std('time')
GFDL_ANN = GFDL_TCO.groupby('time.year').mean('time')
GFDL_std = GFDL_TCO.groupby('time.year').std('time')
GISS_ANN = GISS_TCO.groupby('time.year').mean('time')
GISS_std = GISS_TCO.groupby('time.year').std('time')
MRI_ANN = MRI_TCO.groupby('time.year').mean('time')
MRI_std = MRI_TCO.groupby('time.year').std('time')
UKESM_ANN = UKESM_TCO.groupby('time.year').mean('time')
UKESM_std = UKESM_TCO.groupby('time.year').std('time')
E3SM_ANN = E3SM_xr.groupby('time.year').mean('time')
E3SM_std = E3SM_xr.groupby('time.year').std('time')

#----- plotting -----
fig = plt.figure(figsize=(10,5))
plt.plot(time_range_year,CESM_ANN, label='CESM', linewidth = 1)
plt.fill_between(time_range_year,CESM_ANN+CESM_std,CESM_ANN-CESM_std, alpha=.5, linewidth=0)
plt.plot(time_range_year,GFDL_ANN, label='GFDL',linewidth = 1)
plt.fill_between(time_range_year,GFDL_ANN+GFDL_std,GFDL_ANN-GFDL_std, alpha=.5, linewidth=0)
plt.plot(time_range_year,GISS_ANN, label='GISS', linewidth = 1)
plt.fill_between(time_range_year,GISS_ANN+GISS_std,GISS_ANN-GISS_std, alpha=.5, linewidth=0)
plt.plot(time_range_year,MRI_ANN, label='MRI', linewidth = 1)
plt.fill_between(time_range_year,MRI_ANN+MRI_std,MRI_ANN-MRI_std, alpha=.5, linewidth=0)
plt.plot(time_range_year,UKESM_ANN, label='UKESM', linewidth = 1)
plt.fill_between(time_range_year,UKESM_ANN+UKESM_std,UKESM_ANN-UKESM_std, alpha=.5, linewidth=0)
plt.plot(time_range_year[0:t_end_year],E3SM_ANN, label='E3SM', linewidth = 2, color='k')
plt.fill_between(time_range_year[0:t_end_year],E3SM_ANN+E3SM_std,E3SM_ANN-E3SM_std, alpha=.5, linewidth=0, color='gray')

plt.title('Tropospheric-ozone burden (TCO)')
plt.legend(loc='upper left')
#plt.ylim(200,1)
#plt.axes().set_yscale("log")
plt.xlabel('Time')
plt.ylabel('Tg')

pylab.savefig(pathout+'CMIP_TCO_comparison.png', dpi=600)

EOF

# Run diagnostics
command="python -u CMIP6_TCO_comparison.py"
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
f=${www}/${case}/e3sm_chem_diags_${Y1}_${Y2}/plots/
mkdir -p ${f}
if [ $? != 0 ]; then
  cd ..
  echo 'ERROR (2)' > {{ prefix }}.status
  exit 1
fi

# Copy files
mv *.png ${f}
if [[ "${ncfile_save}" == "true" ]]; then
   mv *.nc ${results_dir}
fi

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

