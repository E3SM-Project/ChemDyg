#!/bin/bash
{% include 'inclusions/slurm_header.bash' %}
{{ environment_commands }}

# To load custom E3SM Diags environment, comment out line above using {# ... #}
# and uncomment lines below

#module load anaconda3/2019.03
#source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh
#conda activate e3sm_diags_env_dev
module load ncl

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
# diagnostics_base_path is set by zppy using the mache package
obsDir="{{ diagnostics_base_path }}/observations/Atm/ChemDyg_inputs"
ncfile_save="{{ ncfile_save }}"
if [[ "${ncfile_save}" == "true" ]]; then
   results_dir={{ output }}/post/atm/chemdyg_output
   mkdir -p ${results_dir}
fi

# Create temporary workdir
workdir=`mktemp -d tmp.${id}.XXXX`
cd ${workdir}

# Create local links to input climo files
tsDir={{ output }}/post/atm/{{ grid }}/ts/monthly/{{ '%dyr' % (ypf) }}
mkdir -p ts
mkdir -p figs
mkdir -p data
ln -s ${obsDir}/O3_hole/*_obs.nc ./ts

# Create symbolic links to input files
input={{ input }}/{{ input_subdir }}
eamfile={{ input_files }}
for (( year=${y1}; year<=${y2}; year++ ))
do
  YYYY=`printf "%04d" ${year}`
  for file in ${input}/${case}.${eamfile}.${YYYY}-*.nc
  do
    ln -s ${file} ./ts
  done
done

ln -s ${input}/${case}.eam.h0.${Y1}-01.nc ./ts

{%- if frequency != 'monthly' %}
# For non-monthly input files, need to add the last file of the previous year
year={{ year1 - 1 }}
YYYY=`printf "%04d" ${year}`
mapfile -t files < <( ls ${input}/{{ case }}.${eamfile}.${YYYY}-*.nc 2> /dev/null )
{% raw -%}
if [ ${#files[@]} -ne 0 ]
then
  ln -s ${files[-1]} ./ts
fi
{%- endraw %}
# as well as first file of next year to ensure that first and last years are complete
year={{ year2 + 1 }}
YYYY=`printf "%04d" ${year}`
mapfile -t files < <( ls ${input}/{{ case }}.${eamfile}.${YYYY}-*.nc 2> /dev/null )
{% raw -%}
if [ ${#files[@]} -ne 0 ]
then
  ln -s ${files[0]} ./ts
fi
{%- endraw %}
{%- endif %}


# Run E3SM chem Diags
echo
echo ===== RUN E3SM CHEM DIAGS  =====
echo

# Prepare configuration file
cat > O3hole_diags_columns.py << EOF
#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
import pylab
from matplotlib.dates import DateFormatter

filename = '${short}'
path = './ts/'
pathout = './'

refer = xr.open_dataset(path+filename+'.eam.h0.${Y1}-01.nc')
file_in = xr.open_mfdataset(path+filename+'.${eamfile}.*')

obs_min_in = xr.open_dataset(path + 'to3mins_obs.nc')
obs_area_in = xr.open_dataset(path + 'to3areas_obs.nc')
area_avg = obs_area_in['area_avg']
area_std = obs_area_in['area_std']
oz_avg = obs_min_in['oz_avg']
oz_std = obs_min_in['oz_std']

O3_thrd = 220.0  # DU

startyear = ${y1}
endyear = ${y2}
years = endyear - startyear + 1
lat = refer['lat']
rearth = 6.37122e6
AREA_in = refer['area'] * rearth * rearth

TOZ_sel = file_in['TCO'].where((lat < 0), drop=True)+file_in['SCO'].where((lat < 0), drop=True)
AREA_sel = AREA_in.where((lat < 0), drop=True)
TOZ_min = TOZ_sel.min(axis=1)

t_period = len(TOZ_sel)
startdate = str(2000)+'-'+str(np.array(file_in['time'].dt.month[0]))+'-'+str(np.array(file_in['time'].dt.day[0]))
time_range = pd.date_range(start=startdate, periods=t_period, freq='D')
startindex = time_range.is_year_start.tolist().index(True)
endindex = startindex + years*365

O3_area = TOZ_min.copy()

for i in range(startindex,endindex):
    d = {'TOZ': np.array(TOZ_sel[i]), 'AREA': np.array(AREA_sel)}
    df = pd.DataFrame(data=d)

    df_sort = df.sort_values(by=['TOZ'])
    TOZ_sort = np.array(df_sort['TOZ'])
    AREA_sort = np.array(df_sort['AREA'])

    if TOZ_sort.min() < O3_thrd:
        result = np.where(df_sort < O3_thrd)
        TOZ_index = result[0].max()
        O3_area[i] = AREA_sort[0:TOZ_index].sum()
    else:
        O3_area[i] = 0.

O3_area_time = O3_area[startindex:endindex].sel(time=O3_area[startindex:endindex].time.dt.month.isin([7, 8, 9, 10, 11, 12]))
TOZ_min_time = TOZ_min[startindex:endindex].sel(time=TOZ_min[startindex:endindex].time.dt.month.isin([7, 8, 9, 10, 11, 12]))

O3_array = O3_area_time.values.reshape((years,184))
O3_mean = O3_array.mean(axis=0) *1.e-12
O3_std = O3_array.std(axis=0) *1.e-12
TOZ_array = TOZ_min_time.values.reshape((years,184))
TOZ_mean = TOZ_array.mean(axis=0)
TOZ_std = TOZ_array.std(axis=0)

# ----- writiing ncfile -----
oz_avg_xr = oz_avg.to_dataset(name='oz_avg')
oz_std_xr = oz_std.to_dataset(name='oz_std')
area_avg_xr = area_avg.to_dataset(name='area_avg')
area_std_xr = area_std.to_dataset(name='area_std')
TOZ_mean_xr = xr.DataArray(TOZ_mean, coords=[oz_avg['time']], dims=["time"])
TOZ_mean_xr = TOZ_mean_xr.assign_attrs(units="DU", description='Climatological mean of daily SH minimum E3SM total column ozone')
TOZ_std_xr = xr.DataArray(TOZ_std, coords=[oz_avg['time']], dims=["time"])
TOZ_std_xr = TOZ_std_xr.assign_attrs(units="DU", description='Standard deviation of daily SH minimum E3SM total column ozone')
TOZ_mean_xr = TOZ_mean_xr.to_dataset(name='TOZ_mean')
TOZ_std_xr = TOZ_std_xr.to_dataset(name='TOZ_std')
O3_mean_xr = xr.DataArray(O3_mean, coords=[oz_avg['time']], dims=["time"])
O3_mean_xr = O3_mean_xr.assign_attrs(units="million of km2", description='Climatological mean of daily mean E3SM O3 hole area')
O3_std_xr = xr.DataArray(O3_std, coords=[oz_avg['time']], dims=["time"])
O3_std_xr = O3_std_xr.assign_attrs(units="million of km2", description='Standard deviation of daily mean O3 hole area')
O3_mean_xr = O3_mean_xr.to_dataset(name='O3_mean')
O3_std_xr = O3_std_xr.to_dataset(name='O3_std')
ds = xr.merge([oz_avg_xr, oz_std_xr, TOZ_mean_xr, TOZ_std_xr, 
               area_avg_xr, area_std_xr, O3_mean_xr, O3_std_xr]) 
ds.to_netcdf(pathout+'E3SM_O3_hole_${y1}-${y2}.nc')

#-------- climo plot ---------------
npdate = np.array(time_range[182:366])
fig, ax = plt.subplots(figsize=(10, 5))
plt.plot(npdate,oz_avg, label ='Obs.')
plt.fill_between(npdate,oz_avg+oz_std,oz_avg-oz_std, alpha=.5, linewidth=0)
plt.plot(npdate,TOZ_mean, label ='E3SM')
plt.fill_between(npdate,TOZ_mean+TOZ_std,TOZ_mean-TOZ_std, alpha=.5, linewidth=0)
date_form = DateFormatter("%b-%d")
ax.xaxis.set_major_formatter(date_form)
plt.legend(loc = 'upper left')
plt.title('SH minimum total O3 ('+str(startyear)+' - '+str(endyear)+')')
plt.xlabel('Date')
plt.ylabel('O3 conc. (DU)',fontsize='large')
pylab.savefig(pathout+'to3mins_O3hole.png', dpi=600)

fig, ax = plt.subplots(figsize=(10, 5))
plt.plot(npdate,area_avg, label ='Obs.')
plt.fill_between(npdate,area_avg+area_std,area_avg-area_std, alpha=.5, linewidth=0)
plt.plot(npdate,O3_mean, label ='E3SM')
plt.fill_between(npdate,O3_mean+O3_std,O3_mean-O3_std, alpha=.5, linewidth=0)
date_form = DateFormatter("%b-%d")
ax.xaxis.set_major_formatter(date_form)
plt.legend(loc = 'upper left')
plt.title('O3 hole area (million of km2) ('+str(startyear)+' - '+str(endyear)+')')
plt.xlabel('Date')
plt.ylabel('Area (million of km2)',fontsize='large')
pylab.savefig(pathout+'to3areas_O3hole.png', dpi=600)

EOF

# Run diagnostics
command="python O3hole_diags_columns.py"
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
if [ -d "${f}" ]; then
   mv ./*.png ${f}
fi
if [[ "${ncfile_save}" == "true" ]]; then
   mv *.nc ${results_dir}
fi

# Change file permissions
chmod -R go+rX,go-w ${f}

if [ $? != 0 ]; then
  cd ..
  echo 'ERROR (2)' > {{ prefix }}.status
  exit 1
fi

# Copy files
#cp ./data/*.nc ${tsDir}
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

