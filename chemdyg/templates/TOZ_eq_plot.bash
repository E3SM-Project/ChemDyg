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
#tsDir={{ output }}/post/atm/{{ grid }}/ts/monthly/{{ '%dyr' % (ypf) }}
mkdir -p ts
#mkdir -p figs
#mkdir -p data
#ln -s ${obsDir}/*_obs.nc ./ts
#cd ts
#ln -s ${tsDir}/TCO*.nc ./ts
#ln -s ${cmipDir}/*TCO.nc ./ts
#cd ..
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

ln -s ${input}/${case}.eam.h0.${y1}-01.nc ./ts

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
echo ===== RUN O3 EQ PLOT  =====
echo

# Prepare configuration file
cat > TOZ_eq_plot.py << EOF
#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
import pylab

filename = '${short}'
path = './ts/'
pathout = './'
refer = xr.open_dataset(path+filename+'.eam.h0.${Y1}-01.nc')
file_in = xr.open_mfdataset(path+filename+'.${eamfile}.*')
years = ${y2} - ${y1} +1 

lat = file_in['lat']

TOZ_sel = file_in['TCO'].where((lat < -60), drop=True)+file_in['SCO'].where((lat < -60), drop=True)
AREA_sel = refer['AREA'].where((lat < -60), drop=True)
AREA_64S = refer['AREA'].where((lat < -64), drop=True).sum()
AREA_64S = np.array(AREA_64S)
TOZ_min = TOZ_sel.min(axis=1)

t_period = len(TOZ_sel)
startdate = str(np.array(file_in['time'].dt.year[0]))+'-'+str(np.array(file_in['time'].dt.month[0]))+'-'+str(np.array(file_in['time'].dt.day[0]))
time_range = pd.date_range(start=startdate, periods=t_period, freq='D')
startindex = time_range.is_year_start.tolist().index(True)
endindex = startindex + years*365

O3_64S = np.zeros(years*365)

for i in range(startindex,endindex):
    ii = i - startindex
    d = {'TOZ': np.array(TOZ_sel[i]), 'AREA': np.array(AREA_sel[0])}
    df = pd.DataFrame(data=d)

    df_sort = df.sort_values(by=['TOZ'])
    TOZ_sort = np.array(df_sort['TOZ'])
    AREA_sort = np.array(df_sort['AREA'])

    AREA_sum = np.zeros(len(TOZ_sort))
    AREA_sum[0]=AREA_sort[0]
    for k in range(1,len(TOZ_sort)):
        AREA_sum[k] = AREA_sum[k-1]+AREA_sort[k]

    result = np.where(AREA_sum > AREA_64S)
    TOZ_index = result[0].min()
    O3_64S[ii] = (TOZ_sort[0:TOZ_index-1]*AREA_sort[0:TOZ_index-1]).mean()/(AREA_sort[0:TOZ_index-1].mean())

fig = plt.figure(figsize=(10,5))
plt.plot(time_range[startindex:endindex],O3_64S)
plt.title('Total column ozone conc. with equivalent latitude (64S)')
plt.xlabel('Time')
plt.ylabel('O3 conc. (DU)',fontsize='large')
pylab.savefig(pathout+'TOZ_PDF_timeseries.png', dpi=300)

#-------- climo plot ---------------
O3_array = O3_64S.reshape((years,365))
O3_mean = O3_array.mean(axis=0)
O3_std = O3_array.std(axis=0)
TOZ_array = TOZ_min[startindex:endindex].values.reshape((years,365))
TOZ_mean = TOZ_array.mean(axis=0)
TOZ_std = TOZ_array.std(axis=0)

fig = plt.figure(figsize=(10,5))
plt.plot(time_range[0:365],TOZ_mean, label ='E3SMv2')
plt.fill_between(time_range[0:365],TOZ_mean+TOZ_std,TOZ_mean-TOZ_std, alpha=.5, linewidth=0)
plt.plot(time_range[0:365],O3_mean, label ='TOZ(64S)')
plt.fill_between(time_range[0:365],O3_mean+O3_std,O3_mean-O3_std, alpha=.5, linewidth=0)
plt.xlim(time_range[181],time_range[365])
plt.legend(loc = 'upper left')
plt.title('SH minimum total O3')
plt.xlabel('Time')
plt.ylabel('O3 conc. (DU)',fontsize='large')
pylab.savefig(pathout+'TOZ_PDF_climo.png', dpi=300)

EOF

# Run diagnostics
command="python -u TOZ_eq_plot.py"
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
if [ -d "${f}" ]; then
   mv ./*.png ${www}/${case}/e3sm_chem_diags/plots/
fi

if [ $? != 0 ]; then
  cd ..
  echo 'ERROR (2)' > {{ prefix }}.status
  exit 1
fi

# Copy files
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

