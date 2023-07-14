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

ln -s ${input}/${case}.eam.h0.${Y1}-01.nc ./ts

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
## as well as first file of next year to ensure that first and last years are complete
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


# Run E3SM chem Diags
echo
echo ===== RUN Temperature EQ PLOT  =====
echo

# Prepare configuration file
cat > Temp_eq_plot.py << EOF
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

lat = file_in['lat']
TOZ_h = file_in['T'][:,13].where((lat < -55), drop=True)
TOZ_m = file_in['T'][:,21].where((lat < -55), drop=True)
TOZ_l = file_in['T'][:,29].where((lat < -55), drop=True)
AREA_sel = refer['AREA'].where((lat < -55), drop=True)
TOZ_h_sel = TOZ_h.sel(time=TOZ_h.time.dt.month.isin([6,7,8,9,10,11,12]))
TOZ_m_sel = TOZ_m.sel(time=TOZ_m.time.dt.month.isin([6,7,8,9,10,11,12]))
TOZ_l_sel = TOZ_l.sel(time=TOZ_l.time.dt.month.isin([6,7,8,9,10,11,12]))

temp_lat_h = np.zeros(10)
temp_lat_m = np.zeros(10)
temp_lat_l = np.zeros(10)
for i in range(0,20,2):
    ilat = i + 60
    AREA_refer = refer['AREA'].where((lat < -ilat), drop=True).sum()
    AREA_refer = np.array(AREA_refer)

    t_period = len(TOZ_h_sel['time'])
    time_range = pd.date_range(start='1/1/2000', periods=t_period, freq='D')
    temp_h = np.zeros(t_period) # height from layer 13 to 29
    temp_m = np.zeros(t_period) # height from layer 13 to 29
    temp_l = np.zeros(t_period) # height from layer 13 to 29

    for t in range(t_period):
        dh = {'TOZ': np.array(TOZ_h_sel[t]), 'AREA': np.array(AREA_sel[0])}
        dfh = pd.DataFrame(data=dh)

        dfh_sort = dfh.sort_values(by=['TOZ'])
        TOZ_h_sort = np.array(dfh_sort['TOZ'])
        AREA_h_sort = np.array(dfh_sort['AREA'])

        AREA_h_sum = np.zeros(len(TOZ_h_sort))
        AREA_h_sum[0]=AREA_h_sort[0]
        for k in range(1,len(TOZ_h_sort)):
            AREA_h_sum[k] = AREA_h_sum[k-1]+AREA_h_sort[k]

        result_h = np.where(AREA_h_sum > AREA_refer)
        TOZ_h_index = result_h[0].min()
        temp_h[t] = (TOZ_h_sort[0:TOZ_h_index-1]*AREA_h_sort[0:TOZ_h_index-1]).mean()/(AREA_h_sort[0:TOZ_h_index-1].mean())
        # ---------
        dm = {'TOZ': np.array(TOZ_m_sel[t]), 'AREA': np.array(AREA_sel[0])}
        dfm = pd.DataFrame(data=dm)

        dfm_sort = dfm.sort_values(by=['TOZ'])
        TOZ_m_sort = np.array(dfm_sort['TOZ'])
        AREA_m_sort = np.array(dfm_sort['AREA'])

        AREA_m_sum = np.zeros(len(TOZ_m_sort))
        AREA_m_sum[0]=AREA_m_sort[0]
        for k in range(1,len(TOZ_m_sort)):
            AREA_m_sum[k] = AREA_m_sum[k-1]+AREA_m_sort[k]

        result_m = np.where(AREA_m_sum > AREA_refer)
        TOZ_m_index = result_m[0].min()
        temp_m[t] = (TOZ_m_sort[0:TOZ_m_index-1]*AREA_m_sort[0:TOZ_m_index-1]).mean()/(AREA_m_sort[0:TOZ_m_index-1].mean())

        # ---------
        dl = {'TOZ': np.array(TOZ_l_sel[t]), 'AREA': np.array(AREA_sel[0])}
        dfl = pd.DataFrame(data=dl)

        dfl_sort = dfl.sort_values(by=['TOZ'])
        TOZ_l_sort = np.array(dfl_sort['TOZ'])
        AREA_l_sort = np.array(dfl_sort['AREA'])

        AREA_l_sum = np.zeros(len(TOZ_l_sort))
        AREA_l_sum[0]=AREA_l_sort[0]
        for k in range(1,len(TOZ_l_sort)):
            AREA_l_sum[k] = AREA_l_sum[k-1]+AREA_l_sort[k]

        result_l = np.where(AREA_l_sum > AREA_refer)
        TOZ_l_index = result_l[0].min()
        temp_l[t] = (TOZ_l_sort[0:TOZ_l_index-1]*AREA_l_sort[0:TOZ_l_index-1]).mean()/(AREA_l_sort[0:TOZ_l_index-1].mean())

    ii = int(i/2)
    temp_lat_h[ii] = temp_h.mean()
    temp_lat_m[ii] = temp_m.mean()
    temp_lat_l[ii] = temp_l.mean()

latlist = np.arange(60,80,2)
fig = plt.figure(figsize=(10,5))
plt.plot(latlist,temp_lat_h)
plt.plot(latlist,temp_lat_m)
plt.plot(latlist,temp_lat_l)
plt.title('Mean temp.(Jul. to Dec.) with equivalent latitude ${Y1}~${Y2}')
plt.xlabel('Lat')
plt.ylabel('Temperature (k)',fontsize='large')
plt.legend( ['14 km alt.','20 km alt.','25 km alt.'])
pylab.savefig(pathout+'temp_PDF_climo.png', dpi=300)

EOF

# Run diagnostics
command="python -u Temp_eq_plot.py"
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

# Change file permissions
chmod -R go+rX,go-w ${f}

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

