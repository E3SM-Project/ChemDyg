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
# diagnostics_base_path is set by zppy using the mache package
noaaDir="{{ diagnostics_base_path }}/observations/Atm/ChemDyg_inputs"
results_dir=${tag}_${Y1}-${Y2}

# Create temporary workdir
workdir=`mktemp -d tmp.${id}.XXXX`
cd ${workdir}

# Create local links to input climo files
tsDir={{ output }}/post/atm/{{ grid }}/ts/monthly/{{ '%dyr' % (ypf) }}
mkdir -p ts
mkdir -p ts/TANG
#cd ts
ln -s ${tsDir}/CO_SRF*.nc ./ts
ln -s ${noaaDir}/CO_NOAA/stations_met_selected.txt ./ts
ln -s ${noaaDir}/CO_NOAA/TANG/*.nc ./ts/TANG
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
cat > NOAA_CO_comparison.py << EOF
#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from calendar import monthrange
from scipy.stats import linregress
import pandas as pd
import pylab

path = './ts/'
pathout = './'

short_name = '${short}'
startyear = '${Y1}'
endyear = '${Y2}'

e3smdataname = 'CO_SRF_'+startyear+'01_'+endyear+'12.nc'
file_in = xr.open_dataset(path+e3smdataname)

timeperiod = len(file_in['time'])
startdate = str(np.array(file_in['time'].dt.year[0]))+'-'+str(np.array(file_in['time'].dt.month[0]))+'-01'
time_range_month = pd.date_range(startdate,  periods=timeperiod, freq='M')
file_in['time'] = time_range_month

CO = file_in['CO_SRF']*1.e9
lat = file_in['lat']
lon = file_in['lon']

#read filenames from NOAA
fid= open(path+'stations_met_selected.txt','r')

count = 0
Sta = []
lat0 = []
lon0 = []
alt0 = []
mnoaa = []
file2 = []
for line in fid:
    columns = line.split()
    Sta.append(columns[0])
    lat0.append(columns[2])
    lon0.append(columns[3])
    alt0.append(columns[4])
    mnoaa.append(columns[5])
    file2.append(columns[6])

for n in range(len(Sta)):
    filename = file2[n]
    file = xr.open_dataset(path+filename)
    timeperiod_noaa = len(file['time'])
    startdate_noaa = str(np.array(file['time'].dt.year[0]))+'-'+str(np.array(file['time'].dt.month[0]))+'-01'
    time_range_noaa = pd.date_range(startdate_noaa,  periods=timeperiod_noaa, freq='M')
    file['time'] = time_range_noaa

    CO_noaa = file['CO']
    lat_noaa = file['lat']
    lon_noaa = file['lon']
    site_code = file['site_code']

    new_lon = lon_noaa
    if lon_noaa < 0:
        new_lon = lon_noaa+360

    CO_sel_1D = np.zeros(len(time_range_noaa))
    for t in range(len(time_range_noaa)):
        if time_range_noaa[t] >= time_range_month[0]:
            CO_sel_1D[t] = CO.sel(time = time_range_noaa[t], lat=lat_noaa,lon=new_lon,method ='nearest')#(lev,lat,lon)
            if time_range_noaa[t]> time_range_month[-1]:
                CO_sel_1D[t::] = 'nan'
                break
        else:
            CO_sel_1D[t] = 'nan'

    nmonth_noaa = np.arange(0,len(CO_noaa),1)
    mask_noaa = ~np.isnan(CO_noaa)
    slope_noaa, intercept, r_value, p_value, std_err = linregress(nmonth_noaa[mask_noaa], CO_noaa[mask_noaa])
    lin_noaa = nmonth_noaa*slope_noaa+intercept
    lin_noaa_xa = xr.DataArray(lin_noaa, coords=[time_range_noaa], dims=["time"])
    diff_noaa = CO_noaa - lin_noaa_xa

    nmonth = np.arange(0,len(CO_sel_1D),1)
    mask = ~np.isnan(CO_sel_1D)
    if (len(nmonth[mask]) != 0):
        slope_e3sm, intercept, r_value, p_value, std_err = linregress(nmonth[mask], CO_sel_1D[mask])
        lin_e3sm = nmonth*slope_e3sm+intercept
        lin_e3sm_xa = xr.DataArray(lin_e3sm, coords=[time_range_noaa], dims=["time"])
        diff = CO_sel_1D - lin_e3sm_xa

        # plotting
        fig, (ax1,ax2) = plt.subplots(2, 1,figsize=(10, 5))
        ax1.plot(CO_noaa['time'],CO_sel_1D,'k')
        ax1.plot(CO_noaa['time'][mask],lin_e3sm_xa[mask],'k--')
        ax1.plot(CO_noaa['time'],CO_noaa,'r')
        ax1.plot(CO_noaa['time'][mask_noaa],lin_noaa_xa[mask_noaa],'r--')
        ax1.set_title('Surface CO at '+Sta[n]+' (Lat '+str(lat_noaa[0].values)+', Lon '+str(lon_noaa[0].values)+')')
        line1 = 'E3SM mean:'+str(np.round(CO_sel_1D[mask].mean(),2))
        line2 = 'E3SM trend:'+str(np.round(slope_e3sm*12,2))+' ppb/yr'
        line3 = 'NOAA mean:'+str(np.round(CO_noaa[mask_noaa].mean().values,2))
        line4 = 'NOAA trend:'+str(np.round(slope_noaa*12,2))+' ppb/yr'
        ax1.legend([line1,line2, line3,line4])

        ax2.plot(CO_noaa['time'][mask],diff[mask],'k')
        ax2.plot(CO_noaa['time'][mask_noaa],diff_noaa[mask_noaa],'r')
        ax2.set_xlabel('time')
        #ax2.set_ylim(-60,60)
        ax2.set_ylabel('Anomalies')
        pylab.savefig(pathout+'NOAA_CO_'+Sta[n]+'.png', dpi=600)
    else:
        # plotting
        fig, (ax1,ax2) = plt.subplots(2, 1,figsize=(10, 5))
        ax1.plot(CO_noaa['time'],CO_noaa,'r')
        ax1.plot(CO_noaa['time'][mask_noaa],lin_noaa_xa[mask_noaa],'r--')
        ax1.set_title('Surface CO at '+Sta[n]+' (Lat '+str(lat_noaa[0].values)+', Lon '+str(lon_noaa[0].values)+')')
        line3 = 'NOAA mean:'+str(np.round(CO_noaa[mask_noaa].mean().values,2))
        line4 = 'NOAA trend:'+str(np.round(slope_noaa*12,2))+' ppb/yr'
        ax1.legend([line3,line4])
        print('E3SM data is not avaiable within this time period')
        ax2.plot(CO_noaa['time'][mask_noaa],diff_noaa[mask_noaa],'r')
        ax2.set_xlabel('time')
        #ax2.set_ylim(-60,60)
        ax2.set_ylabel('Anomalies')
        pylab.savefig(pathout+'NOAA_CO_'+Sta[n]+'.png', dpi=600)

EOF

# Run diagnostics
command="python -u NOAA_CO_comparison.py"
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
mv *.png ${www}/${case}/e3sm_chem_diags/plots/
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

