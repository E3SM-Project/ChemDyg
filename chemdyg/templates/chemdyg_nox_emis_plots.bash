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
Y1="{{ '%04d' % (year1) }}"
Y2="{{ '%04d' % (year2) }}"
run_type="{{ run_type }}"
ncfile_save="{{ ncfile_save }}"
if [[ "${ncfile_save}" == "true" ]]; then
   results_dir={{ output }}/post/atm/chemdygfiles
   mkdir -p ${results_dir}
fi

# Create temporary workdir
workdir=`mktemp -d tmp.${id}.XXXX`
cd ${workdir}

# Create local links to input climo files
tsDir={{ output }}/post/atm/{{ grid }}/clim/{{ '%dyr' % (ypf) }}
mkdir -p climo
#cd ts
ln -s ${tsDir}/*ANN*.nc ./climo
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
cat > e3sm_nox_emis_plots.py << EOF
#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
import pylab
from scipy import interpolate
import cartopy.crs as ccrs

path = './climo/'
pathout = './'

short_name = '${short}'
startyear = '${Y1}'
endyear = '${Y2}'

filename = short_name+'_ANN_'+startyear+'01_'+endyear+'12_climo.nc'

file_in = xr.open_dataset(path+filename)

rearth = 6.37122e6
AREA = file_in['area'] * rearth * rearth #m2
lat = file_in['lat']
lon = file_in['lon']
lev = file_in['lev']


NOx_acf = file_in['NO2_TDAcf'][0,:,:,:] *1.e-9*3600*24*365 #kg N/m2/sec -> Tg N/m2/year
NOx_lgt = file_in['NO_TDLgt'][0,:,:,:] *1.e-9*3600*24*365  #kg N/m2/sec -> Tg N/m2/yr

NOx_acf_3d = NOx_acf.copy()
NOx_lgt_3d = NOx_lgt.copy()

for i in range(len(lev)):
    NOx_acf_3d[i,:,:] = NOx_acf[i,:,:] * AREA # Tg N/m2/year ->Tg N/year
    NOx_lgt_3d[i,:,:] = NOx_lgt[i,:,:] * AREA

NOx_acf_1d = NOx_acf_3d.sum(axis=1).sum(axis=1)
NOx_acf_2d = NOx_acf_3d.sum(axis=0)
NOx_lgt_1d = NOx_lgt_3d.sum(axis=1).sum(axis=1)
NOx_lgt_2d = NOx_lgt_3d.sum(axis=0)

# ----- writing ncfile -----
NOx_acf_1d = NOx_acf_1d.assign_attrs(units="Tg N/year", description="NOx aircraft emission")
NOx_acf_2d = NOx_acf_2d.assign_attrs(units="Tg N/year", description="NOx aircraft emission")
NOx_lgt_1d = NOx_lgt_1d.assign_attrs(units="Tg N/year", description="NOx lightning emission")
NOx_lgt_2d = NOx_lgt_2d.assign_attrs(units="Tg N/year", description="NOx lightning emission")
ds1 = NOx_acf_1d.to_dataset(name='NOx_acf_1d')
ds2 = NOx_acf_2d.to_dataset(name='NOx_acf_2d')
ds3 = NOx_lgt_1d.to_dataset(name='NOx_lgt_1d')
ds4 = NOx_lgt_2d.to_dataset(name='NOx_lgt_2d')
ds = xr.merge([ds1, ds2, ds3, ds4])
ds.to_netcdf(pathout+'E3SM_NOx_emission_'+startyear+'-'+endyear+'.nc')

# ----- plotting -----
fig = plt.figure(figsize=(18,12))
ax1 = fig.add_subplot(221)
plt.plot(NOx_lgt_1d,lev) #levels=np.arange(0, 1000, 20),
plt.title('NOx lightning emission (Tg N/yr) \n '+filename)
plt.ylim(1000,0)
plt.xlabel('Emission (Tg N/yr)')
plt.ylabel('Pressure (hPa)')

ax2 = fig.add_subplot(222, projection=ccrs.PlateCarree(central_longitude=180))
crossplot = plt.contourf(lon, lat, NOx_lgt_2d, #levels=np.arange(0, 1000, 20),
                 cmap = 'jet',extend='both',transform=ccrs.PlateCarree())
ax2.coastlines(color = "grey")
plt.title('NOx lightning column emission (Tg N/yr)')
plt.xlabel('Lon')
plt.ylabel('Lat')
cross_colorbar = fig.colorbar(crossplot)
gl = ax2.gridlines(draw_labels=True,crs=ccrs.PlateCarree())
gl.xlines = False
gl.ylines = False
gl.top_labels = False
gl.right_labels = False

ax3 = fig.add_subplot(223)
plt.plot(NOx_acf_1d,lev) #levels=np.arange(0, 1000, 20),
plt.title('NOx aircraft emission (Tg N/yr) \n '+filename)
plt.ylim(1000,0)
plt.xlabel('Emission (Tg N/yr)')
plt.ylabel('Pressure (hPa)')

ax4 = fig.add_subplot(224, projection=ccrs.PlateCarree(central_longitude=180))
crossplot = plt.contourf(lon, lat, NOx_acf_2d, #levels=np.arange(0, 1000, 20),
                   cmap = 'jet',extend='both',transform=ccrs.PlateCarree())
ax4.coastlines(color = "grey")
plt.title('NOx aircraft column emission (Tg N/yr)')
plt.xlabel('Lon')
plt.ylabel('Lat')
cross_colorbar = fig.colorbar(crossplot)
gl2 = ax4.gridlines(draw_labels=True,crs=ccrs.PlateCarree())
gl2.xlines = False
gl2.ylines = False
gl2.top_labels = False
gl2.right_labels = False

pylab.savefig(pathout+'nox_emis_plot.png', dpi=300)
plt.close()

EOF

# Run diagnostics
command="python -u e3sm_nox_emis_plots.py"
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

