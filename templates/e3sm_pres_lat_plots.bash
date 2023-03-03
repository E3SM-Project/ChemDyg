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
referDir="{{reference_data_path}}"
results_dir=${tag}_${Y1}-${Y2}

# Create temporary workdir
workdir=`mktemp -d tmp.${id}.XXXX`
cd ${workdir}

# Create local links to input climo files
tsDir={{ output }}/post/atm/{{ grid }}/clim/{{ '%dyr' % (ypf) }}
mkdir -p climo
#cd ts
ln -s ${tsDir}/*ANN*.nc ./climo
ln -s ${referDir}/v2.LR.amip_0101_ANN_198501_201412_climo.nc ./climo
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
cat > e3sm_pres_lat_plots.py << EOF
#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from calendar import monthrange
import pandas as pd
import pylab
from scipy import interpolate

path = './climo/'
pathout = './'

short_name = '${short}'
startyear = '${Y1}'
endyear = '${Y2}'

filename = short_name+'_ANN_'+startyear+'01_'+endyear+'12_climo.nc'
refername = 'v2.LR.amip_0101_ANN_198501_201412_climo.nc'

file_in = xr.open_dataset(path+filename)
refer_in = xr.open_dataset(path+refername)

tpp = file_in['TROP_P'][0,:,:].mean(axis=1)
tpp_3d = file_in['TROPE3D_P'][0,:,:].mean(axis=1)
tpp_refer = refer_in['TROP_P'][0,:,:].mean(axis=1)

lev = file_in['hyam']*100000+file_in['hybm']*file_in['PS'][0,:,:]
lev_refer = refer_in['hyam']*100000+refer_in['hybm']*refer_in['PS'][0,:,:]

lev_2d = 0.01*lev.mean(axis=2)
lev_refer_2d = 0.01*lev_refer.mean(axis=2)

o3 = file_in['O3'][0,:,:,:]
o3_refer = refer_in['O3'][0,:,:,:]
o3_2d = o3.mean(axis=2)
o3_refer_2d = o3_refer.mean(axis=2)
o3_new = o3_2d.copy()
o3_refer_new = o3_refer_2d.copy()

Q = file_in['Q'][0,:,:,:]*28.96/18
Q_refer = refer_in['Q'][0,:,:,:]*28.96/18
Q_2d = Q.mean(axis=2)
Q_refer_2d = Q_refer.mean(axis=2)
Q_new = Q_2d.copy()
Q_refer_new = Q_refer_2d.copy()

T = file_in['T'][0,:,:,:]
T_refer = refer_in['T'][0,:,:,:]
T_2d = T.mean(axis=2)
T_refer_2d = T_refer.mean(axis=2)
T_new = T_2d.copy()
T_refer_new = T_refer_2d.copy()

theda = T*(100000/lev)**0.286
theda_refer = T_refer*(100000/lev)**0.286
theda_2d = theda.mean(axis=2)
theda_refer_2d = theda_refer.mean(axis=2)
theda_new = theda_2d.copy()
theda_refer_new = theda_refer_2d.copy()

for i in range(180):
    f = interpolate.interp1d(lev_2d[:,i], o3_2d[:,i])
    f1 = interpolate.interp1d(lev_refer_2d[:,i], o3_refer_2d[:,i])
    Qf = interpolate.interp1d(lev_2d[:,i], Q_2d[:,i])
    Qf1 = interpolate.interp1d(lev_refer_2d[:,i], Q_refer_2d[:,i])
    Tf = interpolate.interp1d(lev_2d[:,i], T_2d[:,i])
    Tf1 = interpolate.interp1d(lev_refer_2d[:,i], T_refer_2d[:,i])
    thedaf = interpolate.interp1d(lev_2d[:,i], theda_2d[:,i])
    thedaf1 = interpolate.interp1d(lev_refer_2d[:,i], theda_refer_2d[:,i])

    for k in range(72):
        if o3['lev'][k] < lev_2d[0,i] or o3['lev'][k] > lev_2d[-1,i]:
            o3_new[k,i] = 'nan'
            Q_new[k,i] = 'nan'
            T_new[k,i] = 'nan'
            theda_new[k,i] = 'nan'
        else:
            o3_new[k,i] = f(o3['lev'][k])
            Q_new[k,i] = Qf(o3['lev'][k])
            T_new[k,i] = Tf(o3['lev'][k])
            theda_new[k,i] = thedaf(o3['lev'][k])

        if o3['lev'][k] < lev_refer_2d[0,i] or o3['lev'][k] > lev_refer_2d[-1,i]:
            o3_refer_new[k,i] = 'nan'
            Q_refer_new[k,i] = 'nan'
            T_refer_new[k,i] = 'nan'
            theda_refer_new[k,i] = 'nan'
        else:
            o3_refer_new[k,i] = f1(o3['lev'][k])
            Q_refer_new[k,i] = Qf1(o3['lev'][k])
            T_refer_new[k,i] = Tf1(o3['lev'][k])
            theda_refer_new[k,i] = thedaf1(o3['lev'][k])
diff = o3_new - o3_refer_new
diff_relate = diff/o3_refer_new
Qdiff = Q_new - Q_refer_new
Qdiff_relate = Qdiff/Q_refer_new
Tdiff = T_new - T_refer_new
Tdiff_relate = Tdiff/T_refer_new
thedadiff = theda_new - theda_refer_new
thedadiff_relate = thedadiff/theda_refer_new

# plotting

fig = plt.figure(figsize=(18,12))
plt.subplot(2, 2, 1)
crossplot = plt.contourf(o3['lat'], o3['lev'], o3_new*1e6, levels=np.arange(0, 14, 1),
                        cmap = 'jet',extend='both')
CS = plt.contour(T['lat'], T['lev'], theda_new, linewidths=0.5, colors = 'w',levels=np.arange(200, 500, 20))
plt.clabel(CS, inline=True, fontsize=7)
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m', label = 'TROP3D_chemUCI')
plt.title('O3 conc. (ppm) and potential temp. (K)\n '+filename)
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
plt.legend(loc = 'upper left')
cross_colorbar = fig.colorbar(crossplot)

plt.subplot(2, 2, 2)
crossplot = plt.contourf(o3['lat'], o3['lev'], o3_refer_new*1e6, levels=np.arange(0, 14, 1),
                        cmap = 'jet',extend='both')
CS = plt.contour(T['lat'], T['lev'], theda_refer_new, linewidths=0.5,colors = 'w', levels=np.arange(200, 500, 20))
plt.clabel(CS, inline=True, fontsize=7)
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('O3 conc. (ppm) and potential temp. (K)\n v2.LR.amip_0101')
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')
#plt.legend( ['TOZ'])
#pylab.savefig(pathout+'O3_v2.LR.amip_0101.png', dpi=300)

#fig = plt.figure(figsize=(10,5))
plt.subplot(2, 2, 3)
crossplot = plt.contourf(o3['lat'], o3['lev'], diff*1e6, levels=np.arange(-1, 1, 0.1),
                        cmap = 'bwr',extend='both')
CS = plt.contour(T['lat'], T['lev'], thedadiff, linewidths=0.5, colors = 'grey',levels=np.arange(-5, 5, 1))
plt.clabel(CS, inline=True, fontsize=7)
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m',label = 'TROP3D_chemUCI')
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('O3 conc. difference (ppm) and potential temp. difference (K)')
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')

plt.subplot(2, 2, 4)
crossplot = plt.contourf(o3['lat'], o3['lev'], diff_relate, levels=np.arange(-1, 1, 0.1),
                        cmap = 'bwr',extend='both')
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m',label = 'TROP3D_chemUCI')
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('O3 conc. relative difference')
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')
pylab.savefig(pathout+'O3_pres_lat_plot.png', dpi=300)

# plot2
fig = plt.figure(figsize=(18,12))
plt.subplot(2, 2, 1)
crossplot = plt.contourf(o3['lat'], o3['lev'], o3_new*1e9, levels=np.arange(0, 500, 25),
                        cmap = 'jet',extend='both')
CS = plt.contour(T['lat'], T['lev'], theda_new, linewidths=0.5, colors = 'w',levels=np.arange(200, 500, 20))
plt.clabel(CS, inline=True, fontsize=7)
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m', label = 'TROP3D_chemUCI')
plt.title('O3 conc. (ppb) \n '+filename)
plt.ylim(1000,100)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
plt.legend(loc = 'upper left')
cross_colorbar = fig.colorbar(crossplot)

plt.subplot(2, 2, 2)
crossplot = plt.contourf(o3['lat'], o3['lev'], o3_refer_new*1e9, levels=np.arange(0, 500, 25),
                        cmap = 'jet',extend='both')
CS = plt.contour(T['lat'], T['lev'], theda_refer_new, linewidths=0.5,colors = 'w', levels=np.arange(200, 500, 20))
plt.clabel(CS, inline=True, fontsize=7)
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('O3 conc. (ppb) \n v2.LR.amip_0101')
plt.ylim(1000,100)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')

plt.subplot(2, 2, 3)
crossplot = plt.contourf(o3['lat'], o3['lev'], diff*1e9, levels=np.arange(-100, 100, 10),
                        cmap = 'bwr',extend='both')
CS = plt.contour(T['lat'], T['lev'], thedadiff, linewidths=0.5, colors = 'grey',levels=np.arange(-5, 5, 1))
plt.clabel(CS, inline=True, fontsize=7)
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m',label = 'TROP3D_chemUCI')
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('O3 conc. difference (ppb)')
plt.ylim(1000,100)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')

plt.subplot(2, 2, 4)
crossplot = plt.contourf(o3['lat'], o3['lev'], diff_relate, levels=np.arange(-1, 1, 0.1),
                        cmap = 'bwr',extend='both')
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m',label = 'TROP3D_chemUCI')
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('O3 conc. relative difference')
plt.ylim(1000,100)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')
pylab.savefig(pathout+'O3_pres_lat_plot_mb.png', dpi=300)

# plot 3
fig = plt.figure(figsize=(18,12))
plt.subplot(2, 2, 1)
crossplot = plt.contourf(o3['lat'], o3['lev'], T_new, levels=np.arange(190, 290, 5),
                        cmap = 'jet',extend='both')
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m', label = 'TROP3D_chemUCI')
plt.title('Temperature (K) \n '+filename)
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
plt.legend(loc = 'upper left')
cross_colorbar = fig.colorbar(crossplot)

plt.subplot(2, 2, 2)
crossplot = plt.contourf(o3['lat'], o3['lev'], T_refer_new, levels=np.arange(190, 290, 5),
                        cmap = 'jet',extend='both')
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('Temperature (K) \n v2.LR.amip_0101')
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')

plt.subplot(2, 2, 3)
crossplot = plt.contourf(o3['lat'], o3['lev'], Tdiff, levels=np.arange(-5, 5, 0.5),
                        cmap = 'bwr',extend='both')
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m',label = 'TROP3D_chemUCI')
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('Temperature difference (K)')
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')

plt.subplot(2, 2, 4)
crossplot = plt.contourf(o3['lat'], o3['lev'], Tdiff_relate, levels=np.arange(-1, 1, 0.1),
                        cmap = 'bwr',extend='both')
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m',label = 'TROP3D_chemUCI')
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('Temperature relative difference')
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')
pylab.savefig(pathout+'T_pres_lat_plot.png', dpi=300)

# plot 4
fig = plt.figure(figsize=(18,12))
plt.subplot(2, 2, 1)
crossplot = plt.contourf(o3['lat'], o3['lev'], Q_new*1.e6, levels=np.arange(1, 10, 0.2),
                        cmap = 'jet',extend='both')
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m', label = 'TROP3D_chemUCI')
plt.title('Specific humidity (ppm) \n '+filename)
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
plt.legend(loc = 'upper left')
cross_colorbar = fig.colorbar(crossplot)

plt.subplot(2, 2, 2)
crossplot = plt.contourf(o3['lat'], o3['lev'], Q_refer_new*1.e6, levels=np.arange(1, 10, 0.2),
                        cmap = 'jet',extend='both')
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('Specific humidity (ppm) \n v2.LR.amip_0101')
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')

plt.subplot(2, 2, 3)
crossplot = plt.contourf(o3['lat'], o3['lev'], Qdiff*1.e6, levels=np.arange(-5, 5, 0.5),
                        cmap = 'bwr',extend='both')
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m',label = 'TROP3D_chemUCI')
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('Specific humidity difference (ppm)')
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')

plt.subplot(2, 2, 4)
crossplot = plt.contourf(o3['lat'], o3['lev'], Qdiff_relate, levels=np.arange(-1, 1, 0.1),
                        cmap = 'bwr',extend='both')
plt.plot(o3['lat'],tpp*0.01,'-y',label = 'TROP_chemUCI')
plt.plot(o3['lat'],tpp_3d*0.01,'-m',label = 'TROP3D_chemUCI')
plt.plot(o3['lat'],tpp_refer*0.01,'-c',label = 'TROP_v2')
plt.title('Specific humidity relative difference')
plt.ylim(1000,1)
plt.yscale('log')
plt.xlabel('Lat')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(crossplot)
plt.legend(loc = 'upper left')
pylab.savefig(pathout+'Q_pres_lat_plot.png', dpi=300)
EOF

# Run diagnostics
command="python -u e3sm_pres_lat_plots.py"
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

