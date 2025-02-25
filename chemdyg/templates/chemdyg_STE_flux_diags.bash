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
ncfile_save="{{ ncfile_save }}"
if [[ "${ncfile_save}" == "true" ]]; then
   results_dir={{ output }}/post/atm/chemdyg_output
   mkdir -p ${results_dir}
fi

# Create temporary workdir
workdir=`mktemp -d tmp.${id}.XXXX`
cd ${workdir}

# Create local links to input climo files
#tsDir={{ output }}/post/atm/{{ grid }}/ts/monthly/{{ '%dyr' % (ypf) }}
mkdir -p ts

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

#mapDir="{{ diagnostics_base_path }}/observations/Atm/ChemDyg_inputs"
mapDir={{ mapping_file }}
ln -s ${mapDir} ./


# Run E3SM chem Diags
echo
echo ===== RUN E3SM CHEM DIAGS  =====
echo

# Prepare configuration file
cat > STE_flux_diags.py << EOF
#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from calendar import monthrange
import pandas as pd
import pylab 
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point
import os

path = './ts/'
pathout = './'

short_name = '${short}'
startyear = ${y1}
endyear = ${y2}
nyears = endyear - startyear + 1

filename = short_name+'.eam.h0.*.nc'
filenameh1 = short_name+'.${eamfile}.*.nc'

varname = ["O3"]

h0_in = xr.open_mfdataset(path+filename)
h1_in = xr.open_mfdataset(path+filenameh1)

timeperiod = len(h0_in['time'])
startdate = str(np.array(h0_in['time'].dt.year[0]))+'-01-01'

time_range_month = pd.date_range(startdate,  periods=timeperiod, freq='ME')
h0_in['time'] = time_range_month
h1_in['time'] = time_range_month
ncol = h0_in['ncol']
rearth = 6.37122e6 # Earth radius: m

area_rad = h0_in['area'][0]         # radian (ncol)
area = area_rad * rearth * rearth  # m2
lat = h0_in['lat'][0]
NH = area.where(lat >= 0)
SH = area.where(lat < 0)

time = h0_in['time']
year = np.array(time.dt.year)
month = np.array(time.dt.month)

STE_2D = h1_in['O3_2DMSD_trop'][0:-1,:].copy()
STE_time = []
STE_NH_time = []
STE_SH_time = []

for var in range(len(varname)):
  
    for i in range(len(time)):

        dt = monthrange(2001,month[i])[1]*3600*24 #considered different days in different months

        MSD = h1_in[varname[var]+'_2DMSD_trop'][i,:] #kg/m2

        TDB = h0_in[varname[var]+'_2DTDB_trop'][i,:]
        TDD = h0_in[varname[var]+'_2DTDD_trop'][i,:]
        TDE = h0_in[varname[var]+'_2DTDE_trop'][i,:]
        TDI = h0_in[varname[var]+'_2DTDI_trop'][i,:]
        TDA = h0_in[varname[var]+'_2DTDA_trop'][i,:]
        TDL = h0_in[varname[var]+'_2DTDL_trop'][i,:]
        TDN = h0_in[varname[var]+'_2DTDN_trop'][i,:]
        TDO = h0_in[varname[var]+'_2DTDO_trop'][i,:]
        TDS = h0_in[varname[var]+'_2DTDS_trop'][i,:]
        TDU = h0_in[varname[var]+'_2DTDU_trop'][i,:]

        total_td = (TDO+TDE+TDI+TDA+TDL+TDN+TDU+TDB+TDS+TDD)

        MSD_total = (MSD*area).sum()
        td_temp = total_td*dt
        TTD_total = (td_temp*area).sum()

        if i == 0:
            STE = 'nan'
        else:
            temp = MSD_old+td_temp
            STE = ((MSD-temp)*area).sum()
            STE_NH = ((MSD-temp)*NH).sum()
            STE_SH = ((MSD-temp)*SH).sum()
            STE_time.append(STE)
            STE_NH_time.append(STE_NH)
            STE_SH_time.append(STE_SH)
            STE_2D[i-1,:] = MSD-temp
        MSD_old = MSD

y = xr.concat(STE_time,dim="time")
yN = xr.concat(STE_NH_time,dim="time")
yS = xr.concat(STE_SH_time,dim="time")

y_ax = np.array(y[11::])
yN_ax = np.array(yN[11::])
yS_ax = np.array(yS[11::])
y_ann = np.zeros(12)
yN_ann = np.zeros(12)
yS_ann = np.zeros(12)
y_std = np.zeros(12)
yN_std = np.zeros(12)
yS_std = np.zeros(12)

if nyears > 1:
    yann = y_ax.reshape((nyears-1,12))
    yNann = yN_ax.reshape((nyears-1,12))
    ySann = yS_ax.reshape((nyears-1,12))
    for i in range(12):
        if i == 0:
           y_ann[i] = np.mean(yann[:,i])*1.E-9*12
           y_std[i] = np.std(yann[:,i])*1.E-9*12
           yN_ann[i] = np.mean(yNann[:,i])*1.E-9*12
           yN_std[i] = np.std(yNann[:,i])*1.E-9*12
           yS_ann[i] = ySann[:,i].mean()*1.E-9*12
           yS_std[i] = ySann[:,i].std()*1.E-9*12
        else:
           temp = np.append(yann[:,i],np.array(y[i-1]))
           temp1 = np.append(yNann[:,i],np.array(yN[i-1]))
           temp2 = np.append(ySann[:,i],np.array(yS[i-1]))
           y_ann[i] = np.mean(temp)*1.E-9*12
           y_std[i] = np.std(temp)*1.E-9*12
           yN_ann[i] = np.mean(temp1)*1.E-9*12
           yN_std[i] = np.std(temp1)*1.E-9*12
           yS_ann[i] = np.mean(temp2)*1.E-9*12
           yS_std[i] = np.std(temp2)*1.E-9*12
else:
     y_ann[0] = 'nan'
     y_ann[1:12] = y[0:11]*1.E-9*12
     yN_ann[0] = 'nan'
     yN_ann[1:12] = y[0:11]*1.E-9*12
     yS_ann[0] = 'nan'
     yS_ann[1:12] = y[0:11]*1.E-9*12

STE_month = xr.DataArray(np.array(STE_2D), coords=[month[1::],ncol], dims=["month","ncol"])
STE_month_mean= np.zeros((12,len(ncol)))
for i in range(12):
    if len(time) == 12:
        STE_month_mean[0,:] = 0
        STE_month_mean[1::,:] = 1.e3*12*STE_month
        break
    else:
        m = i+1
        STE_month_mean[i,:] = 1.e3*12*STE_month.sel(month=m).mean(axis=0)

# ----- writing ncfile -----
y_xr = xr.DataArray(y, coords=[time_range_month[1::]], dims=["time"])
yN_xr = xr.DataArray(yN, coords=[time_range_month[1::]], dims=["time"])
yS_xr = xr.DataArray(yS, coords=[time_range_month[1::]], dims=["time"])
y_xr = y_xr.assign_attrs(_FillValue = 1.0e36,units="kg/month", description="O3 flux")
yN_xr = yN_xr.assign_attrs(_FillValue = 1.0e36, units="kg/month", description="NH O3 flux")
yS_xr = yS_xr.assign_attrs(_FillValue = 1.0e36, units="kg/month", description="SH O3 flux")
ds1 = y_xr.to_dataset(name='y')
ds2 = yN_xr.to_dataset(name='yN')
ds3 = yS_xr.to_dataset(name='yS')
y_ann_xr = xr.DataArray(y_ann, coords=[np.arange(1,13)], dims=["month"])
y_std_xr = xr.DataArray(y_std, coords=[np.arange(1,13)], dims=["month"])
y_ann_xr = y_ann_xr.assign_attrs(_FillValue = 1.0e36, units="Tg/year", description='O3 flux mean')
y_std_xr = y_std_xr.assign_attrs(_FillValue = 1.0e36, units="Tg/year", description='O3 flux standard deviation')
yN_ann_xr = xr.DataArray(yN_ann, coords=[np.arange(1,13)], dims=["month"])
yN_std_xr = xr.DataArray(yN_std, coords=[np.arange(1,13)], dims=["month"])
yN_ann_xr = yN_ann_xr.assign_attrs(_FillValue = 1.0e36, units="Tg/year", description='NH O3 flux mean')
yN_std_xr = yN_std_xr.assign_attrs(_FillValue = 1.0e36, units="Tg/year", description='NH O3 flux standard deviation')
yS_ann_xr = xr.DataArray(yS_ann, coords=[np.arange(1,13)], dims=["month"])
yS_std_xr = xr.DataArray(yS_std, coords=[np.arange(1,13)], dims=["month"])
yS_ann_xr = yS_ann_xr.assign_attrs(_FillValue = 1.0e36, units="Tg/year", description='SH O3 flux mean')
yS_std_xr = yS_std_xr.assign_attrs(_FillValue = 1.0e36, units="Tg/year", description='SH O3 flux standard deviation')
ds4 = y_ann_xr.to_dataset(name='y_ann')
ds5 = yN_ann_xr.to_dataset(name='yN_ann')
ds6 = yS_ann_xr.to_dataset(name='yS_ann')
ds7 = y_std_xr.to_dataset(name='y_std')
ds8 = yN_std_xr.to_dataset(name='yN_std')
ds9 = yS_std_xr.to_dataset(name='yS_std')
ds = xr.merge([ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9])
ds.to_netcdf(pathout+'E3SM_O3_STE_${y1}-${y2}.nc')

STE_xr = xr.DataArray(STE_month_mean, name='STE',coords=[np.arange(1,13),ncol], dims=["month","ncol"],
        attrs=dict(units="g/m2/year", description='Ozone STE flux', _FillValue = 1.0e36))
STE_xr.to_netcdf(pathout+'E3SM_O3_STE_month_mean.nc')

STE_time_xr = xr.DataArray(np.array(STE_2D), name='STE_2D',coords=[time_range_month[1::],ncol], dims=["time","ncol"],attrs=dict(units="kg/m2/month", description='Ozone STE flux', _FillValue = 1.0e36))
STE_time_xr.to_netcdf(pathout+'E3SM_O3_STE2D_${y1}-${y2}.nc')

# ----- time series plot ----- 
y_mean = 1.E-9*12*np.array(y.mean())
yN_mean = 1.E-9*12*np.array(yN.mean())
yS_mean = 1.E-9*12*np.array(yS.mean())

fig = plt.figure(figsize=(10,5))
plt.plot(time_range_month[1::],y*1.E-9*12)
plt.plot(time_range_month[1::],yN*1.E-9*12)
plt.plot(time_range_month[1::],yS*1.E-9*12)

plt.title('O3 STE flux (Tg/year)')
plt.xlabel('Time')
#plt.ylabel('Tg')
line1 = 'Global mean:'+str(np.round(y_mean,2))
line2 = 'NH mean:'+str(np.round(yN_mean,2))
line3 = 'SH mean:'+str(np.round(yS_mean,2))
plt.legend( [line1,line2,line3])
pylab.savefig(pathout+'O3_STE_flux.png', dpi=300)

# ----- annual plot -----
fig = plt.figure(figsize=(10,5))
month_txt = np.arange(1,13)
plt.plot(month_txt,y_ann, linewidth = 1)
plt.fill_between(month_txt,y_ann+y_std,y_ann-y_std, alpha=.5, linewidth=0)
plt.plot(month_txt,yN_ann, linewidth = 1)
plt.fill_between(month_txt,yN_ann+yN_std,yN_ann-yN_std, alpha=.5, linewidth=0)
plt.plot(month_txt,yS_ann, linewidth = 1)
plt.fill_between(month_txt,yS_ann+yS_std,yS_ann-yS_std, alpha=.5, linewidth=0)

plt.title('O3 STE flux (Tg/year)')
#plt.legend(loc='upper left')
#plt.xlim(1,12,1)
#plt.axes().set_yscale("log")
plt.xlabel('Month')
#plt.ylabel('Tg')
line1 = 'Global mean:'+str(np.round(y_mean,2))
line2 = 'NH mean:'+str(np.round(yN_mean,2))
line3 = 'SH mean:'+str(np.round(yS_mean,2))
plt.legend( [line1,'Global std',line2,'NH std',line3,'SH std'])
pylab.savefig(pathout+'O3_STE_flux_monthly.png', dpi=300)

#----- 12 month lat-lon plot---------
cmd = "ncks --map=${mapDir} E3SM_O3_STE_month_mean.nc E3SM_O3_STE_month_latlon.nc"

remap = os.system(cmd)  
print('remapping')

## plot lat-lon STE flux
filename = './E3SM_O3_STE_month_latlon.nc'
file_in = xr.open_dataset(filename)

STE = file_in['STE'] #g/m2/year
nmonth = ['January','February','March','April','May','June','July','August','September','October','November','December']

fig, axs = plt.subplots(nrows=3,ncols=4,
                        subplot_kw={'projection': ccrs.PlateCarree()},
                        figsize=(18,12))
axs=axs.flatten()

#Loop over all of the models
for i in range(12):
    cs=axs[i].contourf(file_in['lon'], file_in['lat'], STE[i,:,:], levels=np.arange(-6.5, 7.5, 1),
        cmap="bwr",extend='both',transform=ccrs.PlateCarree())
    axs[i].set_title(nmonth[i] +' mean: '+str(np.round(y_ann[i],2))+' Tg/yr')
    axs[i].coastlines()
    # Longitude labels
    axs[i].set_xticks(np.arange(-180,181,60), crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    axs[i].xaxis.set_major_formatter(lon_formatter)

    # Latitude labels
    axs[i].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    axs[i].yaxis.set_major_formatter(lat_formatter)
        # Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.25, top=0.9, left=0.05, right=0.95,
                    wspace=0.15, hspace=0.3)

# Add a colorbar axis at the bottom of the graph
cbar_ax = fig.add_axes([0.2, 0.15, 0.6, 0.02])

# Draw the colorbar
cbar=fig.colorbar(cs, ticks=np.arange(-6.5, 7.5, 1), cax=cbar_ax, orientation='horizontal')
cbar.ax.tick_params(labelsize=14)
#fig.tight_layout()
# Add a big title at the top
plt.suptitle('Ozone STE flux (g/m2/year) ${y1}-${y2} ', fontsize=24)

pylab.savefig('STE_lat_lon_12month.png', dpi=600)

EOF

# Run diagnostics
command="python -u STE_flux_diags.py"
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
mv *.png ${f}
rm map_*.nc 
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

