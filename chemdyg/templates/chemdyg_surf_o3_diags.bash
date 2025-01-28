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
# diagnostics_base_path is set by zppy using the mache package
obsDir="{{ diagnostics_base_path }}/observations/Atm/ChemDyg_inputs"
ncfile_save="{{ ncfile_save }}"
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
tsDir1={{ output }}/post/atm/{{ grid1 }}/ts/hourly/{{ '%dyr' % (ypf) }}
tsDir2={{ output }}/post/atm/{{ grid2 }}/ts/hourly/{{ '%dyr' % (ypf) }}
#tsDir={{ output }}/post/atm/{{ grid }}
mkdir -p ts
ln -s ${obsDir}/surfO3/mda8.surfO3.*.nc ./ts
#cd ts
ln -s ${tsDir1}/O3_SRF_${Y1}*.nc ./ts/O3_SRF_{{ grid1 }}.nc
ln -s ${tsDir2}/O3_SRF_${Y1}*.nc ./ts/O3_SRF_{{ grid2 }}.nc
#ln -s ${cmipDir}/*TCO.nc ./ts
#cd ..
# Create symbolic links to input files
#input={{ input }}/{{ input_subdir }}
#for (( year=${y1}; year<=${y2}; year++ ))
#do
#  YYYY=`printf "%04d" ${year}`
#  for file in ${input}/${case}.{{ input_filesh3 }}.${YYYY}-*.nc
#  do
#    ln -s ${file} ./ts
#  done
  #for file in ${input}/${case}.{{ input_files2 }}.${YYYY}-*.nc
  #do
  #  ln -s ${file} .
  #done
#done

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
echo ===== RUN E3SM CHEM DIAGS  =====
echo

# Prepare configuration file
cat > surf_O3_diags.py << EOF
#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import xarray as xr
import pylab
import math

pathin = './ts/'
pathout = './'

short_name = '${short}'

EUdatain = xr.open_dataset(pathin+'O3_SRF_MDA8EU1.0x1.0_nco.nc')
USdatain = xr.open_dataset(pathin+'O3_SRF_MDA8US1.0x1.0_nco.nc')
EUmaskin = xr.open_dataset(pathin+'mda8.surfO3.EU.2000.2009.nc')
USmaskin = xr.open_dataset(pathin+'mda8.surfO3.US.2000.2009.nc')

# ## UTC to local time
rearth = 6.37122e6 # Earth radius: m
# US region
lat_US = USdatain['lat']
lon_US = USdatain['lon']
nlat_US = len(lat_US)
nlon_US = len(lon_US)

time_US = USdatain['time']
ntime_US = len(time_US)

O3_US = USdatain['O3_SRF']

arearad_US = USdatain['area']# radian
area_US = arearad_US * rearth * rearth  # m2

# EU region
lat_EU = EUdatain['lat']
lon_EU = EUdatain['lon']
nlat_EU = len(lat_EU)
nlon_EU = len(lon_EU)

time_EU = EUdatain['time']
ntime_EU = len(time_EU)

O3_EU = EUdatain['O3_SRF']

arearad_EU = EUdatain['area']# radian
area_EU = arearad_EU * rearth * rearth  # m2

#Time zone 0 is between 360-7.5 and 360+7.5
O3local_US = np.zeros((ntime_US,nlat_US,nlon_US))
O3local_EU = np.zeros((ntime_EU,nlat_EU,nlon_EU))
#US longitudes: 234 to 293 degree
#check time zone boundaries
for i in range(0,9):
    tmp_US = 360. - 7.5 - 15.*i
    #print("Time zone -",i,": ",tmp_US," to ",tmp_US+15.0)

#EU longitudes: -12 to 33 degree
#check time zone boundaries
for i in range(-2,4):
    tmp_EU = - 7.5 + 15.*i
    #print("Time zone -",i,": ",tmp_EU," to ",tmp_EU+15.0)

for i in range(nlon_US):
    #calculate the local time zone for each longitude
    Tzone_US = math.floor((lon_US[i] - (360.-7.5))/15.)
    #print("i=",i,", lon(i)=",lon_US[i].values,", Tzone=",Tzone_US)
    ts = 0 - Tzone_US
    O3local_US[0:Tzone_US,:,i] = O3_US[ts::,:,i].copy()

for i in range(nlon_EU):
    #calculate the local time zone for each longitude
    Tzone_EU = math.floor((lon_EU[i] - (-7.5))/15.)
    #print("i=",i,", lon(i)=",lon_EU[i].values,", Tzone=", Tzone_EU)
    if Tzone_EU < 0:
        ts = 0 - Tzone_EU
        O3local_EU[0:Tzone_EU,:,i] = O3_EU[ts::,:,i].copy()
    else:
        te = ntime_EU-Tzone_EU
        O3local_EU[Tzone_EU::,:,i] = O3_EU[0:te,:,i].copy()

O3local_EU[O3local_EU == 0.] = 'nan'
O3local_US[O3local_US == 0.] = 'nan'

# ## calculate diurnal cycle
O3EU_xr = xr.DataArray(O3local_EU, coords=[time_EU,lat_EU,lon_EU], dims=["time","lat","lon"])
O3US_xr = xr.DataArray(O3local_US, coords=[time_US,lat_US,lon_US], dims=["time","lat","lon"])
O3EU_JJA = O3EU_xr.sel(time=O3EU_xr.time.dt.month.isin([6, 7, 8]))
O3EU_DJF = O3EU_xr.sel(time=O3EU_xr.time.dt.month.isin([12, 1, 2]))
O3US_JJA = O3US_xr.sel(time=O3US_xr.time.dt.month.isin([6, 7, 8]))
O3US_DJF = O3US_xr.sel(time=O3US_xr.time.dt.month.isin([12, 1, 2]))

EUmask = np.transpose(EUmaskin['MDA8_SurfO3'][0,:,:].values, (1,0))
USmask = np.transpose(USmaskin['MDA8_SurfO3'][0,:,:].values, (1,0))

mask_EU = np.ma.masked_invalid(EUmask).mask
mask_US = np.ma.masked_invalid(USmask).mask

for i in range(len(O3EU_JJA)):
    tmp = O3EU_JJA[i,:,:].values
    tmp[mask_EU] = 'nan'
    O3EU_JJA[i,:,:] = tmp
    tmp2 = O3US_JJA[i,:,:].values
    tmp2[mask_US] = 'nan'
    O3US_JJA[i,:,:] = tmp2

for i in range(len(O3EU_DJF)):
    tmp = O3EU_DJF[i,:,:].values
    tmp[mask_EU] = 'nan'
    O3EU_DJF[i,:,:] = tmp
    tmp2 = O3US_DJF[i,:,:].values
    tmp2[mask_US] = 'nan'
    O3US_DJF[i,:,:] = tmp2

tmp = area_US.values
tmp[mask_US] = 'nan'
area_US[:,:] = tmp
tmp = area_EU.values
tmp[mask_EU] = 'nan'
area_EU[:,:] = tmp

# US subregion: Western (WNA) and eastern (ENA) North America is split by 264°E
lon1 = 234
lon2 = 263
lon3 = 293
O3US_JJA_WNA = O3US_JJA.sel(lon=slice(lon1,lon2))
O3US_JJA_ENA = O3US_JJA.sel(lon=slice(lon2+1,lon3))
O3US_DJF_WNA = O3US_DJF.sel(lon=slice(lon1,lon2))
O3US_DJF_ENA = O3US_DJF.sel(lon=slice(lon2+1,lon3))

area_WNA = area_US.sel(lon=slice(lon1,lon2))
area_ENA = area_US.sel(lon=slice(lon2+1,lon3))

# EU subregion: Southern (SEU) and northern (NEU) Europe is split by 53°N. (Schnell et al., ACP 2015)
lat1 = 35
lat2 = 52
lat3 = 71
O3EU_JJA_SEU = O3EU_JJA.sel(lat=slice(lat1,lat2))
O3EU_JJA_NEU = O3EU_JJA.sel(lat=slice(lat2+1,lat3))
O3EU_DJF_SEU = O3EU_DJF.sel(lat=slice(lat1,lat2))
O3EU_DJF_NEU = O3EU_DJF.sel(lat=slice(lat2+1,lat3))

area_SEU = area_EU.sel(lat=slice(lat1,lat2))
area_NEU = area_EU.sel(lat=slice(lat2+1,lat3))

# calculate area-weighted averages over the region
O3_ENA_JJA = np.zeros(len(O3EU_JJA))
O3_WNA_JJA = np.zeros(len(O3EU_JJA))
O3_NEU_JJA = np.zeros(len(O3EU_JJA))
O3_SEU_JJA = np.zeros(len(O3EU_JJA))
for i in range(len(O3EU_JJA)):
    O3_ENA_JJA[i] = (O3US_JJA_ENA[i,:,:]*area_ENA[:,:]).sum()/area_ENA[:,:].sum()
    O3_WNA_JJA[i] = (O3US_JJA_WNA[i,:,:]*area_WNA[:,:]).sum()/area_WNA[:,:].sum()
    O3_NEU_JJA[i] = (O3EU_JJA_NEU[i,:,:]*area_NEU[:,:]).sum()/area_NEU[:,:].sum()
    O3_SEU_JJA[i] = (O3EU_JJA_SEU[i,:,:]*area_SEU[:,:]).sum()/area_SEU[:,:].sum()

# calculate area-weighted averages over the region
O3_ENA_DJF = np.zeros(len(O3EU_DJF))
O3_WNA_DJF = np.zeros(len(O3EU_DJF))
O3_NEU_DJF = np.zeros(len(O3EU_DJF))
O3_SEU_DJF = np.zeros(len(O3EU_DJF))
for i in range(len(O3EU_DJF)):
    O3_ENA_DJF[i] = (O3US_DJF_ENA[i,:,:]*area_ENA[:,:]).sum()/area_ENA[:,:].sum()
    O3_WNA_DJF[i] = (O3US_DJF_WNA[i,:,:]*area_WNA[:,:]).sum()/area_WNA[:,:].sum()
    O3_NEU_DJF[i] = (O3EU_DJF_NEU[i,:,:]*area_NEU[:,:]).sum()/area_NEU[:,:].sum()
    O3_SEU_DJF[i] = (O3EU_DJF_SEU[i,:,:]*area_SEU[:,:]).sum()/area_SEU[:,:].sum()

ndays = int(len(O3EU_JJA)/24)
n_d = int(ndays*24)
O3_ENA_JJA_24h = 1.e9*O3_ENA_JJA[0:n_d].reshape((ndays,24)).mean(axis=0)
O3_WNA_JJA_24h = 1.e9*O3_WNA_JJA[0:n_d].reshape((ndays,24)).mean(axis=0)
O3_NEU_JJA_24h = 1.e9*O3_NEU_JJA[0:n_d].reshape((ndays,24)).mean(axis=0)
O3_SEU_JJA_24h = 1.e9*O3_SEU_JJA[0:n_d].reshape((ndays,24)).mean(axis=0)

ndays = int(len(O3EU_DJF)/24)
n_d = int(ndays*24)
O3_ENA_DJF_24h = 1.e9*O3_ENA_DJF[0:n_d].reshape((ndays,24)).mean(axis=0)
O3_WNA_DJF_24h = 1.e9*O3_WNA_DJF[0:n_d].reshape((ndays,24)).mean(axis=0)
O3_NEU_DJF_24h = 1.e9*O3_NEU_DJF[0:n_d].reshape((ndays,24)).mean(axis=0)
O3_SEU_DJF_24h = 1.e9*O3_SEU_DJF[0:n_d].reshape((ndays,24)).mean(axis=0)

# ----- writing ncfile -----
O3_ENA_JJA_24h_xr = xr.DataArray(O3_ENA_JJA_24h, name= 'O3_ENA_JJA_24h', coords=[np.arange(0,24)], dims=["hour"], 
                       attrs=dict(units="ppb", description="JJA O3 abundance in ENA") )
O3_WNA_JJA_24h_xr = xr.DataArray(O3_WNA_JJA_24h, name= 'O3_WNA_JJA_24h', coords=[np.arange(0,24)], dims=["hour"],  
		       attrs=dict(units="ppb", description="JJA O3 abundance in WNA") )
O3_NEU_JJA_24h_xr = xr.DataArray(O3_NEU_JJA_24h, name= 'O3_NEU_JJA_24h', coords=[np.arange(0,24)], dims=["hour"], 
		       attrs=dict(units="ppb", description="JJA O3 abundance in NEU") )
O3_SEU_JJA_24h_xr = xr.DataArray(O3_SEU_JJA_24h, name= 'O3_SEU_JJA_24h', coords=[np.arange(0,24)], dims=["hour"], 
		       attrs=dict(units="ppb", description="JJA O3 abundance in SEU") )
O3_ENA_DJF_24h_xr = xr.DataArray(O3_ENA_DJF_24h, name= 'O3_ENA_DJF_24h', coords=[np.arange(0,24)], dims=["hour"], 
		       attrs=dict(units="ppb", description="DJF O3 abundance in ENA") )
O3_WNA_DJF_24h_xr = xr.DataArray(O3_WNA_DJF_24h, name= 'O3_WNA_DJF_24h', coords=[np.arange(0,24)], dims=["hour"], 
		       attrs=dict(units="ppb", description="DJF O3 abundance in WNA") )
O3_NEU_DJF_24h_xr = xr.DataArray(O3_NEU_DJF_24h, name= 'O3_NEU_DJF_24h', coords=[np.arange(0,24)], dims=["hour"], 
		       attrs=dict(units="ppb", description="DJF O3 abundance in NEU") )
O3_SEU_DJF_24h_xr = xr.DataArray(O3_SEU_DJF_24h, name= 'O3_SEU_DJF_24h', coords=[np.arange(0,24)], dims=["hour"], 
                       attrs=dict(units="ppb", description="DJF O3 abundance in SEU") )
ds1 = O3_ENA_JJA_24h_xr.to_dataset()
ds2 = O3_WNA_JJA_24h_xr.to_dataset()
ds3 = O3_NEU_JJA_24h_xr.to_dataset()
ds4 = O3_SEU_JJA_24h_xr.to_dataset()
ds5 = O3_ENA_DJF_24h_xr.to_dataset()
ds6 = O3_WNA_DJF_24h_xr.to_dataset()
ds7 = O3_NEU_DJF_24h_xr.to_dataset()
ds8 = O3_SEU_DJF_24h_xr.to_dataset()

# ## read observations
obs_WNA_JJA = [30.5009, 29.4010, 28.0499, 26.9680, 25.3009, 23.5642, 24.1381, 28.3088, 33.9598, 39.3651,
               43.9309, 47.3310, 49.5168, 50.7556, 51.2856, 51.2519, 50.6595, 49.2926, 46.2216, 41.3574,
               37.2174, 34.2102, 32.9063, 31.5910]
model_WNA_JJA = [
   [38.3237 , 46.6769 , 33.3072 , 36.3390 , 34.2388 , 36.3785 , 40.5814 , 48.6485 , 43.9095 ],
   [35.4660 , 45.0604 , 31.6952 , 34.3149 , 32.8997 , 35.0905 , 39.9923 , 47.6264 , 41.9392 ],
   [32.9176 , 43.4621 , 30.1940 , 32.9813 , 31.6681 , 33.9253 , 39.4184 , 46.6544 , 40.1781 ],
   [30.5005 , 41.5636 , 28.7727 , 31.8824 , 30.5436 , 32.8657 , 38.7692 , 45.7193 , 38.5836 ],
   [28.2665 , 38.5992 , 27.4133 , 30.5863 , 29.9481 , 32.1295 , 37.8178 , 44.7974 , 38.7750 ],
   [26.2499 , 35.6021 , 26.3357 , 29.8903 , 29.7641 , 32.9888 , 37.3223 , 43.4690 , 43.9065 ],
   [25.3259 , 35.0645 , 26.4127 , 31.9010 , 29.9329 , 34.4405 , 38.2364 , 41.5245 , 53.1372 ],
   [26.7237 , 38.0902 , 27.4703 , 36.1526 , 32.4335 , 38.2055 , 39.7206 , 41.1645 , 63.7333 ],
   [32.6813 , 43.3365 , 29.6613 , 40.5435 , 37.0242 , 43.1281 , 41.5703 , 45.7494 , 71.6495 ],
   [37.7971 , 48.4279 , 33.6659 , 44.6653 , 41.8538 , 46.1664 , 43.8257 , 52.2338 , 76.4492 ],
   [46.9608 , 51.8775 , 38.7573 , 47.9082 , 45.3092 , 47.9730 , 46.0343 , 56.1869 , 77.3953 ],
   [52.2351 , 53.9511 , 42.2944 , 50.1719 , 47.5056 , 48.9980 , 47.8039 , 58.5587 , 76.1436 ],
   [53.4974 , 55.1033 , 44.3418 , 51.7482 , 48.8586 , 49.5905 , 49.0500 , 60.0239 , 74.1259 ],
   [54.3193 , 55.6967 , 45.4762 , 52.8330 , 49.6813 , 49.9542 , 49.8348 , 60.8858 , 71.8295 ],
   [55.1993 , 56.0431 , 46.1103 , 53.5286 , 50.1650 , 50.1129 , 50.1570 , 61.3386 , 69.4160 ],
   [55.5535 , 56.1782 , 46.4396 , 53.9071 , 50.3337 , 50.3484 , 49.9777 , 61.4782 , 66.8310 ],
   [55.1026 , 55.9553 , 46.5237 , 53.9179 , 50.0054 , 50.2433 , 49.2330 , 60.2584 , 64.0391 ],
   [54.7270 , 55.2110 , 46.3001 , 52.6201 , 48.4556 , 49.3351 , 47.9109 , 59.2620 , 60.9097 ],
   [54.4865 , 53.9889 , 45.5469 , 49.0965 , 45.6867 , 47.5959 , 46.1285 , 59.4183 , 58.1459 ],
   [52.8066 , 52.7962 , 43.9853 , 45.9742 , 43.1892 , 45.2548 , 44.4524 , 56.5113 , 55.8285 ],
   [50.7549 , 51.9591 , 41.6280 , 43.8229 , 40.9599 , 42.9891 , 43.3776 , 53.5707 , 52.8816 ],
   [47.7057 , 51.0427 , 39.2061 , 42.5959 , 38.9826 , 41.0359 , 42.5583 , 52.1389 , 50.7435 ],
   [44.7201 , 49.7600 , 36.9913 , 40.7815 , 37.2328 , 39.2998 , 41.8349 , 50.9485 , 48.3976 ],
   [41.4339 , 48.2594 , 35.0207 , 38.6263 , 35.6598 , 37.7502 , 41.1780 , 49.7985 , 45.9767 ]]

obs_ENA_JJA = [23.1592, 22.1130, 20.7995, 19.8345, 18.4364, 17.2065, 18.0162, 22.4712, 28.5121, 34.3992,
               39.2628, 42.6399, 44.7002, 45.8349, 46.2733, 46.1549, 45.4153, 43.6565, 40.0355, 35.0141,
               30.6742, 27.7447, 25.7583, 24.2590]
model_ENA_JJA = [
   [48.8895 , 50.8043 , 28.3385 , 37.8247 , 29.7490 , 42.4351 , 50.4484 , 40.8915 , 46.3619 ],
   [46.6944 , 48.2353 , 26.5267 , 36.4186 , 28.3330 , 40.7719 , 49.9278 , 39.7626 , 43.5495 ],
   [44.4365 , 45.8051 , 24.8381 , 35.1164 , 27.0132 , 39.2339 , 49.3943 , 38.6428 , 40.9052 ],
   [42.0379 , 43.3292 , 23.2469 , 33.9348 , 25.8025 , 37.8102 , 48.7067 , 37.5302 , 38.7375 ],
   [39.1714 , 40.2385 , 21.7405 , 32.8863 , 25.0111 , 36.7614 , 47.1522 , 36.4171 , 39.3111 ],
   [36.2254 , 37.6483 , 20.7514 , 32.1033 , 24.9443 , 36.5831 , 45.2737 , 35.0323 , 48.0262 ],
   [34.1064 , 37.9502 , 21.7027 , 32.6430 , 24.6709 , 39.7320 , 45.1236 , 33.1200 , 63.3745 ],
   [33.8522 , 42.7065 , 23.7482 , 35.7091 , 25.7706 , 44.7980 , 46.5591 , 32.6801 , 80.2375 ],
   [36.2156 , 49.7151 , 26.1726 , 39.9301 , 28.7145 , 49.0839 , 48.9931 , 37.3894 , 95.3729 ],
   [47.1667 , 56.7697 , 29.2899 , 43.9923 , 33.1352 , 53.2481 , 51.9010 , 44.2167 , 99.8032 ],
   [52.8289 , 62.1164 , 33.4579 , 47.7633 , 37.6277 , 55.9573 , 54.5911 , 48.4616 , 99.0580 ],
   [55.7811 , 65.2774 , 37.0141 , 50.7368 , 41.3107 , 57.3174 , 56.7351 , 50.9411 , 97.4854 ],
   [58.7134 , 67.0403 , 39.5735 , 52.9906 , 44.0544 , 58.1812 , 58.3162 , 52.4757 , 93.4153 ],
   [61.5334 , 68.1118 , 41.3670 , 54.6536 , 45.8857 , 58.8289 , 59.3850 , 53.3922 , 89.0636 ],
   [64.0902 , 68.7105 , 42.5161 , 55.6860 , 46.9467 , 59.4893 , 59.9324 , 53.8869 , 85.4602 ],
   [66.0143 , 68.8763 , 43.1347 , 56.1466 , 47.2759 , 59.5150 , 59.8831 , 54.0421 , 81.7373 ],
   [66.4447 , 68.4353 , 43.2927 , 55.5991 , 46.5578 , 59.1011 , 59.1185 , 53.7700 , 77.7250 ],
   [64.3781 , 67.0596 , 42.9015 , 53.6976 , 44.2865 , 58.6813 , 57.5079 , 52.7311 , 73.2767 ],
   [60.6763 , 64.8216 , 41.6713 , 50.5860 , 41.5141 , 57.0172 , 55.1527 , 49.2939 , 68.4142 ],
   [59.2571 , 62.4000 , 39.5028 , 46.8842 , 38.9994 , 54.1518 , 53.1187 , 47.4939 , 64.1545 ],
   [57.5984 , 60.4663 , 37.0758 , 44.3140 , 36.7158 , 51.1110 , 52.3739 , 45.7012 , 61.0925 ],
   [55.1786 , 58.5433 , 34.6898 , 42.5148 , 34.6887 , 48.4859 , 51.9245 , 44.4039 , 56.8863 ],
   [52.9626 , 56.1338 , 32.3907 , 40.7864 , 32.8740 , 46.2348 , 51.4450 , 43.2494 , 52.8789 ],
   [50.8452 , 53.4756 , 30.2515 , 39.2372 , 31.2313 , 44.2427 , 50.9508 , 42.0906 , 49.4451 ]]

obs_SEU_JJA = [28.5354, 27.5896, 26.6149, 25.3877, 24.3031, 23.3113, 23.7851, 26.3779, 30.4694, 35.2119,
               39.5751, 43.1464, 45.6961, 47.2971, 48.0672, 48.2060, 47.6975, 46.4048, 44.0225, 40.4302,
               36.5313, 33.3716, 31.1859, 29.6564]
model_SEU_JJA = [
   [45.8509 , 49.8064 , 35.5119 , 44.3918 , 30.8168 , 42.1552 , 52.8522 , 45.3775 , 50.9056 ],
   [43.4181 , 47.9464 , 33.8779 , 42.8349 , 29.6878 , 41.0210 , 52.2460 , 43.8402 , 48.7774 ],
   [41.3361 , 46.0617 , 32.3541 , 41.6208 , 28.6572 , 39.9690 , 51.6524 , 44.0605 , 46.7820 ],
   [38.9971 , 43.7289 , 30.9169 , 40.5988 , 27.7854 , 39.0258 , 50.8805 , 43.0971 , 45.3570 ],
   [36.5378 , 40.8346 , 29.6072 , 39.9300 , 27.2985 , 38.5717 , 49.9221 , 42.0253 , 46.3516 ],
   [34.5306 , 38.3432 , 28.9848 , 40.2761 , 27.1049 , 38.4303 , 49.2494 , 40.4805 , 51.5448 ],
   [33.9239 , 37.8515 , 29.7125 , 42.6135 , 27.4929 , 39.3362 , 49.5402 , 38.7669 , 60.2597 ],
   [35.7116 , 40.6475 , 31.4116 , 46.0215 , 29.2570 , 42.8379 , 50.7018 , 39.0679 , 70.0997 ],
   [39.5601 , 45.7889 , 33.8633 , 49.8705 , 32.2381 , 45.8683 , 52.7102 , 43.1166 , 78.8727 ],
   [49.1002 , 50.8744 , 37.4024 , 53.1438 , 35.4948 , 48.2512 , 55.3462 , 47.6808 , 84.5826 ],
   [58.6189 , 54.9370 , 41.0160 , 55.8905 , 38.1773 , 50.1808 , 58.0158 , 50.7346 , 85.1817 ],
   [61.4119 , 57.7052 , 43.7816 , 58.0314 , 40.1618 , 51.4548 , 60.2695 , 52.7666 , 84.3519 ],
   [63.3153 , 59.3662 , 45.6775 , 59.6462 , 41.6055 , 52.3877 , 61.9435 , 54.1858 , 82.7721 ],
   [64.5463 , 60.3795 , 46.8708 , 60.8194 , 42.6484 , 53.0331 , 63.0285 , 55.1570 , 80.4048 ],
   [65.6956 , 60.9430 , 47.5725 , 61.5044 , 43.3248 , 53.5796 , 63.5151 , 55.7726 , 78.2270 ],
   [66.4642 , 61.1004 , 47.9367 , 61.8182 , 43.6067 , 53.8625 , 63.3748 , 56.0587 , 75.9783 ],
   [66.8266 , 60.8157 , 48.0263 , 61.3938 , 43.3100 , 53.4102 , 62.5738 , 56.0014 , 73.4653 ],
   [66.1517 , 59.9688 , 47.7944 , 59.6405 , 41.9692 , 52.5369 , 61.1060 , 55.5710 , 70.6770 ],
   [64.9405 , 58.6604 , 47.0430 , 57.0468 , 40.1062 , 51.3998 , 59.0939 , 54.5094 , 67.6858 ],
   [62.6415 , 57.1977 , 45.5328 , 53.6973 , 38.2669 , 49.6792 , 57.1486 , 52.3196 , 64.2201 ],
   [59.0751 , 55.9545 , 43.4405 , 50.9968 , 36.4945 , 47.8134 , 55.8932 , 49.9532 , 61.7322 ],
   [55.6176 , 54.8138 , 41.2388 , 49.3174 , 34.8527 , 46.1281 , 55.0266 , 48.7264 , 59.0353 ],
   [52.1152 , 53.3611 , 39.1515 , 47.7282 , 33.3502 , 44.6482 , 54.2570 , 47.8453 , 55.9023 ],
   [48.8085 , 51.6215 , 37.2250 , 46.0139 , 31.9823 , 43.3129 , 53.5452 , 46.8698 , 53.2624 ]]

obs_NEU_JJA = [23.9228, 23.1932, 22.6196, 22.1234, 21.8984, 22.1519, 23.2769, 25.2495, 27.5878, 29.8173,
               31.6843, 33.1849, 34.2819, 35.0222, 35.4672, 35.5625, 35.3184, 34.6435, 33.3237, 31.3428,
               29.0870, 27.0676, 25.5739, 24.586]
model_NEU_JJA = [
   [43.9291 , 38.4176 , 31.2814 , 38.8202 , 27.3655 , 34.5674 , 39.8432 , 30.0276 , 46.6842 ],
   [42.5314 , 37.3926 , 30.3225 , 38.1292 , 26.7544 , 33.8045 , 39.3851 , 29.1113 , 45.6737 ],
   [40.7333 , 36.1475 , 29.4172 , 37.4437 , 26.2118 , 33.1527 , 38.8276 , 29.4778 , 44.8565 ],
   [38.7311 , 34.6773 , 28.5900 , 36.7957 , 25.7625 , 32.5996 , 38.1683 , 28.8997 , 44.5593 ],
   [36.7872 , 33.3325 , 27.9644 , 36.4385 , 25.4006 , 32.5455 , 37.6038 , 28.1111 , 45.0750 ],
   [35.3732 , 32.7271 , 27.8138 , 36.9359 , 25.1531 , 32.5950 , 37.3371 , 27.2506 , 46.2769 ],
   [35.0384 , 33.0551 , 28.2189 , 38.1200 , 25.2773 , 33.1410 , 37.4249 , 26.9421 , 48.1974 ],
   [37.4308 , 34.2724 , 29.0923 , 39.7286 , 25.9288 , 34.7751 , 37.8787 , 28.0746 , 50.9247 ],
   [38.7595 , 36.2517 , 30.5190 , 41.6047 , 27.0711 , 36.2639 , 38.7691 , 30.2256 , 53.0194 ],
   [41.7194 , 38.2690 , 32.2827 , 43.1839 , 28.4992 , 37.4706 , 40.0073 , 32.1092 , 55.0903 ],
   [46.7225 , 40.1875 , 33.8881 , 44.5593 , 29.9359 , 38.6342 , 41.3656 , 33.4660 , 56.5556 ],
   [48.2588 , 41.9254 , 35.1711 , 45.7732 , 31.2012 , 39.6176 , 42.6350 , 34.5055 , 57.2717 ],
   [49.6109 , 43.2260 , 36.1802 , 46.8021 , 32.2353 , 40.4168 , 43.7127 , 35.3312 , 57.6838 ],
   [50.7119 , 44.1719 , 36.9583 , 47.5884 , 33.0282 , 40.9969 , 44.5253 , 35.9696 , 57.7889 ],
   [51.4144 , 44.7909 , 37.5161 , 48.0575 , 33.5376 , 41.4874 , 45.0305 , 36.4160 , 57.6090 ],
   [51.9118 , 45.0437 , 37.8588 , 48.2585 , 33.7225 , 41.7639 , 45.2077 , 36.6468 , 57.2478 ],
   [52.2343 , 44.9616 , 37.9830 , 47.9122 , 33.5035 , 41.4911 , 45.0577 , 36.6101 , 56.6113 ],
   [51.8936 , 44.4775 , 37.8731 , 46.6261 , 32.8314 , 40.9259 , 44.5392 , 36.1686 , 55.8507 ],
   [51.3694 , 43.6728 , 37.4693 , 45.0751 , 31.9842 , 40.2564 , 43.6519 , 35.1796 , 54.9572 ],
   [50.8400 , 42.7068 , 36.6983 , 43.5012 , 31.1896 , 39.3037 , 42.5741 , 33.8295 , 53.3520 ],
   [49.5833 , 41.7482 , 35.6219 , 42.0819 , 30.3585 , 38.2664 , 41.6779 , 32.5456 , 52.0973 ],
   [48.1928 , 40.9145 , 34.4348 , 40.9887 , 29.5201 , 37.2079 , 41.0624 , 31.6561 , 50.8349 ],
   [46.7454 , 40.1436 , 33.2925 , 40.1723 , 28.7190 , 36.2123 , 40.6416 , 31.1198 , 49.2191 ],
   [45.2930 , 39.3056 , 32.2140 , 39.4881 , 27.9765 , 35.3000 , 40.2959 , 30.6817 , 47.8586 ]]

# plot JJA diurnal figures
fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=4, sharex=True,
                                    figsize=(24, 5))
ax0.plot(obs_WNA_JJA,'ko-')
ax0.plot(model_WNA_JJA)
ax0.plot(O3_WNA_JJA_24h,color='blue',marker='o')
ax0.set_ylim(0,100)
ax0.set_ylabel('JJA O3 abundance (ppb)',fontsize='x-large')
ax0.set_xlabel('WNA local time (hr)',fontsize='x-large')

ax1.plot(obs_ENA_JJA,'ko-')
ax1.plot(model_ENA_JJA)
ax1.plot(O3_ENA_JJA_24h,color='blue',marker='o')
ax1.set_ylim(0,100)
#ax1.set_ylabel('JJA O3 abundance (ppb)',fontsize='x-large')
ax1.set_xlabel('ENA local time (hr)',fontsize='x-large')

ax2.plot(obs_SEU_JJA,'ko-')
ax2.plot(model_SEU_JJA)
ax2.plot(O3_SEU_JJA_24h,color='blue',marker='o')
ax2.set_ylim(0,100)
#ax2.set_ylabel('JJA O3 abundance (ppb)',fontsize='x-large')
ax2.set_xlabel('SEU local time (hr)',fontsize='x-large')

ax3.plot(obs_NEU_JJA,'ko-')
ax3.plot(model_NEU_JJA)
ax3.plot(O3_NEU_JJA_24h,color='blue',marker='o')
ax3.set_ylim(0,100)
#ax3.set_ylabel('JJA O3 abundance (ppb)',fontsize='x-large')
ax3.set_xlabel('NEU local time (hr)',fontsize='x-large')
ax3.legend(["OBS","MOCAGE","GFDL-AM3","CESM-CAM-SF","UM-CAM","CMAM","MIROC-CHEM","GISS-E2-R","GEOSCCM","UCI CTM","E3SM"],fontsize='x-large',loc='center left', bbox_to_anchor=(1, 0.5))

fig.suptitle(short_name,fontsize='x-large')
pylab.savefig(pathout+'surfO3_hour_cycle_JJA.png',dpi=300)

obs_WNA_DJF = [24.1957, 24.1874, 24.0635, 24.1840, 23.9255, 23.2425, 22.4406, 21.6357, 22.5392, 25.0628,
               27.7583, 30.1011, 31.9906, 33.2991, 33.9168, 33.5080, 31.7000, 28.6241, 26.3518, 25.4080,
               24.9390, 24.0683, 24.3993, 24.3896]
model_WNA_DJF = [
   [30.8906 , 42.2942 , 38.2967 , 30.5533 , 26.2589 , 30.1657 , 48.4584 , 38.8602 , 31.6468 ],
   [30.3224 , 41.8535 , 37.8419 , 30.3875 , 25.9975 , 29.7841 , 48.4007 , 38.4916 , 30.7267 ],
   [29.9192 , 41.4536 , 37.4474 , 30.2730 , 25.7748 , 29.4418 , 48.3597 , 38.1360 , 29.9185 ],
   [29.6233 , 41.0880 , 37.1043 , 30.1022 , 25.5822 , 29.1255 , 48.3335 , 37.7915 , 29.1365 ],
   [29.3642 , 40.7605 , 36.8038 , 29.9171 , 25.4168 , 28.8377 , 48.3153 , 37.4575 , 28.4306 ],
   [29.1161 , 40.4437 , 36.5414 , 29.7796 , 25.2694 , 28.5732 , 48.2982 , 37.1314 , 27.8423 ],
   [28.7370 , 39.9525 , 36.3115 , 29.7046 , 25.3228 , 28.4048 , 48.3690 , 36.8116 , 27.7363 ],
   [28.4672 , 39.4300 , 36.3538 , 29.9431 , 25.8069 , 28.6109 , 48.9071 , 36.4505 , 28.5393 ],
   [28.3117 , 39.2529 , 37.3982 , 31.2232 , 26.5090 , 29.6426 , 49.5428 , 36.3022 , 29.5874 ],
   [28.4780 , 39.7906 , 39.0171 , 33.0796 , 27.6811 , 31.1405 , 50.0404 , 36.7934 , 30.4602 ],
   [29.6541 , 40.9002 , 40.7650 , 34.5525 , 29.1957 , 32.7952 , 50.6239 , 38.2046 , 32.5765 ],
   [33.5074 , 42.3999 , 42.4240 , 35.7913 , 30.6433 , 34.4909 , 51.2866 , 40.4737 , 36.1352 ],
   [34.6429 , 43.8831 , 43.8894 , 36.8583 , 31.7482 , 35.8942 , 51.9113 , 42.3389 , 38.0136 ],
   [38.1907 , 45.0025 , 45.0147 , 37.5899 , 32.3593 , 36.7801 , 52.3123 , 43.4897 , 39.9807 ],
   [40.1866 , 45.6245 , 45.7037 , 37.8378 , 32.3131 , 37.0376 , 52.3577 , 44.0800 , 41.2406 ],
   [40.1239 , 45.7377 , 45.8673 , 37.5279 , 31.5376 , 36.9937 , 51.9574 , 44.1365 , 41.4336 ],
   [39.4511 , 45.5614 , 45.4155 , 36.4835 , 30.5013 , 36.3588 , 51.1493 , 42.8263 , 41.1922 ],
   [39.2938 , 45.4492 , 44.3243 , 34.7820 , 29.6631 , 35.1346 , 50.3087 , 41.3389 , 40.3328 ],
   [39.1879 , 45.3343 , 43.1264 , 33.6550 , 28.9488 , 33.9829 , 49.6789 , 41.2661 , 39.6803 ],
   [37.8116 , 44.9897 , 42.0126 , 32.9195 , 28.3382 , 33.0562 , 49.2503 , 40.8625 , 38.2890 ],
   [36.0192 , 44.5059 , 41.0450 , 32.3614 , 27.8164 , 32.3079 , 48.9518 , 40.4779 , 36.5053 ],
   [34.3126 , 43.9694 , 40.2155 , 31.7777 , 27.3686 , 31.6799 , 48.7429 , 40.0821 , 35.1289 ],
   [32.8892 , 43.4308 , 39.5028 , 31.2511 , 26.9793 , 31.1342 , 48.5935 , 39.6895 , 33.8864 ],
   [31.7936 , 42.9177 , 38.8887 , 30.8558 , 26.6345 , 30.6584 , 48.4874 , 39.3026 , 32.7984 ]]

obs_ENA_DJF = [21.4797, 21.3283, 21.0275, 20.9446, 20.5944, 20.0972, 19.3408, 18.8132, 20.0051, 22.5848,
               25.0611, 27.0421, 28.5685, 29.6214, 30.0936, 29.7857, 28.4903, 26.1960, 24.2754, 23.2598,
               22.7159, 22.2734, 21.9671, 21.7115]
model_ENA_DJF = [
   [34.0494 , 35.4419 , 30.5446 , 19.8759 , 17.0052 , 25.8504 , 52.8639 , 26.9949 , 25.3125 ],
   [33.6443 , 34.8739 , 29.9762 , 19.5858 , 16.8184 , 25.4586 , 52.7544 , 26.3763 , 24.4939 ],
   [33.2834 , 34.3351 , 29.4433 , 19.3425 , 16.6485 , 25.0752 , 52.6448 , 25.7841 , 23.7028 ],
   [32.9313 , 33.8364 , 28.9384 , 19.0941 , 16.4973 , 24.7123 , 52.5406 , 25.2257 , 23.0574 ],
   [32.5155 , 33.3790 , 28.4581 , 18.8408 , 16.3662 , 24.3721 , 52.4426 , 24.6966 , 22.4627 ],
   [31.9952 , 32.9396 , 28.0033 , 18.6055 , 16.2527 , 24.0479 , 52.3341 , 24.1966 , 21.8880 ],
   [31.4075 , 32.5304 , 27.5787 , 18.4039 , 16.4263 , 23.9004 , 52.1957 , 23.7210 , 22.0387 ],
   [30.8314 , 32.7041 , 27.5282 , 18.4057 , 17.5344 , 24.4309 , 52.6354 , 23.2803 , 23.4021 ],
   [30.7125 , 33.6702 , 29.5417 , 19.2545 , 18.7575 , 24.9678 , 53.3253 , 23.7530 , 24.7733 ],
   [31.0489 , 34.9184 , 32.6072 , 21.3294 , 19.8707 , 26.4866 , 53.6770 , 25.9485 , 26.7031 ],
   [31.6187 , 36.4373 , 35.0844 , 23.2601 , 20.9444 , 27.9685 , 54.0236 , 28.6509 , 28.7573 ],
   [32.5619 , 37.9825 , 36.9457 , 24.7542 , 21.8835 , 29.2238 , 54.4251 , 31.2706 , 30.0873 ],
   [35.6292 , 39.3255 , 38.4187 , 26.0106 , 22.5828 , 30.3299 , 54.8234 , 33.3469 , 32.8077 ],
   [36.7542 , 40.3092 , 39.4442 , 26.9238 , 22.9177 , 31.1415 , 55.1155 , 34.8307 , 34.6648 ],
   [37.6755 , 40.7811 , 39.9807 , 27.2880 , 22.7191 , 31.8225 , 55.1800 , 35.7146 , 35.1745 ],
   [38.7917 , 40.6988 , 39.9352 , 27.0762 , 21.7399 , 31.7088 , 54.9646 , 35.9184 , 35.0632 ],
   [38.8149 , 40.1707 , 39.0881 , 25.8557 , 20.2740 , 31.0246 , 54.4716 , 35.2153 , 34.5410 ],
   [37.9710 , 39.7022 , 37.5345 , 24.0951 , 19.3666 , 30.0231 , 53.9541 , 33.4553 , 34.3214 ],
   [36.7399 , 39.3448 , 36.0042 , 22.8748 , 18.8588 , 29.0976 , 53.7087 , 30.9430 , 32.8499 ],
   [36.5732 , 38.8463 , 34.7424 , 22.3103 , 18.4436 , 28.3671 , 53.5322 , 30.3781 , 31.1275 ],
   [36.1027 , 38.2174 , 33.6535 , 21.8013 , 18.0910 , 27.7634 , 53.3820 , 29.9538 , 29.8135 ],
   [35.5568 , 37.5442 , 32.7309 , 21.3031 , 17.7881 , 27.2417 , 53.2484 , 29.1633 , 28.5064 ],
   [35.0448 , 36.8706 , 31.9417 , 20.7886 , 17.5288 , 26.7799 , 53.1258 , 28.4224 , 27.3316 ],
   [34.5675 , 36.2160 , 31.2494 , 20.3176 , 17.3027 , 26.3513 , 53.0083 , 27.7294 , 26.2840 ]]

obs_SEU_DJF = [19.6392, 19.7349, 19.8588, 19.7813, 19.7579, 19.2165, 18.4030, 17.5446, 17.4590, 18.4974,
               20.3319, 22.3379, 24.1436, 25.3241, 25.7407, 25.1759, 23.6083, 21.6232, 20.0912, 19.4052,
               19.1933, 19.1862, 19.3279, 19.4989]
model_SEU_DJF = [
   [28.0471 , 37.0980 , 36.8089 , 27.1745 , 24.7889 , 28.4865 , 47.1324 , 31.0927 , 30.1085 ],
   [27.6573 , 36.8356 , 36.4688 , 27.0560 , 24.7109 , 28.2880 , 47.1248 , 30.6110 , 29.4763 ],
   [27.3316 , 36.5839 , 36.1524 , 26.9662 , 24.6358 , 28.1013 , 47.0662 , 31.0122 , 28.8450 ],
   [27.0789 , 36.3460 , 35.8553 , 26.8612 , 24.5713 , 27.9248 , 47.0271 , 30.7542 , 28.2745 ],
   [26.8553 , 36.1321 , 35.5820 , 26.7679 , 24.5123 , 27.7732 , 46.9906 , 30.5101 , 27.7831 ],
   [26.6023 , 35.9223 , 35.3336 , 26.7058 , 24.4561 , 27.6213 , 46.9447 , 30.2773 , 27.3039 ],
   [26.3662 , 35.6697 , 35.1097 , 26.6377 , 24.4752 , 27.4905 , 46.8942 , 30.0542 , 27.1788 ],
   [26.0749 , 35.4712 , 35.0005 , 26.6614 , 24.7582 , 27.6361 , 47.0235 , 29.8178 , 27.9081 ],
   [25.9516 , 35.5064 , 35.4801 , 27.0744 , 25.2188 , 27.9876 , 47.2893 , 29.7149 , 28.9209 ],
   [26.0859 , 35.7568 , 36.6831 , 28.0795 , 25.7571 , 28.6913 , 47.3963 , 30.2210 , 30.0269 ],
   [26.7340 , 36.3683 , 38.1624 , 29.1937 , 26.3731 , 29.6750 , 47.6041 , 31.3180 , 32.1594 ],
   [27.6093 , 37.2239 , 39.6521 , 30.1423 , 26.9665 , 30.5101 , 47.9791 , 32.9015 , 33.4838 ],
   [29.4430 , 38.0624 , 40.9892 , 30.9262 , 27.3800 , 31.2129 , 48.4066 , 34.5065 , 35.6948 ],
   [31.2168 , 38.7185 , 42.0059 , 31.4231 , 27.4970 , 31.6231 , 48.6941 , 35.6115 , 38.0204 ],
   [31.9080 , 39.0506 , 42.5595 , 31.4625 , 27.2271 , 31.8886 , 48.7340 , 36.0984 , 38.6818 ],
   [32.6057 , 39.0845 , 42.5633 , 31.0600 , 26.6848 , 31.8387 , 48.4960 , 35.9312 , 38.6484 ],
   [32.7250 , 39.0171 , 41.9624 , 30.1540 , 26.1909 , 31.3539 , 48.0878 , 35.1253 , 38.0722 ],
   [32.2427 , 39.0037 , 40.9974 , 29.2540 , 25.8751 , 30.7663 , 47.7558 , 34.0786 , 37.7065 ],
   [31.6070 , 38.9210 , 40.0616 , 28.8569 , 25.6351 , 30.2545 , 47.5739 , 33.4809 , 36.7442 ],
   [30.9452 , 38.6969 , 39.2874 , 28.5880 , 25.4333 , 29.8646 , 47.4496 , 33.1682 , 35.1503 ],
   [30.2821 , 38.3970 , 38.6639 , 28.2984 , 25.2633 , 29.5418 , 47.3557 , 32.8475 , 33.9446 ],
   [29.6803 , 38.0800 , 38.1444 , 27.9887 , 25.1171 , 29.2597 , 47.2876 , 32.5232 , 32.8452 ],
   [29.1630 , 37.7748 , 37.6928 , 27.6856 , 24.9964 , 29.0099 , 47.2359 , 32.2096 , 31.8690 ],
   [28.6999 , 37.4776 , 37.2883 , 27.4484 , 24.8931 , 28.7763 , 47.1956 , 31.8288 , 31.0031 ]]

obs_NEU_DJF = [25.2331, 25.3669, 25.5026, 25.5580, 25.5336, 25.3219, 24.8233, 24.2005, 23.9101, 24.2024,
               24.7741, 25.3970, 25.8655, 26.0926, 25.9669, 25.5360, 25.0080, 24.6475, 24.5610, 24.6473,
               24.7610, 24.8682, 25.0038, 25.1624]
model_NEU_DJF = [
   [22.7130 , 38.7028 , 39.9019 , 26.8923 , 23.1448 , 24.6907 , 51.6187 , 30.4082 , 32.4980 ],
   [22.5934 , 38.7201 , 39.8668 , 26.8670 , 23.1329 , 24.6565 , 51.6396 , 29.9719 , 32.4558 ],
   [22.4652 , 38.7343 , 39.8200 , 26.8360 , 23.1240 , 24.6276 , 51.6486 , 30.6960 , 32.3399 ],
   [22.3933 , 38.7404 , 39.7767 , 26.8021 , 23.1201 , 24.6026 , 51.6591 , 30.6517 , 32.2317 ],
   [22.3268 , 38.7563 , 39.7380 , 26.7768 , 23.1172 , 24.5805 , 51.6722 , 30.6123 , 32.1568 ],
   [22.2656 , 38.7688 , 39.7041 , 26.7689 , 23.1153 , 24.5606 , 51.6790 , 30.5754 , 32.0412 ],
   [22.2194 , 38.7691 , 39.6730 , 26.7585 , 23.1418 , 24.5433 , 51.6784 , 30.5395 , 31.9963 ],
   [22.1625 , 38.7983 , 39.6639 , 26.7474 , 23.2656 , 24.6066 , 51.7147 , 30.5130 , 32.1335 ],
   [22.1360 , 38.8552 , 39.7915 , 26.8069 , 23.4943 , 24.6874 , 51.8504 , 30.5843 , 32.3346 ],
   [22.2064 , 38.8374 , 40.1667 , 27.0811 , 23.7540 , 24.8140 , 51.9996 , 30.8407 , 32.5796 ],
   [22.3763 , 38.7280 , 40.6811 , 27.5199 , 23.9755 , 25.0586 , 52.0512 , 31.2226 , 33.2124 ],
   [22.6140 , 38.5642 , 41.1333 , 27.8859 , 24.1295 , 25.2392 , 52.0758 , 31.6283 , 33.5116 ],
   [22.9540 , 38.4256 , 41.4250 , 28.1324 , 24.1796 , 25.3951 , 52.0885 , 31.9791 , 33.8006 ],
   [23.4249 , 38.3716 , 41.5139 , 28.2040 , 24.1045 , 25.3995 , 52.0629 , 32.1937 , 34.4469 ],
   [23.6178 , 38.3611 , 41.3876 , 28.0621 , 23.9155 , 25.4269 , 51.9767 , 32.1980 , 34.5943 ],
   [23.7735 , 38.4088 , 41.1041 , 27.7602 , 23.7034 , 25.3574 , 51.8403 , 31.9772 , 34.6280 ],
   [23.8124 , 38.5321 , 40.7844 , 27.4247 , 23.5545 , 25.2262 , 51.7232 , 31.6419 , 34.3897 ],
   [23.7271 , 38.6517 , 40.5620 , 27.2672 , 23.4806 , 25.1358 , 51.6865 , 31.3873 , 34.2219 ],
   [23.5828 , 38.7191 , 40.4271 , 27.2373 , 23.4294 , 25.0608 , 51.6892 , 31.2648 , 34.0091 ],
   [23.4141 , 38.7567 , 40.3320 , 27.1937 , 23.3809 , 24.9986 , 51.6948 , 31.1820 , 33.6325 ],
   [23.2504 , 38.7771 , 40.2559 , 27.1203 , 23.3378 , 24.9428 , 51.7003 , 31.1027 , 33.3447 ],
   [23.1130 , 38.7881 , 40.1899 , 27.0413 , 23.3004 , 24.8901 , 51.7071 , 31.0264 , 33.1173 ],
   [22.9973 , 38.8074 , 40.1276 , 26.9820 , 23.2702 , 24.8418 , 51.7141 , 30.9540 , 32.9389 ],
   [22.8943 , 38.8233 , 40.0622 , 26.9575 , 23.2468 , 24.7968 , 51.7227 , 30.8463 , 32.7497 ]]

# plot DJF diurnal figures
fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=4, sharex=True,
                                    figsize=(24, 5))
ax0.plot(obs_WNA_DJF,'ko-')
ax0.plot(model_WNA_DJF)
ax0.plot(O3_WNA_DJF_24h,color='blue',marker='o')
ax0.set_ylim(0,100)
ax0.set_ylabel('DJF O3 abundance (ppb)',fontsize='x-large')
ax0.set_xlabel('WNA local time (hr)',fontsize='x-large')

ax1.plot(obs_ENA_DJF,'ko-')
ax1.plot(model_ENA_DJF)
ax1.plot(O3_ENA_DJF_24h,color='blue',marker='o')
ax1.set_ylim(0,100)
#ax1.set_ylabel('DJF O3 abundance (ppb)',fontsize='x-large')
ax1.set_xlabel('ENA local time (hr)',fontsize='x-large')

ax2.plot(obs_SEU_DJF,'ko-')
ax2.plot(model_SEU_DJF)
ax2.plot(O3_SEU_DJF_24h,color='blue',marker='o')
ax2.set_ylim(0,100)
#ax2.set_ylabel('DJF O3 abundance (ppb)',fontsize='x-large')
ax2.set_xlabel('SEU local time (hr)',fontsize='x-large')

ax3.plot(obs_NEU_DJF,'ko-')
ax3.plot(model_NEU_DJF)
ax3.plot(O3_NEU_DJF_24h,color='blue',marker='o')
ax3.set_ylim(0,100)
#ax3.set_ylabel('DJF O3 abundance (ppb)',fontsize='x-large')
ax3.set_xlabel('NEU local time (hr)',fontsize='x-large')
ax3.legend(["OBS","MOCAGE","GFDL-AM3","CESM-CAM-SF","UM-CAM","CMAM","MIROC-CHEM","GISS-E2-R","GEOSCCM","UCI CTM","E3SM"],fontsize='x-large',loc='center left', bbox_to_anchor=(1, 0.5))

fig.suptitle(short_name,fontsize='x-large')
pylab.savefig(pathout+'surfO3_hour_cycle_DJF.png',dpi=300)

# ## Calculate model daily time series of maximum daily 8 hour average (MDA8) on grid boxs
day_timearray = O3US_xr[:,0,0].sel(time=O3US_xr.time.dt.hour.isin([1]))
ndays = int(len(O3local_US)/24)
O3US_navg = np.zeros((ndays,17,nlat_US,nlon_US))
for d in range(ndays):
    nds = d*24
    nde = d*24+23
    tmp = O3local_US[nds:nde+1,:,:]
    for i in range(17):
        iend = i+8
        O3US_navg[d,i,:,:] = tmp[i:iend,:,:].mean(axis=0)

ndays = int(len(O3local_EU)/24)
O3EU_navg = np.zeros((ndays,17,nlat_EU,nlon_EU))
for d in range(ndays):
    nds = d*24
    nde = d*24+23
    tmp = O3local_EU[nds:nde+1,:,:]
    for i in range(17):
        iend = i+8
        O3EU_navg[d,i,:,:] = tmp[i:iend,:,:].mean(axis=0)

O3US_MDA8 = O3US_navg.max(axis=1)
O3EU_MDA8 = O3EU_navg.max(axis=1)
O3EU_MDA8_xr = xr.DataArray(O3EU_MDA8, coords=[day_timearray['time'][0:ndays],lat_EU,lon_EU], dims=["time","lat","lon"])
O3US_MDA8_xr = xr.DataArray(O3US_MDA8, coords=[day_timearray['time'][0:ndays],lat_US,lon_US], dims=["time","lat","lon"])

# US subregion: Western (WNA) and eastern (ENA) North America is split by 264°E
lon1 = 234
lon2 = 263
lon3 = 293
O3_MDA8_WNA = O3US_MDA8_xr.sel(lon=slice(lon1,lon2))
O3_MDA8_ENA = O3US_MDA8_xr.sel(lon=slice(lon2+1,lon3))

# EU subregion: Southern (SEU) and northern (NEU) Europe is split by 53°N. (Schnell et al., ACP 2015)
lat1 = 35
lat2 = 52
lat3 = 71
O3_MDA8_SEU = O3EU_MDA8_xr.sel(lat=slice(lat1,lat2))
O3_MDA8_NEU = O3EU_MDA8_xr.sel(lat=slice(lat2+1,lat3))

# calculate area-weighted averages over the region
O3_MDA8_ENA_1D = np.zeros(ndays)
O3_MDA8_WNA_1D = np.zeros(ndays)
O3_MDA8_NEU_1D = np.zeros(ndays)
O3_MDA8_SEU_1D = np.zeros(ndays)
for i in range(ndays):
    O3_MDA8_ENA_1D[i] = (O3_MDA8_ENA[i,:,:]*area_ENA[:,:]).sum()/area_ENA[:,:].sum()
    O3_MDA8_WNA_1D[i] = (O3_MDA8_WNA[i,:,:]*area_WNA[:,:]).sum()/area_WNA[:,:].sum()
    O3_MDA8_NEU_1D[i] = (O3_MDA8_NEU[i,:,:]*area_NEU[:,:]).sum()/area_NEU[:,:].sum()
    O3_MDA8_SEU_1D[i] = (O3_MDA8_SEU[i,:,:]*area_SEU[:,:]).sum()/area_SEU[:,:].sum()

O3_MDA8_ENA_1D = xr.DataArray(O3_MDA8_ENA_1D, coords=[day_timearray['time'][0:ndays]], dims=["time"])
O3_MDA8_WNA_1D = xr.DataArray(O3_MDA8_WNA_1D, coords=[day_timearray['time'][0:ndays]], dims=["time"])
O3_MDA8_NEU_1D = xr.DataArray(O3_MDA8_NEU_1D, coords=[day_timearray['time'][0:ndays]], dims=["time"])
O3_MDA8_SEU_1D = xr.DataArray(O3_MDA8_SEU_1D, coords=[day_timearray['time'][0:ndays]], dims=["time"])

# calculate monthly mean
O3_MDA8_ENA_month = np.zeros(12)
O3_MDA8_WNA_month = np.zeros(12)
O3_MDA8_SEU_month = np.zeros(12)
O3_MDA8_NEU_month = np.zeros(12)
for n in range(12):
    O3_MDA8_ENA_month[n] = 1.e9*O3_MDA8_ENA_1D.sel(time=O3_MDA8_ENA_1D.time.dt.month.isin([n+1])).mean()
    O3_MDA8_WNA_month[n] = 1.e9*O3_MDA8_WNA_1D.sel(time=O3_MDA8_WNA_1D.time.dt.month.isin([n+1])).mean()
    O3_MDA8_NEU_month[n] = 1.e9*O3_MDA8_NEU_1D.sel(time=O3_MDA8_NEU_1D.time.dt.month.isin([n+1])).mean()
    O3_MDA8_SEU_month[n] = 1.e9*O3_MDA8_SEU_1D.sel(time=O3_MDA8_SEU_1D.time.dt.month.isin([n+1])).mean()

# ----- writing ncfile -----
O3_MDA8_ENA_month_xr = xr.DataArray(O3_MDA8_ENA_month, name= 'O3_MDA8_ENA_month',coords=[np.arange(1,13)], dims=["month"], 
		       attrs=dict(units="ppb", description="Mean MDA8 O3 in ENA") )
O3_MDA8_WNA_month_xr = xr.DataArray(O3_MDA8_WNA_month, name= 'O3_MDA8_WNA_month',coords=[np.arange(1,13)], dims=["month"], 
		       attrs=dict(units="ppb", description="Mean MDA8 O3 in WNA") )
O3_MDA8_NEU_month_xr = xr.DataArray(O3_MDA8_NEU_month, name= 'O3_MDA8_NEU_month',coords=[np.arange(1,13)], dims=["month"], 
		       attrs=dict(units="ppb", description="Mean MDA8 O3 in NEU") )
O3_MDA8_SEU_month_xr = xr.DataArray(O3_MDA8_SEU_month, name= 'O3_MDA8_SEU_month',coords=[np.arange(1,13)], dims=["month"], 
		       attrs=dict(units="ppb", description="Mean MDA8 O3 in SEU") )
M_ds1 = O3_MDA8_ENA_month_xr.to_dataset()
M_ds2 = O3_MDA8_WNA_month_xr.to_dataset()
M_ds3 = O3_MDA8_NEU_month_xr.to_dataset()
M_ds4 = O3_MDA8_SEU_month_xr.to_dataset()
ds = xr.merge([ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, M_ds1, M_ds2, M_ds3, M_ds4])
ds.to_netcdf(pathout+'E3SM_surf_O3_${y1}-${y2}.nc')

obs_MDA8_WNA = [31.0984, 36.6661, 43.0965, 48.3994, 49.6384, 49.6071, 50.3713, 49.5312, 44.6972, 36.8627, 31.6061, 29.5710]
model_MDA8_WNA = [
   [43.1680 , 47.0648 , 44.5557 , 35.7887 , 31.2343 , 35.3478 , 54.0288 , 42.5445 , 38.9795 ],
   [46.2735 , 53.8325 , 49.9738 , 40.4311 , 36.1273 , 40.3459 , 55.5681 , 46.6724 , 46.9572 ],
   [51.6867 , 57.7473 , 59.0850 , 48.4947 , 42.7767 , 45.8694 , 53.2984 , 50.4897 , 56.6458 ],
   [55.3723 , 57.1608 , 59.5253 , 54.3579 , 46.1207 , 48.5654 , 52.6441 , 52.3024 , 63.3639 ],
   [55.5295 , 56.2589 , 52.8305 , 56.3682 , 48.2531 , 49.5094 , 50.8611 , 53.8727 , 66.4470 ],
   [54.2433 , 55.8219 , 47.5407 , 54.5623 , 49.8092 , 50.2604 , 50.6447 , 57.2156 , 69.5429 ],
   [57.4092 , 55.7206 , 44.5332 , 52.4228 , 49.2546 , 50.7349 , 49.2466 , 60.8051 , 78.0097 ],
   [57.0190 , 56.1212 , 45.4025 , 50.6025 , 48.5467 , 48.8057 , 50.1063 , 63.5242 , 73.1782 ],
   [53.4050 , 53.6400 , 45.0664 , 47.9220 , 43.8233 , 45.9149 , 50.3046 , 59.9753 , 64.1024 ],
   [49.1912 , 48.4553 , 45.1410 , 42.4838 , 37.4595 , 42.3178 , 49.4363 , 53.2331 , 51.6985 ],
   [40.0962 , 43.9325 , 45.4502 , 36.8730 , 31.2290 , 37.3144 , 51.2268 , 45.0805 , 41.3357 ],
   [36.3502 , 42.1482 , 44.3250 , 35.7777 , 29.1316 , 33.8621 , 52.7842 , 42.2301 , 36.4849 ]]

obs_MDA8_ENA = [28.1407, 34.1113, 40.8175, 46.6844, 46.6356, 46.4522, 44.8304, 43.6285, 39.8481, 32.8362, 28.7033, 26.3418]
model_MDA8_ENA = [
   [41.6565 , 41.4892 , 39.3432 , 24.9715 , 22.7300 , 29.6039 , 57.9775 , 33.4723 , 32.9278 ],
   [44.5603 , 49.0856 , 44.4851 , 31.0694 , 26.4763 , 36.6313 , 59.6039 , 39.1950 , 40.8158 ],
   [57.7220 , 56.2229 , 54.9875 , 42.4649 , 34.4897 , 44.8599 , 58.3185 , 46.3104 , 53.5632 ],
   [63.9104 , 56.7349 , 58.1684 , 50.8131 , 39.5676 , 50.8916 , 57.9459 , 50.4146 , 65.4485 ],
   [67.0825 , 59.7765 , 51.3450 , 55.2346 , 42.0875 , 55.4860 , 56.2301 , 51.9356 , 76.5619 ],
   [65.6163 , 66.3732 , 44.8331 , 55.3549 , 44.3847 , 58.6632 , 60.7842 , 51.6397 , 89.4350 ],
   [67.5176 , 69.4285 , 39.6722 , 54.4518 , 45.8832 , 60.1055 , 61.1931 , 52.6916 , 98.5627 ],
   [67.1127 , 69.6566 , 41.7218 , 53.9301 , 45.5292 , 59.5167 , 59.8131 , 55.7459 , 94.9392 ],
   [60.8915 , 61.9487 , 42.4120 , 48.8709 , 37.4028 , 54.5349 , 60.0613 , 54.8240 , 76.9298 ],
   [56.6932 , 50.5007 , 44.3311 , 38.7635 , 27.1349 , 43.5275 , 58.1663 , 46.7174 , 55.6003 ],
   [42.8081 , 40.2881 , 41.6558 , 29.2906 , 21.8218 , 34.3545 , 54.2578 , 38.2824 , 40.5866 ],
   [36.5460 , 36.5013 , 37.3797 , 24.3641 , 20.8365 , 28.2270 , 55.0493 , 33.7089 , 31.4236 ]]

obs_MDA8_SEU = [24.0730, 30.1891, 38.3144, 44.5721, 46.0596, 47.3043, 47.0860, 45.7252, 38.3107, 28.8196, 23.2101, 21.3718]
model_MDA8_SEU = [
   [32.8006 , 40.5661 , 42.4594 , 30.5010 , 28.2834 , 31.2122 , 51.5169 , 35.0653 , 36.6899 ],
   [43.8016 , 45.9531 , 47.5480 , 35.8148 , 31.8737 , 35.3949 , 52.7614 , 40.0192 , 45.0964 ],
   [50.8908 , 52.5812 , 57.1668 , 45.4018 , 37.3531 , 42.9385 , 54.2644 , 45.9386 , 56.4627 ],
   [59.4038 , 57.4834 , 60.3536 , 53.7361 , 41.5996 , 47.1585 , 54.3462 , 49.1161 , 64.9394 ],
   [64.0267 , 59.4886 , 55.1465 , 59.1822 , 44.0810 , 53.9897 , 55.9760 , 50.5648 , 73.1877 ],
   [66.6392 , 60.3116 , 49.5541 , 61.7676 , 44.0707 , 55.7243 , 60.8568 , 51.6056 , 84.4542 ],
   [70.1350 , 61.8520 , 46.2231 , 61.6699 , 43.0099 , 54.1853 , 62.8203 , 55.0551 , 82.1707 ],
   [66.3729 , 59.3456 , 46.4685 , 58.2335 , 41.0323 , 50.3617 , 65.8021 , 59.5235 , 80.9377 ],
   [57.3504 , 53.0131 , 46.1174 , 49.7043 , 37.8429 , 44.6161 , 59.4199 , 54.1577 , 60.7815 ],
   [46.6855 , 45.8750 , 43.7910 , 40.2267 , 32.0525 , 38.7502 , 54.3286 , 43.7350 , 48.3137 ],
   [34.9456 , 40.9941 , 42.7763 , 32.4176 , 28.9690 , 33.6047 , 52.1443 , 37.2586 , 39.3603 ],
   [29.0324 , 39.2547 , 39.9866 , 29.5678 , 26.9581 , 30.6259 , 51.1459 , 34.4884 , 33.6556 ]]

obs_MDA8_NEU =[28.4584, 32.0270, 38.9426, 43.5293, 42.8207, 38.2410, 34.8256, 32.6627, 29.5661, 25.7159, 24.8785, 25.4187]
model_MDA8_NEU = [
   [23.5919 , 42.4909 , 42.9792 , 29.5189 , 25.9425 , 26.3222 , 54.5146 , 33.1945 , 36.1156 ],
   [29.5744 , 46.5573 , 47.4245 , 31.9172 , 29.2139 , 29.1291 , 56.5881 , 36.2453 , 40.0708 ],
   [41.5829 , 52.1085 , 57.7614 , 39.8228 , 33.2750 , 35.0527 , 58.9328 , 41.0161 , 48.8154 ],
   [52.2116 , 54.5242 , 60.5114 , 48.8900 , 37.5822 , 40.9922 , 56.6173 , 44.7975 , 55.7258 ],
   [56.3413 , 51.3267 , 50.8769 , 49.3614 , 37.2386 , 44.7930 , 49.9556 , 44.0420 , 58.5395 ],
   [55.6993 , 47.1878 , 41.2864 , 51.3845 , 36.6013 , 44.1233 , 46.0945 , 38.5119 , 58.0179 ],
   [56.0045 , 45.3446 , 36.7440 , 47.5873 , 33.6336 , 42.5226 , 44.1008 , 36.1021 , 60.9459 ],
   [51.1850 , 43.9449 , 36.0218 , 45.1807 , 30.4656 , 39.0970 , 48.8618 , 35.3262 , 55.6645 ],
   [44.6227 , 42.0186 , 37.3661 , 40.5090 , 28.7220 , 34.2051 , 50.3438 , 35.0858 , 48.7166 ],
   [35.5239 , 38.9883 , 40.8686 , 35.0013 , 27.4796 , 30.5452 , 51.0688 , 34.6483 , 42.1752 ],
   [29.5563 , 38.0949 , 42.5559 , 30.3209 , 24.8675 , 27.3245 , 53.9359 , 35.7548 , 36.8644 ],
   [29.3354 , 39.1637 , 42.4841 , 28.9785 , 25.3637 , 25.6507 , 53.4103 , 35.9839 , 34.4121 ]]

# plot monthly figures
fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=4, sharex=True,
                                    figsize=(24, 5))
monthlist = np.arange(1,13)
ax0.plot(monthlist,obs_MDA8_WNA,'ko-')
ax0.plot(monthlist,model_MDA8_WNA)
ax0.plot(monthlist,O3_MDA8_WNA_month,color='blue',marker='o')
ax0.set_ylim(0,100)
ax0.set_ylabel('Mean MDA8 O3 (ppb)',fontsize='x-large')
ax0.set_xlabel('WNA (month)',fontsize='x-large')

ax1.plot(monthlist,obs_MDA8_ENA,'ko-')
ax1.plot(monthlist,model_MDA8_ENA)
ax1.plot(monthlist,O3_MDA8_ENA_month,color='blue',marker='o')
ax1.set_ylim(0,100)
#ax1.set_ylabel('DJF O3 abundance (ppb)',fontsize='x-large')
ax1.set_xlabel('ENA (month)',fontsize='x-large')

ax2.plot(monthlist,obs_MDA8_SEU,'ko-')
ax2.plot(monthlist,model_MDA8_SEU)
ax2.plot(monthlist,O3_MDA8_SEU_month,color='blue',marker='o')
ax2.set_ylim(0,100)
#ax2.set_ylabel('DJF O3 abundance (ppb)',fontsize='x-large')
ax2.set_xlabel('SEU (month))',fontsize='x-large')

ax3.plot(monthlist,obs_MDA8_NEU,'ko-')
ax3.plot(monthlist,model_MDA8_NEU)
ax3.plot(monthlist,O3_MDA8_NEU_month,color='blue',marker='o')
ax3.set_ylim(0,100)
#ax3.set_ylabel('DJF O3 abundance (ppb)',fontsize='x-large')
ax3.set_xlabel('NEU (month)',fontsize='x-large')
ax3.legend(["OBS","MOCAGE","GFDL-AM3","CESM-CAM-SF","UM-CAM","CMAM","MIROC-CHEM","GISS-E2-R","GEOSCCM","UCI CTM","E3SM"],fontsize='x-large',loc='center left', bbox_to_anchor=(1, 0.5))

fig.suptitle(short_name,fontsize='x-large')
pylab.savefig(pathout+'MDA8_month_cycle.png',dpi=300)

EOF

# Run diagnostics
command="python -u surf_O3_diags.py"
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

