#!/bin/bash
{% include 'inclusions/slurm_header.bash' %}
{{ environment_commands }}

# To load custom E3SM Diags environment, comment out line above using {# ... #}
# and uncomment lines below

#module load anaconda3/2019.03
#source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh
#conda activate e3sm_diags_env_dev
#module load ncl

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

# Basic definitions
case="{{ case }}"
short="{{ short_name }}"
www="{{ www }}"
y1={{ year1 }}
y2={{ year2 }}
Y1="{{ '%04d' % (year1) }}"
Y2="{{ '%04d' % (year2) }}"
run_type="{{ run_type }}"
obsDir="{{ diagnostics_base_path }}/observations/Atm/ChemDyg_inputs"
ncfile_save="{{ ncfile_save }}"
if [[ "${ncfile_save}" == "true" ]]; then
   results_dir={{ output }}/post/atm/chemdygfiles
   mkdir -p ${results_dir}
fi

# Create temporary workdir
workdir=`mktemp -d tmp.${id}.XXXX`
cd ${workdir}

# Create local links to input climo files
tsDir={{ output }}/post/atm/{{ grid }}/ts/monthly/{{ '%dyr' % (ypf) }}
mkdir -p ts
#
ln -s ${obsDir}/QBOmetrics/CMZM_198410-202212.nc ./ts
ln -s ${obsDir}/QBOmetrics/ERA5_T_197901-202012.nc ./ts
ln -s ${obsDir}/QBOmetrics/MSR-ASSIM_197901-202212.nc ./ts
ln -s ${obsDir}/QBOmetrics/linoz_1850-2500_CMIP6_Hist_SSP370_10deg_58km_c20210202.nc ./ts
#
ln -s ${obsDir}/QBOmetrics/pca_index.nc ./ts
#sim
ln -s ${tsDir}/TCO_*.nc ./ts
ln -s ${tsDir}/SCO_*.nc ./ts
ln -s ${tsDir}/T_*.nc  ./ts
ln -s ${tsDir}/O3_*.nc  ./ts

# Run E3SM chem Diags
echo
echo ===== RUN E3SM CHEM DIAGS  =====
echo

# Prepare configuration file
 cat > QBO_metrics.py << EOF
#!/usr/bin/env python
# coding: utf-8
#
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab
##=====================
#1 functions for use
##=====================
#########################
def pick_trans(angle):
    angle_ind=np.copy(angle)
    angle_ind[:]=0
    ##
    for j in range(0,len(angle)-1-1,1):
        ##
        if (angle[j] > 0. and angle[j+1] < 0.):
            angle_ind[j]=1
        elif (angle[j] < 0. and angle[j+1] > 0.):
            angle_ind[j]=-1
        ##
    return angle_ind
#########################
def qbo_avg(input1,angle_ind,lat1,lat2,index,index2):
    ##
    dims_nam=input1.dims
    if (index2 == 0):
        input2=input1.rename({dims_nam[1]:"lat",dims_nam[2]:"lev"})
    else           : 
        input2=input1.copy(deep=True)
    ##
    neg=np.where(angle_ind==-1)[0]
    pos=np.where(angle_ind==1)[0]
    ##
    if (index == 0): ##;;for negative
            phase=neg
    else           : ##;;for positive
            phase=pos
    ##
    dims=input1.shape
    aa_tim_med=np.zeros([len(phase),29,dims[2]],dtype=float)*np.nan
    ##
    for j in range(len(phase)-1):
        ##
        if ((phase[j]-14 < 0) or \
            (phase[j]+14+1 > dims[0])):
                        continue
        ##
        aa_280_med=input2[phase[j]-14:phase[j]+14+1,:,:].sel(lat=slice(lat1,lat2))
        rad    = 4.0*np.arctan(1.0)/180.0
        clat   = np.cos(aa_280_med.lat*rad)
        aa_280=aa_280_med*clat
        aa_2801= aa_280.sum(dim="lat")/sum(clat)
        aa_tim_med[j,:,:]=aa_2801[:,:]
        ##
    ##calculate and give back metadata
    output=np.nanmean(aa_tim_med,axis=0)
    output=aa_2801.copy(deep=True,data=output)
    ##
    output["time"]=range(-14,14+1,1)
    return output
#######################
#######################
def qbo_avg_srf(input1,angle_ind,index,index2):
    ##  
    dims_nam=input1.dims
    if (index2 == 1):
        input2=input1.rename({dims_nam[1]:"lat"})
    else	    :
        input2=input1.copy(deep=True)
    ##  
    neg=np.where(angle_ind==-1)[0]
    pos=np.where(angle_ind==1)[0]
    ##
    if (index == 0): ##;;for negative
            phase=neg
    else           : ##;;for positive
            phase=pos
    ##
    dims=input1.shape
    ##
    if (index2 == 1):
        aa_tim_med=np.zeros([len(phase),29,dims[1]],dtype=float)*np.nan
    else	    :
        aa_tim_med=np.zeros([len(phase),29,dims[1]],dtype=np.float64)*np.nan
    ##
    for j in range(len(phase)-1):
        #print(str(phase[j]-14)+" "+str((phase[j]+14)))
        ##
        if ((phase[j]-14 < 0) or \
            (phase[j]+14+1 > dims[0])):
                        continue
        aa_280_med=input2[phase[j]-14:phase[j]+14+1,:]
        aa_tim_med[j,:,:]=aa_280_med
        ##
    ##calculate and give back metadata
    output=np.nanmean(aa_tim_med,axis=0)
    output=aa_280_med.copy(deep=True,data=output)
    ##
    output["time"]=range(-14,14+1,1)
    return output
##=============
def formula(T,o3col,var_comb):
        ##
        o3_clim    =var_comb.o3_clim
        t_clim     =var_comb.t_clim
        o3col_clim =var_comb.o3col_clim
        PmL_clim   =var_comb.PmL_clim
        dPmL_dO3   =var_comb.dPmL_dO3
        dPmL_dT    =var_comb.dPmL_dT
        dPmL_dO3col=var_comb.dPmL_dO3col
        ##
        tau= -1./(dPmL_dO3+1e-30)
        ##
        fss= o3_clim + (PmL_clim + dPmL_dT * (T.values-t_clim) + dPmL_dO3col * ( o3col - o3col_clim ) ) * tau
        xr.where(fss<0.,0.,fss)
        ##
        return fss,tau
#######################
def get_anom(input1):
    ##
    output=input1.copy(deep=True)
    ##
    for i in range(0,11+1):
        output[i::12,:,:]=input1[i::12,:,:]-input1[i::12,:,:].mean(dim="time")
    return(output)
#######################
def get_anom2d(input1):
    ##
    output=input1.copy(deep=True)
    ##
    for i in range(0,11+1):
        output[i::12,:]=input1[i::12,:]-input1[i::12,:].mean(dim="time")
    return(output)
#####################
#input variable extension
def datvar_ext(fil,yr1,yr2):
    vars=["o3_clim", "t_clim", "o3col_clim" ,"PmL_clim" ,"dPmL_dO3" , "dPmL_dT"  ,"dPmL_dO3col"]
    ##
    ff1=xr.open_dataset(fil,decode_times=False)
    ##
    dat=ff1.date[312:431+1]
    dat_ful=ff1.date[0:504].copy(deep=True)
    #
    ilev27=ff1.ilev
    ##(/7,120,27,20/)
    var=ff1[vars[:]].isel(time=slice(312,431+1)).copy(deep=True)
    ##(/7,504,27,20/)
    var_ful=ff1[vars[:]].isel(time=slice(0,504)).copy(deep=True)
    ##
    ##use 5 yrs interval data to fill the data
    for i in range(0,6+1):
        #1979 assign 1975
        dat_ful[0:11+1]=dat[0:11+1].values+40000
        var_ful[vars[i]][0:11+1,:,:]=var[vars[i]][0:11+1,:,:].values
        #1980-2019 assign ispan(1980,2020,5)
        for j in range(0,7+1):
            for k in range(0,4+1):
                dat_ful[12+12*k+60*j:23+12*k+60*j+1]=dat[12+12*j:23+12*j+1].values+k*10000
                var_ful[vars[i]][12+12*k+60*j:23+12*k+60*j+1,:,:]=\
                    var[vars[i]][12+12*j:23+12*j+1,:,:].values
        #2020 assign 2020
        j=8
        dat_ful[492:503+1]=dat[12+12*j:23+12*j+1].values
        var_ful[vars[i]][492:503+1,:,:]=\
            var[vars[i]][12+12*j:23+12*j+1,:,:].values

    ####
    y1i=(yr1-1979)*12
    y2i=(yr2-1979+1)*12
    ####
    return dat_ful[y1i:y2i],var_ful.isel(time=slice(y1i,y2i)),ilev27
    #return dat_ful,var_ful,ilev27
    #return dat_ful[y1i:y2i],var_ful[:][y1i:y2i,:,:],ilev27
    ##
#######################
#ozone iterate
def iterate(T,o3col,var_comb,dp):
        ##constants
        Av = 6.022e23
        grav = 9.8
        mair = 28.97
        ##
        o3_out=var_comb.o3_clim.copy(deep=True)
        o3_out[:]=0.
        ##
        o3col_out=o3_out.copy(deep=True)
        o3col_out[:]=0.
        ##
        tau_out=o3col_out.copy(deep=True)
        tau_out[:]=0.
        ##o3col is the overhead ozone above the current level
        sp=T.shape
        ##
        for i in range(0,sp[1]-2):
        ##no need to calculate bottom
            ##derive o3 mixing ratio of next level
            o3_ss=formula(T.isel(lev=i),o3col[:,i,:],var_comb.isel(lev=i))
            ##transition of o3 mix ratio to o3col
            o3_out[:,i,:]=o3_ss[0][:,:]
            tau_out[:,i,:]=o3_ss[1][:,:]
            ##
            product=o3_out[:,0:i+1,:]*dp[0:i+1,np.newaxis]
            ##transfer to DU
            o3col_out[:,i+1,:]=\
            product.sum()*(1.e3*Av/grav/mair/2.687e20)
            ##update the o3col after each loop
            o3col[:,i+1+1,:]=o3col_out[:,i+1+1,:].values
            ##
        return o3_out,tau_out
        ##
#######################
def cal_steady_o3(T,var_comb,ff1):
    ##
    p=ff1.lev.copy(deep=True)
    sp=p.shape
    ##get dp
    dp=p.isel(lev=slice(1,sp[0]-1+1)).values-p.isel(lev=slice(0,sp[0]-2+1)).values
    dp=dp*100.
    #
    Tsp=T.shape
    ##overhead ozone
    o3col_qboi=np.zeros((Tsp[0],sp[0],20))*np.nan  ##overhead ozone
    qboi=iterate(T3_low,o3col_qboi,var_comb,dp)
    o3_qboi=qboi[0].copy(deep=True)
    ##
    o3_anom=get_anom(o3_qboi)
    T_qboi=get_anom(T3_low)
    ss=o3_anom-o3_qboi
    T_full=T3_low-273.15
    return o3_anom,T_qboi,qboi,T_full
##########################
##==============
##2 read in data
##==============
#check if model data has sufficient time
if (${y2}-${y1}<2):
	print("Insufficent data for QBO diagnostics, recommend at least 3 yrs, exiting")
	quit()
#
path = './ts/'
pathout = './'

#read in QBO index
fil1="pca_index.nc"
pca_ind=xr.open_dataset(path+fil1,decode_times=False)
#pick QBOi transition phase
pca5=pca_ind.angl[:,5].values
phase=pick_trans(pca5)
##read in obs
fil2="MSR-ASSIM_197901-202212.nc"
tco=xr.open_dataset(path+fil2,decode_times=False)
tco_anom=get_anom2d(tco.total_ozone_column)
tco_lat=tco_anom.mean(dim="longitude")
##
fil3="CMZM_198410-202212.nc"
ozo=xr.open_dataset(path+fil3,decode_times=False)
ozo2=xr.open_dataset(path+fil3,decode_times=False)
##use 1985-2020
ozo1=ozo.merged_ozone_concentration[3:75,:,:].copy(deep=True)
ozo1["time"]=ozo1["time"]-600
ozo1[:]=np.nan
ozo2=ozo.merged_ozone_concentration[3:435,:,:]
#connect to form 1979-2020, 1979-1984 is missing value
ozo_comb=xr.concat([ozo1,ozo2],dim="time")
ozo_lev=ozo_comb 
#
#change level dimension in CMZM from altitude(km) to pressure(hPa), reference calculation using 1976 US Standard Atmosphere from NCL function stdatmus_z2tdp
p_ozo=np.array(\
[264.3627,226.3206,193.304,165.1041,141.018,120.4457,102.8746,87.86682,75.04845,64.10007,54.74889,
46.77886,39.9979,34.22434,29.30492,25.11023,21.53094,18.47457,15.8629,13.62965,11.71867,10.08232,
8.680187,7.482283,6.461222,5.589235,4.843167,4.203671,3.654547,3.182204,2.775216,2.423955,2.120299,
1.857378,1.629374,1.431348,1.259103,1.109063,0.9775451,0.8616231,0.7594477])
##
#read in E3SM
#E3SM tco
fil2="TCO_${y1}01_${y2}12.nc"
tco_sim1=xr.open_dataset(path+fil2,decode_times=False)
fil2="SCO_${y1}01_${y2}12.nc"
tco_sim2=xr.open_dataset(path+fil2,decode_times=False)
tco_sim=tco_sim1.TCO.mean(dim="lon")+tco_sim2.SCO.mean(dim="lon")
tco_anom_sim=get_anom2d(tco_sim)
tco_lat_sim=tco_anom_sim
#E3SM ozone
fil="O3_${y1}01_${y2}12.nc"
ozo3=xr.open_dataset(path+fil,decode_times=False)
ozo_lev_sim=get_anom(ozo3.O3.mean(dim="lon"))
#E3SM T
filT="T_${y1}01_${y2}12.nc"
T3=xr.open_dataset(path+filT,decode_times=False)
T_lev_sim=get_anom(T3.T.mean(dim="lon"))
#
#read in linoz forcing
fil_linoz="linoz_1850-2500_CMIP6_Hist_SSP370_10deg_58km_c20210202.nc"
linoz=xr.open_dataset(path+fil_linoz,decode_times=False)
ss1=datvar_ext(path+fil_linoz,${y1},${y2})
var_comb=ss1[1]
##==============
##3 process data
##==============
##
#fig1 MSR composite
tco_lat=qbo_avg_srf(tco_lat,phase,0,1)
tco_lat_sim=qbo_avg_srf(tco_lat_sim,phase,0,0)
##fig2 O3 composite
ozo_lev=get_anom(ozo_lev)
#tropical (15S-15N)
lat1=-15
lat2=15
o3_lev_trop=qbo_avg(ozo_lev,phase,lat1,lat2,0,0)
o3_lev_trop["lev"]=p_ozo
o3_lev_trop_sim=qbo_avg(ozo_lev_sim.transpose("time","lat","lev"),phase,lat1,lat2,0,1)
#extraopical
lat1=30
lat2=60
o3_lev_Next=qbo_avg(ozo_lev,phase,lat1,lat2,0,0)
o3_lev_Next["lev"]=p_ozo
o3_lev_Next_sim=qbo_avg(ozo_lev_sim.transpose("time","lat","lev"),phase,lat1,lat2,0,1)
lat1=-60
lat2=-30
o3_lev_Sext=qbo_avg(ozo_lev,phase,lat1,lat2,0,0)
o3_lev_Sext["lev"]=p_ozo
o3_lev_Sext_sim=qbo_avg(ozo_lev_sim.transpose("time","lat","lev"),phase,lat1,lat2,0,1)
##
##fig3 steady state ozone
##
#return full var_comb and input to calculation
#
##interpolate
T3_low=T3.T.mean(dim="lon").interp(lat=linoz.lat,lev=linoz.lev,method="linear") 
#linoz get to 1979-2020
#
ss_o3=cal_steady_o3(T3_low,var_comb,linoz) #,${y1},${y2}) 
#
lat1=-15
lat2=15
ss_o3_lev_trop=qbo_avg(ss_o3[0].transpose("time","lat","lev"),phase,lat1,lat2,0,1)
ss_t_lev_trop=qbo_avg(ss_o3[1].transpose("time","lat","lev"),phase,lat1,lat2,0,1)
#
lat1=30
lat2=60
ss_o3_lev_Next=qbo_avg(ss_o3[0].transpose("time","lat","lev"),phase,lat1,lat2,0,1)
ss_t_lev_Next=qbo_avg(ss_o3[1].transpose("time","lat","lev"),phase,lat1,lat2,0,1)
#
lat1=-60
lat2=-30
ss_o3_lev_Sext=qbo_avg(ss_o3[0].transpose("time","lat","lev"),phase,lat1,lat2,0,1)
ss_t_lev_Sext=qbo_avg(ss_o3[1].transpose("time","lat","lev"),phase,lat1,lat2,0,1)
##turn to ppm for ozone
o3_lev_trop         =o3_lev_trop*1e+6
o3_lev_trop_sim     =o3_lev_trop_sim*1e+6
o3_lev_Next         =o3_lev_Next*1e+6
o3_lev_Next_sim     =o3_lev_Next_sim*1e+6
o3_lev_Sext         =o3_lev_Sext*1e+6
o3_lev_Sext_sim     =o3_lev_Sext_sim*1e+6
##
ss_o3_lev_trop         =ss_o3_lev_trop*1e+6
ss_o3_lev_Next         =ss_o3_lev_Next*1e+6
ss_o3_lev_Sext         =ss_o3_lev_Sext*1e+6
##================================================
#4 writing ncfile
##================================================
tco_lat_xr=tco_lat.to_dataset(name='tco_lat')
tco_lat_sim_xr=tco_lat_sim.to_dataset(name='tco_lat_sim')
#
o3_lev_trop_xr=o3_lev_trop.to_dataset(name='o3_lev_trop')
o3_lev_Next_xr=o3_lev_Next.to_dataset(name='o3_lev_Next')
o3_lev_Sext_xr=o3_lev_Sext.to_dataset(name='o3_lev_Sext')
#
ss_o3_lev_trop_xr=ss_o3_lev_trop.to_dataset(name='ss_o3_lev_trop')
ss_o3_lev_Next_xr=ss_o3_lev_Next.to_dataset(name='ss_o3_lev_Next')
ss_o3_lev_Sext_xr=ss_o3_lev_Sext.to_dataset(name='ss_o3_lev_Sext')
#
ss_t_lev_trop_xr=ss_t_lev_trop.to_dataset(name='ss_t_lev_trop')
ss_t_lev_Next_xr=ss_t_lev_Next.to_dataset(name='ss_t_lev_Next')
ss_t_lev_Sext_xr=ss_t_lev_Sext.to_dataset(name='ss_t_lev_Sext')
ds = xr.merge([tco_lat_xr,tco_lat_sim_xr,o3_lev_trop_xr,o3_lev_Next_xr,o3_lev_Sext_xr, ss_o3_lev_trop_xr,
    ss_o3_lev_Next_xr,ss_o3_lev_Sext_xr,ss_t_lev_trop_xr,ss_t_lev_Next_xr,ss_t_lev_Sext_xr])
ds.to_netcdf(pathout+'E3SM_QBO_Metric_${y1}-${y2}.nc')
##================================================
#5 draw plot
##================================================
fig = plt.figure(figsize=(18,12))
##
plt.subplot(2, 2, 1)
plot1 = plt.contourf(tco_lat['time'], tco_lat['lat'], tco_lat.transpose("lat","time"),levels=np.linspace(-22,22,12),cmap = 'RdBu_r',extend='both')
plt.title('Latitude-time plot of MSR TCO anomaly (DU) (1979-2020)\n ')
plt.ylim(-60,60)
plt.xlim(-14,14)
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Latitude')
cross_colorbar = fig.colorbar(plot1,ticks=np.linspace(-22,22,12))
##
plt.subplot(2, 2, 2)
plot2 = plt.contourf(tco_lat_sim['time'], tco_lat_sim['lat'], tco_lat_sim.transpose("lat","time"),levels=np.linspace(-22,22,12),cmap = 'RdBu_r',extend='both')
plt.title('Latitude-time plot of E3SM TCO anomalous (DU) (${y1}-${y2})\n ')
plt.ylim(-60,60)
plt.xlim(-14,14)
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Latitude')
cross_colorbar = fig.colorbar(plot2, ticks=np.linspace(-22,22,12))
plt.savefig("QBO_fig1.png", bbox_inches='tight', dpi=600)
##
fig = plt.figure(figsize=(18,24))
##
plt.subplot(4, 2, 1)
plot3 = plt.contourf(o3_lev_trop['time'], o3_lev_trop['lev'], o3_lev_trop.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of CMZM O3 anomaly (ppm) (15S-15N) (1985-2020)\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot3, ticks=np.linspace(-0.4,0.4,11))
##
plt.subplot(4, 2, 2)
plot4 = plt.contourf(o3_lev_trop_sim['time'], o3_lev_trop_sim['lev'], o3_lev_trop_sim.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of E3SM O3 anomaly (ppm) (15S-15N) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot4, ticks=np.linspace(-0.4,0.4,11))
##
plt.subplot(4, 2, 3)
plot3 = plt.contourf(o3_lev_Next['time'], o3_lev_Next['lev'], o3_lev_Next.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of CMZM O3 anomaly (ppm) (30N-60N) (1985-2020)\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot3, ticks=np.linspace(-0.4,0.4,11))
##
plt.subplot(4, 2, 4)
plot3 = plt.contourf(o3_lev_Next_sim['time'], o3_lev_Next_sim['lev'], o3_lev_Next_sim.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of E3SM O3 anomaly (ppm) (30N-60N) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot3, ticks=np.linspace(-0.4,0.4,11))
##
plt.subplot(4, 2, 5)
plot3 = plt.contourf(o3_lev_Sext['time'], o3_lev_Sext['lev'], o3_lev_Sext.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of CMZM O3 anomaly (ppm) (30S-60S) (1985-2020)\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot3, ticks=np.linspace(-0.4,0.4,11))
###
plt.subplot(4, 2, 6)
plot3 = plt.contourf(o3_lev_Sext_sim['time'], o3_lev_Sext_sim['lev'], o3_lev_Sext_sim.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of E3SM O3 anomaly (ppm) (30S-60S) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot3, ticks=np.linspace(-0.4,0.4,11))
##
plt.savefig("QBO_fig2.png", bbox_inches='tight',  dpi=600)
##
fig = plt.figure(figsize=(32,24))
##
plt.subplot(3, 3, 1)
plot1 = plt.contourf(ss_t_lev_trop['time'], ss_t_lev_trop['lev'], ss_t_lev_trop.transpose("lev","time"),levels=np.linspace(-2,2,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of E3SM temperature anomaly (K) (15S-15N) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot1, ticks=np.linspace(-2,2,11))
##
plt.subplot(3, 3, 2)
plot2 = plt.contourf(ss_o3_lev_trop['time'], ss_o3_lev_trop['lev'], ss_o3_lev_trop.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of E3SM Steady State O3 anomaly (ppm) (15S-15N) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot2, ticks=np.linspace(-0.4,0.4,11))
##
plt.subplot(3, 3, 3)
plot3 = plt.contourf(o3_lev_trop_sim['time'], o3_lev_trop_sim['lev'], o3_lev_trop_sim.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title(' Pressure-time plot of E3SM O3 anomaly (ppm) (30N-60N) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot3, ticks=np.linspace(-0.4,0.4,11))
##
plt.subplot(3, 3, 4)
plot4 = plt.contourf(ss_t_lev_Next['time'], ss_t_lev_Next['lev'], ss_t_lev_Next.transpose("lev","time"),levels=np.linspace(-2,2,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of E3SM temperature anomaly (K) (30N-60N) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot4, ticks=np.linspace(-2,2,11))
##
plt.subplot(3, 3, 5)
plot5 = plt.contourf(ss_o3_lev_Next['time'], ss_o3_lev_Next['lev'], ss_o3_lev_Next.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of E3SM Steady State O3 anomaly (ppm) (30N-60N) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot5, ticks=np.linspace(-0.4,0.4,11))
##

plt.subplot(3, 3, 6)
plot6 = plt.contourf(o3_lev_Next_sim['time'], o3_lev_Next_sim['lev'], o3_lev_Next_sim.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of E3SM O3 anomaly (ppm) (30N-60N) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot6, ticks=np.linspace(-0.4,0.4,11))
##
plt.subplot(3, 3, 7)
plot7 = plt.contourf(ss_t_lev_Sext['time'], ss_t_lev_Sext['lev'], ss_t_lev_Sext.transpose("lev","time"),levels=np.linspace(-2,2,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of E3SM temperature anomaly  (K) (30S-60S) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot7, ticks=np.linspace(-2,2,11))
##
plt.subplot(3, 3, 8)
plot8 = plt.contourf(ss_o3_lev_Sext['time'], ss_o3_lev_Sext['lev'], ss_o3_lev_Sext.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of E3SM Steady State O3 anomaly (ppm) (30S-60S) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot8, ticks=np.linspace(-0.4,0.4,11))
##
plt.subplot(3, 3, 9)
plot9 = plt.contourf(o3_lev_Sext_sim['time'], o3_lev_Sext_sim['lev'], o3_lev_Sext_sim.transpose("lev","time"),levels=np.linspace(-0.4,0.4,11),cmap = 'RdBu_r',extend='both')
plt.title('Pressure-time plot of E3SM O3 anomaly (ppm) (30S-60S) (${y1}-${y2})\n ')
plt.ylim(100,1)
plt.xlim(-14,14)
plt.yscale('log')
plt.xlabel('Time (month, QBOE->QBOW)')
plt.ylabel('Pressure (hPa)')
cross_colorbar = fig.colorbar(plot9, ticks=np.linspace(-0.4,0.4,11))
plt.savefig("QBO_fig3.png", bbox_inches='tight',  dpi=600)
##

EOF

# Run diagnostics
command="python QBO_metrics.py"
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


