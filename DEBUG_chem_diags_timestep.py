#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from calendar import monthrange
import pandas as pd
from time import process_time
import os

### add your input infor here ###
# File name
short_name = '20221103.v2.LR.amip.NGD_v3atm.chrysalis'
# E3SM h0 and h1 directory
path = '/lcrc/group/e3sm/ac.mwu/archive/20221103.v2.LR.amip.NGD_v3atm.chrysalis/archive/atm/hist/'
# desired output path
pathout = '/lcrc/group/e3sm/public_html/diagnostic_output/ac.lee1061'

filename = short_name+'.eam.h0.1985*.nc'
filenameh1 = short_name+'.eam.h1.1985*.nc'
################################

htmlpath = pathout+'/'+short_name+'/html'
textpath = pathout+'/'+short_name+'/text'

if os.path.exists(htmlpath) == False:
    os.mkdir(htmlpath)
if os.path.exists(textpath) == False:
    os.mkdir(textpath)

varname = ["O3","OH","HO2","H2O2","CH2O","CH3O2","CH3OOH","NO","NO2","NO3","N2O5",
           "HNO3","HO2NO2","PAN","CO","C2H6","C3H8","C2H4","ROHO2","CH3COCH3","C2H5O2",
           "C2H5OOH","CH3CHO","CH3CO3","ISOP","ISOPO2","MVKMACR",
           "MVKO2","E90","N2OLNZ","NOYLNZ","CH4LNZ", "H2OLNZ","DMS","SO2","H2SO4","SOAG"]
mol_mass = [   
   47.9982000000000,        17.0068000000000,        33.0062000000000,
   34.0136000000000,        30.0252000000000,        47.0320000000000,
   48.0394000000000,        30.0061400000000,        46.0055400000000,
   62.0049400000000,        108.010480000000,        63.0123400000000,
   79.0117400000000,        121.047940000000,        28.0104000000000,
   30.0664000000000,        44.0922000000000,        28.0516000000000,
   77.0572000000000,        58.0768000000000,        61.0578000000000,
   62.0652000000000,        44.0510000000000,        75.0424000000000,
   68.1142000000000,        117.119800000000,        70.0878000000000,
   119.093400000000,        47.9982000000000,        28.0134800000000,
   14.0067400000000,        16.0406000000000,        18.0142000000000,
   62.1324000000000,        64.0648000000000,        98.0784000000000,
   12.0110000000000 ]

WD_list = ["C2H5OOH","CH2O","CH3CHO","CH3OOH","H2O2","H2SO4","HNO3","HO2NO2","SO2"]
trop_list =['O3','O3LNZ','N2OLNZ','CH4LNZ']
layer = ['','_L1','_L2','_L3','_L4','_trop']

refer_in = xr.open_mfdataset(path+filename)
h1_in = xr.open_mfdataset(path+filenameh1)
variablelist = list(refer_in.keys())

seasons = refer_in['time']

O3_thrd = 150.e-9  # ppb
p_thrd = 100.e0  # hPa
rearth = 6.37122e6 # Earth radius: m

area_rad = refer_in['area']         # radian (ncol)
area = area_rad * rearth * rearth  # m2

for ss in range(len(seasons)):
    mass = refer_in['MASS']
    mass_L1 = (mass[ss,0:26,:].sum()/mass[ss,:,:].sum()) * 100
    mass_L2 = (mass[ss,26:38,:].sum()/mass[ss,:,:].sum()) * 100
    mass_L3 = (mass[ss,38:58,:].sum()/mass[ss,:,:].sum()) * 100
    mass_L4 = (mass[ss,58::,:].sum()/mass[ss,:,:].sum()) * 100

    linehead = '<h> E3SM chem tendency check (MSD units: Tg; TD? units: Tg/year)</h>'
    linehead = linehead + '<pre> (WD:wet deposition; TD:tendency; TR: tendency after the stratosphere reset) </pre>'
    linehead = linehead + '<pre> (O:processes outside of chemistry; I:implicit solver; R: reaction rate reset; E:explicit solver; A:aero_model_gasaerexch; L:Linoz; N:reset negative values to zero; U:setting upper boundary values; B:setting lower boundary values; S:surface emission; D:dry deposition) </pre>'
    linehead = linehead + '<pre> (L1: top_of_model to 100 hPa (air mass '+"{0:.1f}".format(np.array(mass_L1))+'%); L2: 100 to 267 hPa (air mass '+"{0:.1f}".format(np.array(mass_L2))+'%); L3: 267 to 856 hPa (air mass '+"{0:.1f}".format(np.array(mass_L3))+'%); L4: 856 hPa to surface (air mass '+"{0:.1f}".format(np.array(mass_L4))+'%)) </pre>'
    linehead = linehead + '<pre>'+short_name+'</pre>'

    line_ann = linehead + '<p> Time: '+str(ss)+' </p>'
    line_ann = line_ann + '<pre> Chemistry           TDO            TDI            TRI            TDE            TRE            TDA            TDL            TDN            TDU            TDB            TDS             TDD            WD          total_TD         MSD            VMR    </pre>'
    line_txt = ' Chemistry           TDO            TDI            TRI            TDE            TRE            TDA            TDL            TDN            TDU            TDB            TDS             TDD            WD          total_TD         MSD            VMR     \n'

    fileout_ann = open(htmlpath+'/chem_clim_t'+str(ss)+'.html',"w")
    fileout_txt = open(textpath+'/chem_clim_t'+str(ss)+'.txt',"w")

    h0_in = refer_in

    for var in range(len(varname)):
        if varname[var] in trop_list:
            total_layer = len(layer)
        else:
            total_layer = len(layer)-1

        if varname[var]+'_2DTDI' in variablelist:
            print(varname[var])
        else:
            print(varname[var]+' not exist')
            continue

        for ll in range(total_layer):
            MSD = h1_in[varname[var]+'_2DMSD'+layer[ll]][ss,:] #kg/m2

            TDB = h0_in[varname[var]+'_2DTDB'+layer[ll]][ss,:] #kg/m2/sec
            TDD = h0_in[varname[var]+'_2DTDD'+layer[ll]][ss,:]
            TDE = h0_in[varname[var]+'_2DTDE'+layer[ll]][ss,:]
            TDI = h0_in[varname[var]+'_2DTDI'+layer[ll]][ss,:]
            TRI = h0_in[varname[var]+'_2DTRI'+layer[ll]][ss,:]
            TRE = h0_in[varname[var]+'_2DTRE'+layer[ll]][ss,:]
            TDA = h0_in[varname[var]+'_2DTDA'+layer[ll]][ss,:]
            TDL = h0_in[varname[var]+'_2DTDL'+layer[ll]][ss,:]
            TDN = h0_in[varname[var]+'_2DTDN'+layer[ll]][ss,:]
            TDO = h0_in[varname[var]+'_2DTDO'+layer[ll]][ss,:]
            TDS = h0_in[varname[var]+'_2DTDS'+layer[ll]][ss,:]
            TDU = h0_in[varname[var]+'_2DTDU'+layer[ll]][ss,:]
    
            if varname[var] in WD_list:
                WD = h0_in['WD'+layer[ll]+'_'+varname[var]]#kg/m2/sec
                total_td = (WD+TDO+TDE+TDI+TRI+TRE+TDA+TDL+TDN+TDU+TDB+TDS+TDD)
                WD_total = (WD*area).sum()
            else:
                total_td = (TDO+TDE+TDI+TRI+TRE+TDA+TDL+TDN+TDU+TDB+TDS+TDD)
    
    # annunal mean
            TDO_total = ((TDO*area).sum()) #kg
            TDE_total = ((TDE*area).sum())
            TDI_total = ((TDI*area).sum())
            TRI_total = ((TRI*area).sum())
            TRE_total = ((TRE*area).sum())
            TDA_total = ((TDA*area).sum())
            TDL_total = ((TDL*area).sum())
            TDN_total = ((TDN*area).sum())
            TDU_total = ((TDU*area).sum())
            TDB_total = ((TDB*area).sum())
            TDS_total = ((TDS*area).sum())
            TDD_total = ((TDD*area).sum())
            TD_total = ((total_td*area).sum())
            MSD_total = ((MSD*area).sum())
            if ll == 0:
                MSD_vmr = (MSD_total/mol_mass[var]) / (mass[ss,:,:].sum()/28.96)
            elif ll == 1:
                MSD_vmr = (MSD_total/mol_mass[var]) / (mass[ss,0:26,:].sum()/28.96)
            elif ll == 2:
                MSD_vmr = (MSD_total/mol_mass[var]) / (mass[ss,26:38,:].sum()/28.96)
            elif ll == 3:
                MSD_vmr = (MSD_total/mol_mass[var]) / (mass[ss,38:58,:].sum()/28.96)
            elif ll == 4:
                MSD_vmr = (MSD_total/mol_mass[var]) / (mass[ss,58::,:].sum()/28.96)


    # write out annual chem tendency
            unit_covet = 1.e-9*365*24*3600 # kg/sec -> Tg/year
            line_ann = line_ann + '<pre> '+ format(varname[var]+layer[ll],'12s')
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDO_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDI_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TRI_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDE_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TRE_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDA_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDL_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDN_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDU_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDB_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDS_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TDD_total)*unit_covet)
            if varname[var] in WD_list:
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(WD_total)*unit_covet)
            else:
                line_ann = line_ann + '     ----------'
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(TD_total)*unit_covet)
            line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(MSD_total)*1.e-9) 
            if ll == 5:
                line_ann = line_ann + '     ----------' +'</pre>'
            else:
                line_ann = line_ann + '     '+"{0:+.3e}".format(np.array(MSD_vmr)) + '</pre>'

            line_txt = line_txt + '     '+ format(varname[var]+layer[ll],'12s')
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDO_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDI_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TRI_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDE_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TRE_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDA_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDL_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDN_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDU_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDB_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDS_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TDD_total)*unit_covet)
            if varname[var] in WD_list:
                line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(WD_total)*unit_covet)
            else:
                line_txt = line_txt + '     ----------'
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(TD_total)*unit_covet)
            line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(MSD_total)*1.e-9)
            if ll == 5:
                line_txt = line_txt + '     ----------' +'\n'
            else:
                line_txt = line_txt + '     '+"{0:+.3e}".format(np.array(MSD_vmr)) + '\n'
            
    
    fileout_ann.write(line_ann)
    fileout_txt.write(line_txt)
    fileout_ann.close()
    fileout_txt.close()
    h0_in.close()

