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
tsDir={{ output }}/post/atm/{{ grid }}/clim/{{ '%dyr' % (ypf) }}
mkdir -p climo
#cd climo
ln -s ${tsDir}/*.nc ./climo
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
cat > e3sm_chem_index.py << EOF
#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from calendar import monthrange
import pandas as pd

pathout = './'

indexfile = open(pathout+'index.html',"w")
index = '<h><b> E3SM chemistry diagnostics package (ChemDyg) </b></h>'
index = index + '<pre> Test: ${short} </pre>'
index = index + '<pre> Reference: Observations and Reanalysis </pre>'
index = index + '<pre>   </pre>'
index = index + '<h><b> Pressure-Latitude plots </b></h>'
index = index + '<pre> <a href="O3_pres_lat_plot.png">O3</a>  <a href="O3_pres_lat_plot_mb.png">O3 (Trop)</a>  <a href="Q_pres_lat_plot.png">Q</a>  <a href="T_pres_lat_plot.png">T</a> </pre>'
index = index + '<h><b> Latitude-Longitude plots </b></h>'
index = index + '<pre> <a href="TMQ_lat_lon_plot.png">TMQ</a> </pre>'
index = index + '<h><b> NOx aircraft and lightning emission plots </b></h>'
index = index + '<pre> <a href="nox_emis_plot.png">NOx emission</a> </pre>'
index = index + '<h><b> Tropospheric O3 burden comparison with CMIP6 </b></h>'
index = index + '<pre> <a href="CMIP_tropO3_comparison.png">plot</a> </pre>'
index = index + '<h><b> Surface ozone diurnal cycle comparison with CMIP6 </b></h>'
index = index + '<pre> <a href="surfO3_hour_cycle_DJF.png">DJF</a>   <a href="surfO3_hour_cycle_JJA.png">JJA</a> </pre>'
index = index + '<h><b> Surface ozone daily maximum 8-hour average (MDA8) annual cycle comparison with CMIP6 </b></h>'
index = index + '<pre> <a href="MDA8_month_cycle.png">plot</a> </pre>'
index = index + '<h><b> Surface CO comparison with NOAA observations </b></h>'
index = index + '<pre> <a href="NOAA_CO_BRW.png"> BRW </a>  <a href="NOAA_CO_CGO.png"> CGO </a>  <a href="NOAA_CO_ICE.png"> ICE </a> <a href="NOAA_CO_KUM.png"> KUM </a>  <a href="NOAA_CO_MHD.png"> MHD </a>  <a href="NOAA_CO_MID.png"> MID </a>  <a href="NOAA_CO_PSA.png"> PSA </a>  <a href="NOAA_CO_RPB.png"> RPB </a>  <a href="NOAA_CO_SMO.png"> SMO </a>  <a href="NOAA_CO_SYO.png"> SYO </a>  <a href="NOAA_CO_WIS.png"> WIS </a>  <a href="NOAA_CO_ZEP.png"> ZEP </a>  </pre>'
index = index + '<h><b> Ozone hole </b></h>'
index = index + '<pre> <a href="to3areas_O3hole.png">areas</a>  <a href="to3mins_O3hole.png">mins</a> </pre>'
index = index + '<h><b> Ozone hole with equivalent latitude </b></h>'
index = index + '<pre> <a href="TOZ_PDF_timeseries.png">time_series</a>  <a href="TOZ_PDF_climo.png">climate</a> </pre>'
index = index + '<h><b> Ozone STE flux </b></h>'
index = index + '<pre> <a href="O3_STE_flux.png">time series</a> <a href="O3_STE_flux_monthly.png">monthly</a> <a href="STE_lat_lon_12month.png">lat-lon plot</a> </pre>'
index = index + '<h><b> QBO diagnostics </b></h>'
index = index + '<pre> <a href="QBO_fig1.png">Lat-time plot for TCO anom vs QBO phase</a> </pre>'
index = index + '<pre>Total column ozone (TCO) anomaly (DU, relative to 1979-2020 mean) as a function of QBO phase for (left) Multi-Sensor Reanalysis version 2, \n(right) E3SMv2 ObsQBO simulation. 0 is centered on the month when QBO transits from QBOE to QBOW. </pre>'
index = index + '<pre> <a href="QBO_fig2.png">Pres-time plot for O3 anom vs QBO phase</a> </pre>'
index = index + '<pre>Pressure-time plot of the anomalous ozone concentration (ppm, relative to 1985-2020 mean) as a function of QBO phase for (left row) CMZM, \n(right row) E3SMv2 ObsQBO simulation. 0 is centered on the month when QBO transits from QBOE to QBOW. </pre>'
index = index + '<pre> <a href="QBO_fig3.png">Pres-time plot for T/O3 anom vs QBO phase</a> </pre>'
index = index + '<pre>Pressure-time anomalous (left row) temperature (K), (middle row)  steady state ozone (SSO, ppm), and (right row) WSTAR (m/s) for E3SMv2. \n0 is centered on the month when QBO transits from QBOE to QBOW. SSO is derived using a linearized ozone chemistry model (Linoz). WSTAR is \nthe transformed Eulerian mean variable for vertical transport. </pre>'
index = index + '<h><b> Temperature with equivalent latitude </b></h>'
index = index + '<pre> <a href="temp_PDF_climo.png">plot</a> </pre>'
index = index + '<h><b> Chemistry tendency table </b></h>'
index = index + '<pre> <a href="chem_clim_ANN.html">ANN</a><a href="chem_clim_ANN.txt">(txt)</a>   <a href="chem_clim_DJF.html">DJF</a><a href="chem_clim_DJF.txt">(txt)</a>   <a href="chem_clim_MAM.html">MAM</a><a href="chem_clim_MAM.txt">(txt)</a>   <a href="chem_clim_JJA.html">JJA</a><a href="chem_clim_JJA.txt">(txt)</a>   <a href="chem_clim_SON.html">SON</a><a href="chem_clim_SON.txt">(txt)</a>   </pre>'
index = index +'<pre> <a href="chem_budget.html">closure check and burden</a><a href="chem_budget.txt">(txt)</a> </pre>'
index = index + '<h><b> Chemistry production/loss tendency table </b></h>'
index = index + '<pre> '+' <a href="chem_prodloss_ANN.html"> ANN </a>   <a href="chem_prodloss_DJF.html"> DJF </a>   <a href="chem_prodloss_MAM.html"> MAM </a>   <a href="chem_prodloss_JJA.html"> JJA </a>   <a href="chem_prodloss_SON.html"> SON </a>   </pre>'
index = index + '<h><b> Chem high-level summary table </b></h>'
index = index + '<pre> <a href="chem_summary_table.html">table</a> </pre>'

indexfile.write(index)
indexfile.close()

EOF

# Run diagnostics
command="python -u e3sm_chem_index.py"
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
cp index.html ${f}

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

