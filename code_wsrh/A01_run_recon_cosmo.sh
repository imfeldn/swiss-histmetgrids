#!/bin/bash

set -o errexit
set -o pipefail

# set definitions for the reconstruction 
date=2025-04-02
start=1961-01-01
end=2020-10-31
scen=scen_hist2_wind
window=80
dist=Gower
stnval=all
calonly=T
climoff=T
L=1000
Z=500

echo "Calculate new analog days"
Rscript A04_masterscript_cosmo.R $scen $stnval $dist $start $end $calonly $climoff $date
 

echo "Extract analogs"
scene=${scen}
scen=${scen}_${window}
echo ${scen}

. A10_extract_netcdf_loop.sh ${date} ${start} ${end} ${scen} ${dist}
 
echo "Do EnKF for wind speed"
scen=${scene}
. A20_EnKF_loop.sh  ${scen} ${date} ${start} ${end} ${L} ${Z} ${climoff} ${dist} ${window} 
