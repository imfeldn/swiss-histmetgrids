#!/bin/bash

set -o errexit
set -o pipefail

date=2025-02-18
start=1763-01-01
end=2020-12-31
scen=all_hist
window=45
dist=Gower
stnval=all
calonly=T
climoff=T
L=1000
Z=500
hist=hist1

echo "Calculate new analog days"
Rscript A04_masterscript_txtnsd.R $scen $stnval $dist $start $end $calonly $climoff $date


echo "Extract analogs"
scene=${scen}
scen=${scen}_${window}
echo ${scen}

. A10_extract_netcdf_loop.sh ${date} ${start} ${end} ${scen} ${dist}

echo "Do EnKF for wind speed"
scen=${scene}
. A20_EnKF_tmaxtmin_loop.sh  ${scen} ${date} ${start} ${end} ${L} ${Z} ${climoff} ${dist} ${window} ${hist}
