#!/bin/bash

set -o errexit
set -o pipefail

## calculate observation errors 
Rscript A21_obs_error_tmaxtmin.R ${scen} ${start} ${end}

## calculate biases between grid and obs
Rscript A22_prep_tmean_tmax_tmin.R

start2=${start}
end2=${end}
climoff=${climoff}
years=($(seq ${start2:0:4} ${end2:0:4}))
len=${#years[@]}
  
  
echo "start EnKF"
echo "start is  "${start2}
echo "scen is   "${scen}
echo "length is "${len}
echo "offset is "${climoff}
  
n=0
maxjobs=5
  
## this was just to try out when the code was not working
for (( i = 0; i < ${len}; i++ ));do
  
  echo ${years[${i}]} &
    ## R-script of EnKF - select tmin or tmax for assimilation
    Rscript A23_EnKF_tmeanstns_and_tmaxtmingrid.R ${years[${i}]} ${scen} ${date} ${start} ${end} ${L} ${Z} ${climoff} ${dist} ${window} tmin ${hist} &
    
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
  wait # wait until all have finished (not optimal, but most times good enough)
  echo $n wait
fi
  
done
  
  
  
  
  