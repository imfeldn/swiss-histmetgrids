#!/bin/bash

set -o errexit
set -o pipefail

#### calculate input for EnKF separately can be run "manually"" ####
# calculate observation errors 
Rscript A21_obs_error.R 

# calculate bias between grid and observation for pressure
Rscript A22_EnKF_help_pressure.R 

# calculate bias between grid and observation for temp
Rscript A23_EnKF_help_temp.R

# calculate covariance matrix based on weather types
Rscript A24_calcPH.R


start2=${start} ## change if you dont want to run all years
end2=${end} ## change if you dont want to run all years
climoff=${climoff}
years=($(seq ${start2:0:4} ${end2:0:4}))
len=${#years[@]}

echo "start EnKF"
echo "start is  "${start}
echo "scen is   "${scen}
echo "length is "${len}

n=0
maxjobs=4


for (( i = 0; i < ${len}; i++ ));do
	 	
	echo ${years[${i}]} &
	## R-script of EnKF
	Rscript A25_EnKF_uvwind_all_pblend.R ${years[${i}]} ${scen} ${date} ${start} ${end} ${L} ${Z} ${climoff} ${dist} U ${window} &
	Rscript A25_EnKF_uvwind_all_pblend.R ${years[${i}]} ${scen} ${date} ${start} ${end} ${L} ${Z} ${climoff} ${dist} V ${window} &
  
  if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
	
 done



