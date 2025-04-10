#!/bin/bash

set -o errexit
set -o pipefail

n=0
maxjobs=4

date=${date}
start=${start}
end=${end}
scen=${scen}
dist=${dist}

echo ${date}

years=($(seq ${start:0:4} ${end:0:4}))
len=${#years[@]}

conda activate /home/anaconda3/envs/my_env

for (( i = 0; i < ${len}; i++ ));do
	 	
	echo ${years[${i}]} &
	 python3 A1_extract_netcdf_relhum_cosmo.py output/ARM_run_${date}/${scen}/ ${years[${i}]} ${start} ${dist} ${scen} ${date} &
   python3 A1_extract_netcdf_minrelhum_cosmo.py output/ARM_run_${date}/${scen}/ ${years[${i}]} ${start} ${dist} ${scen} ${date} &
  	 
  if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished 
        echo $n wait
    fi
	
 done

conda deactivate
echo "done"
