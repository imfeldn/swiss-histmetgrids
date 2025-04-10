#!/bin/bash

set -o errexit
set -o pipefail

n=0
maxjobs=4

date=${date}
start=1864
end=2020
scen=${scen}
dist=${dist}

echo ${date}

years=($(seq ${start:0:4} ${end:0:4}))
len=${#years[@]}

conda activate /home/anaconda3/envs/my_env

for (( i = 0; i < ${len}; i++ ));do
	 	
	echo ${years[${i}]} &
   python3 A11_extract_netcdf_sunshineduration.py output/ARM_run_${date}/${scen}/ ${years[${i}]} ${start} ${dist} ${scen} ${date} &
  	 
  if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished 
        echo $n wait
    fi
	
 done

conda deactivate
echo "done"
