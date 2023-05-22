#!/bin/bash

set -e

echo "Uncertainty Propagation"

cp ../SPM_input.nc ./SPM_input_ori.nc
cp ../SP_output.nc ./SP_output_ori.nc
cp ../SP_check.chp ./SP_check_ori.chp

echo "Getting input parameters and generating data"
samps=$(python3 generate_inp_params.py)

echo "Preparing database store"
if [ ! -d "data_store_up" ] 
then
    mkdir data_store_up
fi
cp data.csv ./data_store_up/inputs.csv

echo "Generating Input files"
i=1
while [ $i -le $samps ]
do
    echo $i"/"$samps
    python3 make_input_file.py $((i-1))
    mv SPM_input.nc ../
    (cd ../ && make -s exe)
    mv ../SP_output.nc ./data_store_up/SP_output_$((i-1)).nc
    ((i++))
done

echo "Visualising Results"
python3 visual_up_data.py $samps

cp ./SPM_input_ori.nc ../SPM_input.nc
cp ./SP_output_ori.nc ../SP_output.nc
cp ./SP_check_ori.chp ../SP_check.chp
