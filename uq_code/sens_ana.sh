#!/bin/bash

set -e

echo "Sensitivity Analysis"
echo "Input Perturbation of parameters:"

mv ../SPM_input.nc ./SPM_input_ori.nc
mv ../SP_output.nc ./SP_output_ori.nc 2> /dev/null || true
mv ../SP_check.chp ./SP_check_ori.chp 2> /dev/null || true

echo "Getting input Parameters and generating data"
python3 generate_inp_params_sens.py

echo "Preparing database store"
if [ ! -d "data_store_sens" ] 
then
    mkdir data_store_sens
fi

echo "Generating input files"
for i in {0..9}
do
echo $((i+1))"/10"
python3 make_input_file.py $i
mv SPM_input.nc ../
(cd ../ && make -s exe)
mv ../SP_output.nc ./data_store_sens/SP_output_$i.nc
done

mv data.csv ./data_store_sens/inputs.csv

echo "Visualising Results"
python3 visual_uq_res.py

mv ./sens_data.csv ./data_store_sens/
mv ./SPM_input_ori.nc ../SPM_input.nc
mv ./SP_output_ori.nc ../SP_output.nc 2> /dev/null || true
mv ./SP_check_ori.chp ../SP_check.chp 2> /dev/null || true
