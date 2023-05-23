#!/bin/bash

set -e

echo "SENSITIVITY ANALYSIS"
echo ""

echo "Preparing database store"
if [ ! -d "data_store_sens" ] 
then
    mkdir data_store_sens
fi

mv ../SPM_input.nc ./SPM_input_ori.nc
mv ../SP_output.nc ./SP_output_ori.nc 2> /dev/null || true
mv ../SP_check.chp ./SP_check_ori.chp 2> /dev/null || true

echo "Getting input Parameters and generating data"
python3 generate_inp_params_sens.py

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
mv ./SPM_input_ori.nc ./data_store_sens/
mv ./SP_output_ori.nc ./data_store_sens/ 2> /dev/null || true
mv ./SP_check_ori.chp ./data_store_sens/ 2> /dev/null || true

echo "Visualising Results"
python3 visual_uq_res.py

mv ./sens_data.csv ./data_store_sens/
cp ./data_store_sens/SPM_input_ori.nc ../SPM_input.nc
cp ./data_store_sens/SP_output_ori.nc ../SP_output.nc 2> /dev/null || true
cp ./data_store_sens/SP_check_ori.chp ../SP_check.chp 2> /dev/null || true
