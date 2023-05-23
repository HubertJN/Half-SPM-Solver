#!/bin/bash

set -e

echo "UNCERTAINTY PROPAGATION"
echo ""

mv ../SPM_input.nc ./SPM_input_ori.nc
mv ../SP_output.nc ./SP_output_ori.nc 2> /dev/null || true
mv ../SP_check.chp ./SP_check_ori.chp 2> /dev/null || true

echo "Getting input parameters and generating data"
samps=$(python3 generate_inp_params.py $1)
echo $samps
echo "Preparing database store"
if [ ! -d "data_store_up" ] 
then
    mkdir data_store_up
fi

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

mv ./data.csv ./data_store_up/inputs.csv
mv ./SPM_input_ori.nc ./data_store_up/
mv ./SP_output_ori.nc ./data_store_up/ 2> /dev/null || true
mv ./SP_check_ori.chp ./data_store_up/ 2> /dev/null || true

echo "Visualising Results"
python3 visual_up_data.py $1

mv ./voltage_confidence_up.csv ./data_store_up/
mv ./std_V_dat.csv ./data_store_up/ 2> /dev/null || true
cp ./data_store_up/SPM_input_ori.nc ../SPM_input.nc
cp ./data_store_up/SP_output_ori.nc ../SP_output.nc 2> /dev/null || true
cp ./data_store_up/SP_check_ori.chp ../SP_check.chp 2> /dev/null || true
