#!/bin/bash

set -e

echo "Sensitivity Analysis"
echo "Input Perturbation of parameters:"
echo "(If you are unsure, 1e-6 is suitable for the default parameters)"
read eps

#echo "Do you want to visualise when complete?"
#echo "y/n"
#read visual

cp ../SPM_input.nc ./SPM_input_ori.nc

echo "Getting input Parameters and generating data"
python3 generate_inp_params_sens.py $eps

echo "Preparing database store"
if [ ! -d "data_store_sens" ] 
then
    mkdir data_store_sens
fi

echo Generating input file
for i in {0..9}
do
echo $i"/9"
python3 make_input_file.py $i
mv SPM_input.nc ../
(cd ../ && make -s exe)
mv ../SP_output.nc ./data_store_sens/SP_output_$i.nc
#mv SP_output.nc ./uq_code/data_store_sens/SP_output_$i.nc
#cd uq_code
done

echo "Visualising Results"
python3 visual_uq_res.py $eps

cp SPM_input_ori.nc ../SPM_input.nc
