#!/bin/bash
#Script to perform sensitivity analysis for the pde solver

#Step when errors occur
set -e

echo "SENSITIVITY ANALYSIS"
echo ""

#If the data base directory doesnt exist, then create it
echo "Preparing database store"
if [ ! -d "data_store_sens" ] 
then
    mkdir data_store_sens
fi

#Move the original input, output, and checkpoint files into the uq_code directory...
#...to clear main directory of files that will be overwritten
mv ../SPM_input.nc ./SPM_input_ori.nc
mv ../SP_output.nc ./SP_output_ori.nc 2> /dev/null || true
mv ../SP_check.chp ./SP_check_ori.chp 2> /dev/null || true

#Call the generation script to get a .csv of input parameters
echo "Getting input Parameters and generating data"
python3 generate_inp_params_sens.py

#Create a loop that:
#1-Generates an input file using the input parameters in the ith row of the .csv file
#2-Moves the input file from uq_code to the main directory
#3-Run the main program
#4-Moves the generated output file into the data store directory
echo "Generating input files"
for i in {0..9}
do
echo $((i+1))"/10"
python3 make_input_file.py $i
mv SPM_input.nc ../
(cd ../ && make -s exe)
mv ../SP_output.nc ./data_store_sens/SP_output_$i.nc
done

#Move the original files, and the sample .csv file, into the data store directory
mv data.csv ./data_store_sens/inputs.csv
mv ./SPM_input_ori.nc ./data_store_sens/
mv ./SP_output_ori.nc ./data_store_sens/ 2> /dev/null || true
mv ./SP_check_ori.chp ./data_store_sens/ 2> /dev/null || true

#Call the visualisation script
echo "Visualising Results"
python3 visual_uq_res.py

#Move the generated data into the data store file...
#...and move the original data files back into the main directory
mv ./sens_data.csv ./data_store_sens/
cp ./data_store_sens/SPM_input_ori.nc ../SPM_input.nc
cp ./data_store_sens/SP_output_ori.nc ../SP_output.nc 2> /dev/null || true
cp ./data_store_sens/SP_check_ori.chp ../SP_check.chp 2> /dev/null || true
