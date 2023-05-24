#!/bin/bash
#Script to perform uncertainty propagation on the pde solver

#Stop when errors occur
set -e

echo "UNCERTAINTY PROPAGATION"
echo ""

#Move the original data files into the uq_code directory
mv ../SPM_input.nc ./SPM_input_ori.nc
mv ../SP_output.nc ./SP_output_ori.nc 2> /dev/null || true
mv ../SP_check.chp ./SP_check_ori.chp 2> /dev/null || true

#Generate the input parameters and save the number of samples to the samps variable
#$1 tells whether the sensitivity analysis data is available...
#...(see generate_inp_params.py for more details)
echo "Getting input parameters and generating data"
samps=$(python3 generate_inp_params.py $1)

#Create the data storage directory, if it doesn't already exist
echo "Preparing database store"
if [ ! -d "data_store_up" ] 
then
    mkdir data_store_up
fi

#This loop will:
#1-make an input file using the ith sample in the .csv file
#2-move the input file to the main directory
#3-run the main program
#4-move the output file back to data storage directory
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

#Move the sample data and original input, output, checkpoint files into the data store directory
mv ./data.csv ./data_store_up/inputs.csv
mv ./SPM_input_ori.nc ./data_store_up/
mv ./SP_output_ori.nc ./data_store_up/ 2> /dev/null || true
mv ./SP_check_ori.chp ./data_store_up/ 2> /dev/null || true

#Call the visualisation code, using $1 to know if sensitivity analysis data is present
echo "Visualising Results"
python3 visual_up_data.py $1

#Move the visualisation data and sensitivity analysis uncertainty data into the data store directory...
#...and move the original input input, output, and checkpoint file back to the main directory
mv ./voltage_confidence_up.csv ./data_store_up/
mv ./std_V_dat.csv ./data_store_up/ 2> /dev/null || true
cp ./data_store_up/SPM_input_ori.nc ../SPM_input.nc
cp ./data_store_up/SP_output_ori.nc ../SP_output.nc 2> /dev/null || true
cp ./data_store_up/SP_check_ori.chp ../SP_check.chp 2> /dev/null || true
