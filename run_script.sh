#!/bin/bash

num_dat=9
echo Generating Samples
python3 generate_inp_params_sens.py
#python3 generate_inp_params.py $num_dat

echo Generating input file
for i in {0..9}
do
echo $i
python3 make_input_file.py $i
make exe
mv ./SP_output.nc ./data_store/SP_output_$i.nc
done

