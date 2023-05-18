#!/bin/bash

num_dat=8
echo Generating Samples
python3 generate_inp_params_sens.py $num_dat

echo Generating input file
for i in {0..8}
do
echo $i
python3 make_input_file.py $i
make exe
mv ./SP_output.nc ./data_store/SP_output_$i.nc
done

