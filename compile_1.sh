#!/bin/bash
#chmod u+r+x compile.sh if permission denied

python3 -m venv venv
source venv/bin/activate
vim requirements.txt - names of mods
pip install requirements.txt 
pip install pandas
pip install matplotlib
pip freeze list >> requirements.txt - add all existing modules to requirements
#deactivate



f2pyfiles="datafitPde.f90"
#rm *.o *.mod
gfortran -c -W $f2pyfiles -llapack 
f2py -c -m f90_src $f2pyfiles -llapack #--f2cmap .kind_map


