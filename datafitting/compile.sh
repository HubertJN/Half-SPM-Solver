#!/bin/bash

#make sure vevnv is activated and requirements are installed!

f2pyfiles="datafitPde.f90"
#rm *.o *.mod
gfortran -c -W $f2pyfiles -llapack 
f2py -c -m f90_src $f2pyfiles -llapack #--f2cmap .kind_map


