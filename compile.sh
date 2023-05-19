#!/bin/bash
#chmod u+r+x compile.sh if permission denied
f2pyfiles="pde.f90"
#rm *.o *.mod
gfortran -c -Wall $f2pyfiles -llapack 
f2py -c -m f90_src $f2pyfiles -llapack 


