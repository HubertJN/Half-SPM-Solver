#This script generates input files to run the pde solver given a data file of input parameters and the sample number

#Import relevant packages
import netCDF4 as NC
from netCDF4 import Dataset
import pandas as pd
import numpy as np
import sys

#Get the sample number
row_num = int(sys.argv[1])

#Open the .csv file with the input parameters and extract the correct row of samples
dat = pd.read_csv('data.csv')

Temp = dat['temp'][row_num]
Rad = dat['rad'][row_num]
Thick = dat['thick'][row_num]
Rr_coef = dat['rr_coef'][row_num]
Dif_coef = dat['dif_coef'][row_num]
Init_c = dat['init_c'][row_num]
Max_c = dat['max_c'][row_num]
Vol_per = dat['vol_per'][row_num]
Iapp = dat['iapp'][row_num]

#Read other unchanged parameters from the original input file
dat_inp = NC.Dataset("SPM_input_ori.nc", "r", format="NETCDF4")
Sim_steps = dat_inp['sim_steps'][:][0]
Dt = dat_inp['dt'][:][0]
Out_steps = dat_inp['out_steps'][:][0]
Space_steps = dat_inp['space_steps'][:][0]
Volt_do = dat_inp['volt_do'][:][0]
Checkpoint = dat_inp['checkpoint'][:][0]

#Create a netCDF input file and save all the parameters
#(since this input file is temporary it performs only mandatory commands in constructing the input file)
rootgrp = Dataset('SPM_input.nc', 'w', format='NETCDF4')

#creating dimensions for vector holding all the input variables
temp_dim = rootgrp.createDimension('temp_dim', 1)
rad_dim = rootgrp.createDimension('rad_dim', 1)
thick_dim = rootgrp.createDimension('thick_dim', 1)
rr_coef_dim = rootgrp.createDimension('rr_coef_dim', 1)
dif_coef_dim = rootgrp.createDimension('dif_coef_dim', 1)
max_c_dim = rootgrp.createDimension('max_c_dim', 1)
init_c_dim = rootgrp.createDimension('init_c_dim', 1)
iapp_dim = rootgrp.createDimension('iapp_dim', 1)
vol_per_dim = rootgrp.createDimension('vol_per_dim', 1)
sim_steps_dim = rootgrp.createDimension('sim_steps_dim', 1)
dt_dim = rootgrp.createDimension('dt_dim', 1)
out_steps_dim = rootgrp.createDimension('out_steps_dim', 1)
space_steps_dim = rootgrp.createDimension('space_steps_dim', 1)
volt_do_dim = rootgrp.createDimension('volt_do_dim', 1)
checkpoint_dim = rootgrp.createDimension('checkpoint_dim', 1)

#createing variable 
temp = rootgrp.createVariable('temp', 'f8', ('temp_dim',))
rad = rootgrp.createVariable('rad', 'f8', ('rad_dim',))
thick = rootgrp.createVariable('thick', 'f8', ('thick_dim',))
rr_coef = rootgrp.createVariable('rr_coef', 'f8', ('rr_coef_dim',))
dif_coef = rootgrp.createVariable('dif_coef', 'f8', ('dif_coef_dim',))
max_c = rootgrp.createVariable('max_c', 'f8', ('max_c_dim',))
init_c = rootgrp.createVariable('init_c', 'f8', ('init_c_dim',))
iapp = rootgrp.createVariable('iapp', 'f8', ('iapp_dim',))
vol_per = rootgrp.createVariable('vol_per', 'f8', ('vol_per_dim',))
sim_steps = rootgrp.createVariable('sim_steps', 'i4', ('sim_steps_dim',))
dt = rootgrp.createVariable('dt', 'f8', ('dt_dim',))
out_steps = rootgrp.createVariable('out_steps', 'i4', ('out_steps_dim',))
space_steps = rootgrp.createVariable('space_steps', 'i4', ('space_steps_dim',))
volt_do = rootgrp.createVariable('volt_do', 'i4', ('volt_do_dim',))
checkpoint = rootgrp.createVariable('checkpoint', 'i4', ('checkpoint_dim',))

# writing data to input_parameters variable
#Purposefully disable checkpoint since it's a feature that could be implemented for sampling...
#...but currently isnt
temp[0] = Temp
rad[0] = Rad
thick[0] = Thick
rr_coef[0] = Rr_coef
dif_coef[0] = Dif_coef
max_c[0] = Max_c
init_c[0] = Init_c
iapp[0] = Iapp
vol_per[0] = Vol_per
sim_steps[0] = Sim_steps
dt[0] = Dt
out_steps[0] = Out_steps
space_steps[0] = Space_steps
volt_do[0] = Volt_do
checkpoint[0] = 0

#closing the NetCDF file
rootgrp.close()
