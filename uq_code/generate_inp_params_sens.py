#Script to generate input parameters to perform sensitivity analysis of voltage with respect to input parameters

#Import relevant packages
import numpy as np
import pandas as pd
import sys
import netCDF4 as NC

#Set the perturbation percentage to 1e-4%
eps = 1e-6

#Import initial parameters
dat_inp = NC.Dataset("SPM_input_ori.nc", "r", format="NETCDF4")

temp = dat_inp['temp'][:][0]
rad = dat_inp['rad'][:][0]
thick = dat_inp['thick'][:][0]
rr_coef = dat_inp['rr_coef'][:][0]
dif_coef = dat_inp['dif_coef'][:][0]
init_c = dat_inp['init_c'][:][0]
max_c = dat_inp['max_c'][:][0]
vol_per = dat_inp['vol_per'][:][0]
iapp = dat_inp['iapp'][:][0]

#Create a row vector of data to store the mean values
dat = np.array([[temp, rad, thick, rr_coef, dif_coef, init_c, max_c, vol_per, iapp]])

#Create an array to save the inputs with each variable perturbed individually
a = np.tile(dat, (9, 1))

#Perturb the diagonal elements to create a matrix of input parameters
for i in range(9):
    a[i,i] = a[i,i] + a[i,i]*eps

#Attach the mean values to the top of the matrix
A = np.append(dat, a, axis=0)

#Save the data to a data frame and export to .csv
dat_fram = {'temp': A[:, 0],
            'rad': A[:, 1],
            'thick': A[:, 2],
            'rr_coef': A[:, 3],
            'dif_coef': A[:, 4],
            'init_c': A[:, 5],
            'max_c': A[:, 6],
            'vol_per': A[:, 7],
            'iapp': A[:, 8],
            }
df = pd.DataFrame(dat_fram, columns=['temp', 'rad', 'thick', 'rr_coef', 'dif_coef', 'init_c', 'max_c', 'vol_per', 'iapp'])
df.to_csv('data.csv', index=False)
