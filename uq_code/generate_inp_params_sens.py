import numpy as np
import pandas as pd
import sys
import netCDF4 as NC

eps = 1e-6 #float(sys.argv[1])

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

#temp, rad, thick, rr_coef, dif_coef, init_c, max_c, vol_per, iapp
dat = np.array([[temp, rad, thick, rr_coef, dif_coef, init_c, max_c, vol_per, iapp]])

a = np.tile(dat, (9, 1))

for i in range(9):
    a[i,i] = a[i,i] + a[i,i]*eps
A = np.append(dat, a, axis=0)

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
