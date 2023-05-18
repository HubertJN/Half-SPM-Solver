import numpy as np
import pandas as pd

eps = 1e-6

#temp, rad, thick, rr_coef, dif_coef, init_c, max_c, vol_per, iapp
dat = np.array([[294.15, 5.22e-6, 75.6e-6, 3.42, 1.48e-15, 47023.326, 51765.0, 66.5, 48.685491723466406]])

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
