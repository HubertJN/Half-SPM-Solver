import numpy as np
import pandas as pd
import scipy.stats as st
from scipy.stats import qmc
import sys

eps = 1e-6
ent = 1
temp = st.norm(294.15, ent)
rad  = st.norm(5.22e-6, ent)
thick = st.norm(75.6e-6, ent)
rr_coef = st.norm(3.42, ent)
dif_coef = st.norm(1.48e-15, ent)
init_c = st.norm(47023.326, ent)
max_c = st.norm(51765.0, ent)
vol_per = st.norm(66.5, ent)
iapp = st.norm(48.685491723466406, ent)

dat = [temp.mean(), rad.mean(), thick.mean(), rr_coef.mean(), dif_coef.mean(), init_c.mean(), max_c.mean(), vol_per.mean(), iapp.mean()]

perturbs = np.ones(9)
perturbs_vars_dat = np.ones((1, 9))

a = np.tile(dat, (9, 1))

for i in range(9):
    a[i,i] = a[i,i] + a[i,i]*eps

dat_fram = {'temp': a[:, 0],
            'rad': a[:, 1],
            'thick': a[:, 2],
            'rr_coef': a[:, 3],
            'dif_coef': a[:, 4],
            'init_c': a[:, 5],
            'max_c': a[:, 6],
            'vol_per': a[:, 7],
            'iapp': a[:, 8],
            }
df = pd.DataFrame(dat_fram, columns=['temp', 'rad', 'thick', 'rr_coef', 'dif_coef', 'init_c', 'max_c', 'vol_per', 'iapp'])
df.to_csv('data.csv', index=False)
