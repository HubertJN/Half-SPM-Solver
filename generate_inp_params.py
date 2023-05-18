import numpy as np
import pandas as pd
import scipy.stats as st
from scipy.stats import qmc
import sys

num_dat = int(sys.argv[1])
num_vars = 9

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

sampler = qmc.LatinHypercube(d=num_vars)
sample = sampler.random(n=num_dat)
dat = np.empty((num_dat, num_vars))

dat[:, 0] = temp.ppf(sample[:, 0])
dat[:, 1] = rad.ppf(sample[:, 1])
dat[:, 2] = thick.ppf(sample[:, 2])
dat[:, 3] = rr_coef.ppf(sample[:, 3])
dat[:, 4] = dif_coef.ppf(sample[:, 4])
dat[:, 5] = init_c.ppf(sample[:, 5])
dat[:, 6] = max_c.ppf(sample[:, 6])
dat[:, 7] = vol_per.ppf(sample[:, 7])
dat[:, 8] = iapp.ppf(sample[:, 8])

dat[0, :] = [temp.mean(), rad.mean(), thick.mean(), rr_coef.mean(), dif_coef.mean(), init_c.mean(), max_c.mean(), vol_per.mean(), iapp.mean()]

dat_fram = {'temp': dat[:, 0],
            'rad': dat[:, 1],
            'thick': dat[:, 2],
            'rr_coef': dat[:, 3],
            'dif_coef': dat[:, 4],
            'init_c': dat[:, 5],
            'max_c': dat[:, 6],
            'vol_per': dat[:, 7],
            'iapp': dat[:, 8],
            }
df = pd.DataFrame(dat_fram, columns=['temp', 'rad', 'thick', 'rr_coef', 'dif_coef', 'init_c', 'max_c', 'vol_per', 'iapp'])
df.to_csv('data.csv', index=False)
