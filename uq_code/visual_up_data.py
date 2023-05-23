import matplotlib.pyplot as plt
import numpy as np
import sys
import netCDF4 as nc
import pandas as pd

num_samps = int(sys.argv[1])

dat_mu = nc.Dataset(f"data_store_up/SP_output_0.nc", "r", format="NETCDF4")
volt_dat_mu = np.array(dat_mu['volt'][:][:,0])
t_steps = dat_mu['sim_steps'][:][0]
x = time_axis = np.linspace(0, t_steps, t_steps).astype(int)

volt_array = np.empty((t_steps, num_samps))

for i in range(1, num_samps-1, 1):
    dat = nc.Dataset(f"data_store_up/SP_output_{i}.nc", "r", format="NETCDF4")
    volt_array[:, i] = np.array(dat['volt'][:][:,0])
    
dat = nc.Dataset(f"data_store_up/SP_output_{num_samps-1}.nc", "r", format="NETCDF4")
volt_array[:, (num_samps-1)] = np.array(dat['volt'][:][:,0])

ana = np.percentile(volt_array, [2.5, 97.5], axis=1)
plt.fill_between(x, ana[0, :], ana[1, :], alpha=0.2, label='95% Confidence (random samps)', color='C0')
plt.plot(x, ana[0, :], color='C0', alpha=0.5)
plt.plot(x, ana[1, :], color='C0', alpha=0.5)

dat_fram = {'mean volt': volt_dat_mu,
            '5th percentile': ana[0,:],
            '95th percentile': ana[1,:]}

df = pd.DataFrame(dat_fram, columns=['mean', '5th percentile', '95th percentile'])
df.to_csv('voltage_confidence_up.csv', index=False)


plt.plot(x, volt_dat_mu, color='black', label='Mean Voltage')
plt.title('Uncertainty in Voltage Calculation', pad=15.0)
plt.xlabel('Step Number')
plt.ylabel('Voltage [V]')
plt.grid()

if (sys.argv[2] == 'True'):
    std_V_dat = pd.read_csv('./std_V_dat.csv')

    std_V = np.array(std_V_dat['std_V'][:])
    low_b = volt_dat_mu-(2*std_V)
    up_b = volt_dat_mu+(2*std_V)
    plt.fill_between(x, low_b, up_b, alpha=0.2, label='95% Confidence (SDs)', color='C3')
    plt.plot(x, low_b, alpha=0.5, color='C3')
    plt.plot(x, up_b, alpha=0.5, color='C3')

plt.legend()
plt.show()

