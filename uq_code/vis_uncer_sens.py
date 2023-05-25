#Script to visualise the uncertainty from the sensitivity analysis data

#Import relevant packages
import matplotlib.pyplot as plt
import numpy as np
import sys
import netCDF4 as nc
import pandas as pd

#Get the results from the mean parameters and save to an array
dat_mu = nc.Dataset(f"data_store_sens/SP_output_0.nc", "r", format="NETCDF4")
volt_dat_mu = np.array(dat_mu['volt'][:][:,0])

#Get the number of time steps and the step length in order to create a time x axis in seconds
t_steps = dat_mu['sim_steps'][:][0]
dt = dat_mu['dt'][:][0]
x = (np.linspace(0, t_steps, t_steps).astype(int))*dt

#Plot the mean voltage
plt.plot(x, volt_dat_mu, color='black', label='Mean Voltage')
plt.title('Uncertainty in Voltage Calculation', pad=15.0)
plt.xlabel('Time [s]')
plt.ylabel('Voltage [V]')
plt.grid()

#Import the standard deviations calculated earlier from sensitivity analysis
std_V_dat = pd.read_csv('./data_store_sens/std_V_dat.csv')

#Save data to an array, create 95% confidence range and plot
std_V = np.array(std_V_dat['std_V'][:])
low_b = volt_dat_mu-(2*std_V)
up_b = volt_dat_mu+(2*std_V)
plt.fill_between(x, low_b, up_b, alpha=0.2, label='95% Confidence (SDs)', color='C3')
plt.plot(x, low_b, alpha=0.5, color='C3')
plt.plot(x, up_b, alpha=0.5, color='C3')

#Display the plot
plt.legend()

#Close file
dat_mu.close()

#Display
plt.show()
