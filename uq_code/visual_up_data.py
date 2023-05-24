#This script visualises the results from an uncertainty propagation calculation

#Import relevant packages
import matplotlib.pyplot as plt
import numpy as np
import sys
import netCDF4 as nc
import pandas as pd

#Read the input parameters from the original input file
inp_dat = nc.Dataset(f"data_store_up/SPM_input_ori.nc", "r", format="NETCDF4")
num_samps = inp_dat['no_samples'][:][0]

#Get the results from the mean parameters and save to an array
dat_mu = nc.Dataset(f"data_store_up/SP_output_0.nc", "r", format="NETCDF4")
volt_dat_mu = np.array(dat_mu['volt'][:][:,0])

#Get the number of time steps and the step length in order to create a time x axis in seconds
t_steps = dat_mu['sim_steps'][:][0]
dt = dat_mu['dt'][:][0]
x = (np.linspace(0, t_steps, t_steps).astype(int))*dt

#Create an array to store all the voltage data from each run
volt_array = np.empty((t_steps, (num_samps-1)))

#Open the files from the data base, extract the voltage and save to the array
for i in range(1, num_samps, 1):
    dat = nc.Dataset(f"data_store_up/SP_output_{i}.nc", "r", format="NETCDF4")
    volt_array[:, (i-1)] = np.array(dat['volt'][:][:,0])
    dat.close()

#Calculate the 2.5th and 97.5th percentile of the voltage at each time step to get a 95% confidence
ana = np.percentile(volt_array, [2.5, 97.5], axis=1)

#Plot the range onto a graph
plt.fill_between(x, ana[0, :], ana[1, :], alpha=0.2, label='95% Confidence (random samps)', color='C0')
plt.plot(x, ana[0, :], color='C0', alpha=0.5)
plt.plot(x, ana[1, :], color='C0', alpha=0.5)

#Save the data to a data frame, and a .csv file for the user to be able to move to a different plotting software
dat_fram = {'mean volt': volt_dat_mu,
            '5th percentile': ana[0,:],
            '95th percentile': ana[1,:]}
df = pd.DataFrame(dat_fram, columns=['mean', '5th percentile', '95th percentile'])
df.to_csv('voltage_confidence_up.csv', index=False)

#Plot the mean voltage
plt.plot(x, volt_dat_mu, color='black', label='Mean Voltage')
plt.title('Uncertainty in Voltage Calculation', pad=15.0)
plt.xlabel('Time [s]')
plt.ylabel('Voltage [V]')
plt.grid()

#If the user has performed the sensitivity analysis, and calculated uncertainty using dV_dx (see generate_inp_params.py)...
#...They can plot this data on top of the other uncertainty propagation data 
if (sys.argv[1] == 'True'):
    std_V_dat = pd.read_csv('./data_store_sens/std_V_dat.csv')

    std_V = np.array(std_V_dat['std_V'][:])
    low_b = volt_dat_mu-(2*std_V)
    up_b = volt_dat_mu+(2*std_V)
    plt.fill_between(x, low_b, up_b, alpha=0.2, label='95% Confidence (SDs)', color='C3')
    plt.plot(x, low_b, alpha=0.5, color='C3')
    plt.plot(x, up_b, alpha=0.5, color='C3')

#Display the plot
plt.legend()

inp_dat.close()
dat_mu.close()

plt.show()

