#Script to generate a sensitivity analysis plot over time

#Import relevant packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import netCDF4 as nc
import sys

#SENSITIVITY ANALYSIS
#Assume independence between variables so that first order sensitivities are enough

#set variable means x_i
#run simulation to obtain quantity of interest Q
#get the mean Q_0
#give nominated variances to variables - (1e-4%) as var_i, sigma_i = sqrt(var_i)
#run simulation (num of variables) times with perturbed variables by h = eps*x_i
#get first order sensitivity (Q(x_i)-Q_0)/h
#obtain vector of sensitivities dQ_dx
#Calculate covariance matrix (only diagonal terms here) Sigma
#Obtain variance in Q via var_Q = dQ_dx.T @ Sigma @ dQ_dx
#Obtain uncertainty estimate via sigma_Q = sqrt(var_Q)

#Set mean of variables here
params = ['temp','rad','thick','rr_coef','dif_coef','init_c','max_c','vol_per','iapp']

#Import the original parameters from the original input file and save to a vector
dat_inp = nc.Dataset("data_store_sens/SPM_input_ori.nc", "r", format="NETCDF4")

mu = np.array([dat_inp['temp'][:][0],
               dat_inp['rad'][:][0],
               dat_inp['thick'][:][0],
               dat_inp['rr_coef'][:][0],
               dat_inp['dif_coef'][:][0],
               dat_inp['init_c'][:][0],
               dat_inp['max_c'][:][0],
               dat_inp['vol_per'][:][0],
               dat_inp['iapp'][:][0]])

sim_steps = dat_inp['sim_steps'][:][0]
dt = dat_inp['dt'][:][0]

#Auxilliary 2D array of variable means, identical columns
Mu_s = np.tile(mu.reshape((9,1)),(1,sim_steps))

#Set value of perturbation (eps) here
#Set h as a matrix for element-wise division, identical columns with elements = h_i
eps = 1e-6
h_i = (eps*mu).reshape((9,1))
h = np.tile(h_i,(1, sim_steps))

#Obtain mean voltage curve V_0 here, each row identical
#9 rows representing variables
#columns represent time steps
V_0 = np.zeros((9,sim_steps))
dat_mu = nc.Dataset("data_store_sens/SP_output_0.nc", "r", format="NETCDF4")
volt_0 = np.array(dat_mu['volt'][:][:,0])
V_0 = np.tile(volt_0, (9, 1))

#Obtain perturbed voltages here
#9 rows representing variables
#Columns representing timesteps
Vs = np.zeros((9,sim_steps))

for i in range(1, 10, 1):
    dat = nc.Dataset(f"data_store_sens/SP_output_{i}.nc", "r", format="NETCDF4")
    Vs[i-1, :] = np.array(dat['volt'][:][:,0])
    dat.close()

#Compute dV/dx as a function of time
dV_dx = np.divide((Vs - V_0),h)

#Compute absolute scaled sensitivity using element-wise multiplication
#each column i of the matrix represents the sensitivities for each paramter at time t=i
V_first_sensitivities = np.abs(Mu_s*dV_dx)

#Put the data into a data frame, and save to .csv
dat_fram = {'temp_sens': V_first_sensitivities[0, :],
            'rad_sens': V_first_sensitivities[1, :],
            'thick_sens': V_first_sensitivities[2, :],
            'rr_coef_sens': V_first_sensitivities[3, :],
            'dif_coef_sens': V_first_sensitivities[4, :],
            'init_c_sens': V_first_sensitivities[5, :],
            'max_c_sens': V_first_sensitivities[6, :],
            'vol_per_sens': V_first_sensitivities[7, :],
            'iapp_sens': V_first_sensitivities[8, :],
            'temp_dvdx': dV_dx[0, :],
            'rad_dvdx': dV_dx[1, :],
            'thick_dvdx': dV_dx[2, :],
            'rr_coef_dvdx': dV_dx[3, :],
            'dif_coef_dvdx': dV_dx[4, :],
            'init_c_dvdx': dV_dx[5, :],
            'max_c_dvdx': dV_dx[6, :],
            'vol_per_dvdx': dV_dx[7, :],
            'iapp_dvdx': dV_dx[8, :]}

cols = ['temp_sens', 'rad_sens', 'thick_sens', 'rr_coef_sens', 'dif_coef_sens', 'init_c_sens', 'max_c_sens', 'vol_per_sens', 'iapp_sens',
        'temp_dvdx', 'rad_dvdx', 'thick_dvdx', 'rr_coef_dvdx', 'dif_coef_dvdx', 'init_c_dvdx', 'max_c_dvdx', 'vol_per_dvdx', 'iapp_dvdx']
df = pd.DataFrame(dat_fram, columns=cols)
df.to_csv('sens_data.csv', index=False)


#plot as animation
fig, ax = plt.subplots(figsize=(10,6))

x = np.arange(0, 9, 1)

#time between frames
intervaltime = 0.5

#plot first frame
line, = plt.plot(x, V_first_sensitivities[:, 0], 'o-')

#definition of animation
def animate(i):
    line.set_ydata(V_first_sensitivities[:, i])
    time.set_text(('T=')+str(i*dt)+(' s'))
    return line, time,

#Design plot
ax.set_xticks(x)
ax.set_xticklabels(params, rotation=45)
ax.set_ylabel('Absolute scaled sensitivity of $V$')
ax.set_xlabel('Parameter')
ax.set_ylim(np.min(V_first_sensitivities)+10**(-9)-0.1, np.max(V_first_sensitivities)+1)
ax.grid()
ax.set_title('First Order Sensitivities Over Time')

time = ax.text(0.1,0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5}, transform=ax.transAxes)

#plot animation
animate_sensitivities = FuncAnimation(fig, animate,interval=10, frames=range(1,sim_steps), blit=True)

dat_inp.close()
dat_mu.close()

plt.show()
