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
#give nominated variances to variables - 1%? as var_i, sigma_i = sqrt(var_i)
#run simulation (num of variables) times with perturbed variables by h = eps*x_i
#get first order sensitivity (Q(x_i)-Q_0)/h
#obtain vector of sensitivities dQ_dx
#Calculate covariance matrix (only diagonal terms here) Sigma
#Obtain variance in Q via var_Q = dQ_dx.T @ Sigma @ dQ_dx
#Obtain uncertainty estimate via sigma_Q = sqrt(var_Q)

#Set mean of variables here
#as well as initial dataframe
#mu should be a 9x1 array
params = ['temp','rad','thick','rr_coef','dif_coef','init_c','max_c','vol_per','iapp']
dat_inp = nc.Dataset("SPM_input_ori.nc", "r", format="NETCDF4")
#mu = np.array([294.15, 5.22e-6, 75.6e-6, 3.42, 1.48e-15, 47023.326, 51765.0, 66.5, 48.685491723466406])
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
volts_df = pd.DataFrame(dict(parameter=params,mu=mu))

#Auxilliary 2D array of variable means, identical columns
Mu_s = np.tile(mu.reshape((9,1)),(1,sim_steps))

#Set value of perturbation here
#h_i = eps*mu_i
#Set h as a matrix for element-wise division, identical columns with elements = h_i
eps = float(sys.argv[1]) #1e-6
h_i = (eps*mu).reshape((9,1))
h = np.tile(h_i,(1, sim_steps))

#Obtain mean voltage curve V_0 here, each row identical
#9 rows representing variables
#1000 columns representing timesteps
V_0 = np.zeros((9,sim_steps))
dat_mu = nc.Dataset("data_store_sens/SP_output_0.nc", "r", format="NETCDF4")
volt_0 = np.array(dat_mu['volt'][:][:,0])
V_0 = np.tile(volt_0, (9, 1))

#Obtain perturbed voltages here
#9 rows representing variables
#1000 columns representing timesteps
Vs = np.zeros((9,sim_steps))

for i in range(1, 10, 1):
    dat = nc.Dataset(f"data_store_sens/SP_output_{i}.nc", "r", format="NETCDF4")
    Vs[i-1, :] = np.array(dat['volt'][:][:,0])

#Compute dV/dx as a function of time
dV_dx = np.divide((Vs - V_0),h)

#Compute absolute scaled sensitivity using element-wise multiplication
V_first_sensitivities = np.abs(Mu_s*dV_dx)

#each column i of the matrix represents the sensitivities for each paramter at time t=i
#so now set covariance matrix
#here is a placeholder where the standard deviations are 1%
std_devs = 0.01*mu
Sigma = np.zeros((len(mu),len(mu)))
for i in range(len(std_devs)):
    Sigma[i,i] = std_devs[i]**2
    
#Now calculate covariance matrix of V in time
#The diagonal entries are the actual variance of V at each time
#Off-diagonal entries are correlations - not important here
var_V = dV_dx.T @ Sigma @ dV_dx

#extract standard deviations
std_V = np.zeros(sim_steps)
for i in range(sim_steps):
    std_V[i] = np.sqrt(var_V[i,i])


#plot as animation
fig, ax = plt.subplots(figsize=(10,6))

Y = []
for i in range (sim_steps):
    V_df = volts_df.copy()
    V_df['scaled_sensitivity'] = V_first_sensitivities[:, i]
    y = V_df['scaled_sensitivity'].abs().to_numpy()
    Y.append(y)

x = np.arange(0, len(V_df.index), 1)

#time between frames
intervaltime = 0.5

#plot first frame
line, = plt.plot(x, Y[0], 'o-')

#definition of animation
def animate(i):
    line.set_ydata(Y[i])
    time.set_text(('T=')+str(i))
    return line, time,

ax.set_xticks(x)
ax.set_xticklabels(volts_df.parameter, rotation=45)
ax.set_ylabel('Absolute scaled sensitivity of $V$')
ax.set_xlabel('Parameter')
ax.set_ylim(np.min(Y)+10**(-9), np.max(Y)+1)
ax.set_yscale('log')
ax.grid()
ax.set_title('First Order Sensitivities Over Time')

time = ax.text(0.1,0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5}, transform=ax.transAxes)

#plot animation
animate_sensitivities = FuncAnimation(fig, animate,interval=10, frames=range(1,sim_steps), blit=True)
#animate_sensitivities.save('sensitivities.mp4')
plt.show()
