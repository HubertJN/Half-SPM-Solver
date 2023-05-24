import netCDF4 as NC
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import (StrMethodFormatter, AutoMinorLocator)
from matplotlib.animation import FuncAnimation

# reading in the data
SP_sim1 = NC.Dataset("./data/SP_sim1.nc", "r", format="NETCDF4")
SP_sim2 = NC.Dataset("./data/SP_sim2.nc", "r", format="NETCDF4")
SP_sim3 = NC.Dataset("./data/SP_sim3.nc", "r", format="NETCDF4")
SP_sim4 = NC.Dataset("./data/SP_sim4.nc", "r", format="NETCDF4")
SP_sim5 = NC.Dataset("./data/SP_sim5.nc", "r", format="NETCDF4")

pybamm_sim1 = NC.Dataset("./data/PyBaMM_sim1.nc", "r", format="NETCDF4")
pybamm_sim2 = NC.Dataset("./data/PyBaMM_sim2.nc", "r", format="NETCDF4")
pybamm_sim3 = NC.Dataset(".//data/PyBaMM_sim3.nc", "r", format="NETCDF4")
pybamm_sim4 = NC.Dataset("./data/PyBaMM_sim4.nc", "r", format="NETCDF4")
pybamm_sim5 = NC.Dataset("./data/PyBaMM_sim5.nc", "r", format="NETCDF4")

SP_dis025 = NC.Dataset("./data/SP_dis025.nc", "r", format="NETCDF4")
SP_dis05 = NC.Dataset("./data/SP_dis05.nc", "r", format="NETCDF4")
SP_dis2 = NC.Dataset("./data/SP_dis2.nc", "r", format="NETCDF4")
SP_dis4 = NC.Dataset("./data/SP_dis4.nc", "r", format="NETCDF4")
SP_dis8 = NC.Dataset("./data/SP_dis8.nc", "r", format="NETCDF4")

# calculate relative absolute error and root mean square error of the 5 simulations
cp1 = pybamm_sim1['conc'][:]
cs1 = SP_sim1['conc'][:].T
r1 = pybamm_sim1['rad'][:]*10**(6)
err_sim1 = abs(cp1 - cs1)/cp1
if np.min(cp1) < np.min(cs1):
    c1_min = np.min(cp1)
else:
    c1_min = np.min(cs1)  
if np.max(cp1) < np.max(cs1):
    c1_max = np.max(cp1)
else:
    c1_max = np.max(cs1)
rmse1 = np.sqrt(np.sum((cs1-cp1)**2)/cp1.size)
print('The root mean suqaure error for charging with default values is {x:=5.3f}'.format(x=rmse1))

cp2 = pybamm_sim2['conc'][:]
cs2 = SP_sim2['conc'][:].T
err_sim2 = abs(cp2 - cs2)/cp2
r2 = pybamm_sim2['rad'][:]*10**(6)
if np.min(cp2) < np.min(cs2):
    c2_min = np.min(cp2)
else:
    c2_min = np.min(cs2)  
if np.max(cp2) < np.max(cs2):
    c2_max = np.max(cp2)
else:
    c2_max = np.max(cs2)  
rmse2 = np.sqrt(np.sum((cs2-cp2)**2)/cp2.size)
print('The root mean suqaure error for discharging with default values is {x:=5.3f}'.format(x=rmse2))

# PyBamm only evaluates every 60th timestep - reshape our own concentration data to have the same shape to get RMSE
cp3_all = pybamm_sim3['conc'][:]
r3 = pybamm_sim3['rad'][:]*10**(6)
t3, idx = np.unique(pybamm_sim3['time'][:].astype(int), return_index=True)
cp3 = cp3_all[:,idx]
cs3_all = SP_sim3['conc'][:].T
cs3 = np.zeros_like(cp3)
for i in range (cs3.shape[0]):
    for j in range (cs3.shape[1]-1):
        cs3[i,j] = cs3_all[i,j*60]
cs3[:,-1] = cs3_all[:,-1]
err_sim3 = abs(cp3 - cs3)/cp3
if np.min(cp3) < np.min(cs3):
    c3_min = np.min(cp3)
else:
    c3_min = np.min(cs3)  
if np.max(cp1) < np.max(cs3):
    c3_max = np.max(cp3)
else:
    c3_max = np.max(cs3)
rmse3 = np.sqrt(np.sum((cs3-cp3)**2)/cp3.size)
print('The root mean suqaure error when running a GITT cycle is {x:=5.3f}'.format(x=rmse3))


cp4 = pybamm_sim4['conc'][:]
cs4 = SP_sim4['conc'][:].T
err_sim4 = abs(cp4 - cs4)/cp4
r4 = pybamm_sim4['rad'][:]*10**(6)
if np.min(cp4) < np.min(cs4):
    c4_min = np.min(cp4)
else:
    c4_min = np.min(cs4)  
if np.max(cp4) < np.max(cs4):
    c4_max = np.max(cp4)
else:
    c4_max = np.max(cs4)  
rmse4 = np.sqrt(np.sum((cs4-cp4)**2)/cp4.size)
print('The root mean suqaure error when increasing the diffusion coefficient is {x:=5.3f}'.format(x=rmse4))

cp5 = pybamm_sim5['conc'][:]
cs5 = SP_sim5['conc'][:].T
err_sim5 = abs(cp5 - cs5)/cp5
r5 = pybamm_sim5['rad'][:]*10**(6)
if np.min(cp5) < np.min(cs5):
    c5_min = np.min(cp5)
else:
    c5_min = np.min(cs5)  
if np.max(cp5) < np.max(cs5):
    c5_max = np.max(cp5)
else:
    c5_max = np.max(cs5)
rmse5 = np.sqrt(np.sum((cs5-cp5)**2)/cp5.size)
print('The root mean suqaure error when increasing the radius is {x:=5.3f}'.format(x=rmse5))

# plot with subplots showing the evolution of concentration over time from both PyBaMM and our
# results on the left and the evolution of the relative absolute error on the right
fig, axs = plt.subplots(5,2, figsize=(10,12))

t_steps = np.shape(cp1)[1]
time_axis = np.linspace(0, t_steps-1, t_steps-1).astype(int)
intervaltime = 10

# setting up the first line for each subplot
conc_cp1, = axs[0,0].plot(r1, cp1[:,0], label='PyBaMM')
conc_cs1, = axs[0,0].plot(r1, cs1[:,0], label='own SPM model')
err1, = axs[0,1].plot(r1, err_sim1[:,0])

conc_cp2, = axs[1,0].plot(r2, cp2[:,0], label='PyBaMM')
conc_cs2, = axs[1,0].plot(r2, cs2[:,0], label='own SPM model')
err2, = axs[1,1].plot(r2, err_sim2[:,0])

conc_cp4, = axs[2,0].plot(r4, cp4[:,0], label='PyBaMM')
conc_cs4, = axs[2,0].plot(r4, cs4[:,0], label='own SPM model')
err4, = axs[2,1].plot(r4, err_sim4[:,0])

conc_cp5, = axs[3,0].plot(r5, cp5[:,0], label='PyBaMM')
conc_cs5, = axs[3,0].plot(r5, cs5[:,0], label='own SPM model')
err5, = axs[3,1].plot(r5, err_sim5[:,0])

conc_cp3, = axs[4,0].plot(r3, cp3[:,0], label='PyBaMM')
conc_cs3, = axs[4,0].plot(r3, cs3[:,0], label='own SPM model')
err3, = axs[4,1].plot(r3, err_sim3[:,0])

# animtion of each subplots
def animate(i):
    conc_cp1.set_ydata(cp1[:,i])
    conc_cs1.set_ydata(cs1[:,i])
    err1.set_ydata(err_sim1[:,i])
    conc_cp2.set_ydata(cp2[:,i])
    conc_cs2.set_ydata(cs2[:,i])
    err2.set_ydata(err_sim2[:,i])
    conc_cp4.set_ydata(cp4[:,i])
    conc_cs4.set_ydata(cs4[:,i])
    err4.set_ydata(err_sim4[:,i])
    conc_cp5.set_ydata(cp5[:,i])
    conc_cs5.set_ydata(cs5[:,i])
    err5.set_ydata(err_sim5[:,i])
    return conc_cp1, conc_cs1, err1, conc_cp2, conc_cs2, err2, conc_cp4, conc_cs4, err4, conc_cp5, conc_cs5, err5,
animation = FuncAnimation(fig, animate, interval=intervaltime, frames=time_axis, blit=True)

# separate animation for the GITT experiment since the timeframe is different to the other animations
def animate_GITT(i):
    conc_cp3.set_ydata(cp3[:,i])
    conc_cs3.set_ydata(cs3[:,i])
    err3.set_ydata(err_sim3[:,i])
    return conc_cp3, conc_cs3, err3,
animation_GITT = FuncAnimation(fig, animate_GITT, interval=intervaltime, frames=np.arange(0,len(t3),1), blit=True)

#customisign the various subplots
axs[0,0].set_xlabel('Distance from particle centre [$\mu m$]', size=8)
axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].tick_params(axis='x', labelsize=7)
axs[0,0].set_ylabel('Lithium concentration $mol*m^{-3}$', size=8)
axs[0,0].set_ylim(c1_min-1000, c1_max+1000)
axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].tick_params(axis='y', labelsize=7)
axs[0,0].ticklabel_format(axis='both', style="sci", useMathText=True)
axs[0,0].xaxis.offsetText.set_fontsize(7)
axs[0,0].yaxis.offsetText.set_fontsize(7)
axs[0,0].set_title('Lithium Concentration during charging', size=10, pad=15.0)
axs[0,0].legend(loc='lower left')

axs[0,1].set_xlabel('Distance from particle centre [$\mu m$]', size=8)
axs[0,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,1].tick_params(axis='x', labelsize=7)
axs[0,1].set_ylabel('Absolute error', size=8)
axs[0,1].set_ylim(np.min(err_sim1)-0.05, np.max(err_sim1)+0.05)
axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,1].tick_params(axis='y', labelsize=7)
axs[0,1].ticklabel_format(axis='both', style="sci", useMathText=True)
axs[0,1].xaxis.offsetText.set_fontsize(7)
axs[0,1].yaxis.offsetText.set_fontsize(7)
axs[0,1].set_title('Relative absolute error during charging', size=10, pad=15.0)
axs[0,1].text(0, np.max(err_sim1)+0.02, r"RMSE = {x:=5.3f} $mol*m^{{-3}}$".format(x=rmse1), fontsize=8)

axs[1,0].set_xlabel('Distance from particle centre [$\mu m$]', size=8)
axs[1,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[1,0].tick_params(axis='x', labelsize=7)
axs[1,0].set_ylabel('Lithium concentration $mol*m^{-3}$', size=8)
axs[1,0].set_ylim(c2_min-1000, c2_max+1000)
axs[1,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[1,0].tick_params(axis='y', labelsize=7)
axs[1,0].ticklabel_format(axis='both', style="sci", useMathText=True)
axs[1,0].xaxis.offsetText.set_fontsize(7)
axs[1,0].yaxis.offsetText.set_fontsize(7)
axs[1,0].set_title('Lithium Concentration during discharging', size=10, pad=15.0)
axs[1,0].legend(loc='lower left')

axs[1,1].set_xlabel('Distance from particle centre [$\mu m$]', size=8)
axs[1,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[1,1].tick_params(axis='x', labelsize=7)
axs[1,1].set_ylabel('Absolute error', size=8)
axs[1,1].set_ylim(np.min(err_sim2)-0.05, np.max(err_sim2)+0.05)
axs[1,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[1,1].tick_params(axis='y', labelsize=7)
axs[1,1].ticklabel_format(axis='both', style="sci", useMathText=True)
axs[1,1].xaxis.offsetText.set_fontsize(7)
axs[1,1].yaxis.offsetText.set_fontsize(7)
axs[1,1].set_title('Relatove absolute error during discharging', size=10, pad=15.0)
axs[1,1].text(0, np.max(err_sim2)+0.02, r"RMSE = {x:=5.3f} $mol*m^{{-3}}$".format(x=rmse2), fontsize=8)

axs[2,0].set_xlabel('Distance from particle centre [$\mu m$]', size=8)
axs[2,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[2,0].tick_params(axis='x', labelsize=7)
axs[2,0].set_ylabel('Lithium concentration $mol*m^{-3}$', size=8)
axs[2,0].set_ylim(c4_min-8000, c4_max+1000)
axs[2,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[2,0].tick_params(axis='y', labelsize=7)
axs[2,0].ticklabel_format(axis='both', style="sci", useMathText=True)
axs[2,0].xaxis.offsetText.set_fontsize(7)
axs[2,0].yaxis.offsetText.set_fontsize(7)
axs[2,0].set_title('Lithium Concentration with larger diffusion coefficient', size=10, pad=15.0)
axs[2,0].legend(loc='lower left')

axs[2,1].set_xlabel('Distance from particle centre [$\mu m$]', size=8)
axs[2,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[2,1].tick_params(axis='x', labelsize=7)
axs[2,1].set_ylabel('Absolute error', size=8)
axs[2,1].set_ylim(np.min(err_sim4)-0.05, np.max(err_sim4)+0.05)
axs[2,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[2,1].tick_params(axis='y', labelsize=7)
axs[2,1].ticklabel_format(axis='both', style="sci", useMathText=True)
axs[2,1].xaxis.offsetText.set_fontsize(7)
axs[2,1].yaxis.offsetText.set_fontsize(7)
axs[2,1].set_title('Relative absolute error with larger diffusion coefficient', size=10, pad=15.0)
axs[2,1].text(0, np.max(err_sim4)+0.03, r'RMSE = {x:=5.3f} $mol*m^{{-3}}$'.format(x=rmse4), fontsize=8)

axs[3,0].set_xlabel('Distance from particle centre [$\mu m$]', size=8)
axs[3,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[3,0].tick_params(axis='x', labelsize=7)
axs[3,0].set_ylabel('Lithium concentration $mol*m^{-3}$', size=8)
axs[3,0].set_ylim(c5_min-1000, c5_max+1000)
axs[3,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[3,0].tick_params(axis='y', labelsize=7)
axs[3,0].ticklabel_format(axis='both', style="sci", useMathText=True)
axs[3,0].xaxis.offsetText.set_fontsize(7)
axs[3,0].yaxis.offsetText.set_fontsize(7)
axs[3,0].set_title('Lithium Concentration with larger diffusion coefficient', size=10, pad=15.0)
axs[3,0].legend(loc='lower left')

axs[3,1].set_xlabel('Distance from particle centre [$\mu m$]', size=8)
axs[3,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[3,1].tick_params(axis='x', labelsize=7)
axs[3,1].set_ylabel('Absolute error', size=8)
axs[3,1].set_ylim(np.min(err_sim5)-0.05, np.max(err_sim5)+0.05)
axs[3,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[3,1].tick_params(axis='y', labelsize=7)
axs[3,1].ticklabel_format(axis='both', style="sci", useMathText=True)
axs[3,1].xaxis.offsetText.set_fontsize(7)
axs[3,1].yaxis.offsetText.set_fontsize(7)
axs[3,1].set_title('Relative absolute error with larger radius', size=10, pad=15.0)
axs[3,1].text(0, np.max(err_sim5)+0.01, r'RMSE = {x:=5.3f} $mol*m^{{-3}}$'.format(x=rmse5), fontsize=8)

axs[4,0].set_xlabel('Distance from particle centre [$\mu m$]', size=8)
axs[4,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[4,0].tick_params(axis='x', labelsize=7)
axs[4,0].set_ylabel('Lithium concentration $mol*m^{-3}$', size=8)
axs[4,0].set_ylim(c3_min-2000, c3_max+10)
axs[4,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[4,0].tick_params(axis='y', labelsize=7)
axs[4,0].ticklabel_format(axis='both', style="sci", useMathText=True)
axs[4,0].xaxis.offsetText.set_fontsize(7)
axs[4,0].yaxis.offsetText.set_fontsize(7)
axs[4,0].set_title('Lithium Concentration during GITT cycle', size=10, pad=15.0)
axs[4,0].legend(loc='lower left')

axs[4,1].set_xlabel('Distance from particle centre [$\mu m$]', size=8)
axs[4,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[4,1].tick_params(axis='x', labelsize=7)
axs[4,1].set_ylabel('Absolute error', size=8)
axs[4,1].set_ylim(np.min(err_sim3)-0.05, np.max(err_sim3)+0.05)
axs[4,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[4,1].tick_params(axis='y', labelsize=7)
axs[4,1].ticklabel_format(axis='both', style="sci", useMathText=True)
axs[4,1].xaxis.offsetText.set_fontsize(7)
axs[4,1].yaxis.offsetText.set_fontsize(7)
axs[4,1].set_title('Relative absolute error during GITT cycle', size=10, pad=15.0)
axs[4,1].text(0, np.max(err_sim3)+0.02, r"RMSE = {x:=5.3f} $mol*m^{{-3}}$".format(x=rmse3), fontsize=8)

plt.tight_layout()
plt.show()
plt.close()

# calculating root mean square error for discretisation convergence
c_dis025_all = SP_dis025['conc'][:].T
c_dis025 = np.zeros_like(cp1)
for i in range (c_dis025.shape[0]):
    for j in range (c_dis025.shape[1]-1):
        c_dis025[i,j] = c_dis025_all[i,j*4]
c_dis025[:,-1] = c_dis025_all[:,-1]
rmse_dis025 = np.sqrt(np.sum((c_dis025-cp1)**2)/cp1.size).astype(float)
print(rmse_dis025)

c_dis05_all = SP_dis05['conc'][:].T
c_dis05 = np.zeros_like(cp1)
for i in range (c_dis05.shape[0]):
    for j in range (c_dis05.shape[1]-1):
        c_dis05[i,j] = c_dis05_all[i,j*2]
c_dis05[:,-1] = c_dis05_all[:,-1]
rmse_dis05 = np.sqrt(np.sum((c_dis05-cp1)**2)/cp1.size).astype(float)
print(rmse_dis05)

rmse_dis1 = np.sqrt(np.sum((cs1-cp1)**2)/cp1.size).astype(float)
print(rmse_dis1)

c_dis2 = SP_dis2['conc'][:].T
cp1_dis2 = np.zeros_like(c_dis2)
for i in range (cp1_dis2.shape[0]):
    for j in range (cp1_dis2.shape[1]-1):
        cp1_dis2[i,j] = cp1[i,j*2]
cp1_dis2[:,-1] = cp1[:,-1]
rmse_dis2 = np.sqrt(np.sum((c_dis2-cp1_dis2)**2)/cp1_dis2.size).astype(float)
print(rmse_dis2)

c_dis4 = SP_dis4['conc'][:].T
cp1_dis4 = np.zeros_like(c_dis4)
for i in range (cp1_dis4.shape[0]):
    for j in range (cp1_dis4.shape[1]-1):
        cp1_dis4[i,j] = cp1[i,j*4]
cp1_dis4[:,-1] = cp1[:,-1]
rmse_dis4 = np.sqrt(np.sum((c_dis4-cp1_dis4)**2)/cp1_dis4.size).astype(float)
print(rmse_dis4)

c_dis8 = SP_dis8['conc'][:].T
cp1_dis8 = np.zeros_like(c_dis8)
for i in range (cp1_dis8.shape[0]):
    for j in range (cp1_dis8.shape[1]-1):
        cp1_dis8[i,j] = cp1[i,j*8]
cp1_dis8[:,-1] = cp1[:,-1]
rmse_dis8 = np.sqrt(np.sum((c_dis8-cp1_dis8)**2)/cp1_dis8.size).astype(float)
print(rmse_dis8)

x = np.array([0.25, 0.5, 1.0, 2.0, 4.0 ,8.0])
y = np.array([rmse_dis025, rmse_dis05, rmse_dis1, rmse_dis2, rmse_dis4, rmse_dis8])

figure, ax = plt.subplots()
ax.plot(x,y)

plt.show()
plt.close()

SP_sim1.close()
SP_sim2.close()
SP_sim3.close()
SP_sim4.close()
SP_sim5.close()
pybamm_sim1.close()
pybamm_sim2.close()
pybamm_sim3.close()
pybamm_sim4.close()
pybamm_sim5.close()
SP_dis025.close()
SP_dis05.close()
SP_dis2.close()
SP_dis4.close()
SP_dis8.close()
