import netCDF4 as NC
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import (StrMethodFormatter, AutoMinorLocator)
from matplotlib.animation import FuncAnimation

# reading in the NetCDF data from SP_output.nc
dat = NC.Dataset("SP_output.nc", "r", format="NETCDF4")

#creating figure and axes
fig = plt.figure()
ax1 = plt.subplot2grid((2,2),(0,0))
ax2 = plt.subplot2grid((2,2),(0,1))
ax3 = plt.subplot2grid((2,2),(1,0), colspan=2)

# animation of the concentration data 
# the main parts of the code were taken from:
# https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/


# Subplot 1 - Animtion of lithium concentration 
# accessing the 2d concentration data and creating a time and space axis
c = dat['conc'][:]
t_steps = dat['sim_steps'][:][0]
r_steps = dat['space_steps'][:][0]
rad = dat['rad'][:][0]*10**(6)
time_axis = np.linspace(0, t_steps-1, t_steps-1).astype(int)
x_axis = np.linspace(0, rad, r_steps, endpoint=True)

intervaltime = dat['dt'][:][0]*10**3   #time between each animations step in ms

# plotting the first graph which is a 2D line of lithium ion concentration
# with respect to the particle radius
line = ax1.plot(x_axis, c[0,:])[0]

# function evolving the animation, which is called each frame 
# the x-axis is the particle radius, so only the x-axis data
# i.e. the lithium concantration changes with simulation time
def animate(t):
    line.set_ydata(c[t,:])

# creating the animation using matplotlib's FunAnimation
# the arguments to be passed:
# figure of the graph  - fig from fix,ax = plt.subplots())
# function evolving the animation - animate(t)
# interval time between each timesetep - intervaltime
# frames, the array to iterate over - time_axis
animation = FuncAnimation(fig, animate, interval=intervaltime, frames=time_axis)

# customise the graph with axes labels, major and minor ticks, labels etc.
ax1.set_xlabel('particle radius [$\mu m$]', size=8)
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='x', labelsize=7)
ax1.set_ylabel('Lithium concentration $mol*m^{-3}$', size=8)
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='y', labelsize=7)
ax1.set_title('Lithium concentration', size=10)

# Subplot 2 - pcolorplot of concentration data

# a lot of the data is the same as in subplot 1
dr = rad/(r_steps-1)
 
#rainbow colormap, values over the range are black, under goes white
colour = cm.get_cmap('gist_rainbow').copy()
colour.set_under(color='w')
colour.set_over(color='k')
max_colour = np.max(c)
min_colour = np.min(c)
ticklist = np.linspace(c.min(), c.max(), 6)

#create finer grid for plotting
gridsize = 200
x = np.linspace(-rad,rad,gridsize,endpoint=True)
y = np.linspace(-rad,rad,gridsize,endpoint=True)

#form 2D array of radii
z = np.array([np.sqrt(i*i+j*j) for j in y for i in x])

#form concentration values; give negative concentration for colourmap if out of range - goes white
for idx, r in enumerate(z):
	if (r <= rad):
		ri = np.rint(r/dr).astype(np.int32)
		z[idx] = c[0,:][ri]
	else:
		z[idx] = -1

# reshape to form 2D array for pcolor
Z = z.reshape(gridsize,gridsize)

# form meshgrid and plot using pcolor
X, Y = np.meshgrid(x,y)
cplot = ax2.pcolor(X,Y,Z, cmap=colour, vmin=min_colour, vmax=max_colour, edgecolors='face')
ax2.set_xlabel('particle radius [$\mu m$]', size=8)
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.tick_params(axis='x', labelsize=7)
ax2.set_ylabel('particle radius [$\mu m$]', size=8)
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.tick_params(axis='y', labelsize=7)
ax2.set_title('concentration coloured contour plot', size=10)
cbar = ax2.figure.colorbar(cplot, ax=ax2, cmap=colour, ticks=ticklist)
cbar.ax.tick_params(labelsize=7)
cbar.ax.set_ylabel('Lithium concentratoin', size=8)
cbar.ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.1f}"))

# animation function for the pcolor plot animation
# works generally the same as for subplot 1
# the concentration data for each timestep is created as in the pcolor plot above
# while iterating through the time dimension of the concentratoin data
# the resulting 2D array is flattened in the end to a 1D array that can be used by FuncAnimation
# the inteval and frames are the same as in subplot 1
def animate_pcol(t):
    z =  np.array([np.sqrt(i*i+j*j) for j in y for i in x])
    for idx, r in enumerate(z):
	    if (r <= rad):
		    ri = np.rint(r/dr).astype(np.int32)
		    z[idx] = c[t,:][ri]
	    else:
		    z[idx] = -1
    Z = z.reshape(gridsize,gridsize)
    cplot.set_array(Z.flatten())
animation_pcol = FuncAnimation(fig, animate_pcol, interval=intervaltime, frames=time_axis)


# Subplot 3 - plot of voltage output

# accessing the voltage data and setting the number of time steps as x-axis
vol = dat['volt'][:][:,0]
time = np.linspace(0, t_steps, t_steps).astype(int)

# plotting the 2D graph and customising it
ax3.plot(vol, time)
ax3.set_xlabel('time [s]', size=8)
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.tick_params(axis='x', labelsize=7)
ax3.set_ylabel('Voltage [V]', size=8)
ax3.yaxis.set_minor_locator(AutoMinorLocator())
ax3.tick_params(axis='y', labelsize=7)
ax3.set_title('Voltage output', size=10)

plt.draw()
plt.tight_layout()
plt.show()