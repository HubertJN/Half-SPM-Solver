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
do_volt = dat['volt_do'][:][0]
dt = dat['dt'][:][0]

#creating figure and axes dependent on whether voltage data is to be written or not
if do_volt == 1:
	fig = plt.figure(figsize=(10,8))
	ax1 = plt.subplot2grid((2,2),(0,0))
	ax2 = plt.subplot2grid((2,2),(0,1))
	ax3 = plt.subplot2grid((2,2),(1,0), colspan=2)
else:
	fig = plt.figure(figsize=(10,4))
	ax1 = plt.subplot2grid((1,2),(0,0))
	ax2 = plt.subplot2grid((1,2),(0,1))


# Subplot 1 - Animtion of lithium concentration 
# accessing the 2d concentration data and creating a time and space axis
c = dat['conc'][:]

t_steps = np.shape(c)[0]
r_steps = dat['space_steps'][:][0]
rad = dat['rad'][:][0]*10**(6)
time_axis = np.linspace(0, t_steps-1, t_steps-1).astype(int)
x_axis = np.linspace(0, rad, r_steps, endpoint=True)

# animation of the concentration data 
# the main parts of the  animation code were taken from:
# https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/

intervaltime = 0.5   #time between each animations step in ms

# plotting the first graph which is a 2D line of lithium ion concentration
# with respect to the particle radius
line, = ax1.plot(x_axis, c[0])

# function evolving the animation, which is called each frame 
# the x-axis is the particle radius, so only the x-axis data
# i.e. the lithium concantration changes with simulation time
def animate(t):
    line.set_ydata(c[t,:])
    return line,

# creating the animation using matplotlib's FunAnimation
# the arguments to be passed:
# figure of the graph  - fig from fix,ax = plt.subplots())
# function evolving the animation - animate(t)
# interval time between each timesetep - intervaltime
# frames, the array to iterate over - time_axis
# blit - blitting set to true
animation = FuncAnimation(fig, animate, interval=intervaltime, frames=time_axis, blit=True)

# customise the graph with axes labels, major and minor ticks, labels etc.
ax1.set_xlabel('Distance from particle centre [$\mu m$]', size=8)
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='x', labelsize=7)
ax1.set_ylabel('Lithium concentration $mol*m^{-3}$', size=8)
ax1.set_ylim(np.min(c)-1*10**(-9), np.max(c)+1*10**(-9))
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='y', labelsize=7)
ax1.ticklabel_format(axis='both', style="sci", useMathText=True)
ax1.xaxis.offsetText.set_fontsize(7)
ax1.yaxis.offsetText.set_fontsize(7)
ax1.set_title('Lithium Concentration Across Particle Radius', size=10, pad=15.0)

# Subplot 2 - pcolorplot of concentration data

# most of the required data is the same as in subplot 1
# get the discretisation steps of the radius
dr = rad/(r_steps-1)
 
# reversed virdis colormap, values over the range are black, under goes white
colour = cm.get_cmap('viridis_r').copy()
colour.set_under(color='w')
colour.set_over(color='k')
max_colour = np.max(c)
min_colour = np.min(c)
colour_range = max_colour - min_colour
ticklist = np.linspace(min_colour, max_colour, 6)

#create grid for plotting
#invert arrays to make sure smallest circle is on top
rs = np.linspace(0,rad,r_steps,endpoint=True)
rs = rs[::-1]
cinv = c[:,::-1]

#draw patches of circles and add them to collection
#zorder to ensure smallest circle is on top
patches = []
for i, r in enumerate(rs):
    circle = mpl.patches.Circle((0,0), r, zorder=i)
    patches.append(circle)
    
p = mpl.collections.PatchCollection(patches, cmap=colour)
p.set_clim([min_colour,max_colour])
colours = np.array(c[0,:])
p.set_array(colours)
ax2.add_collection(p)

#set figure labels
ax2.set_xlim(-rad,rad)
ax2.set_ylim(-rad,rad)
ax2.set_xlabel('Distance from particle centre [$\mu m$]', size=8)
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.tick_params(axis='x', labelsize=7)
ax2.set_ylabel('Distance from particle centre [$\mu m$]', size=8)
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.tick_params(axis='y', labelsize=7)
ax2.ticklabel_format(axis='both', style="sci", useMathText=True)
ax2.xaxis.offsetText.set_fontsize(7)
ax2.yaxis.offsetText.set_fontsize(7)
ax2.set_title('Contour Plot of Concentration inside Particle', size=10, pad=15.0)
cbar = ax2.figure.colorbar(p, ax=ax2, cmap=colour, ticks=ticklist)
cbar.ax.tick_params(labelsize=7)
cbar.ax.set_ylabel('Lithium concentration $mol*m^{-3}$', size=8)
cbar.ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.10f}"))

# animation function for the circle plot animation
#only changes the colours of the circles without redrawing them
# the concentration data for each timestep is feeded into a colour map
# the inteval and frames are the same as in subplot 1 and blitting is set to true
def animate_pcol(t):
    colours = np.array(c[t,:])
    p.set_array(colours)
    
    return p,

animation_pcol = FuncAnimation(fig, animate_pcol, interval=intervaltime, frames=time_axis, blit=True)

# Subplot 3 - plot of voltage output if do_volt set to true
if do_volt == 1:

	# accessing the voltage data and setting the number of time steps as x-axis
	volt = dat['volt'][:][:,0]
	time = np.linspace(0, t_steps, t_steps)

	# plotting the 2D graph and customising it
	ax3.plot(time*dt, volt)
	ax3.set_xlabel('Time [s]', size=8)
	ax3.xaxis.set_minor_locator(AutoMinorLocator())
	ax3.tick_params(axis='x', labelsize=7)
	ax3.set_ylabel('Voltage [V]', size=8)
	ax3.yaxis.set_minor_locator(AutoMinorLocator())
	ax3.tick_params(axis='y', labelsize=7)
	ax3.ticklabel_format(axis='both', style="sci", useMathText=True)
	ax3.xaxis.offsetText.set_fontsize(7)
	ax3.yaxis.offsetText.set_fontsize(7)
	ax3.set_title('Voltage output', size=10, pad=15.0)

	# animation of an curser that moves with the time axis with the same interval and frames
	# as the other two animations
	vl = ax3.axvline(time[0]*dt, color='black', linestyle=':')
	ax3.set_xlim()

	def animate_t_bar(t):
		vl.set_xdata(time[t]*dt)
		return vl,

	animation_t_bar = FuncAnimation(fig, animate_t_bar, interval=intervaltime, frames=time_axis, blit=True)

plt.draw()
plt.tight_layout()
plt.show()

dat.close()
