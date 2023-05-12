import netCDF4 as NC
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

#import data here
#some assigning variables here

#dummy arrays for plotting purposes
size = 50
rad = 2
dr = rad/(size-1)
r_dummy = np.linspace(0,rad,size,endpoint=True)
c_dummy = np.linspace(1,5,size,endpoint=True)

#rainbow colormap, values over the range are black, under goes white
colour = cm.get_cmap('gist_rainbow').copy()
colour.set_under(color='w')
colour.set_over(color='k')
max_colour = np.max(c_dummy)
min_colour = np.min(c_dummy)
#range_colour = max_colour - min_colour

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
		z[idx] = c_dummy[ri]
	else:
		z[idx] = -1

#reshape to form 2D array for pcolor
Z = z.reshape(gridsize,gridsize)

#form meshgrid and plot using pcolor
X, Y = np.meshgrid(x,y)
fig, ax = plt.subplots()
cplot = ax.pcolor(X,Y,Z, cmap=colour, vmin=min_colour, vmax=max_colour, edgecolors='face')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('concentration coloured contour plot')
fig.colorbar(cplot, ax=ax)


plt.show()
