import netCDF4 as NC
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import (StrMethodFormatter, AutoMinorLocator)
from matplotlib.animation import FuncAnimation
import os
cwd = os.getcwd()
print(cwd)

# reading in the data
SP_sim1 = NC.Dataset("./benchmarking/data/SP_sim1.nc", "r", format="NETCDF4")
SP_sim2 = NC.Dataset("SP_sim2.nc", "r", format="NETCDF4")
SP_sim3 = NC.Dataset("SP_sim3.nc", "r", format="NETCDF4")
SP_sim4 = NC.Dataset("SP_sim4.nc", "r", format="NETCDF4")
SP_sim5 = NC.Dataset("SP_sim5.nc", "r", format="NETCDF4")

pybamm_sim1 = NC.Dataset("PyBaMM_sim1.nc", "r", format="NETCDF4")
pybamm_sim2 = NC.Dataset("PyBaMM_sim2.nc", "r", format="NETCDF4")
pybamm_sim3 = NC.Dataset("PyBaMM_sim3.nc", "r", format="NETCDF4")
pybamm_sim4 = NC.Dataset("PyBaMM_sim4.nc", "r", format="NETCDF4")
pybamm_sim5 = NC.Dataset("PyBaMM_sim5.nc", "r", format="NETCDF4")

SP_dis025 = NC.Dataset("SP_dis025.nc", "r", format="NETCDF4")
SP_dis05 = NC.Dataset("SP_dis05.nc", "r", format="NETCDF4")
SP_dis1 = NC.Dataset("SP_dis1.nc", "r", format="NETCDF4")
SP_dis4 = NC.Dataset("SP_dis4.nc", "r", format="NETCDF4")
SP_dis6 = NC.Dataset("SP_dis6.nc", "r", format="NETCDF4")
SP_dis8 = NC.Dataset("SP_dis8.nc", "r", format="NETCDF4")

