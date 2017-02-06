import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm
from cmath import *
import colorsys
from scipy import signal
from netCDF4 import Dataset
import datetime
import pandas as pd

model = 'HADGEM3ES'
type = '*'

model_fname = '/work/home/db876/modelling/%s/%s_SURFACE_O3_2005_2010_*_*_*_D_%s.nc'%(model,model,type)

root_grp = Dataset(model_fname)
var = root_grp.variables['o3'][:]
var = var*1e9
date = root_grp.variables['date'][:]
time = root_grp.variables['time'][:]
lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
grid_size = root_grp.variables['grid_size'][:]
grid_size = grid_size[0]

print lat_c
print lon_c

z = var[10,:,:]

#------------------------------------------
#set up plot
fig =plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')
 
#setup basemap projection
m = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c')

m.drawcoastlines()
m.drawmapboundary()
parallels = np.arange(-90,91,15)
meridians = np.arange(-180,151,30)

pl = m.pcolor(lon_e,lat_e,z, vmin=np.min(z), vmax=np.max(z),linewidth=0.5,cmap=plt.cm.coolwarm,picker = 5)
cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')

plt.show()
