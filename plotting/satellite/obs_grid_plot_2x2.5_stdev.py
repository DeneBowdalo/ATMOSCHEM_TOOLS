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

#set up plot
fig =plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')

model_f = 'OMI_O3_trop_ave_2x2.5_stdev.npy'

# load in model values
values = np.load(model_f)

#lat lon edges for 2x5 grid
lat_e = np.arange(-59.,60,2)
lon_e = np.arange(-181.25,179,2.5)
lat_c = np.arange(-58.,59,2)
lon_c = np.arange(-180,178,2.5)

#get size of grid
grid_dim_c = len(lat_c)*len(lon_c)

lat = np.arange(59)
lon = np.arange(144)

ave_values = np.empty((59,144))
#average times
for i in lat:
    for j in lon:
	ave_values[i,j] = np.average(values[:,i,j])
 
#setup basemap projection
m = Basemap(projection='cyl',llcrnrlat=-59,urcrnrlat=59,\
                   llcrnrlon=lon_e[0],\
                    urcrnrlon=lon_e[-1],\
                   resolution='c')


m.drawcoastlines()
m.drawmapboundary()
parallels = np.arange(-50,51,20)
meridians = np.arange(-180,151,30)
plt.xticks(meridians)
plt.yticks(parallels)
m.drawparallels(parallels)
m.drawmeridians(meridians)

#plot model gridboxes
poly = m.pcolor(lon_e, lat_e, ave_values,vmin=0, vmax=10,cmap = plt.cm.coolwarm)
cb = plt.colorbar(poly, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
cb.set_label('Stdev', fontsize = 16)
plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)
cb.ax.tick_params(labelsize=16)
plt.title('Tropospheric O3 OMI 2x2.5 Stdevs Oct 2004 - Jan 2014 ', fontsize = 18)
plt.show()



