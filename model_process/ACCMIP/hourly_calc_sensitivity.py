import modules
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

data_type = raw_input('amp, ph or ave?\n')
if data_type != 'ave':
    period = raw_input('\nd or s?\n')
    if period == 'd':
        min_ph = 0
        max_ph = 24
    elif period == 's':
        min_ph = 0
        max_ph = 12
else:
    period = ''

typ = raw_input('\nabs, pc or ave?\n')


data = Dataset('ACCMIP_hourly_%s%s_std.nc'%(period,data_type))
ave = data.variables['average'][:]
abs_std = data.variables['absolute_std'][:]
frac_std = data.variables['fractional_std'][:]
lat_c = data.variables['lat_centre'][:]
lon_c = data.variables['lon_centre'][:]
lat_e = data.variables['lat_edges'][:]
lon_e = data.variables['lon_edges'][:]

if typ == 'ave':
    z = ave
if typ == 'abs':
    z = abs_std
elif typ == 'pc':
    z = frac_std

#-------------------------------------------
#set up plot
fig =plt.figure(figsize=(18,13))
fig.patch.set_facecolor('white')
 
#setup basemap projection
m = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c')

m.drawcoastlines()
m.drawmapboundary()
parallels = np.arange(-90,91,15)
meridians = np.arange(-180,151,30)

linear_lats = []
linear_lons=[]
lat_i = 0
lon_i = 0

n_boxes = len(lat_c) * len(lon_c)

for i in range(n_boxes):
    current_lat = lat_c[lat_i]
    current_lon = lon_c[lon_i]

    linear_lats.append(current_lat)
    linear_lons.append(current_lon)
    
    if lon_i == (len(lon_c)-1):
        lat_i+=1
        lon_i=0
    else:
        lon_i+=1

if (data_type == 'ph') & (typ == 'ave'): 
    pl = m.pcolor(lon_e,lat_e,z, vmin=min_ph, vmax=max_ph,linewidth=0.5,cmap=plt.cm.hsv,picker = 5)
else:
    pl = m.pcolor(lon_e,lat_e,z, vmin=np.min(z), vmax=np.max(z),linewidth=0.5,cmap=plt.cm.rainbow,picker = 5)
cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')   

plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)

cb.ax.tick_params(labelsize=16)

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500+100")

plt.show()
