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
from scipy.stats import pearsonr

model_f = 'corr_grid_pearson.nc'
root_grp = Dataset(model_f)
co_data = root_grp.variables['co'][:]
no_data = root_grp.variables['no'][:]
no2_data = root_grp.variables['no2'][:]
psfc_data = root_grp.variables['gmao_psfc'][:]
temp_data = root_grp.variables['gmao_temp'][:]
ws_data = root_grp.variables['wind_speed'][:]
wd_data = root_grp.variables['wind_direction'][:]

lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
n_boxes = len(lat_c)*len(lon_c)

var = raw_input('CO, NO, NO2, PSFC, TEMP, WIND_SPEED or WIND_DIR?\n')

if var == 'CO':
    z = co_data
if var == 'NO':
    z = no_data
if var == 'NO2':
    z = no2_data
if var == 'PSFC':
    z = psfc_data
if var == 'TEMP':
    z = temp_data
if var == 'WIND_SPEED':
    z = ws_data
if var == 'WIND_DIR':
    z = wd_data

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
plt.xticks(meridians)
plt.yticks(parallels)
m.drawparallels(parallels)
m.drawmeridians(meridians)

pl = m.pcolor(lon_e,lat_e,z, vmin=-1, vmax=1,linewidth=0.5,cmap=plt.cm.coolwarm,picker = 5)
cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
cb.set_label("Pearson's Correlation Statistic", fontsize = 16)   

plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)
plt.title("Pearson's Correlation between Surface O3 and Surface %s"%(var),fontsize=20)


cb.ax.tick_params(labelsize=16)

plt.show()