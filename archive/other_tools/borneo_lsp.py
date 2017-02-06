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
import modules
from bpch import bpch
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker

model_name = 'GEOS_CHEM'
species = 'ISOP'
model_range = '2005_2010'

#put together model file name

model_f = model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+'v90103_2x2.5_GEOS5.nc'

root_grp = Dataset(model_f)
model_var = root_grp.variables[species.lower()][:]
model_date = root_grp.variables['date'][:]
model_time = root_grp.variables['time'][:]
lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
grid_size = root_grp.variables['grid_size'][:]
grid_size = grid_size[0]
model_var_mask = np.ma.masked_where(model_var<0,model_var)

gridbox_count = len(lat_c)*len(lon_c)

lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,1.,114.)

y = model_var[:,lat_n,lon_n]

y = y*1e9

#set up plot
fig =plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1) 

x = modules.date_process(model_date,model_time,2005)
 
ofac=4 
 
model_periods,model_mag,model_ph,model_fr,model_fi,amp_corr = modules.take_lomb(x,y,ofac,1./24)

def form2(x, pos):
    """ This function returns a string with 3 decimal places, given the input x"""
    return '%.2f' % x   

def form5(x, pos):
    """ This function returns a string with 3 decimal places, given the input x"""
    return '%.6f' % x

ax.loglog(model_periods,model_mag,color='red',markersize=20) 
plt.grid(True)
xformatter = FuncFormatter(form2)
yformatter = FuncFormatter(form5)
ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
plt.ylabel('Amplitude (ppbv)',fontsize=22)
plt.xlabel('Period (Days)',fontsize=22)
plt.title('Lomb Scargle Periodogram of Borneo Isoprene (ppbv), for GEOS_CHEM v90103 2x2.5',fontsize=24,y=1.04)
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)

plt.savefig('BORNEO_LSP.png')

plt.show()




