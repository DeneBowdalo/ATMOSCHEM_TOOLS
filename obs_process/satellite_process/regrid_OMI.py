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

#load in model values
values = np.load('OMI_O3_trop_ave.npy')

print values.shape

#mask values array
values = np.ma.masked_where(values == -99999,values)

#lat lon edges for 1x1.25 grid
lat_e = np.arange(-60.,60.5,1)
lon_e = np.arange(-180,181,1.25)
lat_c = np.arange(-59.5,60,1)
lon_c = np.arange(-179.375,180,1.25)

lat_c_regrid = np.arange(-88.,89.,2)
lat_c_regrid = np.insert(lat_c_regrid,0,-89.5)
lat_c_regrid = np.append(lat_c_regrid,89.5)
lon_c_regrid = np.arange(-180,178,2.5)
lat_e_regrid = np.arange(-89.,90,2)
lat_e_regrid = np.insert(lat_e_regrid,0,-90.)
lat_e_regrid = np.append(lat_e_regrid,90.)
lon_e_regrid = np.arange(-181.25,179,2.5)

print lat_e_regrid.shape
print lon_e_regrid.shape
print lat_c_regrid.shape
print lon_c_regrid.shape

#regrid array

lat = np.arange(60)
lon = np.arange(144)
time = np.arange(111) 
 
new_values = np.empty((111,60,144))
stdev_array = np.empty((111,60,144))
 
lon_inds_start = np.arange(1,287,2)
#lon_inds_start = np.arange(1,287,2)
lon_inds_start = np.insert(lon_inds_start,0,287)
lon_inds_end = np.arange(0,288,2) 

lat_inds_start = np.arange(1,120,2)
lat_inds_end = np.arange(2,120,2)

lat = np.arange(59)
lon = np.arange(144)
time = np.arange(111)

print lat_inds_start.shape
print lat_inds_end.shape
print lon_inds_start.shape
print lon_inds_end.shape

print lon_inds_start

for i in lat:
    for j in lon:
        for k in time:
            print k,i,j
            val_1 = values[k,lat_inds_start[i],lon_inds_start[j]]
            val_2 = values[k,lat_inds_end[i],lon_inds_start[j]]
            val_3 = values[k,lat_inds_start[i],lon_inds_end[j]]
            val_4 = values[k,lat_inds_end[i],lon_inds_end[j]]
            print val_1,val_2,val_3,val_4
            stdev_array[k,i,j] = np.ma.std((val_1,val_2,val_3,val_4))
            new_values[k,i,j] = np.ma.average((val_1,val_2,val_3,val_4))
 
np.save('OMI_O3_trop_ave_2x2.5.npy',new_values)
np.save('OMI_O3_trop_ave_2x2.5_stdev.npy',stdev_array)


