import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from cmath import *
from math import *
import csv
import datetime
import lomb_phase
from scipy import signal
import multiprocessing
import datetime
import time
import modules
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, shiftgrid
from agpy import AG_image_tools

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.
start_year = 2005

species_list = ['O3','CO','NO','NO2']
#species_list = ['C2H6','C3H8','ISOP']
#species_list = ['NO2']


for spec in species_list:
    #read model netcdf file
    model_root_grp = Dataset('/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_%s_2005_2010_v90103_2x2.5_GEOS5_H_*.nc'%(spec))

    model_data = model_root_grp.variables[spec.lower()][:]*1e9
    model_date = model_root_grp.variables['date'][:]
    model_time = model_root_grp.variables['time'][:]
    lat_c = model_root_grp.variables['lat_centre'][:]
    lon_c = model_root_grp.variables['lon_centre'][:]
    lat_e = model_root_grp.variables['lat_edges'][:]
    lon_e = model_root_grp.variables['lon_edges'][:]

    n_boxes = len(lat_c)*len(lon_c)

    model_data = np.average(model_data,axis=0)
    print model_data

    out = np.fft.fft2(model_data)    
    out = np.fft.fftshift(out)#[-46:,-72:]
    out[:46,:72] = 0                                                                                                                                                                                                                          
    out[47:,73:] = 0
    #out[:,:76] = 0                                                                                                                                                                                                                          
    #out[:,77:] = 0   

    #cut out column
    #out[1:,:] = 0

    fr = out.real
    fi = out.imag
    ph = np.arctan2(fi,fr)
    mag = np.abs(out)
    #remove dc component 
    #out[0, 0] = 0
    #psd1D = AG_image_tools.radialprofile.azimuthalAverage(mag) 
    #plt.semilogy(psd1D,label=spec)


    #fig, (ax) = plt.subplots(1,figsize=(23,12))
    #fig.patch.set_facecolor('white')

    #cb = plt.colorbar(all,shrink=0.8,orientation = 'horizontal', format='%.2f',pad = 0.05)

    #cb.ax.tick_params(labelsize=21)

    #plt.title(spec,fontsize=21)
    #plt.show()

    #take ifft2
    #all = plt.imshow(mag,interpolation='nearest')
    #plt.show()
    
    out = np.fft.ifftshift(out)
    out2 = np.abs(np.fft.ifft2(out))
    out2 = np.flipud(out2)    
 
    print out2
    all = plt.imshow(out2,interpolation='nearest')
    plt.show()

end = datetime.datetime.now() 
diff = end - start
print 'time = ', diff.seconds
