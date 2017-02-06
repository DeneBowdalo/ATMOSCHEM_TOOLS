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

#read in o3 data
model_f = '/work/home/db876/plotting_tools/model_files/GEOS_CHEM_SURFACE_O3_2005_2010_v90103_2x2.5_GEOS5.nc'
root_grp = Dataset(model_f)
o3_data = root_grp.variables['o3'][:]
lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
n_boxes = len(lat_c)*len(lon_c)

variables = ['CO','NO','NO2','GMAO_PSFC','GMAO_TEMP','WIND_DIRECTION','WIND_SPEED']

#save out corr to netcdf
root_grp_grid = Dataset('corr_grids.nc', 'w')
root_grp_grid.description = 'Pearsons Correlation stat for other species and O3 - Program written by Dene Bowdalo'

root_grp_grid.createDimension('lat_centre', len(lat_c))
root_grp_grid.createDimension('lon_centre', len(lon_c))
root_grp_grid.createDimension('lat_edges', len(lat_e))
root_grp_grid.createDimension('lon_edges', len(lon_e))

co_array = root_grp_grid.createVariable('co', 'f8', ('lat_centre','lon_centre'))
no_array = root_grp_grid.createVariable('no', 'f8', ('lat_centre','lon_centre'))
no2_array = root_grp_grid.createVariable('no2', 'f8', ('lat_centre','lon_centre'))
psfc_array = root_grp_grid.createVariable('gmao_psfc', 'f8', ('lat_centre','lon_centre'))
temp_array = root_grp_grid.createVariable('gmao_temp', 'f8', ('lat_centre','lon_centre'))
wd_array = root_grp_grid.createVariable('wind_direction', 'f8', ('lat_centre','lon_centre'))
ws_array = root_grp_grid.createVariable('wind_speed', 'f8', ('lat_centre','lon_centre'))
lat_centre = root_grp_grid.createVariable('lat_centre', 'f8', ('lat_centre',))
lon_centre = root_grp_grid.createVariable('lon_centre', 'f8', ('lon_centre',))
lat_edge = root_grp_grid.createVariable('lat_edges', 'f8', ('lat_edges',))
lon_edge = root_grp_grid.createVariable('lon_edges', 'f8', ('lon_edges',))

lat_centre[:] = lat_c
lon_centre[:] = lon_c
lat_edge[:] = lat_e     
lon_edge[:] = lon_e 

netcdf_list = [co_array,no_array,no2_array,psfc_array,temp_array,wd_array,ws_array]


#read in other species data
for i in range(len(variables)):
    species = variables[i]
    netcdf_var = netcdf_list[i]
    
    print species
    
    corr_big = np.empty((len(lat_c),len(lon_c)))

    model_f = '/work/home/db876/plotting_tools/model_files/GEOS_CHEM_SURFACE_'+species+'_2005_2010_v90103_2x2.5_GEOS5.nc'
    root_grp = Dataset(model_f)
    other_data = root_grp.variables[species.lower()][:]

    lon_i = 0
    lat_i = 0
    
    for siten in range(n_boxes):
        o3_cut = o3_data[:,lat_i,lon_i]
        other_cut = other_data[:,lat_i,lon_i]
        
        #normalise data
        o3_cut = (o3_cut - np.ma.average(o3_cut)) / (np.ma.max(o3_cut) - np.ma.min(o3_cut))
        other_cut = (other_cut - np.ma.average(other_cut)) / (np.ma.max(other_cut) - np.ma.min(other_cut))
        
        #do pearsons correlation on whole grid
        #corr = pearsonr(other_cut,o3_cut)
        #corr = corr[0]
        print len(o3_cut)
        corr = np.correlate(other_cut,o3_cut,mode='full')
        print len(corr)
        
        lag_ind = np.argmax(corr)
        corr = corr[lag_ind]        
        offsets = -lags*distancePerLag
        
        
        corr_big[lat_i,lon_i] = corr
        
        if lon_i == (len(lon_c)-1):
            lat_i+=1
            lon_i=0
        else:
            lon_i+=1
    
    #save out grid to netcdf
    netcdf_var[:] = corr_big
            
