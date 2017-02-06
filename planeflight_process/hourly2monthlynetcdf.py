import datetime
import glob
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import os

species_list  = ['O3','CO','NO','NO2','ISOP']#,'C2H6','C3H8','OH','PAN','RO2','TRA_46','GMAO_ABSH','GMAO_PSFC','GMAO_TEMP','GMAO_UWND','GMAO_VWND']

model_vers = os.getcwd().split('/')[-5]
grid_res = os.getcwd().split('/')[-4]

#start_year = int(raw_input('Start Year?\n'))
#end_year = int(raw_input('\nEnd Year?\n'))

start_year = 2005
end_year = 2009

year_range = range(start_year,end_year+1)

start = datetime.datetime(start_year,1,1,0,0)
end = datetime.datetime(end_year+1,1,1,0,0)

dt_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='M')]
dt_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='M')]

n_months = len(dt_dates)

time_pd = pd.date_range(start = start,end = end, freq = 'H')[:-1]

#set grid
if grid_res == '4x5':
    lat_c = np.arange(-86.,87.,4)
    lat_c = np.insert(lat_c,0,-89)
    lat_c = np.append(lat_c,89)
    lon_c = np.arange(-180,176,5)
    lat_e = np.arange(-88.,90,4)
    lat_e = np.insert(lat_e,0,-90.)
    lat_e = np.append(lat_e,90.)
    lon_e = np.arange(-182.5,178,5)
if grid_res == '2x2.5':
    lat_c = np.arange(-88.,89.,2)
    lat_c = np.insert(lat_c,0,-89.5)
    lat_c = np.append(lat_c,89.5)
    lon_c = np.arange(-180,178,2.5)
    lat_e = np.arange(-89.,90,2)
    lat_e = np.insert(lat_e,0,-90.)
    lat_e = np.append(lat_e,90.)
    lon_e = np.arange(-181.25,179,2.5)

n_lats = len(lat_c)
n_lons = len(lon_c)
grid_size = n_lats * n_lons

key = lambda x: pd.Period(str(x.year)+'-'+str(x.month))

full_data = np.empty((n_months,n_lats,n_lons))

for species in species_list:
    print species
    
    g = Dataset('GEOSCHEM_SURFACE_%s_%s_%s_%s_%s_GEOS5_H_*.nc'%(species,start_year,end_year+1,model_vers,grid_res))
    hourly_data = g.variables[species.lower()][:]
    
    for lat_i in range(n_lats):
       for lon_i in range(n_lons): 
            print lat_i,lon_i
            var_pd = pd.Series(hourly_data[:,lat_i,lon_i], index=time_pd)
            group=var_pd.groupby(key)
            group_ave = group.mean()
            full_data[:,lat_i,lon_i] = np.array(group_ave.values.tolist())
                
            
    #get monthly average in each gridbox

    root_grp = Dataset('GEOSCHEM_SURFACE_%s_%s_%s_%s_%s_GEOS5_M_*.nc'%(species,start_year,end_year+1,model_vers,grid_res), 'w')
    root_grp.description = 'Gridded MONTHLY Surface %s in ppb - Program written by Dene Bowdalo'%(species)

    # dimensions
    #set variables
    root_grp.createDimension('%s'%(species.lower()), n_months)
    root_grp.createDimension('dates', n_months)
    root_grp.createDimension('times', n_months)
    root_grp.createDimension('lat_centre', len(lat_c))
    root_grp.createDimension('lon_centre', len(lon_c))
    root_grp.createDimension('lat_edges', len(lat_e))
    root_grp.createDimension('lon_edges', len(lon_e))
    root_grp.createDimension('grid_s', 1)
 
    #set variables
    #set variables
    species = root_grp.createVariable('%s'%(species.lower()), 'f8', ('%s'%(species.lower()),'lat_centre','lon_centre'))
    d = root_grp.createVariable('date', 'i8', ('dates',))
    t = root_grp.createVariable('time', 'i8', ('times',))
    lat_centre = root_grp.createVariable('lat_centre', 'f8', ('lat_centre',))
    lon_centre = root_grp.createVariable('lon_centre', 'f8', ('lon_centre',))
    lat_edges = root_grp.createVariable('lat_edges', 'f8', ('lat_edges',))
    lon_edges = root_grp.createVariable('lon_edges', 'f8', ('lon_edges',))
    gs = root_grp.createVariable('grid_size', str, ('grid_s',))
 
    species[:] = full_data
    d[:] = dt_dates
    t[:] = dt_times
    lat_centre[:] = lat_c
    lon_centre[:] = lon_c
    lat_edges[:] = lat_e
    lon_edges[:] = lon_e
    gs[0] = grid_res

