import datetime
import glob
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import os
from netCDF4 import num2date, date2num

#set model version and grid resolution
model_vers = 'v1001'
grid_res = '4x5'

#get start year and end year from raw input (last year is last year data is in)
start_year = int(raw_input('Start Year?\n'))
end_year = int(raw_input('\nEnd Year?\n'))
year_range = range(start_year,end_year+1)

#setup datetime array
start = datetime.datetime(start_year,1,1,0,0,0)
end = datetime.datetime(end_year+1,1,1,0,0,0)
datetime_array = pd.date_range(start,end,freq='H')[:-1].tolist()
n_hours = len(datetime_array)

#get list of planeflight logs
all_files = []
for year in year_range:
    files = glob.glob('plane.log.%s*'%(year))
    all_files=np.append(all_files,files)
all_files.sort()

#set format and read file
read = np.genfromtxt(all_files[0],names=True)

float_array = (len(read.dtype.names)-4)*',f4'
species_list = list(read.dtype.names)

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
 
big_data = np.empty((n_hours,n_lats,n_lons))
big_data_flat = np.empty(grid_size*n_hours)

#setup netcdf output file

root_grp = Dataset('GEOSCHEM_SURFACE_%s_%s_%s_%s_GEOS5_H_*.nc'%(start_year,end_year+1,model_vers,grid_res), 'w')
root_grp.description = 'Gridded Hourly Surface Species in ppb - Program written by Dene Bowdalo'

#set dimensions
root_grp.createDimension('time', n_hours)
root_grp.createDimension('lat_centre', len(lat_c))
root_grp.createDimension('lon_centre', len(lon_c))
root_grp.createDimension('lat_edges', len(lat_e))
root_grp.createDimension('lon_edges', len(lon_e))
root_grp.createDimension('grid_s', 1)

#set variables
t = root_grp.createVariable('time', 'f8', ('time',))
lat_centre = root_grp.createVariable('lat_centre', 'f4', ('lat_centre',))
lon_centre = root_grp.createVariable('lon_centre', 'f4', ('lon_centre',))
lat_edges = root_grp.createVariable('lat_edges', 'f4', ('lat_edges',))
lon_edges = root_grp.createVariable('lon_edges', 'f4', ('lon_edges',))
gs = root_grp.createVariable('grid_size', str, ('grid_s',))

#set units
t.units = 'hours since 0001-01-01 00:00:00'
t.calendar = 'gregorian'

#process times to save out (gregorian from 1 AD)
times = date2num(datetime_array, units=t.units, calendar=t.calendar)

#write out time and grid to netcdf
t[:] = times
lat_centre[:] = lat_c
lon_centre[:] = lon_c
lat_edges[:] = lat_e
lon_edges[:] = lon_e
gs[0] = grid_res

#do not read these species from planeflight logs to output
species_2_rm = ['POINT','TYPE','LAT','LON','PRESS','YYYYMMDD','HHMM','N2O5','PPN','PMN','R4N2','H2O2','MP','CH2O','MO2','ETO2','PRPE','ALK4','ACET','ALD2','MEK','RCHO','MVK',
                'SO2','DMS','SO4','REA_N2O5','REA_307','REA_323','REA_324','REA_325','REA_326','MSA','TRA_1','TRA_2','TRA_3','TRA_4','TRA_5','TRA_6','TRA_7',
                'TRA_8','TRA_9','TRA_10','TRA_11','TRA_12','TRA_13','TRA_14','TRA_15','TRA_16','TRA_17','TRA_18','TRA_19','TRA_20','TRA_21','TRA_22','TRA_23',
                'TRA_24','TRA_25','TRA_26','TRA_27','TRA_28','TRA_29','TRA_30','TRA_31','TRA_32','TRA_33','TRA_34','TRA_35','TRA_36','TRA_37','TRA_38','TRA_39',                           
                'TRA_40','TRA_41','TRA_42','TRA_43','REA_327','REA_328','REA_329','NO3','HNO4','HNO3','HNO2','HO2','TRA_44','TRA_45','TRA_46','TRA_47','TRA_48',
                'TRA_49','TRA_50','TRA_51','TRA_52','TRA_53','GMAO_ABSH','GMAO_SURF','GMAO_HFLUX'] 
for s in species_2_rm:
    if s in species_list:
        species_list.remove(s)

#create output grid 
big_data = np.empty((len(species_list),n_hours,n_lats,n_lons))

#iterate through files, reading and putting data into output grid
start=0
end=24
for num in range(len(all_files)):
    print all_files[num]
    data = np.genfromtxt(all_files[num],names=True,dtype='i10,S3,S8,S4%s'%(float_array))

    lat = 0
    lon = 0
    for i in range(grid_size):
        for j in range(len(species_list)):
            big_data[j,start:end,lat,lon] = data[species_list[j]][i::grid_size]
        if lon == (n_lons-1):
            lat+=1
            lon = 0
        else:
            lon+=1
    start+=24
    end+=24

for i in range(len(species_list)):
    current_species = species_list[i].lower()
    species = root_grp.createVariable('%s'%(current_species), 'f4', ('time','lat_centre','lon_centre'))
    species[:] = big_data[i,:,:,:]

root_grp.close()
