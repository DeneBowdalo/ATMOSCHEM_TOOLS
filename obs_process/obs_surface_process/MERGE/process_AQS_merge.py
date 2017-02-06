import glob
import numpy as np
from netCDF4 import Dataset
import os
import datetime
import pandas as pd
from netCDF4 import num2date, date2num
import modules

#data range
start_year = 1970
end_year = 2014  

fname_species = os.getcwd().split('/')[-2]
if fname_species[:3] == 'NO2':
    no2_type = fname_species[4:]
    species = 'NO2'
else:
    no2_type = 'na'
    species = fname_species
    
#OUTPUT RESOLUTION CAN BE HOURLY (H), HOURLY & DAILY (HD), HOURLY, DAILY AND MONTHLY (HDM)
output_res = 'HDM'

#create grid time array
start = datetime.datetime(year = start_year, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = end_year+1, month = 1, day = 1, hour = 0, minute = 0) 
grid_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]
grid_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='H')][:-1]

n_points = len(grid_dates)

#make output dates and times
dt = pd.date_range(start,end,freq='H')[:-1].tolist()
output_res_times = date2num(dt, units='hours since 0001-01-01 00:00:00', calendar='gregorian')
obs_time_pd = pd.date_range(start = start,end = end, freq = 'H')[:-1]

#setup netcdf output
root_grp = Dataset('AQS_SURFACE_%s_%s_%s_%s.nc'%(fname_species,start_year,end_year+1,output_res), 'w')
root_grp.description = 'Hourly Surface %s at AQS sites - Program written by Dene Bowdalo'%(fname_species)

# dimensions
root_grp.createDimension('flex', n_points)
root_grp.createDimension('flex4', None)

#save time out
times = root_grp.createVariable('time', 'f8', ('flex',))
times[:] = output_res_times

#set counts at 0 
alli1 = 0
alli2 = 0
alli3 = 0
alli4 = 0
alli5 = 0

alln1 = 0
alln2 = 0
alln3 = 0
alln4 = 0
alln5 = 0
alln6 = 0
alln7 = 0
alln8 = 0

#set empty lists
unknown_mm = []
unknown_mm_refs = []
unknown_tz_refs = []
exit_nometa_refs = []
exit_nometa_lats = []
exit_nometa_lons = []
exit_nometa_pg = []
exit_anyvaliddata_refs = []
exit_anyvaliddata_lats = []
exit_anyvaliddata_lons = []
exit_anyvaliddata_pg = []
exit_nokeymeta_refs = []
exit_nokeymeta_lats = []
exit_nokeymeta_lons = []
exit_nokeymeta_pg = []
exit_resolution_refs = []
exit_resolution_lats = []
exit_resolution_lons = []
exit_resolution_pg = []
exit_badmeasurementmethod_refs = []
exit_badmeasurementmethod_lats = []
exit_badmeasurementmethod_lons = []
exit_badmeasurementmethod_pg = []

#read in a chunk at a time and write data to joint netcdf
chunk_files = glob.glob('AQS_split/AQS**%s.nc'%(output_res))

print chunk_files

for cf in chunk_files:
    chunk_read = Dataset(cf)
    valid_refs = [str(i) for i in chunk_read.groups]

    for site_ref in valid_refs:

        chunk_group = chunk_read.groups[site_ref]

        data = chunk_group.variables[species.lower()][:]
        all_mm = chunk_group.variables['measurement method'][:]
        n_dup = chunk_group.variables['n_duplicate_array'][:]
        site_name = chunk_group.site_name
        country = chunk_group.country
        contact = chunk_group.data_contact
        lat = chunk_group.latitude
        lon = chunk_group.longitude
        alt = chunk_group.altitude
        data_tz = chunk_group.data_timezone
        local_tz = chunk_group.local_timezone
        raw_class_name = chunk_group.raw_site_class
        unit = chunk_group.raw_units
        p_unit = chunk_group.processed_units 
        file_res = chunk_group.native_resolution
        flask_flag = chunk_group.flask_flag
    
        if flask_flag == 'Yes':
            all_st = [-1]*(len(data)-1)
            all_st.append(3)
        else:
            all_st = [-1]*(len(data))

        meta = [lat,lon,alt,raw_class_name,file_res,unit,p_unit,data_tz,local_tz,site_name,country,contact]

        modules.write_out_data(site_ref,'EPA AQS',root_grp,species,data,all_st,all_mm,True,meta,n_dup)

    #add counts from each file
    alli1 = alli1 + chunk_read.variables['invalid_nometa_count'][0]
    alli2 = alli2 + chunk_read.variables['invalid_anyvaliddata_count'][0]
    alli3 = alli3 + chunk_read.variables['invalid_nokeymeta_count'][0]
    alli4 = alli4 + chunk_read.variables['invalid_resolution_count'][0]
    alli5 = alli5 + chunk_read.variables['invalid_badmeasurementmethod_count'][0]

    alln1 = alln1 + chunk_read.variables['n_obs_all'][0]
    alln2 = alln2 + chunk_read.variables['n_obs_after_nometa'][0]
    alln3 = alln3 + chunk_read.variables['n_obs_after_flagsandlod'][0]
    alln4 = alln4 + chunk_read.variables['n_obs_after_duplicate'][0]
    alln5 = alln5 + chunk_read.variables['n_obs_after_anyvaliddata'][0]
    alln6 = alln6 + chunk_read.variables['n_obs_after_nokeymeta'][0]
    alln7 = alln7 + chunk_read.variables['n_obs_after_resolution'][0]
    alln8 = alln8 + chunk_read.variables['n_obs_after_badmeasurementmethod'][0]
    
    #append exit refs etc...
    unknown_mm = np.append(unknown_mm,chunk_read.variables['UNKNOWN_MM'][:])
    unknown_mm_refs = np.append(unknown_mm_refs,chunk_read.variables['UNKNOWN_MM_REFS'][:])
    unknown_tz_refs = np.append(unknown_tz_refs,chunk_read.variables['UNKNOWN_TZ_REFS'][:])
    exit_nometa_refs = np.append(exit_nometa_refs,chunk_read.variables['exit_nometa_refs'][:])
    exit_nometa_lats = np.append(exit_nometa_lats,chunk_read.variables['exit_nometa_lats'][:])
    exit_nometa_lons = np.append(exit_nometa_lons,chunk_read.variables['exit_nometa_lons'][:])
    exit_nometa_pg = np.append(exit_nometa_pg,chunk_read.variables['exit_nometa_pg'][:])
    exit_anyvaliddata_refs = np.append(exit_anyvaliddata_refs,chunk_read.variables['exit_anyvaliddata_refs'][:])
    exit_anyvaliddata_lats = np.append(exit_anyvaliddata_lats,chunk_read.variables['exit_anyvaliddata_lats'][:])
    exit_anyvaliddata_lons = np.append(exit_anyvaliddata_lons,chunk_read.variables['exit_anyvaliddata_lons'][:])
    exit_anyvaliddata_pg = np.append(exit_anyvaliddata_pg,chunk_read.variables['exit_anyvaliddata_pg'][:])
    exit_nokeymeta_refs = np.append(exit_nokeymeta_refs,chunk_read.variables['exit_nokeymeta_refs'][:])
    exit_nokeymeta_lats = np.append(exit_nokeymeta_lats,chunk_read.variables['exit_nokeymeta_lats'][:])
    exit_nokeymeta_lons = np.append(exit_nokeymeta_lons,chunk_read.variables['exit_nokeymeta_lons'][:])
    exit_nokeymeta_pg = np.append(exit_nokeymeta_pg,chunk_read.variables['exit_nokeymeta_pg'][:])
    exit_resolution_refs = np.append(exit_resolution_refs,chunk_read.variables['exit_resolution_refs'][:])
    exit_resolution_lats = np.append(exit_resolution_lats,chunk_read.variables['exit_resolution_lats'][:])
    exit_resolution_lons = np.append(exit_resolution_lons,chunk_read.variables['exit_resolution_lons'][:])
    exit_resolution_pg = np.append(exit_resolution_pg,chunk_read.variables['exit_resolution_pg'][:])
    exit_badmeasurementmethod_refs = np.append(exit_badmeasurementmethod_refs,chunk_read.variables['exit_badmeasurementmethod_refs'][:])
    exit_badmeasurementmethod_lats = np.append(exit_badmeasurementmethod_lats,chunk_read.variables['exit_badmeasurementmethod_lats'][:])
    exit_badmeasurementmethod_lons = np.append(exit_badmeasurementmethod_lons,chunk_read.variables['exit_badmeasurementmethod_lons'][:])
    exit_badmeasurementmethod_pg = np.append(exit_badmeasurementmethod_pg,chunk_read.variables['exit_badmeasurementmethod_pg'][:])

    del chunk_read

#write out total counts
dim = root_grp.createDimension('flex2', 1)

i1 = root_grp.createVariable('invalid_nometa_count','int', ('flex2',))
i2 = root_grp.createVariable('invalid_anyvaliddata_count','int', ('flex2',))
i3 = root_grp.createVariable('invalid_nokeymeta_count','int', ('flex2',))
i4 = root_grp.createVariable('invalid_resolution_count','int', ('flex2',))
i5 = root_grp.createVariable('invalid_badmeasurementmethod_count','int', ('flex2',))

n1 = root_grp.createVariable('n_obs_all','int', ('flex2',))
n2 = root_grp.createVariable('n_obs_after_nometa','int', ('flex2',))
n3 = root_grp.createVariable('n_obs_after_flagsandlod','int', ('flex2',))
n4 = root_grp.createVariable('n_obs_after_duplicate','int', ('flex2',))
n5 = root_grp.createVariable('n_obs_after_anyvaliddata','int', ('flex2',))
n6 = root_grp.createVariable('n_obs_after_nokeymeta','int', ('flex2',))
n7 = root_grp.createVariable('n_obs_after_resolution','int', ('flex2',))
n8 = root_grp.createVariable('n_obs_after_badmeasurementmethod','int', ('flex2',))

#site number exit counts
i1[:] = alli1
i2[:] = alli2
i3[:] = alli3
i4[:] = alli4
i5[:] = alli5

#n obs counts after checks
n1[:] = alln1
n2[:] = alln2
n3[:] = alln3
n4[:] = alln4
n5[:] = alln5
n6[:] = alln6
n7[:] = alln7
n8[:] = alln8

#write out exit refs etc..
dim = root_grp.createDimension('flex3', None)
o1 = root_grp.createVariable('UNKNOWN_MM',str, ('flex3',))
o2 = root_grp.createVariable('UNKNOWN_MM_REFS',str, ('flex3',))
o3 = root_grp.createVariable('UNKNOWN_TZ_REFS',str, ('flex3',))

e1a = root_grp.createVariable('exit_nometa_refs',str, ('flex3',))
e1b = root_grp.createVariable('exit_nometa_lats',str, ('flex3',))
e1c = root_grp.createVariable('exit_nometa_lons',str, ('flex3',))
e1d = root_grp.createVariable('exit_nometa_pg',str, ('flex3',))
e2a = root_grp.createVariable('exit_anyvaliddata_refs',str, ('flex3',))
e2b = root_grp.createVariable('exit_anyvaliddata_lats',str, ('flex3',))
e2c = root_grp.createVariable('exit_anyvaliddata_lons',str, ('flex3',))
e2d = root_grp.createVariable('exit_anyvaliddata_pg',str, ('flex3',))
e3a = root_grp.createVariable('exit_nokeymeta_refs',str, ('flex3',))
e3b = root_grp.createVariable('exit_nokeymeta_lats',str, ('flex3',))
e3c = root_grp.createVariable('exit_nokeymeta_lons',str, ('flex3',))
e3d = root_grp.createVariable('exit_nokeymeta_pg',str, ('flex3',))
e4a = root_grp.createVariable('exit_resolution_refs',str, ('flex3',))
e4b = root_grp.createVariable('exit_resolution_lats',str, ('flex3',))
e4c = root_grp.createVariable('exit_resolution_lons',str, ('flex3',))
e4d = root_grp.createVariable('exit_resolution_pg',str, ('flex3',))
e5a = root_grp.createVariable('exit_badmeasurementmethod_refs',str, ('flex3',))
e5b = root_grp.createVariable('exit_badmeasurementmethod_lats',str, ('flex3',))
e5c = root_grp.createVariable('exit_badmeasurementmethod_lons',str, ('flex3',))
e5d = root_grp.createVariable('exit_badmeasurementmethod_pg',str, ('flex3',))

o1[:] = np.array(list(set(unknown_mm)),dtype='object')
o2[:] = np.array(unknown_mm_refs,dtype='object')
o3[:] = np.array(unknown_tz_refs,dtype='object')

e1a[:] = np.array(exit_nometa_refs,dtype='object')
e1b[:] = np.array(exit_nometa_lats,dtype='object')
e1c[:] = np.array(exit_nometa_lons,dtype='object')
e1d[:] = np.array(exit_nometa_pg,dtype='object')
e2a[:] = np.array(exit_anyvaliddata_refs,dtype='object')
e2b[:] = np.array(exit_anyvaliddata_lats,dtype='object')
e2c[:] = np.array(exit_anyvaliddata_lons,dtype='object')
e2d[:] = np.array(exit_anyvaliddata_pg,dtype='object')
e3a[:] = np.array(exit_nokeymeta_refs,dtype='object')
e3b[:] = np.array(exit_nokeymeta_lats,dtype='object')
e3c[:] = np.array(exit_nokeymeta_lons,dtype='object')
e3d[:] = np.array(exit_nokeymeta_pg,dtype='object')
e4a[:] = np.array(exit_resolution_refs,dtype='object')
e4b[:] = np.array(exit_resolution_lats,dtype='object')
e4c[:] = np.array(exit_resolution_lons,dtype='object')
e4d[:] = np.array(exit_resolution_pg,dtype='object')
e5a[:] = np.array(exit_badmeasurementmethod_refs,dtype='object')
e5b[:] = np.array(exit_badmeasurementmethod_lats,dtype='object')
e5c[:] = np.array(exit_badmeasurementmethod_lons,dtype='object')
e5d[:] = np.array(exit_badmeasurementmethod_pg,dtype='object')

root_grp.close()
