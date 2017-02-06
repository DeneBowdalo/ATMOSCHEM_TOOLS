import numpy as np
import glob
from collections import Counter
import datetime
import csv
from netCDF4 import Dataset
import math
import modules
from time import sleep
import pandas as pd
import os
import multiprocessing
from netCDF4 import num2date, date2num
from pyexcel_xls import get_data
from tzwhere import tzwhere
import pytz
tz_root = tzwhere.tzwhere(shapely=True,forceTZ=True)

process_group = 'NAPS'

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

#set run type as serial or parallel
run_type = 'parallel'

#get molecular mass and required concentration resolution of species
data_resolution,mol_mass,aqs_code,airbase_code = modules.get_spec_specs(species)

#SET EXIT ERROR COUNTS AT 0
inv_nometa_count = 0
inv_anyvaliddata_count = 0
inv_nokeymeta_count = 0
inv_resolution_count = 0
inv_badmeasurementmethod_count = 0

exit_counts = np.array([inv_nometa_count,inv_anyvaliddata_count,inv_nokeymeta_count,inv_resolution_count,inv_badmeasurementmethod_count])

#SET N OBS COUNTS
n_obs_all = 0
n_obs_after_nometa = 0
n_obs_after_flagsandlod = 0
n_obs_after_duplicate = 0
n_obs_after_anyvaliddata = 0
n_obs_after_nokeymeta = 0
n_obs_after_resolution = 0
n_obs_after_badmeasurementmethod = 0

n_obs_counts = np.array([n_obs_all,n_obs_after_nometa,n_obs_after_flagsandlod,n_obs_after_duplicate,n_obs_after_anyvaliddata,n_obs_after_nokeymeta,n_obs_after_resolution,n_obs_after_badmeasurementmethod])

#set exit code output 
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

unknown_mm = []
unknown_mm_refs = []
unknown_local_tz = []

#read obs files
all_files = glob.glob('/work/home/db876/observations/surface/%s/CANADANAPS/*.hly'%(fname_species))
all_files = sorted(all_files)

year_array = np.arange(start_year,end_year+1)

if species != 'ISOP':
    valid_files = []
    for i in year_array: 
        for f in all_files:
            if str(i) in f:
                valid_files.append(f)
                
#read obs metadata file
meta_refs = []
meta_sitenames = []
meta_class = []
meta_tz = []
meta_lats = []
meta_lons = []
meta_alts = []
meta_contacts = []
meta_countries = []

with open('NAPS_META.csv', 'rb') as f:
        reader = csv.reader(f,delimiter=',')
        count = 0
        for row in reader:
            if count == 0:
                sitename_i = row.index('MAIN_STATIONS_STATION_NAME')
                class_i = row.index('Type')
                tz_i = row.index('TimeZone')
                lat_i = row.index('Lat_Decimal')
                lon_i = row.index('Long_Decimal')
                alt_i = row.index('Elevation_m')
                contact_i = row.index('Contact')
                country_i = row.index('COUNTRY')
            
            if count != 0:
                current_ref = row[0]
                if len(current_ref) == 5:
                    current_ref = '0'+current_ref
                
                #if tz is blank then do not add metdata
                if (row[16] == ''):
                    continue
                else:
                    meta_refs.append(current_ref)
                    meta_sitenames.append(row[sitename_i])
                    meta_class.append(row[class_i])
                    meta_tz.append(row[tz_i])
                    meta_lats=np.append(meta_lats,row[lat_i])
                    meta_lons=np.append(meta_lons,row[lon_i])
                    meta_alts=np.append(meta_alts,row[alt_i])
                    meta_contacts = np.append(meta_contacts,row[contact_i])
                    meta_countries = np.append(meta_countries,row[country_i])
                
            count+=1
            
#find n_hours and n_days between start and end date
d0 = datetime.date(start_year, 1, 1)
d1 = datetime.date(end_year+1, 1, 1)
delta = d1 - d0
n_months = ((end_year+1)-start_year)*12
n_days = delta.days
n_hours = n_days*24
n_points = n_hours

#setup netcdf output
root_grp = Dataset('NAPS_SURFACE_%s_%s_%s_%s.nc'%(fname_species,start_year,end_year+1,output_res), 'w')
root_grp.description = 'Hourly Surface %s at NAPS sites - Program written by Dene Bowdalo'%(fname_species)

# dimensions
root_grp.createDimension('flex', n_points)
root_grp.createDimension('flex4', None)

#create grid time array
start = datetime.datetime(year = start_year, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = end_year+1, month = 1, day = 1, hour = 0, minute = 0)
grid_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]
grid_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='H')][:-1]

#make output dates and times
dt = pd.date_range(start,end,freq='H')[:-1].tolist()
output_res_times = date2num(dt, units='hours since 0001-01-01 00:00:00', calendar='gregorian')
obs_time_pd = pd.date_range(start = start,end = end, freq = 'H')[:-1]

#save time out
times = root_grp.createVariable('time', 'f8', ('flex',))
times[:] = output_res_times

all_refs = []
refs = []
valid_refs  = []
yyyymmdd = []
hhmm = []

vals = []
hour_array = ['0000','0100','0200','0300','0400','0500','0600','0700','0800','0900','1000',
              '1100','1200','1300','1400','1500','1600','1700','1800','1900','2000','2100',
              '2200','2300']


if species != 'ISOP':
    for file_i in range(len(valid_files)):
        year_refs = []
        with open(valid_files[file_i], 'rb') as f:
            reader = csv.reader(f,delimiter='?')
            print valid_files[file_i]
            for row in reader:
                row = row[0]
                year_refs.append(row[3:9])
                for i in range(24):
                    all_refs.append(row[3:9])
                for i in range(24):
                    yyyymmdd.append(row[9:13]+row[13:15]+row[15:17])
            
                start_n = 29
                stop_n = 33
                for i in range(24):
                    vals.append(row[start_n:stop_n])
                    hhmm.append(hour_array[i])
                    start_n+=4
                    stop_n+=4

            year_refs = set(year_refs)
            refs = np.append(refs,[i for i in year_refs])

elif species == 'ISOP':
    all_years = np.array([i.strip('../CANADANAPS/')[3:] for i in glob.glob('../CANADANAPS/*')]).astype(int)
    all_years = sorted(all_years)
    for y in all_years:
        files = glob.glob('../CANADANAPS/VOC%s/*'%(y))
        files = [i.replace('../CANADANAPS/VOC%s/'%(y),'') for i in files]
        year_refs = np.array([i.split('_')[0][1:] for i in files])
        refs = np.append(refs,[i for i in year_refs])
        for i in range(len(refs)):
            if len(refs[i]) == 5:
                refs[i] = '0'+refs[i]
  
valid_refs = set(refs)
valid_refs = sorted([i for i in valid_refs])
print 'all refs len = ', len(valid_refs)

all_refs = np.array(all_refs)
yyyymmdd = np.array(yyyymmdd)
hhmm = np.array(hhmm)
vals = np.array(vals)

def site_iter_process(valid_refs,c):
    #set local counts
    inv_nometa = 0
    inv_anyvaliddata = 0
    inv_nokeymeta = 0
    inv_resolution = 0
    inv_badmeasurementmethod = 0
    n_all = 0
    n_after_nometa = 0
    n_after_flagsandlod = 0
    n_after_duplicate = 0
    n_after_anyvaliddata = 0
    n_after_nokeymeta = 0
    n_after_resolution = 0
    n_after_badmeasurementmethod = 0

    #set local unknown lists
    unknown_mm_list = []
    unknown_mm_refs_list = []
    unknown_local_tz_list = []
    
    site_resolutions = []

    site_ref = valid_refs[c]
 
    data_valid = True
    print 'ref = ',site_ref,c
    
    if species != 'ISOP':
        site_test = all_refs == site_ref
        site_yyyymmdd = yyyymmdd[site_test]
        site_hhmm = hhmm[site_test]
        site_vals = vals[site_test]
        n_dup_array = np.array([0]*len(site_vals))
    else:
        if site_ref[0] == '0':
            site_ref = site_ref[1:]
        files = []
        site_yyyymmdd = []
        site_hhmm = []
        site_vals = []
        n_dup_array = []
        for y in all_years:
            try:
                files.append(glob.glob('../CANADANAPS/VOC%s/S%s*'%(y,site_ref)))
            except:
                pass
        files = [item for sublist in files for item in sublist]
        for f in files:
            print f
            all_data = get_data(f)
            all_data = all_data.values()
            test_header_range = range(0,10)
            for x in test_header_range:
                headers = all_data[0][x]
                if 'Isoprene' in headers:
                    header_x = x
                    break
            data_cut = all_data[0][header_x+1:]
            var_i = headers.index('Isoprene')
            #date_i = headers.index('Sample Date')
            date_i = headers.index('Compounds')
            time_i = headers.index('START TIME')
            duration_i = headers.index('DURATION')
            
            for i in range(len(data_cut)):
                row_cut = data_cut[i]
                
                try:
                    dur = float(row_cut[duration_i])
                    if dur.is_integer() == False:
                        dur = round(dur,0)
                except:
                    #round to nearest hour if necessary
                    if float(row_cut[duration_i].strftime("%M")) != 0:
                        if dur >= 30:
                            dur = float(row_cut[duration_i].strftime("%H"))+1
                        else:
                            dur = float(row_cut[duration_i].strftime("%H"))
                    else:
                        dur = float(row_cut[duration_i].strftime("%H"))
                    
                if dur.is_integer() == False:
                    print 'duration is float'
                    1+'a'
    
                try:
                    val =np.float64(row_cut[var_i])
                except:
                    val = -99999
                
                if dur == 1:
                    site_resolutions.append('H')
                    
                    #if (val >= 0.01):
                    #    site_vals.append([val])
                    #else:
                    #    site_vals.append([-99999]) 
                    site_vals.append([val])
                    
                    n_dup_array.append([0])
                    site_yyyymmdd.append([row_cut[date_i].strftime("%Y%m%d")])
                    try:
                        site_hhmm.append([row_cut[time_i][:2]+row_cut[time_i][3:5]])
                    except:
                        #round to nearest hour if necessary
                        ti = row_cut[time_i].strftime("%H%M")
                        if float(row_cut[time_i].strftime("%M")) != 0:
                            print 'non whole time = ', row_cut[time_i]
                            if float(row_cut[time_i].strftime("%M")) >= 30:
                                site_hhmm.append([datetime.time(hour=int(ti[:2])+1,minute=0).strftime("%H%M")])
                            else:
                                site_hhmm.append([datetime.time(hour=int(ti[:2]),minute=0).strftime("%H%M")])
                                
                        else:
                            site_hhmm.append([row_cut[time_i].strftime("%H%M")])
                #deal with sample lens > 1 hour
                else:
                    if output_res == 'H':
                        continue
                    else:
                        site_resolutions.append('D')
                        
                        #if (val >= 0.01):
                        #    site_vals.append([val])
                        #else:
                        #    site_vals.append([-99999]) 
                        site_vals.append([val])
                    
                        n_dup_array.append([0])
                        
                        try:
                            site_yyyymmdd.append([row_cut[date_i].strftime("%Y%m%d")])
                        except:
                            print row_cut[date_i]
                            1+'a'
                        try:
                            site_hhmm.append([row_cut[time_i][:2]+row_cut[time_i][3:5]])
                        except:
                            #round to nearest hour if necessary
                            ti = row_cut[time_i].strftime("%H%M")
                            if float(row_cut[time_i].strftime("%M")) != 0:
                                print 'non whole time = ', row_cut[time_i]
                                if float(row_cut[time_i].strftime("%M")) >= 30:
                                    site_hhmm.append([datetime.time(hour=int(ti[:2])+1,minute=0).strftime("%H%M")])
                                else:
                                    site_hhmm.append([datetime.time(hour=int(ti[:2]),minute=0).strftime("%H%M")])
                                
                            else:
                                site_hhmm.append([row_cut[time_i].strftime("%H%M")])
                            
                        current_year = int(site_yyyymmdd[-1][0][:4])
                        current_month = int(site_yyyymmdd[-1][0][4:6])
                        current_day = int(site_yyyymmdd[-1][0][6:])
                        current_hh = int(site_hhmm[-1][0][:2])
                        current_mm = int(site_hhmm[-1][0][2:])
                        
                        s = datetime.datetime(year = current_year, month = current_month, day = current_day, hour = current_hh, minute = current_mm)
                        e = s+datetime.timedelta(hours=dur)
                        day_dates = [d.strftime('%Y%m%d') for d in pd.date_range(s,e,freq='H')][1:-1]
                        day_hours = [d.strftime('%H%M') for d in pd.date_range(s,e,freq='H')][1:-1]
                        
                        site_yyyymmdd.append(day_dates)
                        site_hhmm.append(day_hours)
                        site_vals.append([site_vals[-1][0]]*len(day_dates))

                        #append to n duplicated array
                        n_dup_array.append([0])
                        n_dup_array.append([1]*len(day_dates))
    
    
    if species == 'ISOP':
        site_yyyymmdd = [item for sublist in site_yyyymmdd for item in sublist]
        site_hhmm = [item for sublist in site_hhmm for item in sublist]
        site_vals = [item for sublist in site_vals for item in sublist]
        n_dup_array = np.array([item for sublist in n_dup_array for item in sublist])
        if len(site_ref) == 5:
            site_ref = '0'+site_ref
    
    site_vals = np.float64(site_vals)
    
    #add val to total obs count
    n_all += len(site_vals)
    
    #test if site_ref in meta_refs, if not then exit
    if site_ref not in meta_refs:
        print site_ref
        inv_nometa+=1
        print 'Site Invalid. No Metadata for ref'
        if no2_type == 'MOLYBDENUM':
            n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod = 0,0,0,0,0,0,0,0,0,0,0,0,0
        exit_c_list = np.array([inv_nometa,inv_anyvaliddata,inv_nokeymeta,inv_resolution,inv_badmeasurementmethod])
        n_c_list = np.array([n_all,n_after_nometa,n_after_flagsandlod,n_after_duplicate,n_after_anyvaliddata,n_after_nokeymeta,n_after_resolution,n_after_badmeasurementmethod])
        unknown_list = [unknown_mm_list,unknown_mm_refs_list,unknown_local_tz_list]
        meta = ['na','na','na','na','na','na','na','na','na','na','na','na']
        exit_r = 'nometa'
        return c,['na'],['na'],['na'],False,meta,exit_c_list,n_c_list,unknown_list,exit_r,np.zeros(1)
    n_after_nometa += len(site_vals)
 
    #convert all invalids to -99999
    test_inv = site_vals < 0
    site_vals[test_inv] = -99999

    #create max possible grid
    full_data = np.empty(n_hours)
    full_data_after_flagsandlod = np.empty(n_hours)
    big_n_dup_array = np.zeros(n_hours)
    full_data[:] = -99999
    full_data_after_flagsandlod[:] = -99999
    
    #get meta
    meta_index = meta_refs.index(site_ref)
    data_tz = np.float32(meta_tz[meta_index])
    all_tz = [data_tz]
    try:
        lat = np.float32(meta_lats[meta_index])
    except:
        lat = 'na'
    try:
        lon = np.float32(meta_lons[meta_index])
    except:
        lon = 'na'
    try:
        alt = np.float32(meta_alts[meta_index])
    except:
        alt = 'na'
    raw_class_name = meta_class[meta_index]
    site_name = meta_sitenames[meta_index]
    unit = 'na'
    contact = meta_contacts[meta_index]
    country = meta_countries[meta_index]

    #adjust dates and times if tz is not equal to 0
    tz = int(data_tz)
    if tz != 0:
        for i in range(len(site_yyyymmdd)):
            #create datetime
            dt = datetime.datetime(int(site_yyyymmdd[i][:4]),int(site_yyyymmdd[i][4:6]),int(site_yyyymmdd[i][6:]),int(site_hhmm[i][:2]),int(site_hhmm[i][2:]))
            if tz > 0:
                dt  = dt - datetime.timedelta(hours = int(tz))
            elif tz < 0:
                dt  = dt + datetime.timedelta(hours = np.abs(int(tz)))
            site_yyyymmdd[i] = dt.strftime("%Y%m%d")
            site_hhmm[i] = dt.strftime("%H%M")

    #put vals into full grid
    date_con = np.array(site_yyyymmdd).astype(int)
    time_con = np.array(site_hhmm).astype(int)
    
    #remove data < 1970 and >= 2015
    test_inds = (date_con >= 19700101) & (date_con < 20150101)
    date_con = date_con[test_inds]
    time_con = time_con[test_inds]
    site_vals = site_vals[test_inds]
    n_dup_array = n_dup_array[test_inds]
    
    #set st_big and mm_big
    st_big = ['continuous']*len(site_vals)
        
    if species == 'O3':
        mm_big = ['ultraviolet photometry']*len(site_vals)
    elif species == 'NO':
        mm_big = ['chemiluminescence']*len(site_vals)
    elif species == 'NO2':
        mm_big = ['chemiluminescence (conversion-molybdenum)']*len(site_vals)
    elif species == 'CO':
        mm_big = ['non-dispersive infrared spectrometry']*len(site_vals)
    elif species == 'ISOP':
        mm_big = ['gas chromatography mass selective detection']*len(site_vals)
    
    #get obs valid after flagsandlod
    test = site_vals != -99999
    valid_hours_dup = np.sum(n_dup_array[test])
    n_obs_valid = len(site_vals[test]-valid_hours_dup)
    n_after_flagsandlod += n_obs_valid
    
    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
    converted_time = modules.date_process(date_con,time_con,start_year)
    converted_time = np.round(converted_time,decimals=5)
    syn_grid_time = np.arange(0,n_days,1./24)
    syn_grid_time = np.round(syn_grid_time,decimals=5)
    raw_indices = np.searchsorted(syn_grid_time, converted_time, side='left')
    site_vals = np.array(site_vals)
    full_data_after_flagsandlod[raw_indices] = site_vals
    raw_st = np.copy(st_big)
    raw_mm = np.copy(mm_big)
    
    # test and remove duplicate and overlap points
    converted_time,site_vals,mm_big,st_big,n_dup_array = modules.remove_duplicate_points(site_ref,converted_time,site_vals,mm_big,st_big,n_dup_array,output_res)
    test = site_vals != -99999
    valid_hours_dup = np.sum(n_dup_array[test])
    n_obs_valid = int(len(site_vals[test]) - valid_hours_dup)
    print 'n obs valid = ',n_obs_valid
    n_after_duplicate += n_obs_valid
    
    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
    indices = np.searchsorted(syn_grid_time, converted_time, side='left')
    full_data[indices] = site_vals
    big_n_dup_array[indices] = n_dup_array
        
    #if species is CO then convert units from ppmv to ppbv
    if species == 'CO':
        valid_inds = full_data != -99999 
        full_data[valid_inds] = full_data[valid_inds]*1e3   
        
    #if species is ISOP then connvert units from mg/m3 to ppbv
    if species == 'ISOP':
        #calculate conversion factor from mg/m3 assuming 25 degC and 1 atm 
        #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
        conv_fact = 8.3144/mol_mass*(273.15+25)/(1013.25/10)
        valid_inds = full_data != -99999 
        full_data[valid_inds] = full_data[valid_inds]*conv_fact
    
    key_meta = [lat,lon,alt]
    
    #set site file resolution
    if (species == 'O3') or (species == 'CO') or(species == 'NO') or (species == 'NO2'):
        file_res = 'H'
    else:
        # if no valid data then site res does not matter
        if len(site_resolutions) == 0:
            file_res = 'na'
        else:
            #if all site resolutions are same continue then take first file_res
            all_same = all(x == site_resolutions[0] for x in site_resolutions)
            if all_same == True:
                file_res = site_resolutions[0]
            else:
            #otherwise take highest frequency res as file_res 
                if 'M' in site_resolutions:
                    file_res = 'M'
                elif 'D' in site_resolutions:
                    file_res = 'D'
                else:
                    file_res = 'H'
        
    #get sampling/instrument grids
    raw_st_grid,p_st_grid,p_st_grid_after_flagsandlod,raw_mm_grid,p_mm_grid,p_mm_grid_after_flagsandlod,unknown_mm_list,unknown_mm_refs_list = modules.sampling_and_instruments_by_processgroup(site_ref,process_group,species,raw_st,raw_mm,full_data_after_flagsandlod,full_data,raw_indices,unknown_mm_list,unknown_mm_refs_list,no2_type)
    
    #do quality checks  
    data_valid,full_data,valid_hours_dup,p_st_grid,p_mm_grid,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_after_anyvaliddata,inv_nokeymeta,n_after_nokeymeta,inv_resolution,n_after_resolution,inv_badmeasurementmethod,n_after_badmeasurementmethod,exit_r = modules.primary_quality_control(site_ref,species,file_res,no2_type,grid_dates,full_data,big_n_dup_array,valid_hours_dup,raw_st_grid,p_st_grid,p_st_grid_after_flagsandlod,raw_mm_grid,p_mm_grid,p_mm_grid_after_flagsandlod,data_resolution,n_obs_valid,key_meta,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_after_anyvaliddata,inv_nokeymeta,n_after_nokeymeta,inv_resolution,n_after_resolution,inv_badmeasurementmethod,n_after_badmeasurementmethod)
    if data_valid == False:
        exit_c_list = np.array([inv_nometa,inv_anyvaliddata,inv_nokeymeta,inv_resolution,inv_badmeasurementmethod])
        n_c_list = np.array([n_all,n_after_nometa,n_after_flagsandlod,n_after_duplicate,n_after_anyvaliddata,n_after_nokeymeta,n_after_resolution,n_after_badmeasurementmethod])
        unknown_list = [unknown_mm_list,unknown_mm_refs_list,unknown_local_tz_list]
        meta = [lat,lon,alt,'na','na','na','na','na','na','na','na','na']
        return c,['na'],['na'],['na'],False,meta,exit_c_list,n_c_list,unknown_list,exit_r,np.zeros(1)
    
    #make tz int after checks
    data_tz = np.float32(data_tz)
    
    #set processed unit
    p_unit = 'pbbv'
            
    #get local timezone
    try:
        local_tz_name = tz_root.tzNameAt(lat,lon,forceTZ=True)
        pytz_obj = pytz.timezone(local_tz_name)
        datetime_offset = pytz_obj.utcoffset(datetime.datetime(2000,1,1))
        if datetime_offset < datetime.timedelta(0):
            local_tz = -(24-int(datetime_offset.seconds/60/60))
        else:
            local_tz = int(datetime_offset.seconds/60/60)
    except:
        local_tz = 'na'
        print 'TIMEZONE NOT KNOWN, SITE IS %s'%(site_ref)
        unknown_local_tz_list.append(site_ref)

    #pack meta
    meta = [lat,lon,alt,raw_class_name,file_res,unit,p_unit,data_tz,local_tz,site_name,country,contact]

    #if blank strings in meta then convert to 'na'
    for i in range(len(meta)):
        try:
            if meta[i].strip() == '':
                meta[i] = 'na'
        except:
            pass

    print set(raw_st_grid)
    print set(raw_mm_grid)
    print set(p_st_grid)
    print set(p_mm_grid)
    print meta

    exit_c_list = np.array([inv_nometa,inv_anyvaliddata,inv_nokeymeta,inv_resolution,inv_badmeasurementmethod])
    n_c_list = np.array([n_all,n_after_nometa,n_after_flagsandlod,n_after_duplicate,n_after_anyvaliddata,n_after_nokeymeta,n_after_resolution,n_after_badmeasurementmethod])
    print 'exit counts = ', exit_c_list
    print 'n obs counts = ', n_c_list

    unknown_list = [unknown_mm_list,unknown_mm_refs_list,unknown_local_tz_list]

    return c,full_data,p_st_grid,p_mm_grid,data_valid,meta,exit_c_list,n_c_list,unknown_list,'na',big_n_dup_array
                

if run_type == 'serial':
    for c in range(len(valid_refs)):
        c,full_data,p_st_grid,p_mm_grid,data_valid,meta,ec,nc,unknown,exit_code,n_dup = site_iter_process(valid_refs,c)
        
        #change ref name
        valid_refs[c] = 'canaps'+valid_refs[c]
        #add counts up
        exit_counts = exit_counts + ec
        n_obs_counts = n_obs_counts + nc
        #append to unknown lists
        unknown_mm=np.append(unknown_mm,unknown[0])
        unknown_mm_refs=np.append(unknown_mm_refs,unknown[1])
        unknown_local_tz=np.append(unknown_local_tz,unknown[2])
        
        #append exit refs
        if no2_type != 'MOLYBDENUM':
            if exit_code != 'na':
                if exit_code == 'nometa':
                    exit_nometa_refs.append(valid_refs[c])
                    exit_nometa_lats.append(str(meta[0]))
                    exit_nometa_lons.append(str(meta[1]))
                    exit_nometa_pg.append(process_group)
                elif exit_code == 'anyvaliddata':
                    exit_anyvaliddata_refs.append(valid_refs[c])
                    exit_anyvaliddata_lats.append(str(meta[0]))
                    exit_anyvaliddata_lons.append(str(meta[1]))
                    exit_anyvaliddata_pg.append(process_group)
                elif exit_code == 'nokeymeta':
                    exit_nokeymeta_refs.append(valid_refs[c])
                    exit_nokeymeta_lats.append(str(meta[0]))
                    exit_nokeymeta_lons.append(str(meta[1]))
                    exit_nokeymeta_pg.append(process_group)
                elif exit_code == 'resolution':
                    exit_resolution_refs.append(valid_refs[c])
                    exit_resolution_lats.append(str(meta[0]))
                    exit_resolution_lons.append(str(meta[1]))
                    exit_resolution_pg.append(process_group)
                elif exit_code == 'badmeasurementmethod':
                    exit_badmeasurementmethod_refs.append(valid_refs[c])
                    exit_badmeasurementmethod_lats.append(str(meta[0]))
                    exit_badmeasurementmethod_lons.append(str(meta[1]))
                    exit_badmeasurementmethod_pg.append(process_group)
        
        modules.write_out_data(valid_refs[c],process_group,root_grp,species,full_data,p_st_grid,p_mm_grid,data_valid,meta,n_dup)
            
elif run_type == 'parallel':
    if __name__ == '__main__':
        pool = multiprocessing.Pool(processes=16)
        results = [pool.apply_async(site_iter_process, (valid_refs,c)) for c in range(len(valid_refs))]
        big_array = [r.get() for r in results]
        pool.terminate()
    
    indices_array = []
    full_data_array = []
    p_st_grid_array = []
    p_mm_grid_array = []
    data_valid_array = []
    meta_array = []
    exit_counts_array = []
    n_obs_counts_array = []
    unknown_array = []
    exit_code_array = []
    n_dup_array = []
    
    for i in range(len(big_array)):
        cut = big_array[i]
        indices_array.append(cut[0])
        full_data_array.append(cut[1])
        p_st_grid_array.append(cut[2])
        p_mm_grid_array.append(cut[3])
        data_valid_array.append(cut[4])
        meta_array.append(cut[5])
        exit_counts_array.append(cut[6])
        n_obs_counts_array.append(cut[7])
        unknown_array.append(cut[8])
        exit_code_array.append(cut[9])
        n_dup_array.append(cut[10])

    indices_array = np.array(indices_array)
    full_data_array = np.array(full_data_array)
    p_st_grid_array = np.array(p_st_grid_array)
    p_mm_grid_array = np.array(p_mm_grid_array)
    data_valid_array = np.array(data_valid_array)
    meta_array = np.array(meta_array)
    exit_counts_array = np.array(exit_counts_array)
    n_obs_counts_array = np.array(n_obs_counts_array)
    unknown_array = np.array(unknown_array)
    exit_code_array = np.array(exit_code_array)
    n_dup_array = np.array(n_dup_array)
    
    #sort arrays by indices array for sanity
    full_data_array = full_data_array[indices_array]
    p_st_grid_array = p_st_grid_array[indices_array]
    p_mm_grid_array = p_mm_grid_array[indices_array]
    data_valid_array = data_valid_array[indices_array]
    meta_array = meta_array[indices_array]
    exit_counts_array = exit_counts_array[indices_array]
    n_obs_counts_array = n_obs_counts_array[indices_array]
    unknown_array = unknown_array[indices_array]
    exit_code_array = exit_code_array[indices_array]
    n_dup_array = n_dup_array[indices_array]
    
    for c in range(len(valid_refs)):    
        #change ref name
        valid_refs[c] = 'canaps'+valid_refs[c]
        
        #add counts up
        exit_counts = exit_counts + exit_counts_array[c]
        n_obs_counts = n_obs_counts + n_obs_counts_array[c]
        
        #append to unknown lists
        unknown_mm=np.append(unknown_mm,unknown_array[c][0])
        unknown_mm_refs=np.append(unknown_mm_refs,unknown_array[c][1])
        unknown_local_tz=np.append(unknown_local_tz,unknown_array[c][2])
        
        #append exit refs
        if no2_type != 'MOLYBDENUM':
            if exit_code_array[c] != 'na':
                if exit_code_array[c] == 'nometa':
                    exit_nometa_refs.append(valid_refs[c])
                    exit_nometa_lats.append(str(meta_array[c][0]))
                    exit_nometa_lons.append(str(meta_array[c][1]))
                    exit_nometa_pg.append(process_group)
                elif exit_code_array[c] == 'anyvaliddata':
                    exit_anyvaliddata_refs.append(valid_refs[c])
                    exit_anyvaliddata_lats.append(str(meta_array[c][0]))
                    exit_anyvaliddata_lons.append(str(meta_array[c][1]))
                    exit_anyvaliddata_pg.append(process_group)
                elif exit_code_array[c] == 'nokeymeta':
                    exit_nokeymeta_refs.append(valid_refs[c])
                    exit_nokeymeta_lats.append(str(meta_array[c][0]))
                    exit_nokeymeta_lons.append(str(meta_array[c][1]))
                    exit_nokeymeta_pg.append(process_group)
                elif exit_code_array[c] == 'resolution':
                    exit_resolution_refs.append(valid_refs[c])
                    exit_resolution_lats.append(str(meta_array[c][0]))
                    exit_resolution_lons.append(str(meta_array[c][1]))
                    exit_resolution_pg.append(process_group)
                elif exit_code_array[c] == 'badmeasurementmethod':
                    exit_badmeasurementmethod_refs.append(valid_refs[c])
                    exit_badmeasurementmethod_lats.append(str(meta_array[c][0]))
                    exit_badmeasurementmethod_lons.append(str(meta_array[c][1]))
                    exit_badmeasurementmethod_pg.append(process_group)
    
        modules.write_out_data(valid_refs[c],process_group,root_grp,species,full_data_array[c],p_st_grid_array[c],p_mm_grid_array[c],data_valid_array[c],meta_array[c],n_dup_array[c]) 

#save out processing stats
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
i1[:] = exit_counts[0]
i2[:] = exit_counts[1]
i3[:] = exit_counts[2]
i4[:] = exit_counts[3]
i5[:] = exit_counts[4]

print exit_counts

#n obs counts after checks
n1[:] = n_obs_counts[0]
n2[:] = n_obs_counts[1]
n3[:] = n_obs_counts[2]
n4[:] = n_obs_counts[3]
n5[:] = n_obs_counts[4]
n6[:] = n_obs_counts[5]
n7[:] = n_obs_counts[6]
n8[:] = n_obs_counts[7]


print n_obs_counts

print 'Unknown mm = ', set(unknown_mm)
print 'Unknown mm Sites = ', unknown_mm_refs 
print 'Unknown local tz sites = ', unknown_local_tz

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
o3[:] = np.array(unknown_local_tz,dtype='object')

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

print 'NOMETA_EXIT_REFS = ', exit_nometa_refs
print 'ANYVALIDDATA_EXIT_REFS = ', exit_anyvaliddata_refs
print 'NOKEYMETA_EXIT_REFS = ', exit_nokeymeta_refs
print 'RESOLUTION_EXIT_REFS = ', exit_resolution_refs
print 'BADMEASUREMENTMETHOD_EXIT_REFS = ', exit_badmeasurementmethod_refs

root_grp.close()
