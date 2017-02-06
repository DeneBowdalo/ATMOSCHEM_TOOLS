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
from tzwhere import tzwhere
import pytz
tz_root = tzwhere.tzwhere(shapely=True,forceTZ=True)

process_group = 'CASTNET'

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
all_files = glob.glob('/work/home/db876/observations/surface/%s/CASTNET/*.csv'%(fname_species))
all_files = sorted(all_files)

year_array = np.arange(start_year,end_year+1)

valid_files = []
for i in year_array: 
    for f in all_files:
        if str(i) in f:
            valid_files.append(f)
            
#read obs metadata file
meta_refs = []
meta_sitename = []
meta_class = []
meta_tz = []
meta_lats = []
meta_lons = []
meta_alts = []

with open('CASTNET_META.csv', 'rU') as f:
        reader = csv.reader(f,delimiter=',')
        count = 0
        for row in reader:
            if count != 0:
                meta_refs.append(row[0])
                meta_sitename.append(row[2])
                meta_tz.append(row[8])
                meta_lats=np.append(meta_lats,row[9])
                meta_lons=np.append(meta_lons,row[10])
                meta_alts=np.append(meta_alts,row[11])
                meta_class.append(row[12])
            else:
                pass
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
root_grp = Dataset('CASTNET_SURFACE_%s_%s_%s_%s.nc'%(fname_species,start_year,end_year+1,output_res), 'w')
root_grp.description = 'Hourly Surface %s at CASTNET sites - Program written by Dene Bowdalo'%(fname_species)

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

if species == 'O3':
    for file_i in range(len(valid_files)):
        year_refs = []
        with open(valid_files[file_i], 'rb') as f:
            reader = csv.reader(f,delimiter=',')
            print valid_files[file_i]
            count=0
            for row in reader:
                if count == 0:
                    ref_index = row.index('SITE_ID')
                    datetime_index = row.index('DATE_TIME')
                    ozone_index = row.index('OZONE')
                    ozone_flag_index = row.index('OZONE_F')
                    qa_code_index = row.index('QA_CODE')
            
                elif count != 0:
                    year_refs.append(row[ref_index])
                    all_refs.append(row[ref_index])
                    dt = row[datetime_index]
                    yyyymmdd.append(dt[0:4] + dt[5:7] + dt[8:10])
                    hhmm.append(dt[11:13] + dt[14:16])
                    #DO QUALITY CHECKS TO ENSURE DATA IS VALID
                    if (row[qa_code_index] != 'B') & (row[qa_code_index] != 'C') & (row[qa_code_index] != 'D') & (row[qa_code_index] != 'E') & (row[qa_code_index] != 'F') & (row[qa_code_index] != 'GS') & (row[qa_code_index] != 'I') & (row[qa_code_index] != 'K') & (row[qa_code_index] != 'M') & (row[qa_code_index] != 'P') & (row[qa_code_index] != 'S') & (row[qa_code_index] != 'U') & (row[qa_code_index] != 'X') & (row[qa_code_index] != '1X') & (row[qa_code_index] != '<') & (row[qa_code_index] != '^') & (row[qa_code_index] != '_'):
                        try:
                            vals.append(np.float64(row[ozone_index]))
                        except:
                            vals.append(-99999) 
                    else:
                        vals.append(-99999)
                count+=1
            year_refs = set(year_refs)
            refs = np.append(refs,[i for i in year_refs])

if (species == 'CO') or (species == 'NO'):
    for file_i in range(len(valid_files)):
        year_refs = []
        with open(valid_files[file_i], 'rb') as f:
            reader = csv.reader(f,delimiter=',')
            print valid_files[file_i]
            count=0
            for row in reader:
                if count == 0:
                    ref_index = row.index('SITE_ID')
                    datetime_index = row.index('DATE_TIME')
                    value_index = row.index('VALUE')
                    value_flag_index = row.index('VALUE_F')
                    param_index = row.index('PARAMETER')
                    qa_code_index = row.index('QA_CODE')
            
                elif count != 0:
                    if row[param_index] == species:
                        year_refs.append(row[ref_index])
                        all_refs.append(row[ref_index])
                        dt = row[datetime_index]
                        yyyymmdd.append(dt[0:4] + dt[5:7] + dt[8:10])
                        hhmm.append(dt[11:13] + dt[14:16])
                        #DO QUALITY CHECKS TO ENSURE DATA IS VALID
                        if (row[qa_code_index] != 'B') & (row[qa_code_index] != 'C') & (row[qa_code_index] != 'D') & (row[qa_code_index] != 'E') & (row[qa_code_index] != 'F') & (row[qa_code_index] != 'GS') & (row[qa_code_index] != 'I') & (row[qa_code_index] != 'K') & (row[qa_code_index] != 'M') & (row[qa_code_index] != 'P') & (row[qa_code_index] != 'S') & (row[qa_code_index] != 'U') & (row[qa_code_index] != 'X') & (row[qa_code_index] != '1X') & (row[qa_code_index] != '<') & (row[qa_code_index] != '^') & (row[qa_code_index] != '_'):
                            try:
                                vals.append(np.float64(row[value_index]))
                            except:
                                vals.append(-99999) 
                        else:
                            vals.append(-99999)
                count+=1
            year_refs = set(year_refs)
            refs = np.append(refs,[i for i in year_refs]) 
   
refs_count = Counter(refs)
refs = set(refs)

valid_refs = [i for i in refs]
print 'all refs len = ', len(valid_refs)
valid_refs = sorted(valid_refs)      

all_refs = np.array(all_refs)
yyyymmdd = np.array(yyyymmdd)
hhmm = np.array(hhmm)
vals = np.array(vals)

tz_dict = {'EA':-5,'EAno':-5,'CE':-6,'MO':-7,'MOno':-7,'PA':-8,'AL':-9,'AT':-4,'HAno':-10}
#convert timezones to numbers
for i in range(len(meta_tz)):
    tz = meta_tz[i]
    tz_num = tz_dict[tz]
    meta_tz[i] = tz_num
    
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

    data_valid = True
    site_ref = valid_refs[c]
    print 'ref = ',site_ref,c
    site_test = all_refs == site_ref
    
    site_yyyymmdd = yyyymmdd[site_test]
    site_hhmm = hhmm[site_test]
    site_vals = vals[site_test]
    site_vals = np.array(site_vals)
 
    #add val to total obs count
    n_all += len(site_vals)
 
    #test if site_ref in meta_refs, if not then exit
    if site_ref not in meta_refs:
        inv_nometa += 1
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
    
    #convert blank values to -99999
    test_inv = site_vals == ''
    site_vals[test_inv] = -99999
    
    #convert number invalids to -99999
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
    unit = 'na'
    raw_class_name = meta_class[meta_index]
    site_name = meta_sitename[meta_index]
    country = 'United States'
    contact = 'puchalski.melissa@epa.gov'
    
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
    
    #set st_big and mm_big
    st_big = ['continuous']*len(site_vals)
    if species == 'O3':
        mm_big = ['ultraviolet photometry']*len(site_vals)
    elif (species == 'NO'):
        mm_big = ['chemiluminescence']*len(site_vals)
    elif (species == 'CO'):
        mm_big = ['non-dispersive infrared spectroscopy']*len(site_vals)
        
    #get obs valid
    test = site_vals >= 0
    n_obs_valid = len(site_vals[test])
    n_after_flagsandlod += n_obs_valid
    print site_vals,n_after_flagsandlod
    
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
    converted_time,site_vals,mm_big,st_big,na = modules.remove_duplicate_points(site_ref,converted_time,site_vals,mm_big,st_big,'blank',output_res)
    test = site_vals >= 0
    n_obs_valid = int(len(site_vals[test]))
    n_after_duplicate += n_obs_valid
    
    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
    indices = np.searchsorted(syn_grid_time, converted_time, side='left')
    full_data[indices] = site_vals 
    
    key_meta = [lat,lon,alt]
    
    #set site file resolution
    file_res = 'H'
    
    #get sampling/instrument grids
    raw_st_grid,p_st_grid,p_st_grid_after_flagsandlod,raw_mm_grid,p_mm_grid,p_mm_grid_after_flagsandlod,unknown_mm_list,unknown_mm_refs_list = modules.sampling_and_instruments_by_processgroup(site_ref,process_group,species,raw_st,raw_mm,full_data_after_flagsandlod,full_data,raw_indices,unknown_mm_list,unknown_mm_refs_list,no2_type)

    #do quality checks                                                                                                                                                                                                                                                                                                     
    data_valid,full_data,valid_hours_dup,p_st_grid,p_mm_grid,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_after_anyvaliddata,inv_nokeymeta,n_after_nokeymeta,inv_resolution,n_after_resolution,inv_badmeasurementmethod,n_after_badmeasurementmethod,exit_r = modules.primary_quality_control(site_ref,species,file_res,no2_type,grid_dates,full_data,big_n_dup_array,0,raw_st_grid,p_st_grid,p_st_grid_after_flagsandlod,raw_mm_grid,p_mm_grid,p_mm_grid_after_flagsandlod,data_resolution,n_obs_valid,key_meta,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_after_anyvaliddata,inv_nokeymeta,n_after_nokeymeta,inv_resolution,n_after_resolution,inv_badmeasurementmethod,n_after_badmeasurementmethod)
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

