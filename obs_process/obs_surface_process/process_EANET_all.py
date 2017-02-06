import numpy as np
import glob
from collections import Counter
import datetime
import csv
from netCDF4 import Dataset
import modules
from time import sleep
import pandas as pd
from calendar import monthrange
import os
import multiprocessing
from netCDF4 import num2date, date2num
from dateutil.relativedelta import relativedelta
from tzwhere import tzwhere
import pytz
tz_root = tzwhere.tzwhere(shapely=True,forceTZ=True)

def insensitive_glob(pattern):
    def either(c):
        return '[%s%s]'%(c.lower(),c.upper()) if c.isalpha() else c
    return glob.glob(''.join(map(either,pattern)))

process_group = 'EANET'

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
data_resolution,species_mw,aqs_code,airbase_code = modules.get_spec_specs(species)

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

year_array = range(start_year,end_year+1)
n_years = (end_year-start_year)+1

#read obs files
all_files = glob.glob('/work/home/db876/observations/surface/%s/EANET/AT*'%(fname_species))
all_files = modules.natsorted(all_files)

ref_list = []
valid_refs = []

for i in range(len(all_files)):
    f = all_files[i].replace(".csv","")
    f = f.replace("/work/home/db876/observations/surface/%s/EANET/"%(fname_species),"")
    f = f[4:]
    f = f.lower()
    ref_list.append(f)
 
refs = set(ref_list)
refs = [i for i in refs]
refs = sorted(refs)
valid_refs = np.array(refs)
print 'all refs len = ', len(valid_refs)

#find n_hours and n_days between start and end date
d0 = datetime.date(start_year, 1, 1)
d1 = datetime.date(end_year+1, 1, 1)
delta = d1 - d0
n_months = ((end_year+1)-start_year)*12
n_days = delta.days
n_hours = n_days*24
n_points = n_hours

#setup netcdf output
root_grp = Dataset('EANET_SURFACE_%s_%s_%s_%s.nc'%(fname_species,start_year,end_year+1,output_res), 'w')
root_grp.description = 'Hourly Surface %s at EANET sites- Program written by Dene Bowdalo'%(fname_species)

# dimensions
root_grp.createDimension('flex', n_points)
root_grp.createDimension('flex4', None)

#read in metadata
count = 0

met_refs = []
met_sitenames = []
file_refs = []
met_lats_degs = []
met_lats_mins = []
met_lats_dirs = []
met_lons_degs = []
met_lons_mins = []
met_lons_dirs = []
met_alts = []
met_class = []
met_tz = []
met_aq_check = []
met_species = []
met_instruments = []
met_country = []

with open('EANET_COREMETA.txt', 'U') as f:
    reader = csv.reader(f,delimiter=',')
    for x,row in enumerate(reader):
        if count == 0:
            met_sitenames.append(row[0].strip())
            line = row[0].strip().lower().replace(" ", "")
            met_refs.append(line)
            file_refs.append('EA'+line[0:2].upper()+'%03i'%x)
            
        #read lat
        elif count == 2:
            deg = row[0]
            deg = deg[10:]
            deg = deg.strip()
            
            min = row[1]
            
            dir = row[2].strip()

            met_lats_degs.append(deg)
            met_lats_mins.append(min)
            met_lats_dirs.append(dir)

        #read_lon
        elif count == 3:
            deg = row[0]
            deg = deg[11:]
            deg = deg.strip()
            
            min = row[1]
            
            dir = row[2].strip()

            met_lons_degs.append(deg)
            met_lons_mins.append(min)
            met_lons_dirs.append(dir)
            
        #read alt
        elif count == 4:
            line = row[0].strip()
            line = line[24:-1]
            line = line.strip()
            met_alts.append(line)
            
        #read site classification
        elif count  == 5:
            line = row[0].strip()
            line = line[26:]
            line = line.strip()
            met_class.append(line)

        #read timezone
        elif count  == 6:
            line = row[0]
            line = line[10:]
            line = line.strip()
            met_tz.append(line)
 
        #read AQ CHECK
        elif count  == 7:
            line = row[0]
            line = line[4:]
            line = line.strip()
            met_aq_check.append(line)
            
        #read species and instrument
        elif count  == 8:
            line = row[0]
            line = line[16:]
            line_list = line.split(';')
            line_list = np.array(line_list)
            
            if len(line_list) > 1:
                len_list = len(line_list)
                species_inds = range(0,len_list,2)
                instrument_inds = range(1,len_list,2)
                
                met_species.append(line_list[species_inds])
                met_instruments.append(line_list[instrument_inds])
            else:
                met_species.append('na')
                met_instruments.append('na')
            
        #read country
        elif count  == 9:
            line = row[0]
            line = line[9:]
            line = line.strip()
            met_country.append(line)
            
        elif count == 10:
            count = -1
            
        count+=1

#convert lats, lons to -90:90,-180,180 from degree,minutes,seconds
def conversion(all_degs,all_mins,all_dirs):
    conv_array = []
    direction = {'N':1, 'S':-1, 'E': 1, 'W':-1}
    for i in range(len(all_degs)):
        deg = all_degs[i]
        min = all_mins[i]
        sec = 0
        current_dir = all_dirs[i]
        dir = direction[current_dir]
        
        conv =  (np.int64(deg)+np.int64(min)/60.0+sec/3600.0) * dir
        conv = np.around(conv,5)
        conv_array=np.append(conv_array,conv)
    
    return conv_array
        
met_lats = conversion(met_lats_degs,met_lats_mins,met_lats_dirs)
met_lons = conversion(met_lons_degs,met_lons_mins,met_lons_dirs)  

#check sites that have aq data
del_list = np.where(met_aq_check == 'False')
del_list = del_list[0]
alt_meta_refs = np.delete(met_refs,del_list)
valid_refs = [x for x in valid_refs if x in alt_meta_refs]
print 'n refs after aq present check = ',len(valid_refs)

#check sites that have species data for
del_list = []
for i in range(len(met_refs)):
    if species not in met_species[i]:
        del_list.append(i)
alt_meta_refs = np.delete(met_refs,del_list)
valid_refs = [x for x in valid_refs if x in alt_meta_refs]
print 'n refs after species check = ',len(valid_refs)

#remove refs that have no datafiles for species
del_list = []
for i in range(len(valid_refs)):
    ref = valid_refs[i]
    s_files = insensitive_glob('/work/home/db876/observations/surface/%s/EANET/*%s.csv'%(fname_species,ref))
    site_files = []
    for y in year_array:
        for f in s_files:
            if str(y)[-2:] in f:
                site_files.append(f)
    if len(site_files) == 0:
        del_list.append(i)
valid_refs = np.delete(valid_refs,del_list)
valid_refs = sorted(valid_refs)
print 'n refs after checking have files for refs = ',len(valid_refs)

#change refs to be created ref not sitename
new_refs = []
for v in valid_refs:
    new_refs.append(file_refs[met_refs.index(v)])
valid_refs = new_refs

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

    ref = valid_refs[c]
    print 'ref = ',ref,c
    
    #get site instrument for species
    met_i = file_refs.index(ref)
    file_name = met_refs[met_i]
    site_name = met_sitenames[met_i]
    print site_name
    site_species = list(met_species[met_i])
    print site_species
    site_instruments = list(met_instruments[met_i])
    m_method = site_instruments[site_species.index(species)]

    
    site_resolutions = []
    data_valid = True

    s_files = insensitive_glob('/work/home/db876/observations/surface/%s/EANET/*%s.csv'%(fname_species,file_name))
    site_files = []
    for y in year_array:
        for f in s_files:
            if str(y)[-2:] in f:
                site_files.append(f)
                           
    site_files = modules.natsorted(site_files)
    
    years = []
    months = []
    days = []
    hours = []

    vals = []
    yyyymmdd = []
    hhmm = []
    
    n_dup_array = []
    
    last_year_index = len(site_files)
    for y in year_array:
        got_year = False
        for file in site_files:
            last_file_split = file.split('/')[-1]
            if str(y)[2:] in last_file_split:
                got_year = True
                break
        if got_year == False:
            timedelta_diff = datetime.date(y+1, 1, 1) - datetime.date(y, 1, 1)
            ndays_missing = timedelta_diff.days
            continue
            
        print file
        
        valid = True
        with open(file, 'rb') as f:
            reader = csv.reader(f,delimiter=',')
            counter = 0
            
            #get resolution
            for row in reader:
                if counter == 0:
                    all_units = row
            
                elif counter == 1:   
                    file_res = 'H'
            
                    try:
                        hour_index = row.index('Hour')
                    except:
                        file_res = 'D'
                    try:
                        day_index = row.index('Day')
                    except:
                        file_res = 'M'
                    month_index = row.index('Month')
                    year_index = row.index('Year')
                    
                    try:
                        spec_index = row.index(species.upper())
                        unit = all_units[spec_index] 
                    except:
                        valid = False
                        break
                    
                    #make sure each year units are ppb
                    if unit != 'ppb':
                        print 'Units not ppb!'
                        1+'a'
                        
                if counter == 2:
                    if file_res == 'H':
                        yyyy = row[year_index]
                        mm = row[month_index]
                        dd = row[day_index]
                        hh = row[hour_index]
                    elif file_res == 'D':
                        yyyy = row[year_index]
                        mm = row[month_index]
                        dd = row[day_index]
                        hh = 1
                    elif file_res == 'M':
                        yyyy = row[year_index]
                        mm = row[month_index]
                        dd = 1
                        hh = 1
        
                    start_datetime = datetime.datetime(int(yyyy),int(mm),int(dd),int(hh))
                
                if counter == 3:
                    if file_res == 'H':
                        yyyy = row[year_index]
                        mm = row[month_index]
                        dd = row[day_index]
                        hh = row[hour_index]
                    elif file_res == 'D':
                        yyyy = row[year_index]
                        mm = row[month_index]
                        dd = row[day_index]
                        hh = 1
                    elif file_res == 'M':
                        yyyy = row[year_index]
                        mm = row[month_index]
                        dd = 1
                        hh = 1
            
                    present_datetime = datetime.datetime(int(yyyy),int(mm),int(dd),int(hh))
                
                    time_delt = present_datetime-start_datetime
                    hour_delt = datetime.timedelta(hours=1)
                    day_delt = datetime.timedelta(hours=24)
                    week_delt = datetime.timedelta(hours=24*7)
                    month_delt = datetime.timedelta(hours=24*28)
                
                    print time_delt
            
                    if (time_delt < day_delt):
                        print 'Hourly Data'
                        file_res = 'H'
                        site_resolutions.append(file_res)
                
                    elif (time_delt > hour_delt) & (time_delt < week_delt):
                        print 'Daily Data'
                        file_res = 'D'
                        site_resolutions.append(file_res)
                
                    elif (time_delt > week_delt):
                        print 'Monthly Data'
                        file_res = 'M'
                        site_resolutions.append(file_res)       
                                    
                counter+=1
  
        #READ IN DATA   
        if valid == True:     
            #limit to sites with hourly date files for, if required
            if output_res == 'H':
                if file_res != 'H':
                    print 'Not processing as only want hourly files'
                    continue
            if output_res == 'HD':
                if file_res == 'M':
                    print 'Not processing as only want hourly and daily files'
                    continue
            with open(file, 'rb') as f:       
                reader = csv.reader(f,delimiter=',')
                counter = 0
                val_count = 0
                for row in reader:
            
                    if counter >= 2:
                        yyyy = row[year_index]
                        mm = row[month_index]   
                        
                        #add to n_obs_all
                        n_all+=1
                        n_after_nometa+=1
                        
                        if file_res == 'H':
                            try:
                                vals=np.append(vals,np.float64(row[spec_index]))
                            except:
                                vals=np.append(vals,-99999)
                            
                            current_datetime = present_datetime+relativedelta(hours=val_count)
                            yyyymmdd.append(current_datetime.strftime("%Y%m%d"))
                            hhmm.append(current_datetime.strftime("%H%M"))
                            n_dup_array=np.append(n_dup_array,0)
                            
                        elif file_res == 'D':
                            try:
                                vals=np.append(vals,[np.float64(row[spec_index])]*24)
                            except:
                                vals=np.append(vals,[-99999]*24)
                            
                            current_datetime = present_datetime+relativedelta(days=val_count)
                            next_datetime = present_datetime+relativedelta(days=val_count+1)
                            all_datetimes = pd.date_range(current_datetime,next_datetime,freq='H')[:-1]
                            for d in all_datetimes:
                                yyyymmdd.append(d.strftime("%Y%m%d"))
                                hhmm.append(d.strftime("%H%M"))  
                            n_dup_array=np.append(n_dup_array,0)
                            n_dup_array=np.append(n_dup_array,[1]*23)

                    
                        elif file_res == 'M':
                            month_days = monthrange(int(yyyy), int(mm))[1]
                            try:
                                vals=np.append(vals,[np.float64(row[spec_index])]*(month_days*24))
                            except:
                                vals=np.append(vals,[-99999]*(month_days*24))
                            
                            current_datetime = present_datetime+relativedelta(months=int(mm)-1)
                            next_datetime = present_datetime+relativedelta(months=int(mm))
                            all_datetimes = pd.date_range(current_datetime,next_datetime,freq='H')[:-1]
                            for d in all_datetimes:
                                yyyymmdd.append(d.strftime("%Y%m%d"))
                                hhmm.append(d.strftime("%H%M"))
                            n_dup_array=np.append(n_dup_array,0)
                            n_dup_array=np.append(n_dup_array,[1]*((month_days*24)-1))
                        
                        val_count+=1    
                    counter+=1
                        
        else:
            print 'Species is not in file header. Skipping Year'
            timedelta_diff = datetime.date(y+1, 1, 1) - datetime.date(y, 1, 1)
            ndays_missing = timedelta_diff.days
            print 'ndays missing = ', ndays_missing          
    
    #test if have no data due to not required time resolution, if so exit
    if len(vals) == 0:
        n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod = 0,0,0,0,0,0,0,0,0,0,0,0,0
        exit_c_list = np.array([inv_nometa,inv_anyvaliddata,inv_nokeymeta,inv_resolution,inv_badmeasurementmethod])
        n_c_list = np.array([n_all,n_after_nometa,n_after_flagsandlod,n_after_duplicate,n_after_anyvaliddata,n_after_nokeymeta,n_after_resolution,n_after_badmeasurementmethod])
        unknown_list = [unknown_mm_list,unknown_mm_refs_list,unknown_local_tz_list]
        meta = ['na','na','na','na','na','na','na','na','na','na','na','na']
        return c,['na'],['na'],['na'],False,meta,exit_c_list,n_c_list,unknown_list,'nothourly',np.zeros(0)
    
    #create max possible grid
    full_data = np.empty(n_hours)
    full_data_after_flagsandlod = np.empty(n_hours)
    big_n_dup_array = np.zeros(n_hours)
    full_data[:] = -99999
    full_data_after_flagsandlod[:] = -99999
    
    #convert blank values to -99999
    test_inv = vals == ''
    vals[test_inv] = -99999
    
    #convert number invalids to -99999
    test_inv = vals < 0
    vals[test_inv] = -99999
    
    #if all site resolutions are same continue then take first file_res
    all_same = all(x == site_resolutions[0] for x in site_resolutions)
    if all_same == True:
        file_res = site_resolutions[0]
    else:
    #otherwise take lowest frequency res as file_res 
        if 'M' in site_resolutions:
            file_res = 'M'
        elif 'D' in site_resolutions:
            file_res = 'D'
        else:
            file_res = 'H'
    
    #get meta
    i_ref = file_refs.index(ref)
    site_ref = ref
    data_tz = np.float32(met_tz[i_ref])
    all_tz = [data_tz] 
    lat = np.float32(met_lats[i_ref])
    lon = np.float32(met_lons[i_ref])
    alt = np.float32(met_alts[i_ref])
    raw_class_name = met_class[i_ref]
    country = met_country[i_ref]
    unit = str(unit)
    contact = 'Ayako Aoyagi, Asia Center for Air Pollution Research, eanetdata@acap.asia'
    
    
    #adjust dates and times if tz is not equal to 0
    tz = int(data_tz)
    if tz != 0:
        for i in range(len(yyyymmdd)):
            #create datetime
            dt = datetime.datetime(int(yyyymmdd[i][:4]),int(yyyymmdd[i][4:6]),int(yyyymmdd[i][6:]),int(hhmm[i][:2]),int(hhmm[i][2:]))
            if tz > 0:
                dt  = dt - datetime.timedelta(hours = int(tz))
            elif tz < 0:
                dt  = dt + datetime.timedelta(hours = np.abs(int(tz)))
            yyyymmdd[i] = dt.strftime("%Y%m%d")
            hhmm[i] = dt.strftime("%H%M")
        
    #put vals into full grid
    date_con = np.array(yyyymmdd).astype(int)
    time_con = np.array(hhmm).astype(int)
    
    #remove data < 1970 and >= 2015
    test_inds = (date_con >= 19700101) & (date_con < 20150101)
    date_con = date_con[test_inds]
    time_con = time_con[test_inds]
    vals = vals[test_inds]
    n_dup_array = n_dup_array[test_inds]
    
    #set st_big and mm_big
    st_big = ['continuous']*len(vals)
    mm_big = [m_method]*len(vals)
    
    #get obs valid  
    test = vals >= 0
    valid_hours_dup = np.sum(n_dup_array[test])
    n_obs_valid = int(len(vals[test]) - valid_hours_dup)
    n_after_flagsandlod += n_obs_valid
    
    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
    converted_time = modules.date_process(date_con,time_con,start_year)
    converted_time = np.round(converted_time,decimals=5)
    syn_grid_time = np.arange(0,n_days,1./24)
    syn_grid_time = np.round(syn_grid_time,decimals=5)
    raw_indices = np.searchsorted(syn_grid_time, converted_time, side='left')
    vals = np.array(vals)
    full_data_after_flagsandlod[raw_indices] = vals
    raw_st = np.copy(st_big)
    raw_mm = np.copy(mm_big)
    
    # test and remove duplicate and overlap points
    converted_time,vals,mm_big,st_big,n_dup_array = modules.remove_duplicate_points(site_ref,converted_time,vals,mm_big,st_big,n_dup_array,output_res)
    test = vals >= 0
    valid_hours_dup = np.sum(n_dup_array[test])
    n_obs_valid = int(len(vals[test]) - valid_hours_dup)
    n_after_duplicate += n_obs_valid
    
    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
    indices = np.searchsorted(syn_grid_time, converted_time, side='left')
    full_data[indices] = vals 
    big_n_dup_array[indices] = n_dup_array
        
    key_meta = [lat,lon,alt]
    
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

