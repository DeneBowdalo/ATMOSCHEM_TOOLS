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

#data range
start_year = int(raw_input('Start Year?\n'))
end_year = int(raw_input('\nEnd Year?\n'))
  
species = raw_input('\nSpecies?\n')
def date_process(date,time,first_year):
    year=(date//10000)
    month=((date-year*10000)//100)
    day=(date-year*10000-month*100)

    hour=time//100
    min=(time-hour*100)

    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(first_year,1,1,0,0,0) \
              for i in range(len(year))]

    processed_dates=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
    return processed_dates

#read obs files
all_files = glob.glob('/work/home/db876/observations/surface/%s/EMEP/*'%(species))
all_files = modules.natsorted(all_files)

#get all refs
ref_list = []
valid_refs=[]

for i in range(len(all_files)):
    f = all_files[i]
    f = f.replace("/work/home/db876/observations/surface/%s/EMEP/"%(species), "")
    f = f[:7]
    ref_list.append(f)

refs = set(ref_list)
refs = sorted([i for i in refs])
refs = np.array(refs)
print 'all refs len = ', len(refs)

#check each refs has all years in range, with no gaps
year_array = np.arange(start_year,end_year+1)

for ref in refs:
    years = []
    years_there = []
    for file in all_files:
        if ref in file:
            file = file.replace("/work/home/db876/observations/surface/%s/EMEP/"%(species), "")
            year = int(file[8:12])
            years=np.append(years,year)
    for n in range(len(year_array)):
        if year_array[n] in years:  
            years_there = np.append(years_there,year_array[n])
    if np.array_equal(years_there, year_array) == True:
        valid_refs.append(ref)

print 'n refs after year rm = ',len(valid_refs)


#find n_hours and n_days between start and end date
d0 = datetime.date(start_year, 1, 1)
d1 = datetime.date(end_year+1, 1, 1)
delta = d1 - d0
n_days = delta.days
n_hours = n_days*24

#setup netcdf output
root_grp = Dataset('EMEP_SURFACE_%s_%s_%s.nc'%(species,start_year,end_year+1), 'w')
root_grp.description = 'Hourly Surface %s at EMEP sites in ppb - Program written by Dene Bowdalo'%(species)

# dimensions
root_grp.createDimension('date', n_hours)
root_grp.createDimension('time', n_hours)    
root_grp.createDimension('species', n_hours)

#create grid time array
start = datetime.datetime(year = start_year, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = end_year+1, month = 1, day = 1, hour = 0, minute = 0)
            
grid_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]
grid_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='H')][:-1]

met_refs = []
met_lats_degs = []
met_lats_mins = []
met_lats_secs = []
met_lats_dirs = []
met_lons_degs = []
met_lons_mins = []
met_lons_secs = []
met_lons_dirs = []
met_alts = []

#read in meta data
meta_f = '/work/home/db876/observations/surface/%s/process/EMEP_META.csv'%(species)
with open(meta_f, 'rU') as f:
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        if row[0] in valid_refs:
            met_refs.append(row[0])
            
            met_lats_degs.append(row[2])
            met_lats_mins.append(row[3])
            met_lats_secs.append(row[4])
            met_lats_dirs.append(row[5])
            
            met_lons_degs.append(row[6])
            met_lons_mins.append(row[7])
            met_lons_secs.append(row[8])
            met_lons_dirs.append(row[9])
            
            met_alts.append(row[10])

#convert lats, lons to -90:90,-180,180 from degree,minutes,seconds
def conversion(all_degs,all_mins,all_secs,all_dirs):
    conv_array = []
    direction = {'N':1, 'S':-1, 'E': 1, 'W':-1}
    for i in range(len(all_degs)):
        deg = all_degs[i]
        min = all_mins[i]
        sec = all_secs[i]
        current_dir = all_dirs[i]
        dir = direction[current_dir]
        
        conv =  (np.int64(deg)+np.int64(min)/60.0+np.int64(sec)/3600.0) * dir
        conv = np.around(conv,5)
        conv_array=np.append(conv_array,conv)
    
    return conv_array
    
met_lats = conversion(met_lats_degs,met_lats_mins,met_lats_secs,met_lats_dirs)
met_lons = conversion(met_lons_degs,met_lons_mins,met_lons_secs,met_lons_dirs)

#check site is not urban using anthrome map from 2000
anthfile = '/work/home/db876/plotting_tools/core_tools/anthro2_a2000.nc'
anthload = Dataset(anthfile)
class_result,class_name = modules.anthrome_classify(anthload,met_lats.astype('float64'),met_lons.astype('float64'))
del_list = np.where(class_result == 'invalid')
del_list = del_list[0]
alt_meta_refs = np.delete(met_refs,del_list)
valid_refs = [x for x in valid_refs if x in alt_meta_refs]
print 'n refs after class remove = ',len(valid_refs)

#read files site at a time
for ref_i in range(len(valid_refs)):

    site_ref = valid_refs[ref_i]
    print 'Current Ref is = ', valid_refs[ref_i]
    #find if sites have full valid range from start year and finishing in end year
    s_files = glob.glob('/work/home/db876/observations/surface/%s/EMEP/%s*'%(species,site_ref))
    
    year_files = [file.replace("/work/home/db876/observations/surface/%s/EMEP/"%(species), "") for file in s_files]
    year_files = [file[8:12] for file in year_files]
    site_files = []
    for y in year_array:
        for i in range(len(s_files)):
            if str(y) in year_files[i]:
                site_files.append(s_files[i])
                           
    #site_files = modules.natsorted(site_files)
    site_files = modules.natsorted(site_files)

    #get lat,lon, altitude, anthrome class
    meta_index = met_refs.index(site_ref)
    latitude = met_lats[meta_index]
    longitude = met_lons[meta_index]
    altitude = met_alts[meta_index]
    anth_class_name = class_name[meta_index]
    
    full_exit = False

    yymmdd = []
    hhmm = []

    vals = []
    flags = []

    #create max possible o3 grid
    full_data = np.empty(n_hours)
    full_data[:] = -99999

    for file_i in range(len(site_files)):
        if full_exit == False:
            data_start = 9999999
            count = 0
            start_read = False
   
            with open(site_files[file_i], 'rb') as f:
                reader = csv.reader(f,delimiter=' ')
                print site_files[file_i]
                for row in reader:
                    try:
                        row = filter(lambda a: a != '', row)
                    except:
                        pass
                    try:
                        row = filter(lambda a: a != ',', row)
                    except:
                        pass
                    
                    #print row
                    
                    #get start date of file
                    if row[0] == 'Startdate:':
                        data = row[1]
                        s_yyyy = data[:4]
                        s_mm = data[4:6]
                        s_dd = data[6:8]
                        
                    #get timezone  
                    #if row[0] == 'Timezone:':
                    #    timezone = row[1]
                    #    print timezone
            
                    #get latitude
                    #try:
                    #    if row[1] == 'latitude:':
                    #        latitude = row[2]
                    #except:
                    #    pass
                    
                    #get longitude
                    #try:
                    #    if row[1] == 'longitude:':
                    #        longitude = row[2]
                    #except:
                    #    pass
                                   
                    #get altitude
                    #try:
                    #    if row[1] == 'altitude:':
                    #        altitude = row[2]
                    #        altitude = altitude[:-1]
                    #except:
                    #    pass
                         
                    #get unit
                    if row[0] == 'Unit:':
                        unit = row[1]
                
                    #get data
                    if start_read == True:
                        #calc dates, times, and take o3 vals
                        start_datetime = datetime.datetime(int(s_yyyy),1,1,0)
                            
                        time_since_start = np.float64(row[0])
                        days_since_start = math.trunc(time_since_start)
                        remainder = time_since_start - days_since_start
                        unrounded_hour = remainder*24
                        hour = np.round(unrounded_hour)
                        time_delta = datetime.timedelta(days = days_since_start,hours = hour)
                        calc_datetime = start_datetime + time_delta
                
                        calc_year = calc_datetime.strftime("%Y")
                        calc_yymmdd = calc_datetime.strftime("%Y%m%d") 
                        calc_hhmm = calc_datetime.strftime("%H%M")
                
                        yymmdd.append(calc_yymmdd)
                        hhmm.append(calc_hhmm)
                        vals = np.append(vals,np.float64(row[2]))
                        flags = np.append(flags,np.float64(row[3]))
            
                    
                    if row[0] == 'start_time':
                        start_read = True
                    
                    count+=1
            
               #do checks to make sure data goes fully from start year to end year
               # if file_i == 0:
               #     if int(calc_year) < end_year:
               #         print 'Data ends before final year'
               #         full_exit = True
        
        
                #if file_i == len(site_files) - 1:
                #if file_i == 1:    
                #    if int(s_yyyy) > start_year:
                #        print 'Data starts after first year'
                #        full_exit = True
                
                if (full_exit == False) & (file_i == len(site_files)-1):    
                    #check units for each site file:
                    if unit == 'ug/m3':
                        #print 'converting units, temp = 20degC'
                        #calculate conversion factor from mg/m3 assuming 20 degC and 1 atm - default for O3 instruments
                        #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                        conv_fact = 8.3144/48*(273.15+20)/(1013.25/10)
                        vals = conv_fact*vals
                        print 'Converting Units from ug/m3 20degC to ppbv'
            
                    #convert all invalids by flags to -99999
                    test_inv = flags != 0
                    
                    if len(test_inv) != 0:
                        vals[test_inv] = -99999
                    
                    #do additional invalid test, as flags not always correct
                    test_inv_2 = vals > 300
                    vals[test_inv_2] = -99999
            
                    #put o3 vals into full grid
                    date_con = np.array(yymmdd).astype(int)
                    time_con = np.array(hhmm).astype(int)
                    
                    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
                    converted_time = date_process(date_con,time_con,start_year)
                    converted_time = np.round(converted_time,decimals=5)
                    syn_grid_time = np.arange(0,n_days,1./24)
                    syn_grid_time = np.round(syn_grid_time,decimals=5)
                    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
                    indices = np.searchsorted(syn_grid_time, converted_time, side='left')
                    vals = np.array(vals)
                    #make sure no data is past end year
                    index_test = indices < len(full_data)
                    indices = indices[index_test]
                    vals = vals[index_test]
                    full_data[indices] = vals
                    
                    #check amount of valid data
                    a_test = full_data != -99999
                    valid = full_data[a_test]
                    valid_line = len(full_data) /2
                    if len(valid) < valid_line:
                        full_exit = True
                        print 'Site Invalid, more than 50% of data missing.' 
                    data_complete = (100./len(full_data)) * len(valid)
                        
                    #check for data gap > 1 year
                    count = 0
                    max_count = 0
                    for i in range(len(full_data)):
                        if full_data[i] == -99999:
                            count+=1
                            #print grid_dates[i]
                        else:
                            count = 0
                        
                        if count > max_count:
                            max_count = count
                            
                    if max_count >= 8750:
                        print 'Site Invalid, gaps greater than 365 days in data.'
                        full_exit = True
                        
                    #check there are no persistent gaps > 1 month, in at least half of years of files
                    year_inv_count = 0
                    year_range = range(start_year,end_year+1)
                    all_years = [i[:4] for i in grid_dates]
                    all_years = np.array(all_years)
                    for i in range(len(year_range)):
                        year_test = all_years == str(year_range[i])
                        year_vals = full_data[year_test]
                        max_count = 0
                        count = 0
                        for j in year_vals:
                            if j == -99999:
                                count+=1
                            else:   
                                count = 0
            
                            if count > max_count:
                                max_count = count
        
                        # if there is a gap > 28 days, add 1 to year_inv_count
                        if max_count > 672:
                            year_inv_count+=1
    
                    if year_inv_count >= 3:
                        full_exit = True
                        print 'Persisent Data gap in at least half of years > 28 days'
                    
                    
                    #check that data resolution is at least to 1ppb
                    #do by checking checking differences with 1st valid point and all other valid points
                    #exclude 0 differences, if minimum difference is > 1.1ppbv exclude site as data resolution is too poor
    
                    test = full_data != -99999 
                    val_data = full_data[test]  
                    min_diff = np.abs(np.diff(val_data))
                    test = min_diff != 0
                    min_diff = min_diff[test]
                    try:
                        min_diff = np.min(min_diff)    
                        if min_diff > 1.1:
                            full_exit = True
                            print 'Data resolution is not sufficient. Needs to be at least to 1ppbv.'
                            print 'min diff = ', min_diff
                    except:
                        full_exit = True
                        print 'Data resolution is not sufficient. Needs to be at least to 1ppbv.'
                    
                    #check timezone
                    #if timezone != 'UTC':
                    #    print 'Timezone is not UTC. Is %s'%(timezone)
                    #    fail = 1+'a' 
                    
                    #Remove sites above 1000m from sea level
                    if np.float64(altitude) >= 1000:
                        full_exit = True
                        print 'Site is over 1000m from sea level.'

    if full_exit == False:
        print 'Site is Valid'
        country = ''
        while country == '':
            try:
                country = modules.get_country(np.float64(latitude),np.float64(longitude))                       
            except:
                pass
    
        
    #save out netcdf file
        ref = root_grp.createGroup('%s'%(site_ref.lower()))

    #set variables
        dates = ref.createVariable('date', 'i8', ('date',))
        times = ref.createVariable('time', 'i8', ('time',))
        spec = ref.createVariable(species.lower(), 'f8', ('species',))

        #set group attributes
        ref.latitude = np.float64(latitude)
        ref.longitude = np.float64(longitude)
        ref.altitude = np.float64(altitude)
        ref.process_group = 'EMEP'
        ref.country = country
        ref.data_completeness = data_complete
        ref.anthrome_site_class = anth_class_name
        ref.raw_site_class = 'na'
        ref.unit = 'ppbv'

        dates[:] = grid_dates
        times[:] = grid_times
        spec[:] = full_data
    else:
        print 'Site is Invalid'
    
