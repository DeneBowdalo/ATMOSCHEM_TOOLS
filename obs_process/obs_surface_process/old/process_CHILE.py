import urllib
import json
import numpy as np
import glob
import os
import datetime
from netCDF4 import Dataset
import modules
import csv
import matplotlib.pyplot as plt
from time import sleep

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

n_years = (end_year-start_year)+1

meta_data = np.genfromtxt('chile_meta.txt',delimiter=',',dtype=[('ref','S5'),('lat',float),('lon',float),('alt',float)],skip_header=2)

met_refs = meta_data['ref']
met_refs = met_refs.tolist()
met_lats = meta_data['lat']
met_lons = meta_data['lon']
met_alts = meta_data['alt']

# #get alts from lat,lon
# for lat,lon in zip(lats,lons):
#     ELEVATION_BASE_URL = 'http://maps.googleapis.com/maps/api/elevation/json?'
#     URL_PARAMS = "locations=%s,%s&sensor=%s" % (lat, lon, "false")
#     url=ELEVATION_BASE_URL + URL_PARAMS
#     f = urllib.urlopen(url)
#     response = json.loads(f.read().decode())
#     result = response["results"][0]
#     alt = float(result["elevation"])
#     alt = np.around(alt,2)
#     alts=np.append(alts,alt)
# data = np.vstack((refs,lats,lons,alts))
# data = np.transpose(data)
# np.savetxt('chile_meta_pro.txt',data,delimiter=',',fmt = '%s')

#rename files
# files = glob.glob('dat*')
# 
# count = 1
# for filename in os.listdir("."):
#     if filename.startswith("dat"):
#         if count < 10:
#             count = '0'+str(count)
#         os.rename(filename, 'CH0%s'%(count))
#         count = int(count)
#         count+=1

#find n_hours and n_days between start and end date
d0 = datetime.datetime(start_year, 1, 1,0)
d1 = datetime.datetime(end_year+1, 1, 1,0)
delta = d1 - d0
n_days = delta.days
n_hours = n_days*24

#setup netcdf output
root_grp = Dataset('CHILE_SURFACE_%s_%s_%s.nc'%(species,start_year,end_year+1), 'w')
root_grp.description = 'Hourly Surface %s at EANET sites in ppb - Program written by Dene Bowdalo'%(species)

# dimensions
root_grp.createDimension('date', n_hours)
root_grp.createDimension('time', n_hours)    
root_grp.createDimension('species', n_hours)

#check site is not urban using anthrome map from 2000
anthfile = '/work/home/db876/plotting_tools/core_tools/anthro2_a2000.nc'
print anthfile
anthload = Dataset(anthfile)
class_result,class_name = modules.anthrome_classify(anthload,met_lats.astype('float64'),met_lons.astype('float64'))
del_list = np.where(class_result == 'invalid')
del_list = del_list[0]
valid_refs = np.delete(met_refs,del_list)
print 'n refs after class remove = ',len(valid_refs)
print valid_refs

#create grid time array
datetime_obj = datetime.datetime(year = start_year, month = 1, day = 1, hour = 0, minute = 0) 
grid_dates = []
grid_times = []
all_datetimes = []
for i in range(n_hours):
    one_date = datetime_obj.strftime("%Y%m%d")
    one_time = datetime_obj.strftime("%H%M")
    grid_dates = np.append(grid_dates,one_date)
    grid_times = np.append(grid_times,one_time)
    datetime_obj = datetime_obj+datetime.timedelta(hours=1)
    all_datetimes.append(datetime_obj)

for ref in valid_refs:
    valid_flag = True
    valid_write = True
    vals = []
    date_con = []
    time_con = []
    date_range_count = 'na'
    start_read_count = 'na'
    print ref
    file = glob.glob('/work/home/db876/observations/surface/%s/CHILE/%s'%(species,ref))
    file = file[0]
    print file
    
    while (valid_flag == True):
        with open(file, 'rb') as f:
            reader = csv.reader(f,delimiter=',')
            counter = 0
            for row in reader:
                if counter == 0:
                    current_site = row[0][1:]
                try:
                    if row[0][2:10] == 'Database':
                        if row[0][-5:] != 'Horas':
                            valid_flag = False
                            valid_write = False
                        date_range_count = counter+2
                except:
                    pass

                if counter == date_range_count:
                    start_date_str = row[0][2:4] + row[0][4:6] + row[0][6:8] + row[0][9:11]
                    end_date_str = row[0][14:16] + row[0][16:18] + row[0][18:20] + row[0][21:23]

                try:
                    if row[0][:] == 'EOH':
                        print 'start reading'
                        start_read_count = counter+2
                except:
                    pass
            
                try:
                    if row[0][:] == 'EOF':
                        print 'end reading'
                        valid_flag = False
                        counter = -999
                except:
                    pass
            
            
                if counter >= start_read_count:
                    date = row[0]
                    time = row[1][1:]
                
                    #if time is 2400, change it to 0000, and add a day to datetime
                    if time == '2400':
                        time = '0000'
                        add_day = True
                    else:
                        add_day = False
                                
                    present_date_str = date+time
                    present_datetime = datetime.datetime.strptime(present_date_str,'%y%m%d%H%M')
                    if add_day == True:
                        present_datetime = present_datetime + datetime.timedelta(days=1)
                
                    date = present_datetime.strftime("%Y%m%d")
                
                    #print 'current datetime = ',present_datetime
                    #print 'start datetime = ',d0
                
                    if (present_datetime >= d0) & (present_datetime < d1):
                        #print 'yes'
                
                        date_con.append(date)
                        time_con.append(time)
                
                        off_data = row[2].strip()
                        prelim_data = row[3].strip()
                
                        if off_data != '':
                            vals.append(np.float64(off_data))
                        elif prelim_data != '':
                            vals.append(np.float64(prelim_data))
                        else:
                            vals.append(-99999)
            
                counter+=1
            
    #put vals into full grid
    date_con = np.array(date_con).astype(int)
    time_con = np.array(time_con).astype(int)

    #create max possible o3 grid
    full_data = np.empty(n_hours)
    full_data[:] = -99999
        
    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
    converted_time = date_process(date_con,time_con,start_year)
    converted_time = np.round(converted_time,decimals=5)
    syn_grid_time = np.arange(0,n_days,1./24)
    syn_grid_time = np.round(syn_grid_time,decimals=5)
    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
    indices = np.searchsorted(syn_grid_time, converted_time, side='left')
    vals = np.array(vals)
    full_data[indices] = vals
        
    #correct timezone to UTC
    tz = -4
    if tz < 0:
        #get rid of values at start and append -99999's at end
        cut = full_data[:tz]
        for num in range(np.abs(tz)):
            cut = np.insert(cut,0, -99999)
        full_data = cut
    elif tz > 0:
        #put -99999's at start and get rid of values at end
        cut = full_data[tz:]
        for num in range(tz):
            cut = np.append(cut, -99999)
        full_data = cut
                
    #check amount of valid data, if less than 50 %, skip.
    a_test = full_data != -99999
    valid = full_data[a_test]
    valid_line = len(full_data) /2
    if len(valid) < valid_line:
        valid_write = False
        print 'Site Invalid, more than 50% of data missing.'
    data_complete = (100./len(full_data)) * len(valid)

    #check there are no gaps > 365 days
    count = 0
    max_count = 0
    for i in full_data:
        if i == -99999:
            count+=1
        else:
            count = 0

        if count > max_count:
            max_count = count

    if max_count > 8750:
        valid_write = False
        print 'Site Invalid, gaps greater than 365 days in data.'
    
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
        valid_write = False
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
            valid_write = False
            print 'Data resolution is not sufficient. Needs to be at least to 1ppbv.'
            print 'min diff = ', min_diff
    except:
        valid_write = False
        print 'Data resolution is not sufficient. Needs to be at least to 1ppbv.'
        
    meta_index = met_refs.index(ref)
    lat = met_lats[meta_index]  
    lon = met_lons[meta_index]
    alt = met_alts[meta_index]
    raw_class = 'na'
    anth_class_name = class_name[meta_index]
                 
    #Remove sites above 1000m from sea level
    if np.float64(alt) >= 1000:
        valid_write = False
        print 'Site is over 1000m from sea level.'
                            
    if valid_write == True:
        country = ''
        while country == '':
            try:
                country = modules.get_country(np.float64(lat),np.float64(lon))                       
            except:
                pass
    
        #save out netcdf file
        net_ref = root_grp.createGroup('%s'%(ref.lower()))

        #set variables
        dates = net_ref.createVariable('date', 'i8', ('date',))
        times = net_ref.createVariable('time', 'i8', ('time',))
        spec = net_ref.createVariable(species.lower(), 'f8', ('species',))

        #set group attributes
        net_ref.latitude = np.float64(lat)
        net_ref.longitude = np.float64(lon)
        net_ref.altitude = np.float64(alt)
        net_ref.process_group = 'CHILE'
        net_ref.country = country
        net_ref.data_completeness = data_complete
        net_ref.anthrome_site_class = anth_class_name
        net_ref.raw_site_class = raw_class
        net_ref.unit = 'ppbv'

        dates[:] = grid_dates
        times[:] = grid_times
        spec[:] = full_data
        print 'Site is Valid'
        valid_write=False 
    else:
        print 'Site is Invalid'
                
                
