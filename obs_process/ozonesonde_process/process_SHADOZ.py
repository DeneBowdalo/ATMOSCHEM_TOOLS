import glob
import numpy as np
import csv
import datetime
from collections import OrderedDict
from netCDF4 import Dataset
from collections import Counter
import modules
from time import sleep
import pandas as pd
from collections import defaultdict
from netCDF4 import num2date, date2num

def date_process(date,time,start_year):
    year=(date//10000)
    month=((date-year*10000)//100)
    day=(date-year*10000-month*100)

    hour=(time//10000)
    min=((time-hour*10000)//100)
    sec=(time-hour*10000-min*100)

    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),np.int(sec[i]))- \
              datetime.datetime(int(start_year),1,1,0,0,0) \
              for i in range(len(year))]

    processed_dates=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
    
    return np.array(processed_dates)

species = 'O3'

#data range
start_year = 2009
end_year = 2011
start_year_file = 2005
end_year_file = 2012

#load in model pressure edges
data = Dataset('/work/home/db876/plotting_tools/model_files/GEOSCHEM_VERT_%s_%s_%s_v902_2x2.5_GEOS5_D_*.nc'%(species,start_year_file,end_year_file)) 
model_press_edges = data.variables['pressure_edges'][:]
model_time = data.variables['time'][:]
model_lat_centres = data.variables['lat_centre'][:]
model_lat_edges = data.variables['lat_edges'][:]
model_lon_centres = data.variables['lon_centre'][:]
model_lon_edges = data.variables['lon_edges'][:]

#cut time
start_time = date2num(datetime.datetime(int(start_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
end_time = date2num(datetime.datetime(int(end_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
start_ind = np.searchsorted(model_time,start_time)
end_ind = np.searchsorted(model_time,end_time)
if (start_time < model_time[0]) & (end_time > model_time[-1]):
    model_time = model_time[:]
    model_press_edges = model_press_edges[:,:,:,:]
elif start_time < model_time[0]:
    model_time = model_time[:end_ind]
    model_press_edges[:end_ind,:,:,:]
elif end_time > model_time[-1]:
    model_time = model_time[start_ind:]
    model_press_edges[start_ind:,:,:,:]
else:
    model_time = model_time[start_ind:end_ind]
    model_press_edges[start_ind:end_ind,:,:,:]

ave_model_press_edges = np.average(model_press_edges,axis=0)
ave_model_press_centres = np.empty((47,len(model_lat_centres),len(model_lon_centres)))
for lat_i in range(len(model_lat_centres)):
    for lon_i in range(len(model_lon_centres)):
        for alt_i in range(47):
            ave_model_press_centres[alt_i,lat_i,lon_i] = (ave_model_press_edges[alt_i,lat_i,lon_i] + ave_model_press_edges[alt_i+1,lat_i,lon_i]) / 2.

all_val_years_int = range(start_year,end_year)
all_val_years_int = np.array(all_val_years_int)
all_val_years = [str(i) for i in all_val_years_int]

year_range = (end_year) - start_year

#setup netcdf output
root_grp = Dataset('%s_RADIOSONDES_SHADOZ_%s_%s.nc'%(species,start_year,end_year), 'w')
root_grp.description = 'SHADOZ Radiosondes of %s at sites in ppb - Program written by Dene Bowdalo'%(species)

site_count = 0

#---------------------------------------------------------------------------------
#process SHADOZ data

print '\nProcessing SHADOZ data\n'

files = glob.glob('/work/home/db876/observations/ozonesonde/SHADOZ/*')
files=modules.natsorted(files)

site_names = []
#separate out sites into each location
for i in files:
    i = i.replace('/work/home/db876/observations/ozonesonde/SHADOZ/','')
    split = i.split('_')
    site_names.append(split[0])
site_names = np.sort(list(set(site_names)))

for site_name in site_names:
    print '\n'
    print site_name
    
    data_valid = True
    
    file_n = 0
    
    all_var = []
    all_press = []
    datetimes = []
    yyyymmdd = []
    hhmmss = []
    
    site_files = []
    files = sorted(glob.glob('/work/home/db876/observations/ozonesonde/SHADOZ/%s*'%(site_name)))
    for y in all_val_years_int:
        for f in files:
            a = f.split('_')[1][:4] 
            if str(y) == a:
                site_files.append(f)
        

    for j in site_files:
        skip_n = -1
        
        prev_time = -1
        current_time = -1
        
        file_var = []
        file_press = []
        file_time = []
        
        with open(j, 'rb') as f:
            reader = csv.reader(f,delimiter=' ')
            counter = 0
            for row in reader:
                
                row = [var for var in row if var]
                
                if counter == 0:
                    skip_n = int(row[0])
                
                if counter > 0:
                    if (row[0] == 'Latitude') & (file_n == 0):
                        current_lat = row[-1]
                        if current_lat[0] == '+':
                            current_lat = current_lat[1:]
                        print 'lat = %s'%(current_lat)
                    if (row[0] == 'Longitude') & (file_n == 0):
                        current_lon = row[-1]
                        if current_lon[0] == '+':
                            current_lon = current_lon[1:]
                        print 'lon = %s'%(current_lon) 
                    if (row[0] == 'Launch') & (row[1] == 'Date'):
                        current_launch_date = row[-1]
                        #print current_launch_date
                    if (row[0] == 'Launch') & (row[1] == 'Time'):
                        current_launch_time = row[-1]
                        #print current_launch_time
                        s_yyyy = current_launch_date[:4]
                        s_mm = current_launch_date[4:6]
                        s_dd = current_launch_date[6:8]
                        
                        try:
                            s_hh = current_launch_time[0:2]
                            s_min = current_launch_time[3:5]
                            test = int(s_hh)
                        except:
                            s_hh = current_launch_time[0]
                            s_min = current_launch_time[2:]
                        
                        #test if dates/times contained are actually contained.
                        if (s_yyyy == '9000') or (s_yyyy == '9999') or (s_hh == '90') or (s_hh == '99') :
                            print 'file date/times are not valid in %s'%(j)
                            break
                        
                        start_datetime = datetime.datetime(int(s_yyyy),int(s_mm),int(s_dd),int(s_hh),int(s_min))
                                
        
                    if counter == skip_n-1:
                        #get number of params there should be on row
                        param_n = len(row)
                        time_ind = row.index('sec')
                        press_ind = row.index('hPa')
                        var_ind = row.index('ppmv')
                          
                    if counter >= skip_n:
                        #make sure that row contains numbers
                        try:
                            valid = float(row[0])
                        except:
                            #print 'first col not number'
                            continue
                        
                        #make sure number of columns is right
                        if len(row) != param_n:
                            continue
                        
                        prev_time = current_time
                        
                        current_time = row[time_ind]
                        current_var = row[var_ind]
                        current_press = row[press_ind]                        

                        #convert invalid numbers to -99999
                        if (current_time == '9000') or (current_time == ''):
                            current_time = -99999
                        if (current_var == '9000') or (current_var == ''):
                            current_var = -99999
                        if (current_press == '9000') or (current_press == ''):
                            current_var = -99999                        

                        if current_var != -99999: 
                            current_time = float(current_time)
                            current_var = float(current_var)
                            current_press = float(current_press)                            

                            #convert from ppmv to ppbv
                            current_var = current_var*1e3
                            
                            #if current time is less than 0 skip row
                            if (current_time < 0) & (current_time != -99999):
                                current_time = prev_time
                                continue
                            
                            #if current var is negative then skip row
                            if current_var < 0:
                                continue
                        
                            file_var.append(current_var)
                            file_time.append(current_time)
                            file_press.append(current_press)                        

                        else:
                            current_time = prev_time 
                                                
                    
                counter+=1     
    
        #check if there is valid data for file
        if len(file_var) > 0:
            file_n+=1
        
            file_var = np.array(file_var)
            file_time = np.array(file_time)
            file_press = np.array(file_press)
            
            #if all times are same values then must all be -99999, make assumption that all times go from base time in 1s intervals
            if all(x == file_time[0] for x in file_time) == True:    
                file_time = np.arange(0.,len(file_time),1.)
            #if there are some times in array that are -99999, but not all, 
            #then get indices of these rows and delete respective rows for all variables
            else:
                inv_test = file_time != -99999
                file_var = file_var[inv_test]
                file_time = file_time[inv_test]
                file_press = file_press[inv_test]
        
            #sort arrays ascending by time
            sort_inds = sorted(range(len(file_time)), key=lambda k: file_time[k])
            file_var = file_var[sort_inds]
            file_time = file_time[sort_inds]
            file_press = file_press[sort_inds]                

            #check for time duplicates( including data less than a second in resolution from previous point)
            #then remove duplicates
            dup_time = np.diff(file_time)
            #add a 1 to start of array to make length right
            dup_time = np.insert(dup_time,0,1.)
            dup_test = dup_time >= 1.
            file_var = file_var[dup_test] 
            file_time = file_time[dup_test]
            file_press = file_press[dup_test]
                
            #process times
            for time in file_time:
                calc_datetime = start_datetime + datetime.timedelta(seconds=time)
                
                datetimes.append(calc_datetime)
                yyyymmdd.append(calc_datetime.strftime("%Y%m%d"))
                hhmmss.append(calc_datetime.strftime("%H%M%S"))
                
            all_var.append(file_var)
            all_press.append(file_press)
             
        else:
            print 'File %s not valid. No valid data'%(j)    
    
    #make list out of list of lists
    all_var = [val for sublist in all_var for val in sublist]
    all_press = [val for sublist in all_press for val in sublist]
    
    #if there is at least some site data then write out otherwise don't.
    if len(all_var) == 0:
        print 'Site is not valid. No valid data.'
    else:
        all_var = np.array(all_var)
        all_press = np.array(all_press)        

        calc_time = date_process(np.array(yyyymmdd).astype('float'),np.array(hhmmss).astype('float'),start_year)
        
        #check that none of times are duplicated and all times go forward
        time_gaps = np.diff(calc_time)
        for i in range(len(time_gaps)):
            if time_gaps[i] == 0:
                print 'Dulicate times for site'
                #print sorted([k for k,v in Counter(calc_time).items() if v>1])
                c = 1 + 'a'
            if time_gaps[i] < 0:
                print 'Time not increasing for site'
                print i
                print calc_time[i-1], calc_time[i],calc_time[i+1] 
                print yyyymmdd[i-1],hhmmss[i-1],yyyymmdd[i],hhmmss[i],yyyymmdd[i+1],hhmmss[i+1]
                c = 1 + 'a'
        
        #check that data resolution is at least to 1ppb
        #do by checking checking differences with 1st valid point and all other valid points
        #exclude 0 differences, if minimum difference is > 1.1ppbv exclude site as data resolution is too poor

        min_diff = np.abs(np.diff(all_var))
        test = min_diff != 0
        min_diff = min_diff[test]
        try:
            min_diff = np.min(min_diff)    
            if min_diff > 1.1:
                data_valid = False
                print 'Data resolution is not sufficient. Needs to be at least to 1ppbv.'
                print 'min diff = ', min_diff
        except:
            data_valid = False
            print 'Data resolution is not sufficient. Needs to be at least to 1ppbv.'

        #check there are no persistent gaps > 3 months, more than half n of files
        #calc 3 months in days ~ 90 days
        #add start year and end year times to check gaps at start and end of array
        #if there is gap > year also make this invalid
        
        start_year_time = date_process(np.array([int('%s0101'%(start_year))]),np.array([000000]),start_year)
        end_year_time = date_process(np.array([int('%s1231'%(end_year-1))]),np.array([235959]),start_year)
        adj_time = np.insert(calc_time,0,0.0)
        adj_time = np.append(adj_time,end_year_time)
        
        time_gaps = np.diff(adj_time)
        
        inv_count = 0
        
        max_count = round(len(all_val_years_int)/2.)
        
        for i in time_gaps:
            if i > 90:
                inv_count+=1
                if inv_count >= max_count:
                    data_valid = False
                    print 'Persisent Data gap > 3 months'
                    break
            if i > 365:
                    data_valid = False
                    print 'Data gap > 1 Year'
                    break

        #calc gregorian time                                                                                                                                                  
        obs_time = date2num(datetimes, units='hours since 0001-01-01 00:00:00', calendar='gregorian')
                    
        if data_valid == True:  
            #split datasets into model pressure bands by daily timestep
           
            #find lat and lon indices of site
            lat_i,lon_i = modules.obs_model_gridbox(model_lat_edges,model_lon_edges,np.float64(current_lat),np.float64(current_lon))
           
            bad_inds = []
            level_inds = []

            for i in range(len(obs_time)): 
                obs_t = obs_time[i]
                obs_press = all_press[i]
                obs_var = all_var[i]

                time_ind = np.searchsorted(model_time,obs_t)
                if obs_t > model_time[-1]:
                    bad_inds.append(i)
                    continue
                else:
                    time_model_press_edges = model_press_edges[time_ind,:,:,:]            

                #find obs pressures between model pressure edges(pressures are reversed from highest to lowest so need to correct index to account for this)
                level_ind = np.searchsorted(time_model_press_edges[:,lat_i,lon_i][::-1],obs_press) 
                level_ind = 47 - level_ind 

                #if obs_pressure greater than top pressure edge then do not process this point
                if obs_press < time_model_press_edges[-1,lat_i,lon_i]:                                                                                                                                                                                    
                    bad_inds.append(i)
                    continue
                else:
                    if level_ind == -1:
                        level_ind = 0
                    level_inds.append(level_ind)

            #remove bad inds
            obs_time = np.delete(obs_time,bad_inds)
            all_var = np.delete(all_var,bad_inds)
    
            #calc gregorian time
            dt_times = date2num(datetimes, units='hours since 0001-01-01 00:00:00', calendar='gregorian')

            #put data into big array
            all_data = np.empty((len(all_var),47))
            all_data[:,:] = -99999

            for i,l in zip(range(len(all_var)),level_inds):
                all_data[i,l] = all_var[i] 
           
            country = ''
            while country == '':
                try:
                    country = modules.get_country(np.float64(current_lat),np.float64(current_lon))                       
                except:
                    pass

            # dimensions
            root_grp.createDimension('time%s'%(site_count), len(obs_time))    
            root_grp.createDimension('species%s'%(site_count), len(all_var))
            root_grp.createDimension('press%s'%(site_count), 47)

            #save out netcdf file
            ref = root_grp.createGroup('%s'%(site_name.lower()))

            #set variables
            times = ref.createVariable('time', 'f8', ('time%s'%(site_count),))
            spec = ref.createVariable(species.lower(), 'f4', ('species%s'%(site_count),'press%s'%(site_count)))
            p_c = ref.createVariable('pressure_centres', 'f4', ('press%s'%(site_count),))
 
            #set group attributes
            ref.latitude = np.float64(current_lat)
            ref.longitude = np.float64(current_lon)
            ref.process_group = 'SHADOZ'
            ref.country = country
            ref.unit = 'ppbv'
            print 'Measurement Type = ECC'
            ref.measurement_type = 'ECC'

            times[:] = obs_time
            spec[:] = all_data
            p_c[:] = ave_model_press_centres[:,lat_i,lon_i]
    
            print 'Site is Valid'
            
            site_count+=1
        else:
            print 'Site is Invalid'
    

