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
from scipy import stats

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
start_year = 2005
end_year = 2010
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
root_grp = Dataset('%s_RADIOSONDES_NDACC_%s_%s.nc'%(species,start_year,end_year), 'w')
root_grp.description = 'NDACC Radiosondes of %s at sites in ppb - Program written by Dene Bowdalo'%(species)

site_count = 0

#---------------------------------------------------------------------------------
#process NDACC data

print '\nProcessing NDACC data\n'

files = glob.glob('/work/home/db876/observations/ozonesonde/NDACC/*')
files=modules.natsorted(files)

sites = []
for f in files:
    sites.append(f.split('/')[-1])

valid_sites = []
for site in sites:
    site_files = glob.glob('/work/home/db876/observations/ozonesonde/NDACC/%s/ames/o3sonde/*'%(site))
    if len(site_files) > 0:
        valid_sites.append(site)
        
print valid_sites  
#valid_sites = valid_sites[4:]
valid_sites = ['natal']
for site in valid_sites:
    print '\n'
    print site
    
    data_valid = True
        
    all_var = []
    all_press = []
    datetimes = []
    yyyymmdd = []
    hhmmss = []
    
    lats = []
    lons = []

    site_files = glob.glob('/work/home/db876/observations/ozonesonde/NDACC/%s/ames/o3sonde/*'%(site))
    
    #sort files
    years = []
    time_joint = []
    for s in site_files:
        fname = (s.split('/')[-1])
        year = fname[2:4]
        month = fname[4:6]
        day =  fname[6:8]
        if int(year) < 50:
            year = '20'+year
        elif int(year) > 50:
            year = '19'+year
        
        time_joint.append(year+month+day)
    site_files = [y for (x,y) in sorted(zip(time_joint,site_files))]
    print site_files
    
    for s in site_files:
        print s
        
        fname = (s.split('/')[-1])
        year = fname[2:4]
        month = fname[4:6]
        day =  fname[6:8]
        if int(year) < 50:
            year = '20'+year
        elif int(year) > 50:
            year = '19'+year
        
        print year,month,day
        if year in all_val_years:
            comment_n = 1
            ozonesonde_n = 4
                
            skip_n = 999999
            start_n = 999999
            start_var_n = 999999
            press_ind = 999999
            time_ind = 999999
            var_ind = 999999
            meta_ind = 999999
            lat_ind = 999999
            lon_ind = 999999
            launchtime_ind = 999999
            format_code = 999999
            ind_var_n = 999999
            n_var_n = 999999
            inv_n = 999999
            
            got_press = False
            got_time = False
            got_var = False
            got_lat = False
            got_lon = False
            got_launchtime = False
            got_launchtime_hour = False
            got_launchtime_min = False
            got_meta = False
        
            minutes_flag = False
            
            inv_time = 'na'
            inv_press = 'na'
            inv_var = 'na'
            
            prev_time = -1
            current_time = -1
            
            file_var = []
            file_press = []
            file_time = []
            
            read = False
            
            with open(s, 'rb') as f:
                reader = csv.reader(f,delimiter=' ')
                counter = 0
                for row in reader:
                    row = [var for var in row if var]
                    print counter
                    print row
                
                    #get skip row n and nasa ames format
                    if counter == comment_n:
                        skip_n = int(row[0])
                        start_n = (counter+skip_n)+1
                        format_code = row[1]
                    
                    #get inds for nasa ames file format
                    if format_code == '2110':
                        ind_var_n = 9
                    elif format_code == '2160':
                        ind_var_n = 10
                    if format_code == '2110':
                        n_var_n = 11
                    elif format_code == '2160':
                        n_var_n = 12
                    if format_code == '2110': 
                        inv_n = 13
                    elif format_code == '2160':
                        inv_n = 14
                    
                    #check measurment type is ECC
                    if counter == ozonesonde_n:
                        if 'ecc' not in ''.join(row).lower():
                            #if line just says ozonesonde then assume ECC
                            if ''.join(row).lower() != 'ozonesonde':
                                1+'a'   
                    
                    #get independent variable numbers
                    #if (site == 'lauder') & (format_code == '2110'):
                    #    n_independent = 2
                    #else:
                    n_independent = 1
                    
                    if counter == ind_var_n:
                        if 'time' in ''.join(row).lower():
                            print 'time ind',0
                            print ''.join(row).lower()
                            1+'a'
                            time_ind = 0
                            got_time = True
                            inv_time = 'na'
                            #check time is in seconds
                            if ('seconds' not in ''.join(row).lower()) & ('[s]' not in ''.join(row).lower()) & ('(s)' not in ''.join(row).lower()):
                                if 'minutes' in ''.join(row).lower():
                                    minutes_flag = True
                                else:
                                    1+'a'
                        if 'pressure' in ''.join(row).lower():
                            print 'press ind',0
                            press_ind = 0
                            got_press = True
                            inv_press = 'na'
                            #check pressure is in hPa
                            if 'hpa' not in ''.join(row).lower():
                                1+'a'
                    
                    #get number of dependent variables in file
                    if counter == n_var_n:
                        n_dependent = int(row[0])
                        param_n = n_independent+n_dependent
                        var_start = counter+3
                    
                    #get invalid row
                    if counter == inv_n:
                        inv_row = row
                                               
                    #get dependent variable inds
                    
                    #pressure
                    if (''.join(row).lower() == 'pressure[hpa]') & (got_press == False):
                        print 'press ind',1
                        press_ind = (counter - var_start)+n_independent
                        got_press = True
                        inv_press = inv_row[press_ind-n_independent]
                    if (''.join(row).lower() == 'pressure(hpa)') & (got_press == False):
                        print 'press ind',2
                        press_ind = (counter - var_start)+n_independent 
                        got_press = True
                        inv_press = inv_row[press_ind-n_independent]
                    if (''.join(row).lower() == 'pressureatobservation(hpa),measuredusingaradiosonde,precision0.5hpa') & (got_press == False):
                        print 'press ind',3
                        press_ind = (counter - var_start)+n_independent 
                        got_press = True
                        inv_press = inv_row[press_ind-n_independent]
                    if (''.join(row).lower() == '(pres)pressure(hectopascals)') & (got_press == False):
                        print 'press ind',4
                        press_ind = (counter - var_start)+n_independent 
                        got_press = True
                        inv_press = inv_row[press_ind-n_independent]
                    
                    #time
                    if (''.join(row).lower() == 'timeafterlaunch(s)') & (got_time == False):
                        print 'time ind',1
                        time_ind = (counter - var_start)+n_independent
                        got_time = True
                        inv_time = inv_row[time_ind-n_independent]
                    
                    #ozone partial pressure
                    if (''.join(row).lower() == 'ozonepartialpressure(mpa)') & (got_var == False):
                        print 'var ind',0
                        var_ind = (counter - var_start)+n_independent
                        got_var = True
                        inv_var = inv_row[var_ind-n_independent]
                      
                    if (''.join(row).lower() == 'ozonepartialpressure[mpa]') & (got_var == False):
                        print 'var ind',1
                        var_ind = (counter - var_start)+n_independent  
                        got_var = True
                        inv_var = inv_row[var_ind-n_independent]
                    
                    if (''.join(row).lower() == 'ozonepartialpressure(mpa),derived') & (got_var == False):
                        print 'var ind',2
                        var_ind = (counter - var_start)+n_independent  
                        got_var = True
                        inv_var = inv_row[var_ind-n_independent]
                        
                    if (''.join(row).lower() == '(ozpp)ozonepartialpressure(nanobars)') & (got_var == False):
                        print 'var ind',3
                        var_ind = (counter - var_start)+n_independent  
                        got_var = True
                        inv_var = inv_row[var_ind-n_independent]
                    
                    
                    #get meta inds
                    
                    #meta start
                    if (''.join(row).lower() == 'numberoflevels') & (got_meta == False):
                        print 'start meta', 0
                        start_var_n = counter
                        got_meta = True
                        
                    if (''.join(row).lower() == 'numberofpressurelevels') & (got_meta == False):
                        print 'start meta', 1
                        start_var_n = counter
                        got_meta = True    
                    
                    if (''.join(row).lower() == 'numberoflevelsreported') & (got_meta == False):    
                        print 'start meta', 2
                        start_var_n = counter
                        got_meta = True
                        
                    if (''.join(row).lower() == '(nlev)numberofpressurelevels(x1)') & (got_meta == False):
                        print 'start meta', 3
                        start_var_n = counter
                        got_meta = True
                        
                    #latitude
                    if (''.join(row).lower() == 'latitudeofstation') & (got_lat == False):
                        print 'lat ind', 0
                        lat_ind = counter-start_var_n
                        got_lat = True
                    
                    if (''.join(row).lower() == 'latitudeofstation(degrees)') & (got_lat == False):  
                        print 'lat ind', 1
                        lat_ind = counter-start_var_n
                        got_lat = True
                    
                    if (''.join(row).lower() == 'latitudeofstation(decimaldegrees)') & (got_lat == False):
                        print 'lat ind', 2
                        lat_ind = counter-start_var_n
                        got_lat = True
                    
                    if (''.join(row).lower() == 'latitudeofstation(nispositive)(decimaldegrees)') & (got_lat == False):
                        print 'lat ind', 3
                        lat_ind = counter-start_var_n
                        got_lat = True
                
                    if (''.join(row).lower() == 'stationlatitude[decimaldegreesn]') & (got_lat == False):
                        print 'lat ind', 4
                        lat_ind = counter-start_var_n
                        got_lat = True
                      
                    if (''.join(row).lower() == 'latitudeofobservingstation(degrees)') & (got_lat == False):
                        print 'lat ind', 5
                        lat_ind = counter-start_var_n
                        got_lat = True  
                     
                    if (''.join(row).lower() == '(slat)stationlatitude(degreesn)') & (got_lat == False):
                        print 'lat ind', 6
                        lat_ind = counter-start_var_n
                        got_lat = True  
                        
                    
                        
                    #longitude
                    
                    if (''.join(row).lower() == 'stationlongitude[decimaldegreese](range:0.00-359.99)') & (got_lon == False):
                        print 'lon ind', 0
                        lon_ind = counter-start_var_n
                        got_lon = True
                        
                    if (''.join(row).lower() == 'longitudeofstation(eispositive)(decimaldegrees)') & (got_lon == False):    
                        print 'lon ind', 1
                        lon_ind = counter-start_var_n
                        got_lon = True
                        
                    if (''.join(row).lower() == 'eastlongitudeofstation(degrees)') & (got_lon == False):
                        print 'lon ind', 2
                        lon_ind = counter-start_var_n
                        got_lon = True
                    
                    if (''.join(row).lower() == 'eastlongitudeofstation(decimaldegrees)') & (got_lon == False):
                        print 'lon ind', 3
                        lon_ind = counter-start_var_n
                        got_lon = True
                        
                    if (''.join(row).lower() == 'longitudeofobservingstation(degrees)') & (got_lon == False):
                        print 'lon ind', 4
                        lon_ind = counter-start_var_n
                        got_lon = True  
                    
                    if (''.join(row).lower() == '(slon)stationlongitude(degreese)') & (got_lon == False):
                        print 'lon ind', 5
                        lon_ind = counter-start_var_n
                        got_lon = True      
                    
                    
                    #starttime
                 
                    if (''.join(row).lower() == 'launchtime(utdecimalhoursfrom0hoursondaygivenbydate)') & (got_launchtime == False):      
                        print 'launch time ind', 0
                        launchtime_ind = counter-start_var_n
                        got_launchtime = True
                        
                    if (''.join(row).lower() == 'launchtime(decimaluthoursfrom0hoursondaygivenbydate)') & (got_launchtime == False): 
                        print 'launch time ind', 1
                        launchtime_ind = counter-start_var_n
                        got_launchtime = True
                    
                    if (''.join(row).lower() == 'launchtime[decimaluthoursfrom0hoursondaygivenbydate]') & (got_launchtime == False): 
                        print 'launch time ind', 2
                        launchtime_ind = counter-start_var_n
                        got_launchtime = True 
                    
                    if (''.join(row).lower() == 'launchtime(uthoursfrom0hoursondaygivenbydate)') & (got_launchtime == False): 
                        print 'launch time ind', 3
                        launchtime_ind = counter-start_var_n
                        got_launchtime = True
                        
                    if (''.join(row).lower() == 'hour-launch(ut)') & (got_launchtime == False):
                        print 'launch time hour ind', 0
                        launchtime_hour_ind = counter-start_var_n
                        got_launchtime_hour = True
                       
                    if (''.join(row).lower() == 'minute-launch(ut)') & (got_launchtime == False): 
                        print 'launch time min ind', 0
                        launchtime_min_ind = counter-start_var_n
                        got_launchtime_min = True

                    
                    #get file sitename, metadata is line after this one
                    if format_code == '2110':
                        if 'note' in ''.join(row).lower():
                            print row,counter
                            meta_ind = counter+1
                            lat_ind = lat_ind+1
                            lon_ind =  lon_ind+1
                            launchtime_hour_ind = launchtime_hour_ind+1
                            launchtime_min_ind = launchtime_min_ind+1
                            
                    elif format_code == '2160':
                        row_site = ''.join(row)
                        row_site = row_site.lower().replace('-','')
                        row_site = row_site.replace('(','')
                        row_site = row_site.replace(')','')
                        row_site = row_site.replace("'",'')
                        row_site = row_site.replace('.','')
                        row_site = row_site.replace(',','')
                        row_site = ''.join(i for i in row_site if not i.isdigit())
                        if (counter > start_var_n) & (site.lower() == row_site):
                            print row,counter
                            meta_ind = counter+1
                        
                    #read meta using inds
                    if counter == meta_ind:
                        print row,counter
                        if got_launchtime == True:
                            starttime = row[launchtime_ind]
                        if (got_launchtime_hour == True) & (got_launchtime_min == True):
                            starttime_hour = row[launchtime_hour_ind]
                            starttime_min = row[launchtime_min_ind]
                            starttime = float(starttime_hour) + (float(starttime_min)*(1./60.))
                            got_launchtime = True
                            
                        lon = row[lon_ind]  
                        lat = row[lat_ind]
                        print 'Latitude = %s'%(lat)
                        print 'Longitude = %s'%(lon)
                        lons.append(lon)
                        lats.append(lat)
                        
                        print 'decimal hours from start = ',starttime
                    
                        start_datetime = datetime.datetime(int(year),int(month),int(day))+datetime.timedelta(hours=float(starttime))
                        print start_datetime
                        
                    #try to read datalines
                    if counter > (start_n):
                    
                        #make sure number of columns is right
                        if len(row) == param_n:
                            #make sure that all row are numbers
                            try:
                                valid = [float(i) for i in row]
                                #make sure all meta has been got
                                if (got_press==True) & (got_time==True) & (got_var==True) & (got_meta==True) & (got_lat==True) & (got_lon==True) & (got_launchtime==True):
                                    read = True
                                else:
                                    #print 'not all flags triggered'
                                    read = False    
                            except:
                                print 'conversion failed'
                                read=False
                        else:
                            print 'n variables not right'
                            read=False
                            
                        #if read is False and there is already data wrote out then delete what has been wrote
                        if (read == False) & (len(file_var) > 0):
                            print 'deleting appended data, started reading too soon'
                            file_var = []
                            file_time = []
                            file_press = []
                            prev_time = -1
                            current_time = -1
                    
                        #start reading data
                        if read == True:
                            if len(file_var) == 0:
                                print 'reading'
                                print press_ind,time_ind,var_ind
                                print inv_press,inv_time,inv_var
                                print row
                                
                            
                            prev_time = current_time
                    
                            current_time = row[time_ind]
                            current_var = row[var_ind]
                            current_press = row[press_ind]                        

                            #convert invalid numbers to -99999
                            if (current_time == inv_time) or (current_time == ''):
                                current_time = -99999
                            if (current_var == inv_var) or (current_var == '') or (float(current_press) == 0.0):
                                current_var = -99999
                            if (current_press == inv_press) or (current_press == '') or (float(current_press) == 0.0):
                                current_var = -99999                        

                            if current_var != -99999: 
                                current_time = float(current_time)
                                current_var = float(current_var)
                                current_press = float(current_press)                            
                        
                                #convert var from partial pressure(millipascals) to ppbv mixing ratio
                                #convert partial presure into hectopascals, 1 millipascal = 1 x 10-5 hpa
                                current_var = current_var*1e-5
                                #convert partial pressure into mixing ratio, partial pressure = mixing_ratio * pressure
                                #so mixing ratio = partial pressure/ pressure
                                current_var = current_var / current_press
                                #convert mixing ratio into ppbv
                                current_var = current_var*1e9
                        
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
        else:
            continue
    
        #check if there is valid data for file
        if len(file_var) > 0:
            file_var = np.array(file_var)
            file_time = np.array(file_time)
            file_press = np.array(file_press)
        
            inv_test = file_time != -99999
            file_var = file_var[inv_test]
            file_time = file_time[inv_test]
            file_press = file_press[inv_test]
            
            if minutes_flag == True:
                #get into seconds
                file_time = np.array(file_time)*60.
                
            #correct for issue where scaling fact of 100 applied to data (for no good reason), everything except time affected
            if float(file_press[0]) > 10000:
                print 'correcting 100x scale factor'
                file_var = np.array(file_var)/100.
                file_press = np.array(file_press)/100.
            
            #process times
            for time in file_time:
                calc_datetime = start_datetime + datetime.timedelta(seconds=time)
                
                datetimes.append(calc_datetime)
                yyyymmdd.append(calc_datetime.strftime("%Y%m%d"))
                hhmmss.append(calc_datetime.strftime("%H%M%S"))
            
            all_var.append(file_var)
            all_press.append(file_press)
         
        else:
            print 'File %s not valid. No valid data'%(s)    

    #make list out of list of lists
    all_var = [val for sublist in all_var for val in sublist]
    all_press = [val for sublist in all_press for val in sublist]
    
    #if there is at least some site data then write out otherwise don't.
    if len(all_var) == 0:
        print 'Site is not valid. No valid data.'
    else:
        all_var = np.array(all_var)
        all_press = np.array(all_press)  
        datetimes = np.array(datetimes)  
        
        #get mode of lat/lons
        lat = stats.mode(lats)[0][0]
        lon = stats.mode(lons)[0][0]

        calc_time = date_process(np.array(yyyymmdd).astype('float'),np.array(hhmmss).astype('float'),start_year)
    
        #sort arrays ascending by time
        sort_inds = sorted(range(len(calc_time)), key=lambda k: calc_time[k])
        all_var = all_var[sort_inds]
        calc_time = calc_time[sort_inds]
        all_press = all_press[sort_inds]                

        #check for time duplicates( including data less than a second in resolution from previous point)
        #then remove duplicates
        dup_time = np.diff(calc_time)
        #add a 1 to start of array to make length right
        dup_time = np.insert(dup_time,0,1.)
        dup_test = dup_time >= 1.
        all_var = all_var[dup_test] 
        calc_time = calc_time[dup_test]
        all_press = all_press[dup_test]
        datetimes = datetimes[dup_test]
    
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
            lat_i,lon_i = modules.obs_model_gridbox(model_lat_edges,model_lon_edges,np.float64(lat),np.float64(lon))
       
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
                    country = modules.get_country(np.float64(lat),np.float64(lon))                       
                except:
                    pass

            # dimensions
            root_grp.createDimension('time%s'%(site_count), len(obs_time))    
            root_grp.createDimension('species%s'%(site_count), len(all_var))
            root_grp.createDimension('press%s'%(site_count), 47)

            #save out netcdf file
            ref = root_grp.createGroup('%s'%(site.lower()))

            #set variables
            times = ref.createVariable('time', 'f8', ('time%s'%(site_count),))
            spec = ref.createVariable(species.lower(), 'f4', ('species%s'%(site_count),'press%s'%(site_count)))
            p_c = ref.createVariable('pressure_centres', 'f4', ('press%s'%(site_count),))

            #set group attributes
            ref.latitude = np.float64(lat)
            ref.longitude = np.float64(lon)
            ref.process_group = 'NDACC'
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

                        
                        