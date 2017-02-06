from netCDF4 import Dataset
import numpy as np
import glob
import datetime
import pandas as pd
from netCDF4 import num2date, date2num

species_list = ['O3']
#species_list = ['O3']

#HOURLY HISTORICAL MODELS 2000 - 2015 (with 5 species)
models = ['GISSE2R']
#models = ['CESMCAM','CMAM','GEOSCCM','GFDLAM3','GISSE2R','MIROCCHEM','MOCAGE','UMCAM']
#start_year = 2000
#end_year = 2011

#HOURLY RCP2.6 & 8.5 MODELS 2095-2111
#models = ['CESMCAM','CMAM','GFDLAM3','GISSE2R','MIROCCHEM','MOCAGE','UMCAM']
start_year = 2095
end_year = 2110

#HOURLY RCP8.5em2000 MODELS 2095-2111
#models = ['CESMCAM','GFDLAM3','GISSE2R','MIROCCHEM','MOCAGE','UMCAM']
#start_year = 2095
#end_year = 2110

#choose  version of processing. hist, rcp26, rcp85 or rcp85em2000.
vers = 'rcp85em2000'

year_range = range(start_year,end_year+1)

if vers == 'hist':
    vers = 'historical'
    tag = '*'
elif vers == 'rcp26':
    tag = 'rcp26'
elif vers == 'rcp85':
    tag = 'rcp85'
elif vers == 'rcp85em2000':
    tag = 'rcp85em2000'  
#---------------------------------------------------
leap_years = range(1852,2200,4)
#If the year can be evenly divided by 100, it is NOT a leap year, unless the year is also evenly divisible by 400. Then it is a leap year
#Thus remove follwing years
leap_years.remove(1900)
leap_years.remove(2100)

#find n_hours and n_days between start and end date
d0 = datetime.datetime(start_year, 1, 1, 0, 0, 0)
d1 = datetime.datetime(end_year+1, 1, 1, 0, 0, 0)
delta = d1 - d0
n_days = delta.days
n_hours = n_days*24

#---------------------------------------------------  
for i in range(len(models)):
    model = models[i]
    print '\n------------------------------------'
    print model

    for y in range(len(species_list)):
        print species_list[y]
    
        valid_years = np.arange(start_year,end_year+1)

        files = glob.glob('%s/%s/%s/H/*'%(model,species_list[y],vers))
        files = sorted(files)
        
        if len(files) == 0:
            continue
        
        file_startdate = []
        for t in files:
            file_startdate.append(t.split('_')[-1][:4])

        valid_files = []
        for i in range(len(files)):
            for j in valid_years:
                if str(j) in file_startdate[i]:
                    if str(j-1) not in file_startdate[i]:            
                        valid_files.append(files[i])
                        print i

        leap_add_0 = -1
        leap_add_1 = 0
        leap_active = False
        has_leaps = False
        year_time = False
    
        #get start ref
        if (model == 'GFDLAM3') or (model == 'HADGEM2') or (model == 'NCARCAM3.5'):
            if vers == 'historical':
                start_ref = 2000
            if vers == 'rcp26':
                start_ref = 2100
            if vers == 'rcp85':
                start_ref = 2100
            if vers == 'rcp85em2000':
                start_ref = 2100
        elif (model == 'GISSE2R'):
            if vers == 'rcp85em2000':
                start_ref = 2050
            else:
                start_ref = 1850
        else:
            start_ref = 1850

        #---------------------------------------------------  
        #check if array has leap years and set flag accordingly
        print valid_files
        for count in range(len(valid_files)):
            root_grp = Dataset(valid_files[count])
            data = root_grp.variables['vmr%s'%(species_list[y].lower())][:]
            time = root_grp.variables['time'][:]
        
            #if first time is less than 365, time is referenced form start of year
            if time[0] < 365:
                year_time = True
            else:
                year_time = False
             
            #check if array has leap years and set flag accordingly 
            if has_leaps != True:
                if len(time) > 8760:
                    has_leaps = True
                    #calc start ref offset from start year
                    diff = start_year-start_ref
                    start_dt = datetime.datetime(start_ref,1,1,0,0)
                    end_dt = datetime.datetime(start_year,1,1,0,0)
                    diff_dt = end_dt-start_dt
                    year_diff = diff_dt.days
        
            #determine if time is on hour on on half hour
            if count == 0:
                dec = time[0] % 1  
                dec = np.around(dec,decimals=3)
                hour_interval = 1./24.
                hour_interval = np.around(hour_interval,decimals=3)
                if dec - hour_interval != 0:
                    hour = 'mid'            
                    print 'Time is on Half Hour'
                else:
                    hour = 'start'        
                    print 'Time is on Hour'
    
        if year_time == True:
            print '\nTime is referenced from start of year'
        else:
            print '\nTime is not referenced from start of year'
    
        if has_leaps == False:
            print 'Does not have leap years'
        else:
            print 'Has leap years'

        #---------------------------------------------------  
        #get data
    
        for count in range(len(valid_files)):

            root_grp = Dataset(valid_files[count])
            data = root_grp.variables['vmr%s'%(species_list[y].lower())][:]
            time = root_grp.variables['time'][:]
            
            #issue for GISSE2R r1i1p3 version sims, where 1st val of year is blank, thus change to -99999
            if (model == 'GISSE2R') & ('r1i1p3' in valid_files[count]):
                data[0,:,:] = -99999
        
            #get current year from file that are reading
            #current_year = year_range[count] 
            current_year = int(valid_files[count].split('_')[-1][:4])
    
            #set current days to be difference between start ref time and start file year time
            if count == 0:
                start_dt = datetime.datetime(start_year,1,1,0,0)
                end_dt = datetime.datetime(current_year,1,1,0,0)
                diff_dt = end_dt-start_dt
                current_days = diff_dt.days
    
            print '\nprocessing %s'%(valid_files[count])
            print 'len of year vals = ', len(time) 
        
            #print time
            #---------------------------------
            #Do Time corrections
        
            #process time year difference if files have no leap years (therefore time ref has no leap years)
            if has_leaps == False:
                year_diff = ((start_year-start_ref)*365)
        
            #if time is from historical start ref (1850) (or other time), take offset to make start year 0 and count from that, otherwise time is for each year from 0 
            if year_time == False:
                time = time-year_diff
                print 'historical time orig'
            else:    
                time = current_days + time
                print 'year time orig'
        
            #get current days processed  
            print current_year,current_days
            print time[0],time[-1]
            d0 = datetime.datetime(current_year, 1, 1, 0, 0)
            d1 = datetime.datetime(current_year+1, 1, 1, 0, 0)
            delta = d1 - d0
            current_days = current_days + delta.days

            #add leap years into year time array if they are not there
            if has_leaps == False:
                if current_year in leap_years:
                    print current_year, 'is leap year, altering year times accordingly, adding in day'
                    leap_add_0+=1
                    leap_add_1+=1
                    #59 is when days of year would equal Feb 29th on a leap year
                    test0 = time < 59
                    test1 = time >= 59
                    time[test0]=time[test0]+leap_add_0
                    time[test1]=time[test1]+leap_add_1
                    leap_active = True
                elif leap_active == True:
                    time=time+leap_add_1
                
            #print time
        
            #-----------------------------------
            #process grid
            if count == 0:
                lat_centre = root_grp.variables['lat'][:]
                lon_centre = root_grp.variables['lon'][:]

                if (model != 'GFDLAM3') & (model != 'OSLOCTM2') & (model != 'UMCAM'):
                    lat_bnds = root_grp.variables['lat_bnds'][:]
                    lon_bnds = root_grp.variables['lon_bnds'][:]

                    lat_edges = lat_bnds[:,0]
                    lon_edges = lon_bnds[:,0]              
        
                elif model == 'GFDLAM3':
                    lat_edges = np.arange(-90,89,2)                                                                                                                                                                                                               
                    lon_edges = np.arange(0,359,2.5) 

                elif model == 'OSLOCTM2':
                    lat_edges = [-90]
                    for c in range(len(lat_centre)-1):
                        lat_edges.append((lat_centre[c]+lat_centre[c+1])/2.)                                                                                                                                                                
                
                    lon_edges = [-1.40625]
                    for c in range(len(lon_centre)-1):
                        lon_edges.append((lon_centre[c]+lon_centre[c+1])/2.)
            
                elif model == 'UMCAM':
                    lat_bnds = root_grp.variables['lat_bnds'][:]
                    lon_bnds = root_grp.variables['lon_bnds'][:]

                    lat_edges = lat_bnds[:,0]
                    lon_edges = lon_bnds[:,0] 
                
                    lat_edges[0] = -90


                #add 90 to lat_band_eges
                lat_edges = np.append(lat_edges,90)

                #add wrapped lon band edge
                lon_edges = np.append(lon_edges,(360-np.abs(lon_edges[0])))
                ind = np.where(lon_edges==0)[0]
        
                lon_centre = lon_centre-180
                lon_edges = lon_edges-180
        
                # Change extreme centre values to be correct
                if lat_centre[0] == -90:
                    lat_centre[0] = (lat_edges[0] + lat_edges[1])/2. 

                if lat_centre[-1] == 90:
                    lat_centre[-1] = (lat_edges[-2] + lat_edges[-1])/2.

                #correct for lon edges having error for GEOS CCM, model centres are constant but lon edges has error near 0. Change 0 edge to -1.25
                if model == 'GEOS_CCM':
                    ind = np.where(lon_edges==0)[0]
                    lon_edges[ind] = -1.25
                
                #lat centres are not correct at edges so correct them
                if (model == 'EMAC') or  (model == 'MIROCCHEM') or (model == 'OSLOCTM2') or (model == 'CMAM'):
                    lat_centre[0] = lat_edges[0] + ((np.abs(lat_edges[0]) - np.abs(lat_edges[1]))/2.)
                    lat_centre[-1] = lat_edges[-2] + ((lat_edges[-1] - lat_edges[-2])/2.)

                all_o3 = data
                all_time = time
            
            elif count > 0:
                all_o3 = np.vstack((all_o3,data))
                all_time = np.append(all_time,time) 

        #----------------------------------------------------------
        #make fullgrid of -99999's
        if y == 0:
            all_full_grid = np.empty((len(species_list),n_hours,len(lat_centre),len(lon_centre)))
        full_grid = np.empty((n_hours,len(lat_centre),len(lon_centre)))
        full_grid[:,:,:] = -99999
    

        lat_n = 0
        lon_n = 0

        grid_size = len(lon_centre)*len(lat_centre)

        if hour == 'start':
            syn_grid_time = np.arange(0,n_days,1./24)
        elif hour == 'mid':
            syn_grid_time = np.arange((1./24.)/2.,n_days,1./24)

        test_grid_time = np.round(syn_grid_time,decimals=4)
        syn_grid_time = np.round(syn_grid_time,decimals=5)
        all_time = np.round(all_time,decimals=4)

        #--------------------------------------------------------
        #if model == GFDLAM3, and vers is rcp26 or rcp85 or rcp85em2000 then remove a day from time as it includes a leap year in 2100.
        if model == 'GFDLAM3':
            if (vers == 'rcp26') or (vers == 'rcp85') or (vers == 'rcp85em2000'):
                all_time = all_time - 1. 
    
    
        #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
        val_test =  all_time <= test_grid_time[-1]
        all_time = all_time[val_test]
        try:
            all_o3 = all_o3[val_test,0,:,:]
        except:
            all_o3 = all_o3[val_test,:,:]

        print 'grid time=', test_grid_time
        print 'all time=', all_time
        indices = np.searchsorted(test_grid_time, all_time, side='left')
    
        for box in range(grid_size):
            #create max possible grid
            full_data = np.empty(n_hours)
            full_data[:] = -99999
    
            if (lat_n == 0) & (lon_n == 0):
                print full_data
            
            full_data[indices] = all_o3[:,lat_n,lon_n]
            if (lat_n == 0) & (lon_n == 0):
                print full_data

            full_grid[:,lat_n,lon_n] = full_data

            if (lat_n == 0) & (lon_n == 0):
                print full_data
                print full_grid[:,lat_n,lon_n]

            if lon_n == (len(lon_centre)-1):
                lat_n+=1
                lon_n = 0
            else:
                lon_n+=1

        #---------------------------------------------------------
        #change lon grids from 0 - 360 to -180:180
        lon_mid = len(lon_centre)/2
        cut1 = full_grid[:,:,:lon_mid]
        cut2 = full_grid[:,:,lon_mid:]
        full_grid = np.concatenate((cut2,cut1),axis=2)

        #--------------------------------------------------------
        #convert time to date and time components
        if hour == 'start':
            m = 0
        if hour == 'mid':
            m = 30
        start = datetime.datetime(start_year,1,1,0,m)
        end = datetime.datetime(end_year+1,1,1,0,m)
        datetime_array = pd.date_range(start,end,freq='H')[:-1].tolist()

        l_lat_c = len(lat_centre)/2
        l_lon_c = len(lon_centre)/2

        lat_size = lat_centre[l_lat_c-0] - lat_centre[l_lat_c-1]
        lon_size = lon_centre[l_lon_c-0] - lon_centre[l_lon_c-1]

        lat_size = np.around(lat_size,decimals=4)  
        lon_size = np.around(lon_size,decimals=4)  

        if (lat_size % 1) == 0:
            lat_size = int(lat_size)

        if (lon_size % 1) == 0:
            lon_size = int(lon_size)

        print '\nlat size = ', lat_size
        print 'lon size = ', lon_size, '\n' 
        grid_size = '%sx%s'%(lat_size,lon_size)
        grid_size = str(grid_size)

        #convert all negative numbers from raw data to -99999, don't know why they are there,ask.
        invs = np.where(full_grid <= 0)
        full_grid[invs] = -99999
        
        valid = np.where(full_grid >= 0)
        full_grid[valid] = full_grid[valid]*1e9

        all_full_grid[y,:,:,:] = full_grid

    #save out to netcdf
    #-------------------------------------------
    #setup netcdf output
    root_grp = Dataset('%s/%s_SURFACE_%s_%s_*_*_*_H_%s.nc'%(model,model,start_year,end_year+1,tag), 'w')
    root_grp.description = 'Gridded Hourly Surface species in ppb - Program written by Dene Bowdalo'
    
    root_grp.createDimension('time', len(datetime_array))
    root_grp.createDimension('lat_centre', len(lat_centre))
    root_grp.createDimension('lon_centre', len(lon_centre))
    root_grp.createDimension('lat_edges', len(lat_edges))
    root_grp.createDimension('lon_edges', len(lon_edges))
    root_grp.createDimension('grid_s', 1)
    
    t = root_grp.createVariable('time', 'f8', ('time',))
    lat_c = root_grp.createVariable('lat_centre', 'f4', ('lat_centre',))
    lon_c = root_grp.createVariable('lon_centre', 'f4', ('lon_centre',))
    lat_e = root_grp.createVariable('lat_edges', 'f4', ('lat_edges',))
    lon_e = root_grp.createVariable('lon_edges', 'f4', ('lon_edges',))
    gs = root_grp.createVariable('grid_size', str, ('grid_s',))
    
    #set units
    t.units = 'hours since 0001-01-01 00:00:00'
    t.calendar = 'gregorian'

    #process times to save out (gregorian)
    times = date2num(datetime_array, units=t.units, calendar=t.calendar)
    
    t[:] = times
    lat_c[:] = lat_centre
    lon_c[:] = lon_centre
    lat_e[:] = lat_edges
    lon_e[:] = lon_edges
    gs[0] = grid_size

    for y in range(len(species_list)):
        species_data = all_full_grid[y,:,:,:]
        
        test = np.all(species_data==-99999)
        
        if test == False:
            root_grp.createDimension(species_list[y].lower(), len(datetime_array))
            sp = root_grp.createVariable(species_list[y].lower(), 'f4', (species_list[y].lower(),'lat_centre','lon_centre'))
            sp[:] = species_data
    
    del all_full_grid
    del full_grid
    del species_data

    root_grp.close()
