from netCDF4 import Dataset
import numpy as np
import glob
import datetime
import pandas as pd

start_year = 2005
end_year = 2009

year_range = range(start_year,end_year+1)

#historical or historical_sd
vers = 'historical_sd'

if vers == 'historical':
    #DAILY 2005-2010 MODELS HINDCAST
    models = ['HADGEM3ES','IPSL','MIROC3','SOCOL3','UKCA']
    tag = '*'
    
if vers == 'historical_sd':
    #DAILY 2005-2010 MODELS HINDCAST SPECIFIED DYNAMICS
    models = ['IPSL','MIROC3','MRI']
    tag = 'SD'

#---------------------------------------------------
leap_years = range(1960,2200,4)

#create valid time range
d0 = datetime.date(start_year, 1, 1)
d1 = datetime.date(end_year, 12, 31)
delta = d1 - d0

start = datetime.datetime(start_year,1,1,0,0)
end = datetime.datetime(end_year+1,1,1,0,0)
dt_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='M')]
dt_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='M')]

#datetime and offsets from start ref
start_ref = 1960
d_ref = datetime.date(start_ref, 1, 1)
deltastart = d0 - d_ref
deltaend = d1 - d_ref
valid_start = deltastart.days  
valid_end = deltaend.days

#---------------------------------------------------  
for i in range(len(models)):
    model = models[i]
    print model

    files = glob.glob('%s/%s/D/*'%(model,vers))
    files = sorted(files)
    
    time = []
    ccmi_year = []
    ccmi_month = []
    ccmi_day = []
    ccmi_hour = []
    ccmi_minute = []
    
    for count in range(len(files)):
        root_grp = Dataset(files[count])
        time = np.append(time,root_grp.variables['time'][:])
        if model != 'SOCOL3':
            try:
                data = np.concatenate((data,root_grp.variables['vmro3'][:,0,:,:]),axis=0)
            except:
                data = root_grp.variables['vmro3'][:,0,:,:]
        else:
            try:
                data = np.concatenate((data,root_grp.variables['sfvmro3'][:,:,:]),axis=0)
            except:
                data = root_grp.variables['sfvmro3'][:,:,:]
        try:
            ccmi_year=np.append(ccmi_year,root_grp.variables['ccmi_year'][:])
            ccmi_month=np.append(ccmi_month,root_grp.variables['ccmi_month'][:])
            ccmi_day=np.append(ccmi_day,root_grp.variables['ccmi_day'][:])
            have_ccmi_time = True
        except:
            have_ccmi_time = False
            
                
        if count == 0:
            if model != 'SOCOL3':
                 try:
                     print root_grp.variables['lev'] 
                     print root_grp.variables['lev'][:]
                     print root_grp.variables['lev'][0]
                 except:
                     print root_grp.variables['plev']            
                     print root_grp.variables['plev'][:]
                     print root_grp.variables['plev'][0]

            if (model == 'HADGEM3ES') or (model == 'IPSL') or (model == 'MRI') or (model == 'SOCOL3') or (model == 'UKCA'):
                file_startdate = files[count][-20:-12]
                file_enddate = files[count][-11:-3]
                file_startyear = file_startdate[:4]
                file_startmonth = file_startdate[4:6]
                file_startday = file_startdate[6:]
                file_starthour = 0
                file_startmin = 0
        
            if (model == 'MIROC3'):
                file_startdate = files[count][-16:-10]
                file_enddate = files[count][-9:-3]
                file_startyear = file_startdate[:4]
                file_startmonth = file_startdate[4:]
                file_startday = 0
                file_starthour = 0
                file_startmin = 0
        
            lat_centre = np.sort(root_grp.variables['lat'][:])
            lon_centre = np.sort(root_grp.variables['lon'][:])

            if model == 'MIROC3':
                #make lat edges
                lat_edges = [-90]
                for i in range(len(lat_centre)-1):
                    offset = (lat_centre[i+1] - lat_centre[i])/2.
                    lat_edges=np.append(lat_edges,lat_centre[i]+offset)
                lat_edges=np.append(lat_edges,90)
                
                #make lon edges
                offset = ((lon_centre[0]+360) - lon_centre[-1])/2.
                lon_edges = [0-offset]
                for i in range(len(lon_centre)-1):
                    offset = (lon_centre[i+1] - lon_centre[i])/2.
                    lon_edges=np.append(lon_edges,lon_centre[i]+offset)
                lon_edges=np.append(lon_edges,360-offset)
                
                #correct extreme lat centres
                offset = (lat_edges[-1] - lat_edges[-2])/2.
                lat_centre[0] = lat_edges[0] + offset
                lat_centre[-1] = lat_edges[-1] - offset
            
            else:
                lat_bnds = root_grp.variables['lat_bnds'][:]
                lon_bnds = root_grp.variables['lon_bnds'][:]

                lat_edges = np.sort(lat_bnds[:,0])
                lon_edges = np.sort(lon_bnds[:,0])

                #add 360 to lon band edges
                lon_edges = np.append(lon_edges,(360-np.abs(lon_edges[0])))
                #ind = np.where(lon_edges==0)[0]

            #print 'orig lat centre', lat_centre
            #print 'orig lat edges', lat_edges   
            #print 'orig lon centre', lon_centre
            #print 'orig lon edges', lon_edges
            
            #if model equal to SOCOL3 then extreme lat edges not equal to -90/90 but -87.86379884/87.86379884. (as it has no pole?)
            #change lat edges to be equal to -90/90 for uniformity, and also change end lat centres
            if model == 'SOCOL3':
                lat_edges[0] = -90
                lat_edges = np.append(lat_edges,90)
                lat_centre[0] = lat_edges[0]+((lat_edges[1] - lat_edges[0])/2.)
                lat_centre[-1] = lat_edges[-1]-((lat_edges[-1] - lat_edges[-2])/2.)
                
            #if model equal to UKCA correct lon edges at extremes as not correct
            if model == 'UKCA':
                diff = lon_edges[2] - lon_edges[1]
                lon_edges[0] = lon_edges[1] - diff
                lon_edges[-1] = lon_edges[-2]+diff
            
            #correct extreme edges to be correct
            if (lat_edges[0] != -90):
                lat_edges = np.insert(lat_edges,0,-90)
            
            if (lat_edges[-1] != 90):
                lat_edges = np.append(lat_edges,90)
            
            # Change extreme lat centre values to be correct
            if (lat_centre[0] == -90):
                gap = np.abs(lat_edges[0]) - np.abs(lat_edges[1])
                new_centre = (gap/2.) + lat_edges[0]
                lat_centre[0] = new_centre

            if (lat_centre[-1] == 90):
                gap = lat_edges[-1] - lat_edges[-2]
                new_centre = lat_edges[-2] + (gap/2.) 
                lat_centre[-1] = new_centre
            
            #if model equal to MRI, then correct extreme lat centres to be between 86.405 and 90, as default is not
            if model == 'MRI':
                diff = (lat_edges[-1] - lat_edges[-2])/2.
                lat_centre[0] = lat_edges[0] + offset
                lat_centre[-1] = lat_edges[-1] - offset
                
            if lon_centre[0] > -100:
                lon_centre = lon_centre-180
                lon_edges = lon_edges-180
                
            #print 'corrected lat centre',lat_centre
            #print 'corrected lat edges', lat_edges
            #print 'corrected lon centre',lon_centre
            #print 'corrected lon edges', lon_edges
            
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

            print 'lat size = ', lat_size
            print 'lon size = ', lon_size 
            grid_size = '%sx%s'%(lat_size,lon_size)
            grid_size = str(grid_size)
                
            #make fullgrid of -99999's
            full_grid = np.empty((len(dt_dates),len(lat_centre),len(lon_centre)))
            full_grid[:,:,:] = -99999
    
    #if have ccmi time fill in empty array (deals with leap years)
    if have_ccmi_time == True:
        timeres = ccmi_day[1] - ccmi_day[0]
        ccmi_dt = []
        for i in range(len(ccmi_year)):
            ccmi_dt.append(datetime.datetime(int(ccmi_year[i]),int(ccmi_month[i]),int(ccmi_day[i]),0,0))
        
    else:
        timeres = time[1] - time[0]
        #print timeres
        if timeres == 1:
        #currently only IPSL and MRI, both leap years so processing is straightforward
            year = int(file_startyear)
            month = int(file_startmonth)
            day = int(file_startday)
            hour = int(file_starthour)
            min = int(file_startmin)     
            base_datetime = datetime.datetime(year,month,day,hour,min)
            
            time = time - time[0]
            
            ccmi_dt = []
            for i in range(len(time)):
                ccmi_dt.append(base_datetime+datetime.timedelta(days=time[i]))
            
        #if timeres equal to 10, make datetime array from times, with days on 1,11,21,1 days etc...
        if timeres == 10:
        
            year = int(file_startyear)
            month = int(file_startmonth)
            day = int(file_startday)
            hour = int(file_starthour)
            min = int(file_startmin)
            
            ccmi_dt = []
            for i in range(len(time)):
                ccmi_dt.append(datetime.datetime(year,month,day,hour,min))
                if day == 21:
                    day = 1
                    month = month +1
                    if month == 13:
                        month = 1
                        year = year+1
                else:
                    day= day+10

    start_ind = np.searchsorted(ccmi_dt,valid_dt_range_list[0])
    end_ind = np.searchsorted(ccmi_dt,valid_dt_range_list[-1],side='right')
    ccmi_dt = ccmi_dt[start_ind:end_ind]
    data = data[start_ind:end_ind,:,:]
    
    #get indices for ccmi time array in full time array 
    inds = np.searchsorted(valid_dt_range_list,ccmi_dt)
    full_grid[inds,:,:] = data[:,:,:]
    
    #delete data for future iterations
    del data

#---------------------------------------------------------
    #change lon grids from 0 - 360 to -180:180
    lon_mid = len(lon_centre)/2
    cut1 = full_grid[:,:,:lon_mid]
    cut2 = full_grid[:,:,lon_mid:]
    full_grid = np.concatenate((cut2,cut1),axis=2)


    #convert all negative numbers from raw data to -99999, don't know why they are there,ask.
    invs = np.where(full_grid < 0)
    full_grid[invs] = -99999
    
    #convert very big numbers from raw data to -99999,
    invs = np.where(full_grid >= 1e10)
    full_grid[invs] = -99999

    #save out to netcdf
    #-------------------------------------------

    #setup netcdf output
    root_grp = Dataset('%s/%s_SURFACE_O3_%s_%s_*_*_*_D_%s.nc'%(model,model,start_year,end_year+1,tag), 'w')
    root_grp.description = 'Gridded Hourly Surface O3 in ppb - Program written by Dene Bowdalo'

    # dimensions
    #set variables
    root_grp.createDimension('o3', len(dt_dates))
    root_grp.createDimension('dates', len(dt_dates))
    root_grp.createDimension('times', len(dt_dates))
    root_grp.createDimension('lat_centre', len(lat_centre))
    root_grp.createDimension('lon_centre', len(lon_centre))
    root_grp.createDimension('lat_edges', len(lat_edges))
    root_grp.createDimension('lon_edges', len(lon_edges))
    root_grp.createDimension('grid_s', 1)

    #set variables
    species = root_grp.createVariable('o3', 'f8', ('o3','lat_centre','lon_centre'))
    d = root_grp.createVariable('date', 'i8', ('dates',))
    t = root_grp.createVariable('time', 'i8', ('times',))
    lat_c = root_grp.createVariable('lat_centre', 'f8', ('lat_centre',))
    lon_c = root_grp.createVariable('lon_centre', 'f8', ('lon_centre',))
    lat_e = root_grp.createVariable('lat_edges', 'f8', ('lat_edges',))
    lon_e = root_grp.createVariable('lon_edges', 'f8', ('lon_edges',))
    gs = root_grp.createVariable('grid_size', str, ('grid_s',))

    species[:] = full_grid
    d[:] = dt_dates
    t[:] = dt_times
    lat_c[:] = lat_centre
    lon_c[:] = lon_centre
    lat_e[:] = lat_edges
    lon_e[:] = lon_edges
    gs[0] = grid_size
