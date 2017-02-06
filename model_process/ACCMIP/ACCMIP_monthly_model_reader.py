from netCDF4 import Dataset
import numpy as np
import glob
import datetime
import pandas as pd
from netCDF4 import num2date, date2num

#species_list = ['O3','NO','NO2','CO','ISOP']
species_list = ['O3','NO','NO2','CO','ISOP','CH4','OH','PAN','TEMP','PRECIP','EMINOX','EMINH3','EMIVOC']

#MONTHLY HISTORICAL MODELS 2000 - 2015 (with 5 species)
models = ['GFDLAM3']
#models = ['CESMCAM','CMAM','EMAC','GEOSCCM','GFDLAM3','GISSE2R','HADGEM2','LMDZORINCA','MIROCCHEM','MOCAGE','NCARCAM3.5','OSLOCTM2','STOCHADAM3','UMCAM']
#start_year = 2000
#end_year = 2011

#MONTHLY RCP2.6 MODELS 2095-2111
#models = ['CESMCAM','CMAM','GFDLAM3','GISSE2R','HADGEM2','LMDZORINCA','MIROCCHEM','MOCAGE','NCARCAM3.5','OSLOCTM2','STOCHADAM3','UMCAM']

#MONTHLY RCP8.5 MODELS 2095-2111
#models = ['CESMCAM','CMAM','EMAC','GFDLAM3','GISSE2R','HADGEM2','LMDZORINCA','MIROCCHEM','MOCAGE','NCARCAM3.5','OSLOCTM2','STOCHADAM3','UMCAM']
start_year = 2095
end_year = 2110

#choose  version of processing. hist, rcp26 or rcp85 or rcp85em2000.
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
for i in range(len(models)):
    model = models[i]
    print '\n------------------------------------'
    print model
    
    for y in range(len(species_list)):
        print species_list[y]
    
        #create max possible grid
        n_months = ((end_year+1)-start_year)*12
        full_data = np.empty(n_months)
        full_data[:] = -99999

        valid_years = np.arange(start_year,end_year+1)

        files = glob.glob('%s/%s/%s/M/*'%(model,species_list[y],vers))
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

        print valid_files
        #---------------------------------------------------  
        #get data
        first_flag = True
        start = 0
        end = 12
        for count in range(len(year_range)):
            for f in range(len(valid_files)):
                if str(year_range[count]) in valid_files[f]:
                    root_grp = Dataset(valid_files[f])
                    if ('EMI' not in species_list[y]) & (species_list[y] != 'TEMP') & (species_list[y] != 'PRECIP'):
                        data = root_grp.variables['vmr%s'%(species_list[y].lower())][:]
                    else:
                        data = root_grp.variables['%s'%(species_list[y].lower())][:]
                    time = root_grp.variables['time'][:]

                    if first_flag == True:
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

                        #if first lat edge is not -90 change it to be
                        if lat_edges[0] != -90:
                            lat_edges[0] = -90 

                        #correct for lon edges having error for GEOS CCM, model centres are constant but lon edges has error near 0. Change 0 edge to -1.25
                        if model == 'GEOS_CCM':
                            ind = np.where(lon_edges==0)[0]
                            lon_edges[ind] = -1.25
        
                        #lat centres are not correct at edges so correct them
                        if (model == 'EMAC') or  (model == 'MIROCCHEM') or (model == 'OSLOCTM2') or (model == 'CMAM'):
                            lat_centre[0] = lat_edges[0] + ((np.abs(lat_edges[0]) - np.abs(lat_edges[1]))/2.)
                            lat_centre[-1] = lat_edges[-2] + ((lat_edges[-1] - lat_edges[-2])/2.)

                        #make fullgrid of -99999's
                        if y == 0:
                            all_full_grid = np.empty((len(species_list),n_months,len(lat_centre),len(lon_centre)))
                        full_grid = np.empty((n_months,len(lat_centre),len(lon_centre)))
                        full_grid[:,:,:] = -99999
                        
                        first_flag = False
                    
                    if ('EMI' not in species_list[y]) & ('PRECIP' not in species_list[y]):
                        full_grid[start:end,:,:] = data[:,0,:,:]
                    else:
                        full_grid[start:end,:,:] = data[:,:,:]
                    break
             
            start+=12
            end+=12      

        grid_size = len(lon_centre)*len(lat_centre)
        #---------------------------------------------------------
        #change lon grids from 0 - 360 to -180:180
        lon_mid = len(lon_centre)/2
        cut1 = full_grid[:,:,:lon_mid]
        cut2 = full_grid[:,:,lon_mid:]
        full_grid = np.concatenate((cut2,cut1),axis=2)

        #--------------------------------------------------------
        #convert time to date and time components
        start = datetime.datetime(start_year,1,1,0,0)
        end = datetime.datetime(end_year+1,1,1,0,0)
        datetime_array = pd.date_range(start,end,freq='M').tolist()

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

        #convert all negative numbers from raw data to -99999, don't know why they are there,ask.
        invs = np.where(full_grid < 0)
        full_grid[invs] = -99999
        
        if ('EMI' not in species_list[y]) & (species_list[y] != 'TEMP') & (species_list[y] != 'PRECIP'):
            valid = np.where(full_grid >= 0)
            full_grid[valid] = full_grid[valid]*1e9
    
        all_full_grid[y,:,:,:] = full_grid

    #save out to netcdf
    #------------------------------------------
    #setup netcdf output
    root_grp = Dataset('%s/%s_SURFACE_%s_%s_*_*_*_M_%s.nc'%(model,model,start_year,end_year+1,tag), 'w')
    root_grp.description = 'Gridded Monthly Surface species in ppb - Program written by Dene Bowdalo'

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
