import numpy as np
from netCDF4 import Dataset
import datetime
import pandas as pd
import modules
import os
from netCDF4 import num2date, date2num
import matplotlib.pyplot as plt
#import seaborn as sns
from collections import Counter
from scipy.stats import mode

species = os.getcwd().split('/')[-2]

#All netcdf files for different datasets
start_year_full = 1970
end_year_full = 2015

#SET FILE RES, CAN BE HOURLY(H), DAILY(D), MONTHLY(M) OR ANNUAL(Y) 
file_res = 'H'

#OUTPUT RESOLUTION CAN BE HOURLY (H), HOURLY & DAILY (HD), HOURLY, DAILY AND MONTHLY (HDM)
output_res = 'HDM'

#OUTPUT SET CAN BE STANDARD(S), PERIODIC (P), TREND (T), NO RESTRAINTS (N)
#if contains N (No-restraint): (No raw class check, no anthrome class check, no altitude class check, no manual extreme data check)
#if contains P (Periodic): (remove uneven timezone data, test for data gaps > 1 year)
#If contains S (Standard): (standard set with all checks except periodic)
#If contains T (Trend): (same as standard but does not remove partial years)
#Can combine NP to have no restraint and periodic data
#*Data representativeness is checked for all data
#*Median absolute deviation used to remove outliers for all data

output_set = 'S'

#IF OUTPUT TYPE IS PERIODIC THEN SET STEP AND GAP REQUIRED and START AND END YEAR SET ,ELSE LEAVE STEP AS 1 AND GAP AS 45 AND YEARS AS 1970 AND 2015
gap = 45
step = 1
start_year_set = 1970
end_year_set = 2015

start_years = np.arange(start_year_set,(end_year_set+1)-gap,step)
end_years  = np.arange(start_year_set+gap,end_year_set+1,step)

print start_years
print end_years


if species == 'O3':
    AIRBASE = 'AIRBASE_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    AQS = 'AQS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    CAPMON = 'CAPMON_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    CASTNET = 'CASTNET_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    EANET = 'EANET_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    EMEP = 'EMEP_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    GAW = 'GAW_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    NAPS = 'NAPS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    SEARCH = 'SEARCH_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    all_files = [AIRBASE,AQS,CAPMON,CASTNET,EANET,EMEP,GAW,NAPS,SEARCH]
    #all_files = [CAPMON,CASTNET,EANET,EMEP,GAW,NAPS,SEARCH]
    #all_files = [AIRBASE]
elif species == 'CO':
    AIRBASE = 'AIRBASE_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    AQS = 'AQS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    CASTNET = 'CASTNET_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    EMEP = 'EMEP_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    GAW = 'GAW_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    NAPS = 'NAPS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    SEARCH = 'SEARCH_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    all_files = [AIRBASE,AQS,CASTNET,EMEP,GAW,NAPS,SEARCH]
    #all_files = [GAW,NAPS]
elif species == 'NO':
    AIRBASE = 'AIRBASE_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    AQS = 'AQS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    CASTNET = 'CASTNET_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    EANET = 'EANET_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    EMEP = 'EMEP_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    GAW = 'GAW_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    NAPS = 'NAPS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    SEARCH = 'SEARCH_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    all_files = [AIRBASE,AQS,CASTNET,EANET,EMEP,GAW,NAPS,SEARCH]
    #all_files = [GAW]
elif species == 'NO2-MOLYBDENUM':
    AIRBASE = 'AIRBASE_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    AQS = 'AQS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)  
    EANET = 'EANET_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    EMEP = 'EMEP_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    GAW = 'GAW_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    NAPS = 'NAPS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    SEARCH = 'SEARCH_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    all_files = [AIRBASE,AQS,EANET,EMEP,GAW,NAPS,SEARCH]
elif species == 'NO2-PHOTOLYTIC':
    AIRBASE = 'AIRBASE_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    AQS = 'AQS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)  
    EANET = 'EANET_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    EMEP = 'EMEP_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    GAW = 'GAW_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    NAPS = 'NAPS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    SEARCH = 'SEARCH_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    all_files = [AIRBASE,AQS,EANET,EMEP,GAW,NAPS,SEARCH]
elif species == 'ISOP':
    AIRBASE = 'AIRBASE_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    AQS = 'AQS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    EMEP = 'EMEP_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    GAW = 'GAW_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res) 
    NAPS = 'NAPS_SURFACE_%s_%s_%s_%s.nc'%(species,start_year_full,end_year_full,output_res)
    all_files = [AIRBASE,AQS,EMEP,GAW,NAPS]

#SET EXIT ERROR COUNTS AT 0
inv_nometa_count = 0
inv_anyvaliddata_count = 0
inv_nokeymeta_count = 0
inv_resolution_count = 0
inv_badmeasurementmethod_count = 0
inv_duplicatesites_count = 0
inv_rawclass_count = 0
inv_anthromeclass_count = 0
inv_altitude_count = 0
inv_night_count = 0
inv_representativeness_count = 0
inv_extreme_count = 0
inv_partialyear_count = 0
inv_timezone_count = 0
inv_bigdatagap_count = 0
n_final = 0

#SET N OBS COUNTS
n_obs_all = 0
n_obs_after_nometa = 0
n_obs_after_flagsandlod = 0
n_obs_after_duplicatepoints = 0
n_obs_after_anyvaliddata = 0
n_obs_after_nokeymeta = 0
n_obs_after_resolution = 0
n_obs_after_badmeasurementmethod = 0
n_obs_after_duplicatesites = 0
n_obs_after_rawclass = 0
n_obs_after_anthromeclass = 0
n_obs_after_altitude = 0
n_obs_after_night = 0
n_obs_after_representativeness = 0
n_obs_after_extreme = 0
n_obs_after_partialyear = 0
n_obs_after_timezone = 0
n_obs_after_bigdatagap = 0

#set stage output data
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
exit_duplicatesites_refs = []
exit_duplicatesites_lats = []
exit_duplicatesites_lons = []
exit_duplicatesites_pg = []
exit_rawclass_refs = []
exit_rawclass_lats = []
exit_rawclass_lons = []
exit_rawclass_pg = []
exit_anthromeclass_refs = []
exit_anthromeclass_lats = []
exit_anthromeclass_lons = []
exit_anthromeclass_pg = []
exit_altitude_refs = []
exit_altitude_lats = []
exit_altitude_lons = []
exit_altitude_pg = []
exit_night_refs = []
exit_night_lats = []
exit_night_lons = []
exit_night_pg = []
exit_representativeness_refs = []
exit_representativeness_lats = []
exit_representativeness_lons = []
exit_representativeness_pg = []
exit_extreme_refs = []
exit_extreme_lats = []
exit_extreme_lons = []
exit_extreme_pg = []
exit_partialyear_refs = []
exit_partialyear_lats = []
exit_partialyear_lons = []
exit_partialyear_pg = []
exit_timezone_refs = []
exit_timezone_lats = []
exit_timezone_lons = []
exit_timezone_pg = []
exit_bigdatagap_refs = []
exit_bigdatagap_lats = []
exit_bigdatagap_lons = []
exit_bigdatagap_pg = []


#load all data into big lists
all_refs = []
all_data = []
all_dup = []
all_instru = []
all_lats = []
all_lons = []
all_alts = []
all_sitenames = []
all_contacts = []
all_processgroups = []
all_countries = []
all_rawclass = []
all_runit = []
all_punit = []
all_nativeres = []
all_datatz = []
all_localtz = []
all_flaskflag = []

for f in all_files:
    print f
    root_grp_read = Dataset(f)
    
    all_time = root_grp_read.variables['time'][:]
    
    #iterate primary counts
    n_obs_all+=root_grp_read.variables['n_obs_all'][0]
    n_obs_after_nometa+=root_grp_read.variables['n_obs_after_nometa'][0]
    n_obs_after_flagsandlod+=root_grp_read.variables['n_obs_after_flagsandlod'][0]
    n_obs_after_duplicatepoints+=root_grp_read.variables['n_obs_after_duplicate'][0]
    n_obs_after_anyvaliddata+=root_grp_read.variables['n_obs_after_anyvaliddata'][0]
    n_obs_after_nokeymeta+=root_grp_read.variables['n_obs_after_nokeymeta'][0] 
    n_obs_after_resolution+=root_grp_read.variables['n_obs_after_resolution'][0] 
    n_obs_after_badmeasurementmethod+=root_grp_read.variables['n_obs_after_badmeasurementmethod'][0] 
    inv_nometa_count+=root_grp_read.variables['invalid_nometa_count'][0]
    inv_anyvaliddata_count+=root_grp_read.variables['invalid_anyvaliddata_count'][0] 
    inv_nokeymeta_count+=root_grp_read.variables['invalid_nokeymeta_count'][0] 
    inv_resolution_count+=root_grp_read.variables['invalid_resolution_count'][0] 
    inv_badmeasurementmethod_count+=root_grp_read.variables['invalid_badmeasurementmethod_count'][0] 
    exit_nometa_refs = np.append(exit_nometa_refs,root_grp_read.variables['exit_nometa_refs'][:])
    exit_nometa_lats = np.append(exit_nometa_lats,root_grp_read.variables['exit_nometa_lats'][:])
    exit_nometa_lons = np.append(exit_nometa_lons,root_grp_read.variables['exit_nometa_lons'][:])
    exit_nometa_pg = np.append(exit_nometa_pg,root_grp_read.variables['exit_nometa_pg'][:])
    exit_anyvaliddata_refs = np.append(exit_anyvaliddata_refs,root_grp_read.variables['exit_anyvaliddata_refs'][:])
    exit_anyvaliddata_lats = np.append(exit_anyvaliddata_lats,root_grp_read.variables['exit_anyvaliddata_lats'][:])
    exit_anyvaliddata_lons = np.append(exit_anyvaliddata_lons,root_grp_read.variables['exit_anyvaliddata_lons'][:])
    exit_anyvaliddata_pg = np.append(exit_anyvaliddata_pg,root_grp_read.variables['exit_anyvaliddata_pg'][:])
    exit_nokeymeta_refs = np.append(exit_nokeymeta_refs,root_grp_read.variables['exit_nokeymeta_refs'][:])
    exit_nokeymeta_lats = np.append(exit_nokeymeta_lats,root_grp_read.variables['exit_nokeymeta_lats'][:])
    exit_nokeymeta_lons = np.append(exit_nokeymeta_lons,root_grp_read.variables['exit_nokeymeta_lons'][:])
    exit_nokeymeta_pg = np.append(exit_nokeymeta_pg,root_grp_read.variables['exit_nokeymeta_pg'][:])
    exit_resolution_refs = np.append(exit_resolution_refs,root_grp_read.variables['exit_resolution_refs'][:])
    exit_resolution_lats = np.append(exit_resolution_lats,root_grp_read.variables['exit_resolution_lats'][:])
    exit_resolution_lons = np.append(exit_resolution_lons,root_grp_read.variables['exit_resolution_lons'][:])
    exit_resolution_pg = np.append(exit_resolution_pg,root_grp_read.variables['exit_resolution_pg'][:])
    exit_badmeasurementmethod_refs = np.append(exit_badmeasurementmethod_refs,root_grp_read.variables['exit_badmeasurementmethod_refs'][:])
    exit_badmeasurementmethod_lats = np.append(exit_badmeasurementmethod_lats,root_grp_read.variables['exit_badmeasurementmethod_lats'][:])
    exit_badmeasurementmethod_lons = np.append(exit_badmeasurementmethod_lons,root_grp_read.variables['exit_badmeasurementmethod_lons'][:])
    exit_badmeasurementmethod_pg = np.append(exit_badmeasurementmethod_pg,root_grp_read.variables['exit_badmeasurementmethod_pg'][:])
    
    valid_refs_dict = root_grp_read.groups
    valid_refs = []
    for i in valid_refs_dict.keys():
        valid_refs = np.append(valid_refs,i)
    #valid_refs = valid_refs_dict.keys()
    print valid_refs
    #valid_refs = np.array(['aqs270750005','aqs270750005','aqs290390001','aqs292110001','aqs300650004','aqs300830001','aqs380130001','aqs440030002','aqs460330132','aqs550410007'])
    #valid_refs = np.array(['canaps030701','canaps092401'])

    for site_ref in valid_refs:
        site_group = root_grp_read.groups[site_ref]
        
        all_refs.append(site_ref)
        all_data.append(site_group.variables[species.split('-')[0].lower()][:])   
        all_instru.append(site_group.variables['measurement method'][:])
        all_dup.append(site_group.variables['n_duplicate_array'][:])
        all_lats.append(site_group.latitude)
        all_lons.append(site_group.longitude)
        all_alts.append(site_group.altitude)
        all_sitenames.append(site_group.site_name)
        all_contacts.append(site_group.data_contact)
        all_processgroups.append(site_group.process_group)
        all_countries.append(site_group.country)
        all_rawclass.append(site_group.raw_site_class)
        all_runit.append(site_group.raw_units)
        all_punit.append(site_group.processed_units)
        all_nativeres.append(site_group.native_resolution)
        all_datatz.append(site_group.data_timezone)
        all_localtz.append(site_group.local_timezone)
        all_flaskflag.append(site_group.flask_flag)
        
#convert to numpy arrays
all_refs = np.array(all_refs)
all_data = np.array(all_data)
all_instru = np.array(all_instru)
all_dup = np.array(all_dup)
all_lats = np.array(all_lats).astype('float32')
all_lons = np.array(all_lons).astype('float32')
all_alts = np.array(all_alts).astype('float32')
all_sitenames = np.array(all_sitenames)
all_contacts = np.array(all_contacts)
all_processgroups = np.array(all_processgroups)
all_countries = np.array(all_countries)
all_rawclass = np.array(all_rawclass)
all_runit = np.array(all_runit)
all_punit = np.array(all_punit)
all_nativeres = np.array(all_nativeres)
all_datatz = np.array(all_datatz)
all_localtz = np.array(all_localtz)
all_flaskflag= np.array(all_flaskflag)

print 'N refs before Duplicate removal = ', len(all_refs)

#get unique site instruments
unique_instru = []
for l in all_instru:
    l = l[(l != -1) & (l != 0) & (l != 99)]
    unique_instru.append(list(set(l)))
    
unique_flat = [item for sublist in unique_instru for item in sublist]
unique_set = list(set(unique_flat))
print 'Unique Instruments = ', unique_set

#remove duplicate sites
del_inds = modules.remove_duplicate_sites(all_refs,all_lats,all_lons,all_alts,all_data,all_nativeres,unique_instru)
inv_duplicatesites_count+=len(del_inds)
exit_duplicatesites_refs = np.append(exit_duplicatesites_refs,all_refs[del_inds])
exit_duplicatesites_lats = np.append(exit_duplicatesites_lats,all_lats[del_inds].astype('str'))
exit_duplicatesites_lons = np.append(exit_duplicatesites_lons,all_lons[del_inds].astype('str'))
exit_duplicatesites_pg = np.append(exit_duplicatesites_pg,all_processgroups[del_inds])

all_count = np.arange(len(all_refs))
keep_inds = [int(i) for i in all_count if int(i) not in del_inds] 

print 'N refs after Duplicate removal = ', len(keep_inds)

all_refs = all_refs[keep_inds]
all_data = all_data[keep_inds]
all_instru = all_instru[keep_inds]
all_dup = all_dup[keep_inds]
all_lats = all_lats[keep_inds]
all_lons = all_lons[keep_inds]
all_alts = all_alts[keep_inds]
all_sitenames = all_sitenames[keep_inds]
all_contacts = all_contacts[keep_inds]
all_processgroups = all_processgroups[keep_inds]
all_countries = all_countries[keep_inds]
all_rawclass = all_rawclass[keep_inds]
all_runit = all_runit[keep_inds]
all_punit = all_punit[keep_inds]
all_nativeres = all_nativeres[keep_inds]
all_datatz = all_datatz[keep_inds]
all_localtz = all_localtz[keep_inds]
all_flaskflag = all_flaskflag[keep_inds]

#iterate through steps
for start_year,end_year in zip(start_years,end_years):
    print start_year,end_year
    #find n_hours and n_days between start and end date
    d0 = datetime.date(start_year, 1, 1)
    d1 = datetime.date(end_year, 1, 1)
    delta = d1 - d0
    n_years = end_year-start_year
    n_months = ((end_year)-start_year)*12
    n_days = delta.days
    n_hours = n_days*24

    #SET TIME OUTPUT
    if file_res == 'H':
        n_points = n_hours 
    if file_res == 'D':
        n_points = n_days
    if file_res == 'M':
        n_points = n_months
    if file_res == 'Y':
        n_points = n_years

    #setup netcdf output
    root_grp = Dataset('GLOBAL_SURFACE_%s_%s_%s_%s_%s.nc'%(species,start_year,end_year,file_res,output_res+output_set), 'w')
    root_grp.description = '%s Surface %s at sites in ppb - Program written by Dene Bowdalo'%(file_res,species)

    # dimensions
    root_grp.createDimension('time', n_points) 
    times = root_grp.createVariable('time', 'f8', ('time',))
    root_grp.createDimension(species.lower(), n_points)

    #create grid time array
    start = datetime.datetime(year = start_year, month = 1, day = 1, hour = 0, minute = 0)
    end = datetime.datetime(year = end_year, month = 1, day = 1, hour = 0, minute = 0)
    grid_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]

    if (file_res == 'H') or (file_res == 'D'):
        dt = pd.date_range(start,end,freq=file_res)[:-1].tolist()
        obs_time_pd = pd.date_range(start = start,end = end, freq = 'H')[:-1]
        grid_times = date2num(dt, units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    elif (file_res == 'M'):
        dt = pd.date_range(start,end,freq='MS')[:-1].tolist()
        obs_time_pd = pd.date_range(start = start,end = end, freq = 'MS')[:-1]
        grid_times = date2num(dt, units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    elif (file_res == 'Y'):
        dt = pd.date_range(start,end,freq='AS')[:-1].tolist()
        obs_time_pd = pd.date_range(start = start,end = end, freq = 'AS')[:-1]
        grid_times = date2num(dt, units='hours since 0001-01-01 00:00:00', calendar='gregorian')

    #save out time
    times[:] = grid_times
    
    #get datetime cut
    datetime_cut = num2date(grid_times,units='hours since 0001-01-01 00:00:00', calendar='gregorian')

    #start write out 
    for x in range(len(all_refs)):
        site_ref = all_refs[x]
        print site_ref

        #GET DATA FOR REF
        vals = all_data[x]
        instru_data = all_instru[x]
        flask_flag = all_flaskflag[x]
        lat = all_lats[x]
        lon = all_lons[x]
        alt = all_alts[x]
        site_name = all_sitenames[x]
        contact = all_contacts[x]
        process_group = all_processgroups[x]
        country = all_countries[x]
        raw_class_name = all_rawclass[x]
        r_unit = all_runit[x]
        p_unit = all_punit[x]
        native_res = all_nativeres[x]
        data_tz = all_datatz[x]
        local_tz = all_localtz[x]
        n_dup_array = all_dup[x]
        flask_flag = all_flaskflag[x]

        #DO SECONDARY CONTROL CHECKS
        test = vals >= 0
        valid_hours_dup = np.sum(n_dup_array[test])
        n_obs_valid = int(len(vals[test]) - valid_hours_dup)
        n_obs_after_duplicatesites+=n_obs_valid
        meta = [lat,lon,alt,raw_class_name,data_tz]
        
        #IF OUTPUT RES IS PERIODIC THEN CUT TO TIME SPAN WANTED
        if (start_year > start_year_full) or (end_year < end_year_full):
            indices = np.searchsorted(all_time,grid_times,side='left')
            vals = vals[indices]
            instru_data = instru_data[indices]
 
        #t = vals < 0
        #vals[t] = np.nan
        #print len(vals), len(vals[~np.isnan(vals)])
        #plt.plot(grid_times,vals,color='blue',marker='o',linestyle='None')
        #plt.show()
 
        vals,data_valid,anthrome_class_name,instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap = modules.secondary_quality_control(site_ref,species,n_dup_array,flask_flag,instru_data,meta,vals,n_obs_valid,output_set,process_group,start_year,end_year,grid_dates,datetime_cut,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap)
        if exit_r == 'rawclass':
            exit_rawclass_refs = np.append(exit_rawclass_refs,site_ref)
            exit_rawclass_lats = np.append(exit_rawclass_lats,str(lat))
            exit_rawclass_lons = np.append(exit_rawclass_lons,str(lon))
            exit_rawclass_pg = np.append(exit_rawclass_pg,process_group)
        elif exit_r == 'anthromeclass':
            exit_anthromeclass_refs = np.append(exit_anthromeclass_refs,site_ref)
            exit_anthromeclass_lats = np.append(exit_anthromeclass_lats,str(lat))
            exit_anthromeclass_lons = np.append(exit_anthromeclass_lons,str(lon))
            exit_anthromeclass_pg = np.append(exit_anthromeclass_pg,process_group)
        elif exit_r == 'altitude':
            exit_altitude_refs = np.append(exit_altitude_refs,site_ref)
            exit_altitude_lats = np.append(exit_altitude_lats,str(lat))
            exit_altitude_lons = np.append(exit_altitude_lons,str(lon))
            exit_altitude_pg = np.append(exit_altitude_pg,process_group)
        elif exit_r == 'night':
            exit_night_refs = np.append(exit_night_refs,site_ref)
            exit_night_lats = np.append(exit_night_lats,str(lat))
            exit_night_lons = np.append(exit_night_lons,str(lon))
            exit_night_pg = np.append(exit_night_pg,process_group)    
        elif exit_r == 'representativeness':
            exit_representativeness_refs = np.append(exit_representativeness_refs,site_ref)
            exit_representativeness_lats = np.append(exit_representativeness_lats,str(lat))
            exit_representativeness_lons = np.append(exit_representativeness_lons,str(lon))
            exit_representativeness_pg = np.append(exit_representativeness_pg,process_group)
        elif exit_r == 'extreme':
            exit_extreme_refs = np.append(exit_extreme_refs,site_ref)
            exit_extreme_lats = np.append(exit_extreme_lats,str(lat))
            exit_extreme_lons = np.append(exit_extreme_lons,str(lon))
            exit_extreme_pg = np.append(exit_extreme_pg,process_group)
        elif exit_r == 'partialyear':
            exit_partialyear_refs = np.append(exit_partialyear_refs,site_ref)
            exit_partialyear_lats = np.append(exit_partialyear_lats,str(lat))
            exit_partialyear_lons = np.append(exit_partialyear_lons,str(lon))
            exit_partialyear_pg = np.append(exit_partialyear_pg,process_group)
        elif exit_r == 'timezone':
            exit_timezone_refs = np.append(exit_timezone_refs,site_ref)
            exit_timezone_lats = np.append(exit_timezone_lats,str(lat))
            exit_timezone_lons = np.append(exit_timezone_lons,str(lon))
            exit_timezone_pg = np.append(exit_timezone_pg,process_group)
        elif exit_r == 'bigdatagap':
            exit_bigdatagap_refs = np.append(exit_bigdatagap_refs,site_ref)
            exit_bigdatagap_lats = np.append(exit_bigdatagap_lats,str(lat))
            exit_bigdatagap_lons = np.append(exit_bigdatagap_lons,str(lon))
            exit_bigdatagap_pg = np.append(exit_bigdatagap_pg,process_group)
    
        if data_valid == True:
            
            #print len(vals), len(vals[vals >= 0])
            #t = vals < 0
            #vals[t] = np.nan
            #plt.plot(vals,color='green')
            #plt.show()

        
            full_data = np.empty(n_hours)
            full_data[:] = -99999
            full_data[:] = vals


            #TAKE AVERAGE IF NEEDED
            if file_res != 'H':
                full_data = modules.take_average(full_data,obs_time_pd,file_res)

            try:
                ref = root_grp.createGroup('%s'%(site_ref))
            except:
                ref = root_grp.createGroup('%s_%s'%(site_ref,process_group))

            #set variables
            big_data = ref.createVariable(species.lower(), 'f8', (species.lower(),))
            instru = ref.createVariable('measurement method', 'int', (species.lower(),))
            big_data[:] = full_data
            instru[:] = instru_data
            
            #set core meta attributes
            ref.latitude = lat
            ref.longitude = lon
            ref.altitude = alt
            ref.site_name = site_name
            ref.country = country
            ref.data_contact = contact
            ref.process_group = process_group
            ref.data_timezone = data_tz
            ref.local_timezone = local_tz
            ref.raw_site_class = raw_class_name
            ref.anthrome_site_class = anthrome_class_name
            #set final class, must be urban or rural
            final_class_name = 'rural'
            #raw check
            if (process_group == 'AirBase') or (process_group == 'EPA AQS') or (process_group == 'CAPMON') or (process_group == 'CASTNET') or (process_group == 'EANET') or (process_group == 'SEARCH'):
                if ('urban' in raw_class_name.lower()) or ('traffic' in raw_class_name.lower()) or ('industrial' in raw_class_name.lower()):
                    final_class_name = 'urban'
            elif (process_group == 'NAPS'):
                if ('i' == raw_class_name.lower()) or ('c' == raw_class_name.lower()):
                  final_class_name = 'urban'  
            #secondary check
            if anthrome_class_name == 'Dense Settlement':
                final_class_name = 'urban' 
            #manual check
            inv_refs = ['abat2ka71','abat2kl17','abat2m226','abat2sv24','abat2vk26','abat30402','abat30406','abat30603','abat30901','abat30992','abat4s184','abat4s401','abat4s406','abat4s412','abat4s413','abat4s414','abat60103','abat60104','abat60106','abat60114','abat60115','abat60119','abat60139','abat60140','abat60141','abat60145','abat60160','abat60161','abat60162','abat60172','abat60197','abat72126','abat82707','abba0034a','abba0035a',
                        'abbg0019a','abbg0056a','abbg0057a','abbg0058a','abbg0074a','abch0023a','abch0036a','abch0040a','abch0041a','abcz0asan','abcz0ljnm','abcz0molj','abcz0molo','abcz0skls','abcz0tovk','abcz0tozr','abcz0uchm','abcz0udcm','abdebb016','abdebb037','abdeby009','abdeby067','abdeby069','abdehe002','abdehh054','abdeni026','abdenw066','abderp042','abdest006','abdeth012','abdeth017','abee0021a','abes1180a','abes1364a', 'abes1400a','abes1538a','abes1851a','abfr0037a','abfr02005','abfr02037','abfr03014','abfr03038','abfr03047','abfr03062','abfr03066','abfr05070','abfr06018','abfr07002','abfr07003','abfr07035','abfr08003','abfr08711','abfr10007','abfr10014','abfr10016','abfr10026','abfr10041','abfr1054a','abfr1056a','abfr12003','abfr1250a','abfr13007','abfr14003','abfr14005','abfr14006','abfr14007','abfr14012','abfr15002','abfr15007','abfr16054','abfr16056','abfr16057','abfr18053','abfr22004','abfr22016','abfr23058','abfr23065','abfr23093','abfr23174','abfr23179','abfr23188','abfr24030','abfr25028','abfr25048','abfr26003','abfr26006','abfr27001','abfr28001','abfr28003','abfr28011','abfr28131','abfr32014','abfr33305','abfr41001','abgb0036r','abgb0040r','abgb0051a','abgb0962a','abgr0033a','abgr0042a','abgr0110r','abie0132a','abit0440a','abit0441a','abit0447a','abit0558a','abit0659a','abit0663a','abit0684a',                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                        'abit0692a','abit0709a','abit0782a','abit0838a','abit0862a','abit0906a','abit0982a','abit1042a','abit1061a','abit1065a','abit1103a','abit1156a','abit1222a','abit1225a','abit1238a','abit1292a','abit1328a','abit1340a','abit1439a','abit1459a','abit1468a','abit1557a','abit1571a','abit1670a','abit1725a','abit1845a','abit1868a','abit1919a',                                                                                                                                             
                        'abit1925a','abit1944a','abit2032a','abit2051a','ablt00001','ablt00012','ablt0021a','abmk0030a','abmk0048a','abnl00107','abnl00442','abno0075a','abpl0011a','abpl0021a','abpl0033a','abpl0039a','abpl0043a','abpl0051a','abpl0126a','abpl0115a','abpl0119a','abpl0122a','abpl0126a','abpl0134a','abpl0136a','abpl0144a','abpl0148a','abpl0162a','abpl0198a','abpl0221a','abpl0256a','abpl0276a','abpl0321a','abpt03101','abpt03103','abro0074a','abro0075a','abro0077a','abro0079a','abro0084a','abro0191a','abro0208a','abrs0030a','abrs0031a','abrs0033a','abrs0034a','abrs0039a','abrs0040a','abrs0043a','abrs0044a','abrs0045a','abrs0046a','abrs0048a','abse0053a','absk0003a','absk0004a','absk0008a','absk0015a','absk0028a','absk0047a','absk0051a','absk0267a','aqs060290011','aqs060651016','aqs180830004','aqs181290001','aqs320030007','aqs400219002','aqs400719003','aqs420110001','canaps051802','canaps060807','canaps061004','canaps065001','canaps090402','canaps090601','canaps101101','canaps101701','canaps103205','canaps104003','canaps104003','canaps104301','bsc','bur','car','log','pld','plv','sof','vrn',
                        'abat60180','abbg0071a','abes1811a','abfr01014','abfr03032','abfr07038','abfr12052','abfr14022','abfr14042','abfr14051','abfr15018','abfr17016','abfr19010','abfr19015','abfr20061','abfr24033','abfr26019','abfr29439','abfr36019','abfr37003','abfr39007','abfr39008','abfr39009','abit0824a','aqs320030023','canaps040302','canaps060809','canaps092901','canaps101501','canaps102401','canaps102701','aqs291830008','ablv0rim1','abse0054a','abse0055a',
                        'abbg0026a','abbg0038a','abdebb010','abdest074','abdeub020','abfi00066','abfi00424','abfi00582','abfr02026','abfr05086','abfr06134','abfr08401','abfr12053','abfr15001','abfr18049','abfr22022','abfr23181','abfr23182','abfr34034','abfr34046','abfr41002','abgb0929a','abgb0995a','abgb0999a','abgb1005a','abgr0120a','abie0140a','abit1307a','abit1663a','abit1729a','abit2054a','abpl0151a','abpl0152a','abpl0187a','abpl0210a','abro0126a','abse0030a','abse0080a','absi0033a','absk0031a','aqs050350005','aqs350451233','canaps060610','canaps061104','canaps061201','canaps061802','canaps062001','canaps090502','canaps102102','canaps102801','canaps107100','canaps119004','canaps010401','abit1311a','ch0031r','canaps060602','canaps060806','snl','abis0015a','canaps030201','canaps094301','abpl0212a'] 
            if site_ref.lower() in inv_refs:
                final_class_name = 'urban'
            ref.final_class = final_class_name
            ref.raw_units = r_unit
            ref.processed_units = p_unit
            ref.native_resolution = native_res  
            
            n_final+=1

root_grp.close()

#SETUP STATS OUTPUT  

if (output_res == 'HDM') & (file_res == 'H'):
    root_grp_stats = Dataset('STATS_%s.nc'%(output_res+output_set), 'w')
    root_grp_stats.createDimension('count', 1)
    root_grp_stats.createDimension('exits', None)

    i1 = root_grp_stats.createVariable('invalid_nometa_count','int', ('count',))
    i2 = root_grp_stats.createVariable('invalid_anyvaliddata_count','int', ('count',))
    i3 = root_grp_stats.createVariable('invalid_nokeymeta_count','int', ('count',))
    i4 = root_grp_stats.createVariable('invalid_resolution_count','int', ('count',))
    i5 = root_grp_stats.createVariable('invalid_badmeasurementmethod_count','int', ('count',))
    i6 = root_grp_stats.createVariable('invalid_duplicatesites_count','int', ('count',))
    i7 = root_grp_stats.createVariable('invalid_rawclass_count','int', ('count',))
    i8 = root_grp_stats.createVariable('invalid_anthromeclass_count','int', ('count',))
    i9 = root_grp_stats.createVariable('invalid_altitude_count','int', ('count',))
    i10 = root_grp_stats.createVariable('invalid_night_count','int', ('count',))
    i11 = root_grp_stats.createVariable('invalid_representativeness_count','int', ('count',))
    i12 = root_grp_stats.createVariable('invalid_extreme_count','int', ('count',))
    i13 = root_grp_stats.createVariable('invalid_partialyear_count','int', ('count',))
    i14 = root_grp_stats.createVariable('invalid_timezone_count','int', ('count',))
    i15 = root_grp_stats.createVariable('invalid_bigdatagap_count','int', ('count',))
    i16 = root_grp_stats.createVariable('n_final_count','int', ('count',))

    i1[:] = inv_nometa_count
    i2[:] = inv_anyvaliddata_count 
    i3[:] = inv_nokeymeta_count 
    i4[:] = inv_resolution_count
    i5[:] = inv_badmeasurementmethod_count 
    i6[:] = inv_duplicatesites_count
    i7[:] = inv_rawclass_count
    i8[:] = inv_anthromeclass_count
    i9[:] = inv_altitude_count
    i10[:] = inv_night_count
    i11[:] = inv_representativeness_count
    i12[:] = inv_extreme_count
    i13[:] = inv_partialyear_count
    i14[:] = inv_timezone_count
    i15[:] = inv_bigdatagap_count
    i16[:] = n_final

    n1 = root_grp_stats.createVariable('n_obs_all','int', ('count',))
    n2 = root_grp_stats.createVariable('n_obs_after_nometa','int', ('count',))
    n3 = root_grp_stats.createVariable('n_obs_after_flagsandlod','int', ('count',))
    n4 = root_grp_stats.createVariable('n_obs_after_duplicatepoints','int', ('count',))
    n5 = root_grp_stats.createVariable('n_obs_after_anyvaliddata','int', ('count',))
    n6 = root_grp_stats.createVariable('n_obs_after_nokeymeta','int', ('count',))
    n7 = root_grp_stats.createVariable('n_obs_after_resolution','int', ('count',))
    n8 = root_grp_stats.createVariable('n_obs_after_badmeasurementmethod','int', ('count',))
    n9 = root_grp_stats.createVariable('n_obs_after_duplicatesites','int', ('count',))
    n10 = root_grp_stats.createVariable('n_obs_after_rawclass','int', ('count',))
    n11 = root_grp_stats.createVariable('n_obs_after_anthromeclass','int', ('count',))
    n12 = root_grp_stats.createVariable('n_obs_after_altitude','int', ('count',))
    n13 = root_grp_stats.createVariable('n_obs_after_night','int', ('count',))
    n14 = root_grp_stats.createVariable('n_obs_after_representativeness','int', ('count',))
    n15 = root_grp_stats.createVariable('n_obs_after_extreme','int', ('count',))
    n16 = root_grp_stats.createVariable('n_obs_after_partialyear','int', ('count',))
    n17 = root_grp_stats.createVariable('n_obs_after_timezone','int', ('count',))
    n18 = root_grp_stats.createVariable('n_obs_after_bigdatagap','int', ('count',))


    n1[:] = n_obs_all
    n2[:] = n_obs_after_nometa
    n3[:] = n_obs_after_flagsandlod
    n4[:] = n_obs_after_duplicatepoints
    n5[:] = n_obs_after_anyvaliddata
    n6[:] = n_obs_after_nokeymeta
    n7[:] = n_obs_after_resolution
    n8[:] = n_obs_after_badmeasurementmethod
    n9[:] = n_obs_after_duplicatesites
    n10[:] = n_obs_after_rawclass
    n11[:] = n_obs_after_anthromeclass
    n12[:] = n_obs_after_altitude
    n13[:] = n_obs_after_night
    n14[:] = n_obs_after_representativeness
    n15[:] = n_obs_after_extreme
    n16[:] = n_obs_after_partialyear
    n17[:] = n_obs_after_timezone
    n18[:] = n_obs_after_bigdatagap


    e1 = root_grp_stats.createVariable('exit_nometa_refs',str, ('exits',))
    e2 = root_grp_stats.createVariable('exit_nometa_lats',str, ('exits',))
    e3 = root_grp_stats.createVariable('exit_nometa_lons',str, ('exits',))
    e4 = root_grp_stats.createVariable('exit_nometa_pg',str, ('exits',))
    e5 = root_grp_stats.createVariable('exit_anyvaliddata_refs',str, ('exits',))
    e6 = root_grp_stats.createVariable('exit_anyvaliddata_lats',str, ('exits',))
    e7 = root_grp_stats.createVariable('exit_anyvaliddata_lons',str, ('exits',))
    e8 = root_grp_stats.createVariable('exit_anyvaliddata_pg',str, ('exits',))
    e9 = root_grp_stats.createVariable('exit_nokeymeta_refs',str, ('exits',))
    e10 = root_grp_stats.createVariable('exit_nokeymeta_lats',str, ('exits',))
    e11 = root_grp_stats.createVariable('exit_nokeymeta_lons',str, ('exits',))
    e12 = root_grp_stats.createVariable('exit_nokeymeta_pg',str, ('exits',))
    e13 = root_grp_stats.createVariable('exit_resolution_refs',str, ('exits',))
    e14 = root_grp_stats.createVariable('exit_resolution_lats',str, ('exits',))
    e15 = root_grp_stats.createVariable('exit_resolution_lons',str, ('exits',))
    e16 = root_grp_stats.createVariable('exit_resolution_pg',str, ('exits',)) 
    e17 = root_grp_stats.createVariable('exit_badmeasurementmethod_refs',str, ('exits',))
    e18 = root_grp_stats.createVariable('exit_badmeasurementmethod_lats',str, ('exits',))
    e19 = root_grp_stats.createVariable('exit_badmeasurementmethod_lons',str, ('exits',))
    e20 = root_grp_stats.createVariable('exit_badmeasurementmethod_pg',str, ('exits',))
    e21 = root_grp_stats.createVariable('exit_duplicatesites_refs',str, ('exits',))
    e22 = root_grp_stats.createVariable('exit_duplicatesites_lats',str, ('exits',))
    e23 = root_grp_stats.createVariable('exit_duplicatesites_lons',str, ('exits',))
    e24 = root_grp_stats.createVariable('exit_duplicatesites_pg',str, ('exits',))
    e25 = root_grp_stats.createVariable('exit_rawclass_refs',str, ('exits',))
    e26 = root_grp_stats.createVariable('exit_rawclass_lats',str, ('exits',))
    e27 = root_grp_stats.createVariable('exit_rawclass_lons',str, ('exits',))
    e28 = root_grp_stats.createVariable('exit_rawclass_pg',str, ('exits',))
    e29 = root_grp_stats.createVariable('exit_anthromeclass_refs',str, ('exits',))
    e30 = root_grp_stats.createVariable('exit_anthromeclass_lats',str, ('exits',))
    e31 = root_grp_stats.createVariable('exit_anthromeclass_lons',str, ('exits',))
    e32 = root_grp_stats.createVariable('exit_anthromeclass_pg',str, ('exits',))
    e33 = root_grp_stats.createVariable('exit_altitude_refs',str, ('exits',))
    e34 = root_grp_stats.createVariable('exit_altitude_lats',str, ('exits',))
    e35 = root_grp_stats.createVariable('exit_altitude_lons',str, ('exits',))
    e36 = root_grp_stats.createVariable('exit_altitude_pg',str, ('exits',))
    e37 = root_grp_stats.createVariable('exit_night_refs',str, ('exits',))
    e38 = root_grp_stats.createVariable('exit_night_lats',str, ('exits',))
    e39 = root_grp_stats.createVariable('exit_night_lons',str, ('exits',))
    e40 = root_grp_stats.createVariable('exit_night_pg',str, ('exits',))
    e41 = root_grp_stats.createVariable('exit_representativeness_refs',str, ('exits',))
    e42 = root_grp_stats.createVariable('exit_representativeness_lats',str, ('exits',))
    e43 = root_grp_stats.createVariable('exit_representativeness_lons',str, ('exits',))
    e44 = root_grp_stats.createVariable('exit_representativeness_pg',str, ('exits',)) 
    e45 = root_grp_stats.createVariable('exit_extreme_refs',str, ('exits',))
    e46 = root_grp_stats.createVariable('exit_extreme_lats',str, ('exits',))
    e47 = root_grp_stats.createVariable('exit_extreme_lons',str, ('exits',))
    e48 = root_grp_stats.createVariable('exit_extreme_pg',str, ('exits',))
    e49 = root_grp_stats.createVariable('exit_partialyear_refs',str, ('exits',))
    e50 = root_grp_stats.createVariable('exit_partialyear_lats',str, ('exits',))
    e51 = root_grp_stats.createVariable('exit_partialyear_lons',str, ('exits',))
    e52 = root_grp_stats.createVariable('exit_partialyear_pg',str, ('exits',))


    e1[:] = np.array(exit_nometa_refs,dtype='object')
    e2[:] = np.array(exit_nometa_lats,dtype='object') 
    e3[:] = np.array(exit_nometa_lons,dtype='object')
    e4[:] = np.array(exit_nometa_pg,dtype='object')
    e5[:] = np.array(exit_anyvaliddata_refs,dtype='object')
    e6[:] = np.array(exit_anyvaliddata_lats,dtype='object')
    e7[:] = np.array(exit_anyvaliddata_lons,dtype='object')
    e8[:] = np.array(exit_anyvaliddata_pg,dtype='object')
    e9[:] = np.array(exit_nokeymeta_refs,dtype='object')
    e10[:] = np.array(exit_nokeymeta_lats,dtype='object')
    e11[:] = np.array(exit_nokeymeta_lons,dtype='object')
    e12[:] = np.array(exit_nokeymeta_pg,dtype='object')
    e13[:] = np.array(exit_resolution_refs,dtype='object')
    e14[:] = np.array(exit_resolution_lats,dtype='object')
    e15[:] = np.array(exit_resolution_lons,dtype='object')
    e16[:] = np.array(exit_resolution_pg,dtype='object')
    e17[:] = np.array(exit_badmeasurementmethod_refs,dtype='object')
    e18[:] = np.array(exit_badmeasurementmethod_lats,dtype='object')
    e19[:] = np.array(exit_badmeasurementmethod_lons,dtype='object')
    e20[:] = np.array(exit_badmeasurementmethod_pg,dtype='object')
    e21[:] = np.array(exit_duplicatesites_refs,dtype='object')
    e22[:] = np.array(exit_duplicatesites_lats,dtype='object')
    e23[:] = np.array(exit_duplicatesites_lons,dtype='object')
    e24[:] = np.array(exit_duplicatesites_pg,dtype='object')
    e25[:] = np.array(exit_rawclass_refs,dtype='object')
    e26[:] = np.array(exit_rawclass_lats,dtype='object')
    e27[:] = np.array(exit_rawclass_lons,dtype='object')
    e28[:] = np.array(exit_rawclass_pg,dtype='object')
    e29[:] = np.array(exit_anthromeclass_refs,dtype='object')
    e30[:] = np.array(exit_anthromeclass_lats,dtype='object')
    e31[:] = np.array(exit_anthromeclass_lons,dtype='object')
    e32[:] = np.array(exit_anthromeclass_pg,dtype='object')
    e33[:] = np.array(exit_altitude_refs,dtype='object')
    e34[:] = np.array(exit_altitude_lats,dtype='object')
    e35[:] = np.array(exit_altitude_lons,dtype='object')
    e36[:] = np.array(exit_altitude_pg,dtype='object')
    e37[:] = np.array(exit_night_refs,dtype='object')
    e38[:] = np.array(exit_night_lats,dtype='object')
    e39[:] = np.array(exit_night_lons,dtype='object')
    e40[:] = np.array(exit_night_pg,dtype='object')
    e41[:] = np.array(exit_representativeness_refs,dtype='object')
    e42[:] = np.array(exit_representativeness_lats,dtype='object')
    e43[:] = np.array(exit_representativeness_lons,dtype='object')
    e44[:] = np.array(exit_representativeness_pg,dtype='object')
    e45[:] = np.array(exit_extreme_refs,dtype='object')
    e46[:] = np.array(exit_extreme_lats,dtype='object')
    e47[:] = np.array(exit_extreme_lons,dtype='object')
    e48[:] = np.array(exit_extreme_pg,dtype='object')
    e49[:] = np.array(exit_partialyear_refs,dtype='object')
    e50[:] = np.array(exit_partialyear_lats,dtype='object')
    e51[:] = np.array(exit_partialyear_lons,dtype='object')
    e52[:] = np.array(exit_partialyear_pg,dtype='object')


    root_grp_stats.close() 
