# -*- coding: utf-8 -*-
"""Routines for processing AirBase data.
Dene Bowdalo
"""

import os, sys,glob
import numpy as np
import datetime
from itertools import izip
import re   #regular expressions
import time #handles time and date calculations
import calendar
import pytz
import csv
from numbers import Number
try:
    from isounidecode import unidecode
except:
    print "No isounidecode; Processing GAW, AirBase may have problems"
#from Helper.Helper import unique
import pandas as pd
import logging
import modules
from netCDF4 import Dataset
from scipy import stats
import multiprocessing
from netCDF4 import num2date, date2num
from tzwhere import tzwhere
import pytz
tz_root = tzwhere.tzwhere(shapely=True,forceTZ=True)

process_group = 'AirBase'

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
output_res = 'H'

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

start = datetime.datetime(year = start_year, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = end_year+1, month = 1, day = 1, hour = 0, minute = 0)
end_pandas = datetime.datetime(year = end_year, month = 12, day = 31, hour = 23, minute = 0)
 
"""
EU AirBase Ozone Data
"""
###################################
# EU AirBase filename handling routines
###################################
def ingest_AIRBASE(species,code,species_mw):
    path = '/work/home/db876/observations/surface/O3/AIRBASE/'
    files = sorted(glob.glob(path+'*_v8/???????%s?????hour.*'%code))
    daily_files = sorted(glob.glob(path+'*_v8/???????%s?????day.*'%code))
    monthly_files = sorted(glob.glob(path+'*_v8/???????%s?????month.*'%code))
    site_files = sorted(glob.glob(path+'*_v8/*_stations.csv'))
    meas_conf_files = sorted(glob.glob(path+'*_v8/*_measurement_configurations.csv'))
    xml_files = sorted(glob.glob(path+'AirBase_v8_xmldata/*.xml'))
    
    #process to get all refs for hourly, daily and monthly resolutions. Preference for hourly then daily and monthly, 
    #with refs unable to overlap between resolutions. However, can be mulitple files for ref in same resolution, so
    #file process ahead carefully avoids removing these.
    
    refs = [f.split('/')[-1][:7] for f in files]
    group_codes = [f.split('/')[-1][12:17] for f in files]
    monthly_refs = [f.split('/')[-1][:7] for f in monthly_files]
    monthly_group_codes = [f.split('/')[-1][12:17] for f in monthly_files]

    data_resolutions = ['hr']*len(refs)
    
    #get files and refs based on output_res
    daily_refs = [f.split('/')[-1][:7] for f in daily_files]
    daily_group_codes = [f.split('/')[-1][12:17] for f in daily_files]
    
    #limit to sites with hourly date files for, if required add in extra resolutions
    if (output_res == 'HD') or (output_res == 'HDM'):
        refs1 = refs[:]
        refs2 = refs[:]
        #get daily files
        for i in range(len(daily_refs)):
            r_day = daily_refs[i]
            if r_day not in refs1:
                refs2.append(r_day)
                data_resolutions.append('da')
                files.append(daily_files[i])
                group_codes.append(daily_group_codes[i])

        if output_res == 'HDM':
            #get monthly files
            refs3 = refs2[:]
            for i in range(len(monthly_refs)):
                r_month = monthly_refs[i]
                if r_month not in refs2:
                    refs3.append(r_month) 
                    data_resolutions.append('mo') 
                    files.append(monthly_files[i])
                    group_codes.append(monthly_group_codes[i])

            refs = refs3[:]
        
    #refs = daily_refs
    #files = daily_files
    #data_resolutions = ['da']*len(refs)
    #group_codes = daily_group_codes
    
    tz_meta = parse_AIRBASE_xml(xml_files)
    metadata = parse_AIRBASE_header(site_files,meas_conf_files,tz_meta,species)
    O3all = [parse_AIRBASE_data(fi, metadata,species,species_mw,res) for fi,res in zip(files,data_resolutions)]
    O3all = filter(None, O3all)
    #O3all = [item for sublist in O3all for item in sublist]
    return O3all,group_codes,data_resolutions

###################################
# EU AirBase file parsing
###################################

def parse_AIRBASE_xml(xml_files):
    import untangle
    timezone_meta = {}
    # loop over countries
    for fi in xml_files:
        print fi
        cmeta = untangle.parse(fi)
        # loop over stations in each file
        for i in xrange(len(cmeta.airbase.country.station)):
            # get station name
            station = 'AB'+unidecode(cmeta.airbase.country.station[i].station_european_code.cdata)
            # get time zone information.  If there is none, assume UTC
            try:
                TZ = unidecode(cmeta.airbase.country.station[i].network_info.network_time_reference_basis.cdata)
            except:
                TZ = 'UTC'
            timezone_meta[station] = TZ
            tz_meta = pd.Series(timezone_meta)
            
    return tz_meta


def parse_AIRBASE_header(site_files, meas_conf_files,tz_meta,species):
    columns = ['SHORT NAME', 'STATION NAME','TIME ZONE', 'LONGITUDE', 'LATITUDE', 'ALTITUDE', 'COUNTRY/TERRITORY', 'COUNTY', 'LAND_USE', 'MEASUREMENT UNIT', 'PARAMETER', 'CALIBRATION STANDARD ID', 'MEASUREMENT METHOD', 'TIME INTERVAL', 'STATION CATEGORY','IS EMEP','MEASUREMENT AUTOMATIC']
    metadata_all = pd.DataFrame(columns=columns)
    
    for site_file,meas_conf_file in zip(site_files,meas_conf_files):
        # read metadata by country
        site_meta = pd.read_csv(site_file, sep='\t')
        site_meta.set_index(site_meta.station_european_code,inplace=True)
        meas_meta = pd.read_csv(meas_conf_file, sep='\t')
        meas_meta.set_index(meas_meta.station_european_code,inplace=True)
        metadata = pd.DataFrame(columns=columns, index=site_meta.station_european_code)

        metadata.ix[:,'SHORT NAME'] = site_meta.station_european_code
        metadata.ix[:,'COUNTRY/TERRITORY'] = site_meta.country_name
        metadata.ix[:,'STATION NAME'] = site_meta.station_name.apply(unidecode)
        metadata.ix[:,'DATE_START'] = site_meta.station_start_date
        metadata.ix[:,'STATION CATEGORY'] = site_meta.station_type_of_area+', '+site_meta.station_ozone_classification
        metadata.ix[:,'LAND_USE'] = site_meta.type_of_station
        metadata.ix[:,'LONGITUDE'] = site_meta.station_longitude_deg
        metadata.ix[:,'LATITUDE'] = site_meta.station_latitude_deg
        metadata.ix[:,'ALTITUDE'] = site_meta.station_altitude
        metadata.ix[:,'COUNTY'] = site_meta.station_city
        metadata.ix[:,'IS EMEP'] = site_meta.EMEP_station

        if species != 'ISOP':
            ozone_meas_meta = meas_meta.ix[meas_meta.component_caption==species,:]
        else:
            ozone_meas_meta = meas_meta.ix[meas_meta.component_caption=='CH2=CH-C(CH3)=CH2',:]
            
        # reduce metadata to just the sites that have data
        has_ozone = metadata.ix[:,'SHORT NAME'].apply(lambda x: x in ozone_meas_meta.station_european_code)
        metadata = metadata.ix[has_ozone,:]
        if all(ozone_meas_meta.measurement_unit == '\xc2\xb5g/m3'):
            metadata.ix[ozone_meas_meta.station_european_code,'MEASUREMENT UNIT'] = 'ug/m3'
        else:
            metadata.ix[ozone_meas_meta.station_european_code,'MEASUREMENT UNIT'] = ozone_meas_meta.measurement_unit

        metadata.ix[ozone_meas_meta.station_european_code,'CALIBRATION STANDARD ID'] = ozone_meas_meta.calibration_method
        metadata.ix[ozone_meas_meta.station_european_code,'PARAMETER'] = ozone_meas_meta.component_name
        metadata.ix[ozone_meas_meta.station_european_code,'MEASUREMENT METHOD'] = ozone_meas_meta.measurement_technique_principle +ozone_meas_meta.measurement_equipment
        metadata.ix[ozone_meas_meta.station_european_code,'TIME INTERVAL'] = ozone_meas_meta.integration_time_frequency.astype(str) + ozone_meas_meta.integration_time_unit
        metadata.ix[ozone_meas_meta.station_european_code,'MEASUREMENT AUTOMATIC'] = ozone_meas_meta.measurement_automatic
        
        metadata.ix[:,'SHORT NAME'] = metadata.ix[:,'SHORT NAME'].apply(lambda x: 'AB'+x)
        metadata.index = metadata.ix[:,'SHORT NAME']
        
        metadata.ix[:,'TIME ZONE'] = tz_meta[metadata.index]
        
        metadata_all = metadata_all.append(metadata)   
    return metadata_all
    
def parse_AIRBASE_data(datafile,metadata,species,species_mw,res):
    print os.path.split(datafile)[-1]
    SHORTNAME = 'AB'+os.path.split(datafile)[-1][0:7]
    RawData = pd.read_csv(datafile,sep='\t', header = None)
    RawData.set_index(0,drop=True, inplace=True)
    
    #if len(RawData.columns)<48:
    #    return
    
    # use pivot table to get into series?.  Need to deal with bad data points.
    if (res == 'hr') or (res == 'da'):
        Data = RawData.loc[:,np.arange(1,RawData.columns[-1]+1,2)]
        Flags = RawData.loc[:,np.arange(2,RawData.columns[-1]+1,2)]
    
    elif res == 'mo':
        Data = RawData.loc[:,[1,2]]
        Flags = RawData.loc[:,3]
    
    if res == 'hr':
        Data.columns =np.arange(24)
        Flags.columns = np.arange(24)
    elif res == 'da':
        Data.columns =np.arange(31)
        Flags.columns = np.arange(31)
    elif res == 'mo':
        Data.columns =np.arange(2)
        Flags.columns = np.arange(1)
    
    if (res == 'hr') or (res == 'da'):
        Data = Data.stack()
        Flags = Flags.stack()
    
    ntimes = len(Data)
    tzd = {'UTC' : 0, 'CET' : 1, 'EET' : 2}
    if res == 'hr':
        Data.index = (pd.DatetimeIndex(pd.to_datetime(Data.index.get_level_values(0)+' '+Data.index.get_level_values(1).map('{:02d}:00:00'.format))) - pd.Timedelta(tzd[metadata.ix[SHORTNAME,'TIME ZONE']], unit='h')).tz_localize('UTC')
    elif (res == 'da'):
        date = Data.index.get_level_values(0)
        days = Data.index.get_level_values(1)
        dt_array = []
        for i in range(len(date)):  
            dt_array.append(pd.to_datetime(date[i]) + pd.Timedelta(int(days[i]),unit='d'))
        Data.index = pd.DatetimeIndex(dt_array).tz_localize('UTC')
    elif (res == 'mo'):
        Data.index = (pd.DatetimeIndex(pd.to_datetime(Data.index.get_level_values(0)))).tz_localize('UTC')
    
    Flags.index = Data.index
    
    Data1 = pd.concat([Data,Flags], axis=1)
    
    if (res == 'hr') or (res == 'da'):
        Data1.columns = [species, 'Flags']
    elif res == 'mo':
        Data1.columns = ['EndDate', species, 'Flags']

    #remove any lines where data is equal to zero, as these are lines which are artefact of having to have 31 days for each month.
    inv_test = Data1.ix[:,species] != 0   
    Data1 = Data1.ix[inv_test,:]

    """
    Measurement values 1...24 and quality flags 1...24 in columns 2-49 
    correspond to the hours of that date.
    Value 1 is understood to represent the measurements during the first 
    hour of the day (i.e.: 00:00-00:59), and so on.

    Quality flags in raw data
    Flag values indicate the quality of the preceding measurement value.
    A quality flag value > 0 indicates valid measurement data. 
    A quality flag <= 0 indicates invalid or missing data.
    from http://ftp.eea.europa.eu/www/AirBase_v8/airbase_v8_products.pdf"""


    # screen bad data
    baddata = Data1['Flags']<=0
    Data1.ix[baddata,species] = np.nan
    Data1.ix[Data1[species]< 0.0,species] = np.nan

    """ convert unit
    Note: official documentation for EU air quality states T,P of 293 and 1013 hPa
    http://eur-lex.europa.eu/legal-content/EN/TXT/PDF/?uri=CELEX:32002L0003&from=EN
    """
    if metadata.ix[SHORTNAME,'MEASUREMENT UNIT'] == 'ug/m3':
        Data1.ix[:,species] = Data1.ix[:,species]*convert_ugm3_ppbv(MW=species_mw)
    
    if metadata.ix[SHORTNAME,'MEASUREMENT UNIT'] == 'mg/m3':
        Data1.ix[:,species] = (Data1.ix[:,species]*convert_mgm3_ppbv(MW=species_mw))*1e3
    
    if (res == 'hr') or (res == 'da'):
        Data = Data1[species]
    elif res == 'mo':
        Data = Data1[['EndDate',species]]
    
    try:
        return metadata.ix[SHORTNAME,:].to_dict(), Data.ix[start:end_pandas] 
    except:
        return metadata.ix[SHORTNAME,:].to_dict(), Data
    
def convert_ugm3_ppbv(MW=species_mw):
    """Conversion factor from ug/m3 to ppbv for gases based on a given temperature, pressure, and molecular weight.
    """
    return (1.0/MW)*8.3144*(293.)/(1013./10.0)
    
def convert_mgm3_ppbv(MW=species_mw):
    """Conversion factor from mg/m3 to ppbv for gases based on a given temperature, pressure, and molecular weight.
    """
    return ((1.0/MW)*8.3144*(293.)/(1013./10.0))
       

#find n_hours and n_days between start and end date
d0 = datetime.date(start_year, 1, 1)
d1 = datetime.date(end_year+1, 1, 1)
delta = d1 - d0
n_months = ((end_year+1)-start_year)*12
n_days = delta.days
n_hours = n_days*24
n_points = n_hours
 
data,group_codes,resolutions = ingest_AIRBASE(species,airbase_code, species_mw)
n_refs = len(data)


#MAKE SURE THERE ARE NO DUPLICATE REFS
#IF SO GO THROUGH PROCEDURE TO DETERMINE WHY THEY ARE DUPLICATE, IF CANNOT UNDERSTAND WHY THEN APPEND DATA TOGETHER
valid_refs =[]
for i in range(n_refs):
    valid_refs.append(data[i][0]['SHORT NAME'])
    try:
        if data[i][1].index.min() < min_date:
            min_date = data[i][1].index.min()
            
    except:
        min_date = data[i][1].index.min()

del_list = []
unique_refs = [valid_refs[i] for i in sorted(np.unique(valid_refs, return_index=True)[1])]
dup_list = []

for i in range(len(unique_refs)):   
    indices = [j for j, x in enumerate(valid_refs) if x == unique_refs[i]]
    dup_list.append(indices)
    
for i in range(len(dup_list)):
    current_inds = dup_list[i]
    current_ref = unique_refs[i]
    if len(current_inds) > 1:
        print '%s is duplicated'%(current_ref)
        iteration = 0
        for x in range(len(current_inds)-1):
            current_lat = data[current_inds[x]][0]['LATITUDE']
            current_lon = data[current_inds[x]][0]['LONGITUDE']
            current_alt = data[current_inds[x]][0]['ALTITUDE']
            current_mm = data[current_inds[x]][0]['MEASUREMENT METHOD']

            next_lat = data[current_inds[x+1]][0]['LATITUDE']
            next_lon = data[current_inds[x+1]][0]['LONGITUDE']
            next_alt = data[current_inds[x+1]][0]['ALTITUDE']
            next_mm = data[current_inds[x+1]][0]['MEASUREMENT METHOD']
                
            #test if lat/lon/altitude are different 
            if (np.around(current_lat,decimals=2) != np.around(next_lat,decimals=2)) or (np.around(current_lon,decimals=2) != np.around(next_lon,decimals=2)) or (np.abs(np.diff([current_alt,next_alt]))[0] >= 50):
                print 'Next Site Ref is Duplicate, but different location. Changing next ref name to stop duplication.'
                valid_refs[current_inds[x+1]] = current_ref+str(iteration) 
                print '**************',current_ref,valid_refs[current_inds[x+1]]
                data[current_inds[x+1]][0]['SHORT NAME'] = valid_refs[current_inds[x+1]]
                iteration+=1
            else: 
                #test if appended ref has different measurement method to previous
                #if so then make new ref
                if current_mm != next_mm:
                    print 'Next site ref is Duplicate, but with different measurement methods. Changing next ref name to stop duplication.'
                    valid_refs[current_inds[x+1]] = current_ref+str(iteration) 
                    print '**************',current_ref,valid_refs[current_inds[x+1]]
                    data[current_inds[x+1]][0]['SHORT NAME'] = valid_refs[current_inds[x+1]]
                    iteration+=1
                else:
                    print 'Next site ref is Duplicate, with same location and measurment method - Appending' 
                    data[current_inds[x]][1].append(data[current_inds[x+1]][1])
                    del_list.append(current_inds[x+1])

#if there are indices in del list then remove slices from valid_refs, data, group_codes and resolutions arrays
n_deleted = 0
if len(del_list) > 0:
    for i in del_list:
        valid_refs.pop(i-n_deleted)
        data.pop(i-n_deleted) 
        resolutions.pop(i-n_deleted)
        n_deleted+=1
    
grid_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]
grid_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='H')][:-1]

#make output dates and times
dt = pd.date_range(start,end,freq='H')[:-1].tolist()
output_res_times = date2num(dt, units='hours since 0001-01-01 00:00:00', calendar='gregorian')
obs_time_pd = pd.date_range(start = start,end = end, freq = 'H')[:-1]
    
syn_grid_time = np.arange(0,len(grid_dates)/24,1./24)
syn_grid_time = np.round(syn_grid_time,decimals=5)

#setup netcdf output
root_grp = Dataset('AIRBASE_SURFACE_%s_%s_%s_%s.nc'%(fname_species,start_year,end_year+1,output_res), 'w')
root_grp.description = 'Hourly Surface %s at sites - Program written by Dene Bowdalo'%(fname_species)

# dimensions
root_grp.createDimension('flex', n_points)
root_grp.createDimension('flex4', None)

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

    data_valid = True

    site_data = data[c]
    site_meta = site_data[0]    
    file_res = resolutions[c]    
    
    #get data and metadata
    try:
        lat = np.float32(site_meta['LATITUDE'])
    except:
        lat = 'na'
    try:
        lon = np.float32(site_meta['LONGITUDE'])
    except:
        lon = 'na'
    try:
        alt = np.float32(site_meta['ALTITUDE'])
    except:
        alt = 'na'
    land_use_class = site_meta['LAND_USE']
    if pd.isnull(land_use_class) == True:
        land_use_class = 'na'
    station_class = site_meta['STATION CATEGORY']
    if pd.isnull(station_class) == True:
        station_class = 'na'
    raw_class_name = land_use_class+' '+station_class
    mm = site_meta['MEASUREMENT METHOD']
    if pd.isnull(mm) == True:
        mm = ''
    country = site_meta['COUNTRY/TERRITORY']
    if pd.isnull(country) == True:
        country = 'na'
    site_name = site_meta['STATION NAME']
    if pd.isnull(site_name) == True:
        site_name = 'na'
    continuous_check = site_meta['MEASUREMENT AUTOMATIC']
    if pd.isnull(continuous_check) == True:
        continuous_check = 'na'
    unit = site_meta['MEASUREMENT UNIT']
    #integration_time = site_meta['TIME INTERVAL']
    tz = site_meta['TIME ZONE']
    contact = 'david.simoens@eea.europa.eu'
    #convert timezone from str to int
    tzd = {'UTC' : 0, 'CET' : 1, 'EET' : 2}
    data_tz = tzd[tz]
    all_tz = [data_tz]
        
    if (file_res == 'hr') or (file_res == 'da'):
        var = np.array(site_data[1].values.tolist())
    elif file_res == 'mo':
        all_var = np.array(site_data[1].values.tolist())
        var = np.array(all_var[:,1]).astype('float64')
        end_times = all_var[:,0]
        end_date_con = [d[:4]+d[5:7]+d[8:10] for d in end_times]
        end_time_con = [d[11:13]+d[14:] for d in end_times]
        
    times = site_data[1].index
    date_con = [d.strftime('%Y%m%d') for d in times]
    time_con = [d.strftime('%H%M') for d in times]
    
    #get ref
    site_ref = valid_refs[c]
    site_group = group_codes[c]
    
    print 'ref == %s, %s'%(site_ref,c) 
    print 'res = ',file_res
    
    #add var to total obs count
    n_all += len(var)
    n_after_nometa += len(var)
    
    #if file resolution is daily or monthly then replicate times after point, to fill hourly data array.
    count=0
    if file_res == 'hr':
        n_dup_array = np.zeros(len(var))
    
    elif file_res == 'da':
        n_dup_array = []
        file_hours = len(date_con)
        for i in range(file_hours):
            current_hh = int(time_con[count][:2])
            current_mm = int(time_con[count][2:])
            s = datetime.datetime(year = start_year, month = 1, day = 1, hour = current_hh, minute = current_mm)
            e = datetime.datetime(year = start_year, month = 1, day = 2, hour = current_hh, minute = current_mm)
            day_hours = [d.strftime('%H%M') for d in pd.date_range(s,e,freq='H')][1:-1]

            date_con = np.insert(date_con,count+1,[date_con[count]]*23)
            time_con = np.insert(time_con,count+1,day_hours)
            var = np.insert(var,count+1,[var[count]]*23)
            
            #append to n duplicated array
            n_dup_array=np.append(n_dup_array,0)
            n_dup_array=np.append(n_dup_array,[1]*23)
       
            count +=24

    elif file_res == 'mo':
        n_dup_array = []
        file_hours = len(date_con)
    
        for i in range(file_hours):
            current_year = int(date_con[count][:4])
            current_month = int(date_con[count][4:6])
            current_day = int(date_con[count][6:])
            current_hour = int(time_con[count][:2])
            current_min = int(time_con[count][2:])
        
            next_year = int(end_date_con[i][:4])
            next_month = int(end_date_con[i][4:6])
            next_day = int(end_date_con[i][6:])
            next_hour = int(end_time_con[i][:2])
            next_min = int(end_time_con[i][2:])
        
            s = datetime.datetime(year = current_year, month = current_month, day = current_day, hour = current_hour, minute = current_min)
            e = datetime.datetime(year = next_year, month = next_month, day = next_day, hour = next_hour, minute = next_min)
        
            day_date = [d.strftime('%Y%m%d') for d in pd.date_range(s,e,freq='H')][1:-1]
            day_hour = [d.strftime('%H%M') for d in pd.date_range(s,e,freq='H')][1:-1]
            date_con = np.insert(date_con,count+1,day_date)
            time_con = np.insert(time_con,count+1,day_hour)
            var = np.insert(var,count+1,[var[count]]*len(day_date))
            
            #append to n duplicated array
            n_dup_array=np.append(n_dup_array,0)
            n_dup_array=np.append(n_dup_array,[1]*len(day_date))
            
            count += (len(day_date)+1)

    date_con = np.array(date_con).astype(int)
    time_con = np.array(time_con).astype(int)
    
    #remove data < 1970 and >= 2015
    test_inds = (date_con >= 19700101) & (date_con < 20150101)
    date_con = date_con[test_inds]
    time_con = time_con[test_inds]
    var = var[test_inds]
    n_dup_array = n_dup_array[test_inds]
    
    #convert nans to -99999's
    nan_inds = np.isnan(var)
    var[nan_inds] = -99999
    
    if continuous_check == 'yes':
        st_big = ['continuous']*len(var)
    else:
        st_big = ['filter']*len(var)    
    mm_big = [mm]*len(var)
    
    #get obs valid
    test = var >= 0
    valid_hours_dup = np.sum(n_dup_array[test])
    n_obs_valid = int(len(var[test]) - valid_hours_dup)
    n_after_flagsandlod += n_obs_valid
    
    #create max possible grid
    full_data = np.empty(len(grid_dates))
    full_data_after_flagsandlod = np.empty(n_hours)
    big_n_dup_array = np.zeros(n_hours)
    full_data[:] = -99999
    full_data_after_flagsandlod[:] = -99999
    
    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
    converted_time = modules.date_process(date_con,time_con,start_year)
    converted_time = np.round(converted_time,decimals=5)
    raw_indices = np.searchsorted(syn_grid_time, converted_time, side='left')
    var = np.array(var)
    full_data_after_flagsandlod[raw_indices] = var
    raw_st = np.copy(st_big)
    raw_mm = np.copy(mm_big)
    
    #test and remove duplicate and overlap points
    converted_time,var,mm_big,st_big,n_dup_array = modules.remove_duplicate_points(site_ref,converted_time,var,mm_big,st_big,n_dup_array,output_res)
    test = var >= 0
    valid_hours_dup = np.sum(n_dup_array[test])
    n_obs_valid = int(len(var[test]) - valid_hours_dup)
    n_after_duplicate += n_obs_valid
    
    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
    indices = np.searchsorted(syn_grid_time, converted_time, side='left')
    full_data[indices] = var
    big_n_dup_array[indices] = n_dup_array
            
    key_meta = [lat,lon,alt]
    
    #convert file res to standard format
    if file_res == 'hr':
        file_res = 'H'
    elif file_res == 'da':
        file_res = 'D'
    elif file_res == 'mo':
        file_res = 'M'
    
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
e5d= root_grp.createVariable('exit_badmeasurementmethod_pg',str, ('flex3',))

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
