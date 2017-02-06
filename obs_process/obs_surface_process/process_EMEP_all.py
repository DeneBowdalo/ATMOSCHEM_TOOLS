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
import re
from dateutil import relativedelta
from collections import Counter
from operator import itemgetter
from calendar import monthrange
from scipy import stats
import os
import multiprocessing
from netCDF4 import num2date, date2num
from tzwhere import tzwhere
import pytz
tz_root = tzwhere.tzwhere(shapely=True,forceTZ=True)

process_group = 'EMEP'

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

n_years = int(end_year-start_year)

#get molecular mass and required concentration resolution of species
data_resolution,mol_mass,aqs_code,airbase_code = modules.get_spec_specs(species)

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

#EMEP COUNTRIES AND TIMEZONE, BY COUNTRY CODE
def EMEP_COUNTRIES(c_code):
    countries = {'AL':'Albania','AM':'Armenia','AT':'Austria','AZ':'Azerbaijan','BY':'Belarus','BE':'Belgium','BA':'Bosnia and Herzegovina'
                ,'BG':'Bulgaria','CA':'Canada','HR':'Croatia','CS':'Serbia and Montenegro','CY':'Cyprus','CZ':'Czech Republic','DK':'Denmark','EE':'Estonia'
                ,'FI':'Finland','FR':'France','DE':'Germany','GE':'Georgia','GL':'Greenland','GR':'Greece','HU':'Hungary','KE':'Kenya','ID':'Indonesia','IS':'Iceland','IE':'Ireland'
                ,'IT':'Italy','KG':'Kyrgyzstan','LI':'Liechtenstein','LV':'Latvia','LT':'Lithuania','LU':'Luxembourg','NL':'Netherlands'
                ,'NO':'Norway','PL':'Poland','PT':'Portugal','MA':'Malta','MD':'Moldova','ME':'Montenegro','MT':'Monaco','RO':'Romania'
                ,'RU':'Russia','RS':'Serbia','SK':'Slovakia','SI':'Slovenia','ES':'Spain','SE':'Sweden'
                ,'CH':'Switzerland','MK':'Macedonia','TM':'Turkmenistan','TR':'Turkey','UA':'Ukraine','US':'United States','UZ':'Uzbekistan','TJ':'Tajikistan'
                ,'GB':'United Kingdom','YU':'Yugoslavia'}
    return countries[c_code]

#read obs files
all_files = glob.glob('/work/home/db876/observations/surface/%s/EMEP/*'%(fname_species))
all_files = modules.natsorted(all_files)

#get all refs
ref_list = []
valid_refs=[]

for i in range(len(all_files)):
    f = all_files[i]
    f = f.replace("/work/home/db876/observations/surface/%s/EMEP/"%(fname_species), "")
    f = f[:7]
    ref_list.append(f)

refs = set(ref_list)
refs = sorted([i for i in refs])
valid_refs = np.array(refs)
print 'all refs len = ', len(refs)

#limit to sites with hourly date files
if output_res == 'H':
    del_i = []
    for i in range(len(valid_refs)):
        s_files = glob.glob('/work/home/db876/observations/surface/%s/EMEP/%s*'%(fname_species,valid_refs[i]))
        flag = False
        for f in s_files:
            if '.1h.' in f:
                flag = True       
        if flag == False:
            del_i.append(i)
        else:
            pass
    valid_refs = np.delete(valid_refs,del_i)

#limit to sites with hourly and daily date files
if output_res == 'HD':
    del_i = []
    for i in range(len(valid_refs)):
        s_files = glob.glob('/work/home/db876/observations/surface/%s/EMEP/%s*'%(fname_species,valid_refs[i]))
        flag = False
        for f in s_files:
            if ('.1h.' in f) or ('.1d.' in f):
                flag = True       
        if flag == False:
            del_i.append(i)
        else:
            pass
    valid_refs = np.delete(valid_refs,del_i)

year_array = np.arange(start_year,end_year+1)

#find n_hours and n_days between start and end date
d0 = datetime.date(start_year, 1, 1)
d1 = datetime.date(end_year+1, 1, 1)
delta = d1 - d0
n_months = ((end_year+1)-start_year)*12
n_days = delta.days
n_hours = n_days*24
n_points = n_hours

#setup netcdf output
root_grp = Dataset('EMEP_SURFACE_%s_%s_%s_%s.nc'%(fname_species,start_year,end_year+1,output_res), 'w')
root_grp.description = 'Hourly Surface %s at EMEP sites- Program written by Dene Bowdalo'%(fname_species)

# dimensions
root_grp.createDimension('flex', n_points)
root_grp.createDimension('flex4', None)

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

#read files site at a time
#for ref_i in range(len(valid_refs)):
    site_ref = valid_refs[c]

    all_latitudes = []
    all_longitudes = []
    all_altitudes = []
    all_unit = []
    all_site_name = []
    all_country = []
    all_contact = []
    mm_big = []
    meta_valid_list = []

    data_valid = True

    print 'Current Ref is = ', site_ref,c
    #find if sites have full valid range from start year and finishing in end year
    s_files = glob.glob('/work/home/db876/observations/surface/%s/EMEP/%s*'%(fname_species,site_ref))
    year_files = [file.replace("/work/home/db876/observations/surface/%s/EMEP/"%(fname_species), "") for file in s_files]
    cut_year_files = [file[8:12] for file in year_files]
    site_files = []
    for y in year_array:
        for i in range(len(s_files)):
            if str(y) in cut_year_files[i]:
                site_files.append(s_files[i])
                  
    site_files = modules.natsorted(site_files)
    
    #test for duplicate file years, if duplicates break processing
    file_years = []
    for file in site_files:
        last_file_split = file.split('/')[-1]
        file_years=np.append(file_years,last_file_split[8:12])
    for y in year_array:
        test = file_years == str(y)
        if len(file_years[test]) > 1:
            print 'Site has duplicate files for %s. Breaking processing'%(y)
            1+'a'

    if site_files == []:
        print 'No valid files for site\n'
        return
    
    #remove daily/monthly files if necessary
    if output_res == 'H':
        del_i = []
        for i in range(len(site_files)):
            if '.1d.' in site_files[i]:
                del_i.append(i)
            elif '.1mo.' in site_files[i]:
                del_i.append(i)
        site_files=np.delete(site_files,del_i)
    elif output_res == 'HD':
        del_i = []
        for i in range(len(site_files)):
            if '.1mo.' in site_files[i]:
                del_i.append(i)
        site_files=np.delete(site_files,del_i)
    
    for y in year_array:
        bad_meta = False
        got_year = False
        for file in site_files:
            last_file_split = file.split('/')[-1]
            if str(y) in last_file_split[8:12]:
                got_year = True
                break
        if got_year == False:
            #fill in data for missing year
            timedelta_diff = datetime.date(y+1, 1, 1) - datetime.date(y, 1, 1)
            ndays_missing = timedelta_diff.days       
            continue
    
        count = 0
        with open(file, 'rb') as f:
            reader = csv.reader(f,delimiter=' ')
            print file
            for row in reader:
                try:
                    row = filter(lambda a: a != '', row)
                except:
                    pass
                try:
                    row = filter(lambda a: a != ',', row)
                except:
                    pass
                                
                #get start date of file
                if row[0] == 'Startdate:':
                    data = row[1]
                    s_yyyy = data[:4]
                    s_mm = data[4:6]
                    s_dd = data[6:8]
                    s_hh = data[8:10]
                    s_min = data[10:12]
                    start_datetime = datetime.datetime(int(s_yyyy),1,1,0,0)
                
                #get unit
                if row[0] == 'Unit:':
                    try:
                        if len(row) == 3:
                            unit_part1 = row[1]
                            unit_part2 = row[2]
                            unit = unit_part1+'_'+unit_part2
                        
                        elif len(row) == 2:
                            unit = row[1] 
                        all_unit.append(unit)
                    except:
                        bad_meta = True
        
                #get resolution
                if row[0] == 'Resolution':
                    if row[1] == 'code:':
                        file_res = row[2]
                        print 'Resolution = %s'%file_res
                
                #get latitude
                if row[0] == 'Station':
                    if row[1] == 'latitude:':
                        latitude = row[2]
                        all_latitudes.append(latitude)
            
                #get longitude
                if row[0] == 'Station':
                    if row[1] == 'longitude:':
                        longitude = row[2]
                        all_longitudes.append(longitude)
                    
                #get altitude
                if row[0] == 'Station':
                    if row[1] == 'altitude:':
                        altitude = row[2][:-1]
                        all_altitudes.append(altitude)
                        
                #get site name
                if row[0] == 'Station':
                    if row[1] == 'name:':
                        site_name = row[2]
                        all_site_name.append(site_name)
            
                #get period
                if row[0] == 'Period':
                    period_code = row[2]
                
                #get stats method
                if row[0] == 'Statistics:':
                    try:
                        st = row[1] + row[2]
                        if st != 'arithmeticmean':
                            print 'Not Arithmetic Mean!'
                            print row[1]
                            print 1+'a'  
                    except:
                        print 'Not Arithmetic Mean!'
                        print row[1]
                        print 1+'a'
            
                #get instrument method and name
                if row[0] == 'Instrument':
                    if row[1] == 'type:':
                        mm_list = row[2:]
                        if len(mm_list) > 1:
                            site_mm = ''
                            for x in range(len(mm_list)):
                                site_mm = site_mm+mm_list[x]+' '
                            site_mm = site_mm.strip()
                        else:
                            site_mm = mm_list[0]
                
                    if row[1] == 'name:':
                        mn_list = row[2:]
                        if len(mn_list) > 1:
                            site_mn = ''
                            for x in range(len(mn_list)):
                                site_mn = site_mn+mn_list[x]+' '
                            site_mn = site_mn.strip()
                        else:
                            site_mn = mn_list[0]
                
                #get method ref
                if row[0] == 'Method':
                    if row[1] == 'ref:':
                        try:
                            mf_list = row[2:]
                            if len(mf_list) > 1:
                                site_mf = ''
                                for x in range(len(mf_list)):
                                    site_mf = site_mf+mf_list[x]+' '
                                site_mf = site_mf.strip()
                            else:
                                site_mf = mf_list[0]
                        except:
                            site_mf = ''
                
                    #put together intrument type+instrument_name+method_ref
                    mm = site_mm+site_mn+site_mf
                
                #get contact
                if row[0] == 'Originator:':
                    try:
                        contact_list = row[1:]
                        if len(contact_list) > 1:
                            site_contact = ''
                            for x in range(len(mf_list)):
                                site_contact = site_contact+contact_list[x]+' '
                            site_contact = site_contact.strip()
                        else:
                            site_contact = site_contact[0]
                    except:
                        site_contact = ''
                    all_contact.append(site_contact)
                
                #get country
                site_country = EMEP_COUNTRIES(file.split('/')[-1][:2])
                all_country.append(site_country)
                
                if row[0] == 'starttime':
                    skip_n = count+1
                    if species == 'ISOP':
                        spec_ind = row.index('C5H8')
                        try:
                            flag_ind = row.index('flag_C5H8')
                        except:
                            flag_ind = row.index('flag')
                    else:
                        spec_ind = row.index(species)
                        try:
                            flag_ind = row.index('flag_'+species)
                        except:
                            flag_ind = row.index('flag')
                    
                count+=1
            
        read = np.loadtxt(file,dtype="f8,f8,f8,f8",skiprows=skip_n,usecols=(0,1,spec_ind,flag_ind),unpack=True)
        read = np.array(read)
        times_since_start = read[0,:]
        endtimes_since_start = read[1,:]
        conc = read[2,:]
        conc = np.array(conc).astype('float64')
        flags = read[3,:]

        dates = []
        times = []
        enddates = []
        endtimes = []
        times_since_start = np.float64(times_since_start)   
        endtimes_since_start = np.float64(endtimes_since_start)  
        for x in range(len(times_since_start)):
            days_since_start = math.trunc(times_since_start[x])
            enddays_since_start = math.trunc(endtimes_since_start[x])
            remainder = times_since_start[x] - days_since_start
            remainder_end = endtimes_since_start[x] - enddays_since_start
            unrounded_hour = remainder*24
            unrounded_hour_end = remainder_end*24
            hour = np.round(unrounded_hour)
            hour_end = np.round(unrounded_hour_end)
            time_delta = datetime.timedelta(days = days_since_start,hours = hour)
            time_delta_end = datetime.timedelta(days = enddays_since_start,hours = hour_end)
            calc_datetime = start_datetime + time_delta
            calc_datetime_end = start_datetime + time_delta_end
            calc_yyyymmdd = calc_datetime.strftime("%Y%m%d") 
            calc_hhmm = calc_datetime.strftime("%H%M")  
            end_calc_yyyymmdd = calc_datetime_end.strftime("%Y%m%d") 
            end_calc_hhmm = calc_datetime_end.strftime("%H%M")
            dates.append(calc_yyyymmdd)
            times.append(calc_hhmm)
            enddates.append(end_calc_yyyymmdd)
            endtimes.append(end_calc_hhmm)
            
        conc = np.float64(conc)
        flags = np.float64(flags)
        
        #add to n_obs_all
        n_all += len(conc)
        
        #IF bad_meta == True then set all file vals as nans
        if bad_meta == True:
            conc[:] = np.NaN
        meta_valid_list.append(bad_meta)
        
        #DO INLINE INVALID AND FLAG CONVERT to NaN
        test = conc < 0
        conc[test] = np.NaN
        
        test = flags != 0
        conc[test] = np.NaN
            
        #convert units by line (only if value is >= than 0
        try:
            if (unit.lower() != 'ppb') & (unit.lower() != 'ppbv'):
                if unit == 'ug/m3':
                    #calculate conversion factor from mg/m3 assuming 293K and 1013 hPa - in EU LAW
                    #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                    conv_fact = 8.3144/mol_mass*(293.)/(1013./10.)
                    conc = conv_fact*conc
                elif unit == 'ug_N/m3':
                    conv_fact = 8.3144/14.00674*(293.)/(1013./10.)
                    conc = conv_fact*conc
                elif (unit == 'ppm') or (unit == 'ppmv'):
                    conc = conc*1e3
                    #print 'Converting Units from ppmv to ppbv'
                elif (unit == 'ppt') or (unit == 'pptv'):
                    conc = conc/1e3
                    #print 'Converting Units from pptv to ppbv'
                else:
                    print 'Unknown Unit'
                    1+'a'
        except:
            pass
        
        #remove 9.999 from ISOP dataset
        if species == 'ISOP':
            test = conc == 9.999
            conc[test] = np.NaN
        
        #if file resolution is daily or monthly then replicate times after point, to fill hourly data array.
        count=0
        if file_res == '1h':
            n_dups = np.zeros(len(conc))
        elif file_res == '1d':
            n_dups = []
            #if measurement method is flask, then put leave flask measurement in as hourly measurement, the first hour of month 
            file_hours = len(dates)
            for i in range(file_hours):
                current_year = int(dates[count][:4])
                current_month = int(dates[count][4:6])
                current_day = int(dates[count][6:])
                current_hh = int(times[count][:2])
                current_mm = int(times[count][2:])
        
                next_year = int(enddates[i][:4])
                next_month = int(enddates[i][4:6])
                next_day = int(enddates[i][6:])
                next_hh = int(endtimes[i][:2])
                next_mm =  int(endtimes[i][2:])
                
                s = datetime.datetime(year = current_year, month = current_month, day = current_day, hour = current_hh, minute = current_mm)
                e = datetime.datetime(year = next_year, month = next_month, day = next_day, hour = next_hh, minute = next_mm)
                day_dates = [d.strftime('%Y%m%d') for d in pd.date_range(s,e,freq='H')][1:-1]
                day_hours = [d.strftime('%H%M') for d in pd.date_range(s,e,freq='H')][1:-1]

                dates = np.insert(dates,count+1,day_dates)
                times = np.insert(times,count+1,day_hours)
                conc = np.insert(conc,count+1,[conc[count]]*len(day_dates))

                #append to n duplicated array
                n_dups=np.append(n_dups,0)
                n_dups=np.append(n_dups,[1]*len(day_dates))

                count +=(len(day_dates)+1)
        
        elif file_res == '1mo':
            n_dups = []
            #if measurement method is flask, then put leave flask measurement in as hourly measurement, the first hour of month 
            file_hours = len(dates)
            for i in range(file_hours):
                current_year = int(dates[count][:4])
                current_month = int(dates[count][4:6])
                current_day = int(dates[count][6:])
                current_hh = int(times[count][:2])
                current_mm = int(times[count][2:])
    
                next_year = int(enddates[i][:4])
                next_month = int(enddates[i][4:6])
                next_day = int(enddates[i][6:])
                next_hh = int(endtimes[i][:2])
                next_mm =  int(endtimes[i][2:])
    
                s = datetime.datetime(year = current_year, month = current_month, day = current_day, hour = current_hh, minute = current_mm)
                e = datetime.datetime(year = next_year, month = next_month, day = next_day, hour = next_hh, minute = next_mm)
        
                day_dates = [d.strftime('%Y%m%d') for d in pd.date_range(s,e,freq='H')][1:-1]
                day_hours = [d.strftime('%H%M') for d in pd.date_range(s,e,freq='H')][1:-1]
                dates = np.insert(dates,count+1,day_dates)
                times = np.insert(times,count+1,day_hours)
                conc = np.insert(conc,count+1,[conc[count]]*len(day_dates))
                
                #append to n duplicated array
                n_dups=np.append(n_dups,0)
                n_dups=np.append(n_dups,[1]*len(day_dates))
                
                count += (len(day_dates)+1)
        
        data = [dates,times,conc,n_dups]
        
        #put measurnement methods and into big list len of times
        mm_big=np.append(mm_big,[mm]*len(dates))
      
        try:
            big_list = np.hstack((big_list,data))
        except:
            big_list = np.array(data)
                
    if (y == year_array[-1]):    

        #get dates and times
        date_con = big_list[0,:]
        time_con = big_list[1,:]
          
        #get vals
        vals = np.array(big_list[2,:]).astype('float64')
        
        #get n dup array
        n_dup_array = np.array(big_list[3,:]).astype(float).astype(int)

        #if all files have missing key meta then exit
        if all(i == True for i in meta_valid_list) == True:
            inv_nometa += 1
            print 'Site Invalid. No Metadata for ref'
            if no2_type == 'MOLYBDENUM':
                n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod = 0,0,0,0,0,0,0,0,0,0,0,0,0
            exit_c_list = np.array([inv_nometa,inv_anyvaliddata,inv_nokeymeta,inv_resolution,inv_badmeasurementmethod])
            n_c_list = np.array([n_all,n_after_nometa,n_after_flagsandlod,n_after_duplicate,n_after_anyvaliddata,n_after_nokeymeta,n_after_resolution,n_after_badmeasurementmethod])
            unknown_list = [unknown_mm_list,unknown_mm_refs_list,unknown_local_tz_list]
            meta = ['na','na','na','na','na','na','na','na','na','na','na','na']
            exit_r = 'nometa'
            return c,['na'],['na'],['na'],False,meta,exit_c_list,n_c_list,unknown_list,exit_r,np.zeros(1)
        valid_hours_dup = np.sum(n_dup_array)
        n_after_nometa += (len(vals)-valid_hours_dup)

        #delete big list
        del big_list

        date_con = np.array(date_con).astype(int)
        time_con = np.array(time_con).astype(int)
        
        #remove data < 1970 and >= 2015
        test_inds = (date_con >= 19700101) & (date_con < 20150101)
        date_con = date_con[test_inds]
        time_con = time_con[test_inds]
        vals = vals[test_inds]
        mm_big = mm_big[test_inds]
        n_dup_array = n_dup_array[test_inds]
        
        #set st_big as 'continuous'
        st_big = ['continuous']*len(vals)
        
        #convert all Nans back to -99999
        test = np.isnan(vals)
        vals[test] = -99999
        
        #get obs valid
        test = vals >= 0
        valid_hours_dup = np.sum(n_dup_array[test])
        n_obs_valid = int(len(vals[test]) - valid_hours_dup)
        n_after_flagsandlod += n_obs_valid
        
        #create max possible species grid, measurement method and sampling type grids
        full_data = np.empty(n_hours)
        full_data_after_flagsandlod = np.empty(n_hours)
        big_n_dup_array = np.zeros(n_hours)
        full_data[:] = -99999
        full_data_after_flagsandlod[:] = -99999
        
        #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
        converted_time = modules.date_process(date_con,time_con,start_year)
        converted_time = np.round(converted_time,decimals=5)
        syn_grid_time = np.arange(0,n_days,1./24)
        syn_grid_time = np.round(syn_grid_time,decimals=5)
        #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
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
    
        #get mode of metadata
        try:
            lat = np.float32(stats.mode(all_latitudes)[0][0]) 
        except:
            lat = 'na'
        try:
            lon = np.float32(stats.mode(all_longitudes)[0][0])  
        except:
            lon = 'na'
        try:
            alt = np.float32(stats.mode(all_altitudes)[0][0]) 
        except:
            alt = 'na'
        unit = stats.mode(all_unit)[0][0]
        #remove empty strings from extra meta before mode test
        try:
            site_name = stats.mode(filter(None, all_site_name))[0][0]
        except:
            site_name = 'na'
        try:
            country = stats.mode(filter(None, all_country))[0][0]
        except:
            country = 'na'
        try:
            contact = stats.mode(filter(None, all_contact))[0][0] 
        except:
            contact = 'na'
    
        #set data tz - all EMEP times are UTC
        data_tz = 0
        all_tz = [data_tz]
    
        key_meta = [lat,lon,alt]
        
        #convert file res to standard format
        if file_res == '1h':
            file_res = 'H'
        elif file_res == '1d':
            file_res = 'D'
        elif file_res == '1mo':
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

        #set metadata not available as na
        raw_class_name = 'na'
    
        #set processed unit
        p_unit = 'ppbv'
    
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
