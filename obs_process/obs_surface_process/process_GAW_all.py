import numpy as np
import matplotlib.pyplot as plt
import glob
import lomb_phase
import modules
import datetime
import csv
from netCDF4 import Dataset
from time import sleep
import pandas as pd
from scipy import stats
import os
import multiprocessing
from netCDF4 import num2date, date2num
from tzwhere import tzwhere
import pytz
tz_root = tzwhere.tzwhere(shapely=True,forceTZ=True)

process_group = 'WMO GAW'

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

#read obs files
hourly_files = glob.glob('/work/home/db876/observations/surface/%s/GAW/*.hr*.dat'%(fname_species))
hourly_files = modules.natsorted(hourly_files)
daily_files = glob.glob('/work/home/db876/observations/surface/%s/GAW/*.da.*'%(fname_species))
daily_files = modules.natsorted(daily_files)
monthly_files = glob.glob('/work/home/db876/observations/surface/%s/GAW/*.mo.*'%(fname_species))
monthly_files = modules.natsorted(monthly_files)

#get all refs
ref_list_hourly = []
ref_list_daily = []
ref_list_monthly = []

valid_refs=[]

#get refs of different resolution
for i in range(len(hourly_files)):
    f = hourly_files[i]
    f = f.replace("/work/home/db876/observations/surface/%s/GAW/"%(fname_species), "")
    re = f[:3]
    if (re == 'poc') or (re == 'scs'):
        re = f[:7]
    ref_list_hourly.append(re)
    
for i in range(len(daily_files)):
    f = daily_files[i]
    f = f.replace("/work/home/db876/observations/surface/%s/GAW/"%(fname_species), "")
    re = f[:3]
    if (re == 'poc') or (re == 'scs'):
        re = f[:7]
    ref_list_daily.append(re)

for i in range(len(monthly_files)):
    f = monthly_files[i]
    f = f.replace("/work/home/db876/observations/surface/%s/GAW/"%(fname_species), "")
    re = f[:3]
    if (re == 'poc') or (re == 'scs'):
        re = f[:7]
    ref_list_monthly.append(re)  

#add hourly resolution files
refs = set(ref_list_hourly)
refs = [i for i in refs]
data_resolutions = ['hr']*len(refs)
    
if (output_res == 'HD') or (output_res == 'HDM'):
    #add daily resolution files
    ref_list_daily = set(ref_list_daily)
    ref_list_daily = [i for i in ref_list_daily]
    #remove refs if in ref list already
    for r in ref_list_daily:
        if r not in refs:
            refs.append(r) 
            data_resolutions.append('da')

if (output_res == 'HDM'):   
    #add monthly resolution files
    ref_list_monthly = set(ref_list_monthly)
    ref_list_monthly = [i for i in ref_list_monthly]
    # #remove refs if in ref list already
    for r in ref_list_monthly:
        if r not in refs:
            refs.append(r) 
            data_resolutions.append('mo')

#add in separate refs for simultaneous NO2 measurements 
if species == 'NO2':
    site_i = refs.index('rig')
    refs[site_i] = 'rig_molyb'
    refs.append('rig_photo')
    data_resolutions.append('hr')

refs, data_resolutions = zip(*sorted(zip(refs, data_resolutions)))

refs = np.array(refs)
print 'all refs len = ', len(refs)

valid_refs = np.copy(refs)

print valid_refs

#find n_hours and n_days between start and end date
d0 = datetime.date(start_year, 1, 1)
d1 = datetime.date(end_year+1, 1, 1)
delta = d1 - d0
n_years = (end_year+1)-start_year
n_months = n_years*12
n_days = delta.days
n_hours = n_days*24
n_points = n_hours

#create grid time array
start = datetime.datetime(year = start_year, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = end_year+1, month = 1, day = 1, hour = 0, minute = 0)
grid_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]
grid_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='H')][:-1]

#make output dates and times
output_res_dates_strings = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]
dt = pd.date_range(start,end,freq='H')[:-1].tolist()
output_res_times = date2num(dt, units='hours since 0001-01-01 00:00:00', calendar='gregorian')
obs_time_pd = pd.date_range(start = start,end = end, freq = 'H')[:-1]

#setup netcdf output
root_grp = Dataset('GAW_SURFACE_%s_%s_%s_%s.nc'%(fname_species,start_year,end_year+1,output_res), 'w')
root_grp.description = 'Hourly Surface %s at GAW sites - Program written by Dene Bowdalo'%(fname_species)

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

    #for each valid location process
    #limit obs data due for each site in valid_obs_site_names
    #for c in range(len(valid_refs)):
    
    all_lat = []
    all_lon = []
    all_alt = []
    all_tz = []
    all_unit = []
    all_site_name = []
    all_country = []
    all_contact = []
    mm_big = []
    st_big = []
    site_years = []

    site_ref = valid_refs[c]
    
    print 'ref = ',site_ref,c

    file_valid = True
    data_valid = True

    file_res = data_resolutions[c]

    #read files for each valid site
    #do some special handling for NO2
    if (species == 'NO2') & (site_ref[:3].lower() == 'rig'):
        m_type = site_ref[4:]
        s_files = sorted(glob.glob('/work/home/db876/observations/surface/%s/GAW/%s**%s.%s**.dat'%(fname_species,site_ref.lower()[:3],m_type,file_res))) 
    else:
        s_files = sorted(glob.glob('/work/home/db876/observations/surface/%s/GAW/%s**.%s**.dat'%(fname_species,site_ref.lower(),file_res))) 
    
    if file_res == 'hr':
        site_files = sorted(s_files, key = lambda x: x.split(".hr")[1])

    else:
        site_files = sorted(s_files)

    delete_inds = []
    if file_res == 'hr':
        #limit site files before and after year limit
        
        for i in range(len(site_files)):
            f = site_files[i]
            year = f.split(".hr")[1][:4]
            site_years.append(year)
            if int(year) < int(start_year):
                delete_inds.append(i)
            if int(year) > int(end_year):
                delete_inds.append(i)

        site_files = np.delete(site_files,delete_inds)
        print site_files

    site_file_len = len(site_files)
    s_count = 0
    start_ind = 0
    end_ind = 0
    for f in site_files:
    
        #get metadata
        units = [] 
        mycsv = csv.reader(open(f))
        row_count = 0
        for row in mycsv:
            if row_count == 11:
                val = " ".join(row)
                lat = val.replace(" ", "")
                lat = lat[12:]
                try:
                    lat = float(lat)
                except:
                    lat = ''
                all_lat.append(lat)
            # get lon
            if row_count == 12:
                val = " ".join(row)
                lon = val.replace(" ", "")
                lon = lon[13:]
                try:
                    lon = float(lon)
                except:
                   lon = ''
                all_lon.append(lon)
            # get altitude
            if row_count == 13:
                val = " ".join(row)
                alt = val.replace(" ", "")
                alt = alt[12:] 
                try:
                    alt = float(alt) 
                except:
                    alt = ''
                all_alt.append(alt)
            # get units
            if row_count == 20:
                val = " ".join(row)
                unit = val.replace(" ", "")
                unit = unit[19:]      
                all_unit.append(unit)
            # get measurement method
            if row_count == 21:
                val = " ".join(row)
                mm = val.replace(" ", "")
                mm = mm[21:]  
            # get sampling type
            if row_count == 22:
                val = " ".join(row)
                st = val.replace(" ", "")
                st = st[16:]  
            #get timezone
            if row_count == 23:
                val = " ".join(row)
                tz = val.replace(" ", "")
                tz = tz[12:]  
            #get site name
            if row_count == 6:
                val = " ".join(row)
                site_name = val.replace(" ", "")
                site_name = site_name.split(':')[-1]
                all_site_name.append(site_name)
            #get country
            if row_count == 9:
                val = " ".join(row)
                country = val.replace(" ", "")
                country = country.split(':')[-1]
                all_country.append(country)
            #get contact    
            if row_count == 16:
                val = " ".join(row)
                contact = val.replace(" ", "")
                contact = contact.split(':')[-1]
                all_contact.append(contact)
            #get flags index
            if " ".join(row).replace(" ", "")[3:11] == 'DATETIME':
                row_list = row[0].split()
                date_index = 0
                time_index = 1
                enddate_index = 2
                endtime_index = 3
                spec_index = 4
                try:
                    flag_index = row_list.index('F')-1
                except:
                    flag_index = -99999

            row_count+=1   
        
        print f
        #separate read lines if file has comments or does not.
        if flag_index != -99999:
            try:
                read = np.loadtxt(f,dtype="S10,S5,S10,S5,f8,f8",comments='C',usecols=(date_index,time_index,enddate_index,endtime_index,spec_index,flag_index),unpack =True) 	
            except:
                read = np.loadtxt(f,dtype="S10,S5,S10,S5,f8,S6",comments='C',usecols=(date_index,time_index,enddate_index,endtime_index,spec_index,flag_index),unpack =True)
        else:
            read = np.loadtxt(f,dtype="S10,S5,S10,S5,f8",comments='C',usecols=(date_index,time_index,enddate_index,endtime_index,spec_index),unpack =True)
        read = np.array(read)
    
        dates = read[0,:]
        times = read[1,:]
        enddates = read[2,:]
        endtimes = read[3,:]
        conc = read[4,:]
        conc = np.array(conc).astype('float64')
    
        #add to n_obs_all
        n_all += len(conc)
        n_after_nometa += len(conc)
    
        #change all vals < 0 to np.NaN
        inv_test = conc < 0
        conc[inv_test] = np.NaN
    
        #if have flags then remove invalid by flags not equal to zero
        if flag_index != -99999:
            flags = read[5,:]  
            try:
                flags = np.array(flags).astype('float64')
                inv_test = flags >= 200 
                conc[inv_test] = np.NaN
            except:
                flags = np.array(flags)
                inv_test = flags != 'V00'
                conc[inv_test] = np.NaN
            
        start_ind = end_ind
        end_ind+=len(conc)
    
        s_count+=1
     
        # test if units are in ppb for each file - if not convert
        if (unit != 'ppb') & (unit != 'ppbv'):
            if (unit == 'ug/m3') or (unit == 'ug/m3-20C'): 
                print 'converting units, temp = 20degC'
                #calculate conversion factor from mg/m3 assuming 20 degC and 1 atm - default for GAW site O3 instruments
                #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                conv_fact = 8.3144/mol_mass*(273.15+20)/(1013.25/10)
                conc = conv_fact*conc
            elif (unit == 'ugN/m3') or (unit == 'ugN/m3-20C'):
                print 'converting units, temp = 20degC, molecular mass N'
                #calculate conversion factor from mg/m3 assuming 20 degC and 1 atm - default for GAW site O3 instruments
                #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                conv_fact = 8.3144/14.00674*(273.15+20)/(1013.25/10)
                conc = conv_fact*conc
            elif (unit == 'ugN/m3-25C'):
                print 'converting units, temp = 25degC, molecular mass N'
                #calculate conversion factor from mg/m3 assuming 25 degC and 1 atm 
                #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                conv_fact = 8.3144/14.00674*(273.15+25)/(1013.25/10)
                conc = conv_fact*conc
            elif (unit == 'mgN/m3-20C'):
                print 'converting units, temp = 20degC, molecular mass N'
                #calculate conversion factor from mg/m3 assuming 20 degC and 1 atm - default for GAW site O3 instruments
                #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                conv_fact = 8.3144/14.00674*(273.15+20)/(1013.25/10)
                conc = (conv_fact*conc)*1e3
            elif (unit == 'mgN/m3-25C'):
                print 'converting units, temp = 25degC, molecular mass N'
                #calculate conversion factor from mg/m3 assuming 25 degC and 1 atm
                #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                conv_fact = 8.3144/14.00674*(273.15+25)/(1013.25/10)
                conc = (conv_fact*conc)*1e3
            elif (unit == 'ug/m3-25C') or (unit == 'ug/m3at25C'):
                print 'converting units, temp = 25degC'
                #calculate conversion factor from mg/m3 assuming 25 degC and 1 atm
                #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                conv_fact = 8.3144/mol_mass*(273.15+25)/(1013.25/10)
                conc = conv_fact*conc
            elif (unit == 'mg/m3-20C'):
                print 'converting units, temp = 20degC'
                #calculate conversion factor from mg/m3 assuming 20 degC and 1 atm
                #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                conv_fact = 8.3144/mol_mass*(273.15+20)/(1013.25/10)
                conc = (conv_fact*conc)*1e3
            elif (unit == 'mg/m3-25C'):
                print 'converting units, temp = 25degC'
                #calculate conversion factor from mg/m3 assuming 25 degC and 1 atm
                #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                conv_fact = 8.3144/mol_mass*(273.15+25)/(1013.25/10)
                conc = (conv_fact*conc)*1e3
            elif (unit == 'ppm') or (unit == 'ppmv'):
                conc = conc*1.e3
            elif (unit == 'ppt') or (unit == 'pptv'):
                conc = conc/1.e3
        
            else:
                print 'Unknown Unit'
                print unit
                1+'a'
                break
            
        #do some manual changes of timezones for GAW sites that do not have timezone metadata
        if tz != 'UTC':
            if tz == '':
                if site_ref.lower() in ['plm']:
                    tz = -5
        
                if site_ref.lower() in ['kos','edm','vdl','nwr']:
                    tz = 0

                if site_ref.lower() in ['jfj','kps','rig','pay','glh','cmn','zep','dig','hhe','ktb','stp','ivn','jcz','kam','lzp','snz','zbl','kmw','don','mhn','nia','roq','spm']: 
                    tz = 1

                if site_ref.lower() in ['rcv','aht','oul','uto','vir','fdt','sem','stn']:
                    tz = 2
                
                if site_ref.lower() in ['dak']:
                    tz = 3
                
                if site_ref.lower() in ['shp']:
                    tz = 4
                    
                if site_ref.lower() in ['isk']:
                    tz = 5
    
                if site_ref.lower() in ['hkg']:
                    tz = 8

                if site_ref.lower() in ['cgo']:
                    tz = 10
            else:        
                tz = tz.replace('LocaltimeUTC', '')
                tz = tz.replace('OtherUTC', '')
                tz = tz.replace('Localtime', '')
                tz = tz.replace(':', '.')
        
                try:
                    before, sep, after = tz.rpartiton('.')
                    after = int(after)
                    conv = (100./60) * after
                    tz = before+sep+str(conv)
                except:
                    pass 
        
        else: 
            tz = 0
    
        tz = float(tz)
        all_tz.append(tz)
    
        #process dates from date, time to days since start year
        dates = [s.replace('-', '') for s in dates]			
        times = [s.replace(':', '') for s in times]
        enddates = [s.replace('-', '') for s in enddates]			
        endtimes = [s.replace(':', '') for s in endtimes]
    
        if file_res == 'hr':
            #some times go from 0100 to 2400, assume this is when sites report ave for hour previous. Thus all times should have hour minused
            for i in range(len(times)):
                if times[i] == '2400':
                    current_date = dates[i]
                    test = np.array(dates) == current_date
                    indices = [i for i, x in enumerate(test) if x]
                    for x in indices:
                        current_time = times[x]
                        if current_time == '2400':
                            current_time = '0000'
                        date_datetime = datetime.datetime(int(current_date[0:4]),int(current_date[4:6]),int(current_date[6:]),int(current_time[:2]),int(current_time[2:]))
                        date_datetime = date_datetime - datetime.timedelta(hours = 1)
                        times[x] = date_datetime.strftime("%H%M")
    
            #adjust dates and times if tz is not equal to 0
            if tz != 0:
                for i in range(len(dates)):
                    #create datetime
                    dt = datetime.datetime(int(dates[i][:4]),int(dates[i][4:6]),int(dates[i][6:]),int(times[i][:2]),int(times[i][2:]))
                    if tz > 0:
                        dt  = dt - datetime.timedelta(hours = int(tz))
                    elif tz < 0:
                        dt  = dt + datetime.timedelta(hours = np.abs(int(tz)))
                    dates[i] = dt.strftime("%Y%m%d")
                    times[i] = dt.strftime("%H%M")
            
        data = [dates,times,enddates,endtimes,conc]
        
        #put measurnement methods and sampling types into big lists len of times
        mm_big=np.append(mm_big,[mm]*len(dates))
        st_big=np.append(st_big,[st]*len(dates))
        
        try:
            big_list = np.hstack((big_list,data))
        except:
            big_list = np.array(data)    
            
        #start write out process 
        if (s_count == site_file_len):	
        
            if len(site_years)!=len(set(site_years)):
                print 'Site has duplicate years. Breaking'
                1+'a'
            
            #get dates and times
            date_con = big_list[0,:]
            time_con = big_list[1,:]
            enddate_con = big_list[2,:]
            endtime_con = big_list[3,:]
          
            #get vals
            vals = np.array(big_list[4,:]).astype(float) 

            #delete big list
            del big_list
        
            #if file resolution is daily or monthly then replicate times after point, to fill hourly data array.
            count=0
            if file_res == 'hr':
                n_dup_array = np.zeros(len(vals))
            elif file_res == 'da':
            #if measurement method is flask, then put leave flask measurement in as hourly measurement, the first hour of day 
                if st == 'flask':
                    n_dup_array = np.zeros(len(vals))
                else:
                    n_dup_array = []
                    file_hours = len(date_con)
                    for i in range(file_hours):
                        current_year = int(date_con[count][:4])
                        current_month = int(date_con[count][4:6])
                        current_day = int(date_con[count][6:])
                        current_hh = int(time_con[count][:2])
                        current_mm = int(time_con[count][2:])
                        
                        next_year = int(enddate_con[i][:4])
                        next_month = int(enddate_con[i][4:6])
                        next_day = int(enddate_con[i][6:])
                        next_hh = int(endtime_con[i][:2])
                        next_mm =  int(endtime_con[i][2:])
                        
                        if (next_year == 9999) or (next_month == 99) or (next_day == 99) or (next_hh == 99) or (next_mm == 99):
                            next_dt = datetime.datetime(current_year,current_month,current_day,current_hh,current_mm)+datetime.timedelta(hours=24)
                            next_year = int(next_dt.strftime('%Y%m%d')[:4])
                            next_month = int(next_dt.strftime('%Y%m%d')[4:6])
                            next_day = int(next_dt.strftime('%Y%m%d')[6:])
                            next_hh = int(next_dt.strftime('%H%M')[:2])
                            next_mm = int(next_dt.strftime('%H%M')[2:])
                        
                        s = datetime.datetime(year = current_year, month = current_month, day = current_day, hour = current_hh, minute = current_mm)
                        e = datetime.datetime(year = next_year, month = next_month, day = next_day, hour = next_hh, minute = next_mm)
                        
                        day_dates = [d.strftime('%Y%m%d') for d in pd.date_range(s,e,freq='H')][1:-1]
                        day_hours = [d.strftime('%H%M') for d in pd.date_range(s,e,freq='H')][1:-1]
                        date_con = np.insert(date_con,count+1,day_dates)
                        time_con = np.insert(time_con,count+1,day_hours)
                        vals = np.insert(vals,count+1,[vals[count]]*len(day_dates))
                        mm_big = np.insert(mm_big,count+1,[mm_big[count]]*len(day_dates))
                        st_big = np.insert(st_big,count+1,[st_big[count]]*len(day_dates))
           
                        #append to n duplicated array
                        n_dup_array=np.append(n_dup_array,0)
                        n_dup_array=np.append(n_dup_array,[1]*len(day_dates))
           
                        count += (len(day_dates)+1)
        
            elif file_res == 'mo':
                #if measurement method is flask, then put leave flask measurement in as hourly measurement, the first hour of month 
                if st == 'flask':
                    n_dup_array = np.zeros(len(vals))
                else:
                    n_dup_array = []
                    file_hours = len(date_con)
                    for i in range(file_hours):
                        current_year = int(date_con[count][:4])
                        current_month = int(date_con[count][4:6])
                        current_day = int(date_con[count][6:])
                        current_hh = int(time_con[count][:2])
                        current_mm = int(time_con[count][2:])
                        
                        next_year = int(enddate_con[i][:4])
                        next_month = int(enddate_con[i][4:6])
                        next_day = int(enddate_con[i][6:])
                        next_hh = int(endtime_con[i][:2])
                        next_mm =  int(endtime_con[i][2:])
                        
                        if (next_year == 9999) or (next_month == 99) or (next_day == 99) or (next_hh == 99) or (next_mm == 99):
                            next_month = current_month+1
                            if next_month > 12:
                                next_month = 1
                                next_year = current_year+1
                            else:
                                next_year = current_year 
                            next_day = 1
                            next_hh = 0
                            next_mm = 0
                            
                        s = datetime.datetime(year = current_year, month = current_month, day = current_day, hour = current_hh, minute = current_mm)
                        e = datetime.datetime(year = next_year, month = next_month, day = next_day, hour = next_hh, minute = next_mm)
            
                        day_dates = [d.strftime('%Y%m%d') for d in pd.date_range(s,e,freq='H')][1:-1]
                        day_hours = [d.strftime('%H%M') for d in pd.date_range(s,e,freq='H')][1:-1]
                        date_con = np.insert(date_con,count+1,day_dates)
                        time_con = np.insert(time_con,count+1,day_hours)
                        vals = np.insert(vals,count+1,[vals[count]]*len(day_dates))
                        mm_big = np.insert(mm_big,count+1,[mm_big[count]]*len(day_dates))
                        st_big = np.insert(st_big,count+1,[st_big[count]]*len(day_dates))
                        
                        #append to n duplicated array
                        n_dup_array=np.append(n_dup_array,0)
                        n_dup_array=np.append(n_dup_array,[1]*len(day_dates))
                        
                        count += (len(day_dates)+1)
    
            date_con = np.array(date_con).astype(int)
            time_con = np.array(time_con).astype(int)
            
            #remove data < 1970 and >= 2015
            test_inds = (date_con >= 19700101) & (date_con < 20150101)
            date_con = date_con[test_inds]
            time_con = time_con[test_inds]
            vals = vals[test_inds]
            st_big = st_big[test_inds]
            mm_big = mm_big[test_inds]
            n_dup_array = n_dup_array[test_inds]
            
            #by species, make sure mm is correct for sites with missing data (instruments taken from gawsis)
            if species == 'O3':
                if site_ref.upper() in ['AHT','BHD','CGO','CMN','DAK','DIG','DON','EDM','FDT','GLH','HKG','ICE','ISK','JFJ','KOS','KPS','MHN','MKN','NIA','NMY','NWR','OUL','PAY','RCV','RIG','ROQ','SHP','SPM','TAR','USH','UTO','VDL','VIR','ZEP']:
                    mm_big = ['ultraviolet photometry']*len(vals)
                elif site_ref.upper() in ['BRW','MLO','SPO']:
                    test1 = date_con < 19760101
                    test2 = date_con >= 19760101
                    mm_big_1 = ['electrochemical concentration cell']*len(vals[test1])
                    mm_big_2 = ['ultraviolet photometry']*len(vals[test2])
                    mm_big = np.append(mm_big_1,mm_big_2)
                elif site_ref.upper() in ['HPB']:
                    test1 = date_con < 19870101
                    test2 = (date_con >= 19870101) & (date_con < 19880601)
                    test3 = date_con >= 19880601
                    mm_big_1 = ['potassium iodide']*len(vals[test1])
                    mm_big_2 = ['ethylene chemiluminescence']*len(vals[test2])
                    mm_big_3 = ['ultraviolet photometry']*len(vals[test3])
                    mm_big = np.append(mm_big_1,mm_big_2)
                    mm_big = np.append(mm_big,mm_big_3)
                elif site_ref.upper() in ['NMY']:
                    test1 = date_con < 19950101
                    test2 = date_con >= 19950101
                    mm_big_1 = ['electrochemical concentration cell']*len(vals[test1])
                    mm_big_2 = ['ultraviolet photometry']*len(vals[test2])
                    mm_big = np.append(mm_big_1,mm_big_2)  
                elif site_ref.upper() in ['SMO']:
                    test1 = date_con < 19810101
                    test2 = date_con >= 19810101
                    mm_big_1 = ['electrochemical concentration cell']*len(vals[test1])
                    mm_big_2 = ['ultraviolet photometry']*len(vals[test2])
                    mm_big = np.append(mm_big_1,mm_big_2) 
                    
            elif species == 'NO':
                if site_ref.upper() in ['HHE','JFJ','KOS','KTB','PAY','RIG','SSL','STP','ZSF','ZUG']:
                    mm_big = ['chemiluminescence']*len(vals)
            
            elif species == 'NO2':
                if site_ref.upper() in ['BEO','DON','GLH','HHE','IVN','KMW','KOS','KTB','PAY','RIG_MOLYB','SNB','STP']:
                    mm_big = ['chemiluminescence (conversion-molybdenum)']*len(vals)
                elif site_ref.upper() in ['CMN','CVO','HPB','JFJ','RIG_PHOTO','SSL','ZSF']:
                    mm_big = ['chemiluminescence (conversion-photolysis)']*len(vals)
                elif site_ref.upper() in ['CAR','FDT','LOG','LZP']:
                    mm_big = ['colorimetric']*len(vals)
                elif site_ref.upper() in ['BUR','DIG','IRB','JCZ','KAM','KPS','LEB','PLD','PLM','PLV','RCV','SNZ','SOF','SWL','VRN','ZBL','ZSN']:
                    mm_big = ['spectrophotometry']*len(vals)
                elif site_ref.upper() in ['BKT','JKR']:
                    mm_big = ['ion chromatography']*len(vals)
                elif site_ref.upper() in ['SEM','STN']:  
                    mm_big = ['diffusive sampler']*len(vals)
                elif site_ref.upper() in ['MHN','NIA','ROQ','SPM']:
                    test1 = date_con < 20040201
                    test2 = date_con >= 20040201
                    mm_big_1 = ['colorimetric']*len(vals[test1])
                    mm_big_2 = ['chemiluminescence (conversion-molybdenum)']*len(vals[test2])
                    mm_big = np.append(mm_big_1,mm_big_2)
                elif site_ref.upper() in ['AHT','OUL']:
                    test1 = date_con < 19970101
                    test2 = date_con >= 19970101
                    mm_big_1 = ['colorimetric']*len(vals[test1])
                    mm_big_2 = ['chemiluminescence (conversion-molybdenum)']*len(vals[test2])
                    mm_big = np.append(mm_big_1,mm_big_2)
                elif site_ref.upper() in ['UTO']:
                    test1 = date_con < 19940101
                    test2 = date_con >= 19940101
                    mm_big_1 = ['colorimetric']*len(vals[test1])
                    mm_big_2 = ['chemiluminescence (conversion-molybdenum)']*len(vals[test2])
                    mm_big = np.append(mm_big_1,mm_big_2)
                elif site_ref.upper() in ['VIR']:
                    test1 = date_con < 19950201
                    test2 = date_con >= 19950201
                    mm_big_1 = ['colorimetric']*len(vals[test1])
                    mm_big_2 = ['chemiluminescence (conversion-molybdenum)']*len(vals[test2])
                    mm_big = np.append(mm_big_1,mm_big_2)
            
            elif species == 'CO':
                if site_ref.upper() in ['HKG','KTB','MKN','ZSF']:
                    mm_big = ['non-dispersive infrared spectroscopy']*len(vals)
                elif site_ref.upper() in ['POC935N','POC935S']:
                    mm_big = ['gas chromatography reduction gas detection']*len(vals)
                elif site_ref.upper() in ['SSL']:    
                    mm_big = ['vacuum ultraviolet fluorescence']*len(vals)
            
            elif species == 'ISOP':
                #hpb only samples once at midday every day, so set as bad mm by making instrumentation blank.
                if site_ref.upper() in ['HPB']:
                    mm_big = ['']*len(vals)
                    
            #convert nans to -99999's
            nan_inds = np.isnan(vals)
            vals[nan_inds] = -99999
            
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
                lat = np.float32(stats.mode(all_lat)[0][0]) 
            except:
                lat = 'na'
            try:
                lon = np.float32(stats.mode(all_lon)[0][0])  
            except:
                lon = 'na'
            try:
                alt = np.float32(stats.mode(all_alt)[0][0]) 
            except:
                alt = 'na'
            unit = stats.mode(all_unit)[0][0]
            data_tz = np.float32(stats.mode(all_tz)[0][0])
            #remove empty strings from extra meta before mode test
            site_name = stats.mode(filter(None, all_site_name))[0][0]
            country = stats.mode(filter(None, all_country))[0][0]
            contact = stats.mode(filter(None, all_contact))[0][0]    
            
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
    
            #make tz int after checks
            data_tz = np.float32(data_tz)
            
            #set metadata not available as na
            raw_class_name = 'na'

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

