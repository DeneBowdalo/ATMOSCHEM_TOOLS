import numpy as np
import matplotlib.pyplot as plt
import glob
import lomb_phase
import modules
import datetime
import csv
from netCDF4 import Dataset

#process obs
#----------------------------------------
#data range
first_year = 2006
last_year = 2011

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

#read in obs data
obs_data = np.genfromtxt('GAW_Site_Data.dat',delimiter=',',names=True,dtype=None)
obs_refs = obs_data['REF']
obs_site_names = obs_data['SITE_NAME']
obs_lats = obs_data['LAT']
obs_lons = obs_data['LON']

#find valid obs sites for data range
delete_array = []
for index in range(len(obs_refs)):
    obs_files = sorted(glob.glob('../O3_obs_2006_2011/%s*'%(obs_refs[index].lower())))

    #test to remove sites that do not have files for every year in data range

    #get years of files
    years = []    
    for i in obs_files:
        years.append(int(i[-8:-4]))

    #check if no files at all for site
    if len(years) == 0:
        delete_array = np.append(delete_array,index)
        continue
    if years[0] != first_year:
        delete_array = np.append(delete_array,index)
        continue
    if years[-1] != last_year:
        delete_array = np.append(delete_array,index)
        continue
           
    for current_year in years:    
        try:
            if current_year - previous_year != 1:
                delete_array = np.append(delete_array,index)
                break    
            previous_year = current_year
        except:
            previous_year = current_year
    del previous_year
               
valid_obs_site_names = np.delete(obs_site_names,delete_array)

#find n_hours and n_days between start and end date
d0 = datetime.date(first_year, 1, 1)
d1 = datetime.date(last_year+1, 1, 1)
delta = d1 - d0
n_days = delta.days
n_hours = n_days*24

#setup netcdf output
root_grp = Dataset('GAW_SURFACE_O3_%s_%s.nc'%(first_year,last_year+1), 'w')
root_grp.description = 'Hourly Surface O3 at GAW sites in ppb - Program written by Dene Bowdalo'

# dimensions
root_grp.createDimension('date', n_hours)
root_grp.createDimension('time', n_hours)    
root_grp.createDimension('o3', n_hours)

#for each valid location process
#limit obs data due for each site in valid_obs_site_names
start_array = [0,8760,17520,26304,35064,43824]
end_array = [8760,17520,26304,35064,43824,52584]

for site in valid_obs_site_names:
    site_ref = obs_refs[obs_site_names == site]
    site_ref = site_ref[0]
    
    print site_ref
    
    all_tz = []

    
    #read files for each valid site
    site_files = sorted(glob.glob('../O3_obs_2006_2011/%s*'%(site_ref.lower()))) 
    site_file_len = len(site_files)
    s_count = 0
    start_ind = 0
    end_ind = 0
    for f in site_files:
        print f
        read = np.loadtxt(f,dtype="S10,S5,f8",comments='C',usecols=(0,1,4),unpack =True) 	
        read = np.array(read)
                
        if s_count == 0:
            big_list = read
            end_ind = len(read[2,:])
        else:
            big_list = np.hstack((big_list,read))
            start_ind = end_ind
            end_ind+=len(read[2,:])
        s_count+=1
        units = [] 
        mycsv = csv.reader(open(f))
        row_count = 0
        for row in mycsv:
            # get units
            if row_count == 20:
                val = " ".join(row)
                unit = val.replace(" ", "")
                unit = unit[19:]
            row_count+=1
        print 'unit = ', unit    
        
        
        # test if units are in ppb for each file - if not convert from microg/m3	
        if (unit != 'ppb') & (unit != 'ppbv'):
            if unit == 'ug/m3': 
                print 'converting units, temp = 20degC'
                conc = read[2,:]
                conc = np.array(conc)
                conc = conc.astype(float)
                print conc[0:100]
                #calculate conversion factor from mg/m3 assuming 20 degC and 1 atm - default for GAW site O3 instruments
                #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                conv_fact = 8.3144/48*(273.15+20)/(1013.25/10)
                print conv_fact
                convert = conv_fact*conc
                print convert[0:100]
                big_list[2,start_ind:end_ind] = convert
            elif unit == 'ug/m3-25C':
                print 'converting units, temp = 25degC'
                conc = read[2,:]
                conc = np.array(conc)
                conc = conc.astype(float)
                print conc[0:100]
                #calculate conversion factor from mg/m3 assuming 25 degC and 1 atm
                #R/MW*(TEMP0C(K)*TEMP(degC)/P(hPa)/10
                conv_fact = 8.3144/48*(273.15+25)/(1013.25/10)
                print conv_fact
                convert = conv_fact*conc
                print convert[0:100]
                big_list[2,start_ind:end_ind] = convert
            
            else:
                print 'Unknown Unit'
                break
        
    # for each file, extract site info 
        mycsv = csv.reader(open(f))
        row_count = 0
        for row in mycsv:
            #get site name
            if row_count == 6:
                val = " ".join(row)
                site_name = val.replace(" ", "")
                site_name = site_name[15:]
            # get lat	
            if row_count == 11:
                val = " ".join(row)
                lat = val.replace(" ", "")
                lat = lat[12:]
                lat = float(lat)
            # get lon
            if row_count == 12:
                val = " ".join(row)
                lon = val.replace(" ", "")
                lon = lon[13:]
                lon = float(lon)
            # get altitude
            if row_count == 13:
                val = " ".join(row)
                alt = val.replace(" ", "")
                alt = alt[12:] 
                alt = float(alt)
            #get time zone
            if row_count == 23:
                val = " ".join(row)
                tz = val.replace(" ", "")
                tz = tz[12:]
            row_count+=1
    
        #Add each yearly time zone to array
        if tz != 'UTC':
            print site_ref
            if site_ref == 'PAY':
                tz = 1
            
            if site_ref == 'RCV':
                tz = 2
            
            if site_ref == 'RIG':
                tz = 1
            
            if site_ref == 'JFJ':
                tz = 1
            
            if site_ref == 'CMN':
                tz = 1
            
            if (tz == 'Localtime	-1') & (lon > -7.5) & (lon < 30):
                tz = 1
            
            if tz == '':
                tz = 0
                
            if type(tz) != int :    
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
                    1+1 
                tz = float(tz)
                tz = round(tz)
                tz = int(tz)
           
        else:
            tz = 0
        
        all_tz.append(tz)
        print all_tz
          
        if s_count == site_file_len:    
               
            #process dates from date, time to days since start year
            date_con = [s.replace('-', '') for s in big_list[0]]			
            #date_con = np.int32(date_con)
            time_con = [s.replace(':', '') for s in big_list[1]]
            #time_con = np.int32(time_con)
            
            print len(date_con)
            #some times go from 0100 to 2400, as opposed to 0000 to 2300, need to test if this is case, and if so correct
            for i in range(len(time_con)):
                if time_con[i] == '2400':
                    time_con[i] = '0000'
                    str_date = date_con[i]
                    date_datetime = datetime.datetime(int(str_date[0:4]),int(str_date[4:6]),int(str_date[6:]))
                    date_datetime = date_datetime + datetime.timedelta(days = 1)
                    str_date = date_datetime.strftime("%Y%m%d")
                    date_con[i] = int(str_date)
            
            last_date_val = int('%s1231'%(last_year))
            if int(date_con[-1]) > int(last_date_val):
                date_con = date_con[:-1]
                time_con = time_con[:-1]
            
            date_con = np.array(date_con).astype(int)
            time_con = np.array(time_con).astype(int)
            
            #create grid time array
            datetime_obj = datetime.datetime(year = first_year, month = 1, day = 1, hour = 0, minute = 0) 
            grid_dates = []
            grid_times = []
            print n_hours
            for i in range(n_hours):
                one_date = datetime_obj.strftime("%Y%m%d")
                one_time = datetime_obj.strftime("%H%M")
                grid_dates = np.append(grid_dates,one_date)
                grid_times = np.append(grid_times,one_time)
                datetime_obj = datetime_obj+datetime.timedelta(hours=1)
            
            #create max possible o3 grid
            o3_data = np.empty(n_hours)
            o3_data[:] = -99999
            
            print o3_data.shape
            print date_con[-25:]
            print time_con[-25:]
            
            #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
            converted_time = date_process(date_con,time_con,first_year)
            converted_time = np.round(converted_time,decimals=5)
            syn_grid_time = np.arange(0,n_days,1./24)
            syn_grid_time = np.round(syn_grid_time,decimals=5)
			#find matching times between actual times and grid of times, return big array of indices of matched indices in grid
            print converted_time[-25:]
            print syn_grid_time[-25:]
            indices = np.searchsorted(syn_grid_time, converted_time, side='left')
            o3_data[indices] = np.array(big_list[2]).astype(float) 
            for i in range(len(date_con)):
                obs_date = date_con[i]
                obs_time = time_con[i]
            
            
            #make timezone corrections    
            for count in range(len(all_tz)):
                tz = all_tz[count]
                start_ind = start_array[count]
                end_ind = end_array[count]
                o3_year = o3_data[start_ind:end_ind]
                print 'tz = ',tz
                
                print 'before = ', o3_year[-15:]
                
                if tz < 0:
                #get rid of values at start and append -99999's at end
                    cut = o3_year[:tz]
                    for num in range(np.abs(tz)):
                        cut = np.insert(cut,0, -99999)
                    try:
                        other_start_ind = start_ind-(np.abs(tz))
                        other_end_ind = start_ind
                        cut[:tz] = o3_data[other_start_ind:other_end_ind]
                    except:
                        pass
                    o3_data[start_ind:end_ind] = cut
                
                    
                elif tz > 0:
                #put -99999's at start and get rid of values at end
                    cut = o3_year[tz:]
                    print 'before = ', cut[-15:]
                    for num in range(tz):
                        cut = np.append(cut, -99999)
                    print 'tz = ',tz
                    print 'after = ', cut[-15:]
                    try:
                        other_start_ind = end_ind
                        other_end_ind = end_ind+tz
                        cut[-tz:] = o3_data[other_start_ind:other_end_ind]
                    except:
                        pass
                    o3_data[start_ind:end_ind] = cut
      


            #make sure all invalids set to -99999
            invalid_test = o3_data < 0
            o3_data[invalid_test] = -99999

            #save out netcdf file
            ref = root_grp.createGroup('%s'%(site_ref))
            
            #set variables
            dates = ref.createVariable('date', 'i8', ('date',))
            times = ref.createVariable('time', 'i8', ('time',))
            o3 = ref.createVariable('o3', 'f8', ('o3',))
            
            #set group attributes
            ref.site_name = site_name
            ref.latitude = lat
            ref.longitude = lon
            ref.altitude = alt
            ref.timezone = tz
    
            dates[:] = grid_dates
            times[:] = grid_times
            o3[:] = o3_data


root_grp.close()
            

