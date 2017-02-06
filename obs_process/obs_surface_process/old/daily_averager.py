import numpy as np
from netCDF4 import Dataset
import datetime
import pandas as pd

species = 'O3'
start_year = 2005
end_year = 2009

root_grp = Dataset('GLOBAL_SURFACE_%s_%s_%s_D.nc'%(species,start_year,end_year+1), 'w')
root_grp.description = 'Observational Daily Surface %s in ppb - Program written by Dene Bowdalo'%(species)

d0 = datetime.date(start_year, 1, 1)
d1 = datetime.date(end_year+1, 1, 1)
delta = d1 - d0
n_days = delta.days

root_grp.createDimension('date', n_days)
root_grp.createDimension('time', n_days)    
root_grp.createDimension('o3', n_days)

print n_days

#setup datetimes
start_dt = datetime.datetime(start_year,1,1,0,0)
end_dt = datetime.datetime(end_year,12,31,0,0)
dt_range = pd.date_range(start = start_dt,end = end_dt, freq = 'D')
dt_dates = []
dt_times = []
for i in dt_range:
    dt_dates.append(i.strftime("%Y%m%d"))
    dt_times.append(i.strftime("%H%M")) 

print len(dt_dates)

fname = 'GLOBAL_SURFACE_%s_%s_%s_H.nc'%(species,start_year,end_year+1)

data = Dataset(fname)
g = data.groups['mhd']

obs_date = g.variables['date'][:]
obs_time = g.variables['time'][:]
#----------------------------------------
valid_refs_dict = data.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

#----------------------------------------
#process obs dates and obs times to datetimes, then process pandas objects
year_val = []
month_val = []
day_val = []
hour_val = []
minute_val = []

#change obs times to datetimes
obs_date = obs_date.astype('str')
obs_time = obs_time.astype('str')

for date in obs_date:
    year_val.append(int(date[0:4]))
    month_val.append(int(date[4:6]))
    day_val.append(int(date[6:8]))

for time in obs_time:
    if np.float64(time) == 0:
        hour_val.append(0)
        minute_val.append(0)
    elif np.float64(time) == 30:
        hour_val.append(0)
        minute_val.append(30)
    else:
        hour_val.append(int(time[0:-2]))
        minute_val.append(int(time[-2:]))

obs_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]
obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
obs_key = lambda x: pd.Period(str(x.year)+'-'+str(x.month)+'-'+str(x.day))

for ref in valid_refs:
    print ref
    g = data.groups[ref]
    obs_var = g.variables['o3'][:]
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_var_pd = pd.Series(obs_var_mask, index=obs_time_pd)
    obs_lat = g.latitude
    obs_lon = g.longitude
    obs_alt = g.altitude
    obs_process = g.process_group
    obs_country = g.country
    obs_completeness = g.completeness
    obs_anthrome_class = g.anthrome_class
    obs_raw_class = g.raw_class
    obs_unit = g.unit

    #group the data
    group=obs_var_pd.groupby(obs_key) 

    #calculate the fraction of hours in the month with good data
    datafrac=group.count().div(group.size(),axis='rows')
    # calculate the mean
    obs_mean = group.mean()
    # now, zero out the mean data where datafrac<0.5
    obs_mean[datafrac<0.5]= -99999

    # calculate the monthly mean                                                                                                                                                 
    obs_monthly = group.mean()
    obs_monthly = np.array(obs_monthly)

    r = root_grp.createGroup('%s'%(ref))
    #set variables
    dates = r.createVariable('date', 'i8', ('date',))                                                                                                                          
    times = r.createVariable('time', 'i8', ('time',))
    o3 = r.createVariable('o3', 'f8', ('o3',))

    #set group attributes
    r.latitude = obs_lat
    r.longitude = obs_lon
    r.altitude = obs_alt
    r.process_group = obs_process
    r.country = obs_country
    r.completeness = obs_completeness
    r.anthrome_class = obs_anthrome_class
    r.raw_class = obs_raw_class
    r.unit = obs_unit

    dates[:] = dt_dates
    times[:] = dt_times
    o3[:] = obs_monthly
