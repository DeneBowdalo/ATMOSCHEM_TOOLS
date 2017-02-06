import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import lomb_phase
import modules
import datetime
from netCDF4 import Dataset
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
from scipy.stats import norm
import scipy.signal
import matplotlib.dates as dt
import os
import pandas as pd

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-2]
print '\nSpecies is %s\n'%(species)

#process model
#----------------------------------------
#read in model data

start_year = raw_input('Start Year?\n')
end_year = raw_input('\nEnd Year?\n')

#find valid_models
model_files = glob.glob('/work/home/db876/plotting_tools/model_files/*%s_%s_%s*'%(species,start_year,end_year))
model = []
vres = []
version = []
hres = []
met = []
timeres = []
additional = []

for f in model_files:
    f = f.replace("/work/home/db876/plotting_tools/model_files/", "")
    data_split = f.split('_')
    model = np.append(model,data_split[0])
    vres = np.append(vres,data_split[1])        
    version = np.append(version,data_split[5])        
    hres = np.append(hres,data_split[6])
    met = np.append(met,data_split[7])
    timeres = np.append(timeres,data_split[8])    
    additional = np.append(additional,data_split[9])

model_unique = np.unique(model)

if len(model_unique) == 1:
    model_chose = model_unique[0]
else:
    model_chose = raw_input('\nChoose Model\n%s\n'%('   '.join(i for i in model_unique))) 

inds = np.where((model == model_chose)) 

vres = vres[inds]
version = version[inds]
hres = hres[inds]
met = met[inds]
timeres = timeres[inds]
additional = additional[inds]

vres = np.unique(vres)
version = np.unique(version)
hres = np.unique(hres)
met = np.unique(met)
timeres = np.unique(timeres)
additional = np.unique(additional)

if len(vres) == 1:
    vres = vres[0]
else:
    vres = raw_input('\nChoose Vertical Resolution\n%s\n'%('   '.join(i for i in vres)))

if len(version) == 1:
    version = version[0]
else:
    version = raw_input('\nChoose Version\n%s\n'%('   '.join(i for i in version)))     

if len(hres) == 1:
    hres = hres[0]
else:
    hres = raw_input('\nChoose Horizontal Resolution\n%s\n'%('   '.join(i for i in hres))) 

if len(met) == 1:
    met = met[0]
else:
    met = raw_input('\nChoose Meteorology\n%s\n'%('   '.join(i for i in met))) 

if len(timeres) == 1:
    timeres = timeres[0]
else:
    timeres = raw_input('\nChoose Time Resolution\n%s\n'%('   '.join(i for i in timeres)))

if len(additional) == 1:
    additional = additional[0]
else:
    additional = raw_input('\nChoose Additional Detail\n%s\n'%('   '.join(i for i in additional)))

print '/work/home/db876/plotting_tools/model_files/'+model_chose+'_'+vres+'_'+species+'_'+start_year+'_'+end_year+'_'+version+'_'+hres+'_'+met+'_'+timeres+'_'+additional
model_fname = '/work/home/db876/plotting_tools/model_files/'+model_chose+'_'+vres+'_'+species+'_'+start_year+'_'+end_year+'_'+version+'_'+hres+'_'+met+'_'+timeres+'_'+additional

root_grp = Dataset(model_fname)
model_var = root_grp.variables[species.lower()][:]
model_date = root_grp.variables['date'][:]
model_time = root_grp.variables['time'][:]
lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
grid_size = root_grp.variables['grid_size'][:]
grid_size = grid_size[0]

gridbox_count = len(lat_c)*len(lon_c)

#process obs
#----------------------------------------
#data range

#read in obs data
root_grp = Dataset('../process/GLOBAL_%s_%s_1970_2015_%s_ALL.nc'%(vres,species,timeres))

valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

site_ref = raw_input('\nChoose site from list. Sites with full set of yearly files between %s & %s are:\n%s\n'%(start_year,end_year,'   '.join(i for i in valid_refs)))

#read in specific site data
site_group = root_grp.groups[site_ref]

#read in variables for site
obs_var = site_group.variables[species.lower()][:]
obs_date = site_group.variables['date'][:]
obs_time = site_group.variables['time'][:]
obs_lat = site_group.latitude
obs_lon = site_group.longitude
obs_alt = site_group.altitude
obs_group = site_group.process_group

obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)

#----------------------------------------
#find model data gridbox to compare with obs.

#get model gridbox for obs site
lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)

#print 'Obs Lat & Lon = ', obs_lat,obs_lon
#print 'Lat Centre = ',lat_c
#print 'Lon Centre = ',lon_c
#print lat_n,lon_n
#print lat_c[lat_n], lon_c[lon_n]
#print model_var.shape

model_var = model_var[:,lat_n,lon_n]
model_var = model_var*1e9

model_var_mask = np.ma.masked_where(model_var<=0,model_var)

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
obs_var_pd = pd.Series(obs_var_mask, index=obs_time_pd)


#process model dates and model times to datetimes, then process pandas objects

year_val = []
month_val = []
day_val = []
hour_val = []
minute_val = []

#change model times to datetimes
model_date = model_date.astype('str')
model_time = model_time.astype('str')

for date in model_date:
    year_val.append(int(date[0:4]))
    month_val.append(int(date[4:6]))
    day_val.append(int(date[6:8]))

for time in model_time:
    if np.float64(time) == 0:
        hour_val.append(0)
        minute_val.append(0)
    elif np.float64(time) == 30:
        hour_val.append(0)
        minute_val.append(30)
    else:
        hour_val.append(int(time[0:-2]))
        minute_val.append(int(time[-2:]))

model_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]

model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
model_var_pd = pd.Series(model_var_mask, index=model_time_pd)

#----------------------------------------
#set plotting area & background to white
fig=plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)

ax.grid(True)
#ax.xaxis.set_major_formatter(dt.DateFormatter('%d/%m/%Y'))
ax.plot_date(obs_time_pd.to_pydatetime(), obs_var_pd, color='black', markersize = 3, label = 'Observations')
ax.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='red', markersize = 3, label = '%s %s %s %s'%(model_chose,version,grid_size,met))

leg=ax.legend(loc=1, prop={'size':20})
leg.get_frame().set_alpha(0.4)
ax.set_xlabel('Time',fontsize = 15)
ax.set_ylabel('Concentration (ppb)',fontsize = 15,labelpad=40)
ax.set_title('Time Series of Surface %s at %s, for Observations & %s %s-%s \n'%(species,site_ref,model_chose,start_year,end_year),fontsize=18)
for tick in ax.get_xaxis().get_major_ticks():
    tick.set_pad(11.)

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500+100")
plt.show()

