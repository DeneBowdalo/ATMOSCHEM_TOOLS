import numpy as np
import matplotlib.pyplot as plt
import glob
import lomb_phase
import modules
import datetime
from netCDF4 import Dataset
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
from scipy.stats import norm
import random

#process obs
#----------------------------------------
#data range
first_year = 2006
last_year = 2011

#read in obs data
root_grp = Dataset('GAW_SURFACE_O3_2006_2012.nc')
valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)
    
#read gaw site data file
obs_data = np.genfromtxt('GAW_Site_Data.dat',delimiter=',',names=True,dtype=None)
obs_refs = obs_data['REF']
obs_site_names = obs_data['SITE_NAME']

#get valid site names form refs from netcdf file
for i in range(len(obs_refs)):
    for j in valid_refs:
        if obs_refs[i] == j:
            try:
                ref_test.append(i)
            except:
                ref_test = [i]   
valid_obs_site_names = obs_site_names[ref_test]
               
site = raw_input('Choose site from list. Sites with full set of yearly files between %i & %i are:\n%s\n'%(first_year,last_year+1,'   '.join(i for i in valid_obs_site_names)))

#limit obs data due for site
site_ref = obs_refs[obs_site_names == site]
site_ref = site_ref[0]

#read in specific site data
site_group = root_grp.groups[site_ref]

#read in variables for site
obs_var = site_group.variables['o3'][:]
obs_date = site_group.variables['date'][:]
obs_time = site_group.variables['time'][:]
obs_lat = site_group.variables['lat'][:]
obs_lon = site_group.variables['lon'][:]
obs_alt = site_group.variables['alt'][:]
obs_timezone = site_group.variables['timezone'][:]

#process obs dates and obs_times to datetimes
#year_val = []
#month_val = []
#day_val = []
#hour_val = []
#minute_val = []

#for date in dates:
#    year_val.append(int(date[0:4]))
#    month_val.append(int(date[4:6]))
#    day_val.append(int(date[6:8]))

#for time in times:
#    hour_val.append(int(time[0:2]))
#    minute_val.append(int(time[2:4]))

#model_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]

#obs_datetimes = model_datetimes[valids]
#ax1.xaxis.set_major_formatter(dt.DateFormatter('%d/%m\n%H:%M:%S'))

#cut out invalid obs data
valids = obs_var > 0
obs_var = obs_var[valids]
obs_date = obs_date[valids]
obs_time = obs_time[valids]

obs_time = modules.date_process(obs_date,obs_time)

          
#process model
#----------------------------------------

#read in model data
model_dict = {'4x5':'binary_logs/4x5_GRID_O3.npy','2x2.5':'binary_logs/2x2.5_GRID_O3.npy'}
model_version = raw_input('\nChoose Model Version.\n%s\n'%('   '.join(i for i in model_dict)))
model_f = model_dict[model_version]

model_data = read = np.load(model_f)
if model_version == '2x2.5':
    model_time = np.arange(0,2191,1./24)
else:
    model_time = np.arange(0,2190,1./24)

#get model grid dims. for sim. type
lat_c,lat_e,lon_c,lon_e = modules.model_grids(model_version)
gridbox_count = len(lat_c)*len(lon_c)

#get model gridbox for obs site
gridbox_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)

model_var = model_data[gridbox_n::gridbox_count]
model_var = model_var*1e9

#----------------------------------------

#align obs and model periods
#obs_time,model_time,obs_var,model_var = modules.align_periods(obs_time,model_time,obs_var,model_var)

#take lomb for obs and model
ofac = raw_input('\nChoose oversampling factor. (Typically 4.)\n')

#obs lomb
obs_periods,obs_mag,obs_ph,obs_sig,obs_i_freq = modules.take_lomb(obs_time,obs_var,ofac)

obs_noise_count = [0]*len(obs_mag)

for i in range(10):
    #randomly mix arrays 
    obs_rand_i = random.sample(range(len(obs_var)), len(obs_var))
    model_rand_i = random.sample(range(len(model_var)), len(model_var))

    obs_var_rand = obs_var[obs_rand_i]
    model_var_rand = model_var[model_rand_i]

    #obs random lomb
    obs_periods_rand,obs_mag_rand,obs_ph_rand,obs_sig_rand,obs_i_freq_rand = modules.take_lomb(obs_time,obs_var_rand,ofac)
    for i in range(len(obs_mag)):
        if obs_mag_rand[i] > obs_mag[i]:
            obs_noise_count[i]=obs_noise_count[i]+1
    print obs_noise_count[0]

fap = [] 
for i in obs_noise_count:
    fap.append(i/10.)

a = fap == 0
b = fap[b]
print len(b)


#set plotting area & background to white
fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)

plt.semilogx(obs_periods,fap)

#plt.loglog(obs_periods,obs_mag, color='black', label = 'Obs.')
#plt.loglog(model_periods,model_mag, color='red', alpha=0.5, label = 'GEOS %s'%(model_version))

plt.show()
