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
import scipy.signal
import matplotlib.dates as dt

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
obs_name = site_group.site_name
obs_lat = site_group.latitude
obs_lon = site_group.longitude
obs_alt = site_group.altitude
obs_timezone = site_group.timezone

#cut out invalid obs data
valids = obs_var > 0
obs_var = obs_var[valids]

#process model dates and model times to datetimes 
year_val = []
month_val = []
day_val = []
hour_val = []
minute_val = []

#change model times to datetimes
obs_date = obs_date.astype('str')
obs_time = obs_time.astype('str')

for date in obs_date:
    year_val.append(int(date[0:4]))
    month_val.append(int(date[4:6]))
    day_val.append(int(date[6:8]))

for time in obs_time:
    if time == '0':
        hour_val.append(0)
        minute_val.append(0)
    else:
        hour_val.append(int(time[0:-2]))
        minute_val.append(int(time[-2:]))

model_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]

#change obs times to datetimes

year_val = []
month_val = []
day_val = []
hour_val = []
minute_val = []

obs_date = obs_date[valids]
obs_time = obs_time[valids]

for date in obs_date:
    year_val.append(int(date[0:4]))
    month_val.append(int(date[4:6]))
    day_val.append(int(date[6:8]))

for time in obs_time:
    if time == '0':
        hour_val.append(0)
        minute_val.append(0)
    else:
        hour_val.append(int(time[0:-2]))
        minute_val.append(int(time[-2:]))

obs_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]

#obs_time = modules.date_process(obs_date,obs_time)

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

print gridbox_n

model_var = model_data[gridbox_n::gridbox_count]
model_var = model_var*1e9

#----------------------------------------
#set plotting area & background to white
#fig=plt.figure())

fig, (ax, ax2) = plt.subplots(2, sharex=True,figsize=(20,12))
fig.patch.set_facecolor('white')

ax.grid(True)
ax.xaxis.set_major_formatter(dt.DateFormatter('%d/%m/%Y'))
ax.plot(obs_datetimes,obs_var, color='black',label = 'Obs.')
ax.plot(model_datetimes,model_var, color='red',label = 'GEOS %s'%(model_version))
leg=ax.legend(loc=1, prop={'size':21})
leg.get_frame().set_alpha(0.4)
#ax.set_xlabel('Time',fontsize = 15)
ax.set_ylabel('Concentration (ppb)',fontsize = 21)

#ax2 = fig.add_subplot(2,1,2)

all_time = np.arange(0,365.25*4,1./24)
all_time = np.append(all_time,np.arange(365.25*4,(365.25*4)+365*2,1./24))

obs_annual_amplitude = 5.41153541631
obs_hannual_amplitude = 3.70229249913
obs_daily_amplitude = 1.18745139921
obs_annual_phase = 1.35696303697
obs_hannual_phase = -2.49892975397
obs_daily_phase = -1.94122014668

model_annual_amplitude = 2.79727030284
model_hannual_amplitude = 2.26023742293
model_daily_amplitude = 1.65568478096
model_annual_phase = 1.24157192294
model_hannual_phase = -2.5757850734
model_daily_phase = -2.1170157038

obs_annual_wave = obs_annual_amplitude*(np.cos(((2*np.pi)*all_time/365.25)-(obs_annual_phase)))
obs_hannual_wave = obs_hannual_amplitude*(np.cos(((2*np.pi)*all_time/182.625)-(obs_hannual_phase)))
obs_daily_wave = obs_daily_amplitude*(np.cos(((2*np.pi)*all_time/1.)-(obs_daily_phase)))

model_annual_wave = model_annual_amplitude*(np.cos(((2*np.pi)*all_time/365.25)-(model_annual_phase)))
model_hannual_wave = model_hannual_amplitude*(np.cos(((2*np.pi)*all_time/182.625)-(model_hannual_phase)))
model_daily_wave = model_daily_amplitude*(np.cos(((2*np.pi)*all_time/1.)-(model_daily_phase)))

obs_wave = np.mean(obs_var) + (obs_daily_wave)
model_wave = np.mean(model_var) + (model_daily_wave) 

#obs_wave = np.mean(obs_var) + (obs_annual_wave + obs_hannual_wave + obs_daily_wave)
#model_wave = np.mean(model_var) + (model_annual_wave + model_hannual_wave + model_daily_wave)
ax2.xaxis.set_major_formatter(dt.DateFormatter('%d/%m/%Y'))
ax2.plot(obs_datetimes,obs_var, color='black',alpha=0.3)
ax2.plot(model_datetimes,model_var, color='red', alpha=0.3)
ax2.plot(model_datetimes,obs_wave,color='black',label='Obs. Sinusoidal Waveform',linewidth=2)
ax2.plot(model_datetimes,model_wave,color='red',label='GEOS 2x2.5 Sinusoidal Waveform',linewidth=2)

#convert phase to time
obs_phase = modules.convert_phase_units_actual_single(obs_annual_phase,12)
model_phase = modules.convert_phase_units_actual_single(model_annual_phase,12)

print obs_phase
print model_phase

ax2.grid(True)
leg=ax2.legend(loc=1, prop={'size':21})
leg.get_frame().set_alpha(0.4)

ax2.set_xlabel('Time',fontsize = 21)
ax2.set_ylabel('Concentration (ppb)',fontsize = 21)
#ax.set_title(r'Time Series of Surface $O_3$ at %s, for Obs. & Model'%(site),fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=18)
ax2.tick_params(axis='both', which='major', labelsize=18)
#for tick in ax.get_xaxis().get_major_ticks():
#    tick.set_pad(11.)
#    #tick.label1 = tick._get_text1()
for tick in ax2.get_xaxis().get_major_ticks():
    tick.set_pad(11.)
    
ax.text(0.025, 0.9,'a', fontsize = 34,fontweight='bold',ha='center', va='center',transform=ax.transAxes)
ax2.text(0.025, 0.9,'b', fontsize = 34, fontweight='bold',ha='center', va='center',transform=ax2.transAxes)

fig.tight_layout()

plt.savefig('ts_macehead_sig_sinusoids.png')
plt.show()
