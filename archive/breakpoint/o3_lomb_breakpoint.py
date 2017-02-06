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

#process obs
#----------------------------------------
#data range
first_year = 2006
last_year = 2011

#read in obs data
root_grp = Dataset('/work/home/db876/plotting_tools/GAW_SURFACE_O3_2006_2012.nc')
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
model_dict = {'4x5':'/work/home/db876/plotting_tools/binary_logs/4x5_GRID_O3.npy','2x2.5':'/work/home/db876/plotting_tools/binary_logs/2x2.5_GRID_O3.npy'}
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

gap_obs = raw_input('\nDo you wish to impose Observational Gaps on Model data? y or n?\n') 

if gap_obs == 'n':
    #align obs and model periods
    obs_time,model_time,obs_var,model_var = modules.align_periods(obs_time,model_time,obs_var,model_var)
else:
    if model_version == '4x5':
        valids = valids[24:]
        model_time = model_time[valids]
        model_var = model_var[valids]
    else:
        model_time = model_time[valids]
        model_var = model_var[valids]


#take lomb for obs and model
ofac = raw_input('\nChoose oversampling factor. (Typically 4.)\n')

#obs lomb
obs_periods,obs_mag,obs_ph,obs_fr,obs_fi = modules.take_lomb(obs_time,obs_var,ofac,1./24)

#model lomb
model_periods,model_mag,model_ph,model_fr,model_fi = modules.take_lomb(model_time,model_var,ofac,1./24)

obs_time = np.array(obs_time)
obs_var = np.array(obs_var)

#set plotting area & background to white
fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)

obs_periods,obs_mag, obs_breakpoint = modules.find_breakpoint(obs_periods,obs_mag)
model_periods,model_mag, model_breakpoint = modules.find_breakpoint(model_periods,model_mag)

plt.loglog(obs_periods,obs_mag, color='black', label = 'Obs.')
plt.axvline(x=obs_breakpoint, color = 'blue', linestyle = '--')
plt.loglog(model_periods,model_mag, color='red', label = 'GEOS %s'%(model_version))
plt.axvline(x=model_breakpoint, color = 'green', linestyle = '--')

def form2(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.2f' % x	

def form5(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.5f' % x

xformatter = FuncFormatter(form2)
yformatter = FuncFormatter(form5)

plt.grid(True)
leg=plt.legend(loc=4, prop={'size':24})
leg.get_frame().set_alpha(0.4)
ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
plt.xlabel('Period (Days)',fontsize = 21)
plt.ylabel('Amplitude (ppb)',fontsize = 21)
#plt.title(r'Lomb-Scargle Periodogram of Surface $O_3$ at %s, for Obs. & Model'%(site),fontsize=22)

ax.tick_params(axis='both', which='major', labelsize=18)
for tick in ax.get_xaxis().get_major_ticks():
    tick.set_pad(11.)
fig.tight_layout()
#plt.savefig('mace_head_spectra.png')
plt.show()
