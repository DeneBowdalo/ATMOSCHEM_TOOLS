import numpy as np
import matplotlib.pyplot as plt
import glob
import lomb_phase
import modules
import datetime
from netCDF4 import Dataset
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
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
obs_time = np.array(obs_time)
          
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
#obs_periods,obs_mag = modules.take_lomb(obs_time,obs_var,ofac)

#remove data randomly from model data for tests
num_99_9 = int((len(model_time)/100.)*99.9)
num_99 = int((len(model_time)/100.)*99)
num_95 = int((len(model_time)/100.)*95)
num_90 = int((len(model_time)/100.)*90)
num_80 = int((len(model_time)/100.)*80)
num_60 = int((len(model_time)/100.)*60)
num_40 = int((len(model_time)/100.)*40)
num_20 = int((len(model_time)/100.)*20)
num_5 = int((len(model_time)/100.)*5)
num_1 = int((len(model_time)/100.)*1)

inds_99_9= random.sample(range(len(model_time)), num_99_9)
inds_99= random.sample(range(len(model_time)), num_99)
inds_95= random.sample(range(len(model_time)), num_95)
inds_90= random.sample(range(len(model_time)), num_90)
inds_80= random.sample(range(len(model_time)), num_80)
inds_60= random.sample(range(len(model_time)), num_60)
inds_40= random.sample(range(len(model_time)), num_40)
inds_20= random.sample(range(len(model_time)), num_20)
inds_5= random.sample(range(len(model_time)), num_5)
inds_1= random.sample(range(len(model_time)), num_1)

inds_99_9 = sorted(inds_99_9)
inds_99 = sorted(inds_99)
inds_95 = sorted(inds_95)
inds_90 = sorted(inds_90)
inds_80 = sorted(inds_80)
inds_60 = sorted(inds_60)
inds_40 = sorted(inds_40)
inds_20 = sorted(inds_20)
inds_5 = sorted(inds_5)
inds_1 = sorted(inds_1)


##cut out 1 year
year_inds = np.arange(0,8736,1)
year_inds_2 = np.arange(17496,52560,1)

year_inds = np.append(year_inds,year_inds_2)

#cut out a year
year_inds_b = np.arange(0,8736,1)
year_inds_2_b = np.arange(17496,26280,1)
year_inds_3_b = np.arange(26280,35040,1)
year_inds_4_b = np.arange(43800,52560,1)

year_inds_b = np.append(year_inds_b,year_inds_2_b)
year_inds_b = np.append(year_inds_b,year_inds_3_b)
year_inds_b = np.append(year_inds_b,year_inds_4_b)

model_time_99_9 = model_time[inds_99_9]
model_time_99 = model_time[inds_99]
model_time_95 = model_time[inds_95]
model_time_90 = model_time[inds_90]
model_time_80 = model_time[inds_80]
model_time_60 = model_time[inds_60]
model_time_40 = model_time[inds_40]
model_time_20 = model_time[inds_20]
model_time_5 = model_time[inds_5]
model_time_1 = model_time[inds_1]
model_time_yg = model_time[year_inds]
model_time_2yg = model_time[year_inds_b]

model_var_99_9 = model_var[inds_99_9]
model_var_99 = model_var[inds_99]
model_var_95 = model_var[inds_95]
model_var_90 = model_var[inds_90]
model_var_80 = model_var[inds_80]
model_var_60 = model_var[inds_60]
model_var_40 = model_var[inds_40]
model_var_20 = model_var[inds_20]
model_var_5 = model_var[inds_5]
model_var_1 = model_var[inds_1]
model_var_yg = model_var[year_inds]
model_var_2yg = model_var[year_inds_b]

previous = -1
for i in inds_99:
    if i - previous == 2:
        print 'yes', i-1
    previous+=1

#model lomb
model_periods,model_mag= modules.take_lomb(model_time,model_var,ofac)
model_periods_99_9,model_mag_99_9 = modules.take_lomb(model_time_99_9,model_var_99_9,ofac)
model_periods_99,model_mag_99 = modules.take_lomb(model_time_99,model_var_99,ofac)
model_periods_95,model_mag_95 = modules.take_lomb(model_time_95,model_var_95,ofac)
model_periods_90,model_mag_90 = modules.take_lomb(model_time_90,model_var_90,ofac)
model_periods_80,model_mag_80 = modules.take_lomb(model_time_80,model_var_80,ofac)
model_periods_60,model_mag_60 = modules.take_lomb(model_time_60,model_var_60,ofac)
model_periods_40,model_mag_40 = modules.take_lomb(model_time_40,model_var_40,ofac)
model_periods_20,model_mag_20 = modules.take_lomb(model_time_20,model_var_20,ofac)
model_periods_5,model_mag_5 = modules.take_lomb(model_time_5,model_var_5,ofac)
model_periods_1,model_mag_1 = modules.take_lomb(model_time_1,model_var_1,ofac)
model_periods_yg,model_mag_yg = modules.take_lomb(model_time_yg,model_var_yg,ofac)
model_periods_2yg,model_mag_2yg = modules.take_lomb(model_time_2yg,model_var_2yg,ofac)

#set plotting area & background to white
fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)

#plt.loglog(obs_periods,obs_mag, color='black', alpha = 0.8, label = 'Obs.')
plt.loglog(model_periods,model_mag, color='blue',alpha=0.6, label = 'GEOS %s 100pc data'%(model_version))
plt.loglog(model_periods_99_9,model_mag_99_9, color='black', alpha=0.6, label = 'GEOS %s 99.9 pc data'%(model_version))
plt.loglog(model_periods_99,model_mag_99, color='green', alpha=0.6, label = 'GEOS %s 99 pc data'%(model_version))
plt.loglog(model_periods_95,model_mag_95, color='red', alpha=0.6, label = 'GEOS %s 95pc data'%(model_version))
plt.loglog(model_periods_90,model_mag_90, color='purple', alpha=0.6, label = 'GEOS %s 90pc data'%(model_version))
plt.loglog(model_periods_80,model_mag_80, color='orange', alpha=0.6, label = 'GEOS %s 80pc data'%(model_version))
plt.loglog(model_periods_60,model_mag_60, color='cyan', alpha=0.6, label = 'GEOS %s 60pc data'%(model_version))
plt.loglog(model_periods_40,model_mag_40, color='blue', alpha=0.6,linestyle='--', label = 'GEOS %s 40pc data'%(model_version))
plt.loglog(model_periods_20,model_mag_20, color='black', alpha=0.6,linestyle='--', label = 'GEOS %s 20pc data'%(model_version))
plt.loglog(model_periods_5,model_mag_5, color='green', alpha=0.6,linestyle='--', label = 'GEOS %s 5pc data'%(model_version))
plt.loglog(model_periods_1,model_mag_1, color='red', alpha=0.6,linestyle='--', label = 'GEOS %s 1pc data'%(model_version))
plt.loglog(model_periods_yg,model_mag_yg, color='purple', alpha=0.6,linestyle='--', label = 'GEOS %s Year Gap data'%(model_version))
plt.loglog(model_periods_2yg,model_mag_2yg, color='orange', alpha=0.6,linestyle='--', label = 'GEOS %s 2 Year Gap data'%(model_version))

def form2(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.2f' % x	

def form6(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.6f' % x

xformatter = FuncFormatter(form6)
yformatter = FuncFormatter(form6)

#plt.xlim(0.99,1.01)
plt.grid(True)
leg=plt.legend(loc=4, prop={'size':21})
leg.get_frame().set_alpha(0.4)
ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
plt.xlabel('Period (Days)',fontsize = 21)
plt.ylabel('Magnitude (ppb)',fontsize = 21)
plt.title(r'Lomb-Scargle Periodogram of Surface $O_3$  at %s, for Obs. & Model'%(site),fontsize=22)

#p = [0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6,1e-10,1e-11,1e-12]
#for k in range(len(p)):
#    z = -np.log( 1 - (1-p[k])**(1/model_i_freq)) 
#    plt.axhline(y=z, label = p)
    
plt.show()
