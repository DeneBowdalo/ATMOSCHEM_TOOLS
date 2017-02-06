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
root_grp = Dataset('GLOBAL_SURFACE_O3_2006_2012.nc')
valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)
                   
site_ref = raw_input('Choose site from list. Sites with full set of yearly files between %i & %i are:\n%s\n'%(first_year,last_year+1,'   '.join(i for i in valid_refs)))

#read in specific site data
site_group = root_grp.groups[site_ref]

#read in variables for site
obs_var = site_group.variables['o3'][:]
obs_date = site_group.variables['date'][:]
obs_time = site_group.variables['time'][:]
obs_lat = site_group.latitude
obs_lon = site_group.longitude
obs_alt = site_group.altitude
obs_group = site_group.process_group

print obs_var
print len(obs_var)
print np.max(obs_var)

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
model_dict = {'4x5':'/work/home/db876/plotting_tools/binary_logs/v90103_4x5_GRID_O3.npy','2x2.5':'/work/home/db876/plotting_tools/binary_logs/v90103_2x2.5_GRID_O3.npy'}
model_version = raw_input('\nChoose Model Version.\n%s\n'%('   '.join(i for i in model_dict)))
model_f = model_dict[model_version]

model_data = read = np.load(model_f)
model_time = np.arange(0,2191,1./24)

#get model grid dims. for sim. type
lat_c,lat_e,lon_c,lon_e = modules.model_grids(model_version)
gridbox_count = len(lat_c)*len(lon_c)

#get model gridbox for obs site
gridbox_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)

model_var = model_data[gridbox_n::gridbox_count]
model_var = model_var*1e9

#----------------------------------------
model_time = model_time[valids]
model_var = model_var[valids]


#take lomb for obs and model
ofac = raw_input('\nChoose oversampling factor. (Typically 4.)\n')

#shift times to start from 0
obs_time = np.array(obs_time)
model_time = np.array(model_time)
obs_time_from0 = modules.phase_start_correct(obs_time)
model_time_from0 = modules.phase_start_correct(model_time)

#obs lomb
obs_periods,obs_mag,obs_ph,obs_fr,obs_fi,obs_amp_corr = modules.take_lomb(obs_time_from0,obs_var,ofac,1./24)

#model lomb
model_periods,model_mag,model_ph,model_fr,model_fi,model_amp_corr = modules.take_lomb(model_time_from0,model_var,ofac,1./24)

obs_time = np.array(obs_time)
obs_var = np.array(obs_var)

lon_step_time  = 24./360.

#convert site_lon to 0 to 360 degs
if obs_lon < 0:
    obs_lon = 360-np.abs(obs_lon)
    
#transform from UTC time to solar time 
sun_time = lon_step_time*obs_lon
print 'sun time = ', sun_time
time_diff = sun_time - 0
if time_diff > 12:
    time_diff = time_diff-24
print 'time diff =', time_diff


#correct magnitude and phase for spectral leakage 
zoomfact = 1000
daily_obs_mag,daily_obs_phase = modules.periodic_interp(obs_fr,obs_fi,zoomfact,obs_periods,1.,len(obs_var),obs_amp_corr)
daily_model_mag,daily_model_phase = modules.periodic_interp(model_fr,model_fi,zoomfact,model_periods,1.,len(model_var),model_amp_corr)
ha_obs_mag,ha_obs_phase = modules.periodic_interp(obs_fr,obs_fi,zoomfact,obs_periods,(365.25/2.),len(obs_var),obs_amp_corr)
ha_model_mag,ha_model_phase = modules.periodic_interp(model_fr,model_fi,zoomfact,model_periods,(365.25/2.),len(model_var),model_amp_corr)
annual_obs_mag,annual_obs_phase = modules.periodic_interp(obs_fr,obs_fi,zoomfact,obs_periods,365.25,len(obs_var),obs_amp_corr)
annual_model_mag,annual_model_phase = modules.periodic_interp(model_fr,model_fi,zoomfact,model_periods,365.25,len(model_var),model_amp_corr)

#correct for phase shift from sites where raw times do not start from 0
daily_obs_phase = modules.phase_start_point_correct(1.,daily_obs_phase,obs_time)
daily_model_phase = modules.phase_start_point_correct(1.,daily_model_phase,model_time)
ha_obs_phase = modules.phase_start_point_correct((365.25/2.),ha_obs_phase,obs_time)
ha_model_phase = modules.phase_start_point_correct((365.25/2.),ha_model_phase,model_time)
annual_obs_phase = modules.phase_start_point_correct(365.25,annual_obs_phase,obs_time)
annual_model_phase = modules.phase_start_point_correct(365.25,annual_model_phase,model_time)

#convert phase to time
daily_obs_phase = modules.convert_phase_units_actual_single(daily_obs_phase,24)
daily_model_phase = modules.convert_phase_units_actual_single(daily_model_phase,24)
ha_obs_phase = modules.convert_phase_units_actual_single(ha_obs_phase,6)
ha_model_phase = modules.convert_phase_units_actual_single(ha_model_phase,6)
annual_obs_phase = modules.convert_phase_units_actual_single(annual_obs_phase,12)
annual_model_phase = modules.convert_phase_units_actual_single(annual_model_phase,12)



#Maske daily phase SolaR time from UTC
daily_obs_phase = daily_obs_phase + time_diff
daily_model_phase = daily_model_phase + time_diff

if daily_obs_phase >= 24:
	daily_obs_phase = daily_obs_phase-24
if daily_obs_phase < 0:
	daily_obs_phase = 24-np.abs(daily_obs_phase)

if daily_model_phase >= 24:
	daily_model_phase = daily_model_phase-24
if daily_model_phase < 0:
	daily_model_phase = 24-np.abs(daily_model_phase)


print 'Obs Daily Amp = ', daily_obs_mag
print 'Model Daily Amp = ',daily_model_mag
print 'Obs Daily Phase = ', daily_obs_phase
print 'Model Daily Phase = ',daily_model_phase, '\n'
print 'Obs Half-Annual Amp = ', ha_obs_mag
print 'Model Half-Annual Amp = ',ha_model_mag
print 'Obs Half-Annual Phase = ', ha_obs_phase                                                                                                                                                                                                 
print 'Model Half-Annual Phase = ',ha_model_phase, '\n'
print 'Obs Annual Amp = ', annual_obs_mag
print 'Model Annual Amp = ',annual_model_mag
print 'Obs Annual Phase = ', annual_obs_phase                                                                                                                                                                                             
print 'Model Annual Phase = ',annual_model_phase




