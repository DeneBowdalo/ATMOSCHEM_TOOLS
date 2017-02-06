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
full_obs_len = len(obs_var)
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

#calculate model gap indices from obs.
if model_version == '4x5':
    valids = valids[24:]
if full_obs_len < len(model_time):
    print "Can't estimate model gaps"
print len(valids)
print len(model_time)
model_time_gap = model_time[valids]
model_var_gap = model_var[valids]


#----------------------------------------

#align obs and model periods
obs_time,model_time,obs_var,model_var = modules.align_periods(obs_time,model_time,obs_var,model_var)

#take lomb for obs and model
ofac = raw_input('\nChoose oversampling factor. (Typically 4.)\n')

#obs lomb
obs_periods,obs_mag,obs_ph,obs_sig,obs_i_freq = modules.take_lomb(obs_time,obs_var,ofac)

#model lomb
model_periods,model_mag,model_ph,model_sig,model_i_freq = modules.take_lomb(model_time,model_var,ofac)

#calculate high freq. noise due to gaps - apply obs. gaps to model array 
model_periods_gap,model_mag_gap,model_ph_gap,model_sig_gap,model_i_freq_gap = modules.take_lomb(model_time_gap,model_var_gap,ofac)

#calculate lomb with instrumental error for O3 implied on model sims which are gapped
err_model_var = []
for i in range(len(model_var_gap)):
    #value of 0.0166 uncertainty for each point determined by standard error of (1 ppbv / sqrt(60*60))= 1/60. = 0.0166
    err_model_var = np.append(err_model_var,model_var_gap[i]+random.normalvariate(0,0.0166))
model_periods_err,model_mag_err,model_ph_err,model_sig_err,model_i_freq_err = modules.take_lomb(model_time_gap,err_model_var,ofac)

#remove significant peaks from spectrums


#do log bin smoothing on all lombs
smoothed_obs_periods, smoothed_obs_mag = modules.log_bin_smooth(len(obs_periods),20,obs_periods,obs_mag)
smoothed_model_periods, smoothed_model_mag = modules.log_bin_smooth(len(model_periods),20,model_periods,model_mag)
smoothed_model_periods_gap, smoothed_model_mag_gap = modules.log_bin_smooth(len(model_periods_gap),20,model_periods_gap,model_mag_gap)
smoothed_model_periods_err, smoothed_model_mag_err = modules.log_bin_smooth(len(model_periods_err),20,model_periods_gap,model_mag_err)

#calculate offset in noise at high frequency between model and model with gaps. Shows noise due to gaps
gap_noise = smoothed_model_mag_gap[-1]-smoothed_model_mag[-1]
print 'gap noise = ', gap_noise

#calculate offset in noise due to instrument error
instrument_noise = smoothed_model_mag_err[-1] - (smoothed_model_mag_gap[-1]-smoothed_model_mag[-1])
print 'instrument noise = ', instrument_noise


print 'total obs. offset = ', gap_noise+instrument_noise
#calculate difference between model and obs. at high frequencies, after removing gap & instrument noise
high_freq_diff = smoothed_obs_mag[-1] - smoothed_model_mag[-1]
print 'obs. - model noise = ', high_freq_diff
#high_freq_diff = high_freq_diff - gap_noise 
#high_freq_diff = high_freq_diff - instrument_noise
#print 'obs. - model noise = ', high_freq_diff
print 'obs. offset percent of freq. diff. = ', ((gap_noise+instrument_noise)/high_freq_diff)*100


#model_sig = sorted(model_sig, reverse=True)
#model_sig = np.array(model_sig)
#test = model_sig == 100

#cut_periods = model_periods[test]
#cut_mag = model_mag[test]

#set plotting area & background to white
fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)

plt.loglog(smoothed_obs_periods, smoothed_obs_mag, color = 'green', marker = 'x', alpha = 1,markersize=2,label = 'Obs')
plt.loglog(smoothed_model_periods, smoothed_model_mag, color = 'black', marker = 'x', alpha = 1,markersize=2,label = 'Model')
plt.loglog(smoothed_model_periods_gap, smoothed_model_mag_gap, color = 'red', marker = 'x', alpha = 1,markersize=2,label = 'Gap Model')
plt.loglog(smoothed_model_periods_err, smoothed_model_mag_err, color = 'blue', marker = 'x', alpha = 1,markersize=2,label = 'Error Model')
#plt.loglog(obs_periods,obs_mag, color='black', label = 'Obs.')
#plt.loglog(model_periods,model_mag, color='red', alpha=0.5, label = 'GEOS %s'%(model_version))

def form2(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.2f' % x	

def form6(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.6f' % x

#xformatter = FuncFormatter(form2)
#yformatter = FuncFormatter(form6)

#plt.grid(True)
leg=plt.legend(loc=0, prop={'size':21})
leg.get_frame().set_alpha(0.4)
#ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
#plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
#plt.xlabel('Period (Days)',fontsize = 21)
#plt.ylabel('Magnitude (ppb)',fontsize = 21)
#plt.title(r'Lomb-Scargle Periodogram of Surface $O_3$  at %s, for Obs. & Model'%(site),fontsize=22)

p = [0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6,1e-10,1e-11,1e-12]

#baliunas method
#for k in range(len(p)):
#    z = -np.log( 1 - (1-p[k])**(1/model_i_freq)) 
#    plt.axhline(y=z, label = p)

#hocke method
#for i in range(len(p)):
#    z=-np.log(1-(1.-p[i])**(1./len(model_time)))
#   plt.axhline(y=z, label = p)

plt.show()
