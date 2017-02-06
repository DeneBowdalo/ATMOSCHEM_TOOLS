import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from cmath import *
from math import *
import csv
import datetime
import lomb_phase
from scipy import signal
import multiprocessing
import datetime
import time
import modules
from netCDF4 import Dataset
from scipy import stats
from scipy import interpolate

invalid_sites = []

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres,run_type = modules.get_obs_info(present_dir)
if run_type == 'model':
    model_fname,species,start_year,end_year,vres,timeres = modules.get_model_info(present_dir)

def run_LSP(vals,run_type,x):

    print obs_refs[x]
    
    #check obs vals are valid
    valid = vals >= 0
    vals = vals[valid]
    
    if run_type == 'obs':
        site_lon = obs_lons[x]
        valid_times = obs_ref_time[valid]
    
    if run_type == 'model':
        lat_i = lat_indices[x]
        lon_i = lon_indices[x]
        current_lat = lat_c[lat_i]
        current_lon = lon_c[lon_i]
        site_lon = lon_c[lon_i]
        valid_times = model_ref_time[valid]
    
    #make time start from 0    
    valid_times_from0 = modules.time_from0(valid_times)
    
    samp_step = 1./24
    f = interpolate.interp1d(valid_times_from0, vals)
    valid_times_from0 = np.arange(np.min(valid_times_from0),np.max(valid_times_from0),samp_step) 
    vals =  f(valid_times_from0)
      
    #convert site_lon to 0 to 360 degs
    if site_lon < 0:
        site_lon = 360-np.abs(site_lon)
    
    #transform factor for conversion from UTC time to solar time
    sun_time = lon_step_time*site_lon
    time_diff = sun_time - 0
    if time_diff > 12:
        time_diff = time_diff-24

    #take lsp
    ofac = 1
    periodic_periods = [1./4.,1./3.,1./2.,1.,365.25/4.,365.25/3.,365.25/2.,365.25]
    periods,mag,ph,fr,fi,amp_corr = modules.lomb_spectra(valid_times_from0,vals,ofac,samp_step,w=True,kp=periodic_periods) 
    
    lower_data_limit = 0.0
    upper_data_limit = 100.0
    bp_lower_limit = 5.0
    bp_upper_limit = 20.0
    
    #get info of met. regimes through model fit.
    grad1,grad2,bp1,line1_periods,line1_mag,line2_periods,line2_mag,ave1,ave2,med1,med2,sum1,sum2,line1_s,line1_e,line2_s,line2_e = modules.spectra_fit_free_piecewise(periods,mag,ofac,lower_data_limit,upper_data_limit,bp_lower_limit,bp_upper_limit)

    #correct all phases for start point (not actually being from 0 - just corrected to be)
    ph = modules.phase_start_time_relative_correct(periods,ph,valid_times)

    #convert phase to time(days)
    ph = modules.radians_to_time_spectra(ph,periods)
        
    return (x,periods,mag,ph,grad1,grad2,bp1,line1_periods,line1_mag,line2_periods,line2_mag,ave1,ave2,med1,med2,sum1,sum2,line1_s,line1_e,line2_s,line2_e)
        
obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_raw_class,obs_anthrome_class,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)

if run_type == 'obs':
    var = obs_std_var

if run_type == 'model':
    model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

    lat_indices = []
    lon_indices = []
    var = []

    for obs_lat,obs_lon in zip(obs_lats, obs_lons):    
        lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    
        var.append(model_std_var[:,lat_n,lon_n])
    
        lat_indices.append(lat_n)
        lon_indices.append(lon_n)

if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=16)
    results = [pool.apply_async(run_LSP, (var[x],run_type,x)) for x in range(len(obs_refs))]
    big_array = [r.get() for r in results]

    pool.terminate()
    
indices_array = []
periods_array = []
mag_array = []
phase_array = []
grad1_array = []
grad2_array = []
bp1_array = []
line1_periods = []
line1_mag = []
line2_periods = []
line2_mag = []
ave1_array = []
ave2_array = []
med1_array = []
med2_array = []
sum1_array = []
sum2_array = []
line1_s = []
line1_e = []
line2_s = []
line2_e = []


for i in range(len(big_array)):
    cut = big_array[i]
    
    indices_array.append(cut[0])
    periods_array.append(cut[1])
    mag_array.append(cut[2])
    phase_array.append(cut[3])
    grad1_array.append(cut[4])
    grad2_array.append(cut[5])
    bp1_array.append(cut[6])
    line1_periods.append(cut[7])
    line1_mag.append(cut[8])
    line2_periods.append(cut[9])
    line2_mag.append(cut[10])
    ave1_array.append(cut[11])
    ave2_array.append(cut[12])
    med1_array.append(cut[13])
    med2_array.append(cut[14])
    sum1_array.append(cut[15])
    sum2_array.append(cut[16])
    line1_s.append(cut[17])
    line1_e.append(cut[18])
    line2_s.append(cut[19])
    line2_e.append(cut[20])

    
periods_array = np.array(periods_array)
mag_array = np.array(mag_array)
phase_array = np.array(phase_array) 
grad1_array = np.array(grad1_array)
grad2_array = np.array(grad2_array)
bp1_array = np.array(bp1_array)
line1_periods = np.array(line1_periods)
line1_mag = np.array(line1_mag)
line2_periods = np.array(line2_periods)
line2_mag = np.array(line2_mag)
ave1_array = np.array(ave1_array)
ave2_array = np.array(ave2_array)
med1_array = np.array(med1_array)
med2_array = np.array(med2_array)
sum1_array = np.array(sum1_array)
sum2_array = np.array(sum2_array)
line1_s = np.array(line1_s)
line1_e = np.array(line1_e)
line2_s = np.array(line2_s)
line2_e = np.array(line2_e)

#sort arrays by indices array for sanity
periods_array = periods_array[indices_array]
mag_array = mag_array[indices_array]
phase_array = phase_array[indices_array]
grad1_array = grad1_array[indices_array]
grad2_array = grad2_array[indices_array]
bp1_array = bp1_array[indices_array]
line1_periods = line1_periods[indices_array]
line1_mag = line1_mag[indices_array]
line2_periods = line2_periods[indices_array]
line2_mag = line2_mag[indices_array]
ave1_array = ave1_array[indices_array]
ave2_array = ave2_array[indices_array]
med1_array = med1_array[indices_array]
med2_array = med2_array[indices_array]
sum1_array = sum1_array[indices_array]
sum2_array = sum2_array[indices_array]
line1_s = line1_s[indices_array]
line1_e = line1_e[indices_array]
line2_s = line2_s[indices_array]
line2_e = line2_e[indices_array]


#save out spectra to netcdf
root_grp_spec = Dataset('LSP_spectra.nc', 'w')
root_grp_spec.createDimension('site_n', None)

for i in range(len(obs_refs)):
    site_ref = obs_refs[i]
    
    ref = root_grp_spec.createGroup('%s'%(site_ref.lower()))
    amplitude = ref.createVariable('amplitude', 'f4', ('site_n',))
    phase = ref.createVariable('phase', 'f4', ('site_n',))
    period = ref.createVariable('period', 'f4', ('site_n',))
    w_period = ref.createVariable('weather_period', 'f4', ('site_n',))
    w_amp = ref.createVariable('weather_amplitude', 'f4', ('site_n',))
    mw_period = ref.createVariable('macroweather_period', 'f4', ('site_n',))
    mw_amp = ref.createVariable('macroweather_amplitude', 'f4', ('site_n',))
    
    ref.weather_gradient = grad1_array[i]
    ref.macroweather_gradient = grad2_array[i]
    ref.weather_macroweather_transition = bp1_array[i]
    ref.weather_median = med1_array[i]
    ref.macroweather_median = med2_array[i]
    ref.weather_average = ave1_array[i]
    ref.macroweather_average = ave2_array[i]
    ref.weather_sum = sum1_array[i]
    ref.macroweather_sum = sum2_array[i]
    ref.weather_start_amplitude = line1_s[i]
    ref.weather_end_amplitude = line1_e[i]
    ref.macroweather_start_amplitude = line2_s[i]
    ref.macroweather_end_amplitude = line2_e[i]
    
    amplitude[:] = mag_array[i]
    phase[:] = phase_array[i]
    period[:] = periods_array[i]
    w_period[:] = line1_periods[i]
    w_amp[:] = line1_mag[i]
    mw_period[:] = line2_periods[i]
    mw_amp[:] = line2_mag[i]
    
    
root_grp_spec.close()


end = datetime.datetime.now() 
diff = end - start
print 'time = ', diff.seconds




