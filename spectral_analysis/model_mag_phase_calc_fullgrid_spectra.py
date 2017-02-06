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

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

#get species from current directory
present_dir = os.getcwd()
model_fname,species,start_year,end_year,vres,timeres = modules.get_model_info(present_dir)

def run_LSP(vals,x):

    lat_i = lat_indices[x]
    lon_i = lon_indices[x]

    print lat_i,lon_i

    current_lat = lat_c[lat_i]
    current_lon = lon_c[lon_i]
    site_lon = lon_c[lon_i]
 
    #check model vals are valid
    valid = vals >= 0
    vals = vals[valid]
    valid_times = model_ref_time[valid]
    
    #convert site_lon to 0 to 360 degs
    if site_lon < 0:
        site_lon = 360-np.abs(site_lon)
    
    #transform factor for conversion from UTC time to solar time
    sun_time = lon_step_time*site_lon
    time_diff = sun_time - 0
    if time_diff > 12:
        time_diff = time_diff-24
             
    #make time start from 0    
    valid_times_from0 = modules.time_from0(valid_times)

    ofac = 1
    samp_step = 1./24
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
 
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

lat_i = 0
lon_i = 0

linear_data = []
for siten in range(n_boxes):
    linear_data.append(model_std_var[:,lat_i,lon_i])
   
    lat_indices.append(lat_i)
    lon_indices.append(lon_i)
   
    if lon_i == (len(lon_c)-1):
        lat_i+=1
        lon_i=0
    else:
        lon_i+=1

if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=16)
    results = [pool.apply_async(run_LSP, (linear_data[x],x)) for x in range(gridbox_count)]
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

grad1_array = np.reshape(grad1_array,(len(lat_c),len(lon_c)))
grad2_array = np.reshape(grad2_array,(len(lat_c),len(lon_c)))
bp1_array = np.reshape(bp1_array,(len(lat_c),len(lon_c)))
ave1_array = np.reshape(ave1_array,(len(lat_c),len(lon_c)))
ave2_array = np.reshape(ave2_array,(len(lat_c),len(lon_c)))
med1_array = np.reshape(med1_array,(len(lat_c),len(lon_c)))
med2_array = np.reshape(med2_array,(len(lat_c),len(lon_c)))
sum1_array = np.reshape(sum1_array,(len(lat_c),len(lon_c)))
sum2_array = np.reshape(sum2_array,(len(lat_c),len(lon_c)))
line1_s = np.reshape(line1_s,(len(lat_c),len(lon_c)))
line1_e = np.reshape(line1_e,(len(lat_c),len(lon_c)))
line2_s = np.reshape(line2_s,(len(lat_c),len(lon_c)))
line2_e = np.reshape(line2_e,(len(lat_c),len(lon_c)))

line1_periods = np.reshape(line1_periods,(-1,len(lat_c),len(lon_c)))
line2_periods = np.reshape(line2_periods,(-1,len(lat_c),len(lon_c)))
line1_mag = np.reshape(line1_mag,(-1,len(lat_c),len(lon_c)))
line2_mag = np.reshape(line2_mag,(-1,len(lat_c),len(lon_c)))

spectra_len = len(periods_array[0])
period_rs = np.empty((spectra_len,len(lat_c),len(lon_c)))
mag_rs = np.empty((spectra_len,len(lat_c),len(lon_c)))
phase_rs = np.empty((spectra_len,len(lat_c),len(lon_c)))

lat_i = 0
lon_i = 0
for i in range(len(periods_array)):
    current_period = periods_array[i]
    current_mag = mag_array[i]
    current_phase = phase_array[i]
    
    period_rs[:,lat_i,lon_i] = current_period
    mag_rs[:,lat_i,lon_i] = current_mag
    phase_rs[:,lat_i,lon_i] = current_phase
    
    lon_i+=1
    
    if lon_i == len(lon_c):
        lat_i+=1
        lon_i = 0
    
periods_array = np.copy(period_rs)
mag_array = np.copy(mag_rs)
phase_array = np.copy(phase_rs)
    
#save out sig periods to netcdf
root_grp_period = Dataset('LSP_spectra.nc', 'w')

root_grp_period.createDimension('lat_centre', len(lat_c))
root_grp_period.createDimension('lon_centre', len(lon_c))
root_grp_period.createDimension('lat_edges', len(lat_e))
root_grp_period.createDimension('lon_edges', len(lon_e))
root_grp_period.createDimension('spectra_len', None)

period = root_grp_period.createVariable('period', 'f4', ('spectra_len','lat_centre','lon_centre'))      
amp = root_grp_period.createVariable('amplitude', 'f4', ('spectra_len','lat_centre','lon_centre'))
ph = root_grp_period.createVariable('phase', 'f4', ('spectra_len','lat_centre','lon_centre'))   
weather_period = root_grp_period.createVariable('weather_period', 'f4', ('spectra_len','lat_centre','lon_centre'))
weather_amp = root_grp_period.createVariable('weather_amplitude', 'f4', ('spectra_len','lat_centre','lon_centre'))
mw_period = root_grp_period.createVariable('macroweather_period', 'f4', ('spectra_len','lat_centre','lon_centre'))
mw_amp = root_grp_period.createVariable('macroweather_amplitude', 'f4', ('spectra_len','lat_centre','lon_centre'))
   
g1 = root_grp_period.createVariable('weather_gradient', 'f4', ('lat_centre','lon_centre'))
g2 = root_grp_period.createVariable('macroweather_gradient', 'f4', ('lat_centre','lon_centre'))
b1 = root_grp_period.createVariable('weather_macroweather_transition', 'f4', ('lat_centre','lon_centre'))
weather_ave = root_grp_period.createVariable('weather_average', 'f4', ('lat_centre','lon_centre'))
macroweather_ave = root_grp_period.createVariable('macroweather_average', 'f4', ('lat_centre','lon_centre'))
weather_med = root_grp_period.createVariable('weather_median', 'f4', ('lat_centre','lon_centre'))
macroweather_med = root_grp_period.createVariable('macroweather_median', 'f4', ('lat_centre','lon_centre'))
weather_sum = root_grp_period.createVariable('weather_sum', 'f4', ('lat_centre','lon_centre'))
macroweather_sum = root_grp_period.createVariable('macroweather_sum', 'f4', ('lat_centre','lon_centre'))
weather_s_a = root_grp_period.createVariable('weather_start_amplitude', 'f4', ('lat_centre','lon_centre'))
macroweather_s_a = root_grp_period.createVariable('macroweather_start_amplitude', 'f4', ('lat_centre','lon_centre'))
weather_e_a = root_grp_period.createVariable('weather_end_amplitude', 'f4', ('lat_centre','lon_centre'))
macroweather_e_a = root_grp_period.createVariable('macroweather_end_amplitude', 'f4', ('lat_centre','lon_centre'))

lat_centre = root_grp_period.createVariable('lat_centre', 'f4', ('lat_centre',))
lon_centre = root_grp_period.createVariable('lon_centre', 'f4', ('lon_centre',))
lat_edge = root_grp_period.createVariable('lat_edges', 'f4', ('lat_edges',))
lon_edge = root_grp_period.createVariable('lon_edges', 'f4', ('lon_edges',))

period[:] = periods_array
weather_period[:] = line1_periods
mw_period[:] = line2_periods
amp[:] = mag_array
weather_amp[:] = line1_mag
mw_amp[:] = line2_mag
ph[:] = phase_array
g1[:] = grad1_array
g2[:] = grad2_array
b1[:] = bp1_array
weather_ave[:] = ave1_array
macroweather_ave[:] = ave2_array
weather_med[:] =  med1_array
macroweather_med[:] = med2_array
weather_sum[:] = sum1_array
macroweather_sum[:] = sum2_array
weather_s_a[:] = line1_s
macroweather_s_a[:] = line2_s
weather_e_a[:] = line1_e
macroweather_e_a[:] = line2_e

lat_centre[:] = lat_c
lon_centre[:] = lon_c
lat_edge[:] = lat_e     
lon_edge[:] = lon_e          
        
root_grp_period.close()

end = datetime.datetime.now() 
diff = end - start
print 'time = ', diff.seconds
