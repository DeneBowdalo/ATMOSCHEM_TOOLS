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

invalid_sites = []

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-1]

print species

start_year = 1983
end_year = 2008
nyears = 5
vres = 'SURFACE'
timeres = 'H'

#obs_fname,species,start_year,end_year,vres,timeres = modules.get_obs_info(present_dir)

def run_LSP(model_data, x):

    print obs_refs[x]
    
    vals = model_data
    
    #check obs vals are valid
    valid = vals >= 0
    vals = vals[valid]
    model_time_val = model_time[valid]
    model_date_val = model_date[valid]

    full_times = modules.date_process(model_date,model_time,start_year)
    if timeres == 'M':
        full_times_year = full_times[:12]
    else:
        full_times_year = full_times[:8766]
    full_times_day = full_times[:24]

    valid_times = modules.date_process(model_date_val,model_time_val,start_year)  
      
    site_lon = obs_lons[x]

    #convert site_lon to 0 to 360 degs
    if site_lon < 0:
        site_lon = 360-np.abs(site_lon)
    
    #transform from UTC time to solar time 
    sun_time = lon_step_time*site_lon
    time_diff = sun_time - 0
    if time_diff > 12:
        time_diff = time_diff-24

    #make time start from 0    
    valid_times_from0 = modules.phase_start_correct(valid_times)

    periodic_periods = [1./4.,1./3.,1./2.,1.,365.25/4.,365.25/3.,365.25/2.,365.25]
    periods,mag,ph,fr,fi = modules.take_lomb_spec(valid_times_from0,vals,w=True,key_periods=periodic_periods)
    
    #get mean of values
    mean_array = np.average(vals)
    
    #correct all phases for start point (not actually being from 0 - just corrected to be)
    ph = modules.phase_start_point_correct_all(periodic_periods,ph,valid_times)

    key_diurnal_periods = [1./4.,1./3.,1./2.,1.]
    key_seasonal_periods = [365.25/4.,365.25/3.,365.25/2.,365.25]

    diurnal_mags = mag[:4]
    seasonal_mags = mag[4:]
    seasonal_phs = ph[4:]

    #get individual mags and phases
    daily_h3_mag = mag[0]
    daily_h2_mag = mag[1]
    daily_h1_mag = mag[2]
    orig_daily_mag = mag[3]
    daily_h3_ph = ph[0]
    daily_h2_ph = ph[1]
    daily_h1_ph = ph[2]
    orig_daily_ph = ph[3]
    
    seasonal_h3_mag = mag[4]
    seasonal_h2_mag = mag[5]
    seasonal_h1_mag = mag[6]
    annual_mag = mag[7]
    seasonal_h3_ph = ph[4]
    seasonal_h2_ph = ph[5]
    seasonal_h1_ph = ph[6]
    annual_ph = ph[7]

    #convert sub diurnal phases from UTC to solar time
    daily_h3_ph = modules.solar_time_phase_corrector(daily_h3_ph,6,time_diff)
    daily_h2_ph = modules.solar_time_phase_corrector(daily_h2_ph,24./3.,time_diff)
    daily_h1_ph = modules.solar_time_phase_corrector(daily_h1_ph,12,time_diff)
    orig_daily_ph = modules.solar_time_phase_corrector(orig_daily_ph,24,time_diff)
    diurnal_phs = [daily_h3_ph,daily_h2_ph,daily_h1_ph,orig_daily_ph]

    #convolve annual cycle and harmonics to seasonal waveform for 1 year
    seasonal_mag,seasonal_min_ph,seasonal_max_ph,seasonal_waveform,seasonal_ff = modules.period_convolution(key_seasonal_periods,full_times_year,seasonal_mags,seasonal_phs,mean_array)

    #convolve diurnal cycle and harmonics to diurnal waveform for 1 day
    diurnal_mag,diurnal_min_ph,diurnal_max_ph,diurnal_waveform,diurnal_ff = modules.period_convolution(key_diurnal_periods,full_times_day,diurnal_mags,diurnal_phs,mean_array)
    
    #convolve all 
    full_mag,full_min_ph,full_max_ph,full_waveform,full_ff = modules.period_convolution(periodic_periods,full_times,mag,ph,mean_array)

    #convert phase to time
    daily_h3_ph = modules.convert_phase_units_actual_single(daily_h3_ph,6.)
    daily_h2_ph = modules.convert_phase_units_actual_single(daily_h2_ph,24./3.)
    daily_h1_ph = modules.convert_phase_units_actual_single(daily_h1_ph,12.)
    orig_daily_ph = modules.convert_phase_units_actual_single(orig_daily_ph,24.)
    diurnal_min_ph = modules.convert_phase_units_actual_single(diurnal_min_ph,24.)
    diurnal_max_ph = modules.convert_phase_units_actual_single(diurnal_max_ph,24.)
    seasonal_h3_ph = modules.convert_phase_units_actual_single(seasonal_h3_ph,3.)
    seasonal_h2_ph = modules.convert_phase_units_actual_single(seasonal_h2_ph,4.)
    seasonal_h1_ph = modules.convert_phase_units_actual_single(seasonal_h1_ph,6.)
    annual_ph = modules.convert_phase_units_actual_single(annual_ph,12.)
    seasonal_min_ph = modules.convert_phase_units_actual_single(seasonal_min_ph,12.)
    seasonal_max_ph = modules.convert_phase_units_actual_single(seasonal_max_ph,12.)

    return (x,daily_h3_mag,daily_h3_ph,daily_h2_mag,daily_h2_ph,daily_h1_mag,daily_h1_ph,orig_daily_mag,orig_daily_ph,diurnal_mag,diurnal_min_ph,diurnal_max_ph,seasonal_h3_mag,seasonal_h3_ph,seasonal_h2_mag,seasonal_h2_ph,seasonal_h1_mag,seasonal_h1_ph,annual_mag,annual_ph,seasonal_mag,seasonal_min_ph,seasonal_max_ph,mean_array,diurnal_waveform,seasonal_waveform,full_waveform)

start_years = range(start_year,(end_year-nyears)+1)
end_years = range(start_year+nyears,end_year+1)

obs_refs = ['cpt','arh','mbi','mcm','nmy','syo','alt','ca0420g','dk0010g','no0042g','129401','hedo','ijira','ochiishi','oki','rishiri','sado-seki','tappi','yusuhara','dmv','ryo','tkb','yon','ice','abis0007a','abfr38011','abfr38012','150010006','780200001','vii423','ogasawara','bmw','cvo','mnm','rpb','smo','bhd','cgo','lau','sja']
root_grp = Dataset('/work/home/db876/observations/surface/%s/process/GLOBAL_SURFACE_%s_1970_2015_H_ALL.nc'%(species,species))

obs_lats = []
obs_lons = []

for ref in obs_refs:
    site_group = root_grp.groups[ref] 
    obs_lats = np.append(obs_lats,site_group.latitude)
    obs_lons = np.append(obs_lons,site_group.longitude)

for c in range(len(start_years)):

    sy = start_years[c]
    ey = end_years[c]

    #read in obs names
    model_root_grp = Dataset('/work/home/db876/plotting_tools/model_files/GFDLAM3_SURFACE_%s_%s_%s_*_*_*_H_nudged.nc'%(species,sy,ey))

    model_data = model_root_grp.variables[species.lower()][:]*1e9
    model_date = model_root_grp.variables['date'][:]
    model_time = model_root_grp.variables['time'][:]
    lat_c = model_root_grp.variables['lat_centre'][:]
    lon_c = model_root_grp.variables['lon_centre'][:]
    lat_e = model_root_grp.variables['lat_edges'][:]
    lon_e = model_root_grp.variables['lon_edges'][:]

    lat_indices = []
    lon_indices = []
    gridbox_num = []
    linear_data = []

    for obs_lat,obs_lon in zip(obs_lats, obs_lons):    
        lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    
        linear_data.append(model_data[:,lat_n,lon_n])
    
        lat_indices.append(lat_n)
        lon_indices.append(lon_n)

    # obs_refs = obs_refs[:5]
    # obs_lats = obs_lats[:5]
    # obs_lons = obs_lons[:5]
    # obs_alt = obs_alt[:5]
    # linear_data = linear_data[:5]

    if __name__ == '__main__':
        pool = multiprocessing.Pool(processes=16)
        results = [pool.apply_async(run_LSP, (linear_data[x],x)) for x in range(len(obs_refs))]
        big_array = [r.get() for r in results]

        pool.terminate()

    indices_array = []
    daily_h3_mag_array = []
    daily_h3_phase_array = []
    daily_h2_mag_array = []
    daily_h2_phase_array = []
    daily_h1_mag_array = []
    daily_h1_phase_array = []
    orig_daily_mag_array = []
    orig_daily_phase_array = []
    daily_mag_array = []
    daily_phase_min_array = []
    daily_phase_max_array = []
    seasonal_h3_mag_array = []
    seasonal_h3_phase_array = []
    seasonal_h2_mag_array = []
    seasonal_h2_phase_array = []
    seasonal_h1_mag_array = []
    seasonal_h1_phase_array = []
    annual_mag_array = []
    annual_phase_array = [] 
    seasonal_mag_array = []
    seasonal_phase_min_array = []
    seasonal_phase_max_array = []
    mean_array = []
    daily_waveform = []
    seasonal_waveform = []
    full_waveform = []

    for i in range(len(big_array)):
        cut = big_array[i]
    
        indices_array.append(cut[0])
        daily_h3_mag_array.append(cut[1])
        daily_h3_phase_array.append(cut[2])
        daily_h2_mag_array.append(cut[3])
        daily_h2_phase_array.append(cut[4])
        daily_h1_mag_array.append(cut[5])
        daily_h1_phase_array.append(cut[6])
        orig_daily_mag_array.append(cut[7])
        orig_daily_phase_array.append(cut[8])
        daily_mag_array.append(cut[9])
        daily_phase_min_array.append(cut[10])
        daily_phase_max_array.append(cut[11])
        seasonal_h3_mag_array.append(cut[12])
        seasonal_h3_phase_array.append(cut[13])
        seasonal_h2_mag_array.append(cut[14])
        seasonal_h2_phase_array.append(cut[15])
        seasonal_h1_mag_array.append(cut[16])
        seasonal_h1_phase_array.append(cut[17])
        annual_mag_array.append(cut[18])
        annual_phase_array.append(cut[19]) 
        seasonal_mag_array.append(cut[20])
        seasonal_phase_min_array.append(cut[21])
        seasonal_phase_max_array.append(cut[22])
        mean_array.append(cut[23])
        daily_waveform.append(cut[24])
        seasonal_waveform.append(cut[25])
        full_waveform.append(cut[26])

    daily_h3_mag_array = np.array(daily_h3_mag_array)
    daily_h3_phase_array = np.array(daily_h3_phase_array)
    daily_h2_mag_array = np.array(daily_h2_mag_array)
    daily_h2_phase_array = np.array(daily_h2_phase_array)
    daily_h1_mag_array = np.array(daily_h1_mag_array)
    daily_h1_phase_array = np.array(daily_h1_phase_array)
    orig_daily_mag_array = np.array(orig_daily_mag_array)
    orig_daily_phase_array = np.array(orig_daily_phase_array)
    daily_mag_array = np.array(daily_mag_array)
    daily_phase_min_array = np.array(daily_phase_min_array)  
    daily_phase_max_array = np.array(daily_phase_max_array)
    seasonal_h3_mag_array = np.array(seasonal_h3_mag_array)
    seasonal_h3_phase_array = np.array(seasonal_h3_phase_array)
    seasonal_h2_mag_array = np.array(seasonal_h2_mag_array)
    seasonal_h2_phase_array = np.array(seasonal_h2_phase_array)
    seasonal_h1_mag_array = np.array(seasonal_h1_mag_array)
    seasonal_h1_phase_array = np.array(seasonal_h1_phase_array)
    annual_mag_array = np.array(annual_mag_array)
    annual_phase_array = np.array(annual_phase_array)
    seasonal_mag_array = np.array(seasonal_mag_array) 
    seasonal_phase_min_array = np.array(seasonal_phase_min_array) 
    seasonal_phase_max_array = np.array(seasonal_phase_max_array)
    mean_array = np.array(mean_array)
    daily_waveform = np.array(daily_waveform)
    seasonal_waveform = np.array(seasonal_waveform)
    full_waveform = np.array(full_waveform)

    #sort arrays by indices array for sanity
    daily_h3_mag_array = daily_h3_mag_array[indices_array]
    daily_h3_phase_array = daily_h3_phase_array[indices_array]
    daily_h2_mag_array = daily_h2_mag_array[indices_array]
    daily_h2_phase_array = daily_h2_phase_array[indices_array]
    daily_h1_mag_array = daily_h1_mag_array[indices_array]
    daily_h1_phase_array = daily_h1_phase_array[indices_array]
    orig_daily_mag_array = orig_daily_mag_array[indices_array]
    orig_daily_phase_array = orig_daily_phase_array[indices_array]
    daily_mag_array = daily_mag_array[indices_array]
    daily_phase_min_array = daily_phase_min_array[indices_array]
    daily_phase_max_array = daily_phase_max_array[indices_array]
    seasonal_h3_mag_array = seasonal_h3_mag_array[indices_array]
    seasonal_h3_phase_array = seasonal_h3_phase_array[indices_array]
    seasonal_h2_mag_array = seasonal_h2_mag_array[indices_array]
    seasonal_h2_phase_array = seasonal_h2_phase_array[indices_array]
    seasonal_h1_mag_array = seasonal_h1_mag_array[indices_array]
    seasonal_h1_phase_array = seasonal_h1_phase_array[indices_array]
    annual_mag_array = annual_mag_array[indices_array]
    annual_phase_array = annual_phase_array[indices_array]
    seasonal_mag_array = seasonal_mag_array[indices_array]
    seasonal_phase_min_array = seasonal_phase_min_array[indices_array]
    seasonal_phase_max_array = seasonal_phase_max_array[indices_array]
    mean_array = mean_array[indices_array]
    daily_waveform = daily_waveform[indices_array]
    seasonal_waveform = seasonal_waveform[indices_array]
    full_waveform = full_waveform[indices_array]

    #obs_refs = obs_refs.tolist()
    
    #save out sig periods to netcdf
    root_grp_period = Dataset('GFDLAM3_nudged_sig_periods_%s_%s.nc'%(sy,ey), 'w')
    root_grp_period.description = 'Amplitudes and Phases for key periods from Lomb-Scargle Periodogram analysis of GFDLAM3 - Program written by Dene Bowdalo'
    root_grp_period.createDimension('diurnal', 24)
    root_grp_period.createDimension('seasonal', 8766)
    root_grp_period.createDimension('all', len(full_waveform[0]))

    for i in range(len(obs_refs)):
        # dimensions
        site_ref = obs_refs[i]
    
        ref_period = root_grp_period.createGroup('%s'%(site_ref.lower()))
        ref_period.daily_harmonic3_amplitude = daily_h3_mag_array[i]
        ref_period.daily_harmonic2_amplitude = daily_h2_mag_array[i]
        ref_period.daily_harmonic1_amplitude = daily_h1_mag_array[i]
        ref_period.original_daily_amplitude = orig_daily_mag_array[i]
        ref_period.daily_amplitude = daily_mag_array[i]
        ref_period.seasonal_harmonic3_amplitude = seasonal_h3_mag_array[i]
        ref_period.seasonal_harmonic2_amplitude = seasonal_h2_mag_array[i]
        ref_period.seasonal_harmonic1_amplitude = seasonal_h1_mag_array[i]
        ref_period.annual_amplitude = annual_mag_array[i]
        ref_period.seasonal_amplitude = seasonal_mag_array[i]
        ref_period.daily_harmonic3_phase = daily_h3_phase_array[i]
        ref_period.daily_harmonic2_phase = daily_h2_phase_array[i]
        ref_period.daily_harmonic1_phase = daily_h1_phase_array[i]
        ref_period.original_daily_phase = orig_daily_phase_array[i]
        ref_period.daily_phase = daily_phase_max_array[i]
        ref_period.daily_phase_min = daily_phase_min_array[i]
        ref_period.seasonal_harmonic3_phase = seasonal_h3_phase_array[i]
        ref_period.seasonal_harmonic2_phase = seasonal_h2_phase_array[i]
        ref_period.seasonal_harmonic1_phase = seasonal_h1_phase_array[i]
        ref_period.annual_phase = annual_phase_array[i]
        ref_period.seasonal_phase = seasonal_phase_max_array[i]
        ref_period.seasonal_phase_min = seasonal_phase_min_array[i]
        ref_period.average = mean_array[i]
        d_w = ref_period.createVariable('daily_waveform', 'f8', ('diurnal',))
        s_w = ref_period.createVariable('seasonal_waveform', 'f8', ('seasonal',))
        a_w = ref_period.createVariable('all_waveform', 'f8', ('all',))
        d_w[:] = daily_waveform[i]
        s_w[:] = seasonal_waveform[i]
        a_w[:] = full_waveform[i]
        
    root_grp_period.close()

end = datetime.datetime.now() 
diff = end - start
print 'time = ', diff.seconds

