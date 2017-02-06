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
paths = present_dir.split('/')
species = paths[-2]
timeres = 'H'

#get all obs files
files = glob.glob('*HN.nc')


def run_LSP(vals, x):

    print obs_refs[x]
    
    #check obs vals are valid
    valid = vals >= 0
    vals = vals[valid]
    valid_times = obs_ref_time[valid]
    valid_datetimes = obs_datetime_time[valid]

    if timeres == 'H':
        full_times_year = obs_ref_time[:8766]
    elif timeres == 'D':
        full_times_year = obs_ref_time[:365]
    elif timeres == 'M':
        full_times_year = obs_ref_time[:12]
    full_times_day = obs_ref_time[:24] 
      
    site_lon = obs_lons[x]
    site_lat = obs_lats[x]

    #convert site_lon to 0 to 360 degs
    if site_lon < 0:
        site_lon = 360-np.abs(site_lon)
    
    #transform from UTC time to solar time 
    sun_time = lon_step_time*site_lon
    time_diff = sun_time - 0
    if time_diff > 12:
        time_diff = time_diff-24
        
    #cut vals into seasons
    valid_times_spring,vals_spring,valid_times_summer,vals_summer,valid_times_autumn,vals_autumn,valid_times_winter,vals_winter = modules.cut_season(valid_datetimes,valid_times,vals)
    
    #cut vals into day/night
    valid_times_day,vals_day,valid_times_night,vals_night = modules.cut_daynight(valid_datetimes,valid_times,vals,site_lat,site_lon,timeres)
    
    
    #make time start from 0    
    valid_times_from0 = modules.time_from0(valid_times)
    valid_times_from0_spring = modules.time_from0(valid_times_spring)
    valid_times_from0_summer = modules.time_from0(valid_times_summer)
    valid_times_from0_autumn = modules.time_from0(valid_times_autumn)
    valid_times_from0_winter = modules.time_from0(valid_times_winter)
    valid_times_from0_day = modules.time_from0(valid_times_day)
    valid_times_from0_night = modules.time_from0(valid_times_night)
    
    key_diurnal_periods = [1./4.,1./3.,1./2.,1.]
    key_seasonal_periods = [365.25/4.,365.25/3.,365.25/2.,365.25]
    periodic_periods = [1./4.,1./3.,1./2.,1.,365.25/4.,365.25/3.,365.25/2.,365.25]
    periods,mag,ph,fr,fi = modules.take_lomb_spec(valid_times_from0,vals,w=True,key_periods=periodic_periods)
    periods_spring,mag_spring,ph_spring,fr_spring,fi_spring = modules.take_lomb_spec(valid_times_from0_spring,vals_spring,w=True,key_periods=key_diurnal_periods)
    periods_summer,mag_summer,ph_summer,fr_summer,fi_summer = modules.take_lomb_spec(valid_times_from0_summer,vals_summer,w=True,key_periods=key_diurnal_periods)
    periods_autumn,mag_autumn,ph_autumn,fr_autumn,fi_autumn = modules.take_lomb_spec(valid_times_from0_autumn,vals_autumn,w=True,key_periods=key_diurnal_periods)
    periods_winter,mag_winter,ph_winter,fr_winter,fi_winter = modules.take_lomb_spec(valid_times_from0_winter,vals_winter,w=True,key_periods=key_diurnal_periods)
    periods_day,mag_day,ph_day,fr_day,fi_day = modules.take_lomb_spec(valid_times_from0_day,vals_day,w=True,key_periods=key_seasonal_periods)
    periods_night,mag_night,ph_night,fr_night,fi_night = modules.take_lomb_spec(valid_times_from0_night,vals_night,w=True,key_periods=key_seasonal_periods)
    
    
    #get mean of values
    mean = np.average(vals)
    mean_spring = np.average(vals_spring)
    mean_summer = np.average(vals_summer)
    mean_autumn = np.average(vals_autumn)
    mean_winter = np.average(vals_winter)
    mean_day = np.average(vals_day)
    mean_night = np.average(vals_night)
    
    #correct all phases for start point (not actually being from 0 - just corrected to be)
    ph = modules.phase_start_time_relative_correct(periodic_periods,ph,valid_times)
    ph_spring = modules.phase_start_time_relative_correct(key_diurnal_periods,ph_spring,valid_times_spring)
    ph_summer = modules.phase_start_time_relative_correct(key_diurnal_periods,ph_summer,valid_times_summer)
    ph_autumn = modules.phase_start_time_relative_correct(key_diurnal_periods,ph_autumn,valid_times_autumn)
    ph_winter = modules.phase_start_time_relative_correct(key_diurnal_periods,ph_winter,valid_times_winter)
    ph_day = modules.phase_start_time_relative_correct(key_seasonal_periods,ph_day,valid_times_day)
    ph_night = modules.phase_start_time_relative_correct(key_seasonal_periods,ph_night,valid_times_night)
    
    diurnal_mags = mag[:4]
    diurnal_mags_spring = mag_spring[:]
    diurnal_mags_summer = mag_summer[:]
    diurnal_mags_autumn = mag_autumn[:]
    diurnal_mags_winter = mag_winter[:]
    seasonal_mags = mag[4:]
    seasonal_mags_day = mag_day[:]
    seasonal_mags_night = mag_night[:]
    seasonal_phs = ph[4:]
    seasonal_phs_day = ph_day[:]
    seasonal_phs_night = ph_night[:]

    #get individual mags and phases
    daily_h3_mag = mag[0]
    daily_h2_mag = mag[1]
    daily_h1_mag = mag[2]
    orig_daily_mag = mag[3]
    daily_h3_ph = ph[0]
    daily_h2_ph = ph[1]
    daily_h1_ph = ph[2]
    orig_daily_ph = ph[3]
    daily_h3_mag_spring = mag_spring[0]
    daily_h2_mag_spring = mag_spring[1]
    daily_h1_mag_spring = mag_spring[2]
    orig_daily_mag_spring = mag_spring[3]
    daily_h3_ph_spring = ph_spring[0]
    daily_h2_ph_spring = ph_spring[1]
    daily_h1_ph_spring = ph_spring[2]
    orig_daily_ph_spring = ph_spring[3]
    daily_h3_mag_summer = mag_summer[0]
    daily_h2_mag_summer = mag_summer[1]
    daily_h1_mag_summer = mag_summer[2]
    orig_daily_mag_summer = mag_summer[3]
    daily_h3_ph_summer = ph_summer[0]
    daily_h2_ph_summer = ph_summer[1]
    daily_h1_ph_summer = ph_summer[2]
    orig_daily_ph_summer = ph_summer[3]
    daily_h3_mag_autumn = mag_autumn[0]
    daily_h2_mag_autumn = mag_autumn[1]
    daily_h1_mag_autumn = mag_autumn[2]
    orig_daily_mag_autumn = mag_autumn[3]
    daily_h3_ph_autumn = ph_autumn[0]
    daily_h2_ph_autumn = ph_autumn[1]
    daily_h1_ph_autumn = ph_autumn[2]
    orig_daily_ph_autumn = ph_autumn[3]
    daily_h3_mag_winter = mag_winter[0]
    daily_h2_mag_winter = mag_winter[1]
    daily_h1_mag_winter = mag_winter[2]
    orig_daily_mag_winter = mag_winter[3]
    daily_h3_ph_winter = ph_winter[0]
    daily_h2_ph_winter = ph_winter[1]
    daily_h1_ph_winter = ph_winter[2]
    orig_daily_ph_winter = ph_winter[3]
    seasonal_h3_mag = mag[4]
    seasonal_h2_mag = mag[5]
    seasonal_h1_mag = mag[6]
    annual_mag = mag[7]
    seasonal_h3_ph = ph[4]
    seasonal_h2_ph = ph[5]
    seasonal_h1_ph = ph[6]
    annual_ph = ph[7]
    seasonal_h3_mag_day = mag_day[0]
    seasonal_h2_mag_day = mag_day[1]
    seasonal_h1_mag_day = mag_day[2]
    annual_mag_day = mag_day[3]
    seasonal_h3_ph_day = ph_day[0]
    seasonal_h2_ph_day = ph_day[1]
    seasonal_h1_ph_day = ph_day[2]
    annual_ph_day = ph_day[3]
    seasonal_h3_mag_night = mag_night[0]
    seasonal_h2_mag_night = mag_night[1]
    seasonal_h1_mag_night = mag_night[2]
    annual_mag_night = mag_night[3]
    seasonal_h3_ph_night = ph_night[0]
    seasonal_h2_ph_night = ph_night[1]
    seasonal_h1_ph_night = ph_night[2]
    annual_ph_night = ph_night[3]

    #convert sub diurnal phases from UTC to solar time
    daily_h3_ph = modules.solar_time_phase_corrector(daily_h3_ph,6,time_diff)
    daily_h2_ph = modules.solar_time_phase_corrector(daily_h2_ph,24./3.,time_diff)
    daily_h1_ph = modules.solar_time_phase_corrector(daily_h1_ph,12,time_diff)
    orig_daily_ph = modules.solar_time_phase_corrector(orig_daily_ph,24,time_diff)
    diurnal_phs = [daily_h3_ph,daily_h2_ph,daily_h1_ph,orig_daily_ph]
    daily_h3_ph_spring = modules.solar_time_phase_corrector(daily_h3_ph_spring,6,time_diff)
    daily_h2_ph_spring = modules.solar_time_phase_corrector(daily_h2_ph_spring,24./3.,time_diff)
    daily_h1_ph_spring = modules.solar_time_phase_corrector(daily_h1_ph_spring,12,time_diff)
    orig_daily_ph_spring = modules.solar_time_phase_corrector(orig_daily_ph_spring,24,time_diff)
    diurnal_phs_spring = [daily_h3_ph_spring,daily_h2_ph_spring,daily_h1_ph_spring,orig_daily_ph_spring]
    daily_h3_ph_summer = modules.solar_time_phase_corrector(daily_h3_ph_summer,6,time_diff)
    daily_h2_ph_summer = modules.solar_time_phase_corrector(daily_h2_ph_summer,24./3.,time_diff)
    daily_h1_ph_summer = modules.solar_time_phase_corrector(daily_h1_ph_summer,12,time_diff)
    orig_daily_ph_summer = modules.solar_time_phase_corrector(orig_daily_ph_summer,24,time_diff)
    diurnal_phs_summer = [daily_h3_ph_summer,daily_h2_ph_summer,daily_h1_ph_summer,orig_daily_ph_summer]
    daily_h3_ph_autumn = modules.solar_time_phase_corrector(daily_h3_ph_autumn,6,time_diff)
    daily_h2_ph_autumn = modules.solar_time_phase_corrector(daily_h2_ph_autumn,24./3.,time_diff)
    daily_h1_ph_autumn = modules.solar_time_phase_corrector(daily_h1_ph_autumn,12,time_diff)
    orig_daily_ph_autumn = modules.solar_time_phase_corrector(orig_daily_ph_autumn,24,time_diff)
    diurnal_phs_autumn = [daily_h3_ph_autumn,daily_h2_ph_autumn,daily_h1_ph_autumn,orig_daily_ph_autumn]
    daily_h3_ph_winter = modules.solar_time_phase_corrector(daily_h3_ph_winter,6,time_diff)
    daily_h2_ph_winter = modules.solar_time_phase_corrector(daily_h2_ph_winter,24./3.,time_diff)
    daily_h1_ph_winter = modules.solar_time_phase_corrector(daily_h1_ph_winter,12,time_diff)
    orig_daily_ph_winter = modules.solar_time_phase_corrector(orig_daily_ph_winter,24,time_diff)
    diurnal_phs_winter = [daily_h3_ph_winter,daily_h2_ph_winter,daily_h1_ph_winter,orig_daily_ph_winter]

    #convolve annual cycle and harmonics to seasonal waveform for 1 year
    seasonal_mag,seasonal_min_ph,seasonal_max_ph,seasonal_waveform,seasonal_ff = modules.period_convolution(key_seasonal_periods,full_times_year,seasonal_mags,seasonal_phs,mean)
    seasonal_mag_day,seasonal_min_ph_day,seasonal_max_ph_day,seasonal_waveform_day,seasonal_ff_day = modules.period_convolution(key_seasonal_periods,full_times_year,seasonal_mags_day,seasonal_phs_day,mean_day)
    seasonal_mag_night,seasonal_min_ph_night,seasonal_max_ph_night,seasonal_waveform_night,seasonal_ff_night = modules.period_convolution(key_seasonal_periods,full_times_year,seasonal_mags_night,seasonal_phs_night,mean_night)
    na,na,na,seasonal_waveform_extended,na = modules.period_convolution(key_seasonal_periods,obs_ref_time,seasonal_mags,seasonal_phs,mean)
    na,na,na,seasonal_waveform_day_extended,na = modules.period_convolution(key_seasonal_periods,obs_ref_time,seasonal_mags_day,seasonal_phs_day,mean_day)
    na,na,na,seasonal_waveform_night_extended,na = modules.period_convolution(key_seasonal_periods,obs_ref_time,seasonal_mags_night,seasonal_phs_night,mean_night)

    #convolve diurnal cycle and harmonics to diurnal waveform for 1 day
    diurnal_mag,diurnal_min_ph,diurnal_max_ph,diurnal_waveform,diurnal_ff = modules.period_convolution(key_diurnal_periods,full_times_day,diurnal_mags,diurnal_phs,mean)
    diurnal_mag_spring,diurnal_min_ph_spring,diurnal_max_ph_spring,diurnal_waveform_spring,diurnal_ff_spring = modules.period_convolution(key_diurnal_periods,full_times_day,diurnal_mags_spring,diurnal_phs_spring,mean_spring)
    diurnal_mag_summer,diurnal_min_ph_summer,diurnal_max_ph_summer,diurnal_waveform_summer,diurnal_ff_summer = modules.period_convolution(key_diurnal_periods,full_times_day,diurnal_mags_summer,diurnal_phs_summer,mean_summer)
    diurnal_mag_autumn,diurnal_min_ph_autumn,diurnal_max_ph_autumn,diurnal_waveform_autumn,diurnal_ff_autumn = modules.period_convolution(key_diurnal_periods,full_times_day,diurnal_mags_autumn,diurnal_phs_autumn,mean_autumn)
    diurnal_mag_winter,diurnal_min_ph_winter,diurnal_max_ph_winter,diurnal_waveform_winter,diurnal_ff_winter = modules.period_convolution(key_diurnal_periods,full_times_day,diurnal_mags_winter,diurnal_phs_winter,mean_winter)
    na,na,na,diurnal_waveform_extended,na = modules.period_convolution(key_diurnal_periods,obs_ref_time,diurnal_mags,diurnal_phs,mean)
    na,na,na,diurnal_waveform_spring_extended,na = modules.period_convolution(key_diurnal_periods,obs_ref_time,diurnal_mags_spring,diurnal_phs_spring,mean_spring)
    na,na,na,diurnal_waveform_summer_extended,na = modules.period_convolution(key_diurnal_periods,obs_ref_time,diurnal_mags_summer,diurnal_phs_summer,mean_summer)
    na,na,na,diurnal_waveform_autumn_extended,na = modules.period_convolution(key_diurnal_periods,obs_ref_time,diurnal_mags_autumn,diurnal_phs_autumn,mean_autumn)
    na,na,na,diurnal_waveform_winter_extended,na = modules.period_convolution(key_diurnal_periods,obs_ref_time,diurnal_mags_winter,diurnal_phs_winter,mean_winter)
    
    #convolve all 
    full_mag,full_min_ph,full_max_ph,full_waveform,full_ff = modules.period_convolution(periodic_periods,obs_ref_time,mag,ph,mean)

    #convert phase to time
    daily_h3_ph = modules.radians_to_time(daily_h3_ph,6.)
    daily_h2_ph = modules.radians_to_time(daily_h2_ph,24./3.)
    daily_h1_ph = modules.radians_to_time(daily_h1_ph,12.)
    orig_daily_ph = modules.radians_to_time(orig_daily_ph,24.)
    diurnal_min_ph = modules.radians_to_time(diurnal_min_ph,24.)
    diurnal_max_ph = modules.radians_to_time(diurnal_max_ph,24.)
    daily_h3_ph_spring = modules.radians_to_time(daily_h3_ph_spring,6.)
    daily_h2_ph_spring = modules.radians_to_time(daily_h2_ph_spring,24./3.)
    daily_h1_ph_spring = modules.radians_to_time(daily_h1_ph_spring,12.)
    orig_daily_ph_spring = modules.radians_to_time(orig_daily_ph_spring,24.)
    diurnal_min_ph_spring = modules.radians_to_time(diurnal_min_ph_spring,24.)
    diurnal_max_ph_spring = modules.radians_to_time(diurnal_max_ph_spring,24.)
    daily_h3_ph_summer = modules.radians_to_time(daily_h3_ph_summer,6.)
    daily_h2_ph_summer = modules.radians_to_time(daily_h2_ph_summer,24./3.)
    daily_h1_ph_summer = modules.radians_to_time(daily_h1_ph_summer,12.)
    orig_daily_ph_summer = modules.radians_to_time(orig_daily_ph_summer,24.)
    diurnal_min_ph_summer = modules.radians_to_time(diurnal_min_ph_summer,24.)
    diurnal_max_ph_summer = modules.radians_to_time(diurnal_max_ph_summer,24.)
    daily_h3_ph_autumn = modules.radians_to_time(daily_h3_ph_autumn,6.)
    daily_h2_ph_autumn = modules.radians_to_time(daily_h2_ph_autumn,24./3.)
    daily_h1_ph_autumn = modules.radians_to_time(daily_h1_ph_autumn,12.)
    orig_daily_ph_autumn = modules.radians_to_time(orig_daily_ph_autumn,24.)
    diurnal_min_ph_autumn = modules.radians_to_time(diurnal_min_ph_autumn,24.)
    diurnal_max_ph_autumn = modules.radians_to_time(diurnal_max_ph_autumn,24.)
    daily_h3_ph_winter = modules.radians_to_time(daily_h3_ph_winter,6.)
    daily_h2_ph_winter = modules.radians_to_time(daily_h2_ph_winter,24./3.)
    daily_h1_ph_winter = modules.radians_to_time(daily_h1_ph_winter,12.)
    orig_daily_ph_winter = modules.radians_to_time(orig_daily_ph_winter,24.)
    diurnal_min_ph_winter = modules.radians_to_time(diurnal_min_ph_winter,24.)
    diurnal_max_ph_winter = modules.radians_to_time(diurnal_max_ph_winter,24.)
    seasonal_h3_ph = modules.radians_to_time(seasonal_h3_ph,3.)
    seasonal_h2_ph = modules.radians_to_time(seasonal_h2_ph,4.)
    seasonal_h1_ph = modules.radians_to_time(seasonal_h1_ph,6.)
    annual_ph = modules.radians_to_time(annual_ph,12.)
    seasonal_min_ph = modules.radians_to_time(seasonal_min_ph,12.)
    seasonal_max_ph = modules.radians_to_time(seasonal_max_ph,12.)
    seasonal_h3_ph_day = modules.radians_to_time(seasonal_h3_ph_day,3.)
    seasonal_h2_ph_day = modules.radians_to_time(seasonal_h2_ph_day,4.)
    seasonal_h1_ph_day = modules.radians_to_time(seasonal_h1_ph_day,6.)
    annual_ph_day = modules.radians_to_time(annual_ph_day,12.)
    seasonal_min_ph_day = modules.radians_to_time(seasonal_min_ph_day,12.)
    seasonal_max_ph_day = modules.radians_to_time(seasonal_max_ph_day,12.)
    seasonal_h3_ph_night = modules.radians_to_time(seasonal_h3_ph_night,3.)
    seasonal_h2_ph_night = modules.radians_to_time(seasonal_h2_ph_night,4.)
    seasonal_h1_ph_night = modules.radians_to_time(seasonal_h1_ph_night,6.)
    annual_ph_night = modules.radians_to_time(annual_ph_night,12.)
    seasonal_min_ph_night = modules.radians_to_time(seasonal_min_ph_night,12.)
    seasonal_max_ph_night = modules.radians_to_time(seasonal_max_ph_night,12.)

    #get % variability of periodic waveform with raw time series
    pc_var_daily = modules.periodic_variance_percent(vals,full_waveform,mag,mag[:4],valid)
    pc_var_daily_spring = modules.periodic_variance_percent(vals,full_waveform,mag,mag_spring,valid)
    pc_var_daily_summer = modules.periodic_variance_percent(vals,full_waveform,mag,mag_summer,valid)
    pc_var_daily_autumn = modules.periodic_variance_percent(vals,full_waveform,mag,mag_autumn,valid)
    pc_var_daily_winter = modules.periodic_variance_percent(vals,full_waveform,mag,mag_winter,valid)
    pc_var_seasonal = modules.periodic_variance_percent(vals,full_waveform,mag,mag[4:],valid)
    pc_var_seasonal_day = modules.periodic_variance_percent(vals,full_waveform,mag,mag_day,valid)
    pc_var_seasonal_night = modules.periodic_variance_percent(vals,full_waveform,mag,mag_night,valid)
    pc_var_full = modules.periodic_variance_percent(vals,full_waveform,mag,mag,valid)

    return (x,daily_h3_mag,daily_h3_ph,daily_h2_mag,daily_h2_ph,daily_h1_mag,daily_h1_ph,orig_daily_mag,orig_daily_ph,diurnal_mag,diurnal_min_ph,diurnal_max_ph,
            seasonal_h3_mag,seasonal_h3_ph,seasonal_h2_mag,seasonal_h2_ph,seasonal_h1_mag,seasonal_h1_ph,annual_mag,annual_ph,seasonal_mag,seasonal_min_ph,seasonal_max_ph,
            mean,diurnal_waveform,seasonal_waveform,full_waveform,diurnal_ff,seasonal_ff,full_ff,daily_h3_mag_spring,daily_h3_ph_spring,daily_h2_mag_spring,daily_h2_ph_spring,
            daily_h1_mag_spring,daily_h1_ph_spring,orig_daily_mag_spring,orig_daily_ph_spring,diurnal_mag_spring,diurnal_min_ph_spring,diurnal_max_ph_spring,daily_h3_mag_summer,
            daily_h3_ph_summer,daily_h2_mag_summer,daily_h2_ph_summer,daily_h1_mag_summer,daily_h1_ph_summer,orig_daily_mag_summer,orig_daily_ph_summer,diurnal_mag_summer,
            diurnal_min_ph_summer,diurnal_max_ph_summer,daily_h3_mag_autumn,daily_h3_ph_autumn,daily_h2_mag_autumn,daily_h2_ph_autumn,daily_h1_mag_autumn,daily_h1_ph_autumn,
            orig_daily_mag_autumn,orig_daily_ph_autumn,diurnal_mag_autumn,diurnal_min_ph_autumn,diurnal_max_ph_autumn,daily_h3_mag_winter,daily_h3_ph_winter,daily_h2_mag_winter,
            daily_h2_ph_winter,daily_h1_mag_winter,daily_h1_ph_winter,orig_daily_mag_winter,orig_daily_ph_winter,diurnal_mag_winter,diurnal_min_ph_winter,diurnal_max_ph_winter,
            mean_spring,mean_summer,mean_autumn,mean_winter,diurnal_ff_spring,diurnal_ff_summer,diurnal_ff_autumn,diurnal_ff_winter,diurnal_waveform_spring,diurnal_waveform_summer,
            diurnal_waveform_autumn,diurnal_waveform_winter,seasonal_h3_mag_day,seasonal_h3_ph_day,seasonal_h2_mag_day,seasonal_h2_ph_day,seasonal_h1_mag_day,seasonal_h1_ph_day,
            annual_mag_day,annual_ph_day,seasonal_mag_day,seasonal_min_ph_day,seasonal_max_ph_day,mean_day,seasonal_ff_day,seasonal_waveform_day,seasonal_h3_mag_night,seasonal_h3_ph_night,
            seasonal_h2_mag_night,seasonal_h2_ph_night,seasonal_h1_mag_night,seasonal_h1_ph_night,annual_mag_night,annual_ph_night,seasonal_mag_night,seasonal_min_ph_night,seasonal_max_ph_night,
            mean_night,seasonal_ff_night,seasonal_waveform_night,pc_var_daily,pc_var_daily_spring,pc_var_daily_summer,pc_var_daily_autumn,pc_var_daily_winter,pc_var_seasonal,pc_var_seasonal_day,
            pc_var_seasonal_night,pc_var_full)

for f in range(len(files)):
    obs_fname = files[f]
    last_split = obs_fname.split('/')[-1]
    start_year = last_split.split('_')[3] 
    end_year = last_split.split('_')[4]
    
    print obs_fname,start_year,end_year

    obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)

    #if no data for year range got to next range
    if len(obs_refs) == 0:
        continue

    if __name__ == '__main__':
        pool = multiprocessing.Pool(processes=16)
        results = [pool.apply_async(run_LSP, (obs_std_var[x],x)) for x in range(len(obs_refs))]
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
    daily_ff = []
    seasonal_ff = []
    full_ff = []
    daily_h3_mag_array_spring = []
    daily_h3_phase_array_spring = []
    daily_h2_mag_array_spring = []
    daily_h2_phase_array_spring = []
    daily_h1_mag_array_spring = []
    daily_h1_phase_array_spring = []
    orig_daily_mag_array_spring = []
    orig_daily_phase_array_spring = []
    daily_mag_array_spring = []
    daily_phase_min_array_spring = []
    daily_phase_max_array_spring = []
    daily_h3_mag_array_summer = []
    daily_h3_phase_array_summer = []
    daily_h2_mag_array_summer = []
    daily_h2_phase_array_summer = []
    daily_h1_mag_array_summer = []
    daily_h1_phase_array_summer = []
    orig_daily_mag_array_summer = []
    orig_daily_phase_array_summer = []
    daily_mag_array_summer = []
    daily_phase_min_array_summer = []
    daily_phase_max_array_summer = []
    daily_h3_mag_array_autumn = []
    daily_h3_phase_array_autumn = []
    daily_h2_mag_array_autumn = []
    daily_h2_phase_array_autumn = []
    daily_h1_mag_array_autumn = []
    daily_h1_phase_array_autumn = []
    orig_daily_mag_array_autumn = []
    orig_daily_phase_array_autumn = []
    daily_mag_array_autumn = []
    daily_phase_min_array_autumn = []
    daily_phase_max_array_autumn = []
    daily_h3_mag_array_winter = []
    daily_h3_phase_array_winter = []
    daily_h2_mag_array_winter = []
    daily_h2_phase_array_winter = []
    daily_h1_mag_array_winter = []
    daily_h1_phase_array_winter = []
    orig_daily_mag_array_winter = []
    orig_daily_phase_array_winter = []
    daily_mag_array_winter = []
    daily_phase_min_array_winter = []
    daily_phase_max_array_winter = []
    mean_spring_array = []
    mean_summer_array = []
    mean_autumn_array = []
    mean_winter_array = []
    daily_ff_spring = []
    daily_ff_summer = []
    daily_ff_autumn = []
    daily_ff_winter = []
    daily_waveform_spring = []
    daily_waveform_summer = []
    daily_waveform_autumn = []
    daily_waveform_winter = []
    seasonal_h3_mag_day = []
    seasonal_h3_ph_day = []
    seasonal_h2_mag_day = []
    seasonal_h2_ph_day = []
    seasonal_h1_mag_day = []
    seasonal_h1_ph_day = []
    annual_mag_day = []
    annual_ph_day = []
    seasonal_mag_day = []
    seasonal_min_ph_day = []
    seasonal_max_ph_day = []
    mean_day = []
    seasonal_ff_day = []
    seasonal_waveform_day = []
    seasonal_h3_mag_night = []
    seasonal_h3_ph_night = []
    seasonal_h2_mag_night = []
    seasonal_h2_ph_night = []
    seasonal_h1_mag_night = []
    seasonal_h1_ph_night = []
    annual_mag_night = []
    annual_ph_night = []
    seasonal_mag_night = []
    seasonal_min_ph_night = []
    seasonal_max_ph_night = []
    mean_night = []
    seasonal_ff_night = []
    seasonal_waveform_night = []
    pc_var_daily = []
    pc_var_daily_spring = []
    pc_var_daily_summer = []
    pc_var_daily_autumn = []
    pc_var_daily_winter = []
    pc_var_seasonal = []
    pc_var_seasonal_day = []
    pc_var_seasonal_night = []
    pc_var_full = []

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
        daily_ff.append(cut[27])
        seasonal_ff.append(cut[28])
        full_ff.append(cut[29])
        daily_h3_mag_array_spring.append(cut[30])
        daily_h3_phase_array_spring.append(cut[31])
        daily_h2_mag_array_spring.append(cut[32])
        daily_h2_phase_array_spring.append(cut[33])
        daily_h1_mag_array_spring.append(cut[34])
        daily_h1_phase_array_spring.append(cut[35])
        orig_daily_mag_array_spring.append(cut[36])
        orig_daily_phase_array_spring.append(cut[37])
        daily_mag_array_spring.append(cut[38])
        daily_phase_min_array_spring.append(cut[39])
        daily_phase_max_array_spring.append(cut[40])
        daily_h3_mag_array_summer.append(cut[41])
        daily_h3_phase_array_summer.append(cut[42])
        daily_h2_mag_array_summer.append(cut[43])
        daily_h2_phase_array_summer.append(cut[44])
        daily_h1_mag_array_summer.append(cut[45])
        daily_h1_phase_array_summer.append(cut[46])
        orig_daily_mag_array_summer.append(cut[47])
        orig_daily_phase_array_summer.append(cut[48])
        daily_mag_array_summer.append(cut[49])
        daily_phase_min_array_summer.append(cut[50])
        daily_phase_max_array_summer.append(cut[51])
        daily_h3_mag_array_autumn.append(cut[52])
        daily_h3_phase_array_autumn.append(cut[53])
        daily_h2_mag_array_autumn.append(cut[54])
        daily_h2_phase_array_autumn.append(cut[55])
        daily_h1_mag_array_autumn.append(cut[56])
        daily_h1_phase_array_autumn.append(cut[57])
        orig_daily_mag_array_autumn.append(cut[58])
        orig_daily_phase_array_autumn.append(cut[59])
        daily_mag_array_autumn.append(cut[60])
        daily_phase_min_array_autumn.append(cut[61])
        daily_phase_max_array_autumn.append(cut[62])
        daily_h3_mag_array_winter.append(cut[63])
        daily_h3_phase_array_winter.append(cut[64])
        daily_h2_mag_array_winter.append(cut[65])
        daily_h2_phase_array_winter.append(cut[66])
        daily_h1_mag_array_winter.append(cut[67])
        daily_h1_phase_array_winter.append(cut[68])
        orig_daily_mag_array_winter.append(cut[69])
        orig_daily_phase_array_winter.append(cut[70])
        daily_mag_array_winter.append(cut[71])
        daily_phase_min_array_winter.append(cut[72])
        daily_phase_max_array_winter.append(cut[73])
        mean_spring_array.append(cut[74])
        mean_summer_array.append(cut[75])
        mean_autumn_array.append(cut[76])
        mean_winter_array.append(cut[77])
        daily_ff_spring.append(cut[78])
        daily_ff_summer.append(cut[79])
        daily_ff_autumn.append(cut[80])
        daily_ff_winter.append(cut[81])
        daily_waveform_spring.append(cut[82])
        daily_waveform_summer.append(cut[83])
        daily_waveform_autumn.append(cut[84])
        daily_waveform_winter.append(cut[85])
        seasonal_h3_mag_day.append(cut[86])
        seasonal_h3_ph_day.append(cut[87]) 
        seasonal_h2_mag_day.append(cut[88]) 
        seasonal_h2_ph_day.append(cut[89]) 
        seasonal_h1_mag_day.append(cut[90]) 
        seasonal_h1_ph_day.append(cut[91]) 
        annual_mag_day.append(cut[92]) 
        annual_ph_day.append(cut[93]) 
        seasonal_mag_day.append(cut[94]) 
        seasonal_min_ph_day.append(cut[95])
        seasonal_max_ph_day.append(cut[96]) 
        mean_day.append(cut[97]) 
        seasonal_ff_day.append(cut[98]) 
        seasonal_waveform_day.append(cut[99]) 
        seasonal_h3_mag_night.append(cut[100]) 
        seasonal_h3_ph_night.append(cut[101]) 
        seasonal_h2_mag_night.append(cut[102]) 
        seasonal_h2_ph_night.append(cut[103]) 
        seasonal_h1_mag_night.append(cut[104]) 
        seasonal_h1_ph_night.append(cut[105]) 
        annual_mag_night.append(cut[106]) 
        annual_ph_night.append(cut[107]) 
        seasonal_mag_night.append(cut[108]) 
        seasonal_min_ph_night.append(cut[109]) 
        seasonal_max_ph_night.append(cut[110]) 
        mean_night.append(cut[111]) 
        seasonal_ff_night.append(cut[112]) 
        seasonal_waveform_night.append(cut[113]) 
        pc_var_daily.append(cut[114]) 
        pc_var_daily_spring.append(cut[115]) 
        pc_var_daily_summer.append(cut[116]) 
        pc_var_daily_autumn.append(cut[117]) 
        pc_var_daily_winter.append(cut[118]) 
        pc_var_seasonal.append(cut[119]) 
        pc_var_seasonal_day.append(cut[120]) 
        pc_var_seasonal_night.append(cut[121]) 
        pc_var_full.append(cut[122]) 

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
    daily_ff = np.array(daily_ff)
    seasonal_ff = np.array(seasonal_ff)
    full_ff = np.array(full_ff)
    daily_h3_mag_array_spring = np.array(daily_h3_mag_array_spring)
    daily_h3_phase_array_spring = np.array(daily_h3_phase_array_spring)
    daily_h2_mag_array_spring = np.array(daily_h2_mag_array_spring)
    daily_h2_phase_array_spring = np.array(daily_h2_phase_array_spring)
    daily_h1_mag_array_spring = np.array(daily_h1_mag_array_spring)
    daily_h1_phase_array_spring = np.array(daily_h1_phase_array_spring)
    orig_daily_mag_array_spring = np.array(orig_daily_mag_array_spring)
    orig_daily_phase_array_spring = np.array(orig_daily_phase_array_spring)
    daily_mag_array_spring = np.array(daily_mag_array_spring)
    daily_phase_min_array_spring = np.array(daily_phase_min_array_spring)
    daily_phase_max_array_spring = np.array(daily_phase_max_array_spring)
    daily_h3_mag_array_summer = np.array(daily_h3_mag_array_summer)
    daily_h3_phase_array_summer = np.array(daily_h3_phase_array_summer)
    daily_h2_mag_array_summer = np.array(daily_h2_mag_array_summer)
    daily_h2_phase_array_summer = np.array(daily_h2_phase_array_summer)
    daily_h1_mag_array_summer = np.array(daily_h1_mag_array_summer)
    daily_h1_phase_array_summer = np.array(daily_h1_phase_array_summer)
    orig_daily_mag_array_summer = np.array(orig_daily_mag_array_summer)
    orig_daily_phase_array_summer = np.array(orig_daily_phase_array_summer)
    daily_mag_array_summer = np.array(daily_mag_array_summer)
    daily_phase_min_array_summer = np.array(daily_phase_min_array_summer)
    daily_phase_max_array_summer = np.array(daily_phase_max_array_summer) 
    daily_h3_mag_array_autumn = np.array(daily_h3_mag_array_autumn)
    daily_h3_phase_array_autumn = np.array(daily_h3_phase_array_autumn)
    daily_h2_mag_array_autumn = np.array(daily_h2_mag_array_autumn)
    daily_h2_phase_array_autumn = np.array(daily_h2_phase_array_autumn)
    daily_h1_mag_array_autumn = np.array(daily_h1_mag_array_autumn)
    daily_h1_phase_array_autumn = np.array(daily_h1_phase_array_autumn)
    orig_daily_mag_array_autumn = np.array(orig_daily_mag_array_autumn)
    orig_daily_phase_array_autumn = np.array(orig_daily_phase_array_autumn)
    daily_mag_array_autumn = np.array(daily_mag_array_autumn)
    daily_phase_min_array_autumn = np.array(daily_phase_min_array_autumn)
    daily_phase_max_array_autumn = np.array(daily_phase_max_array_autumn)
    daily_h3_mag_array_winter = np.array(daily_h3_mag_array_winter)
    daily_h3_phase_array_winter = np.array(daily_h3_phase_array_winter)
    daily_h2_mag_array_winter = np.array(daily_h2_mag_array_winter)
    daily_h2_phase_array_winter = np.array(daily_h2_phase_array_winter)
    daily_h1_mag_array_winter = np.array(daily_h1_mag_array_winter)
    daily_h1_phase_array_winter = np.array(daily_h1_phase_array_winter)
    orig_daily_mag_array_winter = np.array(orig_daily_mag_array_winter)
    orig_daily_phase_array_winter = np.array(orig_daily_phase_array_winter)
    daily_mag_array_winter = np.array(daily_mag_array_winter)
    daily_phase_min_array_winter = np.array(daily_phase_min_array_winter)
    daily_phase_max_array_winter = np.array(daily_phase_max_array_winter)
    mean_spring_array = np.array(mean_spring_array)
    mean_summer_array = np.array(mean_summer_array)
    mean_autumn_array = np.array(mean_autumn_array)
    mean_winter_array = np.array(mean_winter_array)
    daily_ff_spring = np.array(daily_ff_spring)
    daily_ff_summer = np.array(daily_ff_summer)
    daily_ff_autumn = np.array(daily_ff_autumn)
    daily_ff_winter = np.array(daily_ff_winter)
    daily_waveform_spring = np.array(daily_waveform_spring)
    daily_waveform_summer = np.array(daily_waveform_summer)
    daily_waveform_autumn = np.array(daily_waveform_autumn)
    daily_waveform_winter = np.array(daily_waveform_winter)
    seasonal_h3_mag_day = np.array(seasonal_h3_mag_day)
    seasonal_h3_ph_day = np.array(seasonal_h3_ph_day)
    seasonal_h2_mag_day = np.array(seasonal_h2_mag_day)
    seasonal_h2_ph_day = np.array(seasonal_h2_ph_day)
    seasonal_h1_mag_day = np.array(seasonal_h1_mag_day)
    seasonal_h1_ph_day = np.array(seasonal_h1_ph_day)
    annual_mag_day = np.array(annual_mag_day)
    annual_ph_day = np.array(annual_ph_day)
    seasonal_mag_day = np.array(seasonal_mag_day)
    seasonal_min_ph_day = np.array(seasonal_min_ph_day)
    seasonal_max_ph_day = np.array(seasonal_max_ph_day)
    mean_day = np.array(mean_day)
    seasonal_ff_day = np.array(seasonal_ff_day)
    seasonal_waveform_day = np.array(seasonal_waveform_day)
    seasonal_h3_mag_night = np.array(seasonal_h3_mag_night)
    seasonal_h3_ph_night = np.array(seasonal_h3_ph_night)
    seasonal_h2_mag_night = np.array(seasonal_h2_mag_night)
    seasonal_h2_ph_night = np.array(seasonal_h2_ph_night)
    seasonal_h1_mag_night = np.array(seasonal_h1_mag_night)
    seasonal_h1_ph_night = np.array(seasonal_h1_ph_night)
    annual_mag_night = np.array(annual_mag_night)
    annual_ph_night = np.array(annual_ph_night)
    seasonal_mag_night = np.array(seasonal_mag_night)
    seasonal_min_ph_night = np.array(seasonal_min_ph_night)
    seasonal_max_ph_night = np.array(seasonal_max_ph_night)
    mean_night = np.array(mean_night)
    seasonal_ff_night = np.array(seasonal_ff_night)
    seasonal_waveform_night = np.array(seasonal_waveform_night)
    pc_var_daily = np.array(pc_var_daily)
    pc_var_daily_spring = np.array(pc_var_daily_spring)
    pc_var_daily_summer = np.array(pc_var_daily_summer)
    pc_var_daily_autumn = np.array(pc_var_daily_autumn)
    pc_var_daily_winter = np.array(pc_var_daily_winter)
    pc_var_seasonal = np.array(pc_var_seasonal)
    pc_var_seasonal_day = np.array(pc_var_seasonal_day)
    pc_var_seasonal_night = np.array(pc_var_seasonal_night)
    pc_var_full = np.array(pc_var_full)


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
    daily_ff = daily_ff[indices_array]
    seasonal_ff = seasonal_ff[indices_array]
    full_ff = full_ff[indices_array]
    daily_h3_mag_array_spring = daily_h3_mag_array_spring[indices_array]
    daily_h3_phase_array_spring = daily_h3_phase_array_spring[indices_array]
    daily_h2_mag_array_spring = daily_h2_mag_array_spring[indices_array]
    daily_h2_phase_array_spring = daily_h2_phase_array_spring[indices_array]
    daily_h1_mag_array_spring = daily_h1_mag_array_spring[indices_array]
    daily_h1_phase_array_spring = daily_h1_phase_array_spring[indices_array]
    orig_daily_mag_array_spring = orig_daily_mag_array_spring[indices_array]
    orig_daily_phase_array_spring = orig_daily_phase_array_spring[indices_array]
    daily_mag_array_spring = daily_mag_array_spring[indices_array]
    daily_phase_min_array_spring = daily_phase_min_array_spring[indices_array]
    daily_phase_max_array_spring = daily_phase_max_array_spring[indices_array]
    daily_h3_mag_array_summer = daily_h3_mag_array_summer[indices_array]
    daily_h3_phase_array_summer = daily_h3_phase_array_summer[indices_array]
    daily_h2_mag_array_summer = daily_h2_mag_array_summer[indices_array]
    daily_h2_phase_array_summer = daily_h2_phase_array_summer[indices_array]
    daily_h1_mag_array_summer = daily_h1_mag_array_summer[indices_array]
    daily_h1_phase_array_summer = daily_h1_phase_array_summer[indices_array]
    orig_daily_mag_array_summer = orig_daily_mag_array_summer[indices_array]
    orig_daily_phase_array_summer = orig_daily_phase_array_summer[indices_array]
    daily_mag_array_summer = daily_mag_array_summer[indices_array]
    daily_phase_min_array_summer = daily_phase_min_array_summer[indices_array]
    daily_phase_max_array_summer = daily_phase_max_array_summer[indices_array]
    daily_h3_mag_array_autumn = daily_h3_mag_array_autumn[indices_array]
    daily_h3_phase_array_autumn = daily_h3_phase_array_autumn[indices_array]
    daily_h2_mag_array_autumn = daily_h2_mag_array_autumn[indices_array]
    daily_h2_phase_array_autumn = daily_h2_phase_array_autumn[indices_array]
    daily_h1_mag_array_autumn = daily_h1_mag_array_autumn[indices_array]
    daily_h1_phase_array_autumn = daily_h1_phase_array_autumn[indices_array]
    orig_daily_mag_array_autumn = orig_daily_mag_array_autumn[indices_array]
    orig_daily_phase_array_autumn = orig_daily_phase_array_autumn[indices_array]
    daily_mag_array_autumn = daily_mag_array_autumn[indices_array]
    daily_phase_min_array_autumn = daily_phase_min_array_autumn[indices_array]
    daily_phase_max_array_autumn = daily_phase_max_array_autumn[indices_array]
    daily_h3_mag_array_winter = daily_h3_mag_array_winter[indices_array]
    daily_h3_phase_array_winter = daily_h3_phase_array_winter[indices_array]
    daily_h2_mag_array_winter = daily_h2_mag_array_winter[indices_array]
    daily_h2_phase_array_winter = daily_h2_phase_array_winter[indices_array]
    daily_h1_mag_array_winter = daily_h1_mag_array_winter[indices_array]
    daily_h1_phase_array_winter = daily_h1_phase_array_winter[indices_array]
    orig_daily_mag_array_winter = orig_daily_mag_array_winter[indices_array]
    orig_daily_phase_array_winter = orig_daily_phase_array_winter[indices_array]
    daily_mag_array_winter = daily_mag_array_winter[indices_array]
    daily_phase_min_array_winter = daily_phase_min_array_winter[indices_array]
    daily_phase_max_array_winter = daily_phase_max_array_winter[indices_array]
    mean_spring_array = mean_spring_array[indices_array]
    mean_summer_array = mean_summer_array[indices_array]
    mean_autumn_array = mean_autumn_array[indices_array]
    mean_winter_array = mean_winter_array[indices_array]
    daily_ff_spring = daily_ff_spring[indices_array]
    daily_ff_summer = daily_ff_summer[indices_array]
    daily_ff_autumn = daily_ff_autumn[indices_array]
    daily_ff_winter = daily_ff_winter[indices_array]
    daily_waveform_spring = daily_waveform_spring[indices_array]
    daily_waveform_summer = daily_waveform_summer[indices_array]
    daily_waveform_autumn = daily_waveform_autumn[indices_array]
    daily_waveform_winter = daily_waveform_winter[indices_array]
    seasonal_h3_mag_day = seasonal_h3_mag_day[indices_array]
    seasonal_h3_ph_day = seasonal_h3_ph_day[indices_array]
    seasonal_h2_mag_day = seasonal_h2_mag_day[indices_array]
    seasonal_h2_ph_day = seasonal_h2_ph_day[indices_array]
    seasonal_h1_mag_day = seasonal_h1_mag_day[indices_array]
    seasonal_h1_ph_day = seasonal_h1_ph_day[indices_array]
    annual_mag_day = annual_mag_day[indices_array]
    annual_ph_day = annual_ph_day[indices_array]
    seasonal_mag_day = seasonal_mag_day[indices_array]
    seasonal_min_ph_day = seasonal_min_ph_day[indices_array]
    seasonal_max_ph_day = seasonal_max_ph_day[indices_array]
    mean_day = mean_day[indices_array]
    seasonal_ff_day = seasonal_ff_day[indices_array]
    seasonal_waveform_day = seasonal_waveform_day[indices_array]
    seasonal_h3_mag_night = seasonal_h3_mag_night[indices_array]
    seasonal_h3_ph_night = seasonal_h3_ph_night[indices_array]
    seasonal_h2_mag_night = seasonal_h2_mag_night[indices_array]
    seasonal_h2_ph_night = seasonal_h2_ph_night[indices_array]
    seasonal_h1_mag_night = seasonal_h1_mag_night[indices_array]
    seasonal_h1_ph_night = seasonal_h1_ph_night[indices_array]
    annual_mag_night = annual_mag_night[indices_array]
    annual_ph_night = annual_ph_night[indices_array]
    seasonal_mag_night = seasonal_mag_night[indices_array]
    seasonal_min_ph_night = seasonal_min_ph_night[indices_array]
    seasonal_max_ph_night = seasonal_max_ph_night[indices_array]
    mean_night = mean_night[indices_array]
    seasonal_ff_night = seasonal_ff_night[indices_array]
    seasonal_waveform_night = seasonal_waveform_night[indices_array]
    pc_var_daily = pc_var_daily[indices_array]
    pc_var_daily_spring = pc_var_daily_spring[indices_array]
    pc_var_daily_summer = pc_var_daily_summer[indices_array]
    pc_var_daily_autumn = pc_var_daily_autumn[indices_array]
    pc_var_daily_winter = pc_var_daily_winter[indices_array]
    pc_var_seasonal = pc_var_seasonal[indices_array]
    pc_var_seasonal_day = pc_var_seasonal_day[indices_array]
    pc_var_seasonal_night = pc_var_seasonal_night[indices_array]
    pc_var_full = pc_var_full[indices_array]
    
    #save out sig periods to netcdf
    root_grp_period = Dataset('obs_sig_periods_%s_%s.nc'%(start_year,end_year), 'w')
    root_grp_period.description = 'Amplitudes and Phases for key periods from Lomb-Scargle Periodogram analysis of Observations - Program written by Dene Bowdalo'
    root_grp_period.createDimension('diurnal', 24)

    if timeres == 'H':
        root_grp_period.createDimension('seasonal', 8766)
    elif timeres == 'D':
        root_grp_period.createDimension('seasonal',365)
    elif timeres == 'M':
        root_grp_period.createDimension('seasonal',12)

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
        ref_period.average_spring = mean_spring_array[i]
        ref_period.average_summer = mean_summer_array[i]
        ref_period.average_autumn = mean_autumn_array[i]
        ref_period.average_winter = mean_winter_array[i]
        ref_period.average_day = mean_day[i]
        ref_period.average_night = mean_night[i]
        ref_period.daily_harmonic3_amplitude_spring = daily_h3_mag_array_spring[i]
        ref_period.daily_harmonic2_amplitude_spring = daily_h2_mag_array_spring[i]
        ref_period.daily_harmonic1_amplitude_spring = daily_h1_mag_array_spring[i]
        ref_period.original_daily_amplitude_spring = orig_daily_mag_array_spring[i]
        ref_period.daily_amplitude_spring = daily_mag_array_spring[i]
        ref_period.daily_harmonic3_phase_spring = daily_h3_phase_array_spring[i]
        ref_period.daily_harmonic2_phase_spring = daily_h2_phase_array_spring[i]
        ref_period.daily_harmonic1_phase_spring = daily_h1_phase_array_spring[i]
        ref_period.original_daily_phase_spring = orig_daily_phase_array_spring[i]
        ref_period.daily_phase_spring = daily_phase_max_array_spring[i]
        ref_period.daily_phase_min_spring = daily_phase_min_array_spring[i]
        ref_period.daily_harmonic3_amplitude_summer = daily_h3_mag_array_summer[i]
        ref_period.daily_harmonic2_amplitude_summer = daily_h2_mag_array_summer[i]
        ref_period.daily_harmonic1_amplitude_summer = daily_h1_mag_array_summer[i]
        ref_period.original_daily_amplitude_summer = orig_daily_mag_array_summer[i]
        ref_period.daily_amplitude_summer = daily_mag_array_summer[i]
        ref_period.daily_harmonic3_phase_summer = daily_h3_phase_array_summer[i]
        ref_period.daily_harmonic2_phase_summer = daily_h2_phase_array_summer[i]
        ref_period.daily_harmonic1_phase_summer = daily_h1_phase_array_summer[i]
        ref_period.original_daily_phase_summer = orig_daily_phase_array_summer[i]
        ref_period.daily_phase_summer = daily_phase_max_array_summer[i]
        ref_period.daily_phase_min_summer = daily_phase_min_array_summer[i]
        ref_period.daily_harmonic3_amplitude_autumn = daily_h3_mag_array_autumn[i]
        ref_period.daily_harmonic2_amplitude_autumn = daily_h2_mag_array_autumn[i]
        ref_period.daily_harmonic1_amplitude_autumn = daily_h1_mag_array_autumn[i]
        ref_period.original_daily_amplitude_autumn = orig_daily_mag_array_autumn[i]
        ref_period.daily_amplitude_autumn = daily_mag_array_autumn[i]
        ref_period.daily_harmonic3_phase_autumn = daily_h3_phase_array_autumn[i]
        ref_period.daily_harmonic2_phase_autumn = daily_h2_phase_array_autumn[i]
        ref_period.daily_harmonic1_phase_autumn = daily_h1_phase_array_autumn[i]
        ref_period.original_daily_phase_autumn = orig_daily_phase_array_autumn[i]
        ref_period.daily_phase_autumn = daily_phase_max_array_autumn[i]
        ref_period.daily_phase_min_autumn = daily_phase_min_array_autumn[i]
        ref_period.daily_harmonic3_amplitude_winter = daily_h3_mag_array_winter[i]
        ref_period.daily_harmonic2_amplitude_winter = daily_h2_mag_array_winter[i]
        ref_period.daily_harmonic1_amplitude_winter = daily_h1_mag_array_winter[i]
        ref_period.original_daily_amplitude_winter = orig_daily_mag_array_winter[i]
        ref_period.daily_amplitude_winter = daily_mag_array_winter[i]
        ref_period.daily_harmonic3_phase_winter = daily_h3_phase_array_winter[i]
        ref_period.daily_harmonic2_phase_winter = daily_h2_phase_array_winter[i]
        ref_period.daily_harmonic1_phase_winter = daily_h1_phase_array_winter[i]
        ref_period.original_daily_phase_winter = orig_daily_phase_array_winter[i]
        ref_period.daily_phase_winter = daily_phase_max_array_winter[i]
        ref_period.daily_phase_min_winter = daily_phase_min_array_winter[i]
        ref_period.seasonal_harmonic3_amplitude_day = seasonal_h3_mag_day[i]
        ref_period.seasonal_harmonic2_amplitude_day = seasonal_h2_mag_day[i]
        ref_period.seasonal_harmonic1_amplitude_day = seasonal_h1_mag_day[i]
        ref_period.annual_amplitude_day = annual_mag_day[i]
        ref_period.seasonal_amplitude_day = seasonal_mag_day[i]
        ref_period.seasonal_harmonic3_phase_day = seasonal_h3_ph_day[i]
        ref_period.seasonal_harmonic2_phase_day = seasonal_h2_ph_day[i]
        ref_period.seasonal_harmonic1_phase_day = seasonal_h1_ph_day[i]
        ref_period.annual_phase_day = annual_ph_day[i]
        ref_period.seasonal_phase_day = seasonal_max_ph_day[i]
        ref_period.seasonal_phase_min_day = seasonal_min_ph_day[i]
        ref_period.seasonal_harmonic3_amplitude_night = seasonal_h3_mag_night[i]
        ref_period.seasonal_harmonic2_amplitude_night = seasonal_h2_mag_night[i]
        ref_period.seasonal_harmonic1_amplitude_night = seasonal_h1_mag_night[i]
        ref_period.annual_amplitude_night = annual_mag_night[i]
        ref_period.seasonal_amplitude_night = seasonal_mag_night[i]
        ref_period.seasonal_harmonic3_phase_night = seasonal_h3_ph_night[i]
        ref_period.seasonal_harmonic2_phase_night = seasonal_h2_ph_night[i]
        ref_period.seasonal_harmonic1_phase_night = seasonal_h1_ph_night[i]
        ref_period.annual_phase_night = annual_ph_night[i]
        ref_period.seasonal_phase_night = seasonal_max_ph_night[i]
        ref_period.seasonal_phase_min_night = seasonal_min_ph_night[i]
        d_w = ref_period.createVariable('daily_waveform', 'f4', ('diurnal',))
        d_w_spring = ref_period.createVariable('daily_waveform_spring', 'f4', ('diurnal',))
        d_w_summer = ref_period.createVariable('daily_waveform_summer', 'f4', ('diurnal',))
        d_w_autumn = ref_period.createVariable('daily_waveform_autumn', 'f4', ('diurnal',))
        d_w_winter = ref_period.createVariable('daily_waveform_winter', 'f4', ('diurnal',))
        s_w = ref_period.createVariable('seasonal_waveform', 'f4', ('seasonal',))
        s_w_day = ref_period.createVariable('seasonal_waveform_day', 'f4', ('seasonal',))
        s_w_night = ref_period.createVariable('seasonal_waveform_night', 'f4', ('seasonal',))
        a_w = ref_period.createVariable('all_waveform', 'f4', ('all',))
        d_w[:] = daily_waveform[i]
        d_w_spring[:] = daily_waveform_spring[i]
        d_w_summer[:] = daily_waveform_summer[i]
        d_w_autumn[:] = daily_waveform_autumn[i]
        d_w_winter[:] = daily_waveform_winter[i]
        s_w[:] = seasonal_waveform[i]
        s_w_day[:] = seasonal_waveform_day[i]
        s_w_night[:] = seasonal_waveform_night[i]
        a_w[:] = full_waveform[i]
        ref_period.daily_ff = daily_ff[i]
        ref_period.daily_ff_spring = daily_ff_spring[i]
        ref_period.daily_ff_summer = daily_ff_summer[i]
        ref_period.daily_ff_autumn = daily_ff_autumn[i]
        ref_period.daily_ff_winter = daily_ff_winter[i]
        ref_period.seasonal_ff = seasonal_ff[i]
        ref_period.seasonal_ff_day = seasonal_ff_day[i]
        ref_period.seasonal_ff_night = seasonal_ff_night[i]
        ref_period.full_ff = full_ff[i]
        ref_period.periodic_variance_daily = pc_var_daily[i]
        ref_period.periodic_variance_daily_spring = pc_var_daily_spring[i]
        ref_period.periodic_variance_daily_summer = pc_var_daily_summer[i]
        ref_period.periodic_variance_daily_autumn = pc_var_daily_autumn[i]
        ref_period.periodic_variance_daily_winter = pc_var_daily_winter[i]
        ref_period.periodic_variance_seasonal = pc_var_seasonal[i]
        ref_period.periodic_variance_seasonal_day = pc_var_seasonal_day[i]
        ref_period.periodic_variance_seasonal_night = pc_var_seasonal_night[i]
        ref_period.periodic_variance_all = pc_var_full[i]
    
    root_grp_period.close()

end = datetime.datetime.now() 
diff = end - start
print 'time = ', diff.seconds

