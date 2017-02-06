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

    site_lon = lon_c[lon_i]
    site_lat = lat_c[lat_i]

    #check model vals are valid
    valid = vals >= 0
    vals = vals[valid]
    orig_time = model_ref_time
    valid_times = model_ref_time[valid]
    orig_datetimes = model_datetime_time
    valid_datetimes = model_datetime_time[valid]
    
    if timeres == 'H':
        full_times_year = model_ref_time[:8766]
    elif timeres == 'D':
        full_times_year = model_ref_time[:365]
    elif timeres == 'M':
        full_times_year = model_ref_time[:12]
    full_times_day = model_ref_time[:24]

    #convert site_lon to 0 to 360 degs
    #keep -180 to 180 lon as orig_lon
    orig_lon = site_lon
    if site_lon < 0:
        site_lon = 360-np.abs(site_lon)
    
    #transform factor for conversion from UTC time to solar time 
    time_diff = lon_step_time*site_lon
    if time_diff > 12:
        time_diff = time_diff-24
        
    #cut vals into seasons
    valid_times_spring,vals_spring,valid_inds_spring,valid_times_summer,vals_summer,valid_inds_summer,valid_times_autumn,vals_autumn,valid_inds_autumn,valid_times_winter,vals_winter,valid_inds_winter = modules.cut_season(valid_datetimes,valid_times,vals)
             
    #cut vals into day/night (use -180 to 180 lon)
    valid_times_day,vals_day,valid_times_night,vals_night = modules.cut_daynight(valid_datetimes,valid_times,vals,site_lat,orig_lon,timeres)
    
    #make time start from 0    
    valid_times_from0 = modules.time_from0(valid_times)
    valid_times_from0_spring = modules.time_from0(valid_times_spring)
    valid_times_from0_summer = modules.time_from0(valid_times_summer)
    valid_times_from0_autumn = modules.time_from0(valid_times_autumn)
    valid_times_from0_winter = modules.time_from0(valid_times_winter)
    valid_times_from0_day = modules.time_from0(valid_times_day)
    valid_times_from0_night = modules.time_from0(valid_times_night)

    key_diurnal_periods = [1./12.,1./11.,1./10.,1./9.,1./8.,1./7.,1./6.,1./5.,1./4.,1./3.,1./2.,1.]
    #key_seasonal_periods = [365.25/12.,365.25/11.,365.25/10.,365.25/9.,365.25/8.,365.25/7.,365.25/6.,365.25/5.,365.25/4.,365.25/3.,365.25/2.,365.25]
    key_seasonal_periods = [365.25/4.,365.25/3.,365.25/2.,365.25]
    #periodic_periods = [1./12.,1./11.,1./10.,1./9.,1./8.,1./7.,1./6.,1./5.,1./4.,1./3.,1./2.,1.,365.25/12.,365.25/11.,365.25/10.,365.25/9.,365.25/8.,365.25/7.,365.25/6.,365.25/5.,365.25/4.,365.25/3.,365.25/2.,365.25]
    periodic_periods = [1./12.,1./11.,1./10.,1./9.,1./8.,1./7.,1./6.,1./5.,1./4.,1./3.,1./2.,1.,365.25/4.,365.25/3.,365.25/2.,365.25]
    periods,mag,ph,fr,fi = modules.lomb_specific(valid_times_from0,vals,w=True,key_periods=periodic_periods)
    periods_spring,mag_spring,ph_spring,fr_spring,fi_spring = modules.lomb_specific(valid_times_from0_spring,vals_spring,w=True,key_periods=key_diurnal_periods)
    periods_summer,mag_summer,ph_summer,fr_summer,fi_summer = modules.lomb_specific(valid_times_from0_summer,vals_summer,w=True,key_periods=key_diurnal_periods)
    periods_autumn,mag_autumn,ph_autumn,fr_autumn,fi_autumn = modules.lomb_specific(valid_times_from0_autumn,vals_autumn,w=True,key_periods=key_diurnal_periods)
    periods_winter,mag_winter,ph_winter,fr_winter,fi_winter = modules.lomb_specific(valid_times_from0_winter,vals_winter,w=True,key_periods=key_diurnal_periods)
    periods_day,mag_day,ph_day,fr_day,fi_day = modules.lomb_specific(valid_times_from0_day,vals_day,w=True,key_periods=key_seasonal_periods)
    periods_night,mag_night,ph_night,fr_night,fi_night = modules.lomb_specific(valid_times_from0_night,vals_night,w=True,key_periods=key_seasonal_periods)

    #process RAW STATS
    #STANDARD
    mean,p1,p5,p25,p50,p75,p95,p99 = modules.raw_stats(vals)
    #SEASON
    mean_spring,p1_spring,p5_spring,p25_spring,p50_spring,p75_spring,p95_spring,p99_spring = modules.raw_stats(vals_spring)
    mean_summer,p1_summer,p5_summer,p25_summer,p50_summer,p75_summer,p95_summer,p99_summer = modules.raw_stats(vals_summer)
    mean_autumn,p1_autumn,p5_autumn,p25_autumn,p50_autumn,p75_autumn,p95_autumn,p99_autumn = modules.raw_stats(vals_autumn)
    mean_winter,p1_winter,p5_winter,p25_winter,p50_winter,p75_winter,p95_winter,p99_winter = modules.raw_stats(vals_winter)
    #DAY/NIGHT
    mean_day,p1_day,p5_day,p25_day,p50_day,p75_day,p95_day,p99_day = modules.raw_stats(vals_day)
    mean_night,p1_night,p5_night,p25_night,p50_night,p75_night,p95_night,p99_night = modules.raw_stats(vals_night)

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
    
    #put mags and phases into appropriate arrays
    diurnal_mags = mag[:12]
    diurnal_phs = ph[:12]
    diurnal_mags_spring = mag_spring[:]
    diurnal_phs_spring = ph_spring[:]
    diurnal_mags_summer = mag_summer[:]
    diurnal_phs_summer = ph_summer[:]
    diurnal_mags_autumn = mag_autumn[:]
    diurnal_phs_autumn = ph_autumn[:]
    diurnal_mags_winter = mag_winter[:]
    diurnal_phs_winter = ph_winter[:]
    seasonal_mags = mag[12:]
    seasonal_mags_day = mag_day[:]
    seasonal_mags_night = mag_night[:]
    seasonal_phs = ph[12:]
    seasonal_phs_day = ph_day[:]
    seasonal_phs_night = ph_night[:]
    all_mags_spring = np.concatenate((diurnal_mags_spring,seasonal_mags))
    all_mags_summer = np.concatenate((diurnal_mags_summer,seasonal_mags))
    all_mags_autumn = np.concatenate((diurnal_mags_autumn,seasonal_mags))
    all_mags_winter = np.concatenate((diurnal_mags_winter,seasonal_mags))
    all_phs_spring = np.concatenate((diurnal_phs_spring,seasonal_phs))
    all_phs_summer = np.concatenate((diurnal_phs_summer,seasonal_phs))
    all_phs_autumn = np.concatenate((diurnal_phs_autumn,seasonal_phs))
    all_phs_winter = np.concatenate((diurnal_phs_winter,seasonal_phs))
    all_mags_day = np.concatenate((diurnal_mags,seasonal_mags_day))
    all_mags_night = np.concatenate((diurnal_mags,seasonal_mags_night))

    #set key individual mags and phases
    daily_h11_mag = mag[0]
    daily_h10_mag = mag[1]
    daily_h9_mag = mag[2]
    daily_h8_mag = mag[3]
    daily_h7_mag = mag[4]
    daily_h6_mag = mag[5]
    daily_h5_mag = mag[6]
    daily_h4_mag = mag[7]
    daily_h3_mag = mag[8]
    daily_h2_mag = mag[9]
    daily_h1_mag = mag[10]
    orig_daily_mag = mag[11]
    
    #seasonal_h11_mag = mag[12]
    #seasonal_h10_mag = mag[13]
    #seasonal_h9_mag = mag[14]
    #seasonal_h8_mag = mag[15]
    #seasonal_h7_mag = mag[16]
    #seasonal_h6_mag = mag[17]
    #seasonal_h5_mag = mag[18]
    #seasonal_h4_mag = mag[19]
    seasonal_h3_mag = mag[12]
    seasonal_h2_mag = mag[13]
    seasonal_h1_mag = mag[14]
    annual_mag = mag[15]

    #convolve annual cycle and harmonics to seasonal waveform for 1 year
    seasonal_mag,seasonal_ph,seasonal_waveform = modules.period_convolution(key_seasonal_periods,full_times_year,seasonal_mags,seasonal_phs,mean)
    seasonal_mag_day,seasonal_ph_day,seasonal_waveform_day = modules.period_convolution(key_seasonal_periods,full_times_year,seasonal_mags_day,seasonal_phs_day,mean_day)
    seasonal_mag_night,seasonal_ph_night,seasonal_waveform_night = modules.period_convolution(key_seasonal_periods,full_times_year,seasonal_mags_night,seasonal_phs_night,mean_night)
    remainder = np.mod(len(orig_time),len(seasonal_waveform))
    if remainder > 0:
        seasonal_waveform_extended = np.append(list(seasonal_waveform)*(len(orig_time)//len(seasonal_waveform)),seasonal_waveform[-remainder:])
    else:
        seasonal_waveform_extended = np.array(list(seasonal_waveform)*(len(orig_time)//len(seasonal_waveform)))

    #convolve diurnal cycle and harmonics to diurnal waveform for 1 day
    diurnal_mag,diurnal_ph,diurnal_ave_waveform = modules.period_convolution(key_diurnal_periods,full_times_day,diurnal_mags,diurnal_phs,mean)
    diurnal_mag_spring,diurnal_ph_spring,diurnal_waveform_spring = modules.period_convolution(key_diurnal_periods,full_times_day,diurnal_mags_spring,diurnal_phs_spring,mean_spring)
    diurnal_mag_summer,diurnal_ph_summer,diurnal_waveform_summer = modules.period_convolution(key_diurnal_periods,full_times_day,diurnal_mags_summer,diurnal_phs_summer,mean_summer)
    diurnal_mag_autumn,diurnal_ph_autumn,diurnal_waveform_autumn = modules.period_convolution(key_diurnal_periods,full_times_day,diurnal_mags_autumn,diurnal_phs_autumn,mean_autumn)
    diurnal_mag_winter,diurnal_ph_winter,diurnal_waveform_winter = modules.period_convolution(key_diurnal_periods,full_times_day,diurnal_mags_winter,diurnal_phs_winter,mean_winter)
    remainder = np.mod(len(orig_time),len(diurnal_ave_waveform))
    if remainder > 0:
        diurnal_ave_waveform_extended = np.append(list(diurnal_ave_waveform)*(len(orig_time)//len(diurnal_ave_waveform)),diurnal_ave_waveform[-remainder:])
    else:
        diurnal_ave_waveform_extended = np.array(list(diurnal_ave_waveform)*(len(orig_time)//len(diurnal_ave_waveform)))
    
    #convert diurnal phase and waveforms from UTC to local solar time (in radians)
    diurnal_ph = modules.solar_time_phase_corrector(diurnal_ph,24.,time_diff)
    diurnal_ph_spring = modules.solar_time_phase_corrector(diurnal_ph_spring,24.,time_diff)
    diurnal_ph_summer = modules.solar_time_phase_corrector(diurnal_ph_summer,24.,time_diff)
    diurnal_ph_autumn = modules.solar_time_phase_corrector(diurnal_ph_autumn,24.,time_diff)
    diurnal_ph_winter = modules.solar_time_phase_corrector(diurnal_ph_winter,24.,time_diff)
    diurnal_ave_waveform = np.append(diurnal_ave_waveform[-time_diff:],diurnal_ave_waveform[:-time_diff])
    diurnal_waveform_spring = np.append(diurnal_waveform_spring[-time_diff:],diurnal_waveform_spring[:-time_diff])
    diurnal_waveform_summer = np.append(diurnal_waveform_summer[-time_diff:],diurnal_waveform_summer[:-time_diff])
    diurnal_waveform_autumn = np.append(diurnal_waveform_autumn[-time_diff:],diurnal_waveform_autumn[:-time_diff])
    diurnal_waveform_winter = np.append(diurnal_waveform_winter[-time_diff:],diurnal_waveform_winter[:-time_diff])
    diurnal_ave_waveform_extended = np.append(diurnal_ave_waveform_extended[-time_diff:],diurnal_ave_waveform_extended[:-time_diff])

    #convert phase from radians to time
    #make sure diurnal phases are on hour by rounding - some numerical rounding issues previous
    diurnal_ph = round(modules.radians_to_time(diurnal_ph,24.))
    diurnal_ph_spring = round(modules.radians_to_time(diurnal_ph_spring,24.))
    diurnal_ph_summer = round(modules.radians_to_time(diurnal_ph_summer,24.))
    diurnal_ph_autumn = round(modules.radians_to_time(diurnal_ph_autumn,24.))
    diurnal_ph_winter = round(modules.radians_to_time(diurnal_ph_winter,24.))
    seasonal_ph = modules.radians_to_time(seasonal_ph,12.)
    seasonal_ph_day = modules.radians_to_time(seasonal_ph_day,12.)
    seasonal_ph_night = modules.radians_to_time(seasonal_ph_night,12.)
    
    #make extended waveform for diurnal cycle putting together diurnal waveforms from seasons
    #also make full waveforms combining both seasonal varying diurnal cycles: full_season_waveform & average diurnal cycles: full_ave_waveform
    #take mean from seasonal_waveform_extended
    seasonal_waveform_extended_nomean = seasonal_waveform_extended - mean
    if timeres == 'H':
        diurnal_season_waveform_extended = []
        full_ave_waveform = []
        full_season_waveform = []
        n_days = len(orig_datetimes)/24
        start_n = 0
        end_n = 24
        for i in range(n_days):
            if (int(orig_datetimes[start_n].month) == 3) or (int(orig_datetimes[start_n].month) == 4) or (int(orig_datetimes[start_n].month) == 5):
                if site_lat >= 0:
                    diurnal_season_waveform_extended=np.append(diurnal_season_waveform_extended,diurnal_waveform_spring)
                    full_ave_waveform=np.append(full_ave_waveform,((diurnal_ave_waveform-mean)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean)
                    full_season_waveform=np.append(full_season_waveform,((diurnal_waveform_spring-mean_spring)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean_spring)
                else:
                    diurnal_season_waveform_extended=np.append(diurnal_season_waveform_extended,diurnal_waveform_autumn)
                    full_ave_waveform=np.append(full_ave_waveform,((diurnal_ave_waveform-mean)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean)
                    full_season_waveform=np.append(full_season_waveform,((diurnal_waveform_autumn-mean_autumn)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean_autumn)
            elif (int(orig_datetimes[start_n].month) == 6) or (int(orig_datetimes[start_n].month) == 7) or (int(orig_datetimes[start_n].month) == 8):
                if site_lat >= 0:
                    diurnal_season_waveform_extended=np.append(diurnal_season_waveform_extended,diurnal_waveform_summer)
                    full_ave_waveform=np.append(full_ave_waveform,((diurnal_ave_waveform-mean)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean)
                    full_season_waveform=np.append(full_season_waveform,((diurnal_waveform_summer-mean_summer)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean_summer)
                else:
                    diurnal_season_waveform_extended=np.append(diurnal_season_waveform_extended,diurnal_waveform_winter)
                    full_ave_waveform=np.append(full_ave_waveform,((diurnal_ave_waveform-mean)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean)
                    full_season_waveform=np.append(full_season_waveform,((diurnal_waveform_winter-mean_winter)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean_winter)
            elif (int(orig_datetimes[start_n].month) == 9) or (int(orig_datetimes[start_n].month) == 10) or (int(orig_datetimes[start_n].month) == 11):
                if site_lat >= 0:
                    diurnal_season_waveform_extended=np.append(diurnal_season_waveform_extended,diurnal_waveform_autumn)
                    full_ave_waveform=np.append(full_ave_waveform,((diurnal_ave_waveform-mean)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean)
                    full_season_waveform=np.append(full_season_waveform,((diurnal_waveform_autumn-mean_autumn)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean_autumn)
                else:
                    diurnal_season_waveform_extended=np.append(diurnal_season_waveform_extended,diurnal_waveform_spring)
                    full_ave_waveform=np.append(full_ave_waveform,((diurnal_ave_waveform-mean)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean)
                    full_season_waveform=np.append(full_season_waveform,((diurnal_waveform_spring-mean_spring)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean_spring)
            elif (int(orig_datetimes[start_n].month) == 12) or (int(orig_datetimes[start_n].month) == 1) or (int(orig_datetimes[start_n].month) == 2):
                if site_lat >= 0:
                    diurnal_season_waveform_extended=np.append(diurnal_season_waveform_extended,diurnal_waveform_winter)
                    full_ave_waveform=np.append(full_ave_waveform,((diurnal_ave_waveform-mean)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean)
                    full_season_waveform=np.append(full_season_waveform,((diurnal_waveform_winter-mean_winter)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean_winter)
                else:
                    diurnal_season_waveform_extended=np.append(diurnal_season_waveform_extended,diurnal_waveform_summer)
                    full_ave_waveform=np.append(full_ave_waveform,((diurnal_ave_waveform-mean)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean)
                    full_season_waveform=np.append(full_season_waveform,((diurnal_waveform_summer-mean_summer)+seasonal_waveform_extended_nomean[start_n:end_n]))#+mean_summer)        
            start_n+=24
            end_n+=24
        
        full_ave_waveform = full_ave_waveform + mean
        full_season_waveform = full_season_waveform + mean
        
    elif (timeres == 'D') or (timeres == 'M'):
        diurnal_season_waveform_extended = diurnal_ave_waveform_extended
        full_ave_waveform = seasonal_waveform_extended
        full_season_waveform = seasonal_waveform_extended
    
    #get % variability of periodic waveform with raw time series
    #also get total variance and noise variance %
    pc_var_daily,na,na = modules.periodic_variance_percent(vals,full_ave_waveform,mag,mag[:12],valid)
    pc_var_seasonal,na,na = modules.periodic_variance_percent(vals,full_ave_waveform,mag,mag[12:],valid)
    pc_var_full,pc_var_noise,total_var = modules.periodic_variance_percent(vals,full_ave_waveform,mag,mag,valid)

    #remove mean from standard diurnal and seasonal waveforms
    #diurnal_ave_waveform = diurnal_ave_waveform - mean
    #seasonal_waveform = seasonal_waveform - mean
    #diurnal_waveform_spring = diurnal_waveform_spring - mean_spring
    #diurnal_waveform_summer = diurnal_waveform_summer - mean_summer
    #diurnal_waveform_autummn = diurnal_waveform_autumn - mean_autumn
    #diurnal_waveform_winter = diurnal_waveform_winter - mean_winter
    #seasonal_waveform_day = seasonal_waveform_day - mean_day
    #seasonal_waveform_night = seasonal_waveform_night - mean_night

    return (x,diurnal_mag,diurnal_ph,seasonal_mag,seasonal_ph,mean,p1,p5,p25,p50,p75,p95,p99,diurnal_ave_waveform,seasonal_waveform,full_ave_waveform,pc_var_daily,pc_var_seasonal,pc_var_full,pc_var_noise,total_var,
              diurnal_mag_spring,diurnal_ph_spring,mean_spring,p1_spring,p5_spring,p25_spring,p50_spring,p75_spring,p95_spring,p99_spring,diurnal_waveform_spring,
              diurnal_mag_summer,diurnal_ph_summer,mean_summer,p1_summer,p5_summer,p25_summer,p50_summer,p75_summer,p95_summer,p99_summer,diurnal_waveform_summer,
              diurnal_mag_autumn,diurnal_ph_autumn,mean_autumn,p1_autumn,p5_autumn,p25_autumn,p50_autumn,p75_autumn,p95_autumn,p99_autumn,diurnal_waveform_autumn,
              diurnal_mag_winter,diurnal_ph_winter,mean_winter,p1_winter,p5_winter,p25_winter,p50_winter,p75_winter,p95_winter,p99_winter,diurnal_waveform_winter,
              seasonal_mag_day,seasonal_ph_day,mean_day,p1_day,p5_day,p25_day,p50_day,p75_day,p95_day,p99_day,seasonal_waveform_day,
              seasonal_mag_night,seasonal_ph_night,mean_night,p1_night,p5_night,p25_night,p50_night,p75_night,p95_night,p99_night,seasonal_waveform_night,
              daily_h11_mag,daily_h10_mag,daily_h9_mag,daily_h8_mag,daily_h7_mag,daily_h6_mag,daily_h5_mag,daily_h4_mag,daily_h3_mag,daily_h2_mag,daily_h1_mag,orig_daily_mag,seasonal_h3_mag,seasonal_h2_mag,seasonal_h1_mag,annual_mag)
    
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

lat_indices = []
lon_indices = []
linear_data = []    

lat_i = 0
lon_i = 0

print gridbox_count

for siten in range(gridbox_count):
    linear_data.append(model_std_var[:,lat_i,lon_i])
   
    lat_indices.append(lat_i)
    lon_indices.append(lon_i)
   
    if lon_i == (len(lon_c)-1):
        lat_i+=1
        lon_i=0
    else:
        lon_i+=1

#linear_data = linear_data[:5]
#lat_indices = lat_indices[:5]
#lon_indices = lon_indices[:5]
#lat_c = [lat_c[0]]
#lon_c = lon_c[:5]

#for x in range(len(linear_data)):
#    big = run_LSP(linear_data[x],x)

if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=16)
    results = [pool.apply_async(run_LSP, (linear_data[x],x)) for x in range(len(linear_data))]
    big_array = [r.get() for r in results]
    pool.terminate()

indices_array = []

diurnal_mag = []
diurnal_ph = []
seasonal_mag = []
seasonal_ph = []
mean = []
p1 = []
p5 = []
p25 = []
p50 = []
p75 = []
p95 = []
p99 = []
diurnal_ave_waveform = []
#diurnal_ave_waveform_extended = []
#diurnal_season_waveform_extended = []
seasonal_waveform = []
#seasonal_waveform_extended = []
full_ave_waveform = []
#full_season_waveform = []
pc_var_daily = []
pc_var_seasonal = []
pc_var_full = []
pc_var_noise = []
total_var = []

diurnal_mag_spring = []
diurnal_ph_spring = []
mean_spring = []
p1_spring = []
p5_spring = []
p25_spring = []
p50_spring = []
p75_spring = []
p95_spring = []
p99_spring = []
diurnal_waveform_spring = []

diurnal_mag_summer = []
diurnal_ph_summer = []
mean_summer = []
p1_summer = []
p5_summer = []
p25_summer = []
p50_summer = []
p75_summer = []
p95_summer = []
p99_summer = []
diurnal_waveform_summer = []

diurnal_mag_autumn = []
diurnal_ph_autumn = []
mean_autumn = []
p1_autumn = []
p5_autumn = []
p25_autumn = []
p50_autumn = []
p75_autumn = []
p95_autumn = []
p99_autumn = []
diurnal_waveform_autumn = []

diurnal_mag_winter = []
diurnal_ph_winter = []
mean_winter = []
p1_winter = []
p5_winter = []
p25_winter = []
p50_winter = []
p75_winter = []
p95_winter = []
p99_winter = []
diurnal_waveform_winter = []

seasonal_mag_day = []
seasonal_ph_day = []
mean_day = []
p1_day = []
p5_day = []
p25_day = []
p50_day = []
p75_day = []
p95_day = []
p99_day = []
seasonal_waveform_day = []

seasonal_mag_night = []
seasonal_ph_night = []
mean_night = []
p1_night = []
p5_night = []
p25_night = []
p50_night = []
p75_night = []
p95_night = []
p99_night = []
seasonal_waveform_night = []

daily_h11_mag = []
daily_h10_mag = []
daily_h9_mag = []
daily_h8_mag = []
daily_h7_mag = []
daily_h6_mag = []
daily_h5_mag = []
daily_h4_mag = []
daily_h3_mag = []
daily_h2_mag = []
daily_h1_mag = []
daily_mag = []

annual_h3_mag = []
annual_h2_mag = []
annual_h1_mag = []
annual_mag = []

all_output = [diurnal_mag,diurnal_ph,seasonal_mag,seasonal_ph,mean,p1,p5,p25,p50,p75,p95,p99,diurnal_ave_waveform,seasonal_waveform,full_ave_waveform,pc_var_daily,pc_var_seasonal,pc_var_full,pc_var_noise,total_var,
              diurnal_mag_spring,diurnal_ph_spring,mean_spring,p1_spring,p5_spring,p25_spring,p50_spring,p75_spring,p95_spring,p99_spring,diurnal_waveform_spring,
              diurnal_mag_summer,diurnal_ph_summer,mean_summer,p1_summer,p5_summer,p25_summer,p50_summer,p75_summer,p95_summer,p99_summer,diurnal_waveform_summer,
              diurnal_mag_autumn,diurnal_ph_autumn,mean_autumn,p1_autumn,p5_autumn,p25_autumn,p50_autumn,p75_autumn,p95_autumn,p99_autumn,diurnal_waveform_autumn,
              diurnal_mag_winter,diurnal_ph_winter,mean_winter,p1_winter,p5_winter,p25_winter,p50_winter,p75_winter,p95_winter,p99_winter,diurnal_waveform_winter,
              seasonal_mag_day,seasonal_ph_day,mean_day,p1_day,p5_day,p25_day,p50_day,p75_day,p95_day,p99_day,seasonal_waveform_day,
              seasonal_mag_night,seasonal_ph_night,mean_night,p1_night,p5_night,p25_night,p50_night,p75_night,p95_night,p99_night,seasonal_waveform_night,
              daily_h11_mag,daily_h10_mag,daily_h9_mag,daily_h8_mag,daily_h7_mag,daily_h6_mag,daily_h5_mag,daily_h4_mag,daily_h3_mag,daily_h2_mag,daily_h1_mag,daily_mag,annual_h3_mag,annual_h2_mag,annual_h1_mag,annual_mag]

all_output_str = ['diurnal_mag','diurnal_ph','seasonal_mag','seasonal_ph','mean','p1','p5','p25','p50','p75','p95','p99','diurnal_ave_waveform','seasonal_waveform','full_ave_waveform','pc_var_daily','pc_var_seasonal','pc_var_full','pc_var_noise','total_var',
              'diurnal_mag_spring','diurnal_ph_spring','mean_spring','p1_spring','p5_spring','p25_spring','p50_spring','p75_spring','p95_spring','p99_spring','diurnal_waveform_spring',
              'diurnal_mag_summer','diurnal_ph_summer','mean_summer','p1_summer','p5_summer','p25_summer','p50_summer','p75_summer','p95_summer','p99_summer','diurnal_waveform_summer',
              'diurnal_mag_autumn','diurnal_ph_autumn','mean_autumn','p1_autumn','p5_autumn','p25_autumn','p50_autumn','p75_autumn','p95_autumn','p99_autumn','diurnal_waveform_autumn',
              'diurnal_mag_winter','diurnal_ph_winter','mean_winter','p1_winter','p5_winter','p25_winter','p50_winter','p75_winter','p95_winter','p99_winter','diurnal_waveform_winter',
              'seasonal_mag_day','seasonal_ph_day','mean_day','p1_day','p5_day','p25_day','p50_day','p75_day','p95_day','p99_day','seasonal_waveform_day',
              'seasonal_mag_night','seasonal_ph_night','mean_night','p1_night','p5_night','p25_night','p50_night','p75_night','p95_night','p99_night','seasonal_waveform_night',
              'daily_h11_mag','daily_h10_mag','daily_h9_mag','daily_h8_mag','daily_h7_mag','daily_h6_mag','daily_h5_mag','daily_h4_mag','daily_h3_mag','daily_h2_mag','daily_h1_mag','daily_mag','annual_h3_mag','annual_h2_mag','annual_h1_mag','annual_mag']


#write data into arrays
for i in range(len(big_array)):
    cut = big_array[i]
    
    indices_array.append(cut[0])
    
    for j in range(len(all_output)):
        all_output[j].append(cut[j+1])

#convert all to numpy arrays
indices_array = np.array(indices_array)
for i in range(len(all_output)):
    all_output[i] = np.array(all_output[i])

#sort arrays by indices array for sanity
for i in range(len(all_output)):
    all_output[i] = all_output[i][indices_array]
    
#reshape arrays (those which need to be lat x lon) into model grid
for i in range(len(all_output)):
    if all_output_str[i] not in ['diurnal_ave_waveform','diurnal_ave_waveform_extended','diurnal_season_waveform_extended','diurnal_waveform_spring','diurnal_waveform_summer','diurnal_waveform_autumn','diurnal_waveform_winter','seasonal_waveform','seasonal_waveform_extended','seasonal_waveform_day','seasonal_waveform_night','full_ave_waveform','full_season_waveform']:
        all_output[i] = np.reshape(all_output[i],(len(lat_c),len(lon_c)))

#reshape waveform arrays into model grid
d_rs = np.empty((24,len(lat_c),len(lon_c)))
d_rs_spring = np.empty((24,len(lat_c),len(lon_c)))
d_rs_summer = np.empty((24,len(lat_c),len(lon_c)))
d_rs_autumn = np.empty((24,len(lat_c),len(lon_c)))
d_rs_winter = np.empty((24,len(lat_c),len(lon_c)))
if timeres == 'H':
    s_rs = np.empty((8766,len(lat_c),len(lon_c)))
    s_rs_day = np.empty((8766,len(lat_c),len(lon_c)))
    s_rs_night = np.empty((8766,len(lat_c),len(lon_c)))
elif timeres == 'D':
    s_rs = np.empty((365,len(lat_c),len(lon_c)))
    s_rs_day = np.empty((365,len(lat_c),len(lon_c)))
    s_rs_night = np.empty((365,len(lat_c),len(lon_c)))
elif timeres == 'M':
    s_rs = np.empty((12,len(lat_c),len(lon_c)))
    s_rs_day = np.empty((12,len(lat_c),len(lon_c)))
    s_rs_night = np.empty((12,len(lat_c),len(lon_c)))
fw_len = len(full_ave_waveform[0])
a_rs = np.empty((fw_len,len(lat_c),len(lon_c)))
#as_rs = np.empty((fw_len,len(lat_c),len(lon_c)))
#d_e_rs = np.empty((fw_len,len(lat_c),len(lon_c)))
#ds_e_rs = np.empty((fw_len,len(lat_c),len(lon_c)))
#s_e_rs = np.empty((fw_len,len(lat_c),len(lon_c)))

lat_i = 0
lon_i = 0
for i in range(len(diurnal_ave_waveform)):
    current_d = diurnal_ave_waveform[i]
    current_d_spring = diurnal_waveform_spring[i]
    current_d_summer = diurnal_waveform_summer[i]
    current_d_autumn = diurnal_waveform_autumn[i]
    current_d_winter = diurnal_waveform_winter[i]
    current_s = seasonal_waveform[i]
    current_s_day = seasonal_waveform_day[i]
    current_s_night = seasonal_waveform_night[i]
    current_a = full_ave_waveform[i]
    #current_as = full_season_waveform[i]
    #current_d_e = diurnal_ave_waveform_extended[i]
    #current_ds_e = diurnal_season_waveform_extended[i]
    #current_s_e = seasonal_waveform_extended[i]
    
    d_rs[:,lat_i,lon_i] = current_d
    d_rs_spring[:,lat_i,lon_i] = current_d_spring
    d_rs_summer[:,lat_i,lon_i] = current_d_summer
    d_rs_autumn[:,lat_i,lon_i] = current_d_autumn
    d_rs_winter[:,lat_i,lon_i] = current_d_winter
    s_rs[:,lat_i,lon_i] = current_s
    s_rs_day[:,lat_i,lon_i] = current_s_day
    s_rs_night[:,lat_i,lon_i] = current_s_night
    a_rs[:,lat_i,lon_i] = current_a
    #as_rs[:,lat_i,lon_i] = current_as
    #d_e_rs[:,lat_i,lon_i] = current_d_e
    #ds_e_rs[:,lat_i,lon_i] = current_ds_e
    #s_e_rs[:,lat_i,lon_i] = current_s_e
    
    lon_i+=1
    
    if lon_i == len(lon_c):
        lat_i+=1
        lon_i = 0
        
diurnal_ave_waveform = np.copy(d_rs)
diurnal_waveform_spring = np.copy(d_rs_spring)
diurnal_waveform_summer = np.copy(d_rs_summer)
diurnal_waveform_autumn = np.copy(d_rs_autumn)
diurnal_waveform_winter = np.copy(d_rs_winter)
seasonal_waveform = np.copy(s_rs)
seasonal_waveform_day = np.copy(s_rs_day)
seasonal_waveform_night = np.copy(s_rs_night)
full_ave_waveform = np.copy(a_rs)
#full_season_waveform = np.copy(as_rs)
#diurnal_ave_waveform_extended = np.copy(d_e_rs)
#diurnal_season_waveform_extended = np.copy(ds_e_rs)
#seasonal_waveform_extended = np.copy(s_e_rs)

#save out sig periods to netcdf
root_grp_period = Dataset('LSP_stats.nc', 'w')

root_grp_period.createDimension('diurnal', 24)
if timeres == 'H':
    root_grp_period.createDimension('seasonal', 8766)
elif timeres == 'D':
    root_grp_period.createDimension('seasonal', 365)
elif timeres == 'M':    
    root_grp_period.createDimension('seasonal', 12)     
root_grp_period.createDimension('full', fw_len)
root_grp_period.createDimension('lat_centre', len(lat_c))
root_grp_period.createDimension('lon_centre', len(lon_c))
root_grp_period.createDimension('lat_edges', len(lat_e))
root_grp_period.createDimension('lon_edges', len(lon_e))

#set variables
d_amp = root_grp_period.createVariable('diurnal_amplitude', 'f4', ('lat_centre','lon_centre'))
d_ph = root_grp_period.createVariable('diurnal_phase', 'f4', ('lat_centre','lon_centre'))
s_amp = root_grp_period.createVariable('seasonal_amplitude', 'f4', ('lat_centre','lon_centre'))
s_ph = root_grp_period.createVariable('seasonal_phase', 'f4', ('lat_centre','lon_centre')) 
ave = root_grp_period.createVariable('average', 'f4', ('lat_centre','lon_centre')) 
per1 = root_grp_period.createVariable('percentile1', 'f4', ('lat_centre','lon_centre')) 
per5 = root_grp_period.createVariable('percentile5', 'f4', ('lat_centre','lon_centre')) 
per25 = root_grp_period.createVariable('percentile25', 'f4', ('lat_centre','lon_centre')) 
per50 = root_grp_period.createVariable('median', 'f4', ('lat_centre','lon_centre')) 
per75 = root_grp_period.createVariable('percentile75', 'f4', ('lat_centre','lon_centre')) 
per95 = root_grp_period.createVariable('percentile95', 'f4', ('lat_centre','lon_centre')) 
per99 = root_grp_period.createVariable('percentile99', 'f4', ('lat_centre','lon_centre')) 
dw = root_grp_period.createVariable('diurnal_average_waveform', 'f4', ('diurnal','lat_centre','lon_centre')) 
#dw_e = root_grp_period.createVariable('diurnal_average_waveform_extended', 'f4', ('full','lat_centre','lon_centre')) 
#dws_e = root_grp_period.createVariable('diurnal_season_waveform_extended', 'f4', ('full','lat_centre','lon_centre')) 
sw = root_grp_period.createVariable('seasonal_waveform', 'f4', ('seasonal','lat_centre','lon_centre'))
#sw_e = root_grp_period.createVariable('seasonal_waveform_extended', 'f4', ('full','lat_centre','lon_centre')) 
fw = root_grp_period.createVariable('full_average_waveform', 'f4', ('full','lat_centre','lon_centre'))
#fws = root_grp_period.createVariable('full_season_waveform', 'f4', ('full','lat_centre','lon_centre'))
pv_d = root_grp_period.createVariable('diurnal_periodic_variability', 'f4', ('lat_centre','lon_centre')) 
pv_s = root_grp_period.createVariable('seasonal_periodic_variability', 'f4', ('lat_centre','lon_centre')) 
pv_f = root_grp_period.createVariable('full_periodic_variability', 'f4', ('lat_centre','lon_centre')) 
pv_n = root_grp_period.createVariable('noise_periodic_variability', 'f4', ('lat_centre','lon_centre')) 
vt = root_grp_period.createVariable('total_variability', 'f4', ('lat_centre','lon_centre'))

d_amp_sp = root_grp_period.createVariable('diurnal_amplitude_spring', 'f4', ('lat_centre','lon_centre'))
d_ph_sp = root_grp_period.createVariable('diurnal_phase_spring', 'f4', ('lat_centre','lon_centre'))
ave_sp = root_grp_period.createVariable('average_spring', 'f4', ('lat_centre','lon_centre')) 
per1_sp = root_grp_period.createVariable('percentile1_spring', 'f4', ('lat_centre','lon_centre')) 
per5_sp = root_grp_period.createVariable('percentile5_spring', 'f4', ('lat_centre','lon_centre')) 
per25_sp = root_grp_period.createVariable('percentile25_spring', 'f4', ('lat_centre','lon_centre')) 
per50_sp = root_grp_period.createVariable('median_spring', 'f4', ('lat_centre','lon_centre')) 
per75_sp = root_grp_period.createVariable('percentile75_spring', 'f4', ('lat_centre','lon_centre')) 
per95_sp = root_grp_period.createVariable('percentile95_spring', 'f4', ('lat_centre','lon_centre')) 
per99_sp = root_grp_period.createVariable('percentile99_spring', 'f4', ('lat_centre','lon_centre')) 
dw_sp = root_grp_period.createVariable('diurnal_waveform_spring', 'f4', ('diurnal','lat_centre','lon_centre')) 

d_amp_su = root_grp_period.createVariable('diurnal_amplitude_summer', 'f4', ('lat_centre','lon_centre'))
d_ph_su = root_grp_period.createVariable('diurnal_phase_summer', 'f4', ('lat_centre','lon_centre'))
ave_su = root_grp_period.createVariable('average_summer', 'f4', ('lat_centre','lon_centre')) 
per1_su = root_grp_period.createVariable('percentile1_summer', 'f4', ('lat_centre','lon_centre')) 
per5_su = root_grp_period.createVariable('percentile5_summer', 'f4', ('lat_centre','lon_centre')) 
per25_su = root_grp_period.createVariable('percentile25_summer', 'f4', ('lat_centre','lon_centre')) 
per50_su = root_grp_period.createVariable('median_summer', 'f4', ('lat_centre','lon_centre')) 
per75_su = root_grp_period.createVariable('percentile75_summer', 'f4', ('lat_centre','lon_centre')) 
per95_su = root_grp_period.createVariable('percentile95_summer', 'f4', ('lat_centre','lon_centre')) 
per99_su = root_grp_period.createVariable('percentile99_summer', 'f4', ('lat_centre','lon_centre')) 
dw_su = root_grp_period.createVariable('diurnal_waveform_summer', 'f4', ('diurnal','lat_centre','lon_centre')) 

d_amp_au = root_grp_period.createVariable('diurnal_amplitude_autumn', 'f4', ('lat_centre','lon_centre'))
d_ph_au = root_grp_period.createVariable('diurnal_phase_autumn', 'f4', ('lat_centre','lon_centre'))
ave_au = root_grp_period.createVariable('average_autumn', 'f4', ('lat_centre','lon_centre')) 
per1_au = root_grp_period.createVariable('percentile1_autumn', 'f4', ('lat_centre','lon_centre')) 
per5_au = root_grp_period.createVariable('percentile5_autumn', 'f4', ('lat_centre','lon_centre')) 
per25_au = root_grp_period.createVariable('percentile25_autumn', 'f4', ('lat_centre','lon_centre')) 
per50_au = root_grp_period.createVariable('median_autumn', 'f4', ('lat_centre','lon_centre')) 
per75_au = root_grp_period.createVariable('percentile75_autumn', 'f4', ('lat_centre','lon_centre')) 
per95_au = root_grp_period.createVariable('percentile95_autumn', 'f4', ('lat_centre','lon_centre')) 
per99_au = root_grp_period.createVariable('percentile99_autumn', 'f4', ('lat_centre','lon_centre')) 
dw_au = root_grp_period.createVariable('diurnal_waveform_autumn', 'f4', ('diurnal','lat_centre','lon_centre')) 

d_amp_wi = root_grp_period.createVariable('diurnal_amplitude_winter', 'f4', ('lat_centre','lon_centre'))
d_ph_wi = root_grp_period.createVariable('diurnal_phase_winter', 'f4', ('lat_centre','lon_centre'))
ave_wi = root_grp_period.createVariable('average_winter', 'f4', ('lat_centre','lon_centre')) 
per1_wi = root_grp_period.createVariable('percentile1_winter', 'f4', ('lat_centre','lon_centre')) 
per5_wi = root_grp_period.createVariable('percentile5_winter', 'f4', ('lat_centre','lon_centre')) 
per25_wi = root_grp_period.createVariable('percentile25_winter', 'f4', ('lat_centre','lon_centre')) 
per50_wi = root_grp_period.createVariable('median_winter', 'f4', ('lat_centre','lon_centre')) 
per75_wi = root_grp_period.createVariable('percentile75_winter', 'f4', ('lat_centre','lon_centre')) 
per95_wi = root_grp_period.createVariable('percentile95_winter', 'f4', ('lat_centre','lon_centre')) 
per99_wi = root_grp_period.createVariable('percentile99_winter', 'f4', ('lat_centre','lon_centre')) 
dw_wi = root_grp_period.createVariable('diurnal_waveform_winter', 'f4', ('diurnal','lat_centre','lon_centre')) 
   
s_amp_da = root_grp_period.createVariable('seasonal_amplitude_day', 'f4', ('lat_centre','lon_centre'))
s_ph_da = root_grp_period.createVariable('seasonal_phase_day', 'f4', ('lat_centre','lon_centre'))
ave_da = root_grp_period.createVariable('average_day', 'f4', ('lat_centre','lon_centre')) 
per1_da = root_grp_period.createVariable('percentile1_day', 'f4', ('lat_centre','lon_centre')) 
per5_da = root_grp_period.createVariable('percentile5_day', 'f4', ('lat_centre','lon_centre')) 
per25_da = root_grp_period.createVariable('percentile25_day', 'f4', ('lat_centre','lon_centre')) 
per50_da = root_grp_period.createVariable('median_day', 'f4', ('lat_centre','lon_centre')) 
per75_da = root_grp_period.createVariable('percentile75_day', 'f4', ('lat_centre','lon_centre')) 
per95_da = root_grp_period.createVariable('percentile95_day', 'f4', ('lat_centre','lon_centre')) 
per99_da = root_grp_period.createVariable('percentile99_day', 'f4', ('lat_centre','lon_centre')) 
sw_da = root_grp_period.createVariable('seasonal_waveform_day', 'f4', ('seasonal','lat_centre','lon_centre')) 

s_amp_ni = root_grp_period.createVariable('seasonal_amplitude_night', 'f4', ('lat_centre','lon_centre'))
s_ph_ni = root_grp_period.createVariable('seasonal_phase_night', 'f4', ('lat_centre','lon_centre'))
ave_ni = root_grp_period.createVariable('average_night', 'f4', ('lat_centre','lon_centre')) 
per1_ni = root_grp_period.createVariable('percentile1_night', 'f4', ('lat_centre','lon_centre')) 
per5_ni = root_grp_period.createVariable('percentile5_night', 'f4', ('lat_centre','lon_centre')) 
per25_ni = root_grp_period.createVariable('percentile25_night', 'f4', ('lat_centre','lon_centre')) 
per50_ni = root_grp_period.createVariable('median_night', 'f4', ('lat_centre','lon_centre')) 
per75_ni = root_grp_period.createVariable('percentile75_night', 'f4', ('lat_centre','lon_centre')) 
per95_ni = root_grp_period.createVariable('percentile95_night', 'f4', ('lat_centre','lon_centre')) 
per99_ni = root_grp_period.createVariable('percentile99_night', 'f4', ('lat_centre','lon_centre')) 
sw_ni = root_grp_period.createVariable('seasonal_waveform_night', 'f4', ('seasonal','lat_centre','lon_centre')) 

d_h12_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic12', 'f4', ('lat_centre','lon_centre'))
d_h11_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic11', 'f4', ('lat_centre','lon_centre'))
d_h10_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic10', 'f4', ('lat_centre','lon_centre'))
d_h9_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic9', 'f4', ('lat_centre','lon_centre'))
d_h8_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic8', 'f4', ('lat_centre','lon_centre'))
d_h7_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic7', 'f4', ('lat_centre','lon_centre'))
d_h6_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic6', 'f4', ('lat_centre','lon_centre'))
d_h5_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic5', 'f4', ('lat_centre','lon_centre'))
d_h4_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic4', 'f4', ('lat_centre','lon_centre'))
d_h3_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic3', 'f4', ('lat_centre','lon_centre'))
d_h2_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic2', 'f4', ('lat_centre','lon_centre'))
d_h1_amp = root_grp_period.createVariable('diurnal_amplitude_harmonic1', 'f4', ('lat_centre','lon_centre'))
#s_h12_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic12', 'f4', ('lat_centre','lon_centre'))
#s_h11_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic11', 'f4', ('lat_centre','lon_centre'))
#s_h10_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic10', 'f4', ('lat_centre','lon_centre'))
#s_h9_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic9', 'f4', ('lat_centre','lon_centre'))
#s_h8_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic8', 'f4', ('lat_centre','lon_centre'))
#s_h7_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic7', 'f4', ('lat_centre','lon_centre'))
#s_h6_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic6', 'f4', ('lat_centre','lon_centre'))
#s_h5_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic5', 'f4', ('lat_centre','lon_centre'))
s_h4_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic4', 'f4', ('lat_centre','lon_centre'))
s_h3_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic3', 'f4', ('lat_centre','lon_centre'))
s_h2_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic2', 'f4', ('lat_centre','lon_centre'))
s_h1_amp = root_grp_period.createVariable('seasonal_amplitude_harmonic1', 'f4', ('lat_centre','lon_centre'))

lat_centre = root_grp_period.createVariable('lat_centre', 'f4', ('lat_centre',))
lon_centre = root_grp_period.createVariable('lon_centre', 'f4', ('lon_centre',))
lat_edge = root_grp_period.createVariable('lat_edges', 'f4', ('lat_edges',))
lon_edge = root_grp_period.createVariable('lon_edges', 'f4', ('lon_edges',))

#save out
d_amp[:] = diurnal_mag
d_ph[:] = diurnal_ph
s_amp[:] = seasonal_mag
s_ph[:] = seasonal_ph
ave[:] = mean
per1[:] = p1
per5[:] = p5
per25[:] = p25
per50[:] = p50
per75[:] = p75
per95[:] =  p95
per99[:] = p99
dw[:] =  diurnal_ave_waveform
#dw_e[:] = diurnal_ave_waveform_extended
#dws_e[:] = diurnal_season_waveform_extended
sw[:] = seasonal_waveform
#sw_e[:] = seasonal_waveform_extended
fw[:] = full_ave_waveform
#fws[:] = full_season_waveform
pv_d[:] = pc_var_daily 
pv_s[:] = pc_var_seasonal
pv_f[:] = pc_var_full
pv_n[:] = pc_var_noise
vt[:] = total_var

d_amp_sp[:] = diurnal_mag_spring
d_ph_sp[:] = diurnal_ph_spring
ave_sp[:] = mean_spring
per1_sp[:] = p1_spring
per5_sp[:] = p5_spring
per25_sp[:] = p25_spring
per50_sp[:] = p50_spring
per75_sp[:] = p75_spring
per95_sp[:] = p95_spring
per99_sp[:] = p99_spring
dw_sp[:] = diurnal_waveform_spring

d_amp_su[:] = diurnal_mag_summer
d_ph_su[:] = diurnal_ph_summer
ave_su[:] = mean_summer
per1_su[:] = p1_summer 
per5_su[:] = p5_summer 
per25_su[:] = p25_summer 
per50_su[:] = p50_summer  
per75_su[:] = p75_summer 
per95_su[:] = p95_summer  
per99_su[:] = p99_summer 
dw_su[:] = diurnal_waveform_summer

d_amp_au[:] = diurnal_mag_autumn
d_ph_au[:] = diurnal_ph_autumn
ave_au[:] = mean_autumn
per1_au[:] = p1_autumn
per5_au[:] = p5_autumn
per25_au[:] = p25_autumn
per50_au[:] = p50_autumn
per75_au[:] = p75_autumn 
per95_au[:] = p95_autumn
per99_au[:] = p99_autumn
dw_au[:] = diurnal_waveform_autumn

d_amp_wi[:] = diurnal_mag_winter
d_ph_wi[:] = diurnal_ph_winter
ave_wi[:] = mean_winter
per1_wi[:] = p1_winter
per5_wi[:] = p5_winter
per25_wi[:] = p25_winter
per50_wi[:] = p50_winter
per75_wi[:] = p75_winter
per95_wi[:] = p95_winter
per99_wi[:] = p99_winter
dw_wi[:] = diurnal_waveform_winter
   
s_amp_da[:] = seasonal_mag_day
s_ph_da[:] = seasonal_ph_day
ave_da[:] = mean_day
per1_da[:] = p1_day
per5_da[:] = p5_day
per25_da[:] = p25_day
per50_da[:] = p50_day
per75_da[:] = p75_day
per95_da[:] = p95_day
per99_da[:] = p99_day
sw_da[:] = seasonal_waveform_day

s_amp_ni[:] = seasonal_mag_night
s_ph_ni[:] = seasonal_ph_night
ave_ni[:] = mean_night
per1_ni[:] = p1_night
per5_ni[:] = p5_night
per25_ni[:] = p25_night 
per50_ni[:] = p50_night
per75_ni[:] = p75_night 
per95_ni[:] = p95_night 
per99_ni[:] = p99_night
sw_ni[:] = seasonal_waveform_night

d_h12_amp[:] = daily_h11_mag
d_h11_amp[:] = daily_h10_mag
d_h10_amp[:] = daily_h9_mag
d_h9_amp[:] = daily_h8_mag
d_h8_amp[:] = daily_h7_mag
d_h7_amp[:] = daily_h6_mag
d_h6_amp[:] = daily_h5_mag
d_h5_amp[:] = daily_h4_mag
d_h4_amp[:] = daily_h3_mag
d_h3_amp[:] = daily_h2_mag
d_h2_amp[:] = daily_h1_mag
d_h1_amp[:] = daily_mag
#s_h12_amp[:] = annual_h11_mag
#s_h11_amp[:] = annual_h10_mag
#s_h10_amp[:] = annual_h9_mag
#s_h9_amp[:] = annual_h8_mag
#s_h8_amp[:] = annual_h7_mag
#s_h7_amp[:] = annual_h6_mag
#s_h6_amp[:] = annual_h5_mag
#s_h5_amp[:] = annual_h4_mag
s_h4_amp[:] = annual_h3_mag
s_h3_amp[:] = annual_h2_mag
s_h2_amp[:] = annual_h1_mag
s_h1_amp[:] = annual_mag

lat_centre[:] = lat_c
lon_centre[:] = lon_c
lat_edge[:] = lat_e     
lon_edge[:] = lon_e          
        
root_grp_period.close()

end = datetime.datetime.now() 
diff = end - start
print 'time = ', diff.seconds
