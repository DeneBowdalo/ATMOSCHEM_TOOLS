import modules
import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm
from cmath import *
import colorsys
from scipy import signal
from netCDF4 import Dataset
import datetime
import pandas as pd

models = ['CESMCAM', 'CMAM', 'GFDLAM3','GISSE2R', 'MIROCCHEM', 'MOCAGE', 'UMCAM']

for model in models:
    if (model == 'CESMCAM') or (model == 'CMAM') or (model == 'MIROCCHEM') or (model == 'MOCAGE'):
        year2000 = '2000_2003'
        year2100 = '2100_2103'

    if (model == 'GFDLAM3'):
        year2000 = '2001_2004'
        year2100 = '2101_2104'

    if (model == 'GISSE2R'):
        year2000 = '2000_2003'
        year2100 = '2102_2105'

    if (model == 'UMCAM'):
        year2000 = '2000_2003'
        year2100 = '2097_2100'        

    #read in 2000 period data
    f2000 = '/work/home/db876/grid/O3/2000_2012/%s/%s_SURFACE_*_*_*_H_*/LSP_stats.nc'%(year2000,model)

    #read in 2100 period data
    f2100 = '/work/home/db876/grid/O3/2095_2111/%s/%s_SURFACE_*_*_*_H_rcp85/LSP_stats.nc'%(year2100,model)

    #read in 2100 model fixed emissions period data
    f2100e = '/work/home/db876/grid/O3/2095_2111/%s/%s_SURFACE_*_*_*_H_rcp85em2000/LSP_stats.nc'%(year2100,model)

    root_2000 = Dataset(f2000)
    lat_e = root_2000.variables['lat_edges'][:]
    lon_e = root_2000.variables['lon_edges'][:]
    lat_c = root_2000.variables['lat_centre'][:]
    lon_c = root_2000.variables['lon_centre'][:]

    (diurnal_mag_2000,diurnal_ph_2000,seasonal_mag_2000,seasonal_ph_2000,mean_2000,p1_2000,p5_2000,p25_2000,p50_2000,p75_2000,p95_2000,p99_2000,diurnal_ave_waveform_2000,seasonal_waveform_2000,full_ave_waveform_2000,pc_var_daily_2000,pc_var_seasonal_2000,pc_var_full_2000,pc_var_noise_2000,total_var_2000,
    diurnal_mag_spring_2000,diurnal_ph_spring_2000,mean_spring_2000,p1_spring_2000,p5_spring_2000,p25_spring_2000,p50_spring_2000,p75_spring_2000,p95_spring_2000,p99_spring_2000,diurnal_waveform_spring_2000,
    diurnal_mag_summer_2000,diurnal_ph_summer_2000,mean_summer_2000,p1_summer_2000,p5_summer_2000,p25_summer_2000,p50_summer_2000,p75_summer_2000,p95_summer_2000,p99_summer_2000,diurnal_waveform_summer_2000,
    diurnal_mag_autumn_2000,diurnal_ph_autumn_2000,mean_autumn_2000,p1_autumn_2000,p5_autumn_2000,p25_autumn_2000,p50_autumn_2000,p75_autumn_2000,p95_autumn_2000,p99_autumn_2000,diurnal_waveform_autumn_2000,
    diurnal_mag_winter_2000,diurnal_ph_winter_2000,mean_winter_2000,p1_winter_2000,p5_winter_2000,p25_winter_2000,p50_winter_2000,p75_winter_2000,p95_winter_2000,p99_winter_2000,diurnal_waveform_winter_2000,
    seasonal_mag_day_2000,seasonal_ph_day_2000,mean_day_2000,p1_day_2000,p5_day_2000,p25_day_2000,p50_day_2000,p75_day_2000,p95_day_2000,p99_day_2000,seasonal_waveform_day_2000,
    seasonal_mag_night_2000,seasonal_ph_night_2000,mean_night_2000,p1_night_2000,p5_night_2000,p25_night_2000,p50_night_2000,p75_night_2000,p95_night_2000,p99_night_2000,seasonal_waveform_night_2000,
    daily_h3_mag_2000,daily_h2_mag_2000,daily_h1_mag_2000,daily_mag_2000,annual_h3_mag_2000,annual_h2_mag_2000,annual_h1_mag_2000,annual_mag_2000) = modules.get_periodic_specific('na',f2000)  

    if model != 'GISSE2R':
        (diurnal_mag_2100,diurnal_ph_2100,seasonal_mag_2100,seasonal_ph_2100,mean_2100,p1_2100,p5_2100,p25_2100,p50_2100,p75_2100,p95_2100,p99_2100,diurnal_ave_waveform_2100,seasonal_waveform_2100,full_ave_waveform_2100,pc_var_daily_2100,pc_var_seasonal_2100,pc_var_full_2100,pc_var_noise_2100,total_var_2100,
        diurnal_mag_spring_2100,diurnal_ph_spring_2100,mean_spring_2100,p1_spring_2100,p5_spring_2100,p25_spring_2100,p50_spring_2100,p75_spring_2100,p95_spring_2100,p99_spring_2100,diurnal_waveform_spring_2100,
        diurnal_mag_summer_2100,diurnal_ph_summer_2100,mean_summer_2100,p1_summer_2100,p5_summer_2100,p25_summer_2100,p50_summer_2100,p75_summer_2100,p95_summer_2100,p99_summer_2100,diurnal_waveform_summer_2100,
        diurnal_mag_autumn_2100,diurnal_ph_autumn_2100,mean_autumn_2100,p1_autumn_2100,p5_autumn_2100,p25_autumn_2100,p50_autumn_2100,p75_autumn_2100,p95_autumn_2100,p99_autumn_2100,diurnal_waveform_autumn_2100,
        diurnal_mag_winter_2100,diurnal_ph_winter_2100,mean_winter_2100,p1_winter_2100,p5_winter_2100,p25_winter_2100,p50_winter_2100,p75_winter_2100,p95_winter_2100,p99_winter_2100,diurnal_waveform_winter_2100,
        seasonal_mag_day_2100,seasonal_ph_day_2100,mean_day_2100,p1_day_2100,p5_day_2100,p25_day_2100,p50_day_2100,p75_day_2100,p95_day_2100,p99_day_2100,seasonal_waveform_day_2100,
        seasonal_mag_night_2100,seasonal_ph_night_2100,mean_night_2100,p1_night_2100,p5_night_2100,p25_night_2100,p50_night_2100,p75_night_2100,p95_night_2100,p99_night_2100,seasonal_waveform_night_2100,
        daily_h3_mag_2100,daily_h2_mag_2100,daily_h1_mag_2100,daily_mag_2100,annual_h3_mag_2100,annual_h2_mag_2100,annual_h1_mag_2100,annual_mag_2100) = modules.get_periodic_specific('na',f2100)  

    if model != 'CMAM':
        (diurnal_mag_2100e,diurnal_ph_2100e,seasonal_mag_2100e,seasonal_ph_2100e,mean_2100e,p1_2100e,p5_2100e,p25_2100e,p50_2100e,p75_2100e,p95_2100e,p99_2100e,diurnal_ave_waveform_2100e,seasonal_waveform_2100e,full_ave_waveform_2100e,pc_var_daily_2100e,pc_var_seasonal_2100e,pc_var_full_2100e,pc_var_noise_2100e,total_var_2100e,
        diurnal_mag_spring_2100e,diurnal_ph_spring_2100e,mean_spring_2100e,p1_spring_2100e,p5_spring_2100e,p25_spring_2100e,p50_spring_2100e,p75_spring_2100e,p95_spring_2100e,p99_spring_2100e,diurnal_waveform_spring_2100e,
        diurnal_mag_summer_2100e,diurnal_ph_summer_2100e,mean_summer_2100e,p1_summer_2100e,p5_summer_2100e,p25_summer_2100e,p50_summer_2100e,p75_summer_2100e,p95_summer_2100e,p99_summer_2100e,diurnal_waveform_summer_2100e,
        diurnal_mag_autumn_2100e,diurnal_ph_autumn_2100e,mean_autumn_2100e,p1_autumn_2100e,p5_autumn_2100e,p25_autumn_2100e,p50_autumn_2100e,p75_autumn_2100e,p95_autumn_2100e,p99_autumn_2100e,diurnal_waveform_autumn_2100e,
        diurnal_mag_winter_2100e,diurnal_ph_winter_2100e,mean_winter_2100e,p1_winter_2100e,p5_winter_2100e,p25_winter_2100e,p50_winter_2100e,p75_winter_2100e,p95_winter_2100e,p99_winter_2100e,diurnal_waveform_winter_2100e,
        seasonal_mag_day_2100e,seasonal_ph_day_2100e,mean_day_2100e,p1_day_2100e,p5_day_2100e,p25_day_2100e,p50_day_2100e,p75_day_2100e,p95_day_2100e,p99_day_2100e,seasonal_waveform_day_2100e,
        seasonal_mag_night_2100e,seasonal_ph_night_2100e,mean_night_2100e,p1_night_2100e,p5_night_2100e,p25_night_2100e,p50_night_2100e,p75_night_2100e,p95_night_2100e,p99_night_2100e,seasonal_waveform_night_2100e,
        daily_h3_mag_2100e,daily_h2_mag_2100e,daily_h1_mag_2100e,daily_mag_2100e,annual_h3_mag_2100e,annual_h2_mag_2100e,annual_h1_mag_2100e,annual_mag_2100e) = modules.get_periodic_specific('na',f2100e)  
    
    #params

    params_2000 = [diurnal_mag_2000,diurnal_ph_2000,seasonal_mag_2000,seasonal_ph_2000,mean_2000,p1_2000,p5_2000,p25_2000,p50_2000,p75_2000,p95_2000,p99_2000,pc_var_daily_2000,pc_var_seasonal_2000,pc_var_full_2000,pc_var_noise_2000,total_var_2000,
    diurnal_mag_spring_2000,diurnal_ph_spring_2000,mean_spring_2000,p1_spring_2000,p5_spring_2000,p25_spring_2000,p50_spring_2000,p75_spring_2000,p95_spring_2000,p99_spring_2000,
    diurnal_mag_summer_2000,diurnal_ph_summer_2000,mean_summer_2000,p1_summer_2000,p5_summer_2000,p25_summer_2000,p50_summer_2000,p75_summer_2000,p95_summer_2000,p99_summer_2000,
    diurnal_mag_autumn_2000,diurnal_ph_autumn_2000,mean_autumn_2000,p1_autumn_2000,p5_autumn_2000,p25_autumn_2000,p50_autumn_2000,p75_autumn_2000,p95_autumn_2000,p99_autumn_2000,
    diurnal_mag_winter_2000,diurnal_ph_winter_2000,mean_winter_2000,p1_winter_2000,p5_winter_2000,p25_winter_2000,p50_winter_2000,p75_winter_2000,p95_winter_2000,p99_winter_2000,
    daily_h3_mag_2000,daily_h2_mag_2000,daily_h1_mag_2000,daily_mag_2000,annual_h3_mag_2000,annual_h2_mag_2000,annual_h1_mag_2000,annual_mag_2000]

    if model != 'GISSE2R':
        params_2100 = [diurnal_mag_2100,diurnal_ph_2100,seasonal_mag_2100,seasonal_ph_2100,mean_2100,p1_2100,p5_2100,p25_2100,p50_2100,p75_2100,p95_2100,p99_2100,pc_var_daily_2100,pc_var_seasonal_2100,pc_var_full_2100,pc_var_noise_2100,total_var_2100,
        diurnal_mag_spring_2100,diurnal_ph_spring_2100,mean_spring_2100,p1_spring_2100,p5_spring_2100,p25_spring_2100,p50_spring_2100,p75_spring_2100,p95_spring_2100,p99_spring_2100,
        diurnal_mag_summer_2100,diurnal_ph_summer_2100,mean_summer_2100,p1_summer_2100,p5_summer_2100,p25_summer_2100,p50_summer_2100,p75_summer_2100,p95_summer_2100,p99_summer_2100,
        diurnal_mag_autumn_2100,diurnal_ph_autumn_2100,mean_autumn_2100,p1_autumn_2100,p5_autumn_2100,p25_autumn_2100,p50_autumn_2100,p75_autumn_2100,p95_autumn_2100,p99_autumn_2100,
        diurnal_mag_winter_2100,diurnal_ph_winter_2100,mean_winter_2100,p1_winter_2100,p5_winter_2100,p25_winter_2100,p50_winter_2100,p75_winter_2100,p95_winter_2100,p99_winter_2100,
        daily_h3_mag_2100,daily_h2_mag_2100,daily_h1_mag_2100,daily_mag_2100,annual_h3_mag_2100,annual_h2_mag_2100,annual_h1_mag_2100,annual_mag_2100]

    if model != 'CMAM':
        params_2100e = [diurnal_mag_2100e,diurnal_ph_2100e,seasonal_mag_2100e,seasonal_ph_2100e,mean_2100e,p1_2100e,p5_2100e,p25_2100e,p50_2100e,p75_2100e,p95_2100e,p99_2100e,pc_var_daily_2100e,pc_var_seasonal_2100e,pc_var_full_2100e,pc_var_noise_2100e,total_var_2100e,
        diurnal_mag_spring_2100e,diurnal_ph_spring_2100e,mean_spring_2100e,p1_spring_2100e,p5_spring_2100e,p25_spring_2100e,p50_spring_2100e,p75_spring_2100e,p95_spring_2100e,p99_spring_2100e,
        diurnal_mag_summer_2100e,diurnal_ph_summer_2100e,mean_summer_2100e,p1_summer_2100e,p5_summer_2100e,p25_summer_2100e,p50_summer_2100e,p75_summer_2100e,p95_summer_2100e,p99_summer_2100e,
        diurnal_mag_autumn_2100e,diurnal_ph_autumn_2100e,mean_autumn_2100e,p1_autumn_2100e,p5_autumn_2100e,p25_autumn_2100e,p50_autumn_2100e,p75_autumn_2100e,p95_autumn_2100e,p99_autumn_2100e,
        diurnal_mag_winter_2100e,diurnal_ph_winter_2100e,mean_winter_2100e,p1_winter_2100e,p5_winter_2100e,p25_winter_2100e,p50_winter_2100e,p75_winter_2100e,p95_winter_2100e,p99_winter_2100e,
        daily_h3_mag_2100e,daily_h2_mag_2100e,daily_h1_mag_2100e,daily_mag_2100e,annual_h3_mag_2100e,annual_h2_mag_2100e,annual_h1_mag_2100e,annual_mag_2100e]

    params_str_2000 = ['diurnal_mag_2000','diurnal_ph_2000','seasonal_mag_2000','seasonal_ph_2000','mean_2000','p1_2000','p5_2000','p25_2000','p50_2000','p75_2000','p95_2000','p99_2000','pc_var_daily_2000','pc_var_seasonal_2000','pc_var_full_2000','pc_var_noise_2000','total_var_2000',
    'diurnal_mag_spring_2000','diurnal_ph_spring_2000','mean_spring_2000','p1_spring_2000','p5_spring_2000','p25_spring_2000','p50_spring_2000','p75_spring_2000','p95_spring_2000','p99_spring_2000',
    'diurnal_mag_summer_2000','diurnal_ph_summer_2000','mean_summer_2000','p1_summer_2000','p5_summer_2000','p25_summer_2000','p50_summer_2000','p75_summer_2000','p95_summer_2000','p99_summer_2000',
    'diurnal_mag_autumn_2000','diurnal_ph_autumn_2000','mean_autumn_2000','p1_autumn_2000','p5_autumn_2000','p25_autumn_2000','p50_autumn_2000','p75_autumn_2000','p95_autumn_2000','p99_autumn_2000',
    'diurnal_mag_winter_2000','diurnal_ph_winter_2000','mean_winter_2000','p1_winter_2000','p5_winter_2000','p25_winter_2000','p50_winter_2000','p75_winter_2000','p95_winter_2000','p99_winter_2000',
    'daily_h3_mag_2000','daily_h2_mag_2000','daily_h1_mag_2000','daily_mag_2000','annual_h3_mag_2000','annual_h2_mag_2000','annual_h1_mag_2000','annual_mag_2000']

    params_str_2100 = ['diurnal_mag_2100','diurnal_ph_2100','seasonal_mag_2100','seasonal_ph_2100','mean_2100','p1_2100','p5_2100','p25_2100','p50_2100','p75_2100','p95_2100','p99_2100','pc_var_daily_2100','pc_var_seasonal_2100','pc_var_full_2100','pc_var_noise_2100','total_var_2100',
    'diurnal_mag_spring_2100','diurnal_ph_spring_2100','mean_spring_2100','p1_spring_2100','p5_spring_2100','p25_spring_2100','p50_spring_2100','p75_spring_2100','p95_spring_2100','p99_spring_2100',
    'diurnal_mag_summer_2100','diurnal_ph_summer_2100','mean_summer_2100','p1_summer_2100','p5_summer_2100','p25_summer_2100','p50_summer_2100','p75_summer_2100','p95_summer_2100','p99_summer_2100',
    'diurnal_mag_autumn_2100','diurnal_ph_autumn_2100','mean_autumn_2100','p1_autumn_2100','p5_autumn_2100','p25_autumn_2100','p50_autumn_2100','p75_autumn_2100','p95_autumn_2100','p99_autumn_2100',
    'diurnal_mag_winter_2100','diurnal_ph_winter_2100','mean_winter_2100','p1_winter_2100','p5_winter_2100','p25_winter_2100','p50_winter_2100','p75_winter_2100','p95_winter_2100','p99_winter_2100',
    'daily_h3_mag_2100','daily_h2_mag_2100','daily_h1_mag_2100','daily_mag_2100','annual_h3_mag_2100','annual_h2_mag_2100','annual_h1_mag_2100','annual_mag_2100']

    params_str_2100e = ['diurnal_mag_2100e','diurnal_ph_2100e','seasonal_mag_2100e','seasonal_ph_2100e','mean_2100e','p1_2100e','p5_2100e','p25_2100e','p50_2100e','p75_2100e','p95_2100e','p99_2100e','pc_var_daily_2100e','pc_var_seasonal_2100e','pc_var_full_2100e','pc_var_noise_2100e','total_var_2100e',
    'diurnal_mag_spring_2100e','diurnal_ph_spring_2100e','mean_spring_2100e','p1_spring_2100e','p5_spring_2100e','p25_spring_2100e','p50_spring_2100e','p75_spring_2100e','p95_spring_2100e','p99_spring_2100e',
    'diurnal_mag_summer_2100e','diurnal_ph_summer_2100e','mean_summer_2100e','p1_summer_2100','p5_summer_2100e','p25_summer_2100e','p50_summer_2100e','p75_summer_2100e','p95_summer_2100e','p99_summer_2100e',
    'diurnal_mag_autumn_2100e','diurnal_ph_autumn_2100e','mean_autumn_2100e','p1_autumn_2100','p5_autumn_2100e','p25_autumn_2100e','p50_autumn_2100e','p75_autumn_2100e','p95_autumn_2100e','p99_autumn_2100e',
    'diurnal_mag_winter_2100e','diurnal_ph_winter_2100e','mean_winter_2100e','p1_winter_2100','p5_winter_2100e','p25_winter_2100e','p50_winter_2100e','p75_winter_2100e','p95_winter_2100e','p99_winter_2100e',
    'daily_h3_mag_2100e','daily_h2_mag_2100e','daily_h1_mag_2100e','daily_mag_2100e','annual_h3_mag_2100e','annual_h2_mag_2100e','annual_h1_mag_2100e','annual_mag_2100e']


    for p in range(len(params_2000)):
        param_2000 = params_2000[p]
        param_str_2000 = params_str_2000[p]
        if model != 'GISSE2R':
            param_2100 = params_2100[p]
            param_str_2100 = params_str_2100[p]
        if model != 'CMAM':
            param_2100e = params_2100e[p]
            param_str_2100e = params_str_2100e[p]

        #set up plot
        fig =plt.figure(figsize=(17,9.5))
        fig.patch.set_facecolor('white')
        
        ax1 = plt.subplot2grid((2,2), (0,0))
        ax2 = plt.subplot2grid((2,2), (1,0))
        ax3 = plt.subplot2grid((2,2), (0,1))

        #setup basemap projection
        m = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c',ax=ax1)
        m2 = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c',ax=ax2)
        m3 = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c',ax=ax3)

        m.drawcoastlines()
        m.drawmapboundary()
        m2.drawcoastlines()
        m2.drawmapboundary()
        m3.drawcoastlines()
        m3.drawmapboundary()

        #take diff
        if model != 'CMAM':
            diff = param_2100e - param_2000
        if model != 'GISSE2R':
            diff2 = param_2100 - param_2000
        if (model != 'CMAM') & (model != 'GISSE2R'):
            diff3 = param_2100 - param_2100e

        if 'diurnal_ph' in param_str_2000:
            cb = plt.cm.hsv
            cb_l = 'hour'
            
            if model != 'CMAM':
                test = diff < -12
                diff[test] = 12 - (np.abs(diff[test]) - 12) 
                test = diff > 12
                diff[test] = -12 + (np.abs(diff[test]) - 12)
            
            if model != 'GISSE2R':
                test = diff2 < -12
                diff2[test] = 12 - (np.abs(diff2[test]) - 12) 
                test = diff2 > 12
                diff2[test] = -12 + (np.abs(diff2[test]) - 12)
            
            if (model != 'CMAM') & (model != 'GISSE2R'):
                test = diff3 < -12
                diff3[test] = 12 - (np.abs(diff3[test]) - 12) 
                test = diff3 > 12
                diff3[test] = -12 + (np.abs(diff3[test]) - 12)
            
            minv = -12
            maxv = 12

        elif 'seasonal_ph' in param_str_2000:
            cb = plt.cm.hsv
            cb_l = 'month'
            
            if model != 'CMAM':
                test = diff < -6
                diff[test] = 6 - (np.abs(diff[test]) - 6) 
                test = diff > 6
                diff[test] = -6 + (np.abs(diff[test]) - 6)
            
            if model != 'GISSE2R':
                test = diff2 < -6
                diff2[test] = 6 - (np.abs(diff2[test]) - 6) 
                test = diff2 > 6
                diff2[test] = -6 + (np.abs(diff2[test]) - 6)
            
            if (model != 'CMAM') & (model != 'GISSE2R'):
                test = diff3 < -6
                diff3[test] = 6 - (np.abs(diff3[test]) - 6) 
                test = diff > 6
                diff3[test] = -6 + (np.abs(diff3[test]) - 6)
            
            minv = -6
            maxv = 6
    
        else:
            if 'pc_var' in param_str_2000:
                cb_l = '% of variance'
            elif 'total_var' in param_str_2000:
                cb_l = 'variance'
            else:
                cb_l = 'ppb'  
            cb = plt.cm.coolwarm
            
            if model == 'CMAM':
                if (np.abs(np.min(diff2)) > np.max(diff2)):
                    minv = np.min(diff2)
                    maxv = np.abs(minv)
                else:
                    minv = -np.max(diff2)
                    maxv = np.max(diff2)
                    
            elif model == 'GISSE2R':
                if (np.abs(np.min(diff)) > np.max(diff)):
                    minv = np.min(diff)
                    maxv = np.abs(minv)
                else:
                    minv = -np.max(diff)
                    maxv = np.max(diff)
            else:
                if (np.abs(np.min(diff)) > np.abs(np.max(diff))) & (np.abs(np.min(diff)) > np.abs(np.min(diff2))) & (np.abs(np.min(diff)) > np.abs(np.max(diff2))) & (np.abs(np.min(diff)) > np.abs(np.min(diff3))) & (np.abs(np.min(diff)) > np.abs(np.max(diff3))):
                    minv = np.min(diff)
                    maxv = np.abs(minv)
                elif (np.abs(np.min(diff2)) > np.abs(np.max(diff2))) & (np.abs(np.min(diff2)) > np.abs(np.min(diff))) & (np.abs(np.min(diff2)) > np.abs(np.max(diff))) & (np.abs(np.min(diff2)) > np.abs(np.min(diff3))) & (np.abs(np.min(diff2)) > np.abs(np.max(diff3))):
                    minv = np.min(diff2)
                    maxv = np.abs(minv)
                elif (np.abs(np.min(diff3)) > np.abs(np.max(diff3))) & (np.abs(np.min(diff3)) > np.abs(np.min(diff))) & (np.abs(np.min(diff3)) > np.abs(np.max(diff))) & (np.abs(np.min(diff3)) > np.abs(np.min(diff2))) & (np.abs(np.min(diff3)) > np.abs(np.max(diff2))):
                    minv = np.min(diff3)
                    maxv = np.abs(minv)
                elif (np.abs(np.max(diff)) > np.abs(np.min(diff))) & (np.abs(np.max(diff)) > np.abs(np.max(diff2))) & (np.abs(np.max(diff)) > np.abs(np.min(diff2))) & (np.abs(np.max(diff)) > np.abs(np.max(diff3))) & (np.abs(np.max(diff)) > np.abs(np.min(diff3))):
                    minv = -np.max(diff)
                    maxv = np.abs(minv)
                elif (np.abs(np.max(diff2)) > np.abs(np.min(diff2))) & (np.abs(np.max(diff2)) > np.abs(np.max(diff))) & (np.abs(np.max(diff2)) > np.abs(np.min(diff))) & (np.abs(np.max(diff2)) > np.abs(np.max(diff3))) & (np.abs(np.max(diff2)) > np.abs(np.min(diff3))):
                    minv = -np.max(diff2)
                    maxv = np.abs(minv)
                elif (np.abs(np.max(diff3)) > np.abs(np.min(diff3))) & (np.abs(np.max(diff3)) > np.abs(np.max(diff))) & (np.abs(np.max(diff3)) > np.abs(np.min(diff))) & (np.abs(np.max(diff3)) > np.abs(np.max(diff2))) & (np.abs(np.max(diff3)) > np.abs(np.min(diff2))):
                    minv = -np.max(diff3)
                    maxv = np.abs(minv)    
                
         
        if model != 'CMAM':
            pl = m.pcolor(lon_e,lat_e,diff, vmin=minv, vmax=maxv,linewidth=0.5,cmap=cb)
        if model != 'GISSE2R':
            pl = m3.pcolor(lon_e,lat_e,diff2, vmin=minv, vmax=maxv,linewidth=0.5,cmap=cb)
        if (model != 'CMAM') & (model != 'GISSE2R'):
            pl = m2.pcolor(lon_e,lat_e,diff3, vmin=minv, vmax=maxv,linewidth=0.5,cmap=cb)
        plt.tight_layout()
        
        cbar_ax = fig.add_axes([0.53, 0.25, 0.44, 0.07])
        cb = fig.colorbar(pl,orientation='horizontal',format='%.1f',cax=cbar_ax)
        cb.ax.tick_params(labelsize=18)
        cb.set_label('%s change'%(cb_l), fontsize = 23)
        
        if model != 'CMAM':
            ax1.set_title('%s - %s'%(param_str_2100e,param_str_2000),fontsize=20) 
        if model != 'GISSE2R':
            ax3.set_title('%s - %s'%(param_str_2100,param_str_2000),fontsize=20) 
        if (model != 'CMAM') & (model != 'GISSE2R'):
            ax2.set_title('%s - %s'%(param_str_2100,param_str_2100e),fontsize=20) 
        
        plt.savefig('plots/diff_%s_%s.png'%(model,param_str_2000))
        #plt.show()




