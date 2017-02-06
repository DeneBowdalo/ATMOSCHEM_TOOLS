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
import matplotlib.dates as dates

plot_type = 's'

#models = ['CESMCAM', 'CMAM', 'GFDLAM3','GISSE2R', 'MIROCCHEM', 'MOCAGE', 'UMCAM']
models = ['GFDLAM3']

areas = modules.model_areas()
area_boundaries,area_labels = modules.model_area_dicts()

#read in obs
fobs = '/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_O3_2000_2003_H_HP.nc'
obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_raw_class,obs_anthrome_class,obs_gap_inds = modules.read_obs_all(fobs,'O3',2000,2003)

if plot_type == 'd':
    dt= obs_datetime_time[:24]
if plot_type == 's':
    dt = obs_datetime_time[:8766]

for i in range(len(models)):

    fig, axes = plt.subplots(nrows=5, ncols=5,figsize=(19,13))
    fig.patch.set_facecolor('white')

    model = models[i]
    print model
    
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

    print 'a'

    #read in 2000 model period data
    f2000 = '/work/home/db876/grid/O3/2000_2012/%s/%s_SURFACE_*_*_*_H_*/LSP_stats.nc'%(year2000,model)

    #read in 2100 model period data
    f2100 = '/work/home/db876/grid/O3/2095_2111/%s/%s_SURFACE_*_*_*_H_rcp85/LSP_stats.nc'%(year2100,model)

    #read in 2100 model fixed emissions period data
    f2100e = '/work/home/db876/grid/O3/2095_2111/%s/%s_SURFACE_*_*_*_H_rcp85em2000/LSP_stats.nc'%(year2100,model)

    root_2000 = Dataset(f2000) 
    lat_e = root_2000.variables['lat_edges'][:]
    lon_e = root_2000.variables['lon_edges'][:]
    lat_c = root_2000.variables['lat_centre'][:]
    lon_c = root_2000.variables['lon_centre'][:]
    size = len(lat_c)*len(lon_c)

    m = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c')

    print 'b'

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
    
    print 'c'

    if model != 'CMAM':
        (diurnal_mag_2100e,diurnal_ph_2100e,seasonal_mag_2100e,seasonal_ph_2100e,mean_2100e,p1_2100e,p5_2100e,p25_2100e,p50_2100e,p75_2100e,p95_2100e,p99_2100e,diurnal_ave_waveform_2100e,seasonal_waveform_2100e,full_ave_waveform_2100e,pc_var_daily_2100e,pc_var_seasonal_2100e,pc_var_full_2100e,pc_var_noise_2100e,total_var_2100e,
        diurnal_mag_spring_2100e,diurnal_ph_spring_2100e,mean_spring_2100e,p1_spring_2100e,p5_spring_2100e,p25_spring_2100e,p50_spring_2100e,p75_spring_2100e,p95_spring_2100e,p99_spring_2100e,diurnal_waveform_spring_2100e,
        diurnal_mag_summer_2100e,diurnal_ph_summer_2100e,mean_summer_2100e,p1_summer_2100e,p5_summer_2100e,p25_summer_2100e,p50_summer_2100e,p75_summer_2100e,p95_summer_2100e,p99_summer_2100e,diurnal_waveform_summer_2100e,
        diurnal_mag_autumn_2100e,diurnal_ph_autumn_2100e,mean_autumn_2100e,p1_autumn_2100e,p5_autumn_2100e,p25_autumn_2100e,p50_autumn_2100e,p75_autumn_2100e,p95_autumn_2100e,p99_autumn_2100e,diurnal_waveform_autumn_2100e,
        diurnal_mag_winter_2100e,diurnal_ph_winter_2100e,mean_winter_2100e,p1_winter_2100e,p5_winter_2100e,p25_winter_2100e,p50_winter_2100e,p75_winter_2100e,p95_winter_2100e,p99_winter_2100e,diurnal_waveform_winter_2100e,
        seasonal_mag_day_2100e,seasonal_ph_day_2100e,mean_day_2100e,p1_day_2100e,p5_day_2100e,p25_day_2100e,p50_day_2100e,p75_day_2100e,p95_day_2100e,p99_day_2100e,seasonal_waveform_day_2100e,
        seasonal_mag_night_2100e,seasonal_ph_night_2100e,mean_night_2100e,p1_night_2100e,p5_night_2100e,p25_night_2100e,p50_night_2100e,p75_night_2100e,p95_night_2100e,p99_night_2100e,seasonal_waveform_night_2100e,
        daily_h3_mag_2100e,daily_h2_mag_2100e,daily_h1_mag_2100e,daily_mag_2100e,annual_h3_mag_2100e,annual_h2_mag_2100e,annual_h1_mag_2100e,annual_mag_2100e) = modules.get_periodic_specific('na',f2100e)  

    print 'd'

    count = 0 
    for ax in axes.flat:
        try:
            area = areas[count]
        except:
            ax.axis('off')
            continue

        area_lat_inds = []
        area_lon_inds = []

        print area

        area_grid = area_boundaries[area]
        area_label = area_labels[area]

        lat_i = 0
        lon_i = 0
        for i in range(size):                                                                                                                            
            #check centre of gridbox within area limits
            if (lat_c[lat_i] >= area_grid[0]) & (lat_c[lat_i] < area_grid[1]) & (lon_c[lon_i] >= area_grid[2]) & (lon_c[lon_i] < area_grid[3]):
                #check centre of gridbox is over land
                if m.is_land(lon_c[lon_i],lat_c[lat_i]) == True:
                    area_lat_inds.append(lat_i)
                    area_lon_inds.append(lon_i)
        
            if lon_i == (len(lon_c)-1):                                                                                                                  
                lat_i+=1
                lon_i=0
            else:
                lon_i+=1

        area_lat_inds = np.array(area_lat_inds)    
        area_lon_inds = np.array(area_lon_inds)

        print seasonal_waveform_2000.shape
        if plot_type == 'd':
            model2000_w = np.nanmean(diurnal_ave_waveform_2000[:,area_lat_inds,area_lon_inds],axis=1)
            if model != 'GISSE2R':
                model2100_w = np.nanmean(diurnal_ave_waveform_2100[:,area_lat_inds,area_lon_inds],axis=1)
            if model != 'CMAM':
                model2100e_w = np.nanmean(diurnal_ave_waveform_2100e[:,area_lat_inds,area_lon_inds],axis=1)
        
        if plot_type == 's':
            model2000_w = np.nanmean(seasonal_waveform_2000[:,area_lat_inds,area_lon_inds],axis=1)
            if model != 'GISSE2R':
                model2100_w = np.nanmean(seasonal_waveform_2100[:,area_lat_inds,area_lon_inds],axis=1)
            if model != 'CMAM':
                model2100e_w = np.nanmean(seasonal_waveform_2100e[:,area_lat_inds,area_lon_inds],axis=1)
        
        print model2000_w.shape
        ax.annotate(area,xy=(0.01,0.88),xycoords='axes fraction',alpha=5)
        ax.plot_date(dt,model2000_w,color = 'blue',linestyle='-',linewidth=1,markersize=1,markeredgecolor='None')
        if model != 'GISSE2R':
            ax.plot_date(dt,model2100_w,color = 'red',linestyle='-',linewidth=1,markersize=1,markeredgecolor='None')
        if model != 'CMAM':
            ax.plot_date(dt,model2100e_w,color = 'green',linestyle='-',linewidth=1,markersize=1,markeredgecolor='None')
    
        if plot_type == 'd':
            ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        if plot_type == 's':
            ax.xaxis.set_major_formatter(dates.DateFormatter('%m'))
            
        count+=1
        
    h1, = ax.plot([1,1],color='blue',marker='o',linestyle='None',markersize=10)
    h2, = ax.plot([1,1],color='red',marker='o',linestyle='None',markersize=10)
    h3, = ax.plot([1,1],color='green',marker='o',linestyle='None',markersize=10)

    plt.legend((h1,h2,h3),['%s 2000'%(model),'%s 2100 RCP8.5'%(model),'%s 2100 RCP8.5 2000 Emissions'%(model)],loc='lower left',prop={'size':13},fancybox=True,ncol=1,markerscale=1,bbox_to_anchor=(-1.0,-0.1))
    h1.set_visible(False)
    h2.set_visible(False)
    h3.set_visible(False)

    #plt.savefig('plots/area_%s_%s.png'%(plot_type,model))

    plt.show()
