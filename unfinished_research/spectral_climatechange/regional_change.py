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

#models = ['CESMCAM', 'CMAM', 'GFDLAM3','GISSE2R', 'MIROCCHEM', 'MOCAGE', 'UMCAM']
models = ['GFDLAM3']

areas = modules.model_areas()
#areas = ['CE_NA', 'SE_NA', 'S_NA', 'SW_NA', 'C_NA', 'NE_NA', 'NW_NA', 'N_EU', 'C_EU', 'SW_EU', 'S_EU', 'E_EU', 'NW_EU', 'N_O', 'S_O']

print areas

area_boundaries,area_labels = modules.model_area_dicts()
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
    size = len(lat_c)*len(lon_c)

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
   
    #read in monthly emissions and species for year 2000
    monthly_2000 = Dataset('/work/home/db876/modelling/ACCMIP/GFDLAM3/GFDLAM3_SURFACE_2000_2012_*_*_*_M_*.nc')
    monthly_2100 = Dataset('/work/home/db876/modelling/ACCMIP/GFDLAM3/GFDLAM3_SURFACE_2095_2111_*_*_*_M_rcp85.nc')
    monthly_2100e = Dataset('/work/home/db876/modelling/ACCMIP/GFDLAM3/GFDLAM3_SURFACE_2095_2111_*_*_*_M_rcp85em2000.nc') 
    monthly_2000_ch4 = np.sum(monthly_2000.variables['ch4'][12:24,:,:],axis=0)
    monthly_2000_co = np.sum(monthly_2000.variables['co'][12:24,:,:],axis=0)
    monthly_2000_no = np.sum(monthly_2000.variables['no'][12:24,:,:],axis=0)
    monthly_2000_no2 = np.sum(monthly_2000.variables['no2'][12:24,:,:],axis=0)    
    monthly_2000_nox = monthly_2000_no + monthly_2000_no2
    monthly_2000_isop = np.sum(monthly_2000.variables['isop'][12:24,:,:],axis=0)
    monthly_2000_eminox = np.sum(monthly_2000.variables['eminox'][12:24,:,:],axis=0)
    monthly_2000_emivoc = np.sum(monthly_2000.variables['emivoc'][12:24,:,:],axis=0)
    monthly_2000_eminh3 = np.sum(monthly_2000.variables['eminh3'][12:24,:,:],axis=0)
    monthly_2100_ch4 = np.sum(monthly_2100.variables['ch4'][72:84,:,:],axis=0)
    monthly_2100_co = np.sum(monthly_2100.variables['co'][72:84,:,:],axis=0)
    monthly_2100_no = np.sum(monthly_2100.variables['no'][72:84,:,:],axis=0)
    monthly_2100_no2 = np.sum(monthly_2100.variables['no2'][72:84,:,:],axis=0)    
    monthly_2100_nox = monthly_2100_no + monthly_2100_no2
    monthly_2100_isop = np.sum(monthly_2100.variables['isop'][72:84,:,:],axis=0)
    monthly_2100_eminox = np.sum(monthly_2100.variables['eminox'][72:84,:,:],axis=0)
    monthly_2100_emivoc = np.sum(monthly_2100.variables['emivoc'][72:84,:,:],axis=0)
    monthly_2100_eminh3 = np.sum(monthly_2100.variables['eminh3'][72:84,:,:],axis=0)
    monthly_2100e_ch4 = np.sum(monthly_2100e.variables['ch4'][72:84,:,:],axis=0)
    monthly_2100e_co = np.sum(monthly_2100e.variables['co'][72:84,:,:],axis=0)
    monthly_2100e_no = np.sum(monthly_2100e.variables['no'][72:84,:,:],axis=0)
    monthly_2100e_no2 = np.sum(monthly_2100e.variables['no2'][72:84,:,:],axis=0)                                                                                            
    monthly_2100e_nox = monthly_2100e_no + monthly_2100e_no2
    monthly_2100e_isop = np.sum(monthly_2100e.variables['isop'][72:84,:,:],axis=0)
    monthly_2100e_eminox = np.sum(monthly_2100e.variables['eminox'][72:84,:,:],axis=0)
    monthly_2100e_emivoc = np.sum(monthly_2100e.variables['emivoc'][72:84,:,:],axis=0)
    monthly_2100e_eminh3 = np.sum(monthly_2100e.variables['eminh3'][72:84,:,:],axis=0)

    #set up plot
    fig, axes = plt.subplots(nrows=10, ncols=10,figsize=(19,13))
    fig.patch.set_facecolor('white')
        
    #setup basemap projection
    m = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c')

    count = 0
    r1 = [0,2,4,6,8,20,22,24,26,28,40,42,44,46,48,60,62,64,66,68,80,82,84,86]
    r2 = [1,3,5,7,9,21,23,25,27,29,41,43,45,47,49,61,63,65,67,69,81,83,85,87]
    r3 = [10,12,14,16,18,30,32,34,36,38,50,52,54,56,58,70,72,74,76,78,90,92,94,96]
    r4 = [11,13,15,17,19,31,33,35,37,39,51,53,55,57,59,71,73,75,77,79,91,93,95,97]
    for area in areas:
        print area
        print count
        area_lat_inds = []
        area_lon_inds = []         

        #area_color = cmaplist[count]
        area_grid = area_boundaries[area]
        #area_tag = area_tags[area]
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
        
        area_mean_2000 = np.average(mean_2000[area_lat_inds,area_lon_inds])
        area_mean_2100e = np.average(mean_2100e[area_lat_inds,area_lon_inds])    
        area_mean_2100 = np.average(mean_2100[area_lat_inds,area_lon_inds])

        area_winter_mean_2000 = np.average(mean_winter_2000[area_lat_inds,area_lon_inds])
        area_winter_mean_2100e = np.average(mean_winter_2100e[area_lat_inds,area_lon_inds])    
        area_winter_mean_2100 = np.average(mean_winter_2100[area_lat_inds,area_lon_inds])

        area_summer_mean_2000 = np.average(mean_summer_2000[area_lat_inds,area_lon_inds])
        area_summer_mean_2100e = np.average(mean_summer_2100e[area_lat_inds,area_lon_inds])    
        area_summer_mean_2100 = np.average(mean_summer_2100[area_lat_inds,area_lon_inds])

        area_seasonal_mag_2000 = np.average(seasonal_mag_2000[area_lat_inds,area_lon_inds])
        area_seasonal_mag_2100e = np.average(seasonal_mag_2100e[area_lat_inds,area_lon_inds])    
        area_seasonal_mag_2100 = np.average(seasonal_mag_2100[area_lat_inds,area_lon_inds])

        area_pc_var_daily_2000 = np.average(pc_var_daily_2000[area_lat_inds,area_lon_inds])
        area_pc_var_daily_2100e = np.average(pc_var_daily_2100e[area_lat_inds,area_lon_inds])    
        area_pc_var_daily_2100 = np.average(pc_var_daily_2100[area_lat_inds,area_lon_inds])

        area_pc_var_seasonal_2000 = np.average(pc_var_seasonal_2000[area_lat_inds,area_lon_inds])                                                              
        area_pc_var_seasonal_2100e = np.average(pc_var_seasonal_2100e[area_lat_inds,area_lon_inds])    
        area_pc_var_seasonal_2100 = np.average(pc_var_seasonal_2100[area_lat_inds,area_lon_inds])

        area_pc_var_full_2000 = np.average(pc_var_full_2000[area_lat_inds,area_lon_inds])
        area_pc_var_full_2100e = np.average(pc_var_full_2100e[area_lat_inds,area_lon_inds])    
        area_pc_var_full_2100 = np.average(pc_var_full_2100[area_lat_inds,area_lon_inds])

        area_pc_var_noise_2000 = np.average(pc_var_noise_2000[area_lat_inds,area_lon_inds])
        area_pc_var_noise_2100e = np.average(pc_var_noise_2100e[area_lat_inds,area_lon_inds])    
        area_pc_var_noise_2100 = np.average(pc_var_noise_2100[area_lat_inds,area_lon_inds])

        area_total_var_2000 = np.average(total_var_2000[area_lat_inds,area_lon_inds])
        area_total_var_2100e = np.average(total_var_2100e[area_lat_inds,area_lon_inds])    
        area_total_var_2100 = np.average(total_var_2100[area_lat_inds,area_lon_inds])

        all_var = [area_total_var_2000,area_total_var_2100,area_total_var_2100e]
        pie2000 = [area_pc_var_daily_2000,area_pc_var_seasonal_2000,area_pc_var_noise_2000]
        pie2100 = [area_pc_var_daily_2100,area_pc_var_seasonal_2100,area_pc_var_noise_2100]
        pie2100e = [area_pc_var_daily_2100e,area_pc_var_seasonal_2100e,area_pc_var_noise_2100e]
        colors = ['yellowgreen', 'red', 'lightskyblue']

        ax = axes.flat[r1[count]]
        ax.annotate(area,xy=(0.20,0.94),xycoords='axes fraction',alpha=5)
        p1,t1 = ax.pie(pie2000,colors=colors,radius=(all_var[0]/np.max(all_var))*1.02)    
        ax = axes.flat[r2[count]]
        p2,t2 = ax.pie(pie2100,colors=colors,radius=(all_var[1]/np.max(all_var))*1.02)
        ax = axes.flat[r3[count]]
        p3,t3 = ax.pie(pie2100e,colors=colors,radius=(all_var[2]/np.max(all_var))*1.02)
        ax = axes.flat[r4[count]]
        ax.axis('off')
        
        width = 1
        ch4 = [np.average(monthly_2000_ch4[area_lat_inds,area_lon_inds]),np.average(monthly_2100_ch4[area_lat_inds,area_lon_inds]),np.average(monthly_2100e_ch4[area_lat_inds,area_lon_inds])]
        co = [np.average(monthly_2000_co[area_lat_inds,area_lon_inds]),np.average(monthly_2100_co[area_lat_inds,area_lon_inds]),np.average(monthly_2100e_co[area_lat_inds,area_lon_inds])]
        nox = [np.average(monthly_2000_nox[area_lat_inds,area_lon_inds]),np.average(monthly_2100_nox[area_lat_inds,area_lon_inds]),np.average(monthly_2100e_nox[area_lat_inds,area_lon_inds])]
        isop = [np.average(monthly_2000_isop[area_lat_inds,area_lon_inds]),np.average(monthly_2100_isop[area_lat_inds,area_lon_inds]),np.average(monthly_2100e_isop[area_lat_inds,area_lon_inds])]
        
        #inds = np.arange(3)
        #ax.plot(inds, nox, color='b')
        #ax.plot(inds, isop, color='g')
        #ax.plot(inds, co, color='r')        
        #ax.plot(inds, ch4, color='y') 

        count+=1
   
    axes.flat[88].axis('off')
    axes.flat[89].axis('off')
    axes.flat[98].axis('off')
    axes.flat[99].axis('off')
    #cbar_ax = fig.add_axes([0.53, 0.25, 0.44, 0.07])
    #cb = fig.colorbar(pl,orientation='horizontal',format='%.1f',cax=cbar_ax)
    #cb.ax.tick_params(labelsize=18)
    #cb.set_label('%s'%(cb_l), fontsize = 23)
        
    #ax1.set_title('%s'%(param_str_2000),fontsize=20) 

    plt.subplots_adjust(wspace=0.01,hspace=0.01)
    plt.show()
