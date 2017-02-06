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

areas = ['SW_NA','NW_NA','NE_NA','CE_NA','SE_NA','C_NA','S_EU','SW_EU','CW_EU','CS_EU','C_EU','E_EU','NW_EU','N_EU']
cmaplist = plt.cm.get_cmap('Set1')                                                                                                                       
array = np.linspace(0,1,len(areas))
cmaplist = [cmaplist(i) for i in array]

area_boundaries,area_tags,area_labels = modules.area_dicts()
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
    monthly_2000_ch4 = np.average(monthly_2000.variables['ch4'][12:24,:,:],axis=0)
    monthly_2000_co = np.average(monthly_2000.variables['co'][12:24,:,:],axis=0)
    monthly_2000_no = np.average(monthly_2000.variables['no'][12:24,:,:],axis=0)
    monthly_2000_no2 = np.average(monthly_2000.variables['no2'][12:24,:,:],axis=0)    
    monthly_2000_nox = monthly_2000_no + monthly_2000_no2
    monthly_2000_isop = np.average(monthly_2000.variables['isop'][12:24,:,:],axis=0)
    monthly_2000_eminox = np.sum(monthly_2000.variables['eminox'][12:24,:,:],axis=0)
    monthly_2000_emivoc = np.sum(monthly_2000.variables['emivoc'][12:24,:,:],axis=0)
    monthly_2000_eminh3 = np.sum(monthly_2000.variables['eminh3'][12:24,:,:],axis=0)
    monthly_2100_ch4 = np.average(monthly_2100.variables['ch4'][72:84,:,:],axis=0)
    monthly_2100_co = np.average(monthly_2100.variables['co'][72:84,:,:],axis=0)
    monthly_2100_no = np.average(monthly_2100.variables['no'][72:84,:,:],axis=0)
    monthly_2100_no2 = np.average(monthly_2100.variables['no2'][72:84,:,:],axis=0)    
    monthly_2100_nox = monthly_2100_no + monthly_2100_no2
    monthly_2100_isop = np.average(monthly_2100.variables['isop'][72:84,:,:],axis=0)
    monthly_2100_eminox = np.sum(monthly_2100.variables['eminox'][72:84,:,:],axis=0)
    monthly_2100_emivoc = np.sum(monthly_2100.variables['emivoc'][72:84,:,:],axis=0)
    monthly_2100_eminh3 = np.sum(monthly_2100.variables['eminh3'][72:84,:,:],axis=0)
    monthly_2100e_ch4 = np.average(monthly_2100e.variables['ch4'][72:84,:,:],axis=0)
    monthly_2100e_co = np.average(monthly_2100e.variables['co'][72:84,:,:],axis=0)
    monthly_2100e_no = np.average(monthly_2100e.variables['no'][72:84,:,:],axis=0)
    monthly_2100e_no2 = np.average(monthly_2100e.variables['no2'][72:84,:,:],axis=0)                                                                                            
    monthly_2100e_nox = monthly_2100e_no + monthly_2100e_no2
    monthly_2100e_isop = np.average(monthly_2100e.variables['isop'][72:84,:,:],axis=0)
    monthly_2100e_eminox = np.sum(monthly_2100e.variables['eminox'][72:84,:,:],axis=0)
    monthly_2100e_emivoc = np.sum(monthly_2100e.variables['emivoc'][72:84,:,:],axis=0)
    monthly_2100e_eminh3 = np.sum(monthly_2100e.variables['eminh3'][72:84,:,:],axis=0)

    vocnox_2000 = (monthly_2000_co+monthly_2000_ch4+monthly_2000_isop)/monthly_2000_nox
    vocnox_2100 = (monthly_2100_co+monthly_2100_ch4+monthly_2100_isop)/monthly_2100_nox
    vocnox_2100e = (monthly_2100e_co+monthly_2100e_ch4+monthly_2100e_isop)/monthly_2100e_nox

    diff_vocnox_a = vocnox_2100e - vocnox_2000
    diff_vocnox_b = vocnox_2100 - vocnox_2000
    diff_vocnox_c = vocnox_2100e - vocnox_2100
    print np.min(diff_vocnox_a),np.max(diff_vocnox_a)

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

    lim = 10000
    
    #pl = m.pcolor(lon_e,lat_e,diff_vocnox_a, vmin=-lim, vmax=lim,linewidth=0.5,cmap=plt.cm.coolwarm)
    #pl = m2.pcolor(lon_e,lat_e,diff_vocnox_b, vmin=-lim, vmax=lim,linewidth=0.5,cmap=plt.cm.coolwarm)
    #pl = m3.pcolor(lon_e,lat_e,diff_vocnox_c, vmin=-lim, vmax=lim,linewidth=0.5,cmap=plt.cm.coolwarm)
    #pl = m.pcolor(lon_e,lat_e,vocnox_2000, vmin=0, vmax=1000,linewidth=0.5,cmap=plt.cm.gist_earth)
    #pl = m2.pcolor(lon_e,lat_e,vocnox_2100, vmin=0, vmax=1000,linewidth=0.5,cmap=plt.cm.gist_earth)
    #pl = m3.pcolor(lon_e,lat_e,vocnox_2100e, vmin=0, vmax=1000,linewidth=0.5,cmap=plt.cm.gist_earth)
    pl = m.pcolor(lon_e,lat_e,monthly_2000_isop, vmin=0, vmax=5,linewidth=0.5,cmap=plt.cm.gist_earth)
    pl = m2.pcolor(lon_e,lat_e,monthly_2100_isop, vmin=0, vmax=5,linewidth=0.5,cmap=plt.cm.gist_earth)
    pl = m3.pcolor(lon_e,lat_e,monthly_2100e_isop, vmin=0, vmax=5,linewidth=0.5,cmap=plt.cm.gist_earth)
    plt.tight_layout()

    plt.show()
