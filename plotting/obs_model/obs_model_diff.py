#import matplotlib
#matplotlib.use('Agg')
import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm
from cmath import *
import colorsys
from scipy import stats
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
import modules
from collections import OrderedDict
import operator
from netCDF4 import Dataset
from scipy.odr import Model, RealData, ODR, Data
from pylab import *
import pandas as pd

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-3]
date_range = paths[-2]

start_year = date_range[:4]

def form2(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.2f' % x	

def form5(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.7f' % x

xformatter = FuncFormatter(form2)
yformatter = FuncFormatter(form5)

def onpick(event):
    global pl
    
    global ind
    global fig2
    
    ind = event.ind
    ind = ind[0]
    #x_data = event.xdata
    #y_data = event.ydata

    #find ind of closest lat/lon
    #ind = modules.find_nearest_point_index(obs_lons,obs_lats,x_data,y_data)
    
    try:
        for i in range(len(pl)):
            pl.pop(0).remove()
            first_run = False  
        
    except:
        first_run = True
        pass
    
    pl = m.plot([X[ind]], [Y[ind]], 'o', ms=12, alpha=0.6, color='yellow',zorder=20)
    
    #get model timeseries for site clicked
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lats[ind],obs_lons[ind])
    model_var_pick = model_var[:,lat_n,lon_n]
    model_var_pick = model_var_pick*1e9
    model_var_mask = np.ma.masked_where(model_var_pick<=0,model_var_pick)
   
    if model_name == 'MACC':
        model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
        count = 0
        valids = []
        for i in range(len(model_time_pd)):
            if count == 0:
                valids.append(i)
                count+=1
            elif count == 2:
                count = 0
            else:
                count+=1
        model_time_pd = model_time_pd[valids]
        model_var_pd = pd.Series(model_var_mask, index=model_time_pd) 
    else:    
        model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
        model_var_pd = pd.Series(model_var_mask, index=model_time_pd)
    
    #get obs timeseries for site clicked
    ref = obs_refs[ind]
    obs_ts_group = obs_root_grp.groups[ref]
    obs_var = obs_ts_group.variables[species.lower()][:]
    group = obs_ts_group.process_group
    lat = obs_ts_group.latitude
    lon = obs_ts_group.longitude
    lon = obs_ts_group.longitude
    alt = obs_ts_group.altitude
    complete = obs_ts_group.completeness
    a_class = obs_ts_group.anthrome_class
    r_class = obs_ts_group.raw_class
    continent = loc_dict[tags[ind]]
    country = obs_ts_group.country
    
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
    obs_var_pd = pd.Series(obs_var_mask, index=obs_time_pd)
    
    #create sine wave from amp/phase
    obs_date_l = obs_date.astype(int)
    obs_time_l = obs_time.astype(int)
    obs_times = modules.date_process(obs_date_l,obs_time_l,start_year)
    obs_times = np.array(obs_times)
    pi2 = np.pi*2
    
    #convert phases to radians
    calc = pi2/6.
    obs_ha_phase_r = obs_ha_phase[ind] * calc
    calc = pi2/12.
    obs_annual_phase_r = obs_annual_phase[ind] * calc
    
    ha_obs_wave = obs_ha_mag[ind]*(np.cos((pi2*obs_times/(365.25/2.))-(obs_ha_phase_r)))
    annual_obs_wave = obs_annual_mag[ind]*(np.cos((pi2*obs_times/(365.25))-(obs_annual_phase_r)))
    seasonal_obs_wave = (ha_obs_wave+annual_obs_wave)+obs_ave[ind]
    obs_seasonal_wave_pd = pd.Series(seasonal_obs_wave, index=obs_time_pd)
    
    #create sine wave from amp/phase
    mod_date_l = model_date.astype(int)
    mod_time_l = model_time.astype(int)
    mod_times = modules.date_process(mod_date_l,mod_time_l,start_year)
    mod_times = np.array(mod_times)
    pi2 = np.pi*2
    
    #convert phases to radians
    calc = pi2/6.
    model_ha_phase_r = model_ha_phase[ind] * calc
    calc = pi2/12.
    model_annual_phase_r = model_annual_phase[ind] * calc
    
    ha_model_wave = model_ha_mag[ind]*(np.cos((pi2*mod_times/(365.25/2.))-(model_ha_phase_r)))
    annual_model_wave = model_annual_mag[ind]*(np.cos((pi2*mod_times/(365.25))-(model_annual_phase_r)))
    seasonal_model_wave = (ha_model_wave+annual_model_wave)+model_ave[ind]
    model_seasonal_wave_pd = pd.Series(seasonal_model_wave, index=model_time_pd)

    
    #get spectra data
    site_group_obs = root_grp_obs_spec.groups[ref]
    site_group_mod = root_grp_mod_spec.groups[ref]
    
    obs_period = site_group_obs.variables['period'][:]
    mod_period = site_group_mod.variables['period'][:]
    
    obs_amp = site_group_obs.variables['amplitude'][:]
    mod_amp = site_group_mod.variables['amplitude'][:]
    
    fig.canvas.draw()
        
    if first_run == False:
        plt.close(fig2)
        fig2, (axo,axo2) = plt.subplots(2,figsize=(24,12))
        fig2.patch.set_facecolor('white')
        
        #fig2 = plt.figure()
        
        axo.plot_date(obs_time_pd.to_pydatetime(), obs_var_pd, color='black', markersize = 3, label = 'Observations')
        axo.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='red',alpha=0.5, markersize = 3, label = '%s %s %s %s'%(model_name,version,grid_size,met),markeredgecolor='None')
        axo.plot_date(obs_time_pd.to_pydatetime(), obs_seasonal_wave_pd, color='yellow', markersize = 3, label = 'Obs Seasonal Waveform',markeredgecolor='None')
        axo.plot_date(model_time_pd.to_pydatetime(), model_seasonal_wave_pd, color='green', markersize = 3, label = 'Model Seasonal Waveform',markeredgecolor='None')
        
        axo2.loglog(obs_period,obs_amp,color='black',label='Obs')
        axo2.loglog(mod_period,mod_amp,color='red',label = '%s %s %s %s'%(model_name,version,grid_size,met))

        axo2.text(0.01, 0.95, 'Obs D Amp = %8.2f ppb'%(obs_daily_mag[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.92, 'Model D Amp = %8.2f ppb'%(model_daily_mag[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.text(0.01, 0.85, 'Obs HA Amp = %8.2f ppb'%(obs_ha_mag[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.82, 'Model HA Amp = %8.2f ppb'%(model_ha_mag[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.text(0.01, 0.75, 'Obs A Amp = %8.2f ppb'%(obs_annual_mag[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.72, 'Model A Amp = %8.2f ppb'%(model_annual_mag[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.text(0.01, 0.55, 'Obs D Phase = %8.2f'%(obs_daily_phase[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.52, 'Model D Phase = %8.2f'%(model_daily_phase[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.text(0.01, 0.45, 'Obs HA Phase = %8.2f'%(obs_ha_phase[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.42, 'Model HA Phase = %8.2f'%(model_ha_phase[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        obs_a_ph = obs_annual_phase[ind]
        mod_a_ph = model_annual_phase[ind]
        
        if obs_a_ph > 12:
            obs_a_ph = obs_a_ph-12.
        if mod_a_ph > 12:
            mod_a_ph = mod_a_ph-12.
        
        axo2.text(0.01, 0.35, 'Obs A Phase = %8.2f'%(obs_a_ph),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.32, 'Model A Phase = %8.2f'%(mod_a_ph),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.text(0.01, 0.15, 'Obs Ave = %8.2f ppb'%(obs_ave[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.12, 'Model Ave = %8.2f ppb'%(model_ave[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.axvline(1.,ymin=0,ymax=1,color='blue',linestyle='--')
        axo2.axvline(182.625,ymin=0,ymax=1,color='blue',linestyle='--')
        axo2.axvline(365.25,ymin=0,ymax=1,color='blue',linestyle='--')
        
        axo2.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axo2.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
        plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
        
        
        axo.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s, Grid Index = %s,%s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,lat_n,lon_n)) 
        
        plt.legend(loc = 'lower right')
        plt.tight_layout()
        axo.grid()
        axo2.grid()
        
        plt.show()
    else:
        #fig2 = plt.figure()
        fig2, (axo,axo2) = plt.subplots(2,figsize=(24,12))
        fig2.patch.set_facecolor('white')
        
        axo.plot_date(obs_time_pd.to_pydatetime(), obs_var_pd, color='black', markersize = 3, label = 'Observations')
        axo.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='red', markersize = 3,alpha=0.5, label = '%s %s %s %s'%(model_name,version,grid_size,met),markeredgecolor='None')
        axo.plot_date(obs_time_pd.to_pydatetime(), obs_seasonal_wave_pd, color='yellow', markersize = 3, label = 'Obs Seasonal Waveform',markeredgecolor='None')
        axo.plot_date(model_time_pd.to_pydatetime(), model_seasonal_wave_pd, color='green', markersize = 3, label = 'Model Seasonal Waveform',markeredgecolor='None')
        
        axo2.loglog(obs_period,obs_amp,color='black',label='Obs')
        axo2.loglog(mod_period,mod_amp,color='red', label = '%s %s %s %s'%(model_name,version,grid_size,met))

        axo2.text(0.01, 0.95, 'Obs D Amp = %8.2f ppb'%(obs_daily_mag[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.92, 'Model D Amp = %8.2f ppb'%(model_daily_mag[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.text(0.01, 0.85, 'Obs HA Amp = %8.2f ppb'%(obs_ha_mag[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.82, 'Model HA Amp = %8.2f ppb'%(model_ha_mag[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.text(0.01, 0.75, 'Obs A Amp = %8.2f ppb'%(obs_annual_mag[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.72, 'Model A Amp = %8.2f ppb'%(model_annual_mag[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.text(0.01, 0.55, 'Obs D Phase = %8.2f'%(obs_daily_phase[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.52, 'Model D Phase = %8.2f'%(model_daily_phase[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.text(0.01, 0.45, 'Obs HA Phase = %8.2f'%(obs_ha_phase[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.42, 'Model HA Phase = %8.2f'%(model_ha_phase[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        obs_a_ph = obs_annual_phase[ind]
        mod_a_ph = model_annual_phase[ind]
        
        if obs_a_ph > 12:
            obs_a_ph = obs_a_ph-12.
        if mod_a_ph > 12:
            mod_a_ph = mod_a_ph-12.
        
        axo2.text(0.01, 0.35, 'Obs A Phase = %8.2f'%(obs_a_ph),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.32, 'Model A Phase = %8.2f'%(mod_a_ph),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.text(0.01, 0.15, 'Obs Ave = %8.2f ppb'%(obs_ave[ind]),transform=axo2.transAxes,fontweight='bold')
        axo2.text(0.01, 0.12, 'Model Ave = %8.2f ppb'%(model_ave[ind]),transform=axo2.transAxes,fontweight='bold',color='red')
        
        axo2.axvline(1.,ymin=0,ymax=1,color='blue',linestyle='--')
        axo2.axvline(182.625,ymin=0,ymax=1,color='blue',linestyle='--')
        axo2.axvline(365.25,ymin=0,ymax=1,color='blue',linestyle='--')
        
        axo2.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axo2.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
        plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
        
        
        axo.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s, Grid Index = %s,%s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,lat_n,lon_n))
        
        plt.legend(loc = 'lower right')
        plt.tight_layout()
        axo.grid()
        axo2.grid()
        
        plt.show()

#-----------------------------------------------------
#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-3]
model_range = paths[-2]
print '\nSpecies is %s\n'%(species)

model_details = paths[-1]
print model_details

model_split = model_details.split('_SFC')
model_name = model_split[0]
model_other = model_split[1]
model_other_split = model_other.split('_')
#model_other_split = model_other_split.remove('')

#read in model time series data

if len(model_other_split) == 1:
    version = ''
    grid = ''
    met = ''
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'.nc'

elif len(model_other_split) == 2:
    version = model_other_split[1] 
    grid = ''
    met = ''
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'.nc'

elif len(model_other_split) == 3:
    version = model_other_split[1] 
    grid = model_other_split[2]
    met = ''
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'_'+grid+'.nc'

elif len(model_other_split) == 4:
    version = model_other_split[1] 
    grid = model_other_split[2]
    met = model_other_split[3]
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'_'+grid+'_'+met+'.nc'
    

root_grp = Dataset(model_f)
model_var = root_grp.variables[species.lower()][:]
model_date = root_grp.variables['date'][:]
model_time = root_grp.variables['time'][:]
lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
grid_size = root_grp.variables['grid_size'][:]
grid_size = grid_size[0]

#process model dates and model times to datetimes, then process pandas objects

year_val = []
month_val = []
day_val = []
hour_val = []
minute_val = []

#change model times to datetimes
model_date = model_date.astype('str')
model_time = model_time.astype('str')

for date in model_date:
    year_val.append(int(date[0:4]))
    month_val.append(int(date[4:6]))
    day_val.append(int(date[6:8]))

for time in model_time:
    if np.float64(time) == 0:
        hour_val.append(0)
        minute_val.append(0)
    elif np.float64(time) == 30:
        hour_val.append(0)
        minute_val.append(30)
    else:
        hour_val.append(int(time[0:-2]))
        minute_val.append(int(time[-2:]))

model_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]

#------------------------------------------------------------
#read in obs time series data
obs_root_grp = Dataset('/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_O3_%s.nc'%(model_range))
obs_refs_dict = obs_root_grp.groups

obs_refs = []
obs_lats = []
obs_lons = []
obs_alt = []

for i in obs_refs_dict.keys():
    i = i.encode('ascii')
    obs_refs = np.append(obs_refs,i)

for ref in obs_refs:
    obs_ts_group = obs_root_grp.groups[ref] 
    obs_lats = np.append(obs_lats,obs_ts_group.latitude)
    obs_lons = np.append(obs_lons,obs_ts_group.longitude)
    obs_alt = np.append(obs_alt,obs_ts_group.altitude)
    obs_date = obs_ts_group.variables['date'][:]
    obs_time = obs_ts_group.variables['time'][:]
    
for i in range(len(obs_refs)):
    obs_refs[i] = obs_refs[i].lower()
    
#process obs dates and obs times to datetimes, then process pandas objects
year_val = []
month_val = []
day_val = []
hour_val = []
minute_val = []

#change obs times to datetimes
obs_date = obs_date.astype('str')
obs_time = obs_time.astype('str')

for date in obs_date:
    year_val.append(int(date[0:4]))
    month_val.append(int(date[4:6]))
    day_val.append(int(date[6:8]))

for time in obs_time:
    if np.float64(time) == 0:
        hour_val.append(0)
        minute_val.append(0)
    elif np.float64(time) == 30:
        hour_val.append(0)
        minute_val.append(30)
    else:
        hour_val.append(int(time[0:-2]))
        minute_val.append(int(time[-2:]))

obs_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]

#get obs lat_lon grid central points
obs_lats_centre, obs_lons_centre, model_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_lats,obs_lons) 

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

#--------------------------------------------------------
#load in periodic lsp data
root_grp_obs = Dataset('../obs/obs_sig_periods.nc')
root_grp_mod = Dataset('model_sig_periods.nc')

model_daily_mag = []
model_daily_phase = []
obs_daily_mag = []
obs_daily_phase = []
model_ave = []
obs_ave = []
model_seasonal_mag = []
model_seasonal_phase = []
obs_seasonal_mag = []
obs_seasonal_phase = []

for ref in obs_refs:
    site_group_obs = root_grp_obs.groups[ref]
    site_group_mod = root_grp_mod.groups[ref]
    
    model_daily_mag = np.append(model_daily_mag,site_group_mod.daily_amplitude)
    model_daily_phase = np.append(model_daily_phase,site_group_mod.daily_phase)
    model_ave = np.append(model_ave,site_group_mod.average)
    model_seasonal_mag = np.append(model_seasonal_mag,site_group_mod.seasonal_amplitude)
    model_seasonal_phase = np.append(model_seasonal_phase,site_group_mod.seasonal_phase)
    
    obs_daily_mag = np.append(obs_daily_mag,site_group_obs.daily_amplitude)
    obs_daily_phase = np.append(obs_daily_phase,site_group_obs.daily_phase)
    obs_ave = np.append(obs_ave,site_group_obs.average)
    obs_seasonal_mag = np.append(obs_seasonal_mag,site_group_obs.seasonal_amplitude)
    obs_seasonal_phase = np.append(obs_seasonal_phase,site_group_obs.seasonal_phase)

tags = modules.get_tags(obs_refs)
loc_dict = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic Sites','NA':'North America','SA':'South America'}

param = raw_input('\nmag, ph or ave?\n')

if param != 'ave':
    period_type = raw_input('d or s?\n')


if param == 'mag':
    if period_type == 'd':
        diff = model_daily_mag - obs_daily_mag
    if period_type == 's':
        diff = model_seasonal_mag - obs_seasonal_mag

if param == 'ph':
    if period_type == 'd':
        diff = model_daily_phase - obs_daily_phase
    if period_type == 's':
        diff = model_seasonal_phase - obs_seasonal_phase
        
if param == 'ave':
    diff = model_ave - obs_ave

#set up plot

latlower_setup = [24,30,20]
latupper_setup = [72,72,50]
lonwest_setup = [-170,-15,120]
loneast_setup = [-50,35,157]
label_out = ['us','eu','asia']

for x in range(3):
    fig, (ax) = plt.subplots(1,figsize=(23,12))
    fig.patch.set_facecolor('white')
 
    #setup basemap projection
    m = Basemap(projection='cyl',llcrnrlat=latlower_setup[x],urcrnrlat=latupper_setup[x],llcrnrlon=lonwest_setup[x],urcrnrlon=loneast_setup[x],resolution='c')


    m.drawcoastlines()
    m.drawmapboundary()
    #parallels = np.arange(-90,91,15)
    #meridians = np.arange(-180,181,30)
    #plt.xticks(meridians)
    #plt.yticks(parallels)
    #m.drawparallels(parallels)
    #m.drawmeridians(meridians)

    X,Y = m(obs_lons,obs_lats)


    m_size = 350

    for i in range(len(obs_lons)):
        if param == 'mag':
            max_diff = np.max(abs(diff))
            all = m.scatter(X[i],Y[i],c=diff[i], s=m_size, vmin = -max_diff,vmax = max_diff, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.coolwarm,zorder=10)
        if param == 'ph':
            if (period_type == 'd'):
                all = m.scatter(X[i],Y[i],c=diff[i], s=m_size, vmin = -12,vmax = 12, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.hsv,zorder=10)
            if (period_type == 'ha'):
                all = m.scatter(X[i],Y[i],c=diff[i], s=m_size, vmin = -3,vmax = 3, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.hsv,zorder=10)    
            if (period_type == 's') or (period_type == 'a') :
                all = m.scatter(X[i],Y[i],c=diff[i], s=m_size, vmin = -6,vmax = 6, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.hsv,zorder=10)
        if param == 'ave':
            max_diff = np.max(abs(diff))
            all = m.scatter(X[i],Y[i],c=diff[i], s=m_size, vmin = -max_diff,vmax = max_diff, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.coolwarm,zorder=10)
    t = m.plot(X,Y, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

    cb = plt.colorbar(all, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f',pad = 0.05)
    
    if (param == 'mag') or (param == 'ave') : 
        cb.set_label('Concentration (ppb)', fontsize = 26)    
    if param == 'ph':
        cb.set_label('Months', fontsize = 26)
    cb.ax.tick_params(labelsize=21)
    #plt.title('Observational average surface %s between %s'%(species,date_range),fontsize=20)
    #cb.ax.tick_params(labelsize=16)

    plt.tight_layout(pad = 3.08)

    mng = plt.get_current_fig_manager()
    mng.window.wm_geometry("+2500+2000")

    fig.canvas.draw()

    fig.canvas.mpl_connect('pick_event', onpick)

    #plt.savefig('/work/home/db876/plotting_tools/plots/AGU_plots/seasonal_%s_diff_%s.png'%(param,label_out[x]),format='png')

    plt.show()
