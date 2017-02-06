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

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

#get species from current directory
present_dir = os.getcwd()

model_fname,species,start_year,end_year = modules.get_model_info(present_dir)

def interactive(event):
    all_periods = np.array([daily_h3_amp,daily_h2_amp,daily_h1_amp,orig_daily_amp,daily_amp,daily_h3_ph,daily_h2_ph,daily_h1_ph,orig_daily_ph,daily_ph,seasonal_h3_amp,seasonal_h2_amp,seasonal_h1_amp,annual_amp,seasonal_amp,seasonal_h3_ph,seasonal_h2_ph,seasonal_h1_ph,annual_ph,seasonal_ph,ave])
    modules.clicker_interactive_map_model(event,plot_type_2,lat_e,lon_e,linear_lats,linear_lons,date,time,datetimes,start_year,all_periods,var,fig,m)
    
#-------------------------------------------------------
#read in model time series data
    
root_grp = Dataset(model_fname)
var = root_grp.variables[species.lower()][:]
var = var*1e9
date = root_grp.variables['date'][:]
time = root_grp.variables['time'][:]
lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
grid_size = root_grp.variables['grid_size'][:]
grid_size = grid_size[0]

n_boxes = len(lon_c)*len(lat_c)

#process model dates and model times to datetimes, then process pandas objects

year_val = []
month_val = []
day_val = []
hour_val = []
minute_val = []

#change model times to datetimes
date = date.astype('str')
time = time.astype('str')

for d in date:
    year_val.append(int(d[0:4]))
    month_val.append(int(d[4:6]))
    day_val.append(int(d[6:8]))

for t in time:
    if np.float64(t) == 0:
        hour_val.append(0)
        minute_val.append(0)
    elif np.float64(t) == 30:
        hour_val.append(0)
        minute_val.append(30)
    else:
        hour_val.append(int(t[0:-2]))
        minute_val.append(int(t[-2:]))

datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]

#--------------------------------------------
#
site_group = Dataset('model_sig_periods.nc')

daily_h3_amp = site_group.variables['daily_harmonic3_amplitude'][:]
daily_h2_amp = site_group.variables['daily_harmonic2_amplitude'][:]
daily_h1_amp = site_group.variables['daily_harmonic1_amplitude'][:]
orig_daily_amp = site_group.variables['original_daily_amplitude'][:]
daily_amp = site_group.variables['daily_amplitude'][:]
daily_h3_ph = site_group.variables['daily_harmonic3_phase'][:]
daily_h2_ph = site_group.variables['daily_harmonic2_phase'][:]
daily_h1_ph = site_group.variables['daily_harmonic1_phase'][:]
orig_daily_ph = site_group.variables['original_daily_phase'][:]
daily_ph = site_group.variables['daily_phase'][:]
seasonal_h3_amp = site_group.variables['seasonal_harmonic3_amplitude'][:]
seasonal_h2_amp = site_group.variables['seasonal_harmonic2_amplitude'][:]
seasonal_h1_amp = site_group.variables['seasonal_harmonic1_amplitude'][:]
annual_amp = site_group.variables['annual_amplitude'][:]
seasonal_amp = site_group.variables['seasonal_amplitude'][:]
seasonal_h3_ph = site_group.variables['seasonal_harmonic3_phase'][:]
seasonal_h2_ph = site_group.variables['seasonal_harmonic2_phase'][:]
seasonal_h1_ph = site_group.variables['seasonal_harmonic1_phase'][:]
annual_ph = site_group.variables['annual_phase'][:]
seasonal_ph = site_group.variables['seasonal_phase'][:]
ave = site_group.variables['average'][:]

lat_e = site_group.variables['lat_edges'][:]
lon_e = site_group.variables['lon_edges'][:]
lat_c = site_group.variables['lat_centre'][:]
lon_c = site_group.variables['lon_centre'][:]

#amp or phase
data_type = raw_input('\namp, ph or ave?\n')
if data_type == 'amp':
    type_label = 'Amplitude'
if data_type == 'ph':
    type_label = 'Phase'
if data_type == 'ave':
    type_label = 'Average'

if (data_type != 'ave'):
    period = raw_input('\nd or s?\n')

if data_type == 'amp':
    if period == 'd':
        z = daily_amp
    if period == 's':
        z = seasonal_amp
        phase_type = 'Month'
        
    period = period.title()
        
if data_type == 'ph':
    if period == 'd':
        z = daily_ph
        phase_type = 'Hour'
        ph_min = 0
        ph_max = 24
    if period == 's':
        z = seasonal_ph
        phase_type = 'Month'
        ph_min = 0
        ph_max = 12

    period = period.title()

    #fundamental_min = np.min(annual_amp)
    #annual_amp = annual_amp-all_min
    #fundamental_max = np.max(annual_amp)
    #color_ratio = int(255./all_max)
    #r_colors = np.empty((len(lat_c),len(lon_c)))
    #for i in range(len(lat_c)):
    #    for j in range(len(lon_c)):
    #        r_colors[i,j] = annual_amp[i,j]*color_ratio

    #harmonic1_min = np.min(seasonal_h1_amp)
    #seasonal_h1_amp = seasonal_h1_amp-all_min
    #harmonic1_max = np.max(seasonal_h1_amp)
    #color_ratio = int(255./all_max)
    #g_colors = np.empty((len(lat_c),len(lon_c)))
    #for i in range(len(lat_c)):
    #    for j in range(len(lon_c)):
    #        g_colors[i,j] = seasonal_h1_amp[i,j]*color_ratio

    #harmonic2_min = np.min(seasonal_h2_amp)
    #seasonal_h2_amp = seasonal_h2_amp-all_min 
    #harmonic2_max = np.max(seasonal_h2_amp) 
    #color_ratio = int(255./all_max)
    #b_colors = np.empty((len(lat_c),len(lon_c)))
    #for i in range(len(lat_c)):
    #    for j in range(len(lon_c)):
    #        b_colors[i,j] = seasonal_h2_amp[i,j]*color_ratio    

    #rgb_mix = np.empty((len(lat_c),len(lon_c),3))
    #for i in range(len(lat_c)):
    #    for j in range(len(lon_c)):
    #        rgb_mix[i,j,0] = r_colors[i,j]
    #        rgb_mix[i,j,1] = g_colors[i,j]
    #        rgb_mix[i,j,2] = b_colors[i,j]
    
#z = percent_amps[:,:,0]

#for i in range(len(lon_c)):
    #print percent_amps[0,i,:]
      
    #harmonic3_min = np.min(seasonal_h3_amp)
    #harmonic3_max = np.max(seasonal_h3_amp) 

    #r_type = raw_input('\nh1, h2 or h3?\n')
    #if period == 'seasonal':
    #    if r_type == 'h1':
    #        z = seasonal_h1_amp/annual_amp
    #    if r_type == 'h2':
    #        z = seasonal_h2_amp/annual_amp
    #    if r_type == 'h3':
    #        z = seasonal_h3_amp/annual_amp
    #if period == 'daily':
    #    if r_type == 'h1':
    #        z = daily_h1_amp/orig_daily_amp 
    #    if r_type == 'h2':
    #        z = daily_h2_amp/orig_daily_amp 
    #    if r_type == 'h3':
    #        z = daily_h3_amp/orig_daily_amp


if data_type == 'ave':
    z = ave
   
#z[:46,:] = z[:46,:]-6
#test = z < 0
#z[test] =  12 - np.abs(z[test])

#------------------------------------------
#set up plot
fig =plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')
 
#setup basemap projection
m = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c')

m.drawcoastlines()
m.drawmapboundary()
parallels = np.arange(-90,91,15)
meridians = np.arange(-180,151,30)
#plt.xticks(meridians)
#plt.yticks(parallels)
#m.drawparallels(parallels)
#m.drawmeridians(meridians)

linear_lats = []
linear_lons=[]
lat_i = 0
lon_i = 0

for i in range(n_boxes):
    current_lat = lat_c[lat_i]
    current_lon = lon_c[lon_i]

    linear_lats.append(current_lat)
    linear_lons.append(current_lon)
    
    if lon_i == (len(lon_c)-1):
        lat_i+=1
        lon_i=0
    else:
        lon_i+=1


if (data_type == 'amp') or (data_type == 'ave'):
    pl = m.pcolor(lon_e,lat_e,z, vmin=np.min(z), vmax=np.max(z),linewidth=0.5,cmap=plt.cm.coolwarm,picker = 5)
    cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
    cb.set_label('Concentration (ppb)', fontsize = 16)    
if data_type == 'ph':
    pl = m.pcolor(lon_e,lat_e,z, vmin=ph_min, vmax=ph_max,linewidth=0.5,cmap=plt.cm.hsv,picker = 5)
    cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
    cb.set_label('%s'%(phase_type), fontsize = 16)

plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)

if (data_type != 'ave'):
    plt.title('%s %s for Surface %s between %s-%s'%(period,type_label,species,start_year,end_year),fontsize=20)
else:
    plt.title('%s for Surface %s between %s-%s'%(type_label,species,start_year,end_year),fontsize=20)

cb.ax.tick_params(labelsize=16)

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500+100")


fig.canvas.mpl_connect('pick_event', interactive)

plt.show()
