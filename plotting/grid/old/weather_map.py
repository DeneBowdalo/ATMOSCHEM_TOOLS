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

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

#get species from current directory
present_dir = os.getcwd()

model_fname,species,start_year,end_year = modules.get_model_info(present_dir)

def grid_interactive(event):
    global pl 
    global ind
    global fig2
    
    ind = event.ind
    ind = ind[0]

    #find ind of closest lat/lon
    #ind = modules.find_nearest_point_index(obs_lons,obs_lats,x_data,y_data)
    
    try:
        for i in range(len(pl)):
            pl.pop(0).remove()
            first_run = False  
        
    except:
        first_run = True
        pass
    
    pl = m.plot([linear_lons[ind]], [linear_lats[ind]], 's', ms=20, alpha=0.6, color='yellow',zorder=20)

    #get model spectra for site clicked
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,linear_lats[ind],linear_lons[ind])
    period = site_group.variables['period'][:,lat_n,lon_n]
    spectrum = site_group.variables['amplitude'][:,lat_n,lon_n]
    #period = periods[:,lat_n,lon_n]
    #spectrum = amplitude[:,lat_n,lon_n]
    
    fig.canvas.draw()
        
    if first_run == False:
        plt.close(fig2)
        fig2, (ax1) = plt.subplots(1, 1, figsize =(24,12))
        fig2.patch.set_facecolor('white')

        ax1.plot(period, spectrum, color='black', markersize = 3)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        
        #plt.legend(loc = 'lower right')
        plt.tight_layout()
        ax1.grid()
        
        plt.show()
    else:
        fig2, (ax1) = plt.subplots(1, 1, figsize =(24,12))
        fig2.patch.set_facecolor('white')
        
        ax1.plot(period, spectrum, color='black', markersize = 3)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        
        #plt.legend(loc = 'lower right')
        plt.tight_layout()
        ax1.grid()
        
        plt.show()


#--------------------------------------------
#
site_group = Dataset('model_spectra.nc')

#read in variables for site
#period = site_group.variables['period'][0]
#print period.shape
#period = period[:,0,0]


#amplitude = site_group.variables['amplitude'][:]
#phase = site_group.variables['phase'][:]
tur_grad = site_group.variables['turbulent_gradient'][:]
syn_grad = site_group.variables['synoptic_gradient'][:]
mw_grad = site_group.variables['macroweather_gradient'][:]
ts_bp = site_group.variables['turbulence_synoptic_transition'][:]
smw_bp = site_group.variables['synoptic_macroweather_transition'][:]
lat_e = site_group.variables['lat_edges'][:]
lon_e = site_group.variables['lon_edges'][:]
lat_c = site_group.variables['lat_centre'][:]
lon_c = site_group.variables['lon_centre'][:]

n_boxes = len(lon_c)*len(lat_c)


data_type = raw_input('grad, grad_f, bp or p?\n')

if data_type == 'grad':
    period = raw_input('\nt, s or mw?\n')
    if period == 't':
        z = tur_grad
    if period == 's':
        z = syn_grad
    if period == 'mw':   
        z = mw_grad
        
if data_type == 'grad_f':
    period = raw_input('\nt, s or mw?\n')
    if period == 't':
        z = tur_grad_f
    if period == 's':
        z = syn_grad_f
    if period == 'mw':   
        z = mw_grad_f
        
if data_type == 'bp':
    period = raw_input('\nts or smw?\n')
    if period == 'ts':
        z = ts_bp
    if period == 'smw':
        z = smw_bp
        
if data_type == 'p':
    period = raw_input('\nt, s or mw?\n')
    if period == 't':
        z = p1
    if period == 's':
        z = p2
    if period == 'mw':   
        z = p3


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
plt.xticks(meridians)
plt.yticks(parallels)
m.drawparallels(parallels)
m.drawmeridians(meridians)

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

if data_type == 'grad':
    if period == 't':
        min = np.min(z)
        max = np.max(z)
    if period == 's':
        min = np.min(z)
        max = np.max(z)
    if period == 'mw':
        min = np.min(z)
        max = np.max(z)
    pl = m.pcolor(lon_e,lat_e,z, vmin=min, vmax=max,linewidth=0.5,picker = 5)
    cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
    cb.set_label('Gradient', fontsize = 24)

if data_type == 'grad_f':
    if period == 't':
        min = np.min(z)
        max = np.max(z)
    if period == 's':
        min = np.min(z)
        max = np.max(z)
    if period == 'mw':
        min = np.min(z)
        max = np.max(z)
    pl = m.pcolor(lon_e,lat_e,z, vmin=min, vmax=max,linewidth=0.5,picker = 5)
    cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
    cb.set_label('Gradient', fontsize = 24)
    
if data_type == 'bp':
    if period == 'ts':
        min = 0.2
        max = 1
    if period == 'smw':
        min = 5
        max = 25
    pl = m.pcolor(lon_e,lat_e,z, vmin=min, vmax=max,linewidth=0.5,picker = 5)
    cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
    cb.set_label('Transtion Point', fontsize = 24)

if data_type == 'p':
    if period == 't':
        min = np.min(z)
        max = np.max(z)
    if period == 's':
        min = np.min(z)
        max = np.max(z)
    if period == 'mw':
        min = np.min(z)
        max = np.max(z)
    pl = m.pcolor(lon_e,lat_e,z, vmin=min, vmax=max,linewidth=0.5,picker = 5)
    cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
    cb.set_label('Average Amplitude (ppb)', fontsize = 24)

plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)
#plt.title('%s Day %s for Surface %s between %s_%s'%(period_chose,type_label,species,start_year,end_year),fontsize=20)
cb.ax.tick_params(labelsize=16)

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500+100")

fig.canvas.mpl_connect('pick_event', grid_interactive)

plt.show()
