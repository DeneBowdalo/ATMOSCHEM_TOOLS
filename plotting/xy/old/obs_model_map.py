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

#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres = modules.get_obs_info(present_dir)
model_fname,species,start_year,end_year = modules.get_model_info(present_dir)

obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

root_obs_ts = Dataset(obs_fname)

model_root_grp = Dataset('/work/home/db876/global_grid/%s/%s_%s/%s/model_sig_periods.nc'%(species,start_year,end_year,last_dir))
lat_c = model_root_grp.variables['lat_centre'][:]
lon_c = model_root_grp.variables['lon_centre'][:]
lat_e = model_root_grp.variables['lat_edges'][:]
lon_e = model_root_grp.variables['lon_edges'][:]
model_daily_mag = model_root_grp.variables['daily_amplitude'][:]
model_daily_phase = model_root_grp.variables['daily_phase'][:]
model_ave = model_root_grp.variables['average'][:]
model_seasonal_mag = model_root_grp.variables['seasonal_amplitude'][:]
model_seasonal_phase = model_root_grp.variables['seasonal_phase'][:]

#--------------------------------------------------------
#load in periodic lsp data
root_grp_obs = Dataset('../obs_%s_%s/obs_sig_periods.nc'%(vres,timeres))

tags = modules.get_tags(obs_refs)

obs_daily_mag = []
obs_daily_phase = []
obs_ave = []
obs_seasonal_mag = []
obs_seasonal_phase = []

for ref in obs_refs:
    site_group_obs = root_grp_obs.groups[ref.lower()]

    obs_daily_mag = np.append(obs_daily_mag,site_group_obs.daily_amplitude)
    obs_daily_phase = np.append(obs_daily_phase,site_group_obs.daily_phase)
    obs_ave = np.append(obs_ave,site_group_obs.average)
    obs_seasonal_mag = np.append(obs_seasonal_mag,site_group_obs.seasonal_amplitude)
    obs_seasonal_phase = np.append(obs_seasonal_phase,site_group_obs.seasonal_phase)
    
n_boxes = len(lat_c)*len(lon_c)

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
        z_obs = obs_daily_mag
        z_model = model_daily_mag
        min = 0
        max = np.max(z_model)
    if period == 's':
        z_obs = obs_seasonal_mag
        z_model = model_seasonal_mag
        min = 0
        max = np.max(z_model)
        
        phase_type = 'Month'
        
    period = period.title()
        
if data_type == 'ph':
    if period == 'd':
        z_obs = obs_daily_phase
        z_model = model_daily_phase
        min = 0
        max = 24
        phase_type = 'Hour'
    if period == 's':
        min = 0
        max = 12
        z_obs = obs_seasonal_phase
        z_model = model_seasonal_phase
        
        #change obs SH seasonal phase to match up with NH
        #for i in range(len(obs_lats)):
        #    if obs_lats[i] < 0 :
        #        now = z_obs[i] - 6.
        #        if now < 0:
        #            now = 12. - np.abs(now)
        #        z_obs[i] = now
        
        #change model SH seasonal phase to match up with NH
        #out = z_model[:45,:] - 6.
        #test = out < 0
        #out[test] = 12. - np.abs(out[test])
        #z_model[:45,:] = out
        
        phase_type = 'Month'

if data_type == 'ave':
    z_obs = obs_ave
    z_model = model_ave
    min = np.min(z_model)
    max = np.max(z_model)
    
#--------------------------
#set up plot
fig =plt.figure(figsize=(11,9.5))
fig.patch.set_facecolor('white')
ax2 = plt.subplot2grid((4,2), (0,0), colspan=1)
ax3 = plt.subplot2grid((4,2), (0,1), colspan=1)
ax4 = plt.subplot2grid((4,2), (1,0), colspan=1)
ax5 = plt.subplot2grid((4,2), (1,1), colspan=1)
ax = plt.subplot2grid((4,2), (2,0), colspan=2, rowspan = 2)

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

#plot model absolute
#--------------------
if (data_type == 'amp') or (data_type == 'ave'):
    if data_type == 'amp':
        pl = m.pcolor(lon_e,lat_e,z_model, vmin=min, vmax=max,linewidth=0.5,cmap=plt.cm.coolwarm)
    if data_type == 'ave':
        pl = m.pcolor(lon_e,lat_e,z_model, vmin=min, vmax=max,linewidth=0.5,cmap=plt.cm.coolwarm)
 
if data_type == 'ph':
    if period == 'd':
        pl = m.pcolor(lon_e,lat_e,z_model, vmin=0, vmax=24,linewidth=0.5,cmap=plt.cm.hsv)
    if period == 's':
        pl = m.pcolor(lon_e,lat_e,z_model, vmin=0, vmax=12,linewidth=0.5,cmap=plt.cm.hsv)

#ax.set_xlabel('Longitude',fontsize = 24)
#ax.set_ylabel('Latitude',fontsize = 24)

#-------------------------------------------
#plot obs absolute
tags = modules.get_tags(obs_refs)

latlower_setup = [20,30,20,lat_e[0]]
latupper_setup = [80,72,55,lat_e[-1]]
lonwest_setup = [-170,-15,115,lon_e[0]]
loneast_setup = [-50,35,155,lon_e[-1]]
label_out = ['NA','EU','AS','zz']
label = ['NA','EU','AS','ROW']
size_dict = {'NA':30,'EU':30,'AS':40,'AF':60,'O':80,'OC':80,'SA':80,'ARC':80,'ANT':80}

obs_lons = np.array(obs_lons)
obs_lats = np.array(obs_lats)
tags = np.array(tags)

ax_list = [ax2,ax3,ax4,ax5]
count = 0 

for ax in ax_list:
    #setup basemap projection
    m = Basemap(projection='cyl',llcrnrlat=latlower_setup[count],urcrnrlat=latupper_setup[count],llcrnrlon=lonwest_setup[count],urcrnrlon=loneast_setup[count],resolution='c',ax = ax)

    m.drawcoastlines()
    m.drawmapboundary()
    
    if count == 3:
        test = []
        other_tags = ['AF','ANT','ARC','OC','O','SA']
        for i in tags:
            if i in other_tags:
                test.append(True)
            else:
                test.append(False)
        test = np.array(test)
    else:
        current_tag = label_out[count] 
        test = tags == current_tag
    
    
    current_z = z_obs[test]
    current_lons = obs_lons[test]
    current_lats = obs_lats[test]
    
    X,Y = m(current_lons,current_lats)

    if count == 0:
        m_size= 25
    if count == 1:
        m_size = 15
    if count == 2:
        m_size = 100
    if count == 3:
        m_size = 100

    for i in range(len(current_lons)):
        if data_type == 'amp':
            all = m.scatter(X[i],Y[i],c=current_z[i], s=m_size, vmin = min,vmax = max, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.coolwarm,zorder=10)
        if data_type == 'ph':
            if (period == 'd'):
                max_diff = 12
                all = m.scatter(X[i],Y[i],c=current_z[i], s=m_size, vmin = 0,vmax = 24, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.hsv,zorder=10)  
            if (period == 's') :
                max_diff = 6
                all = m.scatter(X[i],Y[i],c=current_z[i], s=m_size, vmin = 0,vmax = 12, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.hsv,zorder=10)
        if data_type == 'ave':
            all = m.scatter(X[i],Y[i],c=current_z[i], s=m_size, vmin = min,vmax = max, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.coolwarm,zorder=10)
    
    count+=1

plt.tight_layout(pad = 1.08)
fig.subplots_adjust(wspace=0.01)
fig.subplots_adjust(hspace=0.01)

#plot colorbar
fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.82, 0.13, 0.05, 0.7])

if (data_type == 'amp') or (data_type == 'ave'):
    cb = fig.colorbar(pl, cax=cbar_ax)
    cb.set_label('Concentration (ppb)', fontsize = 30)

elif (data_type == 'ph'):
    if period == 'd':
        cb = fig.colorbar(pl, cax=cbar_ax,ticks=[0,2,4,6,8,10,12,14,16,18,20,22])
        cb.ax.set_yticklabels(['0','2','4','6','8','10','12','14','16','18','20','22'])
        cb.set_label('Hour', fontsize = 30)

    if period == 's':
        month_lengths = [0,31,28.25,31,30,31,30,31,31,30,31,30]
        annual_ratio = 12./365.25
        current_days = 0
        month_ticks = []
        for i in month_lengths:
            current_days = current_days+i
            month_ticks.append(current_days * annual_ratio)
        cb = fig.colorbar(pl, cax=cbar_ax,ticks=month_ticks)
        cb.ax.set_yticklabels(['1','2','3','4','5','6','7','8','9','10','11','12'])
        cb.set_label('Month', fontsize = 30)

else:
    cb = fig.colorbar(pl, cax=cbar_ax)

cb.ax.tick_params(labelsize=24)
ax.tick_params(axis='both', which='major', labelsize=30)
ax.tick_params(axis='both', which='minor', labelsize=30)


plt.annotate('NA',xy=(-15.8,0.90), xycoords='axes fraction', fontsize=30)
plt.annotate('AS',xy=(-11,0.54), xycoords='axes fraction', fontsize=30)
plt.annotate('EU',xy=(-6.5,1.16), xycoords='axes fraction', fontsize=30)
plt.annotate('ROW',xy=(-8.0,0.55), xycoords='axes fraction', fontsize=30)


for tick in ax.get_xaxis().get_major_ticks():
    tick.set_pad(11.)

plt.show()
