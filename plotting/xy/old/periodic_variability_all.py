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

typ = raw_input('\nobs or model?\n')

def interactive(event):
    modules.clicker_interactive_map_obsmodel(event,species,lat_e,lon_e,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_daily_waveform,obs_seasonal_waveform,obs_full_waveform,model_daily_waveform,model_seasonal_waveform,model_full_waveform,fig,all_m)

#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres = modules.get_obs_info(present_dir)
model_fname,species,start_year,end_year = modules.get_model_info(present_dir)

obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

#get obs lat_lon grid central points
obs_lats_centre, obs_lons_centre, model_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_lats,obs_lons) 

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

#--------------------------------------------------------
#load in periodic lsp data
obs_period_grp = Dataset('../obs_%s_%s/obs_sig_periods.nc'%(vres,timeres))
model_period_grp = Dataset('model_sig_periods.nc')

model_daily_mag = []
model_daily_phase = []
model_daily_phase_min = []
obs_daily_mag = []
obs_daily_phase = []
obs_daily_phase_min = []
model_ave = []
obs_ave = []
model_seasonal_mag = []
model_seasonal_phase = []
model_seasonal_phase_min = []
obs_seasonal_mag = []
obs_seasonal_phase = []
obs_seasonal_phase_min = [] 
model_daily_ff = []
model_seasonal_ff = []
obs_daily_ff = []
obs_seasonal_ff = []
obs_daily_waveform = []
model_daily_waveform = []
obs_seasonal_waveform = []                                                                                                                                                                                                                       
model_seasonal_waveform = []
obs_full_waveform = []                                                                                                                                                                                                                       
model_full_waveform = []

for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    site_group_mod = model_period_grp.groups[ref]
    
    model_daily_mag = np.append(model_daily_mag,site_group_mod.daily_amplitude)
    model_daily_phase = np.append(model_daily_phase,site_group_mod.daily_phase)
    model_daily_phase_min = np.append(model_daily_phase_min,site_group_mod.daily_phase_min)
    model_ave = np.append(model_ave,site_group_mod.average)
    model_seasonal_mag = np.append(model_seasonal_mag,site_group_mod.seasonal_amplitude)
    model_seasonal_phase = np.append(model_seasonal_phase,site_group_mod.seasonal_phase)
    model_seasonal_phase_min = np.append(model_seasonal_phase_min,site_group_mod.seasonal_phase_min)   
    model_daily_ff = np.append(model_daily_ff,site_group_mod.daily_ff) 
    model_seasonal_ff = np.append(model_seasonal_ff,site_group_mod.seasonal_ff)
    model_daily_waveform.append(site_group_mod.variables['daily_waveform'][:])
    model_seasonal_waveform.append(site_group_mod.variables['seasonal_waveform'][:]) 
    model_full_waveform.append(site_group_mod.variables['all_waveform'][:]) 
 
    obs_daily_mag = np.append(obs_daily_mag,site_group_obs.daily_amplitude)
    obs_daily_phase = np.append(obs_daily_phase,site_group_obs.daily_phase)
    obs_daily_phase_min = np.append(obs_daily_phase_min,site_group_obs.daily_phase_min)
    obs_ave = np.append(obs_ave,site_group_obs.average)
    obs_seasonal_mag = np.append(obs_seasonal_mag,site_group_obs.seasonal_amplitude)
    obs_seasonal_phase = np.append(obs_seasonal_phase,site_group_obs.seasonal_phase)
    obs_seasonal_phase_min = np.append(obs_seasonal_phase_min,site_group_obs.seasonal_phase_min)
    obs_daily_ff = np.append(obs_daily_ff,site_group_obs.daily_ff) 
    obs_seasonal_ff = np.append(obs_seasonal_ff,site_group_obs.seasonal_ff)
    obs_daily_waveform.append(site_group_obs.variables['daily_waveform'][:])
    obs_seasonal_waveform.append(site_group_obs.variables['seasonal_waveform'][:]) 
    obs_full_waveform.append(site_group_obs.variables['all_waveform'][:])

#calculate periodic variance
n_sites = len(obs_refs)
year_range = int(end_year) - int(start_year)
d0 = datetime.date(int(start_year), 1, 1)
d1 = datetime.date(int(end_year), 1, 1)
delta = d1 - d0
n_days = delta.days
n_years = year_range
n_hours = n_days*24
obs_daily_waveforms_long = []
obs_seasonal_waveforms_long = []
model_daily_waveforms_long = []
model_seasonal_waveforms_long = []

obs_waveform_len = len(obs_full_waveform[0])
model_waveform_len = len(model_full_waveform[0])

for i in range(n_sites):
    obs_daily_waveforms_long.append(np.repeat(obs_daily_waveform[i],n_days))
    mod_seasonal_obs = np.repeat(obs_seasonal_waveform[i],n_years)
    obs_seasonal_waveforms_long.append(mod_seasonal_obs[:obs_waveform_len])
    
    model_daily_waveforms_long.append(np.repeat(model_daily_waveform[i],n_days))
    mod_seasonal_model = np.repeat(model_seasonal_waveform[i],n_years)
    model_seasonal_waveforms_long.append(mod_seasonal_model[:model_waveform_len])
    
obs_daily_waveforms = obs_daily_waveforms_long[:]
obs_seasonal_waveforms = obs_seasonal_waveforms_long[:]
model_daily_waveforms = model_daily_waveforms_long[:]
model_seasonal_waveforms = model_seasonal_waveforms_long[:]

total_obs_var = []
diurnal_obs_var = []
seasonal_obs_var = []
full_obs_var = []
total_model_var = []
diurnal_model_var = []
seasonal_model_var = []
full_model_var = []

for i in range(len(obs_var)):
    total_obs_var.append(np.var(obs_var[i]))
    total_model_var.append(np.var(model_var[i]))
    
    d_w_obs = np.array(obs_daily_waveforms[i])[tests[i]]
    s_w_obs = np.array(obs_seasonal_waveforms[i])[tests[i]]
    f_w_obs = np.array(obs_full_waveform[i])[tests[i]]
    d_w_model = np.array(model_daily_waveforms[i])[tests[i]]
    s_w_model = np.array(model_seasonal_waveforms[i])[tests[i]]
    f_w_model = np.array(model_full_waveform[i])[tests[i]]
    
    diurnal_obs_var.append(np.var(d_w_obs))
    seasonal_obs_var.append(np.var(s_w_obs))
    full_obs_var.append(np.var(f_w_obs))
    diurnal_model_var.append(np.var(d_w_model))
    seasonal_model_var.append(np.var(s_w_model))
    full_model_var.append(np.var(f_w_model))
    
total_obs_var = np.array(total_obs_var)
diurnal_obs_var = np.array(diurnal_obs_var)
seasonal_obs_var = np.array(seasonal_obs_var)
full_obs_var = np.array(full_obs_var)
total_model_var = np.array(total_model_var)
diurnal_model_var = np.array(diurnal_model_var)
seasonal_model_var = np.array(seasonal_model_var)
full_model_var = np.array(full_model_var)

ratio_obs = 100./total_obs_var
ratio_model = 100./total_model_var

diurnal_obs_frac = ratio_obs*diurnal_obs_var 
seasonal_obs_frac = ratio_obs*seasonal_obs_var 
all_obs_frac = ratio_obs*full_obs_var
diurnal_model_frac = ratio_model*diurnal_model_var 
seasonal_model_frac = ratio_model*seasonal_model_var 
all_model_frac = ratio_model*full_model_var

if typ == 'obs':
    diff1 = np.copy(diurnal_obs_frac)
    diff2 = np.copy(seasonal_obs_frac)
elif typ == 'model':
    diff1 = np.copy(diurnal_model_frac)
    diff2 = np.copy(seasonal_model_frac)

#-------------------------
#set up plot

latlower_setup = [20,30,20,-90]
latupper_setup = [80,72,55,90]
lonwest_setup = [-170,-15,115,-181.25]
loneast_setup = [-50,35,155,178.75]
label_out = ['NA','EU','AS','ZZ']
label = ['NA','EU','AS','ROW']

obs_lons = np.array(obs_lons)
obs_lats = np.array(obs_lats)
tags = np.array(tags)

#--------------------------------
fig =plt.figure(figsize=(11.5,12))
fig.patch.set_facecolor('white')

ax = plt.subplot2grid((4,2), (0,0), colspan=2, rowspan=4)
ax.axhline(y=0.5,color='red')

ax1 = plt.subplot2grid((4,2), (0,0), colspan=1)
ax2 = plt.subplot2grid((4,2), (0,1), colspan=1)
ax3 = plt.subplot2grid((4,2), (1,0), colspan=1)
ax4 = plt.subplot2grid((4,2), (1,1), colspan=1)
ax5 = plt.subplot2grid((4,2), (2,0), colspan=1)
ax6 = plt.subplot2grid((4,2), (2,1), colspan=1)
ax7 = plt.subplot2grid((4,2), (3,0), colspan=1)
ax8 = plt.subplot2grid((4,2), (3,1), colspan=1)

fig.patch.set_facecolor('white')

max_diff = 100
min_diff = 0

#first plot
#---------------------------------------
count = 0

ax_list = [ax1,ax2,ax3,ax4]

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
    
    
    current_diff = diff1[test]
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
        #max_diff = np.max(abs(current_diff))
        all1 = m.scatter(X[i],Y[i],c=current_diff[i], s=m_size, vmin = min_diff,vmax = max_diff, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.rainbow,zorder=10)
           
    if (count == 0) or (count == 3):
        ax.text(0.03, 0.18, label[count], transform=ax.transAxes,fontsize=30, va='top')
    if (count == 1) or (count == 2):
        ax.text(0.03, 0.97, label[count], transform=ax.transAxes,fontsize=30, va='top')

    count+=1

#second plot
#---------------------------------------
count = 0

ax_list = [ax5,ax6,ax7,ax8]

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
    
    
    current_diff = diff2[test]
    current_lons = obs_lons[test]
    current_lats = obs_lats[test]
    
    X,Y = m(current_lons,current_lats)

    if count == 0:
        m_size= 30
    if count == 1:
        m_size = 15
    if count == 2:
        m_size = 100
    if count == 3:
        m_size = 100

    for i in range(len(current_lons)):
        all2 = m.scatter(X[i],Y[i],c=current_diff[i], s=m_size, vmin =  min_diff,vmax = max_diff, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.rainbow,zorder=10)
    
    if (count == 0) or (count == 3):
        ax.text(0.03, 0.18, label[count], transform=ax.transAxes,fontsize=30, va='top')
    if (count == 1) or (count == 2):
        ax.text(0.03, 0.97, label[count], transform=ax.transAxes,fontsize=30, va='top')

    count+=1

plt.tight_layout(pad = 1.08)
fig.subplots_adjust(wspace=0.01)
fig.subplots_adjust(hspace=0.01)

#plot colorbar
fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.82, 0.52, 0.05, 0.35])
cbar_ax2 = fig.add_axes([0.82, 0.02, 0.05, 0.35])

t = np.arange(0,101,20)

cb = fig.colorbar(all1, cax=cbar_ax,ticks=t)
cb2 = fig.colorbar(all2, cax=cbar_ax2,ticks=t)

cb.ax.tick_params(labelsize=24)
cb2.ax.tick_params(labelsize=24)

cb.set_label('% of Total Variance', fontsize = 25)
cb2.set_label('% of Total Variance', fontsize = 25)

plt.annotate('Diurnal',(-0.8,2.7),xycoords='axes fraction',fontsize=30, va='top')
plt.annotate('Seasonal',(-0.8,1.3),xycoords='axes fraction',fontsize=30, va='top')

plt.show()

