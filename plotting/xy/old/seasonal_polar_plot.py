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
import matplotlib.lines as mlines
import operator
from collections import OrderedDict

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres = modules.get_obs_info(present_dir)

obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)

#--------------------------------------------------------
#load in periodic lsp data
root_grp_obs = Dataset('../obs_%s_%s/obs_sig_periods.nc'%(vres,timeres))
root_grp_model = Dataset('model_sig_periods.nc')

tags = modules.get_tags(obs_refs)

obs_daily_mag = []
obs_daily_phase = []
obs_ave = []
obs_seasonal_mag = []
obs_seasonal_phase = []

model_daily_mag = []
model_daily_phase = []
model_ave = []
model_seasonal_mag = []
model_seasonal_phase = []


for ref in obs_refs:
    site_group_obs = root_grp_obs.groups[ref.lower()]
    site_group_model = root_grp_model.groups[ref.lower()]
    
    obs_daily_mag = np.append(obs_daily_mag,site_group_obs.daily_amplitude)
    obs_daily_phase = np.append(obs_daily_phase,site_group_obs.daily_phase)
    obs_ave = np.append(obs_ave,site_group_obs.average)
    obs_seasonal_mag = np.append(obs_seasonal_mag,site_group_obs.seasonal_amplitude)
    obs_seasonal_phase = np.append(obs_seasonal_phase,site_group_obs.seasonal_phase)
    
    model_daily_mag = np.append(model_daily_mag,site_group_model.daily_amplitude)
    model_daily_phase = np.append(model_daily_phase,site_group_model.daily_phase)
    model_ave = np.append(model_ave,site_group_model.average)
    model_seasonal_mag = np.append(model_seasonal_mag,site_group_model.seasonal_amplitude)
    model_seasonal_phase = np.append(model_seasonal_phase,site_group_model.seasonal_phase)

    
        
#change obs SH seasonal phase to match up with NH
for i in range(len(obs_lats)):
    if obs_lats[i] < 0 :
        now_obs = obs_seasonal_phase[i] - 6.
        now_model = model_seasonal_phase[i] - 6.
        
        if now_obs < 0:
            now_obs = 12. - np.abs(now_obs)
            obs_seasonal_phase[i] = now_obs
            
        if now_model < 0:
            now_model = 12. - np.abs(now_model)
            model_seasonal_phase[i] = now_model
        
    
        
#------------------------------------------
tags = modules.get_tags(obs_refs)


loc_colors = {'ANT':'magenta','ARC':'purple','AF':'orange','AS':'blue','EU':'red','OC':'black','O':'yellow','NA':'green','SA':'pink'}
all_locs = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic Sites','NA':'North America','SA':'South America'}
order_dict = {'ANT':8,'ARC':9,'AF':6,'AS':3,'EU':2,'OC':5,'O':4,'NA':1,'SA':7}

order_tags = []
for tag in tags:
    n = order_dict[tag]
    order_tags.append(n)

tags_argsort = np.argsort(order_tags)

tags = tags[tags_argsort]
obs_seasonal_mag = obs_seasonal_mag[tags_argsort]
obs_seasonal_phase = obs_seasonal_phase[tags_argsort]
model_seasonal_mag = model_seasonal_mag[tags_argsort]
model_seasonal_phase = model_seasonal_phase[tags_argsort]

print np.max(model_seasonal_phase)

tag_colors = [loc_colors[i] for i in tags]
tag_locs = [all_locs[i] for i in tags]

#set up plot
fig =plt.figure(figsize=(17,12))
fig.patch.set_facecolor('white')

pi2 = np.pi*2.
ratio = pi2/12.
obs_seasonal_phase_rads = [i*ratio for i in obs_seasonal_phase]
model_seasonal_phase_rads = [i*ratio for i in model_seasonal_phase]

ax1 = fig.add_subplot(121, polar=True)
ax1.set_theta_zero_location('N')
ax1.set_theta_direction(-1) 
for i in range(len(obs_seasonal_mag)):
    ax1.scatter(obs_seasonal_phase_rads[i],obs_seasonal_mag[i],c=tag_colors[i],s=60,label = tag_locs[i])

month_lengths = [0,31,28.25,31,30,31,30,31,31,30,31,30]
pi2 = np.pi*2.
ratio = pi2/365.25
current_days = 0
deg_months = []
for i in month_lengths:
    current_days = current_days + i
    val = current_days*ratio
    deg_months.append(val)
ax1.set_xticks(deg_months)
ax1.set_xticklabels(['1','2','3','4','5','6','7','8','9','10','11','12'])
ax1.set_rmax(np.max(model_seasonal_mag)+2)
ax1.set_title("Observations", fontsize=25,y=1.10,weight='bold')

ax1.tick_params(axis='both', which='major', labelsize=22)

ax2 = fig.add_subplot(122, polar=True)
ax2.set_theta_zero_location('N')
ax2.set_theta_direction(-1) 
for i in range(len(model_seasonal_mag)):
    ax2.scatter(model_seasonal_phase_rads[i],model_seasonal_mag[i],c=tag_colors[i],s=60,label=tag_locs[i])
ax2.set_xticks(deg_months)
ax2.set_xticklabels(['1','2','3','4','5','6','7','8','9','10','11','12'])
ax2.set_rmax(np.max(model_seasonal_mag)+2)
ax2.set_title("Model", fontsize=25,y=1.10,weight='bold')



ax2.tick_params(axis='both', which='major', labelsize=22)

ax1.text(0.30, 0.90, 'Amplitude (ppb)', transform=ax1.transAxes,fontsize=20, va='top')
ax1.text(0.85, 0.90, 'Phase (Month)', transform=ax1.transAxes,fontsize=20, va='top')


plt.tight_layout(pad = 3.08)

all_handles, all_labels = ax1.get_legend_handles_labels()
hl = sorted(zip(all_handles, all_labels), key=operator.itemgetter(1))
all_handles, all_labels = zip(*hl)
by_label = OrderedDict(zip(all_labels, all_handles))
plt.legend(by_label.values(),by_label.keys(), loc='lower left',prop={'size':22},fancybox=True,ncol=4,markerscale=2,bbox_to_anchor=(-1.0,-0.25))


#plt.savefig('polar_seasonal_distribution.png')
plt.show()
