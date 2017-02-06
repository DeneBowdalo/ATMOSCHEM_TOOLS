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
obs_period_grp = Dataset('../obs_%s_%s/obs_sig_periods.nc'%(vres,timeres))
model_period_grp = Dataset('model_sig_periods.nc')

tags = modules.get_tags(obs_refs)

obs_diurnal_mag = []
obs_diurnal_phase = []

model_diurnal_mag = []
model_diurnal_phase = []


for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref.lower()]
    site_group_model = model_period_grp.groups[ref.lower()]
    
    obs_diurnal_mag = np.append(obs_diurnal_mag,site_group_obs.daily_amplitude)
    obs_diurnal_phase = np.append(obs_diurnal_phase,site_group_obs.daily_phase)
    
    model_diurnal_mag = np.append(model_diurnal_mag,site_group_model.daily_amplitude)
    model_diurnal_phase = np.append(model_diurnal_phase,site_group_model.daily_phase)

#------------------------------------------

loc_colors = {'ANT':'magenta','ARC':'purple','AF':'orange','AS':'blue','EU':'red','OC':'black','O':'yellow','NA':'green','SA':'pink'}
all_locs = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic Sites','NA':'North America','SA':'South America'}
order_dict = {'ANT':8,'ARC':9,'AF':6,'AS':3,'EU':2,'OC':5,'O':4,'NA':1,'SA':7}

order_tags = []
for tag in tags:
    n = order_dict[tag]
    order_tags.append(n)

tags_argsort = np.argsort(order_tags)

tags = tags[tags_argsort]
obs_diurnal_mag = obs_diurnal_mag[tags_argsort]
obs_diurnal_phase = obs_diurnal_phase[tags_argsort]
model_diurnal_mag = model_diurnal_mag[tags_argsort]
model_diurnal_phase = model_diurnal_phase[tags_argsort]

tag_colors = [loc_colors[i] for i in tags]
tag_locs = [all_locs[i] for i in tags]

#set up plot
fig =plt.figure(figsize=(17,12))
fig.patch.set_facecolor('white')

pi2 = np.pi*2.
ratio = pi2/24.
obs_diurnal_phase_rads = [i*ratio for i in obs_diurnal_phase]
model_diurnal_phase_rads = [i*ratio for i in model_diurnal_phase]

ax1 = fig.add_subplot(121, polar=True)
ax1.set_theta_zero_location('N')
ax1.set_theta_direction(-1) 
for i in range(len(obs_diurnal_mag)):
    ax1.scatter(obs_diurnal_phase_rads[i],obs_diurnal_mag[i],c=tag_colors[i],s=60,label = tag_locs[i])
ax1.set_xticks(ratio*np.arange(24))
ax1.set_xticklabels(['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'])
ax1.set_rmax(np.max(model_diurnal_mag)+2)
ax1.set_title("Observations", fontsize=25,y=1.10,weight='bold')

ax1.tick_params(axis='both', which='major', labelsize=22)

ax2 = fig.add_subplot(122, polar=True)
ax2.set_theta_zero_location('N')
ax2.set_theta_direction(-1) 
for i in range(len(model_diurnal_mag)):
    ax2.scatter(model_diurnal_phase_rads[i],model_diurnal_mag[i],c=tag_colors[i],s=60,label=tag_locs[i])
ax2.set_xticks(ratio*np.arange(24))
ax2.set_xticklabels(['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'])
ax2.set_rmax(np.max(model_diurnal_mag)+2)
ax2.set_title("Model", fontsize=25,y=1.10,weight='bold')

ax2.tick_params(axis='both', which='major', labelsize=22)

ax1.text(0.30, 0.90, 'Amplitude (ppb)', transform=ax1.transAxes,fontsize=20, va='top')
ax1.text(0.85, 0.97, 'Phase (Hour)', transform=ax1.transAxes,fontsize=20, va='top')

plt.tight_layout(pad = 3.08)

all_handles, all_labels = ax1.get_legend_handles_labels()
hl = sorted(zip(all_handles, all_labels), key=operator.itemgetter(1))
all_handles, all_labels = zip(*hl)
by_label = OrderedDict(zip(all_labels, all_handles))
plt.legend(by_label.values(),by_label.keys(), loc='lower left',prop={'size':22},fancybox=True,ncol=4,markerscale=2,bbox_to_anchor=(-1.0,-0.25))

plt.show()
