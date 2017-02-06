#import matplotlib
#matplotlib.use('Agg')
import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
#from mpl_toolkits.basemap import Basemap, shiftgrid
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


def interactive(event):
    modules.clicker_interactive_xy_obsmodel_single(event,species,lat_e,lon_e,obs_lats,obs_lons,obs_date,obs_time,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_chosen,model_chosen,obs_d_waveform,obs_s_waveform,obs_all_waveform,model_d_waveform,model_s_waveform,model_all_waveform,fig,ax)

#-----------------------------------------------------
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
obs_d_waveform = []
obs_s_waveform = []
obs_all_waveform = []
model_d_waveform = []
model_s_waveform = []
model_all_waveform = []
model_daily_ff = []
model_seasonal_ff = []
obs_daily_ff = []
obs_seasonal_ff = []

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
    model_d_waveform.append(site_group_mod.variables['daily_waveform'][:])
    model_s_waveform.append(site_group_mod.variables['seasonal_waveform'][:])
    model_all_waveform.append(site_group_mod.variables['all_waveform'][:])
    model_daily_ff = np.append(model_daily_ff,site_group_mod.daily_ff) 
    model_seasonal_ff = np.append(model_seasonal_ff,site_group_mod.seasonal_ff)
    
    obs_daily_mag = np.append(obs_daily_mag,site_group_obs.daily_amplitude)
    obs_daily_phase = np.append(obs_daily_phase,site_group_obs.daily_phase)
    obs_daily_phase_min = np.append(obs_daily_phase_min,site_group_obs.daily_phase_min)
    obs_ave = np.append(obs_ave,site_group_obs.average)
    obs_seasonal_mag = np.append(obs_seasonal_mag,site_group_obs.seasonal_amplitude)
    obs_seasonal_phase = np.append(obs_seasonal_phase,site_group_obs.seasonal_phase)
    obs_seasonal_phase_min = np.append(obs_seasonal_phase_min,site_group_obs.seasonal_phase_min)
    obs_d_waveform.append(site_group_obs.variables['daily_waveform'][:])
    obs_s_waveform.append(site_group_obs.variables['seasonal_waveform'][:])
    obs_all_waveform.append(site_group_obs.variables['all_waveform'][:])
    obs_daily_ff = np.append(obs_daily_ff,site_group_obs.daily_ff) 
    obs_seasonal_ff = np.append(obs_seasonal_ff,site_group_obs.seasonal_ff)

month_lengths = [0,31,59.25,90.25,120.25,151.25,181.25,212.25,243.25,273.25,304.25,334.25,365.25]
#month_lengths = [273.25,304.25,334.25,365.25,396.25,424.5,455.5,485.5,516.5,546.5,577.5,608.5]
year_fact = 12./365.25
month_array  = [i*year_fact for i in month_lengths]
month_strings = ["JAN","","MAR","","MAY","","JUL","","SEP","","NOV",""]
#month_strings = ["OCT","","DEC","","FEB","","APR","","JUN","","AUG",""]
axis_nums = month_array
		
#--------------------------------------------------------------------------
#rearrange data so that NA plotted 1, EU plotted 2, AS plotted 3, O plotted 4, OC plotted 5, AF plotted 6, SA plotted 7, ANT plotted 8, ARC plotted 9
loc_dict = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic Sites','NA':'North America','SA':'South America'}
order_dict = {'ANT':8,'ARC':9,'AF':6,'AS':3,'EU':2,'OC':5,'O':4,'NA':1,'SA':7}

order_tags = []
for tag in tags:
    n = order_dict[tag]
    order_tags.append(n)

tags_argsort = np.argsort(order_tags)

tags = tags[tags_argsort]
obs_refs = obs_refs[tags_argsort]
obs_lats = obs_lats[tags_argsort]
obs_lons = obs_lons[tags_argsort]
obs_alt = obs_alt[tags_argsort]
obs_daily_mag = obs_daily_mag[tags_argsort]
model_daily_mag = model_daily_mag[tags_argsort]
obs_daily_phase = obs_daily_phase[tags_argsort]
model_daily_phase = model_daily_phase[tags_argsort]
obs_ave = obs_ave[tags_argsort]
model_ave = model_ave[tags_argsort]
model_seasonal_mag = model_seasonal_mag[tags_argsort]
model_seasonal_phase = model_seasonal_phase[tags_argsort]
obs_seasonal_mag = obs_seasonal_mag[tags_argsort]
obs_seasonal_phase = obs_seasonal_phase[tags_argsort]
obs_d_waveform = np.array(obs_d_waveform)[tags_argsort]
obs_s_waveform = np.array(obs_s_waveform)[tags_argsort]
obs_all_waveform = np.array(obs_all_waveform)[tags_argsort]
model_d_waveform = np.array(model_d_waveform)[tags_argsort]
model_s_waveform = np.array(model_s_waveform)[tags_argsort]
model_all_waveform = np.array(model_all_waveform)[tags_argsort]
obs_daily_ff = obs_daily_ff[tags_argsort]
obs_seasonal_ff = obs_seasonal_ff[tags_argsort]
model_daily_ff = model_daily_ff[tags_argsort]
model_seasonal_ff = model_seasonal_ff[tags_argsort]


#---------------------------------------------------------------------------

#setup shape dict for plotting
loc_colors = {'ANT':'magenta','ARC':'purple','AF':'orange','AS':'blue','EU':'red','OC':'black','O':'yellow','NA':'green','SA':'pink'}
loc_shapes = {'ANT':'o','ARC':'+','AF':'*','AS':'d','EU':'s','OC':'>','O':'v','NA':'<','SA':'^'}

typ = raw_input('amp, ph_min, ph_max, ave or ff?\n')
if typ != 'ave':
    period = raw_input('\nd or s?\n')
else:
    period = ''
#--------------------------------------------------------------------
#set up plot
fig, (ax) = plt.subplots(1,figsize=(15,13))
fig.patch.set_facecolor('white')

marker_size = 15
font_size = 25
pad_size = 5

#get colors and labels for site points
plot_colors = []

for t in range(len(obs_refs)):
    plot_colors.append(loc_colors[tags[t]])
    
#----------------------------------------------------------------------
#Daily Amplitude Plot

if typ == 'amp':
    if period == 'd':
        obs_chosen = obs_daily_mag
        model_chosen = model_daily_mag
        
        for t in range(len(obs_refs)):
            line, = ax.plot(obs_daily_mag[t], model_daily_mag[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

        ax.plot(obs_daily_mag, model_daily_mag, color = 'black', marker = 'o', markersize = 1, linestyle= 'None',zorder=1, picker = 5)

        obs_max = np.max(obs_daily_mag)+2
        model_max = np.max(model_daily_mag)+2
        if np.min(obs_daily_mag) < np.min(model_daily_mag):
            min_val = np.min(obs_daily_mag)
        else:
            min_val = np.min(model_daily_mag)
        if obs_max > model_max:
            ax.set_xlim(min_val,obs_max)
            ax.set_ylim(min_val,obs_max)
        else:
            ax.set_xlim(min_val,model_max)
            ax.set_ylim(min_val,model_max)
        x = np.arange(0,1000,0.005)
        y = np.arange(0,1000,0.005)

        handles, labels = ax.get_legend_handles_labels()
        hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
        handles, labels = zip(*hl)
        by_label = OrderedDict(zip(labels, handles))
        ax.plot(x,y,'k--',alpha = 0.5)
        ax.set_ylabel('Model (ppb)',fontsize = font_size)
        ax.set_xlabel('Observations (ppb)',fontsize = font_size)
			
#----------------------------------------------------------------------
#Seasonal Amplitude Plot
 
if typ == 'amp':
    if period == 's': 
        obs_chosen = obs_seasonal_mag
        model_chosen = model_seasonal_mag
        
        for t in range(len(obs_refs)):
            ax.plot(obs_seasonal_mag[t], model_seasonal_mag[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

        ax.plot(obs_seasonal_mag, model_seasonal_mag, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

        obs_max = np.max(obs_seasonal_mag)+1
        model_max = np.max(model_seasonal_mag)+1
        if np.min(obs_seasonal_mag) < np.min(model_seasonal_mag):
            min_val = np.min(obs_seasonal_mag)
        else:
            min_val = np.min(model_seasonal_mag)
        if obs_max > model_max:
            ax.set_xlim(min_val,obs_max)
            ax.set_ylim(min_val,obs_max)
        else:
            ax.set_xlim(min_val,model_max)
            ax.set_ylim(min_val,model_max)
        x = np.arange(0,1000,1)
        y = np.arange(0,1000,1)

        ax.plot(x,y,'k--',alpha = 0.5)
        ax.set_ylabel('Model (ppb)',fontsize = font_size)
        ax.set_xlabel('Observations (ppb)',fontsize = font_size)

#----------------------------------------------------------------------
#Daily max Phase Plot

if typ == 'ph_max':
    if period == 'd': 
        obs_chosen = obs_daily_phase
        model_chosen = model_daily_phase
        
        for t in range(len(obs_refs)):
            ax.plot(obs_daily_phase[t], model_daily_phase[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

        ax.plot(obs_daily_phase, model_daily_phase, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

        ax.set_xlim(0,24)
        ax.set_ylim(0,24)
        x = np.arange(0,50,1)
        y = np.arange(0,50,1)

        ax.set_xticks(range(0,25,2))
        ax.set_yticks(range(0,25,2))

        ax.set_ylabel('Model (Hours)',fontsize = font_size)
        ax.set_xlabel('Observations (Hours)',fontsize = font_size)
        ax.plot(x,y,'k--',alpha = 0.5)
        for tick in ax.get_xaxis().get_major_ticks():
            tick.set_pad(pad_size)
        for tick in ax.get_yaxis().get_major_ticks():
            tick.set_pad(pad_size)
            
#----------------------------------------------------------------------
#Daily min Phase Plot

if typ == 'ph_min':
    if period == 'd': 
        obs_chosen = obs_daily_phase_min
        model_chosen = model_daily_phase_min
        
        for t in range(len(obs_refs)):
            ax.plot(obs_daily_phase_min[t], model_daily_phase_min[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

        ax.plot(obs_daily_phase_min, model_daily_phase_min, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

        ax.set_xlim(0,24)
        ax.set_ylim(0,24)
        x = np.arange(0,50,1)
        y = np.arange(0,50,1)

        ax.set_xticks(range(0,25,2))
        ax.set_yticks(range(0,25,2))

        ax.set_ylabel('Model (Hours)',fontsize = font_size)
        ax.set_xlabel('Observations (Hours)',fontsize = font_size)
        ax.plot(x,y,'k--',alpha = 0.5)
        for tick in ax.get_xaxis().get_major_ticks():
            tick.set_pad(pad_size)
        for tick in ax.get_yaxis().get_major_ticks():
            tick.set_pad(pad_size)
            
#----------------------------------------------------------------------
#Seasonal max Phase Plot
if typ == 'ph_max':
    if period == 's': 
        obs_chosen = obs_seasonal_phase
        model_chosen = model_seasonal_phase
        #change phase to run from sep to sep, instead of jan to jan
        #obs_seasonal_phase_corr, model_seasonal_phase_corr = modules.annual_phase_shift(obs_seasonal_phase,model_seasonal_phase)

        for t in range(len(obs_refs)):
            ax.plot(obs_seasonal_phase[t], model_seasonal_phase[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

        ax.plot(obs_seasonal_phase, model_seasonal_phase, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)
        ax.set_xlim(0,11)
        ax.set_ylim(0,11)
        x = np.arange(0,50,1)
        y = np.arange(0,50,1)

        ax.set_xticks(month_array)
        ax.set_yticks(month_array)
        ax.set_xticklabels(month_strings,fontsize=9.2)
        ax.set_yticklabels(month_strings,fontsize=9.2)
        ax.set_ylabel('Model (Months)',fontsize = font_size)
        ax.set_xlabel('Observations (Months)',fontsize = font_size)
        ax.plot(x,y,'k--',alpha = 0.5)
        for tick in ax.get_xaxis().get_major_ticks():
            tick.set_pad(pad_size)
        for tick in ax.get_yaxis().get_major_ticks():
            tick.set_pad(pad_size)
            
#----------------------------------------------------------------------
#Seasonal min Phase Plot
if typ == 'ph_min':
    if period == 's': 
        obs_chosen = obs_seasonal_phase_min
        model_chosen = model_seasonal_phase_min
        #change phase to run from sep to sep, instead of jan to jan
        #obs_seasonal_phase_corr, model_seasonal_phase_corr = modules.annual_phase_shift(obs_seasonal_phase,model_seasonal_phase)

        for t in range(len(obs_refs)):
            ax.plot(obs_seasonal_phase_min[t], model_seasonal_phase_min[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

        ax.plot(obs_seasonal_phase_min, model_seasonal_phase_min, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)
        ax.set_xlim(0,11)
        ax.set_ylim(0,11)
        x = np.arange(0,50,1)
        y = np.arange(0,50,1)

        ax.set_xticks(month_array)
        ax.set_yticks(month_array)
        ax.set_xticklabels(month_strings,fontsize=9.2)
        ax.set_yticklabels(month_strings,fontsize=9.2)
        ax.set_ylabel('Model (Months)',fontsize = font_size)
        ax.set_xlabel('Observations (Months)',fontsize = font_size)
        ax.plot(x,y,'k--',alpha = 0.5)
        for tick in ax.get_xaxis().get_major_ticks():
            tick.set_pad(pad_size)
        for tick in ax.get_yaxis().get_major_ticks():
            tick.set_pad(pad_size)

#----------------------------------------------------------------------
#Average Plot   

if typ == 'ave':
    obs_chosen = obs_ave
    model_chosen = model_ave

    for t in range(len(obs_refs)):
        ax.plot(obs_ave[t], model_ave[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

    ax.plot(obs_ave, model_ave, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

    obs_max = np.max(obs_ave)+1
    model_max = np.max(model_ave)+1
    if np.min(obs_ave) < np.min(model_ave):
        min_val = np.min(obs_ave)-1
    else:
        min_val = np.min(model_ave)-1
    if obs_max > model_max:
        ax.set_xlim(min_val,obs_max)
        ax.set_ylim(min_val,obs_max)
    else:
        ax.set_xlim(min_val,model_max)
        ax.set_ylim(min_val,model_max)
    x = np.arange(0,1000,1)
    y = np.arange(0,1000,1)

    ax.plot(x,y,'k--')

    ax.set_ylabel('Model (ppb)',fontsize = font_size)
    ax.set_xlabel('Observations (ppb)',fontsize = font_size)

#----------------------------------------------------------------------
#Daily Form Factor Plot   
if typ == 'ff':
    if period == 'd':
        obs_chosen = obs_daily_ff
        model_chosen = model_daily_ff

        for t in range(len(obs_refs)):
            ax.plot(obs_daily_ff[t], model_daily_ff[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

        ax.plot(obs_daily_ff, model_daily_ff, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

        obs_max = np.max(obs_daily_ff)+0.1
        model_max = np.max(model_daily_ff)+0.1
        if np.min(obs_daily_ff) < np.min(model_daily_ff):
            min_val = np.min(obs_daily_ff)-0.1
        else:
            min_val = np.min(model_daily_ff)-0.1
        if obs_max > model_max:
            ax.set_xlim(min_val,obs_max)
            ax.set_ylim(min_val,obs_max)
        else:
            ax.set_xlim(min_val,model_max)
            ax.set_ylim(min_val,model_max)
        x = np.arange(0,1000,1)
        y = np.arange(0,1000,1)

        ax.plot(x,y,'k--')

        ax.set_ylabel('Model Form Factor',fontsize = font_size)
        ax.set_xlabel('Observations Form Factor',fontsize = font_size)
        
#----------------------------------------------------------------------
#Seasonal Form Factor Plot   
if typ == 'ff':
    if period == 's':
        obs_chosen = obs_seasonal_ff
        model_chosen = model_seasonal_ff

        for t in range(len(obs_refs)):
            ax.plot(obs_seasonal_ff[t], model_seasonal_ff[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

        ax.plot(obs_seasonal_ff, model_seasonal_ff, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

        obs_max = np.max(obs_seasonal_ff)+0.1
        model_max = np.max(model_seasonal_ff)+0.1
        if np.min(obs_seasonal_ff) < np.min(model_seasonal_ff):
            min_val = np.min(obs_seasonal_ff)-0.1
        else:
            min_val = np.min(model_seasonal_ff)-0.1
        if obs_max > model_max:
            ax.set_xlim(min_val,obs_max)
            ax.set_ylim(min_val,obs_max)
        else:
            ax.set_xlim(min_val,model_max)
            ax.set_ylim(min_val,model_max)
        x = np.arange(0,1000,1)
        y = np.arange(0,1000,1)

        ax.plot(x,y,'k--')

        ax.set_ylabel('Model Form Factor',fontsize = font_size)
        ax.set_xlabel('Observations Form Factor',fontsize = font_size)

#--------------------------------------------------
#make axis labels tight to plots
ax.yaxis.labelpad = 12
ax.xaxis.labelpad = 12

ax.tick_params(axis='both', which='major', labelsize=24)
ax.tick_params(axis='both', which='minor', labelsize=24)

ax.grid(True)

x = np.arange(0,1000,1)
y = np.arange(0,1000,1)

ax.plot(x,y,'k--',zorder=20,linewidth=3)

#plot legend
handles, labels = ax.get_legend_handles_labels()
hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
handles, labels = zip(*hl)
by_label = OrderedDict(zip(labels, handles))
leg = ax.legend(by_label.values(), by_label.keys(), loc = 'upper center', bbox_to_anchor=(0.95,0.4),fancybox=True,ncol=1,prop={'size':22})
leg.set_zorder(100) 

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500-800")
fig.canvas.mpl_connect('pick_event', interactive)

plt.savefig('seasonal_amp_xy.png')

plt.show()
