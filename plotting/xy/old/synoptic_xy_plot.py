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

def form2(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.2f' % x	

def form5(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.7f' % x

xformatter = FuncFormatter(form2)
yformatter = FuncFormatter(form5)

#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres = modules.get_obs_info(present_dir)
model_fname,species,start_year,end_year = modules.get_model_info(present_dir)
    
root_grp = Dataset(model_fname)
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
obs_root_grp = Dataset(obs_fname)
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
#load in spectra lsp data
root_grp_obs = Dataset('../obs_%s_%s/obs_spectra.nc'%(vres,timeres))
root_grp_mod = Dataset('model_spectra.nc')

obs_periods = []
model_periods = []
obs_gradient = []
model_gradient = []
obs_bp = []
model_bp = []
obs_med = []
model_med = []
obs_ave = []
model_ave = []
obs_ys = []
model_ys = []
obs_ye = []
model_ye = []
obs_syn_period = []
obs_mw_period = []
model_syn_period = []
model_mw_period = []
obs_syn_amp = []
obs_mw_amp = []
model_syn_amp = []
model_mw_amp = []

for ref in obs_refs:
    site_group_obs = root_grp_obs.groups[ref]
    site_group_mod = root_grp_mod.groups[ref]   
    obs_syn_period.append(site_group_obs.variables['synoptic_period'][:])
    model_syn_period.append(site_group_mod.variables['synoptic_period'][:])
    obs_mw_period.append(site_group_obs.variables['macroweather_period'][:])
    model_mw_period.append(site_group_mod.variables['macroweather_period'][:])
    obs_syn_amp.append(site_group_obs.variables['synoptic_amplitude'][:])
    model_syn_amp.append(site_group_mod.variables['synoptic_amplitude'][:])
    obs_mw_amp.append(site_group_obs.variables['macroweather_amplitude'][:])
    model_mw_amp.append(site_group_mod.variables['macroweather_amplitude'][:])
  
    model_gradient = np.append(model_gradient,site_group_mod.synoptic_gradient)
    model_bp = np.append(model_bp,site_group_mod.synoptic_macroweather_transition)
    model_med = np.append(model_med,site_group_mod.synoptic_median)    
    model_ave =  np.append(model_ave,site_group_mod.synoptic_average)
    model_ys = np.append(model_ys,site_group_mod.synoptic_start_amplitude)
    model_ye = np.append(model_ye,site_group_mod.synoptic_end_amplitude)

    obs_gradient = np.append(obs_gradient,site_group_obs.synoptic_gradient)
    obs_bp = np.append(obs_bp,site_group_obs.synoptic_macroweather_transition)
    obs_med = np.append(obs_med,site_group_obs.synoptic_median)  
    obs_ave = np.append(obs_ave,site_group_obs.synoptic_average)
    obs_ys = np.append(obs_ys,site_group_obs.synoptic_start_amplitude)
    obs_ye = np.append(obs_ye,site_group_obs.synoptic_end_amplitude)


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
obs_gradient = obs_gradient[tags_argsort]
model_gradient = model_gradient[tags_argsort]
obs_bp = obs_bp[tags_argsort]
model_bp = model_bp[tags_argsort]
obs_med = obs_med[tags_argsort]
model_med = model_med[tags_argsort]
obs_ave = obs_ave[tags_argsort]
model_ave = model_ave[tags_argsort]
obs_ys = obs_ys[tags_argsort]
model_ys = model_ys[tags_argsort]
obs_ye = obs_ye[tags_argsort]
model_ye = model_ye[tags_argsort]


#setup shape dict for plotting
loc_colors = {'ANT':'magenta','ARC':'purple','AF':'orange','AS':'blue','EU':'red','OC':'black','O':'yellow','NA':'green','SA':'pink'}
loc_shapes = {'ANT':'o','ARC':'+','AF':'*','AS':'d','EU':'s','OC':'>','O':'v','NA':'<','SA':'^'}

#--------------------------------------------------------------------
#set up plot
fig, ((ax, ax2), (ax3, ax4),(ax5, ax6)) = plt.subplots(3, 2,figsize=(13,16))
fig.patch.set_facecolor('white')

marker_size = 6
font_size = 10
pad_size = 5

#get colors and labels for site points
plot_colors = []

for t in range(len(obs_refs)):
    plot_colors.append(loc_colors[tags[t]])
    
#----------------------------------------------------------------------
#Gradient

for t in range(len(obs_refs)):
    line, = ax.plot(obs_gradient[t], model_gradient[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

ax.plot(obs_gradient, model_gradient, color = 'black', marker = 'o', markersize = 1, linestyle= 'None',zorder=1, picker = 5)

obs_max = np.max(obs_gradient)+0.01
model_max = np.max(model_gradient)+0.01
if np.min(obs_gradient) < np.min(model_gradient):
    min_val = np.min(obs_gradient)-0.01
else:
    min_val = np.min(model_gradient)-0.01
if obs_max > model_max:
    ax.set_xlim(min_val,obs_max)
    ax.set_ylim(min_val,obs_max)
else:
    ax.set_xlim(min_val,model_max)
    ax.set_ylim(min_val,model_max)
x = np.arange(-1000,1000,0.005)
y = np.arange(-1000,1000,0.005)

handles, labels = ax.get_legend_handles_labels()
hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
handles, labels = zip(*hl)
by_label = OrderedDict(zip(labels, handles))
ax.plot(x,y,'k--',alpha = 0.5)

ax.tick_params(axis='both', which='major', labelsize=20)
ax.tick_params(axis='both', which='minor', labelsize=20)

for tick in ax.get_xaxis().get_major_ticks():
    tick.set_pad(7.5)
for tick in ax.get_yaxis().get_major_ticks():
    tick.set_pad(7.5)

ax.text(0.65, 0.97, 'Gradient', transform=ax.transAxes,fontsize=19, fontweight='bold', va='top')
		
#----------------------------------------------------------------------
#Median
obs_med = np.log(obs_med)
model_med = np.log(model_med)

for t in range(len(obs_refs)):
    ax2.plot(obs_med[t], model_med[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

ax2.plot(obs_med, model_med, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

obs_max = np.max(obs_med)+0.01
model_max = np.max(model_med)+0.01
if np.min(obs_med) < np.min(model_med):
    min_val = np.min(obs_med)-0.01
else:
    min_val = np.min(model_med)-0.01
if obs_max > model_max:
    ax2.set_xlim(min_val,obs_max)
    ax2.set_ylim(min_val,obs_max)
else:
    ax2.set_xlim(min_val,model_max)
    ax2.set_ylim(min_val,model_max)

ax2.plot(x,y,'k--',alpha = 0.5)

ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='minor', labelsize=20)

for tick in ax2.get_xaxis().get_major_ticks():
    tick.set_pad(7.5)
for tick in ax2.get_yaxis().get_major_ticks():
    tick.set_pad(7.5)
    
ax2.text(0.23, 0.97, 'Log Average Amplitude', transform=ax2.transAxes,fontsize=19, fontweight='bold', va='top')

#----------------------------------------------------------------------
#Y 
obs_ys = np.log(obs_ys)
model_ys = np.log(model_ys)

for t in range(len(obs_refs)):
    ax3.plot(obs_ys[t], model_ys[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

ax3.plot(obs_ys, model_ys, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

obs_max = np.max(obs_ys)+0.01
model_max = np.max(model_ys)+0.01
if np.min(obs_ys) < np.min(model_ys):
    min_val = np.min(obs_ys)-0.01
else:
    min_val = np.min(model_ys)-0.01
if obs_max > model_max:
    ax3.set_xlim(min_val,obs_max)
    ax3.set_ylim(min_val,obs_max)
else:
    ax3.set_xlim(min_val,model_max)
    ax3.set_ylim(min_val,model_max)

ax3.plot(x,y,'k--',alpha = 0.5)

ax3.tick_params(axis='both', which='major', labelsize=20)
ax3.tick_params(axis='both', which='minor', labelsize=20)

for tick in ax3.get_xaxis().get_major_ticks():
    tick.set_pad(7.5)
for tick in ax3.get_yaxis().get_major_ticks():
    tick.set_pad(7.5)
    
ax3.text(0.25, 0.97, 'Log Lowest Amplitude', transform=ax3.transAxes,fontsize=19, fontweight='bold', va='top')

#----------------------------------------------------------------------                                                                                                           
#Y 
obs_ye = np.log(obs_ye)
model_ye = np.log(model_ye)

for t in range(len(obs_refs)):
    ax4.plot(obs_ye[t], model_ye[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)

ax4.plot(obs_ye, model_ye, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

obs_max = np.max(obs_ye)+0.01
model_max = np.max(model_ye)+0.01
if np.min(obs_ye) < np.min(model_ye):
    min_val = np.min(obs_ye)-0.01
else:
    min_val = np.min(model_ye)-0.01
if obs_max > model_max:
    ax4.set_xlim(min_val,obs_max)
    ax4.set_ylim(min_val,obs_max)
else:
    ax4.set_xlim(min_val,model_max)
    ax4.set_ylim(min_val,model_max)

ax4.plot(x,y,'k--',alpha = 0.5)

ax4.tick_params(axis='both', which='major', labelsize=20)
ax4.tick_params(axis='both', which='minor', labelsize=20)

for tick in ax4.get_xaxis().get_major_ticks():
    tick.set_pad(7.5)
for tick in ax4.get_yaxis().get_major_ticks():
    tick.set_pad(7.5)
    
ax4.text(0.22, 0.97, 'Log Highest Amplitude', transform=ax4.transAxes,fontsize=19, fontweight='bold', va='top')
#--------------------------------------------------------------------------------
#remove ax 5 and ax6 for table

#ax4.set_frame_on(False)
ax5.set_frame_on(False)
ax6.set_frame_on(False)
#ax4.axes.get_yaxis().set_visible(False)
ax5.axes.get_yaxis().set_visible(False)
ax6.axes.get_yaxis().set_visible(False)
#ax4.axes.get_xaxis().set_visible(False) 
ax5.axes.get_xaxis().set_visible(False)
ax6.axes.get_xaxis().set_visible(False)

#----------------------------------------------------------------------------------
#plot legend
handles, labels = ax.get_legend_handles_labels()
hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
handles, labels = zip(*hl)
by_label = OrderedDict(zip(labels, handles))
leg = ax6.legend(by_label.values(),by_label.keys(), loc='lower right',prop={'size':22},fancybox=True,ncol=1)

#-------------------------------------------------------------------------------
#do orthagonal regression and output stats

def f(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x + B[1]


def orth_regression(obs,model):
    linear = Model(f)
    mydata = RealData(obs, model)
    myodr = ODR(mydata, linear, beta0=[1., 0.])
    myoutput = myodr.run()
    params = myoutput.beta
    gradient = params[0]
    y_intercept = params[1]
    res_var = myoutput.res_var
    return np.around(gradient,2), np.around(y_intercept,2), np.around(res_var,2)

#gradient_grad,gradient_y,gradient_res = orth_regression(obs_gradient,model_gradient)
#ave_grad,ave_y,ave_res = orth_regression(obs_ave,model_ave)
#y_first_grad,y_first_y,y_first_res = orth_regression(y_obs_first,y_model_first)
#y_last_grad,y_last_y,y_last_res = orth_regression(y_obs_last,y_model_last) 

#gradient_res_sum,gradient_res_ave = modules.orthogonal_1_1_res_var(obs_gradient,model_gradient)
#ave_res_sum,ave_res_ave = modules.orthogonal_1_1_res_var(obs_ave,model_ave)
#y_first_res_sum,y_first_res_ave = modules.orthogonal_1_1_res_var(y_obs_first,y_model_first)
#y_last_res_sum,y_last_res_ave = modules.orthogonal_1_1_res_var(y_obs_last,y_model_last)

#----------------------------------------------------------------------
# Add a table at the right of the axes
#col_labels=['Gradient','Y-Intercept','Ave Residual Variance']
#row_labels=['Gradient','Log Average Amplitude','Log Lowest Amplitude','Log Highest Amplitude']
#table = ax5.table(cellText=[[gradient_grad,gradient_y,gradient_res_ave],[ave_grad,ave_y,ave_res_ave],[y_first_grad,y_first_y,y_first_res_ave],[y_last_grad,y_last_y,y_last_res_ave]],colWidths=[0.27,0.32,0.53],rowLabels=row_labels,colLabels=col_labels,loc='bottom',bbox=[0.8, 0, 1.286,1])

plt.tight_layout(pad = 3.08)

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500-800")

fig.canvas.mpl_connect('pick_event', onpick)


plt.show()
