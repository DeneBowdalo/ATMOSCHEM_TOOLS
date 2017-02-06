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

print 'yes'

def form2(x, pos):
    """ This function returns a string with 3 decimal places, given the input x"""
    return '%.2f' % x

def form6(x, pos):
    """ This function returns a string with 3 decimal places, given the input x"""
    return '%.6f' % x

#set up plot
fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)


period = 'daily'

plot_type = 'phase'

third_var = 'yes'

if  period == 'daily':
	model_phase_file = 'model_phases/GAW_model_daily_phases.npy'
	obs_phase_file = 'obs_phases/GAW_obs_daily_phases.npy' 
	model_mag_file = 'model_magnitudes/GAW_model_daily_magnitudes.npy'
	obs_mag_file = 'obs_magnitudes/GAW_obs_daily_magnitudes.npy' 
	title = 'Daily'
	con_num = 24
	label = 'Time (Hours)'
	phase_min = 0
	phase_max = 24
	mnames = range(phase_min,phase_max)
	axis_nums = range(phase_min,phase_max)


if period == 'annual': 
	model_phase_file = 'model_phases/GAW_model_annual_phases.npy'
	obs_phase_file = 'obs_phases/GAW_obs_annual_phases.npy'
	model_mag_file = 'model_magnitudes/GAW_model_annual_magnitudes.npy'
	obs_mag_file = 'obs_magnitudes/GAW_obs_annual_magnitudes.npy'
	title = 'Annual'
	con_num = 12
	label = ' Time (Months)'
	phase_min = 9
	phase_max = 21
	mnames=["SEP","OCT","NOV","DEC","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG"]
	axis_nums = range(phase_min,phase_max)
	
	
# load in model values
model_mag = np.load(model_mag_file)
model_phase = np.load(model_phase_file)

#read in obs values
obs_mag = np.load(obs_mag_file)
obs_phase = np.load(obs_phase_file)

#lat lon centres & edges for 2x2.5 grid 
lat_c = np.arange(-88.,89.,2)
lat_c = np.insert(lat_c,0,-89.5)
lat_c = np.append(lat_c,89.5)

lon_c = np.arange(-180,178,2.5)

lat_e = np.arange(-89.,90,2)
lat_e = np.insert(lat_e,0,-90.)
lat_e = np.append(lat_e,90.)

lon_e = np.arange(-181.25,179,2.5)

#read in obs lats & lons
obs_refs,obs_locs,obs_lats,obs_lons,obs_alt,obs_number = np.loadtxt('GAW_site_indices',dtype='S20,S20,f5,f5,f5,i5',unpack=True)

#get obs lat_lon grid central points
obs_lats_centre, obs_lons_centre, model_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_lats,obs_lons) 

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)
loc_dict = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic Sites','NA':'North America','SA':'South America'}

# remove sites over 1000m
valid_indices = []
high_sites = []
for num in range(len(obs_alt)):	
	if obs_alt[num] > 1000:
		high_sites.append(num)
	else:
		valid_indices.append(num)

tags[high_sites] = 'H'

#setup shape dict for plotting
loc_colors = {'ANT':'red','ARC':'pink','AF':'black','AS':'orange','EU':'purple','OC':'blue','O':'yellow','NA':'magenta','SA':'green'}
loc_shapes = {'ANT':'o','ARC':'+','AF':'*','AS':'d','EU':'s','OC':'>','O':'v','NA':'<','SA':'^'}

#limit obs and model data to valid data
tags_test = tags != 'H'
tags = tags[tags_test]
obs_refs = obs_refs[valid_indices]
obs_alt = obs_alt[valid_indices]
obs_mag = obs_mag[valid_indices]
model_mag = model_mag[valid_indices]
obs_phase = obs_phase[valid_indices]
model_phase = model_phase[valid_indices]

# create way of formatting number of decimal places
xformatter = FuncFormatter(form2)
yformatter = FuncFormatter(form2)

#plot data

#when plotting magnitude x-y data, colour points by phase(average of obs & model phase), and give shape on location
if plot_type == 'mag': 
    if period == 'daily':
        for t in range(len(obs_refs)):
            if third_var == 'yes':
                plt.scatter(model_mag[t], obs_mag[t], c = obs_phase[t],  cmap=plt.cm.hsv, marker = loc_shapes[tags[t]], s = 100, vmin = phase_min, vmax= phase_max, label = loc_dict[tags[t]])	
            else:
                plt.plot(model_mag[t], obs_mag[t], color = loc_colors[tags[t]], marker = 'o', markersize = 15, label = loc_dict[tags[t]])
		#plot site labels
        for i, txt in enumerate(obs_refs):
            plt.annotate(txt, (model_mag[i]+(model_mag[i]/15),obs_mag[i] -(obs_mag[i]/25) ), fontweight='bold')
		
        ax.set_xscale('log')
        ax.set_yscale('log')
        if third_var == 'yes':
            cb = plt.colorbar(ticks=range(phase_min,phase_max+1))	
            cb_ax = cb.ax
            cb_ax.text(-1,1.025,'Phase (Hours)', fontsize = 18)
        obs_max = np.max(obs_mag)+2
        model_max = np.max(model_mag)+2
        if np.min(obs_mag) < np.min(model_mag):
            min_val = np.min(obs_mag)-0.001
        else:
            min_val = np.min(model_mag)-0.001
        if obs_max > model_max:
            plt.xlim(min_val,obs_max)
            plt.ylim(min_val,obs_max)
        else:
            plt.xlim(min_val,model_max)
            plt.ylim(min_val,model_max)
        x = np.arange(0,1000,0.005)
        y = np.arange(0,1000,0.005)
		
        handles, labels = plt.gca().get_legend_handles_labels()
        hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
        handles, labels = zip(*hl)
        by_label = OrderedDict(zip(labels, handles))
        leg = ax.legend(by_label.values(), by_label.keys(), loc = 0, prop={'size':20})
		#lines=leg.get_children()
		#print lines
		#lines[0].set_color('red')
		#print lines
        plt.plot(x,y,'k--',alpha = 0.5)
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())	
        plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
        plt.xlabel('GEOS-Chem 4x5 Magnitude (ppbv)',fontsize = 20)
        plt.ylabel('Obs. Magnitude (ppbv)',fontsize = 20)
        plt.title(r'%s Surface $O_3$ Magnitude, X-Y Plot of Model v Obs.'%title, fontsize= 22)
			
    elif period == 'annual':
        for t in range(len(obs_refs)):
            if third_var == 'yes':
                plt.scatter(model_mag[t], obs_mag[t], c = obs_phase[t],  cmap=plt.cm.hsv, marker = loc_shapes[tags[t]], s = 100, vmin = 0, vmax= 12, label = loc_dict[tags[t]])
            else:
                plt.plot(model_mag[t], obs_mag[t], color = loc_colors[tags[t]], marker = 'o', markersize = 15, label = loc_dict[tags[t]])
		#plot site labels
        for i, txt in enumerate(obs_alt):
            x_pos = (model_mag[i]+0.15)
            y_pos = (obs_mag[i]-0.075)
            plt.annotate(txt, (x_pos,y_pos) , fontweight='bold')
		
        if third_var == 'yes':
            cb = plt.colorbar(ticks=range(0,13))
            cb.ax.set_yticklabels(["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC","JAN"])
            cb_ax = cb.ax
            cb_ax.text(-1,1.025,'Phase (Months)', fontsize = 18)
        obs_max = np.max(obs_mag)+1
        model_max = np.max(model_mag)+1
        if np.min(obs_mag) < np.min(model_mag):
            min_val = np.min(obs_mag)-1
        else:
            min_val = np.min(model_mag)-1
        if obs_max > model_max:
            plt.xlim(min_val,obs_max)
            plt.ylim(min_val,obs_max)
        else:
            plt.xlim(min_val,model_max)
            plt.ylim(min_val,model_max)
        x = np.arange(0,1000,1)
        y = np.arange(0,1000,1)
		
        handles, labels = plt.gca().get_legend_handles_labels()
        hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
        handles, labels = zip(*hl)
        by_label = OrderedDict(zip(labels, handles))
        leg = ax.legend(by_label.values(), by_label.keys(), loc = 0, prop={'size':20})
        plt.plot(x,y,'k--',alpha = 0.5)
        plt.xlabel('GEOS-Chem 4x5 Magnitude (ppbv)',fontsize = 20)
        plt.ylabel('Obs. Magnitude (ppbv)',fontsize = 20)
        plt.title(r'%s Surface $O_3$ Magnitude, X-Y Plot of Model v Obs.'%title, fontsize= 22)

if plot_type == 'phase':
	if period == 'daily':
		for t in range(len(obs_refs)):
			if third_var == 'yes':
				plt.scatter(model_phase[t], obs_phase[t], c = obs_mag[t], cmap=plt.cm.coolwarm, marker = loc_shapes[tags[t]], s = 100, norm=LogNorm(vmin = np.min(obs_mag), vmax= np.max(obs_mag)), label = loc_dict[tags[t]])	
			else:
				plt.plot(model_phase[t], obs_phase[t], color = loc_colors[tags[t]], marker = 'o', markersize = 15, label = loc_dict[tags[t]])
		#plot site labels
		for i, txt in enumerate(obs_refs):
			x_pos = (model_phase[i]+0.25)
			y_pos = (obs_phase[i]-0.13)
			plt.annotate(txt, (x_pos,y_pos) , fontweight='bold')
		
		if third_var == 'yes':
			cb = plt.colorbar(ticks =[np.min(obs_mag),0.1,0.25,0.5,1,2,4,8,16,np.max(obs_mag)],format='%.2f')	
			cb_ax = cb.ax
			cb_ax.text(-1,1.025,'Magnitude (ppbv)', fontsize = 18)

		plt.xlim(phase_min,phase_max)
		plt.ylim(phase_min,phase_max)
		x = np.arange(0,50,1)
		y = np.arange(0,50,1)
		
		handles, labels = plt.gca().get_legend_handles_labels()
		hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
		handles, labels = zip(*hl)
		by_label = OrderedDict(zip(labels, handles))
		leg = ax.legend(by_label.values(), by_label.keys(), loc = 6, prop={'size':20})
		plt.xticks(range(phase_min,phase_max+1))
		plt.yticks(range(phase_min,phase_max+1))
		plt.plot(x,y,'k--',alpha = 0.5)
		plt.xlabel('GEOS-Chem 4x5 Phase (Hours)',fontsize = 20)
		plt.ylabel('Obs. Phase (Hours)',fontsize = 20)
		plt.title(r'%s Surface $O_3$ Phase, X-Y Plot of Model v Obs.'%title, fontsize= 22)
		
	elif period == 'annual':
		#change phase to run from sep to sep, instead of jan to jan
		obs_phase, model_phase = modules.annual_phase_shift(obs_phase,model_phase)
	
		for t in range(len(obs_refs)):
			if third_var == 'yes':
				plt.scatter(model_phase[t], obs_phase[t], c = obs_mag[t], cmap=plt.cm.coolwarm, marker = loc_shapes[tags[t]], s = 100, vmin = np.min(obs_mag), vmax= np.max(obs_mag), label = loc_dict[tags[t]])	
			else:
				plt.plot(model_phase[t], obs_phase[t], color = loc_colors[tags[t]], marker = 'o', markersize = 15, label = loc_dict[tags[t]])
		#plot site labels
		for i, txt in enumerate(obs_refs):
			x_pos = (model_phase[i]+0.18)
			y_pos = (obs_phase[i]-0.06)
			plt.annotate(txt, (x_pos,y_pos) , fontweight='bold')
		
		if third_var == 'yes':
			cb = plt.colorbar(ticks =[np.min(obs_mag),1,2.5,5,7.5,10,12.5,15,np.max(obs_mag)],format='%.2f')	
			cb_ax = cb.ax
			cb_ax.text(-1,1.025,'Magnitude (ppbv)', fontsize = 18)

		plt.xlim(phase_min,phase_max)
		plt.ylim(phase_min,phase_max)
		x = np.arange(0,50,1)
		y = np.arange(0,50,1)
		
		handles, labels = plt.gca().get_legend_handles_labels()
		hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
		handles, labels = zip(*hl)
		by_label = OrderedDict(zip(labels, handles))
		leg = ax.legend(by_label.values(), by_label.keys(), loc = 0, prop={'size':20})
		plt.xticks(range(phase_min,phase_max+1), ["SEP","OCT","NOV","DEC","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP"])
		plt.yticks(range(phase_min,phase_max+1), ["SEP","OCT","NOV","DEC","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP"])
        #plt.xticks(range(15,phase_max+1), ["MAR","APR","MAY","JUN","JUL","AUG","SEP"])
        #plt.yticks(range(15,phase_max+1), ["MAR","APR","MAY","JUN","JUL","AUG","SEP"])
        plt.plot(x,y,'k--',alpha = 0.5)
        plt.xlabel('GEOS-Chem 4x5 Phase (Months)',fontsize = 20)
        plt.ylabel('Obs. Phase (Months)',fontsize = 20)
        plt.title(r'%s Surface $O_3$ Phase, X-Y Plot of Model v Obs.'%title, fontsize= 22)		

if plot_type == 'mag':
	slope, intercept, r_value, p_value, std_err = stats.linregress(model_mag,obs_mag)
	plt.annotate('Gradient = %.4f'%slope, xy=(.78,.01), xycoords='axes fraction', fontsize =20)

if plot_type == 'phase':
	slope, intercept, r_value, p_value, std_err = stats.linregress(model_phase,obs_phase)
	plt.annotate('Gradient = %.4f'%slope, xy=(.78,.01), xycoords='axes fraction', fontsize =20)
	line = slope*model_phase+intercept
	plt.plot(model_phase,line)
plt.show()
