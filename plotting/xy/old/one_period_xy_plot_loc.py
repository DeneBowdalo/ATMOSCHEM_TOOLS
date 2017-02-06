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


#loc group names
loc_names = ['NA','EU','AS','ROW']

def interactive(event):
    modules.clicker_interactive_map_obsmodel(event,species,lat_e,lon_e,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_daily_waveform,obs_seasonal_waveform,obs_full_waveform,model_daily_waveform,model_seasonal_waveform,model_full_waveform,fig,all_m)

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
model_seasonal_mag = []
model_seasonal_phase = []
model_seasonal_phase_min = []
obs_daily_mag = []
obs_daily_phase = []
obs_daily_phase_min = []
obs_seasonal_mag = []
obs_seasonal_phase = []
obs_seasonal_phase_min = []
model_ave = []
obs_ave = []
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
    model_seasonal_mag = np.append(model_seasonal_mag,site_group_mod.seasonal_amplitude)
    model_seasonal_phase = np.append(model_seasonal_phase,site_group_mod.seasonal_phase)
    model_seasonal_phase_min = np.append(model_seasonal_phase_min,site_group_mod.seasonal_phase_min)
    model_ave = np.append(model_ave,site_group_mod.average)
    model_daily_ff = np.append(model_daily_ff,site_group_mod.daily_ff) 
    model_seasonal_ff = np.append(model_seasonal_ff,site_group_mod.seasonal_ff)
    model_daily_waveform.append(site_group_mod.variables['daily_waveform'][:])
    model_seasonal_waveform.append(site_group_mod.variables['seasonal_waveform'][:]) 
    model_full_waveform.append(site_group_mod.variables['all_waveform'][:]) 
    
    obs_daily_mag = np.append(obs_daily_mag,site_group_obs.daily_amplitude)
    obs_daily_phase = np.append(obs_daily_phase,site_group_obs.daily_phase)
    obs_daily_phase_min = np.append(obs_daily_phase_min,site_group_obs.daily_phase_min)
    obs_seasonal_mag = np.append(obs_seasonal_mag,site_group_obs.seasonal_amplitude)
    obs_seasonal_phase = np.append(obs_seasonal_phase,site_group_obs.seasonal_phase)
    obs_seasonal_phase_min = np.append(obs_seasonal_phase_min,site_group_obs.seasonal_phase_min)
    obs_ave = np.append(obs_ave,site_group_obs.average)
    obs_daily_ff = np.append(obs_daily_ff,site_group_obs.daily_ff) 
    obs_seasonal_ff = np.append(obs_seasonal_ff,site_group_obs.seasonal_ff)
    obs_daily_waveform.append(site_group_obs.variables['daily_waveform'][:])
    obs_seasonal_waveform.append(site_group_obs.variables['seasonal_waveform'][:]) 
    obs_full_waveform.append(site_group_obs.variables['all_waveform'][:])

param = raw_input('\nmag, ph_min, ph_max, ave or ff?\n')

obs = []
model = []


if (param != 'ave'):
    period = raw_input('d or s?\n')

        
if param == 'mag':
    if period == 'd':
        current_obs = np.copy(obs_daily_mag)
        current_model =  np.copy(model_daily_mag)
        min = [0,0,0,0]
        max = [20,20,20,5]
    if period == 's':
        current_obs = np.copy(obs_seasonal_mag)
        current_model =  np.copy(model_seasonal_mag)
        min = [0,0,0,0]
        max = [30,30,30,20]

if param == 'ph_max':
    if period == 'd':
        current_obs = np.copy(obs_daily_phase)
        current_model =  np.copy(model_daily_phase)
        min = 0
        max = 24
    if period == 's':
        current_obs = np.copy(obs_seasonal_phase)
        current_model =  np.copy(model_seasonal_phase)
        min = 0
        max = 12
        #change SH seasonal phase to match up with NH
        #test = obs_lats < 0
        #current_obs[test] = current_obs[test] - 6. 
        #current_model[test] = current_model[test] - 6. 
        #test1 = current_obs < 0
        #test2 = current_model < 0
        #current_obs[test1] = 12. - np.abs(current_obs[test1])
        #current_model[test2] = 12. - np.abs(current_model[test2])
        
if param == 'ph_min':
    if period == 'd':
        current_obs = np.copy(obs_daily_phase_min)
        current_model =  np.copy(model_daily_phase_min)
        min = 0
        max = 24
    if period == 's':
        current_obs = np.copy(obs_seasonal_phase_min)
        current_model =  np.copy(model_seasonal_phase_min)
        min = 0
        max = 12
        #change SH seasonal phase to match up with NH
        #test = obs_lats < 0
        #current_obs[test] = current_obs[test] - 6. 
        #current_model[test] = current_model[test] - 6. 
        #test1 = current_obs < 0
        #test2 = current_model < 0
        #current_obs[test1] = 12. - np.abs(current_obs[test1])
        #current_model[test2] = 12. - np.abs(current_model[test2])

if param == 'ff':
    if period == 'd':
        current_obs = np.copy(obs_daily_ff)
        current_model =  np.copy(model_daily_ff)
        min = [1,1,1,1]
        max = [1.3,1.4,1.3,1.3]
    if period == 's':
        current_obs = np.copy(obs_seasonal_ff)
        current_model =  np.copy(model_seasonal_ff)
        min = [1,1,1,1]
        max = [1.5,1.4,1.4,1.3]
      
if param == 'ave':
    current_obs = np.copy(obs_ave)
    current_model = np.copy(model_ave)
    min = [15,15,15,10]
    max = [50,50,50,50]

month_lengths = [0,31,28.25,31,30,31,30,31,31,30,31,30,31]
annual_ratio = 12./365.25
current_days = 0
actual_month_ratio = []
for i in month_lengths:
    current_days = current_days+i
    actual_month_ratio.append(current_days * annual_ratio)


#ints = [0,2,4,6,8,10]
#actual_month_ratio = np.array(actual_month_ratio)
#actual_month_ratio = actual_month_ratio[ints]
month_strings = ["1","2","3","4","5","6","7","8","9","10","11","12","1"]

#---------------------------------------------------------------------------

#setup shape dict for plotting
loc_colors = {'ANT':'magenta','ARC':'purple','AF':'orange','AS':'blue','EU':'red','OC':'black','O':'yellow','NA':'green','SA':'pink'}
loc_shapes = {'ANT':'o','ARC':'+','AF':'*','AS':'d','EU':'s','OC':'>','O':'v','NA':'<','SA':'^'}
loc_dict = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic Sites','NA':'North America','SA':'South America'}
all_locs = ['Antarctica','Arctic','Africa','Asia','Europe','Oceania','Oceanic Sites','North America','South America']

point_sizes = [10,10,15,15]

#--------------------------------------------------------------------
fig, axes = plt.subplots(nrows=3, ncols=2,figsize=(13,16))
fig.patch.set_facecolor('white')
t = 0
actual_t = 0
all_handles = []
all_labels = []

all_m = []

for ax in axes.flat:
    if (t != 4) and (t != 5):
        current_tag = loc_names[actual_t]
    
        if current_tag != 'ROW':
            test = tags == current_tag
        else:
            test = []
            other_tags = ['AF','ANT','ARC','OC','O','SA']
            for i in tags:
                if i in other_tags:
                    test.append(True)
                else:
                    test.append(False)
            test = np.array(test)
    
        cut_obs = current_obs[test]
        cut_model = current_model[test]
        cut_tags = tags[test]
        unique_tags = np.unique(cut_tags)
    
        for tag in unique_tags:
            test = cut_tags == tag
            group_obs = cut_obs[test]
            group_model = cut_model[test]
            line, = ax.plot(group_obs, group_model, color = loc_colors[tag],linestyle = 'None', marker = 'o', markersize = point_sizes[actual_t])
            ax.plot([1],[1],color = loc_colors[tag],linestyle = 'None',marker = 'o', markersize = 0.000001,label = loc_dict[tag],picker = 5)
    
        ax.text(0.03, 0.97, loc_names[actual_t], transform=ax.transAxes,fontsize=28, fontweight='bold', va='top')
        ax.grid(True)

        if param != 'ph':
            ax.set_xlim(min[actual_t],max[actual_t])
            ax.set_ylim(min[actual_t],max[actual_t])
        else:
            ax.set_xlim(min,max)
            ax.set_ylim(min,max)

        if param == 'ph':
            if period == 'd':
                ax.set_xticks([0,4,8,12,16,20,24])
                ax.set_xticklabels(['0','4','8','12','16','20','0'])
                ax.set_yticks([0,4,8,12,16,20,24])
                ax.set_yticklabels(['0','4','8','12','16','20','0'])
            if period == 's':
                ax.set_xticks(actual_month_ratio)
                ax.set_xticklabels(month_strings)
                ax.set_yticks(actual_month_ratio)
                ax.set_yticklabels(month_strings)
            

        x = np.arange(0,1000,0.005)
        y = np.arange(0,1000,0.005)

        ax.plot(x,y,'k--',alpha = 0.5)
    
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.tick_params(axis='both', which='minor', labelsize=20)
    
        for tick in ax.get_xaxis().get_major_ticks():
            tick.set_pad(7.5)
        for tick in ax.get_yaxis().get_major_ticks():
            tick.set_pad(7.5)
    

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


        def orth_regression(model,obs):
            linear = Model(f)
            mydata = RealData(obs, model)
            myodr = ODR(mydata, linear, beta0=[1., 0.])
            myoutput = myodr.run()
            params = myoutput.beta
            gradient = params[0]
            y_intercept = params[1]
            res_var = myoutput.res_var
            return np.around(gradient,2), np.around(y_intercept,2), np.around(res_var,2)

        current_grad,current_y,current_res = orth_regression(cut_obs,cut_model)

        res_sum,res_ave = modules.orthogonal_1_1_res_var(cut_obs,cut_model)

        if actual_t == 0:
            na_grad = np.copy(current_grad)
            na_y = np.copy(current_y)
            na_res = np.copy(res_ave)
        if actual_t == 1:
            eu_grad = np.copy(current_grad)
            eu_y = np.copy(current_y)
            eu_res = np.copy(res_ave)
        if actual_t == 2:
            as_grad = np.copy(current_grad)
            as_y = np.copy(current_y)
            as_res = np.copy(res_ave)
        if actual_t == 3:
            row_grad = np.copy(current_grad)
            row_y = np.copy(current_y)
            row_res = np.copy(res_ave)
        

        handles, labels = ax.get_legend_handles_labels()
        all_handles = np.append(all_handles,handles)
        all_labels = np.append(all_labels,labels)
        
        actual_t+=1    

    if t == 4:
       #remove ax 
        ax.set_frame_on(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_xaxis().set_visible(False)
        
        row_labels=['Gradient','Y-Intercept','Average Residuals']
        col_labels=['NA','EU','AS','ROW']
        table = ax.table(cellText=[[na_grad,eu_grad,as_grad,row_grad],[na_y,eu_y,as_y,row_y],[na_res,eu_res,as_res,row_res]],rowLabels=row_labels,colLabels=col_labels,loc='bottom',bbox = [0.55,0,0.95,0.95])
        table_props=table.properties()
        table_cells=table_props['child_artists']
        for cell in table_cells:
            cell.set_fontsize(22)
        
    if t == 5:
       #remove ax 
        ax.set_frame_on(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_xaxis().set_visible(False)
        
        all_handles = np.ravel(all_handles)
        all_labels = np.ravel(all_labels)
        hl = sorted(zip(all_handles, all_labels), key=operator.itemgetter(1))
        all_handles, all_labels = zip(*hl)
        by_label = OrderedDict(zip(all_labels, all_handles))
        ax.legend(by_label.values(),by_label.keys(), loc='center right',prop={'size':22},fancybox=True,ncol=1,markerscale=20000000)
        

    t+=1 

    all_m.append(ax)
        
    #----------------------------------------------------------------------

plt.tight_layout(pad = 3.08)

if (param == 'ph_max') or (param == 'ph_min'):
    if period == 'd':
        plt.figtext(0.39, 0.315, 'Observational Hour', fontsize=24)
        plt.figtext(0.01, 0.71, 'Model Hour', fontsize=24,rotation=90)
    if period == 's':
        plt.figtext(0.36, 0.315, 'Observational Month', fontsize=24)
        plt.figtext(0.01, 0.71, 'Model Month', fontsize=24,rotation=90)

if param == 'ff':
    plt.figtext(0.39, 0.315, 'Observational Form Factor', fontsize=24)
    plt.figtext(0.01, 0.71, 'Model Form Factor', fontsize=24,rotation=90)
    
else:
    plt.figtext(0.275, 0.315, 'Observational Concentration (ppb)', fontsize=24)
    plt.figtext(0.01, 0.78, 'Model Concentration (ppb)', fontsize=24,rotation=90)

#fig.subplots_adjust(bottom=0.15)

plt.show()


