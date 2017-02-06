fimport numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
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
import pandas as pd
from bokeh.plotting import figure,output_notebook,show,ColumnDataSource
from bokeh.models import HoverTool,PanTool,WheelZoomTool,BoxZoomTool,ResetTool,ResizeTool,TapTool,CustomJS,HBox,VBox,Select
from bokeh.models.widgets import Dropdown
from bokeh.io import vform
from bokeh.models.glyphs import Text

#set up writing to ipython notebook
output_notebook()

#-----------------------------------------------------
#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres = modules.get_obs_info(present_dir)
model_fname,species,start_year,end_year = modules.get_model_info(present_dir)
    
#read in obs_ts_data
obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_raw_class,obs_anthrome_class,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

#get obs lat_lon grid central points
obs_lats_centre, obs_lons_centre, model_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_lats,obs_lons) 

#get observational location tags 
tags = modules.get_tags(obs_refs)

#--------------------------------------------------------
#load in periodic lsp data
obs_periodic_fname = '../obs_%s_%s/obs_sig_periods.nc'%(vres,timeres)
model_periodic_fname = 'model_sig_periods.nc'

obs_diurnal_mag ,obs_diurnal_ph,obs_seasonal_mag,obs_seasonal_ph,obs_mean,obs_p1,
obs_p5,obs_p25,obs_p50,obs_p75,obs_p95,obs_p99,obs_diurnal_ave_waveform,
obs_diurnal_ave_waveform_extended,obs_diurnal_season_waveform_extended,obs_seasonal_waveform,
obs_seasonal_waveform_extended,obs_full_ave_waveform,obs_full_season_waveform,obs_pc_var_daily, 
obs_pc_var_seasonal,obs_pc_var_full,obs_diurnal_mag_spring,obs_diurnal_ph_spring, 
obs_mean_spring,obs_p1_spring,obs_p5_spring,obs_p25_spring,obs_p50_spring,obs_p75_spring, 
obs_p95_spring,obs_p99_spring,obs_diurnal_waveform_spring,obs_pc_var_daily_spring,obs_diurnal_mag_summer, 
obs_diurnal_ph_summer,obs_mean_summer,obs_p1_summer,obs_p5_summer,obs_p25_summer,obs_p50_summer, 
obs_p75_summer,obs_p95_summer,obs_p99_summer,obs_diurnal_waveform_summer,obs_pc_var_daily_summer, 
obs_diurnal_mag_autumn,obs_diurnal_ph_autumn,obs_mean_autumn,obs_p1_autumn,obs_p5_autumn, 
obs_p25_autumn,obs_p50_autumn,obs_p75_autumn,obs_p95_autumn,obs_p99_autumn,obs_diurnal_waveform_autumn, 
obs_pc_var_daily_autumn,obs_diurnal_mag_winter,obs_diurnal_ph_winter,obs_mean_winter, 
obs_p1_winter,obs_p5_winter,obs_p25_winter,obs_p50_winter,obs_p75_winter,obs_p95_winter,           
obs_p99_winter,obs_diurnal_waveform_winter,obs_pc_var_daily_winter,obs_seasonal_mag_day, 
obs_seasonal_ph_day,v_mean_day,obs_p1_day,obs_p5_day,obs_p25_day,obs_p50_day, 
obs_p75_day,obs_p95_day,obs_p99_day,obs_seasonal_waveform_day,obs_pc_var_seasonal_day,obs_seasonal_mag_night,
obs_seasonal_ph_night,obs_mean_night,obs_p1_night,obs_p5_night,obs_p25_night,obs_p50_night,obs_p75_night,
obs_p95_night,obs_p99_night,obs_seasonal_waveform_night,obs_pc_var_seasonal_night,obs_daily_h3_mag,
obs_daily_h2_mag,obs_daily_h1_mag,obs_daily_mag,obs_annual_h3_mag,obs_annual_h2_mag,obs_annual_h1_mag,obs_annual_mag,
model_diurnal_mag ,model_diurnal_ph,model_seasonal_mag,model_seasonal_ph,model_mean,model_p1,
model_p5,model_p25,model_p50,model_p75,model_p95,model_p99,model_diurnal_ave_waveform,
model_diurnal_ave_waveform_extended,model_diurnal_season_waveform_extended,model_seasonal_waveform,
model_seasonal_waveform_extended,model_full_ave_waveform,model_full_season_waveform,model_pc_var_daily, 
model_pc_var_seasonal,model_pc_var_full,model_diurnal_mag_spring,model_diurnal_ph_spring, 
model_mean_spring,model_p1_spring,model_p5_spring,model_p25_spring,model_p50_spring,model_p75_spring, 
model_p95_spring,model_p99_spring,model_diurnal_waveform_spring,model_pc_var_daily_spring,model_diurnal_mag_summer, 
model_diurnal_ph_summer,model_mean_summer,model_p1_summer,model_p5_summer,model_p25_summer,model_p50_summer, 
model_p75_summer,model_p95_summer,model_p99_summer,model_diurnal_waveform_summer,model_pc_var_daily_summer, 
model_diurnal_mag_autumn,model_diurnal_ph_autumn,model_mean_autumn,model_p1_autumn,model_p5_autumn, 
model_p25_autumn,model_p50_autumn,model_p75_autumn,model_p95_autumn,model_p99_autumn,model_diurnal_waveform_autumn, 
model_pc_var_daily_autumn,model_diurnal_mag_winter,model_diurnal_ph_winter,model_mean_winter, 
model_p1_winter,model_p5_winter,model_p25_winter,model_p50_winter,model_p75_winter,model_p95_winter,           
model_p99_winter,model_diurnal_waveform_winter,model_pc_var_daily_winter,model_seasonal_mag_day, 
model_seasonal_ph_day,model_mean_day,model_p1_day,model_p5_day,model_p25_day,model_p50_day, 
model_p75_day,model_p95_day,model_p99_day,model_seasonal_waveform_day,model_pc_var_seasonal_day,model_seasonal_mag_night,
model_seasonal_ph_night,model_mean_night,model_p1_night,model_p5_night,model_p25_night,model_p50_night,model_p75_night,
model_p95_night,model_p99_night,model_seasonal_waveform_night,model_pc_var_seasonal_night,model_daily_h3_mag,
model_daily_h2_mag,model_daily_h1_mag,model_daily_mag,model_annual_h3_mag,model_annual_h2_mag,model_annual_h1_mag,model_annual_mag = modules.get_periodic_specific(obs_periodic_fname,model_periodic_fname)


#group data

#xy_stats
obs_diurnal_mag ,obs_diurnal_ph,obs_seasonal_mag,obs_seasonal_ph,obs_mean,obs_p1,obs_p5,obs_p25,obs_p50,obs_p75,obs_p95,obs_p99,obs_pc_var_daily,obs_pc_var_seasonal,obs_pc_var_full,
obs_diurnal_mag_spring,obs_diurnal_ph_spring,obs_mean_spring,obs_p1_spring,obs_p5_spring,obs_p25_spring,obs_p50_spring,obs_p75_spring, obs_p95_spring,obs_p99_spring,obs_pc_var_daily_spring,
obs_diurnal_mag_summer, obs_diurnal_ph_summer,obs_mean_summer,obs_p1_summer,obs_p5_summer,obs_p25_summer,obs_p50_summer, obs_p75_summer,obs_p95_summer,obs_p99_summer,obs_pc_var_daily_summer,
obs_diurnal_mag_autumn,obs_diurnal_ph_autumn,obs_mean_autumn,obs_p1_autumn,obs_p5_autumn, obs_p25_autumn,obs_p50_autumn,obs_p75_autumn,obs_p95_autumn,obs_p99_autumn,obs_pc_var_daily_autumn,
obs_diurnal_mag_winter,obs_diurnal_ph_winter,obs_mean_winter, obs_p1_winter,obs_p5_winter,obs_p25_winter,obs_p50_winter,obs_p75_winter,obs_p95_winter,obs_p99_winter,obs_pc_var_daily_winter,
obs_seasonal_mag_day,obs_seasonal_ph_day,v_mean_day,obs_p1_day,obs_p5_day,obs_p25_day,obs_p50_day,obs_p75_day,obs_p95_day,obs_p99_day,obs_pc_var_seasonal_night
#set up xy plot
#s1 = figure(width=1000, plot_height=1000, title=None)
#s1.circle(x, y0, size=20, color="navy", alpha=0.5)

#month_lengths = [0,31,59.25,90.25,120.25,151.25,181.25,212.25,243.25,273.25,304.25,334.25,365.25]
#month_lengths = [273.25,304.25,334.25,365.25,396.25,424.5,455.5,485.5,516.5,546.5,577.5,608.5]
#year_fact = 12./365.25
#month_array  = [i*year_fact for i in month_lengths]
#month_strings = ["JAN","","MAR","","MAY","","JUL","","SEP","","NOV",""]
#axis_nums = month_array
		
#--------------------------------------------------------------------------
#rearrange data so that NA plotted 1, EU plotted 2, AS plotted 3, O plotted 4, OC plotted 5, AF plotted 6, SA plotted 7, ANT plotted 8, ARC plotted 9
# loc_dict = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic Sites','NA':'North America','SA':'South America'}
# order_dict = {'ANT':8,'ARC':9,'AF':6,'AS':3,'EU':2,'OC':5,'O':4,'NA':1,'SA':7}
# 
# order_tags = []
# for tag in tags:
#     n = order_dict[tag]
#     order_tags.append(n)
# 
# tags_argsort = np.argsort(order_tags)
# 
# tags = tags[tags_argsort]
# obs_refs = obs_refs[tags_argsort]
# obs_lats = obs_lats[tags_argsort]
# obs_lons = obs_lons[tags_argsort]
# obs_alt = obs_alt[tags_argsort]
# obs_daily_mag = obs_daily_mag[tags_argsort]
# model_daily_mag = model_daily_mag[tags_argsort]
# obs_daily_phase = obs_daily_phase[tags_argsort]
# model_daily_phase = model_daily_phase[tags_argsort]
# obs_ave = obs_ave[tags_argsort]
# model_ave = model_ave[tags_argsort]
# model_seasonal_mag = model_seasonal_mag[tags_argsort]
# model_seasonal_phase = model_seasonal_phase[tags_argsort]
# obs_seasonal_mag = obs_seasonal_mag[tags_argsort]
# obs_seasonal_phase = obs_seasonal_phase[tags_argsort]
# obs_daily_ff = obs_daily_ff[tags_argsort]
# obs_seasonal_ff = obs_seasonal_ff[tags_argsort]
# model_daily_ff = model_daily_ff[tags_argsort]
# model_seasonal_ff = model_seasonal_ff[tags_argsort]
# obs_daily_waveform = np.array(obs_daily_waveform)[tags_argsort]
# obs_seasonal_waveform = np.array(obs_seasonal_waveform)[tags_argsort]
# obs_full_waveform = np.array(obs_full_waveform)[tags_argsort]
# model_daily_waveform = np.array(model_daily_waveform)[tags_argsort]
# model_seasonal_waveform = np.array(model_seasonal_waveform)[tags_argsort]
# model_full_waveform = np.array(model_full_waveform)[tags_argsort]
# #---------------------------------------------------------------------------
# 
# #setup shape dict for plotting
# loc_colors = {'ANT':'magenta','ARC':'purple','AF':'orange','AS':'blue','EU':'red','OC':'black','O':'yellow','NA':'green','SA':'pink'}
# loc_shapes = {'ANT':'o','ARC':'+','AF':'*','AS':'d','EU':'s','OC':'>','O':'v','NA':'<','SA':'^'}
# 
# #--------------------------------------------------------------------
# #set up plot
# fig, ((ax, ax2,ax3), (ax4, ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3, 3,figsize=(14,12))
# fig.patch.set_facecolor('white')
# fig.subplots_adjust(left=0.18,bottom=0.12)
# 
# marker_size = 6
# font_size = 10
# pad_size = 5
# 
# #get colors and labels for site points
# plot_colors = []
# 
# for t in range(len(obs_refs)):
#     plot_colors.append(loc_colors[tags[t]])
#     
# #----------------------------------------------------------------------
# #Daily Amplitude Plot
# 
# for t in range(len(obs_refs)):
#     line, = ax.plot(obs_daily_mag[t], model_daily_mag[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
# 
# ax.plot(obs_daily_mag, model_daily_mag, color = 'black', marker = 'o', markersize = 1, linestyle= 'None',zorder=1, picker = 5)
# 
# obs_max = np.max(obs_daily_mag)+2
# model_max = np.max(model_daily_mag)+2
# if np.min(obs_daily_mag) < np.min(model_daily_mag):
#     min_val = np.min(obs_daily_mag)-0.001
# else:
#     min_val = np.min(model_daily_mag)-0.001
# if obs_max > model_max:
#     ax.set_xlim(min_val,obs_max)
#     ax.set_ylim(min_val,obs_max)
# else:
#     ax.set_xlim(min_val,model_max)
#     ax.set_ylim(min_val,model_max)
# x = np.arange(0,1000,0.005)
# y = np.arange(0,1000,0.005)
# 
# handles, labels = ax.get_legend_handles_labels()
# hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
# handles, labels = zip(*hl)
# by_label = OrderedDict(zip(labels, handles))
# ax.plot(x,y,'k--',alpha = 0.5)
# ax.set_ylabel('Model (ppb)',fontsize = font_size)
# ax.set_xlabel('Observations (ppb)',fontsize = font_size)
# 
# #----------------------------------------------------------------------
# #Seasonal Amplitude Plot
#  
# for t in range(len(obs_refs)):
#     ax2.plot(obs_seasonal_mag[t], model_seasonal_mag[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
# 
# ax2.plot(obs_seasonal_mag, model_seasonal_mag, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)
# 
# obs_max = np.max(obs_seasonal_mag)+1
# model_max = np.max(model_seasonal_mag)+1
# if np.min(obs_seasonal_mag) < np.min(model_seasonal_mag):
#     min_val = np.min(obs_seasonal_mag)-1
# else:
#     min_val = np.min(model_seasonal_mag)-1
# if obs_max > model_max:
#     ax2.set_xlim(min_val,obs_max)
#     ax2.set_ylim(min_val,obs_max)
# else:
#     ax2.set_xlim(min_val,model_max)
#     ax2.set_ylim(min_val,model_max)
# x = np.arange(0,1000,1)
# y = np.arange(0,1000,1)
# 
# ax2.plot(x,y,'k--',alpha = 0.5)
# ax2.set_ylabel('Model (ppb)',fontsize = font_size)
# ax2.set_xlabel('Observations (ppb)',fontsize = font_size)
# 
# #----------------------------------------------------------------------
# #Daily Phase Plot
# 
# for t in range(len(obs_refs)):
#     ax4.plot(obs_daily_phase[t], model_daily_phase[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
# 
# ax4.plot(obs_daily_phase, model_daily_phase, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)
# 
# ax4.set_xlim(0,24)
# ax4.set_ylim(0,24)
# x = np.arange(0,50,1)
# y = np.arange(0,50,1)
# 
# ax4.set_xticks(range(0,25,2))
# ax4.set_yticks(range(0,25,2))
# 
# ax4.set_ylabel('Model (ppb)',fontsize = font_size)
# ax4.set_xlabel('Observations (Hours)',fontsize = font_size)
# ax4.plot(x,y,'k--',alpha = 0.5)
# for tick in ax4.get_xaxis().get_major_ticks():
#     tick.set_pad(pad_size)
# for tick in ax4.get_yaxis().get_major_ticks():
#     tick.set_pad(pad_size)
# 		
# #----------------------------------------------------------------------
# #Seasonal Phase Plot
# 
# #change phase to run from sep to sep, instead of jan to jan
# #obs_seasonal_phase_corr, model_seasonal_phase_corr = modules.annual_phase_shift(obs_seasonal_phase,model_seasonal_phase)
# 
# for t in range(len(obs_refs)):
#     ax5.plot(obs_seasonal_phase[t], model_seasonal_phase[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
# 
# ax5.plot(obs_seasonal_phase, model_seasonal_phase, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)
# 
# ax5.set_xlim(0,11)
# ax5.set_ylim(0,11)
# #ax5.set_xlim(9,21)
# #ax5.set_ylim(9,21)
# x = np.arange(0,50,1)
# y = np.arange(0,50,1)
# 
# ax5.set_xticks(month_array)
# ax5.set_yticks(month_array)
# ax5.set_xticklabels(month_strings,fontsize=9.2)
# ax5.set_yticklabels(month_strings,fontsize=9.2)
# ax5.set_ylabel('Model (ppb)',fontsize = font_size)
# ax5.set_xlabel('Observations (Months)',fontsize = font_size)
# ax5.plot(x,y,'k--',alpha = 0.5)
# for tick in ax5.get_xaxis().get_major_ticks():
#     tick.set_pad(pad_size)
# for tick in ax5.get_yaxis().get_major_ticks():
#     tick.set_pad(pad_size)		
# 
# #----------------------------------------------------------------------
# #Average Plot   
# 
# for t in range(len(obs_refs)):
#     ax3.plot(obs_ave[t], model_ave[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
# 
# ax3.plot(obs_ave, model_ave, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)
# 
# obs_max = np.max(obs_ave)+1
# model_max = np.max(model_ave)+1
# if np.min(obs_ave) < np.min(model_ave):
#     min_val = np.min(obs_ave)-1
# else:
#     min_val = np.min(model_ave)-1
# if obs_max > model_max:
#     ax3.set_xlim(min_val,obs_max)
#     ax3.set_ylim(min_val,obs_max)
# else:
#     ax3.set_xlim(min_val,model_max)
#     ax3.set_ylim(min_val,model_max)
# x = np.arange(0,1000,1)
# y = np.arange(0,1000,1)
# 
# ax3.plot(x,y,'k--',alpha = 0.5)
# ax3.set_ylabel('Model (ppb)',fontsize = font_size)
# ax3.set_xlabel('Observations (ppb)',fontsize = font_size)
# 
# #----------------------------------------------------------------------
# #Daily Form Factor Plot   
# 
# for t in range(len(obs_refs)):
#     ax7.plot(obs_daily_ff[t], model_daily_ff[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
# 
# ax7.plot(obs_daily_ff, model_daily_ff, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)
# 
# obs_max = np.max(obs_daily_ff)+0.1
# model_max = np.max(model_daily_ff)+0.1
# if np.min(obs_daily_ff) < np.min(model_daily_ff):
#     min_val = np.min(obs_daily_ff)-0.1
# else:
#     min_val = np.min(model_daily_ff)-0.1
# if obs_max > model_max:
#     ax7.set_xlim(min_val,obs_max)
#     ax7.set_ylim(min_val,obs_max)
# else:
#     ax7.set_xlim(min_val,model_max)
#     ax7.set_ylim(min_val,model_max)
# x = np.arange(0,1000,1)
# y = np.arange(0,1000,1)
# 
# ax7.plot(x,y,'k--')
# 
# ax7.set_ylabel('Model Form Factor',fontsize = font_size)
# ax7.set_xlabel('Observational Form Factor',fontsize = font_size)
#         
# #----------------------------------------------------------------------
# #Seasonal Form Factor Plot   
# 
# for t in range(len(obs_refs)):
#     ax8.plot(obs_seasonal_ff[t], model_seasonal_ff[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
# 
# ax8.plot(obs_seasonal_ff, model_seasonal_ff, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)
# 
# obs_max = np.max(obs_seasonal_ff)+0.1
# model_max = np.max(model_seasonal_ff)+0.1
# if np.min(obs_seasonal_ff) < np.min(model_seasonal_ff):
#     min_val = np.min(obs_seasonal_ff)-0.1
# else:
#     min_val = np.min(model_seasonal_ff)-0.1
# if obs_max > model_max:
#     ax8.set_xlim(min_val,obs_max)
#     ax8.set_ylim(min_val,obs_max)
# else:
#     ax8.set_xlim(min_val,model_max)
#     ax8.set_ylim(min_val,model_max)
# x = np.arange(0,1000,1)
# y = np.arange(0,1000,1)
# 
# ax8.plot(x,y,'k--')
# 
# ax8.set_ylabel('Model Form Factor',fontsize = font_size)
# ax8.set_xlabel('Observational Form Factor',fontsize = font_size)
# 
# #------------------------------------------------
# #plot big labels
# plt.figtext(0.24, 0.92, 'Diurnal', fontsize=22)
# plt.figtext(0.50, 0.92, 'Seasonal', fontsize=22)
# plt.figtext(0.76, 0.92, 'Average', fontsize=22)
# 
# plt.figtext(0.01, 0.77, 'Amplitude', fontsize=22)
# plt.figtext(0.01, 0.50, 'Phase', fontsize=22)
# plt.figtext(0.01, 0.23, 'Form Factor', fontsize=22)
# 
# #--------------------------------------------------
# 
# #make axis labels tight to plots
# ax.yaxis.labelpad = 0 
# ax2.yaxis.labelpad = 0
# ax3.yaxis.labelpad = 0
# ax4.yaxis.labelpad = 0
# ax5.yaxis.labelpad = 0
# ax7.yaxis.labelpad = 0
# ax8.yaxis.labelpad = 0
# 
# ax.xaxis.labelpad = 0 
# ax2.xaxis.labelpad = 0
# ax3.xaxis.labelpad = 0
# ax4.xaxis.labelpad = 0
# ax5.xaxis.labelpad = 0
# ax7.xaxis.labelpad = 0
# ax8.xaxis.labelpad = 0
# 
# #plot legend
# handles, labels = ax5.get_legend_handles_labels()
# hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
# handles, labels = zip(*hl)
# by_label = OrderedDict(zip(labels, handles))
# leg = ax.legend(by_label.values(), by_label.keys(), loc = 'upper center', bbox_to_anchor=(1.7,-2.55),fancybox=True,ncol=4)
# 
# #handles, labels = ax5.get_legend_handles_labels()
# #lgd = ax5.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.45,-1.35),fancybox=True,ncol=5)
# 
# #--------------------------------------------------------------------------------
# #remove ax 8 and ax9 for table
# 
# #ax3.set_frame_on(False)
# ax6.set_frame_on(False)
# #ax8.set_frame_on(False)
# ax9.set_frame_on(False)
# #ax3.axes.get_yaxis().set_visible(False)
# ax6.axes.get_yaxis().set_visible(False)
# #ax8.axes.get_yaxis().set_visible(False)
# ax9.axes.get_yaxis().set_visible(False)
# #ax3.axes.get_xaxis().set_visible(False)
# ax6.axes.get_xaxis().set_visible(False)
# #ax8.axes.get_xaxis().set_visible(False)
# ax9.axes.get_xaxis().set_visible(False)
# 
# 
# mng = plt.get_current_fig_manager()
# mng.window.wm_geometry("+2500-800")
# 
# fig.canvas.mpl_connect('pick_event', interactive)
# 
# plt.show()


