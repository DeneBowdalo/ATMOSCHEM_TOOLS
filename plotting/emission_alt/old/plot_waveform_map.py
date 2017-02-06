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
import matplotlib.dates as dates
import psutil
import gc
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from matplotlib import dates

pt = 'seasonal'

species = 'NO2-MOLYBDENUM'

start_year = 2009
end_year = 2011

all_areas = ['SW_NA','C_NA','NW_NA','NE_NA','CE_NA','SE_NA','S_NA','NW_EU','C_EU','N_EU','E_EU','S_EU','SW_EU','NE_AS']
areas = [['SW_NA','C_NA','NW_NA','NE_NA','CE_NA','SE_NA','S_NA'],['NW_EU','C_EU','N_EU','E_EU','S_EU','SW_EU'],['NE_AS']]
area_boundaries,area_tags,area_labels = modules.area_dicts()

color_dict = {'SW_NA':0,'C_NA':2,'NW_NA':4,'NE_NA':6,'CE_NA':8,'SE_NA':10,'S_NA':12,'NW_EU':1,'C_EU':3,'N_EU':13,'E_EU':7,'S_EU':9,'SW_EU':11,'NE_AS':5}

#read in obs ts data
obs_fname = '/work/home/db876/observations/surface/%s/process/GLOBAL_SURFACE_%s_2009_2011_H_HP.nc'%(species,species)
obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_var,obs_lats,obs_lons,obs_alt,obs_groups,obs_raw_class,obs_anthrome_class,gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)

#read in std model data                                                                                                                                  
model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_2009_2011_v1001_4x5_GEOS5_H_STD.nc'
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

obs_grp = Dataset('/work/home/db876/xy/%s/2009_2011/2009_2011/obs_SURFACE_H/LSP_stats.nc'%(species))
area_boundaries,area_tags,area_labels = modules.area_dicts()
model_grp = Dataset('/work/home/db876/xy/%s/2009_2011/2009_2011/GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_STD/LSP_stats.nc'%(species))

if pt == 'seasonal':
    obs_waveforms = obs_grp.variables['seasonal_waveform'][:]
    model_waveforms = model_grp.variables['seasonal_waveform'][:]
    obs_wf = np.empty((len(all_areas),8766))
    model_wf = np.empty((len(all_areas),8766))

elif pt == 'diurnal':
    obs_waveforms = obs_grp.variables['diurnal_waveform_winter'][:]
    model_waveforms = model_grp.variables['diurnal_waveform_winter'][:] 
    obs_wf = np.empty((len(all_areas),24))
    model_wf = np.empty((len(all_areas),24))

#iterate through areas
for a in range(len(all_areas)):

    area = all_areas[a]
    print area

    area_grid = area_boundaries[area]
    area_tag = area_tags[area]
    area_label = area_labels[area]

    cut_test = modules.area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag)
    #if dont have any obs for area, then put nans into arrays
    if np.all(cut_test == False):
        print 'none'
        obs_wf[a,:] = np.NaN
        model_wf[a,:] = np.NaN

    else:
        #get cut obs/model waveforms for each area
        obs_sites = obs_waveforms[cut_test,:]
        obs_wf[a,:] = np.nanmean(obs_sites,axis=0)
    
        model_sites = model_waveforms[cut_test,:]
        model_wf[a,:] = np.nanmean(model_sites,axis=0)


cmap = cm.get_cmap('Set1',13)
ratio = 1./13
lw = 4

fig =plt.figure(figsize=(17.5,13.5))
fig.patch.set_facecolor('white')
gs1 = gridspec.GridSpec(1, 1)
gs1.update(top=0.92,bottom=0.35,left=0.02,right=0.48,wspace=0,hspace=0)
ax1 = plt.subplot(gs1[0, 0])
gs2 = gridspec.GridSpec(1, 1)
gs2.update(top=0.92,bottom=0.35,left=0.52,right=0.98,wspace=0,hspace=0)
ax2 = plt.subplot(gs2[0, 0])
gs3 = gridspec.GridSpec(1, 1)
gs3.update(top=0.35,bottom=0.05,left=0.02,right=0.24,wspace=0,hspace=0)
ax3 = plt.subplot(gs3[0, 0])

x_dif = 0.16
y_dif = 0.13

#SW NA
area_i = color_dict['SW_NA']
x1 = 0.02
y1 = 0.36
ax4 = plt.axes([x1, y1, x_dif, y_dif])
ax4.spines['left'].set_color(cmap(area_i))
ax4.spines['left'].set_linewidth(lw)
ax4.spines['bottom'].set_color(cmap(area_i))
ax4.spines['bottom'].set_linewidth(lw)
ax4.spines['right'].set_color(cmap(area_i))
ax4.spines['right'].set_linewidth(lw)
ax4.spines['top'].set_color(cmap(area_i))
ax4.spines['top'].set_linewidth(lw)

#C NA
area_i = color_dict['C_NA']
x1 = 0.20
y1 = 0.79
ax5 = plt.axes([x1, y1, x_dif, y_dif])
ax5.spines['left'].set_color(cmap(area_i))
ax5.spines['left'].set_linewidth(lw)
ax5.spines['bottom'].set_color(cmap(area_i))
ax5.spines['bottom'].set_linewidth(lw)
ax5.spines['right'].set_color(cmap(area_i))
ax5.spines['right'].set_linewidth(lw)
ax5.spines['top'].set_color(cmap(area_i))
ax5.spines['top'].set_linewidth(lw)

#NW NA
area_i = color_dict['NW_NA']
x1 = 0.02
y1 = 0.86
ax6 = plt.axes([x1, y1, x_dif, y_dif])
ax6.spines['left'].set_color(cmap(area_i))
ax6.spines['left'].set_linewidth(lw)
ax6.spines['bottom'].set_color(cmap(area_i))
ax6.spines['bottom'].set_linewidth(lw)
ax6.spines['right'].set_color(cmap(area_i))
ax6.spines['right'].set_linewidth(lw)
ax6.spines['top'].set_color(cmap(area_i))
ax6.spines['top'].set_linewidth(lw)


#NE NA
area_i = color_dict['NE_NA']
x1 = 0.38
y1 = 0.75
ax7 = plt.axes([x1, y1, x_dif, y_dif])
ax7.spines['left'].set_color(cmap(area_i))
ax7.spines['left'].set_linewidth(lw)
ax7.spines['bottom'].set_color(cmap(area_i))
ax7.spines['bottom'].set_linewidth(lw)
ax7.spines['right'].set_color(cmap(area_i))
ax7.spines['right'].set_linewidth(lw)
ax7.spines['top'].set_color(cmap(area_i))
ax7.spines['top'].set_linewidth(lw)


#CE NA
area_i = color_dict['CE_NA']
x1 = 0.415
y1 = 0.47
ax8 = plt.axes([x1, y1, x_dif, y_dif])
ax8.spines['left'].set_color(cmap(area_i))
ax8.spines['left'].set_linewidth(lw)
ax8.spines['bottom'].set_color(cmap(area_i))
ax8.spines['bottom'].set_linewidth(lw)
ax8.spines['right'].set_color(cmap(area_i))
ax8.spines['right'].set_linewidth(lw)
ax8.spines['top'].set_color(cmap(area_i))
ax8.spines['top'].set_linewidth(lw)


#SE NA
area_i = color_dict['SE_NA']
x1 = 0.38
y1 = 0.315
ax9 = plt.axes([x1, y1, x_dif, y_dif])
ax9.spines['left'].set_color(cmap(area_i))
ax9.spines['left'].set_linewidth(lw)
ax9.spines['bottom'].set_color(cmap(area_i))
ax9.spines['bottom'].set_linewidth(lw)
ax9.spines['right'].set_color(cmap(area_i))
ax9.spines['right'].set_linewidth(lw)
ax9.spines['top'].set_color(cmap(area_i))
ax9.spines['top'].set_linewidth(lw)


#S NA
area_i = color_dict['S_NA']
x1 = 0.20
y1 = 0.29
ax10 = plt.axes([x1, y1, x_dif, y_dif])
ax10.spines['left'].set_color(cmap(area_i))
ax10.spines['left'].set_linewidth(lw)
ax10.spines['bottom'].set_color(cmap(area_i))
ax10.spines['bottom'].set_linewidth(lw)
ax10.spines['right'].set_color(cmap(area_i))
ax10.spines['right'].set_linewidth(lw)
ax10.spines['top'].set_color(cmap(area_i))
ax10.spines['top'].set_linewidth(lw)


#NW EU
area_i = color_dict['NW_EU']
x1 = 0.56
y1 = 0.78
ax11 = plt.axes([x1, y1, x_dif, y_dif])
ax11.spines['left'].set_color(cmap(area_i))
ax11.spines['left'].set_linewidth(lw)
ax11.spines['bottom'].set_color(cmap(area_i))
ax11.spines['bottom'].set_linewidth(lw)
ax11.spines['right'].set_color(cmap(area_i))
ax11.spines['right'].set_linewidth(lw)
ax11.spines['top'].set_color(cmap(area_i))
ax11.spines['top'].set_linewidth(lw)

#C EU
area_i = color_dict['C_EU']
x1 = 0.74
y1 = 0.34
ax12 = plt.axes([x1, y1, x_dif, y_dif])
ax12.spines['left'].set_color(cmap(area_i))
ax12.spines['left'].set_linewidth(lw)
ax12.spines['bottom'].set_color(cmap(area_i))
ax12.spines['bottom'].set_linewidth(lw)
ax12.spines['right'].set_color(cmap(area_i))
ax12.spines['right'].set_linewidth(lw)
ax12.spines['top'].set_color(cmap(area_i))
ax12.spines['top'].set_linewidth(lw)


#N EU
area_i = color_dict['N_EU']
x1 = 0.74
y1 = 0.86
ax13 = plt.axes([x1, y1, x_dif, y_dif])
ax13.spines['left'].set_color(cmap(area_i))
ax13.spines['left'].set_linewidth(lw)
ax13.spines['bottom'].set_color(cmap(area_i))
ax13.spines['bottom'].set_linewidth(lw)
ax13.spines['right'].set_color(cmap(area_i))
ax13.spines['right'].set_linewidth(lw)
ax13.spines['top'].set_color(cmap(area_i))
ax13.spines['top'].set_linewidth(lw)


#E EU
area_i = color_dict['E_EU']
x1 = 0.84
y1 = 0.53
ax14 = plt.axes([x1, y1, x_dif, y_dif])
ax14.spines['left'].set_color(cmap(area_i))
ax14.spines['left'].set_linewidth(lw)
ax14.spines['bottom'].set_color(cmap(area_i))
ax14.spines['bottom'].set_linewidth(lw)
ax14.spines['right'].set_color(cmap(area_i))
ax14.spines['right'].set_linewidth(lw)
ax14.spines['top'].set_color(cmap(area_i))
ax14.spines['top'].set_linewidth(lw)


#S EU
area_i = color_dict['S_EU']
x1 = 0.79
y1 = 0.18
ax15 = plt.axes([x1, y1, x_dif, y_dif])
ax15.spines['left'].set_color(cmap(area_i))
ax15.spines['left'].set_linewidth(lw)
ax15.spines['bottom'].set_color(cmap(area_i))
ax15.spines['bottom'].set_linewidth(lw)
ax15.spines['right'].set_color(cmap(area_i))
ax15.spines['right'].set_linewidth(lw)
ax15.spines['top'].set_color(cmap(area_i))
ax15.spines['top'].set_linewidth(lw)

#SW EU
area_i = color_dict['SW_EU']
x1 = 0.56
y1 = 0.26
ax16 = plt.axes([x1, y1, x_dif, y_dif])
ax16.spines['left'].set_color(cmap(area_i))
ax16.spines['left'].set_linewidth(lw)
ax16.spines['bottom'].set_color(cmap(area_i))
ax16.spines['bottom'].set_linewidth(lw)
ax16.spines['right'].set_color(cmap(area_i))
ax16.spines['right'].set_linewidth(lw)
ax16.spines['top'].set_color(cmap(area_i))
ax16.spines['top'].set_linewidth(lw)


#NE_AS
area_i = color_dict['NE_AS']
x1 = 0.18
y1 = 0.08
ax17 = plt.axes([x1, y1, x_dif, y_dif])
ax17.spines['left'].set_color(cmap(area_i))
ax17.spines['left'].set_linewidth(lw)
ax17.spines['bottom'].set_color(cmap(area_i))
ax17.spines['bottom'].set_linewidth(lw)
ax17.spines['right'].set_color(cmap(area_i))
ax17.spines['right'].set_linewidth(lw)
ax17.spines['top'].set_color(cmap(area_i))
ax17.spines['top'].set_linewidth(lw)


axes = [ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14,ax15,ax16,ax17]

#----------------------------------------------
#plot map areas and points
ax_list = [ax1,ax2,ax3]

latlower_setup = [23,32.5,30]
latupper_setup = [68,72,51.5]
lonwest_setup = [-135,-15,126]
loneast_setup = [-56,35,151]
sects = ['NA','EU','EA']

area_count = 0
for count,ax in enumerate(ax_list):
    area_list = areas[count]

    m = Basemap(projection='cyl',llcrnrlat=latlower_setup[count],urcrnrlat=latupper_setup[count],llcrnrlon=lonwest_setup[count],urcrnrlon=loneast_setup[count],resolution='c',ax = ax)                               

    m.shadedrelief()

    for x,area in enumerate(area_list):
        area_grid = area_boundaries[area]
        area_tag = area_tags[area]
        area_label = area_labels[area]
        cut_test = modules.area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag)

        current_lons = obs_lons[cut_test]
        current_lats = obs_lats[cut_test]
        m.plot(current_lons,current_lats,marker='o',color=cmap(color_dict[area]),linestyle='None')
        area_count+=1

for ax_count,area in enumerate(all_areas):
    area_label = area_labels[area]
    
    #if any NaNs for a area, dont plot
    #turn axis off also
    skip_area = False
    for yy in obs_wf[ax_count,:]:
        if np.isnan(yy) == True:
            skip_area = True
    if skip_area == True:
        axes[ax_count].axis('off')
    else:
        #plot obs and v10 model
        if pt == 'seasonal':
            m_dt = pd.date_range(start = obs_datetime_time[0],end = obs_datetime_time[8766], freq = 'H')[:-1]
        elif pt == 'diurnal':
            m_dt = pd.date_range(start = obs_datetime_time[0],end = obs_datetime_time[24], freq = 'H')[:-1]

        axes[ax_count].plot(m_dt,obs_wf[ax_count,:],color='black')
        axes[ax_count].plot(m_dt,model_wf[ax_count,:],color='red')
    
        if pt == 'seasonal':
            period = dates.MonthLocator()
            dfmt = dates.DateFormatter('%m')
            axes[ax_count].xaxis.set_major_locator(period)
            axes[ax_count].xaxis.set_major_formatter(dfmt)
        elif pt == 'diurnal':
            period = dates.HourLocator()
            dfmt = dates.DateFormatter('%H')
            axes[ax_count].xaxis.set_major_locator(period)
            axes[ax_count].xaxis.set_major_formatter(dfmt)

x_dif = 0.2
y_dif = 0.2
x1 = 0.41
y1 = 0.11
ax18 = plt.axes([x1, y1, x_dif, y_dif])
a = ax18.scatter([0],[0],s=0.001,c='black')
b = ax18.scatter([0],[0],s=0.001,c='red')
ax18.axis('off')
ax18.axes.get_xaxis().set_visible(False)
ax18.axes.get_yaxis().set_visible(False)
leg = ax18.legend([a,b],['Observations','GEOS-Chem 4x5 v1001 GEOS5'], loc = 'center', prop={'size':20},markerscale=400,ncol=1,frameon=False)
leg.get_frame().set_facecolor('white')

plt.show()
