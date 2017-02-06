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
import matplotlib.dates as dates

present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-4]

start_year = 2009
end_year = 2011

#read in obs ts data
obs_fname = '/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_%s_2009_2011_H_PERIODIC.nc'%(species)
obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)

#read in std model data
model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_2009_2011_v1001_4x5_GEOS5_H_STD.nc'
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

#--------------------------------------------------------
#load in periodic lsp data


obs_grp = Dataset('../obs_SURFACE_H/obs_sig_periods.nc')
std_model_grp = Dataset('../GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_STD/model_sig_periods.nc')

alt_model_dirs = ['GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5NMVOC1.0',
                  'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ANMVOC0.25',
                  'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ANMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ANMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ANMVOC4.0',
                  'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ANMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ANMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0BNMVOC0.25',
                  'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0BNMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25BNMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0BNMVOC4.0',
                  'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25BNMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0BNMVOC0.25']

#nox_array = [0.25,0.5,2.0,4.0,1.0,1.0,0.25,4.0,0.25,4.0,1.0,1.0,0.25,4.0,0.25,4.0]
#anmvoc_array = [1.0,1.0,1.0,1.0,0.25,4.0,0.25,4.0,4.0,0.25]
#bnmvoc_array = [0.25,4.0,0.25,4.0,4.0,0.25]

nox_array = [0,1,3,4,2,2,0,4,0,4,2,2,0,4,0,4]
anmvoc_array = [2,2,2,2,0,4,0,4,4,0]
bnmvoc_array = [0,4,0,4,4,0]

#nmodels = len(alt_model_dirs)
nmodels_linspace = np.linspace(0, 1, 5)

obs_daily_amp = []
obs_daily_ph = []
obs_daily_ff = []
obs_seasonal_amp = []
obs_seasonal_ph = []
obs_seasonal_ff = []
obs_ave = []
std_model_daily_amp = []
std_model_daily_ph = []
std_model_daily_ff = []
std_model_seasonal_amp = []
std_model_seasonal_ph = []
std_model_seasonal_ff = []
std_model_ave = []

for ref in obs_refs:
    site_group_obs = obs_grp.groups[ref]
    site_group_std_model = std_model_grp.groups[ref]    

    obs_daily_amp.append(site_group_obs.daily_amplitude)
    obs_daily_ph.append(site_group_obs.daily_phase)
    obs_daily_ff.append(site_group_obs.daily_ff)
    obs_seasonal_amp.append(site_group_obs.seasonal_amplitude)
    obs_seasonal_ph.append(site_group_obs.seasonal_phase)
    obs_seasonal_ff.append(site_group_obs.seasonal_ff)
    obs_ave.append(site_group_obs.average)
    
    std_model_daily_amp.append(site_group_std_model.daily_amplitude)
    std_model_daily_ph.append(site_group_std_model.daily_phase)
    std_model_daily_ff.append(site_group_std_model.daily_ff)
    std_model_seasonal_amp.append(site_group_std_model.seasonal_amplitude)
    std_model_seasonal_ph.append(site_group_std_model.seasonal_phase)
    std_model_seasonal_ff.append(site_group_std_model.seasonal_ff)
    std_model_ave.append(site_group_std_model.average)

obs_daily_amp = np.array(obs_daily_amp)
obs_daily_ph = np.array(obs_daily_ph)
obs_daily_ff = np.array(obs_daily_ff)
obs_seasonal_amp = np.array(obs_seasonal_amp)
obs_seasonal_ph = np.array(obs_seasonal_ph)
obs_seasonal_ff = np.array(obs_seasonal_ff)
obs_ave = np.array(obs_ave)
std_model_daily_amp = np.array(std_model_daily_amp)
std_model_daily_ph = np.array(std_model_daily_ph)
std_model_daily_ff = np.array(std_model_daily_ff)
std_model_seasonal_amp = np.array(std_model_seasonal_amp)
std_model_seasonal_ph = np.array(std_model_seasonal_ph)
std_model_seasonal_ff = np.array(std_model_seasonal_ff)
std_model_ave = np.array(std_model_ave)

print obs_daily_amp.shape
print std_model_daily_amp.shape
print std_model_seasonal_amp.shape
print std_model_ave.shape

#-----------------------------------
#get area
areas = ['ANT','S_O','OC','AF','SA','NE_NA','CE_NA','SE_NA','S_NA','SW_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','NE_AS','SE_AS','C_AS','S_AS','N_O','AL','ARC']

plot_type = raw_input('\nd, s or ave?\n')

if plot_type != 'ave':
    plot_type_2 = raw_input('\namp, ph or ff?\n')

#if plot_type == 'd':
#    obs_datetimes = obs_datetime_time[:24]
#    model_datetimes = model_datetime_time[:24]
#if plot_type == 's':
#    obs_datetimes = obs_datetime_time[:8766]
#    model_datetimes = model_datetime_time[:8766]

#obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
#model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')

area_boundaries,area_tags,area_labels = modules.area_dicts()

fig, axes = plt.subplots(nrows=5, ncols=5,figsize=(19,13))
fig.patch.set_facecolor('white')

diff_array = np.empty((len(areas),len(alt_model_dirs)+1))

anmvoc_ind = 0
bnmvoc_ind = 0
mix_color_array = []
mix_color_array.append('yellow')

for m in range(len(alt_model_dirs)):
    print m
    print '../%s/model_sig_periods.nc'%(alt_model_dirs[m])
    alt_model_grp = Dataset('../%s/model_sig_periods.nc'%(alt_model_dirs[m]))

    alt_model_daily_amp = []
    alt_model_daily_ph = []
    alt_model_daily_ff = []
    alt_model_seasonal_amp = []
    alt_model_seasonal_ph = []
    alt_model_seasonal_ff = []
    alt_model_ave = []

    for ref in obs_refs:
        site_group_alt_model = alt_model_grp.groups[ref]
    
        alt_model_daily_amp.append(site_group_alt_model.daily_amplitude)
        alt_model_daily_ph.append(site_group_alt_model.daily_phase)
        alt_model_daily_ff.append(site_group_alt_model.daily_ff)
        alt_model_seasonal_amp.append(site_group_alt_model.seasonal_amplitude)
        alt_model_seasonal_ph.append(site_group_alt_model.seasonal_phase)
        alt_model_seasonal_ff.append(site_group_alt_model.seasonal_ff)
        alt_model_ave.append(site_group_alt_model.average)

    alt_model_daily_amp = np.array(alt_model_daily_amp)
    alt_model_daily_ph = np.array(alt_model_daily_ph)
    alt_model_daily_ff = np.array(alt_model_daily_ff)
    alt_model_seasonal_amp = np.array(alt_model_seasonal_amp)
    alt_model_seasonal_ph = np.array(alt_model_seasonal_ph)
    alt_model_seasonal_ff = np.array(alt_model_seasonal_ff)
    alt_model_ave = np.array(alt_model_ave)

    #get mix color 
    nox_color = plt.cm.RdBu_r(nmodels_linspace[nox_array[m]])
    print alt_model_dirs[m]
    if ('ANMVOC' in alt_model_dirs[m]):
        anmvoc_color = plt.cm.PRGn(nmodels_linspace[anmvoc_array[anmvoc_ind]])  
        mix_color = ((nox_color[0]+anmvoc_color[0])/2, (nox_color[1]+anmvoc_color[1])/2, (nox_color[2]+anmvoc_color[2])/2, (nox_color[3]+anmvoc_color[3])/2)
        anmvoc_ind+=1
    elif 'BNMVOC' in alt_model_dirs[m]:
        bnmvoc_color = plt.cm.BrBG_r(nmodels_linspace[bnmvoc_array[bnmvoc_ind]])
        mix_color = ((nox_color[0]+bnmvoc_color[0])/2, (nox_color[1]+bnmvoc_color[1])/2, (nox_color[2]+bnmvoc_color[2])/2, (nox_color[3]+bnmvoc_color[3])/2)
        bnmvoc_ind+=1
    #NMVOC
    else:
        anmvoc_color = plt.cm.PRGn(nmodels_linspace[anmvoc_array[anmvoc_ind]])  
        mix_color = ((nox_color[0]+anmvoc_color[0])/2, (nox_color[1]+anmvoc_color[1])/2, (nox_color[2]+anmvoc_color[2])/2, (nox_color[3]+anmvoc_color[3])/2)
        anmvoc_ind+=1    
    
    mix_color_array.append(mix_color)

    count = 0
    for a in range(len(areas)):
    
        area = areas[a]

        area_grid = area_boundaries[area]
        area_tag = area_tags[area]
        area_label = area_labels[area]

        cut_test = modules.area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag)
        
        obs_d_amp = np.nanmean(obs_daily_amp[cut_test])
        obs_d_ph = np.nanmean(obs_daily_ph[cut_test])
        obs_d_ff = np.nanmean(obs_daily_ff[cut_test])
        obs_s_amp = np.nanmean(obs_seasonal_amp[cut_test])
        obs_s_ph = np.nanmean(obs_seasonal_ph[cut_test])
        obs_s_ff = np.nanmean(obs_seasonal_ff[cut_test])
        obs_a = np.nanmean(obs_ave[cut_test])

        std_model_d_amp = np.nanmean(std_model_daily_amp[cut_test])
        std_model_d_ph = np.nanmean(std_model_daily_ph[cut_test])
        std_model_d_ff = np.nanmean(std_model_daily_ff[cut_test])
        std_model_s_amp = np.nanmean(std_model_seasonal_amp[cut_test])
        std_model_s_ph = np.nanmean(std_model_seasonal_ph[cut_test])
        std_model_s_ff = np.nanmean(std_model_seasonal_ff[cut_test])
        std_model_a = np.nanmean(std_model_ave[cut_test])

        alt_model_d_amp = np.nanmean(alt_model_daily_amp[cut_test])
        alt_model_d_ph = np.nanmean(alt_model_daily_ph[cut_test])
        alt_model_d_ff = np.nanmean(alt_model_daily_ff[cut_test])
        alt_model_s_amp = np.nanmean(alt_model_seasonal_amp[cut_test])
        alt_model_s_ph = np.nanmean(alt_model_seasonal_ph[cut_test])
        alt_model_s_ff = np.nanmean(alt_model_seasonal_ff[cut_test])
        alt_model_a = np.nanmean(alt_model_ave[cut_test])

        if plot_type == 'ave':
            ave_obs_param = obs_a
            ave_std_model_param = std_model_a
            ave_alt_model_param = alt_model_a

        if plot_type == 'd':
            if plot_type_2 == 'amp':
                ave_obs_param = obs_d_amp
                ave_std_model_param = std_model_d_amp
                ave_alt_model_param = alt_model_d_amp
            elif plot_type_2 == 'ph':
                ave_obs_param = obs_d_ph
                ave_std_model_param = std_model_d_ph
                ave_alt_model_param = alt_model_d_ph
            elif plot_type_2 == 'ff':
                ave_obs_param = obs_d_ff
                ave_std_model_param = std_model_d_ff
                ave_alt_model_param = alt_model_d_ff
                
            ave_obs_cs = np.average([obs_d_amp,max_phmax_seasonal_diff,max_ff_seasonal_diff],axis=0)
                
        if plot_type == 's':
            if plot_type_2 == 'amp':
                ave_obs_param = obs_s_amp
                ave_std_model_param = std_model_s_amp
                ave_alt_model_param = alt_model_s_amp
            elif plot_type_2 == 'ph':
                ave_obs_param = obs_s_ph
                ave_std_model_param = std_model_s_ph
                ave_alt_model_param = alt_model_s_ph
            elif plot_type_2 == 'ff':
                ave_obs_param = obs_s_ff
                ave_std_model_param = std_model_s_ff
                ave_alt_model_param = alt_model_s_ff
                
        if m == 0:
            diff_array[a,0] = ave_std_model_param - ave_obs_param
        diff_array[a,m+1] = ave_alt_model_param - ave_obs_param 
        
        count+=1



count = 0
for ax in axes.flat:

    try:
        area = areas[count]
    except:
        ax.axis('off')
        continue
        
    area_label = area_labels[area]
    
    model_types = ['STD','NOX0.25 NMVOC1','NOX0.5 NMVOC1','NOX2.0 NMVOC1','NOX4.0 NMVOC1','NOX1.0 ANMVOC0.25','NOX1.0 ANMVOC4.0','NOX0.25 ANMVOC0.25','NOX4.0 ANMVOC4.0','NOX0.25 ANMVOC4.0','NOX4.0 ANMVOC0.25','NOX1.0 BNMVOC0.25','NOX1.0 BNMVOC4.0','NOX0.25 BNMVOC0.25','NOX4.0 BNMVOC4.0','NOX0.25 BNMVOC4.0','NOX4.0 BNMVOC0.25']
    
    diff_cut = list(diff_array[count,:])  
    diffs = [abs(np.diff((0,i))[0]) for i in diff_cut]
    joined_list_1 = zip(diffs,diff_cut)    
    joined_list_2 = zip(diffs,mix_color_array)
    joined_list_3 = zip(diffs,model_types)
    joined_list_1.sort()
    joined_list_2.sort()
    joined_list_3.sort()
    model_diff_sorted = [y for x, y in joined_list_1]
    color_sorted = [y for x, y in joined_list_2]  
    model_type_sorted = [y for x, y in joined_list_3]
    
    width = 1
    ind = np.arange(len(model_diff_sorted))
    for t in range(len(ind)):
        ax.bar(ind[t]+width, model_diff_sorted[t], width, color=color_sorted[t])
    
    for tic in ax.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
        tic.label1On = tic.label2On = False
            
    ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)
    ax.text(0.1, -0.1, '1. %s'%(model_type_sorted[0]), ha='center', va='center',transform=ax.transAxes,fontsize=9)
    ax.text(0.75, -0.1, '2. %s'%(model_type_sorted[1]), ha='center', va='center',transform=ax.transAxes,fontsize=9)
    ax.text(0.1, -0.21, '3. %s'%(model_type_sorted[2]), ha='center', va='center',transform=ax.transAxes,fontsize=9)
    ax.text(0.75, -0.21, '4. %s'%(model_type_sorted[3]), ha='center', va='center',transform=ax.transAxes,fontsize=9)
    
    count+=1

plt.tight_layout(pad = 3.08)

h2, = ax.plot([1,1],color=mix_color_array[0],marker='o',linestyle='None',markersize=10)
h3, = ax.plot([1,1],color=mix_color_array[1],marker='o',linestyle='None',markersize=10)
h4, = ax.plot([1,1],color=mix_color_array[2],marker='o',linestyle='None',markersize=10)
h5, = ax.plot([1,1],color=mix_color_array[3],marker='o',linestyle='None',markersize=10)
h6, = ax.plot([1,1],color=mix_color_array[4],marker='o',linestyle='None',markersize=10)
h7, = ax.plot([1,1],color=mix_color_array[5],marker='o',linestyle='None',markersize=10)
h8, = ax.plot([1,1],color=mix_color_array[6],marker='o',linestyle='None',markersize=10)
h9, = ax.plot([1,1],color=mix_color_array[7],marker='o',linestyle='None',markersize=10)
h10, = ax.plot([1,1],color=mix_color_array[8],marker='o',linestyle='None',markersize=10)
h11, = ax.plot([1,1],color=mix_color_array[9],marker='o',linestyle='None',markersize=10)
h12, = ax.plot([1,1],color=mix_color_array[10],marker='o',linestyle='None',markersize=10)
h13, = ax.plot([1,1],color=mix_color_array[11],marker='o',linestyle='None',markersize=10)
h14, = ax.plot([1,1],color=mix_color_array[12],marker='o',linestyle='None',markersize=10)
h15, = ax.plot([1,1],color=mix_color_array[13],marker='o',linestyle='None',markersize=10)
h16, = ax.plot([1,1],color=mix_color_array[14],marker='o',linestyle='None',markersize=10)
h17, = ax.plot([1,1],color=mix_color_array[15],marker='o',linestyle='None',markersize=10)
h18, = ax.plot([1,1],color=mix_color_array[16],marker='o',linestyle='None',markersize=10)


plt.legend((h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18),['STD','NOX0.25 NMVOC1','NOX0.5 NMVOC1','NOX2.0 NMVOC1','NOX4.0 NMVOC1','NOX1.0 ANMVOC0.25','NOX1.0 ANMVOC4.0','NOX0.25 ANMVOC0.25','NOX4.0 ANMVOC4.0','NOX0.25 ANMVOC4.0','NOX4.0 ANMVOC0.25','NOX1.0 BNMVOC0.25','NOX1.0 BNMVOC4.0','NOX0.25 BNMVOC0.25','NOX4.0 BNMVOC4.0','NOX0.25 BNMVOC4.0','NOX4.0 BNMVOC0.25'],loc='lower left',prop={'size':7},fancybox=True,ncol=2,markerscale=1,bbox_to_anchor=(-0.2,0))
h2.set_visible(False)
h3.set_visible(False)
h4.set_visible(False)
h5.set_visible(False)
h6.set_visible(False)
h7.set_visible(False)
h8.set_visible(False)
h9.set_visible(False)
h10.set_visible(False)
h11.set_visible(False)
h12.set_visible(False)
h13.set_visible(False)
h14.set_visible(False)
h15.set_visible(False)
h16.set_visible(False)
h17.set_visible(False)
h18.set_visible(False)

plt.show()
