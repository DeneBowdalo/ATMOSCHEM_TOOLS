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

start_year = 2005
end_year = 2010

obs_fname = '/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_%s_2005_2010_H_PERIODIC.nc'%(species)
model_fname = '/work/home/db876/plotting_tools/model_files/MIROCCHEM_SURFACE_2000_2012_*_*_*_H_*.nc'

obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

#--------------------------------------------------------
#load in periodic lsp data
nmodels = 4
nmodels_linspace = np.linspace(0, 1, nmodels)

obs_period_grp = Dataset('../obs_SURFACE_H/obs_sig_periods.nc')
geoschemv90103_period_grp = Dataset('../../../2005_2012/2005_2010/GEOSCHEM_SURFACE_v90103_2x2.5_GEOS5_H_*/model_sig_periods.nc')
geoschemv901034x5_period_grp = Dataset('../../../2005_2012/2005_2010/GEOSCHEM_SURFACE_v90103_4x5_GEOS5_H_*/model_sig_periods.nc')
geoschemv902_period_grp = Dataset('../../../2005_2012/2005_2010/GEOSCHEM_SURFACE_v902_2x2.5_GEOS5_H_*/model_sig_periods.nc')
geoschemv1001_period_grp = Dataset('../../../2005_2012/2005_2010/GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_*/model_sig_periods.nc')

obs_daily_waveforms = []
obs_seasonal_waveforms = []
geoschemv90103_daily_waveforms = []
geoschemv90103_seasonal_waveforms = []
geoschemv901034x5_daily_waveforms = []
geoschemv901034x5_seasonal_waveforms = []
geoschemv902_daily_waveforms = []
geoschemv902_seasonal_waveforms = []
geoschemv1001_daily_waveforms = []
geoschemv1001_seasonal_waveforms = []

for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    site_group_geoschemv90103 = geoschemv90103_period_grp.groups[ref]
    site_group_geoschemv901034x5 = geoschemv901034x5_period_grp.groups[ref]
    site_group_geoschemv902 = geoschemv902_period_grp.groups[ref]
    site_group_geoschemv1001 = geoschemv1001_period_grp.groups[ref]

    geoschemv90103_daily_waveforms.append(site_group_geoschemv90103.variables['daily_waveform'][:])                                                                                        
    geoschemv90103_seasonal_waveforms.append(site_group_geoschemv90103.variables['seasonal_waveform'][:])

    geoschemv901034x5_daily_waveforms.append(site_group_geoschemv901034x5.variables['daily_waveform'][:])                                                                               
    geoschemv901034x5_seasonal_waveforms.append(site_group_geoschemv901034x5.variables['seasonal_waveform'][:])

    geoschemv902_daily_waveforms.append(site_group_geoschemv902.variables['daily_waveform'][:])                                                                                
    geoschemv902_seasonal_waveforms.append(site_group_geoschemv902.variables['seasonal_waveform'][:])

    geoschemv1001_daily_waveforms.append(site_group_geoschemv1001.variables['daily_waveform'][:])                                                                             
    geoschemv1001_seasonal_waveforms.append(site_group_geoschemv1001.variables['seasonal_waveform'][:])

    obs_daily_waveforms.append(site_group_obs.variables['daily_waveform'][:])
    obs_seasonal_waveforms.append(site_group_obs.variables['seasonal_waveform'][:])

obs_daily_waveforms = np.array(obs_daily_waveforms)
obs_seasonal_waveforms = np.array(obs_seasonal_waveforms)
geoschemv90103_daily_waveforms = np.array(geoschemv90103_daily_waveforms)                                                                                                                  
geoschemv90103_seasonal_waveforms = np.array(geoschemv90103_seasonal_waveforms)                                                                                                            
geoschemv901034x5_daily_waveforms = np.array(geoschemv901034x5_daily_waveforms)                                                                                                         
geoschemv901034x5_seasonal_waveforms = np.array(geoschemv901034x5_seasonal_waveforms)
geoschemv902_daily_waveforms = np.array(geoschemv902_daily_waveforms)                                                                                                          
geoschemv902_seasonal_waveforms = np.array(geoschemv902_seasonal_waveforms)                                                                                                    
geoschemv1001_daily_waveforms = np.array(geoschemv1001_daily_waveforms)                                                                                                       
geoschemv1001_seasonal_waveforms = np.array(geoschemv1001_seasonal_waveforms)

test = obs_daily_waveforms < 0
obs_daily_waveforms[test] = np.nan
test = obs_seasonal_waveforms < 0
obs_seasonal_waveforms[test] = np.nan

#-----------------------------------
#get area
areas = ['ANT','S_O','OC','AF','SA','NE_NA','CE_NA','SE_NA','S_NA','SW_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','NE_AS','SE_AS','C_AS','S_AS','N_O','AL','ARC']

plot_type = raw_input('\nd or s?\n')

if plot_type == 'd':
    obs_datetimes = obs_datetime_time[:24]
    model_datetimes = model_datetime_time[:24]
if plot_type == 's':
    obs_datetimes = obs_datetime_time[:8766]
    model_datetimes = model_datetime_time[:8766]

obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')

area_boundaries,area_tags,area_labels = modules.area_dicts()

fig, axes = plt.subplots(nrows=5, ncols=5,figsize=(19,13))
fig.patch.set_facecolor('white')
count = 0
for ax in axes.flat:
    try:
        area = areas[count]
    except:
        ax.axis('off')
        continue
    
    print area

    area_grid = area_boundaries[area]
    area_tag = area_tags[area]
    area_label = area_labels[area]

    cut_test = modules.area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag)

    obs_d_w = np.nanmean(obs_daily_waveforms[cut_test,:],axis=0)
    obs_s_w = np.nanmean(obs_seasonal_waveforms[cut_test,:],axis=0)

    geoschemv90103_d_w = np.average(geoschemv90103_daily_waveforms[cut_test,:],axis=0)                                                                                                     
    geoschemv90103_s_w = np.average(geoschemv90103_seasonal_waveforms[cut_test,:],axis=0)

    geoschemv901034x5_d_w = np.average(geoschemv901034x5_daily_waveforms[cut_test,:],axis=0)
    geoschemv901034x5_s_w = np.average(geoschemv901034x5_seasonal_waveforms[cut_test,:],axis=0)

    geoschemv902_d_w = np.average(geoschemv902_daily_waveforms[cut_test,:],axis=0)                                                                                             
    geoschemv902_s_w = np.average(geoschemv902_seasonal_waveforms[cut_test,:],axis=0)

    geoschemv1001_d_w = np.average(geoschemv1001_daily_waveforms[cut_test,:],axis=0)                                                                                          
    geoschemv1001_s_w = np.average(geoschemv1001_seasonal_waveforms[cut_test,:],axis=0)

    if plot_type == 'd':
        ave_obs_waveform = obs_d_w
        ave_geoschemv90103_waveform = geoschemv90103_d_w 
        ave_geoschemv901034x5_waveform = geoschemv901034x5_d_w 
        ave_geoschemv902_waveform = geoschemv902_d_w
        ave_geoschemv1001_waveform = geoschemv1001_d_w

    if plot_type == 's':
        ave_obs_waveform = obs_s_w
        ave_geoschemv90103_waveform = geoschemv90103_s_w 
        ave_geoschemv901034x5_waveform = geoschemv901034x5_s_w 
        ave_geoschemv902_waveform = geoschemv902_s_w 
        ave_geoschemv1001_waveform = geoschemv1001_s_w

    ax.plot_date(obs_time_pd.to_pydatetime(),ave_obs_waveform,color = 'black',linestyle='--',linewidth=1,markersize=5,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_geoschemv90103_waveform,color = plt.cm.jet(nmodels_linspace[0]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_geoschemv901034x5_waveform,color = plt.cm.jet(nmodels_linspace[1]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_geoschemv902_waveform,color = plt.cm.jet(nmodels_linspace[2]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_geoschemv1001_waveform,color = plt.cm.jet(nmodels_linspace[3]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')

    ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)
    if plot_type == 'd':
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
    if plot_type == 's':
        ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))

    tst = False
    for t in range(len(ave_obs_waveform)):
        b = np.isnan(ave_obs_waveform[t])
        if b == False:
            tst= True
    if tst == False:
        ax.plot_date(obs_time_pd.to_pydatetime(),len(obs_time_pd)*[1],markersize=0.000001)
    
    count+=1
    
plt.tight_layout(pad = 3.08)

h1, = ax.plot([1,1],color='black',marker='o',linestyle='None',markersize=10)
h2, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[0]),marker='o',linestyle='None',markersize=10)
h3, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[1]),marker='o',linestyle='None',markersize=10)
h4, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[2]),marker='o',linestyle='None',markersize=10)
h5, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[3]),marker='o',linestyle='None',markersize=10)

plt.legend((h1,h2,h3,h4,h5),['Observations','GEOSChemv90103 2x2.5','GEOSChemv90103 4x5','GEOSChemv902 2x2.5','GEOSChemv1001 4x5'],loc='lower left',prop={'size':10},fancybox=True,ncol=1,markerscale=1,bbox_to_anchor=(-0.2,0))
h1.set_visible(False)
h2.set_visible(False)
h3.set_visible(False)
h4.set_visible(False)
h5.set_visible(False)

plt.show()
