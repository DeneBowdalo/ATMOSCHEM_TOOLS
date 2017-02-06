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

model_fname = '/work/home/db876/plotting_tools/model_files/SOCOL3_SURFACE_O3_2005_2010_*_*_*_D_*.nc'
model_ts_grp = Dataset(model_fname)
model_var_mirocchem = model_ts_grp.variables['o3'][:]
model_date = model_ts_grp.variables['date'][:]
model_time = model_ts_grp.variables['time'][:]

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
obs_fname = '/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_O3_2005_2010_D.nc'
obs_ts_grp = Dataset(obs_fname)
obs_refs_dict = obs_ts_grp.groups

obs_refs = []
obs_lats = []
obs_lons = []
obs_alt = []
obs_country = []

for i in obs_refs_dict.keys():
    i = i.encode('ascii')
    obs_refs = np.append(obs_refs,i)

for ref in obs_refs:
    obs_site_group = obs_ts_grp.groups[ref] 
    obs_country = np.append(obs_country,obs_site_group.country)
    obs_lats = np.append(obs_lats,obs_site_group.latitude)
    obs_lons = np.append(obs_lons,obs_site_group.longitude)
    obs_alt = np.append(obs_alt,obs_site_group.altitude)
    obs_date = obs_site_group.variables['date'][:]
    obs_time = obs_site_group.variables['time'][:]
    
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

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

#--------------------------------------------------------
#load in periodic lsp data
nmodels = 10
nmodels_linspace = np.linspace(0, 1, nmodels)

obs_period_grp = Dataset('../obs_SURFACE_D/obs_sig_periods.nc')
#hadgem3es_period_grp = Dataset('../HADGEM3ES_SURFACE_*_*_*_D_*/model_sig_periods.nc')
miroc3_period_grp = Dataset('../MIROC3_SURFACE_*_*_*_D_*/model_sig_periods.nc')
ipsl_period_grp = Dataset('../IPSL_SURFACE_*_*_*_D_*/model_sig_periods.nc')
geoschemv90103_period_grp = Dataset('../GEOSCHEM_SURFACE_v90103_2x2.5_GEOS5_D_*/model_sig_periods.nc')
geoschemv902_period_grp = Dataset('../GEOSCHEM_SURFACE_v902_2x2.5_GEOS5_D_*/model_sig_periods.nc')
socol3_period_grp = Dataset('../SOCOL3_SURFACE_*_*_*_D_*/model_sig_periods.nc')
ukca_period_grp = Dataset('../UKCA_SURFACE_*_*_*_D_*/model_sig_periods.nc')
miroc3sd_period_grp = Dataset('../MIROC3_SURFACE_*_*_*_D_SD/model_sig_periods.nc')
ipslsd_period_grp = Dataset('../IPSL_SURFACE_*_*_*_D_SD/model_sig_periods.nc')
mrisd_period_grp = Dataset('../MRI_SURFACE_*_*_*_D_SD/model_sig_periods.nc')


obs_seasonal_waveforms = []
obs_full_waveforms = []
hadgem3es_seasonal_waveforms = []
hadgem3es_full_waveforms = []
miroc3_seasonal_waveforms = []
miroc3_full_waveforms = []
ipsl_seasonal_waveforms = []
ipsl_full_waveforms = []
geoschemv90103_seasonal_waveforms = []
geoschemv90103_full_waveforms = [] 
geoschemv902_seasonal_waveforms = []
geoschemv902_full_waveforms = [] 
socol3_seasonal_waveforms = []
socol3_full_waveforms = [] 
ukca_seasonal_waveforms = []
ukca_full_waveforms = []
miroc3sd_seasonal_waveforms = []
miroc3sd_full_waveforms = []
ipslsd_seasonal_waveforms = []
ipslsd_full_waveforms = []
mrisd_seasonal_waveforms = []
mrisd_full_waveforms = [] 

for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    #site_group_hadgem3es = hadgem3es_period_grp.groups[ref]
    site_group_miroc3 = miroc3_period_grp.groups[ref]    
    site_group_ipsl = ipsl_period_grp.groups[ref]
    site_group_geoschemv90103 = geoschemv90103_period_grp.groups[ref]
    site_group_geoschemv902 = geoschemv902_period_grp.groups[ref]
    site_group_socol3 = socol3_period_grp.groups[ref]
    site_group_ukca = ukca_period_grp.groups[ref]
    site_group_miroc3sd = miroc3sd_period_grp.groups[ref]
    site_group_ipslsd = ipslsd_period_grp.groups[ref]
    site_group_mrisd = mrisd_period_grp.groups[ref]

    #hadgem3es_seasonal_waveforms.append(site_group_hadgem3es.variables['seasonal_waveform'][:])
    #hadgem3es_full_waveforms.append(site_group_hadgem3es.variables['all_waveform'][:])

    miroc3_seasonal_waveforms.append(site_group_miroc3.variables['seasonal_waveform'][:])
    miroc3_full_waveforms.append(site_group_miroc3.variables['all_waveform'][:]) 
    
    ipsl_seasonal_waveforms.append(site_group_ipsl.variables['seasonal_waveform'][:])
    ipsl_full_waveforms.append(site_group_ipsl.variables['all_waveform'][:]) 

    geoschemv90103_seasonal_waveforms.append(site_group_geoschemv90103.variables['seasonal_waveform'][:])
    geoschemv90103_full_waveforms.append(site_group_geoschemv90103.variables['all_waveform'][:]) 

    geoschemv902_seasonal_waveforms.append(site_group_geoschemv902.variables['seasonal_waveform'][:])
    geoschemv902_full_waveforms.append(site_group_geoschemv902.variables['all_waveform'][:]) 

    socol3_seasonal_waveforms.append(site_group_socol3.variables['seasonal_waveform'][:])
    socol3_full_waveforms.append(site_group_socol3.variables['all_waveform'][:])

    ukca_seasonal_waveforms.append(site_group_ukca.variables['seasonal_waveform'][:])
    ukca_full_waveforms.append(site_group_ukca.variables['all_waveform'][:])

    miroc3sd_seasonal_waveforms.append(site_group_miroc3sd.variables['seasonal_waveform'][:])
    miroc3sd_full_waveforms.append(site_group_miroc3sd.variables['all_waveform'][:])

    ipslsd_seasonal_waveforms.append(site_group_ipslsd.variables['seasonal_waveform'][:])
    ipslsd_full_waveforms.append(site_group_ipslsd.variables['all_waveform'][:])

    mrisd_seasonal_waveforms.append(site_group_mrisd.variables['seasonal_waveform'][:])
    mrisd_full_waveforms.append(site_group_mrisd.variables['all_waveform'][:])

    obs_seasonal_waveforms.append(site_group_obs.variables['seasonal_waveform'][:])
    obs_full_waveforms.append(site_group_obs.variables['all_waveform'][:])


obs_seasonal_waveforms = np.array(obs_seasonal_waveforms)
obs_full_waveforms = np.array(obs_full_waveforms)	
#hadgem3es_seasonal_waveforms = np.array(hadgem3es_seasonal_waveforms)
#hadgem3es_full_waveforms = np.array(hadgem3es_full_waveforms)
miroc3_seasonal_waveforms = np.array(miroc3_seasonal_waveforms)
miroc3_full_waveforms = np.array(miroc3_full_waveforms)
ipsl_seasonal_waveforms = np.array(ipsl_seasonal_waveforms)
ipsl_full_waveforms = np.array(ipsl_full_waveforms)
geoschemv90103_seasonal_waveforms = np.array(geoschemv90103_seasonal_waveforms)                                                                                                            
geoschemv90103_full_waveforms = np.array(geoschemv90103_full_waveforms)
geoschemv902_seasonal_waveforms = np.array(geoschemv902_seasonal_waveforms)                                                                                                    
geoschemv902_full_waveforms = np.array(geoschemv902_full_waveforms) 
socol3_seasonal_waveforms = np.array(socol3_seasonal_waveforms)                                                                                                          
socol3_full_waveforms = np.array(socol3_full_waveforms)
ukca_seasonal_waveforms = np.array(ukca_seasonal_waveforms)
ukca_full_waveforms = np.array(ukca_full_waveforms)
miroc3sd_seasonal_waveforms = np.array(miroc3sd_seasonal_waveforms)                                                                                                          
miroc3sd_full_waveforms = np.array(miroc3sd_full_waveforms)
ipslsd_seasonal_waveforms = np.array(ipslsd_seasonal_waveforms)                                                                                                          
ipslsd_full_waveforms = np.array(ipslsd_full_waveforms)
mrisd_seasonal_waveforms = np.array(mrisd_seasonal_waveforms)                                                                                                          
mrisd_full_waveforms = np.array(mrisd_full_waveforms)

#test = obs_daily_waveforms < 0
#obs_daily_waveforms[test] = np.nan
#test = obs_seasonal_waveforms < 0
#obs_seasonal_waveforms[test] = np.nan
#test = obs_full_waveforms < 0
#obs_full_waveforms[test] = np.nan

#-----------------------------------
#get area
areas = ['ANT','OC','S_O','AF','SE_US','S_US','W_US','N_US','NE_US','W_CAN','E_CAN','S_EU','C_EU','NW_EU','N_EU','E_EU','AS','N_O','ARC']

plot_type = 's'
#plot_type = raw_input('\nd, s or full?\n')

obs_datetimes = obs_datetimes[:365]
model_datetimes = model_datetimes[:365]

obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'D')
model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'D')

area_boundaries = {'NE_US':[37,50,-90,-60],'SE_US':[25,37,-90,-60],'S_US':[25,37,-105,-90],'W_US':[25,50,-130,-105],'N_US':[37,60,-105,-90],'N_EU':[57,80,5,45],'C_EU':[47,57,5,20],'S_EU':[30,47,-10,40],'E_EU':[47,57,20,40],'NW_EU':[47,70,-15,5],'AS':[0,0,0,0],'N_O':[0,0,0,0],'S_O':[0,0,0,0],'OC':[0,0,0,0],'AF':[0,0,0,0],'E_CAN':[40,80,-95,-50],'W_CAN':[40,80,-150,-95],'ANT':[0,0,0,0],'ARC':[0,0,0,0]}
area_countries = {'NE_US':'United States','SE_US':'United States','S_US':'United States','W_US':'United States','N_US':'United States','C_EU':'EU','N_EU':'EU','S_EU':'EU','E_EU':'EU','NW_EU':'EU','AS':'AS','N_O':'O','S_O':'O','OC':'OC','AF':'AF','E_CAN':'Canada','W_CAN':'Canada','ANT':'ANT','ARC':'ARC'}
area_labels = {'NE_US':'NE US','SE_US':'SE US','S_US':'S US','W_US':'W US','N_US':'N US','N_EU':'N EU','NW_EU':'NW EU','C_EU':'C EU','E_EU':'E EU','S_EU':'S EU','AS':'Asia','N_O':'NH Oceanic','S_O':'SH Oceanic','OC':'Oceania','AF':'Africa','E_CAN':'East Canada','W_CAN':'West Canada','ANT':'Antarctica','ARC':'Arctic'}


fig, axes = plt.subplots(nrows=4, ncols=5,figsize=(19,13))
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
    area_country = area_countries[area]
    area_label = area_labels[area]

    if (area == 'C_EU') or (area == 'E_EU') or (area == 'S_EU') or (area == 'NW_EU') or (area == 'N_EU'):
        cut_test = (obs_lats >= area_grid[0]) & (obs_lats < area_grid[1]) & (obs_lons >= area_grid[2]) & (obs_lons < area_grid[3]) & (tags == area_country)
    elif (area == 'AS') or (area == 'OC') or (area == 'AF') or (area == 'ANT') or (area == 'ARC'):
        cut_test = tags == area_country
    elif (area == 'N_O'):
        cut_test = (tags == area_country) & (obs_lats >= 0)
    elif (area == 'S_O'):
        cut_test = (tags == area_country) & (obs_lats < 0)
    else:
        cut_test = (obs_lats >= area_grid[0]) & (obs_lats < area_grid[1]) & (obs_lons >= area_grid[2]) & (obs_lons < area_grid[3]) & (obs_country == area_country)

    obs_s_w = np.nanmean(obs_seasonal_waveforms[cut_test,:],axis=0)
    obs_full_w = np.nanmean(obs_full_waveforms[cut_test,:],axis=0)

    #hadgem3es_s_w = np.average(hadgem3es_seasonal_waveforms[cut_test,:],axis=0)
    #hadgem3es_full_w = np.average(hadgem3es_full_waveforms[cut_test,:],axis=0)

    miroc3_s_w = np.average(miroc3_seasonal_waveforms[cut_test,:],axis=0)
    miroc3_full_w = np.average(miroc3_full_waveforms[cut_test,:],axis=0)

    ipsl_s_w = np.average(ipsl_seasonal_waveforms[cut_test,:],axis=0)
    ipsl_full_w = np.average(ipsl_full_waveforms[cut_test,:],axis=0)

    geoschemv90103_s_w = np.average(geoschemv90103_seasonal_waveforms[cut_test,:],axis=0)
    geoschemv90103_full_w = np.average(geoschemv90103_full_waveforms[cut_test,:],axis=0)

    geoschemv902_s_w = np.average(geoschemv902_seasonal_waveforms[cut_test,:],axis=0)
    geoschemv902_full_w = np.average(geoschemv902_full_waveforms[cut_test,:],axis=0) 

    socol3_s_w = np.average(socol3_seasonal_waveforms[cut_test,:],axis=0)
    socol3_full_w = np.average(socol3_full_waveforms[cut_test,:],axis=0)

    ukca_s_w = np.average(ukca_seasonal_waveforms[cut_test,:],axis=0)
    ukca_full_w = np.average(ukca_full_waveforms[cut_test,:],axis=0)

    miroc3sd_s_w = np.average(miroc3sd_seasonal_waveforms[cut_test,:],axis=0)
    miroc3sd_full_w = np.average(miroc3sd_full_waveforms[cut_test,:],axis=0)

    ipslsd_s_w = np.average(ipslsd_seasonal_waveforms[cut_test,:],axis=0)
    ipslsd_full_w = np.average(ipslsd_full_waveforms[cut_test,:],axis=0) 

    mrisd_s_w = np.average(mrisd_seasonal_waveforms[cut_test,:],axis=0)                                                                                                    
    mrisd_full_w = np.average(mrisd_full_waveforms[cut_test,:],axis=0) 

    if plot_type == 's':
        ave_obs_waveform = obs_s_w
        #ave_hadgem3es_waveform = hadgem3es_s_w
        ave_miroc3_waveform = miroc3_s_w
        ave_ipsl_waveform = ipsl_s_w    
        ave_geoschemv90103_waveform = geoschemv90103_s_w 
        ave_geoschemv902_waveform = geoschemv902_s_w 
        ave_socol3_waveform = socol3_s_w 
        ave_ukca_waveform = ukca_s_w 
        ave_miroc3sd_waveform = miroc3sd_s_w
        ave_ipslsd_waveform = ipslsd_s_w
        ave_mrisd_waveform = mrisd_s_w

    if plot_type == 'full':
        ave_obs_waveform = obs_full_w
        #ave_hadgem3es_waveform = hadgem3es_full_w
        ave_miroc3_waveform = miroc3_full_w
        ave_ipsl_waveform = ipsl_full_w    
        ave_geoschemv90103_waveform = geoschemv90103_full_w 
        ave_geoschemv902_waveform = geoschemv902_full_w 
        ave_socol3_waveform = socol3_full_w 
        ave_ukca_waveform = ukca_full_w 
        ave_miroc3sd_waveform = miroc3sd_full_w
        ave_ipslsd_waveform = ipslsd_full_w
        ave_mrisd_waveform = mrisd_full_w

    ax.plot_date(obs_time_pd.to_pydatetime(),ave_obs_waveform,color = 'black',linestyle='_',linewidth=3,markeredgecolor='None')
    #ax.plot_date(model_time_pd.to_pydatetime(),ave_hadgem3es_waveform,color = plt.cm.jet(nmodels_linspace[0]),linestyle='_',linewidth=1,marker='o',markersize=1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_miroc3_waveform,color = plt.cm.jet(nmodels_linspace[1]),linestyle='_',linewidth=1,marker='o',markersize=1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_ipsl_waveform,color = plt.cm.jet(nmodels_linspace[2]),linestyle='_',linewidth=1,marker='o',markersize=1,markeredgecolor='None') 
    ax.plot_date(model_time_pd.to_pydatetime(),ave_geoschemv90103_waveform,color = plt.cm.jet(nmodels_linspace[3]),linestyle='_',linewidth=1,marker='o',markersize=1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_geoschemv902_waveform,color = plt.cm.jet(nmodels_linspace[4]),linestyle='_',linewidth=1,marker='o',markersize=1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_socol3_waveform,color = plt.cm.jet(nmodels_linspace[5]),linestyle='_',linewidth=1,marker='o',markersize=1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_ukca_waveform,color = plt.cm.jet(nmodels_linspace[6]),linestyle='_',linewidth=1,marker='o',markersize=1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_miroc3sd_waveform,color = plt.cm.jet(nmodels_linspace[7]),linestyle='_',linewidth=1,marker='o',markersize=1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_ipslsd_waveform,color = plt.cm.jet(nmodels_linspace[8]),linestyle='_',linewidth=1,marker='o',markersize=1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_mrisd_waveform,color = plt.cm.jet(nmodels_linspace[9]),linestyle='_',linewidth=1,marker='o',markersize=1,markeredgecolor='None')

    ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)
    if plot_type == 's':
        ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
    
    count+=1
    
plt.tight_layout(pad = 3.08)

h1, = ax.plot([1,1],color='black',marker='o',linestyle='None',markersize=10)
#h2, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[0]),marker='o',linestyle='None',markersize=10)
h3, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[1]),marker='o',linestyle='None',markersize=10)
h4, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[2]),marker='o',linestyle='None',markersize=10)
h5, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[3]),marker='o',linestyle='None',markersize=10)
h6, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[4]),marker='o',linestyle='None',markersize=10)
h7, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[5]),marker='o',linestyle='None',markersize=10)
h8, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[6]),marker='o',linestyle='None',markersize=10)
h9, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[7]),marker='o',linestyle='None',markersize=10)
h10, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[8]),marker='o',linestyle='None',markersize=10)
h11, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[9]),marker='o',linestyle='None',markersize=10)

plt.legend((h1,h3,h4,h5,h6,h7,h8,h9,h10,h11),['Observations','MIROC3','IPSL','GEOSChemv90103','GEOSChemv902','SOCOL3','UKCA','MIROC3-SD','IPSL-SD','MRI-SD'],loc='lower left',prop={'size':11},fancybox=True,ncol=1,markerscale=1,bbox_to_anchor=(-0.1,-0.2))
h1.set_visible(False)
#h2.set_visible(False)
h3.set_visible(False)
h4.set_visible(False)
h5.set_visible(False)
h6.set_visible(False)
h7.set_visible(False)
h8.set_visible(False)
h9.set_visible(False)
h10.set_visible(False)
h11.set_visible(False)

plt.show()
