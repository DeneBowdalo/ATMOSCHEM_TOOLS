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
species = paths[-3]

model_fname = '/work/home/db876/plotting_tools/model_files/MIROCCHEM_SURFACE_%s_2005_2010_*_*_*_M_*.nc'%(species)
model_ts_grp = Dataset(model_fname)
model_var_mirocchem = model_ts_grp.variables[species.lower()][:]
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
obs_fname = '/work/home/db876/observations/surface/%s/process/GLOBAL_SURFACE_%s_2005_2010_M_PERIODIC.nc'%(species,species)
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
nmodels = 7
nmodels_linspace = np.linspace(0, 1, nmodels)

obs_period_grp = Dataset('../obs_SURFACE_M/obs_sig_periods.nc')
cesmcam_period_grp = Dataset('../CESMCAM_SURFACE_*_*_*_M_*/model_sig_periods.nc')
cmam_period_grp = Dataset('../CMAM_SURFACE_*_*_*_M_*/model_sig_periods.nc')
geosccm_period_grp = Dataset('../GEOSCCM_SURFACE_*_*_*_M_*/model_sig_periods.nc')
#geoschemv90103_period_grp = Dataset('../GEOSCHEM_SURFACE_v90103_2x2.5_GEOS5_M_*/model_sig_periods.nc')
geoschemv902_period_grp = Dataset('../GEOSCHEM_SURFACE_v902_2x2.5_GEOS5_M_*/model_sig_periods.nc')
#geoschemv1001_period_grp = Dataset('../GEOSCHEM_SURFACE_v1001_2x2.5_GEOS5_M_*/model_sig_periods.nc')
gfdlam3_period_grp = Dataset('../GFDLAM3_SURFACE_*_*_*_M_*/model_sig_periods.nc')
gisse2r_period_grp = Dataset('../GISSE2R_SURFACE_*_*_*_M_*/model_sig_periods.nc')
mirocchem_period_grp = Dataset('../MIROCCHEM_SURFACE_*_*_*_M_*/model_sig_periods.nc') 
#mrisd_period_grp = Dataset('../MRI_SURFACE_*_*_*_D_SD/model_sig_periods.nc')
#socol3_period_grp = Dataset('../SOCOL3_SURFACE_*_*_*_D_*/model_sig_periods.nc') 
#ukca_period_grp = Dataset('../UKCA_SURFACE_*_*_*_D_*/model_sig_periods.nc')

obs_seasonal_waveforms = []
cesmcam_seasonal_waveforms = []
cmam_seasonal_waveforms = []
geosccm_seasonal_waveforms = []
geoschemv90103_seasonal_waveforms = []
geoschemv902_seasonal_waveforms = []
geoschemv1001_seasonal_waveforms = []
gfdlam3_seasonal_waveforms = []
gisse2r_seasonal_waveforms = []
mirocchem_seasonal_waveforms = []
mrisd_seasonal_waveforms = []
socol3_seasonal_waveforms = []
ukca_seasonal_waveforms = []

for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    site_group_cesmcam = cesmcam_period_grp.groups[ref]
    site_group_cmam = cmam_period_grp.groups[ref]    
    site_group_geosccm = geosccm_period_grp.groups[ref]
    #site_group_geoschemv90103 = geoschemv90103_period_grp.groups[ref]
    site_group_geoschemv902 = geoschemv902_period_grp.groups[ref]
    #site_group_geoschemv1001 = geoschemv1001_period_grp.groups[ref]
    site_group_gfdlam3 = gfdlam3_period_grp.groups[ref]
    site_group_gisse2r = gisse2r_period_grp.groups[ref]
    site_group_mirocchem = mirocchem_period_grp.groups[ref]
    #site_group_mrisd = mrisd_period_grp.groups[ref]
    #site_group_socol3 = socol3_period_grp.groups[ref]     
    #site_group_ukca = ukca_period_grp.groups[ref]

    cesmcam_seasonal_waveforms.append(site_group_cesmcam.variables['seasonal_waveform'][:])

    cmam_seasonal_waveforms.append(site_group_cmam.variables['seasonal_waveform'][:])
    
    geosccm_seasonal_waveforms.append(site_group_geosccm.variables['seasonal_waveform'][:])

    #geoschemv90103_seasonal_waveforms.append(site_group_geoschemv90103.variables['seasonal_waveform'][:])

    geoschemv902_seasonal_waveforms.append(site_group_geoschemv902.variables['seasonal_waveform'][:])

    #geoschemv1001_seasonal_waveforms.append(site_group_geoschemv1001.variables['seasonal_waveform'][:])

    gfdlam3_seasonal_waveforms.append(site_group_gfdlam3.variables['seasonal_waveform'][:])

    gisse2r_seasonal_waveforms.append(site_group_gisse2r.variables['seasonal_waveform'][:])

    mirocchem_seasonal_waveforms.append(site_group_mirocchem.variables['seasonal_waveform'][:])

    #mrisd_seasonal_waveforms.append(site_group_mrisd.variables['seasonal_waveform'][:])

    #socol3_seasonal_waveforms.append(site_group_socol3.variables['seasonal_waveform'][:])

    #ukca_seasonal_waveforms.append(site_group_ukca.variables['seasonal_waveform'][:])

    obs_seasonal_waveforms.append(site_group_obs.variables['seasonal_waveform'][:])

obs_seasonal_waveforms = np.array(obs_seasonal_waveforms)
cesmcam_seasonal_waveforms = np.array(cesmcam_seasonal_waveforms)
cmam_seasonal_waveforms = np.array(cesmcam_seasonal_waveforms)
geosccm_seasonal_waveforms = np.array(geosccm_seasonal_waveforms)
#geoschemv90103_seasonal_waveforms = np.array(geoschemv90103_seasonal_waveforms)                                                                                                            
geoschemv902_seasonal_waveforms = np.array(geoschemv902_seasonal_waveforms)                                                                                                    
#geoschemv1001_seasonal_waveforms = np.array(geoschemv1001_seasonal_waveforms)
gfdlam3_seasonal_waveforms = np.array(gfdlam3_seasonal_waveforms)                                                                                                          
gisse2r_seasonal_waveforms = np.array(gisse2r_seasonal_waveforms)
mirocchem_seasonal_waveforms = np.array(mirocchem_seasonal_waveforms)                                                                                                          
#mrisd_seasonal_waveforms = np.array(mrisd_seasonal_waveforms) 
#socol3_seasonal_waveforms = np.array(socol3_seasonal_waveforms) 
#ukca_seasonal_waveforms = np.array(ukca_seasonal_waveforms)

test = obs_seasonal_waveforms < 0
obs_seasonal_waveforms[test] = np.nan

#-----------------------------------
#get area
areas = ['ANT','S_O','OC','AF','SA','NE_NA','CE_NA','SE_NA','S_NA','SW_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','NE_AS','SE_AS','C_AS','S_AS','N_O','AL','ARC']

plot_type = 's'

obs_time_pd = []
model_time_pd = []

for m in range(1,13):
    obs_time_pd.append(datetime.datetime(2005,m,1,0,0))
    model_time_pd.append(datetime.datetime(2005,m,1,0,0))

#obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'M')
#model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'M')

obs_time_pd = pd.DatetimeIndex(obs_time_pd)
model_time_pd = pd.DatetimeIndex(model_time_pd)

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

    obs_s_w = np.nanmean(obs_seasonal_waveforms[cut_test,:],axis=0)
    cesmcam_s_w = np.average(cesmcam_seasonal_waveforms[cut_test,:],axis=0)
    cmam_s_w = np.average(cmam_seasonal_waveforms[cut_test,:],axis=0)
    geosccm_s_w = np.average(geosccm_seasonal_waveforms[cut_test,:],axis=0)
    #geoschemv90103_s_w = np.average(geoschemv90103_seasonal_waveforms[cut_test,:],axis=0)
    geoschemv902_s_w = np.average(geoschemv902_seasonal_waveforms[cut_test,:],axis=0)
    #geoschemv1001_s_w = np.average(geoschemv1001_seasonal_waveforms[cut_test,:],axis=0)
    gfdlam3_s_w = np.average(gfdlam3_seasonal_waveforms[cut_test,:],axis=0)
    gisse2r_s_w = np.average(gisse2r_seasonal_waveforms[cut_test,:],axis=0)
    mirocchem_s_w = np.average(mirocchem_seasonal_waveforms[cut_test,:],axis=0)
    #mrisd_s_w = np.average(mrisd_seasonal_waveforms[cut_test,:],axis=0)
    #socol3_s_w = np.average(socol3_seasonal_waveforms[cut_test,:],axis=0) 
    #ukca_s_w = np.average(ukca_seasonal_waveforms[cut_test,:],axis=0) 

    if plot_type == 's':
        ave_obs_waveform = obs_s_w
        ave_cesmcam_waveform = cesmcam_s_w
        ave_cmam_waveform = cmam_s_w
        ave_geosccm_waveform = geosccm_s_w    
        #ave_geoschemv90103_waveform = geoschemv90103_s_w 
        ave_geoschemv902_waveform = geoschemv902_s_w 
        #ave_geoschemv1001_waveform = geoschemv1001_s_w
        ave_gfdlam3_waveform = gfdlam3_s_w 
        ave_gisse2r_waveform = gisse2r_s_w 
        ave_mirocchem_waveform = mirocchem_s_w
        #ave_mrisd_waveform = mrisd_s_w
        #ave_socol3_waveform = socol3_s_w
        #ave_ukca_waveform = ukca_s_w

    ax.plot_date(obs_time_pd.to_pydatetime(),ave_obs_waveform,color = 'black',linestyle='--',linewidth=1,marker='o',markersize = 5,markeredgecolor='None')#,linewidth=1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_cesmcam_waveform,color = plt.cm.jet(nmodels_linspace[0]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_cmam_waveform,color = plt.cm.jet(nmodels_linspace[1]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_geosccm_waveform,color = plt.cm.jet(nmodels_linspace[2]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None') 
    #ax.plot_date(model_time_pd.to_pydatetime(),ave_geoschemv90103_waveform,color = plt.cm.jet(nmodels_linspace[3]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_geoschemv902_waveform,color = plt.cm.jet(nmodels_linspace[3]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    #ax.plot_date(model_time_pd.to_pydatetime(),ave_geoschemv1001_waveform,color = plt.cm.jet(nmodels_linspace[5]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_gfdlam3_waveform,color = plt.cm.jet(nmodels_linspace[4]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_gisse2r_waveform,color = plt.cm.jet(nmodels_linspace[5]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_mirocchem_waveform,color = plt.cm.jet(nmodels_linspace[6]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    #ax.plot_date(model_time_pd.to_pydatetime(),ave_mrisd_waveform,color = plt.cm.jet(nmodels_linspace[9]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    #ax.plot_date(model_time_pd.to_pydatetime(),ave_socol3_waveform,color = plt.cm.jet(nmodels_linspace[10]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')
    #ax.plot_date(model_time_pd.to_pydatetime(),ave_ukca_waveform,color = plt.cm.jet(nmodels_linspace[11]),linestyle='--',linewidth=1,marker='o',markersize = 1,markeredgecolor='None')

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
h6, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[4]),marker='o',linestyle='None',markersize=10)
h7, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[5]),marker='o',linestyle='None',markersize=10)
h8, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[6]),marker='o',linestyle='None',markersize=10)
#h9, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[7]),marker='o',linestyle='None',markersize=10)
#h10, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[8]),marker='o',linestyle='None',markersize=10)
#h11, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[9]),marker='o',linestyle='None',markersize=10)
#h12, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[10]),marker='o',linestyle='None',markersize=10)
#h13, = ax.plot([1,1],color=plt.cm.jet(nmodels_linspace[11]),marker='o',linestyle='None',markersize=10)

plt.legend((h1,h2,h3,h4,h5,h6,h7,h8),['Observations','CESMCAM','CMAM','GEOSCCM','GEOSChemv902','GFDLAM3','GISSE2R','MIROCCHEM'],loc='lower left',prop={'size':10},fancybox=True,ncol=2,markerscale=1,bbox_to_anchor=(-0.2,0))
h1.set_visible(False)
h2.set_visible(False)
h3.set_visible(False)
h4.set_visible(False)
h5.set_visible(False)
h6.set_visible(False)
h7.set_visible(False)
h8.set_visible(False)
#h9.set_visible(False)
#h10.set_visible(False)
#h11.set_visible(False)
#h12.set_visible(False)
#h13.set_visible(False)

plt.show()
