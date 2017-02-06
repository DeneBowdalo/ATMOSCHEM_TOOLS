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
from bpch import bpch
from MChem_tools import *

species = 'O3'

param = raw_input('NOX, ISOP or METH?\n')

#-----------------------------
#read in model data
model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_O3_2005_2010_v90103_4x5_GEOS5_H_*.nc'
model_ts_grp = Dataset(model_fname)
model_var = model_ts_grp.variables[species.lower()][:]
model_date = model_ts_grp.variables['date'][:]
model_time = model_ts_grp.variables['time'][:]
lat_e = model_ts_grp.variables['lat_edges'][:]
lon_e = model_ts_grp.variables['lon_edges'][:]
lat_c = model_ts_grp.variables['lat_centre'][:]
lon_c = model_ts_grp.variables['lon_centre'][:]
grid_size = model_ts_grp.variables['grid_size'][:]
grid_size = grid_size[0]

#------------------------------------------------------------
#read in obs time series data
obs_fname = '/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_O3_2005_2010_M_PERIODIC.nc'

obs_ts_grp = Dataset(obs_fname)
obs_refs_dict = obs_ts_grp.groups

obs_data = []
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
    obs_data.append(obs_site_group.variables[species.lower()][:])
    obs_country = np.append(obs_country,obs_site_group.country)
    obs_lats = np.append(obs_lats,obs_site_group.latitude)
    obs_lons = np.append(obs_lons,obs_site_group.longitude)
    obs_alt = np.append(obs_alt,obs_site_group.altitude)
    obs_date = obs_site_group.variables['date'][:]
    obs_time = obs_site_group.variables['time'][:]
    
for i in range(len(obs_refs)):
    obs_refs[i] = obs_refs[i].lower()

print len(obs_refs)
    
obs_data = np.array(obs_data)
test = obs_data == -99999
obs_data[test] = np.NaN

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
#load in model data
nox_a_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150411b/Run/'
nox_b_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150411d/Run/'
nox_c_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150421a/Run/'
nox_d_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150411c/Run/'
nox_e_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150411a/Run/'
isop_a_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150520b/Run/'
isop_b_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150519b/Run/'
isop_c_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150421a/Run/'
isop_d_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150518b/Run/'
isop_e_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150521b/Run/'
meth_a_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150522b/Run/'
meth_b_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150524b/Run/' 
meth_c_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150421a/Run/'
meth_d_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150525b/Run/'
meth_e_file = '/work/home/pe534/Reactivity/PresentDay/GC_20150523b/Run/'

ctm_f_nox_a = open_ctm_bpch(nox_a_file,'ctm.bpch')
nox_a_data = np.transpose(get_gc_data_np(ctm_f_nox_a, species,category="IJ-AVG-$", debug=False))
ctm_f_nox_b = open_ctm_bpch(nox_b_file,'ctm.bpch')
nox_b_data = np.transpose(get_gc_data_np(ctm_f_nox_b, species,category="IJ-AVG-$", debug=False)) 
ctm_f_nox_c = open_ctm_bpch(nox_c_file,'ctm.bpch')
nox_c_data = np.transpose(get_gc_data_np(ctm_f_nox_c, species,category="IJ-AVG-$", debug=False))
ctm_f_nox_d = open_ctm_bpch(nox_d_file,'ctm.bpch')
nox_d_data = np.transpose(get_gc_data_np(ctm_f_nox_d, species,category="IJ-AVG-$", debug=False)) 
ctm_f_nox_e = open_ctm_bpch(nox_e_file,'ctm.bpch')
nox_e_data = np.transpose(get_gc_data_np(ctm_f_nox_e, species,category="IJ-AVG-$", debug=False)) 
ctm_f_isop_a = open_ctm_bpch(isop_a_file,'ctm.bpch')
isop_a_data = np.transpose(get_gc_data_np(ctm_f_isop_a, species,category="IJ-AVG-$", debug=False)) 
ctm_f_isop_b = open_ctm_bpch(isop_b_file,'ctm.bpch')
isop_b_data = np.transpose(get_gc_data_np(ctm_f_isop_b, species,category="IJ-AVG-$", debug=False))
ctm_f_isop_c = open_ctm_bpch(isop_c_file,'ctm.bpch')
isop_c_data = np.transpose(get_gc_data_np(ctm_f_isop_c, species,category="IJ-AVG-$", debug=False))
ctm_f_isop_d = open_ctm_bpch(isop_d_file,'ctm.bpch')
isop_d_data = np.transpose(get_gc_data_np(ctm_f_isop_d, species,category="IJ-AVG-$", debug=False))
ctm_f_isop_e = open_ctm_bpch(isop_e_file,'ctm.bpch')
isop_e_data = np.transpose(get_gc_data_np(ctm_f_isop_e, species,category="IJ-AVG-$", debug=False))

ctm_f_meth_a = open_ctm_bpch(meth_a_file,'ctm.bpch')
meth_a_data = np.transpose(get_gc_data_np(ctm_f_meth_a, species,category="IJ-AVG-$", debug=False))
ctm_f_meth_b = open_ctm_bpch(meth_b_file,'ctm.bpch')
meth_b_data = np.transpose(get_gc_data_np(ctm_f_meth_b, species,category="IJ-AVG-$", debug=False))
ctm_f_meth_c = open_ctm_bpch(meth_c_file,'ctm.bpch')
meth_c_data = np.transpose(get_gc_data_np(ctm_f_meth_c, species,category="IJ-AVG-$", debug=False))
ctm_f_meth_d = open_ctm_bpch(meth_d_file,'ctm.bpch')
meth_d_data = np.transpose(get_gc_data_np(ctm_f_meth_d, species,category="IJ-AVG-$", debug=False))
ctm_f_meth_e = open_ctm_bpch(meth_e_file,'ctm.bpch')
meth_e_data = np.transpose(get_gc_data_np(ctm_f_meth_e, species,category="IJ-AVG-$", debug=False))

cut_lat_n = []
cut_lon_n = []
#cut data into sites for O3 PERIODIC 2005-2010 sites 
for obs_lat,obs_lon in zip(obs_lats, obs_lons):    
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    cut_lat_n.append(lat_n)
    cut_lon_n.append(lon_n)

nox_a_data = nox_a_data[:,0,cut_lat_n,cut_lon_n]*1e9
nox_b_data = nox_b_data[:,0,cut_lat_n,cut_lon_n]*1e9
nox_c_data = nox_c_data[:,0,cut_lat_n,cut_lon_n]*1e9
nox_d_data = nox_d_data[:,0,cut_lat_n,cut_lon_n]*1e9
nox_e_data = nox_e_data[:,0,cut_lat_n,cut_lon_n]*1e9

isop_a_data = isop_a_data[:,0,cut_lat_n,cut_lon_n]*1e9
isop_b_data = isop_b_data[:,0,cut_lat_n,cut_lon_n]*1e9
isop_c_data = isop_c_data[:,0,cut_lat_n,cut_lon_n]*1e9
isop_d_data = isop_d_data[:,0,cut_lat_n,cut_lon_n]*1e9
isop_e_data = isop_e_data[:,0,cut_lat_n,cut_lon_n]*1e9

meth_a_data = meth_a_data[:,0,cut_lat_n,cut_lon_n]*1e9
meth_b_data = meth_b_data[:,0,cut_lat_n,cut_lon_n]*1e9
meth_c_data = meth_c_data[:,0,cut_lat_n,cut_lon_n]*1e9
meth_d_data = meth_d_data[:,0,cut_lat_n,cut_lon_n]*1e9
meth_e_data = meth_e_data[:,0,cut_lat_n,cut_lon_n]*1e9

print obs_data.shape
obs_data = obs_data[:,6:30]



print obs_data.shape
print nox_a_data.shape

#get monthly average for each type

start_dt = datetime.datetime(2005,7,1,0,0)
end_dt = datetime.datetime(2007,7,1,0,0)
time_pd = pd.date_range(start = start_dt,end = end_dt, freq = 'M')

array1 = range(0,12)
array2 = range(12,24)

nox_a_ave = []
nox_b_ave = []
nox_c_ave = []
nox_d_ave = []
nox_e_ave = []
isop_a_ave = []
isop_b_ave = []
isop_c_ave = []
isop_d_ave = []                                                                                                                                                                                                                                
isop_e_ave = []
meth_a_ave = []
meth_b_ave = []
meth_c_ave = []
meth_d_ave = []
meth_e_ave = []
obs_ave = []

for n in range(len(obs_refs)):
    nox_a_year_ave = []
    nox_b_year_ave = []
    nox_c_year_ave = []
    nox_d_year_ave = []
    nox_e_year_ave = []
    isop_a_year_ave = []
    isop_b_year_ave = []
    isop_c_year_ave = []
    isop_d_year_ave = []
    isop_e_year_ave = [] 
    meth_a_year_ave = []
    meth_b_year_ave = []
    meth_c_year_ave = []
    meth_d_year_ave = []
    meth_e_year_ave = []
    obs_year_ave = []
    for x in range(12):
        nox_a_year_ave.append(np.average([nox_a_data[array1[x],n],nox_a_data[array2[x],n]]))
        nox_b_year_ave.append(np.average([nox_b_data[array1[x],n],nox_b_data[array2[x],n]]))
        nox_c_year_ave.append(np.average([nox_c_data[array1[x],n],nox_c_data[array2[x],n]]))
        nox_d_year_ave.append(np.average([nox_d_data[array1[x],n],nox_d_data[array2[x],n]]))
        nox_e_year_ave.append(np.average([nox_e_data[array1[x],n],nox_e_data[array2[x],n]]))
        isop_a_year_ave.append(np.average([isop_a_data[array1[x],n],isop_a_data[array2[x],n]]))
        isop_b_year_ave.append(np.average([isop_b_data[array1[x],n],isop_b_data[array2[x],n]]))
        isop_c_year_ave.append(np.average([isop_c_data[array1[x],n],isop_c_data[array2[x],n]]))
        isop_d_year_ave.append(np.average([isop_d_data[array1[x],n],isop_d_data[array2[x],n]]))                                                                                                                                                      
        isop_e_year_ave.append(np.average([isop_e_data[array1[x],n],isop_e_data[array2[x],n]]))    
        
        meth_a_year_ave.append(np.average([meth_a_data[array1[x],n],meth_a_data[array2[x],n]]))
        meth_b_year_ave.append(np.average([meth_b_data[array1[x],n],meth_b_data[array2[x],n]]))
        meth_c_year_ave.append(np.average([meth_c_data[array1[x],n],meth_c_data[array2[x],n]]))
        meth_d_year_ave.append(np.average([meth_d_data[array1[x],n],meth_d_data[array2[x],n]]))                                                                                   
        meth_e_year_ave.append(np.average([meth_e_data[array1[x],n],meth_e_data[array2[x],n]]))
    
        obs_year_ave.append(np.average([obs_data[n,array1[x]],obs_data[n,array2[x]]]))    

    nox_a_ave.append(nox_a_year_ave)
    nox_b_ave.append(nox_b_year_ave)
    nox_c_ave.append(nox_c_year_ave)
    nox_d_ave.append(nox_d_year_ave)
    nox_e_ave.append(nox_e_year_ave)
    isop_a_ave.append(isop_a_year_ave)
    isop_b_ave.append(isop_b_year_ave)
    isop_c_ave.append(isop_c_year_ave)
    isop_d_ave.append(isop_d_year_ave)
    isop_e_ave.append(isop_e_year_ave)
    meth_a_ave.append(meth_a_year_ave)
    meth_b_ave.append(meth_b_year_ave)
    meth_c_ave.append(meth_c_year_ave)
    meth_d_ave.append(meth_d_year_ave)
    meth_e_ave.append(meth_e_year_ave)
    obs_ave.append(obs_year_ave)

nox_a_ave = np.array(nox_a_ave)
nox_b_ave = np.array(nox_b_ave)
nox_c_ave = np.array(nox_c_ave)
nox_d_ave = np.array(nox_d_ave)
nox_e_ave = np.array(nox_e_ave)
isop_a_ave = np.array(isop_a_ave)
isop_b_ave = np.array(isop_b_ave)
isop_c_ave = np.array(isop_c_ave)
isop_d_ave = np.array(isop_d_ave)
isop_e_ave = np.array(isop_e_ave)
meth_a_ave = np.array(meth_a_ave)
meth_b_ave = np.array(meth_b_ave)
meth_c_ave = np.array(meth_c_ave)
meth_d_ave = np.array(meth_d_ave)
meth_e_ave = np.array(meth_e_ave)
obs_ave = np.array(obs_ave)

start_dt = datetime.datetime(2005,1,1,0,0)
end_dt = datetime.datetime(2006,1,1,0,0)                                                                                                                                                                                                      
time_pd = pd.date_range(start = start_dt,end = end_dt, freq = 'M')

print obs_data.shape

#-----------------------------------
#get area
areas = ['ANT','S_O','OC','AF','SA','NE_NA','CE_NA','SE_NA','S_NA','SW_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','NE_AS','SE_AS','C_AS','S_AS','N_O','AL','ARC']

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

    area_grid = area_boundaries[area]
    area_tag = area_tags[area]
    area_label = area_labels[area]

    cut_test = modules.area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag)

    nox_a = np.nanmean(nox_a_ave[cut_test,:],axis=0)
    nox_b = np.nanmean(nox_b_ave[cut_test,:],axis=0)
    nox_c = np.nanmean(nox_c_ave[cut_test,:],axis=0)
    nox_d = np.nanmean(nox_d_ave[cut_test,:],axis=0)
    nox_e = np.nanmean(nox_e_ave[cut_test,:],axis=0)
    isop_a = np.nanmean(isop_a_ave[cut_test,:],axis=0)
    isop_b = np.nanmean(isop_b_ave[cut_test,:],axis=0)
    isop_c = np.nanmean(isop_c_ave[cut_test,:],axis=0)
    isop_d = np.nanmean(isop_d_ave[cut_test,:],axis=0)
    isop_e = np.nanmean(isop_e_ave[cut_test,:],axis=0)
    meth_a = np.nanmean(meth_a_ave[cut_test,:],axis=0)
    meth_b = np.nanmean(meth_b_ave[cut_test,:],axis=0)
    meth_c = np.nanmean(meth_c_ave[cut_test,:],axis=0)
    meth_d = np.nanmean(meth_d_ave[cut_test,:],axis=0)
    meth_e = np.nanmean(meth_e_ave[cut_test,:],axis=0)
    obs = np.nanmean(obs_ave[cut_test,:],axis=0)

    nox_a = np.concatenate((nox_a[6:],nox_a[:6]))
    nox_b = np.concatenate((nox_b[6:],nox_b[:6]))
    nox_c = np.concatenate((nox_c[6:],nox_c[:6]))
    nox_d = np.concatenate((nox_d[6:],nox_d[:6]))
    nox_e = np.concatenate((nox_e[6:],nox_e[:6]))
    isop_a = np.concatenate((isop_a[6:],isop_a[:6]))
    isop_b = np.concatenate((isop_b[6:],isop_b[:6]))
    isop_c = np.concatenate((isop_c[6:],isop_c[:6]))
    isop_d = np.concatenate((isop_d[6:],isop_d[:6]))
    isop_e = np.concatenate((isop_e[6:],isop_e[:6]))
    meth_a = np.concatenate((meth_a[6:],meth_a[:6]))
    meth_b = np.concatenate((meth_b[6:],meth_b[:6]))
    meth_c = np.concatenate((meth_c[6:],meth_c[:6]))
    meth_d = np.concatenate((meth_d[6:],meth_d[:6]))
    meth_e = np.concatenate((meth_e[6:],meth_e[:6]))
    obs = np.concatenate((obs[6:],obs[:6]))

    if param == 'NOX':
        ax.plot_date(time_pd.to_pydatetime(),obs,color = 'black',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),nox_a,color = 'blue',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),nox_b,color = 'purple',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')    
        ax.plot_date(time_pd.to_pydatetime(),nox_c,color = 'green',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),nox_d,color = 'red',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),nox_e,color = 'orange',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
    elif param == 'ISOP':        
        ax.plot_date(time_pd.to_pydatetime(),obs,color = 'black',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),isop_a,color = 'blue',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),isop_b,color = 'purple',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')    
        ax.plot_date(time_pd.to_pydatetime(),isop_c,color = 'green',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),isop_d,color = 'red',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),isop_e,color = 'orange',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
    elif param == 'METH':        
        ax.plot_date(time_pd.to_pydatetime(),obs,color = 'black',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),meth_a,color = 'blue',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),meth_b,color = 'purple',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),meth_c,color = 'green',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),meth_d,color = 'red',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')
        ax.plot_date(time_pd.to_pydatetime(),meth_e,color = 'orange',linestyle='--',linewidth=1,markersize=3,markeredgecolor='None')


    #ax.set_xlabel('Time',fontsize=25)
    #ax.set_ylabel('Concentration (ppb)',fontsize=25)
    ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)
    #ax.fill_between(obs_time_pd.to_pydatetime(),ave_obs_waveform, ave_model_waveform,where=ave_model_waveform>ave_obs_waveform, facecolor='yellow', interpolate=True)
    #ax.fill_between(obs_time_pd.to_pydatetime(),ave_obs_waveform, ave_model_waveform,where=ave_obs_waveform>ave_model_waveform, facecolor='blue', interpolate=True)
    ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
    
    tst = False
    for t in range(len(nox_a)):
        b = np.isnan(nox_a[t])
        if b == False:
            tst= True
    if tst == False:
        ax.plot_date(time_pd.to_pydatetime(),len(time_pd)*[1],markersize=0.000001)
    
    count+=1
    
plt.tight_layout(pad = 3.08)

h1, = ax.plot([1,1],color='black',marker='o',linestyle='None',markersize=10)
h2, = ax.plot([1,1],color='blue',marker='o',linestyle='None',markersize=10)
h3, = ax.plot([1,1],color='purple',marker='o',linestyle='None',markersize=10)
h4, = ax.plot([1,1],color='green',marker='o',linestyle='None',markersize=10)
h5, = ax.plot([1,1],color='red',marker='o',linestyle='None',markersize=10)
h6, = ax.plot([1,1],color='orange',marker='o',linestyle='None',markersize=10)

plt.legend((h1,h2,h3,h4,h5,h6),['Observations','%s 0.5'%(param),'%s 0.75'%(param),'%s 1'%(param),'%s 1.5'%(param),'%s 2'%(param)],loc='lower left',prop={'size':14},fancybox=True,ncol=1,markerscale=1,bbox_to_anchor=(-0.25,-0.25))
h1.set_visible(False)
h2.set_visible(False)

plt.show()
