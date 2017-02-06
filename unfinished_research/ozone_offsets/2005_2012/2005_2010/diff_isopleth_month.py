#mport matplotlib
#matplotlib.use('Agg')
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
from netCDF4 import Dataset
from scipy.odr import Model, RealData, ODR, Data
from pylab import *
import pandas as pd
from scipy import stats
import matplotlib.dates as dates

start_year = 2005
end_year = 2010

start = datetime.datetime(year = start_year, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = end_year, month = 1, day = 1, hour = 0, minute = 0)

full_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]
full_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='H')][:-1]
full_time_pd = pd.date_range(start = start,end = end, freq = 'H')[:-1]

#full_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='M')]
#full_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='M')]
#full_time_pd = pd.date_range(start = start,end = end, freq = 'M')

#------------------------------------------------------------

model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_O3_2005_2012_v902_2x2.5_GEOS5_H_*.nc'
model_ts_grp = Dataset(model_fname)
model_o3_var = model_ts_grp.variables['o3'][:]*1e9
model_date = model_ts_grp.variables['date'][:]
model_time = model_ts_grp.variables['time'][:]
lat_e = model_ts_grp.variables['lat_edges'][:]
lon_e = model_ts_grp.variables['lon_edges'][:]
lat_c = model_ts_grp.variables['lat_centre'][:]
lon_c = model_ts_grp.variables['lon_centre'][:]
grid_size = model_ts_grp.variables['grid_size'][:]
grid_size = grid_size[0]

model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_CO_2005_2012_v902_2x2.5_GEOS5_H_*.nc'
model_ts_grp = Dataset(model_fname)
model_co_var = model_ts_grp.variables['co'][:]*1e9

model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_NO2_2005_2012_v902_2x2.5_GEOS5_H_*.nc'
model_ts_grp = Dataset(model_fname)
model_no2_var = model_ts_grp.variables['no2'][:]*1e9

#------------------------------------------------------------
#read in obs time series data
obs_o3_grp = Dataset('/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_O3_1970_2015_H_ALL.nc')
obs_no2_grp = Dataset('/work/home/db876/observations/surface/NO2/process/GLOBAL_SURFACE_NO2_1970_2015_H_ALL.nc')
obs_co_grp = Dataset('/work/home/db876/observations/surface/CO/process/GLOBAL_SURFACE_CO_1970_2015_H_ALL.nc')
#obs_o3_grp = Dataset('/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_O3_1970_2015_M_ALL.nc')
#obs_no2_grp = Dataset('/work/home/db876/observations/surface/NO2/process/GLOBAL_SURFACE_NO2_1970_2015_M_ALL.nc')
#obs_co_grp = Dataset('/work/home/db876/observations/surface/CO/process/GLOBAL_SURFACE_CO_1970_2015_M_ALL.nc')

obs_o3_refs_dict = obs_o3_grp.groups
obs_no2_refs_dict = obs_no2_grp.groups
obs_co_refs_dict = obs_co_grp.groups

obs_o3_refs = []
obs_o3_lats = []
obs_o3_lons = []
obs_o3_alt = []
obs_o3_data = []
obs_o3_country = []
obs_no2_refs = []
obs_no2_lats = []
obs_no2_lons = []
obs_no2_alt = []
obs_no2_data = []
obs_no2_country = []
obs_co_refs = []
obs_co_lats = []
obs_co_lons = []
obs_co_alt = []
obs_co_data = []
obs_co_country = []

for i in obs_o3_refs_dict.keys():
    i = i.encode('ascii')
    obs_o3_refs = np.append(obs_o3_refs,i)

for i in obs_no2_refs_dict.keys():
    i = i.encode('ascii')
    obs_no2_refs = np.append(obs_no2_refs,i)
    
for i in obs_co_refs_dict.keys():
    i = i.encode('ascii')
    obs_co_refs = np.append(obs_co_refs,i)

for ref in obs_o3_refs:
    obs_site_group = obs_o3_grp.groups[ref] 
    obs_o3_country = np.append(obs_o3_country,obs_site_group.country)
    obs_o3_data.append(obs_site_group.variables['o3'][:])
    obs_o3_lats = np.append(obs_o3_lats,obs_site_group.latitude)
    obs_o3_lons = np.append(obs_o3_lons,obs_site_group.longitude)
    obs_o3_alt = np.append(obs_o3_alt,obs_site_group.altitude)
obs_o3_date = obs_site_group.variables['date'][:]
obs_o3_time = obs_site_group.variables['time'][:]

for ref in obs_no2_refs:
    obs_site_group = obs_no2_grp.groups[ref] 
    obs_no2_country = np.append(obs_no2_country,obs_site_group.country)
    obs_no2_data.append(obs_site_group.variables['no2'][:])
    obs_no2_lats = np.append(obs_no2_lats,obs_site_group.latitude)
    obs_no2_lons = np.append(obs_no2_lons,obs_site_group.longitude)
    obs_no2_alt = np.append(obs_no2_alt,obs_site_group.altitude)
obs_no2_date = obs_site_group.variables['date'][:]
obs_no2_time = obs_site_group.variables['time'][:]

for ref in obs_co_refs:
    obs_site_group = obs_co_grp.groups[ref] 
    obs_co_country = np.append(obs_co_country,obs_site_group.country)
    obs_co_data.append(obs_site_group.variables['co'][:])
    obs_co_lats = np.append(obs_co_lats,obs_site_group.latitude)
    obs_co_lons = np.append(obs_co_lons,obs_site_group.longitude)
    obs_co_alt = np.append(obs_co_alt,obs_site_group.altitude)
obs_co_date = obs_site_group.variables['date'][:]
obs_co_time = obs_site_group.variables['time'][:]

for i in range(len(obs_o3_refs)):
    obs_o3_refs[i] = obs_o3_refs[i].lower()
    
for i in range(len(obs_no2_refs)):
    obs_no2_refs[i] = obs_no2_refs[i].lower()
    
for i in range(len(obs_co_refs)):
    obs_co_refs[i] = obs_co_refs[i].lower()
    
#convert obs_lons to same grid as model 
if lon_e[0] < -180:
    test = obs_o3_lons > lon_e[-1]
    diff = obs_o3_lons[test] - lon_e[-1]
    obs_o3_lons[test] = lon_e[0] + diff
        
    test = obs_no2_lons > lon_e[-1]
    diff = obs_no2_lons[test] - lon_e[-1]
    obs_no2_lons[test] = lon_e[0] + diff
    
    test = obs_co_lons > lon_e[-1]
    diff = obs_co_lons[test] - lon_e[-1]
    obs_co_lons[test] = lon_e[0] + diff

#get obs lat_lon grid central points
obs_o3_lats_centre, obs_o3_lons_centre, model_o3_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_o3_lats,obs_o3_lons) 
obs_no2_lats_centre, obs_no2_lons_centre, model_no2_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_no2_lats,obs_no2_lons) 
obs_co_lats_centre, obs_co_lons_centre, model_co_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_co_lats,obs_co_lons) 

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
obs_o3_tags = modules.get_tags(obs_o3_refs)
obs_no2_tags = modules.get_tags(obs_no2_refs)
obs_co_tags = modules.get_tags(obs_co_refs)

#cut model data
model_o3_data = model_o3_var[:,model_o3_indices[:,0],model_o3_indices[:,1]]
model_no2_data = model_no2_var[:,model_no2_indices[:,0],model_no2_indices[:,1]]
model_co_data = model_co_var[:,model_co_indices[:,0],model_co_indices[:,1]]

obs_o3_data = np.array(obs_o3_data)
obs_no2_data = np.array(obs_no2_data)
obs_co_data = np.array(obs_co_data)

#convert all -99999's to nans
test = obs_o3_data == -99999
obs_o3_data[test] = np.NaN
test = obs_no2_data == -99999
obs_no2_data[test] = np.NaN
test = obs_co_data == -99999
obs_co_data[test] = np.NaN

#--------------------------------------------------
#areas = ['ANT','S_O','OC','AF','NE_NA','CE_NA','SE_NA','S_NA','SW_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','NE_AS','SE_AS','C_AS','S_AS','N_O','AL','ARC'] 
areas = ['NE_NA','CE_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','N_O'] 

co_max = [700,700,1500,700,1100,900,900,500,2600,250]
no2_max = [50,50,30,25,30,30,30,15,40,0.15]

area_boundaries,area_tags,area_labels = modules.area_dicts()

fig, axes = plt.subplots(nrows=3, ncols=4,figsize=(19,13))
fig.patch.set_facecolor('white')
count = 0

start_boolean = (obs_o3_date ==  int(full_dates[0])) & (obs_o3_time ==  int(full_times[0])) 
end_boolean = (obs_o3_date ==  int(full_dates[-1])) & (obs_o3_time ==  int(full_times[-1]))
start_obs_ind = np.where(start_boolean == True)[0][0]
end_obs_ind = (np.where(end_boolean == True)[0][0]) +1

start_boolean = (model_date ==  int(full_dates[0])) & (model_time ==  int(full_times[0])) 
end_boolean = (model_date ==  int(full_dates[-1])) & (model_time ==  int(full_times[-1]))
start_model_ind = np.where(start_boolean == True)[0][0]
end_model_ind = (np.where(end_boolean == True)[0][0]) +1

obs_o3_data = obs_o3_data[:,start_obs_ind:end_obs_ind]
obs_co_data = obs_co_data[:,start_obs_ind:end_obs_ind]
obs_no2_data = obs_no2_data[:,start_obs_ind:end_obs_ind]

model_o3_data = model_o3_data[start_model_ind:end_model_ind,:]
model_co_data = model_co_data[start_model_ind:end_model_ind,:]
model_no2_data = model_no2_data[start_model_ind:end_model_ind,:]

for s in range(len(12)):
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

        cut_test_o3 = modules.area_cut(area,obs_o3_lats,obs_o3_lons,obs_o3_tags,area_grid,area_tag)
        cut_test_no2 = modules.area_cut(area,obs_no2_lats,obs_no2_lons,obs_no2_tags,area_grid,area_tag)
        cut_test_co = modules.area_cut(area,obs_co_lats,obs_co_lons,obs_co_tags,area_grid,area_tag)

        obs_o3_cut = obs_o3_data[cut_test_o3,:]
        obs_no2_cut = obs_no2_data[cut_test_no2,:]
        obs_co_cut = obs_co_data[cut_test_co,:]
    
        model_o3_cut = model_o3_data[:,cut_test_o3]
        model_no2_cut = model_no2_data[:,cut_test_no2]
        model_co_cut = model_co_data[:,cut_test_co]
    
        if obs_o3_cut.shape[0] == 0:
            count+=1
            continue
        if obs_no2_cut.shape[0] == 0:
            count+=1
            continue
        if obs_co_cut.shape[0] == 0:
            count+=1
            continue
    
        obs_o3_cut_ave = np.nanmean(obs_o3_cut,axis=0)
        obs_no2_cut_ave = np.nanmean(obs_no2_cut,axis=0)
        obs_co_cut_ave = np.nanmean(obs_co_cut,axis=0)
        model_o3_cut_ave = np.nanmean(model_o3_cut,axis=1)
        model_no2_cut_ave = np.nanmean(model_no2_cut,axis=1)
        model_co_cut_ave = np.nanmean(model_co_cut,axis=1)
    
        all_nans = np.logical_or(np.logical_or(np.isnan(obs_o3_cut_ave), np.isnan(obs_no2_cut_ave)), np.isnan(obs_co_cut_ave))
    
        obs_o3_cut_ave[all_nans] = np.NaN
        obs_no2_cut_ave[all_nans] = np.NaN
        obs_co_cut_ave[all_nans] = np.NaN 
        model_o3_cut_ave[all_nans] = np.NaN
        model_no2_cut_ave[all_nans] = np.NaN
        model_co_cut_ave[all_nans] = np.NaN 

        diff_o3 = obs_o3_cut_ave - model_o3_cut_ave
        diff_no2 = obs_no2_cut_ave - model_no2_cut_ave
        diff_co = obs_co_cut_ave - model_co_cut_ave
    
        no2_diff_max = np.nanmax(np.abs(diff_no2))
        co_diff_max = np.nanmax(np.abs(diff_co))

        cb = ax.scatter(diff_no2,diff_co,s=5,c=diff_o3,vmin=-40,vmax=40,cmap=plt.cm.coolwarm,edgecolors='none')    
        ax.set_xlabel('NO2 (ppb)')
        ax.set_ylabel('CO (ppb)')
        ax.set_xlim(-no2_diff_max,no2_diff_max)
        ax.set_ylim(-co_diff_max,co_diff_max)

        ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)
    
        ax.axhline(0,linestyle='--',color='black',alpha=0.8)
        ax.axvline(0,linestyle='--',color='black',alpha=0.8)

        count+=1
    
plt.colorbar(cb,orientation='horizontal',label='O3 (ppb)')
plt.tight_layout(pad = 3.08)

plt.show()
