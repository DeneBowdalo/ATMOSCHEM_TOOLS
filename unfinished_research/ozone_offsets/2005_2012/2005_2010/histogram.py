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

species = raw_input('O3, CO, NO or NO2?\n')

start_year = 2005
end_year = 2010

start = datetime.datetime(year = start_year, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = end_year, month = 1, day = 1, hour = 0, minute = 0)

full_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]
full_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='H')][:-1]
full_time_pd = pd.date_range(start = start,end = end, freq = 'H')[:-1]

#------------------------------------------------------------

model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_%s_2005_2012_v902_2x2.5_GEOS5_H_*.nc'%(species)
model_ts_grp = Dataset(model_fname)
model_var = model_ts_grp.variables[species.lower()][:]
model_time = model_ts_grp.variables['time'][:]
lat_e = model_ts_grp.variables['lat_edges'][:]
lon_e = model_ts_grp.variables['lon_edges'][:]
lat_c = model_ts_grp.variables['lat_centre'][:]
lon_c = model_ts_grp.variables['lon_centre'][:]
grid_size = model_ts_grp.variables['grid_size'][:]
grid_size = grid_size[0]


#------------------------------------------------------------
#read in obs time series data
obs_grp = Dataset('/work/home/db876/observations/surface/%s/process/GLOBAL_SURFACE_%s_1970_2015_H_ALL.nc'%(species,species))

obs_refs_dict = obs_grp.groups

obs_refs = []
obs_lats = []
obs_lons = []
obs_alt = []
obs_data = []
obs_country = []

for i in obs_refs_dict.keys():
    i = i.encode('ascii')
    obs_refs = np.append(obs_refs,i)

for ref in obs_refs:
    obs_site_group = obs_grp.groups[ref] 
    obs_country = np.append(obs_country,obs_site_group.country)
    obs_data.append(obs_site_group.variables[species.lower()][:])
    obs_lats = np.append(obs_lats,obs_site_group.latitude)
    obs_lons = np.append(obs_lons,obs_site_group.longitude)
    obs_alt = np.append(obs_alt,obs_site_group.altitude)
obs_date = obs_site_group.variables['date'][:]
obs_time = obs_site_group.variables['time'][:]

for i in range(len(obs_refs)):
    obs_refs[i] = obs_refs[i].lower()
    
#convert obs_lons to same grid as model 
if lon_e[0] < -180:
    test = obs_lons > lon_e[-1]
    diff = obs_lons[test] - lon_e[-1]
    obs_lons[test] = lon_e[0] + diff

#get obs lat_lon grid central points
obs_lats_centre, obs_lons_centre, model_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_lats,obs_lons) 

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
obs_tags = modules.get_tags(obs_refs)

#cut model data
model_data = model_var[:,model_indices[:,0],model_indices[:,1]]

obs_data = np.array(obs_data)

#convert all -99999's to nans
test = obs_data == -99999
obs_data[test] = np.NaN

#--------------------------------------------------
#areas = ['ANT','S_O','OC','AF','NE_NA','CE_NA','SE_NA','S_NA','SW_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','NE_AS','SE_AS','C_AS','S_AS','N_O','AL','ARC'] 
areas = ['NE_NA','CE_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','N_O'] 

if species == 'CO':
    max_array = [700,700,1500,700,1100,900,900,500,2600,500]
if species == 'NO2':
    max_array = [50,50,30,25,30,30,30,15,40,50]
if species == 'NO':
    max_array = [10,10,10,10,10,10,10,10,10,10]  
if species == 'O3':
    max_array = [70,70,70,70,70,70,70,70,70,70] 


area_boundaries,area_tags,area_labels = modules.area_dicts()


start_boolean = (obs_date ==  int(full_dates[0])) & (obs_time ==  int(full_times[0])) 
end_boolean = (obs_date ==  int(full_dates[-1])) & (obs_time ==  int(full_times[-1]))
start_obs_ind = np.where(start_boolean == True)[0][0]
end_obs_ind = (np.where(end_boolean == True)[0][0]) +1

start_boolean = (model_date ==  int(full_dates[0])) & (model_time ==  int(full_times[0])) 
end_boolean = (model_date ==  int(full_dates[-1])) & (model_time ==  int(full_times[-1]))
start_model_ind = np.where(start_boolean == True)[0][0]
end_model_ind = (np.where(end_boolean == True)[0][0]) +1

obs_data = obs_data[:,start_obs_ind:end_obs_ind]
model_data = model_data[start_model_ind:end_model_ind,:]

fig, axes = plt.subplots(nrows=3, ncols=4,figsize=(19,13))
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

    cut_test = modules.area_cut(area,obs_lats,obs_lons,obs_tags,area_grid,area_tag)

    obs_cut = obs_data[cut_test,:]
    
    model_cut = model_data[:,cut_test]
    
    if obs_cut.shape[0] == 0:
        count+=1
        continue
   
    #for b in range(len(obs_data)):
    #    ax.hist(obs_data[b,:],50,range=[0, max_array[count]],histtype='step',color='black',lw=2)
 
    obs_cut_ave = np.nanmean(obs_cut,axis=0)
    model_cut_ave = np.nanmean(model_cut,axis=1)
    
    all_nans = np.isnan(obs_cut_ave)
    
    obs_cut_ave[all_nans] = np.NaN
    model_cut_ave[all_nans] = np.NaN
    
    obs_cut_ave = obs_cut_ave[~np.isnan(obs_cut_ave)]
    model_cut_ave = model_cut_ave[~np.isnan(model_cut_ave)]
    
    weights_obs = np.ones_like(obs_cut_ave)/len(model_cut_ave)
    weights_model = np.ones_like(model_cut_ave)/len(model_cut_ave)
    #ax.hist(obs_cut_ave,50,range=[0, max_array[count]],weights=weights_obs,histtype='step',color='black',lw=2)
    #ax.hist(model_cut_ave,50,range=[0, max_array[count]],weights=weights_model,histtype='step',color='red',lw=2)
    ax.hist(obs_cut_ave,50,range=[0, max_array[count]],histtype='step',color='black',lw=2)
    ax.hist(model_cut_ave,50,range=[0, max_array[count]],histtype='step',color='red',lw=2) 


    ax.set_xlabel('%s (ppb)'%(species))
    ax.set_ylabel('Frac. Count')

    ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)

    count+=1
    
#plt.colorbar(cb,orientation='horizontal',label='O3 (ppb)')
plt.tight_layout(pad = 3.08)

plt.show()
