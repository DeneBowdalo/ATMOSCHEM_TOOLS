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

species = raw_input('O3, CO, NO, or NO2?\n') 

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
model_var = model_ts_grp.variables['%s'%(species.lower())][:]*1e9
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
    obs_data.append(obs_site_group.variables['%s'%(species.lower())][:])
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

#--------------------------------------------------
#areas = ['ANT','S_O','OC','AF','NE_NA','CE_NA','SE_NA','S_NA','SW_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','NE_AS','SE_AS','C_AS','S_AS','N_O','AL','ARC'] 
areas = ['NE_NA','CE_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','N_O'] 

area_boundaries,area_tags,area_labels = modules.area_dicts()

fig, axes = plt.subplots(nrows=3, ncols=4,figsize=(19,13))
fig.patch.set_facecolor('white')
count = 0

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

#convert all -99999's to nans
test = obs_data == -99999
obs_data[test] = np.NaN
model_data[np.transpose(test)] = np.NaN

print obs_data.shape
print model_data.shape

cm = plt.cm.get_cmap('jet')

all_colors = []
all_alpha = []

for m in range(len(full_dates)):
    print m
    current_datetime = datetime.datetime(int(full_dates[m][:4]),int(full_dates[m][4:6]),int(full_dates[m][6:]))
    spring_start = datetime.datetime(int(full_dates[m][:4]),3,20,0,0)
    summer_start = datetime.datetime(int(full_dates[m][:4]),6,21,0,0)
    autumn_start = datetime.datetime(int(full_dates[m][:4]),9,23,0,0)
    winter_start = datetime.datetime(int(full_dates[m][:4]),12,21,0,0)
    if (current_datetime < spring_start) or (current_datetime >= winter_start):
        all_colors.append(0)
        all_alpha.append(1)
        #color = 'blue'
    if (current_datetime >= spring_start) & (current_datetime < summer_start):
        all_colors.append(1)
        all_alpha.append(0.8)
        #color = 'green'
    if (current_datetime >= summer_start) & (current_datetime < autumn_start): 
        all_colors.append(2)
        all_alpha.append(0.6)
        #color = 'yellow'
    if (current_datetime >= autumn_start) & (current_datetime < winter_start): 
        all_colors.append(3)
        all_alpha.append(0.4)
        #color = 'red'

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

    #cut data to just sites in area 
    obs_cut = obs_data[cut_test,:]
    model_cut = model_data[:,cut_test]

    if obs_cut.shape[0] == 0:
        count+=1
        continue
    
    obs_cut_ave = np.nanmean(obs_cut,axis=0)
    model_cut_ave = np.nanmean(model_cut,axis=1)

    obs_min = np.nanmin(obs_cut_ave)
    obs_max = np.nanmax(obs_cut_ave)
    model_min = np.nanmin(model_cut_ave)
    model_max = np.nanmax(model_cut_ave)
    
    if obs_min < model_min:
        current_min = obs_min
    else:
        current_min = model_min

    if obs_max > model_max:
        current_max = obs_max
    else:
        current_max = model_max

    for t in range(len(all_colors)):        
        cb = ax.scatter(obs_cut_ave[t],model_cut_ave[t],s=4,c=all_colors[t],edgecolors='none',vmin=0,vmax=3,alpha=all_alpha[t])    
    ax.set_xlabel('Obs. %s (ppb)'%(species))
    ax.set_ylabel('Model %s (ppb)'%(species))        

    x = np.arange(0,3000,1)
    y = np.arange(0,3000,1)                                                                                                                                                                                                                
    ax.plot(x,y,color='black',linestyle='--')

    ax.set_xlim(current_min,current_max)
    ax.set_ylim(current_min,current_max)
    ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)
    
    count+=1
 
a = ax.plot([-999],[-999],color=cm(0.001),marker='o',markersize=20,label='Winter') 
b = ax.plot([-999],[-999],color=cm(0.333),marker='o',markersize=20,label='Spring') 
c = ax.plot([-999],[-999],color=cm(0.666),marker='o',markersize=20,label='Summer') 
d = ax.plot([-999],[-999],color=cm(0.999),marker='o',markersize=20,label='Autumn') 
l = plt.legend(prop={'size':20})

plt.tight_layout(pad = 3.08)

plt.show()
