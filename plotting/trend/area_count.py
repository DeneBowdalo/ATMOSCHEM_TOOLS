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

species = os.getcwd().split('/')[-1]
print species

start = datetime.datetime(year = 1970, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = 2015, month = 1, day = 1, hour = 0, minute = 0)

full_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='M')]
full_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='M')]

full_time_pd = pd.date_range(start = start,end = end, freq = 'M')

#------------------------------------------------------------

model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_%s_2005_2010_v902_2x2.5_GEOS5_M_*.nc'%(species)

model_ts_grp = Dataset(model_fname)
model_var = model_ts_grp.variables[species.lower()][:]*1e9
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
obs_ts_grp = Dataset('/work/home/db876/observations/surface/%s/process/GLOBAL_SURFACE_%s_1970_2015_M_ALL.nc'%(species,species))
obs_refs_dict = obs_ts_grp.groups

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
    obs_site_group = obs_ts_grp.groups[ref] 
    obs_country = np.append(obs_country,obs_site_group.country)
    obs_data.append(obs_site_group.variables[species.lower()][:])
    obs_lats = np.append(obs_lats,obs_site_group.latitude)
    obs_lons = np.append(obs_lons,obs_site_group.longitude)
    obs_alt = np.append(obs_alt,obs_site_group.altitude)
    obs_date = obs_site_group.variables['date'][:]
    obs_time = obs_site_group.variables['time'][:]
obs_date = obs_site_group.variables['date'][:]
obs_time = obs_site_group.variables['time'][:]

for i in range(len(obs_refs)):
    obs_refs[i] = obs_refs[i].lower()

#get obs lat_lon grid central points
obs_lats_centre, obs_lons_centre, model_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_lats,obs_lons) 

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

#cut model data
model_data = model_var[:,model_indices[:,0],model_indices[:,1]]

#remove -99999's
obs_data = np.array(obs_data)
#test = obs_data < 0
#obs_data[test] = np.nan

#--------------------------------------------------
areas = ['ANT','S_O','OC','AF','SA','NE_NA','CE_NA','SE_NA','S_NA','SW_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','NE_AS','SE_AS','C_AS','S_AS','N_O','AL','ARC']

area_boundaries,area_tags,area_labels = modules.area_dicts()

fig, axes = plt.subplots(nrows=5, ncols=5,figsize=(19,13))
fig.patch.set_facecolor('white')
count = 0

print len(axes.flat)

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

    obs_cut = obs_data[cut_test,:]
    model_cut = model_data[:,cut_test]
    
    for m in range(len(full_dates)):
        obs_test = obs_date ==  int(full_dates[m])
            
        obs_cut_2 = obs_cut[:,obs_test]
        
        valid_test = obs_cut_2 != -99999
        obs_cut_3 = obs_cut_2[valid_test]
            
        ax.plot_date(full_time_pd[m].to_pydatetime(),len(obs_cut_3),color = 'black',linestyle='_',linewidth=1,marker='o',markersize=5,markeredgecolor='None')
    
    ax.set_xlim([datetime.date(1990, 1, 1), datetime.date(2013, 1, 1)])
    ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)
    ax.xaxis.set_major_formatter(dates.DateFormatter('%y'))
    
    count+=1
    
plt.tight_layout(pad = 3.08)

plt.show()
