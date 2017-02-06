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

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
obs_o3_tags = modules.get_tags(obs_o3_refs)
obs_no2_tags = modules.get_tags(obs_no2_refs)
obs_co_tags = modules.get_tags(obs_co_refs)

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

co_max = [700,700,1500,700,1100,900,900,500,2600,500]
no2_max = [50,50,30,25,30,30,30,15,40,50]

#x_ticks = [[0.1,0.25,0.5,1,2.5,5,10,25,50,100],[0.1,0.25,0.5,1,2.5,5,10,25,50,100],[0.1,0.25,0.5,1,2.5,5,10,25,50],[0.1,0.25,0.5,1,2.5,5,10,25,50],[0.5,1,2.5,5,10,25,50],[0.5,1,2.5,5,10,25,50],[0.5,1,2.5,5,10,25,50],[0.1,0.25,0.5,1,2.5,5,10,25],[0.5,1,2.5,5,10,25,50],[0.1,0.25,0.5,1,2.5,5,10,25,50]]
#y_ticks = [[50,100,250,500,1000],[25,50,100,250,500,1000],[50,100,250,500,1000,2500],[50,100,250,500,1000,2500],[50,100,250,500,1000,2500],[50,100,250,500,1000,2500],[50,100,250,500,1000,2500],[50,100,250,500,1000,2500],[10,25,50,100,250,500,1000,2500,5000],[25,50,100,250,500]]

#co_min = [60,50,50]
#no2_min = []

area_boundaries,area_tags,area_labels = modules.area_dicts()

fig, ax = plt.subplots(1,figsize=(19,13))
fig.patch.set_facecolor('white')
count = 0
     
start_boolean = (obs_o3_date ==  int(full_dates[0])) & (obs_o3_time ==  int(full_times[0])) 
end_boolean = (obs_o3_date ==  int(full_dates[-1])) & (obs_o3_time ==  int(full_times[-1]))
start_obs_ind = np.where(start_boolean == True)[0][0]
end_obs_ind = (np.where(end_boolean == True)[0][0]) +1

obs_o3_data = obs_o3_data[:,start_obs_ind:end_obs_ind]
obs_co_data = obs_co_data[:,start_obs_ind:end_obs_ind]
obs_no2_data = obs_no2_data[:,start_obs_ind:end_obs_ind]

obs_o3_data = np.nanmean(obs_o3_data,axis=0)
obs_no2_data = np.nanmean(obs_no2_data,axis=0)
obs_co_data = np.nanmean(obs_co_data,axis=0)
    
all_nans = np.logical_or(np.logical_or(np.isnan(obs_o3_data), np.isnan(obs_no2_data)), np.isnan(obs_co_data))
    
obs_o3_data[all_nans] = np.NaN
obs_no2_data[all_nans] = np.NaN
obs_co_data[all_nans] = np.NaN 
        
no2_max = np.nanmax(obs_no2_data)
co_max = np.nanmax(obs_co_data)

#Average data into boxes
no2_box_edges = np.arange(0,no2_max+0.001,0.1)
co_box_edges = np.arange(0,co_max+0.001,10) 
    
fill_box = np.empty((len(no2_box_edges)-1,len(co_box_edges)-1))
fill_box[:,:] = np.NaN
    
print fill_box.shape
    
for x in range(len(no2_box_edges)-1):
    for y in range(len(co_box_edges)-1):  
          
        lower_no2 = no2_box_edges[x] 
        lower_co = co_box_edges[y]
        upper_no2 = no2_box_edges[x+1]
        upper_co = co_box_edges[y+1]
            
        test = (lower_no2 <= obs_no2_data) & (upper_no2 > obs_no2_data) & (lower_co <= obs_co_data) & (upper_co > obs_co_data)
    
        box_ave_o3 = np.nanmean(obs_o3_data[test])
            
        fill_box[x,y] = box_ave_o3

    #cb = ax.scatter(np.log(obs_no2_cut_ave),np.log(obs_co_cut_ave),s=5,c=obs_o3_cut_ave,vmin=0,vmax=70,edgecolors='none')    
    #cb = ax.scatter(obs_no2_cut_ave,obs_co_cut_ave,s=5,c=obs_o3_cut_ave,vmin=0,vmax=70,edgecolors='none')
    
xx, yy = np.meshgrid(no2_box_edges, co_box_edges)
cb = ax.pcolor(xx,yy,np.transpose(fill_box),vmin=0,vmax=70)
cb.cmap.set_under('white')
ax.set_xlabel('NO2 (ppb)')
ax.set_ylabel('CO (ppb)')
    
plt.colorbar(cb,orientation='horizontal',label='O3 (ppb)')
plt.tight_layout(pad = 3.08)

plt.show()
