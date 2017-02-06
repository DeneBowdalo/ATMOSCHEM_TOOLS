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

start_year = 1970
end_year = 2015

n_years = end_year - start_year
n_months = n_years*12

start = datetime.datetime(year = 1970, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = 2015, month = 1, day = 1, hour = 0, minute = 0)

full_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='M')]
full_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='M')]

dates_year = [full_dates[c] for c in range(0,n_months,12)]

full_time_pd = pd.date_range(start = start,end = end, freq = 'M')

#------------------------------------------------------------

model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_%s_2005_2010_v902_2x2.5_GEOS5_M_*.nc'%(species)
model_ts_grp = Dataset(model_fname)
geoschem_o3 = model_ts_grp.variables[species.lower()][:]*1e9
geoschem_date = model_ts_grp.variables['date'][:]
geoschem_time = model_ts_grp.variables['time'][:]
geoschem_lat_e = model_ts_grp.variables['lat_edges'][:]
geoschem_lon_e = model_ts_grp.variables['lon_edges'][:]
geoschem_lat_c = model_ts_grp.variables['lat_centre'][:]
geoschem_lon_c = model_ts_grp.variables['lon_centre'][:]

test = geoschem_o3 < 0
geoschem_o3[test] = np.NaN

#------------------------------------------------------------
#read in obs time series data
obs_ts_grp = Dataset('/work/home/db876/observations/surface/%s/process/GLOBAL_SURFACE_%s_1970_2015_M_ALL.nc'%(species,species))
obs_refs_dict_o3 = obs_ts_grp.groups

obs_refs_o3 = []
obs_lats_o3 = []
obs_lons_o3 = []
obs_alt_o3 = []
obs_o3 = []
obs_country_o3 = []

for i in obs_refs_dict_o3.keys():
    i = i.encode('ascii')
    obs_refs_o3 = np.append(obs_refs_o3,i)

for ref in obs_refs_o3:
    obs_site_group = obs_ts_grp.groups[ref] 
    obs_country_o3 = np.append(obs_country_o3,obs_site_group.country)
    obs_o3.append(obs_site_group.variables[species.lower()][:])
    obs_lats_o3 = np.append(obs_lats_o3,obs_site_group.latitude)
    obs_lons_o3 = np.append(obs_lons_o3,obs_site_group.longitude)
    obs_alt_o3 = np.append(obs_alt_o3,obs_site_group.altitude)
    obs_date = obs_site_group.variables['date'][:]
    obs_time = obs_site_group.variables['time'][:]
obs_date = obs_site_group.variables['date'][:]
obs_time = obs_site_group.variables['time'][:]

for i in range(len(obs_refs_o3)):
    obs_refs_o3[i] = obs_refs_o3[i].lower()

#get obs lat_lon grid central points
obs_lats_centre, obs_lons_centre, geoschem_model_indices = modules.grid_obs_centre_convergance(geoschem_lat_e,geoschem_lon_e,obs_lats_o3,obs_lons_o3)

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags_o3 = modules.get_tags(obs_refs_o3)

#cut model data
geoschem_o3 = geoschem_o3[:,geoschem_model_indices[:,0],geoschem_model_indices[:,1]]

obs_o3 = np.array(obs_o3)
#--------------------------------------------------
areas = ['ANT','S_O','OC','AF','SA','NE_NA','CE_NA','SE_NA','S_NA','SW_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','NE_AS','SE_AS','C_AS','S_AS','N_O','AL','ARC']

area_boundaries,area_tags,area_labels = modules.area_dicts()

fig, axes = plt.subplots(nrows=5, ncols=5,figsize=(19,13))
fig.patch.set_facecolor('white')
count = 0

all_dates = []
year_count = 1
for i in range(len(full_time_pd)):
    if year_count == 12:
        all_dates=np.append(all_dates,full_time_pd[i-11])
        year_count = 1
    else:    
        year_count+=1

inds = range(0,len(full_time_pd),12)
all_dates = full_time_pd[inds]

for ax in axes.flat:
    try:
        area = areas[count]
    except:
        ax.axis('off')
        continue

    print area
    obs_ave_o3 = []
    obs_twentyfifth_o3 = []
    obs_seventyfifth_o3 = []
    geoschem_ave_o3 = []
    geoschem_twentyfifth_o3 = []
    geoschem_seventyfifth_o3 = [] 

    area_grid = area_boundaries[area]
    area_tag = area_tags[area]
    area_label = area_labels[area]
    
    cut_test_o3 = modules.area_cut(area,obs_lats_o3,obs_lons_o3,tags_o3,area_grid,area_tag)

    obs_cut_o3 = obs_o3[cut_test_o3,:]
    geoschem_cut_o3 = geoschem_o3[:,cut_test_o3]

    n_sites_o3 = obs_cut_o3.shape[0]
    
    print 'nsites = ',n_sites_o3

    if n_sites_o3 > 0:
        for m in range(len(dates_year)):
            obs_area_year = []
            geoschem_area_year = []

            print dates_year[m]
            obs_s = np.where(obs_date ==int(dates_year[m]))[0]
            geoschem_s = np.where(geoschem_date == int(dates_year[m]))[0]                    
            obs_s = obs_s[0]

            if m < (len(dates_year)-1):
                obs_e = np.where(obs_date ==  int(dates_year[m+1]))[0]
                geoschem_e = np.where(geoschem_date == int(dates_year[m+1]))[0]
                obs_e = obs_e[0]

                obs_cut_2 = obs_cut_o3[:,obs_s:obs_e]

                if (len(geoschem_s) == 0) or (len(geoschem_e) == 0):
                    geoschem_cut_2 = np.empty((12,n_sites_o3))
                    geoschem_cut_2[:] = np.NaN
                else:
                    geoschem_s = geoschem_s[0]
                    geoschem_e = geoschem_e[0]
                    geoschem_cut_2 = geoschem_cut_o3[geoschem_s:geoschem_e,:] 
            else:
                obs_cut_2 = obs_cut_o3[:,obs_s:]

                if (len(geoschem_s) == 0):
                    geoschem_cut_2 = np.empty((12,n_sites_o3))
                    geoschem_cut_2[:] = np.NaN
                else:
                    geoschem_s = geoschem_s[0]
                    geoschem_cut_2 = geoschem_cut_o3[geoschem_s:,:] 


            for t in range(n_sites_o3):
                valid_test = obs_cut_2[t,:] != -99999
                obs_year_site = obs_cut_2[t,valid_test]
                geoschem_year_site = geoschem_cut_2[valid_test,t]

                if len(obs_year_site) >= 12:
                    obs_area_year.append(np.average(obs_year_site))
                    geoschem_area_year.append(np.average(geoschem_year_site))

            if len(obs_area_year) > 0:
                obs_ave_o3=np.append(obs_ave_o3,np.average(obs_area_year))
                obs_twentyfifth_o3 = np.append(obs_twentyfifth_o3,np.percentile(obs_area_year, 25))
                obs_seventyfifth_o3 = np.append(obs_seventyfifth_o3,np.percentile(obs_area_year, 75))
                geoschem_ave_o3 = np.append(geoschem_ave_o3,np.average(geoschem_area_year))
                geoschem_twentyfifth_o3 = np.append(geoschem_twentyfifth_o3,np.percentile(geoschem_area_year, 25))
                geoschem_seventyfifth_o3 = np.append(geoschem_seventyfifth_o3,np.percentile(geoschem_area_year, 75))
            else:
                obs_ave_o3 = np.append(obs_ave_o3,np.NaN)
                obs_twentyfifth_o3 = np.append(obs_twentyfifth_o3,np.NaN)
                obs_seventyfifth_o3 = np.append(obs_seventyfifth_o3,np.NaN)
                geoschem_ave_o3 = np.append(geoschem_ave_o3,np.NaN)
                geoschem_twentyfifth_o3 = np.append(geoschem_twentyfifth_o3,np.NaN)
                geoschem_seventyfifth_o3 = np.append(geoschem_seventyfifth_o3,np.NaN)

        ax.plot_date(all_dates.to_pydatetime(),obs_ave_o3,color = 'black',linestyle='--',linewidth=1,marker='o',markersize=5,markeredgecolor='None')                           
        ax.fill_between(all_dates.to_pydatetime(),obs_ave_o3, obs_seventyfifth_o3,where=obs_seventyfifth_o3>obs_ave_o3, facecolor='black', interpolate=True,alpha=0.4)
        ax.fill_between(all_dates.to_pydatetime(),obs_ave_o3, obs_twentyfifth_o3,where=obs_twentyfifth_o3<obs_ave_o3, facecolor='black', interpolate=True,alpha=0.4)
        ax.plot_date(all_dates.to_pydatetime(),geoschem_ave_o3,color = 'red',linestyle='--',linewidth=1,marker='o',markersize=5,markeredgecolor='None')
        ax.fill_between(all_dates.to_pydatetime(),geoschem_ave_o3, geoschem_seventyfifth_o3,where=geoschem_seventyfifth_o3>geoschem_ave_o3, facecolor='red', interpolate=True,alpha=0.4)
        ax.fill_between(all_dates.to_pydatetime(),geoschem_ave_o3, geoschem_twentyfifth_o3,where=geoschem_twentyfifth_o3<geoschem_ave_o3, facecolor='red', interpolate=True,alpha=0.4)
 
    ax.set_xlim([datetime.date(1990, 1, 1), datetime.date(2013, 1, 1)])
    ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)
    ax.xaxis.set_major_formatter(dates.DateFormatter('%y'))
    
    count+=1
    
plt.tight_layout(pad = 3.08)

h1, = ax.plot([1,1],color='black',marker='o',linestyle='None',markersize=10)
h2, = ax.plot([1,1],color='red',marker='o',linestyle='None',markersize=10)

plt.legend((h1,h2),['Observations','GEOS-Chem v902'],loc='lower left',prop={'size':14},fancybox=True,ncol=1,markerscale=2,bbox_to_anchor=(-0.2,0.0))
h1.set_visible(False)
h2.set_visible(False)

plt.show()
