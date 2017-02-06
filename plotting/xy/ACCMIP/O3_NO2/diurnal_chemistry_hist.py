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

model_fname = '/work/home/db876/plotting_tools/model_files/MIROCCHEM_SURFACE_%s_2005_2010_*_*_*_H_*.nc'%(species)
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
obs_fname = '/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_%s_2005_2010_H_PERIODIC.nc'%(species)
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

obs_period_grp = Dataset('../obs_SURFACE_H/obs_sig_periods.nc')
cesmcam_period_grp = Dataset('../CESMCAM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
cmam_period_grp = Dataset('../CMAM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
geosccm_period_grp = Dataset('../GEOSCCM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
geoschemv902_period_grp = Dataset('../GEOSCHEM_SURFACE_v902_2x2.5_GEOS5_H_*/model_sig_periods.nc')
gfdlam3_period_grp = Dataset('../GFDLAM3_SURFACE_*_*_*_H_*/model_sig_periods.nc')
gisse2r_period_grp = Dataset('../GISSE2R_SURFACE_*_*_*_H_*/model_sig_periods.nc')
mirocchem_period_grp = Dataset('../MIROCCHEM_SURFACE_*_*_*_H_*/model_sig_periods.nc') 

obs_daily_amp = []
cesmcam_daily_amp = []
cmam_daily_amp = []
geosccm_daily_amp = []
geoschemv902_daily_amp = []
gfdlam3_daily_amp = []
gisse2r_daily_amp = []
mirocchem_daily_amp = []
obs_ave = []
cesmcam_ave = []
cmam_ave = []
geosccm_ave = []
geoschemv902_ave = []
gfdlam3_ave = []
gisse2r_ave = []
mirocchem_ave = []

for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    site_group_cesmcam = cesmcam_period_grp.groups[ref]
    site_group_cmam = cmam_period_grp.groups[ref]    
    site_group_geosccm = geosccm_period_grp.groups[ref]
    site_group_geoschemv902 = geoschemv902_period_grp.groups[ref]
    site_group_gfdlam3 = gfdlam3_period_grp.groups[ref]
    site_group_gisse2r = gisse2r_period_grp.groups[ref]
    site_group_mirocchem = mirocchem_period_grp.groups[ref]

    cesmcam_daily_amp.append(site_group_cesmcam.daily_amplitude)
    cmam_daily_amp.append(site_group_cmam.daily_amplitude)
    geosccm_daily_amp.append(site_group_geosccm.daily_amplitude)                                                                                              
    geoschemv902_daily_amp.append(site_group_geoschemv902.daily_amplitude)                                                                                
    gfdlam3_daily_amp.append(site_group_gfdlam3.daily_amplitude)
    gisse2r_daily_amp.append(site_group_gisse2r.daily_amplitude)
    mirocchem_daily_amp.append(site_group_mirocchem.daily_amplitude)                                                                                        
    obs_daily_amp.append(site_group_obs.daily_amplitude)
    #cesmcam_daily_amp.append(site_group_cesmcam.seasonal_amplitude)
    #cmam_daily_amp.append(site_group_cmam.seasonal_amplitude)
    #geosccm_daily_amp.append(site_group_geosccm.seasonal_amplitude)                                                                                              
    #geoschemv902_daily_amp.append(site_group_geoschemv902.seasonal_amplitude)                                                                                
    #gfdlam3_daily_amp.append(site_group_gfdlam3.seasonal_amplitude)
    #gisse2r_daily_amp.append(site_group_gisse2r.seasonal_amplitude)
    #mirocchem_daily_amp.append(site_group_mirocchem.seasonal_amplitude)                                                                                        
    #obs_daily_amp.append(site_group_obs.seasonal_amplitude)
    
    cesmcam_ave.append(site_group_cesmcam.average)
    cmam_ave.append(site_group_cmam.average)
    geosccm_ave.append(site_group_geosccm.average)                                                                                              
    geoschemv902_ave.append(site_group_geoschemv902.average)                                                                                
    gfdlam3_ave.append(site_group_gfdlam3.average)
    gisse2r_ave.append(site_group_gisse2r.average)
    mirocchem_ave.append(site_group_mirocchem.average)                                                                                        
    obs_ave.append(site_group_obs.average) 

plt.plot(obs_ave,obs_daily_amp,linestyle='None',marker='x')
plt.show()


obs_daily_amp = np.average(obs_daily_amp)
cesmcam_daily_amp = np.average(cesmcam_daily_amp)
cmam_daily_amp = np.average(cmam_daily_amp)                                                                                                                  
geosccm_daily_amp = np.average(geosccm_daily_amp)                                                                                                                  
geoschemv902_daily_amp = np.average(geoschemv902_daily_amp)                                                                                                          
gfdlam3_daily_amp = np.average(gfdlam3_daily_amp)                                                                                                                
gisse2r_daily_amp = np.average(gisse2r_daily_amp)                                                                                                                
mirocchem_daily_amp = np.average(mirocchem_daily_amp)                                                                                                                

#-----------------------------------
#get area

area_boundaries,area_tags,area_labels = modules.area_dicts()

fig, ax = plt.subplots(1,figsize=(19,13))
fig.patch.set_facecolor('white')

means = [obs_daily_amp,cesmcam_daily_amp,cmam_daily_amp,geosccm_daily_amp,geoschemv902_daily_amp,gfdlam3_daily_amp,gisse2r_daily_amp,mirocchem_daily_amp]

ax.bar(0, obs_daily_amp, 1, color='black')
ax.bar(1, cesmcam_daily_amp, 1, color = plt.cm.jet(nmodels_linspace[0]))
ax.bar(2, cmam_daily_amp, 1, color = plt.cm.jet(nmodels_linspace[1]))
ax.bar(3, geosccm_daily_amp, 1, color = plt.cm.jet(nmodels_linspace[2]))
ax.bar(4, geoschemv902_daily_amp, 1, color = plt.cm.jet(nmodels_linspace[3]))
ax.bar(5, gfdlam3_daily_amp, 1, color = plt.cm.jet(nmodels_linspace[4]))
ax.bar(6, gisse2r_daily_amp, 1, color = plt.cm.jet(nmodels_linspace[5]))
ax.bar(7, mirocchem_daily_amp, 1, color = plt.cm.jet(nmodels_linspace[6]))

ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5])
ax.set_xticklabels( ('Obs.', 'CMAM', 'CESMCAM', 'GEOSCCM', 'GEOSCHEMv902', 'GFDLAM3', 'GISSE2R', 'MIROCCHEM') )
ax.set_ylabel('Average Diurnal Amplitude (ppb)')
#ax.set_ylim(3,8.5)

plt.show()
