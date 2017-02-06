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
from collections import Counter
from itertools import groupby

#-----------------------------------------------------
#get species from current directory
present_dir = os.getcwd()

species = 'O3'
start_year = 1993
end_year = 2012
nyears = 5

start_years = range(start_year,(end_year-nyears)+1)
end_years = range(start_year+nyears,end_year+1)

obsn = len(start_years)
obsn_linspace = np.linspace(0, 1, 16)

plot_type = raw_input('\nd, s or full?\n')

#fig, axes = plt.subplots(nrows=4, ncols=5,figsize=(19,12))
fig, axes = plt.subplots(1,figsize=(19,12))
fig.patch.set_facecolor('white')

#find constant refs in all files
obs_refs = []
all_refs = []
bad_refs = []
obs_dates = []
obs_times = []
obs_country = []
obs_lats = []
obs_lons = []
obs_alt = []
obs_daily_waveforms = []
obs_seasonal_waveforms = []
obs_full_waveforms = []

for c in range(len(start_years)):
    print c
    sy = start_years[c]
    ey = end_years[c]
    #------------------------------------------------------------
    #read in obs time series data
    obs_fname = '/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_O3_%s_%s_H.nc'%(sy,ey)
    obs_ts_grp = Dataset(obs_fname)
    obs_refs_dict = obs_ts_grp.groups
    new_refs = []

    for i in obs_refs_dict.keys():
        i = i.encode('ascii')
        new_refs = np.append(new_refs,i)
    
    all_refs = np.append(all_refs,new_refs)

    if c == 0:
        base_refs = np.copy(new_refs)    

    for ref in new_refs:
        obs_site_group = obs_ts_grp.groups[ref] 
        obs_country = np.append(obs_country,obs_site_group.country)
        obs_lats = np.append(obs_lats,obs_site_group.latitude)
        obs_lons = np.append(obs_lons,obs_site_group.longitude)
        obs_alt = np.append(obs_alt,obs_site_group.altitude)
    obs_dates.append(obs_site_group.variables['date'][:])                                                                                                                    
    obs_times.append(obs_site_group.variables['time'][:])    

    for ref in base_refs:    
        if ref not in new_refs:
            bad_refs = np.append(bad_refs,ref)

print bad_refs

for ref in bad_refs:
    base_refs =  np.delete(base_refs,np.where(base_refs==ref))

print base_refs

obs_refs = np.copy(base_refs)
test = []
for i in all_refs:
    if i in obs_refs:
        test.append(True)
    else:
        test.append(False)
    
test = np.array(test)
obs_country = obs_country[test][:len(obs_refs)]
obs_lats = obs_lats[test][:len(obs_refs)]  
obs_lons = obs_lons[test][:len(obs_refs)]  
obs_alt = obs_alt[test][:len(obs_refs)]  

for c in range(len(start_years)):
    print c
    sy = start_years[c]
    ey = end_years[c]

    #load in periodic lsp data
    obs_period_grp = Dataset('obs_sig_periods_%s_%s.nc'%(sy,ey))

    for i in range(len(obs_refs)):
        obs_refs[i] = obs_refs[i].lower()

    odw = []
    osw = []

    for ref in obs_refs:
        site_group_obs = obs_period_grp.groups[ref]
    
        odw.append(site_group_obs.variables['daily_waveform'][:])
        osw.append(site_group_obs.variables['seasonal_waveform'][:])

    obs_daily_waveforms.append(odw)
    obs_seasonal_waveforms.append(osw)

np.save('obs_refs',obs_refs)
np.save('obs_lats',obs_lats)
np.save('obs_lons',obs_lons)
np.save('obs_country',obs_country)

for c in range(len(start_years)):
    print c
    sy = start_years[c]
    ey = end_years[c]

    obs_daily_waveform = np.array(obs_daily_waveforms[c])
    obs_seasonal_waveform = np.array(obs_seasonal_waveforms[c])

    #process obs dates and obs times to datetimes, then process pandas objects
    year_val = []
    month_val = []
    day_val = []
    hour_val = []
    minute_val = []

    obs_date = obs_dates[c]
    obs_time = obs_times[c]

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
    tags = modules.get_tags(obs_refs)

    #-----------------------------------
    #get area
    areas = ['ANT','AF','SE_US','S_US','W_US','N_US','NE_US','W_CAN','E_CAN','AL','S_EU','C_EU','NW_EU','N_EU','AS','ARC']
    if c == 0:
        if plot_type == 'd':
            obs_datetimes = obs_datetimes[:24]
        if plot_type == 's':
            obs_datetimes = obs_datetimes[:8766]

        obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')

    area_colors,area_boundaries,area_countries,area_labels = modules.area_dicts()
     
    #count = 0
    #for ax in axes.flat:
        #try:
            #area = areas[count]
        #except:
            #ax.axis('off')
            #continue

    count = 0
    for area in areas:
        print area
        area_grid = area_boundaries[area]
        area_country = area_countries[area]
        area_label = area_labels[area]

        cut_test = modules.area_cut(area,obs_lats,obs_lons,tags,area_grid,obs_country,area_country,obs_refs)

        print len(obs_daily_waveform[cut_test,:])
        obs_d_w = np.nanmean(obs_daily_waveform[cut_test,:],axis=0)
        obs_s_w = np.nanmean(obs_seasonal_waveform[cut_test,:],axis=0)

        if plot_type == 'd':
            ave_obs_waveform = obs_d_w
            max_ph = np.argmax(ave_obs_waveform)
            max_amp = (np.max(ave_obs_waveform) - np.min(ave_obs_waveform))/2.
            ave = np.average(ave_obs_waveform)
        if plot_type == 's':
            ave_obs_waveform = obs_s_w
            max_ph = np.argmax(ave_obs_waveform)/24.
            max_amp = (np.max(ave_obs_waveform) - np.min(ave_obs_waveform))/2.
            ave = np.average(ave_obs_waveform)

        #ax.plot(sy,max_ph,color = 'black',marker='x',markersize=10)
        plt.scatter(max_ph,sy,s=50,c=obsn_linspace[count],marker='o')
        
        #ax_dual = ax.twinx()
        #ax_dual.plot(sy,max_amp,color = 'red',marker='o',markersize=10)
        count+=1
        
        #if plot_type == 'd':
        #    ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        #if plot_type == 's':
        #    ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
    
    #count+=1

#h1, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[0]),marker='o',linestyle='None',markersize=8)
#h2, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[1]),marker='o',linestyle='None',markersize=8)
#h3, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[2]),marker='o',linestyle='None',markersize=8)
#h4, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[3]),marker='o',linestyle='None',markersize=8)
#h5, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[4]),marker='o',linestyle='None',markersize=8)
#h6, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[5]),marker='o',linestyle='None',markersize=8)
#h7, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[6]),marker='o',linestyle='None',markersize=8)
#h8, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[7]),marker='o',linestyle='None',markersize=8)
#h9, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[8]),marker='o',linestyle='None',markersize=8)
#h10, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[9]),marker='o',linestyle='None',markersize=8)
#h11, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[10]),marker='o',linestyle='None',markersize=8)
#h12, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[11]),marker='o',linestyle='None',markersize=8)
#h13, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[12]),marker='o',linestyle='None',markersize=8)
#h14, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[13]),marker='o',linestyle='None',markersize=8)
#h15, = plt.plot([1,1],color=plt.cm.cool(obsn_linspace[14]),marker='o',linestyle='None',markersize=8)

#plt.legend((h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15),['93-98','94-99','95-00','96-01','97-02','98-03','99-04','00-05','01-06','02-07','03-08','04-09','05-10','06-11','07-12'],loc='lower left',prop={'size':14},fancybox=True,ncol=2,markerscale=1,bbox_to_anchor=(-0.1,0))
#h1.set_visible(False)
#h2.set_visible(False)
#h3.set_visible(False)
#h4.set_visible(False)
#h5.set_visible(False)
#h6.set_visible(False)
#h7.set_visible(False)
#h8.set_visible(False)
#h9.set_visible(False)
#h10.set_visible(False)
#h11.set_visible(False)
#h12.set_visible(False)
#h13.set_visible(False)
#h14.set_visible(False)
#h15.set_visible(False)

plt.tight_layout(pad = 3.08)

plt.show()
