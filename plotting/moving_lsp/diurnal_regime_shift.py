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
paths = present_dir.split("/")
species = paths[-1]

start_year = 1993
end_year = 2012
nyears = 5

start_years = range(start_year,(end_year-nyears)+1)
end_years = range(start_year+nyears,end_year+1)

obsn = len(start_years)
obsn_linspace = np.linspace(0, 1, obsn)

plot_type = raw_input('\nd or s?\n')

fig, axes = plt.subplots(nrows=5, ncols=5,figsize=(19,12))
fig.patch.set_facecolor('white')

#find constant refs in all files
start_refs = []
end_refs = []

#load in periodic lsp data
start_period_grp = Dataset('obs_sig_periods_1993_1998.nc')
end_period_grp = Dataset('obs_sig_periods_2007_2012.nc')

start_refs_dict = start_period_grp.groups
end_refs_dict = end_period_grp.groups

for i in start_refs_dict.keys():
    i = i.encode('ascii')
    start_refs = np.append(start_refs,i)

for i in range(len(start_refs)):
    start_refs[i] = start_refs[i].lower()

for i in end_refs_dict.keys():
    i = i.encode('ascii')
    end_refs = np.append(end_refs,i)

for i in range(len(end_refs)):
    end_refs[i] = end_refs[i].lower()

obs_refs = list(set(start_refs) & set(end_refs))

site_diff = []

for ref in obs_refs:
    start_group = start_period_grp.groups[ref]
    end_group = end_period_grp.groups[ref]    

    start_diff = np.max(start_group.variables['daily_waveform'][11:18]) / np.max(start_group.variables['daily_waveform'][:6])
    end_diff = np.max(end_group.variables['daily_waveform'][11:18]) / np.max(end_group.variables['daily_waveform'][:6])
    if start_diff < 1.1:
        site_diff.append(np.abs(end_diff - start_diff))
    else:
        site_diff.append(0)

    if ref == 'tkb':
        print np.max(start_group.variables['daily_waveform'][11:18]), np.max(start_group.variables['daily_waveform'][:6])
        print np.max(end_group.variables['daily_waveform'][11:18]), np.max(end_group.variables['daily_waveform'][:6])
        print start_diff
        print end_diff

site_diff, obs_refs = zip(*sorted(zip(site_diff, obs_refs)))

print site_diff[-40:]
print obs_refs[-40:]

1+'a'


#get observational location tags 
tags = modules.get_tags(obs_refs)

areas = ['ANT','S_O','OC','AF','SA','NE_NA','CE_NA','SE_NA','S_NA','SW_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','NE_AS','SE_AS','C_AS','S_AS','N_O','AL','ARC']

np.save('constant_obs_refs',obs_refs)
np.save('constant_obs_lats',obs_lats)
np.save('constant_obs_lons',obs_lons)
np.save('constant_obs_tags',tags)

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

    #-----------------------------------
    #get area
    if c == 0:
        if plot_type == 'd':
            obs_datetimes = obs_datetimes[:24]
        if plot_type == 's':
            obs_datetimes = obs_datetimes[:8766]

        obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')

    area_boundaries,area_tags,area_labels = modules.area_dicts()    
 
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

        obs_d_w = np.nanmean(obs_daily_waveform[cut_test,:],axis=0)
        obs_s_w = np.nanmean(obs_seasonal_waveform[cut_test,:],axis=0)

        if plot_type == 'd':
            ave_obs_waveform = obs_d_w
            ax.plot_date(obs_time_pd.to_pydatetime(),ave_obs_waveform,color = plt.cm.cool(obsn_linspace[c]),linestyle='--',linewidth=1,markeredgecolor='None',markersize=4)
        if plot_type == 's':
            ave_obs_waveform = obs_s_w
            ax.plot_date(obs_time_pd.to_pydatetime(),ave_obs_waveform,color = plt.cm.cool(obsn_linspace[c]),linestyle='_',linewidth=1,markeredgecolor='None',markersize=1)

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

h1, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[0]),marker='o',linestyle='None',markersize=8)
h2, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[1]),marker='o',linestyle='None',markersize=8)
h3, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[2]),marker='o',linestyle='None',markersize=8)
h4, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[3]),marker='o',linestyle='None',markersize=8)
h5, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[4]),marker='o',linestyle='None',markersize=8)
h6, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[5]),marker='o',linestyle='None',markersize=8)
h7, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[6]),marker='o',linestyle='None',markersize=8)
h8, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[7]),marker='o',linestyle='None',markersize=8)
h9, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[8]),marker='o',linestyle='None',markersize=8)
h10, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[9]),marker='o',linestyle='None',markersize=8)
h11, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[10]),marker='o',linestyle='None',markersize=8)
h12, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[11]),marker='o',linestyle='None',markersize=8)
h13, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[12]),marker='o',linestyle='None',markersize=8)
h14, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[13]),marker='o',linestyle='None',markersize=8)
h15, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[14]),marker='o',linestyle='None',markersize=8)

plt.legend((h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15),['93-98','94-99','95-00','96-01','97-02','98-03','99-04','00-05','01-06','02-07','03-08','04-09','05-10','06-11','07-12'],loc='lower left',prop={'size':10},fancybox=True,ncol=2,markerscale=1,bbox_to_anchor=(-0.2,-0.3))
h1.set_visible(False)
h2.set_visible(False)
h3.set_visible(False)
h4.set_visible(False)
h5.set_visible(False)
h6.set_visible(False)
h7.set_visible(False)
h8.set_visible(False)
h9.set_visible(False)
h10.set_visible(False)
h11.set_visible(False)
h12.set_visible(False)
h13.set_visible(False)
h14.set_visible(False)
h15.set_visible(False)

plt.tight_layout(pad = 3.08)

plt.show()
