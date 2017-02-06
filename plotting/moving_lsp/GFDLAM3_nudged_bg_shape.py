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

start_year = 1983
end_year = 2008
nyears = 5

start_years = range(start_year,(end_year-nyears)+1)
end_years = range(start_year+nyears,end_year+1)

obsn = len(start_years)
obsn_linspace = np.linspace(0, 1, 25)

plot_type = raw_input('\nd or s?\n')

fig, axes = plt.subplots(nrows=4, ncols=5,figsize=(19,12))
fig.patch.set_facecolor('white')

#obs_refs = ['cpt','abfr38011','abfr38012','abpt07001','150010006','780200001','vii423','ogasawara','ang','bmw','cvo','mnm','rpb','smo','bhd','cgo','lau','sja']
obs_refs = ['cpt','abfr38011','abfr38012','150010006','780200001','vii423','ogasawara','bmw','cvo','mnm','rpb','smo','bhd','cgo','lau','sja']
#obs_refs = ['060670011', 'alh157', 'bvl130', '450210002', '470930021']
#obs_refs = ['alh157', '450210002', '470930021', 'cad150', 'tkb']
#obs_refs = ['ckt136', 'cpt', '270750005', 'no0039r', 'abcz0lsou']
#obs_refs = ['bmw', 'alt', 'no0042g', 'mnm', 'cpt']
#obs_refs = ['vin140', 'vpi120', 'wes', 'wst109', 'zgt']
#obs_refs = ['brw', 'nmy', 'abdebw031', 'bmw', 'alt']
#obs_refs = ['ie0031r', 'mnm', 'fi0022r', 'no0015r', '470090101']
obs_daily_waveforms = []
obs_seasonal_waveforms = []

for c in range(len(start_years)):
    print c
    sy = start_years[c]
    ey = end_years[c]
    odw = []
    osw = []

    #load in periodic lsp data
    obs_period_grp = Dataset('GFDLAM3_nudged_sig_periods_%s_%s.nc'%(sy,ey))
    a = obs_period_grp.groups

    for ref in obs_refs:
        try:
            site_group_obs = obs_period_grp.groups[ref]
            odw.append(site_group_obs.variables['daily_waveform'][:])
            osw.append(site_group_obs.variables['seasonal_waveform'][:])
        except:
            odw.append(np.array([np.NaN]*24))
            osw.append(np.array([np.NaN]*8766))

    obs_daily_waveforms.append(odw)
    obs_seasonal_waveforms.append(osw)

#get observational location tags 
tags = modules.get_tags(obs_refs)

for c in range(len(start_years)):
    print c
    sy = start_years[c]
    ey = end_years[c]

    obs_daily_waveform = np.array(obs_daily_waveforms[c])
    obs_seasonal_waveform = np.array(obs_seasonal_waveforms[c])
    #-----------------------------------
    #get area
    if c == 0:
        if plot_type == 'd':
            start_dt = datetime.datetime(2008,1,1,0,0,0)
            end_dt = datetime.datetime(2008,1,1,23,0,0)
            obs_time_pd = pd.date_range(start = start_dt,end = end_dt, freq = 'H')
        if plot_type == 's':
            start_dt = datetime.datetime(2008,1,1,0,0,0)
            end_dt = datetime.datetime(2008,12,31,23,0,0)
            obs_time_pd = pd.date_range(start = start_dt,end = end_dt, freq = 'H')[:8766]    

 
    count = 0
    for ax in axes.flat:
        if count <  len(obs_refs):
            obs_d_w = obs_daily_waveform[count,:]
            obs_s_w = obs_seasonal_waveform[count,:]        
            print obs_d_w            

            if plot_type == 'd':
                if np.isnan(obs_d_w).any() == False:
                    ax.plot_date(obs_time_pd.to_pydatetime(),obs_d_w,color = plt.cm.cool(obsn_linspace[c]),linestyle='--',linewidth=1,markeredgecolor='None',markersize=4)
            if plot_type == 's':
                if np.isnan(obs_s_w).any() == False:
                    ax.plot_date(obs_time_pd.to_pydatetime(),obs_s_w,color = plt.cm.cool(obsn_linspace[c]),linestyle='_',linewidth=1,markeredgecolor='None',markersize=1)

            ax.text(0.87, 0.91, obs_refs[count], ha='center', va='center',transform=ax.transAxes,fontsize=15)
            if plot_type == 'd':
                ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
            if plot_type == 's':
                ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
        else:
           ax.axis('off') 

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
h16, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[15]),marker='o',linestyle='None',markersize=8)
h17, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[16]),marker='o',linestyle='None',markersize=8)
h18, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[17]),marker='o',linestyle='None',markersize=8)
h19, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[18]),marker='o',linestyle='None',markersize=8)
h20, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[19]),marker='o',linestyle='None',markersize=8)
#h21, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[20]),marker='o',linestyle='None',markersize=8)
#h22, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[21]),marker='o',linestyle='None',markersize=8)
#h23, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[22]),marker='o',linestyle='None',markersize=8)
#h24, = ax.plot([1,1],color=plt.cm.cool(obsn_linspace[23]),marker='o',linestyle='None',markersize=8)

plt.legend((h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19,h20),['83-88','84-89','85-90','86-91','87-92','88-93','89-94','90-95','91-96','92-97','93-98','94-99','95-00','96-01','97-02','98-03','99-04','00-05','01-06','02-07','03-08','04-09','05-10','06-11'],loc='lower left',prop={'size':14},fancybox=True,ncol=4,markerscale=1,bbox_to_anchor=(-2.5,0.1))
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
h16.set_visible(False)
h17.set_visible(False)
h18.set_visible(False)
h19.set_visible(False)
h20.set_visible(False)
#h21.set_visible(False)
#h22.set_visible(False)
#h23.set_visible(False)
#h24.set_visible(False)

plt.tight_layout(pad = 3.08)

plt.show()
