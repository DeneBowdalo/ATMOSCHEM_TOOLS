import modules
import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm
from cmath import *
import colorsys
from scipy import signal
from netCDF4 import Dataset
import datetime
import pandas as pd
import matplotlib.lines as mlines
import operator
from collections import OrderedDict

obs_fname = '/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_O3_2005_2010_H_PERIODIC.nc'
model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_2005_2012_v90103_2x2.5_GEOS5_H_*.nc'
species = 'O3'
start_year = 2005
end_year = 2010

other_species = 'GMAO_TEMP'

obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)
#model_raw_time,model_ref_time,model_datetime_time,no_model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,'NO',start_year,end_year)
#model_raw_time,model_ref_time,model_datetime_time,no2_model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,'NO2',start_year,end_year)
#nox_model_std_var = no_model_std_var + no2_model_std_var
model_raw_time,model_ref_time,model_datetime_time,nox_model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,other_species,start_year,end_year)


#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

model_time_pd = pd.date_range(start = model_datetime_time[0],end = model_datetime_time[-1], freq = 'H')
key = [model_time_pd.month]
nox_by_site = []
nox_ave = []
#cut all lat/lon indices
for i in range(len(obs_lons)):
    obs_lat = obs_lats[i]
    obs_lon = obs_lons[i]
    lat_i,lon_i = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)

    nox_model_var_pd = pd.Series(nox_model_std_var[:,lat_i,lon_i], index=model_time_pd)
    nox_ave.append(np.average(nox_model_std_var[:,lat_i,lon_i]))
    group=nox_model_var_pd.groupby(key)
    nox_monthly_mean = group.mean()
    
    nox_by_site.append(nox_monthly_mean)
 
nox_by_site = np.array(nox_by_site)
nox_ave = np.array(nox_ave)

#--------------------------------------------------------
#load in periodic lsp data

periodic_obs_fname = '/work/home/db876/xy/O3/2005_2012/2005_2010/obs_SURFACE_H/obs_sig_periods.nc'
periodic_model_fname = '/work/home/db876/xy/O3/2005_2012/2005_2010/GEOSCHEM_SURFACE_v90103_2x2.5_GEOS5_H_*/model_sig_periods.nc'

obs_daily_mag,obs_daily_mag_spring,obs_daily_mag_summer,obs_daily_mag_autumn,obs_daily_mag_winter,obs_daily_phase\
,obs_daily_phase_spring,obs_daily_phase_spring,obs_daily_phase_summer,obs_daily_phase_autumn,obs_daily_phase_winter\
,obs_daily_phase_min,obs_daily_phase_min_spring,obs_daily_phase_min_summer,obs_daily_phase_min_autumn,obs_daily_phase_min_winter\
,obs_ave,obs_ave_spring,obs_ave_summer,obs_ave_autumn,obs_ave_winter,obs_ave_day,obs_ave_night,obs_seasonal_mag,obs_seasonal_mag_day\
,obs_seasonal_mag_night,obs_seasonal_phase,obs_seasonal_phase_day,obs_seasonal_phase_night,obs_seasonal_phase_min\
,obs_seasonal_phase_min_day,obs_seasonal_phase_min_night,obs_daily_ff,obs_daily_ff_spring,obs_daily_ff_summer,obs_daily_ff_autumn\
,obs_daily_ff_winter,obs_seasonal_ff,obs_seasonal_ff_day,obs_seasonal_ff_night,obs_daily_waveform,obs_daily_waveform_spring\
,obs_daily_waveform_summer,obs_daily_waveform_autumn,obs_daily_waveform_winter,obs_seasonal_waveform,obs_seasonal_waveform_day\
,obs_seasonal_waveform_night,obs_full_waveform,obs_periodic_variance_daily,obs_periodic_variance_daily_spring,obs_periodic_variance_daily_summer\
,obs_periodic_variance_daily_autumn,obs_periodic_variance_daily_winter,obs_periodic_variance_seasonal,obs_periodic_variance_seasonal_day\
,obs_periodic_variance_seasonal_night,obs_periodic_variance_all,model_daily_mag,model_daily_mag_spring,model_daily_mag_summer,model_daily_mag_autumn\
,model_daily_mag_winter,model_daily_phase,model_daily_phase_spring,model_daily_phase_spring,model_daily_phase_summer\
,model_daily_phase_autumn,model_daily_phase_winter,model_daily_phase_min,model_daily_phase_min_spring,model_daily_phase_min_summer\
,model_daily_phase_min_autumn,model_daily_phase_min_winter,model_ave,model_ave_spring,model_ave_summer,model_ave_autumn,model_ave_winter\
,model_ave_day,model_ave_night,model_seasonal_mag,model_seasonal_mag_day,model_seasonal_phase,model_seasonal_phase_day,model_seasonal_phase_min\
,model_seasonal_phase_min_day,model_daily_ff,model_daily_ff_spring,model_daily_ff_summer,model_daily_ff_autumn,model_daily_ff_winter,model_seasonal_ff\
,model_seasonal_ff_day,model_seasonal_ff_night,model_daily_waveform,model_daily_waveform_spring,model_daily_waveform_summer,model_daily_waveform_autumn\
,model_daily_waveform_winter,model_seasonal_waveform,model_seasonal_waveform_day,model_seasonal_waveform_night,model_full_waveform,model_periodic_variance_daily\
,model_periodic_variance_daily_spring,model_periodic_variance_daily_summer,model_periodic_variance_daily_autumn,model_periodic_variance_daily_winter,model_periodic_variance_seasonal\
,model_periodic_variance_seasonal_day,model_periodic_variance_seasonal_night,model_periodic_variance_all = modules.get_periodic_specific(periodic_obs_fname,periodic_model_fname,obs_refs)
        
        
#------------------------------------------
tags = modules.get_tags(obs_refs)


loc_colors = {'ANT':'magenta','ARC':'purple','AF':'orange','AS':'blue','EU':'red','OC':'black','O':'yellow','NA':'green','SA':'pink'}
all_locs = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic','NA':'North America','SA':'South America'}
order_dict = {'ANT':8,'ARC':9,'AF':6,'AS':3,'EU':2,'OC':5,'O':4,'NA':1,'SA':7}

order_tags = []
for tag in tags:
    n = order_dict[tag]
    order_tags.append(n)

tags_argsort = np.argsort(order_tags)

tags = tags[tags_argsort]
obs_seasonal_mag = obs_seasonal_mag[tags_argsort]
obs_seasonal_phase = obs_seasonal_phase[tags_argsort]
model_seasonal_mag = model_seasonal_mag[tags_argsort]
model_seasonal_phase = model_seasonal_phase[tags_argsort]
nox_by_site = nox_by_site[tags_argsort]
nox_ave = nox_ave[tags_argsort]
obs_lons = obs_lons[tags_argsort]
obs_lats = obs_lats[tags_argsort]

seasonal_amplitude_diff = model_seasonal_mag - obs_seasonal_mag
seasonal_phase_diff = model_seasonal_phase - obs_seasonal_phase
for d in range(len(seasonal_phase_diff)):
    if seasonal_phase_diff[d] >= 6:
        seasonal_phase_diff[d] = -6.+(seasonal_phase_diff[d]-6.)
    elif seasonal_phase_diff[d] < -6:
        seasonal_phase_diff[d] = 6.-(np.abs(seasonal_phase_diff[d])-6.)

phase_ind = np.round(model_seasonal_phase)

nox_month_ave = []

#cut arrays
#test = seasonal_phase_diff < 1
#obs_lons = obs_lons[test]
#obs_lats = obs_lats[test]
#tags = tags[test]
#seasonal_amplitude_diff = seasonal_amplitude_diff[test]
#seasonal_phase_diff =  seasonal_phase_diff[test]
#nox_ave = nox_ave[test]
#nox_by_site = nox_by_site[test]
#phase_ind = phase_ind[test]

#in month of model seasonal phase get average nox in that month
for j in range(len(nox_by_site)):
    nox_month_ave.append(nox_by_site[j][phase_ind[j]])
    
tag_colors = [loc_colors[i] for i in tags]
tag_locs = [all_locs[i] for i in tags]

#set up daily plot
fig =plt.figure(figsize=(19,14))
fig.patch.set_facecolor('white')

ax1 = plt.subplot2grid((12,12), (0,0),rowspan=6,colspan = 6)
ax2 = plt.subplot2grid((12,12), (0,6),rowspan=6,colspan = 6)

ax1.axhline(0,color='black',zorder=0,linewidth=3,alpha=0.6)
ax2.axhline(0,color='black',zorder=0,linewidth=3,alpha=0.6)

for i in range(len(seasonal_phase_diff)):                                                                                   
    ax1.scatter(nox_ave[i],seasonal_amplitude_diff[i],c=tag_colors[i],s=30,label=tag_locs[i])
    
for i in range(len(seasonal_phase_diff)):                                                                                   
    ax2.scatter(nox_ave[i],seasonal_phase_diff[i],c=tag_colors[i],s=30,label=tag_locs[i])

min_var = np.min(nox_ave)
max_var = np.max(nox_ave)

ax1.set_xlabel(r'GEOS-Chem $%s$ 2005-2010 Average (ppb)'%(other_species),fontsize=14)
ax1.set_ylabel(r'GEOS-Chem $O_3$ Seasonal Amplitude Bias (ppb)',fontsize=14)
ax1.set_xlim(min_var,max_var)

ax2.set_xlabel(r'GEOS-Chem $%s$ 2005-2010 Average (ppb)'%(other_species),fontsize=14)
ax2.set_ylabel(r'GEOS-Chem $O_3$ Seasonal Phase Bias (Months)',fontsize=14)
ax2.set_xlim(min_var,max_var)

ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.tick_params(axis='both', which='minor', labelsize=14)
ax2.tick_params(axis='both', which='major', labelsize=14)
ax2.tick_params(axis='both', which='minor', labelsize=14)

ax3 = plt.subplot2grid((12,12), (6,0),rowspan=3,colspan=3)
ax4 = plt.subplot2grid((12,12), (6,3),rowspan=3,colspan=3)
ax5 = plt.subplot2grid((12,12), (9,0),rowspan=3,colspan=3)
ax6 = plt.subplot2grid((12,12), (9,3),rowspan=3,colspan=3)

latlower_setup = [20,30,20,lat_e[0]]
latupper_setup = [80,72,55,lat_e[-1]]
lonwest_setup = [-170,-15,115,lon_e[0]]
loneast_setup = [-50,35,155,lon_e[-1]]
label_out = ['NA','EU','AS','zz']
label = ['NA','EU','AS','ROW']
size_dict = {'NA':30,'EU':30,'AS':40,'AF':60,'O':80,'OC':80,'SA':80,'ARC':80,'ANT':80}

obs_lons = np.array(obs_lons)
obs_lats = np.array(obs_lats)
tags = np.array(tags)

ax_list = [ax3,ax4,ax5,ax6]
count = 0 

for ax in ax_list:
    #setup basemap projection
    m = Basemap(projection='cyl',llcrnrlat=latlower_setup[count],urcrnrlat=latupper_setup[count],llcrnrlon=lonwest_setup[count],urcrnrlon=loneast_setup[count],resolution='c',ax = ax)

    m.drawcoastlines()
    m.drawmapboundary()
    
    if count == 3:
        test = []
        other_tags = ['AF','ANT','ARC','OC','O','SA']
        for i in tags:
            if i in other_tags:
                test.append(True)
            else:
                test.append(False)
        test = np.array(test)
    else:
        current_tag = label_out[count] 
        test = tags == current_tag
    
    
    current_z = nox_ave[test]
    current_lons = obs_lons[test]
    current_lats = obs_lats[test]
    
    X,Y = m(current_lons,current_lats)

    if count == 0:
        m_size= 15
    if count == 1:
        m_size = 15
    if count == 2:
        m_size = 100
    if count == 3:
        m_size = 100

    for i in range(len(current_lons)):
        all = m.scatter(X[i],Y[i],c=current_z[i], s=m_size, vmin = min_var,vmax = max_var, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.rainbow,zorder=10)

    count+=1

plt.tight_layout(pad = 1.08)

#cb
cbar_ax = fig.add_axes([0.51, 0.05, 0.03, 0.4])
#t = np.arange(min_var,max_var,1)
cb = fig.colorbar(all, orientation = 'vertical', cax=cbar_ax)
cb.ax.tick_params(labelsize=18)
cb.set_label(r'$%s$ 2005-2010'%(other_species)+'\nAverage (ppb)', fontsize = 20)

plt.annotate('NA',xy=(-15,0.65), xycoords='axes fraction', fontsize=30)
plt.annotate('AS',xy=(-14.55,0.34), xycoords='axes fraction', fontsize=30)
plt.annotate('EU',xy=(-6.4,0.97), xycoords='axes fraction', fontsize=30)
plt.annotate('ROW',xy=(-6.9,0.01), xycoords='axes fraction', fontsize=30)

plt.annotate('a',xy=(-14.82,2.25), xycoords='axes fraction', fontsize=36)
plt.annotate('b',xy=(1.4,2.25), xycoords='axes fraction', fontsize=36)
plt.annotate('c',xy=(-14.75,0.99), xycoords='axes fraction', fontsize=36)


all_handles, all_labels = ax1.get_legend_handles_labels()
hl = sorted(zip(all_handles, all_labels), key=operator.itemgetter(1))
all_handles, all_labels = zip(*hl)
by_label = OrderedDict(zip(all_labels, all_handles))
ax2.legend(by_label.values(),by_label.keys(), loc='lower right',prop={'size':13},fancybox=True,ncol=2,markerscale=2)#,bbox_to_anchor=(-0.7,-0.45))

plt.show()
