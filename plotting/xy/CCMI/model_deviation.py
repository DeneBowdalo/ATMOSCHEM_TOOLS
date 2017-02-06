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
from scipy import stats

#read in obs time series data
obs_fname = '/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_O3_2005_2010_H.nc'
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

tags = modules.get_tags(obs_refs)

#--------------------------------------------------------
#load in periodic lsp data
obs_period_grp = Dataset('../obs_SURFACE_H/obs_sig_periods.nc')
cesmcam_period_grp = Dataset('../CESMCAM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
cmam_period_grp = Dataset('../CMAM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
geosccm_period_grp = Dataset('../GEOSCCM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
geoschemv90103_period_grp = Dataset('../GEOSCHEM_SURFACE_v90103_2x2.5_GEOS5_H_*/model_sig_periods.nc')
geoschemv902_period_grp = Dataset('../GEOSCHEM_SURFACE_v902_2x2.5_GEOS5_H_*/model_sig_periods.nc')
gfdlam3_period_grp = Dataset('../GFDLAM3_SURFACE_*_*_*_H_*/model_sig_periods.nc')
gisse2r_period_grp = Dataset('../GISSE2R_SURFACE_*_*_*_H_*/model_sig_periods.nc')
mirocchem_period_grp = Dataset('../MIROCCHEM_SURFACE_*_*_*_H_*/model_sig_periods.nc') 

obs_daily_waveforms = []
obs_seasonal_waveforms = []
obs_full_waveforms = []
obs_ave = []
obs_daily_amp = []
obs_seasonal_amp = []
cesmcam_daily_waveforms = []
cesmcam_seasonal_waveforms = []
cesmcam_full_waveforms = []
cesmcam_ave = []
cesmcam_daily_amp = []
cesmcam_seasonal_amp = []
cmam_daily_waveforms = []
cmam_seasonal_waveforms = []
cmam_full_waveforms = []                                                                                                                                                     
cmam_ave = []
cmam_daily_amp = []
cmam_seasonal_amp = []
geosccm_daily_waveforms = []
geosccm_seasonal_waveforms = []
geosccm_full_waveforms = []
geosccm_ave = []
geosccm_daily_amp = []
geosccm_seasonal_amp = []
geoschemv90103_daily_waveforms = []
geoschemv90103_seasonal_waveforms = []
geoschemv90103_full_waveforms = [] 
geoschemv90103_ave = []
geoschemv90103_daily_amp = []
geoschemv90103_seasonal_amp = []
geoschemv902_daily_waveforms = []
geoschemv902_seasonal_waveforms = []
geoschemv902_full_waveforms = [] 
geoschemv902_ave = []
geoschemv902_daily_amp = []
geoschemv902_seasonal_amp = []
gfdlam3_daily_waveforms = []
gfdlam3_seasonal_waveforms = []
gfdlam3_full_waveforms = [] 
gfdlam3_ave = []
gfdlam3_daily_amp = []
gfdlam3_seasonal_amp = []
gisse2r_daily_waveforms = []
gisse2r_seasonal_waveforms = []
gisse2r_full_waveforms = [] 
gisse2r_ave = []
gisse2r_daily_amp = []
gisse2r_seasonal_amp = []
mirocchem_daily_waveforms = []
mirocchem_seasonal_waveforms = []
mirocchem_full_waveforms = []
mirocchem_ave = []
mirocchem_daily_amp = []
mirocchem_seasonal_amp = []

for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    site_group_cesmcam = cesmcam_period_grp.groups[ref]
    site_group_cmam = cmam_period_grp.groups[ref]    
    site_group_geosccm = geosccm_period_grp.groups[ref]
    site_group_geoschemv90103 = geoschemv90103_period_grp.groups[ref]
    site_group_geoschemv902 = geoschemv902_period_grp.groups[ref]
    site_group_gfdlam3 = gfdlam3_period_grp.groups[ref]
    site_group_gisse2r = gisse2r_period_grp.groups[ref]
    site_group_mirocchem = mirocchem_period_grp.groups[ref]

    cesmcam_daily_waveforms.append(site_group_cesmcam.variables['daily_waveform'][:])
    cesmcam_seasonal_waveforms.append(site_group_cesmcam.variables['seasonal_waveform'][:])
    cesmcam_full_waveforms.append(site_group_cesmcam.variables['all_waveform'][:])
    cesmcam_ave=np.append(cesmcam_ave,site_group_cesmcam.average)
    cesmcam_daily_amp = np.append(cesmcam_daily_amp,site_group_cesmcam.daily_amplitude)
    cesmcam_seasonal_amp = np.append(cesmcam_seasonal_amp,site_group_cesmcam.seasonal_amplitude)

    cmam_daily_waveforms.append(site_group_cmam.variables['daily_waveform'][:])
    cmam_seasonal_waveforms.append(site_group_cmam.variables['seasonal_waveform'][:])
    cmam_full_waveforms.append(site_group_cmam.variables['all_waveform'][:]) 
    cmam_ave=np.append(cmam_ave,site_group_cmam.average)    
    cmam_daily_amp = np.append(cmam_daily_amp,site_group_cmam.daily_amplitude)
    cmam_seasonal_amp = np.append(cmam_seasonal_amp,site_group_cmam.seasonal_amplitude)

    geosccm_daily_waveforms.append(site_group_geosccm.variables['daily_waveform'][:])                                                                                        
    geosccm_seasonal_waveforms.append(site_group_geosccm.variables['seasonal_waveform'][:])
    geosccm_full_waveforms.append(site_group_geosccm.variables['all_waveform'][:]) 
    geosccm_ave=np.append(geosccm_ave,site_group_geosccm.average)
    geosccm_daily_amp = np.append(geosccm_daily_amp,site_group_geosccm.daily_amplitude)
    geosccm_seasonal_amp = np.append(geosccm_seasonal_amp,site_group_geosccm.seasonal_amplitude)

    geoschemv90103_daily_waveforms.append(site_group_geoschemv90103.variables['daily_waveform'][:])                                                                          
    geoschemv90103_seasonal_waveforms.append(site_group_geoschemv90103.variables['seasonal_waveform'][:])
    geoschemv90103_full_waveforms.append(site_group_geoschemv90103.variables['all_waveform'][:]) 
    geoschemv90103_ave=np.append(geoschemv90103_ave,site_group_geoschemv90103.average)
    geoschemv90103_daily_amp = np.append(geoschemv90103_daily_amp,site_group_geoschemv90103.daily_amplitude)
    geoschemv90103_seasonal_amp = np.append(geoschemv90103_seasonal_amp,site_group_geoschemv90103.seasonal_amplitude)

    geoschemv902_daily_waveforms.append(site_group_geoschemv902.variables['daily_waveform'][:])                                                                              
    geoschemv902_seasonal_waveforms.append(site_group_geoschemv902.variables['seasonal_waveform'][:])
    geoschemv902_full_waveforms.append(site_group_geoschemv902.variables['all_waveform'][:]) 
    geoschemv902_ave=np.append(geoschemv902_ave,site_group_geoschemv902.average)
    geoschemv902_daily_amp = np.append(geoschemv902_daily_amp,site_group_geoschemv902.daily_amplitude)
    geoschemv902_seasonal_amp = np.append(geoschemv902_seasonal_amp,site_group_geoschemv902.seasonal_amplitude)

    gfdlam3_daily_waveforms.append(site_group_gfdlam3.variables['daily_waveform'][:])
    gfdlam3_seasonal_waveforms.append(site_group_gfdlam3.variables['seasonal_waveform'][:])
    gfdlam3_full_waveforms.append(site_group_gfdlam3.variables['all_waveform'][:])
    gfdlam3_ave=np.append(gfdlam3_ave,site_group_gfdlam3.average)
    gfdlam3_daily_amp = np.append(gfdlam3_daily_amp,site_group_gfdlam3.daily_amplitude)
    gfdlam3_seasonal_amp = np.append(gfdlam3_seasonal_amp,site_group_gfdlam3.seasonal_amplitude)

    gisse2r_daily_waveforms.append(site_group_gisse2r.variables['daily_waveform'][:])
    gisse2r_seasonal_waveforms.append(site_group_gisse2r.variables['seasonal_waveform'][:])
    gisse2r_full_waveforms.append(site_group_gisse2r.variables['all_waveform'][:])
    gisse2r_ave=np.append(gisse2r_ave,site_group_gisse2r.average)
    gisse2r_daily_amp = np.append(gisse2r_daily_amp,site_group_gisse2r.daily_amplitude)
    gisse2r_seasonal_amp = np.append(gisse2r_seasonal_amp,site_group_gisse2r.seasonal_amplitude)

    mirocchem_daily_waveforms.append(site_group_mirocchem.variables['daily_waveform'][:])                                                                                    
    mirocchem_seasonal_waveforms.append(site_group_mirocchem.variables['seasonal_waveform'][:])
    mirocchem_full_waveforms.append(site_group_mirocchem.variables['all_waveform'][:])
    mirocchem_ave=np.append(mirocchem_ave,site_group_mirocchem.average)
    mirocchem_daily_amp = np.append(mirocchem_daily_amp,site_group_mirocchem.daily_amplitude)
    mirocchem_seasonal_amp = np.append(mirocchem_seasonal_amp,site_group_mirocchem.seasonal_amplitude)

    obs_daily_waveforms.append(site_group_obs.variables['daily_waveform'][:])
    obs_seasonal_waveforms.append(site_group_obs.variables['seasonal_waveform'][:])
    obs_full_waveforms.append(site_group_obs.variables['all_waveform'][:])
    obs_ave = np.append(obs_ave,site_group_obs.average)
    obs_daily_amp = np.append(obs_daily_amp,site_group_obs.daily_amplitude)
    obs_seasonal_amp = np.append(obs_seasonal_amp,site_group_obs.seasonal_amplitude)

obs_daily_waveforms = np.array(obs_daily_waveforms)
obs_seasonal_waveforms = np.array(obs_seasonal_waveforms)
obs_full_waveforms = np.array(obs_full_waveforms)   
cesmcam_daily_waveforms = np.array(cesmcam_daily_waveforms)
cesmcam_seasonal_waveforms = np.array(cesmcam_seasonal_waveforms)
cesmcam_full_waveforms = np.array(cesmcam_full_waveforms)
cmam_daily_waveforms = np.array(cmam_daily_waveforms)                                                                                                                  
cmam_seasonal_waveforms = np.array(cesmcam_seasonal_waveforms)
cmam_full_waveforms = np.array(cmam_full_waveforms)
geosccm_daily_waveforms = np.array(geosccm_daily_waveforms)                                                                                                                  
geosccm_seasonal_waveforms = np.array(geosccm_seasonal_waveforms)
geosccm_full_waveforms = np.array(geosccm_full_waveforms)
geoschemv90103_daily_waveforms = np.array(geoschemv90103_daily_waveforms)                                                                                                    
geoschemv90103_seasonal_waveforms = np.array(geoschemv90103_seasonal_waveforms)                                                                                              
geoschemv90103_full_waveforms = np.array(geoschemv90103_full_waveforms)
geoschemv902_daily_waveforms = np.array(geoschemv902_daily_waveforms)                                                                                                        
geoschemv902_seasonal_waveforms = np.array(geoschemv902_seasonal_waveforms)                                                                                                  
geoschemv902_full_waveforms = np.array(geoschemv902_full_waveforms) 
gfdlam3_daily_waveforms = np.array(gfdlam3_daily_waveforms)                                                                                                                
gfdlam3_seasonal_waveforms = np.array(gfdlam3_seasonal_waveforms)                                                                                                          
gfdlam3_full_waveforms = np.array(gfdlam3_full_waveforms)
gisse2r_daily_waveforms = np.array(gisse2r_daily_waveforms)                                                                                                                
gisse2r_seasonal_waveforms = np.array(gisse2r_seasonal_waveforms)                                                                                                          
gisse2r_full_waveforms = np.array(gisse2r_full_waveforms)
mirocchem_daily_waveforms = np.array(mirocchem_daily_waveforms)                                                                                                              
mirocchem_seasonal_waveforms = np.array(mirocchem_seasonal_waveforms)                                                                                                        
mirocchem_full_waveforms = np.array(mirocchem_full_waveforms)

#calc r2 correlation
cesmcam_daily_r2 = []
cesmcam_seasonal_r2 = []
cesmcam_full_r2 = []
cmam_daily_r2 = []
cmam_seasonal_r2 = []
cmam_full_r2 = []
geosccm_daily_r2 = []
geosccm_seasonal_r2 = []
geosccm_full_r2 = []
geoschemv90103_daily_r2 = []
geoschemv90103_seasonal_r2 = []
geoschemv90103_full_r2 = []
geoschemv902_daily_r2 = []
geoschemv902_seasonal_r2 = []
geoschemv902_full_r2 = [] 
gfdlam3_daily_r2 = []
gfdlam3_seasonal_r2 = []
gfdlam3_full_r2 = [] 
gisse2r_daily_r2 = []
gisse2r_seasonal_r2 = []
gisse2r_full_r2 = []
mirocchem_seasonal_r2 = []
mirocchem_full_r2 = [] 
mirocchem_daily_r2 = []
mirocchem_seasonal_r2 = []
mirocchem_full_r2 = [] 

for i in range(len(obs_daily_waveforms)):
    slope, intercept, r_daily, p_value, std_err = stats.linregress(obs_daily_waveforms[i],cesmcam_daily_waveforms[i])
    slope, intercept, r_seasonal, p_value, std_err = stats.linregress(obs_seasonal_waveforms[i],cesmcam_seasonal_waveforms[i])
    slope, intercept, r_full, p_value, std_err = stats.linregress(obs_full_waveforms[i],cesmcam_full_waveforms[i])
    cesmcam_daily_r2 = np.append(cesmcam_daily_r2,r_daily**2)
    cesmcam_seasonal_r2 = np.append(cesmcam_seasonal_r2,r_seasonal**2)
    cesmcam_full_r2 = np.append(cesmcam_full_r2,r_full**2)

model_ave = [cesmcam_ave,cmam_ave,geosccm_ave,geoschemv90103_ave,geoschemv902_ave,gfdlam3_ave,gisse2r_ave,mirocchem_ave]
model_daily_amp = [cesmcam_daily_amp,cmam_daily_amp,geosccm_daily_amp,geoschemv90103_daily_amp,geoschemv902_daily_amp,gfdlam3_daily_amp,gisse2r_daily_amp,mirocchem_daily_amp]
model_seasonal_amp = [cesmcam_seasonal_amp,cmam_seasonal_amp,geosccm_seasonal_amp,geoschemv90103_seasonal_amp,geoschemv902_seasonal_amp,gfdlam3_seasonal_amp,gisse2r_seasonal_amp,mirocchem_seasonal_amp]

areas = ['ANT','OC','S_O','AF','SE_US','S_US','W_US','N_US','NE_US','W_CAN','E_CAN','S_EU','C_EU','NW_EU','N_EU','E_EU','AS','N_O','ARC']
area_boundaries = {'NE_US':[37,50,-90,-60],'SE_US':[25,37,-90,-60],'S_US':[25,37,-105,-90],'W_US':[25,50,-130,-105],'N_US':[37,60,-105,-90],'N_EU':[57,80,5,45],'C_EU':[47,57,5,20],'S_EU':[30,47,-10,40],'E_EU':[47,57,20,40],'NW_EU':[47,70,-15,5],'AS':[0,0,0,0],'N_O':[0,0,0,0],'S_O':[0,0,0,0],'OC':[0,0,0,0],'AF':[0,0,0,0],'E_CAN':[40,80,-95,-50],'W_CAN':[40,80,-150,-95],'ANT':[0,0,0,0],'ARC':[0,0,0,0]}
area_countries = {'NE_US':'United States','SE_US':'United States','S_US':'United States','W_US':'United States','N_US':'United States','C_EU':'EU','N_EU':'EU','S_EU':'EU','E_EU':'EU','NW_EU':'EU','AS':'AS','N_O':'O','S_O':'O','OC':'OC','AF':'AF','E_CAN':'Canada','W_CAN':'Canada','ANT':'ANT','ARC':'ARC'}
area_labels = {'NE_US':'NE US','SE_US':'SE US','S_US':'S US','W_US':'W US','N_US':'N US','N_EU':'N EU','NW_EU':'NW EU','C_EU':'C EU','E_EU':'E EU','S_EU':'S EU','AS':'Asia','N_O':'NH Oceanic','S_O':'SH Oceanic','OC':'Oceania','AF':'Africa','E_CAN':'East Canada','W_CAN':'West Canada','ANT':'Antarctica','ARC':'Arctic'}

models = ['CESMCAM','CMAM','GEOSCCM','GEOSCHEMV90103','GEOSCHEMV902','GFDL','GISS','MIROCCHEM']
facecolors = ['#7570b3','#771b9e','#1b429e','#9e1b42','#429e1b','#1b9e77','#9e771b','#ffc0cb']
whiskercolors = ['#7570b3','#7570b3','#771b9e','#771b9e','#1b429e','#1b429e','#9e1b42','#9e1b42','#429e1b','#429e1b','#1b9e77','#1b9e77','#9e771b','#9e771b']

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)] 

offset_type = raw_input('ave, amp, ph_max, ph_min or r?\n')
if offset_type != 'ave':
    waveform_type = raw_input('\nd, s or full?\n')

#set up plot
fig =plt.figure(figsize=(19,13))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(111)

plt.axhline(y=0,color='red',linestyle = '--')


count = 0

for i in range(len(models)):
    
    area_offset = []
    
    for area in areas:

        print area

        r, g, b = tableau20[count]  
        tableau20[count] = (r / 255., g / 255., b / 255.) 

        area_grid = area_boundaries[area]
        area_country = area_countries[area]
        area_label = area_labels[area]

        if (area == 'C_EU') or (area == 'E_EU') or (area == 'S_EU') or (area == 'NW_EU') or (area == 'N_EU'):
            cut_test = (obs_lats >= area_grid[0]) & (obs_lats < area_grid[1]) & (obs_lons >= area_grid[2]) & (obs_lons < area_grid[3]) & (tags == area_country)
        elif (area == 'AS') or (area == 'OC') or (area == 'AF') or (area == 'ANT') or (area == 'ARC'):
            cut_test = tags == area_country
        elif (area == 'N_O'):
            cut_test = (tags == area_country) & (obs_lats >= 0)
        elif (area == 'S_O'):
            cut_test = (tags == area_country) & (obs_lats < 0)
        else:
            cut_test = (obs_lats >= area_grid[0]) & (obs_lats < area_grid[1]) & (obs_lons >= area_grid[2]) & (obs_lons < area_grid[3]) & (obs_country == area_country)

        if offset_type == 'ave':
            area_offset.append(np.average(model_ave[i][cut_test] - obs_ave[cut_test]))

        if offset_type == 'amp':
            if waveform_type == 'd':
                area_offset.append(np.average(model_daily_amp[i][cut_test] - obs_daily_amp[cut_test]))
            elif waveform_type == 's':
                area_offset.append(np.average(model_seasonal_amp[i][cut_test] - obs_seasonal_amp[cut_test]))


        #if offset_type = 'ph_max':
            #area_offset = 

        #if offset_type = 'ph_min':
            #area_offset = 

        #if offset_type = 'r':
            #area_offset = 

    ax.plot(np.arange(len(area_offset)),area_offset,marker='o',linestyle='None',markersize=10,color=facecolors[i])

    count+=1

ax.set_xlim([-1,19])
ax.set_xticks(np.arange(0,19))
labels = []
for i in areas:
    labels.append(area_labels[i])
ax.set_xticklabels(labels,rotation=45)

h1, = ax.plot([1,1],color=facecolors[0],marker='o',linestyle='None',markersize=10)
h2, = ax.plot([1,1],color=facecolors[1],marker='o',linestyle='None',markersize=10)
h3, = ax.plot([1,1],color=facecolors[2],marker='o',linestyle='None',markersize=10)
h4, = ax.plot([1,1],color=facecolors[3],marker='o',linestyle='None',markersize=10)
h5, = ax.plot([1,1],color=facecolors[4],marker='o',linestyle='None',markersize=10)
h6, = ax.plot([1,1],color=facecolors[5],marker='o',linestyle='None',markersize=10)
h7, = ax.plot([1,1],color=facecolors[6],marker='o',linestyle='None',markersize=10)
h8, = ax.plot([1,1],color=facecolors[7],marker='o',linestyle='None',markersize=10)

plt.legend((h1,h2,h3,h4,h5,h6,h7,h8),models)
h1.set_visible(False)
h2.set_visible(False)
h3.set_visible(False)
h4.set_visible(False)
h5.set_visible(False)
h6.set_visible(False)
h7.set_visible(False)
h8.set_visible(False)

plt.show()
