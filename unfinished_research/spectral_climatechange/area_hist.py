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
import matplotlib.dates as dates
import seaborn as sns

models = ['CESMCAM', 'CMAM', 'GFDLAM3','GISSE2R', 'MIROCCHEM', 'MOCAGE', 'UMCAM']

#get color set for models
cmaplist = plt.cm.get_cmap('Set1')
array = np.linspace(0,1,len(models))
cmaplist = [cmaplist(i) for i in array]

#read in obs
fobs = '/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_O3_2000_2003_H_HP.nc'
obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_raw_class,obs_anthrome_class,obs_gap_inds = modules.read_obs_all(fobs,'O3',2000,2003)

#obs_std_var[obs_std_var < 0] = np.NaN

tags = modules.get_tags(obs_refs)
areas = ['ANT','S_O','AF','SW_NA','NW_NA','NE_NA','CE_NA','SE_NA','C_NA','S_EU','SW_EU','CW_EU','CS_EU','C_EU','E_EU','NW_EU','N_EU','SE_AS','NE_AS','N_O','AL','ARC']
area_boundaries,area_tags,area_labels = modules.area_dicts()

for i in range(len(models)):

    fig, axes = plt.subplots(nrows=6, ncols=4,figsize=(19,13))
    fig.patch.set_facecolor('white')

    model = models[i]
    print model
    
    if (model == 'CESMCAM') or (model == 'CMAM') or (model == 'MIROCCHEM') or (model == 'MOCAGE'):
        year2000s = 2000
        year2000e = 2003
        year2100s = 2100
        year2100e = 2103

    if (model == 'GFDLAM3'):
        year2000s = 2001
        year2000e = 2004
        year2100s = 2101
        year2100e = 2104

    if (model == 'GISSE2R'):
        year2000s = 2000
        year2000e = 2003
        year2100s = 2102
        year2100e = 2105

    if (model == 'UMCAM'):
        year2000s = 2000
        year2000e = 2003
        year2100s = 2097
        year2100e = 2100


    #read in 2000 model period data
    f2000 = '/work/home/db876/plotting_tools/model_files/%s_SURFACE_2000_2012_*_*_*_H_*.nc'%(model)

    #read in 2100 model period data
    f2100 = '/work/home/db876/plotting_tools/model_files/%s_SURFACE_2095_2111_*_*_*_H_rcp85.nc'%(model)

    #read in 2100 model fixed emissions period data
    f2100e = '/work/home/db876/plotting_tools/model_files/%s_SURFACE_2095_2111_*_*_*_H_rcp85em2000.nc'%(model)

    f2000raw_time,f2000ref_time,f2000datetime_time,f2000std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(f2000,'O3',year2000s,year2000e)
    if model != 'GISSE2R':
        f2100raw_time,f2100ref_time,f2100datetime_time,f2100std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(f2100,'O3',year2100s,year2100e)
    if model != 'CMAM':
        f2100eraw_time,f2100eref_time,f2100edatetime_time,f2100estd_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(f2100e,'O3',year2100s,year2100e)

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
        na,na,spring_inds,na,na,summer_inds,na,na,autumn_inds,na,na,winter_inds = modules.cut_season(obs_datetime_time,obs_ref_time,obs_std_var[0,:])
        data_obs = obs_std_var[cut_test,:]
        #data_obs = data_obs[:,winter_inds]
        data_obs = np.array([item for sublist in data_obs for item in sublist])
        testinv = data_obs > 0                                                                                                                                                                                     
        data_obs = data_obs[testinv]
        bins, edges = np.histogram(data_obs, 34, normed=1)
        ax.step(edges[:-1], bins, where='post',color='black')
    
        na, na, model_indices,lat_indices,lon_indices = modules.obs_model_gridboxes(lat_e,lon_e,obs_lats,obs_lons)
        data_f2000 = f2000std_var[:,lat_indices,lon_indices]
        na,na,spring_inds,na,na,summer_inds,na,na,autumn_inds,na,na,winter_inds = modules.cut_season(f2000datetime_time,f2000ref_time,data_f2000[:,0]) 
        data_f2000 = data_f2000[:,cut_test] 
        #data_f2000 = data_f2000[winter_inds,:]
        data_f2000 = np.array([item for sublist in data_f2000 for item in sublist])       
        testinv = data_f2000 > 0                                                                                                                                                                                     
        data_f2000 = data_f2000[testinv]
        bins, edges = np.histogram(data_f2000, 34, normed=1)
        ax.step(edges[:-1], bins, where='post',color='blue')   
        if model != 'GISSE2R':
            data_f2100 = f2100std_var[:,lat_indices,lon_indices]
            na,na,spring_inds,na,na,summer_inds,na,na,autumn_inds,na,na,winter_inds = modules.cut_season(f2100datetime_time,f2100ref_time,data_f2100[:,0])
            data_f2100 = data_f2100[:,cut_test] 
            #data_f2100 = data_f2100[winter_inds,:]
            data_f2100 = np.array([item for sublist in data_f2100 for item in sublist])
            testinv = data_f2100 > 0                                                                                                                                                                                   
            data_f2100 = data_f2100[testinv]
            bins, edges = np.histogram(data_f2100, 34, normed=1)                                                                                                                                                      
            ax.step(edges[:-1], bins, where='post',color='red')
        if model != 'CMAM':
            data_f2100e = f2100estd_var[:,lat_indices,lon_indices]
            na,na,spring_inds,na,na,summer_inds,na,na,autumn_inds,na,na,winter_inds = modules.cut_season(f2100edatetime_time,f2100eref_time,data_f2100e[:,0])
            data_f2100e = data_f2100e[:,cut_test]
            #data_f2100e = data_f2100e[winter_inds,:] 
            data_f2100e = np.array([item for sublist in data_f2100e for item in sublist])
            testinv = data_f2100e > 0                                                                                                                                                                                   
            data_f2100e = data_f2100e[testinv]
            bins, edges = np.histogram(data_f2100e, 34, normed=1)
            ax.step(edges[:-1], bins, where='post',color='green')        

        ax.annotate(area,xy=(0.01,0.88),xycoords='axes fraction',alpha=5) 

        count+=1
        
    h1, = ax.plot([1,1],color='black',marker='o',linestyle='None',markersize=10)
    h2, = ax.plot([1,1],color='blue',marker='o',linestyle='None',markersize=10)
    h3, = ax.plot([1,1],color='red',marker='o',linestyle='None',markersize=10)
    h4, = ax.plot([1,1],color='green',marker='o',linestyle='None',markersize=10)

    plt.legend((h1,h2,h3,h4),['Observations','%s 2000'%(model),'%s 2100 RCP8.5'%(model),'%s 2100 RCP8.5 2000 Emissions'%(model)],loc='lower left',prop={'size':13},fancybox=True,ncol=1,markerscale=1,bbox_to_anchor=(-1.0,-0.1))
    h1.set_visible(False)
    h2.set_visible(False)
    h3.set_visible(False)
    h4.set_visible(False)

    plt.savefig('plots/hist_%s.png'%(model))

    #plt.show()
