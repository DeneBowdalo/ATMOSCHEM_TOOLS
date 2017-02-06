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
import psutil
import gc
import matplotlib.gridspec as gridspec

present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-4]

start_year = 2009
end_year = 2011

#read in obs ts data
obs_fname = '/work/home/db876/observations/surface/%s/process/GLOBAL_SURFACE_%s_2009_2011_H_HP.nc'%(species,species)
obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_var,obs_lats,obs_lons,obs_alt,obs_groups,obs_raw_class,obs_anthrome_class,gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)

#read in std model data
model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_2009_2011_v1001_4x5_GEOS5_H_STD.nc'
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

#--------------------------------------------------------
#load in periodic lsp data


obs_grp = Dataset('../obs_SURFACE_H/LSP_stats.nc')

alt_model_dirs_a = ['GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ANMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ANMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ANMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ANMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ANMVOC0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ANMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ANMVOC2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ANMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ANMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ANMVOC0.5',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_STD','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ANMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ANMVOC4.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ANMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ANMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ANMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ANMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ANMVOC0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ANMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ANMVOC2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ANMVOC4.0']

alt_model_dirs_b = ['GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25BNMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25BNMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25BNMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25BNMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5BNMVOC0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5BNMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5BNMVOC2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5BNMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0BNMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0BNMVOC0.5',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_STD','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0BNMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0BNMVOC4.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0BNMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0BNMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0BNMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0BNMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0BNMVOC0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0BNMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0BNMVOC2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0BNMVOC4.0']
                    
alt_model_dirs_c = ['GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25DRYDEPO30.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25DRYDEPO30.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25DRYDEPO32.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25DRYDEPO34.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5DRYDEPO30.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5DRYDEPO30.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5DRYDEPO32.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5DRYDEPO34.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0DRYDEPO30.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0DRYDEPO30.5',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_STD','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0DRYDEPO32.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0DRYDEPO34.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0DRYDEPO30.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0DRYDEPO30.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0DRYDEPO32.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0DRYDEPO34.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0DRYDEPO30.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0DRYDEPO30.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0DRYDEPO32.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0DRYDEPO34.0'] 

#if 'rig' in obs_refs:
#    obs_refs[obs_refs.index('rig')] = 'rig_photo'

obs_seasonal_waveforms = obs_grp.variables['seasonal_waveform'][:]

#-----------------------------------
#get area
areas = ['SW_NA','NW_NA','NE_NA','CE_NA','SE_NA','C_NA','S_EU','SW_EU','CW_EU','CS_EU','C_EU','E_EU','NW_EU','N_EU','NE_AS','SE_AS']

area_boundaries,area_tags,area_labels = modules.area_dicts()

diff_wf_s_a_spring = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_summer = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_autumn = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_winter = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_b_spring = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_summer = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_autumn = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_winter = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_c_spring = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_summer = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_autumn = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_winter = np.empty((len(areas),len(alt_model_dirs_c)))

#cut vals into seasons
start = datetime.datetime(year = 2008, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = 2009, month = 1, day = 1, hour = 0, minute = 0)
ref_date_dt = pd.date_range(start,end,freq='H')[:-19]
months = np.array([d.strftime('%m') for d in ref_date_dt]).astype('int')

valid_inds_winter = (months < 3) | (months ==12)
valid_inds_spring = (months >=3) & (months <6)
valid_inds_summer = (months >= 6) & (months <9)
valid_inds_autumn = (months >= 9) & (months <12)

for m in range(len(alt_model_dirs_a)):

    print 'point 1'

    print m
    print '../%s/LSP_stats.nc'%(alt_model_dirs_a[m])
    print '../%s/LSP_stats.nc'%(alt_model_dirs_b[m])
    print '../%s/LSP_stats.nc'%(alt_model_dirs_c[m])
    
    alt_model_grp_a = Dataset('../%s/LSP_stats.nc'%(alt_model_dirs_a[m]))
    alt_model_grp_b = Dataset('../%s/LSP_stats.nc'%(alt_model_dirs_b[m]))
    alt_model_grp_c = Dataset('../%s/LSP_stats.nc'%(alt_model_dirs_c[m]))
        
    alt_model_seasonal_waveforms_a = alt_model_grp_a.variables['seasonal_waveform'][:]
    alt_model_seasonal_waveforms_b = alt_model_grp_b.variables['seasonal_waveform'][:]
    alt_model_seasonal_waveforms_c = alt_model_grp_c.variables['seasonal_waveform'][:]

    day = np.arange(0,24,1)
    year = np.linspace(0,12,8766,endpoint=False)
 
    count = 0
    for a in range(len(areas)):
    
        area = areas[a]
        
        print area

        area_grid = area_boundaries[area]
        area_tag = area_tags[area]
        area_label = area_labels[area]

        cut_test = modules.area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag)
        if np.all(cut_test == False):
            diff_wf_s_a_spring[a,m] = np.NaN
            diff_wf_s_a_summer[a,m] = np.NaN
            diff_wf_s_a_autumn[a,m] = np.NaN
            diff_wf_s_a_winter[a,m] = np.NaN
            diff_wf_s_b_spring[a,m] = np.NaN
            diff_wf_s_b_summer[a,m] = np.NaN
            diff_wf_s_b_autumn[a,m] = np.NaN
            diff_wf_s_b_winter[a,m] = np.NaN
            diff_wf_s_c_spring[a,m] = np.NaN
            diff_wf_s_c_summer[a,m] = np.NaN
            diff_wf_s_c_autumn[a,m] = np.NaN
            diff_wf_s_c_winter[a,m] = np.NaN
        else:
            
            obs_sites = obs_seasonal_waveforms[cut_test,:]
            obs_s_w_spring = np.nanmean(obs_sites[:,valid_inds_spring],axis=0)
            obs_s_w_summer = np.nanmean(obs_sites[:,valid_inds_summer],axis=0)
            obs_s_w_autumn = np.nanmean(obs_sites[:,valid_inds_autumn],axis=0)
            obs_s_w_winter = np.nanmean(obs_sites[:,valid_inds_winter],axis=0)       
     
            model_sites_a = alt_model_seasonal_waveforms_a[cut_test,:]
            model_sites_b = alt_model_seasonal_waveforms_b[cut_test,:]
            model_sites_c = alt_model_seasonal_waveforms_c[cut_test,:]
            alt_model_s_w_a_spring = np.nanmean(model_sites_a[:,valid_inds_spring],axis=0)
            alt_model_s_w_a_summer = np.nanmean(model_sites_a[:,valid_inds_summer],axis=0)
            alt_model_s_w_a_autumn = np.nanmean(model_sites_a[:,valid_inds_autumn],axis=0)
            alt_model_s_w_a_winter = np.nanmean(model_sites_a[:,valid_inds_winter],axis=0)
            alt_model_s_w_b_spring = np.nanmean(model_sites_b[:,valid_inds_spring],axis=0)
            alt_model_s_w_b_summer = np.nanmean(model_sites_b[:,valid_inds_summer],axis=0)
            alt_model_s_w_b_autumn = np.nanmean(model_sites_b[:,valid_inds_autumn],axis=0)
            alt_model_s_w_b_winter = np.nanmean(model_sites_b[:,valid_inds_winter],axis=0)
            alt_model_s_w_c_spring = np.nanmean(model_sites_c[:,valid_inds_spring],axis=0)
            alt_model_s_w_c_summer = np.nanmean(model_sites_c[:,valid_inds_summer],axis=0)
            alt_model_s_w_c_autumn = np.nanmean(model_sites_c[:,valid_inds_autumn],axis=0)
            alt_model_s_w_c_winter = np.nanmean(model_sites_c[:,valid_inds_winter],axis=0)

            diff_wf_s_a_spring[a,m] = np.sum(np.abs(alt_model_s_w_a_spring - obs_s_w_spring))            
            diff_wf_s_a_summer[a,m] = np.sum(np.abs(alt_model_s_w_a_summer - obs_s_w_summer))
            diff_wf_s_a_autumn[a,m] = np.sum(np.abs(alt_model_s_w_a_autumn - obs_s_w_autumn))
            diff_wf_s_a_winter[a,m] = np.sum(np.abs(alt_model_s_w_a_winter - obs_s_w_winter))

            diff_wf_s_b_spring[a,m] = np.sum(np.abs(alt_model_s_w_b_spring - obs_s_w_spring))            
            diff_wf_s_b_summer[a,m] = np.sum(np.abs(alt_model_s_w_b_summer - obs_s_w_summer))
            diff_wf_s_b_autumn[a,m] = np.sum(np.abs(alt_model_s_w_b_autumn - obs_s_w_autumn))
            diff_wf_s_b_winter[a,m] = np.sum(np.abs(alt_model_s_w_b_winter - obs_s_w_winter))

            diff_wf_s_c_spring[a,m] = np.sum(np.abs(alt_model_s_w_c_spring - obs_s_w_spring))            
            diff_wf_s_c_summer[a,m] = np.sum(np.abs(alt_model_s_w_c_summer - obs_s_w_summer))
            diff_wf_s_c_autumn[a,m] = np.sum(np.abs(alt_model_s_w_c_autumn - obs_s_w_autumn))
            diff_wf_s_c_winter[a,m] = np.sum(np.abs(alt_model_s_w_c_winter - obs_s_w_winter))
        
        count+=1
            
    #remove unneeded variables
    try:
        del ave_obs_param
        del ave_alt_model_param_a
        del ave_alt_model_param_b
        del ave_alt_model_param_c
        del alt_model_seasonal_waveforms_a
        del alt_model_seasonal_waveforms_b 
        del alt_model_seasonal_waveforms_c 
        del obs_s_w
        del alt_model_s_w_a 
        del alt_model_s_w_b 
        del alt_model_s_w_c 
        del area_grid
        del area_tag
        del area_label
        del cut_test
    except:
        pass
    gc.collect()   
 
    print '\n'

spring_inds = np.array([0,2,4,6,16,18,20,22,32,34,36,38,48,50,52,54])
summer_inds = np.array([1,3,5,7,17,19,21,23,33,35,37,39,49,51,53,55])
autumn_inds = np.array([8,10,12,14,24,26,28,30,40,42,44,46,56,58,60,62])
winter_inds = np.array([9,11,13,15,25,27,29,31,41,43,45,47,57,59,61,63])

plotter = 'T'
while plotter == 'T':
    fig = plt.figure(figsize = (14,13))
    fig.patch.set_facecolor('white')

    gs1 = gridspec.GridSpec(2, 2)
    gs1.update(top=0.99,bottom=0.79,left=0.04,right=0.27,wspace=0,hspace=0)
    ax1 = plt.subplot(gs1[0, 0])
    ax2 = plt.subplot(gs1[0, 1])
    ax3 = plt.subplot(gs1[1, 0])
    ax4 = plt.subplot(gs1[1, 1])
    gs2 = gridspec.GridSpec(2, 2)
    gs2.update(top=0.99,bottom=0.79,left=0.28, right=0.51,wspace=0,hspace=0)
    ax5 = plt.subplot(gs2[0, 0])
    ax6 = plt.subplot(gs2[0, 1])
    ax7 = plt.subplot(gs2[1, 0])
    ax8 = plt.subplot(gs2[1, 1])
    gs3 = gridspec.GridSpec(2, 2)
    gs3.update(top=0.99,bottom=0.79,left=0.52, right=0.75,wspace=0,hspace=0)
    ax9 = plt.subplot(gs3[0, 0])
    ax10 = plt.subplot(gs3[0, 1])
    ax11 = plt.subplot(gs3[1, 0])
    ax12 = plt.subplot(gs3[1, 1]) 
    gs4 = gridspec.GridSpec(2, 2)
    gs4.update(top=0.99,bottom=0.79,left=0.76, right=0.99,wspace=0,hspace=0)
    ax13 = plt.subplot(gs4[0, 0])
    ax14 = plt.subplot(gs4[0, 1])
    ax15 = plt.subplot(gs4[1, 0])
    ax16 = plt.subplot(gs4[1, 1])
    gs5 = gridspec.GridSpec(2, 2)
    gs5.update(top=0.78,bottom=0.58,left=0.04, right=0.27,wspace=0,hspace=0)
    ax17 = plt.subplot(gs5[0, 0])
    ax18 = plt.subplot(gs5[0, 1])
    ax19 = plt.subplot(gs5[1, 0])
    ax20 = plt.subplot(gs5[1, 1])
    gs6 = gridspec.GridSpec(2, 2)
    gs6.update(top=0.78,bottom=0.58,left=0.28, right=0.51,wspace=0,hspace=0)
    ax21 = plt.subplot(gs6[0, 0])
    ax22 = plt.subplot(gs6[0, 1])
    ax23 = plt.subplot(gs6[1, 0])
    ax24 = plt.subplot(gs6[1, 1])
    gs7 = gridspec.GridSpec(2, 2)
    gs7.update(top=0.78,bottom=0.58,left=0.52, right=0.75,wspace=0,hspace=0)
    ax25 = plt.subplot(gs7[0, 0])
    ax26 = plt.subplot(gs7[0, 1])
    ax27 = plt.subplot(gs7[1, 0])
    ax28 = plt.subplot(gs7[1, 1])
    gs8 = gridspec.GridSpec(2, 2)
    gs8.update(top=0.78,bottom=0.58,left=0.76, right=0.99,wspace=0,hspace=0)
    ax29 = plt.subplot(gs8[0, 0])
    ax30 = plt.subplot(gs8[0, 1])
    ax31 = plt.subplot(gs8[1, 0])
    ax32 = plt.subplot(gs8[1, 1])
    gs9 = gridspec.GridSpec(2, 2)
    gs9.update(top=0.57,bottom=0.37,left=0.04, right=0.27,wspace=0,hspace=0)
    ax33 = plt.subplot(gs9[0, 0])
    ax34 = plt.subplot(gs9[0, 1])
    ax35 = plt.subplot(gs9[1, 0])
    ax36 = plt.subplot(gs9[1, 1])
    gs10 = gridspec.GridSpec(2, 2)    
    gs10.update(top=0.57,bottom=0.37,left=0.28, right=0.51,wspace=0,hspace=0)
    ax37 = plt.subplot(gs10[0, 0])
    ax38 = plt.subplot(gs10[0, 1])
    ax39 = plt.subplot(gs10[1, 0])
    ax40 = plt.subplot(gs10[1, 1])
    gs11 = gridspec.GridSpec(2, 2)
    gs11.update(top=0.57,bottom=0.37,left=0.52, right=0.75,wspace=0,hspace=0)
    ax41 = plt.subplot(gs11[0, 0])
    ax42 = plt.subplot(gs11[0, 1])
    ax43 = plt.subplot(gs11[1, 0])
    ax44 = plt.subplot(gs11[1, 1])
    gs12 = gridspec.GridSpec(2, 2)
    gs12.update(top=0.57,bottom=0.37,left=0.76, right=0.99,wspace=0,hspace=0)
    ax45 = plt.subplot(gs12[0, 0])
    ax46 = plt.subplot(gs12[0, 1])
    ax47 = plt.subplot(gs12[1, 0])
    ax48 = plt.subplot(gs12[1, 1])
    gs13 = gridspec.GridSpec(2, 2)
    gs13.update(top=0.36,bottom=0.16,left=0.04, right=0.27,wspace=0,hspace=0)
    ax49 = plt.subplot(gs13[0, 0])
    ax50 = plt.subplot(gs13[0, 1])
    ax51 = plt.subplot(gs13[1, 0])
    ax52 = plt.subplot(gs13[1, 1])
    gs14 = gridspec.GridSpec(2, 2)
    gs14.update(top=0.36,bottom=0.16,left=0.28, right=0.51,wspace=0,hspace=0)
    ax53 = plt.subplot(gs14[0, 0])
    ax54 = plt.subplot(gs14[0, 1])
    ax55 = plt.subplot(gs14[1, 0])
    ax56 = plt.subplot(gs14[1, 1])
    gs15 = gridspec.GridSpec(2, 2) 
    gs15.update(top=0.36,bottom=0.16,left=0.52, right=0.75,wspace=0,hspace=0)
    ax57 = plt.subplot(gs15[0, 0])
    ax58 = plt.subplot(gs15[0, 1])
    ax59 = plt.subplot(gs15[1, 0])
    ax60 = plt.subplot(gs15[1, 1])
    gs16 = gridspec.GridSpec(2, 2)
    gs16.update(top=0.36,bottom=0.16,left=0.76, right=0.99,wspace=0,hspace=0)                                                                                    
    ax61 = plt.subplot(gs16[0, 0])
    ax62 = plt.subplot(gs16[0, 1])
    ax63 = plt.subplot(gs16[1, 0])
    ax64 = plt.subplot(gs16[1, 1])

    axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14,ax15,ax16,ax17,ax18,ax19,ax20,ax21,ax22,ax23,ax24,ax25,ax26,ax27,ax28,ax29,ax30,ax31,ax32,ax33,ax34,ax35,ax36,ax37,ax38,ax39,ax40,ax41,ax42,ax43,ax44,ax45,ax46,ax47,ax48,ax49,ax50,ax51,ax52,ax53,ax54,ax55,ax56,ax57,ax58,ax59,ax60,ax61,ax62,ax63,ax64]

    set_type = raw_input('\nANMVOC ,BNMVOC or DRYDEPO3?\n')

    plot_grid = raw_input('\nbox or contour?\n')

    area_count = 0
    ax_count = 0
    for area in areas:
        print area
        #ax.axis('off')
        
        area_label = area_labels[area]
    
        #ANMVOC cut
        area_grid_a_spring = diff_wf_s_a_spring[area_count,:]
        area_grid_a_summer = diff_wf_s_a_summer[area_count,:]
        area_grid_a_autumn = diff_wf_s_a_autumn[area_count,:]
        area_grid_a_winter = diff_wf_s_a_winter[area_count,:]
        #BNMVOC cut
        area_grid_b_spring = diff_wf_s_b_spring[area_count,:]
        area_grid_b_summer = diff_wf_s_b_summer[area_count,:]
        area_grid_b_autumn = diff_wf_s_b_autumn[area_count,:]
        area_grid_b_winter = diff_wf_s_b_winter[area_count,:]
        #drydepo3 cut
        area_grid_c_spring = diff_wf_s_c_spring[area_count,:]
        area_grid_c_summer = diff_wf_s_c_summer[area_count,:]
        area_grid_c_autumn = diff_wf_s_c_autumn[area_count,:]
        area_grid_c_winter = diff_wf_s_c_winter[area_count,:]
        
        if set_type == 'ANMVOC':
            area_grid_spring = area_grid_a_spring
            area_grid_summer = area_grid_a_summer
            area_grid_autumn = area_grid_a_autumn
            area_grid_winter = area_grid_a_winter
        elif set_type == 'BNMVOC':
            area_grid_spring = area_grid_b_spring
            area_grid_summer = area_grid_b_summer
            area_grid_autumn = area_grid_b_autumn
            area_grid_winter = area_grid_b_winter
        elif set_type == 'DRYDEPO3':
            area_grid_spring = area_grid_c_spring
            area_grid_summer = area_grid_c_summer
            area_grid_autumn = area_grid_c_autumn
            area_grid_winter = area_grid_c_winter
            
        #set min and max        
        if species == 'O3':
            minval = 0
            maxval = 100000
            t = [0,100000]
            t_str = ['0','+100000']
            cb_label = 'Integrated Seasonal Bias'
            cts = np.linspace(0,100000,40)
        elif species == 'NO':                                                                                                                        
            minval = 0
            maxval = 100000
            t = [0,100000]
            t_str = ['0','+100000']
            cb_label = 'Integrated Seasonal Bias'
            cts = np.linspace(0,100000,40)
        elif species == 'NO2-MOLYBDENUM':                                                                                                                        
            minval = 0
            maxval = 300000
            t = [0,300000]
            t_str = ['0','+300000']
            cb_label = 'Integrated Seasonal Bias'
            cts = np.linspace(0,300000,40)
        elif species == 'CO':
            minval = 0
            maxval = 2000000
            t = [0,2000000]
            t_str = ['0','+2000000']
            cb_label = 'Integrated Seasonal Bias'
            cts = np.linspace(0,2000000,40)

        cmap = matplotlib.cm.jet

        area_grid_spring = np.reshape(area_grid_spring,(5,5))    
        area_grid_summer = np.reshape(area_grid_summer,(5,5))
        area_grid_autumn = np.reshape(area_grid_autumn,(5,5))
        area_grid_winter = np.reshape(area_grid_winter,(5,5))
        masked_array_spring = np.ma.array(area_grid_spring, mask=np.isnan(area_grid_spring))
        masked_array_summer = np.ma.array(area_grid_summer, mask=np.isnan(area_grid_summer))
        masked_array_autumn = np.ma.array(area_grid_autumn, mask=np.isnan(area_grid_autumn))
        masked_array_winter = np.ma.array(area_grid_winter, mask=np.isnan(area_grid_winter))
        cmap.set_bad('w',1.)
        if plot_grid == 'box':
            pl = axes[ax_count].pcolor(masked_array_spring,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+1].pcolor(masked_array_summer,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+2].pcolor(masked_array_autumn,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+3].pcolor(masked_array_winter,vmin = minval,vmax=maxval,cmap =cmap)
        elif plot_grid == 'contour':
            pl = axes[ax_count].contourf(masked_array_spring,cts,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+1].contourf(masked_array_summer,cts,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+2].contourf(masked_array_autumn,cts,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+3].contourf(masked_array_winter,cts,vmin = minval,vmax=maxval,cmap =cmap)
        
        axes[ax_count].set_xticks([1,2,3,4])
        axes[ax_count+1].set_xticks([1,2,3,4])
        axes[ax_count+2].set_xticks([1,2,3,4])
        axes[ax_count+3].set_xticks([1,2,3,4])
        axes[ax_count].set_xticklabels(['','','',''])
        axes[ax_count+1].set_xticklabels(['','','',''])
        axes[ax_count+2].set_xticklabels(['','','',''])
        axes[ax_count+3].set_xticklabels(['','','',''])

        axes[ax_count].set_yticks([1,2,3,4])
        axes[ax_count+1].set_yticks([1,2,3,4])
        axes[ax_count+2].set_yticks([1,2,3,4])
        axes[ax_count+3].set_yticks([1,2,3,4])
        axes[ax_count].set_yticklabels(['','','',''])
        axes[ax_count+1].set_yticklabels(['','','','']) 
        axes[ax_count+2].set_yticklabels(['','','','']) 
        axes[ax_count+3].set_yticklabels(['','','','']) 
        #ax.set_yticklabels(['0.25','0.5','1.0','2.0','4.0'])
        
        #axes[ax_count].axes.get_xaxis().set_visible(False)
        #axes[ax_count].axes.get_yaxis().set_visible(False)
        #axes[ax_count+1].axes.get_xaxis().set_visible(False)
        #axes[ax_count+1].axes.get_yaxis().set_visible(False)
        #axes[ax_count+2].axes.get_xaxis().set_visible(False)
        #axes[ax_count+2].axes.get_yaxis().set_visible(False)
        #axes[ax_count+3].axes.get_xaxis().set_visible(False)
        #axes[ax_count+3].axes.get_yaxis().set_visible(False)

        axes[ax_count].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        axes[ax_count+1].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        axes[ax_count+2].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        axes[ax_count+3].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)

        axes[ax_count+1].text(0.81, 0.85, area_label, ha='center', va='center',transform=axes[ax_count+1].transAxes,fontsize=15)
    
        area_count+=1
        ax_count+=4
        
    #plt.tight_layout(pad = 1.5)
    fig.subplots_adjust(bottom=0.08)
    fig.subplots_adjust(left=0.10)

    fig.text(0.5, 0.12, set_type, ha='center',fontsize=30)
    fig.text(0.01, 0.5, 'NOx', va='center', rotation='vertical',fontsize=30)

    cbar_ax = fig.add_axes([0.58, 0.07, 0.35, 0.06])
    cb = fig.colorbar(pl,orientation='horizontal',cax=cbar_ax,ticks=t)
    cb.set_ticklabels(t_str)
    cb.ax.tick_params(labelsize=20)
    cb.set_label(cb_label,fontsize=20)
    cb.set_clim(minval, maxval)

    plt.show()

    plotter = raw_input('\nAnother plot? T or F?\n')
