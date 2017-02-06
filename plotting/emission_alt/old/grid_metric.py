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

present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-4]

start_year = 2009
end_year = 2011

#read in obs ts data
obs_fname = '/work/home/db876/observations/surface/%s/process/GLOBAL_SURFACE_%s_2009_2011_H_HP.nc'%(species,species)
obs_refs,raw_time,ref_time,datetime_time,std_var,obs_lats,obs_lons,obs_alt,obs_groups,obs_raw_class,obs_anthrome_class,gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)

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

obs_daily_amp = obs_grp.variables['diurnal_amplitude'][:]
obs_daily_ph = obs_grp.variables['diurnal_phase'][:] 
obs_seasonal_amp = obs_grp.variables['seasonal_amplitude'][:] 
obs_seasonal_ph = obs_grp.variables['seasonal_phase'][:] 
obs_ave = obs_grp.variables['average'][:] 
obs_daily_waveforms = obs_grp.variables['diurnal_average_waveform'][:]
obs_seasonal_waveforms = obs_grp.variables['seasonal_waveform'][:]

#-----------------------------------
#get area
areas = ['SW_NA','NW_NA','NE_NA','CE_NA','SE_NA','C_NA','S_EU','SW_EU','CW_EU','CS_EU','C_EU','E_EU','NW_EU','N_EU','NE_AS','SE_AS']

area_boundaries,area_tags,area_labels = modules.area_dicts()

diff_amp_d_a = np.empty((len(areas),len(alt_model_dirs_a)))
diff_ph_d_a = np.empty((len(areas),len(alt_model_dirs_a)))
diff_ave_a = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_d_a = np.empty((len(areas),len(alt_model_dirs_a)))
diff_amp_s_a = np.empty((len(areas),len(alt_model_dirs_a)))
diff_ph_s_a = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a = np.empty((len(areas),len(alt_model_dirs_a)))

diff_amp_d_b = np.empty((len(areas),len(alt_model_dirs_b)))
diff_ph_d_b = np.empty((len(areas),len(alt_model_dirs_b)))
diff_ave_b = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_d_b = np.empty((len(areas),len(alt_model_dirs_b)))                                                                                                                        
diff_amp_s_b = np.empty((len(areas),len(alt_model_dirs_b)))
diff_ph_s_b = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b = np.empty((len(areas),len(alt_model_dirs_b)))

diff_amp_d_c = np.empty((len(areas),len(alt_model_dirs_c)))
diff_ph_d_c = np.empty((len(areas),len(alt_model_dirs_c)))
diff_ave_c = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_d_c = np.empty((len(areas),len(alt_model_dirs_c)))                                                                                                                        
diff_amp_s_c = np.empty((len(areas),len(alt_model_dirs_c)))
diff_ph_s_c = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c = np.empty((len(areas),len(alt_model_dirs_c)))

for m in range(len(alt_model_dirs_a)):

    print 'point 1'

    print m
    print '../%s/LSP_stats.nc'%(alt_model_dirs_a[m])
    print '../%s/LSP_stats.nc'%(alt_model_dirs_b[m])
    print '../%s/LSP_stats.nc'%(alt_model_dirs_c[m])
    
    alt_model_grp_a = Dataset('../%s/LSP_stats.nc'%(alt_model_dirs_a[m]))
    alt_model_grp_b = Dataset('../%s/LSP_stats.nc'%(alt_model_dirs_b[m]))
    alt_model_grp_c = Dataset('../%s/LSP_stats.nc'%(alt_model_dirs_c[m]))
        
    alt_model_daily_amp_a = alt_model_grp_a.variables['diurnal_amplitude'][:]
    alt_model_daily_ph_a = alt_model_grp_a.variables['diurnal_phase'][:]
    alt_model_seasonal_amp_a = alt_model_grp_a.variables['seasonal_amplitude'][:]
    alt_model_seasonal_ph_a = alt_model_grp_a.variables['seasonal_phase'][:]
    alt_model_ave_a = alt_model_grp_a.variables['average'][:]
    alt_model_daily_waveforms_a = alt_model_grp_a.variables['diurnal_average_waveform'][:]
    alt_model_seasonal_waveforms_a = alt_model_grp_a.variables['seasonal_waveform'][:]
            
    alt_model_daily_amp_b = alt_model_grp_b.variables['diurnal_amplitude'][:]
    alt_model_daily_ph_b = alt_model_grp_b.variables['diurnal_phase'][:]
    alt_model_seasonal_amp_b = alt_model_grp_b.variables['seasonal_amplitude'][:]
    alt_model_seasonal_ph_b = alt_model_grp_b.variables['seasonal_phase'][:]
    alt_model_ave_b = alt_model_grp_b.variables['average'][:]
    alt_model_daily_waveforms_b = alt_model_grp_b.variables['diurnal_average_waveform'][:]
    alt_model_seasonal_waveforms_b = alt_model_grp_b.variables['seasonal_waveform'][:]
 
    alt_model_daily_amp_c = alt_model_grp_c.variables['diurnal_amplitude'][:]
    alt_model_daily_ph_c = alt_model_grp_c.variables['diurnal_phase'][:]
    alt_model_seasonal_amp_c = alt_model_grp_c.variables['seasonal_amplitude'][:]
    alt_model_seasonal_ph_c = alt_model_grp_c.variables['seasonal_phase'][:]
    alt_model_ave_c = alt_model_grp_c.variables['average'][:]
    alt_model_daily_waveforms_c = alt_model_grp_c.variables['diurnal_average_waveform'][:]
    alt_model_seasonal_waveforms_c = alt_model_grp_c.variables['seasonal_waveform'][:]

    print 'point 2'

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
            diff_ave_a[a,m] = np.NaN
            diff_ave_b[a,m] = np.NaN
            diff_ave_c[a,m] = np.NaN         
            diff_amp_d_a[a,m] = np.NaN
            diff_amp_d_b[a,m] = np.NaN
            diff_amp_d_c[a,m] = np.NaN
            diff_ph_d_a[a,m] = np.NaN
            diff_ph_d_b[a,m] = np.NaN
            diff_ph_d_c[a,m] = np.NaN
            diff_wf_d_a[a,m] = np.NaN
            diff_wf_d_b[a,m] = np.NaN
            diff_wf_d_c[a,m] = np.NaN
            diff_amp_s_a[a,m] = np.NaN
            diff_amp_s_b[a,m] = np.NaN
            diff_amp_s_c[a,m] = np.NaN
            diff_ph_s_a[a,m] = np.NaN
            diff_ph_s_b[a,m] = np.NaN
            diff_ph_s_c[a,m] = np.NaN
            diff_wf_s_a[a,m] = np.NaN
            diff_wf_s_b[a,m] = np.NaN
            diff_wf_s_c[a,m] = np.NaN
    
        else:
            obs_d_w = np.nanmean(obs_daily_waveforms[cut_test,:],axis=0)
            obs_d_amp = (np.nanmax(obs_d_w)-np.nanmin(obs_d_w))/2.
            obs_d_ph = day[np.nanargmax(obs_d_w)]
            obs_s_w = np.nanmean(obs_seasonal_waveforms[cut_test,:],axis=0)
            obs_s_amp = (np.nanmax(obs_s_w)-np.nanmin(obs_s_w))/2.
            obs_s_ph = year[np.nanargmax(obs_s_w)]
            obs_a = np.nanmean(obs_s_w)

            alt_model_d_w_a = np.nanmean(alt_model_daily_waveforms_a[cut_test,:],axis=0)
            alt_model_d_amp_a = (np.nanmax(alt_model_d_w_a)-np.nanmin(alt_model_d_w_a))/2.
            alt_model_d_ph_a = day[np.nanargmax(alt_model_d_w_a)]
            alt_model_s_w_a = np.nanmean(alt_model_seasonal_waveforms_a[cut_test,:],axis=0)
            alt_model_s_amp_a = (np.nanmax(alt_model_s_w_a)-np.nanmin(alt_model_s_w_a))/2.
            alt_model_s_ph_a = year[np.nanargmax(alt_model_s_w_a)]        
            alt_model_a_a = np.nanmean(alt_model_s_w_a)

            alt_model_d_w_b = np.nanmean(alt_model_daily_waveforms_b[cut_test,:],axis=0)
            alt_model_d_amp_b = (np.nanmax(alt_model_d_w_b)-np.nanmin(alt_model_d_w_b))/2.
            alt_model_d_ph_b = day[np.nanargmax(alt_model_d_w_b)]
            alt_model_s_w_b = np.nanmean(alt_model_seasonal_waveforms_b[cut_test,:],axis=0)
            alt_model_s_amp_b = (np.nanmax(alt_model_s_w_b)-np.nanmin(alt_model_s_w_b))/2.
            alt_model_s_ph_b = year[np.nanargmax(alt_model_s_w_b)]        
            alt_model_a_b = np.nanmean(alt_model_s_w_b)        

            alt_model_d_w_c = np.nanmean(alt_model_daily_waveforms_c[cut_test,:],axis=0)
            alt_model_d_amp_c = (np.nanmax(alt_model_d_w_c)-np.nanmin(alt_model_d_w_c))/2.
            alt_model_d_ph_c = day[np.nanargmax(alt_model_d_w_c)]
            alt_model_s_w_c = np.nanmean(alt_model_seasonal_waveforms_c[cut_test,:],axis=0)
            alt_model_s_amp_c = (np.nanmax(alt_model_s_w_c)-np.nanmin(alt_model_s_w_c))/2.
            alt_model_s_ph_c = year[np.nanargmax(alt_model_s_w_c)]        
            alt_model_a_c = np.nanmean(alt_model_s_w_c)

            ave_obs_param = obs_a
            ave_alt_model_param_a = alt_model_a_a
            ave_alt_model_param_b = alt_model_a_b
            ave_alt_model_param_c = alt_model_a_c
            diff_ave_a[a,m] = ave_alt_model_param_a - ave_obs_param
            diff_ave_b[a,m] = ave_alt_model_param_b - ave_obs_param
            diff_ave_c[a,m] = ave_alt_model_param_c - ave_obs_param

            ave_obs_param = obs_d_amp
            ave_alt_model_param_a = alt_model_d_amp_a
            ave_alt_model_param_b = alt_model_d_amp_b
            ave_alt_model_param_c = alt_model_d_amp_c
            diff_amp_d_a[a,m] = ave_alt_model_param_a - ave_obs_param
            diff_amp_d_b[a,m] = ave_alt_model_param_b - ave_obs_param
            diff_amp_d_c[a,m] = ave_alt_model_param_c - ave_obs_param
                
            ave_obs_param = obs_d_ph
            ave_alt_model_param_a = alt_model_d_ph_a
            ave_alt_model_param_b = alt_model_d_ph_b
            ave_alt_model_param_c = alt_model_d_ph_c
            diff_a = ave_alt_model_param_a - ave_obs_param
            diff_b = ave_alt_model_param_b - ave_obs_param
            diff_c = ave_alt_model_param_c - ave_obs_param
            if diff_a < -12.:
                diff_a = 12. - (np.abs(diff_a)-12.)
            if diff_b < -12.:
                diff_b = 12. - (np.abs(diff_b)-12.)
            if diff_c < -12.:
                diff_c = 12. - (np.abs(diff_c)-12.)   
            if diff_a > 12.:
                diff_a = -12. + (np.abs(diff_a)-12.)
            if diff_b > 12.:
                diff_b = -12. + (np.abs(diff_b)-12.)
            if diff_c > 12.:
                diff_c = -12. + (np.abs(diff_c)-12.)   
            diff_ph_d_a[a,m] = diff_a
            diff_ph_d_b[a,m] = diff_b
            diff_ph_d_c[a,m] = diff_c
                
            ave_obs_param = obs_d_w
            ave_alt_model_param_a = alt_model_d_w_a
            ave_alt_model_param_b = alt_model_d_w_b
            ave_alt_model_param_c = alt_model_d_w_c
            diff_wf_d_a[a,m] = np.sum(np.abs(ave_alt_model_param_a - ave_obs_param))
            diff_wf_d_b[a,m] = np.sum(np.abs(ave_alt_model_param_b - ave_obs_param))
            diff_wf_d_c[a,m] = np.sum(np.abs(ave_alt_model_param_c - ave_obs_param))    
                
            ave_obs_param = obs_s_amp
            ave_alt_model_param_a = alt_model_s_amp_a
            ave_alt_model_param_b = alt_model_s_amp_b
            ave_alt_model_param_c = alt_model_s_amp_c
            diff_amp_s_a[a,m] = ave_alt_model_param_a - ave_obs_param
            diff_amp_s_b[a,m] = ave_alt_model_param_b - ave_obs_param
            diff_amp_s_c[a,m] = ave_alt_model_param_c - ave_obs_param
                
            ave_obs_param = obs_s_ph
            ave_alt_model_param_a = alt_model_s_ph_a
            ave_alt_model_param_b = alt_model_s_ph_b
            ave_alt_model_param_c = alt_model_s_ph_c
            diff_a = ave_alt_model_param_a - ave_obs_param
            diff_b = ave_alt_model_param_b - ave_obs_param
            diff_c = ave_alt_model_param_c - ave_obs_param
            if diff_a < -6.:
                diff_a = 6. - (np.abs(diff_a)-6.)
            if diff_b < -6.:
                diff_b = 6. - (np.abs(diff_b)-6.)
            if diff_c < -6.:
                diff_c = 6. - (np.abs(diff_c)-6.)   
            if diff_a > 6.:
                diff_a = -6. + (np.abs(diff_a)-6.)
            if diff_b > 6.:
                diff_b = -6. + (np.abs(diff_b)-6.)
            if diff_c > 6.:
                diff_c = -6. + (np.abs(diff_c)-6.)    
            diff_ph_s_a[a,m] = diff_a
            diff_ph_s_b[a,m] = diff_b  
            diff_ph_s_c[a,m] = diff_c 
                
            ave_obs_param = obs_s_w
            ave_alt_model_param_a = alt_model_s_w_a
            ave_alt_model_param_b = alt_model_s_w_b
            ave_alt_model_param_c = alt_model_s_w_c
            diff_wf_s_a[a,m] = np.sum(np.abs(ave_alt_model_param_a - ave_obs_param))
            diff_wf_s_b[a,m] = np.sum(np.abs(ave_alt_model_param_b - ave_obs_param))
            diff_wf_s_c[a,m] = np.sum(np.abs(ave_alt_model_param_c - ave_obs_param)) 
        
        count+=1
            
    #remove unneeded variables
    try:
        del ave_obs_param
        del ave_alt_model_param_a
        del ave_alt_model_param_b
        del ave_alt_model_param_c
        del site_group_alt_model_a
        del site_group_alt_model_b
        del site_group_alt_model_c
        del alt_model_grp_a
        del alt_model_grp_b
        del alt_model_grp_c
        del alt_model_daily_amp_a 
        del alt_model_daily_ph_a 
        del alt_model_seasonal_amp_a 
        del alt_model_seasonal_ph_a 
        del alt_model_ave_a 
        del alt_model_daily_waveforms_a 
        del alt_model_seasonal_waveforms_a 
        del alt_model_daily_amp_b 
        del alt_model_daily_ph_b 
        del alt_model_seasonal_amp_b 
        del alt_model_seasonal_ph_b 
        del alt_model_ave_b 
        del alt_model_daily_waveforms_b 
        del alt_model_seasonal_waveforms_b 
        del alt_model_daily_amp_c 
        del alt_model_daily_ph_c 
        del alt_model_seasonal_amp_c 
        del alt_model_seasonal_ph_c 
        del alt_model_ave_c 
        del alt_model_daily_waveforms_c 
        del alt_model_seasonal_waveforms_c 
        del obs_d_amp 
        del obs_d_ph 
        del obs_s_amp 
        del obs_s_ph
        del obs_a 
        del obs_d_w 
        del obs_s_w
        del alt_model_d_amp_a 
        del alt_model_d_ph_a
        del alt_model_s_amp_a 
        del alt_model_s_ph_a
        del alt_model_a_a 
        del alt_model_d_w_a 
        del alt_model_s_w_a 
        del alt_model_d_amp_b 
        del alt_model_d_ph_b
        del alt_model_s_amp_b 
        del alt_model_s_ph_b 
        del alt_model_a_b 
        del alt_model_d_w_b 
        del alt_model_s_w_b 
        del alt_model_d_amp_c 
        del alt_model_d_ph_c 
        del alt_model_s_amp_c 
        del alt_model_s_ph_c
        del alt_model_a_c 
        del alt_model_d_w_c 
        del alt_model_s_w_c 
        del area_grid
        del area_tag
        del area_label
        del cut_test
    except:
        pass
    gc.collect()   
 
    print '\n'

plotter = 'T'
while plotter == 'T':
    fig, axes = plt.subplots(nrows=4, ncols=5,figsize=(17,13))
    fig.patch.set_facecolor('white')

    plot_type = raw_input('\nd, s, all_rank or ave?\n')

    if (plot_type != 'ave') & (plot_type != 'all_rank'):                                                                                                                                                            
        plot_type_2 = raw_input('\namp, ph, wf or rank?\n')
    else:
        plot_type_2 = ''

    set_type = raw_input('\nANMVOC ,BNMVOC or DRYDEPO3?\n')

    plot_grid = raw_input('\nbox or contour?\n')

    count = 0
    for ax in axes.flat:

        try:
            area = areas[count]
        except:
            ax.axis('off')
            continue
        
        area_label = area_labels[area]
    
        #ANMVOC cut
    
        ave_rank = np.argsort(diff_ave_a[count,:])
        amp_d_rank = np.argsort(diff_amp_d_a[count,:])
        ph_d_rank = np.argsort(diff_ph_d_a[count,:])            
        amp_s_rank = np.argsort(diff_amp_s_a[count,:])
        ph_s_rank = np.argsort(diff_ph_s_a[count,:])            
        d_rank = np.sum((ave_rank,amp_d_rank,ph_d_rank),axis=0)
        s_rank = np.sum((ave_rank,amp_s_rank,ph_s_rank),axis=0)

        lowest_s = np.min(s_rank)
        lowest_s_i = np.argmin(s_rank)
        lowest_s_file = alt_model_dirs_a[lowest_s_i]
        
        if plot_type == 'ave':
            area_grid_a = diff_ave_a[count,:]
        elif plot_type == 'd':
            if plot_type_2 == 'amp':
                area_grid_a = diff_amp_d_a[count,:]  
            elif plot_type_2 == 'ph':
                area_grid_a = diff_ph_d_a[count,:]
            elif plot_type_2 == 'wf':
                area_grid_a = diff_wf_d_a[count,:]
            elif plot_type_2 == 'rank':
                area_grid_a = d_rank
        elif plot_type == 's':
            if plot_type_2 == 'amp':
                area_grid_a = diff_amp_s_a[count,:]  
            elif plot_type_2 == 'ph':
                area_grid_a = diff_ph_s_a[count,:]
            elif plot_type_2 == 'wf':
                area_grid_a = diff_wf_s_a[count,:]
            elif plot_type_2 == 'rank':
                area_grid_a = s_rank
        #elif plot_type == 'all_rank':
        #    area_grid_a = all_rank

        #BNMVOC cut

        ave_rank = np.argsort(diff_ave_b[count,:])
        amp_d_rank = np.argsort(diff_amp_d_b[count,:])
        ph_d_rank = np.argsort(diff_ph_d_b[count,:])            
        amp_s_rank = np.argsort(diff_amp_s_b[count,:])
        ph_s_rank = np.argsort(diff_ph_s_b[count,:])            
        d_rank = np.sum((ave_rank,amp_d_rank,ph_d_rank),axis=0)
        s_rank = np.sum((ave_rank,amp_s_rank,ph_s_rank),axis=0)
        #all_rank = np.sum((ave_rank,amp_d_rank,ph_d_rank,ff_d_rank,amp_s_rank,ph_s_rank,ff_s_rank),axis=0)
        
        if np.min(s_rank) < lowest_s:
            lowest_s = np.min(s_rank)
            lowest_s_i = np.argmin(s_rank)
            lowest_s_file = alt_model_dirs_b[lowest_s_i]
            
        
        if plot_type == 'ave':
            area_grid_b = diff_ave_b[count,:]
        elif plot_type == 'd':
            if plot_type_2 == 'amp':
                area_grid_b = diff_amp_d_b[count,:]
            elif plot_type_2 == 'ph':
                area_grid_b = diff_ph_d_b[count,:]
            elif plot_type_2 == 'wf':
                area_grid_b = diff_wf_d_b[count,:]
            elif plot_type_2 == 'rank':
                area_grid_b = d_rank
        elif plot_type == 's':
            if plot_type_2 == 'amp':
                area_grid_b = diff_amp_s_b[count,:]
            elif plot_type_2 == 'ph':
                area_grid_b = diff_ph_s_b[count,:]
            elif plot_type_2 == 'wf':
                area_grid_b = diff_wf_s_b[count,:]
            elif plot_type_2 == 'rank':
                area_grid_b = s_rank
       # elif plot_type == 'all_rank':
       #     area_grid_b = all_rank
            
                
        #drydepo3 cut
        
        ave_rank = np.argsort(diff_ave_c[count,:])
        amp_d_rank = np.argsort(diff_amp_d_c[count,:])
        ph_d_rank = np.argsort(diff_ph_d_c[count,:])            
        amp_s_rank = np.argsort(diff_amp_s_c[count,:])
        ph_s_rank = np.argsort(diff_ph_s_c[count,:])            
        d_rank = np.sum((ave_rank,amp_d_rank,ph_d_rank),axis=0)
        s_rank = np.sum((ave_rank,amp_s_rank,ph_s_rank),axis=0)
        #all_rank = np.sum((ave_rank,amp_d_rank,ph_d_rank,ff_d_rank,amp_s_rank,ph_s_rank,ff_s_rank),axis=0)
        
        if np.min(s_rank) < lowest_s:
            lowest_s = np.min(s_rank)
            lowest_s_i = np.argmin(s_rank)
            lowest_s_file = alt_model_dirs_c[lowest_s_i]
        
        if plot_type == 'ave':
            area_grid_c = diff_ave_c[count,:]
        elif plot_type == 'd':
            if plot_type_2 == 'amp':
                area_grid_c = diff_amp_d_c[count,:]
            elif plot_type_2 == 'ph':
                area_grid_c = diff_ph_d_c[count,:]
            elif plot_type_2 == 'wf':
                area_grid_c = diff_wf_d_c[count,:]
            elif plot_type_2 == 'rank':
                area_grid_c = d_rank
        elif plot_type == 's':
            if plot_type_2 == 'amp':
                area_grid_c = diff_amp_s_c[count,:]
            elif plot_type_2 == 'ph':
                area_grid_c = diff_ph_s_c[count,:]
            elif plot_type_2 == 'wf':
                area_grid_c = diff_wf_s_c[count,:]
            elif plot_type_2 == 'rank':
                area_grid_c = s_rank
        #elif plot_type == 'all_rank':
        #        area_grid_c = all_rank

                
        if set_type == 'ANMVOC':
            area_grid = area_grid_a
        elif set_type == 'BNMVOC':
            area_grid = area_grid_b
        elif set_type == 'DRYDEPO3':
            area_grid = area_grid_c
            
        #set min and max
        
        if (plot_type_2 == 'amp'):
            if plot_type == 's':
                if species == 'O3':
                    minval = -20
                    maxval = 20
                    t = [-20,-10,0,10,20]
                    t_str = ['-20','-10','0','+10','+20']
                if species == 'NO':
                    minval = -10
                    maxval = 10
                    t = [-10,-5,0,5,10]
                    t_str = ['-10','-5','0','+5','+10']               
 
                cb_label = 'Amplitude Bias (ppb)'
                cts = np.linspace(minval,maxval,12)
                
                print np.min(area_grid),np.max(area_grid)
        
        if (plot_type == 'ave'):
            if species == 'O3':
                minval = -20
                maxval = 20
                t = [-20,-10,0,10,20]
                t_str = ['-20','-10','0','+10','+20']
                cts = np.linspace(-20,20,12)
            if species == 'NO':
                minval = -10
                maxval = 10
                t = [-10,-5,0,5,10]
                t_str = ['-10','-5','0','+5','+10']
            
            cb_label = 'Average Bias (ppb)'
        
        elif plot_type == 'all_rank':
            minval = 0
            maxval = 24*7
        elif plot_type_2 == 'rank':
            minval = 0
            maxval = 24*3
            t = [0,24,48,72]
            t_str = ['0','24','48','72']
            cb_label = 'Simulation Rank\n(0 = Best Possible, 72 = Worst Possible)'
            cts = np.linspace(0,72,12)
        elif (plot_type_2 == 'ph'): 
            if plot_type == 'd':
                minval = -12
                maxval = 12
            elif plot_type == 's':
                minval = -6
                maxval = 6   
                t = [-6,-3,0,3,6]
                t_str = ['-6','-3','0','+3','+6']
                cb_label = 'Phase Bias (months)' 
                cts = np.linspace(-6,6,24)

        elif (plot_type_2 == 'wf'):
            if species == 'O3':
                minval = 0
                maxval = 250000
                t = [0,250000]
                t_str = ['0','+250000']
                cb_label = 'Integrated Seasonal Bias'
                cts = np.linspace(0,250000,40)
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

        if (plot_type != 'all_rank') & (plot_type_2 != 'wf') & (plot_type_2 != 'rank'):
            cmap = matplotlib.cm.RdBu_r
        else:
            cmap = matplotlib.cm.cubehelix

        area_grid = np.reshape(area_grid,(5,5))    
        masked_array = area_grid
        masked_array = np.ma.array(area_grid, mask=np.isnan(area_grid))
        cmap.set_bad('w',1.)
        if plot_grid == 'box':
            pl = ax.pcolor(masked_array,vmin = minval,vmax=maxval,cmap =cmap)
        elif plot_grid == 'contour':
            pl = ax.contourf(masked_array,cts,vmin = minval,vmax=maxval,cmap =cmap)
        ax.set_xticks([0,1,2,3,4])
        ax.set_xticklabels(['0.25','0.5','1.0','2.0','4.0'])
        ax.set_yticks([0,1,2,3,4])
        ax.set_yticklabels(['0.25','0.5','1.0','2.0','4.0'])
        
        ax.tick_params(axis='both', which='major', labelsize=18,pad = 7)
            
        ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)
    
        print area, lowest_s_file
    
        count+=1
        
    plt.tight_layout(pad = 1.5)
    fig.subplots_adjust(bottom=0.08)
    fig.subplots_adjust(left=0.10)

    fig.text(0.5, 0.01, set_type, ha='center',fontsize=30)
    fig.text(0.03, 0.5, 'NOx', va='center', rotation='vertical',fontsize=30)

    cbar_ax = fig.add_axes([0.58, 0.15, 0.35, 0.06])
    cb = fig.colorbar(pl,orientation='horizontal',cax=cbar_ax,ticks=t)
    cb.set_ticklabels(t_str)
    cb.ax.tick_params(labelsize=24)
    cb.set_label(cb_label,fontsize=24)
    cb.set_clim(minval, maxval)

    plt.show()

    plotter = raw_input('\nAnother plot? T or F?\n')
