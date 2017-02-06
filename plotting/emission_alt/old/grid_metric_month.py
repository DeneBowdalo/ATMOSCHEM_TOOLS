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

print species

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

alt_model_dirs_d = ['GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ACO0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ACO0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ACO2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ACO4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ACO0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ACO0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ACO2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ACO4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ACO0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ACO0.5',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_STD','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ACO2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ACO4.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ACO0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ACO0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ACO2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ACO4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ACO0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ACO0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ACO2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ACO4.0']


obs_seasonal_waveforms = obs_grp.variables['seasonal_waveform'][:]

#-----------------------------------
#get area
areas = ['SW_NA','NW_NA','NE_NA','CE_NA','SE_NA','C_NA','S_EU','SW_EU','CW_EU','CS_EU','C_EU','E_EU','NW_EU','N_EU','NE_AS','SE_AS']

area_boundaries,area_tags,area_labels = modules.area_dicts()

diff_wf_s_a_jan = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_feb = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_mar = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_apr = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_may = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_jun = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_jul = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_aug = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_sep = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_oct = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_nov = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_a_dec = np.empty((len(areas),len(alt_model_dirs_a)))
diff_wf_s_b_jan = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_feb = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_mar = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_apr = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_may = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_jun = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_jul = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_aug = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_sep = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_oct = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_nov = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_b_dec = np.empty((len(areas),len(alt_model_dirs_b)))
diff_wf_s_c_jan = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_feb = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_mar = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_apr = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_may = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_jun = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_jul = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_aug = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_sep = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_oct = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_nov = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_c_dec = np.empty((len(areas),len(alt_model_dirs_c)))
diff_wf_s_d_jan = np.empty((len(areas),len(alt_model_dirs_d)))
diff_wf_s_d_feb = np.empty((len(areas),len(alt_model_dirs_d)))
diff_wf_s_d_mar = np.empty((len(areas),len(alt_model_dirs_d)))
diff_wf_s_d_apr = np.empty((len(areas),len(alt_model_dirs_d)))
diff_wf_s_d_may = np.empty((len(areas),len(alt_model_dirs_d)))
diff_wf_s_d_jun = np.empty((len(areas),len(alt_model_dirs_d)))
diff_wf_s_d_jul = np.empty((len(areas),len(alt_model_dirs_d)))
diff_wf_s_d_aug = np.empty((len(areas),len(alt_model_dirs_d)))
diff_wf_s_d_sep = np.empty((len(areas),len(alt_model_dirs_d)))
diff_wf_s_d_oct = np.empty((len(areas),len(alt_model_dirs_d)))
diff_wf_s_d_nov = np.empty((len(areas),len(alt_model_dirs_d)))
diff_wf_s_d_dec = np.empty((len(areas),len(alt_model_dirs_d)))

#cut vals into seasons
start = datetime.datetime(year = 2008, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = 2009, month = 1, day = 1, hour = 0, minute = 0)
ref_date_dt = pd.date_range(start,end,freq='H')[:-19]
months = np.array([d.strftime('%m') for d in ref_date_dt]).astype('int')

valid_inds_jan = (months == 1)
valid_inds_feb = (months == 2)
valid_inds_mar = (months == 3)
valid_inds_apr = (months == 4)
valid_inds_may = (months == 5)
valid_inds_jun = (months == 6)
valid_inds_jul = (months == 7)
valid_inds_aug = (months == 8)
valid_inds_sep = (months == 9)
valid_inds_oct = (months == 10)
valid_inds_nov = (months == 11)
valid_inds_dec = (months == 12)

for m in range(len(alt_model_dirs_a)):

    print 'point 1'

    print m
    print '../%s/LSP_stats.nc'%(alt_model_dirs_a[m])
    print '../%s/LSP_stats.nc'%(alt_model_dirs_b[m])
    print '../%s/LSP_stats.nc'%(alt_model_dirs_c[m])
    print '../%s/LSP_stats.nc'%(alt_model_dirs_d[m])
    
    alt_model_grp_a = Dataset('../%s/LSP_stats.nc'%(alt_model_dirs_a[m]))
    alt_model_grp_b = Dataset('../%s/LSP_stats.nc'%(alt_model_dirs_b[m]))
    alt_model_grp_c = Dataset('../%s/LSP_stats.nc'%(alt_model_dirs_c[m]))
    alt_model_grp_d = Dataset('../%s/LSP_stats.nc'%(alt_model_dirs_d[m]))        

    alt_model_seasonal_waveforms_a = alt_model_grp_a.variables['seasonal_waveform'][:]
    alt_model_seasonal_waveforms_b = alt_model_grp_b.variables['seasonal_waveform'][:]
    alt_model_seasonal_waveforms_c = alt_model_grp_c.variables['seasonal_waveform'][:]
    alt_model_seasonal_waveforms_d = alt_model_grp_d.variables['seasonal_waveform'][:]

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
            diff_wf_s_a_jan[a,m] = np.NaN
            diff_wf_s_a_feb[a,m] = np.NaN
            diff_wf_s_a_mar[a,m] = np.NaN
            diff_wf_s_a_apr[a,m] = np.NaN
            diff_wf_s_a_may[a,m] = np.NaN
            diff_wf_s_a_jun[a,m] = np.NaN
            diff_wf_s_a_jul[a,m] = np.NaN
            diff_wf_s_a_aug[a,m] = np.NaN
            diff_wf_s_a_sep[a,m] = np.NaN
            diff_wf_s_a_oct[a,m] = np.NaN
            diff_wf_s_a_nov[a,m] = np.NaN 
            diff_wf_s_a_dec[a,m] = np.NaN
            diff_wf_s_b_jan[a,m] = np.NaN
            diff_wf_s_b_feb[a,m] = np.NaN
            diff_wf_s_b_mar[a,m] = np.NaN
            diff_wf_s_b_apr[a,m] = np.NaN
            diff_wf_s_b_may[a,m] = np.NaN
            diff_wf_s_b_jun[a,m] = np.NaN
            diff_wf_s_b_jul[a,m] = np.NaN
            diff_wf_s_b_aug[a,m] = np.NaN
            diff_wf_s_b_sep[a,m] = np.NaN
            diff_wf_s_b_oct[a,m] = np.NaN
            diff_wf_s_b_nov[a,m] = np.NaN
            diff_wf_s_b_dec[a,m] = np.NaN
            diff_wf_s_c_jan[a,m] = np.NaN
            diff_wf_s_c_feb[a,m] = np.NaN
            diff_wf_s_c_mar[a,m] = np.NaN
            diff_wf_s_c_apr[a,m] = np.NaN
            diff_wf_s_c_may[a,m] = np.NaN
            diff_wf_s_c_jun[a,m] = np.NaN
            diff_wf_s_c_jul[a,m] = np.NaN
            diff_wf_s_c_aug[a,m] = np.NaN
            diff_wf_s_c_sep[a,m] = np.NaN
            diff_wf_s_c_oct[a,m] = np.NaN
            diff_wf_s_c_nov[a,m] = np.NaN
            diff_wf_s_c_dec[a,m] = np.NaN
            diff_wf_s_d_jan[a,m] = np.NaN
            diff_wf_s_d_feb[a,m] = np.NaN
            diff_wf_s_d_mar[a,m] = np.NaN
            diff_wf_s_d_apr[a,m] = np.NaN
            diff_wf_s_d_may[a,m] = np.NaN
            diff_wf_s_d_jun[a,m] = np.NaN
            diff_wf_s_d_jul[a,m] = np.NaN
            diff_wf_s_d_aug[a,m] = np.NaN
            diff_wf_s_d_sep[a,m] = np.NaN
            diff_wf_s_d_oct[a,m] = np.NaN
            diff_wf_s_d_nov[a,m] = np.NaN
            diff_wf_s_d_dec[a,m] = np.NaN
        else:
            obs_sites = obs_seasonal_waveforms[cut_test,:]
            obs_s_w_jan = np.nanmean(obs_sites[:,valid_inds_jan],axis=0)
            obs_s_w_feb = np.nanmean(obs_sites[:,valid_inds_feb],axis=0)
            obs_s_w_mar = np.nanmean(obs_sites[:,valid_inds_mar],axis=0)
            obs_s_w_apr = np.nanmean(obs_sites[:,valid_inds_apr],axis=0)
            obs_s_w_may = np.nanmean(obs_sites[:,valid_inds_may],axis=0)
            obs_s_w_jun = np.nanmean(obs_sites[:,valid_inds_jun],axis=0)
            obs_s_w_jul = np.nanmean(obs_sites[:,valid_inds_jul],axis=0)
            obs_s_w_aug = np.nanmean(obs_sites[:,valid_inds_aug],axis=0)
            obs_s_w_sep = np.nanmean(obs_sites[:,valid_inds_sep],axis=0)
            obs_s_w_oct = np.nanmean(obs_sites[:,valid_inds_oct],axis=0)
            obs_s_w_nov = np.nanmean(obs_sites[:,valid_inds_nov],axis=0)
            obs_s_w_dec = np.nanmean(obs_sites[:,valid_inds_dec],axis=0)    
 
            model_sites_a = alt_model_seasonal_waveforms_a[cut_test,:]
            model_sites_b = alt_model_seasonal_waveforms_b[cut_test,:]
            model_sites_c = alt_model_seasonal_waveforms_c[cut_test,:]
            model_sites_d = alt_model_seasonal_waveforms_d[cut_test,:]
            alt_model_s_w_a_jan = np.nanmean(model_sites_a[:,valid_inds_jan],axis=0)
            alt_model_s_w_a_feb = np.nanmean(model_sites_a[:,valid_inds_feb],axis=0)
            alt_model_s_w_a_mar = np.nanmean(model_sites_a[:,valid_inds_mar],axis=0)
            alt_model_s_w_a_apr = np.nanmean(model_sites_a[:,valid_inds_apr],axis=0)
            alt_model_s_w_a_may = np.nanmean(model_sites_a[:,valid_inds_may],axis=0)
            alt_model_s_w_a_jun = np.nanmean(model_sites_a[:,valid_inds_jun],axis=0)
            alt_model_s_w_a_jul = np.nanmean(model_sites_a[:,valid_inds_jul],axis=0)
            alt_model_s_w_a_aug = np.nanmean(model_sites_a[:,valid_inds_aug],axis=0)
            alt_model_s_w_a_sep = np.nanmean(model_sites_a[:,valid_inds_sep],axis=0)
            alt_model_s_w_a_oct = np.nanmean(model_sites_a[:,valid_inds_oct],axis=0)
            alt_model_s_w_a_nov = np.nanmean(model_sites_a[:,valid_inds_nov],axis=0)
            alt_model_s_w_a_dec = np.nanmean(model_sites_a[:,valid_inds_dec],axis=0)
            alt_model_s_w_b_jan = np.nanmean(model_sites_b[:,valid_inds_jan],axis=0)
            alt_model_s_w_b_feb = np.nanmean(model_sites_b[:,valid_inds_feb],axis=0)                                                                     
            alt_model_s_w_b_mar = np.nanmean(model_sites_b[:,valid_inds_mar],axis=0)
            alt_model_s_w_b_apr = np.nanmean(model_sites_b[:,valid_inds_apr],axis=0)
            alt_model_s_w_b_may = np.nanmean(model_sites_b[:,valid_inds_may],axis=0)
            alt_model_s_w_b_jun = np.nanmean(model_sites_b[:,valid_inds_jun],axis=0)
            alt_model_s_w_b_jul = np.nanmean(model_sites_b[:,valid_inds_jul],axis=0)
            alt_model_s_w_b_aug = np.nanmean(model_sites_b[:,valid_inds_aug],axis=0)
            alt_model_s_w_b_sep = np.nanmean(model_sites_b[:,valid_inds_sep],axis=0)
            alt_model_s_w_b_oct = np.nanmean(model_sites_b[:,valid_inds_oct],axis=0)
            alt_model_s_w_b_nov = np.nanmean(model_sites_b[:,valid_inds_nov],axis=0)
            alt_model_s_w_b_dec = np.nanmean(model_sites_b[:,valid_inds_dec],axis=0)
            alt_model_s_w_c_jan = np.nanmean(model_sites_c[:,valid_inds_jan],axis=0)
            alt_model_s_w_c_feb = np.nanmean(model_sites_c[:,valid_inds_feb],axis=0)                                                                     
            alt_model_s_w_c_mar = np.nanmean(model_sites_c[:,valid_inds_mar],axis=0)
            alt_model_s_w_c_apr = np.nanmean(model_sites_c[:,valid_inds_apr],axis=0)
            alt_model_s_w_c_may = np.nanmean(model_sites_c[:,valid_inds_may],axis=0)
            alt_model_s_w_c_jun = np.nanmean(model_sites_c[:,valid_inds_jun],axis=0)
            alt_model_s_w_c_jul = np.nanmean(model_sites_c[:,valid_inds_jul],axis=0)
            alt_model_s_w_c_aug = np.nanmean(model_sites_c[:,valid_inds_aug],axis=0)
            alt_model_s_w_c_sep = np.nanmean(model_sites_c[:,valid_inds_sep],axis=0)
            alt_model_s_w_c_oct = np.nanmean(model_sites_c[:,valid_inds_oct],axis=0)
            alt_model_s_w_c_nov = np.nanmean(model_sites_c[:,valid_inds_nov],axis=0)
            alt_model_s_w_c_dec = np.nanmean(model_sites_c[:,valid_inds_dec],axis=0)
            alt_model_s_w_d_jan = np.nanmean(model_sites_d[:,valid_inds_jan],axis=0)
            alt_model_s_w_d_feb = np.nanmean(model_sites_d[:,valid_inds_feb],axis=0)
            alt_model_s_w_d_mar = np.nanmean(model_sites_d[:,valid_inds_mar],axis=0)
            alt_model_s_w_d_apr = np.nanmean(model_sites_d[:,valid_inds_apr],axis=0)
            alt_model_s_w_d_may = np.nanmean(model_sites_d[:,valid_inds_may],axis=0)
            alt_model_s_w_d_jun = np.nanmean(model_sites_d[:,valid_inds_jun],axis=0)
            alt_model_s_w_d_jul = np.nanmean(model_sites_d[:,valid_inds_jul],axis=0)
            alt_model_s_w_d_aug = np.nanmean(model_sites_d[:,valid_inds_aug],axis=0)
            alt_model_s_w_d_sep = np.nanmean(model_sites_d[:,valid_inds_sep],axis=0)
            alt_model_s_w_d_oct = np.nanmean(model_sites_d[:,valid_inds_oct],axis=0)
            alt_model_s_w_d_nov = np.nanmean(model_sites_d[:,valid_inds_nov],axis=0)
            alt_model_s_w_d_dec = np.nanmean(model_sites_d[:,valid_inds_dec],axis=0)

            diff_wf_s_a_jan[a,m] = np.sum(np.abs(alt_model_s_w_a_jan - obs_s_w_jan))            
            diff_wf_s_a_feb[a,m] = np.sum(np.abs(alt_model_s_w_a_feb - obs_s_w_feb))
            diff_wf_s_a_mar[a,m] = np.sum(np.abs(alt_model_s_w_a_mar - obs_s_w_mar))
            diff_wf_s_a_apr[a,m] = np.sum(np.abs(alt_model_s_w_a_apr - obs_s_w_apr))
            diff_wf_s_a_may[a,m] = np.sum(np.abs(alt_model_s_w_a_may - obs_s_w_may))
            diff_wf_s_a_jun[a,m] = np.sum(np.abs(alt_model_s_w_a_jun - obs_s_w_jun))
            diff_wf_s_a_jul[a,m] = np.sum(np.abs(alt_model_s_w_a_jul - obs_s_w_jul))
            diff_wf_s_a_aug[a,m] = np.sum(np.abs(alt_model_s_w_a_aug - obs_s_w_aug))
            diff_wf_s_a_sep[a,m] = np.sum(np.abs(alt_model_s_w_a_sep - obs_s_w_sep))       
            diff_wf_s_a_oct[a,m] = np.sum(np.abs(alt_model_s_w_a_oct - obs_s_w_oct))
            diff_wf_s_a_nov[a,m] = np.sum(np.abs(alt_model_s_w_a_nov - obs_s_w_nov))
            diff_wf_s_a_dec[a,m] = np.sum(np.abs(alt_model_s_w_a_dec - obs_s_w_dec))
            diff_wf_s_b_jan[a,m] = np.sum(np.abs(alt_model_s_w_b_jan - obs_s_w_jan))
            diff_wf_s_b_feb[a,m] = np.sum(np.abs(alt_model_s_w_b_feb - obs_s_w_feb))
            diff_wf_s_b_mar[a,m] = np.sum(np.abs(alt_model_s_w_b_mar - obs_s_w_mar))
            diff_wf_s_b_apr[a,m] = np.sum(np.abs(alt_model_s_w_b_apr - obs_s_w_apr))
            diff_wf_s_b_may[a,m] = np.sum(np.abs(alt_model_s_w_b_may - obs_s_w_may))
            diff_wf_s_b_jun[a,m] = np.sum(np.abs(alt_model_s_w_b_jun - obs_s_w_jun))
            diff_wf_s_b_jul[a,m] = np.sum(np.abs(alt_model_s_w_b_jul - obs_s_w_jul))
            diff_wf_s_b_aug[a,m] = np.sum(np.abs(alt_model_s_w_b_aug - obs_s_w_aug))
            diff_wf_s_b_sep[a,m] = np.sum(np.abs(alt_model_s_w_b_sep - obs_s_w_sep))       
            diff_wf_s_b_oct[a,m] = np.sum(np.abs(alt_model_s_w_b_oct - obs_s_w_oct))
            diff_wf_s_b_nov[a,m] = np.sum(np.abs(alt_model_s_w_b_nov - obs_s_w_nov))
            diff_wf_s_b_dec[a,m] = np.sum(np.abs(alt_model_s_w_b_dec - obs_s_w_dec))
            diff_wf_s_c_jan[a,m] = np.sum(np.abs(alt_model_s_w_c_jan - obs_s_w_jan))
            diff_wf_s_c_feb[a,m] = np.sum(np.abs(alt_model_s_w_c_feb - obs_s_w_feb))
            diff_wf_s_c_mar[a,m] = np.sum(np.abs(alt_model_s_w_c_mar - obs_s_w_mar))
            diff_wf_s_c_apr[a,m] = np.sum(np.abs(alt_model_s_w_c_apr - obs_s_w_apr))
            diff_wf_s_c_may[a,m] = np.sum(np.abs(alt_model_s_w_c_may - obs_s_w_may))
            diff_wf_s_c_jun[a,m] = np.sum(np.abs(alt_model_s_w_c_jun - obs_s_w_jun))
            diff_wf_s_c_jul[a,m] = np.sum(np.abs(alt_model_s_w_c_jul - obs_s_w_jul))
            diff_wf_s_c_aug[a,m] = np.sum(np.abs(alt_model_s_w_c_aug - obs_s_w_aug))
            diff_wf_s_c_sep[a,m] = np.sum(np.abs(alt_model_s_w_c_sep - obs_s_w_sep))       
            diff_wf_s_c_oct[a,m] = np.sum(np.abs(alt_model_s_w_c_oct - obs_s_w_oct))
            diff_wf_s_c_nov[a,m] = np.sum(np.abs(alt_model_s_w_c_nov - obs_s_w_nov))
            diff_wf_s_c_dec[a,m] = np.sum(np.abs(alt_model_s_w_c_dec - obs_s_w_dec))
            diff_wf_s_d_jan[a,m] = np.sum(np.abs(alt_model_s_w_d_jan - obs_s_w_jan))
            diff_wf_s_d_feb[a,m] = np.sum(np.abs(alt_model_s_w_d_feb - obs_s_w_feb))
            diff_wf_s_d_mar[a,m] = np.sum(np.abs(alt_model_s_w_d_mar - obs_s_w_mar))
            diff_wf_s_d_apr[a,m] = np.sum(np.abs(alt_model_s_w_d_apr - obs_s_w_apr))
            diff_wf_s_d_may[a,m] = np.sum(np.abs(alt_model_s_w_d_may - obs_s_w_may))
            diff_wf_s_d_jun[a,m] = np.sum(np.abs(alt_model_s_w_d_jun - obs_s_w_jun))
            diff_wf_s_d_jul[a,m] = np.sum(np.abs(alt_model_s_w_d_jul - obs_s_w_jul))
            diff_wf_s_d_aug[a,m] = np.sum(np.abs(alt_model_s_w_d_aug - obs_s_w_aug))
            diff_wf_s_d_sep[a,m] = np.sum(np.abs(alt_model_s_w_d_sep - obs_s_w_sep))
            diff_wf_s_d_oct[a,m] = np.sum(np.abs(alt_model_s_w_d_oct - obs_s_w_oct))
            diff_wf_s_d_nov[a,m] = np.sum(np.abs(alt_model_s_w_d_nov - obs_s_w_nov))
            diff_wf_s_d_dec[a,m] = np.sum(np.abs(alt_model_s_w_d_dec - obs_s_w_dec))


        count+=1
            
    #remove unneeded variables
    try:
        del ave_obs_param
        del ave_alt_model_param_a
        del ave_alt_model_param_b
        del ave_alt_model_param_c
        del ave_alt_model_param_d
        del alt_model_seasonal_waveforms_a
        del alt_model_seasonal_waveforms_b 
        del alt_model_seasonal_waveforms_c
        del alt_model_seasonal_waveforms_d 
        del obs_s_w
        del alt_model_s_w_a 
        del alt_model_s_w_b 
        del alt_model_s_w_c
        del alt_model_s_w_d 
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
    fig = plt.figure(figsize = (14,14))
    fig.patch.set_facecolor('white')

    gs1 = gridspec.GridSpec(4, 3)
    gs1.update(top=0.99,bottom=0.79,left=0.04,right=0.27,wspace=0,hspace=0)
    ax1 = plt.subplot(gs1[0, 0])
    ax2 = plt.subplot(gs1[0, 1])
    ax3 = plt.subplot(gs1[0, 2])
    ax4 = plt.subplot(gs1[1, 0])
    ax5 = plt.subplot(gs1[1, 1])
    ax6 = plt.subplot(gs1[1, 2])
    ax7 = plt.subplot(gs1[2, 0])
    ax8 = plt.subplot(gs1[2, 1])
    ax9 = plt.subplot(gs1[2, 2])
    ax10 = plt.subplot(gs1[3, 0])
    ax11 = plt.subplot(gs1[3, 1])
    ax12 = plt.subplot(gs1[3, 2])
    gs2 = gridspec.GridSpec(4, 3)
    gs2.update(top=0.99,bottom=0.79,left=0.28, right=0.51,wspace=0,hspace=0)
    ax13 = plt.subplot(gs2[0, 0])
    ax14 = plt.subplot(gs2[0, 1])
    ax15 = plt.subplot(gs2[0, 2])
    ax16 = plt.subplot(gs2[1, 0])
    ax17 = plt.subplot(gs2[1, 1])
    ax18 = plt.subplot(gs2[1, 2])
    ax19 = plt.subplot(gs2[2, 0])
    ax20 = plt.subplot(gs2[2, 1])
    ax21 = plt.subplot(gs2[2, 2])
    ax22 = plt.subplot(gs2[3, 0])
    ax23 = plt.subplot(gs2[3, 1])
    ax24 = plt.subplot(gs2[3, 2])    
    gs3 = gridspec.GridSpec(4, 3)
    gs3.update(top=0.99,bottom=0.79,left=0.52, right=0.75,wspace=0,hspace=0)
    ax25 = plt.subplot(gs3[0, 0])
    ax26 = plt.subplot(gs3[0, 1])
    ax27 = plt.subplot(gs3[0, 2])
    ax28 = plt.subplot(gs3[1, 0])
    ax29 = plt.subplot(gs3[1, 1])
    ax30 = plt.subplot(gs3[1, 2])
    ax31 = plt.subplot(gs3[2, 0])
    ax32 = plt.subplot(gs3[2, 1])
    ax33 = plt.subplot(gs3[2, 2])
    ax34 = plt.subplot(gs3[3, 0])
    ax35 = plt.subplot(gs3[3, 1])
    ax36 = plt.subplot(gs3[3, 2])
    gs4 = gridspec.GridSpec(4, 3)
    gs4.update(top=0.99,bottom=0.79,left=0.76, right=0.99,wspace=0,hspace=0)
    ax37 = plt.subplot(gs4[0, 0])
    ax38 = plt.subplot(gs4[0, 1])
    ax39 = plt.subplot(gs4[0, 2])
    ax40 = plt.subplot(gs4[1, 0])
    ax41 = plt.subplot(gs4[1, 1])
    ax42 = plt.subplot(gs4[1, 2])
    ax43 = plt.subplot(gs4[2, 0])
    ax44 = plt.subplot(gs4[2, 1])
    ax45 = plt.subplot(gs4[2, 2])
    ax46 = plt.subplot(gs4[3, 0])
    ax47 = plt.subplot(gs4[3, 1])
    ax48 = plt.subplot(gs4[3, 2])
    gs5 = gridspec.GridSpec(4, 3)
    gs5.update(top=0.78,bottom=0.58,left=0.04, right=0.27,wspace=0,hspace=0)
    ax49 = plt.subplot(gs5[0, 0])
    ax50 = plt.subplot(gs5[0, 1])
    ax51 = plt.subplot(gs5[0, 2])
    ax52 = plt.subplot(gs5[1, 0])
    ax53 = plt.subplot(gs5[1, 1])
    ax54 = plt.subplot(gs5[1, 2])
    ax55 = plt.subplot(gs5[2, 0])
    ax56 = plt.subplot(gs5[2, 1])
    ax57 = plt.subplot(gs5[2, 2])
    ax58 = plt.subplot(gs5[3, 0])
    ax59 = plt.subplot(gs5[3, 1])
    ax60 = plt.subplot(gs5[3, 2])
    gs6 = gridspec.GridSpec(4, 3)
    gs6.update(top=0.78,bottom=0.58,left=0.28, right=0.51,wspace=0,hspace=0)
    ax61 = plt.subplot(gs6[0, 0])
    ax62 = plt.subplot(gs6[0, 1])
    ax63 = plt.subplot(gs6[0, 2])
    ax64 = plt.subplot(gs6[1, 0])
    ax65 = plt.subplot(gs6[1, 1])
    ax66 = plt.subplot(gs6[1, 2])
    ax67 = plt.subplot(gs6[2, 0])
    ax68 = plt.subplot(gs6[2, 1])
    ax69 = plt.subplot(gs6[2, 2])
    ax70 = plt.subplot(gs6[3, 0])
    ax71 = plt.subplot(gs6[3, 1])
    ax72 = plt.subplot(gs6[3, 2])
    gs7 = gridspec.GridSpec(4, 3)
    gs7.update(top=0.78,bottom=0.58,left=0.52, right=0.75,wspace=0,hspace=0)
    ax73 = plt.subplot(gs7[0, 0])
    ax74 = plt.subplot(gs7[0, 1])
    ax75 = plt.subplot(gs7[0, 2])
    ax76 = plt.subplot(gs7[1, 0])
    ax77 = plt.subplot(gs7[1, 1])
    ax78 = plt.subplot(gs7[1, 2])
    ax79 = plt.subplot(gs7[2, 0])
    ax80 = plt.subplot(gs7[2, 1])
    ax81 = plt.subplot(gs7[2, 2])
    ax82 = plt.subplot(gs7[3, 0])
    ax83 = plt.subplot(gs7[3, 1])
    ax84 = plt.subplot(gs7[3, 2])
    gs8 = gridspec.GridSpec(4, 3)
    gs8.update(top=0.78,bottom=0.58,left=0.76, right=0.99,wspace=0,hspace=0)
    ax85 = plt.subplot(gs8[0, 0])
    ax86 = plt.subplot(gs8[0, 1])
    ax87 = plt.subplot(gs8[0, 2])
    ax88 = plt.subplot(gs8[1, 0])
    ax89 = plt.subplot(gs8[1, 1])
    ax90 = plt.subplot(gs8[1, 2])
    ax91 = plt.subplot(gs8[2, 0])
    ax92 = plt.subplot(gs8[2, 1])
    ax93 = plt.subplot(gs8[2, 2])
    ax94 = plt.subplot(gs8[3, 0])
    ax95 = plt.subplot(gs8[3, 1])
    ax96 = plt.subplot(gs8[3, 2])
    gs9 = gridspec.GridSpec(4, 3)
    gs9.update(top=0.57,bottom=0.37,left=0.04, right=0.27,wspace=0,hspace=0)
    ax97 = plt.subplot(gs9[0, 0])
    ax98 = plt.subplot(gs9[0, 1])
    ax99 = plt.subplot(gs9[0, 2])
    ax100 = plt.subplot(gs9[1, 0])
    ax101 = plt.subplot(gs9[1, 1])
    ax102 = plt.subplot(gs9[1, 2])
    ax103 = plt.subplot(gs9[2, 0])
    ax104 = plt.subplot(gs9[2, 1])
    ax105 = plt.subplot(gs9[2, 2])
    ax106 = plt.subplot(gs9[3, 0])
    ax107 = plt.subplot(gs9[3, 1])
    ax108 = plt.subplot(gs9[3, 2])
    gs10 = gridspec.GridSpec(4, 3)    
    gs10.update(top=0.57,bottom=0.37,left=0.28, right=0.51,wspace=0,hspace=0)
    ax109 = plt.subplot(gs10[0, 0])
    ax110 = plt.subplot(gs10[0, 1])
    ax111 = plt.subplot(gs10[0, 2])
    ax112 = plt.subplot(gs10[1, 0])
    ax113 = plt.subplot(gs10[1, 1])
    ax114 = plt.subplot(gs10[1, 2])
    ax115 = plt.subplot(gs10[2, 0])
    ax116 = plt.subplot(gs10[2, 1])
    ax117 = plt.subplot(gs10[2, 2])
    ax118 = plt.subplot(gs10[3, 0])
    ax119 = plt.subplot(gs10[3, 1])
    ax120 = plt.subplot(gs10[3, 2])
    gs11 = gridspec.GridSpec(4, 3)
    gs11.update(top=0.57,bottom=0.37,left=0.52, right=0.75,wspace=0,hspace=0)
    ax121 = plt.subplot(gs11[0, 0])
    ax122 = plt.subplot(gs11[0, 1])
    ax123 = plt.subplot(gs11[0, 2])
    ax124 = plt.subplot(gs11[1, 0])
    ax125 = plt.subplot(gs11[1, 1])
    ax126 = plt.subplot(gs11[1, 2])
    ax127 = plt.subplot(gs11[2, 0])
    ax128 = plt.subplot(gs11[2, 1])
    ax129 = plt.subplot(gs11[2, 2])
    ax130 = plt.subplot(gs11[3, 0])
    ax131 = plt.subplot(gs11[3, 1])
    ax132 = plt.subplot(gs11[3, 2])
    gs12 = gridspec.GridSpec(4, 3)
    gs12.update(top=0.57,bottom=0.37,left=0.76, right=0.99,wspace=0,hspace=0)
    ax133 = plt.subplot(gs12[0, 0])
    ax134 = plt.subplot(gs12[0, 1])
    ax135 = plt.subplot(gs12[0, 2])
    ax136 = plt.subplot(gs12[1, 0])
    ax137 = plt.subplot(gs12[1, 1])
    ax138 = plt.subplot(gs12[1, 2])
    ax139 = plt.subplot(gs12[2, 0])
    ax140 = plt.subplot(gs12[2, 1])
    ax141 = plt.subplot(gs12[2, 2])
    ax142 = plt.subplot(gs12[3, 0])
    ax143 = plt.subplot(gs12[3, 1])
    ax144 = plt.subplot(gs12[3, 2])
    gs13 = gridspec.GridSpec(4, 3)
    gs13.update(top=0.36,bottom=0.16,left=0.04, right=0.27,wspace=0,hspace=0)
    ax145 = plt.subplot(gs13[0, 0])
    ax146 = plt.subplot(gs13[0, 1])
    ax147 = plt.subplot(gs13[0, 2])
    ax148 = plt.subplot(gs13[1, 0])
    ax149 = plt.subplot(gs13[1, 1])
    ax150 = plt.subplot(gs13[1, 2])
    ax151 = plt.subplot(gs13[2, 0])
    ax152 = plt.subplot(gs13[2, 1])
    ax153 = plt.subplot(gs13[2, 2])
    ax154 = plt.subplot(gs13[3, 0])
    ax155 = plt.subplot(gs13[3, 1])
    ax156 = plt.subplot(gs13[3, 2])
    gs14 = gridspec.GridSpec(4, 3)
    gs14.update(top=0.36,bottom=0.16,left=0.28, right=0.51,wspace=0,hspace=0)
    ax157 = plt.subplot(gs14[0, 0])
    ax158 = plt.subplot(gs14[0, 1])
    ax159 = plt.subplot(gs14[0, 2])
    ax160 = plt.subplot(gs14[1, 0])
    ax161 = plt.subplot(gs14[1, 1])
    ax162 = plt.subplot(gs14[1, 2])
    ax163 = plt.subplot(gs14[2, 0])
    ax164 = plt.subplot(gs14[2, 1])
    ax165 = plt.subplot(gs14[2, 2])
    ax166 = plt.subplot(gs14[3, 0])
    ax167 = plt.subplot(gs14[3, 1])
    ax168 = plt.subplot(gs14[3, 2])
    gs15 = gridspec.GridSpec(4, 3) 
    gs15.update(top=0.36,bottom=0.16,left=0.52, right=0.75,wspace=0,hspace=0)
    ax169 = plt.subplot(gs15[0, 0])
    ax170 = plt.subplot(gs15[0, 1])
    ax171 = plt.subplot(gs15[0, 2])
    ax172 = plt.subplot(gs15[1, 0])
    ax173 = plt.subplot(gs15[1, 1])
    ax174 = plt.subplot(gs15[1, 2])
    ax175 = plt.subplot(gs15[2, 0])
    ax176 = plt.subplot(gs15[2, 1])
    ax177 = plt.subplot(gs15[2, 2])
    ax178 = plt.subplot(gs15[3, 0])
    ax179 = plt.subplot(gs15[3, 1])
    ax180 = plt.subplot(gs15[3, 2])
    gs16 = gridspec.GridSpec(4, 3)
    gs16.update(top=0.36,bottom=0.16,left=0.76, right=0.99,wspace=0,hspace=0)                                                                                    
    ax181 = plt.subplot(gs16[0, 0])
    ax182 = plt.subplot(gs16[0, 1])
    ax183 = plt.subplot(gs16[0, 2])
    ax184 = plt.subplot(gs16[1, 0])
    ax185 = plt.subplot(gs16[1, 1])
    ax186 = plt.subplot(gs16[1, 2])
    ax187 = plt.subplot(gs16[2, 0])
    ax188 = plt.subplot(gs16[2, 1])
    ax189 = plt.subplot(gs16[2, 2])
    ax190 = plt.subplot(gs16[3, 0])
    ax191 = plt.subplot(gs16[3, 1])
    ax192 = plt.subplot(gs16[3, 2])

    axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14,ax15,ax16,ax17,ax18,ax19,ax20,ax21,ax22,ax23,ax24,ax25,ax26,ax27,ax28,ax29,ax30,ax31,ax32,ax33,ax34,ax35,ax36,ax37,ax38,ax39,ax40,ax41,ax42,ax43,ax44,ax45,ax46,ax47,ax48,ax49,ax50,ax51,ax52,ax53,ax54,ax55,ax56,ax57,ax58,ax59,ax60,ax61,ax62,ax63,ax64,ax65,ax66,ax67,ax68,ax69,ax70,ax71,ax72,ax73,ax74,ax75,ax76,ax77,ax78,ax79,ax80,ax81,ax82,ax83,ax84,ax85,ax86,ax87,ax88,ax89,ax90,ax91,ax92,ax93,ax94,ax95,ax96,ax97,ax98,ax99,ax100,ax101,ax102,ax103,ax104,ax105,ax106,ax107,ax108,ax109,ax110,ax111,ax112,ax113,ax114,ax115,ax116,ax117,ax118,ax119,ax120,ax121,ax122,ax123,ax124,ax125,ax126,ax127,ax128,ax129,ax130,ax131,ax132,ax133,ax134,ax135,ax136,ax137,ax138,ax139,ax140,ax141,ax142,ax143,ax144,ax145,ax146,ax147,ax148,ax149,ax150,ax151,ax152,ax153,ax154,ax155,ax156,ax157,ax158,ax159,ax160,ax161,ax162,ax163,ax164,ax165,ax166,ax167,ax168,ax169,ax170,ax171,ax172,ax173,ax174,ax175,ax176,ax177,ax178,ax179,ax180,ax181,ax182,ax183,ax184,ax185,ax186,ax187,ax188,ax189,ax190,ax191,ax192]

    set_type = raw_input('\nANMVOC ,BNMVOC, DRYDEPO3 or ACO?\n')

    plot_grid = raw_input('\nbox or contour?\n')

    area_count = 0
    ax_count = 0
    for area in areas:
        print area
        #ax.axis('off')
        
        area_label = area_labels[area]
    
        #ANMVOC cut
        area_grid_a_jan = diff_wf_s_a_jan[area_count,:]
        area_grid_a_feb = diff_wf_s_a_feb[area_count,:]
        area_grid_a_mar = diff_wf_s_a_mar[area_count,:]
        area_grid_a_apr = diff_wf_s_a_apr[area_count,:]
        area_grid_a_may = diff_wf_s_a_may[area_count,:]
        area_grid_a_jun = diff_wf_s_a_jun[area_count,:]
        area_grid_a_jul = diff_wf_s_a_jul[area_count,:]
        area_grid_a_aug = diff_wf_s_a_aug[area_count,:]
        area_grid_a_sep = diff_wf_s_a_sep[area_count,:]
        area_grid_a_oct = diff_wf_s_a_oct[area_count,:]
        area_grid_a_nov = diff_wf_s_a_nov[area_count,:]
        area_grid_a_dec = diff_wf_s_a_dec[area_count,:]
        #BNMVOC cut
        area_grid_b_jan = diff_wf_s_b_jan[area_count,:]
        area_grid_b_feb = diff_wf_s_b_feb[area_count,:]
        area_grid_b_mar = diff_wf_s_b_mar[area_count,:]
        area_grid_b_apr = diff_wf_s_b_apr[area_count,:]
        area_grid_b_may = diff_wf_s_b_may[area_count,:]
        area_grid_b_jun = diff_wf_s_b_jun[area_count,:]
        area_grid_b_jul = diff_wf_s_b_jul[area_count,:]
        area_grid_b_aug = diff_wf_s_b_aug[area_count,:]
        area_grid_b_sep = diff_wf_s_b_sep[area_count,:]
        area_grid_b_oct = diff_wf_s_b_oct[area_count,:]
        area_grid_b_nov = diff_wf_s_b_nov[area_count,:]
        area_grid_b_dec = diff_wf_s_b_dec[area_count,:]
        #drydepo3 cut
        area_grid_c_jan = diff_wf_s_c_jan[area_count,:]
        area_grid_c_feb = diff_wf_s_c_feb[area_count,:]
        area_grid_c_mar = diff_wf_s_c_mar[area_count,:]
        area_grid_c_apr = diff_wf_s_c_apr[area_count,:]
        area_grid_c_may = diff_wf_s_c_may[area_count,:]
        area_grid_c_jun = diff_wf_s_c_jun[area_count,:]
        area_grid_c_jul = diff_wf_s_c_jul[area_count,:]
        area_grid_c_aug = diff_wf_s_c_aug[area_count,:]
        area_grid_c_sep = diff_wf_s_c_sep[area_count,:]
        area_grid_c_oct = diff_wf_s_c_oct[area_count,:]
        area_grid_c_nov = diff_wf_s_c_nov[area_count,:]
        area_grid_c_dec = diff_wf_s_c_dec[area_count,:]
        #ACO cut
        #ACO cut
        area_grid_d_jan = diff_wf_s_d_jan[area_count,:]
        area_grid_d_feb = diff_wf_s_d_feb[area_count,:]
        area_grid_d_mar = diff_wf_s_d_mar[area_count,:]
        area_grid_d_apr = diff_wf_s_d_apr[area_count,:]
        area_grid_d_may = diff_wf_s_d_may[area_count,:]
        area_grid_d_jun = diff_wf_s_d_jun[area_count,:]
        area_grid_d_jul = diff_wf_s_d_jul[area_count,:]
        area_grid_d_aug = diff_wf_s_d_aug[area_count,:]
        area_grid_d_sep = diff_wf_s_d_sep[area_count,:]
        area_grid_d_oct = diff_wf_s_d_oct[area_count,:]
        area_grid_d_nov = diff_wf_s_d_nov[area_count,:]
        area_grid_d_dec = diff_wf_s_d_dec[area_count,:]
        

        if set_type == 'ANMVOC':
            area_grid_jan = area_grid_a_jan
            area_grid_feb = area_grid_a_feb
            area_grid_mar = area_grid_a_mar
            area_grid_apr = area_grid_a_apr
            area_grid_may = area_grid_a_may
            area_grid_jun = area_grid_a_jun
            area_grid_jul = area_grid_a_jul
            area_grid_aug = area_grid_a_aug
            area_grid_sep = area_grid_a_sep
            area_grid_oct = area_grid_a_oct
            area_grid_nov = area_grid_a_nov
            area_grid_dec = area_grid_a_dec
        elif set_type == 'BNMVOC':
            area_grid_jan = area_grid_b_jan
            area_grid_feb = area_grid_b_feb                                                                                                              
            area_grid_mar = area_grid_b_mar
            area_grid_apr = area_grid_b_apr
            area_grid_may = area_grid_b_may
            area_grid_jun = area_grid_b_jun
            area_grid_jul = area_grid_b_jul
            area_grid_aug = area_grid_b_aug
            area_grid_sep = area_grid_b_sep
            area_grid_oct = area_grid_b_oct
            area_grid_nov = area_grid_b_nov
            area_grid_dec = area_grid_b_dec
        elif set_type == 'DRYDEPO3':
            area_grid_jan = area_grid_c_jan
            area_grid_feb = area_grid_c_feb                                                                                                              
            area_grid_mar = area_grid_c_mar
            area_grid_apr = area_grid_c_apr
            area_grid_may = area_grid_c_may
            area_grid_jun = area_grid_c_jun
            area_grid_jul = area_grid_c_jul
            area_grid_aug = area_grid_c_aug
            area_grid_sep = area_grid_c_sep
            area_grid_oct = area_grid_c_oct
            area_grid_nov = area_grid_c_nov
            area_grid_dec = area_grid_c_dec
        elif set_type == 'ACO':            
            area_grid_jan = area_grid_d_jan
            area_grid_feb = area_grid_d_feb
            area_grid_mar = area_grid_d_mar
            area_grid_apr = area_grid_d_apr
            area_grid_may = area_grid_d_may
            area_grid_jun = area_grid_d_jun
            area_grid_jul = area_grid_d_jul
            area_grid_aug = area_grid_d_aug
            area_grid_sep = area_grid_d_sep
            area_grid_oct = area_grid_d_oct
            area_grid_nov = area_grid_d_nov
            area_grid_dec = area_grid_d_dec

        #set min and max        
        if species == 'O3':
            minval = 100
            maxval = 35000
            t = [100,500,1000,2500,5000,10000,15000,20000,35000]
            t_str = ['0','500','1000','2500','5000','10000','15000','20000','35000']
            cb_label = 'Integrated Seasonal Bias'
            #cts = np.linspace(minval,maxval,40)
            cts = [100,500,1000,2500,5000,10000,15000,20000,35000]
        elif species == 'NO':                                                                                                                        
            minval = 1
            maxval = 60000
            t = [1,10,100,500,1000,2500,5000,10000,15000,20000,30000,60000]
            t_str = ['0','10','100','500','1000','2500','5000','10000','15000','20000','30000','60000']
            cb_label = 'Integrated Seasonal Bias'
            #cts = np.linspace(0,20000,40)
            cts = [1,10,100,500,1000,2500,5000,10000,15000,20000,30000,60000]
        elif species == 'NO2-MOLYBDENUM':                                                                                                                        
            minval = 1
            maxval = 60000
            t = [1,10,100,500,1000,2500,5000,10000,15000,20000,30000,60000]
            t_str = ['0','10','100','500','1000','2500','5000','10000','15000','20000','30000','60000']
            cb_label = 'Integrated Seasonal Bias'
            #cts = np.linspace(0,20000,40)
            cts = [1,10,100,500,1000,2500,5000,10000,15000,20000,30000,60000]
        elif species == 'CO':
            minval = 100
            maxval = 250000
            t = [100,500,1000,2500,5000,10000,15000,20000,30000,50000,100000,250000]
            t_str = ['0','500','1000','2500','5000','10000','15000','20000','30000','50000','100000','250000']
            cb_label = 'Integrated Seasonal Bias'
            cts = [100,500,1000,2500,5000,10000,15000,20000,30000,50000,100000,250000]

        cmap = matplotlib.cm.rainbow

        area_grid_jan = np.reshape(area_grid_jan,(5,5))    
        area_grid_feb = np.reshape(area_grid_feb,(5,5))
        area_grid_mar = np.reshape(area_grid_mar,(5,5))
        area_grid_apr = np.reshape(area_grid_apr,(5,5))
        area_grid_may = np.reshape(area_grid_may,(5,5))
        area_grid_jun = np.reshape(area_grid_jun,(5,5))
        area_grid_jul = np.reshape(area_grid_jul,(5,5))
        area_grid_aug = np.reshape(area_grid_aug,(5,5))
        area_grid_sep = np.reshape(area_grid_sep,(5,5))
        area_grid_oct = np.reshape(area_grid_oct,(5,5))
        area_grid_nov = np.reshape(area_grid_nov,(5,5))
        area_grid_dec = np.reshape(area_grid_dec,(5,5))
 
        masked_array_jan = np.ma.array(area_grid_jan, mask=np.isnan(area_grid_jan))
        masked_array_feb = np.ma.array(area_grid_feb, mask=np.isnan(area_grid_feb))
        masked_array_mar = np.ma.array(area_grid_may, mask=np.isnan(area_grid_mar))
        masked_array_apr = np.ma.array(area_grid_apr, mask=np.isnan(area_grid_apr))
        masked_array_may = np.ma.array(area_grid_may, mask=np.isnan(area_grid_may))
        masked_array_jun = np.ma.array(area_grid_jun, mask=np.isnan(area_grid_jun))
        masked_array_jul = np.ma.array(area_grid_jul, mask=np.isnan(area_grid_jul))
        masked_array_aug = np.ma.array(area_grid_aug, mask=np.isnan(area_grid_aug))
        masked_array_sep = np.ma.array(area_grid_sep, mask=np.isnan(area_grid_sep))
        masked_array_oct = np.ma.array(area_grid_oct, mask=np.isnan(area_grid_oct))
        masked_array_nov = np.ma.array(area_grid_nov, mask=np.isnan(area_grid_nov))
        masked_array_dec = np.ma.array(area_grid_dec, mask=np.isnan(area_grid_dec))

        #cmap.set_bad('w',1.)
        if plot_grid == 'box':
            pl = axes[ax_count].pcolor(masked_array_jan,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+1].pcolor(masked_array_feb,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+2].pcolor(masked_array_mar,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+3].pcolor(masked_array_apr,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+4].pcolor(masked_array_may,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+5].pcolor(masked_array_jun,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+6].pcolor(masked_array_jul,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+7].pcolor(masked_array_aug,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+8].pcolor(masked_array_sep,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+9].pcolor(masked_array_oct,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+10].pcolor(masked_array_nov,vmin = minval,vmax=maxval,cmap =cmap)
            pl = axes[ax_count+11].pcolor(masked_array_dec,vmin = minval,vmax=maxval,cmap =cmap)
        elif plot_grid == 'contour':
            pl = axes[ax_count].contourf(masked_array_jan,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)
            pl = axes[ax_count+1].contourf(masked_array_feb,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)
            pl = axes[ax_count+2].contourf(masked_array_mar,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)
            pl = axes[ax_count+3].contourf(masked_array_apr,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)
            pl = axes[ax_count+4].contourf(masked_array_may,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)
            pl = axes[ax_count+5].contourf(masked_array_jun,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)
            pl = axes[ax_count+6].contourf(masked_array_jul,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)
            pl = axes[ax_count+7].contourf(masked_array_aug,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)
            pl = axes[ax_count+8].contourf(masked_array_sep,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)
            pl = axes[ax_count+9].contourf(masked_array_oct,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)
            pl = axes[ax_count+10].contourf(masked_array_nov,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)
            pl = axes[ax_count+11].contourf(masked_array_dec,cts,norm=LogNorm(vmin = minval,vmax=maxval),cmap =cmap)        

        axes[ax_count].set_xticks([1,2,3,4])
        axes[ax_count+1].set_xticks([1,2,3,4])
        axes[ax_count+2].set_xticks([1,2,3,4])
        axes[ax_count+3].set_xticks([1,2,3,4])
        axes[ax_count+4].set_xticks([1,2,3,4])
        axes[ax_count+5].set_xticks([1,2,3,4])
        axes[ax_count+6].set_xticks([1,2,3,4])
        axes[ax_count+7].set_xticks([1,2,3,4])
        axes[ax_count+8].set_xticks([1,2,3,4])
        axes[ax_count+9].set_xticks([1,2,3,4])
        axes[ax_count+10].set_xticks([1,2,3,4])
        axes[ax_count+11].set_xticks([1,2,3,4])

        axes[ax_count].set_xticklabels(['','','',''])
        axes[ax_count+1].set_xticklabels(['','','',''])
        axes[ax_count+2].set_xticklabels(['','','',''])
        axes[ax_count+3].set_xticklabels(['','','',''])
        axes[ax_count+4].set_xticklabels(['','','',''])
        axes[ax_count+5].set_xticklabels(['','','',''])
        axes[ax_count+6].set_xticklabels(['','','',''])
        axes[ax_count+7].set_xticklabels(['','','',''])
        axes[ax_count+8].set_xticklabels(['','','',''])
        axes[ax_count+9].set_xticklabels(['','','',''])
        axes[ax_count+10].set_xticklabels(['','','',''])
        axes[ax_count+11].set_xticklabels(['','','',''])

        axes[ax_count].set_yticks([1,2,3,4])
        axes[ax_count+1].set_yticks([1,2,3,4])
        axes[ax_count+2].set_yticks([1,2,3,4])
        axes[ax_count+3].set_yticks([1,2,3,4])
        axes[ax_count+4].set_yticks([1,2,3,4])
        axes[ax_count+5].set_yticks([1,2,3,4])
        axes[ax_count+6].set_yticks([1,2,3,4])
        axes[ax_count+7].set_yticks([1,2,3,4])
        axes[ax_count+8].set_yticks([1,2,3,4])
        axes[ax_count+9].set_yticks([1,2,3,4])
        axes[ax_count+10].set_yticks([1,2,3,4])
        axes[ax_count+11].set_yticks([1,2,3,4])
        
        axes[ax_count].set_yticklabels(['','','',''])
        axes[ax_count+1].set_yticklabels(['','','','']) 
        axes[ax_count+2].set_yticklabels(['','','','']) 
        axes[ax_count+3].set_yticklabels(['','','','']) 
        axes[ax_count+4].set_yticklabels(['','','',''])
        axes[ax_count+5].set_yticklabels(['','','','']) 
        axes[ax_count+6].set_yticklabels(['','','','']) 
        axes[ax_count+7].set_yticklabels(['','','',''])
        axes[ax_count+8].set_yticklabels(['','','',''])
        axes[ax_count+9].set_yticklabels(['','','','']) 
        axes[ax_count+10].set_yticklabels(['','','','']) 
        axes[ax_count+11].set_yticklabels(['','','',''])        

        #axes[ax_count].axes.get_xaxis().set_visible(False)
        #axes[ax_count].axes.get_yaxis().set_visible(False)
        #axes[ax_count+1].axes.get_xaxis().set_visible(False)
        #axes[ax_count+1].axes.get_yaxis().set_visible(False)
        #axes[ax_count+2].axes.get_xaxis().set_visible(False)
        #axes[ax_count+2].axes.get_yaxis().set_visible(False)
        #axes[ax_count+3].axes.get_xaxis().set_visible(False)
        #axes[ax_count+3].axes.get_yaxis().set_visible(False)

        #axes[ax_count].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        #axes[ax_count+1].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        #axes[ax_count+2].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        #axes[ax_count+3].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        #axes[ax_count+4].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        #axes[ax_count+5].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        #axes[ax_count+6].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        #axes[ax_count+7].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        #axes[ax_count+8].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        #axes[ax_count+9].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        #axes[ax_count+10].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)
        #axes[ax_count+11].grid(b=True,which='major',color='white',linestyle='--',linewidth=0.42)

        axes[ax_count+2].text(0.70, 0.85, area_label, ha='center', va='center',transform=axes[ax_count+2].transAxes,fontsize=15)
    
        area_count+=1
        ax_count+=12
        
    #plt.tight_layout(pad = 1.5)
    fig.subplots_adjust(bottom=0.08)
    fig.subplots_adjust(left=0.10)

    fig.text(0.5, 0.12, set_type, ha='center',fontsize=30)
    fig.text(0.01, 0.5, 'NOx', va='center', rotation='vertical',fontsize=30)

    cbar_ax = fig.add_axes([0.58, 0.07, 0.35, 0.06])
    cb = fig.colorbar(pl,orientation='horizontal',cax=cbar_ax,ticks=t)
    cb.set_ticklabels(t_str)
    cb.ax.tick_params(labelsize=10)
    cb.set_label(cb_label,fontsize=20)
    cb.set_clim(minval, maxval)

    plt.show()

    plotter = raw_input('\nAnother plot? T or F?\n')
