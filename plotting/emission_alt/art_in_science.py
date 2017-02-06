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
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic

#seasonal or diurnal
pt = 'seasonal'
area = 'CE_NA'

species = ['O3']

area_boundaries,area_tags,area_labels = modules.area_dicts()
color_dict = {'SW_NA':0,'C_NA':2,'NW_NA':4,'NE_NA':6,'CE_NA':8,'SE_NA':10,'S_NA':3,'NW_EU':1,'C_EU':3,'N_EU':0,'E_EU':7,'S_EU':9,'SW_EU':11,'NE_AS':6}

#set colormap
cmap = cm.get_cmap('Set1',13)
ratio = 1./13
lw = 4

#cut vals into seasons and hours
start = datetime.datetime(year = 2008, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = 2009, month = 1, day = 1, hour = 0, minute = 0)
ref_date_dt = pd.date_range(start,end,freq='H')[:-19]
months = np.array([d.strftime('%m') for d in ref_date_dt]).astype('int')

if pt == 'seasonal':
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
    cut_inds = [valid_inds_jan,valid_inds_feb,valid_inds_mar,valid_inds_apr,valid_inds_may,valid_inds_jun,valid_inds_jul,valid_inds_aug,valid_inds_sep,valid_inds_oct,valid_inds_nov,valid_inds_dec]

elif pt == 'diurnal':
    cut_inds = [np.array([0,1]),np.array([2,3]),np.array([4,5]),np.array([6,7]),np.array([8,9]),np.array([10,11]),np.array([12,13]),np.array([14,15]),np.array([16,17]),np.array([18,19]),np.array([20,21]),np.array([22,23])]


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

all_models = ['GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ANMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ANMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ANMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ANMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ANMVOC0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ANMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ANMVOC2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ANMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ANMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ANMVOC0.5',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_STD','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ANMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ANMVOC4.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ANMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ANMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ANMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ANMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ANMVOC0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ANMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ANMVOC2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ANMVOC4.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25BNMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25BNMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25BNMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25BNMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5BNMVOC0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5BNMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5BNMVOC2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5BNMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0BNMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0BNMVOC0.5',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_STD','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0BNMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0BNMVOC4.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0BNMVOC0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0BNMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0BNMVOC2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0BNMVOC4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0BNMVOC0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0BNMVOC0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0BNMVOC2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0BNMVOC4.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25DRYDEPO30.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25DRYDEPO30.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25DRYDEPO32.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25DRYDEPO34.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5DRYDEPO30.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5DRYDEPO30.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5DRYDEPO32.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5DRYDEPO34.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0DRYDEPO30.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0DRYDEPO30.5',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_STD','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0DRYDEPO32.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0DRYDEPO34.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0DRYDEPO30.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0DRYDEPO30.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0DRYDEPO32.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0DRYDEPO34.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0DRYDEPO30.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0DRYDEPO30.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0DRYDEPO32.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0DRYDEPO34.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ACO0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ACO0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ACO2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.25ACO4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ACO0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ACO0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ACO2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX0.5ACO4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ACO0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ACO0.5',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_STD','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ACO2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX1.0ACO4.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ACO0.25','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ACO0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0NMVOC1.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ACO2.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX2.0ACO4.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ACO0.25',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ACO0.5','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0NMVOC1.0','GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ACO2.0',
                    'GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_NOX4.0ACO4.0']

#-----------------------------------------------
diff_wf_a = np.empty((len(species),len(alt_model_dirs_a),12))
diff_wf_b = np.empty((len(species),len(alt_model_dirs_b),12))
diff_wf_c = np.empty((len(species),len(alt_model_dirs_c),12))
diff_wf_d = np.empty((len(species),len(alt_model_dirs_d),12))
diff_wf = np.empty((len(species),len(all_models),12))
perdiff_wf = np.empty((len(species),len(all_models),12))
obs_wf = np.empty((len(species),8766))
model_wf = np.empty((len(species),len(all_models),8766))

for x in range(len(species)):

    spec = species[x]

    #read in obs ts data
    obs_fname = '/work/home/db876/observations/surface/%s/process/GLOBAL_SURFACE_%s_2009_2011_H_HP.nc'%(spec,spec)
    obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_var,obs_lats,obs_lons,obs_alt,obs_groups,obs_raw_class,obs_anthrome_class,gap_inds = modules.read_obs_all(obs_fname,spec,2009,2011)
    tags = np.array(modules.get_tags(obs_refs))

    #read in std model data
    model_fname = '/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_2009_2011_v1001_4x5_GEOS5_H_STD.nc'
    model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,spec,2009,2011)

    #load obs lsp data
    obs_grp = Dataset('/work/home/db876/xy/%s/2009_2011/2009_2011/obs_SURFACE_H/LSP_stats.nc'%(spec))
    if pt == 'seasonal':
        obs_waveforms = obs_grp.variables['seasonal_waveform_night'][:]
    elif pt == 'diurnal':
        obs_waveforms = obs_grp.variables['diurnal_average_waveform'][:]
        
    collected_models = []
    all_model_count = 0

    for m in range(len(alt_model_dirs_a)):

        print 'point 1'

        print m
        print '/work/home/db876/xy/%s/2009_2011/2009_2011/%s/LSP_stats.nc'%(spec,alt_model_dirs_a[m])
        print '/work/home/db876/xy/%s/2009_2011/2009_2011/%s/LSP_stats.nc'%(spec,alt_model_dirs_b[m])
        print '/work/home/db876/xy/%s/2009_2011/2009_2011/%s/LSP_stats.nc'%(spec,alt_model_dirs_c[m])
        print '/work/home/db876/xy/%s/2009_2011/2009_2011/%s/LSP_stats.nc'%(spec,alt_model_dirs_d[m])
    
        model_grp_a = Dataset('/work/home/db876/xy/%s/2009_2011/2009_2011/%s/LSP_stats.nc'%(spec,alt_model_dirs_a[m]))
        model_grp_b = Dataset('/work/home/db876/xy/%s/2009_2011/2009_2011/%s/LSP_stats.nc'%(spec,alt_model_dirs_b[m]))
        model_grp_c = Dataset('/work/home/db876/xy/%s/2009_2011/2009_2011/%s/LSP_stats.nc'%(spec,alt_model_dirs_c[m]))
        model_grp_d = Dataset('/work/home/db876/xy/%s/2009_2011/2009_2011/%s/LSP_stats.nc'%(spec,alt_model_dirs_d[m]))         
    
        if pt == 'seasonal':
            model_waveforms_a = model_grp_a.variables['seasonal_waveform_night'][:]
            model_waveforms_b = model_grp_b.variables['seasonal_waveform_night'][:]
            model_waveforms_c = model_grp_c.variables['seasonal_waveform_night'][:]
            model_waveforms_d = model_grp_d.variables['seasonal_waveform_night'][:]
        elif pt == 'diurnal':
            model_waveforms_a = model_grp_a.variables['diurnal_average_waveform'][:]
            model_waveforms_b = model_grp_b.variables['diurnal_average_waveform'][:]
            model_waveforms_c = model_grp_c.variables['diurnal_average_waveform'][:]
            model_waveforms_d = model_grp_d.variables['diurnal_average_waveform'][:]
        

        area_grid = area_boundaries[area]
        area_tag = area_tags[area]
        area_label = area_labels[area]

        cut_test = modules.area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag)
        if np.all(cut_test == False):
            for c in range(12):
                diff_wf_a[x,m,c] = np.NaN
                diff_wf_b[x,m,c] = np.NaN
                diff_wf_c[x,m,c] = np.NaN
                diff_wf_d[x,m,c] = np.NaN
        
        else:
            obs_sites = obs_waveforms[cut_test,:]
            model_sites_a = model_waveforms_a[cut_test,:]
            model_sites_b = model_waveforms_b[cut_test,:]
            model_sites_c = model_waveforms_c[cut_test,:]
            model_sites_d = model_waveforms_d[cut_test,:]
            
            obs_wf[x,:] = np.nanmean(obs_sites,axis=0)
            model_wf[x,all_model_count,:] = np.nanmean(model_sites_a,axis=0)
            model_wf[x,all_model_count+1,:] = np.nanmean(model_sites_b,axis=0)
            model_wf[x,all_model_count+2,:] = np.nanmean(model_sites_c,axis=0)
            model_wf[x,all_model_count+3,:] = np.nanmean(model_sites_d,axis=0)
        
            for c in range(12):
                obs_w_cut = np.nanmean(obs_sites[:,cut_inds[c]],axis=0)
                model_a_w_cut = np.nanmean(model_sites_a[:,cut_inds[c]],axis=0)
                model_b_w_cut = np.nanmean(model_sites_b[:,cut_inds[c]],axis=0)
                model_c_w_cut = np.nanmean(model_sites_c[:,cut_inds[c]],axis=0)
                model_d_w_cut = np.nanmean(model_sites_d[:,cut_inds[c]],axis=0)

                diff_a = np.nanmean(model_a_w_cut) - np.nanmean(obs_w_cut)
                diff_b = np.nanmean(model_b_w_cut) - np.nanmean(obs_w_cut)
                diff_c = np.nanmean(model_c_w_cut) - np.nanmean(obs_w_cut)
                diff_d = np.nanmean(model_d_w_cut) - np.nanmean(obs_w_cut)
            
                diff_wf_a[x,m,c] = diff_a  
                diff_wf_b[x,m,c] = diff_b    
                diff_wf_c[x,m,c] = diff_c  
                diff_wf_d[x,m,c] = diff_d     
                
                diff_wf[x,all_model_count,c] = np.abs(diff_a)
                diff_wf[x,all_model_count+1,c] = np.abs(diff_b)
                diff_wf[x,all_model_count+2,c] = np.abs(diff_c) 
                diff_wf[x,all_model_count+3,c] = np.abs(diff_d)

                perdiff_wf[x,all_model_count,c] = np.abs(np.diff([1.,np.nanmean(model_a_w_cut)/np.nanmean(obs_w_cut)])[0])*100 
                perdiff_wf[x,all_model_count+1,c] = np.abs(np.diff([1.,np.nanmean(model_b_w_cut)/np.nanmean(obs_w_cut)])[0])*100 
                perdiff_wf[x,all_model_count+2,c] = np.abs(np.diff([1.,np.nanmean(model_c_w_cut)/np.nanmean(obs_w_cut)])[0])*100 
                perdiff_wf[x,all_model_count+3,c] = np.abs(np.diff([1.,np.nanmean(model_d_w_cut)/np.nanmean(obs_w_cut)])[0])*100                     
              
        collected_models.append(alt_model_dirs_a[m])
        collected_models.append(alt_model_dirs_b[m])
        collected_models.append(alt_model_dirs_c[m])
        collected_models.append(alt_model_dirs_d[m])
        
        all_model_count+=4
            
        #remove unneeded variables
        try:
            del model_grp_a
            del model_grp_b
            del model_grp_c
            del model_grp_d
            del model_waveforms_a
            del model_waveforms_b
            del model_waveforms_c
            del model_waveforms_d
            del obs_w_cut
            del model_a_w_cut
            del model_b_w_cut
            del model_c_w_cut
            del model_d_w_cut
            del area_grid
            del area_tag
            del area_label
            del cut_test
        except:
            pass
        gc.collect()   
 
        print '\n'


#-----------------------------------------------
#SET MAP AREA AXES
fig =plt.figure(figsize=(13.5,13.5))
fig.patch.set_facecolor('white')
gs1 = gridspec.GridSpec(8, 6)
gs1.update(top=1.0,bottom=0.666,left=0,right=0.333,wspace=0,hspace=0)
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
anmvoc_1 = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12]

ax13 = plt.subplot(gs1[0, 3])
ax14 = plt.subplot(gs1[0, 4])
ax15 = plt.subplot(gs1[0, 5])
ax16 = plt.subplot(gs1[1, 3])
ax17 = plt.subplot(gs1[1, 4])
ax18 = plt.subplot(gs1[1, 5])
ax19 = plt.subplot(gs1[2, 3])
ax20 = plt.subplot(gs1[2, 4])
ax21 = plt.subplot(gs1[2, 5])
ax22 = plt.subplot(gs1[3, 3])
ax23 = plt.subplot(gs1[3, 4])
ax24 = plt.subplot(gs1[3, 5])
bnmvoc_1 = [ax13,ax14,ax15,ax16,ax17,ax18,ax19,ax20,ax21,ax22,ax23,ax24]

ax25 = plt.subplot(gs1[4, 0])
ax26 = plt.subplot(gs1[4, 1])
ax27 = plt.subplot(gs1[4, 2])
ax28 = plt.subplot(gs1[5, 0])
ax29 = plt.subplot(gs1[5, 1])
ax30 = plt.subplot(gs1[5, 2])
ax31 = plt.subplot(gs1[6, 0])
ax32 = plt.subplot(gs1[6, 1])
ax33 = plt.subplot(gs1[6, 2])
ax34 = plt.subplot(gs1[7, 0])
ax35 = plt.subplot(gs1[7, 1])
ax36 = plt.subplot(gs1[7, 2])
co_1 = [ax25,ax26,ax27,ax28,ax29,ax30,ax31,ax32,ax33,ax34,ax35,ax36]

ax37 = plt.subplot(gs1[4, 3])
ax38 = plt.subplot(gs1[4, 4])
ax39 = plt.subplot(gs1[4, 5])
ax40 = plt.subplot(gs1[5, 3])
ax41 = plt.subplot(gs1[5, 4])
ax42 = plt.subplot(gs1[5, 5])
ax43 = plt.subplot(gs1[6, 3])
ax44 = plt.subplot(gs1[6, 4])
ax45 = plt.subplot(gs1[6, 5])
ax46 = plt.subplot(gs1[7, 3])
ax47 = plt.subplot(gs1[7, 4])
ax48 = plt.subplot(gs1[7, 5])
drydep_1 = [ax37,ax38,ax39,ax40,ax41,ax42,ax43,ax44,ax45,ax46,ax47,ax48]


gs2 = gridspec.GridSpec(8, 6)
gs2.update(top=1.0,bottom=0.666,left=0.333,right=0.666,wspace=0,hspace=0)
ax49 = plt.subplot(gs2[0, 0])
ax50 = plt.subplot(gs2[0, 1])
ax51 = plt.subplot(gs2[0, 2])
ax52 = plt.subplot(gs2[1, 0])
ax53 = plt.subplot(gs2[1, 1])
ax54 = plt.subplot(gs2[1, 2])
ax55 = plt.subplot(gs2[2, 0])
ax56 = plt.subplot(gs2[2, 1])
ax57 = plt.subplot(gs2[2, 2])
ax58 = plt.subplot(gs2[3, 0])
ax59 = plt.subplot(gs2[3, 1])
ax60 = plt.subplot(gs2[3, 2])
anmvoc_2 = [ax49,ax50,ax51,ax52,ax53,ax54,ax55,ax56,ax57,ax58,ax59,ax60]

ax61 = plt.subplot(gs2[0, 3])
ax62 = plt.subplot(gs2[0, 4])
ax63 = plt.subplot(gs2[0, 5])
ax64 = plt.subplot(gs2[1, 3])
ax65 = plt.subplot(gs2[1, 4])
ax66 = plt.subplot(gs2[1, 5])
ax67 = plt.subplot(gs2[2, 3])
ax68 = plt.subplot(gs2[2, 4])
ax69 = plt.subplot(gs2[2, 5])
ax70 = plt.subplot(gs2[3, 3])
ax71 = plt.subplot(gs2[3, 4])
ax72 = plt.subplot(gs2[3, 5])
bnmvoc_2 = [ax61,ax62,ax63,ax64,ax65,ax66,ax67,ax68,ax69,ax70,ax71,ax72]

ax73 = plt.subplot(gs2[4, 0])
ax74 = plt.subplot(gs2[4, 1])
ax75 = plt.subplot(gs2[4, 2])
ax76 = plt.subplot(gs2[5, 0])
ax77 = plt.subplot(gs2[5, 1])
ax78 = plt.subplot(gs2[5, 2])
ax79 = plt.subplot(gs2[6, 0])
ax80 = plt.subplot(gs2[6, 1])
ax81 = plt.subplot(gs2[6, 2])
ax82 = plt.subplot(gs2[7, 0])
ax83 = plt.subplot(gs2[7, 1])
ax84 = plt.subplot(gs2[7, 2])
co_2 = [ax73,ax74,ax75,ax76,ax77,ax78,ax79,ax80,ax81,ax82,ax83,ax84]

ax85 = plt.subplot(gs2[4, 3])
ax86 = plt.subplot(gs2[4, 4])
ax87 = plt.subplot(gs2[4, 5])
ax88 = plt.subplot(gs2[5, 3])
ax89 = plt.subplot(gs2[5, 4])
ax90 = plt.subplot(gs2[5, 5])
ax91 = plt.subplot(gs2[6, 3])
ax92 = plt.subplot(gs2[6, 4])
ax93 = plt.subplot(gs2[6, 5])
ax94 = plt.subplot(gs2[7, 3])
ax95 = plt.subplot(gs2[7, 4])
ax96 = plt.subplot(gs2[7, 5])
drydep_2 = [ax85,ax86,ax87,ax88,ax89,ax90,ax91,ax92,ax93,ax94,ax95,ax96]


gs3 = gridspec.GridSpec(8, 6)
gs3.update(top=1.0,bottom=0.666,left=0.666,right=1.0,wspace=0,hspace=0)
ax97 = plt.subplot(gs3[0, 0])
ax98 = plt.subplot(gs3[0, 1])
ax99 = plt.subplot(gs3[0, 2])
ax100 = plt.subplot(gs3[1, 0])
ax101 = plt.subplot(gs3[1, 1])
ax102 = plt.subplot(gs3[1, 2])
ax103 = plt.subplot(gs3[2, 0])
ax104 = plt.subplot(gs3[2, 1])
ax105 = plt.subplot(gs3[2, 2])
ax106 = plt.subplot(gs3[3, 0])
ax107 = plt.subplot(gs3[3, 1])
ax108 = plt.subplot(gs3[3, 2])
anmvoc_3 = [ax97,ax98,ax99,ax100,ax101,ax102,ax103,ax104,ax105,ax106,ax107,ax108]

ax109 = plt.subplot(gs3[0, 3])
ax110 = plt.subplot(gs3[0, 4])
ax111 = plt.subplot(gs3[0, 5])
ax112 = plt.subplot(gs3[1, 3])
ax113 = plt.subplot(gs3[1, 4])
ax114 = plt.subplot(gs3[1, 5])
ax115 = plt.subplot(gs3[2, 3])
ax116 = plt.subplot(gs3[2, 4])
ax117 = plt.subplot(gs3[2, 5])
ax118 = plt.subplot(gs3[3, 3])
ax119 = plt.subplot(gs3[3, 4])
ax120 = plt.subplot(gs3[3, 5])
bnmvoc_3 = [ax109,ax110,ax111,ax112,ax113,ax114,ax115,ax116,ax117,ax118,ax119,ax120]

ax121 = plt.subplot(gs3[4, 0])
ax122 = plt.subplot(gs3[4, 1])
ax123 = plt.subplot(gs3[4, 2])
ax124 = plt.subplot(gs3[5, 0])
ax125 = plt.subplot(gs3[5, 1])
ax126 = plt.subplot(gs3[5, 2])
ax127 = plt.subplot(gs3[6, 0])
ax128 = plt.subplot(gs3[6, 1])
ax129 = plt.subplot(gs3[6, 2])
ax130 = plt.subplot(gs3[7, 0])
ax131 = plt.subplot(gs3[7, 1])
ax132 = plt.subplot(gs3[7, 2])
co_3 = [ax121,ax122,ax123,ax124,ax125,ax126,ax127,ax128,ax129,ax130,ax131,ax132]

ax133 = plt.subplot(gs3[4, 3])
ax134 = plt.subplot(gs3[4, 4])
ax135 = plt.subplot(gs3[4, 5])
ax136 = plt.subplot(gs3[5, 3])
ax137 = plt.subplot(gs3[5, 4])
ax138 = plt.subplot(gs3[5, 5])
ax139 = plt.subplot(gs3[6, 3])
ax140 = plt.subplot(gs3[6, 4])
ax141 = plt.subplot(gs3[6, 5])
ax142 = plt.subplot(gs3[7, 3])
ax143 = plt.subplot(gs3[7, 4])
ax144 = plt.subplot(gs3[7, 5])
drydep_3 = [ax133,ax134,ax135,ax136,ax137,ax138,ax139,ax140,ax141,ax142,ax143,ax144]


gs4 = gridspec.GridSpec(8, 6)
gs4.update(top=0.666,bottom=0.333,left=0,right=0.333,wspace=0,hspace=0)
ax145 = plt.subplot(gs4[0, 0])
ax146 = plt.subplot(gs4[0, 1])
ax147 = plt.subplot(gs4[0, 2])
ax148 = plt.subplot(gs4[1, 0])
ax149 = plt.subplot(gs4[1, 1])
ax150 = plt.subplot(gs4[1, 2])
ax151 = plt.subplot(gs4[2, 0])
ax152 = plt.subplot(gs4[2, 1])
ax153 = plt.subplot(gs4[2, 2])
ax154 = plt.subplot(gs4[3, 0])
ax155 = plt.subplot(gs4[3, 1])
ax156 = plt.subplot(gs4[3, 2])
anmvoc_4 = [ax145,ax146,ax147,ax148,ax149,ax150,ax151,ax152,ax153,ax154,ax155,ax156]

ax157 = plt.subplot(gs4[0, 3])
ax158 = plt.subplot(gs4[0, 4])
ax159 = plt.subplot(gs4[0, 5])
ax160 = plt.subplot(gs4[1, 3])
ax161 = plt.subplot(gs4[1, 4])
ax162 = plt.subplot(gs4[1, 5])
ax163 = plt.subplot(gs4[2, 3])
ax164 = plt.subplot(gs4[2, 4])
ax165 = plt.subplot(gs4[2, 5])
ax166 = plt.subplot(gs4[3, 3])
ax167 = plt.subplot(gs4[3, 4])
ax168 = plt.subplot(gs4[3, 5])
bnmvoc_4 = [ax157,ax158,ax159,ax160,ax161,ax162,ax163,ax164,ax165,ax166,ax167,ax168]

ax169 = plt.subplot(gs4[4, 0])
ax170 = plt.subplot(gs4[4, 1])
ax171 = plt.subplot(gs4[4, 2])
ax172 = plt.subplot(gs4[5, 0])
ax173 = plt.subplot(gs4[5, 1])
ax174 = plt.subplot(gs4[5, 2])
ax175 = plt.subplot(gs4[6, 0])
ax176 = plt.subplot(gs4[6, 1])
ax177 = plt.subplot(gs4[6, 2])
ax178 = plt.subplot(gs4[7, 0])
ax179 = plt.subplot(gs4[7, 1])
ax180 = plt.subplot(gs4[7, 2])
co_4 = [ax169,ax170,ax171,ax172,ax173,ax174,ax175,ax176,ax177,ax178,ax179,ax180]

ax181 = plt.subplot(gs4[4, 3])
ax182 = plt.subplot(gs4[4, 4])
ax183 = plt.subplot(gs4[4, 5])
ax184 = plt.subplot(gs4[5, 3])
ax185 = plt.subplot(gs4[5, 4])
ax186 = plt.subplot(gs4[5, 5])
ax187 = plt.subplot(gs4[6, 3])
ax188 = plt.subplot(gs4[6, 4])
ax189 = plt.subplot(gs4[6, 5])
ax190 = plt.subplot(gs4[7, 3])
ax191 = plt.subplot(gs4[7, 4])
ax192 = plt.subplot(gs4[7, 5])
drydep_4 = [ax181,ax182,ax183,ax184,ax185,ax186,ax187,ax188,ax189,ax190,ax191,ax192]


gs5 = gridspec.GridSpec(8, 6)
gs5.update(top=0.666,bottom=0.333,left=0.333,right=0.666,wspace=0,hspace=0)
ax193 = plt.subplot(gs5[0, 0])
ax194 = plt.subplot(gs5[0, 1])
ax195 = plt.subplot(gs5[0, 2])
ax196 = plt.subplot(gs5[1, 0])
ax197 = plt.subplot(gs5[1, 1])
ax198 = plt.subplot(gs5[1, 2])
ax199 = plt.subplot(gs5[2, 0])
ax200 = plt.subplot(gs5[2, 1])
ax201 = plt.subplot(gs5[2, 2])
ax202 = plt.subplot(gs5[3, 0])
ax203 = plt.subplot(gs5[3, 1])
ax204 = plt.subplot(gs5[3, 2])
anmvoc_5 = [ax193,ax194,ax195,ax196,ax197,ax198,ax199,ax200,ax201,ax202,ax203,ax204]

ax205 = plt.subplot(gs5[0, 3])
ax206 = plt.subplot(gs5[0, 4])
ax207 = plt.subplot(gs5[0, 5])
ax208 = plt.subplot(gs5[1, 3])
ax209 = plt.subplot(gs5[1, 4])
ax210 = plt.subplot(gs5[1, 5])
ax211 = plt.subplot(gs5[2, 3])
ax212 = plt.subplot(gs5[2, 4])
ax213 = plt.subplot(gs5[2, 5])
ax214 = plt.subplot(gs5[3, 3])
ax215 = plt.subplot(gs5[3, 4])
ax216 = plt.subplot(gs5[3, 5])
bnmvoc_5 = [ax205,ax206,ax207,ax208,ax209,ax210,ax211,ax212,ax213,ax214,ax215,ax216]

ax217 = plt.subplot(gs5[4, 0])
ax218 = plt.subplot(gs5[4, 1])
ax219 = plt.subplot(gs5[4, 2])
ax220 = plt.subplot(gs5[5, 0])
ax221 = plt.subplot(gs5[5, 1])
ax222 = plt.subplot(gs5[5, 2])
ax223 = plt.subplot(gs5[6, 0])
ax224 = plt.subplot(gs5[6, 1])
ax225 = plt.subplot(gs5[6, 2])
ax226 = plt.subplot(gs5[7, 0])
ax227 = plt.subplot(gs5[7, 1])
ax228 = plt.subplot(gs5[7, 2])
co_5 = [ax217,ax218,ax219,ax220,ax221,ax222,ax223,ax224,ax225,ax226,ax227,ax228]

ax229 = plt.subplot(gs5[4, 3])
ax230 = plt.subplot(gs5[4, 4])
ax231 = plt.subplot(gs5[4, 5])
ax232 = plt.subplot(gs5[5, 3])
ax233 = plt.subplot(gs5[5, 4])
ax234 = plt.subplot(gs5[5, 5])
ax235 = plt.subplot(gs5[6, 3])
ax236 = plt.subplot(gs5[6, 4])
ax237 = plt.subplot(gs5[6, 5])
ax238 = plt.subplot(gs5[7, 3])
ax239 = plt.subplot(gs5[7, 4])
ax240 = plt.subplot(gs5[7, 5])
drydep_5 = [ax229,ax230,ax231,ax232,ax233,ax234,ax235,ax236,ax237,ax238,ax239,ax240]

gs6 = gridspec.GridSpec(8, 6)
gs6.update(top=0.666,bottom=0.333,left=0.666,right=1.0,wspace=0,hspace=0)
ax241 = plt.subplot(gs6[0, 0])
ax242 = plt.subplot(gs6[0, 1])
ax243 = plt.subplot(gs6[0, 2])
ax244 = plt.subplot(gs6[1, 0])
ax245 = plt.subplot(gs6[1, 1])
ax246 = plt.subplot(gs6[1, 2])
ax247 = plt.subplot(gs6[2, 0])
ax248 = plt.subplot(gs6[2, 1])
ax249 = plt.subplot(gs6[2, 2])
ax250 = plt.subplot(gs6[3, 0])
ax251 = plt.subplot(gs6[3, 1])
ax252 = plt.subplot(gs6[3, 2])
anmvoc_6 = [ax241,ax242,ax243,ax244,ax245,ax246,ax247,ax248,ax249,ax250,ax251,ax252]

ax253 = plt.subplot(gs6[0, 3])
ax254 = plt.subplot(gs6[0, 4])
ax255 = plt.subplot(gs6[0, 5])
ax256 = plt.subplot(gs6[1, 3])
ax257 = plt.subplot(gs6[1, 4])
ax258 = plt.subplot(gs6[1, 5])
ax259 = plt.subplot(gs6[2, 3])
ax260 = plt.subplot(gs6[2, 4])
ax261 = plt.subplot(gs6[2, 5])
ax262 = plt.subplot(gs6[3, 3])
ax263 = plt.subplot(gs6[3, 4])
ax264 = plt.subplot(gs6[3, 5])
bnmvoc_6 = [ax253,ax254,ax255,ax256,ax257,ax258,ax259,ax260,ax261,ax262,ax263,ax264]

ax265 = plt.subplot(gs6[4, 0])
ax266 = plt.subplot(gs6[4, 1])
ax267 = plt.subplot(gs6[4, 2])
ax268 = plt.subplot(gs6[5, 0])
ax269 = plt.subplot(gs6[5, 1])
ax270 = plt.subplot(gs6[5, 2])
ax271 = plt.subplot(gs6[6, 0])
ax272 = plt.subplot(gs6[6, 1])
ax273 = plt.subplot(gs6[6, 2])
ax274 = plt.subplot(gs6[7, 0])
ax275 = plt.subplot(gs6[7, 1])
ax276 = plt.subplot(gs6[7, 2])
co_6 = [ax265,ax266,ax267,ax268,ax269,ax270,ax271,ax272,ax273,ax274,ax275,ax276]

ax277 = plt.subplot(gs6[4, 3])
ax278 = plt.subplot(gs6[4, 4])
ax279 = plt.subplot(gs6[4, 5])
ax280 = plt.subplot(gs6[5, 3])
ax281 = plt.subplot(gs6[5, 4])
ax282 = plt.subplot(gs6[5, 5])
ax283 = plt.subplot(gs6[6, 3])
ax284 = plt.subplot(gs6[6, 4])
ax285 = plt.subplot(gs6[6, 5])
ax286 = plt.subplot(gs6[7, 3])
ax287 = plt.subplot(gs6[7, 4])
ax288 = plt.subplot(gs6[7, 5])
drydep_6 = [ax277,ax278,ax279,ax280,ax281,ax282,ax283,ax284,ax285,ax286,ax287,ax288]

gs7 = gridspec.GridSpec(8, 6)
gs7.update(top=0.333,bottom=0,left=0,right=0.333,wspace=0,hspace=0)
ax289 = plt.subplot(gs7[0, 0])
ax290 = plt.subplot(gs7[0, 1])
ax291 = plt.subplot(gs7[0, 2])
ax292 = plt.subplot(gs7[1, 0])
ax293 = plt.subplot(gs7[1, 1])
ax294 = plt.subplot(gs7[1, 2])
ax295 = plt.subplot(gs7[2, 0])
ax296 = plt.subplot(gs7[2, 1])
ax297 = plt.subplot(gs7[2, 2])
ax298 = plt.subplot(gs7[3, 0])
ax299 = plt.subplot(gs7[3, 1])
ax300 = plt.subplot(gs7[3, 2])
anmvoc_7 = [ax289,ax290,ax291,ax292,ax293,ax294,ax295,ax296,ax297,ax298,ax299,ax300]

ax301 = plt.subplot(gs7[0, 3])
ax302 = plt.subplot(gs7[0, 4])
ax303 = plt.subplot(gs7[0, 5])
ax304 = plt.subplot(gs7[1, 3])
ax305 = plt.subplot(gs7[1, 4])
ax306 = plt.subplot(gs7[1, 5])
ax307 = plt.subplot(gs7[2, 3])
ax308 = plt.subplot(gs7[2, 4])
ax309 = plt.subplot(gs7[2, 5])
ax310 = plt.subplot(gs7[3, 3])
ax311 = plt.subplot(gs7[3, 4])
ax312 = plt.subplot(gs7[3, 5])
bnmvoc_7 = [ax301,ax302,ax303,ax304,ax305,ax306,ax307,ax308,ax309,ax310,ax311,ax312]

ax313 = plt.subplot(gs7[4, 0])
ax314 = plt.subplot(gs7[4, 1])
ax315 = plt.subplot(gs7[4, 2])
ax316 = plt.subplot(gs7[5, 0])
ax317 = plt.subplot(gs7[5, 1])
ax318 = plt.subplot(gs7[5, 2])
ax319 = plt.subplot(gs7[6, 0])
ax320 = plt.subplot(gs7[6, 1])
ax321 = plt.subplot(gs7[6, 2])
ax322 = plt.subplot(gs7[7, 0])
ax323 = plt.subplot(gs7[7, 1])
ax324 = plt.subplot(gs7[7, 2])
co_7 = [ax313,ax314,ax315,ax316,ax317,ax318,ax319,ax320,ax321,ax322,ax323,ax324]

ax325 = plt.subplot(gs7[4, 3])
ax326 = plt.subplot(gs7[4, 4])
ax327 = plt.subplot(gs7[4, 5])
ax328 = plt.subplot(gs7[5, 3])
ax329 = plt.subplot(gs7[5, 4])
ax330 = plt.subplot(gs7[5, 5])
ax331 = plt.subplot(gs7[6, 3])
ax332 = plt.subplot(gs7[6, 4])
ax333 = plt.subplot(gs7[6, 5])
ax334 = plt.subplot(gs7[7, 3])
ax335 = plt.subplot(gs7[7, 4])
ax336 = plt.subplot(gs7[7, 5])
drydep_7 = [ax325,ax326,ax327,ax328,ax329,ax330,ax331,ax332,ax333,ax334,ax335,ax336]

gs8 = gridspec.GridSpec(8, 6)
gs8.update(top=0.333,bottom=0,left=0.333,right=0.666,wspace=0,hspace=0)
ax337 = plt.subplot(gs8[0, 0])
ax338 = plt.subplot(gs8[0, 1])
ax339 = plt.subplot(gs8[0, 2])
ax340 = plt.subplot(gs8[1, 0])
ax341 = plt.subplot(gs8[1, 1])
ax342 = plt.subplot(gs8[1, 2])
ax343 = plt.subplot(gs8[2, 0])
ax344 = plt.subplot(gs8[2, 1])
ax345 = plt.subplot(gs8[2, 2])
ax346 = plt.subplot(gs8[3, 0])
ax347 = plt.subplot(gs8[3, 1])
ax348 = plt.subplot(gs8[3, 2])
anmvoc_8 = [ax337,ax338,ax339,ax340,ax341,ax342,ax343,ax344,ax345,ax346,ax347,ax348]

ax349 = plt.subplot(gs8[0, 3])
ax350 = plt.subplot(gs8[0, 4])
ax351 = plt.subplot(gs8[0, 5])
ax352 = plt.subplot(gs8[1, 3])
ax353 = plt.subplot(gs8[1, 4])
ax354 = plt.subplot(gs8[1, 5])
ax355 = plt.subplot(gs8[2, 3])
ax356 = plt.subplot(gs8[2, 4])
ax357 = plt.subplot(gs8[2, 5])
ax358 = plt.subplot(gs8[3, 3])
ax359 = plt.subplot(gs8[3, 4])
ax360 = plt.subplot(gs8[3, 5])
bnmvoc_8 = [ax349,ax350,ax351,ax352,ax353,ax354,ax355,ax356,ax357,ax358,ax359,ax360]

ax361 = plt.subplot(gs8[4, 0])
ax362 = plt.subplot(gs8[4, 1])
ax363 = plt.subplot(gs8[4, 2])
ax364 = plt.subplot(gs8[5, 0])
ax365 = plt.subplot(gs8[5, 1])
ax366 = plt.subplot(gs8[5, 2])
ax367 = plt.subplot(gs8[6, 0])
ax368 = plt.subplot(gs8[6, 1])
ax369 = plt.subplot(gs8[6, 2])
ax370 = plt.subplot(gs8[7, 0])
ax371 = plt.subplot(gs8[7, 1])
ax372 = plt.subplot(gs8[7, 2])
co_8 = [ax361,ax362,ax363,ax364,ax365,ax366,ax367,ax368,ax369,ax370,ax371,ax372]

ax373 = plt.subplot(gs8[4, 3])
ax374 = plt.subplot(gs8[4, 4])
ax375 = plt.subplot(gs8[4, 5])
ax376 = plt.subplot(gs8[5, 3])
ax377 = plt.subplot(gs8[5, 4])
ax378 = plt.subplot(gs8[5, 5])
ax379 = plt.subplot(gs8[6, 3])
ax380 = plt.subplot(gs8[6, 4])
ax381 = plt.subplot(gs8[6, 5])
ax382 = plt.subplot(gs8[7, 3])
ax383 = plt.subplot(gs8[7, 4])
ax384 = plt.subplot(gs8[7, 5])
drydep_8 = [ax373,ax374,ax375,ax376,ax377,ax378,ax379,ax380,ax381,ax382,ax383,ax384]


gs9 = gridspec.GridSpec(8, 6)
gs9.update(top=0.333,bottom=0,left=0.666,right=1.0,wspace=0,hspace=0)
ax385 = plt.subplot(gs9[0, 0])
ax386 = plt.subplot(gs9[0, 1])
ax387 = plt.subplot(gs9[0, 2])
ax388 = plt.subplot(gs9[1, 0])
ax389 = plt.subplot(gs9[1, 1])
ax390 = plt.subplot(gs9[1, 2])
ax391 = plt.subplot(gs9[2, 0])
ax392 = plt.subplot(gs9[2, 1])
ax393 = plt.subplot(gs9[2, 2])
ax394 = plt.subplot(gs9[3, 0])
ax395 = plt.subplot(gs9[3, 1])
ax396 = plt.subplot(gs9[3, 2])
anmvoc_9 = [ax385,ax386,ax387,ax388,ax389,ax390,ax391,ax392,ax393,ax394,ax395,ax396]

ax397 = plt.subplot(gs9[0, 3])
ax398 = plt.subplot(gs9[0, 4])
ax399 = plt.subplot(gs9[0, 5])
ax400 = plt.subplot(gs9[1, 3])
ax401 = plt.subplot(gs9[1, 4])
ax402 = plt.subplot(gs9[1, 5])
ax403 = plt.subplot(gs9[2, 3])
ax404 = plt.subplot(gs9[2, 4])
ax405 = plt.subplot(gs9[2, 5])
ax406 = plt.subplot(gs9[3, 3])
ax407 = plt.subplot(gs9[3, 4])
ax408 = plt.subplot(gs9[3, 5])
bnmvoc_9 = [ax397,ax398,ax399,ax400,ax401,ax402,ax403,ax404,ax405,ax406,ax407,ax408]

ax409 = plt.subplot(gs9[4, 0])
ax410 = plt.subplot(gs9[4, 1])
ax411 = plt.subplot(gs9[4, 2])
ax412 = plt.subplot(gs9[5, 0])
ax413 = plt.subplot(gs9[5, 1])
ax414 = plt.subplot(gs9[5, 2])
ax415 = plt.subplot(gs9[6, 0])
ax416 = plt.subplot(gs9[6, 1])
ax417 = plt.subplot(gs9[6, 2])
ax418 = plt.subplot(gs9[7, 0])
ax419 = plt.subplot(gs9[7, 1])
ax420 = plt.subplot(gs9[7, 2])
co_9 = [ax409,ax410,ax411,ax412,ax413,ax414,ax415,ax416,ax417,ax418,ax419,ax420]

ax421 = plt.subplot(gs9[4, 3])
ax422 = plt.subplot(gs9[4, 4])
ax423 = plt.subplot(gs9[4, 5])
ax424 = plt.subplot(gs9[5, 3])
ax425 = plt.subplot(gs9[5, 4])
ax426 = plt.subplot(gs9[5, 5])
ax427 = plt.subplot(gs9[6, 3])
ax428 = plt.subplot(gs9[6, 4])
ax429 = plt.subplot(gs9[6, 5])
ax430 = plt.subplot(gs9[7, 3])
ax431 = plt.subplot(gs9[7, 4])
ax432 = plt.subplot(gs9[7, 5])
drydep_9 = [ax421,ax422,ax423,ax424,ax425,ax426,ax427,ax428,ax429,ax430,ax431,ax432]

#----------------------------------------------

#sort diff_wf in order as all_models
si = []
for mm in range(len(all_models)):
    si.append(collected_models.index(all_models[mm]))
diff_wf = diff_wf[:,si,:]
perdiff_wf = perdiff_wf[:,si,:]
model_wf = model_wf[:,si,:]
collected_models = np.array(collected_models)
collected_models = collected_models[si]

#cut diff_wf to remove duplciates
u, indices = np.unique(all_models, return_index=True)
indices = np.sort(indices)
diff_wf = diff_wf[:,indices,:]
perdiff_wf = perdiff_wf[:,indices,:]
model_wf = model_wf[:,indices,:]
collected_models = collected_models[indices]

minval_orig = 0
maxval_orig = 150
cts = np.linspace(minval_orig,maxval_orig,50)
cmaps = [matplotlib.cm.prism,matplotlib.cm.rainbow,matplotlib.cm.Set1,matplotlib.cm.cubehelix,matplotlib.cm.terrain,matplotlib.cm.hsv,matplotlib.cm.Dark2,matplotlib.cm.gist_stern,matplotlib.cm.flag]

area_label = area_labels[area]

#plot species solo
#ANMVOC
for c in range(12):     

    #get inds of all set models - deals with models that are duplicated
    valid_inds = [i for i,x in enumerate(collected_models) if ('ANMVOC' in x) or ('_STD' in x) or ('NMVOC' in x) & ('BNMVOC' not in x) & ('DRYDEP' not in x) & ('ACO' not in x) & ('CH4' not in x)]
    valid_models = [x for i,x in enumerate(collected_models) if ('ANMVOC' in x) or ('_STD' in x) or ('NMVOC' in x) & ('BNMVOC' not in x) & ('DRYDEP' not in x) & ('ACO' not in x) & ('CH4' not in x)]
    #sort valid models by original set order
    si = []
    for mm in range(len(valid_models)):
        si.append(valid_models.index(alt_model_dirs_a[mm]))
    valid_inds = np.array(valid_inds)
    valid_models = np.array(valid_models)
    valid_inds = valid_inds[si]
    valid_models = valid_models[si]
    
    area_grid_1 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_2 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_3 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_4 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_5 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_6 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_7 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_8 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_9 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    
    pl = anmvoc_1[c].contourf(area_grid_1,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[0])
    pl = anmvoc_2[c].contourf(area_grid_2,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[1])
    pl = anmvoc_3[c].contourf(area_grid_3,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[2])
    pl = anmvoc_4[c].contourf(area_grid_4,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[3])
    pl = anmvoc_5[c].contourf(area_grid_5,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[4])
    pl = anmvoc_6[c].contourf(area_grid_6,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[5])
    pl = anmvoc_7[c].contourf(area_grid_7,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[6])
    pl = anmvoc_8[c].contourf(area_grid_8,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[7])
    pl = anmvoc_9[c].contourf(area_grid_9,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[8])
        
    anmvoc_1[c].set_xlim(0,4)
    anmvoc_1[c].set_ylim(0,4)
    anmvoc_1[c].set_xticks([1,2,3,4])
    anmvoc_1[c].set_xticklabels(['','','',''])
    anmvoc_1[c].set_yticks([1,2,3,4])
    anmvoc_1[c].set_yticklabels(['','','','']) 
    anmvoc_2[c].set_xlim(0,4)
    anmvoc_2[c].set_ylim(0,4)
    anmvoc_2[c].set_xticks([1,2,3,4])
    anmvoc_2[c].set_xticklabels(['','','',''])
    anmvoc_2[c].set_yticks([1,2,3,4])
    anmvoc_2[c].set_yticklabels(['','','',''])
    anmvoc_3[c].set_xlim(0,4)
    anmvoc_3[c].set_ylim(0,4)
    anmvoc_3[c].set_xticks([1,2,3,4])
    anmvoc_3[c].set_xticklabels(['','','',''])
    anmvoc_3[c].set_yticks([1,2,3,4])
    anmvoc_3[c].set_yticklabels(['','','',''])
    anmvoc_4[c].set_xlim(0,4)
    anmvoc_4[c].set_ylim(0,4)
    anmvoc_4[c].set_xticks([1,2,3,4])
    anmvoc_4[c].set_xticklabels(['','','',''])
    anmvoc_4[c].set_yticks([1,2,3,4])
    anmvoc_4[c].set_yticklabels(['','','',''])
    anmvoc_5[c].set_xlim(0,4)
    anmvoc_5[c].set_ylim(0,4)
    anmvoc_5[c].set_xticks([1,2,3,4])
    anmvoc_5[c].set_xticklabels(['','','',''])
    anmvoc_5[c].set_yticks([1,2,3,4])
    anmvoc_5[c].set_yticklabels(['','','',''])
    anmvoc_6[c].set_xlim(0,4)
    anmvoc_6[c].set_ylim(0,4)
    anmvoc_6[c].set_xticks([1,2,3,4])
    anmvoc_6[c].set_xticklabels(['','','',''])
    anmvoc_6[c].set_yticks([1,2,3,4])
    anmvoc_6[c].set_yticklabels(['','','',''])
    anmvoc_6[c].set_xlim(0,4)
    anmvoc_6[c].set_ylim(0,4)
    anmvoc_6[c].set_xticks([1,2,3,4])
    anmvoc_6[c].set_xticklabels(['','','',''])
    anmvoc_6[c].set_yticks([1,2,3,4])
    anmvoc_6[c].set_yticklabels(['','','',''])
    anmvoc_7[c].set_xlim(0,4)
    anmvoc_7[c].set_ylim(0,4)
    anmvoc_7[c].set_xticks([1,2,3,4])
    anmvoc_7[c].set_xticklabels(['','','',''])
    anmvoc_7[c].set_yticks([1,2,3,4])
    anmvoc_7[c].set_yticklabels(['','','',''])
    anmvoc_8[c].set_xlim(0,4)
    anmvoc_8[c].set_ylim(0,4)
    anmvoc_8[c].set_xticks([1,2,3,4])
    anmvoc_8[c].set_xticklabels(['','','',''])
    anmvoc_8[c].set_yticks([1,2,3,4])
    anmvoc_8[c].set_yticklabels(['','','',''])
    anmvoc_9[c].set_xlim(0,4)
    anmvoc_9[c].set_ylim(0,4)
    anmvoc_9[c].set_xticks([1,2,3,4])
    anmvoc_9[c].set_xticklabels(['','','',''])
    anmvoc_9[c].set_yticks([1,2,3,4])
    anmvoc_9[c].set_yticklabels(['','','',''])
    
    anmvoc_1[c].set_axis_off()
    anmvoc_2[c].set_axis_off()
    anmvoc_3[c].set_axis_off()
    anmvoc_4[c].set_axis_off()
    anmvoc_5[c].set_axis_off()
    anmvoc_6[c].set_axis_off()
    anmvoc_7[c].set_axis_off()
    anmvoc_8[c].set_axis_off()
    anmvoc_9[c].set_axis_off()

#BNMVOC
for c in range(12):
    #get inds of all set models - deals with models that are duplicated
    valid_inds = [i for i,x in enumerate(collected_models) if ('BNMVOC' in x) or ('_STD' in x) or ('NMVOC' in x) & ('ANMVOC' not in x) & ('DRYDEP' not in x) & ('ACO' not in x) & ('CH4' not in x)]
    valid_models = [x for i,x in enumerate(collected_models) if ('BNMVOC' in x) or ('_STD' in x) or ('NMVOC' in x) & ('ANMVOC' not in x) & ('DRYDEP' not in x) & ('ACO' not in x) & ('CH4' not in x)]
    #sort valid models by original set order
    si = []
    for mm in range(len(valid_models)):
        si.append(valid_models.index(alt_model_dirs_b[mm]))
    valid_inds = np.array(valid_inds)
    valid_models = np.array(valid_models)
    valid_inds = valid_inds[si]
    valid_models = valid_models[si]
    
    area_grid_1 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_2 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_3 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_4 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_5 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_6 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_7 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_8 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_9 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    
    pl = bnmvoc_1[c].contourf(area_grid_1,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[0])
    pl = bnmvoc_2[c].contourf(area_grid_2,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[1])
    pl = bnmvoc_3[c].contourf(area_grid_3,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[2])
    pl = bnmvoc_4[c].contourf(area_grid_4,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[3])
    pl = bnmvoc_5[c].contourf(area_grid_5,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[4])
    pl = bnmvoc_6[c].contourf(area_grid_6,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[5])
    pl = bnmvoc_7[c].contourf(area_grid_7,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[6])
    pl = bnmvoc_8[c].contourf(area_grid_8,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[7])
    pl = bnmvoc_9[c].contourf(area_grid_9,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[8])
    
    bnmvoc_1[c].set_xlim(0,4)
    bnmvoc_1[c].set_ylim(0,4)
    bnmvoc_1[c].set_xticks([1,2,3,4])
    bnmvoc_1[c].set_xticklabels(['','','',''])
    bnmvoc_1[c].set_yticks([1,2,3,4])
    bnmvoc_1[c].set_yticklabels(['','','',''])
    bnmvoc_2[c].set_xlim(0,4)
    bnmvoc_2[c].set_ylim(0,4)
    bnmvoc_2[c].set_xticks([1,2,3,4])
    bnmvoc_2[c].set_xticklabels(['','','',''])
    bnmvoc_2[c].set_yticks([1,2,3,4])
    bnmvoc_2[c].set_yticklabels(['','','',''])
    bnmvoc_3[c].set_xlim(0,4)
    bnmvoc_3[c].set_ylim(0,4)
    bnmvoc_3[c].set_xticks([1,2,3,4])
    bnmvoc_3[c].set_xticklabels(['','','',''])
    bnmvoc_3[c].set_yticks([1,2,3,4])
    bnmvoc_3[c].set_yticklabels(['','','',''])
    bnmvoc_4[c].set_xlim(0,4)
    bnmvoc_4[c].set_ylim(0,4)
    bnmvoc_4[c].set_xticks([1,2,3,4])
    bnmvoc_4[c].set_xticklabels(['','','',''])
    bnmvoc_4[c].set_yticks([1,2,3,4])
    bnmvoc_4[c].set_yticklabels(['','','',''])
    bnmvoc_5[c].set_xlim(0,4)
    bnmvoc_5[c].set_ylim(0,4)
    bnmvoc_5[c].set_xticks([1,2,3,4])
    bnmvoc_5[c].set_xticklabels(['','','',''])
    bnmvoc_5[c].set_yticks([1,2,3,4])
    bnmvoc_5[c].set_yticklabels(['','','',''])
    bnmvoc_6[c].set_xlim(0,4)
    bnmvoc_6[c].set_ylim(0,4)
    bnmvoc_6[c].set_xticks([1,2,3,4])
    bnmvoc_6[c].set_xticklabels(['','','',''])
    bnmvoc_6[c].set_yticks([1,2,3,4])
    bnmvoc_6[c].set_yticklabels(['','','',''])
    bnmvoc_6[c].set_xlim(0,4)
    bnmvoc_6[c].set_ylim(0,4)
    bnmvoc_6[c].set_xticks([1,2,3,4])
    bnmvoc_6[c].set_xticklabels(['','','',''])
    bnmvoc_6[c].set_yticks([1,2,3,4])
    bnmvoc_6[c].set_yticklabels(['','','',''])
    bnmvoc_7[c].set_xlim(0,4)
    bnmvoc_7[c].set_ylim(0,4)
    bnmvoc_7[c].set_xticks([1,2,3,4])
    bnmvoc_7[c].set_xticklabels(['','','',''])
    bnmvoc_7[c].set_yticks([1,2,3,4])
    bnmvoc_7[c].set_yticklabels(['','','',''])
    bnmvoc_8[c].set_xlim(0,4)
    bnmvoc_8[c].set_ylim(0,4)
    bnmvoc_8[c].set_xticks([1,2,3,4])
    bnmvoc_8[c].set_xticklabels(['','','',''])
    bnmvoc_8[c].set_yticks([1,2,3,4])
    bnmvoc_8[c].set_yticklabels(['','','',''])
    bnmvoc_9[c].set_xlim(0,4)
    bnmvoc_9[c].set_ylim(0,4)
    bnmvoc_9[c].set_xticks([1,2,3,4])
    bnmvoc_9[c].set_xticklabels(['','','',''])
    bnmvoc_9[c].set_yticks([1,2,3,4])
    bnmvoc_9[c].set_yticklabels(['','','',''])
    
    bnmvoc_1[c].set_axis_off()
    bnmvoc_2[c].set_axis_off()
    bnmvoc_3[c].set_axis_off()
    bnmvoc_4[c].set_axis_off()
    bnmvoc_5[c].set_axis_off()
    bnmvoc_6[c].set_axis_off()
    bnmvoc_7[c].set_axis_off()
    bnmvoc_8[c].set_axis_off()
    bnmvoc_9[c].set_axis_off()

 
#ACO cut
    #get inds of all set models - deals with models that are duplicated
    valid_inds = [i for i,x in enumerate(collected_models) if ('ACO' in x) or ('_STD' in x) or ('NMVOC' in x) & ('ANMVOC' not in x) & ('BNMVOC' not in x) & ('DRYDEP' not in x) & ('CH4' not in x)]
    valid_models = [x for i,x in enumerate(collected_models) if ('ACO' in x) or ('_STD' in x) or ('NMVOC' in x) & ('ANMVOC' not in x) & ('BNMVOC' not in x) & ('DRYDEP' not in x) & ('CH4' not in x)]
    #sort valid models by original set order
    si = []
    for mm in range(len(valid_models)):
        si.append(valid_models.index(alt_model_dirs_d[mm]))
    valid_inds = np.array(valid_inds)
    valid_models = np.array(valid_models)
    valid_inds = valid_inds[si]
    valid_models = valid_models[si]
    
    area_grid_1 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_2 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_3 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_4 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_5 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_6 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_7 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_8 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_9 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    
    pl = co_1[c].contourf(area_grid_1,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[0])
    pl = co_2[c].contourf(area_grid_2,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[1])
    pl = co_3[c].contourf(area_grid_3,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[2])
    pl = co_4[c].contourf(area_grid_4,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[3])
    pl = co_5[c].contourf(area_grid_5,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[4])
    pl = co_6[c].contourf(area_grid_6,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[5])
    pl = co_7[c].contourf(area_grid_7,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[6])
    pl = co_8[c].contourf(area_grid_8,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[7])
    pl = co_9[c].contourf(area_grid_9,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[8])
    
    co_1[c].set_xlim(0,4)
    co_1[c].set_ylim(0,4)
    co_1[c].set_xticks([1,2,3,4])
    co_1[c].set_xticklabels(['','','',''])
    co_1[c].set_yticks([1,2,3,4])
    co_1[c].set_yticklabels(['','','',''])  
    co_2[c].set_xlim(0,4)
    co_2[c].set_ylim(0,4)
    co_2[c].set_xticks([1,2,3,4])
    co_2[c].set_xticklabels(['','','',''])
    co_2[c].set_yticks([1,2,3,4])
    co_2[c].set_yticklabels(['','','',''])
    co_3[c].set_xlim(0,4)
    co_3[c].set_ylim(0,4)
    co_3[c].set_xticks([1,2,3,4])
    co_3[c].set_xticklabels(['','','',''])
    co_3[c].set_yticks([1,2,3,4])
    co_3[c].set_yticklabels(['','','',''])
    co_4[c].set_xlim(0,4)
    co_4[c].set_ylim(0,4)
    co_4[c].set_xticks([1,2,3,4])
    co_4[c].set_xticklabels(['','','',''])
    co_4[c].set_yticks([1,2,3,4])
    co_4[c].set_yticklabels(['','','',''])
    co_5[c].set_xlim(0,4)
    co_5[c].set_ylim(0,4)
    co_5[c].set_xticks([1,2,3,4])
    co_5[c].set_xticklabels(['','','',''])
    co_5[c].set_yticks([1,2,3,4])
    co_5[c].set_yticklabels(['','','',''])
    co_6[c].set_xlim(0,4)
    co_6[c].set_ylim(0,4)
    co_6[c].set_xticks([1,2,3,4])
    co_6[c].set_xticklabels(['','','',''])
    co_6[c].set_yticks([1,2,3,4])
    co_6[c].set_yticklabels(['','','',''])
    co_6[c].set_xlim(0,4)
    co_6[c].set_ylim(0,4)
    co_6[c].set_xticks([1,2,3,4])
    co_6[c].set_xticklabels(['','','',''])
    co_6[c].set_yticks([1,2,3,4])
    co_6[c].set_yticklabels(['','','',''])
    co_7[c].set_xlim(0,4)
    co_7[c].set_ylim(0,4)
    co_7[c].set_xticks([1,2,3,4])
    co_7[c].set_xticklabels(['','','',''])
    co_7[c].set_yticks([1,2,3,4])
    co_7[c].set_yticklabels(['','','',''])
    co_8[c].set_xlim(0,4)
    co_8[c].set_ylim(0,4)
    co_8[c].set_xticks([1,2,3,4])
    co_8[c].set_xticklabels(['','','',''])
    co_8[c].set_yticks([1,2,3,4])
    co_8[c].set_yticklabels(['','','',''])
    co_9[c].set_xlim(0,4)
    co_9[c].set_ylim(0,4)
    co_9[c].set_xticks([1,2,3,4])
    co_9[c].set_xticklabels(['','','',''])
    co_9[c].set_yticks([1,2,3,4])
    co_9[c].set_yticklabels(['','','',''])
    
    co_1[c].set_axis_off()
    co_2[c].set_axis_off()
    co_3[c].set_axis_off()
    co_4[c].set_axis_off()
    co_5[c].set_axis_off()
    co_6[c].set_axis_off()
    co_7[c].set_axis_off()
    co_8[c].set_axis_off()
    co_9[c].set_axis_off()
    
#drydepo3 cut
for c in range(12):
    #get inds of all set models - deals with models that are duplicated
    valid_inds = [i for i,x in enumerate(collected_models) if ('DRYDEP' in x) or ('_STD' in x) or ('NMVOC' in x) & ('ANMVOC' not in x) & ('BNMVOC' not in x) & ('ACO' not in x) & ('CH4' not in x)]
    valid_models = [x for i,x in enumerate(collected_models) if ('DRYDEP' in x) or ('_STD' in x) or ('NMVOC' in x) & ('ANMVOC' not in x) & ('BNMVOC' not in x) & ('ACO' not in x) & ('CH4' not in x)]
    #sort valid models by original set order
    si = []
    for mm in range(len(valid_models)):
        si.append(valid_models.index(alt_model_dirs_c[mm]))
    valid_inds = np.array(valid_inds)
    valid_models = np.array(valid_models)
    valid_inds = valid_inds[si]
    valid_models = valid_models[si]
    
    area_grid_1 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_2 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_3 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_4 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_5 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_6 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_7 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_8 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    area_grid_9 = np.ma.array(np.reshape(perdiff_wf[0,valid_inds,c],(5,5)),mask=np.isnan(np.reshape(perdiff_wf[0,valid_inds,c],(5,5))))
    
    pl = drydep_1[c].contourf(area_grid_1,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[0])
    pl = drydep_2[c].contourf(area_grid_2,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[1])
    pl = drydep_3[c].contourf(area_grid_3,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[2])
    pl = drydep_4[c].contourf(area_grid_4,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[3])
    pl = drydep_5[c].contourf(area_grid_5,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[4])
    pl = drydep_6[c].contourf(area_grid_6,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[5])
    pl = drydep_7[c].contourf(area_grid_7,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[6])
    pl = drydep_8[c].contourf(area_grid_8,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[7])
    pl = drydep_9[c].contourf(area_grid_9,cts,vmin = minval_orig,vmax=maxval_orig,cmap =cmaps[8])

    drydep_1[c].set_xlim(0,4)
    drydep_1[c].set_ylim(0,4)
    drydep_1[c].set_xticks([1,2,3,4])
    drydep_1[c].set_xticklabels(['','','',''])
    drydep_1[c].set_yticks([1,2,3,4])
    drydep_1[c].set_yticklabels(['','','',''])  
    drydep_2[c].set_xlim(0,4)
    drydep_2[c].set_ylim(0,4)
    drydep_2[c].set_xticks([1,2,3,4])
    drydep_2[c].set_xticklabels(['','','',''])
    drydep_2[c].set_yticks([1,2,3,4])
    drydep_2[c].set_yticklabels(['','','','']) 
    drydep_3[c].set_xlim(0,4)
    drydep_3[c].set_ylim(0,4)
    drydep_3[c].set_xticks([1,2,3,4])
    drydep_3[c].set_xticklabels(['','','',''])
    drydep_3[c].set_yticks([1,2,3,4])
    drydep_3[c].set_yticklabels(['','','','']) 
    drydep_4[c].set_xlim(0,4)
    drydep_4[c].set_ylim(0,4)
    drydep_4[c].set_xticks([1,2,3,4])
    drydep_4[c].set_xticklabels(['','','',''])
    drydep_4[c].set_yticks([1,2,3,4])
    drydep_4[c].set_yticklabels(['','','',''])
    drydep_5[c].set_xlim(0,4)
    drydep_5[c].set_ylim(0,4)
    drydep_5[c].set_xticks([1,2,3,4])
    drydep_5[c].set_xticklabels(['','','',''])
    drydep_5[c].set_yticks([1,2,3,4])
    drydep_5[c].set_yticklabels(['','','',''])
    drydep_6[c].set_xlim(0,4)
    drydep_6[c].set_ylim(0,4)
    drydep_6[c].set_xticks([1,2,3,4])
    drydep_6[c].set_xticklabels(['','','',''])
    drydep_6[c].set_yticks([1,2,3,4])
    drydep_6[c].set_yticklabels(['','','',''])
    drydep_6[c].set_xlim(0,4)
    drydep_6[c].set_ylim(0,4)
    drydep_6[c].set_xticks([1,2,3,4])
    drydep_6[c].set_xticklabels(['','','',''])
    drydep_6[c].set_yticks([1,2,3,4])
    drydep_6[c].set_yticklabels(['','','',''])
    drydep_7[c].set_xlim(0,4)
    drydep_7[c].set_ylim(0,4)
    drydep_7[c].set_xticks([1,2,3,4])
    drydep_7[c].set_xticklabels(['','','',''])
    drydep_7[c].set_yticks([1,2,3,4])
    drydep_7[c].set_yticklabels(['','','',''])
    drydep_8[c].set_xlim(0,4)
    drydep_8[c].set_ylim(0,4)
    drydep_8[c].set_xticks([1,2,3,4])
    drydep_8[c].set_xticklabels(['','','',''])
    drydep_8[c].set_yticks([1,2,3,4])
    drydep_8[c].set_yticklabels(['','','',''])
    drydep_9[c].set_xlim(0,4)
    drydep_9[c].set_ylim(0,4)
    drydep_9[c].set_xticks([1,2,3,4])
    drydep_9[c].set_xticklabels(['','','',''])
    drydep_9[c].set_yticks([1,2,3,4])
    drydep_9[c].set_yticklabels(['','','',''])
    
    drydep_1[c].set_axis_off()
    drydep_2[c].set_axis_off()
    drydep_3[c].set_axis_off()
    drydep_4[c].set_axis_off()
    drydep_5[c].set_axis_off()
    drydep_6[c].set_axis_off()
    drydep_7[c].set_axis_off()
    drydep_8[c].set_axis_off()
    drydep_9[c].set_axis_off()

plt.savefig('art_in_science.png')
#plt.show()
