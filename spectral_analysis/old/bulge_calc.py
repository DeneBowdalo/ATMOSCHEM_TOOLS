#import matplotlib
#matplotlib.use('Agg')
import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from mpl_toolkits.basemap import Basemap, shiftgrid
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

def form2(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.2f' % x	

def form5(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.7f' % x

xformatter = FuncFormatter(form2)
yformatter = FuncFormatter(form5)

#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres = modules.get_obs_info(present_dir)
model_fname,species,start_year,end_year = modules.get_model_info(present_dir)

obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(np.copy(obs_refs))

#--------------------------------------------------------
#load in periodic lsp data
obs_period_grp = Dataset('../obs_%s_%s/obs_sig_periods.nc'%(vres,timeres))
model_period_grp = Dataset('model_sig_periods.nc')

obs_daily_waveforms = []
obs_seasonal_waveforms = []
model_daily_waveforms = []
model_seasonal_waveforms = []
obs_seasonal_amp = []
model_seasonal_amp = []

for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    site_group_mod = model_period_grp.groups[ref]

    obs_daily_waveforms.append(site_group_obs.variables['daily_waveform'][:])
    obs_seasonal_waveforms.append(site_group_obs.variables['seasonal_waveform'][:])
    obs_seasonal_amp.append(site_group_obs.seasonal_amplitude)
    
    model_daily_waveforms.append(site_group_mod.variables['daily_waveform'][:])
    model_seasonal_waveforms.append(site_group_mod.variables['seasonal_waveform'][:])
    model_seasonal_amp.append(site_group_mod.seasonal_amplitude)


obs_seasonal_waveforms = np.array(obs_seasonal_waveforms)
model_seasonal_waveforms = np.array(model_seasonal_waveforms)

#for i in range(len(obs_seasonal_waveforms)):
#    obs_seasonal_waveforms[i] = obs_seasonal_waveforms[i] -  np.average((obs_seasonal_amp[i],model_seasonal_amp[i]))
#    model_seasonal_waveforms[i] = model_seasonal_waveforms[i] -  np.average((obs_seasonal_waveforms[i],model_seasonal_waveforms[i]))

max_ph = 12.
modules.bulge_calc(obs_refs,obs_seasonal_waveforms,model_seasonal_waveforms,max_ph)
    

