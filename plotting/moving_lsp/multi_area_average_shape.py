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

present_dir = os.getcwd()
paths = present_dir.split('/')
species = paths[-2]

ts_files =  glob.glob('*HN.nc')
ts_files = np.array(ts_files[3:])
i = np.array([19,20,21,37])
ts_files = ts_files[i]

fig, axes = plt.subplots(nrows=5, ncols=5,figsize=(24,18))
fig.patch.set_facecolor('white')

obsn_linspace = np.linspace(0, 1, len(ts_files))

start = datetime.datetime(year = 2004, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = 2005, month = 1, day = 1, hour = 0, minute = 0)
grid_dates = pd.date_range(start,end,freq='H')[:8766]


period_type = 's'
output_set = 'rural'

areas = ['ANT','S_O','OC','SA','AF','SW_NA','NW_NA','NE_NA','CE_NA','SE_NA','C_NA','S_EU','SW_EU','CW_EU','CS_EU','C_EU','E_EU','NW_EU','N_EU','SE_AS','NE_AS','N_O','AL']#,'ARC']

area_boundaries,area_tags,area_labels = modules.area_dicts()

for c in range(len(ts_files)):
    #------------------------------------------------------------
    obs_ts_fname = ts_files[c]
    last_split = obs_ts_fname.split('/')[-1]
    start_year = last_split.split('_')[3] 
    end_year = last_split.split('_')[4]
    obs_lsp_fname = 'obs_sig_periods_%s_%s.nc'%(start_year,end_year)
    print ts_files[c]
    print obs_lsp_fname
    print start_year,end_year
    
    obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_raw_class,obs_anthrome_class,obs_gap_inds = modules.read_obs_all(obs_ts_fname,species,start_year,end_year)

    #screen data for urban/high altitude sites
    
    #if type is urban then just keep low altitude urban sites 
    #if type is altitude then just keep rural high altitude sites
    #if type is rural then just keep low altitude rural sites
    #if type is all then do not checks
    
    bad_inds = []
    for x in range(len(obs_process_groups)):
        raw_class = obs_raw_class[x]
        process_group = obs_process_groups[x]
        obs_a = obs_alt[x]
        anthrome_class = obs_anthrome_class[x]
        
        #urban check
        if (output_set == 'rural') or (output_set == 'highalt') or (output_set == 'urban'):
            if (process_group == 'AirBase') or (process_group == 'EPA AQS') or (process_group == 'CAPMON') or (process_group == 'CASTNET') or (process_group == 'EANET') or (process_group == 'SEARCH'):
                if ('urban' in raw_class.lower()) or ('traffic' in raw_class.lower()) or ('industrial' in raw_class.lower()):
                    if (output_set == 'rural') or (output_set == 'highalt'):
                        bad_inds.append(x)
                else:
                    if output_set == 'urban':
                        bad_inds.append(x)
                
            elif (process_group == 'NAPS'):
                if ('i' == raw_class.lower()) or ('c' == raw_class.lower()):
                    if (output_set == 'rural') or (output_set == 'highalt'):
                        bad_inds.append(x)
                    
                else:
                    if output_set == 'urban':
                        bad_inds.append(x) 
        

            if anthrome_class == 'Dense Settlement':
                if (output_set == 'rural') or (output_set == 'highalt'):
                    bad_inds.append(x)
                
            else:
                if output_set == 'urban':
                    bad_inds.append(x)
        
        #altitude check
        if (output_set == 'urban') or (output_set == 'rural'):
            if int(obs_a) >= 1000:
                bad_inds.append(x)

        elif (output_set == 'highalt'):
            if int(obs_a) < 1000:
                bad_inds.append(x)

    bad_inds = np.array(list(set(bad_inds)))
    obs_refs = np.delete(obs_refs,bad_inds)
    obs_lats = np.delete(obs_lats,bad_inds)
    obs_lons = np.delete(obs_lons,bad_inds)

    print len(obs_refs)
    
    if len(obs_refs) == 0:
        continue

    tags = modules.get_tags(obs_refs)

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
    ,obs_periodic_variance_seasonal_night,obs_periodic_variance_all = modules.get_obs_specific(obs_lsp_fname,obs_refs)

    if period_type == 'd':
        waveform = np.array(obs_daily_waveform)
    elif period_type == 's':    
        waveform = np.array(obs_seasonal_waveform)

    #take average of waveform by area
    count = 0
    for ax in axes.flat:
    
        area = areas[count]
        print area
        area_grid = area_boundaries[area]
        area_tag = area_tags[area]
        area_label = area_labels[area]
        
        cut_test = modules.area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag)
        
        if np.all(cut_test == False):
            count+=1
            continue

        ave_waveform = np.nanmean(waveform[cut_test,:],axis=0)

        #plot phase
        ax.plot_date(grid_dates.to_pydatetime(),ave_waveform,color = plt.cm.jet(obsn_linspace[c]),linestyle='_',linewidth=1,markeredgecolor='None',markersize=1)

        if period_type == 'd':
            ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
        if period_type == 's':
            ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))    

        ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=14)
        
        count+=1
        
        
#plt.savefig('changing_%s_%s_%s.png'%(output_type,period_type,output_set))
plt.show()
