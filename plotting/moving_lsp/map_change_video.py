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
import cmath
import matplotlib.animation as animation

#'', d or s
period_type = 's'
#average, amplitude, phase
output_type = 'phase'
#all, rural, urban or highalt
output_set = 'highalt'

if output_type == 'amplitude':
    min = 0
    max = 30
    cm=plt.cm.jet
elif output_type == 'average':
    min = 0
    max = 60
    cm=plt.cm.jet
elif output_type == 'phase':
    if period_type == 'd':
        min = 0
        max = 24
    if period_type == 's':
        min = 0
        max = 12
    cm=plt.cm.hsv

ts_files =  glob.glob('*HN.nc')
ts_files = ts_files[3:]

present_dir = os.getcwd()
paths = present_dir.split('/')
species = paths[-2]

#set up plot
latlower_setup = [20,30,20,-90]
latupper_setup = [80,72,55,90]
lonwest_setup = [-170,-15,115,-180]
loneast_setup = [-50,35,155,180]
label_out = ['NA','EU','EA','ZZ']
label = ['NA','EU','EA','ROW']

#set colorbar 


fig =plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')

#plt.tight_layout(pad = 3.08)
fig.subplots_adjust(bottom=0.16)
cbar_ax = fig.add_axes([0.10, 0.08, 0.80, 0.06])
fig.subplots_adjust(top=0.94)

title = fig.suptitle('',fontsize=30)

def data_gen(f):
    print f

    ax = plt.subplot2grid((2,2), (0,0), colspan=2, rowspan=2)
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=1)
    ax2 = plt.subplot2grid((2,2), (0,1), colspan=1)
    ax3 = plt.subplot2grid((2,2), (1,0), colspan=1)
    ax4 = plt.subplot2grid((2,2), (1,1), colspan=1)
    ax_list = [ax1,ax2,ax3,ax4]

    a = Basemap(projection='cyl',llcrnrlat=latlower_setup[0],urcrnrlat=latupper_setup[0],llcrnrlon=lonwest_setup[0],urcrnrlon=loneast_setup[0],resolution='c',ax = ax1)
    b = Basemap(projection='cyl',llcrnrlat=latlower_setup[1],urcrnrlat=latupper_setup[1],llcrnrlon=lonwest_setup[1],urcrnrlon=loneast_setup[1],resolution='c',ax = ax2)
    c = Basemap(projection='cyl',llcrnrlat=latlower_setup[2],urcrnrlat=latupper_setup[2],llcrnrlon=lonwest_setup[2],urcrnrlon=loneast_setup[2],resolution='c',ax = ax3)
    d = Basemap(projection='cyl',llcrnrlat=latlower_setup[3],urcrnrlat=latupper_setup[3],llcrnrlon=lonwest_setup[3],urcrnrlon=loneast_setup[3],resolution='c',ax = ax4)

    a.drawcoastlines()
    a.drawmapboundary()
    b.drawcoastlines()
    b.drawmapboundary()
    c.drawcoastlines()
    c.drawmapboundary()
    d.drawcoastlines()
    d.drawmapboundary()

    count = 0
    for ax in ax_list:
        if (count == 0) or (count == 3):
            ax.text(0.03, 0.18, label[count], transform=ax.transAxes,fontsize=30, va='top')
        if (count == 1) or (count == 2):
            ax.text(0.03, 0.97, label[count], transform=ax.transAxes,fontsize=30, va='top')
        count+=1

    obs_ts_fname = ts_files[f]
    last_split = obs_ts_fname.split('/')[-1]
    start_year = last_split.split('_')[3] 
    end_year = last_split.split('_')[4]
    obs_lsp_fname = 'obs_sig_periods_%s_%s.nc'%(start_year,end_year)
    print ts_files[f]
    print obs_lsp_fname
    print start_year,end_year
    
    title.set_text('%s - %s'%(start_year,end_year))

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

    #print obs_refs
    #print bad_inds
        
    bad_inds = np.array(list(set(bad_inds)))
    obs_refs = np.delete(obs_refs,bad_inds)
    obs_lats = np.delete(obs_lats,bad_inds)
    obs_lons = np.delete(obs_lons,bad_inds)
    
    #print obs_refs
    
    if len(obs_refs) == 0:
        return

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
    
    if output_type == 'amplitude':
        if period_type == 'd':
            param = obs_daily_mag
        elif period_type == 's':    
            param = obs_seasonal_mag
    if output_type == 'phase':
        if period_type == 'd':
            param = obs_daily_phase
        elif period_type == 's':    
            param = obs_seasonal_phase
    if output_type == 'average':
        param = obs_ave

    #plot
    #---------------------------------------
    count = 0

    for ax in ax_list:

        if count == 3:
            test = []
            other_tags = ['AF','ANT','ARC','OC','O','SA','MA']
            for i in tags:
                if i in other_tags:
                    test.append(True)
                else:
                    test.append(False)
            test = np.array(test)
        else:
            current_tag = label_out[count] 
            test = tags == current_tag
    
        current_lons = obs_lons[test]
        current_lats = obs_lats[test]
        current_refs = obs_refs[test]
        current_param = param[test]

        if count == 0:
            if int(start_year) < 1980:
                m_size = 200
            else:
                m_size= 50
        if count == 1:
            if int(start_year) < 1986:
                m_size = 200
            else:
                m_size = 30
        if count == 2:
            m_size = 250
        if count == 3:
            m_size = 300
        
        if count == 0:
            X,Y = a(current_lons,current_lats)
            for i in range(len(current_lons)):
                all = a.scatter(X[i],Y[i],c=current_param[i], s=m_size, vmin = min,vmax = max, marker = 'o',edgecolor='black',linewidth=0.5,cmap=cm,zorder=10)
        
        elif count == 1:
            X,Y = b(current_lons,current_lats)
            for i in range(len(current_lons)):
                all = b.scatter(X[i],Y[i],c=current_param[i], s=m_size, vmin = min,vmax = max, marker = 'o',edgecolor='black',linewidth=0.5,cmap=cm,zorder=10)
        
        elif count == 2:
            X,Y = c(current_lons,current_lats)
            for i in range(len(current_lons)):
                all = c.scatter(X[i],Y[i],c=current_param[i], s=m_size, vmin = min,vmax = max, marker = 'o',edgecolor='black',linewidth=0.5,cmap=cm,zorder=10)
        
        elif count == 3:
            X,Y = d(current_lons,current_lats)
            for i in range(len(current_lons)):
                all = d.scatter(X[i],Y[i],c=current_param[i], s=m_size, vmin = min,vmax = max, marker = 'o',edgecolor='black',linewidth=0.5,cmap=cm,zorder=10)
        
        count+=1
    
    cb = fig.colorbar(all, cax=cbar_ax,orientation ='horizontal')
    if output_type == 'phase':
        if period_type == 'd':
            cb.set_label('Hour', fontsize = 24)
        elif period_type == 's':
            cb.set_label('Month', fontsize = 24)
    elif output_type == 'amplitude':
        cb.set_label('Amplitude (ppb)', fontsize = 24)
    elif output_type == 'average':
        cb.set_label('Average (ppb)', fontsize = 24)    
        
    cb.ax.tick_params(labelsize=22)
    
ani = animation.FuncAnimation(fig, data_gen, range(len(ts_files)),blit=False)

ani.save('obs_map_phase_change_%s_%s_%s.mp4'%(output_type,period_type,output_set),'ffmpeg',fps=1,extra_args=['-vcodec','libx264', '-preset', 'ultrafast', '-s','hd1080', '-pix_fmt', 'yuv420p'])

#plt.show()
