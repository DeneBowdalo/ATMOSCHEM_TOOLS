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

present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-4]

model_fname = '/work/home/db876/plotting_tools/model_files/MIROCCHEM_SURFACE_2000_2012_*_*_*_H_*.nc'
obs_fname = '/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_%s_2005_2010_H_PERIODIC.nc'%(species)

start_year = 2005
end_year = 2010

obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

#--------------------------------------------------------
#load in periodic lsp data
nmodels = 7
nmodels_linspace = np.linspace(0, 1, nmodels)

obs_period_grp = Dataset('../obs_SURFACE_H/obs_sig_periods.nc')
cesmcam_period_grp = Dataset('../CESMCAM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
cmam_period_grp = Dataset('../CMAM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
geosccm_period_grp = Dataset('../GEOSCCM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
#geoschemv90103_period_grp = Dataset('../../../2005_2012/2005_2010/GEOSCHEM_SURFACE_v90103_2x2.5_GEOS5_H_*/model_sig_periods.nc')
#geoschemv902_period_grp = Dataset('../../../2005_2012/2005_2010/GEOSCHEM_SURFACE_v902_2x2.5_GEOS5_H_*/model_sig_periods.nc')
geoschemv1001_period_grp = Dataset('../../../2005_2012/2005_2010/GEOSCHEM_SURFACE_v1001_4x5_GEOS5_H_*/model_sig_periods.nc')
gfdlam3_period_grp = Dataset('../GFDLAM3_SURFACE_*_*_*_H_*/model_sig_periods.nc')
gisse2r_period_grp = Dataset('../GISSE2R_SURFACE_*_*_*_H_*/model_sig_periods.nc')
mirocchem_period_grp = Dataset('../MIROCCHEM_SURFACE_*_*_*_H_*/model_sig_periods.nc') 
#mrisd_period_grp = Dataset('../MRI_SURFACE_*_*_*_D_SD/model_sig_periods.nc')
#socol3_period_grp = Dataset('../SOCOL3_SURFACE_*_*_*_D_*/model_sig_periods.nc') 
#ukca_period_grp = Dataset('../UKCA_SURFACE_*_*_*_D_*/model_sig_periods.nc')

obs_daily_waveforms = []
obs_seasonal_waveforms = []
cesmcam_daily_waveforms = []
cesmcam_seasonal_waveforms = []
cmam_daily_waveforms = []
cmam_seasonal_waveforms = []
geosccm_daily_waveforms = []
geosccm_seasonal_waveforms = []
geoschemv90103_daily_waveforms = []
geoschemv90103_seasonal_waveforms = []
geoschemv902_daily_waveforms = []
geoschemv902_seasonal_waveforms = []
geoschemv1001_daily_waveforms = []
geoschemv1001_seasonal_waveforms = []
gfdlam3_daily_waveforms = []
gfdlam3_seasonal_waveforms = []
gisse2r_daily_waveforms = []
gisse2r_seasonal_waveforms = []
mirocchem_daily_waveforms = []
mirocchem_seasonal_waveforms = []
mrisd_seasonal_waveforms = []
socol3_seasonal_waveforms = []
ukca_seasonal_waveforms = []

for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    site_group_cesmcam = cesmcam_period_grp.groups[ref]
    site_group_cmam = cmam_period_grp.groups[ref]    
    site_group_geosccm = geosccm_period_grp.groups[ref]
    #site_group_geoschemv90103 = geoschemv90103_period_grp.groups[ref]
    #site_group_geoschemv902 = geoschemv902_period_grp.groups[ref]
    site_group_geoschemv1001 = geoschemv1001_period_grp.groups[ref]
    site_group_gfdlam3 = gfdlam3_period_grp.groups[ref]
    site_group_gisse2r = gisse2r_period_grp.groups[ref]
    site_group_mirocchem = mirocchem_period_grp.groups[ref]
    #site_group_mrisd = mrisd_period_grp.groups[ref]
    #site_group_socol3 = socol3_period_grp.groups[ref]     
    #site_group_ukca = ukca_period_grp.groups[ref]

    cesmcam_daily_waveforms.append(site_group_cesmcam.variables['daily_waveform'][:])
    cesmcam_seasonal_waveforms.append(site_group_cesmcam.variables['seasonal_waveform'][:])

    cmam_daily_waveforms.append(site_group_cmam.variables['daily_waveform'][:])
    cmam_seasonal_waveforms.append(site_group_cmam.variables['seasonal_waveform'][:])
    
    geosccm_daily_waveforms.append(site_group_geosccm.variables['daily_waveform'][:])                                                                                              
    geosccm_seasonal_waveforms.append(site_group_geosccm.variables['seasonal_waveform'][:])

    #geoschemv90103_daily_waveforms.append(site_group_geoschemv90103.variables['daily_waveform'][:])                                                                                        
    #geoschemv90103_seasonal_waveforms.append(site_group_geoschemv90103.variables['seasonal_waveform'][:])

    #geoschemv902_daily_waveforms.append(site_group_geoschemv902.variables['daily_waveform'][:])                                                                                
    #geoschemv902_seasonal_waveforms.append(site_group_geoschemv902.variables['seasonal_waveform'][:])

    geoschemv1001_daily_waveforms.append(site_group_geoschemv1001.variables['daily_waveform'][:])                                                                             
    geoschemv1001_seasonal_waveforms.append(site_group_geoschemv1001.variables['seasonal_waveform'][:])

    gfdlam3_daily_waveforms.append(site_group_gfdlam3.variables['daily_waveform'][:])
    gfdlam3_seasonal_waveforms.append(site_group_gfdlam3.variables['seasonal_waveform'][:])

    gisse2r_daily_waveforms.append(site_group_gisse2r.variables['daily_waveform'][:])
    gisse2r_seasonal_waveforms.append(site_group_gisse2r.variables['seasonal_waveform'][:])

    mirocchem_daily_waveforms.append(site_group_mirocchem.variables['daily_waveform'][:])                                                                                        
    mirocchem_seasonal_waveforms.append(site_group_mirocchem.variables['seasonal_waveform'][:])

    #mrisd_seasonal_waveforms.append(site_group_mrisd.variables['seasonal_waveform'][:])

    #socol3_seasonal_waveforms.append(site_group_socol3.variables['seasonal_waveform'][:])

    #ukca_seasonal_waveforms.append(site_group_ukca.variables['seasonal_waveform'][:])

    obs_daily_waveforms.append(site_group_obs.variables['daily_waveform'][:])
    obs_seasonal_waveforms.append(site_group_obs.variables['seasonal_waveform'][:])

obs_daily_waveforms = np.array(obs_daily_waveforms)
obs_seasonal_waveforms = np.array(obs_seasonal_waveforms)
cesmcam_daily_waveforms = np.array(cesmcam_daily_waveforms)
cesmcam_seasonal_waveforms = np.array(cesmcam_seasonal_waveforms)
cmam_daily_waveforms = np.array(cmam_daily_waveforms)                                                                                                                  
cmam_seasonal_waveforms = np.array(cesmcam_seasonal_waveforms)
geosccm_daily_waveforms = np.array(geosccm_daily_waveforms)                                                                                                                  
geosccm_seasonal_waveforms = np.array(geosccm_seasonal_waveforms)
#geoschemv90103_daily_waveforms = np.array(geoschemv90103_daily_waveforms)                                                                                                                  
#geoschemv90103_seasonal_waveforms = np.array(geoschemv90103_seasonal_waveforms)                                                                                                            
#geoschemv902_daily_waveforms = np.array(geoschemv902_daily_waveforms)                                                                                                          
#geoschemv902_seasonal_waveforms = np.array(geoschemv902_seasonal_waveforms)                                                                                                    
geoschemv1001_daily_waveforms = np.array(geoschemv1001_daily_waveforms)                                                                                                       
geoschemv1001_seasonal_waveforms = np.array(geoschemv1001_seasonal_waveforms)
gfdlam3_daily_waveforms = np.array(gfdlam3_daily_waveforms)                                                                                                                
gfdlam3_seasonal_waveforms = np.array(gfdlam3_seasonal_waveforms)                                                                                                          
gisse2r_daily_waveforms = np.array(gisse2r_daily_waveforms)                                                                                                                
gisse2r_seasonal_waveforms = np.array(gisse2r_seasonal_waveforms)
mirocchem_daily_waveforms = np.array(mirocchem_daily_waveforms)                                                                                                                
mirocchem_seasonal_waveforms = np.array(mirocchem_seasonal_waveforms)                                                                                                          
#mrisd_seasonal_waveforms = np.array(mrisd_seasonal_waveforms) 
#socol3_seasonal_waveforms = np.array(socol3_seasonal_waveforms) 
#ukca_seasonal_waveforms = np.array(ukca_seasonal_waveforms)

test = obs_daily_waveforms < 0
obs_daily_waveforms[test] = np.nan
test = obs_seasonal_waveforms < 0
obs_seasonal_waveforms[test] = np.nan

#-----------------------------------
#get area
areas = ['SW_NA','NW_NA','NE_NA','CE_NA','SE_NA','C_NA','S_EU','SW_EU','CW_EU','CS_EU','C_EU','E_EU','NW_EU','N_EU']

plot_type = raw_input('\nd or s?\n')

if plot_type == 'd':
    obs_datetimes = obs_datetime_time[:24]
    model_datetimes = model_datetime_time[:24]
if plot_type == 's':
    obs_datetimes = obs_datetime_time[:8766]
    model_datetimes = model_datetime_time[:8766]

obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')

area_boundaries,area_tags,area_labels = modules.area_dicts()

#fig, axes = plt.subplots(nrows=5, ncols=5,figsize=(19,13))
fig, axes = plt.subplots(nrows=4, ncols=4,figsize=(19,13))
fig.patch.set_facecolor('white')
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

    obs_d_w = np.nanmean(obs_daily_waveforms[cut_test,:],axis=0)
    obs_s_w = np.nanmean(obs_seasonal_waveforms[cut_test,:],axis=0)

    cesmcam_d_w = np.average(cesmcam_daily_waveforms[cut_test,:],axis=0)
    cesmcam_s_w = np.average(cesmcam_seasonal_waveforms[cut_test,:],axis=0)

    cmam_d_w = np.average(cmam_daily_waveforms[cut_test,:],axis=0)
    cmam_s_w = np.average(cmam_seasonal_waveforms[cut_test,:],axis=0)

    geosccm_d_w = np.average(geosccm_daily_waveforms[cut_test,:],axis=0)                                                                                                           
    geosccm_s_w = np.average(geosccm_seasonal_waveforms[cut_test,:],axis=0)

    #geoschemv90103_d_w = np.average(geoschemv90103_daily_waveforms[cut_test,:],axis=0)                                                                                                     
    #geoschemv90103_s_w = np.average(geoschemv90103_seasonal_waveforms[cut_test,:],axis=0)

    #geoschemv902_d_w = np.average(geoschemv902_daily_waveforms[cut_test,:],axis=0)                                                                                             
    #geoschemv902_s_w = np.average(geoschemv902_seasonal_waveforms[cut_test,:],axis=0)

    geoschemv1001_d_w = np.average(geoschemv1001_daily_waveforms[cut_test,:],axis=0)                                                                                          
    geoschemv1001_s_w = np.average(geoschemv1001_seasonal_waveforms[cut_test,:],axis=0)

    gfdlam3_d_w = np.average(gfdlam3_daily_waveforms[cut_test,:],axis=0)                                                                                                   
    gfdlam3_s_w = np.average(gfdlam3_seasonal_waveforms[cut_test,:],axis=0)

    gisse2r_d_w = np.average(gisse2r_daily_waveforms[cut_test,:],axis=0)                                                                                                   
    gisse2r_s_w = np.average(gisse2r_seasonal_waveforms[cut_test,:],axis=0)

    mirocchem_d_w = np.average(mirocchem_daily_waveforms[cut_test,:],axis=0)                                                                                                   
    mirocchem_s_w = np.average(mirocchem_seasonal_waveforms[cut_test,:],axis=0)

    #mrisd_s_w = np.average(mrisd_seasonal_waveforms[cut_test,:],axis=0)

    #socol3_s_w = np.average(socol3_seasonal_waveforms[cut_test,:],axis=0) 

    #ukca_s_w = np.average(ukca_seasonal_waveforms[cut_test,:],axis=0) 

    if plot_type == 'd':
        ave_obs_waveform = obs_d_w
        ave_cesmcam_waveform = cesmcam_d_w
        ave_cmam_waveform = cmam_d_w
        ave_geosccm_waveform = geosccm_d_w    
        #ave_geoschemv90103_waveform = geoschemv90103_d_w 
        #ave_geoschemv902_waveform = geoschemv902_d_w
        ave_geoschemv1001_waveform = geoschemv1001_d_w
        ave_gfdlam3_waveform = gfdlam3_d_w 
        ave_gisse2r_waveform = gisse2r_d_w 
        ave_mirocchem_waveform = mirocchem_d_w

    if plot_type == 's':
        ave_obs_waveform = obs_s_w
        ave_cesmcam_waveform = cesmcam_s_w
        ave_cmam_waveform = cmam_s_w
        ave_geosccm_waveform = geosccm_s_w    
        #ave_geoschemv90103_waveform = geoschemv90103_s_w 
        #ave_geoschemv902_waveform = geoschemv902_s_w 
        ave_geoschemv1001_waveform = geoschemv1001_s_w
        ave_gfdlam3_waveform = gfdlam3_s_w 
        ave_gisse2r_waveform = gisse2r_s_w 
        ave_mirocchem_waveform = mirocchem_s_w
        #ave_mrisd_waveform = mrisd_s_w
        #ave_socol3_waveform = socol3_s_w
        #ave_ukca_waveform = ukca_s_w

    ave_accmip_waveform = np.nanmean((ave_cesmcam_waveform,ave_cmam_waveform,ave_geosccm_waveform,ave_gisse2r_waveform,ave_mirocchem_waveform),axis=0)   
    
    max_accmip_waveform = []
    min_accmip_waveform = []
    for i in range(len(ave_cesmcam_waveform)):  
        max_accmip_waveform.append(np.max((ave_cesmcam_waveform[i],ave_cmam_waveform[i],ave_geosccm_waveform[i],ave_gisse2r_waveform[i],ave_mirocchem_waveform[i])))
        min_accmip_waveform.append(np.min((ave_cesmcam_waveform[i],ave_cmam_waveform[i],ave_geosccm_waveform[i],ave_gisse2r_waveform[i],ave_mirocchem_waveform[i])))
    min_accmip_waveform = np.array(min_accmip_waveform)
    max_accmip_waveform = np.array(max_accmip_waveform)
 
    ax.fill_between(model_time_pd.to_pydatetime(),min_accmip_waveform, max_accmip_waveform,where=max_accmip_waveform>min_accmip_waveform, facecolor='purple', interpolate=True,alpha=0.1)
    ax.plot_date(obs_time_pd.to_pydatetime(),ave_obs_waveform,color = 'black',linestyle='--',linewidth=1,markersize=5,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_accmip_waveform,color = 'purple',linestyle='--',linewidth=1,marker='o',markersize = 3.5,markeredgecolor='None')
    ax.plot_date(model_time_pd.to_pydatetime(),ave_geoschemv1001_waveform,color = 'red',linestyle='--',linewidth=1,marker='o',markersize = 3.5,markeredgecolor='None')
    
    ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=20)
    if plot_type == 'd':
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
    if plot_type == 's':
        ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
        
    ax.tick_params(axis='both', which='major', labelsize=18,pad = 7)

    tst = False
    for t in range(len(ave_obs_waveform)):
        b = np.isnan(ave_obs_waveform[t])
        if b == False:
            tst= True
    if tst == False:
        ax.plot_date(obs_time_pd.to_pydatetime(),len(obs_time_pd)*[1],markersize=0.000001)
    
    count+=1
    
plt.tight_layout(pad = 1.5)
fig.subplots_adjust(bottom=0.08)
fig.subplots_adjust(left=0.08)

plt.annotate('Time (Months)',(-1.65,-0.2),xycoords='axes fraction',fontsize=30, va='top')
plt.annotate('[ppb]',(-3.80,2.6),xycoords='axes fraction',fontsize=30, va='top',rotation=90)

h1, = ax.plot([1,1],color='black',marker='o',linestyle='None',markersize=10)
h2, = ax.plot([1,1],color='purple',marker='o',linestyle='None',markersize=10)
h3, = ax.plot([1,1],color='red',marker='o',linestyle='None',markersize=10)

plt.legend((h1,h2,h3),['Observations','ACCMIP Model Average','GEOSChem v1001'],loc='lower left',prop={'size':20},fancybox=True,ncol=2,markerscale=2.5,bbox_to_anchor=(-1.2,0.2))
h1.set_visible(False)
h2.set_visible(False)
h3.set_visible(False)

plt.show()
