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
from bpch import bpch



def interactive(event):
    modules.clicker_interactive_xy_obsmodel_multi(event,species,lat_e,lon_e,obs_lats,obs_lons,date,time,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_daily_waveform,obs_seasonal_waveform,obs_full_waveform,model_daily_waveform,model_seasonal_waveform,model_full_waveform,fig,ax,ax2,ax3,ax4,ax5,ax7,ax8)
#-----------------------------------------------------
#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres = modules.get_obs_info(present_dir)
model_fname,species,start_year,end_year = modules.get_model_info(present_dir)
    
obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

#get obs lat_lon grid central points
obs_lats_centre, obs_lons_centre, model_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_lats,obs_lons) 

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

#--------------------------------------------------------
#load in periodic lsp data
obs_period_grp = Dataset('../obs_%s_%s/obs_sig_periods.nc'%(vres,timeres))
model_period_grp = Dataset('model_sig_periods.nc')

model_daily_mag = []
model_daily_phase = []
obs_daily_mag = []
obs_daily_phase = []
model_ave = []
obs_ave = []
model_seasonal_mag = []
model_seasonal_phase = []
obs_seasonal_mag = []
obs_seasonal_phase = []
model_daily_ff = []
model_seasonal_ff = []
obs_daily_ff = []
obs_seasonal_ff = []
obs_daily_waveform = []
model_daily_waveform = []
obs_seasonal_waveform = []                                                                                                                                                                                                                       
model_seasonal_waveform = []
obs_full_waveform = []                                                                                                                                                                                                                       
model_full_waveform = []

for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    site_group_mod = model_period_grp.groups[ref]
    
    model_daily_mag = np.append(model_daily_mag,site_group_mod.daily_amplitude)
    model_daily_phase = np.append(model_daily_phase,site_group_mod.daily_phase)
    model_ave = np.append(model_ave,site_group_mod.average)
    model_seasonal_mag = np.append(model_seasonal_mag,site_group_mod.seasonal_amplitude)
    model_seasonal_phase = np.append(model_seasonal_phase,site_group_mod.seasonal_phase)
    model_daily_ff = np.append(model_daily_ff,site_group_mod.daily_ff) 
    model_seasonal_ff = np.append(model_seasonal_ff,site_group_mod.seasonal_ff)
    model_daily_waveform.append(site_group_mod.variables['daily_waveform'][:])
    model_seasonal_waveform.append(site_group_mod.variables['seasonal_waveform'][:]) 
    model_full_waveform.append(site_group_mod.variables['all_waveform'][:]) 
    
    obs_daily_mag = np.append(obs_daily_mag,site_group_obs.daily_amplitude)
    obs_daily_phase = np.append(obs_daily_phase,site_group_obs.daily_phase)
    obs_ave = np.append(obs_ave,site_group_obs.average)
    obs_seasonal_mag = np.append(obs_seasonal_mag,site_group_obs.seasonal_amplitude)
    obs_seasonal_phase = np.append(obs_seasonal_phase,site_group_obs.seasonal_phase)
    obs_daily_ff = np.append(obs_daily_ff,site_group_obs.daily_ff) 
    obs_seasonal_ff = np.append(obs_seasonal_ff,site_group_obs.seasonal_ff)
    obs_daily_waveform.append(site_group_obs.variables['daily_waveform'][:])
    obs_seasonal_waveform.append(site_group_obs.variables['seasonal_waveform'][:]) 
    obs_full_waveform.append(site_group_obs.variables['all_waveform'][:])

month_lengths = [0,31,59.25,90.25,120.25,151.25,181.25,212.25,243.25,273.25,304.25,334.25,365.25]
#month_lengths = [273.25,304.25,334.25,365.25,396.25,424.5,455.5,485.5,516.5,546.5,577.5,608.5]
year_fact = 12./365.25
month_array  = [i*year_fact for i in month_lengths]
month_strings = ["JAN","","MAR","","MAY","","JUL","","SEP","","NOV",""]
axis_nums = month_array
		
#--------------------------------------------------------------------------
#rearrange data so that NA plotted 1, EU plotted 2, AS plotted 3, O plotted 4, OC plotted 5, AF plotted 6, SA plotted 7, ANT plotted 8, ARC plotted 9
loc_dict = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic Sites','NA':'North America','SA':'South America'}
order_dict = {'ANT':8,'ARC':9,'AF':6,'AS':3,'EU':2,'OC':5,'O':4,'NA':1,'SA':7}

order_tags = []
for tag in tags:
    n = order_dict[tag]
    order_tags.append(n)

tags_argsort = np.argsort(order_tags)

tags = tags[tags_argsort]
obs_refs = obs_refs[tags_argsort]
obs_lats = obs_lats[tags_argsort]
obs_lons = obs_lons[tags_argsort]
obs_alt = obs_alt[tags_argsort]
obs_complete = obs_complete[tags_argsort]
obs_daily_mag = obs_daily_mag[tags_argsort]
model_daily_mag = model_daily_mag[tags_argsort]
obs_daily_phase = obs_daily_phase[tags_argsort]
model_daily_phase = model_daily_phase[tags_argsort]
obs_ave = obs_ave[tags_argsort]
model_ave = model_ave[tags_argsort]
model_seasonal_mag = model_seasonal_mag[tags_argsort]
model_seasonal_phase = model_seasonal_phase[tags_argsort]
obs_seasonal_mag = obs_seasonal_mag[tags_argsort]
obs_seasonal_phase = obs_seasonal_phase[tags_argsort]
obs_daily_ff = obs_daily_ff[tags_argsort]
obs_seasonal_ff = obs_seasonal_ff[tags_argsort]
model_daily_ff = model_daily_ff[tags_argsort]
model_seasonal_ff = model_seasonal_ff[tags_argsort]
obs_daily_ff = obs_daily_ff[tags_argsort]
obs_seasonal_ff = obs_seasonal_ff[tags_argsort]
model_daily_ff = model_daily_ff[tags_argsort]
model_seasonal_ff = model_seasonal_ff[tags_argsort]
obs_daily_waveform = np.array(obs_daily_waveform)[tags_argsort]
obs_seasonal_waveform = np.array(obs_seasonal_waveform)[tags_argsort]
obs_full_waveform = np.array(obs_full_waveform)[tags_argsort]
model_daily_waveform = np.array(model_daily_waveform)[tags_argsort]
model_seasonal_waveform = np.array(model_seasonal_waveform)[tags_argsort]
model_full_waveform = np.array(model_full_waveform)[tags_argsort]

#---------------------------------------------------------------------------

#setup shape dict for plotting
loc_colors = {'ANT':'magenta','ARC':'purple','AF':'orange','AS':'blue','EU':'red','OC':'black','O':'yellow','NA':'green','SA':'pink'}
loc_shapes = {'ANT':'<','ARC':'+','AF':'^','AS':'d','EU':'s','OC':'>','O':'v','NA':'o','SA':'*'}

#--------------------------------------------------------------------
#get model no data
root_grp = Dataset('/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_NO_2005_2010_v90103_2x2.5_GEOS5_H_*.nc')
no = root_grp.variables['no'][:]
no = no*1e9
no = np.ma.average(no,axis=0)
no_ave = []
for obs_lat,obs_lon in zip(obs_lats, obs_lons):
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    no_ave = np.append(no_ave,no[lat_n,lon_n])
#--------------------------------------------------------------------
#get model no2 data
root_grp = Dataset('/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_NO2_2005_2010_v90103_2x2.5_GEOS5_H_*.nc')
no2 = root_grp.variables['no2'][:]
no2 = no2*1e9
no2 = np.ma.average(no2,axis=0)
no2_ave = []
for obs_lat,obs_lon in zip(obs_lats, obs_lons):
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    no2_ave = np.append(no2_ave,no2[lat_n,lon_n])
#--------------------------------------------------------------------
#get model ro2 data
root_grp = Dataset('/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_RO2_2005_2010_v90103_2x2.5_GEOS5_H_*.nc')
ro2 = root_grp.variables['ro2'][:]                                                                                                                                          
ro2 = ro2*1e9
ro2 = np.ma.average(ro2,axis=0)
ro2_ave = []
for obs_lat,obs_lon in zip(obs_lats, obs_lons):
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    ro2_ave = np.append(ro2_ave,ro2[lat_n,lon_n])

#-----------------------------------------------------------------------------
#get all average data
files = glob.glob('/work/home/db876/modelling/GEOS_CHEM/v90103/2x2.5/2x2.5_GEOS5_fullgrid/run/*.ctm.bpch')
files = sorted(files)

o3_prod = np.empty((len(files),len(lat_c),len(lon_c)))
o3_loss = np.empty((len(files),len(lat_c),len(lon_c)))
o3 = np.empty((len(files),len(lat_c),len(lon_c)))
co = np.empty((len(files),len(lat_c),len(lon_c)))
nox = np.empty((len(files),len(lat_c),len(lon_c)))
acet = np.empty((len(files),len(lat_c),len(lon_c)))
isop = np.empty((len(files),len(lat_c),len(lon_c)))
c2h6 = np.empty((len(files),len(lat_c),len(lon_c)))
c3h8 = np.empty((len(files),len(lat_c),len(lon_c))) 
ald2 = np.empty((len(files),len(lat_c),len(lon_c)))
alk4 = np.empty((len(files),len(lat_c),len(lon_c)))
ch2o = np.empty((len(files),len(lat_c),len(lon_c)))
mek = np.empty((len(files),len(lat_c),len(lon_c)))
prpe = np.empty((len(files),len(lat_c),len(lon_c)))

for i in range(len(files)):
    data = bpch(files[i],tracerinfo='/work/home/db876/modelling/GEOS_CHEM/v90103/2x2.5/2x2.5_GEOS5_fullgrid/run/valid_tracerinfo.dat')
    group_pl = data.groups['PORL-L=$']
    po3 = group_pl.variables['PO3']
    lo3 = group_pl.variables['LO3']
    po3 = po3[:,0,:,:]
    lo3 = lo3[:,0,:,:]
    o3_prod[i,:,:] = po3
    o3_loss[i,:,:] = lo3

#need to convert voc species from ppbC to ppbV by dividing by number of Carbon atoms in molecule    
    group_ave = data.groups['IJ-AVG-$']
    o3[i,:,:] = group_ave.variables['O3'][:,0,:,:]
    co[i,:,:] = group_ave.variables['CO'][:,0,:,:]
    nox[i,:,:] = group_ave.variables['NOx'][:,0,:,:] 
    acet[i,:,:] = group_ave.variables['ACET'][:,0,:,:]/3. 
    isop[i,:,:] = group_ave.variables['ISOP'][:,0,:,:]/5. 
    c2h6[i,:,:] = group_ave.variables['C2H6'][:,0,:,:]/2. 
    c3h8[i,:,:] = group_ave.variables['C3H8'][:,0,:,:]/3. 
    ald2[i,:,:] = group_ave.variables['ALD2'][:,0,:,:]/2. 
    alk4[i,:,:] = group_ave.variables['ALK4'][:,0,:,:]/4.
    ch2o[i,:,:] = group_ave.variables['CH2O'][:,0,:,:]
    mek[i,:,:] = group_ave.variables['MEK'][:,0,:,:]/4.
    prpe[i,:,:] = group_ave.variables['PRPE'][:,0,:,:]/3.


print 'nox units = ', group_ave.variables['NOx'].units
print 'acet units = ',group_ave.variables['ACET'].units
print 'isop units = ',group_ave.variables['ISOP'].units
print 'c2h6 units = ',group_ave.variables['C2H6'].units
print 'c3h8 units = ',group_ave.variables['C3H8'].units
print 'ald2 units = ',group_ave.variables['ALD2'].units
print 'alk4 units = ',group_ave.variables['ALK4'].units
print 'ch2o units = ',group_ave.variables['CH2O'].units
print 'mek units = ',group_ave.variables['MEK'].units
print 'prpe units = ',group_ave.variables['PRPE'].units

o3_pl = o3_prod-o3_loss
o3_pl = np.ma.average(o3_pl,axis=0)

o3 = np.ma.average(o3,axis=0)
co = np.ma.average(co,axis=0)
nox = np.ma.average(nox,axis=0)
acet = np.ma.average(acet,axis=0)
isop = np.ma.average(isop,axis=0)
c2h6 = np.ma.average(c2h6,axis=0)
c3h8 = np.ma.average(c3h8,axis=0)
ald2 = np.ma.average(ald2,axis=0)
alk4 = np.ma.average(alk4,axis=0)
ch2o = np.ma.average(ch2o,axis=0)
mek = np.ma.average(mek,axis=0)
prpe = np.ma.average(prpe,axis=0)

o3_pl_data = []
o3_ave = []
co_ave = []
nox_ave = []
acet_ave = []
isop_ave = []
c2h6_ave = []
c3h8_ave = []
ald2_ave = []
alk4_ave = []
ch2o_ave = []
mek_ave = []
prpe_ave = []

for obs_lat,obs_lon in zip(obs_lats, obs_lons):
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    o3_pl_data = np.append(o3_pl_data,o3_pl[lat_n,lon_n])

    o3_ave = np.append(o3_ave,o3[lat_n,lon_n])    
    co_ave = np.append(co_ave,co[lat_n,lon_n])
    nox_ave = np.append(nox_ave,nox[lat_n,lon_n])
    acet_ave = np.append(acet_ave,acet[lat_n,lon_n])
    isop_ave = np.append(isop_ave,isop[lat_n,lon_n])
    c2h6_ave = np.append(c2h6_ave,c2h6[lat_n,lon_n])
    c3h8_ave = np.append(c3h8_ave,c3h8[lat_n,lon_n])
    ald2_ave = np.append(ald2_ave,ald2[lat_n,lon_n])
    alk4_ave = np.append(alk4_ave,alk4[lat_n,lon_n])
    ch2o_ave = np.append(ch2o_ave,ch2o[lat_n,lon_n])
    mek_ave = np.append(mek_ave,mek[lat_n,lon_n])
    prpe_ave = np.append(prpe_ave,prpe[lat_n,lon_n])

voc_ave = acet_ave+ald2_ave+alk4_ave+c2h6_ave+c3h8_ave+ch2o_ave+isop_ave+mek_ave+prpe_ave
no2_no_ratio = no2_ave/no_ave
voc_nox_ratio = voc_ave/nox_ave

#---------------------------------------------------------------------
#get extra parameter to use as colour

extra_param_name = ['obs_complete','alt','lat','lon','model_no','model_no2','model_nox','model_no2_no_ratio','model_co','model_ro2','model_isop','model_acet','model_c2h6','model_c3h8','o3_pl','model_voc_nox_ratio']
extra_param_labels = ['Data Valid Percent','Alitude (m)','Absolute Latitude','Longitude','GEOS Chem v90103 2x2.5 NO (ppb)','GEOS Chem v90103 2x2.5 NO2 (ppb)','GEOS Chem v90103 2x2.5 NOx (ppb)','GEOS Chem v90103 2x2.5 NO2/NO ratio (ppb)','GEOS Chem v90103 2x2.5 CO (ppb)','GEOS Chem v90103 2x2.5 RO2 (ppb)','GEOS Chem v90103 2x2.5 ISOP (ppb)','GEOS Chem v90103 2x2.5 ACET (ppb)','GEOS Chem v90103 2x2.5 C2H6 (ppb)','GEOS Chem v90103 2x2.5 C3H8 (ppb)','GEOS Chem v90103 2x2.5 P/L','GEOS Chem v90103 2x2.5 VOC/NOx ratio (ppb)']
extra_param = [obs_complete,obs_alt,np.abs(obs_lats),obs_lons,no_ave,no2_ave,nox_ave,no2_no_ratio,co_ave,ro2_ave,isop_ave,acet_ave,c2h6_ave,c3h8_ave,o3_pl_data,voc_nox_ratio]

param = raw_input('\nChoose Extra Paramter to Plot.\n%s\n'%('   '.join(str(i) for i in extra_param_name)))
param_index = extra_param_name.index(param)

#--------------------------------------------------------------------
#set up plot
fig, ((ax, ax2,ax3), (ax4, ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3, 3,figsize=(14,12))
fig.patch.set_facecolor('white')
fig.subplots_adjust(left=0.18,bottom=0.12)

marker_size = 25
font_size = 10
pad_size = 5

#get colors and labels for site points
plot_colors = []

for t in range(len(obs_refs)):
    plot_colors.append(loc_colors[tags[t]])
    
#----------------------------------------------------------------------
#Daily Amplitude Plot

#ax = fig.add_subplot(6,1,1)

current_min = np.min(extra_param[param_index])
current_max = np.max(extra_param[param_index])
#current_min = 0
#current_max = 10

for t in range(len(obs_refs)):
    #line, = ax.plot(obs_daily_mag[t], model_daily_mag[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
    #extra = ax.scatter(obs_daily_mag[t],model_daily_mag[t], c = extra_param[param_index][t], marker = loc_shapes[tags[t]], s = marker_size, label = loc_dict[tags[t]],cmap=plt.cm.coolwarm,zorder=10,vmin=np.min(extra_param[param_index]),vmax =np.max(extra_param[param_index]))
    extra = ax.scatter(obs_daily_mag[t],model_daily_mag[t], c = extra_param[param_index][t], marker = loc_shapes[tags[t]], s = marker_size, label = loc_dict[tags[t]],cmap=plt.cm.coolwarm,zorder=10,vmin=current_min,vmax = current_max)

ax.plot(obs_daily_mag, model_daily_mag, color = 'black', marker = 'o', markersize = 1, linestyle= 'None',zorder=1, picker = 5)

#ax.set_xscale('log')
#ax.set_yscale('log')
obs_max = np.max(obs_daily_mag)+2
model_max = np.max(model_daily_mag)+2
if np.min(obs_daily_mag) < np.min(model_daily_mag):
    min_val = np.min(obs_daily_mag)-0.001
else:
    min_val = np.min(model_daily_mag)-0.001
if obs_max > model_max:
    ax.set_xlim(min_val,obs_max)
    ax.set_ylim(min_val,obs_max)
else:
    ax.set_xlim(min_val,model_max)
    ax.set_ylim(min_val,model_max)
x = np.arange(0,1000,0.005)
y = np.arange(0,1000,0.005)

handles, labels = ax.get_legend_handles_labels()
hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
handles, labels = zip(*hl)
by_label = OrderedDict(zip(labels, handles))
ax.plot(x,y,'k--',alpha = 0.5)
ax.set_ylabel('Model (ppb)',fontsize = font_size)
ax.set_xlabel('Observations (ppb)',fontsize = font_size)
			
#----------------------------------------------------------------------
#Seasonal Amplitude Plot

for t in range(len(obs_refs)):
    #ax3.plot(obs_annual_mag[t], model_annual_mag[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
    ax2.scatter(obs_seasonal_mag[t],model_seasonal_mag[t], c=extra_param[param_index][t], marker = loc_shapes[tags[t]], s = marker_size, label = loc_dict[tags[t]],cmap=plt.cm.coolwarm,zorder=10,vmin=current_min,vmax =current_max)

ax2.plot(obs_seasonal_mag, model_seasonal_mag, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

obs_max = np.max(obs_seasonal_mag)+1
model_max = np.max(model_seasonal_mag)+1
if np.min(obs_seasonal_mag) < np.min(model_seasonal_mag):
    min_val = np.min(obs_seasonal_mag)-1
else:
    min_val = np.min(model_seasonal_mag)-1
if obs_max > model_max:
    ax2.set_xlim(min_val,obs_max)
    ax2.set_ylim(min_val,obs_max)
else:
    ax2.set_xlim(min_val,model_max)
    ax2.set_ylim(min_val,model_max)
x = np.arange(0,1000,1)
y = np.arange(0,1000,1)

ax2.plot(x,y,'k--',alpha = 0.5)
ax2.set_ylabel('Model (ppb)',fontsize = font_size)
ax2.set_xlabel('Observations (ppb)',fontsize = font_size)

#----------------------------------------------------------------------
#Daily Phase Plot

for t in range(len(obs_refs)):
    #ax4.plot(obs_daily_phase[t], model_daily_phase[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
    ax4.scatter(obs_daily_phase[t],model_daily_phase[t], c=extra_param[param_index][t], marker = loc_shapes[tags[t]], s = marker_size, label = loc_dict[tags[t]],cmap=plt.cm.coolwarm,zorder=10,vmin=current_min,vmax =current_max)

ax4.plot(obs_daily_phase, model_daily_phase, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

ax4.set_xlim(0,24)
ax4.set_ylim(0,24)
x = np.arange(0,50,1)
y = np.arange(0,50,1)

ax4.set_xticks(range(0,25,2))
ax4.set_yticks(range(0,25,2))

ax4.set_ylabel('Model (ppb)',fontsize = font_size)
ax4.set_xlabel('Observations (Hours)',fontsize = font_size)
ax4.plot(x,y,'k--',alpha = 0.5)
for tick in ax4.get_xaxis().get_major_ticks():
    tick.set_pad(pad_size)
for tick in ax4.get_yaxis().get_major_ticks():
    tick.set_pad(pad_size)
	
#----------------------------------------------------------------------
#Seasonal Phase Plot

#change phase to run from sep to sep, instead of jan to jan
#obs_seasonal_phase_corr, model_seasonal_phase_corr = modules.annual_phase_shift(obs_seasonal_phase,model_seasonal_phase)

for t in range(len(obs_refs)):
    #ax5.plot(obs_annual_phase_corr[t], model_annual_phase_corr[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
     ax5.scatter(obs_seasonal_phase[t],model_seasonal_phase[t], c=extra_param[param_index][t], marker = loc_shapes[tags[t]], s = marker_size, label = loc_dict[tags[t]],cmap=plt.cm.coolwarm,zorder=10,vmin=current_min,vmax =current_max)

ax5.plot(obs_seasonal_phase, model_seasonal_phase, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

ax5.set_xlim(0,11)
ax5.set_ylim(0,11)
#ax5.set_xlim(9,21)
#ax5.set_ylim(9,21)
x = np.arange(0,50,1)
y = np.arange(0,50,1)

ax5.set_xticks(month_array)
ax5.set_yticks(month_array)
ax5.set_xticklabels(month_strings,fontsize=9.2)
ax5.set_yticklabels(month_strings,fontsize=9.2)
ax5.set_ylabel('Model (ppb)',fontsize = font_size)
ax5.set_xlabel('Observations (Months)',fontsize = font_size)
ax5.plot(x,y,'k--',alpha = 0.5)
for tick in ax5.get_xaxis().get_major_ticks():
    tick.set_pad(pad_size)
for tick in ax5.get_yaxis().get_major_ticks():
    tick.set_pad(pad_size)	

#----------------------------------------------------------------------
#Average Plot   
    
for t in range(len(obs_refs)):
    #ax7.plot(obs_ave[t], model_ave[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
    ax3.scatter(obs_ave[t],model_ave[t], c=extra_param[param_index][t], marker = loc_shapes[tags[t]], s = marker_size, label = loc_dict[tags[t]],cmap=plt.cm.coolwarm,zorder=10,vmin=current_min,vmax = current_max)

ax3.plot(obs_ave, model_ave, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

obs_max = np.max(obs_ave)+1
model_max = np.max(model_ave)+1
if np.min(obs_ave) < np.min(model_ave):
    min_val = np.min(obs_ave)-1
else:
    min_val = np.min(model_ave)-1
if obs_max > model_max:
    ax3.set_xlim(min_val,obs_max)
    ax3.set_ylim(min_val,obs_max)
else:
    ax3.set_xlim(min_val,model_max)
    ax3.set_ylim(min_val,model_max)
x = np.arange(0,1000,1)
y = np.arange(0,1000,1)

ax3.plot(x,y,'k--',alpha = 0.5)
ax3.set_ylabel('Model (ppb)',fontsize = font_size)
ax3.set_xlabel('Observations (ppb)',fontsize = font_size)

#----------------------------------------------------------------------
#Daily Form Factor Plot   

for t in range(len(obs_refs)):
    #ax7.plot(obs_daily_ff[t], model_daily_ff[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
    ax7.scatter(obs_daily_ff[t], model_daily_ff[t], c=extra_param[param_index][t], marker = loc_shapes[tags[t]], s = marker_size, label = loc_dict[tags[t]],cmap=plt.cm.coolwarm,zorder=10,vmin=current_min,vmax = current_max)

ax7.plot(obs_daily_ff, model_daily_ff, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

obs_max = np.max(obs_daily_ff)+0.1
model_max = np.max(model_daily_ff)+0.1
if np.min(obs_daily_ff) < np.min(model_daily_ff):
    min_val = np.min(obs_daily_ff)-0.1
else:
    min_val = np.min(model_daily_ff)-0.1
if obs_max > model_max:
    ax7.set_xlim(min_val,obs_max)
    ax7.set_ylim(min_val,obs_max)
else:
    ax7.set_xlim(min_val,model_max)
    ax7.set_ylim(min_val,model_max)
x = np.arange(0,1000,1)
y = np.arange(0,1000,1)

ax7.plot(x,y,'k--')

ax7.set_ylabel('Model Form Factor',fontsize = font_size)
ax7.set_xlabel('Observational Form Factor',fontsize = font_size)
        
#----------------------------------------------------------------------
#Seasonal Form Factor Plot   

for t in range(len(obs_refs)):
    #ax8.plot(obs_seasonal_ff[t], model_seasonal_ff[t], color = loc_colors[tags[t]], marker = 'o', markersize = marker_size, label = loc_dict[tags[t]], zorder=10)
    ax8.scatter(obs_seasonal_ff[t], model_seasonal_ff[t], c=extra_param[param_index][t], marker = loc_shapes[tags[t]], s = marker_size, label = loc_dict[tags[t]],cmap=plt.cm.coolwarm,zorder=10,vmin=current_min,vmax = current_max)

ax8.plot(obs_seasonal_ff, model_seasonal_ff, color = 'black', marker = 'o', markersize = 1, linestyle= 'None', zorder=1, picker = 5)

obs_max = np.max(obs_seasonal_ff)+0.1
model_max = np.max(model_seasonal_ff)+0.1
if np.min(obs_seasonal_ff) < np.min(model_seasonal_ff):
    min_val = np.min(obs_seasonal_ff)-0.1
else:
    min_val = np.min(model_seasonal_ff)-0.1
if obs_max > model_max:
    ax8.set_xlim(min_val,obs_max)
    ax8.set_ylim(min_val,obs_max)
else:
    ax8.set_xlim(min_val,model_max)
    ax8.set_ylim(min_val,model_max)
x = np.arange(0,1000,1)
y = np.arange(0,1000,1)

ax8.plot(x,y,'k--')

ax8.set_ylabel('Model Form Factor',fontsize = font_size)
ax8.set_xlabel('Observational Form Factor',fontsize = font_size)

#------------------------------------------------
#plot big labels
plt.figtext(0.24, 0.92, 'Diurnal', fontsize=22)
plt.figtext(0.50, 0.92, 'Seasonal', fontsize=22)
plt.figtext(0.76, 0.92, 'Average', fontsize=22)

plt.figtext(0.01, 0.77, 'Amplitude', fontsize=22)
plt.figtext(0.01, 0.50, 'Phase', fontsize=22)
plt.figtext(0.01, 0.23, 'Form Factor', fontsize=22)
#--------------------------------------------------

#make axis labels tight to plots
ax.yaxis.labelpad = 0 
ax2.yaxis.labelpad = 0
ax3.yaxis.labelpad = 0
ax4.yaxis.labelpad = 0
ax5.yaxis.labelpad = 0
ax7.yaxis.labelpad = 0
ax8.yaxis.labelpad = 0

ax.xaxis.labelpad = 0 
ax2.xaxis.labelpad = 0
ax3.yaxis.labelpad = 0
ax4.xaxis.labelpad = 0
ax5.xaxis.labelpad = 0
ax7.xaxis.labelpad = 0
ax8.yaxis.labelpad = 0

#plot legend
handles, labels = ax5.get_legend_handles_labels()
hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
handles, labels = zip(*hl)
by_label = OrderedDict(zip(labels, handles))
leg = ax.legend(by_label.values(), by_label.keys(), loc = 'upper center', bbox_to_anchor=(1.7,-2.55),fancybox=True,ncol=4)

cbar = fig.colorbar(extra,orientation='vertical')
cbar.set_label(extra_param_labels[param_index],rotation = 90)

#handles, labels = ax5.get_legend_handles_labels()
#lgd = ax5.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.45,-1.35),fancybox=True,ncol=5)

#--------------------------------------------------------------------------------
#remove ax 8 and ax9 for table

#ax3.set_frame_on(False)
ax6.set_frame_on(False)
#ax8.set_frame_on(False)
ax9.set_frame_on(False)
#ax3.axes.get_yaxis().set_visible(False)
ax6.axes.get_yaxis().set_visible(False)
#ax8.axes.get_yaxis().set_visible(False)
ax9.axes.get_yaxis().set_visible(False)
#ax3.axes.get_xaxis().set_visible(False)
ax6.axes.get_xaxis().set_visible(False)
#ax8.axes.get_xaxis().set_visible(False)
ax9.axes.get_xaxis().set_visible(False)

#plot separation line
#plt.plot([0, 1], [0.355, 0.355],linestyle='--', color='red',lw=2,transform=gcf().transFigure, clip_on=False) 

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500+100")
# 
fig.canvas.mpl_connect('pick_event', interactive)


plt.show()


