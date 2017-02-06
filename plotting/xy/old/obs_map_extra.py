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
from bpch import bpch


def interactive(event):
    modules.clicker_interactive_map_obsmodel(event,species,lat_e,lon_e,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_daily_waveform,obs_seasonal_waveform,obs_full_waveform,model_daily_waveform,model_seasonal_waveform,model_full_waveform,fig,all_m)

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
    model_daily_waveform.append(site_group_mod.variables['daily_waveform'][:])
    model_seasonal_waveform.append(site_group_mod.variables['seasonal_waveform'][:]) 
    model_full_waveform.append(site_group_mod.variables['all_waveform'][:])
    
    obs_daily_mag = np.append(obs_daily_mag,site_group_obs.daily_amplitude)
    obs_daily_phase = np.append(obs_daily_phase,site_group_obs.daily_phase)
    obs_ave = np.append(obs_ave,site_group_obs.average)
    obs_seasonal_mag = np.append(obs_seasonal_mag,site_group_obs.seasonal_amplitude)
    obs_seasonal_phase = np.append(obs_seasonal_phase,site_group_obs.seasonal_phase)
    obs_daily_waveform.append(site_group_obs.variables['daily_waveform'][:])
    obs_seasonal_waveform.append(site_group_obs.variables['seasonal_waveform'][:]) 
    obs_full_waveform.append(site_group_obs.variables['all_waveform'][:])

tags = modules.get_tags(obs_refs)
loc_dict = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic Sites','NA':'North America','SA':'South America'}

#get model no data
root_grp = Dataset('/work/home/db876/plotting_tools/model_files/GEOSCHEM_SURFACE_2005_2012_v90103_2x2.5_GEOS5_H_*.nc')
no = root_grp.variables['no'][:]
print no.shape
no = np.ma.average(no,axis=0)
no_ave = []
for obs_lat,obs_lon in zip(obs_lats, obs_lons):
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    no_ave = np.append(no_ave,no[lat_n,lon_n])
#--------------------------------------------------------------------
#get model no2 data
no2 = root_grp.variables['no2'][:]
no2 = np.ma.average(no2,axis=0)
no2_ave = []
for obs_lat,obs_lon in zip(obs_lats, obs_lons):
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    no2_ave = np.append(no2_ave,no2[lat_n,lon_n])
#--------------------------------------------------------------------
#get model ro2 data
ro2 = root_grp.variables['ro2'][:]                                                                                                                                          
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
param_label = extra_param_labels[param_index]

diff = extra_param[param_index]

#set up plot

latlower_setup = [20,30,20,lat_e[0]]
latupper_setup = [80,72,55,lat_e[-1]]
lonwest_setup = [-170,-15,115,lon_e[0]]
loneast_setup = [-50,35,155,lon_e[-1]]
label_out = ['NA','EU','AS','ZZ']
label = ['NA','EU','AS','ROW']

diff = np.array(diff)
obs_lons = np.array(obs_lons)
obs_lats = np.array(obs_lats)
tags = np.array(tags)

#--------------------------------
fig, axes = plt.subplots(nrows=2, ncols=2,figsize=(23,12.3))
fig.patch.set_facecolor('white')
count = 0

all_m = []

for ax in axes.flat:

    #setup basemap projection
    m = Basemap(projection='cyl',llcrnrlat=latlower_setup[count],urcrnrlat=latupper_setup[count],llcrnrlon=lonwest_setup[count],urcrnrlon=loneast_setup[count],resolution='c',ax = ax)

    m.drawcoastlines()
    m.drawmapboundary()


    if count == 3:
        test = []
        other_tags = ['AF','ANT','ARC','OC','O','SA']
        for i in tags:
            if i in other_tags:
                test.append(True)
            else:
                test.append(False)
        test = np.array(test)
    else:
        current_tag = label_out[count] 
        test = tags == current_tag
    
    
    current_diff = diff[test]
    current_lons = obs_lons[test]
    current_lats = obs_lats[test]
    
    X,Y = m(current_lons,current_lats)

    if count == 0:
        m_size= 100
    if count == 1:
        m_size = 50
    if count == 2:
        m_size = 250
    if count == 3:
        m_size = 300

    for i in range(len(current_lons)):
        min_diff = np.min(abs(diff))
        max_diff = np.max(abs(diff))
        all = m.scatter(X[i],Y[i],c=current_diff[i], s=m_size, vmin = min_diff,vmax = max_diff, marker = 'o',edgecolor='black',linewidth=0.5,cmap=plt.cm.coolwarm,zorder=10)
    
    ax.text(0.03, 0.97, label[count], transform=ax.transAxes,fontsize=34, fontweight='bold', va='top')
    
    m.plot(obs_lons, obs_lats, color = 'black', marker = 'o', markersize = 0.001, linestyle= 'None', zorder=1, picker = 5)
    
    all_m.append(m)
    
    count+=1

plt.tight_layout(pad = 3.08)

fig.subplots_adjust(bottom=0.16)
cbar_ax = fig.add_axes([0.10, 0.08, 0.80, 0.06])
cb = fig.colorbar(all, cax=cbar_ax,orientation ='horizontal')
cb.set_label('%s'%(param_label), fontsize = 24)

cb.ax.tick_params(labelsize=22)

fig.canvas.mpl_connect('pick_event', interactive)

plt.show()
