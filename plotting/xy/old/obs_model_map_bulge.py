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

def interactive(event):
    modules.clicker_interactive_map_obsmodel(event,species,lat_e,lon_e,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_d_waveform,obs_s_waveform,obs_all_waveform,model_d_waveform,model_s_waveform,model_all_waveform,fig,all_m)

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
tags = modules.get_tags(np.copy(obs_refs))

#--------------------------------------------------------
#load in periodic lsp data
obs_period_grp = Dataset('../obs_%s_%s/obs_sig_periods.nc'%(vres,timeres))
model_period_grp = Dataset('model_sig_periods.nc')

obs_d_waveform = []
obs_s_waveform = []
obs_all_waveform = []
model_d_waveform = []
model_s_waveform = []
model_all_waveform = []

#print obs_refs
for ref in obs_refs:
    #print ref
    site_group_obs = obs_period_grp.groups[ref]
    site_group_mod = model_period_grp.groups[ref]

    model_d_waveform.append(site_group_mod.variables['daily_waveform'][:])
    model_s_waveform.append(site_group_mod.variables['seasonal_waveform'][:])
    model_all_waveform.append(site_group_mod.variables['all_waveform'][:])
    
    obs_d_waveform.append(site_group_obs.variables['daily_waveform'][:])
    obs_s_waveform.append(site_group_obs.variables['seasonal_waveform'][:])
    obs_all_waveform.append(site_group_obs.variables['all_waveform'][:])


rms_model_s = []
rms_obs_s = []
mean_model_s = []
mean_obs_s = []
peak_model_s = []
peak_obs_s = []
rms_model_d = []
rms_obs_d = []
mean_model_d = []
mean_obs_d = []
peak_model_d = []
peak_obs_d = []

#obs_s_waveform = obs_s_waveform-np.average(obs_s_waveform)
#model_s_waveform = model_s_waveform-np.average(model_s_waveform)
#obs_d_waveform = obs_d_waveform-np.average(obs_d_waveform)
#model_d_waveform = model_d_waveform-np.average(model_d_waveform)

for i in range(len(model_s_waveform)):
    rms_model_s.append(np.sqrt(np.mean(np.square(model_s_waveform[i]))))
    rms_obs_s.append(np.sqrt(np.mean(np.square(obs_s_waveform[i]))))
    mean_model_s.append(np.average(np.abs(model_s_waveform[i])))
    mean_obs_s.append(np.average(np.abs(obs_s_waveform[i])))
    peak_model_s.append(np.square(np.max(model_s_waveform[i])))
    peak_obs_s.append(np.square(np.max(obs_s_waveform[i])))
    
    rms_model_d.append(np.sqrt(np.mean(np.square(model_d_waveform[i]))))
    rms_obs_d.append(np.sqrt(np.mean(np.square(obs_d_waveform[i]))))
    mean_model_d.append(np.average(np.abs(model_d_waveform[i])))
    mean_obs_d.append(np.average(np.abs(obs_d_waveform[i])))
    peak_model_d.append(np.square(np.max(model_d_waveform[i])))
    peak_obs_d.append(np.square(np.max(obs_d_waveform[i])))

rms_model_s = np.array(rms_model_s)
rms_obs_s = np.array(rms_obs_s)
mean_model_s = np.array(mean_model_s)
mean_obs_s = np.array(mean_obs_s)
peak_model_s = np.array(peak_model_s)
peak_obs_s = np.array(peak_obs_s)
rms_model_d = np.array(rms_model_d)
rms_obs_d = np.array(rms_obs_d)
mean_model_d = np.array(mean_model_d)
mean_obs_d = np.array(mean_obs_d)
peak_model_d = np.array(peak_model_d)
peak_obs_d = np.array(peak_obs_d)

model_ff_s = rms_model_s/mean_model_s
obs_ff_s = rms_obs_s/mean_obs_s
model_cf_s = peak_model_s/rms_model_s
obs_cf_s =  peak_obs_s/rms_obs_s

model_ff_d = rms_model_d/mean_model_d
obs_ff_d = rms_obs_d/mean_obs_d
model_cf_d = peak_model_d/rms_model_d
obs_cf_d =  peak_obs_d/rms_obs_d

sf_diff_s =  model_ff_s - obs_ff_s
cf_diff_s = model_cf_s - obs_cf_s
sf_diff_d =  model_ff_d - obs_ff_d
cf_diff_d = model_cf_d - obs_cf_d

param = raw_input('\nall or bulge?\n')
if param == 'bulge':
    bparam = raw_input('\nall or sign?\n')
    if bparam == 'all':
        bn = raw_input('\nBulge Number?\n')
    if bparam == 'sign':
        bsign = raw_input('\n+ or -?\n')
    
    pt = raw_input('\nChoose plot type. sign, maxoffset, maxoffsetph, aveoffset, totaloffset, startphoffset, endphoffset, totalphoffset, area or time\n')

if param == 'all':
    pt = raw_input('\nChoose plot type. nbulges, n+bulges, n-bulges, aveoffset, ave+offset, total+offset, ave-offset, total-offset, area+, area-, time+ or time-\n')
             
max_ph = 12.
    

#load in bulge data
bulge_data = Dataset('bulge.nc')

nbulges = []
nplusbulges = []
nnegbulges = []
a_aveconcoffset = [] 
a_aveplusconcoffset = []
a_totalplusconcoffset = []
a_avenegconcoffset = []
a_totalnegconcoffset = []
a_pcplusbulgearea = []
a_pcnegbulgearea = []
a_pcplusbulgetime = []
a_pcnegbulgetime = []
b_sign = []
b_maxconcoffset = []
b_maxconcoffsetph = []
b_aveconcoffset = []
b_totalconcoffset = []
b_startphaseoffset = []
b_endphaseoffset = []
b_totalphaseoffset = []
b_pcbulgearea = []
b_pcbulgetime = []


for ref in obs_refs:
    site_grp = bulge_data.groups[ref]
    
    nbulges.append(site_grp.nbulges)
    nplusbulges.append(site_grp.nplusbulges)
    nnegbulges.append(site_grp.nnegbulges)
    a_aveconcoffset.append(site_grp.a_aveconcoffset)
    a_aveplusconcoffset.append(site_grp.a_aveplusconcoffset)
    a_totalplusconcoffset.append(site_grp.a_totalplusconcoffset)
    a_avenegconcoffset.append(site_grp.a_avenegconcoffset)
    a_totalnegconcoffset.append(site_grp.a_totalnegconcoffset)
    a_pcplusbulgearea.append(site_grp.a_pcplusbulgearea)
    a_pcnegbulgearea.append(site_grp.a_pcnegbulgearea)
    a_pcplusbulgetime.append(site_grp.a_pcplusbulgetime)
    a_pcnegbulgetime.append(site_grp.a_pcnegbulgetime)
    
    #bulge_grp = site_grp.groups['bulge1']
    
    if param == 'bulge':
        if bparam == 'all':
            bulge_grp = site_grp.groups['bulge%s'%(bn)]
            
            b_sign.append(bulge_grp.b_sign)
            b_maxconcoffset.append(bulge_grp.b_maxconcoffset)
            b_maxconcoffsetph.append(bulge_grp.b_maxconcoffsetph)
            b_aveconcoffset.append(bulge_grp.b_aveconcoffset)
            b_totalconcoffset.append(bulge_grp.b_totalconcoffset)
            b_startphaseoffset.append(bulge_grp.b_startphaseoffset)
            b_endphaseoffset.append(bulge_grp.b_endphaseoffset)
            b_totalphaseoffset.append(bulge_grp.b_totalphaseoffset)
            b_pcbulgearea.append(bulge_grp.b_pcbulgearea)
            b_pcbulgetime.append(bulge_grp.b_pcbulgetime)
        
        
        if bparam == 'sign':
            found_sign = False
            for n in range(nbulges[-1]):
                bulge_grp = site_grp.groups['bulge%s'%(n+1)]
                sign = bulge_grp.b_sign
            
                if bsign == '+':
                    if sign  == '-':
                        continue
                    else:
                        b_sign.append(bulge_grp.b_sign)
                        b_maxconcoffset.append(bulge_grp.b_maxconcoffset)
                        b_maxconcoffsetph.append(bulge_grp.b_maxconcoffsetph)
                        b_aveconcoffset.append(bulge_grp.b_aveconcoffset)
                        b_totalconcoffset.append(bulge_grp.b_totalconcoffset)
                        b_startphaseoffset.append(bulge_grp.b_startphaseoffset)
                        b_endphaseoffset.append(bulge_grp.b_endphaseoffset)
                        b_totalphaseoffset.append(bulge_grp.b_totalphaseoffset)
                        b_pcbulgearea.append(bulge_grp.b_pcbulgearea)
                        b_pcbulgetime.append(bulge_grp.b_pcbulgetime)
                        found_sign = True
                        break
                
                if bsign == '-':
                    if sign  == '+':
                        continue
                    else:
                        b_sign.append(bulge_grp.b_sign)
                        b_maxconcoffset.append(bulge_grp.b_maxconcoffset)
                        b_maxconcoffsetph.append(bulge_grp.b_maxconcoffsetph)
                        b_aveconcoffset.append(bulge_grp.b_aveconcoffset)
                        b_totalconcoffset.append(bulge_grp.b_totalconcoffset)
                        b_startphaseoffset.append(bulge_grp.b_startphaseoffset)
                        b_endphaseoffset.append(bulge_grp.b_endphaseoffset)
                        b_totalphaseoffset.append(bulge_grp.b_totalphaseoffset)
                        b_pcbulgearea.append(bulge_grp.b_pcbulgearea)
                        b_pcbulgetime.append(bulge_grp.b_pcbulgetime)
                        found_sign = True
                        break
            if found_sign == False:
                b_sign.append(np.nan)
                b_maxconcoffset.append(np.nan)
                b_maxconcoffsetph.append(np.nan)
                b_aveconcoffset.append(np.nan)
                b_totalconcoffset.append(np.nan)
                b_startphaseoffset.append(np.nan)
                b_endphaseoffset.append(np.nan)
                b_totalphaseoffset.append(np.nan)
                b_pcbulgearea.append(np.nan)
                b_pcbulgetime.append(np.nan)

     
if param == 'all':
    if pt == 'nbulges':
        z = nbulges
        min = np.nanmin(z)
        max = np.nanmax(z)
        cm = cm.get_cmap('rainbow', max)
    if pt == 'n+bulges':
        z = nplusbulges
        min = np.nanmin(z)
        max = np.nanmax(z)
        cm = cm.get_cmap('rainbow', max)
    if pt == 'n-bulges':
        z = nnegbulges
        min = np.nanmin(z)
        max = np.nanmax(z)
        cm = cm.get_cmap('rainbow', max)
    if pt == 'aveoffset': 
        #z = sf_diff_s
        z = a_aveconcoffset
        min = np.nanmin(z)
        max = np.nanmax(z)
        
        if np.abs(min) > max:
            max = np.abs(min)
        else:
            min = -max
        
        cm = cm.RdBu_r
    if pt == 'ave+offset': 
        #z = model_ff_s
        z = a_aveplusconcoffset
        min = np.nanmin(z)
        max = np.nanmax(z)
        cm = cm.rainbow
    if pt == 'total+offset': 
        z = a_totalplusconcoffset
        min = np.nanmin(z)
        max = np.nanmax(z)
        cm = cm.rainbow
    if pt == 'ave-offset': 
        z = a_avenegconcoffset
        min = np.nanmin(z)
        max = np.nanmax(z)
        cm = cm.rainbow
    if pt == 'total-offset': 
        z = a_totalnegconcoffset
        min = np.nanmin(z)
        max = np.nanmax(z)
        cm = cm.rainbow
    if pt == 'area+':
        z = a_pcplusbulgearea
        min = 0
        max = 100
        cm = cm.rainbow
    if pt == 'area-':
        z = a_pcnegbulgearea
        min = 0
        max = 100
        cm = cm.rainbow
    if pt == 'time+':
        z = a_pcplusbulgetime
        min = 0
        max = max_ph
        cm = cm.rainbow
    if pt == 'time-':
        z = a_pcplusbulgetime
        min = 0
        max = max_ph
        cm = cm.rainbow
        
if param == 'bulge':
    if pt == 'maxoffset':
        z = b_maxconcoffset
        min = np.nanmin(z)
        max = np.nanmax(z)
        
        if bparam == 'sign':
            if bsign == '+':
                cm = cm.rainbow
            else:
                cm = cm.rainbow
        else:
            cm = cm.rainbow
    
    if pt == 'maxoffsetph':
        z = b_maxconcoffsetph
        min = 0
        max = max_ph
        cm = cm.hsv
    if pt == 'aveoffset': 
        z = b_aveconcoffset
        min = np.nanmin(z)
        max = np.nanmax(z)
        if bparam == 'sign':
            if bsign == '+':
                cm = cm.rainbow
            else:
                cm = cm.rainbow
        else:
            cm = cm.rainbow
    if pt == 'totaloffset': 
        z = b_totalconcoffset
        min = np.nanmin(z)
        max = np.nanmax(z)
        if bparam == 'sign':
            if bsign == '+':
                cm = cm.rainbow
            else:
                cm = cm.rainbow
        else:
            cm = cm.rainbow
    if pt == 'startphoffset': 
        z = b_startphaseoffset
        min = 0
        max = max_ph
        cm = cm.hsv
    if pt == 'endphoffset': 
        z = b_endphaseoffset
        min = 0
        max = max_ph
        cm = cm.hsv
    if pt == 'totalphoffset': 
        z = b_totalphaseoffset
        min = 0
        max = max_ph
        cm = cm.rainbow
    if pt == 'area':
        z = b_pcbulgearea
        min = 0
        max = 100
        if bparam == 'sign':
            if bsign == '+':
                cm = cm.rainbow
            else:
                cm = cm.rainbow
        else:
            cm = cm.rainbow
    if pt == 'time':
        z = b_pcbulgetime
        min = 0
        max = max_ph
        if bparam == 'sign':
            if bsign == '+':
                cm = cm.rainbow
            else:
                cm = cm.rainbow
        else:
            cm = cm.rainbow

    
loc_dict = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic Sites','NA':'North America','SA':'South America'}

#set up plot

latlower_setup = [20,30,20,lat_e[0]]
latupper_setup = [80,72,55,lat_e[-1]]
lonwest_setup = [-170,-15,115,lon_e[0]]
loneast_setup = [-50,35,155,lon_e[-1]]
label_out = ['NA','EU','AS','ZZ']
label = ['NA','EU','AS','ROW']

z = np.array(z)
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
    #parallels = np.arange(-90,91,15)
    #meridians = np.arange(-180,181,30)
    #plt.xticks(meridians)
    #plt.yticks(parallels)
    #m.drawparallels(parallels)
    #m.drawmeridians(meridians)


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
    
    
    current_z = z[test]
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
        all = m.scatter(X[i],Y[i],c=current_z[i], s=m_size, vmin = min,vmax = max, marker = 'o',edgecolor='black',linewidth=0.5,cmap=cm,zorder=10)
            
    ax.text(0.03, 0.97, label[count], transform=ax.transAxes,fontsize=34, fontweight='bold', va='top')

    m.plot(obs_lons, obs_lats, color = 'black', marker = 'o', markersize = 0.001, linestyle= 'None', zorder=1, picker = 5)
    
    all_m.append(m)

    count+=1

plt.tight_layout(pad = 3.08)

fig.subplots_adjust(bottom=0.16)
cbar_ax = fig.add_axes([0.10, 0.08, 0.80, 0.06])
cb = fig.colorbar(all, cax=cbar_ax,orientation ='horizontal')
if param != 'ph':
    cb.set_label('Concentration (ppb)', fontsize = 24)
else:
    if period_type == 'd':
        cb.set_label('Hours', fontsize = 24)
    else:
        cb.set_label('Months', fontsize = 24)

cb.ax.tick_params(labelsize=22)

fig.canvas.mpl_connect('pick_event', interactive)


plt.show()
