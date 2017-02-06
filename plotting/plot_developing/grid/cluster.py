import modules
import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm
from cmath import *
import colorsys
from scipy import signal
from netCDF4 import Dataset
import datetime
import pandas as pd
from scipy.cluster.vq import vq, kmeans, whiten
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn import metrics

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-3]
model_range = paths[-2]
start_year = model_range[:4]

print '\nSpecies is %s\n'%(species)

model_details = paths[-1]
print model_details

model_split = model_details.split('_SFC')
model = model_split[0]
model_other = model_split[1]
other_split = model_other.split('_') 


#define module for interactive plotting
def onpick(event):
    global pl
    
    global ind
    global fig2
    
    ind = event.ind
    print 'ind = ',ind
    ind = ind[0]
    #x_data = event.xdata
    #y_data = event.ydata

    #find ind of closest lat/lon
    #ind = modules.find_nearest_point_index(obs_lons,obs_lats,x_data,y_data)
    
    try:
        for i in range(len(pl)):
            pl.pop(0).remove()
            first_run = False  
        
    except:
        first_run = True
        pass
    
    pl = m.plot([linear_lons[ind]], [linear_lats[ind]], 's', ms=20, alpha=0.6, color='yellow',zorder=20)

    
    #get model timeseries for site clicked
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,linear_lats[ind],linear_lons[ind])
    model_var_pick = model_var[:,lat_n,lon_n]
    model_var_pick = model_var_pick*1e9
    model_var_mask = np.ma.masked_where(model_var_pick<=0,model_var_pick)
    model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
    model_var_pd = pd.Series(model_var_mask, index=model_time_pd)

        
    #create sine wave from amp/phase
    model_date_l = model_date.astype(int)
    model_time_l = model_time.astype(int)
    model_times = modules.date_process(model_date_l,model_time_l,start_year)
    model_times = np.array(model_times)
    pi2 = np.pi*2
    
    #ratio = 100./o3_annual_amp[lat_n,lon_n]
    #ha_percent = ratio*o3_annual_amp[lat_n,lon_n]
    
    #convert phases to radians
    #calc = pi2/24.
    
    #calc = pi2/6.
    #ha_ph_r = ha_ph[lat_n,lon_n] * calc
    #calc = pi2/12.
    #annual_ph_r = annual_ph[lat_n,lon_n] * calc
    
    #ha_model_wave = ha_amp[lat_n,lon_n]*(np.cos((pi2*model_times/(365.25/2.))-(ha_ph_r)))
    #annual_model_wave = annual_amp[lat_n,lon_n]*(np.cos((pi2*model_times/(365.25))-(annual_ph_r)))
    
    #ha_primary = p_ha_ph[lat_n,lon_n]
    #ha_secondary = s_ha_ph[lat_n,lon_n]
    
    #ha_model_wave = ha_model_wave+ave[lat_n,lon_n]
    #annual_model_wave = annual_model_wave+ave[lat_n,lon_n]
    
    #model_ha_wave_pd = pd.Series(ha_model_wave, index=model_time_pd)
    #model_annual_wave_pd = pd.Series(annual_model_wave, index=model_time_pd)
    
    
    fig.canvas.draw()
        
    if first_run == False:
        plt.close(fig2)
        fig2, (axo) = plt.subplots(1,figsize=(24,12))
        fig2.patch.set_facecolor('white')
        
        axo.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='black', markersize = 3, label = 'Observations')
        #axo.plot_date(model_time_pd.to_pydatetime(), model_ha_wave_pd, color='green', markersize = 3, label = 'Ha Waveform',markeredgecolor='None')
        #axo.plot_date(model_time_pd.to_pydatetime(), model_annual_wave_pd, color='red', markersize = 3, label = 'Annual Waveform',markeredgecolor='None')
        
        #axo.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s\nPrimary HA Phase = %s,Primary HA Regime = %s, HA Amp to Annual Amp Percent = %s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,ha_primary,ha_regime,ha_percent)) 
        
        plt.legend(loc = 'lower right')
        plt.tight_layout()
        axo.grid()
        
        plt.show()
    else:
        fig2, (axo) = plt.subplots(1,figsize=(24,12))
        fig2.patch.set_facecolor('white')
        
        axo.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='black', markersize = 3, label = 'Observations')
        #axo.plot_date(model_time_pd.to_pydatetime(), model_ha_wave_pd, color='green', markersize = 3, label = 'Ha Waveform',markeredgecolor='None')
        #axo.plot_date(model_time_pd.to_pydatetime(), model_annual_wave_pd, color='red', markersize = 3, label = 'Annual Waveform',markeredgecolor='None')
        
        #axo.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s\nPrimary HA Phase = %s,Primary HA Regime = %s, HA Amp to Annual Amp Percent = %s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,ha_primary,ha_regime,ha_percent))
        
        plt.legend(loc = 'lower right')
        plt.tight_layout()
        axo.grid()
        
        plt.show()

print '\nSpecies is %s\n'%(species)

#-------------------------------------------------------
#read in model time series data

if len(other_split) == 1:
    version = ''
    grid = ''
    met = ''
    model_f = '/work/home/db876/plotting_tools/model_files/'+model+'_'+'SURFACE'+'_'+species+'_'+model_range+'.nc'

elif len(other_split) == 2:
    version = other_split[1] 
    grid = ''
    met = ''
    model_f = '/work/home/db876/plotting_tools/model_files/'+model+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'.nc'

elif len(other_split) == 3:
    version = other_split[1] 
    grid = other_split[2]
    met = ''
    model_f = '/work/home/db876/plotting_tools/model_files/'+model+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'_'+grid+'.nc'

elif len(other_split) == 4:
    version = other_split[1] 
    grid = other_split[2]
    met = other_split[3]
    model_f = '/work/home/db876/plotting_tools/model_files/'+model+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'_'+grid+'_'+met+'.nc'
    
print model_f
root_grp = Dataset(model_f)
model_var = root_grp.variables[species.lower()][:]
model_date = root_grp.variables['date'][:]
model_time = root_grp.variables['time'][:]
lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
grid_size = root_grp.variables['grid_size'][:]
grid_size = grid_size[0]

n_boxes = len(lon_c)*len(lat_c)

#process model dates and model times to datetimes, then process pandas objects

year_val = []
month_val = []
day_val = []
hour_val = []
minute_val = []

#change model times to datetimes
model_date = model_date.astype('str')
model_time = model_time.astype('str')

for date in model_date:
    year_val.append(int(date[0:4]))
    month_val.append(int(date[4:6]))
    day_val.append(int(date[6:8]))

for time in model_time:
    if np.float64(time) == 0:
        hour_val.append(0)
        minute_val.append(0)
    elif np.float64(time) == 30:
        hour_val.append(0)
        minute_val.append(30)
    else:
        hour_val.append(int(time[0:-2]))
        minute_val.append(int(time[-2:]))

model_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]

#--------------------------------------------
#
o3_site_group = Dataset('/work/home/db876/global_grid/O3/2005_2010/GEOS_CHEM_SFC_v90103_2x2.5_GEOS5/model_sig_periods.nc')
no_site_group = Dataset('/work/home/db876/global_grid/NO/2005_2010/GEOS_CHEM_SFC_v90103_2x2.5_GEOS5/model_sig_periods.nc')
no2_site_group = Dataset('/work/home/db876/global_grid/NO2/2005_2010/GEOS_CHEM_SFC_v90103_2x2.5_GEOS5/model_sig_periods.nc')
co_site_group = Dataset('/work/home/db876/global_grid/CO/2005_2010/GEOS_CHEM_SFC_v90103_2x2.5_GEOS5/model_sig_periods.nc')
isop_site_group = Dataset('/work/home/db876/global_grid/ISOP/2005_2010/GEOS_CHEM_SFC_v90103_2x2.5_GEOS5/model_sig_periods.nc')
temp_site_group = Dataset('/work/home/db876/global_grid/GMAO_TEMP/2005_2010/GEOS_CHEM_SFC_v90103_2x2.5_GEOS5/model_sig_periods.nc')
ws_site_group = Dataset('/work/home/db876/global_grid/WIND_SPEED/2005_2010/GEOS_CHEM_SFC_v90103_2x2.5_GEOS5/model_sig_periods.nc')

# o3_annual_amp = np.ravel(o3_site_group.variables['annual_amplitude'][:])
# no_annual_amp = np.ravel(no_site_group.variables['annual_amplitude'][:])
# no2_annual_amp = np.ravel(no2_site_group.variables['annual_amplitude'][:])
# #co_annual_amp = np.ravel(co_site_group.variables['annual_amplitude'][:])
# #isop_annual_amp = np.ravel(isop_site_group.variables['annual_amplitude'][:])
# temp_annual_amp = np.ravel(temp_site_group.variables['annual_amplitude'][:])
# ws_annual_amp = np.ravel(ws_site_group.variables['annual_amplitude'][:])

o3_annual_amp = np.ravel(o3_site_group.variables['annual_amplitude'][:])
o3_ha_amp = np.ravel(o3_site_group.variables['half_annual_amplitude'][:])
o3_daily_amp = np.ravel(o3_site_group.variables['daily_amplitude'][:])
o3_daily_ph = np.ravel(o3_site_group.variables['daily_phase'][:])
o3_annual_ph = o3_site_group.variables['annual_phase'][:]
o3_ha_ph = o3_site_group.variables['half_annual_phase'][:]

half_lat = len(lat_c)/2
#make sure annual phase in SH and NH are eqivalent - take 6 months off SH Annual phase
o3_ha_ph[:half_lat,:] = o3_ha_ph[:half_lat,:] - 6.
test = o3_ha_ph < 0
o3_ha_ph[test] = 12.- np.abs(o3_ha_ph[test])
o3_ha_ph = np.ravel(o3_ha_ph)


half_lat = len(lat_c)/2
#make sure annual phase in SH and NH are eqivalent - take 6 months off SH Annual phase
o3_annual_ph[:half_lat,:] = o3_annual_ph[:half_lat,:] - 6.
test = o3_annual_ph < 0
o3_annual_ph[test] = 12.- np.abs(o3_annual_ph[test])
o3_annual_ph = np.ravel(o3_annual_ph)

no_annual_amp = np.ravel(no_site_group.variables['annual_amplitude'][:])
no2_annual_amp = np.ravel(no2_site_group.variables['annual_amplitude'][:])
#co_annual_amp = np.ravel(co_site_group.variables['annual_amplitude'][:])
#isop_annual_amp = np.ravel(isop_site_group.variables['annual_amplitude'][:])
temp_annual_amp = np.ravel(temp_site_group.variables['annual_amplitude'][:])
ws_annual_amp = np.ravel(ws_site_group.variables['annual_amplitude'][:])

#o3_annual_amp = (o3_annual_amp-np.average(o3_annual_amp))/np.std(o3_annual_amp)
#no_annual_amp = (no_annual_amp-np.average(no_annual_amp))/np.std(no_annual_amp)
#no2_annual_amp = (no2_annual_amp-np.average(no2_annual_amp))/np.std(no2_annual_amp)
#temp_annual_amp = (temp_annual_amp-np.average(temp_annual_amp))/np.std(temp_annual_amp)

#o3_annual_amp = (o3_annual_amp - np.ma.average(o3_annual_amp)) / (np.ma.max(o3_annual_amp) - np.ma.min(o3_annual_amp))
#no_annual_amp = (no_annual_amp - np.ma.average(no_annual_amp)) / (np.ma.max(no_annual_amp) - np.ma.min(no_annual_amp))
#no2_annual_amp = (no2_annual_amp - np.ma.average(no2_annual_amp)) / (np.ma.max(no2_annual_amp) - np.ma.min(no2_annual_amp))
#co_annual_amp = (co_annual_amp - np.ma.average(co_annual_amp)) / (np.ma.max(co_annual_amp) - np.ma.min(co_annual_amp))
#isop_annual_amp = (isop_annual_amp - np.ma.average(isop_annual_amp)) / (np.ma.max(isop_annual_amp) - np.ma.min(isop_annual_amp))




all_array = np.vstack((o3_annual_amp,o3_annual_ph,o3_daily_amp,o3_daily_ph,o3_ha_ph,o3_ha_amp))
all_array = np.transpose(all_array)


# computing K-Means with K = 4 (4 clusters)
centroids,_ = kmeans(all_array,6)
# assign each sample to a cluster
idx,_ = vq(all_array,centroids)

#centroids,_ = kmeans(all_annual_amp,5)
#db = DBSCAN(eps=0.1, min_samples=5).fit(all_annual_amp)
#core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
#core_samples_mask[db.core_sample_indices_] = True
#labels = db.labels_

print idx

idx = np.reshape(idx,(len(lat_c),len(lon_c)))

#-----------------------------------------
#set up plot
fig =plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')
 
#setup basemap projection
m = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c')


m.drawcoastlines()
m.drawmapboundary()
parallels = np.arange(-90,91,15)
meridians = np.arange(-180,151,30)
plt.xticks(meridians)
plt.yticks(parallels)
m.drawparallels(parallels)
m.drawmeridians(meridians)

linear_lats = []
linear_lons=[]
lat_i = 0
lon_i = 0

for i in range(n_boxes):
    current_lat = lat_c[lat_i]
    current_lon = lon_c[lon_i]

    linear_lats.append(current_lat)
    linear_lons.append(current_lon)
    
    if lon_i == (len(lon_c)-1):
        lat_i+=1
        lon_i=0
    else:
        lon_i+=1


pl = m.pcolor(lon_e,lat_e,idx,linewidth=0.5,cmap=plt.cm.jet,picker = 5)
cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
#cb.set_label('Concentration (ppb)', fontsize = 16)    

plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)

#plt.title('%s %s for Surface %s between %s'%(period,type_label,species,model_range),fontsize=20)


#cb.ax.tick_params(labelsize=16)

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500+100")

fig.canvas.mpl_connect('pick_event', onpick)

plt.show()

# Black removed and is used for noise instead.
unique_labels = set(labels)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    class_member_mask = (labels == k)

    
    
    xy = all_annual_amp[class_member_mask & core_samples_mask]
    print xy[:,0]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)

    #xy = all_annual_amp[class_member_mask & ~core_samples_mask]
    #plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
     #        markeredgecolor='k', markersize=6)

plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.show()

#plt.plot(all_annual_amp[idx==0,0],all_annual_amp[idx==0,1])#,all_annual_amp[idx==0,2],all_annual_amp[idx==0,3],all_annual_amp[idx==0,4],'ob') 
#plt.plot(all_annual_amp[idx==1,0],all_annual_amp[idx==1,1])#,all_annual_amp[idx==1,2],all_annual_amp[idx==1,3],all_annual_amp[idx==1,4],'or')
#plt.plot(all_annual_amp[idx==2,0],all_annual_amp[idx==2,1])#,all_annual_amp[idx==2,2],all_annual_amp[idx==2,3],all_annual_amp[idx==2,4],'og')

#plt.show()

