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
from pylab import hist, show,subplot


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
    
    ratio = 100./annual_amp[lat_n,lon_n]
    ha_percent = ratio*annual_amp[lat_n,lon_n]
    
    #convert phases to radians
    calc = pi2/24.
    
    calc = pi2/6.
    ha_ph_r = ha_ph[lat_n,lon_n] * calc
    calc = pi2/12.
    annual_ph_r = annual_ph[lat_n,lon_n] * calc
    
    ha_model_wave = ha_amp[lat_n,lon_n]*(np.cos((pi2*model_times/(365.25/2.))-(ha_ph_r)))
    annual_model_wave = annual_amp[lat_n,lon_n]*(np.cos((pi2*model_times/(365.25))-(annual_ph_r)))
    
    ha_primary = p_ha_ph[lat_n,lon_n]
    ha_secondary = s_ha_ph[lat_n,lon_n]
    
    ha_model_wave = ha_model_wave+ave[lat_n,lon_n]
    annual_model_wave = annual_model_wave+ave[lat_n,lon_n]
    
    model_ha_wave_pd = pd.Series(ha_model_wave, index=model_time_pd)
    model_annual_wave_pd = pd.Series(annual_model_wave, index=model_time_pd)
    
    
    fig.canvas.draw()
        
    if first_run == False:
        plt.close(fig2)
        fig2, (axo) = plt.subplots(1,figsize=(24,12))
        fig2.patch.set_facecolor('white')
        
        axo.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='black', markersize = 3, label = 'Observations')
        axo.plot_date(model_time_pd.to_pydatetime(), model_ha_wave_pd, color='green', markersize = 3, label = 'Ha Waveform',markeredgecolor='None')
        axo.plot_date(model_time_pd.to_pydatetime(), model_annual_wave_pd, color='red', markersize = 3, label = 'Annual Waveform',markeredgecolor='None')
        
        #axo.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s\nPrimary HA Phase = %s,Primary HA Regime = %s, HA Amp to Annual Amp Percent = %s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,ha_primary,ha_regime,ha_percent)) 
        
        plt.legend(loc = 'lower right')
        plt.tight_layout()
        axo.grid()
        
        plt.show()
    else:
        fig2, (axo) = plt.subplots(1,figsize=(24,12))
        fig2.patch.set_facecolor('white')
        
        axo.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='black', markersize = 3, label = 'Observations')
        axo.plot_date(model_time_pd.to_pydatetime(), model_ha_wave_pd, color='green', markersize = 3, label = 'Ha Waveform',markeredgecolor='None')
        axo.plot_date(model_time_pd.to_pydatetime(), model_annual_wave_pd, color='red', markersize = 3, label = 'Annual Waveform',markeredgecolor='None')
        
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
site_group = Dataset('model_sig_periods.nc')

daily_amp = site_group.variables['daily_amplitude'][:]
daily_ph = site_group.variables['daily_phase'][:]
ha_amp = site_group.variables['half_annual_amplitude'][:]
ha_ph = site_group.variables['half_annual_phase'][:]
annual_amp = site_group.variables['annual_amplitude'][:]
annual_ph = site_group.variables['annual_phase'][:]
ave = site_group.variables['average'][:]
p_ha_ph = site_group.variables['primary_half_annual_phase'][:]
s_ha_ph = site_group.variables['secondary_half_annual_phase'][:]

ra = 100./annual_amp
daily_annual_amp_ratio = ra*daily_amp
ha_annual_amp_ratio = ra*ha_amp

lat_e = site_group.variables['lat_edges'][:]
lon_e = site_group.variables['lon_edges'][:]
lat_c = site_group.variables['lat_centre'][:]
lon_c = site_group.variables['lon_centre'][:]

n_boxes = len(lat_c)*len(lon_c)

z = np.zeros(np.shape(daily_amp))

print len(lat_c)
half_lat = len(lat_c)/2

#make sure annual phase in SH and NH are eqivalent - take 6 months off SH Annual phase
print annual_ph.shape
annual_ph[:half_lat,:] = annual_ph[:half_lat,:] - 6.
test = annual_ph < 0
annual_ph[test] = 12.- np.abs(annual_ph[test])

#p_ha_ph[:half_lat,:] = p_ha_ph[:half_lat,:] - 6.
#test = p_ha_ph < 0
#p_ha_ph[test] = 12.- np.abs(p_ha_ph[test])

start_spring = (12/365.25)*79.25
start_summer = (12/365.25)*172.25
start_autumn = (12/365.25)*265.25
start_winter = (12/365.25)*355.25

mid_spring = start_spring + ((start_summer-start_spring)/2.)
mid_summer = start_summer + ((start_autumn-start_summer)/2.)
mid_autumn = start_autumn + ((start_winter-start_autumn)/2.)

winter_diff = 12.-((12/365.25)*355.25)
winter_diff = (winter_diff+start_spring)/2.
mid_winter = start_spring - winter_diff    

print 'start spring = ', start_spring
print 'mid spring = ', mid_spring
print 'start summer = ', start_summer
print 'mid summer = ', mid_summer
print 'start autumn = ', start_autumn
print 'mid autumn = ', mid_autumn
print 'start winter = ', start_winter
print 'mid winter = ', mid_winter

lat_i = 0
lon_i = 0

for siten in range(n_boxes):

#get continental polluted o3 group

#strong annual cycle, that broadly peaks in summer
#make sure annual phase peak is between mid-spring to mid-autumn
#make sure ha/annual cycle ratio is less than 50%
#make sure daily/annual cycle ratio is greater than 35%
    if (annual_ph[lat_i,lon_i] >= start_summer) & (annual_ph[lat_i,lon_i] < start_autumn):
        if (annual_amp[lat_i,lon_i] >= 10):
            if (daily_ph[lat_i,lon_i] >= 12) & (daily_ph[lat_i,lon_i] < 21):
                if ha_annual_amp_ratio[lat_i,lon_i] < 50:
                    z[lat_i,lon_i] = 1
             
            
#get oceanic o3 group

#Broad annual winter/spring peak
#make sure annual phase peak is between start of winter to end of spring
#make sure ha/annual cycle ratio is less than 50%
#make sure daily/annual cycle ratio is less than 35%    
    if (annual_ph[lat_i,lon_i] >= start_winter) or (annual_ph[lat_i,lon_i] < start_spring): 
        if (daily_amp[lat_i,lon_i] < 5):  
            if(annual_amp[lat_i,lon_i] >= 10):
                if (daily_ph[lat_i,lon_i] >= 21) or (daily_ph[lat_i,lon_i] < 12):
                    if ha_annual_amp_ratio[lat_i,lon_i] < 50:
                        z[lat_i,lon_i] = 2
            
#get continental clean o3 group

#strong annual cycle, that broadly peaks in summer
#make sure annual phase peak is between mid-spring to mid-autumn
#make sure ha/annual cycle ratio is less than 50%
#make sure daily/annual cycle ratio is greater than 35%
    if (annual_ph[lat_i,lon_i] >= mid_winter) or (annual_ph[lat_i,lon_i] < start_summer):
        if (annual_amp[lat_i,lon_i] < 10):
            if ha_annual_amp_ratio[lat_i,lon_i] < 50:
                z[lat_i,lon_i] = 3

#get intermediate o3 group

# V weak annual cycle. Half-annual cycle is dominant
#make sure primary peak half annual phase is between start of spring to mid-summer
#make sure ha/annual cycle ratio is greater than 50%
    if (annual_ph[lat_i,lon_i] >= start_winter) & (annual_ph[lat_i,lon_i] < start_summer):
        if (daily_ph[lat_i,lon_i] >= 12) & (daily_ph[lat_i,lon_i] < 21):
            if ha_annual_amp_ratio[lat_i,lon_i] >= 50:
                z[lat_i,lon_i] = 4




    if lon_i == (len(lon_c)-1):
        lat_i+=1
        lon_i=0
    else:
        lon_i+=1


#z = annual_ph

#plt.hist(np.ravel(annual_amp),bins=30)
#plt.show()

#------------------------------------------
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



pl = m.pcolor(lon_e,lat_e,z, vmin=np.min(z), vmax=np.max(z),linewidth=0.5,picker = 5)
cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
cb.set_label('Group Type', fontsize = 16)    


plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)


#plt.title('%s for Surface %s between %s'%(type_label,species,model_range),fontsize=20)

cb.ax.tick_params(labelsize=16)

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500+100")

fig.canvas.mpl_connect('pick_event', onpick)

plt.show()