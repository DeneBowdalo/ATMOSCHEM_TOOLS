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
import matplotlib.animation as animation

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

#get species from current directory
present_dir = os.getcwd()
model_fname,species,start_year,end_year = modules.get_model_info(present_dir)

model_range = '%s_%s'%(start_year,end_year)

print '\nSpecies is %s\n'%(species)

#-------------------------------------------------------
#read in model time series data

#read model netcdf file

model_root_grp = Dataset(model_fname)

model_var = model_root_grp.variables[species.lower()][:]
model_date = model_root_grp.variables['date'][:]
model_time = model_root_grp.variables['time'][:]
lat_c = model_root_grp.variables['lat_centre'][:]
lon_c = model_root_grp.variables['lon_centre'][:]
lat_e = model_root_grp.variables['lat_edges'][:]
lon_e = model_root_grp.variables['lon_edges'][:]
    
n_boxes = len(lat_c)*len(lon_c)

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

z = model_var[:,45,45]

model_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]

time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
ave_pd = pd.Series(z, index=time_pd)

key = lambda x: pd.Period(str(x.year)+'-'+str(x.month))
group=ave_pd.groupby(key)
monthly_ave = group.mean()

monthly_ave.plot(marker = 'x', markersize = 20, color = 'red')
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
linear_lat_inds = []
linear_lons=[]
linear_lon_inds = []
ave_z = np.empty((60,len(lat_c),len(lon_c)))

lat_i = 0
lon_i = 0

key = lambda x: pd.Period(str(x.year)+'-'+str(x.month))

for i in range(n_boxes):
    print i
    current_lat = lat_c[lat_i]
    current_lon = lon_c[lon_i]
    
    cut = model_var[:,lat_i,lon_i]

    ave_pd = pd.Series(cut, index=time_pd)

    group=ave_pd.groupby(key)
    monthly_ave = group.mean()
    monthly_ave = monthly_ave.tolist()
    
    
    ave_z[:,lat_i,lon_i] = monthly_ave

    linear_lats.append(current_lat)
    linear_lat_inds.append(lat_i)
    linear_lons.append(current_lon)
    linear_lon_inds.append(lon_i)
    
    if lon_i == (len(lon_c)-1):
        lat_i+=1
        lon_i=0
    else:
        lon_i+=1


plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)

def data_gen(i):
    pl = m.pcolor(lon_e,lat_e,z[:,linear_lat_inds[i],linear_lon_inds[i]], vmin=0, vmax=50,linewidth=0.5,cmap=plt.cm.coolwarm) 
    #cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
    #cb.set_label('Concentration (ppb)', fontsize = 16)    

    #plt.title('%s for Surface %s between %s'%(type_label,species,model_range),fontsize=20)

    #cb.ax.tick_params(labelsize=16)

    #mng = plt.get_current_fig_manager()
    #mng.window.wm_geometry("+2500+100")

    return pl

ani = animation.FuncAnimation(fig, data_gen, np.arange(0,n_boxes),blit=False)

ani.save('%s_monthly_ave.mp4'%(species), 'ffmpeg',fps=0.5,extra_args=['-vcodec','libx264', '-preset', 'ultrafast', '-s','hd1080', '-pix_fmt', 'yuv420p'])

plt.show()
