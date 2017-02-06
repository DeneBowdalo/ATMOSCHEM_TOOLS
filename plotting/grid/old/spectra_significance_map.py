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

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-3]
years = paths[-2]
start_year = years[:4]
end_year = years[5:9]
model_path = paths[-1] 


data_split = model_path.split('_SFC') 
print data_split
model = data_split[0]
#model = model[4:]

only_model = False
model_version = False
model_version_grid = False
model_version_grid_met = False

other = data_split[1]
other_split = other.split('_')
other_split.remove('')
if other_split == []:
    print 'only 1 version of model'
    only_model = True
elif len(other_split) == 1:
    version = other_split[0]
    model_version = True
elif len(other_split) == 2:
    version = other_split[0]
    grid = other_split[1]
    model_version_grid = True
elif len(other_split) == 3:
    version = other_split[0]
    grid = other_split[1]    
    met = other_split[2]
    model_version_grid_met = True

print start_year, end_year
print model 


print '\nSpecies is %s\n'%(species)


#--------------------------------------------
#
site_group = Dataset('model_spectra_significance.nc')

#read in variables for site
period = site_group.variables['period'][:]

#get period
period_chose = np.float64(raw_input('Choose Period (in Days)\n'))

the_ind = min(range(len(period)), key=lambda i: abs(period[i]-period_chose))

significance = site_group.variables['chi_significance_points'][:,:,the_ind]
amplitude = site_group.variables['amplitude'][:,:,the_ind]
tau = site_group.variables['tau'][:,:]
lat_e = site_group.variables['lat_edges'][:]
lon_e = site_group.variables['lon_edges'][:]
lat_c = site_group.variables['lat_centre'][:]
lon_c = site_group.variables['lon_centre'][:]

#amp or phase
data_type = raw_input('\nsig or tau?\n')
if data_type == 'sig':
    z = significance
    type_label = 'Chi significance (%)'
else:
    z = tau 
    type_label = 'Tau'


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

if data_type == 'sig':
    pl = m.pcolor(lon_e,lat_e,z, vmin=95, vmax=100,linewidth=0.5,cmap=plt.cm.coolwarm)
    cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
    #cb.set_label('Concentration (ppb)', fontsize = 16)    
else:
    pl = m.pcolor(lon_e,lat_e,z, vmin=0, vmax=20,linewidth=0.5,cmap=plt.cm.coolwarm)
    cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
    #cb.set_label('Day', fontsize = 16)

plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)
plt.title('%s Day %s for Surface %s between %s_%s'%(period_chose,type_label,species,start_year,end_year),fontsize=20)
cb.ax.tick_params(labelsize=16)

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500+100")
plt.show()