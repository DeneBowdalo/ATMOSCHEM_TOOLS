#mport matplotlib
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
from scipy import stats
import matplotlib.dates as dates
from matplotlib import animation

species = raw_input('O3, CO, NO, or NO2?\n') 

start_year = 2005
end_year = 2010

start = datetime.datetime(year = start_year, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = end_year, month = 1, day = 1, hour = 0, minute = 0)

full_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]
full_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='H')][:-1]
full_time_pd = pd.date_range(start = start,end = end, freq = 'H')[:-1] 

#------------------------------------------------------------
#read in obs time series data
obs_grp = Dataset('/work/home/db876/observations/surface/%s/process/GLOBAL_SURFACE_%s_1970_2015_H_ALL.nc'%(species,species))

obs_refs_dict = obs_grp.groups

obs_refs = []
obs_lats = []
obs_lons = []
obs_alt = []
obs_data = []
obs_country = []

for i in obs_refs_dict.keys():
    i = i.encode('ascii')
    obs_refs = np.append(obs_refs,i)

for ref in obs_refs:
    obs_site_group = obs_grp.groups[ref] 
    obs_country = np.append(obs_country,obs_site_group.country)
    obs_data.append(obs_site_group.variables['%s'%(species.lower())][:])
    obs_lats = np.append(obs_lats,obs_site_group.latitude)
    obs_lons = np.append(obs_lons,obs_site_group.longitude)
    obs_alt = np.append(obs_alt,obs_site_group.altitude)
obs_date = obs_site_group.variables['date'][:]
obs_time = obs_site_group.variables['time'][:]

for i in range(len(obs_refs)):
    obs_refs[i] = obs_refs[i].lower()

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
obs_tags = modules.get_tags(obs_refs)

obs_data = np.array(obs_data)

#--------------------------------------------------
areas = ['ANT','S_O','S_NA','S_EU','S_AS','SW_NA','C_EU','SE_AS','SE_NA','E_EU','C_AS','CE_NA','NW_EU','NE_AS','C_NA','N_EU','OC','NW_NA','AF','NE_NA','SA','AL','N_O','ARC','NO GIVEN AREA'] 
#areas = ['NE_NA','CE_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','N_O','NO GIVEN AREA'] 

area_boundaries,area_tags,area_labels = modules.area_dicts()

cm = plt.cm.get_cmap('Set1')

obs_areas = []
obs_color_ratios = []
for i in range(len(obs_lats)):
    obs_areas.append(modules.get_area(area_boundaries,obs_tags[i],obs_lats[i],obs_lons[i],obs_refs[i])) 

area_color_ratios = np.linspace(0.001,0.999,len(areas))



for i in range(len(obs_areas)):
    obs_color_ratios.append(area_color_ratios[areas.index(obs_areas[i])])
obs_color_ratios = np.array(obs_color_ratios)

start_boolean = (obs_date ==  int(full_dates[0])) & (obs_time ==  int(full_times[0])) 
end_boolean = (obs_date ==  int(full_dates[-1])) & (obs_time ==  int(full_times[-1]))
start_obs_ind = np.where(start_boolean == True)[0][0]
end_obs_ind = (np.where(end_boolean == True)[0][0]) +1

obs_data = obs_data[:,start_obs_ind:end_obs_ind]

fig = plt.figure(figsize=[16,9])
ax = fig.add_axes([0.02,0.02,0.96,0.96])
ax.set_aspect('equal')
ax.autoscale_view(False)
fig.patch.set_facecolor('white')

# create basemap
map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
map.ax = ax
mcoa = map.drawcoastlines(linewidth=0.25)
mcou = map.drawcountries(linewidth=0.25)
        
for j in range(len(obs_lats)):
    x, y = np.meshgrid(*map(obs_lons[j], obs_lats[j]))
    map.scatter(x, y, c=cm(obs_color_ratios[j]),s=35)
        
for x in range(len(areas)):
    map.plot([-999],[-999],color=cm(area_color_ratios[x]),marker='o',markersize=9,label=areas[x])
          
l = plt.legend(prop={'size':9},loc=(0.05,0.02))          
ax.set_title('HOURLY %s AREAS'%(species),fontsize=26)

plt.show()


