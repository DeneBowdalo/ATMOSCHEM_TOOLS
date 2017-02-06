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

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-2]
print '\nSpecies is %s\n'%(species)

#read in obs data
obs_files = glob.glob('../process/GLOBAL_SURFACE*')

start_year = raw_input('Start Year?\n')
end_year = raw_input('\nEnd Year?\n')
    
timeres = raw_input('\nH or M?\n')

root_grp = Dataset('../process/GLOBAL_SURFACE_%s_%s_%s_%s_PERIODIC.nc'%(species,start_year,end_year,timeres))

valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

obs_aves = []
obs_lats = []
obs_lons = []
obs_alts = []
obs_groups = []

for ref in valid_refs:  
    #read in specific site data
    site_group = root_grp.groups[ref]

    #read in variables for site
    obs_var = site_group.variables[species.lower()][:]
    obs_lat = site_group.latitude
    obs_lon = site_group.longitude
    obs_alt = site_group.altitude
    obs_group = site_group.process_group
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_ave = np.ma.average(obs_var_mask)
    
    obs_aves.append(obs_ave)
    obs_lats.append(obs_lat)
    obs_lons.append(obs_lon)
    obs_alts.append(obs_alt)
    obs_groups.append(obs_group)

#set up plot
fig =plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')
 
#setup basemap projection
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                   llcrnrlon=-180,\
                    urcrnrlon=180,\
                   resolution='c')


m.drawcoastlines()
m.drawmapboundary()
#parallels = np.arange(-90,91,15)
#meridians = np.arange(-180,181,30)
#plt.xticks(meridians)
#plt.yticks(parallels)
#m.drawparallels(parallels)
#m.drawmeridians(meridians)

X,Y = m(obs_lons,obs_lats)

pl = m.scatter(X,Y,c=obs_aves,s=20, vmin=np.min(obs_aves), vmax=np.max(obs_aves),edgecolor='black',linewidth=0.5,cmap=plt.cm.coolwarm)
cb = plt.colorbar(pl, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')

cb.set_label('Concentration (ppb)', fontsize = 16)    
plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)
plt.title('Observational average surface %s between %s:%s'%(species,start_year,end_year),fontsize=20)
cb.ax.tick_params(labelsize=16)

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500+100")
plt.show()




