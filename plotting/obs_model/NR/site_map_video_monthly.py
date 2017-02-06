from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from collections import OrderedDict 
import operator
import matplotlib as mpl
from matplotlib import animation
import os

species = os.getcwd().split('/')[-2]

#read in obs data
root_grp = Dataset('../process/GLOBAL_SURFACE_%s_1970_2015_M_ALL_NR.nc'%(species))
valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

print valid_refs[0]

site_grp = root_grp.groups[valid_refs[0]]
full_dates = site_grp.variables['date'][:]
full_times = site_grp.variables['time'][:]

fig = plt.figure(figsize=[16,9])
ax = fig.add_axes([0.02,0.02,0.85,0.96])
ax.set_aspect('equal')
ax.autoscale_view(False)
fig.patch.set_facecolor('white')

framerate = 5
nframes = len(full_dates)

cmap = mpl.cm.cubehelix                                                                                                                                                  

if species == 'O3':
    minc = 0
    maxc = 70
elif (species == 'NO') or (species == 'NO2'):
    minc = 0
    maxc = 20
elif species == 'CO':
    minc = 0
    maxc = 500

cnorm = mpl.colors.Normalize(vmin=minc, vmax=maxc)

axcb = fig.add_axes([0.9,0.2,0.04,0.6], frameon=False)
cb=mpl.colorbar.ColorbarBase(axcb, cmap=cmap,norm=cnorm,orientation='vertical')
    
cb.set_label('ppbv', size=18)
cb.ax.tick_params(labelsize=16) 


def animator(i):
    print i
    ax.cla()
    
    # create basemap
    map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
    map.ax = ax
    mcoa = map.drawcoastlines(linewidth=0.25)
    mcou = map.drawcountries(linewidth=0.25)
    #add extra stripe of data to complete sphere
    
    for site_ref in valid_refs:
        site_group = root_grp.groups[site_ref]
        obs_cut = site_group.variables[species.lower()][i]
        obs_lat = site_group.latitude
        obs_lon = site_group.longitude
        obs_group = site_group.process_group
        
        if obs_cut != -99999:        
            x, y = np.meshgrid(*map(obs_lon, obs_lat))
            map.scatter(x, y, c=obs_cut,s=35,cmap=cmap, norm=cnorm)

    ax.set_title('MONTHLY %s %s %s'%(species,full_dates[i],full_times[i]),fontsize=26)
    
print nframes

ani = animation.FuncAnimation(fig, animator, frames=nframes)

ani.save('%s_MONTHLY_video.mp4'%(species), 'ffmpeg',fps=framerate,extra_args=['-vcodec','libx264', '-s','hd1080', '-pix_fmt', 'yuv420p'])
