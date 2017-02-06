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
import matplotlib as mpl
from matplotlib import animation

model = 'IPSL'
type = '*'

if type == 'SD':
    add = ' SD'
else:
    add = ''

model_fname = '/work/home/db876/modelling/%s/%s_SURFACE_O3_2005_2010_*_*_*_D_%s.nc'%(model,model,type)

root_grp = Dataset(model_fname)
var = root_grp.variables['o3'][:]
var = var*1e9
date = root_grp.variables['date'][:]
time = root_grp.variables['time'][:]
lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
grid_size = root_grp.variables['grid_size'][:]
grid_size = grid_size[0]

if (model == 'UKCA') or (model == 'HADGEM3ES') or (model == 'MIROC3'):
    valid = [0,10,20,31,41,51,59,69,79,90,100,110,120,130,140,151,161,171,181,191,201,212,222,232,243,253,263,273,283,293,304,314,324,334,344,354]
    var = var[valid,:,:]
    date = date[valid]
    var = var[:36,:,:]
    framerate = 2
    
else:
    var = var[:365,:,:]
    framerate = 5

nframes = var.shape[0]

#set up plot

fig = plt.figure(figsize=[16,9])
ax = fig.add_axes([0.02,0.02,0.85,0.96])
ax.set_aspect('equal')
ax.autoscale_view(False)
fig.patch.set_facecolor('white')

cmap = mpl.cm.cubehelix                                                                                                                                                  
minc = 0.545
maxc = 0.555
cnorm = mpl.colors.Normalize(vmin=minc, vmax=maxc)

axcb = fig.add_axes([0.9,0.2,0.04,0.6], frameon=False)
cb=mpl.colorbar.ColorbarBase(axcb, cmap=cmap,norm=cnorm,orientation='vertical')
    
cb.set_label('ppbv', size=18)
cb.ax.tick_params(labelsize=16) 

#LON,LAT = np.meshgrid(lon_c,lat_c)
LON,LAT = np.meshgrid(lon_e,lat_e)

def animator(i):
    print i
    ax.cla()
    
    data = var[i,:,:]
    
    # create basemap
    map = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c')
    # make sure map is associated with ax, not axlogo2
    map.ax = ax
    mcoa = map.drawcoastlines(linewidth=0.25)
    mcou = map.drawcountries(linewidth=0.25)
    #add extra stripe of data to complete sphere
    #contour plot data
    pl = map.pcolor(LON,LAT,data, cmap=cmap, norm=cnorm)
    #pl = map.contourf(LON,LAT,data, cmap=cmap, norm=cnorm,antialiased=True)
    ax.set_title('%s %s %s'%(model+add,date[i],time[i]),fontsize=26)
    
    
ani = animation.FuncAnimation(fig, animator, frames=nframes)

ani.save('CCMI_VIDEOS/CCMI%s%s.mp4'%(model,type), 'ffmpeg',fps=framerate,extra_args=['-vcodec','libx264', '-s','hd1080', '-pix_fmt', 'yuv420p'])
