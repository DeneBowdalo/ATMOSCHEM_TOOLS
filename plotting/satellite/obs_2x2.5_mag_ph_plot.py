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

#Type should be Magnitude or Phase
type = raw_input('mag or phase?\n')

#Period should be daily, half_annual, annual
period = raw_input('half_annual or annual?\n')

#Ask if want to save out full image, no borders
full_image = raw_input('Do you want to save out full image? Y or N?\n')
if full_image == 'Y':
    f_name = 'omi_%s_%s.png'%(period,type)

if type == 'mag':
    type_f = 'magnitudes'
    
else:
    type_f = 'phases'

#set up plot
fig =plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')

if full_image == 'Y':
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

model_f = 'obs_%s_2x2.5/obs_%s_%s.npy'%(type_f,period,type_f)

if period == 'half_annual':
    if type_f == 'magnitudes':
        title = 'Half-Annual Magnitude' 
        label = 'Concentration (ppbV)'  
        types = 'mag'
    else:
        title = 'Half-Annual Phase'
        label = 'Time (Months)'
        types = 'phase'
        phase_min = 0
        phase_max = 6

if period == 'annual':
    if type_f == 'magnitudes':
        title = 'Annual Magnitude'
        label = 'Concentration (ppbV)'
        types = 'mag'
    else:
        title = 'Annual Phase'
        label = ' Time (Months)'
        types = 'phase'
        phase_min = 0
        phase_max = 12
	
# load in model values
values = np.load(model_f)
print values.shape

#lat lon edges for 2x5 grid
lat_e = np.arange(-59.,60,2)
lon_e = np.arange(-181.25,179,2.5)
lat_c = np.arange(-58.,59,2)
lon_c = np.arange(-180,178,2.5)
#get size of grid
grid_dim_c = len(lat_c)*len(lon_c)
 
#reshape array to plot on map
start = 0
end = len(lon_c)

for i in range(len(lat_c)):
	new_list = values[start:end]
	new_list = np.array(new_list)
	try:
		z =np.vstack((z,new_list))
	except:
		z = [new_list]
		z=np.array(z)
	start+=len(lon_c)
	end+=len(lon_c)

#setup basemap projection
m = Basemap(projection='cyl',llcrnrlat=-59,urcrnrlat=59,\
                   llcrnrlon=lon_e[0],\
                    urcrnrlon=lon_e[-1],\
                   resolution='c')


if full_image == 'N':
    m.drawcoastlines()
    m.drawmapboundary()
    parallels = np.arange(-50,51,20)
    meridians = np.arange(-180,151,30)
    plt.xticks(meridians)
    plt.yticks(parallels)
    m.drawparallels(parallels)
    m.drawmeridians(meridians)


#if standard dev of gridbox > certain num , mask it by making it black
stdevs = np.load('../process/OMI_O3_trop_ave_2x2.5_stdev.npy')

lat = np.arange(59)
lon = np.arange(144)

for la in lat:
    for lo in lon:
        data = np.average(stdevs[:,int(la),int(lo)])
        if data > 2.5:
            try:
                mask_lats=np.vstack((mask_lats,[lat_e[la],lat_e[la+1]]))
                mask_lons=np.vstack((mask_lons,[lon_e[lo],lon_e[lo+1]]))
            except:
                mask_lats = np.array([lat_e[la],lat_e[la+1]])
                mask_lons = np.array([lon_e[lo],lon_e[lo+1]])



z_stdev = np.zeros((len(mask_lats),1,1))


#plot model gridboxes
if type == 'mag':
	poly = m.pcolor(lon_e, lat_e, z,vmin=0, vmax=np.max(z),cmap = plt.cm.coolwarm)
	#for i in range(len(mask_lons)):
    #        m.pcolor(mask_lons[i], mask_lats[i], z_stdev[i],cmap= plt.cm.gray)
else:
    poly = m.pcolor(lon_e, lat_e, z, vmin=phase_min, vmax=phase_max, cmap=plt.cm.hsv)
    #for i in range(len(mask_lons)):
    #        m.pcolor(mask_lons[i], mask_lats[i], z_stdev[i],cmap= plt.cm.gray)

if full_image == 'N':
    if (period == 'half_annual') & (type == 'phase'):
        cb = plt.colorbar(poly, ticks =[0,1,2,3,4,5],  ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
        cb.ax.set_xticklabels(['M1','M2','M3','M4','M5','M6'])
        
    elif (period == 'annual') & (type == 'phase'):
        cb = plt.colorbar(poly, ticks =[0,1.01848,1.94661,2.96509,3.95072,4.96920,5.95483,6.97331,7.99179,8.97741,9.99589,10.98152],  ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
        cb.ax.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
    else:
        cb = plt.colorbar(poly, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
    cb.set_label('%s'%(label), fontsize = 16)
    plt.xlabel('Longitude',fontsize = 20)
    plt.ylabel('Latitude',fontsize = 20)
    cb.ax.tick_params(labelsize=16)
    plt.title('Tropospheric O3 OMI 2x5 Oct 2004 - Jan 2014 %s'%(title), fontsize = 18)
    plt.show()
else:
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    plt.savefig(f_name,bbox_inches=extent)



