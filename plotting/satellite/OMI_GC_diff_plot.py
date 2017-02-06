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

model_f = 'model_%s/model_%s_%s.npy'%(type_f,period,type_f)
obs_f = 'obs_%s_2x2.5/obs_%s_%s.npy'%(type_f,period,type_f)

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
obs_values = np.load(obs_f)
model_values = np.load(model_f)

#lat lon edges for model 2x5 grid
obs_lat_c = np.arange(-58.,59.,2)
obs_lon_c = np.arange(-180,178,2.5)
obs_lat_e = np.arange(-59.,60,2)
obs_lon_e = np.arange(-181.25,179,2.5)

#lat lon edges for model 2x5 grid
model_lat_c = np.arange(-88.,89.,2)
model_lat_c = np.insert(model_lat_c,0,-89.5)
model_lat_c = np.append(model_lat_c,89.5)
model_lon_c = np.arange(-180,178,2.5)
model_lat_e = np.arange(-89.,90,2)
model_lat_e = np.insert(model_lat_e,0,-90.)
model_lat_e = np.append(model_lat_e,90.)
model_lon_e = np.arange(-181.25,179,2.5)


#reshape obs array to plot on map
start = 0
end = len(obs_lon_c)

for i in range(len(obs_lat_c)):
	new_list = obs_values[start:end]
	new_list = np.array(new_list)
	try:
		z_obs =np.vstack((z_obs,new_list))
	except:
		z_obs = [new_list]
		z_obs=np.array(z_obs)
	start+=len(obs_lon_c)
	end+=len(obs_lon_c)

#reshape model array to plot on map
start = 0
end = len(model_lon_c)

for i in range(len(model_lat_c)):
        new_list = model_values[start:end]
        new_list = np.array(new_list)
        try:
                z_model =np.vstack((z_model,new_list))
        except:
                z_model = [new_list]
                z_model = np.array(z_model)
        start+=len(model_lon_c)
        end+=len(model_lon_c)

#cut off model lats above and below 60
z_model = z_model[16:75,:]

#take differences
z = z_model - z_obs 

#if standard dev of gridbox > certain num , mask it by making it black
stdevs = np.load('../process/OMI_O3_trop_ave_2x2.5_stdev.npy')

lat = np.arange(59)
lon = np.arange(144)

for la in lat:
    for lo in lon:
        data = np.average(stdevs[:,int(la),int(lo)])
        if data > 5:
            try:
                mask_lats=np.vstack((mask_lats,[obs_lat_e[la],obs_lat_e[la+1]]))
                mask_lons=np.vstack((mask_lons,[obs_lon_e[lo],obs_lon_e[lo+1]]))
            except:
                mask_lats = np.array([obs_lat_e[la],obs_lat_e[la+1]])
                mask_lons = np.array([obs_lon_e[lo],obs_lon_e[lo+1]])



z_stdev = np.zeros((len(mask_lats),1,1))

#setup basemap projection
m = Basemap(projection='cyl',llcrnrlat=-59,urcrnrlat=59,\
                   llcrnrlon=obs_lon_e[0],\
                    urcrnrlon=obs_lon_e[-1],\
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


#mask_lons, mask_lats = np.meshgrid(mask_lons, mask_lats)

#print mask_lons
# convert the xs and ys to map coordinates
#xs, ys = m(mask_lons, mask_lats)


#print xs,ys

#plot model gridboxes
if type == 'mag':
    poly = m.pcolor(obs_lon_e, obs_lat_e, z,vmin=-15, vmax=15,cmap = plt.cm.coolwarm)
    #for i in range(len(mask_lons)):
        #m.pcolor(mask_lons[i], mask_lats[i], z_stdev[i],cmap= plt.cm.gray)
else:
    poly = m.pcolor(obs_lon_e, obs_lat_e, z, vmin=np.min(-6), vmax=np.max(6), cmap=plt.cm.coolwarm)
    #for i in range(len(mask_lons)):
        #m.pcolor(mask_lons[i], mask_lats[i], z_stdev[i],cmap= plt.cm.gray)

if full_image == 'N':
    cb = plt.colorbar(poly, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')
    cb.set_label('%s'%(label), fontsize = 16)
    plt.xlabel('Longitude',fontsize = 20)
    plt.ylabel('Latitude',fontsize = 20)
    cb.ax.tick_params(labelsize=16)
    plt.title('Tropospheric O3 OMI 2x2.5 %s - GC 2x2.5 %s '%(title,title), fontsize = 18)
    plt.show()
else:
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    plt.savefig(f_name,bbox_inches=extent)



