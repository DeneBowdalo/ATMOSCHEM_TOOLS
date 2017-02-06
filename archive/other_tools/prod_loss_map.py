import glob
import numpy as np
from bpch import bpch
from matplotlib.colors import LogNorm

files = glob.glob('2005*')

files = sorted(files)

print files

o3_prod = np.empty((len(files),91,144))
o3_loss = np.empty((len(files),91,144))

no2_prod = np.empty((len(files),91,144))
no2_loss = np.empty((len(files),91,144))

ox_prod = np.empty((len(files),91,144))
ox_loss = np.empty((len(files),91,144))

o3_ave = np.empty((len(files),91,144))
no2_ave = np.empty((len(files),91,144))
no_ave = np.empty((len(files),91,144))
isop_ave = np.empty((len(files),91,144))

up_flx = np.empty((len(files),91,144))



for i in range(len(files)):
    data = bpch(files[i],tracerinfo='valid_tracerinfo.dat')
    group_pl = data.groups['PORL-L=$']
    group_ave = data.groups['IJ-AVG-$']
    group_up = data.groups['UP-FLX-$']

    po3 = group_pl.variables['PO3']
    print po3.units
    lo3 = group_pl.variables['LO3']
    pno2 = group_pl.variables['PNO2']
    lno2 = group_pl.variables['LNO2']
    pox = group_pl.variables['POX']
    lox = group_pl.variables['LOX']
    o3_a = group_ave.variables['O3']
    no_a = group_ave.variables['NOx']
    #no2_a = group_ave.variables['NO2']
    isop_a = group_ave.variables['ISOP']
    up_f = group_up.variables['Ox']

    
    po3 = po3[:,0,:,:]
    lo3 = lo3[:,0,:,:]
    pno2 = pno2[:,0,:,:]
    lno2 = lno2[:,0,:,:]
    pox = pox[:,0,:,:]
    lox = lox[:,0,:,:]
    o3_a = o3_a[:,0,:,:]
    no_a = no_a[:,0,:,:]
    isop_a = isop_a[:,0,:,:]
    up_f = up_f[:,0,:,:]
                                                                                                                                                                             
    o3_prod[i,:,:] = po3
    o3_loss[i,:,:] = lo3
    no2_prod[i,:,:] = pno2
    no2_loss[i,:,:] = lno2
    ox_prod[i,:,:] = pox
    ox_loss[i,:,:] = lox
    o3_ave[i,:,:] = o3_a
    no_ave[i,:,:] = no_a
    #no2_ave[i,:,:] = no2_a
    isop_ave[i,:,:] = isop_a    
    up_flx[i,:,:] = up_f
    

o3_pl = o3_prod-o3_loss
no2_pl = no2_prod-no2_loss


o3_pl = np.average(o3_pl,axis=0)
o3_ave = np.average(o3_ave,axis=0)
no_ave = np.average(no_ave,axis=0)
#no2_ave = np.average(no2_ave,axis=0)
isop_ave = np.average(isop_ave,axis=0)



#layer = np.copy(no_ave/isop_ave)
layer = no_ave

#import basemap module
from mpl_toolkits.basemap import Basemap 
import numpy as np

#set up plot                                                                                                                                                                 
import matplotlib.pyplot as plt

#setup plotting window
fig=plt.figure(figsize=(20,12))
ax = fig.add_axes([.05, .1, .9, .8])
fig.patch.set_facecolor('white')

#setup horizontal map
# llcrnrlat is the southern edge of the map, urcrnrlat is the northern edge of the map.
# llcrnrlon is the western edge of the map, urcrnrlon is the eastern edge of the map.

m = Basemap(projection='cyl', llcrnrlat=-90,  urcrnrlat=90, llcrnrlon=-181.25, urcrnrlon=178.75, resolution='c')

#Set Latitudes of grid
lat = np.arange(-89,90,2)
lat = np.insert(lat,0,-90)
lat = np.append(lat,90)

# Set longitudes of grid
lon = np.arange(-181.25,179,2.5)

#Draw map detail
m.drawcoastlines()
m.drawmapboundary()
parallels = np.arange(-90,91,15)
meridians = np.arange(-180,151,30)                                                                                                                                          
plt.xticks(meridians)
plt.yticks(parallels)
m.drawparallels(parallels)
m.drawmeridians(meridians)

#m.parallels()
#m.meridians()

#Plot data on the global map in a 4x5 grid.
poly = m.pcolor(lon, lat, layer)

#Plot a colourbar for the variable plotted.
cb = plt.colorbar(poly, ax = m.ax,shrink=0.8)

plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)
#plt.title('Average Surface O3, Jan 2005',fontsize=20)
plt.text(190,95,'Ozone (ppbv)',fontsize=20)

plt.show()
