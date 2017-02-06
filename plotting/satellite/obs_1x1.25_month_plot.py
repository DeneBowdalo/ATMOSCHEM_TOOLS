import numpy as np
import modules
import glob
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
import matplotlib.pyplot as plt


f = np.load('../process/OMI_O3_trop_ave.npy')

year = raw_input('Year\n')
month = raw_input('\nMonth No\n')

year_dict = {'2004':0,'2005':12,'2006':24,'2007':36,'2008':48,'2009':60,'2010':72,'2011':84,'2012':96,'2013':108}
month_dict = {'1':0,'2':1,'3':2,'4':3,'5':4,'6':5,'7':6,'8':7,'9':8,'10':9,'11':10,'12':11}

start_offset = 9

year = year_dict[year]
month = month_dict[month]


index = (year+month) - start_offset

print index
data = f[index,:,:]
print data[0,0]


values = []

print data.shape

for i in range(len(data[:,0])):
    for j in range(len(data[0,:])):
        values.append(data[i,j])
 

#lat lon edges for 1x2.25 grid
lat_e = np.arange(-60.,60.5,1)
lon_e = np.arange(-180,181,1.25)
lat_c = np.arange(-59.5,60,1)
lon_c = np.arange(-179.375,180,1.25)
 
#get size of grid
grid_dim_c = len(lat_c)*len(lon_c)
  
#reshape array to plot on map
start = 0
end = len(lon_c)
 

for i in range(len(lat_c)):
    new_list = values[start:end]
    #print new_list
    new_list = np.array(new_list)
    #print new_list
    try:
        z =np.vstack((z,new_list))
    except:
        z = [new_list]
        z=np.array(z)
    start+=len(lon_c)
    end+=len(lon_c)
 
print z

#set up plot
fig =plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')

#setup basemap projection
m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=60,\
                    llcrnrlon=lon_e[0],\
                     urcrnrlon=lon_e[-1],\
                    resolution='c')
 
m.drawcoastlines()
m.drawmapboundary()
parallels = np.arange(-60,61,15)
meridians = np.arange(-180,151,30)
plt.xticks(meridians)
plt.yticks(parallels)
m.drawparallels(parallels)
m.drawmeridians(meridians)

 
#plot gridboxes
poly = m.pcolor(lon_e, lat_e, z ,vmin = 0, vmax = 90)
plt.colorbar(poly)

plt.show()
