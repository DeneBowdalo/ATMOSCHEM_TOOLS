from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
import modules
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

lat_e = np.arange(-89.5,90.,1.)
lat_e = np.insert(lat_e,0,-90.0)
lat_e = np.append(lat_e,90.0)
lon_e = np.arange(-180.5,180.,1.)

#lon_e = np.arange(-182.5,177.5+2.5,5)
#lat_e = np.arange(-92,92+2,4.)
#lat_e[0]=-90.
#lat_e[-1]=90.

#lon_e = np.arange(-180.125,179.875+0.125,0.25)
#lat_e = np.arange(-90.125,90.125+0.125,0.25)
#lat_e[0]=-90.
#lat_e[-1]=90.

#lon_e = np.arange(-180.05,179.95+0.05,0.1)
#lat_e = np.arange(-90.05,90.05+0.05,0.1)
#lat_e[0]=-90.
#lat_e[-1]=90.


print lon_e
print lat_e

a = Dataset('NOX_SCALING_1x1v3.nc')
b = a.variables['scale'][6,:,:]

m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c')                               


m.pcolor(lon_e,lat_e,b,linewidth=0.5,cmap=plt.cm.gist_earth,rasterized=True) 
m.drawcoastlines()

plt.colorbar()

plt.show()
