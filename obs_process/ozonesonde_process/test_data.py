import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

a = Dataset('O3_RADIOSONDES_WOUDC_2009_2011.nc')
b = a.groups['stn099']
time  = b.variables['time']
o3  = b.variables['o3']

print o3.shape

#o3 = o3[:,8]

#test = np.where(o3 < 0)
#print test
#print test
#o3[test] = np.NaN

print o3.shape
o3 = np.nanmean(o3,axis=0)

print o3.shape
#test = o3 >= 0
#o3 = o3[test]
#time = time[test]

plt.plot(o3,range(47))
plt.show()
