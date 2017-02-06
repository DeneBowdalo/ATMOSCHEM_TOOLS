import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import seaborn as sns
import os

species = os.getcwd().split('/')[-2]

fig =plt.figure(figsize=(11,9.5))
fig.patch.set_facecolor('white')

ref = 'ee0009r'
root_grp = Dataset('GLOBAL_SURFACE_%s_1970_2015_H_HDMN.nc'%(species))
time = root_grp.variables['time'][:]                                                                                                                     
data = root_grp.groups[ref].variables[species.lower()][:]
data[data <= 0] = np.NaN
plt.plot(range(len(data)),data,marker='o',linestyle='None')
plt.show()


