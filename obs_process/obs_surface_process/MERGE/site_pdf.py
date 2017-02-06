from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from collections import OrderedDict 
import operator
import os
import glob
import seaborn as sns

species = os.getcwd().split('/')[-2]

#SET FILE RES, CAN BE HOURLY(H), DAILY(D), MONTHLY(M) OR ANNUAL(Y) 
file_res = 'H'

#OUTPUT RESOLUTION CAN BE HOURLY (H), HOURLY & DAILY (HD), HOURLY, DAILY AND MONTHLY (HDM)
output_res = 'HDM'

output_set = 'S'

#read in obs data
root_grp = Dataset('GLOBAL_SURFACE_%s_1970_2015_%s_%s.nc'%(species,file_res,output_res+output_set))
valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

#set up plot
fig=plt.figure(figsize=(24,14))
ax = fig.add_axes([.05, .1, .9, .8])
fig.patch.set_facecolor('white')

                                                                                                                                                                                                     
color_dict = {'WMO GAW':'Crimson','EMEP':'pink','CASTNET':'blue','NAPS':'purple','EPA AQS':'orange','CAPMON':'yellow','EANET':'green','AirBase':'black','SEARCH':'grey'}
zorders = {'WMO GAW':4,'EMEP':2,'CASTNET':2,'NAPS':3,'EPA AQS':1,'CAPMON':5,'EANET':6,'AirBase':1,'SEARCH':6}

all_data = []
all_data_ave = []
all_data_max = []
print len(valid_refs)
for site_ref in valid_refs:
    print site_ref
    site_group = root_grp.groups[site_ref]
   
    obs_data = site_group.variables[species.lower()][:]
    valid = obs_data >= 0
    all_data.append(obs_data[valid])
    all_data_ave.append(np.average(obs_data[valid]))    
    all_data_max.append(np.max(obs_data[valid]))
    #obs_lat = site_group.latitude
    #obs_lon = site_group.longitude
    #obs_group = site_group.process_group
    #obs_res = site_group.native_resolution
all_data = np.array([item for sublist in all_data for item in sublist])
all_data_ave = np.array(all_data_ave)

print np.max(all_data)

#weights = np.ones_like(all_data)/float(len(all_data)) 
weights = np.ones_like(all_data_ave)/float(len(all_data_ave))

#ax.hist(all_data,1000,normed=0,color='green',alpha=0.8,weights=weights)
ax.hist(all_data_ave,100,normed=0,color='green',alpha=0.8,weights=weights)

plt.show()
