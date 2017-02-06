from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from collections import OrderedDict 
import operator
import os
import glob

species = os.getcwd().split('/')[-2]

#SET FILE RES, CAN BE HOURLY(H), DAILY(D), MONTHLY(M) OR ANNUAL(Y) 
file_res = 'H'

#OUTPUT RESOLUTION CAN BE HOURLY (H), HOURLY & DAILY (HD), HOURLY, DAILY AND MONTHLY (HDM)
output_res = 'HDM'

#OUTPUT SET CAN BE PERIODIC (P), TREND1 (T), TREND2 (X), NO RESTRAINTS (N)
output_set = 'T'

type_plot = ''

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

#m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
#                      llcrnrlon=-180,\
#                      urcrnrlon=180,\
#                      resolution='c')

m = Basemap(projection='moll',lon_0 = 0,resolution='c')
     
m.fillcontinents(color='#cc9966',lake_color='#99ffff',alpha=0.7,zorder=0.5)
m.drawmapboundary(fill_color='#99ffff')
                                                                                                                                                                                                     
color_dict = {'WMO GAW':'Crimson','EMEP':'pink','CASTNET':'blue','NAPS':'purple','EPA AQS':'orange','CAPMON':'yellow','EANET':'green','AirBase':'black','SEARCH':'grey'}
zorders = {'WMO GAW':4,'EMEP':2,'CASTNET':2,'NAPS':3,'EPA AQS':1,'CAPMON':5,'EANET':6,'AirBase':1,'SEARCH':6}


for site_ref in valid_refs:
    site_group = root_grp.groups[site_ref]
   
    obs_lat = site_group.latitude
    obs_lon = site_group.longitude
    obs_group = site_group.process_group
    obs_res = site_group.native_resolution

    print obs_group
    print obs_lat,obs_lon

    if type_plot == 'res':
        if output_res == 'H':
            continue
        
        
    x, y = np.meshgrid(*m(obs_lon, obs_lat))
    m.scatter(x, y, color=color_dict[obs_group],s=35,label=obs_group,zorder=zorders[obs_group],alpha=1.)



handles, labels = plt.gca().get_legend_handles_labels()
hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
handles, labels = zip(*hl)
by_label = OrderedDict(zip(labels, handles))
leg = ax.legend(by_label.values(), by_label.keys(), loc = 'upper left', prop={'size':23},markerscale=2)
leg.get_frame().set_facecolor('white')

#plt.title('Observational %s Sites %s'%(species,date_range),fontsize = 20)

plt.show()
