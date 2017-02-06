from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from collections import OrderedDict 
import operator
import os
import glob

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-2]
print '\nSpecies is %s\n'%(species)

valid_files = glob.glob('../process/*.nc')

#read in obs data
root_grp = Dataset(valid_files[0])
valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

#set up plot
fig=plt.figure(figsize=(24,14))
ax = fig.add_axes([.05, .1, .9, .8])
fig.patch.set_facecolor('white')

m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                      llcrnrlon=-180,\
                      urcrnrlon=180,\
                      resolution='c')
     
m.fillcontinents(color='#cc9966',lake_color='#99ffff')
m.drawmapboundary(fill_color='#99ffff')
                                                                                                                                                                                                     
color_dict = {'SHADOZ':'red','WOUDC':'black'}
size_dict = {'SHADOZ':100,'WOUDC':175}


for site_ref in valid_refs:
    site_group = root_grp.groups[site_ref]
   
    obs_lat = site_group.latitude
    obs_lon = site_group.longitude
    obs_group = site_group.process_group

    print obs_group
    print obs_lat,obs_lon

    x, y = np.meshgrid(*m(obs_lon, obs_lat))
    m.scatter(x, y, color=color_dict[obs_group],s=size_dict[obs_group],label=obs_group,zorder=20)

handles, labels = plt.gca().get_legend_handles_labels()
hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
handles, labels = zip(*hl)
by_label = OrderedDict(zip(labels, handles))
leg = ax.legend(by_label.values(), by_label.keys(), loc = 3, prop={'size':20})
leg.get_frame().set_facecolor('grey')

plt.show()
