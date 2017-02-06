from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from collections import OrderedDict 
import operator
import os
import glob
import modules

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-2]
print '\nSpecies is %s\n'%(species)


start_year = 2009
end_year = 2011


timeres = 'H'

#read in obs data
root_grp = Dataset('../process/GLOBAL_SURFACE_%s_%s_%s_%s_PERIODIC.nc'%(species,start_year,end_year,timeres))
valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

#set up plot
fig=plt.figure(figsize=(24,14))
ax = fig.add_axes([.05, .1, .9, .8])
fig.patch.set_facecolor('white')

m = Basemap(projection='cyl',llcrnrlat=22,urcrnrlat=75,\
                      llcrnrlon=-143,\
                      urcrnrlon=-51,\
                      resolution='c')

#m = Basemap(projection='moll',lon_0 = 0,resolution='c')
     
m.fillcontinents(color='#cc9966',lake_color='#99ffff',alpha=0.7,zorder=0.5)
m.drawmapboundary(fill_color='#99ffff')
                                                                                                                                                                                                     
#color_dict = {'GAW':'Crimson','EMEP':'pink','CASTNET':'blue','NAPS':'purple','AQS':'orange','CAPMON':'yellow','EANET':'green','AIRBASE':'black'}
#zorders = {'GAW':4,'EMEP':2,'CASTNET':2,'NAPS':3,'AQS':1,'CAPMON':5,'EANET':6,'AIRBASE':1}

areas = ['SW_NA','NW_NA','NE_NA','CE_NA','SE_NA','C_NA']

color_dict = {'SW_NA':'Crimson','NW_NA':'pink','NE_NA':'darkblue','CE_NA':'orange','SE_NA':'green','C_NA':'yellow','S_EU':'black','SW_EU':'magenta','CW_EU':'brown','CS_EU':'orangered','C_EU':'lightskyblue','E_EU':'limegreen','NW_EU':'purple','N_EU':'darkgreen'}

tags = modules.get_tags(valid_refs)
area_boundaries,area_tags,area_labels = modules.area_dicts()

for i in range(len(valid_refs)):
    site_ref = valid_refs[i]
    
    site_group = root_grp.groups[site_ref]
   
    obs_lat = site_group.latitude
    obs_lon = site_group.longitude
    
    tag = tags[i]
    
    obs_area = modules.get_area(area_boundaries,tag,obs_lat,obs_lon,site_ref)

    x, y = np.meshgrid(*m(obs_lon, obs_lat))
    
    try:
        m.scatter(x, y, color=color_dict[obs_area],s=45,alpha=1.)
    except:
        pass


#handles, labels = plt.gca().get_legend_handles_labels()
#hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
#handles, labels = zip(*hl)
#by_label = OrderedDict(zip(labels, handles))
#leg = ax.legend(by_label.values(), by_label.keys(), loc = 'upper left', prop={'size':23},markerscale=2)
#leg.get_frame().set_facecolor('white')

#plt.title('Observational %s Sites %s'%(species,date_range),fontsize = 20)

plt.show()
