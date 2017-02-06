from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from collections import OrderedDict 
import operator
import os
import glob

#SET FILE RES, CAN BE HOURLY(H), DAILY(D), MONTHLY(M) OR ANNUAL(Y) 
file_res = 'H'

#OUTPUT RESOLUTION CAN BE HOURLY (H), HOURLY & DAILY (HD), HOURLY, DAILY AND MONTHLY (HDM)
output_res = 'HDM'

#OUTPUT SET CAN BE PERIODIC (P), TREND1 (T), TREND2 (X), NO RESTRAINTS (N)
output_set = 'T'

species_list = ['O3','NO','NO2-MOLYBDENUM','NO2-PHOTOLYTIC','CO']

color_dict = {'WMO GAW':'Crimson','EMEP':'pink','CASTNET':'blue','NAPS':'purple','EPA AQS':'orange','CAPMON':'yellow','EANET':'green','AirBase':'black','SEARCH':'grey'}
zorders = {'WMO GAW':4,'EMEP':2,'CASTNET':2,'NAPS':3,'EPA AQS':1,'CAPMON':5,'EANET':6,'AirBase':1,'SEARCH':6}


stats_data = Dataset('STATS_%s.nc'%(output_res+output_set))
c1 = stats_data.variables['invalid_nometa_count'][0]
c2 = stats_data.variables['invalid_anyvaliddata_count'][0]
c3 = stats_data.variables['invalid_nokeymeta_count'][0]
c4 = stats_data.variables['invalid_resolution_count'][0]
c5 = stats_data.variables['invalid_badmeasurementmethod_count'][0]
c6 = stats_data.variables['invalid_duplicatesites_count'][0]
c7 = stats_data.variables['invalid_rawclass_count'][0]
c8 = stats_data.variables['invalid_anthromeclass_count'][0]
c9 = stats_data.variables['invalid_altitude_count'][0]
c10 = stats_data.variables['invalid_representativeness_count'][0]
c11 = stats_data.variables['invalid_extreme_count'][0]
c12 = stats_data.variables['invalid_partialyear_count'][0]

e1 = stats_data.variables['exit_nometa_refs'][:]
e2 = stats_data.variables['exit_nometa_lats'][:]
e3 = stats_data.variables['exit_nometa_lons'][:]
e4 = stats_data.variables['exit_nometa_pg'][:]
e5 = stats_data.variables['exit_anyvaliddata_refs'][:]
e6 = stats_data.variables['exit_anyvaliddata_lats'][:]
e7 = stats_data.variables['exit_anyvaliddata_lons'][:]
e8 = stats_data.variables['exit_anyvaliddata_pg'][:]
e9 = stats_data.variables['exit_nokeymeta_refs'][:]
e10 = stats_data.variables['exit_nokeymeta_lats'][:]
e11 = stats_data.variables['exit_nokeymeta_lons'][:]
e12 = stats_data.variables['exit_nokeymeta_pg'][:]
e13 = stats_data.variables['exit_resolution_refs'][:]
e14 = stats_data.variables['exit_resolution_lats'][:]
e15 = stats_data.variables['exit_resolution_lons'][:]
e16 = stats_data.variables['exit_resolution_pg'][:]
e17 = stats_data.variables['exit_badmeasurementmethod_refs'][:]
e18 = stats_data.variables['exit_badmeasurementmethod_lats'][:]
e19 = stats_data.variables['exit_badmeasurementmethod_lons'][:]
e20 = stats_data.variables['exit_badmeasurementmethod_pg'][:]
e21 = stats_data.variables['exit_duplicatesites_refs'][:]
e22 = stats_data.variables['exit_duplicatesites_lats'][:]
e23 = stats_data.variables['exit_duplicatesites_lons'][:]
e24 = stats_data.variables['exit_duplicatesites_pg'][:]
e25 = stats_data.variables['exit_rawclass_refs'][:]
e26 = stats_data.variables['exit_rawclass_lats'][:]
e27 = stats_data.variables['exit_rawclass_lons'][:]
e28 = stats_data.variables['exit_rawclass_pg'][:]
e29 = stats_data.variables['exit_anthromeclass_refs'][:]
e30 = stats_data.variables['exit_anthromeclass_lats'][:]
e31 = stats_data.variables['exit_anthromeclass_lons'][:]
e32 = stats_data.variables['exit_anthromeclass_pg'][:]
e33 = stats_data.variables['exit_altitude_refs'][:]
e34 = stats_data.variables['exit_altitude_lats'][:]
e35 = stats_data.variables['exit_altitude_lons'][:]
e36 = stats_data.variables['exit_altitude_pg'][:]
e37 = stats_data.variables['exit_representativeness_refs'][:]
e38 = stats_data.variables['exit_representativeness_lats'][:]
e39 = stats_data.variables['exit_representativeness_lons'][:]
e40 = stats_data.variables['exit_representativeness_pg'][:]
e41 = stats_data.variables['exit_extreme_refs'][:]
e42 = stats_data.variables['exit_extreme_lats'][:]
e43 = stats_data.variables['exit_extreme_lons'][:]
e44 = stats_data.variables['exit_extreme_pg'][:]
e45 = stats_data.variables['exit_partialyear_refs'][:]
e46 = stats_data.variables['exit_partialyear_lats'][:]
e47 = stats_data.variables['exit_partialyear_lons'][:]
e48 = stats_data.variables['exit_partialyear_pg'][:]

#set up plot
fig =plt.figure(figsize=(16,8))
fig.patch.set_facecolor('white')
ax1 = plt.subplot2grid((4,3), (0,0))
ax2 = plt.subplot2grid((4,3), (0,1))
ax3 = plt.subplot2grid((4,3), (0,2))
ax4 = plt.subplot2grid((4,3), (1,0))
ax5 = plt.subplot2grid((4,3), (1,1))
ax6 = plt.subplot2grid((4,3), (1,2))
ax7 = plt.subplot2grid((4,3), (2,0))
ax8 = plt.subplot2grid((4,3), (2,1))
ax9 = plt.subplot2grid((4,3), (2,2))
ax10 = plt.subplot2grid((4,3), (3,0))
#ax11 = plt.subplot2grid((4,3), (3,1))
#ax12 = plt.subplot2grid((4,3), (3,2))

ax_list = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10]
exit_refs = [e5,e13,e17,e21,e25,e29,e33,e37,e41,e45]
exit_lats = [e6,e14,e18,e22,e26,e30,e34,e38,e42,e46]
exit_lons = [e7,e15,e19,e23,e27,e31,e35,e39,e43,e47]
exit_pg = [e8,e16,e20,e24,e28,e32,e36,e40,e44,e48]

titles = ['No Valid Data','Resolution too coarse','Bad Measurement/Sampling Method','Remove Duplicates','Urban Raw','Urban Anthrome','Altitude > 1500m','All Un-Representative of day','Average too high - Urban','All Partialyears']

ax_count = 0
for ax in ax_list:
    print ax_count

    counter = 0

    m = Basemap(projection='moll',lon_0 = 0,resolution='c',ax = ax)
 
    m.fillcontinents(color='#cc9966',lake_color='#99ffff',alpha=0.7,zorder=0.5)
    m.drawmapboundary(fill_color='#99ffff')
    
    valid_refs = filter(None,exit_refs[ax_count])
    valid_lats = filter(None,exit_lats[ax_count])
    valid_lons = filter(None,exit_lons[ax_count])
    valid_pg = filter(None,exit_pg[ax_count])
    
    print len(valid_refs),len(valid_lats),len(valid_lons),len(valid_pg)
    
    for i in range(len(valid_refs)):
        obs_ref = valid_refs[i]
        print obs_ref
        obs_lat = np.float64(valid_lats[i])
        obs_lon = np.float64(valid_lons[i])
        obs_group = valid_pg[i]

        x, y = np.meshgrid(*m(obs_lon, obs_lat))
        m.scatter(x, y, color=color_dict[obs_group],s=20,label=obs_group,zorder=zorders[obs_group],alpha=1.)

        counter+=1
    
    ax.set_title('%s:%s'%(titles[ax_count],len(valid_refs)))
    
    ax_count+=1

#handles, labels = plt.gca().get_legend_handles_labels()
#hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
#handles, labels = zip(*hl)
#by_label = OrderedDict(zip(labels, handles))
#leg = ax.legend(by_label.values(), by_label.keys(), loc = 'upper left', prop={'size':23},markerscale=2)
#leg.get_frame().set_facecolor('white')

#plt.title('Observational %s Sites %s'%(species,date_range),fontsize = 20)

plt.show()
