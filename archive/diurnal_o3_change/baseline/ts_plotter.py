import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import modules

#Read in Magnitude or phase
typ = raw_input('magnitude or phase?\n')

if typ == 'magnitude':
    f_name = 'magnitudes'
    
else:
    f_name = 'phases'

#read obs data
obs_refs,obs_locs,obs_lats,obs_lons,obs_alt,obs_number = np.loadtxt('GAW_site_indices',dtype='S20,S20,f5,f5,f5,i5',unpack=True)
#get site name from ref

#choose site from ref
#site_name = 'lgb'
site_name = raw_input('What Site?\n%s\n'%('   '.join(str(i) for i in obs_locs)))
all_sims = raw_input('Plot all sims? Y or N?\n')

o_test = obs_locs == site_name
site_ref = obs_refs[o_test]
site_ref = site_ref[0]

#get obs_lon of site
obs_lon = obs_lons[o_test]
#obs_lon = obs_lon[0]

#lat lon edges for 4x5 grid 
lat_e = np.arange(-88,90,4)
lat_e = np.insert(lat_e,0,-90)
lat_e = np.append(lat_e,90)

lon_e = np.arange(-182.5,178,5)

lat_c = np.arange(-86,90,4)
lat_c = np.insert(lat_c,0,-89)
lat_c = np.append(lat_c,89)

lon_c = np.arange(-180,178,5)

def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

#read obs 

obs_files = glob.glob('obs_%s/GAW*'%(f_name))
obs_files = sorted_nicely(obs_files)

#read baseline 4x5 model

base_files = glob.glob('model_%s/GAW*'%(f_name))
base_files = sorted_nicely(base_files)
 
#read convectoff 4x5 model

convectoff_files = glob.glob('/home/db876/diurnal_o3_change/convectoff/model_%s/GAW*'%(f_name))
convectoff_files = sorted_nicely(convectoff_files)


#read drydepoff 4x5 model

drydepoff_files = glob.glob('/home/db876/diurnal_o3_change/drydepoff/model_%s/GAW*'%(f_name))
drydepoff_files = sorted_nicely(drydepoff_files)

#read emissionsoff 4x5 model

emissionsoff_files = glob.glob('/home/db876/diurnal_o3_change/emissionsoff/model_%s/GAW*'%(f_name))
emissionsoff_files = sorted_nicely(emissionsoff_files)

#read transportoff 4x5 model

transportoff_files = glob.glob('/home/db876/diurnal_o3_change/transportoff/model_%s/GAW*'%(f_name))
transportoff_files = sorted_nicely(transportoff_files)


#set arrays to append to

obs_vals = []
base_vals = []
convectoff_vals = []
drydepoff_vals = []
emissionsoff_vals = []
transportoff_vals = []


for f in range(len(obs_files)):
    obs_read = np.load(obs_files[f])
    base_read = np.load(base_files[f])
    convectoff_read = np.load(convectoff_files[f])
    drydepoff_read = np.load(drydepoff_files[f])
    emissionsoff_read = np.load(emissionsoff_files[f])
    transportoff_read = np.load(transportoff_files[f])
    
    #cut to just 1 site
    site_test = obs_read[1,:] == site_ref
    
    obs_read = obs_read[0,site_test]
    base_read = base_read[site_test]
    convectoff_read = convectoff_read[site_test]
    drydepoff_read = drydepoff_read[site_test]
    emissionsoff_read = emissionsoff_read[site_test]
    transportoff_read = transportoff_read[site_test]
    
    
    obs_read = np.float64(obs_read)
   
    obs_vals = np.append(obs_vals,obs_read)
    base_vals = np.append(base_vals,base_read)
    convectoff_vals = np.append(convectoff_vals,convectoff_read)
    drydepoff_vals = np.append(drydepoff_vals,drydepoff_read)
    emissionsoff_vals = np.append(emissionsoff_vals,emissionsoff_read)
    transportoff_vals = np.append(transportoff_vals,transportoff_read)
    
fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)
time = np.arange(0,72)
plt.plot(time, obs_vals, label = 'Observations')
plt.plot(time, base_vals, label = 'GEOS-Chem 4x5')
if all_sims == 'Y':
    plt.plot(time, convectoff_vals, label = 'GEOS-Chem 4x5 Convect off')
    plt.plot(time, drydepoff_vals, label = 'GEOS-Chem 4x5 Drydep off')
    plt.plot(time, emissionsoff_vals, label = 'GEOS-Chem 4x5 Emissions off')
    plt.plot(time, transportoff_vals, label = 'GEOS-Chem 4x5 Transport off')
xticks = [0,12,24,36,48,60,72]
xticklabels = ['2006','2007','2008','2009','2010','2011','2012']
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(xticklabels)
leg = plt.legend(loc=0, prop={'size':21})

plt.xlabel('Time (Year)',fontsize = 21)
if typ == 'magnitude':
    plt.ylabel('Magnitude (ppbv)', fontsize=21)
    plt.title('Monthly Daily Magnitude - %s'%(site_name), fontsize = 21)
else:
    plt.ylabel('Phase (Hours)', fontsize=21)
    plt.title('Monthly Daily Phase - %s'%(site_name), fontsize = 21)
plt.xlim(0,72)

plt.show()    


