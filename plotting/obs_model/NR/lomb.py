import numpy as np
import matplotlib.pyplot as plt
import glob
import lomb_phase
import modules
import datetime
from netCDF4 import Dataset
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
from scipy.stats import norm
import scipy.signal
import os
import numpy.fft

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-2]

#process model
#----------------------------------------
#read in model data

start_year = raw_input('Start Year?\n')
end_year = raw_input('\nEnd Year?\n')

#find valid_models
model_files = glob.glob('/work/home/db876/plotting_tools/model_files/*%s_%s_%s*'%(species,start_year,end_year))
model = []
vres = []
version = []
hres = []
met = []
timeres = []
additional = []

for f in model_files:
    f = f.replace("/work/home/db876/plotting_tools/model_files/", "")
    data_split = f.split('_')
    model = np.append(model,data_split[0])
    vres = np.append(vres,data_split[1])        
    version = np.append(version,data_split[5])        
    hres = np.append(hres,data_split[6])
    met = np.append(met,data_split[7])
    timeres = np.append(timeres,data_split[8])    
    additional = np.append(additional,data_split[9])

model_unique = np.unique(model)

if len(model_unique) == 1:
    model_chose = model_unique[0]
else:
    model_chose = raw_input('\nChoose Model\n%s\n'%('   '.join(i for i in model_unique))) 

inds = np.where((model == model_chose)) 

vres = vres[inds]
version = version[inds]
hres = hres[inds]
met = met[inds]
timeres = timeres[inds]
additional = additional[inds]

vres = np.unique(vres)
version = np.unique(version)
hres = np.unique(hres)
met = np.unique(met)
timeres = np.unique(timeres)
additional = np.unique(additional)

if len(vres) == 1:
    vres = vres[0]
else:
    vres = raw_input('\nChoose Vertical Resolution\n%s\n'%('   '.join(i for i in vres)))

if len(version) == 1:
    version = version[0]
else:
    version = raw_input('\nChoose Version\n%s\n'%('   '.join(i for i in version)))     

if len(hres) == 1:
    hres = hres[0]
else:
    hres = raw_input('\nChoose Horizontal Resolution\n%s\n'%('   '.join(i for i in hres))) 

if len(met) == 1:
    met = met[0]
else:
    met = raw_input('\nChoose Meteorology\n%s\n'%('   '.join(i for i in met))) 

if len(timeres) == 1:
    timeres = timeres[0]
else:
    timeres = raw_input('\nChoose Time Resolution\n%s\n'%('   '.join(i for i in timeres)))

if len(additional) == 1:
    additional = additional[0]
else:
    additional = raw_input('\nChoose Additional Detail\n%s\n'%('   '.join(i for i in additional)))

model_fname = '/work/home/db876/plotting_tools/model_files/'+model_chose+'_'+vres+'_'+species+'_'+start_year+'_'+end_year+'_'+version+'_'+hres+'_'+met+'_'+timeres+'_'+additional

print model_fname
root_grp = Dataset(model_fname)
model_var = root_grp.variables[species.lower()][:]
model_date = root_grp.variables['date'][:]
model_time = root_grp.variables['time'][:]
lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
grid_size = root_grp.variables['grid_size'][:]
grid_size = grid_size[0]

gridbox_count = len(lat_c)*len(lon_c)

#process obs
#----------------------------------------
#data range

#read in obs data
root_grp = Dataset('../process/GLOBAL_%s_%s_1970_2015_%s_ALL_NR.nc'%(vres,species,timeres))

valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

site_ref = raw_input('\nChoose site from list. Sites with full set of yearly files between %s & %s are:\n%s\n'%(start_year,end_year,'   '.join(i for i in valid_refs)))

#read in specific site data
site_group = root_grp.groups[site_ref]

#read in variables for site
obs_var = site_group.variables[species.lower()][:]
obs_date = site_group.variables['date'][:]
obs_time = site_group.variables['time'][:]
obs_lat = site_group.latitude
obs_lon = site_group.longitude
obs_alt = site_group.altitude
obs_group = site_group.process_group

obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
obs_var = obs_var[~np.isnan(obs_var_mask)]

#----------------------------------------
#find model data gridbox to compare with obs.

#get model gridbox for obs site
lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)

model_var = model_var[:,lat_n,lon_n]
model_var = model_var*1e9

model_var_mask = np.ma.masked_where(model_var<=0,model_var)
model_ave = np.ma.average(model_var_mask)
model_var = model_var[~np.isnan(model_var_mask)]

#--------------------------------------------
#get valid data and process time

obs_time = np.array(modules.date_process(obs_date,obs_time,start_year))
model_time = np.array(modules.date_process(model_date,model_time,start_year))

model_test = model_var >= 0
model_var = model_var[model_test]
model_time = model_time[model_test]

obs_test = obs_var >= 0
obs_var = obs_var[obs_test]
obs_time = obs_time[obs_test]

#--------------------------------------------
#take LSP's

#windowing?
wind_set = raw_input('Windowing? Y or N?\n')
if wind_set == 'Y':
    wind_set = True
if wind_set == 'N':
    wind_set = False

#take lomb for obs and model
ofac = int(raw_input('\nChoose oversampling factor. (Typically 4.)\n'))

periodic_periods = [1./10.,1./9.,1./8.,1./7.,1./6.,1./5.,1./4.,1./3.,1./2.,1.,365.25/4.,365.25/3.,365.25/2.,365.25]

samp_step = 1./24

#obs lomb
obs_periods,obs_mag,obs_ph,obs_fr,obs_fi,amp_corr = modules.take_lomb(obs_time,obs_var,ofac,samp_step,w=wind_set,kp=periodic_periods)

#model lomb
model_periods,model_mag,model_ph,model_fr,model_fi,amp_corr = modules.take_lomb(model_time,model_var,ofac,samp_step,w=wind_set,kp=periodic_periods)

obs_time = np.array(obs_time)
obs_var = np.array(obs_var)

#set plotting area & background to white
fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)

plt.loglog(obs_periods,obs_mag, color='black',alpha = 1, label = 'Observations')
plt.loglog(model_periods,model_mag, color='red', alpha= 0.6, label = '%s %s %s %s %s'%(model_chose,version,hres,met,timeres))

def form2(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.2f' % x	

def form5(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.7f' % x

xformatter = FuncFormatter(form2)
yformatter = FuncFormatter(form5)

plt.grid(True)
leg=plt.legend(loc=4, prop={'size':20})
leg.get_frame().set_alpha(0.4)
ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
plt.xlabel('Period (Days)',fontsize = 15)
plt.ylabel('Amplitude (ppb)',fontsize = 15)
plt.title('Lomb-Scargle Periodogram of Surface %s at %s, for Observations & %s %s:%s'%(species,site_ref,model_chose,start_year,end_year),fontsize=22,y=1.02)
#ax.tick_params(axis='both', which='major', labelsize=15)
for tick in ax.get_xaxis().get_major_ticks():
    tick.set_pad(11.)
plt.show()
