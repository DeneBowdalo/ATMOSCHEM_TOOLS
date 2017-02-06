import matplotlib.pyplot as plt
import modules
from netCDF4 import Dataset
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
import numpy as np
import os
#----------------------------------------------------
#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres = modules.get_obs_info(present_dir)
model_fname,species,start_year,end_year = modules.get_model_info(present_dir)

#------------------------------------------------------------------
#load in spectra lsp data
root_grp_obs = Dataset('../obs_%s_%s/obs_spectra.nc')
root_grp_mod = Dataset('model_spectra.nc')

valid_refs_dict = root_grp_obs.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

site_ref = raw_input('\nChoose site from list. Sites with full set of yearly files between %s & %s are:\n%s\n'%(start_year,end_year,'   '.join(i for i in valid_refs)))

site_group_obs = root_grp_obs.groups[site_ref]
site_group_mod = root_grp_mod.groups[site_ref]

obs_amp = site_group_obs.variables['amplitude'][:]
model_amp = site_group_mod.variables['amplitude'][:]
obs_period = site_group_obs.variables['period'][:]
model_period = site_group_mod.variables['period'][:]

#--------------------------------------------------------------------

#set plotting area & background to white
fig=plt.figure(figsize=(25,12))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)

plt.loglog(obs_period,obs_amp, color='black',alpha = 1, label = 'Observations')
plt.loglog(model_period,model_amp, color='red', alpha= 0.6, label = '%s %s %s %s'%(model_name,version,grid,met))

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
plt.title('Lomb-Scargle Periodogram of Surface %s at %s, for Observations & %s %s'%(species,site_ref,model_name,model_range),fontsize=22,y=1.02)
#ax.tick_params(axis='both', which='major', labelsize=15)
for tick in ax.get_xaxis().get_major_ticks():
    tick.set_pad(11.)
plt.tight_layout()
plt.show()
