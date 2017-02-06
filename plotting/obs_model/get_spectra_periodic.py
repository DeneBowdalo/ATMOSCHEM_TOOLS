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
paths = present_dir.split("/")
species = paths[-3]
model_range = paths[-2]
print '\nSpecies is %s\n'%(species)

model_details = paths[-1]
print model_details

model_split = model_details.split('_SFC')
model_name = model_split[0]
model_other = model_split[1]
model_other_split = model_other.split('_')
#model_other_split = model_other_split.remove('')
years = paths[-2]
start_year = years[:4]
end_year = years[5:9] 

if len(model_other_split) == 1:
    version = ''
    grid = ''
    met = ''
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'.nc'

elif len(model_other_split) == 2:
    version = model_other_split[1] 
    grid = ''
    met = ''
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'.nc'

elif len(model_other_split) == 3:
    version = model_other_split[1] 
    grid = model_other_split[2]
    met = ''
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'_'+grid+'.nc'

elif len(model_other_split) == 4:
    version = model_other_split[1] 
    grid = model_other_split[2]
    met = model_other_split[3]
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'_'+grid+'_'+met+'.nc'

#------------------------------------------------------------------
#load in lsp periodic data
root_grp_obs = Dataset('../obs/obs_sig_periods.nc')
root_grp_mod = Dataset('model_sig_periods.nc')

valid_refs_dict = root_grp_obs.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

site_ref = raw_input('\nChoose site from list. Sites with full set of yearly files between %s & %s are:\n%s\n'%(start_year,end_year,'   '.join(i for i in valid_refs)))

site_group_obs = root_grp_obs.groups[site_ref]
site_group_mod = root_grp_mod.groups[site_ref]

obs_daily_amp = site_group_obs.daily_amplitude
obs_ha_amp = site_group_obs.half_annual_amplitude
obs_annual_amp = site_group_obs.annual_amplitude
obs_daily_ph = site_group_obs.daily_phase
obs_ha_ph = site_group_obs.half_annual_phase
obs_annual_ph = site_group_obs.annual_phase
obs_average = site_group_obs.average

mod_daily_amp = site_group_mod.daily_amplitude
mod_ha_amp = site_group_mod.half_annual_amplitude
mod_annual_amp = site_group_mod.annual_amplitude
mod_daily_ph = site_group_mod.daily_phase
mod_ha_ph = site_group_mod.half_annual_phase
mod_annual_ph = site_group_mod.annual_phase
mod_average = site_group_mod.average


print '\nSite is %s\n'%(site_ref) 

print 'Obs. Diurnal Amplitude = ', obs_daily_amp , 'ppb'
print 'Obs. Half-Annual Amplitude = ', obs_ha_amp, 'ppb'
print 'Obs. Annual Amplitude = ', obs_annual_amp, 'ppb'
print 'Obs. Diurnal Phase = ', obs_daily_ph, 'ppb'
print 'Obs. Half-Annual Phase = ', obs_ha_ph, 'ppb'
print 'Obs. Annual Phase = ', obs_annual_ph, 'ppb'
print 'Obs. Average = ', obs_average, 'ppb', '\n'

print 'Model Diurnal Amplitude = ', mod_daily_amp
print 'Model Half-Annual Amplitude = ', mod_ha_amp
print 'Model Annual Amplitude = ', mod_annual_amp
print 'Model Diurnal Phase = ', mod_daily_ph
print 'Model Half-Annual Phase = ', mod_ha_ph
print 'Model Annual Phase = ', mod_annual_ph, '\n'
print 'Model Average = ', mod_average, 'ppb'



