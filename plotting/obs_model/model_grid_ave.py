import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm
from cmath import *
import colorsys
from scipy import signal
from netCDF4 import Dataset
import modules
from bpch import bpch

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-2]
print '\nSpecies is %s\n'%(species)

#process model
#----------------------------------------
#read in model data

model_dict = {}
version_dict = {}
grid_dict = {}
met_dict = {}

#find valid_models
model_files = glob.glob('/work/home/db876/plotting_tools/model_files/*%s*'%(species))
all_models = []
for i in model_files:
    i = i.replace("/work/home/db876/plotting_tools/model_files/", "")
    model_split = i.split('_%s_'%(species))
    part1 = model_split[0]
    part2 = model_split[1]
    start_year = part2[:4]
    end_year = part2[5:9]
    year_range = start_year+'_'+end_year
    model = part1[:-8]
    #check for version, grid_size and metorology
    #if only 1 there is version and grid_size, if 2 there is version and grid_size, if 3 there it is version, grid_size and meterology.
    extra_data = part2[10:-3]
    extra_data_split = extra_data.split('_')
    
    if (len(extra_data_split) == 1) & (extra_data_split != ['']):
        version = extra_data_split[0]
        try:
            key = version_dict[model]
            if version not in key:
                key.append(version) 
                version_dict[model] = key
        except:
            version_dict[model] = [version]
    elif len(extra_data_split) == 2:
        version = extra_data_split[0]
        grid = extra_data_split[1]
        try:
            key = version_dict[model]
            if version not in key:
                key.append(version) 
                version_dict[model] = key
        except:
            version_dict[model] = [version]
        try:
            key = grid_dict[model]
            if grid not in key:
                key.append(grid) 
                grid_dict[model] = key
        except:
            grid_dict[model] = [grid]

    elif len(extra_data_split) == 3:
        version = extra_data_split[0]
        grid = extra_data_split[1]
        met = extra_data_split[2]
        try:
            key = version_dict[model]
            if version not in key:
                key.append(version) 
                version_dict[model] = key
        except:
            version_dict[model] = [version]
        try:
            key = grid_dict[model]
            if grid not in key:
                key.append(grid) 
                grid_dict[model] = key
        except:
            grid_dict[model] = [grid] 
        try:
            key = met_dict[model]
            if met not in key:                                                                                                                                                                                                           
                key.append(met) 
                met_dict[model] = key
        except:
            met_dict[model] = [met]     

    try:
        key = model_dict[model]
        if year_range not in key:
            key.append(year_range) 
            model_dict[model] = key
    except:
        model_dict[model] = [year_range]   
    
#get model and date range
model_name = raw_input('\nChoose Model.\n%s\n'%('   '.join(i for i in model_dict)))
model_range = model_dict[model_name]
if len(model_range) == 1:
    model_range = model_range[0]
else:
    model_range = raw_input('\nChoose Date Range.\n%s\n'%('   '.join(i for i in model_range)))

#get model version
valid_versions = []
try:
    versions = version_dict[model_name]
    for mf in model_files:
        for i in versions:
            test_string = '%s_SURFACE_%s_%s_%s'%(model_name,species,model_range,i)
            if test_string in mf:
                valid_versions.append(i)     
    valid_versions = set(valid_versions)
    valid_versions = [i for i in valid_versions]   
except:
    valid_versions = ''
if valid_versions != '':
    if len(valid_versions) > 1:
        version = raw_input('\nChoose Model Version.\n%s\n'%('   '.join(i for i in valid_versions)))
    else:
        version = valid_versions[0]

#get grid version
valid_grids = []
try:
    grids = grid_dict[model_name]
    for mf in model_files:
        for i in grids:
            test_string = '%s_SURFACE_%s_%s_%s_%s'%(model_name,species,model_range,version,i)
            if test_string in mf:
                valid_grids.append(i)
    valid_grids = set(valid_grids)
    valid_grids = [i for i in valid_grids]
except:
    valid_grids = ''
if valid_grids != '':
    if len(valid_grids) > 1:
        grid = raw_input('\nChoose Model Grid.\n%s\n'%('   '.join(i for i in valid_grids)))
    else:
        grid = valid_grids[0]

#get met version
valid_mets = []
try:
    mets = met_dict[model_name]
    for mf in model_files:
        for i in mets:
            test_string = '%s_SURFACE_%s_%s_%s_%s_%s'%(model_name,species,model_range,version,grid,i)
            if test_string in mf:
                valid_mets.append(i)
    valid_mets = set(valid_mets)
    valid_mets = [i for i in valid_mets]
except:
    valid_mets = ''
if valid_mets != '':
    if len(valid_mets) > 1:
        met = raw_input('\nChoose Model Meteorology.\n%s\n'%('   '.join(i for i in valid_mets)))
    else:
        met = valid_mets[0]

#put together model file name
if (valid_versions == '') & (valid_grids == '') & (valid_mets == ''):
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'.nc'
    version,grid,met = '','',''
elif (valid_versions != '') & (valid_grids == '') & (valid_mets == ''):
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'.nc'
    grid,met = '',''
elif (valid_versions != '') & (valid_grids != '') & (valid_mets == ''):
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'_'+grid+'.nc'
    met = ''
elif (valid_versions != '') & (valid_grids != '') & (valid_mets != ''):
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'_'+grid+'_'+met+'.nc'

root_grp = Dataset(model_f)
model_var = root_grp.variables[species.lower()][:]
model_date = root_grp.variables['date'][:]
model_time = root_grp.variables['time'][:]
lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
grid_size = root_grp.variables['grid_size'][:]
grid_size = grid_size[0]
model_var_mask = np.ma.masked_where(model_var<0,model_var)

#model_var_t = model_var*1e9
test = model_var <= 0
invs =  model_var[test]
print invs

gridbox_count = len(lat_c)*len(lon_c)

#average model gridboxes across time
model_ave = np.ma.average(model_var_mask,axis=0)
model_ave = model_ave*1e9
print model_ave

#process obs
#----------------------------------------
#data range

plot_obs = raw_input('\nOverlay Observation Averages? y or n?\n')

#read in obs data
try:
    root_grp = Dataset('../process/GLOBAL_SURFACE_%s_%s_%s.nc'%(species,model_range[:4],model_range[-4:]))
    note = ''
except:
    all_obs_files = glob.glob('../process/GLOBAL_SURFACE_%s*'%(species)) 
    valid_obs_files = []
    for f in all_obs_files:
        root_grp = Dataset(f)
        valid_refs_dict = root_grp.groups
        valid_refs = []
        for i in valid_refs_dict.keys():
            valid_refs = np.append(valid_refs,i)
        valid_ref = valid_refs[0]
        site_group = root_grp.groups[valid_ref]
        test_obs_var = site_group.variables[species.lower()][:]
        if len(model_var) == len(test_obs_var):
            valid_obs_files.append(f)
    print valid_obs_files
    
    if len(valid_obs_files) == 1:
        m_start_year, m_end_year = start_year,end_year
        file = valid_obs_files[0]
        start_year = file[-12:-8]
        end_year = file[-7:-3]
        note = '- BE AWARE OBS. IS %s to %s,  MODEL IS %s to %s'%(start_year,end_year,m_start_year,m_end_year) 
        print '\nNo equivalent observational data range. Using range between %s and %s that matches in shape.\n'%(start_year,end_year)    
        root_grp = Dataset(file)        

    elif len(valid_obs_files) > 1:
        print '\nNo equivalent observational data range. Choose a range that matches in shape to compare against.\n' 
        s_years = []
        e_years = []
        for i in range(len(valid_files)):
            file = valid_files[i]
            s_years.append(file[-12:-8])
            e_years.append(file[-7:-3])
        ranges = ['%s_%s'%(a,b) for a,b in zip(s_years,e_years)] 
        chosen_range = raw_input('\nChoose Range.\n%s\n'%('   '.join(i for i in ranges)))
        chosen_file = valid_files[ranges.index(chosen_range)]
        m_start_year, m_end_year = start_year,end_year
        start_year = int(chosen_range[:4])
        end_year = int(chosen_range[5:])
        note = '- BE AWARE OBS. IS %s to %s,  MODEL IS %s to %s'%(start_year,end_year,m_start_year,m_end_year)
        root_grp = Dataset(chosen_file)
    else:
        print '\nNo equivalent observational data range. And no observational data range that matches in shape. Adios.\n'
        1+'a'

valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)
       
valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

obs_aves = []
obs_lats = []
obs_lons = []
obs_alts = []
obs_groups = []

for ref in valid_refs:  
    #read in specific site data
    site_group = root_grp.groups[ref]

    #read in variables for site
    obs_var = site_group.variables[species.lower()][:]
    obs_lat = site_group.latitude
    obs_lon = site_group.longitude
    obs_alt = site_group.altitude
    obs_group = site_group.process_group
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_ave = np.ma.average(obs_var_mask)
    
    obs_aves.append(obs_ave)
    obs_lats.append(obs_lat)
    obs_lons.append(obs_lon)
    obs_alts.append(obs_alt)
    obs_groups.append(obs_group)

#---------------------------------------------------------------------------
#check which griboxes each obs_lat & obs_lon convergance point lie in & then return central convergance point of lat and lon in that gridbox

obs_lat_centres = []
obs_lon_centres = []
for obs_lat,obs_lon in zip(obs_lats, obs_lons):
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    obs_lat_centres.append(lat_c[lat_n])
    obs_lon_centres.append(lon_c[lon_n])
 
#---------------------------------------------------------------------------        
#set up plot
fig =plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')
 
print lat_e
print lon_e
print lon_c


#bpch_fnames = glob.glob('/work/home/db876/modelling/GEOS_CHEM/v90103/2x2.5/2x2.5_GEOS5_fullgrid/run/ctm.bpch*')
#bpch_fnames = modules.natsorted(bpch_fnames)
 
#bpch_fnames =  bpch_fnames[2:]
 
#print model_ave.shape
 
#for bpch_fname in bpch_fnames:
#    f = bpch(bpch_fname)
#    group = f.groups['IJ-AVG-$']
#    ozone = group.variables['O3']
#    ozone = ozone[:,0,:,:]
#    try:
#        ozone_add = np.concatenate((ozone_add,ozone),axis=0)
#    except:
#        ozone_add = ozone
#    print ozone_add.shape
#model_ave = np.average(ozone_add,axis=0)
#print model_ave.shape        

#setup basemap projection
m = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],\
                   llcrnrlon=lon_e[0],\
                    urcrnrlon=lon_e[-1],\
                   resolution='c')


m.drawcoastlines()
m.drawmapboundary()
parallels = np.arange(-90,91,15)
meridians = np.arange(-180,151,30)
plt.xticks(meridians)
plt.yticks(parallels)
m.drawparallels(parallels)
m.drawmeridians(meridians)

X,Y = m(obs_lon_centres,obs_lat_centres)

poly = m.pcolor(lon_e, lat_e, model_ave, vmin=np.min(model_ave), vmax=np.max(model_ave),cmap = plt.cm.coolwarm)

if plot_obs == 'y':
    pl = m.scatter(X,Y,c=obs_aves,s=40, vmin=np.min(obs_aves), vmax=np.max(obs_aves),edgecolor='None',linewidth=0.5,cmap=plt.cm.coolwarm)

cb = plt.colorbar(poly, ax = m.ax,shrink=0.8,orientation = 'horizontal', format='%.2f')

cb.set_label('Concentration (ppb)', fontsize = 16)    
plt.xlabel('Longitude',fontsize = 20)
plt.ylabel('Latitude',fontsize = 20)
plt.title('%s average surface %s between %s & %s'%(model_name,species,start_year,end_year),fontsize=20)
cb.ax.tick_params(labelsize=16)

mng = plt.get_current_fig_manager()
mng.window.wm_geometry("+2500+100")

plt.show()




