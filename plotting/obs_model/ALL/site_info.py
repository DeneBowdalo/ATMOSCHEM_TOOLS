from netCDF4 import Dataset
import numpy as np
import os
import glob

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-2]
print '\nSpecies is %s\n'%(species)

timeres = raw_input('\nTime resolution? H, D or M?\n')

#read in obs data
root_grp = Dataset('../process/GLOBAL_SURFACE_%s_1970_2015_%s_ALL.nc'%(species,timeres))
valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

print '%i Sites'%(len(valid_refs))
site_ref = raw_input('Choose site from list.\n%s\n'%('   '.join(i for i in valid_refs)))

#read in specific site data
site_group = root_grp.groups[site_ref]

print site_group
