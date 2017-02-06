import numpy as np
from netCDF4 import Dataset
import modules

data = Dataset('GLOBAL_SURFACE_O3_1970_2015_H_HDMN.nc')
obs_refs = []
diff_array = []
obs_refs_dict = data.groups
for i in obs_refs_dict.keys():
    i = i.encode('ascii')
    obs_refs = np.append(obs_refs,i)

for r in range(len(obs_refs)):
    print r
    obs_ref = obs_refs[r]
    vals = data.groups[obs_ref].variables['o3'][:]
    test_valid = vals >= 0
    min_val_ind = np.min(np.where(test_valid == True)[0])
    max_val_ind = np.max(np.where(test_valid == True)[0]) 
    diff = max_val_ind - min_val_ind
    diff_array = np.append(diff_array,diff)

longest_refs = np.array([x for (y,x) in sorted(zip(diff_array,obs_refs),reverse=True)])
longest_tags = modules.get_tags(longest_refs)
longest_diff = sorted(diff_array,reverse=True)

c = 0
for j in range(len(longest_refs)):
    if longest_tags[j] == 'NA':    
        print longest_refs[j],longest_tags[j],longest_diff[j]/24./365.25
        c+=1
        if c == 20:
            break

print

c = 0
for j in range(len(longest_refs)):
    if longest_tags[j] == 'EU':    
        print longest_refs[j],longest_tags[j],longest_diff[j]/24./365.25
        c+=1
        if c == 20:
            break

print

c = 0
for j in range(len(longest_refs)):
    if longest_tags[j] == 'O':    
        print longest_refs[j],longest_tags[j],longest_diff[j]/24./365.25
        c+=1
        if c == 20:
            break

print

c = 0
for j in range(len(longest_refs)):
    if longest_tags[j] == 'OC':    
        print longest_refs[j],longest_tags[j],longest_diff[j]/24./365.25
        c+=1
        if c == 20:
            break
print

c = 0
for j in range(len(longest_refs)):                                                                                                                            
    if longest_tags[j] == 'AF':
        print longest_refs[j],longest_tags[j],longest_diff[j]/24./365.25
        c+=1
        if c == 20:
            break

print

c = 0
for j in range(len(longest_refs)):                                                                                                                            
    if longest_tags[j] == 'ANT':
        print longest_refs[j],longest_tags[j],longest_diff[j]/24./365.25
        c+=1
        if c == 20:
            break

print

c = 0
for j in range(len(longest_refs)):                                                                                                                            
    if longest_tags[j] == 'ARC':
        print longest_refs[j],longest_tags[j],longest_diff[j]/24./365.25
        c+=1
        if c == 20:
            break

print

c = 0
for j in range(len(longest_refs)):
    if longest_tags[j] == 'EA':                                                                                                                               
        print longest_refs[j],longest_tags[j],longest_diff[j]/24./365.25
        c+=1
        if c == 20:
            break

print

c = 0
for j in range(len(longest_refs)):
    if longest_tags[j] == 'MA':                                                                                                                               
        print longest_refs[j],longest_tags[j],longest_diff[j]/24./365.25
        c+=1
        if c == 20:
            break

