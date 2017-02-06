from netCDF4 import Dataset
import numpy as np

#OUTPUT RESOLUTION CAN BE HOURLY (H), HOURLY & DAILY (HD), HOURLY, DAILY AND MONTHLY (HDM)
output_res = 'HDM'

output_set = 'S'

stats_data = Dataset('STATS_%s.nc'%(output_res+output_set))


c1 = stats_data.variables['invalid_nometa_count'][0] + stats_data.variables['invalid_nokeymeta_count'][0]
c2 = stats_data.variables['invalid_anyvaliddata_count'][0]
c3 = stats_data.variables['invalid_resolution_count'][0]
c4 = stats_data.variables['invalid_badmeasurementmethod_count'][0]
c5 = stats_data.variables['invalid_duplicatesites_count'][0]
c6 = stats_data.variables['invalid_rawclass_count'][0] + stats_data.variables['invalid_anthromeclass_count'][0] 
c7 = stats_data.variables['invalid_altitude_count'][0]
c8 = stats_data.variables['invalid_night_count'][0]
c9 = stats_data.variables['invalid_representativeness_count'][0]
c10 = stats_data.variables['invalid_extreme_count'][0]
c11 = stats_data.variables['invalid_partialyear_count'][0]
c12 = stats_data.variables['n_final_count'][0]

n1 = stats_data.variables['n_obs_all'][0]
n2 = stats_data.variables['n_obs_after_nometa'][0]
n3 = stats_data.variables['n_obs_after_flagsandlod'][0]
n4 = stats_data.variables['n_obs_after_duplicatepoints'][0]
n5 = stats_data.variables['n_obs_after_anyvaliddata'][0]
n6 = stats_data.variables['n_obs_after_nokeymeta'][0]
n7 = stats_data.variables['n_obs_after_resolution'][0]
n8 = stats_data.variables['n_obs_after_badmeasurementmethod'][0]
n9 = stats_data.variables['n_obs_after_duplicatesites'][0]
n10 = stats_data.variables['n_obs_after_rawclass'][0]
n11 = stats_data.variables['n_obs_after_anthromeclass'][0]
n12 = stats_data.variables['n_obs_after_altitude'][0]
n13 = stats_data.variables['n_obs_after_night'][0]
n14 = stats_data.variables['n_obs_after_representativeness'][0]
n15 = stats_data.variables['n_obs_after_extreme'][0]
n16 = stats_data.variables['n_obs_after_partialyear'][0]

total_sites = c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+c11+c12
print 'n sites'
print 'n sites total = %s'%(total_sites) 
print 'n sites after nometa = %s'%(total_sites-c1)
print 'n sites after anyvaliddata = %s'%(total_sites-c1-c2)
print 'n sites after resolution = %s'%(total_sites-c1-c2-c3)
print 'n sites after badmeasurmentmethod = %s'%(total_sites-c1-c2-c3-c4)
print 'n sites after duplicated = %s'%(total_sites-c1-c2-c3-c4-c5)
print 'n sites after urban = %s'%(total_sites-c1-c2-c3-c4-c5-c6)
print 'n sites after altitude = %s'%(total_sites-c1-c2-c3-c4-c5-c6-c7)
print 'n sites after night = %s'%(total_sites-c1-c2-c3-c4-c5-c6-c7-c8)
print 'n sites after representativeness = %s'%(total_sites-c1-c2-c3-c4-c5-c6-c7-c8-c9)
print 'n sites after extreme = %s'%(total_sites-c1-c2-c3-c4-c5-c6-c7-c8-c9-c10)
print 'n sites after partialyear = %s'%(total_sites-c1-c2-c3-c4-c5-c6-c7-c8-c9-c10-c11)
print 'n final = %s'%(c12)

badmeta_diff = n5-n6
nometa = n2-badmeta_diff
flagsandlod = n3-badmeta_diff
duppoints = n4-badmeta_diff
anyvalid = n5-badmeta_diff

print '\n'
print 'counts'
print 'n obs total = %s'%(n1) 
print 'n obs after nometa = %s'%(nometa) 
print 'n obs after flagsandlod = %s'%(flagsandlod)  
print 'n obs after duplicate points = %s'%(duppoints)  
print 'n obs after anyvaliddata = %s'%(anyvalid) 
print 'n obs after resolution = %s'%(n7)
print 'n obs after bad measurement method = %s'%(n8) 
print 'n obs after duplicate sites = %s'%(n9) 
print 'n obs after urban sites = %s'%(n11)
print 'n obs after altitude = %s'%(n12) 
print 'n obs after night = %s'%(n13) 
print 'n obs after representativeness = %s'%(n14) 
print 'n obs after extreme = %s'%(n15) 
print 'n obs after partialyear = %s'%(n16)

