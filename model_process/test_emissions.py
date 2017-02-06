from bpch import bpch
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

std_file = '../standard/HEMCO_diagnostics.200902010000.nc'
alt_file = 'HEMCO_diagnostics.200902010000.nc'

std_root = Dataset(std_file)
alt_root = Dataset(alt_file)

m_set = raw_input('NOX, ANMVOC or BNMVOC?\n')
if m_set == 'NOX':
    spec = raw_input('\nNO or NO2?\n')
elif m_set == 'ANMVOC':
    spec = raw_input('\nACET, ALD2, ALK4, CH2O, C2H6, C3H8, MACR, MEK, PRPE or RCHO?\n')
elif m_set == 'BNMVOC':
    spec = raw_input('\nACET, ALD2, ISOP or PRPE?\n')

#no emissions
if m_set == 'NOX':
    if spec == 'NO':
        std_soil = std_root.variables['NO_SOIL'][0,:,:]
        alt_soil = alt_root.variables['NO_SOIL'][0,:,:]
        std_li = std_root.variables['NO_LIGHTNING'][0,:,:]
        alt_li = alt_root.variables['NO_LIGHTNING'][0,:,:]
        std_plane = std_root.variables['NO_PLANE'][0,:,:]
        alt_plane = alt_root.variables['NO_PLANE'][0,:,:]
        std_biob = std_root.variables['NO_BIOMASS'][0,:,:]
        alt_biob = alt_root.variables['NO_BIOMASS'][0,:,:]
        std_biof = std_root.variables['NO_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['NO_BIOFUEL'][0,:,:]
        std_anth = std_root.variables['NO_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['NO_ANTHRO'][0,:,:]
        std_ship = std_root.variables['NO_SHIP'][0,:,:]
        alt_ship = alt_root.variables['NO_SHIP'][0,:,:]
    elif spec == 'NO2':
        std_anth = std_root.variables['NO2_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['NO2_ANTHRO'][0,:,:]
        std_biof = std_root.variables['NO2_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['NO2_BIOFUEL'][0,:,:]
        std_ship = std_root.variables['NO2_SHIP'][0,:,:]
        alt_ship = alt_root.variables['NO2_SHIP'][0,:,:]   

#no2 emissions

#anth nmvoc emissions
elif m_set == 'ANMVOC':
    if spec == 'ACET':
        std_anth = std_root.variables['ACET_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['ACET_ANTHRO'][0,:,:] 
        std_biof = std_root.variables['ACET_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['ACET_BIOFUEL'][0,:,:]
        std_ship = std_root.variables['ACET_SHIP'][0,:,:]
        alt_ship = alt_root.variables['ACET_SHIP'][0,:,:]
        std_plane = std_root.variables['ACET_PLANE'][0,:,:]
        alt_plane = alt_root.variables['ACET_PLANE'][0,:,:]
        std_biob = std_root.variables['ACET_BIOMASS'][0,:,:]
        alt_biob = alt_root.variables['ACET_BIOMASS'][0,:,:]
    elif spec == 'ALD2':
        std_anth = std_root.variables['ALD2_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['ALD2_ANTHRO'][0,:,:]
        std_biof = std_root.variables['ALD2_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['ALD2_BIOFUEL'][0,:,:]
        std_ship = std_root.variables['ALD2_SHIP'][0,:,:]
        alt_ship = alt_root.variables['ALD2_SHIP'][0,:,:]
        std_plane = std_root.variables['ALD2_PLANE'][0,:,:]
        alt_plane = alt_root.variables['ALD2_PLANE'][0,:,:]
        std_biob = std_root.variables['ALD2_BIOMASS'][0,:,:]
        alt_biob = alt_root.variables['ALD2_BIOMASS'][0,:,:]
    elif spec == 'ALK4':
        std_anth = std_root.variables['ALK4_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['ALK4_ANTHRO'][0,:,:]
        std_biof = std_root.variables['ALK4_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['ALK4_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['ALK4_SHIP'][0,:,:]
        alt_ship = alt_root.variables['ALK4_SHIP'][0,:,:]
        std_plane = std_root.variables['ALK4_PLANE'][0,:,:]
        alt_plane = alt_root.variables['ALK4_PLANE'][0,:,:]
        std_biob = std_root.variables['ALK4_BIOMASS'][0,:,:]
        alt_biob = alt_root.variables['ALK4_BIOMASS'][0,:,:]
    elif spec == 'BENZ':
       # std_anth = std_root.variables['BENZ_ANTHRO'][0,:,:]
        #alt_anth = alt_root.variables['BENZ_ANTHRO'][0,:,:]
        std_biof = std_root.variables['BENZ_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['BENZ_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['BENZ_SHIP'][0,:,:]
        alt_ship = alt_root.variables['BENZ_SHIP'][0,:,:]
    elif spec == 'CH2O':
        std_anth = std_root.variables['CH2O_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['CH2O_ANTHRO'][0,:,:]
        std_biof = std_root.variables['CH2O_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['CH2O_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['CH2O_SHIP'][0,:,:]
        alt_ship = alt_root.variables['CH2O_SHIP'][0,:,:]
        std_plane = std_root.variables['CH2O_PLANE'][0,:,:]
        alt_plane = alt_root.variables['CH2O_PLANE'][0,:,:]
        std_biob = std_root.variables['CH2O_BIOMASS'][0,:,:]
        alt_biob = alt_root.variables['CH2O_BIOMASS'][0,:,:]
    elif spec == 'C2H2':
        std_anth = std_root.variables['C2H2_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['C2H2_ANTHRO'][0,:,:]
        std_biof = std_root.variables['C2H2_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['C2H2_BIOFUEL'][0,:,:] 
    elif spec == 'C2H4':
        std_anth = std_root.variables['C2H4_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['C2H4_ANTHRO'][0,:,:]
        std_biof = std_root.variables['C2H4_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['C2H4_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['C2H4_SHIP'][0,:,:]
        alt_ship = alt_root.variables['C2H4_SHIP'][0,:,:]
    elif spec == 'C2H6':
        std_anth = std_root.variables['C2H6_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['C2H6_ANTHRO'][0,:,:]
        std_biof = std_root.variables['C2H6_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['C2H6_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['C2H6_SHIP'][0,:,:]
        alt_ship = alt_root.variables['C2H6_SHIP'][0,:,:]
        std_plane = std_root.variables['C2H6_PLANE'][0,:,:]
        alt_plane = alt_root.variables['C2H6_PLANE'][0,:,:]
        std_biob = std_root.variables['C2H6_BIOMASS'][0,:,:]
        alt_biob = alt_root.variables['C2H6_BIOMASS'][0,:,:]
    elif spec == 'C3H8':
        std_anth = std_root.variables['C3H8_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['C3H8_ANTHRO'][0,:,:]
        std_biof = std_root.variables['C3H8_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['C3H8_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['C3H8_SHIP'][0,:,:]
        alt_ship = alt_root.variables['C3H8_SHIP'][0,:,:]
        std_plane = std_root.variables['C3H8_PLANE'][0,:,:]
        alt_plane = alt_root.variables['C3H8_PLANE'][0,:,:]
        std_biob = std_root.variables['C3H8_BIOMASS'][0,:,:]
        alt_biob = alt_root.variables['C3H8_BIOMASS'][0,:,:]
    elif spec == 'MACR':
        std_anth = std_root.variables['MACR_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['MACR_ANTHRO'][0,:,:]
        std_biof = std_root.variables['MACR_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['MACR_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['MACR_SHIP'][0,:,:]
        alt_ship = alt_root.variables['MACR_SHIP'][0,:,:]
        std_plane = std_root.variables['MACR_PLANE'][0,:,:]
        alt_plane = alt_root.variables['MACR_PLANE'][0,:,:]
    elif spec == 'MEK':
        std_anth = std_root.variables['MEK_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['MEK_ANTHRO'][0,:,:]
        std_biof = std_root.variables['MEK_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['MEK_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['MEK_SHIP'][0,:,:]
        alt_ship = alt_root.variables['MEK_SHIP'][0,:,:]
        std_biob = std_root.variables['MEK_BIOMASS'][0,:,:]
        alt_biob = alt_root.variables['MEK_BIOMASS'][0,:,:]
    elif spec == 'PRPE':
        std_anth = std_root.variables['PRPE_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['PRPE_ANTHRO'][0,:,:]
        std_biof = std_root.variables['PRPE_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['PRPE_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['PRPE_SHIP'][0,:,:]
        alt_ship = alt_root.variables['PRPE_SHIP'][0,:,:]
        std_plane = std_root.variables['PRPE_PLANE'][0,:,:]
        alt_plane = alt_root.variables['PRPE_PLANE'][0,:,:]
        std_biob = std_root.variables['PRPE_BIOMASS'][0,:,:]
        alt_biob = alt_root.variables['PRPE_BIOMASS'][0,:,:]
    elif spec == 'RCHO':
        std_anth = std_root.variables['RCHO_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['RCHO_ANTHRO'][0,:,:]
        std_biof = std_root.variables['RCHO_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['RCHO_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['RCHO_SHIP'][0,:,:]
        alt_ship = alt_root.variables['RCHO_SHIP'][0,:,:]
        std_plane = std_root.variables['RCHO_PLANE'][0,:,:]
        alt_plane = alt_root.variables['RCHO_PLANE'][0,:,:]
    elif spec == 'TOLU':
        std_anth = std_root.variables['TOLU_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['TOLU_ANTHRO'][0,:,:]
        std_biof = std_root.variables['TOLU_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['TOLU_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['TOLU_SHIP'][0,:,:]
        alt_ship = alt_root.variables['TOLU_SHIP'][0,:,:]
    elif spec == 'XYLE':
        std_anth = std_root.variables['XYLE_ANTHRO'][0,:,:]
        alt_anth = alt_root.variables['XYLE_ANTHRO'][0,:,:]
        std_biof = std_root.variables['XYLE_BIOFUEL'][0,:,:]
        alt_biof = alt_root.variables['XYLE_BIOFUEL'][0,:,:] 
        std_ship = std_root.variables['XYLE_SHIP'][0,:,:]
        alt_ship = alt_root.variables['XYLE_SHIP'][0,:,:]

elif m_set == 'BNMVOC':
    if spec == 'ACET':
        std_nat = std_root.variables['ACET_NATURAL'][0,:,:]
        alt_nat = alt_root.variables['ACET_NATURAL'][0,:,:]
        std_oc = std_root.variables['ACET_OCEANIC'][0,:,:]
        alt_oc = alt_root.variables['ACET_OCEANIC'][0,:,:]
    elif spec == 'ALD2':
        std_nat = std_root.variables['ALD2_NATURAL'][0,:,:]
        alt_nat = alt_root.variables['ALD2_NATURAL'][0,:,:]
    elif spec == 'C2H4':
        std_nat = std_root.variables['C2H4_NATURAL'][0,:,:]
        alt_nat = alt_root.variables['C2H4_NATURAL'][0,:,:]
    elif spec == 'ISOP':
        std_nat = std_root.variables['ISOP_NATURAL'][0,:,:]
        alt_nat = alt_root.variables['ISOP_NATURAL'][0,:,:]
    elif spec == 'PRPE':
        std_nat = std_root.variables['PRPE_NATURAL'][0,:,:]
        alt_nat = alt_root.variables['PRPE_NATURAL'][0,:,:]

if m_set == 'NOX':
    if spec == 'NO':
        param = raw_input('\nANTH, BIOB, BIOF, LI, PLANE, SHIP or SOIL?\n')
    elif spec == 'NO2':
        param = raw_input('\nANTH or BIOF or SHIP?\n')

elif m_set == 'ANMVOC':
    if (spec == 'ACET') or (spec == 'ALD2') or (spec == 'ALK4') or (spec == 'CH2O') or (spec == 'C2H6') or (spec == 'C3H8') or (spec == 'PRPE'):
        param = raw_input('\nANTH, BIOB, BIOF, PLANE or SHIP?\n')
    elif (spec == 'BENZ') or (spec == 'C2H4') or (spec == 'TOLU') or (spec == 'XYLE'):
        param = raw_input('\nANTH, BIOF or SHIP?\n')
    elif spec == 'C2H2':
        param = raw_input('\nANTH or BIOF\n')
    elif spec == 'MEK':
        param = raw_input('\nANTH, BIOB, BIOF or SHIP?\n')
    elif (spec == 'MACR') or (spec == 'RCHO'):
        param = raw_input('\nANTH, BIOF, PLANE or SHIP?\n')

elif m_set == 'BNMVOC':
    if spec == 'ACET':
        param = raw_input('\nOC or NAT?\n')
    elif (spec == 'ALD2') or (spec == 'C2H4') or (spec == 'ISOP') or (spec == 'PRPE'):
        param = 'NAT'

if param == 'SOIL':
    diff_ratio = std_soil/alt_soil
elif param == 'LI':
    diff_ratio = std_li/alt_li
elif param == 'PLANE':
    diff_ratio = std_plane/alt_plane
elif param == 'BIOB':
    diff_ratio = std_biob/alt_biob
elif param == 'BIOF':
    diff_ratio = std_biof/alt_biof
elif param == 'ANTH':
    diff_ratio = std_anth/alt_anth
elif param == 'SHIP':
    diff_ratio = std_ship/alt_ship
elif param == 'NAT':
    diff_ratio = std_nat/alt_nat
elif param == 'OC':
    diff_ratio = std_oc/alt_oc

test = np.isnan(diff_ratio)
diff_ratio[test] = 1

print np.nanmin(diff_ratio)

plt.pcolor(diff_ratio)
plt.colorbar()

plt.show()



