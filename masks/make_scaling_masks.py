import modules
from netCDF4 import Dataset
import numpy as np

vers = 'v6'

#0.1x0.1 
lon_e_01 = np.arange(-180.05,179.95+0.05,0.1)
lat_e_01 = np.arange(-90.05,90.05+0.05,0.1)
lat_e_01[0]=-90.
lat_e_01[-1]=90.

lon_c_01 = np.arange(-180.0,180.0,0.1)
lat_c_01 = np.arange(-90,90+0.05,0.1)
lat_c_01[0]=-89.975
lat_c_01[-1]=89.975

#0.25x0.25 
lon_e_025 = np.arange(-180.125,179.875+0.125,0.25)
lat_e_025 = np.arange(-90.125,90.125+0.125,0.25)
lat_e_025[0]=-90.
lat_e_025[-1]=90.

lon_c_025 = np.arange(-180.0,180.0,0.25)
lat_c_025 = np.arange(-90,90+0.125,0.25)
lat_c_025[0]=-89.9375
lat_c_025[-1]=89.9375

#1x1 
lon_e_1x1 = np.arange(-180.5,179.5+0.5,1.)
lat_e_1x1 = np.arange(-90.5,90.5+0.5,1.)
lat_e_1x1[0]=-90.
lat_e_1x1[-1]=90.

lon_c_1x1 = np.arange(-180.0,180.0,1.)
lat_c_1x1 = np.arange(-90,90+0.5,1.)
lat_c_1x1[0]=-89.75
lat_c_1x1[-1]=89.75

#4x5
lon_e_4x5 = np.arange(-182.5,177.5+2.5,5)
lat_e_4x5 = np.arange(-92,92+2,4.)
lat_e_4x5[0]=-90.
lat_e_4x5[-1]=90.

lon_c_4x5 = np.arange(-180.0,180.0,5.)
lat_c_4x5 = np.arange(-90,90+2,4.)
lat_c_4x5[0]=-89.
lat_c_4x5[-1]=89.


#area boundaries
#stretch some to incorporate additional regions (i.e. china, mexico)
areas = ['SW_NA','C_NA','NW_NA','NE_NA','CE_NA','SE_NA','S_NA','NW_EU','C_EU','N_EU','E_EU','S_EU','SW_EU','NE_AS']
area_boundaries = {'CE_NA':[37.01,43.51,-90.01,-49.99],'SE_NA':[19.99,37.01,-90.01,-64.99],'S_NA':[19.99,40.99,-107.99,-90.1],'SW_NA':[19.99,40.01,-130.01,-107.99],'C_NA':[40.99,80.01,-107.99,-90.01],'NE_NA':[43.51,80.01,-90.01,-49.99],'NW_NA':[40.01,80.01,-170.01,-107.99],'N_EU':[55.01,80.01,3.99,45.01],'C_EU':[46.01,55.01,3.99,12.51],'SW_EU':[29.99,46.01,-10.01,0.01],'S_EU':[29.99,46.01,0.01,45.01],'E_EU':[46.01,55.01,12.51,45.01],'NW_EU':[46.01,70.01,-15.01,3.99],'NE_AS':[10,80,90,160]}


#param area scale factors
#anmvoc_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_NA':[0.5,0.5,1.,1.,1.,1.,1.,1.,1.,0.5,0.5,0.5],'NE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'N_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_EU':[1.,1.,2.,2.,1.,1.,1.,1.,1.,1.,2.,1.],'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'E_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[0.5,0.5,0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5]}
#bnmvoc_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'N_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'E_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]}
#co_scale = {'CE_NA':[0.5,0.5,0.5,1.,1.,1.,1.,1.,1.,1.,0.5,0.5],'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'N_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_EU':[2.,2.,2.,2.,2.,4.,4.,4.,4.,2.,2.,2.],'E_EU':[2.,2.,1.,1.,1.,0.5,1.,1.,1.,2.,2.,2.],'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]}
#drydepo3_scale = {'CE_NA':[1.,1.,1.,1.,2.,4.,4.,4.,2.,2.,1.,1.],'SE_NA':[2.,2.,2.,2.,2.,4.,4.,4.,4.,4.,2.,2.],'S_NA':[2.,2.,2.,2.,2.,2.,4.,4.,2.,2.,2.,2.],'SW_NA':[1.,2.,2.,2.,2.,2.,4.,4.,2.,2.,1.,1.],'C_NA':[1.,1.,1.,2.,2.,2.,4.,4.,4.,1.,1.,1.],'NE_NA':[1.,1.,2.,2.,2.,2.,4.,4.,4.,2.,1.,1.],'NW_NA':[2.,2.,2.,2.,2.,2.,4.,4.,4.,2.,2.,2.],'N_EU':[0.5,0.5,0.5,1.,2.,2.,2.,2.,2.,1.,0.5,0.5],'C_EU':[1.,1.,1.,1.,2.,2.,2.,2.,2.,1.,1.,1.],'SW_EU':[2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.],'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'E_EU':[1.,1.,1.,1.,1.,1.,2.,2.,2.,1.,1.,1.],'NW_EU':[1.,1.,2.,2.,2.,4.,4.,4.,4.,2.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,2.,4.,4.,4.,1.,1.,1.]}
#nox_scale = {'CE_NA':[0.5,0.5,1.,2.,2.,2.,2.,2.,2.,1.,0.5,0.5],'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_NA':[1.,1.,2.,2.,1.,0.5,0.5,1.,2.,2.,2.,1.],'SW_NA':[2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.],'C_NA':[1.,1.,2.,2.,2.,1.,1.,1.,1.,2.,1.,1.],'NE_NA':[2.,2.,2.,2.,2.,1.,2.,2.,2.,2.,2.,2.],'NW_NA':[2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.],'N_EU':[1.,2.,2.,1.,0.5,0.5,0.5,0.5,0.5,2.,2.,1.],'C_EU':[0.5,0.5,1.,2.,2.,2.,2.,2.,2.,1.,1.,0.5],'SW_EU':[2.,2.,2.,2.,1.,1.,1.,1.,1.,1.,1.,2.],'S_EU':[1.,1.,2.,2.,0.25,0.25,0.25,0.25,0.25,2.,2.,1.],'E_EU':[0.5,0.5,1.,0.5,0.5,0.5,0.5,0.5,1.,1.,1.,0.5],'NW_EU':[1.,1.,1.,1.,2.,2.,2.,2.,2.,1.,1.,1.],'NE_AS':[0.5,0.5,0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5]}

#anmvoc_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SE_NA':[0.25,0.25,0.25,1.,1.,1.,1.,1.,1.,1.,0.25,0.25],'S_NA':[0.25,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.25],'SW_NA':[2.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,2.],'C_NA':[0.25,0.25,1.,1.,1.,1.,1.,1.,1.,0.5,0.25,0.25],'NE_NA':[2.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,2.],'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'N_EU':[4.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,4.],'C_EU':[4.,4.,4.,4.,1.,1.,1.,1.,1.,1.,4.,4.],'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'E_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_EU':[2.,2.,4.,4.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[0.25,0.25,0.25,0.5,0.5,1.,1.,1.,1.,1.,0.25,0.25]}
#bnmvoc_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'N_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'E_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]}
#co_scale = {'CE_NA':[0.5,0.5,0.5,1.,1.,1.,1.,1.,1.,0.5,0.5,0.5],'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_NA':[1.,1.,1.,1.,1.,0.5,1.,1.,1.,1.,1.,1.],'NE_NA':[1.,1.,1.,1.,1.,0.5,1.,1.,1.,1.,1.,1.],'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'N_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_EU':[1.,1.,1.,1.,1.,1.,1.,1.,2.,2.,1.,1.],'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_EU':[2.,2.,2.,2.,2.,4.,4.,4.,4.,4.,4.,4.],'E_EU':[2.,2.,1.,1.,1.,0.5,1.,1.,1.,2.,2.,2.],'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]}
#drydepo3_scale = {'CE_NA':[1.,1.,1.,2.,4.,4.,4.,4.,2.,1.,1.,1.],'SE_NA':[1.,1.,1.,2.,2.,4.,4.,4.,4.,4.,1.,1.],'S_NA':[1.,1.,2.,2.,4.,4.,4.,4.,4.,2.,1.,1.],'SW_NA':[1.,2.,2.,2.,2.,2.,4.,4.,2.,2.,1.,1.],'C_NA':[1.,1.,1.,2.,2.,2.,4.,4.,4.,1.,1.,1.],'NE_NA':[1.,1.,2.,2.,2.,1.,4.,4.,4.,2.,1.,1.],'NW_NA':[1.,2.,2.,4.,4.,4.,4.,4.,4.,4.,4.,2.],'N_EU':[1.,1.,1.,1.,2.,2.,2.,2.,2.,1.,1.,1.],'C_EU':[1.,1.,1.,1.,2.,2.,4.,4.,1.,1.,1.,1.],'SW_EU':[2.,2.,2.,4.,4.,4.,4.,4.,4.,4.,4.,2.],'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'E_EU':[1.,1.,1.,1.,1.,1.,2.,2.,2.,1.,1.,1.],'NW_EU':[1.,1.,1.,1.,2.,4.,4.,4.,4.,4.,4.,2.],'NE_AS':[1.,1.,1.,1.,1.,2.,4.,4.,4.,1.,1.,1.]}
#nox_scale = {'CE_NA':[0.5,0.5,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5],'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_NA':[1.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,1.],'SW_NA':[2.,2.,2.,2.,4.,4.,2.,2.,2.,2.,2.,2.],'C_NA':[1.,1.,2.,2.,2.,0.5,1.,1.,1.,2.,1.,1.],'NE_NA':[2.,2.,2.,4.,4.,0.5,4.,4.,2.,2.,2.,2.],'NW_NA':[2.,2.,4.,4.,4.,4.,4.,4.,4.,4.,4.,4.],'N_EU':[2.,2.,2.,1.,0.5,0.5,0.5,0.5,0.5,4.,2.,2.],'C_EU':[1.,1.,2.,2.,2.,2.,4.,4.,2.,1.,1.,1.],'SW_EU':[2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.],'S_EU':[1.,1.,2.,2.,0.25,0.25,0.25,0.25,0.25,2.,2.,2.],'E_EU':[0.5,0.5,1.,0.5,0.5,0.5,0.5,0.5,1.,1.,1.,0.5],'NW_EU':[1.,1.,2.,2.,2.,2.,4.,4.,2.,1.,1.,1.],'NE_AS':[0.5,0.5,0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5]}

#anmvoc_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SE_NA':[0.5,0.5,0.5,1.,1.,1.,1.,1.,1.,1.,1.,0.5],'S_NA':[0.5,0.5,0.5,1.,1.,1.,1.,1.,1.,1.,1.,0.5],'SW_NA':[0.5,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5,0.5],'C_NA':[0.5,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5],'NE_NA':[0.5,0.5,1.,1.,1.,1.,1.,1.,1.,1.,0.5,0.5],'NW_NA':[0.5,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5],'N_EU':[2.,2.,2.,1.,1.,1.,1.,1.,1.,2.,2.,2.],'C_EU':[2.,2.,2.,2.,1.,1.,1.,1.,1.,1.,2.,2.],'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_EU':[0.5,0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5,0.5],'E_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[0.5,0.5,0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5]}
#bnmvoc_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'N_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'E_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]}
#co_scale = {'CE_NA':[0.5,0.5,0.5,1.,1.,1.,1.,1.,1.,1.,0.5,0.5],'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'N_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'E_EU':[2.,2.,1.,1.,1.,1.,1.,1.,1.,2.,2.,1.],'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]}
#drydepo3_scale = {'CE_NA':[1.,1.,1.,2.,2.,4.,4.,4.,2.,2.,1.,1.],'SE_NA':[1.,1.,1.,2.,2.,4.,4.,4.,4.,4.,2.,1.],'S_NA':[1.,1.,1.,2.,2.,4.,4.,4.,4.,2.,2.,1.],'SW_NA':[1.,2.,2.,2.,2.,2.,2.,4.,2.,2.,1.,1.],'C_NA':[1.,1.,1.,1.,2.,2.,2.,2.,2.,1.,1.,1.],'NE_NA':[1.,1.,1.,1.,2.,2.,4.,4.,4.,2.,1.,1.],'NW_NA':[1.,2.,2.,2.,2.,4.,4.,4.,4.,2.,2.,1.],'N_EU':[1.,1.,1.,2.,2.,2.,4.,4.,4.,1.,1.,1.],'C_EU':[1.,1.,1.,1.,2.,2.,2.,2.,2.,1.,1.,1.],'SW_EU':[1.,2.,2.,2.,4.,4.,4.,4.,4.,4.,2.,1.],'S_EU':[1.,1.,1.,1.,2.,2.,2.,2.,2.,1.,1.,1.],'E_EU':[1.,1.,1.,1.,1.,1.,2.,2.,2.,1.,1.,1.],'NW_EU':[1.,1.,1.,1.,2.,4.,4.,4.,4.,2.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,2.,4.,4.,4.,1.,1.,1.]}
#nox_scale = {'CE_NA':[0.5,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5],'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'S_NA':[1.,1.,2.,2.,1.,1.,1.,1.,1.,1.,1.,1.],'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'C_NA':[1.,1.,2.,2.,2.,1.,0.5,0.5,0.5,1.,1.,1.],'NE_NA':[0.5,0.5,0.5,1.,1.,1.,2.,2.,1.,0.5,0.5,0.5],'NW_NA':[2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.],'N_EU':[1.,1.,1.,2.,2.,2.,2.,2.,2.,1.,1.,1.],'C_EU':[0.5,1.,1.,1.,2.,1.,1.,1.,2.,1.,1.,0.5],'SW_EU':[2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.],'S_EU':[0.5,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'E_EU':[0.5,0.5,1.,1.,0.5,0.5,0.5,0.5,1.,1.,1.,0.5],'NW_EU':[1.,1.,1.,2.,2.,2.,2.,2.,2.,2.,2.,2.],'NE_AS':[0.5,0.5,0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5]}

#anmvoc_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'SE_NA':[1.,1.,0.5,0.5,1.,1.,1.,1.,1.,1.,1.,1.],          'S_NA':[0.5,0.5,0.5,0.5,1.,1.,1.,1.,1.,1.,1.,0.5], 'SW_NA':[0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5,0.5,0.5],     'C_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],      'NE_NA':[0.5,0.5,1.,1.,1.,1.,1.,1.,1.,1.,0.5,0.5],  'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],       'N_EU':[2.,2.,2.,2.,1.,1.,1.,1.,1.,1.,2.,2.],            'C_EU':[4.,4.,4.,4.,0.5,1.,1.,1.,1.,0.5,4.,4.],'SW_EU':[0.5,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5],'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5,0.5,1.],      'E_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],    'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[0.5,0.5,0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5]}
#bnmvoc_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'S_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],      'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'C_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],      'NE_NA':[0.5,0.5,1.,1.,1.,1.,1.,1.,1.,1.,0.5,0.5],     'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],       'N_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'C_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],  'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],  'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],        'E_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],    'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]}
#co_scale = {'CE_NA':[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],    'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'S_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],      'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'C_NA':[1.,1.,1.,1.,1.,0.5,1.,1.,1.,1.,1.,1.],     'NE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],      'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],       'N_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'C_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],  'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],  'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],        'E_EU':[1.,1.,2.,2.,2.,1.,1.,1.,2.,2.,2.,1.],    'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]}
#drydepo3_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],          'SE_NA':[2.,2.,1.,1.,2.,2.,4.,4.,4.,2.,2.,2.],            'S_NA':[1.,1.,1.,1.,2.,4.,4.,4.,2.,2.,1.,1.],      'SW_NA':[1.,1.,1.,1.,1.,1.,2.,2.,1.,1.,1.,1.],            'C_NA':[1.,1.,1.,1.,1.,1.,2.,2.,2.,1.,1.,1.],      'NE_NA':[1.,1.,1.,1.,2.,2.,4.,4.,4.,2.,1.,1.],      'NW_NA':[2.,2.,2.,2.,2.,2.,4.,4.,2.,2.,2.,2.],       'N_EU':[1.,1.,1.,1.,1.,2.,2.,2.,2.,1.,1.,1.],            'C_EU':[1.,1.,1.,1.,1.,2.,4.,4.,2.,1.,1.,1.],  'SW_EU':[1.,2.,2.,2.,2.,4.,4.,4.,4.,4.,2.,1.],  'S_EU':[1.,1.,1.,1.,1.,2.,2.,4.,2.,1.,1.,1.],        'E_EU':[1.,1.,1.,1.,1.,2.,2.,2.,1.,1.,1.,1.],    'NW_EU':[1.,1.,1.,1.,2.,4.,4.,4.,4.,4.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,2.,4.,4.,4.,1.,1.,1.]}
#nox_scale = {'CE_NA':[0.5,0.5,0.5,0.5,0.5,0.25,0.25,0.25,0.5,0.5,0.5,0.5],'SE_NA':[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],'S_NA':[0.5,0.5,0.5,1.,1.,0.5,0.5,1.,1.,1.,1.,0.5],'SW_NA':[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],'C_NA':[1.,1.,1.,0.5,0.5,0.5,0.5,0.5,0.5,1.,1.,1.],'NE_NA':[0.5,0.5,0.5,1.,1.,1.,2.,2.,1.,0.5,0.5,0.5],'NW_NA':[0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5,0.5,0.5],'N_EU':[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],'C_EU':[1.,1.,1.,2.,2.,2.,2.,2.,2.,1.,1.,1.],  'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],  'S_EU':[0.5,0.5,0.5,0.5,0.5,1.,2.,2.,1.,0.5,0.5,0.5],'E_EU':[0.5,0.5,1.,1.,2.,2.,2.,2.,2.,1.,0.5,0.5],'NW_EU':[1.,1.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.],'NE_AS':[0.5,0.5,0.5,0.5,0.5,1.,1.,2.,1.,1.,0.5,0.5]}

anmvoc_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],         'SE_NA':[1.,1.,0.5,0.5,1.,1.,1.,1.,1.,1.,1.,1.],          'S_NA':[0.5,0.5,0.5,0.5,1.,1.,1.,1.,1.,1.,1.,0.5], 'SW_NA':[0.5,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5,0.5],'C_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],   'NE_NA':[0.5,0.5,1.,1.,1.,1.,1.,1.,1.,1.,0.5,0.5],  'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],       'N_EU':[2.,2.,2.,2.,1.,1.,1.,1.,1.,1.,2.,2.],            'C_EU':[1.,1.,1.,1.,0.5,1.,1.,1.,1.,1.,1.,1.],   'SW_EU':[0.5,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5],'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,0.5,0.5,1.],      'E_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],    'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[0.5,0.5,0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5]}
bnmvoc_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],         'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'S_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],      'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],   'C_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],   'NE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],      'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],       'N_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'C_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],    'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],  'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],        'E_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],    'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]}
co_scale = {'CE_NA':[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5], 'SE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'S_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],      'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],   'C_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],   'NE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],      'NW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],       'N_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],            'C_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],    'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],  'S_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],        'E_EU':[1.,1.,2.,2.,2.,1.,1.,1.,2.,2.,2.,1.],    'NW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]}
drydepo3_scale = {'CE_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],       'SE_NA':[2.,2.,1.,1.,2.,2.,4.,4.,4.,2.,2.,2.],            'S_NA':[1.,1.,1.,1.,2.,4.,4.,4.,2.,2.,1.,1.],      'SW_NA':[1.,2.,2.,2.,2.,2.,2.,4.,2.,2.,1.,1.],   'C_NA':[1.,1.,1.,1.,1.,2.,2.,2.,2.,1.,1.,1.],   'NE_NA':[1.,1.,1.,1.,2.,2.,4.,4.,4.,2.,1.,1.],      'NW_NA':[2.,2.,2.,2.,2.,2.,4.,4.,2.,2.,2.,2.],       'N_EU':[1.,1.,1.,1.,1.,2.,2.,2.,2.,1.,1.,1.],            'C_EU':[1.,1.,1.,1.,1.,2.,4.,4.,2.,1.,1.,1.],    'SW_EU':[1.,2.,2.,2.,2.,4.,4.,4.,4.,4.,2.,1.],  'S_EU':[1.,1.,1.,1.,1.,2.,2.,4.,2.,1.,1.,1.],        'E_EU':[1.,1.,1.,1.,1.,2.,2.,2.,1.,1.,1.,1.],    'NW_EU':[1.,1.,1.,1.,2.,4.,4.,4.,4.,2.,1.,1.],'NE_AS':[1.,1.,1.,1.,1.,2.,4.,4.,4.,1.,1.,1.]}
nox_scale = {'CE_NA':[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],'SE_NA':[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],'S_NA':[0.5,0.5,0.5,1.,1.,0.5,0.5,1.,1.,1.,1.,0.5],'SW_NA':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],   'C_NA':[1.,1.,1.,1.,1.,1.,0.5,0.5,0.5,1.,1.,1.],'NE_NA':[0.5,0.5,0.5,1.,1.,1.,2.,2.,1.,0.5,0.5,0.5],'NW_NA':[0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5,0.5,0.5],'N_EU':[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],'C_EU':[0.5,0.5,1.,1.,2.,2.,2.,2.,2.,1.,0.5,0.5],'SW_EU':[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],  'S_EU':[0.5,0.5,0.5,0.5,0.5,1.,1.,1.,1.,0.5,0.5,0.5],'E_EU':[0.5,0.5,1.,1.,2.,2.,2.,2.,2.,1.,0.5,0.5],'NW_EU':[1.,1.,1.,1.,2.,2.,2.,2.,2.,1.,2.,2.],'NE_AS':[0.5,0.5,0.5,0.5,0.5,1.,1.,1.,1.,1.,0.5,0.5]}

#loop through scale params, creating netcdfs

def output_scale_netcdf(lons_c,lats_c,lons_e,lats_e,areas,area_boundaries,scale,fname):

	#setup netcdfs
	root_grp = Dataset(fname,'w')

	root_grp.createDimension('time', 12)
	root_grp.createDimension('lon', len(lons_c))
	root_grp.createDimension('lat', len(lats_c))

	t = root_grp.createVariable('time','f8', ('time',))
	lon = root_grp.createVariable('lon','f8', ('lon',))
	lat = root_grp.createVariable('lat','f8', ('lat',))
	scale_v = root_grp.createVariable('scale','f8', ('time','lat','lon',))
	
	t.standard_name = 'time'
	t.units = 'hours since 1985-01-01 00:00:00'
	t.calendar = 'standard'
	
	lon.standard_name = 'longitude'
	lon.long_name = 'longitude'
	lon.units = 'degrees_east'
	lon.axis = 'X'
    
	lat.standard_name = 'latitude'
	lat.long_name = 'latitude'
	lat.units = 'degrees_north'
	lat.axis = 'Y'
	
	scale_v.long_name = 'Scaling ratios'
	scale_v.units = 'unitless'
	scale_v.grid_type = 'gaussian'
	scale_v.gamap_category = 'RATIO_2D'
	
	t[:] = [175320.,176064.,176736.,177480.,178200.,178944.,179664.,180408.,181152.,181872.,182616.,183336.]
	lon[:] = lons_c
	lat[:] = lats_c

	#setup 1x1 deg. grid
	grid = np.ones((12,len(lats_c),len(lons_c)))

	inter_grid = np.empty((len(areas),len(lats_c),len(lons_c)))

	#loop through areas 
	for a in range(len(areas)):
		area = areas[a]
	
		#get boundaries for area
		boundaries = area_boundaries[area]
	
		lat_test = [[]*len(areas)]
		lon_test = [[]*len(areas)]
		
		o = []
		
		XA1 = boundaries[2]
		XA2 = boundaries[3]
		YA1 = boundaries[0]
		YA2 = boundaries[1]
		
		#
	
		for lat_i in range(len(lats_c)):
			for lon_i in range(len(lons_c)):

				#A is area, B is gridbox
				XB1 = lons_e[lon_i]
				XB2 = lons_e[lon_i+1]
				YB1 = lats_e[lat_i]
				YB2 = lats_e[lat_i+1]
				
				Aheight = YA2-YA1
				Bheight = YB2-YB1
				Awidth = XA2-XA1
				Bwidth = XB2-XB1

				#how much of gridbox is covered by area
				#B all in A
				if (XB1 >= XA1) & (XB2 <= XA2):
					overlapfracx = 1. 
				#PARTIAL B IN A
				elif (XB1 >= XA1) & (XB1 < XA2) & (XB2 >= XA2):
					overlapx = XA2 - XB1    
					overlapfracx = overlapx/Bwidth
				elif (XB1 < XA1) & (XB2 > XA1):
					overlapx = XB2 - XA1    
					overlapfracx = overlapx/Bwidth
				else:
					overlapfracx = 0. 

				#B all in A
				if (YB1 >= YA1) & (YB2 <= YA2):
					overlapfracy = 1. 
				#PARTIAL B IN A
				elif (YB1 >= YA1) & (YB1 < YA2) & (YB2 >= YA2):
					overlapy = YA2 - YB1    
					overlapfracy = overlapy/Bheight

				elif (YB1 < YA1) & (YB2 > YA1):
					overlapy = YB2 - YA1    
					overlapfracy = overlapy/Bheight
				else:
					overlapfracy = 0.

				overlap_ratio = overlapfracx*overlapfracy
				
				inter_grid[a,lat_i,lon_i] = overlap_ratio
				
				#if overlap_ratio >= 0.5:
				#lat_test[a].append(lat_i)
				#lon_test[a].append(lon_i)
			
		#lon_test = np.array(lon_test)
		#lat_test = np.array(lat_test)
	
		#get area scaling lists
		#sca = scale[area]

		#loop through months
		#for m in range(12):
		
			#overwrite grid scaling in different areas
			#grid[m,lat_test,lon_test] = sca[m]

	#get index of area with max ratio in lat,lon pairing
	#use this to avoid scaling 2 areas simultaneously
	
	area_max = np.argmax(inter_grid,axis=0)
	
	for a in range(len(areas)):
		area = areas[a]
	
		#get area scaling lists
		sca = scale[area]
	
		for lat_i in range(len(lats_c)):
			for lon_i in range(len(lons_c)):
				if a == area_max[lat_i,lon_i]:
					if inter_grid[a,lat_i,lon_i] > 0:
			
						#loop through months
						for m in range(12):
		
							#overwrite grid scaling in different areas
							grid[m,lat_i,lon_i] = sca[m]
	
	scale_v[:] = grid 
		
	root_grp.close()
	
#process netcdfs 

#0.1X0.1
output_scale_netcdf(lon_c_01,lat_c_01,lon_e_01,lat_e_01,areas,area_boundaries,anmvoc_scale,'ANMVOC_SCALING_0.1x0.1%s.nc'%(vers))
output_scale_netcdf(lon_c_01,lat_c_01,lon_e_01,lat_e_01,areas,area_boundaries,bnmvoc_scale,'BNMVOC_SCALING_0.1x0.1%s.nc'%(vers))
output_scale_netcdf(lon_c_01,lat_c_01,lon_e_01,lat_e_01,areas,area_boundaries,co_scale,'CO_SCALING_0.1x0.1%s.nc'%(vers))
output_scale_netcdf(lon_c_01,lat_c_01,lon_e_01,lat_e_01,areas,area_boundaries,drydepo3_scale,'DRYDEPO3_SCALING_0.1x0.1%s.nc'%(vers))
output_scale_netcdf(lon_c_01,lat_c_01,lon_e_01,lat_e_01,areas,area_boundaries,nox_scale,'NOX_SCALING_0.1x0.1%s.nc'%(vers))

#0.25X0.25
output_scale_netcdf(lon_c_025,lat_c_025,lon_e_025,lat_e_025,areas,area_boundaries,anmvoc_scale,'ANMVOC_SCALING_0.25x0.25%s.nc'%(vers))
output_scale_netcdf(lon_c_025,lat_c_025,lon_e_025,lat_e_025,areas,area_boundaries,bnmvoc_scale,'BNMVOC_SCALING_0.25x0.25%s.nc'%(vers))
output_scale_netcdf(lon_c_025,lat_c_025,lon_e_025,lat_e_025,areas,area_boundaries,co_scale,'CO_SCALING_0.25x0.25%s.nc'%(vers))
output_scale_netcdf(lon_c_025,lat_c_025,lon_e_025,lat_e_025,areas,area_boundaries,drydepo3_scale,'DRYDEPO3_SCALING_0.25x0.25%s.nc'%(vers))
output_scale_netcdf(lon_c_025,lat_c_025,lon_e_025,lat_e_025,areas,area_boundaries,nox_scale,'NOX_SCALING_0.25x0.25%s.nc'%(vers))

#1X1
output_scale_netcdf(lon_c_1x1,lat_c_1x1,lon_e_1x1,lat_e_1x1,areas,area_boundaries,anmvoc_scale,'ANMVOC_SCALING_1x1%s.nc'%(vers))
output_scale_netcdf(lon_c_1x1,lat_c_1x1,lon_e_1x1,lat_e_1x1,areas,area_boundaries,bnmvoc_scale,'BNMVOC_SCALING_1x1%s.nc'%(vers))
output_scale_netcdf(lon_c_1x1,lat_c_1x1,lon_e_1x1,lat_e_1x1,areas,area_boundaries,co_scale,'CO_SCALING_1x1%s.nc'%(vers))
output_scale_netcdf(lon_c_1x1,lat_c_1x1,lon_e_1x1,lat_e_1x1,areas,area_boundaries,drydepo3_scale,'DRYDEPO3_SCALING_1x1%s.nc'%(vers))
output_scale_netcdf(lon_c_1x1,lat_c_1x1,lon_e_1x1,lat_e_1x1,areas,area_boundaries,nox_scale,'NOX_SCALING_1x1%s.nc'%(vers))

#4X5
output_scale_netcdf(lon_c_4x5,lat_c_4x5,lon_e_4x5,lat_e_4x5,areas,area_boundaries,anmvoc_scale,'ANMVOC_SCALING_4x5%s.nc'%(vers))
output_scale_netcdf(lon_c_4x5,lat_c_4x5,lon_e_4x5,lat_e_4x5,areas,area_boundaries,bnmvoc_scale,'BNMVOC_SCALING_4x5%s.nc'%(vers))
output_scale_netcdf(lon_c_4x5,lat_c_4x5,lon_e_4x5,lat_e_4x5,areas,area_boundaries,co_scale,'CO_SCALING_4x5%s.nc'%(vers))
output_scale_netcdf(lon_c_4x5,lat_c_4x5,lon_e_4x5,lat_e_4x5,areas,area_boundaries,drydepo3_scale,'DRYDEPO3_SCALING_4x5%s.nc'%(vers))
output_scale_netcdf(lon_c_4x5,lat_c_4x5,lon_e_4x5,lat_e_4x5,areas,area_boundaries,nox_scale,'NOX_SCALING_4x5%s.nc'%(vers))