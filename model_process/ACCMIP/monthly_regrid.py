import modules
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
import datetime
import pandas as pd

data_type = raw_input('amp, ph or ave?\n')

if data_type != 'ave':
    period = 's'
else:
    period = ''

#read in data
#----------------------------------------------
cesm_cam_root = Dataset('../CESMCAM_SURFACE_*_*_*_M_*/model_sig_periods.nc')
cmam_root = Dataset('../CMAM_SURFACE_*_*_*_M_*/model_sig_periods.nc')
emac_root = Dataset('../EMAC_SURFACE_*_*_*_M_*/model_sig_periods.nc')
geos_ccm_root = Dataset('../GEOSCCM_SURFACE_*_*_*_M_*/model_sig_periods.nc')
geos_chem_root = Dataset('../GEOSCHEM_SURFACE_v90103_2x2.5_GEOS5_M_*/model_sig_periods.nc')
gfdl_root = Dataset('../GFDLAM3_SURFACE_*_*_*_M_*/model_sig_periods.nc')
giss_root = Dataset('../GISSE2R_SURFACE_*_*_*_M_*/model_sig_periods.nc')
hadgem2_root = Dataset('../HADGEM2_SURFACE_*_*_*_M_*/model_sig_periods.nc')
miroc_chem_root = Dataset('../MIROCCHEM_SURFACE_*_*_*_M_*/model_sig_periods.nc')
ncar_root = Dataset('../NCARCAM3.5_SURFACE_*_*_*_M_*/model_sig_periods.nc')
stochad_root = Dataset('../STOCHADAM3_SURFACE_*_*_*_M_*/model_sig_periods.nc')


cesm_cam_sa = cesm_cam_root.variables['seasonal_amplitude'][:]
cesm_cam_sp = cesm_cam_root.variables['seasonal_phase'][:]
cesm_cam_ave = cesm_cam_root.variables['average'][:]
cesm_cam_lat_e = cesm_cam_root.variables['lat_edges'][:]
cesm_cam_lon_e = cesm_cam_root.variables['lon_edges'][:]
cesm_cam_lat_c = cesm_cam_root.variables['lat_centre'][:]
cesm_cam_lon_c = cesm_cam_root.variables['lon_centre'][:]
cesm_cam_n_boxes = len(cesm_cam_lon_c)*len(cesm_cam_lat_c)

cmam_sa = cmam_root.variables['seasonal_amplitude'][:]
cmam_sp = cmam_root.variables['seasonal_phase'][:]
cmam_ave = cmam_root.variables['average'][:]
cmam_lat_e = cmam_root.variables['lat_edges'][:]
cmam_lon_e = cmam_root.variables['lon_edges'][:]
cmam_lat_c = cmam_root.variables['lat_centre'][:]
cmam_lon_c = cmam_root.variables['lon_centre'][:]
cmam_n_boxes = len(cmam_lon_c)*len(cmam_lat_c)

emac_sa = emac_root.variables['seasonal_amplitude'][:]
emac_sp = emac_root.variables['seasonal_phase'][:]
emac_ave = emac_root.variables['average'][:]
emac_lat_e = emac_root.variables['lat_edges'][:]
emac_lon_e = emac_root.variables['lon_edges'][:]
emac_lat_c = emac_root.variables['lat_centre'][:]
emac_lon_c = emac_root.variables['lon_centre'][:]
emac_n_boxes = len(emac_lon_c)*len(emac_lat_c)

geos_ccm_sa = geos_ccm_root.variables['seasonal_amplitude'][:]
geos_ccm_sp = geos_ccm_root.variables['seasonal_phase'][:]
geos_ccm_ave = geos_ccm_root.variables['average'][:]
geos_ccm_lat_e = geos_ccm_root.variables['lat_edges'][:]
geos_ccm_lon_e = geos_ccm_root.variables['lon_edges'][:]
geos_ccm_lat_c = geos_ccm_root.variables['lat_centre'][:]
geos_ccm_lon_c = geos_ccm_root.variables['lon_centre'][:]
geos_ccm_n_boxes = len(geos_ccm_lon_c)*len(geos_ccm_lat_c)

geos_chem_sa = geos_chem_root.variables['seasonal_amplitude'][:]
geos_chem_sp = geos_chem_root.variables['seasonal_phase'][:]
geos_chem_ave = geos_chem_root.variables['average'][:]
geos_chem_lat_e = geos_chem_root.variables['lat_edges'][:]
geos_chem_lon_e = geos_chem_root.variables['lon_edges'][:]
geos_chem_lat_c = geos_chem_root.variables['lat_centre'][:]
geos_chem_lon_c = geos_chem_root.variables['lon_centre'][:]
geos_chem_n_boxes = len(geos_chem_lon_c)*len(geos_chem_lat_c)

gfdl_sa = gfdl_root.variables['seasonal_amplitude'][:]
gfdl_sp = gfdl_root.variables['seasonal_phase'][:]
gfdl_ave = gfdl_root.variables['average'][:]
gfdl_lat_e = gfdl_root.variables['lat_edges'][:]
gfdl_lon_e = gfdl_root.variables['lon_edges'][:]
gfdl_lat_c = gfdl_root.variables['lat_centre'][:]
gfdl_lon_c = gfdl_root.variables['lon_centre'][:]
gfdl_n_boxes = len(gfdl_lon_c)*len(gfdl_lat_c)

giss_sa = giss_root.variables['seasonal_amplitude'][:]
giss_sp = giss_root.variables['seasonal_phase'][:]
giss_ave = giss_root.variables['average'][:]
giss_lat_e = giss_root.variables['lat_edges'][:]
giss_lon_e = giss_root.variables['lon_edges'][:]
giss_lat_c = giss_root.variables['lat_centre'][:]
giss_lon_c = giss_root.variables['lon_centre'][:]
giss_n_boxes = len(giss_lon_c)*len(giss_lat_c)

hadgem2_sa = hadgem2_root.variables['seasonal_amplitude'][:]
hadgem2_sp = hadgem2_root.variables['seasonal_phase'][:]
hadgem2_ave = hadgem2_root.variables['average'][:]
hadgem2_lat_e = hadgem2_root.variables['lat_edges'][:]
hadgem2_lon_e = hadgem2_root.variables['lon_edges'][:]
hadgem2_lat_c = hadgem2_root.variables['lat_centre'][:]
hadgem2_lon_c = hadgem2_root.variables['lon_centre'][:]
hadgem2_n_boxes = len(hadgem2_lon_c)*len(hadgem2_lat_c)

miroc_chem_sa = miroc_chem_root.variables['seasonal_amplitude'][:]
miroc_chem_sp = miroc_chem_root.variables['seasonal_phase'][:]
miroc_chem_ave = miroc_chem_root.variables['average'][:]
miroc_chem_lat_e = miroc_chem_root.variables['lat_edges'][:]
miroc_chem_lon_e = miroc_chem_root.variables['lon_edges'][:]
miroc_chem_lat_c = miroc_chem_root.variables['lat_centre'][:]
miroc_chem_lon_c = miroc_chem_root.variables['lon_centre'][:]
miroc_chem_n_boxes = len(miroc_chem_lon_c)*len(miroc_chem_lat_c)

ncar_sa = ncar_root.variables['seasonal_amplitude'][:]
ncar_sp = ncar_root.variables['seasonal_phase'][:]
ncar_ave = ncar_root.variables['average'][:]
ncar_lat_e = ncar_root.variables['lat_edges'][:]
ncar_lon_e = ncar_root.variables['lon_edges'][:]
ncar_lat_c = ncar_root.variables['lat_centre'][:]
ncar_lon_c = ncar_root.variables['lon_centre'][:]
ncar_n_boxes = len(ncar_lon_c)*len(ncar_lat_c)

stochad_sa = stochad_root.variables['seasonal_amplitude'][:]
stochad_sp = stochad_root.variables['seasonal_phase'][:]
stochad_ave = stochad_root.variables['average'][:]
stochad_lat_e = stochad_root.variables['lat_edges'][:]
stochad_lon_e = stochad_root.variables['lon_edges'][:]
stochad_lat_c = stochad_root.variables['lat_centre'][:]
stochad_lon_c = stochad_root.variables['lon_centre'][:]
stochad_n_boxes = len(stochad_lon_c)*len(stochad_lat_c)

if data_type == 'amp':
    if period == 's':
        cesm_cam_param = cesm_cam_sa    
        cmam_param = cmam_sa
        emac_param = emac_sa
        geos_ccm_param = geos_ccm_sa
        geos_chem_param = geos_chem_sa
        gfdl_param = gfdl_sa
        giss_param = giss_sa
        hadgem2_param = hadgem2_sa
        miroc_chem_param = miroc_chem_sa
        ncar_param = ncar_sa
        stochad_param = stochad_sa
        max_ph = 'na'

elif data_type == 'ph':
    if period == 's': 
        cesm_cam_param = cesm_cam_sp    
        cmam_param = cmam_sp
        emac_param = emac_sp
        geos_ccm_param = geos_ccm_sp
        geos_chem_param = geos_chem_sp
        gfdl_param = gfdl_sp
        giss_param = giss_sp
        hadgem2_param = hadgem2_sp
        miroc_chem_param = miroc_chem_sp
        ncar_param = ncar_sp
        stochad_param = stochad_sp
        max_ph = 12

elif data_type == 'ave':
    cesm_cam_param = cesm_cam_ave    
    cmam_param = cmam_ave
    emac_param = emac_ave
    geos_ccm_param = geos_ccm_ave
    geos_chem_param = geos_chem_ave
    gfdl_param = gfdl_ave
    giss_param = giss_ave
    hadgem2_param = hadgem2_ave
    miroc_chem_param = miroc_chem_ave
    ncar_param = ncar_ave
    stochad_param = stochad_ave
    max_ph = 'na'

#interpolate model grids onto minimum grid size
#----------------------------------------------
cesm_cam_lat_diff =  180./len(cesm_cam_lat_c)
cesm_cam_lon_diff = 360./len(cesm_cam_lon_c)

cmam_lat_diff =  180./len(cmam_lat_c)
cmam_lon_diff =  360./len(cmam_lon_c)

emac_lat_diff =  180./len(emac_lat_c)
emac_lon_diff =  360./len(emac_lon_c)

geos_ccm_lat_diff =  180./len(geos_ccm_lat_c)
geos_ccm_lon_diff =  360./len(geos_ccm_lon_c)

geos_chem_lat_diff =  180./len(geos_chem_lat_c)
geos_chem_lon_diff =  360./len(geos_chem_lon_c)

gfdl_lat_diff =  180./len(gfdl_lat_c)
gfdl_lon_diff =  360./len(gfdl_lon_c)

giss_lat_diff =  180./len(giss_lat_c)
giss_lon_diff =  360./len(giss_lon_c)

hadgem2_lat_diff =  180./len(hadgem2_lat_c)
hadgem2_lon_diff =  360./len(hadgem2_lon_c)

miroc_chem_lat_diff =  180./len(miroc_chem_lat_c)
miroc_chem_lon_diff =  360./len(miroc_chem_lon_c)

ncar_lat_diff =  180./len(ncar_lat_c)
ncar_lon_diff =  360./len(ncar_lon_c)

stochad_lat_diff =  180./len(stochad_lat_c)
stochad_lon_diff =  360./len(stochad_lon_c)

#min_lat = np.min([cesm_cam_lat_diff,cmam_lat_diff,emac_lat_diff,geos_ccm_lat_diff,geos_chem_lat_diff,gfdl_lat_diff,giss_lat_diff,hadgem2_lat_diff,miroc_chem_lat_diff,ncar_lat_diff,stochad_lat_diff])
#min_lon = np.min([cesm_cam_lon_diff,cmam_lon_diff,emac_lon_diff,geos_ccm_lon_diff,geos_chem_lon_diff,gfdl_lon_diff,giss_lon_diff,hadgem2_lon_diff,miroc_chem_lon_diff,ncar_lon_diff,stochad_lon_diff])

min_lat = 1.
min_lon = 1.

new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_cesm_cam = modules.regrid_interp(cesm_cam_param,cesm_cam_lat_e,cesm_cam_lat_c,cesm_cam_lon_e,cesm_cam_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_cmam = modules.regrid_interp(cmam_param,cmam_lat_e,cmam_lat_c,cmam_lon_e,cmam_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_emac = modules.regrid_interp(emac_param,emac_lat_e,emac_lat_c,emac_lon_e,emac_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_geos_ccm = modules.regrid_interp(geos_ccm_param,geos_ccm_lat_e,geos_ccm_lat_c,geos_ccm_lon_e,geos_ccm_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_geos_chem = modules.regrid_interp(geos_chem_param,geos_chem_lat_e,geos_chem_lat_c,geos_chem_lon_e,geos_chem_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_gfdl = modules.regrid_interp(gfdl_param,gfdl_lat_e,gfdl_lat_c,gfdl_lon_e,gfdl_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_giss = modules.regrid_interp(giss_param,giss_lat_e,giss_lat_c,giss_lon_e,giss_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_hadgem2 = modules.regrid_interp(hadgem2_param,hadgem2_lat_e,hadgem2_lat_c,hadgem2_lon_e,hadgem2_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_miroc_chem = modules.regrid_interp(miroc_chem_param,miroc_chem_lat_e,miroc_chem_lat_c,miroc_chem_lon_e,miroc_chem_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_ncar = modules.regrid_interp(ncar_param,ncar_lat_e,ncar_lat_c,ncar_lon_e,ncar_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_stochad = modules.regrid_interp(stochad_param,stochad_lat_e,stochad_lat_c,stochad_lon_e,stochad_lon_c,min_lat,min_lon,data_type,max_ph)

if (data_type == 'amp') or (data_type == 'ave'):
    z_std = np.std((z_cesm_cam,z_cmam,z_emac,z_geos_ccm,z_geos_chem,z_gfdl,z_giss,z_hadgem2,z_miroc_chem,z_ncar,z_stochad),axis=0)
    z_ave = np.average((z_cesm_cam,z_cmam,z_emac,z_geos_ccm,z_geos_chem,z_gfdl,z_giss,z_hadgem2,z_miroc_chem,z_ncar,z_stochad),axis=0)
    z_pc = (z_std/z_ave)*100.
if data_type == 'ph':
    data = np.array([z_cesm_cam,z_cmam,z_emac,z_geos_ccm,z_geos_chem,z_gfdl,z_giss,z_hadgem2,z_miroc_chem,z_ncar,z_stochad])
    z_ave,z_std = modules.yamartino_ave_std(data,max_ph)
    z_pc = (z_std/z_ave)*100.

#save out sig periods to netcdf
root_grp_period = Dataset('ACCMIP_monthly_%s%s_std.nc'%(period,data_type), 'w')
root_grp_period.createDimension('lat_centre', len(new_lat_c))
root_grp_period.createDimension('lon_centre', len(new_lon_c))
root_grp_period.createDimension('lat_edges', len(new_lat_e))
root_grp_period.createDimension('lon_edges', len(new_lon_e))
cesmcam = root_grp_period.createVariable('cesmcam', 'f8', ('lat_centre','lon_centre'))
cmam = root_grp_period.createVariable('cmam', 'f8', ('lat_centre','lon_centre'))
emac = root_grp_period.createVariable('emac', 'f8', ('lat_centre','lon_centre'))  
geosccm = root_grp_period.createVariable('geosccm', 'f8', ('lat_centre','lon_centre'))
geoschem = root_grp_period.createVariable('geoschem', 'f8', ('lat_centre','lon_centre'))
gfdl = root_grp_period.createVariable('gfdl', 'f8', ('lat_centre','lon_centre'))
giss = root_grp_period.createVariable('giss', 'f8', ('lat_centre','lon_centre'))
hadgem2 = root_grp_period.createVariable('hadgem2', 'f8', ('lat_centre','lon_centre'))
mirocchem = root_grp_period.createVariable('mirocchem', 'f8', ('lat_centre','lon_centre'))
ncar = root_grp_period.createVariable('ncar', 'f8', ('lat_centre','lon_centre'))
stochad = root_grp_period.createVariable('stochadam3', 'f8', ('lat_centre','lon_centre'))

ave = root_grp_period.createVariable('average', 'f8', ('lat_centre','lon_centre'))
abs = root_grp_period.createVariable('absolute_std', 'f8', ('lat_centre','lon_centre'))                                                                                                                                                       
frac = root_grp_period.createVariable('fractional_std', 'f8', ('lat_centre','lon_centre'))
lat_centre = root_grp_period.createVariable('lat_centre', 'f8', ('lat_centre',))
lon_centre = root_grp_period.createVariable('lon_centre', 'f8', ('lon_centre',))
lat_edge = root_grp_period.createVariable('lat_edges', 'f8', ('lat_edges',))
lon_edge = root_grp_period.createVariable('lon_edges', 'f8', ('lon_edges',))
lat_centre[:] = new_lat_c
lon_centre[:] = new_lon_c
lat_edge[:] = new_lat_e     
lon_edge[:] = new_lon_e 
ave[:] = z_ave
abs[:] = z_std
frac[:] = z_pc
cesmcam[:]  = z_cesm_cam
cmam[:]  = z_cmam
emac[:] = z_emac
geosccm[:]  = z_geos_ccm
geoschem[:]  = z_geos_chem
gfdl[:]  = z_gfdl
giss[:]  = z_giss
hadgem2[:] = z_hadgem2
mirocchem[:]  = z_miroc_chem 
ncar[:]  = z_ncar
stochad[:] = z_stochad
