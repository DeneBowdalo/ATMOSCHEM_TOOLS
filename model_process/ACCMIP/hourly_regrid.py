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
    period = raw_input('\nd or s?\n')
else:
    period = ''

#read in data
#----------------------------------------------
cesm_cam_root = Dataset('../CESMCAM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
cmam_root = Dataset('../CMAM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
geos_ccm_root = Dataset('../GEOSCCM_SURFACE_*_*_*_H_*/model_sig_periods.nc')
geos_chem_root = Dataset('../GEOSCHEM_SURFACE_v902_2x2.5_GEOS5_H_*/model_sig_periods.nc')
gfdl_root = Dataset('../GFDLAM3_SURFACE_*_*_*_H_*/model_sig_periods.nc')
giss_root = Dataset('../GISSE2R_SURFACE_*_*_*_H_*/model_sig_periods.nc')
miroc_chem_root = Dataset('../MIROCCHEM_SURFACE_*_*_*_H_*/model_sig_periods.nc')

cesm_cam_sa = cesm_cam_root.variables['seasonal_amplitude'][:]
cesm_cam_sp = cesm_cam_root.variables['seasonal_phase'][:]
cesm_cam_da = cesm_cam_root.variables['daily_amplitude'][:]
cesm_cam_dp = cesm_cam_root.variables['daily_phase'][:]
cesm_cam_ave = cesm_cam_root.variables['average'][:]
cesm_cam_lat_e = cesm_cam_root.variables['lat_edges'][:]
cesm_cam_lon_e = cesm_cam_root.variables['lon_edges'][:]
cesm_cam_lat_c = cesm_cam_root.variables['lat_centre'][:]
cesm_cam_lon_c = cesm_cam_root.variables['lon_centre'][:]
cesm_cam_n_boxes = len(cesm_cam_lon_c)*len(cesm_cam_lat_c)

cmam_sa = cmam_root.variables['seasonal_amplitude'][:]
cmam_sp = cmam_root.variables['seasonal_phase'][:]
cmam_da = cmam_root.variables['daily_amplitude'][:]
cmam_dp = cmam_root.variables['daily_phase'][:]
cmam_ave = cmam_root.variables['average'][:]
cmam_lat_e = cmam_root.variables['lat_edges'][:]
cmam_lon_e = cmam_root.variables['lon_edges'][:]
cmam_lat_c = cmam_root.variables['lat_centre'][:]
cmam_lon_c = cmam_root.variables['lon_centre'][:]
cmam_n_boxes = len(cmam_lon_c)*len(cmam_lat_c)

geos_ccm_sa = geos_ccm_root.variables['seasonal_amplitude'][:]
geos_ccm_sp = geos_ccm_root.variables['seasonal_phase'][:]
geos_ccm_da = geos_ccm_root.variables['daily_amplitude'][:]
geos_ccm_dp = geos_ccm_root.variables['daily_phase'][:]
geos_ccm_ave = geos_ccm_root.variables['average'][:]
geos_ccm_lat_e = geos_ccm_root.variables['lat_edges'][:]
geos_ccm_lon_e = geos_ccm_root.variables['lon_edges'][:]
geos_ccm_lat_c = geos_ccm_root.variables['lat_centre'][:]
geos_ccm_lon_c = geos_ccm_root.variables['lon_centre'][:]
geos_ccm_n_boxes = len(geos_ccm_lon_c)*len(geos_ccm_lat_c)

geos_chem_sa = geos_chem_root.variables['seasonal_amplitude'][:]
geos_chem_sp = geos_chem_root.variables['seasonal_phase'][:]
geos_chem_da = geos_chem_root.variables['daily_amplitude'][:]
geos_chem_dp = geos_chem_root.variables['daily_phase'][:]
geos_chem_ave = geos_chem_root.variables['average'][:]
geos_chem_lat_e = geos_chem_root.variables['lat_edges'][:]
geos_chem_lon_e = geos_chem_root.variables['lon_edges'][:]
geos_chem_lat_c = geos_chem_root.variables['lat_centre'][:]
geos_chem_lon_c = geos_chem_root.variables['lon_centre'][:]
geos_chem_n_boxes = len(geos_chem_lon_c)*len(geos_chem_lat_c)

gfdl_sa = gfdl_root.variables['seasonal_amplitude'][:]
gfdl_sp = gfdl_root.variables['seasonal_phase'][:]
gfdl_da = gfdl_root.variables['daily_amplitude'][:]
gfdl_dp = gfdl_root.variables['daily_phase'][:]
gfdl_ave = gfdl_root.variables['average'][:]
gfdl_lat_e = gfdl_root.variables['lat_edges'][:]
gfdl_lon_e = gfdl_root.variables['lon_edges'][:]
gfdl_lat_c = gfdl_root.variables['lat_centre'][:]
gfdl_lon_c = gfdl_root.variables['lon_centre'][:]
gfdl_n_boxes = len(gfdl_lon_c)*len(gfdl_lat_c)

giss_sa = giss_root.variables['seasonal_amplitude'][:]
giss_sp = giss_root.variables['seasonal_phase'][:]
giss_da = giss_root.variables['daily_amplitude'][:]
giss_dp = giss_root.variables['daily_phase'][:]
giss_ave = giss_root.variables['average'][:]
giss_lat_e = giss_root.variables['lat_edges'][:]
giss_lon_e = giss_root.variables['lon_edges'][:]
giss_lat_c = giss_root.variables['lat_centre'][:]
giss_lon_c = giss_root.variables['lon_centre'][:]
giss_n_boxes = len(giss_lon_c)*len(giss_lat_c)

miroc_chem_sa = miroc_chem_root.variables['seasonal_amplitude'][:]
miroc_chem_sp = miroc_chem_root.variables['seasonal_phase'][:]
miroc_chem_da = miroc_chem_root.variables['daily_amplitude'][:]
miroc_chem_dp = miroc_chem_root.variables['daily_phase'][:]
miroc_chem_ave = miroc_chem_root.variables['average'][:]
miroc_chem_lat_e = miroc_chem_root.variables['lat_edges'][:]
miroc_chem_lon_e = miroc_chem_root.variables['lon_edges'][:]
miroc_chem_lat_c = miroc_chem_root.variables['lat_centre'][:]
miroc_chem_lon_c = miroc_chem_root.variables['lon_centre'][:]
miroc_chem_n_boxes = len(miroc_chem_lon_c)*len(miroc_chem_lat_c)

if data_type == 'amp':
    if period == 'd':
        cesm_cam_param = cesm_cam_da    
        cmam_param = cmam_da
        geos_ccm_param = geos_ccm_da
        geos_chem_param = geos_chem_da
        gfdl_param = gfdl_da
        giss_param = giss_da
        miroc_chem_param = miroc_chem_da
        max_ph = 'na'
    elif period == 's':
        cesm_cam_param = cesm_cam_sa    
        cmam_param = cmam_sa
        geos_ccm_param = geos_ccm_sa
        geos_chem_param = geos_chem_sa
        gfdl_param = gfdl_sa
        giss_param = giss_sa
        miroc_chem_param = miroc_chem_sa
        max_ph = 'na'

elif data_type == 'ph':
    if period == 'd':  
        cesm_cam_param = cesm_cam_dp    
        cmam_param = cmam_dp
        geos_ccm_param = geos_ccm_dp
        geos_chem_param = geos_chem_dp
        gfdl_param = gfdl_dp
        giss_param = giss_dp
        miroc_chem_param = miroc_chem_dp
        max_ph = 24
    if period == 's':  
        cesm_cam_param = cesm_cam_sp    
        cmam_param = cmam_sp
        geos_ccm_param = geos_ccm_sp
        geos_chem_param = geos_chem_sp
        gfdl_param = gfdl_sp
        giss_param = giss_sp
        miroc_chem_param = miroc_chem_sp
        max_ph = 12    

elif data_type == 'ave':
    cesm_cam_param = cesm_cam_ave    
    cmam_param = cmam_ave
    geos_ccm_param = geos_ccm_ave
    geos_chem_param = geos_chem_ave
    gfdl_param = gfdl_ave
    giss_param = giss_ave
    miroc_chem_param = miroc_chem_ave
    max_ph = 'na'

#interpolate model grids onto minimum grid size
#----------------------------------------------
cesm_cam_lat_diff =  180./len(cesm_cam_lat_c)
cesm_cam_lon_diff = 360./len(cesm_cam_lon_c)

cmam_lat_diff =  180./len(cmam_lat_c)
cmam_lon_diff =  360./len(cmam_lon_c)

geos_ccm_lat_diff =  180./len(geos_ccm_lat_c)
geos_ccm_lon_diff =  360./len(geos_ccm_lon_c)

geos_chem_lat_diff =  180./len(geos_chem_lat_c)
geos_chem_lon_diff =  360./len(geos_chem_lon_c)

gfdl_lat_diff =  180./len(gfdl_lat_c)
gfdl_lon_diff =  360./len(gfdl_lon_c)

giss_lat_diff =  180./len(giss_lat_c)
giss_lon_diff =  360./len(giss_lon_c)

miroc_chem_lat_diff =  180./len(miroc_chem_lat_c)
miroc_chem_lon_diff =  360./len(miroc_chem_lon_c)


min_lat = np.min([cesm_cam_lat_diff,cmam_lat_diff,geos_ccm_lat_diff,geos_chem_lat_diff,gfdl_lat_diff,giss_lat_diff,miroc_chem_lat_diff])
min_lon = np.min([cesm_cam_lon_diff,cmam_lon_diff,geos_ccm_lon_diff,geos_chem_lon_diff,gfdl_lon_diff,giss_lon_diff,miroc_chem_lon_diff])

min_lat = 1.
min_lon = 1.

new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_cesm_cam = modules.regrid_interp(cesm_cam_param,cesm_cam_lat_e,cesm_cam_lat_c,cesm_cam_lon_e,cesm_cam_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_cmam = modules.regrid_interp(cmam_param,cmam_lat_e,cmam_lat_c,cmam_lon_e,cmam_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_geos_ccm = modules.regrid_interp(geos_ccm_param,geos_ccm_lat_e,geos_ccm_lat_c,geos_ccm_lon_e,geos_ccm_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_geos_chem = modules.regrid_interp(geos_chem_param,geos_chem_lat_e,geos_chem_lat_c,geos_chem_lon_e,geos_chem_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_gfdl = modules.regrid_interp(gfdl_param,gfdl_lat_e,gfdl_lat_c,gfdl_lon_e,gfdl_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_giss = modules.regrid_interp(giss_param,giss_lat_e,giss_lat_c,giss_lon_e,giss_lon_c,min_lat,min_lon,data_type,max_ph)
new_lat_e,new_lat_c,new_lon_e,new_lon_c,z_miroc_chem = modules.regrid_interp(miroc_chem_param,miroc_chem_lat_e,miroc_chem_lat_c,miroc_chem_lon_e,miroc_chem_lon_c,min_lat,min_lon,data_type,max_ph)

if (data_type == 'amp') or (data_type == 'ave'):
    z_std = np.std((z_cesm_cam,z_cmam,z_geos_ccm,z_geos_chem,z_gfdl,z_giss,z_miroc_chem),axis=0)
    z_ave = np.average((z_cesm_cam,z_cmam,z_geos_ccm,z_geos_chem,z_gfdl,z_giss,z_miroc_chem),axis=0)
    z_pc = (z_std/z_ave)*100.
if data_type == 'ph':
    data = np.array([z_cesm_cam,z_cmam,z_geos_ccm,z_geos_chem,z_gfdl,z_giss,z_miroc_chem])
    z_ave,z_std = modules.yamartino_ave_std(data,max_ph)
    z_pc = (z_std/z_ave)*100.

#save out sig periods to netcdf
root_grp_period = Dataset('ACCMIP_hourly_%s%s_std.nc'%(period,data_type), 'w')
root_grp_period.createDimension('lat_centre', len(new_lat_c))
root_grp_period.createDimension('lon_centre', len(new_lon_c))
root_grp_period.createDimension('lat_edges', len(new_lat_e))
root_grp_period.createDimension('lon_edges', len(new_lon_e))
cesmcam = root_grp_period.createVariable('cesmcam', 'f8', ('lat_centre','lon_centre'))
cmam = root_grp_period.createVariable('cmam', 'f8', ('lat_centre','lon_centre'))
geosccm = root_grp_period.createVariable('geosccm', 'f8', ('lat_centre','lon_centre'))
geoschem = root_grp_period.createVariable('geoschem', 'f8', ('lat_centre','lon_centre'))
gfdl = root_grp_period.createVariable('gfdl', 'f8', ('lat_centre','lon_centre'))
giss = root_grp_period.createVariable('giss', 'f8', ('lat_centre','lon_centre'))
mirocchem = root_grp_period.createVariable('mirocchem', 'f8', ('lat_centre','lon_centre'))
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
geosccm[:]  = z_geos_ccm
geoschem[:]  = z_geos_chem
gfdl[:]  = z_gfdl
giss[:]  = z_giss
mirocchem[:]  = z_miroc_chem
