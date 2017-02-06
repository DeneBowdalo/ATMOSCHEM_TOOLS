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

#models = ['CESMCAM', 'CMAM', 'GFDLAM3','GISSE2R', 'MIROCCHEM', 'MOCAGE', 'UMCAM']
models = ['GFDLAM3']
all_spec = ['co','nox','ch4','isop','oh','pan','precip','temp','eminox','eminh3','emivoc']
seasons = ['','winter','spring','summer','autumn']

for model in models:
    if (model == 'CESMCAM') or (model == 'CMAM') or (model == 'MIROCCHEM') or (model == 'MOCAGE'):
        year2000 = '2000_2003'
        year2100 = '2100_2103'
        year2000_s = 0
        year2000_e = 36                                                                                                                                  
        year2100_s = 60
        year2100_e = 96

    if (model == 'GFDLAM3'):
        year2000 = '2001_2004'
        year2100 = '2101_2104'
        year2000_s = 12
        year2000_e = 48
        year2100_s = 72
        year2100_e = 108

    if (model == 'GISSE2R'):
        year2000 = '2000_2003'
        year2100 = '2102_2105'
        year2000_s = 0
        year2000_e = 36                                                                                                                                  
        year2100_s = 84
        year2100_e = 96

    if (model == 'UMCAM'):
        year2000 = '2000_2003'
        year2100 = '2097_2100'        
        year2000_s = 0
        year2000_e = 36                                                                                                                                  
        year2100_s = 24
        year2100_e = 60

    #read in 2000 period data
    f2000 = '/work/home/db876/grid/O3/2000_2012/%s/%s_SURFACE_*_*_*_H_*/LSP_stats.nc'%(year2000,model)

    #read in 2100 period data
    f2100 = '/work/home/db876/grid/O3/2095_2111/%s/%s_SURFACE_*_*_*_H_rcp85/LSP_stats.nc'%(year2100,model)

    #read in 2100 model fixed emissions period data
    f2100e = '/work/home/db876/grid/O3/2095_2111/%s/%s_SURFACE_*_*_*_H_rcp85em2000/LSP_stats.nc'%(year2100,model)

    root_2000 = Dataset(f2000)
    lat_e = root_2000.variables['lat_edges'][:]
    lon_e = root_2000.variables['lon_edges'][:]
    lat_c = root_2000.variables['lat_centre'][:]
    lon_c = root_2000.variables['lon_centre'][:]
    size = len(lat_c)*len(lon_c)

    #read in monthly emissions and species for year 2000
    monthly_2000 = Dataset('/work/home/db876/modelling/ACCMIP/%s/%s_SURFACE_2000_2012_*_*_*_M_*.nc'%(model,model))
    m_2000_ch4 = monthly_2000.variables['ch4'][:,:,:]
    m_2000_co = monthly_2000.variables['co'][:,:,:]
    m_2000_no = monthly_2000.variables['no'][:,:,:]
    m_2000_no2 = monthly_2000.variables['no2'][:,:,:]   
    m_2000_isop = monthly_2000.variables['isop'][:,:,:]
    m_2000_oh = monthly_2000.variables['oh'][:,:,:]
    m_2000_pan = monthly_2000.variables['pan'][:,:,:]
    m_2000_temp = monthly_2000.variables['temp'][:,:,:]
    m_2000_precip = monthly_2000.variables['precip'][:,:,:]
    m_2000_eminox = monthly_2000.variables['eminox'][:,:,:]
    m_2000_emivoc = monthly_2000.variables['emivoc'][:,:,:]
    m_2000_eminh3 = monthly_2000.variables['eminh3'][:,:,:]
    if model != 'GISSE2R':
        monthly_2100 = Dataset('/work/home/db876/modelling/ACCMIP/%s/%s_SURFACE_2095_2111_*_*_*_M_rcp85.nc'%(model,model))
        m_2100_ch4 = monthly_2100.variables['ch4'][:,:,:]
        m_2100_co = monthly_2100.variables['co'][:,:,:]
        m_2100_no = monthly_2100.variables['no'][:,:,:]
        m_2100_no2 = monthly_2100.variables['no2'][:,:,:]    
        m_2100_isop = monthly_2100.variables['isop'][:,:,:]
        m_2100_oh = monthly_2100.variables['oh'][:,:,:]
        m_2100_pan = monthly_2100.variables['pan'][:,:,:]
        m_2100_temp = monthly_2100.variables['temp'][:,:,:]
        m_2100_precip = monthly_2100.variables['precip'][:,:,:]
        m_2100_eminox = monthly_2100.variables['eminox'][:,:,:]
        m_2100_emivoc = monthly_2100.variables['emivoc'][:,:,:]
        m_2100_eminh3 = monthly_2100.variables['eminh3'][:,:,:]
    if model != 'CMAM':
        monthly_2100e = Dataset('/work/home/db876/modelling/ACCMIP/%s/%s_SURFACE_2095_2111_*_*_*_M_rcp85em2000.nc'%(model,model))
        m_2100e_ch4 = monthly_2100e.variables['ch4'][:,:,:]
        m_2100e_co = monthly_2100e.variables['co'][:,:,:]
        m_2100e_no = monthly_2100e.variables['no'][:,:,:]
        m_2100e_no2 = monthly_2100e.variables['no2'][:,:,:]                                                 
        m_2100e_isop = monthly_2100e.variables['isop'][:,:,:]
        m_2100e_oh = monthly_2100e.variables['oh'][:,:,:]
        m_2100e_pan = monthly_2100e.variables['pan'][:,:,:]
        m_2100e_temp = monthly_2100e.variables['temp'][:,:,:]
        m_2100e_precip = monthly_2100e.variables['precip'][:,:,:]
        m_2100e_eminox = monthly_2100e.variables['eminox'][:,:,:]
        m_2100e_emivoc = monthly_2100e.variables['emivoc'][:,:,:]
        m_2100e_eminh3 = monthly_2100e.variables['eminh3'][:,:,:]


    for species in all_spec:
        for s in seasons:
            if s == 'winter':
                inds_2000 = np.sort(range(year2000_s+0,year2000_e,12) + range(year2000_s+1,year2000_e,12) + range(year2000_s+11,year2000_e,12))
                inds_2100 = np.sort(range(year2100_s+0,year2100_e,12) + range(year2100_s+1,year2100_e,12) + range(year2100_s+11,year2100_e,12))
            elif s == 'spring':
                inds_2000 = np.sort(range(year2000_s+2,year2000_e,12) + range(year2000_s+3,year2000_e,12) + range(year2000_s+4,year2000_e,12))                                                                                       
                inds_2100 = np.sort(range(year2100_s+2,year2100_e,12) + range(year2100_s+3,year2100_e,12) + range(year2100_s+4,year2100_e,12))
            elif s == 'summer':
                inds_2000 = np.sort(range(year2000_s+5,year2000_e,12) + range(year2000_s+6,year2000_e,12) + range(year2000_s+7,year2000_e,12))
                inds_2100 = np.sort(range(year2100_s+5,year2100_e,12) + range(year2100_s+6,year2100_e,12) + range(year2100_s+7,year2100_e,12))
            elif s == 'autumn':
                inds_2000 = np.sort(range(year2000_s+8,year2000_e,12) + range(year2000_s+9,year2000_e,12) + range(year2000_s+10,year2000_e,12))
                inds_2100 = np.sort(range(year2100_s+8,year2100_e,12) + range(year2100_s+9,year2100_e,12) + range(year2100_s+10,year2100_e,12))
            else:
                inds_2000 = np.arange(year2000_s,year2000_e)
                inds_2100 = np.arange(year2100_s,year2100_e)

            #read in monthly emissions and species for year 2000
            monthly_2000_ch4 = np.average(m_2000_ch4[inds_2000,:,:],axis=0)
            monthly_2000_co = np.average(m_2000_co[inds_2000,:,:],axis=0)
            monthly_2000_no = np.average(m_2000_no[inds_2000,:,:],axis=0)
            monthly_2000_no2 = np.average(m_2000_no2[inds_2000,:,:],axis=0)
            monthly_2000_nox = monthly_2000_no + monthly_2000_no2
            monthly_2000_isop = np.average(m_2000_isop[inds_2000,:,:],axis=0)
            monthly_2000_oh = np.average(m_2000_oh[inds_2000,:,:],axis=0)
            monthly_2000_pan = np.average(m_2000_pan[inds_2000,:,:],axis=0)
            monthly_2000_temp = np.average(m_2000_temp[inds_2000,:,:],axis=0)
            monthly_2000_precip = np.average(m_2000_precip[inds_2000,:,:],axis=0)
            monthly_2000_eminox = np.sum(m_2000_eminox[inds_2000,:,:],axis=0)
            monthly_2000_emivoc = np.sum(m_2000_emivoc[inds_2000,:,:],axis=0)
            monthly_2000_eminh3 = np.sum(m_2000_eminh3[inds_2000,:,:],axis=0)
            if model != 'GISSE2R':
                monthly_2100_ch4 = np.average(m_2100_ch4[inds_2100,:,:],axis=0)
                monthly_2100_co = np.average(m_2100_co[inds_2100,:,:],axis=0)
                monthly_2100_no = np.average(m_2100_no[inds_2100,:,:],axis=0)
                monthly_2100_no2 = np.average(m_2100_no2[inds_2100,:,:],axis=0)
                monthly_2100_nox = monthly_2100_no + monthly_2100_no2
                monthly_2100_isop = np.average(m_2100_isop[inds_2100,:,:],axis=0)
                monthly_2100_oh = np.average(m_2100_oh[inds_2100,:,:],axis=0)
                monthly_2100_pan = np.average(m_2100_pan[inds_2100,:,:],axis=0)
                monthly_2100_temp = np.average(m_2100_temp[inds_2100,:,:],axis=0)
                monthly_2100_precip = np.average(m_2100_precip[inds_2100,:,:],axis=0)
                monthly_2100_eminox = np.sum(m_2100_eminox[inds_2100,:,:],axis=0)
                monthly_2100_emivoc = np.sum(m_2100_emivoc[inds_2100,:,:],axis=0)
                monthly_2100_eminh3 = np.sum(m_2100_eminh3[inds_2100,:,:],axis=0)
            if model != 'CMAM':
                monthly_2100e_ch4 = np.average(m_2100e_ch4[inds_2100,:,:],axis=0)
                monthly_2100e_co = np.average(m_2100e_co[inds_2100,:,:],axis=0)
                monthly_2100e_no = np.average(m_2100e_no[inds_2100,:,:],axis=0)
                monthly_2100e_no2 = np.average(m_2100e_no2[inds_2100,:,:],axis=0)
                monthly_2100e_nox = monthly_2100e_no + monthly_2100e_no2
                monthly_2100e_isop = np.average(m_2100e_isop[inds_2100,:,:],axis=0)
                monthly_2100e_oh = np.average(m_2100e_oh[inds_2100,:,:],axis=0)
                monthly_2100e_pan = np.average(m_2100e_pan[inds_2100,:,:],axis=0)
                monthly_2100e_temp = np.average(m_2100e_temp[inds_2100,:,:],axis=0)
                monthly_2100e_precip = np.average(m_2100e_precip[inds_2100,:,:],axis=0)
                monthly_2100e_eminox = np.sum(m_2100e_eminox[inds_2100,:,:],axis=0)
                monthly_2100e_emivoc = np.sum(m_2100e_emivoc[inds_2100,:,:],axis=0)
                monthly_2100e_eminh3 = np.sum(m_2100e_eminh3[inds_2100,:,:],axis=0)

            if species == 'nox':
                cb_l = 'NOx (ppb)'
                if model == 'GISSE2R':
                    params = [monthly_2000_nox,monthly_2000_nox,monthly_2100e_nox]
                elif model == 'CMAM':
                    params = [monthly_2000_nox,monthly_2100_nox,monthly_2100_nox]
                else:
                    params = [monthly_2000_nox,monthly_2100_nox,monthly_2100e_nox]
                min = np.min(params)
                max = np.max(params)
            elif species == 'isop':
                cb_l = 'ISOP (ppb)'
                if model == 'GISSE2R':
                    params = [monthly_2000_isop,monthly_2000_isop,monthly_2100e_isop]
                elif model == 'CMAM':
                    params = [monthly_2000_isop,monthly_2100_isop,monthly_2100_isop]
                else:                                                                                                                                        
                    params = [monthly_2000_isop,monthly_2100_isop,monthly_2100e_isop]
                min = np.min(params)
                max = np.max(params)
            elif species == 'co':
                cb_l = 'CO (ppb)'        
                if model == 'GISSE2R':
                    params = [monthly_2000_co,monthly_2000_co,monthly_2100e_co]
                elif model == 'CMAM':
                    params = [monthly_2000_co,monthly_2100_co,monthly_2100_co]
                else:                                                                                                                                        
                    params = [monthly_2000_co,monthly_2100_co,monthly_2100e_co]
                min = np.min(params)
                max = np.max(params)
            elif species == 'oh':
                cb_l = 'OH (ppb)'        
                if model == 'GISSE2R':
                    params = [monthly_2000_oh,monthly_2000_oh,monthly_2100e_oh]
                elif model == 'CMAM':
                    params = [monthly_2000_oh,monthly_2100_oh,monthly_2100_oh]
                else:                                                                                                                                    
                    params = [monthly_2000_oh,monthly_2100_oh,monthly_2100e_oh]
                min = np.min(params)
                max = np.max(params)
            elif species == 'pan':
                cb_l = 'PAN (ppb)'        
                if model == 'GISSE2R':
                    params = [monthly_2000_pan,monthly_2000_pan,monthly_2100e_pan]
                elif model == 'CMAM':
                    params = [monthly_2000_pan,monthly_2100_pan,monthly_2100_pan]
                else:                                                                                                                                    
                    params = [monthly_2000_pan,monthly_2100_pan,monthly_2100e_pan]
                min = np.min(params)                                                                                                                     
                max = np.max(params)
            elif species == 'temp':
                cb_l = 'TEMP (K)'        
                if model == 'GISSE2R':
                    params = [monthly_2000_temp,monthly_2000_temp,monthly_2100e_temp]                                                                       
                elif model == 'CMAM':
                    params = [monthly_2000_temp,monthly_2100_temp,monthly_2100_temp]
                else:                                                                                                                                    
                    params = [monthly_2000_temp,monthly_2100_temp,monthly_2100e_temp]
                min = np.min(params)                                                                                                                     
                max = np.max(params)
            elif species == 'precip':
                cb_l = 'PRECIP'        
                if model == 'GISSE2R':
                    params = [monthly_2000_precip,monthly_2000_precip,monthly_2100e_precip]                                                                    
                elif model == 'CMAM':
                    params = [monthly_2000_precip,monthly_2100_precip,monthly_2100_precip]
                else:                                                                                                                                    
                    params = [monthly_2000_precip,monthly_2100_precip,monthly_2100e_precip]
                min = np.min(params)                                                                                                                     
                max = np.max(params)
            elif species == 'ch4':
                cb_l = 'CH4 (ppb)'
                if model == 'GISSE2R':
                    params = [monthly_2000_ch4,monthly_2000_ch4,monthly_2100e_ch4]
                elif model == 'CMAM':
                    params = [monthly_2000_ch4,monthly_2100_ch4,monthly_2100_ch4]
                else:                                                                                                                                        
                    params = [monthly_2000_ch4,monthly_2100_ch4,monthly_2100e_ch4]
                min = 0
                max= 5000
            elif species == 'eminox':
                cb_l = 'NOx Emissions'
                if model == 'GISSE2R':
                    params = [monthly_2000_eminox,monthly_2000_eminox,monthly_2100e_eminox]
                elif model == 'CMAM':
                    params = [monthly_2000_eminox,monthly_2100_eminox,monthly_2100_eminox]
                else:                                                                                                                                        
                    params = [monthly_2000_eminox,monthly_2100_eminox,monthly_2100e_eminox]
                min = np.min(params)
                max = np.max(params)
            elif species == 'eminh3':
                cb_l = 'NH3 Emissions'                                                                                                                            
                if model == 'GISSE2R':
                    params = [monthly_2000_eminh3,monthly_2000_eminh3,monthly_2100e_eminh3]
                elif model == 'CMAM':
                    params = [monthly_2000_eminh3,monthly_2100_eminh3,monthly_2100_eminh3]
                else:                                                                                                                                        
                    params = [monthly_2000_eminh3,monthly_2100_eminh3,monthly_2100e_eminh3]
                min = np.min(params)
                max = np.max(params)
            elif species == 'emivoc':
                cb_l = 'VOC Emissions'                                                                                                                            
                if model == 'GISSE2R':
                    params = [monthly_2000_emivoc,monthly_2000_emivoc,monthly_2100e_emivoc]
                elif model == 'CMAM':
                    params = [monthly_2000_emivoc,monthly_2100_emivoc,monthly_2100_emivoc]
                else:                                                                                                                                        
                    params = [monthly_2000_emivoc,monthly_2100_emivoc,monthly_2100e_emivoc]
                min = np.min(params)
                max = np.max(params)

            #set up plot
            fig =plt.figure(figsize=(17,9.5))
            fig.patch.set_facecolor('white')
        
            ax1 = plt.subplot2grid((2,2), (0,0))
            ax2 = plt.subplot2grid((2,2), (1,0))
            ax3 = plt.subplot2grid((2,2), (0,1))

            #setup basemap projection
            m = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c',ax=ax1)
            m2 = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c',ax=ax2)
            m3 = Basemap(projection='cyl',llcrnrlat=lat_e[0],urcrnrlat=lat_e[-1],llcrnrlon=lon_e[0],urcrnrlon=lon_e[-1],resolution='c',ax=ax3)

            m.drawcoastlines()
            m.drawmapboundary()
            m2.drawcoastlines()
            m2.drawmapboundary()
            m3.drawcoastlines()
            m3.drawmapboundary()

            pl = m.pcolor(lon_e,lat_e,params[0], vmin=min, vmax=max,linewidth=0.5,cmap=plt.cm.gist_earth)
            if model != 'GISSE2R': 
                pl = m2.pcolor(lon_e,lat_e,params[1], vmin=min, vmax=max,linewidth=0.5,cmap=plt.cm.gist_earth)
            if model != 'CMAM':
                pl = m3.pcolor(lon_e,lat_e,params[2], vmin=min, vmax=max,linewidth=0.5,cmap=plt.cm.gist_earth)
            plt.tight_layout()

            cbar_ax = fig.add_axes([0.53, 0.25, 0.44, 0.07])
            cb = fig.colorbar(pl,orientation='horizontal',format='%.1f',cax=cbar_ax)
            cb.ax.tick_params(labelsize=18)
            cb.set_label('%s'%(cb_l), fontsize = 23)

            ax1.set_title('%s 2000'%(species.upper()),fontsize=20)
            if model != 'GISSE2R':                                                                                                                       
                ax2.set_title('%s 2100'%(species.upper()),fontsize=20)
            if (model != 'CMAM'):
                ax3.set_title('%s 2100e'%(species.upper()),fontsize=20)

            plt.savefig('plots/%s_%s%s.png'%(model,species,s))
            #plt.show()
