import numpy as np
import csv
import glob
import datetime
from datetime import date,time
import matplotlib.pyplot as plt
import logging as log
import scipy.stats as stats
import lomb_phase
import lomb_phase_spec
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft
import collections
from collections import Counter
from collections import defaultdict
from scipy import signal
import numpy.fft as FFT
import scipy.signal
from math import radians, acos, cos, sin, asin, tan,  sqrt, degrees
import math
from mpl_toolkits.basemap import Basemap, shiftgrid, interp, maskoceans
import matplotlib as mpl
from netCDF4 import Dataset
import urllib2
import json
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from incf.countryutils import transformations
from geopy.geocoders import Nominatim
import pandas as pd
from netCDF4 import num2date, date2num
import ephem as ep
#import seaborn as sns

#DENE BOWDALO

#A collection of python routines used during PhD 

#SECTION 1 - OBSERVATIONAL PROCESSING
#SECTION 2 - MODEL PROCESSING
#SECTION 3 - SPECTRAL ANALYSIS
#SECTION 4 - STATISTICS
#SECTION 5 - INTERACTIVE PLOTTING

#-------------------------------------------------------------
#-------------------------------------------------------------
#SECTION 1 - OBSERVATIONAL PROCESSING

#1.01 - get observational filename (from directory structure)
#1.02 - read observational data from netcdf for 1 site
#1.03 - read observational data from netcdf for all sites
#1.04 - quality checks for _ALL tagged observational data
#1.05 - quality checks for _NR tagged observational data
#1.06 - quality checks for _PERIODIC tagged observational data
#1.07 - average observational data into daily/monthly resolution 
#1.08 - write out _ALL/_NR tagged observational data to netcdf
#1.09 - write out _PERIODIC tagged observational data to netcdf
#1.10 - set defined areas for site averaging
#1.11 - set dicts for site area classification
#1.12 - get area cut for sites (using lat,lon and tag)
#1.13 - get area for sites (using lat,lon and tag)
#1.14 - get country (from lat,lon) using google api 
#1.15 - get anthrome classification
#1.16 - get observational continental tags (manually put together)
#1.17 - read obs/model periodic metrics from netcdf
#1.18 - cut time series into 4 seasonal time series


#1.01
#----------

#GET OBS FILE INFO, NEED DIRECTORIES STRUCTURED IN FOLLOWING WAY:
#..../SPECIES/STARTYEAROFFILE_ENDYEAROFFILE/STARTYEARWANTED_ENDYEARWANTED/MODELNAME_VERTRES_VERSION_HORIZRES_MET_TIMERES_ADDITIONAL
#USE * IF NO DETAILS FOR CERTAIN ASPECT, i.e ACCMIP ONLY HAVE 1 VERSION,HORIZRES,MET
def get_obs_info(present_dir):

    paths = present_dir.split("/")
    species = paths[-4]
    years = paths[-2]
    start_year = years[:4]
    end_year = years[5:9]
    model_path = paths[-1] 
    
    data_split = model_path.split('_')
    
    #if reading obs. info from a model directory
    if len(data_split) == 7:
        model = data_split[0]
        vres = data_split[1]
        version = data_split[2]
        hres = data_split[3]
        met = data_split[4]
        timeres = data_split[5]
        additional = data_split[6]
        fname = '/work/home/db876/observations/surface/%s/process/GLOBAL_%s_%s_%s_%s_%s_%sP.nc'%(species,vres,species,start_year,end_year,timeres,timeres)
        dir = 'model'
    
    #if reading obs. info from a obs. directory
    if len(data_split) == 3:
        vres = data_split[1]
        timeres = data_split[2]
        fname = '/work/home/db876/observations/surface/%s/process/GLOBAL_%s_%s_%s_%s_%s_%sP.nc'%(species,vres,species,start_year,end_year,timeres,timeres)
        dir = 'obs'
        
    return fname,species,start_year,end_year,vres,timeres,dir
    
#1.02
#----------

#READ ONE OBS SITE DATA FILE FROM NETCDF
def read_obs_one(obs_file,species,ref,start_year,end_year):
    root_grp = Dataset(obs_file)   
    raw_time = root_grp.variables['time'][:]
    
    ref_group = root_grp.groups[ref]
    std_var = ref_group.variables[species.lower()][:]
    obs_lat = ref_group.latitude
    obs_lon = ref_group.longitude
    obs_alt = ref_group.altitude
    obs_process_group = ref_group.process_group
    obs_raw_class = ref_group.raw_site_class
    obs_anthrome_class = ref_group.raw_site_class
    
    #cut time
    start_time = date2num(datetime.datetime(int(start_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    end_time = date2num(datetime.datetime(int(end_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    start_ind = np.searchsorted(raw_time,start_time)
    end_ind = np.searchsorted(raw_time,end_time)
    if (start_time < raw_time[0]) & (end_time > raw_time[-1]):
        raw_time = raw_time[:]
        std_var = std_var[:]
    elif start_time < raw_time[0]:
        raw_time = raw_time[:end_ind]
        std_var = std_var[:end_ind]
    elif end_time > raw_time[-1]:
        raw_time = raw_time[start_ind:]
        std_var = std_var[start_ind:]
    else:
        raw_time = raw_time[start_ind:end_ind]
        std_var = std_var[start_ind:end_ind]
    
    ref_time = raw_time - raw_time[0]
    
    datetime_time = num2date(raw_time,units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    
    #get raw and ref time in days
    raw_time = raw_time / 24.
    ref_time = ref_time / 24.
    
    std_var_mask = np.ma.masked_where(std_var<=0,std_var)
    pd_var = pd.Series(std_var_mask, index=datetime_time)
    
    return ref,raw_time,ref_time,datetime_time,std_var,pd_var,obs_lat,obs_lon,obs_alt,obs_process_group,obs_raw_class,obs_anthrome_class 
 
#1.03
#----------

#READ ALL OBS SITE DATA FILE FROM NETCDF
def read_obs_all(obs_file,species,start_year,end_year):
    
    root_grp = Dataset(obs_file)
    raw_time = root_grp.variables['time'][:]
    refs_dict = root_grp.groups
    
    refs = []
    obs_lats = []
    obs_lons = []
    obs_alt = []
    std_var = []
    gap_inds = []
    obs_groups = []
    obs_raw_class = []
    obs_anthrome_class = []
    
    for i in refs_dict.keys():
        i = i.encode('ascii')
        refs.append(i)
    
    for ref in refs:
        ref_group = root_grp.groups[ref]
        std_var.append(ref_group.variables[species.lower()][:])
        obs_lats = np.append(obs_lats,ref_group.latitude)
        obs_lons = np.append(obs_lons,ref_group.longitude)
        obs_alt = np.append(obs_alt,ref_group.altitude)
        obs_groups =  np.append(obs_groups,ref_group.process_group)
        obs_raw_class = np.append(obs_raw_class,ref_group.raw_site_class)
        obs_anthrome_class = np.append(obs_anthrome_class,ref_group.anthrome_site_class)
        gap_inds.append(std_var < 0)
    
    start_time = date2num(datetime.datetime(int(start_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    end_time = date2num(datetime.datetime(int(end_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    start_ind = np.searchsorted(raw_time,start_time)
    end_ind = np.searchsorted(raw_time,end_time)
    
    std_var = np.array(std_var)
    
    if (start_time < raw_time[0]) & (end_time > raw_time[-1]):
        raw_time = raw_time[:]
        std_var = std_var[:,:]
    elif start_time < raw_time[0]:
        raw_time = raw_time[:end_ind]
        std_var = std_var[:,:end_ind]
    elif end_time > raw_time[-1]:
        raw_time = raw_time[start_ind:]
        std_var = std_var[:,start_ind:]
    else:
        raw_time = raw_time[start_ind:end_ind]
        std_var = std_var[:,start_ind:end_ind]
    
    ref_time = raw_time - raw_time[0]
    datetime_time = num2date(raw_time,units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    
    #get raw and ref time in days
    raw_time = raw_time / 24.
    ref_time = ref_time / 24.
    
    for i in range(len(refs)):
        refs[i] = refs[i].lower()

    return refs,raw_time,ref_time,datetime_time,np.array(std_var),obs_lats,obs_lons,obs_alt,obs_groups,obs_raw_class,obs_anthrome_class,gap_inds

def sampling_and_instruments_by_processgroup(ref,process_group,species,raw_st,raw_mm,full_data_after_flagsandlod,full_data,raw_indices,unknown_mm,unknown_mm_refs,no2_type):

    n_points = len(full_data)

#CREATE RAW AND PROCESSED SAMPLING/MEASURMENT GRIDS
#RAW IS BEFORE FLAGSANDLOD (NO POINTS REMOVED IF GOT SO FAR)
#FLAGSANDLOD IS AFTER THE POINTS REMOVAL BY THAT TEST
#P IS AFTER DUPLICATE POINT REMOVAL
    raw_mm_grid = np.chararray(n_points,itemsize=200)
    raw_st_grid = np.chararray(n_points,itemsize=200)
    #p_mm_grid = np.chararray(n_points,itemsize=200)
    #p_st_grid = np.chararray(n_points,itemsize=200)
    p_mm_grid = np.empty(n_points)
    p_st_grid = np.empty(n_points)
    
    raw_mm_grid[:] = 'na'
    raw_st_grid[:] = 'na'
    p_mm_grid[:] = -1
    p_st_grid[:] = -1
    
    raw_mm_grid = np.array(raw_mm_grid,dtype='object')
    raw_st_grid = np.array(raw_st_grid,dtype='object')
    #p_mm_grid = np.array(p_mm_grid,dtype='object')
    #p_st_grid = np.array(p_st_grid,dtype='object')
    
#FILL RAW GRIDS WITH RAW MM AND ST (NONE REMOVED BY FLAGSANDLOD OR DUPLICATES)
    raw_st_grid[raw_indices] = raw_st
    raw_mm_grid[raw_indices] = raw_mm

#--------------------------------------------------------------------------------------------------- 
    if process_group == 'AirBase':
        continuous_st_list = ['continuous']
        remote_st_list = []
        flask_st_list = []
        filter_st_list = ['filter']
        
        uv_mm_list = ['Ultraviolet (UV) photometryEnvironnement S.A. Model O331M UV Ozone Analyzer','Ultraviolet (UV) photometryMonitor Labs model 9800','Ultraviolet (UV) photometryThermo model 42 NO/Nox analyser','Ultraviolet (UV) photometryUNKNOWN','Ultraviolet (UV) photometryMCV 48-AUV','Ultraviolet (UV) photometryTeledyne API 400A UV photometric O3 analyser','Ultraviolet (UV) photometryThermo model 48 CO analyser','Ultraviolet (UV) photometryTeledyne API 400E UV photometric O3 analyser','Ultraviolet (UV) photometryHoriba model APOA 300 O3 analyser', \
                      'Ultraviolet (UV) photometry342 M','Ultraviolet (UV) photometryMonitor Labs model 9812 O3 analyser','Ultraviolet (UV) photometryHoriba model APOA 350E O3 analyser','Ultraviolet (UV) photometryENVIRONMENT 1003 AH','Ultraviolet (UV) photometryC.S.I. 3.100','Ultraviolet (UV) photometryDASIBI 1003 O3 analyser','Ultraviolet (UV) photometryMonitor Labs undetermined','Ultraviolet (UV) photometryMonitor Labs model 9810B O3 analyser','Ultraviolet (UV) photometrytoo generic','Ultraviolet (UV) photometryThermo 49 CPS Ozone Primary Standard', \
                      'Ultraviolet (UV) photometryDASIBI','UV fluorescencetoo generic','Ultraviolet (UV) photometryDASIBI 1003-PC O3 analyser','Ultraviolet (UV) photometryThermo model 43 SO2 analyser','Ultraviolet (UV) photometryThermo model 49i O3 analyser','Ultraviolet (UV) photometryDASIBI 1008-PC O3 analyser','Ultraviolet (UV) photometryDASIBI 1008-RS O3 analyser','Ultraviolet (UV) photometryEnvironnement S.A. Model O341M UV Ozone Analyzer','Ultraviolet (UV) photometryISEO Argopol-SAM-XAIR',  \
                      'Ultraviolet (UV) photometryEnvironnement S.A. Model O342M UV Ozone Analyze','Ultraviolet (UV) photometryHoriba model APOA 370 O3 analyser','Ultraviolet (UV) photometryDASIBI 1008-AH O3 analyser','UV fluorescenceThermo 49c','Ultraviolet (UV) photometryPHILIPS K50110/00 UV Photometric O3 analyser','Ultraviolet (UV) photometryMonitor Labs model 8810 O3 analyser','Ultraviolet (UV) photometryPHILIPS K50094 API 400','Ultraviolet (UV) photometryORION','Ultraviolet (UV) photometryThermo model 49w O3 analyser', \
                      'Ultraviolet (UV) photometryMonitor Labs model 9810 O3 analyser','Ultraviolet (UV) photometryCOLUMBIA SCIENTIFIC IC 3100','Ultraviolet (UV) photometry2008A','Ultraviolet (UV) photometryThermo model 43s SO2 analyser','Ultraviolet (UV) photometryMLU','Ultraviolet (UV) photometryThermo model 49 O3 analyser','Ultraviolet (UV) photometryDASIBI 1108 O3 analyser','Ultraviolet (UV) photometryAMIBRACK','Ultraviolet (UV) photometryThermo model 49c O3 analyser','UV fluorescenceUNKNOWN','Ultraviolet (UV) photometryTeledyne API 400 UV photometric O3 analyser', \
                      'UV fluorescenceTeledyne API 400 UV photometric O3 analyser','Ultraviolet (UV) photometryMonitor Labs model 9830 CO analyser','Ultraviolet (UV) photometryDASIBI 5014','Ultraviolet (UV) photometryEnvironics 300/ Environics','Ultraviolet (UV) photometryANALYSIS AUTOMATION Mod. 427','Ultraviolet (UV) photometryANALYSIS AUTOMATION','Ultraviolet (UV) photometryDASIBI 1008 O3 analyser','ultraviolet absorptionORION','Ultraviolet (UV) photometryMonitor Labs model 9811 O3 analyser','Ultraviolet (UV) photometryENVIRONMENT 1003RS', \
                      'UV absorption (ref)UNKNOWN','Ultraviolet (UV) photometryDASIBI 1003-RS O3 analyser','Ultraviolet (UV) photometryHoriba model APOA 350 O3 analyser','Ultraviolet (UV) photometrySFI O342M','UV fluorescenceMonitor Labs undetermined','Ultraviolet (UV) photometryDANI ENVIRONMENT 1003 AH','Ultraviolet (UV) photometryS-5014','Ultraviolet (UV) photometryThermo model 42 NO/Nox analyser','Ultraviolet (UV) photometryUNKNOWN', \
                      'Ultraviolet (UV) photometryHoriba model APNA 360 NOx analyser','Ultraviolet (UV) photometryMonitor Labs undetermined','Ultraviolet (UV) photometryTeledyne API 200A chemiluminescent NOx analyser','UV fluorescenceThermo model 42 NO/Nox analyser','Ultraviolet (UV) photometryContiflo','Ultraviolet (UV) photometryTeledyne API undertermined','UV fluorescenceThermo model 43a SO2 analyser','UV fluorescenceEnvironnement S.A. Model AF21M SO2 Analyzer','UV fluorescenceThermo model 43c SO2 analyser', \
                      'Ultraviolet (UV) photometryTeledyne API undertermined','UV fluorescenceUNKNOWN','UV fluorescenceEnvironnement S.A. Model AF21M SO2 Analyzer','Ultraviolet (UV) photometryUNKNOWN','Ultraviolet (UV) photometryThermo model 43 SO2 analyser','Ultraviolet (UV) photometryMonitor Labs model 9810 O3 analyser','Ultraviolet (UV) photometryHoriba model APNA 360 NOx analyser','UV fluorescenceUNKNOWN','Ultraviolet (UV) photometryTeledyne API 200A chemiluminescent NOx analyser','Ultraviolet (UV) photometryTeledyne API 200A chemiluminescent NOx analyser', \
                      'UV fluorescenceThermo model 43 SO2 analyser','PhotometryFILTER SAMPLING']
    
        vuf_mm_list = []
    
        cld_mol_mm_list = ['ChemiluminescenceHoriba model APNA 300 NOx analyser','ChemiluminescenceHoriba model APNA 300E NOx analyser','ChemiluminescenceHoriba model APNA 350 NOx analyser','ChemiluminescenceHoriba model APNA 350 NOx analyser','ChemiluminescenceHoriba model APNA 350E NOx analyser','ChemiluminescenceHoriba model APNA 350E NOx analyser','ChemiluminescenceHoriba model APNA 350E NOx analyser','ChemiluminescenceHoriba model APHA 360E hydrocarbons analyser','ChemiluminescenceHoriba model APNA 360E NOx analyser','ChemiluminescenceHoriba model APNA 360E NOx analyser', \
                           'ChemiluminescenceHoriba model APHA 360E hydrocarbons analyser','ChemiluminescenceHoriba model APNA 360 NOx analyser','UV fluorescenceHoriba model APNA 360 NOx analyser','ChemiluminescenceHoriba model APNA 360 NOx analyser','ChemiluminescenceHoriba model APNA 360 NOx analyser','ChemiluminescenceHoriba model APNA 370 NOx analyser','chemiluminescenceHORIBA APNA 370','ChemiluminescenceHoriba model APNA 370 NOx analyser', \
                           'ChemiluminescenceEnvironnement S.A. Model AC31M NO2 Analyzer','ChemiluminescenceEnvironnement S.A. Model AC31M NO2 Analyzer','ChemiluminescenceENVIRONMENT','ChemiluminescenceENVIRONMENT','chemiluminescenceENVIRONNEMENT AC 30M','ChemiluminescenceEnvironnement S.A. Model AC30M NO2 Analyzer','ChemiluminescenceEnvironnement S.A. Model AC30M NO2 Analyzer','ChemiluminescenceEnvironnement S.A. Model AC32M NO2 Analyzer','ChemiluminescenceEnvironnement S.A. Model AC32M NO2 Analyzer','ChemiluminescenceSFI AC32M','ChemiluminescenceENVIRONMENT ZC 32M','ChemiluminescenceENVIRONMENT ZC 32M', \
                           'ChemiluminescenceThermo model 42c NO/Nox analyser','ChemiluminescenceThermo model 42C-TL (Trace level Nox)','ChemiluminescenceThermo model 42c NO/Nox analyser','ChemiluminescenceThermo model 42C-TL (Trace level Nox)','ChemiluminescenceThermo model 42i NO/Nox analyser','ChemiluminescenceThermo model 42i-TL (Trace level Nox)','ChemiluminescenceThermo model 42i NO/Nox analyser','ChemiluminescenceThermo model 42i-TL (Trace level Nox)', \
                           'ChemiluminescenceTHERMO ENVIRONMENTAL INSTRUMENTS','ChemiluminescenceTHERMO ELECTRON INSTRUMENTS','ChemiluminescenceTHERMO ELECTRON INSTRUMENTS','ChemiluminescenceThermo model 14B chemiluminescence NO-NO2-Nox','ChemiluminescenceThermo model 14B/E chemiluminescence NO-NO2-Nox','ChemiluminescenceThermo model 14B chemiluminescence NO-NO2-Nox','ChemiluminescenceThermo model 14B/E chemiluminescence NO-NO2-Nox', \
                           'ChemiluminescenceThermo model 42 NO/Nox analyser','ChemiluminescenceThermo model 42 NO/Nox analyser','ChemiluminescenceThermo model 42s NO/Nox analyser','ChemiluminescenceThermo model 42s NO/Nox analyser','ChemiluminescenceThermo model 42w NO/Nox analyser','ChemiluminescenceThermo model 42w NO/Nox analyser', \
                           'ChemiluminescenceMonitor Labs undetermined','ChemiluminescenceMonitor Labs undetermined','ChemiluminescenceMonitor Labs model 8440 NOx analyser','ChemiluminescenceMonitor Labs model 8440 NOx analyser','ChemiluminescenceMonitor Labs model 8840 NOx analyser','ChemiluminescenceMonitor Labs model 8840 NOx analyser','ChemiluminescenceMonitor Labs model 8841 NOx analyser','ChemiluminescenceMonitor Labs model 8841 NOx analyser','ChemiluminescenceMonitor Labs model 8941A NOx analyser','ChemiluminescenceMonitor Labs model 8941A NOx analyser', \
                           'ChemiluminescenceMonitor Labs model 9841A NOx analyser','ChemiluminescenceMonitor Labs model 9841A NOx analyser','ChemiluminescenceMonitor Labs model 9841 NOx analyser','ChemiluminescenceMonitor Labs model 9841 NOx analyser','ChemiluminescenceMonitor Labs model 9841B NOx analyser','ChemiluminescenceMonitor Labs model 9841B NOx analyser','ChemiluminescenceMonitor Labs model 9841B NOx analyser','ChemiluminescenceMonitor Labs model 9841T NOx analyser', \
                           'ChemiluminescenceTeledyne API undertermined','ChemiluminescenceTeledyne API undertermined','ChemiluminescenceTeledyne API undertermined','ChemiluminescenceTeledyne API 200 chemiluminescent NOx analyser','ChemiluminescenceTeledyne API 200 chemiluminescent NOx analyser','ChemiluminescenceTeledyne API 200A chemiluminescent NOx analyser','ChemiluminescenceTeledyne API 200A chemiluminescent NOx analyser','ChemiluminescenceTeledyne API 200A chemiluminescent NOx analyser','ChemiluminescenceTeledyne API 200E chemiluminescent NOx analyser','ChemiluminescenceTeledyne API 200E chemiluminescent NOx analyser', \
                           'ChemiluminescenceColumbia Scientific Industries Models 1600','ChemiluminescenceColumbia Scientific Industries Models 1600', \
                           'ChemiluminescenceAirpointer','ChemiluminescenceAirpointer', \
                           'ChemiluminescenceBENDIX','ChemiluminescenceBendix/Combustion Engineering Model 8101-C Oxides of Nitrogen Analyze', \
                           'ChemiluminescenceECO PHYSICS CLD 700','ChemiluminescenceECO PHYSICS CLD 700','ChemiluminescenceECO PHYSICS CLD 700 AL','ChemiluminescenceECO PHYSICS CLD 700 AL', \
                           'ChemiluminescenceS-5012','ChemiluminescenceS-5012', \
                           'Chemiluminescenceserinus 40 Nox', \
                           'ChemiluminescenceEC9843', \
                            #CAN'T FIND LITERATURE
                           'ChemiluminescenceORION','ChemiluminescenceORION','Beta ray attenuationORION', \
                           'ChemiluminescenceTECAN CLD 502','ChemiluminescenceTECAN CLD 502', \
                           'ChemiluminescenceLAP 884','ChemiluminescenceLAP 884', \
                           'ChemiluminescenceANALYSIS AUTOMATION','ChemiluminescenceANALYSIS AUTOMATION Mod. 447','ChemiluminescenceANALYSIS AUTOMATION Mod. 447', \
                           'ChemiluminescenceDASIBI 2108 NOx analyser','ChemiluminescenceDASIBI 2108 NOx analyser', \
                           'ChemiluminescenceMCV 30-QL','ChemiluminescenceMCV 30-QL', \
                           'ChemiluminescenceMELOY S1600', \
                           'ChemiluminescenceAMBIRACK', \
                            #UNKNOWN INSTRUMENT
                           'Chemiluminescencetoo generic','ChemiluminescenceUNKNOWN','ChemiluminescenceUNKNOWN','ChemiluminescenceUNKNOWN','Chemiluminescencetoo generic','chemiluminescenceUNKNOWN','ChemiluminescenceUNKNOWN','chemiluminescenceUNKNOWN','Chemiluminescencetoo generic', \
                            #method inconsistent with instruments, but trust method
                           'ChemiluminescenceTeledyne API 100A UV Fluorescent SO2 Analyser','ChemiluminescenceTeledyne API 100A UV Fluorescent SO2 Analyser', \
                           'ChemiluminescenceTeledyne API 400 UV photometric O3 analyser','ChemiluminescenceTeledyne API 400 UV photometric O3 analyser', \
                           'ChemiluminescenceThermo model 43 SO2 analyser', \
                           'ChemiluminescenceThermo model 48 CO analyser', \
                           'ChemiluminescenceEnvironnement S.A. Model CO12M CO Analyzer', \
                           'ChemiluminescenceEnvironnement S.A. Model AF21M SO2 Analyzer', \
                           'ChemiluminescenceEnvironnement S.A. Model AF21M SO2 Analyzer', \
                           'ChemiluminescenceEnvironnement S.A. Model AF22M SO2 Analyzer', \
                           'ChemiluminescenceMonitor Labs model 9850 SO2 analyser','ChemiluminescenceMonitor Labs model 9850 SO2 analyser', \
                           'ChemiluminescencePHILIPS K50109/00 Gas Filter Correlation CO analyser','ChemiluminescencePHILIPS K50109/00 Gas Filter Correlation CO analyser','ChemiluminescencePHILIPS 42','ChemiluminescencePHILIPS K50034 API 200A','ChemiluminescencePHILIPS K50034 API 200A','ChemiluminescencePHILIPS K50102 NO','ChemiluminescencePHILIPS K50235/00 NO-NOx-NO2 analyser','ChemiluminescencePHILIPS K50235/00 NO-NOx-NO2 analyser',\
                           'Beta ray attenuationTeledyne API 200A chemiluminescent NOx analyser']
        
        cld_photo_mm_list = []
        cld_eth_mm_list = []
        poti_mm_list = []
        ecc_mm_list = []
            
        ndir_mm_list = ['Non-dispersive infrared spectroscopy (NDIR)Meloy Model SA 700 Fluorescence Sulfur Dioxide Analyze','Non-dispersive infrared spectroscopy (NDIR)Monitor Labs model 9830B CO analyser','Non-dispersive infrared spectroscopy (NDIR)Monitor Labs model 8831 CO analyser','Non-dispersive infrared spectroscopy (NDIR)Thermo model 48 CO analyser','Non-dispersive infrared spectroscopy (NDIR)ORION','Non-dispersive infrared spectroscopy (NDIR)Teledyne API 200A chemiluminescent NOx analyser','Non-dispersive infrared spectroscopy (NDIR)ANALYSIS AUTOMATION','Non-dispersive infrared spectroscopy (NDIR)THERMO ELECTRON INSTRUMENTS', \
                        'Non-dispersive infrared spectroscopy (NDIR)Thermo model 43a SO2 analyser','Non-dispersive infrared spectroscopy (NDIR)Monitor Labs model 8830 CO analyser','Non-dispersive infrared spectroscopy (NDIR)CO ANALAYZER','Non-dispersive infrared spectroscopy (NDIR)Environnement S.A. Model CO12M CO Analyzer','Non-dispersive infrared spectroscopy (NDIR)Thermo model 48i CO analyser','Non-dispersive infrared spectroscopy (NDIR)too generic','Non-dispersive infrared spectroscopy (NDIR)PHILIPS K50093 API 300A','Non-dispersive infrared spectroscopy (NDIR)MLU','Non-dispersive infrared spectroscopy (NDIR)Horiba model APMA 300 CO analyser', \
                        'Non-dispersive infrared spectroscopy (NDIR)MLU 300','Non-dispersive infrared spectroscopy (NDIR)UNKNOWN','Non-dispersive infrared spectroscopy (NDIR)ENVIRONMENT','Non-dispersive infrared spectroscopy (NDIR)Teledyne API 300 gas filter correlation CO analyser','Non-dispersive infrared spectroscopy (NDIR)Thermo model 49 O3 analyser','Non-dispersive infrared spectroscopy (NDIR)Thermo model 48w CO analyser','Non-dispersive infrared spectroscopy (NDIR)Maihak Unor 6N','Non-dispersive infrared spectroscopy (NDIR)Horiba model APMA 360E CO analyser','Non-dispersive infrared spectroscopy (NDIR)Monitor Labs undetermined', \
                        'Non-dispersive infrared spectroscopy (NDIR)Teledyne API 300E gas filter correlation CO analyser','Non-dispersive infrared spectroscopy (NDIR)Teledyne API 100 UV Fluorescent SO2 Analyser','Non-dispersive infrared spectroscopy (NDIR)Environnement S.A. Model CO10M CO Analyzer','Non-dispersive infrared spectroscopy (NDIR)Horiba model APMA 350 CO analyser','Non-dispersive infrared spectroscopy (NDIR)FUJI ZRC','Non-dispersive infrared spectroscopy (NDIR)Teledyne API undertermined','Non-dispersive infrared spectroscopy (NDIR)S-5006','Non-dispersive infrared spectroscopy (NDIR)Horiba model APMA 350E CO analyser', \
                        'Non-dispersive infrared spectroscopy (NDIR)Thermo model 48c CO analyser','Non-dispersive infrared spectroscopy (NDIR)Thermo model 42 NO/Nox analyser','Non-dispersive infrared spectroscopy (NDIR)SFI CO12M','Non-dispersive infrared spectroscopy (NDIR)Horiba model APMA 360CE CO analyser','Non-dispersive infrared spectroscopy (NDIR)PHILIPS 48','Non-dispersive infrared spectroscopy (NDIR)DASIBI 3008 CO analyser','Non-dispersive infrared spectroscopy (NDIR)Teledyne API 300A gas filter correlation CO analyser','Non-dispersive infrared spectroscopy (NDIR)Horiba model APMA 370 CO analyser','Non-dispersive infrared spectroscopy (NDIR)Environnement S.A. Model CO11M CO Analyzer', \
                        'Non-dispersive infrared spectroscopy (NDIR)Horiba model APMA 360 CO analyser','Non-dispersive infrared spectroscopy (NDIR)Monitor Labs model 9841A NOx analyser','Non-dispersive infrared spectroscopy (NDIR)AAL 407','Non-dispersive infrared spectroscopy (NDIR)AMBIRACK','Non-dispersive infrared spectroscopy (NDIR)Monitor Labs model 9830 CO analyser','Non-dispersive infrared spectroscopy (NDIR)Horiba model APMA 300E CO analyser','Non-dispersive infrared spectroscopy (NDIR)PHILIPS K50109/00 Gas Filter Correlation CO analyser','UNKNOWNTeledyne API 300 gas filter correlation CO analyser','UNKNOWNHoriba model APMA 350 CO analyser', \
                        'Infrared gas filter correlationTHERMO ELECTRON INSTRUMENTS 48c','Infrared gas filter correlationHoriba model APMA 360 CO analyser','infrared absorptionUNKNOWN','Infrared gas filter correlationUNKNOWN','Infrared gas filter correlationTeledyne API 300E gas filter correlation CO analyser']
    
        doas_mm_list = ['Differential Optical Absorption Spectroscopy (DOAS)Opsis AR500 Open path monitor','Differential Optical Absorption Spectroscopy (DOAS)UNKNOWN','Ultraviolet (UV) photometryDOAS','Differential Optical Absorption Spectroscopy (DOAS)Environnement S.A. SANOA Multigas Longpath Monitoring System','Differential Optical Absorption Spectroscopy (DOAS)Environnement S.A. SANOA Multigas Longpath Monitoring System']
        crs_mm_list = []
        capss_mm_list = []
        gcfid_mm_list = ['Gas Chromatography (ref)UNKNOWN','chromatographyUNKNOWN','Gas chromatography followed by flame ionization detection (GUNKNOWN','Gas chromatography followed by flame ionization detection (GEnvironnement VOC71M','chromatographyMonitor Labs model 8440 NOx analyser','Gas chromatography (GC) + flame ionisation (GC-FID)UNKNOWN','Gas chromatography followed by flame ionization detection (GAIRMOZONE','Gas chromatography followed by flame ionization detection (GVarian Chrompack','Gas chromotography (GC)UNKNOWN','chromatographyChrompack BTX CP7001 Monitor']
        gcrgd_mm_list = []
        gcms_mm_list = ['Gas chromatography + mass spectrometry (GC-MS)AF 20 M','GAS CHROMATOGRAPHY - MASS SPECTROMETRYUNKNOWN','Gas chromatography + mass spectrometry (GC-MS)UNKNOWN','Gas chromatography + mass spectrometry GC-MS after solvent oMarkes Thermal Desorber + Agilent gas Chromatograph Mass Spectrometer']
        gcpid_mm_list = ['Gas chromatography with photo ionization detectorSYNTECH SPECTRAS GC 955 series undetermined','Gas chromatography with photo ionization detectorUNKNOWN']
        gcftir_mm_list = []
        fid_mm_list = ['Flame ionization detection (FID)Chrompack CP9000']
        flamephoto_mm_list = ['flame photometryThermo model 48 CO analyser','flame photometryTeledyne API 300 gas filter correlation CO analyser']
        ptrms_mm_list = []
        ionchrom_mm_list = ['Ion chromatographyIMPREGNATED FILTER']
        spectrometry_mm_list = ['SpectrometryBUBBLER 24 H','Atomic absorption spectrometry (AAS)UNKNOWN']
        spectrophoto_mm_list = ['SpectrophotometryIMPREGNATED FILTER','SpectrophotometrySequential Air Sampler, Type SS2000. NaI-impregnated glass sinters','SpectrophotometryGlass tubes','Spectrophotometryglass_sinter','SpectrophotometryUNKNOWN',"SpectrophotometryLipinski's aspirator",'SpectrophotometryBUBBLER 24 H','Spectrophotometryglass filter','spectrophotometryUNKNOWN','Spectrophotometryfilter pack','TGS-ANSAFILTER']
        color_mm_list = ["Griess-Saltzman reactionLipinski's aspirator",'Griess-Saltzman reaction101','Griess-Saltzman reactionUNKNOWN',"UNKNOWNLipinski's aspirator",'Griess-Saltzman reactionBUBBLER 24 H',"Griess-Saltzman reactionLipinski's aspirator AGT24",'Griess-Saltzman reactionfilter pack','NEDA Griess-Yloswayaspirator','colorimetryUNKNOWN','Spectrophotometryphotocolorimeter']
        fia_mm_list = []
        coul_mm_list = ['coulometryUNKNOWN']
        diffsamp_mm_list = ['diffusive samplerUNKNOWN','UNKNOWNSEQUENTIAL SAMPLER']
        notclear_mm_list = ['Beta ray attenuationMLU',]

#---------------------------------------------------------------------------------------------------
    elif process_group == 'EPA AQS':
        continuous_st_list = ['continuous']
        remote_st_list = []
        flask_st_list = []
        filter_st_list = []
    
        uv_mm_list = ['INSTRUMENTAL-ULTRAVIOLETABSORPTION','INSTRUMENTAL-ULTRAVIOLET2BMODEL202','INSTRUMENTAL-UVPHOTOMETRIC','INSTRUMENTAL-ULTRAVIOLETRADIATIONABSORBTN','INSTRUMENTAL-ULTRAVIOLET','INSTRUMENTAL-ULTRAVIOLETPHOTOMETRY','INSTRUMENTAL-UVABSORPTIONPHOTOMETRY/UV2BMODEL202AND205','INSTRUMENTAL-ECOTECHSERINUS10','Instrumental-UltravioletAbsorption','Instrumental-UltraVioletPhotometry','Instrumental-UVabsorptionphotometry/UV2Bmodel202and205','Instrumental-EcotechSerinus10','Instrumental-UltraViolet2BModel202']
        vuf_mm_list = []
        cld_mol_mm_list = ['INSTRUMENTAL-CHEMILUMINESCENCE','INSTRUMENTAL-GASPHASECHEMILUMINESCENCE','LOWLEVELNOXINSTRUMENTAL-TECO42SCHEMILUMINESCENCE','INSTRUMENTAL-GAS-PHASECHEMILUMINESCENCE','INSTRUMENTAL-CHEMILUMINESCENCETHERMOELECTRON42C-TL,42I-TL','INSTRUMENTAL-CHEMILUMINESCENCEECOTECHEC9841T','INSTRUMENTAL-CHEMILUMINESCENCEECOTECHEC9843','INSTRUMENTAL-CHEMILUMINESCENCETELEDYNEAPI200EU/501','INSTRUMENTAL-CHEMILUMINESCENCETHERMOELECTRON42C-Y,42I-Y','INSTRUMENTAL-CHEMILUMINESCENCEAPIMODEL265EANDT265','Instrumental-ChemiluminescenceAPIModel265EandT265','Instrumental-ChemiluminescenceThermoElectron42C-TL,42i-TL','Instrumental-ChemiluminescenceEcotechEC9841T','Instrumental-ChemiluminescenceTeledyneAPI200EU/501','Instrumental-Chemiluminescence','Instrumental-ChemiluminescenceEcotechEC9843','Instrumental-GasPhaseChemiluminescence','Instrumental-ChemiluminescenceEcotechEC9841T','Instrumental-ChemiluminescenceTeledyneAPI200EU/501','Instrumental-ChemiluminescenceThermoElectron42C-Y,42i-Y','Instrumental-Chemiluminescence']
        cld_photo_mm_list = ['TELEDYNE-APIMODEL200EUPORT200UP-PHOTOLYTIC-CHEMILUMINESCENCE','INSTRUMENTAL-CHEMILUMINESCENCETELEDYNEAPIT200UPPHOTOLYTIC','Instrumental-ChemiluminescenceTeledyneAPIT200UPphotolytic','Teledyne-APIModel200EUPorT200UP-Photolytic-Chemiluminescence','Teledyne-APIModel200EUPorT200UP-Photolytic-Chemiluminescence']
        cld_eth_mm_list = []
        poti_mm_list = []
        ecc_mm_list = []
        ndir_mm_list = ['INSTRUMENTAL-NONDISPERSIVEINFRARED','INSTRUMENTAL-NONDISPERSIVEINFRAREDPHOTOMETRY','Instrumental-NondispersiveInfreredPhotometry','INSTRUMENTAL-GasFilterCorrelationEcotechEC9830T','INSTRUMENTAL-GasFilterCorrelationTeledyneAPI300EU','INSTRUMENTAL-GASFILTERCORRELATIONCOANALYZER','INSTRUMENTAL-GasFilterCorrelationThermoElectron48C-TL','INSTRUMENTAL-GasFilterCorrelationThermoElectron48i-TLE','INSTRUMENTAL-DUALISOTOPEFLORESCENCE']
        doas_mm_list = ['INSTRUMENTAL-OPENPATHO3ANALYZER','INSTRUMENTAL-OPENPATHNOANALYZER']
        crs_mm_list = []
        capss_mm_list = ['TELEDYNEMODELT500U-CAVITYATTENUATEDPHASESHIFTSPECTROSCOPY','TeledyneModelT500U-CavityAttenuatedPhaseShiftSpectroscopy']
        gcfid_mm_list = ['gas chromatography flame ionisation detection','gas chromatography','gas chromatography flame ionisation detection/mass selective detection']
        gcrgd_mm_list = []
        gcms_mm_list = ['gas chromatography mass spectrometry','gas chromatography mass spectrometry/flame ionisation detection']
        gcpid_mm_list = []
        gcftir_mm_list = ['gas chromatography fourier transform infrared spectroscopy/mass spectrometry']
        fid_mm_list = ['flame ionisation detection']
        flamephoto_mm_list = []
        ptrms_mm_list = []
        ionchrom_mm_list = []
        spectrometry_mm_list = []
        spectrophoto_mm_list = ['INSTRUMENTAL-CHEMILUMINESCENCERHODAMINEBDYE','INSTRUMENTAL-CHEMILUMINESCENCERHODAMINEBDYE']
        color_mm_list = ['INSTRUMENTAL-COLORIMETRIC-GRIESS-SALTZMAN','INSTRUMENTAL-COLORIMETRIC','INSTRUMENTAL-COLORIMETRIC-LYSHKOW(MOD)']
        fia_mm_list = []
        coul_mm_list = ['INSTRUMENTAL-COULOMETRIC']
        diffsamp_mm_list = []
        notclear_mm_list = []
#---------------------------------------------------------------------------------------------------  
    elif process_group == 'EMEP':
        continuous_st_list = ['continuous']
        remote_st_list = []
        flask_st_list = []
        filter_st_list = []
        
        uv_mm_list = ['uv_absuv_abs_35DE03L_uv_abs', 'uv_absuv_12DK01L_uv_abs', 'uv_absuv_abs_1RU01L_uv_abs', 'uv_absuv_abs_02AT03L_ozone1a', 'uv_absuva37FI01L_uv_absorption', 'uv_absuv_abs_3DE03L_uv_abs', 'uv_absuv_abs_18RU01L_uv_abs', \
                      'uv_absuv_abs_6DE03L_uv_abs', 'uv_absuv_abs_14DE03L_uv_abs', 'uv_absuv_abs_uk_0039NO01L_uv_abs', 'uv_absuv_mon_17ES04L_uv_abs', 'uv_absuv_16LV01L_uv_abs', 'uv_absuv_mon_11ES04L_uv_abs', 'uv_absuv_09NL01L_uv_abs', \
                      'uv_absuv_abs_34SE01L_uv_abs', 'uv_absuv_abs_5CZ01L_97_3_3_1j', 'uv_absuv_abs_31DE03L_uv_abs', 'uv_absuv_mon_14ES04L_uv_abs', 'uv_absuv_3PL01L_uv_abs', 'uv_absuv_abs_03CZ01L_uv_abs', 'uv_absuv_abs_uk_0044NO01L_uv_abs', \
                      'uv_absuv_15GB02L_uv_abs', 'uv_absuv_abs_7DE07L_uv_abs', 'uv_absuv_abs_uk_0762NO01L_uv_abs', 'uv_absuv_abs_1CZ01L_uv_abs', 'uv_absO3mon_420CA01L_UV_absorption', 'uv_absuv_7SK01L_uv_abs', 'uv_absuv_abs_41DK01L_uv_abs', \
                      'uv_absuv_01IE01L_uv_abs', 'uv_absuv_6ES04L_uv_abs', 'uv_absuv_5SK01L_uv_abs', 'uv_absuv_abs_17DE03L_uv_abs', 'uv_absuv_abs_08DE08L_uv_abs', 'uv_absTEI49i_061991500FI01L_uv_absorption', 'uv_absuv_abs_7DE03L_uv_abs', \
                      'uv_absuv_mon_13ES04L_uv_abs', 'uv_absuv_abs_03AT03L_ozone1a', 'uv_absuv_2GB02L_uv_abs', 'uv_absuv_abs_31SE01L_uv_abs', 'uv_absuv_abs_uk_0002NO01L_uv_abs', 'uv_absuv_abs_34AT03L_ozone1a', 'uv_absuv_mon_15ES04L_uv_abs', \
                      'uv_absuv_abs_11EE01L_uvabs_AH', 'uv_absuv_abs_02AM01L_uv_abs', 'uv_absuv_abs_40AT03L_ozone1a', 'uv_absES0001R_CM08110023ES04L_uv_abs', 'uv_absuv_abs_uk_0047NO01L_uv_abs', 'uv_absuv_abs_42AT03L_ozone1a', 'uv_absuv_mon_8ES04L_uv_abs', \
                      'uv_absuv_abs_11DE03L_uv_abs', 'uv_absuv_abs_uk_0041NO01L_uv_abs', 'uv_absuv_abs_uk_0488NO01L_uv_abs', 'uv_absO3mon10DK01L_UV_absorption', 'uv_absuv_abs_08SI01L_uv_abs', 'uv_absuv_49GB02L_uv_abs', 'uv_absuv_abs_03CZ01L_97_3_3_1j', \
                      'uv_absuv_abs_3GR04L_uv_abs', 'uv_absThermoUV_LabsMK01L_uv_abs', 'uv_absuv_abs_2CH01L_uv_abs', 'uv_absuv_abs_32SE01L_uv_abs', 'uv_absuv_abs_05CH01L_uv_abs', 'uv_absuv_abs_46DE03L_uv_abs', 'uv_absTEI_49C-54544-300SI01L_ozone', \
                      'uv_absuv_abs_15FR02L_uv_abs', 'uv_absuv_36GB02L_uv_abs', 'uv_absuv_abs_44AT03L_ozone1a', 'uv_absuv_abs_13SE01L_uv_abs', 'uv_absuv_abs_38AT03L_ozone1a', 'uv_absuv_abs_13RU01L_uv_abs', 'uv_absuv_abs_37AT03L_ozone1a', 'uv_absuv_1IT01L_uv_abs', \
                      'uv_absuv_abs_08RO02L_uv_abs', 'uv_absuv_50GB02L_uv_abs', 'uv_absuv_abs_uk_1083NO01L_uv_abs', 'uv_absuv_abs_47DE03L_uv_abs', 'uv_absuv_abs_45AT03L_ozone1a', 'uv_absuv_abs_uk_0056_15mNO01L_uv_abs', 'uv_absuv_abs_04CH01L_uv_abs', 'uv_absuv_abs_41AT03L_ozone1a', \
                      'uv_absuv_abs_uk_0048NO01L_uv_abs', 'uv_absuv_abs_10FR02L_uv_abs', 'uv_absuv_mon_1ES04L_uv_abs', 'uv_absuv_abs_32FR02L_uv_abs', 'uv_absuv_4SK01L_uv_abs', 'uv_absuv_abs_uk_1201NO01L_uv_abs', 'uv_absuv_abs_35SE01L_uv_abs', 'uv_absuv_abs_09EE01L_uvabs', \
                      'uv_absuv_abs_32SI01L_uv_abs', 'uv_absuv_abs_12DE03L_uv_abs', 'uv_absuv_2SK01L_uv_abs', 'uv_absuv_abs_1DE03L_uv_abs', 'uv_absuv_abs_39DE03L_uv_abs', 'uv_absuv_mon_4ES04L_uv_abs', 'uv_absuv_abs_2SE01L_uv_abs', 'uv_absTEI94C_1CH01L_O3', 'uv_absuv_41GB02L_uv_abs', \
                      'uv_absuv_abs_05AT03L_ozone1a', 'uv_absuv_abs_18DE03L_uv_abs', 'uv_absuv_abs_DE43DE09L_uv_abs', 'uv_absuv_abs_3CZ01L_uv_abs', 'uv_absuv_31GB02L_uv_abs', 'uv_absuv_39GB02L_uv_abs', 'uv_absuv_abs_09EE01L_uvabs_AH', 'uv_absuva22FI01L_uv_absorption', 'uv_absO3_FR19FR07L_uv_abs', \
                      'uv_absuv_10LV01L_uv_abs', 'uv_absuv_abs_5SE01L_uv_abs', 'uv_absuv_abs_01CZ01L_97_3_3_1', 'uv_absO3_FR14FR07L_uv_abs', 'uv_absuv_35GB02L_uv_abs', 'uv_absES0017R_49C65415ES04L_uv_abs', 'uv_absuv_abs_5ES04L_uv_abs', 'uv_absuv_abs_1GR01L_uv_abs', 'uv_absO3_FR17FR07L_uv_abs', \
                      'uv_absuv_mon_3ES04L_uv_abs', 'uv_absO3_FR18FR07L_uv_abs', 'uv_absuv_mon_12ES04L_uv_abs', 'uv_absuv_abs_31SI01L_uv_abs', 'uv_absuv_abs_35BE02L_uv_abs', 'uv_absuv_abs_uk_1011NO01L_uv_abs', 'uv_absuv_abs_07NL01L_uv_abs', 'uv_absuv_abs_08FR02L_uv_abs', 'uv_absuv_03PL01L_uv_abs', \
                      'uv_absuv_abs_1IT01L_uv_abs', 'uv_absuv_abs_9FR02L_uv_abs', 'uv_absuv_abs_uk_1007NO01L_uv_abs', 'uv_absO3_FR08FR07L_uv_abs', 'uv_absuva17FI01L_uv_absorption', 'uv_absuv_abs_11NL01L_uv_abs', 'uv_absuv_abs_5DK01L_uv_abs', 'uv_absuvabs_49AT03L_ozone1a', 'uv_absuv_abs_uk_0977NO01L_uv_abs', \
                      'uv_absuv_52GB02L_uv_abs', 'uv_absuv_48GB02L_uv_abs', 'uv_absuv_abs_uk_0052NO01L_uv_abs', 'uv_absuv_abs_18FR02L_uv_abs', 'uv_absuv_38GB02L_uv_abs', 'uv_absuv_abs_uk_0056_25mNO01L_uv_abs', 'uv_absuv_abs_32AT03L_ozone1a', 'uv_absuv_abs_9DE09L_uv_abs', 'uv_absuv_abs_31CH01L_uv_abs', 'uv_absuv_abs_42DE03L_uv_abs', \
                      'uv_absuv_abs_08ES04L_uv_abs', 'uv_absuv_abs_03SE01L_uv_abs', 'uv_absuv_abs_39SE01L_uv_abs', 'uv_absuv_abs_uk_0492NO01L_uv_abs', 'uv_absuv_mon_5ES04L_uv_abs', 'uv_absuv_45GB02L_uv_abs', 'uv_absO3_FR30FR07L_uv_abs', 'uv_absuv_abs_2DE02L_uv_abs', 'uv_absuv_abs_13FR02L_uv_abs', \
                      'uv_absuv_mon_2ES04L_uv_abs', 'uv_absuv_abs_uk_0001NO01L_uv_abs', 'uv_absuv_abs_02CH01L_uv_abs', 'uv_absuv_abs_09NL01L_uv_abs', 'uv_absuv_abs_01MT01L_uv_abs', 'uv_absuv_abs_11SE01L_uv_abs', 'uv_absuv_51GB02L_uv_abs', 'uv_absuv_mon_16ES04L_uv_abs', 'uv_absuv_abs_4IT03L_uv_abs', 'uv_absuv_abs_12FR02L_uv_abs', \
                      'uv_absuv_abs_8DE03L_uv_abs', 'uv_absuv_abs_9DE03L_uv_abs', 'uv_absuv_abs_30AT03L_ozone1a', 'uv_absuv_09_1DK01L_uv_abs', 'uv_absuv_43GB02L_uv_abs', 'uv_abs03_FR30FR04L_uv_abs', 'uv_absuv_10NL01L_uv_abs', 'uv_absuv_abs_2GR02L_uv_abs', 'uv_absO3_FR13FR07L_uv_abs', 'uv_absuv_2CY01L_uv_abs', 'uv_absuv_53GB02L_uv_abs', \
                      'uv_absuv_abs_46AT03L_ozone1a', 'uv_absuva4FI01L_uv_absorption', 'uv_absuv_abs_26DE03L_uv_abs', 'uv_absuv_2PL01L_uv_abs', 'uv_absuv_abs_8_7mFR02L_uv_abs', 'uv_absuv_33GB02L_uv_abs', 'uv_absuv_abs_uk_0058NO01L_uv_abs', 'uv_absuv_abs_12SE01L_uv_abs', 'uv_absuv_abs_48AT03L_ozone1a', 'uv_absuv_abs_16FR02L_uv_abs', 'uv_absuv_abs_31DK01L_uv_abs', \
                      'uv_absuv_abs_01CZ01L_97_3_3_1j', 'uv_absuv_abs_33SE01L_uv_abs', 'uv_absuv_44GB02L_uv_abs', 'uv_absuv_abs_14FR02L_uv_abs', 'uv_absuv_abs_07MK01L_uv_abs', 'uv_absuv_abs_uk_0042NO01L_uv_abs', 'uv_absuv_mon_10ES04L_uv_abs', 'uv_absuv_mon_7ES04L_uv_abs', 'uv_absuv_abs_5DE03L_uv_abs', 'uv_absuv_abs_14SE01L_uv_abs', 'uv_absES0005R_CM08110025ES04L_uv_abs', \
                      'uv_absuv_abs_uk_0055NO01L_uv_abs', 'uv_absuv_abs_uk_0030NO01L_uv_abs', 'uv_absuv_abs_33AT03L_ozone1a', 'uv_absuv_abs_01CH01L_uv_abs', 'uv_absuv_6SK01L_uv_abs', 'uv_absuv_abs_31FR02L_uv_abs', 'uv_absuv_14GB02L_uv_abs', 'uv_abs03_FR19FR05L_uv_abs', 'uv_absuv_abs_2DE03L_uv_abs', 'uv_absuv_abs_47AT03L_ozone1a', 'uv_absuv_abs_uk_0489NO01L_uv_abs', \
                      'uv_absuv_abs_17FR02L_uv_abs', 'uv_absuv_abs_uk_0015NO01L_uv_abs', 'uv_absuv_abs_uk_0043NO01L_uv_abs', 'uv_absuv_abs_15LT01L_uv_abs', 'uv_absuv_02PL01L_uv_abs', 'uv_absuv_abs_1DE01L_uv_abs', 'uv_absuv_abs_2GR03L_uv_abs', 'uv_absuv_abs_01BE02L_uv_abs', 'uv_absuv_abs_91NL01L_uv_abs', 'uv_absuv_abs_4DE03L_uv_abs', 'uv_absAPI_400ASI02L_ozone', \
                      'uv_absuv_abs_08DE03L_uv_abs', 'uv_absuv_abs_03CH01L_uv_abs', 'uv_absuv_abs_4PT01L_uv_abs', 'uv_absuv_abs_38DE03L_uv_abs', 'uv_absuv_abs_32BE02L_uv_abs', 'uv_absuv_13GB02L_uv_abs', 'uv_absuv_abs_uk_0045NO01L_uv_abs', 'uv_absuv_abs_32DK01L_uv_abs', 'uv_absuv_abs_4IT04L_uv_abs', 'uv_absuv_4PL01L_uv_abs', 'uv_absuv_abs_13DE03L_uv_abs', \
                      'uv_absuva09FI01L_uv_absorption', 'uv_absuv_5PL02L_uv_abs', 'uv_absuv_abs_04AT03L_ozone1a', 'uv_absuv_abs_16RU01L_uv_abs', 'uv_absuv_abs_43AT03L_ozone1a', 'uv_absuv_abs_11FR02L_uv_abs', 'uv_absuv_abs_33SI01L_uv_abs', 'uv_absO3_FR09FR07L_uv_abs', 'uv_absuv_abs_2HU01L_uv_abs', 'uv_absuv_mon_9ES04L_uv_abs', 'uv_absuv_abs_45DE03L_uv_abs', \
                      'uv_absO3_FR10FR07L_uv_abs', 'uv_absuv_6GB02L_uv_abs', 'uv_absuv_abs_31GB02L_uv_abs', 'uv_absuv_37GB02L_uv_abs', 'uv_absO3_FR15FR07L_uv_abs', 'uv_absES0008R_CM10060018ES04L_uv_abs', 'uv_absO3_FR16FR07L_uv_abs','uv_absEE01_ua1EE01L_uvabs_AH', 'uv_absCM08110024ES04L_uv_abs', 'uv_absO341MLT01L_ozone1a', 'uv_absuv_abs_h_0001BE01L_uv_abs', \
                      'uv_absuv_abs_644NL01L_uv_abs', 'uv_absCM10060019ES04L_uv_abs', 'uv_absuv_abs_uk_0059NO01L_uv_abs', 'uv_absuv_abs_h_0035BE01L_uv_abs', 'uv_absES0017R_CM10060017ES04L_uv_abs', 'uv_absCM10060020ES04L_uv_abs', 'uv_absCM08110026ES04L_uv_abs', 'uv_absuv_abs_50AT03L_ozone1a', 'uv_absuv_abs_BKTCH01L_uv_abs', 'uv_absuv_abs_53BG01L_uv_abs', \
                      'uv_abs49C-65415-348ES04L_uv_abs', 'uv_absuv_abs_33SI02L_uv_abs', 'uv_absO3MLBPL02L_AO3MONLABS', 'uv_absTEI49I_3CH01L_O3', 'uv_absuv_abs_h_0032BE01L_uv_abs', 'uv_absuv_IE0031GB02L_uv_abs', 'uv_absuv_absNL01L_uv_abs', 'uv_absuv_abs_MKNCH01L_uv_abs', 'uv_absTEI49I_4CH01L_O3', 'uv_absCM10060017ES04L_uv_abs', 'uv_absTEI49I_1CH01L_O3', \
                      'uv_absuv_abs_03RO02L_ub_abs', 'uv_absCM08110023ES04L_uv_abs', 'uv_absTHERMO49CIT04L_UVABS_2', 'uv_absTEI49I_5CH01L_O3', 'uv_absuv_abs_kre_0050CZ06L_uv_abs', 'uv_absCM08110022ES04L_uv_abs', 'uv_absThermo49iPL01L_uv_abs', 'uv_absTEI49I_2CH01L_O3', 'uv_abs49C-65422-348ES04L_uv_abs', 'uv_abs332503151ES04L_uv_abs', \
                      'uv_absTE49WNL01L_LM-W-401-400-402_IL-W-900', 'uv_absuv_abs_10NL01L_uv_abs', 'uv_absO301GB03L_O3_1', 'uv_absO3mon5DK01L_UV_absorption', 'uv_absEE01_ua2EE01L_uvabs_AH', 'uv_absO3mon1DK01L_UV_absorption', 'uv_absA-102ES04L_uv_abs', 'uv_abs56_fluks_15mNO01L_uv_abs']
        
        vuf_mm_list = ['uv_fluorescAL_5001-DE03LDE01L_VUV_Fluorescence','uv_fluorescTE48S_DE43DE09L_uv_fluoroesc','uv_fluorescAL_5001_DE01LDE01L_VUV_Fluorescene','uv_fluorescTEI_48C_TL_MKNCH01L_uv_fluoresc']
        
        cld_mol_mm_list =  ['chemiluminescence_photometerTE42W_91NL01L_LLO_MM_141','chemiluminescence_photometerTE42W_09_dupNL01L_LLO_MM_141','chemiluminescence_photometerThermo_42i-LS_TAD','chemiluminescence_photometerTE42W_09NL01L_LVM_LU_P141','chemiluminescence_photometerEcophysics_CLD770ALppt_1DE09L_chemiluminescence','chemiluminescence_photometerAPI200ENL01L_LVM_LU_P132_133_142','chemiluminescence_photometerTE42W_91NL01L_LVM_LU_P141', \
                            'chemiluminescence_photometerThermo_42i-LS_PYEFR07L_chemiluminescence_Thermo_42i-LS','chemiluminescence_photometerAPI200ENL01L_LM-W-300-301-302-303-304_IL-W-900','chemiluminescence_photometerThermo_M_42SDCZ01L_chemiluminescence_Thermo_42SD','chemiluminescence_photometerECO_Physics_DE43DE09L_no_o3_chemoluminescene','chemiluminescence_photometerTE42W_09NL01L_LLO_MM_141', \
                            'chemiluminescence_photometerTE42WNL01L_LVM_LU_P132_133_141','chemiluminescence_photometerECO_Physics_DE43','chemiluminescence_photometer42ITL-1034346046ES04L_chemilum','chemiluminescence_photometerTE42W_09_dupNL01L_LVM_LU_P141','chemiluminescence_photometer42CTL-0328902535ES04L_chemilum','chemiluminescence_photometer42CTL-0328902535ES04L_chemilum', \
                            'chemiluminescence_photometerThermo_Scientific_Nox_Analyzer_42iSI01L_Mo_converter_chemoluminis_detection','chemiluminescence_photometerTE42W_09NL01L_Chemilum','chemiluminescence_photometerThermo42CMK01L_Nitrogen',  'chemiluminescence_photometerTEI42ITL_Mo_4CH01L_no2', 'chemiluminescence_photometerapna370_Mo_3CH01L_no2','chemiluminescence_photometerapna360_3CH01L_no2', \
                            'chemiluminescence_photometertei42ctl_4CH01L_no2','chemiluminescence_photometertei42ctl_4CH01L_chemilum', \
                            #NO LITERATURE
                            'chemiluminescence_photometerchemilum_08ES04L_chemilum', 'chemiluminescence_photometerchemilum_15GB02L_chemilum', 'chemiluminescence_photometerchemilum_91NL01L_chemilum', 'chemiluminescence_photometerES0005R_1036546009ES04L_chemilum', 'chemiluminescence_photometerchemilum_5ES04L_chemilum', \
                            'chemiluminescence_photometerchemil_11BE01L_AH', 'chemiluminescence_photometerchemilum_15FR02L_chemilum', 'chemiluminescence_photometerchemilum_14ES04L_chemilum', 'chemiluminescence_photometerchemilum_45GB02L_chemilum', 'chemiluminescence_photometerchemilum_09NL01L_chemilum', 'chemiluminescence_photometerchemilum_11ES04L_chemilum', \
                            'chemiluminescence_photometerchemilum_33GB02L_chemilum', 'chemiluminescence_photometerchemilum_36GB02L_chemilum', 'chemiluminescence_photometerchemil_13BE01L_AH','chemiluminescence_photometerchemilum_39GB02L_chemilum', 'chemiluminescence_photometer44N012BE01L_AH', \
                            'chemiluminescence_photometerchemilum_91', 'chemiluminescence_photometer1036546010ES04L_chemilum', 'chemiluminescence_photometerchemilum_44GB02L_chemilum', 'chemiluminescence_photometerES0008R_1036546008ES04L_chemilum', \
                            'chemiluminescence_photometerchemilum_01ES04L_chemilum', 'chemiluminescence_photometerchemilum_9', 'chemiluminescence_photometerchemilum_13GB02L_chemilum', 'chemiluminescence_photometerchemilum_43GB02L_chemilum', \
                            'chemiluminescence_photometerchemilum_9ES04L_chemilum', 'chemiluminescence_photometerchemilum_7_dupES04L_chemilum', 'chemiluminescence_photometerchemilum_51GB02L_chemilum', 'chemiluminescence_photometerchemilum_12ES04L_chemilum', 'chemiluminescence_photometerchem_3', 'chemiluminescence_photometerchemilum_13ES04L_chemilum', \
                            'chemiluminescence_photometerchemilum_53GB02L_chemilum', 'chemiluminescence_photometerchem_3CZ01L_97_3_3', 'chemiluminescence_photometerNO12_dupBE01L_AH', 'chemiluminescence_photometer33290149ES04L_chemilum', 'chemiluminescence_photometerchemilum_2GB02L_chemilum', 'chemiluminescence_photometerchemilum_38GB02L_chemilum', 'chemiluminescence_photometerNO12BE01L_AH', \
                            'chemiluminescence_photometerchemilum_13BE01L_chemilum', 'chemiluminescence_photometerchemilum_16ES04L_chemilum', 'chemiluminescence_photometerchemilum_37GB02L_chemilum', 'chemiluminescence_photometerchemilum_31GB02L_chemilum', 'chemiluminescence_photometerchemilum_13FR02L_chemilum', \
                            'chemiluminescence_photometerchemilum_11BE01L_chemilum', 'chemiluminescence_photometerchemilum_15FR07L_chemilum', 'chemiluminescence_photometerchemilum_50GB02L_chemilum', 'chemiluminescence_photometer328902538ES04L_chemilum', \
                            'chemiluminescence_photometerchemilum_8ES04L_chemilum', 'chemiluminescence_photometerNO29BE01L_AH', 'chemiluminescence_photometerchemilum_644NL01L_chemilum', 'chemiluminescence_photometerAN02GB02L_AN_L', 'chemiluminescence_photometerchemilum_7ES04L_chemilum', 'chemiluminescence_photometer328902537ES04L_chemilum', 'chemiluminescence_photometerchemil_5ES04L_chemilum', \
                            'chemiluminescence_photometerchemilum_11NL01L_chemilum', 'chemiluminescence_photometerchemilum_07NL01L_chemilum', 'chemiluminescence_photometerchemilum_6ES04L_chemilum', 'chemiluminescence_photometerChemilum_02_dupCH01L_Chemilum', 'chemiluminescence_photometerchemilum_17ES04L_chemilum', 'chemiluminescence_photometerchemilum_13FR07L_chemilum', \
                            'chemiluminescence_photometerchemilum_08_7mFR02L_chemilum', 'chemiluminescence_photometerchem_01GR01L_chem01', 'chemiluminescence_photometerAN02GB02L_GB02_AN_L', 'chemiluminescence_photometerchemilum_10ES04L_chemilum', 'chemiluminescence_photometerchemilum_15ES04L_chemilum', 'chemiluminescence_photometerchemilum_14GB02L_chemilum', \
                            'chemiluminescence_photometerES0017R_33290149ES04L_chemilum', 'chemiluminescence_photometerchemilum_10NL01L_chemilum', 'chemiluminescence_photometer44N029BE01L_AH','chemiluminescence_photometerES0017R_33290149ES04L_chemilum', 'chemiluminescence_photometerchemilum_14GB02L_chemilum', 'chemiluminescence_photometerchemilum_08ES04L_chemilum', \
                            'chemiluminescence_photometerTIN-001ES04L_chemilum','chemiluminescence_photometerTIN-009ES04L_chemilum','chemiluminescence_photometerTIN-003ES04L_chemilum','chemiluminescence_photometerTIN-004ES04L_chemilum','chemiluminescence_photometerTIN-010ES04L_chemilum', \
                            'uv_fluorescuv_41DK01L_chemilum','uv_fluorescuv_12_1DK01L_chemilum','uv_fluorescuv_05_dupDK01L_chemilum','uv_fluorescuv_05DK01L_chemilum','uv_fluorescNOxmon1DK01L_Chemiluminescence','uv_fluorescuv_05DK01L_Chemiluminescence','uv_fluorescuv_08DK01L_Chemilum','uv_fluorescuv_08DK01L_chemilum','uv_fluorescuv_12_1DK01L_Chemilum','uv_fluorescuv_09DK01L_chemilum','uv_fluorescuv_nox_12DK01L_Chemiluminescence', \
                            'uv_fluorescuv_08DK01L_Chemiluminescence','uv_fluorescuv_05DK01L_Chemilum','uv_fluorescNOxmon2DK01L_Chemiluminescence','chemiluminescence_photometerChemilum_1PT01L_Chemilum','chemiluminescence_photometerchemilum_07MK01L_chemilum','chemiluminescence_photometerChemilum_9NL01L_Chemilum','chemiluminescence_photometerChemilum_4IT03L_Chemilum', \
                            'chemiluminescence_photometerChemilum_4IT04L_Chemilum','chemiluminescence_photometerorion_08RO02L_chemilum','chemiluminescence_photometerChemilum_8NL01L_Chemilum','chemiluminescence_photometerChemilum_7NL01L_Chemilum','chemiluminescence_photometerChemilum_2NL01L_Chemilum','chemiluminescence_photometerChemilum_10NL01L_Chemilum','chemiluminescence_photometerChemilum_1IT01L_Chemilum', \
                            'chemiluminescence_photometercl37FI01L_chemilum','uv_fluorescuv_05DK01L_chemilium','chemiluminescence_photometerchemilum_35BE02L_chemilum','chemiluminescence_photometerChemilum_03CH01L_Chemilum','chemiluminescence_photometercl1_09EE01L_chemilum','chemiluminescence_photometerChemilum_9EE01L_Chemilum','chemiluminescence_photometerml8841_4CH01L_no2', \
                            'chemiluminescence_photometerchemilum_09FI01L_chemilum','chemiluminescence_photometerchemilum_01BE02L_chemilum','chemiluminescence_photometerchemilum_37FI01L_chemilum','chemiluminescence_photometerChemilum_02CH01L_Chemilum','chemiluminescence_photometerChemilum_11EE01L_Chemilum','chemiluminescence_photometerchemilum_1ES04L_chemilum', \
                            'chemiluminescence_photometerchemilum_02CY01L_chemilum','chemiluminescence_photometerchem_11EE01L_chlum_AD','chemiluminescence_photometerchemilum_32BE02L_chemilum','chemiluminescence_photometerchemilium_02CY01L_chemilum','chemiluminescence_photometerchemilum_06ES04L_chemilum','uv_fluorescuv_05','chemiluminescence_photometerES0009R_TIN-008ES04L_chemilum', \
                            'chemiluminescence_photometerchemilum_4AT03L_chemilum','chemiluminescence_photometerchemilum_22FI01L_chemilum','chemiluminescence_photometerchem_3CZ01L_chemiluminesc', 'chemiluminescence_photometerchemilum_5AT03L_chemilum','chemiluminescence_photometercl2_11','chemiluminescence_photometerEE01_cl2EE01L_chlum_AD', \
                            'chemiluminescence_photometerchemilum_4ES04L_chemilum','chemiluminescence_photometercl22FI01L_chemilum','chemiluminescence_photometerchemilum_36_dupGB02L_chemilum','chemiluminescence_photometerchemilum_96FI01L_chemilum', 'chemiluminescence_photometerchemilum_h_0001BE01L_chemilum','chemiluminescence_photometerchemilum_3ES04L_chemilum', \
                            'chemiluminescence_photometerChemilum_05CH01L_Chemilum','chemiluminescence_photometerml8841_3CH01L_no2','chemiluminescence_photometercl2_09EE01L_chlum_AD','chemiluminescence_photometerchemilum_h_0035BE01L_chemilum','chemiluminescence_photometerchemilum_48AT03L_chemilum','chemiluminescence_photometerchemilum_h_0032BE01L_chemilum', \
                            'chemiluminescence_photometerChemilum_04CH01L_Chemilum','chemiluminescence_photometerchemilum_2AT03L_chemilum','chemiluminescence_photometerchemilum_17FI01L_chemilum','chemiluminescence_photometerchemilum_03CH01L_chemilum','chemiluminescence_photometerChemilum_11EE01L_chlum_AD'] 
        
        cld_photo_mm_list = ['chemiluminescence_photometercranox1CH01L_chemiluminesc','chemiluminescence_photometercranox_1CH01L_NOX','chemiluminescence_photometercranox_01CH01L_NOX','chemiluminescence_photometercranox_1_dupCH01L_NOX','chemiluminescence_photometercranox_1CH01L_chemiluminesc', 'chemiluminescence_photometerbluelight_2CH01L_no2','chemiluminescence_photometerECO_Physics_DE43DE09L_no2_photolytic_conversion', \
                             'chemiluminescence_photometercranox_5CH01L_chemilium','chemiluminescence_photometerTEI42ITL_BLC_2CH01L_no2','chemiluminescence_photometercranox_5','chemiluminescence_photometerbluelight5CH01L_no2','chemiluminescence_photometerTEI42ITL_BLC_5CH01L_no2','chemiluminescence_photometerbluelight_5CH01L_no2']
        
        cld_eth_mm_list = []
        poti_mm_list = []
        ecc_mm_list = []
        
        ndir_mm_list = ['ndirndir_01CH01L_CO','ndirAPMA360_1CH01L_CO', 'ndirAPMA360_01CH01L_CO','ndirHORIBA_AMPA370IT04L_NDIR_3','infrared_absorptionia1_09EE01L_uk','infrared_absorptionin_abs_3CZ01L_97_3_3','infrared_absorptionTE48WNL01L_LVM-LU', 'infrared_absorptionir_abs_51GB02L_ir_abs','infrared_absorptionin_abs_30FR02L_ir_abs','infrared_absorptionCO_FR19FR07L_infrared_absorption','infrared_absorptionCO_FR19FR07L_IR_abs','infrared_absorptionTE48WNL01L_LVM_LU', \
                        'infrared_absorptionir_abs_50GB02L_ir_abs', 'infrared_absorptionCO_FR30FR07L_IR_abs', 'infrared_absorptionTE48WNL01L_LM-W-604-605-606-607_IL-W-900', 'infrared_absorptionEE01_ia1EE01L_inabs_AH','infrared_absorptionorion_08_2RO02L_ir_abs', 'infrared_absorptionCO_FR30FR07L_infrared_absorption','ndirAPMA_360_205003SI01L_carbon_monoxide']
        
        doas_mm_list = ['doasOPSIS5']
        crs_mm_list = ['tracegas_monitorpicarro_G2401_1CH01L_co','tracegas_monitorpicarro1_G2401NO01L_picarro_2014']
        capss_mm_list = []
        
        gcfid_mm_list = ['steel_canistersc_fr8NO01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_05DE08L_GC_FID_Chrompac_VOCAIR_Analyzer', 'steel_canistersc_no42NO01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_lv10NO01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_cz99CZ01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_02DE08L_GC_FID_Chrompac_VOCAIR_Analyzer', \
        'steel_canistersc_no1NO01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_fi96FI01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_ch3CH01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'online_gconline_Varian_3600CXDE09L_GC_FID_Varian_GmbH', 'steel_canistersc_07DE08L_GC_FID_Chrompac_VOCAIR_Analyzer', \
        'steel_canistersc_03DE08L_GC_FID_Chrompac_VOCAIR_Analyzer', 'steel_canistersc_fr8FR02L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_es9ES01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_fr13FR02L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_cz3NO01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'online_gconline_gc_05CH01L_GC_FID_Chrompack_VOCAIR_Analyzer', \
        'steel_canistersc_de2NO01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_fi9FI01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_sk6NO01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_08DE08L_GC_FID_Chrompac_VOCAIR_Analyzer', 'steel_canistersc_it1IT01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_09DE08L_GC_FID_Chrompac_VOCAIR_Analyzer', \
        'steel_canistersc_cz3CZ01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_fr15FR02L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_lv10LV01L_GC_FID_Chrompack_VOCAIR_Analyzer', 'steel_canistersc_sk6SK01L_GC_FID_Chrompack_VOCAIR_Analyzer','steel_canistersteel_canister_01ES01L_gc', 'steel_canisterVOC_FR13FR07L_NMHC_analysis','steel_canister4-8021ES01L_gc', \
        'online_gcVOC01GB03L_voc_1','online_gcSyntech_Spectras_GC955_09NL01L_LM_W_700_701','steel_canisterVOC_FR15FR07L_NMHC_analysis','online_gcVOC02GB02L_voc_1']
        
        gcrgd_mm_list = ['GC-HgORGA3_no42NO01L_RGA3']
        gcms_mm_list = []
        gcpid_mm_list = []
        gcftir_mm_list = []
        fid_mm_list = []
        flamephoto_mm_list = []
        ptrms_mm_list = ['online_ptrIONICON_PTR_TOF_8000_BIRNO01L_internal_proton_transfer_reaction']
        ionchrom_mm_list = ['filter_3packf3p_d_0002NO01L_IC_aero','filter_1packf1p_02CH01L_IC_a']

        spectrometry_mm_list = []
        
        spectrophoto_mm_list = ['abs_solutionabs_03CZ01L_spectro_mod_Jacobs_Hochheiser','abs_solutionabs_01CZ01L_spectro_mod_Jacobs_Hochheiser','abs_solutionabs_30NO01L_TGS_ANSA','abs_solutionabs_d_0047NO01L_TGS_ANSA','abs_solutionabs_30NO01L_ANSA','abs_solutionabs_1GR01L_ANSA','abs_solutionabs_d_0043NO01L_TGS_ANSA','abs_solutionabs_01IE01L_ANSA','abs_solutionabs_d_0015NO01L_ANSA','abs_solutionabs_d_0008NO01L_TGS_ANSA','abs_solutionabs_d_0044NO01L_TGS_ANSA', \
                                'abs_solutionabs_d_0008NO01L_ANSA', 'abs_solutionabs_d_0039NO01L_ANSA','abs_solutionabs_d_0041NO01L_ANSA','abs_solutionabs_d_0001NO01L_ANSA','abs_solutionabs_d_0001NO01L_TGS_ANSA','abs_solutionabs_1GR01L_ANSA','abs_solutionabs_d_0008NO01L_ANSA','abs_solutionabs_d_0001NO01L_ANSA','abs_solutionabs_d_0044NO01L_TGS_ANSA', \
                                'abs_solutionabs_d_0039NO01L_ANSA','abs_solutionabs_01IE01L_ANSA','abs_solutionabs_d_0041NO01L_ANSA','abs_solutionabs_d_0015NO01L_ANSA','abs_solutionabs_d_0001NO01L_TGS_ANSA','abs_solutionabs_30NO01L_ANSA','abs_solutionabs_d_0043NO01L_TGS_ANSA','abs_solutionabs_d_0008NO01L_TGS_ANSA','abs_solutionabs_30NO01L_TGS_ANSA', \
                                'abs_solutionabs_d_0047NO01L_TGS_ANSA','glass_sintergs_16LV33L_indo','glass_sinterglass_2SE01L_NEDA','abs_tubeabst_2SE01L_NEDA', \
                                'abs_solutionabs_d_0916NO01L_NEDA','abs_solutionabs_d_0618NO01L_NEDA','abs_solutionabs_2PL01L_NEDA','glass_sinterglass_13SE01L_NEDA','abs_tubeabst_12SU01L_NEDA','abs_solutionabs_d_0762NO01L_NEDA', 'abs_tubeabst_10SU01L_NEDA','filter_abs_solutionfabs_2HU01L_Trieth-NEDA', \
                                'glass_sinterglass_11SE01L_NEDA','abs_tubeabst_3SU01L_NEDA','abs_solutionabs_d_0796NO01L_NEDA','glass_sinterglass_12SE01L_NEDA','abs_tubeabst_11SE01L_NEDA','abs_solutionabs_d_0797NO01L_NEDA','abs_tubeabst_8SE01L_NEDA','filter_1packf1p_02CZ01L_NEDA', 'filter_abs_solutionfabs_1HU01L_Trieth-NEDA','abs_tubeabst_1SU01L_NEDA', \
                                'abs_tubeabst_14SU01L_NEDA','abs_solutionabs_d_0047NO01L_NEDA', 'abs_tubeabst_5SE01L_NEDA','abs_solutionabs_d_0044NO01L_NEDA','abs_solutionabs_1PL01L_NEDA','abs_tubeabst_5SU01L_NEDA','glass_sinterglass_3SE01L_NEDA','abs_tubeabst_3SE01L_NEDA', 'abs_tubeabst_12SE01L_NEDA','glass_sinterglass_05SE01L_NEDA','glass_sinterglass_08SE01L_NEDA', \
                                'abs_tubeabst_9SU01L_NEDA','glass_sinterglass_05DK01L_NEDA','abs_tubeabst_03DK01L_NEDA','abs_tubeabst_11SU01L_NEDA','abs_tubeabst_4SU01L_NEDA','glass_sinterNILU_SS2000_sequential_air_samplerSI03L_spectrophotometry','glass_sinterglass_d_01IE01L_spectro','abs_solutionabs_sol_08ME01L_spectrophm','glass_sinterNILU_SS2000_sequential_air_samplerSI01L_spectrophotometry']
        
        color_mm_list = ['abs_solutionabs_02PL01L_abs_griess', 'abs_solutionabs_1TR01L_spectro_griess', 'abs_solutionabs_5SE01L_saltzman','abs_solutionabs_07SK01L_saltzman','abs_solutionabs_d_0042NO01L_Griess_air', 'glass_sinterglass_01SE01L_saltzman',  'abs_tubeabst_13RU01L_Griess', 'glass_sinterglass_d_0047NO01L_Griess_air', \
                           'abs_solutionimpingerSK01L_modified_saltzman', 'abs_tubeabst_16RU01L_Griess', 'glass_sinterglass_11SE01L_Griess', 'abs_tubeabst_6SU01L_Griess', 'abs_solutionabs_04SK01L_Griess', 'abs_solutionabs_02CZ01L_Griess', \
                           'glass_sinterglass_05SE01L_Griess', 'glass_sinterglass_d_0044NO01L_Griess_air', 'glass_sinterglass_d_0797NO01L_Griess_air', 'abs_tubeabst_1RU01L_Griess','abs_tubeabst_10LV33L_Griess', 'abs_solutionabs_d_0008NO01L_Griess_air', 'glass_sinterglass_01NO01L_Griess_air', 'abs_tubeabst_5SU01L_Griess', 'glass_sinterglass_d_0002NO01L_Griess_air', \
                           'abs_tubeabst_16LV33L_Griess', 'filter_1packf1p_5PL02L_Griess', 'filter_3packf3p_d_0039NO01L_Griess_air', 'glass_sinterglass_d_0055NO01L_Griess_air', 'glass_sinterglass_12SE01L_Griess', 'abs_solutionabs_01IE01L_Griess', 'glass_sinterglass_d_0039NO01L_Griess_air', \
                           'abs_solutionabs_3PL01L_Griess','filter_abs_solutionfabs_2HU01L_Griess', 'abs_tubeabst_14RU01L_Griess', 'abs_solutionabs_06SK01L_saltzman', 'abs_tubeabst_10LV01L_Griess', 'abs_solutionabs_05SK01L_Griess', 'filter_abs_solutionfabs_12SE01L_saltzman', 'abs_solutionabs_03PL01L_Griess', 'abs_solutionabs_5PL02L_Griess', \
                           'filter_1packf3p_5_1PL02L_f1p_Griess', 'glass_sinterglass_d_1010NO01L_Griess_air', 'glass_sinterglass_d_0796NO01L_Griess_air', 'abs_solutionabs_30NO01L_Griess_air','glass_sinterglass_13SE01L_Griess', 'glass_sinterglass_15LT01L_Griess', 'abs_solutionabs_04PL01L_abs_griess','abs_solutionabs_03PL01L_abs_griess', 'glass_sinterglass_d_0015NO01L_Griess_air', \
                           'abs_solutionabs_06SK01L_Griess', 'glass_sinterglass_d_0001NO01L_Griess_air', 'abs_solutionabs_05CZ01L_Griess','glass_sinterglass_30NO01L_Griess_air', 'abs_solutionabs_02SK01L_saltzman','abs_solutionabs_07SK01L_Griess', 'glass_sinterglass_01MD02L_Griess_air', 'abs_solutionabs_04CZ01L_Griess', 'glass_sinterglass_08SE01L_Griess', \
                           'abs_solutionabs_d_0015NO01L_Griess_air','glass_sinterglass_d_0008NO01L_Griess_air','glass_sinterglass_d_0043NO01L_Griess_air','abs_solutionabs_15LT01L_Griess_air', 'glass_sinterglass_d_0042NO01L_Griess_air','glass_sinterglass_d_0916NO01L_Griess_air','abs_tubeabs_12MD01L_Griess','abs_tubeabst_39NO01L_Griess_air', \
                           'glass_sinterglass_15NO01L_Griess_air', 'abs_solutionabs_2PL01L_abs_griess', 'abs_solutionabs_2PL01L_Griess','abs_solutionabs_d_0041NO01L_Griess_air','abs_solutionabs_4PL01L_Griess','filter_1packf1p_5PL02L_f1p_Griess','glass_sinterglass_d_0041NO01L_Griess_air','glass_sinterglass_d_0030NO01L_Griess_air', \
                           'glass_sinterglass_02SE01L_saltzman', 'abs_solutionabs_1PL01L_Griess','abs_tubeabst_16LV01L_Griess', 'glass_sinterglass_02SE01L_Griess', 'filter_1packf1p_05PL02L_f1p_Griess','abs_solutionabs_04SK01L_saltzman', 'abs_solutionabs_d_0001NO01L_Griess_air','glass_sinterglass_08NO01L_Griess_air','abs_solutionabs_d_0039NO01L_Griess_air', \
                           'glass_sinterglass_1TR01L_spectro_griess','abs_solutionabs_02SK01L_Griess','abs_solutionabs_05SK01L_saltzman','glass_sinterglass_d_0056NO01L_Griess_air', \
                           'continuous_colorimetriccon_col_04FI01L_Saltzman', 'filter_2packf2p_01CZ01L_Griess', 'continuous_colorimetriccon_col_09FI01L_Saltzman', 'filter_abs_solutionfabs_04DE04L_Saltzman','abs_solutionabs_01CZ01L_Griess', 'filter_abs_solutionfabs_09DE03L_Saltzmann', 'filter_abs_solutionfabs_05DE05L_Saltzman', \
                           'glass_sinterA_NaI_impregnated_02AM01L_uv_abs_Griess', 'filter_abs_solutionfabs_5ES01L_Griess', 'continuous_colorimetriccon_col_17FI01L_Saltzman', 'abs_solutionabs_08DE03L_saltzman', 'abs_solutionabs_02DE03L_saltzman', 'filter_abs_solutionfabs_02DE03L_Saltzmann','filter_abs_solutionfabs_7ES01L_Griess', \
                           'continuous_colorimetriccon_col_22FI01L_Saltzman','filter_abs_solutionfabs_07DE03L_Saltzmann','filter_abs_solutionfabs_04DE03L_Saltzman','glass_sinterglass_03DK01L_Griess','filter_abs_solutionfabs_05DE03L_Saltzman','abs_solutionabs_03CZ01L_Griess', 'abs_tubeabst_4BY01L_Griess', \
                           'filter_abs_solutionfabs_07DE07L_Saltzmann', 'abs_solutionabs_07DE03L_saltzman', 'glass_sinterglass_05DK01L_Griess', 'glass_sinterglass_08DK01L_Griess', 'abs_solutionabs_05DE03L_saltzman', 'abs_solutionabs_03DE03L_saltzman','abs_solutionabs_09DE03L_saltzman', 'filter_2packf2p_03CZ01L_Griess', \
                           'glass_sinterglass_9EE01L_Griess', 'abs_solutionabs_04DE03L_saltzman', 'filter_abs_solutionfabs_03DE03L_Saltzman', 'filter_abs_solutionfabs_08DE03L_Saltzmann', 'filter_abs_solutionfabs_01DE03L_Saltzmann','glass_sinterA_NaI_impregnated_02AM01L_spectrometric_Griess','abs_tubeabst_9EE01L_Griess']
        
        fia_mm_list = ['NaJ_solutionNaJ_07DE03L_FIA','NaJ_solutionNaJ_09DE03L_FIA','NaJ_solutionNaJ_05DE03L_FIA', 'NaJ_solutionNaJ_08DE03L_FIA','NaJ_solutionNaJ_03DE03L_NaJ_method','NaJ_solutionNaJ_08DE03L_NaJ_method','NaJ_solutionNaJ_09DE03L_NaJ_fia','NaJ_solutionNaJ_03DE03L_FIA','NaJ_solutionNaJ_09DE03L_NaJ_method','NaJ_solutionNaJ_01DE03L_FIA', \
                       'NaJ_solutionNaJ_04DE03L_FIA','NaJ_solutionNaJ_02DE03L_NaJ_method','NaJ_solutionNaJ_07DE03L_NaJ_method','NaJ_solutionNaJ_02DE03L_FIA','glass_sinterglass_08SE01L_FIA','glass_sinterglass_11SE01L_FIA','glass_sinterglass_02SE01L_FIA','glass_sinterglass_05SE01L_FIA']

        coul_mm_list = []
        diffsamp_mm_list = ['diffusion_tubediffusion_tubeGB02L_UK']
        
        notclear_mm_list = ['abs_solutionabs_3','abs_tubeabst_5','abs_tubeabst_6','abs_tubeabst_3', 'abs_tubeabst_1','abs_tubeabs_08SE01L_AD', 'abs_solutionabs_57NO01L_AW','abs_solutionabs_11', 'abs_solutionabs_12','abs_tubeabs_11_1SE01L_AD','abs_tubeabs_05_1SE01L_AD','abs_tubeimpregnated_filterSE01L_AD','abs_tubeabst_10','abs_tubeabst_12', \
                            'abs_solutionabs_7', 'abs_solutionabs_2','abs_solutionabs_1','abs_solutionabs_8','abs_tubeabs_14SE01L_AD','filter_1packf1p_02', 'filter_3packf3p_d_0002','filter_abs_solutionfabs_5','filter_abs_solutionfabs_8','filter_abs_solutionfabs_3', 'filter_abs_solutionfabs_2','filter_abs_solutionfabs_11','filter_abs_solutionfabs_12', \
                            'filter_3packf3p_2','NANALT01L_NA','glass_sintergs_10LV33L_glass_sinter','NANAIT01L_NA', 'NANANL01L_NA','abs_tubeabst_4','abs_tubeabst_11','abs_tubeabst_9','filter_2packf2p_03','filter_2packf2p_01','filter_abs_solutionfabs_4','filter_abs_solutionfabs_1','filter_abs_solutionfabs_6','filter_3packf3p_1CZ01L_PJ_T21-VS-009', \
                            'filter_3packf3p_08','abs_solutionabs_6', 'abs_solutionabs_4','abs_solutionabs_1CZ01L_PJ_T21-VS-011','NANADE03L_NA','NANAAT01L_CO','NANA']

#--------------------------------------------------------------------------------------------------- 
    elif process_group == 'WMO GAW':
        continuous_st_list = ['continuous']
        remote_st_list = ['continuous(carbondioxide),remotespectroscopicmethod(methaneandsurfaceozone)']
        flask_st_list = ['flask']
        filter_st_list = ['filter']
        
        uv_mm_list = ['ultraviolet photometry','Lightabsorptionanalysis(UV)','Ultravioletabsorptionmethod']
        vuf_mm_list = ['VUVFlourescence','VUVFluorescence','vacuum ultraviolet fluorescence']
        cld_mol_mm_list = ['chemiluminescence','chemiluminescence (conversion-molybdenum)','Chemiluminescence']
        cld_photo_mm_list = ['chemiluminescence (conversion-photolysis)']
        cld_eth_mm_list = ['ethylene chemiluminescence']
        poti_mm_list = ['potassium iodide']
        ecc_mm_list = ['electrochemical concentration cell']
        ndir_mm_list = ['NDIR','Lightabsorptionanalysis(IRexceptNDIR)','1.)NDIR2.)resonancefluorescence','non-dispersive infrared spectroscopy']
        doas_mm_list = []
        crs_mm_list = ['CavityRingdownSpectroscopy','Horibainstruments:non-dispersiveinfrareddetection(NDIR)Picarroinstruments:CavityRingDownSpectroscopyCRDS','Cavityring-downspectroscopy(CRDS)','CavityRing-DownSpectroscopy(CRDS)','Horibainstruments:non-dispersiveinfrareddetection(NDIR);Picarroinstruments:CavityRingDownSpectroscopyCRDS','Horibainstruments:non-dispersiveinfrareddetection(NDIR)Picarroinstrument:CavityRingDownSpectroscopyCRDS']
        capss_mm_list = []
        gcfid_mm_list = ['GasChromatography(FID)','GasChromatography(other)']
        gcrgd_mm_list = ['Gas Chromatography (RGD)','GC-HgO','GasChromatography(RGD)','gas chromatography reduction gas detection']
        gcms_mm_list = []
        gcpid_mm_list = []
        gcftir_mm_list = []
        fid_mm_list = []
        flamephoto_mm_list = []
        ptrms_mm_list = []
        ionchrom_mm_list = ['ion chromatography']
        spectrometry_mm_list = []
        spectrophoto_mm_list = ['spectrophotometry']
        color_mm_list = ['colorimetric']
        fia_mm_list = []
        coul_mm_list = []
        diffsamp_mm_list = ['diffusive sampler']
        notclear_mm_list = []

#---------------------------------------------------------------------------------------------------  
    elif (process_group == 'CAPMON') or (process_group == 'CASTNET') or (process_group == 'EANET') or (process_group == 'NAPS') or (process_group == 'SEARCH'):
        continuous_st_list = ['continuous']
        remote_st_list = []
        flask_st_list = ['flask']
        filter_st_list = []
        
        uv_mm_list = ['ultraviolet photometry']
        vuf_mm_list = []
        cld_mol_mm_list = ['chemiluminescence','chemiluminescence (conversion-molybdenum)']
        cld_photo_mm_list = ['chemiluminescence','chemiluminescence (conversion-photolysis)']
        cld_eth_mm_list = []
        poti_mm_list = []
        ecc_mm_list = []
        ndir_mm_list = ['non-dispersive infrared spectrometry','non-dispersive infrared spectroscopy']
        doas_mm_list = ['differential optical absorption spectrosocopy']
        crs_mm_list = []
        capss_mm_list = []
        gcfid_mm_list = ['gas chromatography flame ionisation detection']
        gcrgd_mm_list = []
        gcms_mm_list = ['gas chromatography mass selective detection']
        gcpid_mm_list = []
        gcftir_mm_list = []
        fid_mm_list = []
        flamephoto_mm_list = []
        ptrms_mm_list = []
        ionchrom_mm_list = []
        spectrometry_mm_list = []
        spectrophoto_mm_list = []
        color_mm_list = []
        fia_mm_list = []
        coul_mm_list = []
        diffsamp_mm_list = []
        notclear_mm_list = []

#---------------------------------------------------------------------------------------------------  
    #process standardised sampling methods into grid
    #GAW has raw sampling types
    
    st_test = np.array([x in continuous_st_list for x in raw_st_grid])
    if True in st_test:
        #p_st_grid[st_test] = 'continuous'
        p_st_grid[st_test] = 1
    
    st_test = np.array([x in remote_st_list for x in raw_st_grid])
    if True in st_test:
        #p_st_grid[st_test] = 'remote'
        p_st_grid[st_test] = 2

    st_test = np.array([x in flask_st_list for x in raw_st_grid])
    if True in st_test:
        #p_st_grid[st_test] = 'flask'
        p_st_grid[st_test] = 3

    st_test = np.array([x in filter_st_list for x in raw_st_grid])
    if True in st_test:
        #p_st_grid[st_test] = 'filter'
        p_st_grid[st_test] = 4
    
    #process standardised measurement methods into grid

    mm_test = np.array([x in uv_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'ultraviolet photometry'
        p_mm_grid[mm_test] = 1
        
    mm_test = np.array([x in vuf_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'vacuum ultraviolet fluorescence'
        p_mm_grid[mm_test] = 2
   
    mm_test = np.array([x in cld_mol_mm_list for x in raw_mm_grid])
    if True in mm_test:
        if (species == 'NO') or (species == 'O3'):
            #p_mm_grid[mm_test] = 'chemiluminescence'
            p_mm_grid[mm_test] = 3
        if species == 'NO2':
            #p_mm_grid[mm_test] = 'chemiluminescence (conversion-molybdenum)'
            p_mm_grid[mm_test] = 4
        
    mm_test = np.array([x in cld_photo_mm_list for x in raw_mm_grid])
    if True in mm_test:
        if (species == 'NO') or (species == 'O3'):
            #p_mm_grid[mm_test] = 'chemiluminescence'
            p_mm_grid[mm_test] = 3
        if species == 'NO2':
            #p_mm_grid[mm_test] = 'chemiluminescence (conversion-photolysis)'
            p_mm_grid[mm_test] = 5
     
    mm_test = np.array([x in cld_eth_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'ethylene chemiluminescence'
        p_mm_grid[mm_test] = 6

    mm_test = np.array([x in poti_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'potassium iodide'
        p_mm_grid[mm_test] = 7

    mm_test = np.array([x in ecc_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'electrochemical concentration cell'
        p_mm_grid[mm_test] = 8
        
    mm_test = np.array([x in ndir_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'non-dispersive infrared spectroscopy'
        p_mm_grid[mm_test] = 9
        
    mm_test = np.array([x in doas_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'differential optical absorption spectrosocopy'
        p_mm_grid[mm_test] = 10
        
    mm_test = np.array([x in crs_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'cavity ringdown spectroscopy'
        p_mm_grid[mm_test] = 11
    
    mm_test = np.array([x in capss_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'cavity attenuated phase shift spectroscopy'
        p_mm_grid[mm_test] = 12
        
    mm_test = np.array([x in gcfid_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'gas chromatography flame ionisation detection'
        p_mm_grid[mm_test] = 13

    mm_test = np.array([x in gcrgd_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'gas chromatography reduction gas detection'
        p_mm_grid[mm_test] = 14
        
    mm_test = np.array([x in gcms_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'gas chromatography mass spectrometry'
        p_mm_grid[mm_test] = 15
        
    mm_test = np.array([x in gcpid_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'gas chromatography photo ionization detection'
        p_mm_grid[mm_test] = 17
        
    mm_test = np.array([x in gcftir_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'gas chromatography fourier transform infrared spectroscopy'
        p_mm_grid[mm_test] = 18
        
    mm_test = np.array([x in fid_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'flame ionisation detection'
        p_mm_grid[mm_test] = 22
        
    mm_test = np.array([x in flamephoto_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'photometric flame photometry'
        p_mm_grid[mm_test] = 23
        
    mm_test = np.array([x in ptrms_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'proton transfer reaction mass spectrometry'
        p_mm_grid[mm_test] = 24
    
    mm_test = np.array([x in ionchrom_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'ion chromatography'
        p_mm_grid[mm_test] = 25
        
    mm_test = np.array([x in spectrometry_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'spectrometry'
        p_mm_grid[mm_test] = 26
        
    mm_test = np.array([x in spectrophoto_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'spectrophotometry'
        p_mm_grid[mm_test] = 27
        
    mm_test = np.array([x in color_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'colorimetric'
        p_mm_grid[mm_test] = 28
        
    mm_test = np.array([x in fia_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'fluid injection analysis'
        p_mm_grid[mm_test] = 29
        
    mm_test = np.array([x in coul_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'coulometry'
        p_mm_grid[mm_test] = 30
    
    mm_test = np.array([x in diffsamp_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'diffusive sampler'
        p_mm_grid[mm_test] = 31
        
    mm_test = np.array([x in notclear_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'nonsensical instrument'
        p_mm_grid[mm_test] = 99
    
    #test for blanks
    blank_mm_list = ['','UNKNOWNUNKNOWN',np.NaN]
    mm_test = np.array([x in blank_mm_list for x in raw_mm_grid])
    if True in mm_test:
        #p_mm_grid[mm_test] = 'blank'
         p_mm_grid[mm_test] = 0
    
    #if there measurement types not in valid measurement list then change then print them to screen and append to lists
    unknown_mm_local = []
    for x in raw_mm_grid:
        if (x not in uv_mm_list) & (x not in vuf_mm_list)  & (x not in cld_mol_mm_list) & (x not in cld_photo_mm_list) & (x not in cld_eth_mm_list) & (x not in poti_mm_list) & (x not in ecc_mm_list) & (x not in ndir_mm_list) & (x not in doas_mm_list) & (x not in crs_mm_list) & (x not in capss_mm_list) & (x not in gcfid_mm_list) & (x not in gcrgd_mm_list) & (x not in gcms_mm_list) & (x not in gcpid_mm_list) & (x not in gcftir_mm_list) & (x not in fid_mm_list) & (x not in flamephoto_mm_list) & (x not in ptrms_mm_list) & (x not in ionchrom_mm_list) & (x not in spectrometry_mm_list) & (x not in spectrophoto_mm_list) & (x not in color_mm_list) & (x not in fia_mm_list) & (x not in coul_mm_list) & (x not in diffsamp_mm_list) & (x not in notclear_mm_list) & (x not in blank_mm_list) & (x != 'na'):
            unknown_mm_local.append(x)
    if len(unknown_mm_local) > 0:
        print 'Unknown mm for ref:', set(unknown_mm_local)
        unknown_mm = np.append(unknown_mm,unknown_mm_local)
        unknown_mm_refs.append(ref)
    
    #FOR PROCESSED GRIDS, WHERE VALUES ARE -99999 THEN make invalid
    #raw grids are before flagsandlod
    #p grids are after duplicate point removal
    invalid_boolean_duplicate = full_data < 0
    invalid_boolean_flagsandlod = full_data_after_flagsandlod < 0
    raw_st_grid = np.copy(p_st_grid)
    raw_mm_grid = np.copy(p_mm_grid)
    p_st_grid_after_flagsandlod = np.copy(p_st_grid)
    p_mm_grid_after_flagsandlod = np.copy(p_mm_grid)
    #p_st_grid[invalid_boolean] = 'na'
    #p_mm_grid[invalid_boolean] = 'na' 
    p_st_grid[invalid_boolean_duplicate] = -1
    p_mm_grid[invalid_boolean_duplicate] = -1 
    p_st_grid_after_flagsandlod[invalid_boolean_flagsandlod] = -1
    p_mm_grid_after_flagsandlod[invalid_boolean_flagsandlod] = -1

    return raw_st_grid,p_st_grid,p_st_grid_after_flagsandlod,raw_mm_grid,p_mm_grid,p_mm_grid_after_flagsandlod,unknown_mm,unknown_mm_refs


#DO FIRST QUALITY CONTROL PROCESS CHECKS ON SURFACE OBSEVRATIONAL DATA - IN PROCESS GROUP PROCESS FILES
def primary_quality_control(site_ref,species,file_res,no2_type,grid_dates,full_data,n_dup_array,valid_hours_dup,raw_st_grid,p_st_grid,p_st_grid_afterflagsandlod,raw_mm_grid,p_mm_grid,p_mm_grid_afterflagsandlod,data_resolution,n_obs_valid,key_meta,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod):
    
    data_valid = True
    
    p_mm_grid_afterduplicate = np.copy(p_mm_grid)
    
#-----------------------------------------------------------
    #CHECK THERE IS AT LEAST 1 VALID DATAPOINT IN DATASET
    valid_data = full_data[np.where(full_data>=0)]
    if len(valid_data) == 0:
        data_valid = False
        inv_anyvaliddata +=1
        print 'Site Invalid, No valid data in dataset'
        #for molybdenum no2 process only save valid data counts
        if no2_type == 'MOLYBDENUM':
            n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod = 0,0,0,0,0,0,0,0,0,0,0,0,0
        return data_valid,full_data,valid_hours_dup,p_st_grid,p_mm_grid,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod,'anyvaliddata'
    n_obs_after_anyvaliddata += n_obs_valid
    
#-----------------------------------------------------------
#CHECK KEY METADATA EXITS FOR EACH SITE IS FLOAT (OR EVEN EXISTS) (LAT,LON,ALT)
#ALSO CHECK LAT AND LON ARE NOT SET AS 0 (MUST BE BUG AS NOTHING THERE BUT OCEAN)
    lat = key_meta[0]
    lon = key_meta[1]
    alt = key_meta[2]
    try:
        lat % 1
        lon % 1
        alt % 1
        
        if (lat == 0.0) & (lon == 0.0):
            1+'a'
    except:
        data_valid = False
        inv_nokeymeta +=1
        print ' Site Invalid,Key metadata is missing for site'
        #for molybdenum no2 process only save valid data counts
        if no2_type == 'MOLYBDENUM':
            n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod = 0,0,0,0,0,0,0,0,0,0,0,0,0
        return data_valid,full_data,valid_hours_dup,p_st_grid,p_mm_grid,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod,'nokeymeta'
        
    n_obs_after_nokeymeta += n_obs_valid
    
#-----------------------------------------------------------
#check data resolution by each year
#do by checking checking differences between valid points either side
#make sure data is more more than 1 day
#exclude 0 differences, if minimum difference is > data_resolution min exclude site as data resolution is too poor
#also test to make sure vals in year are not all the same value
#also test to remove data resolution per year
    if species.lower() == 'o3':
        minres = 4.9
        minres_refs = ['abes1399a','aqs150031004','aqs180970031','aqs180970042','aqs181270024','aqs280490010','aqs280890001','aqs371190034','aqs371191005','aqs371191009','aqs470650028'] 
    elif species.lower() == 'no':
        minres = 1.9
        minres_refs = ['abat31496','abdebb053','abdebb065','abdebb066','abdehe023','abdehe024','abdehe025','abdehe026','abdehe027','abdehe028','abdehe034','abdehe042','abdehe043','abdenw001','abdenw063','abdenw064','abdenw065','abdenw066','abdenw068','abdenw081','abdeth007','abdeth026','abdeth027','abdeth029','abdeth030','abdeth037','abdeth040','abdeth042','abdeth061','abes1654a','abit0988a','abit1174a','abit1519a','181470008','320030023','540250001','550210008','550430003','550591001','551110007',
                       'abat60145','abat80706','abat80807','abat82707','abat82801','abat900ke','abat900za','abat90tab','abat9stef','abbetb006','abbetb011','abbetn043','abbetr001','abbetr012','abbetr201','abbetr501','abbetr701','abbetr731','abbetr740','abbetwol1','abdebb007','abdebb021','abdebb029','abdebb032','abdebb042','abdebb044','abdebb048','abdebb049','abdebb050','abdebb064','abdebb067','abdebb068','abdebb075','abdebb081',
                       'abdebb082','abdebb083','abdebb084','abdehe001','abdehe006','abdehe008','abdehe010','abdehe013','abdehe014','abdehe017','abdehe018','abdehe020','abdehe021','abdehe022','abdehe029','abdehe030','abdehe031','abdehe032','abdehe033','abdehe040','abdehe041','abdenw002','abdenw004','abdenw006','abdenw008','abdenw010','abdenw011','abdenw012','abdenw013','abdenw015','abdenw018','abdenw021','abdenw022','abdenw023','abdenw024','abdenw028','abdenw029','abdenw030','abdenw034','abdenw036',
                       'abdenw037','abdenw038','abdenw039','abdenw040','abdenw042','abdenw043','abdenw047','abdenw048','abdenw050','abdenw051','abdenw053','abdenw058','abdenw059','abdenw060','abdenw062','abdenw067','abdenw071','abdenw072','abdenw077','abdenw078','abdenw079','abdenw080','abdenw082','abdenw094','abdenw095','abdenw096','abdenw097','abdenw098','abdenw099','abdenw100','abdenw101',
                       'abdenw102','abdenw112','abdenw116','abdenw133','abdenw134','abdenw135','abdenw136','abdenw137','abdenw138','abdenw141','abdenw143','abdenw145','abdenw177','abdenw179','abdenw180','abdenw181','abdenw183','abdenw199','abderp001','abderp002','abderp018','abderp019','abderp020','abdeth002','abdeth005','abdeth009','abdeth011','abdeth013','abdeth016','abdeth018','abdeth020','abdeth021','abdeth024','abdeth025','abdeth031','abdeth032','abdeth036','abdeth039','abdeth041','abdeth043','abdeth060','abdeth063','abdeth064','abdeth065','abdeth066','abdeth067',
                       'abdeth072','abdeth073','abdeth075','abdeth080','abdeth081','abdeth084','abes0327a','abes0330a','abes0629a','abes0633a','abes0634a','abes0651a','abes0774a','abes0817a','abes0825a','abes1067a','abes1074a','abes1130a','abes1165a','abes1169a','abes1183a','abes1268a','abes1271a','abes1351a','abes1394a','abes1424a','abes1425a','abes1492a','abes1519a','abes1560a','abes1589a','abes1593a','abes1595a','abes1597a','abes1620a','abes1640a','abes1653a','abes1845a','abgr0002a','abgr0003a','abgr0027a','abgr0028a','abgr0029a','abgr0030a','abgr0031a','abgr0035a','abhu0008a','abit0559a','abit0754a','abit0906a','abit0966a','abit1010a','abit1035a',
                       'abit1524a','abit1692a','abmk0036a','abpt03010','aqs060830011','aqs060850005','aqs090031003','aqs110010043','aqs130890002','aqs180830004','aqs180891016','aqs180970030','aqs180970070','aqs180970073','aqs181571001','aqs181630012','aqs250250040','aqs280490010','aqs340030005','aqs370670022','aqs371190041','aqs440071007','aqs470370011','aqs500070003','aqs500070013','aqs500070014','aqs500210002','aqs540110006','aqs540390004','aqs540690007','aqs550590016','aqs550790080','aqs550890005','hhe','canaps100125']
    elif species.lower() == 'no2':
        if no2_type == 'MOLYBDENUM':
            minres = 1.9
            minres_refs = ['abat31502','abbetn132','abdebb065','abdenw001','abdenw063','abdenw064','abdenw065','abdenw066','abdenw068','abdenw081','abdeth007','abdeth029','abdeth030','abdeth037','abdeth040','abdeth042','abdeth061','abes1436a','abes1441a','abes1542a','abes1543a','abes1648a','abes1654a','abes1670a','abes1671a','abes1689a','abes1690a','abfr28011','abit1548a','132230003','180330002','180770001','180830004','181470002','181470006','540250001','550210008','550430003','550591001','551051002','551051004','551110007','be0032r','065401',
                           'abat80706','abat80807','abat82707','abat82801','abba0029a','abbetb011','abbetmeu1','abbetn066','abbetr001','abbetr012','abbetr512','abdenw002','abdenw004','abdenw006','abdenw008','abdenw012','abdenw013','abdenw015','abdenw018','abdenw021','abdenw022','abdenw023','abdenw029','abdenw030','abdenw036','abdenw038','abdenw039','abdenw040','abdenw042','abdenw047','abdenw050','abdenw051','abdenw058','abdenw060','abdenw062','abdenw067','abdenw072','abdenw077','abdenw078','abdenw079','abdenw080','abdenw094','abdenw095','abdenw096','abdenw098','abdenw100','abdenw116','abdenw141','abdenw143','abdenw177','abdenw179','abdenw180','abdenw181',
                           'abderp002','abdeth002','abdeth005','abdeth009','abdeth011','abdeth013','abdeth016','abdeth017','abdeth018','abdeth020','abdeth021','abdeth024','abdeth025','abdeth031','abdeth032','abdeth036','abdeth041','abdeth060','abdeth062','abdeth063','abdeth065','abdeth067','abdeth072','abdeth073','abdeth080','abdeth084','abes0584a','abes0691a','abes0692a','abes0693a','abes0924a','abes0975a','abes1119a','abes1123a','abes1125a','abes1148a','abes1181a','abes1182a','abes1184a','abes1185a','abes1195a','abes1275a','abes1287a','abes1312a','abes1380a','abes1386a','abes1387a','abes1405a','abes1428a','abes1435a','abes1437a','abes1445a','abes1453a','abes1519a','abes1617a',
                           'abes1619a','abes1620a','abes1623a','abes1624a','abes1625a','abes1675a','abes1687a','abes1688a','abes1691a','abes1751a','abes1765a','abes1826a','abes1834a','abes1849a','abes1877a','abes1884a','abes1885a','abes1886a','abes1911a','abes1912a','abes1914a','abes1915a','abes1926a','abfr28019','aqs130890002','aqs180030006','aqs180970073','aqs181530001','aqs181571001','aqs181630012','aqs212230004','aqs280490010','aqs340030005','aqs340190001','aqs340210005','aqs340230011','aqs370670007','aqs371190041','aqs500070013','aqs500210002','aqs550590016','aqs550790026']
        elif no2_type == 'PHOTOLYTIC':
            minres = 99999
            minres_refs = []
    elif species.lower() == 'co':
        minres_refs1 = ['abat30901','abat30903','abat31401','abat31402','abat31403','abat32301','abat32701','abbetn043','abbetwol1','abdeni048','abes1269a','abes1279a','abes1394a','abes1536a','abes1560a','abes1599a','abes1618a','abes1679a','abes1818a','abes1823a','ablt00003','371190041','abgb0754a','gb0050r','aqs090010010','aqs202090021','aqs250250042','aqs371830014','aqs410510080','aqs060670006','aqs060831008','aqs090010010','aqs090050004','aqs100032004','aqs110010043','aqs202090021','aqs320310016']
        minres_refs2 = ['abbetr201','abbetr501','abbetr512','abes1536a']
        minres_refs3 = ['abdest091','abes1016a','abes1172a','abes1178a','abes1183a','abes1185a','abes1239a','abes1240a','abes1400a','abes1428a','abes1445a','abes1464a','abes1609a','abes1617a','abes1619a','abes1623a','abes1625a','abes1635a','abes1675a','abes1685a','abes1691a','abes1818a','abes1826a','abes1834a','abes1849a','abes1867a','abes1877a','abes1884a','abes1886a','abes1905a','abes1911a','abes1912a','abes1915a','abes1926a','abie0135a','abit1456a']
        if site_ref.lower() in minres_refs1:
            minres = 149
            minres_refs = minres_refs1
        elif site_ref.lower() in minres_refs2:
            minres = 99
            minres_refs = minres_refs2
        elif site_ref.lower() in minres_refs3:
            minres = 49
            minres_refs = minres_refs3
        else:
            minres = 99999
            minres_refs = []
    elif species.lower() == 'isop':
        minres = 99999
        minres_refs = []

    start_year = int(grid_dates[0][:4])
    end_year = int(grid_dates[-1][:4])+1
    year_range = range(start_year,end_year)
    all_years = [i[:4] for i in grid_dates]
    all_years = np.array(all_years)
    all_months = [i[4:6] for i in grid_dates]
    all_months = np.array(all_months)
    for i in range(len(year_range)):
        year_test = all_years == str(year_range[i])
        year_vals = full_data[year_test]
        test = year_vals != -99999
        year_vals = year_vals[test]
        if len(year_vals) >= 24:
            unique_vals = np.unique(year_vals)
            min_diff = np.abs(np.diff(year_vals))
            test = min_diff != 0
            min_diff = min_diff[test]
            if len(min_diff) > 0:
                min_diff = np.min(min_diff)    
                #if min difference > min data resolution then year invalid
                if min_diff > data_resolution:
                    full_data[year_test] = -99999
                    #p_mm_grid[year_test] = 'na'
                    #p_st_grid[year_test] = 'na'
                    p_mm_grid[year_test] = -1
                    p_st_grid[year_test] = -1
            #if only 1 value in year, year is invalid
            if len(unique_vals) == 1:
                full_data[year_test] = -99999 
                p_mm_grid[year_test] = -1
                p_st_grid[year_test] = -1
            #test minimum data resolution per year - if ref in refs manually screened then remove year if min resolution above limit
            if site_ref.lower() in minres_refs:
                if np.min(year_vals) > minres:
                    full_data[year_test] = -99999 
                    p_mm_grid[year_test] = -1
                    p_st_grid[year_test] = -1
                    
            #add additional test for NO - due to high average bias assocciated with poor LOD
            #split into months per year
            #if 6 months or more have min resolution of >= 0.1 ppbv, then remove year
            if (species.lower() == 'no'): 
                inv_count = 0
                cut_months = all_months[year_test]
                for m in ['01','02','03','04','05','06','07','08','09','10','11','12']:
                    month_test = cut_months == m
                    var_month = full_data[year_test][month_test] 
                    var_month = var_month[var_month>=0]
                    try:
                        if np.min(var_month) >= 0.1:
                            inv_count+=1 
                    except:
                        inv_count+=1 
            		
                if inv_count >= 6:
                    full_data[year_test] = -99999 
                    p_mm_grid[year_test] = -1
                    p_st_grid[year_test] = -1
                    
    #CHECK THERE IS AT LEAST 1 VALID DATAPOINT IN DATASET AFTER DATA RESOLUTION REMOVAL
    valid_data = full_data[np.where(full_data>=0)]
    if len(valid_data) == 0:
        data_valid = False
        inv_resolution+=1
        print 'Site Invalid, Data resolution is not sufficient.'
        #for molybdenum no2 process only save valid data counts
        if no2_type == 'MOLYBDENUM':
            n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod = 0,0,0,0,0,0,0,0,0,0,0,0,0
        return data_valid,full_data,valid_hours_dup,p_st_grid,p_mm_grid,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod,'resolution'
        
    test = full_data >= 0
    valid_hours_dup = np.sum(n_dup_array[test])
    updated_n_obs_valid1 = int(len(full_data[test]) - valid_hours_dup)    
    n_obs_after_resolution += updated_n_obs_valid1
    
    p_mm_grid_afterresolution = np.copy(p_mm_grid)
    
#-----------------------------------------------------------
#INSTRUMENT CHECK - INTEGRATING METHODS REMOVED
#FOR EACH SPECIES HAVE INVALID MEASURMENT TYPES
    if species == 'O3':
        #bad_mm = ['ethylene chemiluminescence','potassium iodide','electrochemical concentration cell','gas chromatography flame ionisation detection','photometric flame photometry','chemiluminescence','spectrophotometry','nonsensical instrument','blank'] 
        bad_mm = [3,6,7,8,13,23,27,99,0]

    elif species == 'NO':
        #bad_mm = ['gas chromatography flame ionisation detection','ultraviolet photometry','nonsensical instrument','blank'] 
        bad_mm = [1,9,13,28,99,0]

    elif species == 'NO2':
        if no2_type == 'MOLYBDENUM':
            #bad_mm = ['chemiluminescence (conversion-photolysis)','ultraviolet photometry','ethylene chemiluminescence','potassium iodide','electrochemical concentration cell','vacuum ultraviolet fluorescence','non-dispersive infrared spectroscopy','differential optical absorption spectrosocopy','cavity ringdown spectroscopy','cavity attenuated phase shift spectroscopy','gas chromatography flame ionisation detection','gas chromatography reduction gas detection','gas chromatography mass spectrometry','gas chromatography mass selective detection','gas chromatography photo ionization detection','gas chromatography fourier transform infrared spectroscopy/mass spectrometry','gas chromatography flame ionisation detection/mass selective detection','gas chromatography mass spectrometry/flame ionisation detection','gas chromatography','flame ionisation detection','photometric flame photometry','proton transfer reaction mass spectrometry','ion chromatography','spectrometry','spectrophotometry','colorimetric','fluid injection analysis','coulometry','diffusive sampler','nonsensical instrument','blank']
            bad_mm = [1,2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,99,0]
        elif no2_type == 'PHOTOLYTIC':
            #bad_mm = ['chemiluminescence (conversion-molybdenum)','nonsensical instrument','gas chromatography flame ionisation detection','ultraviolet photometry','photometric flame photometry','blank']
            bad_mm = [1,4,13,23,27,28,29,31,99,0]
    
    elif species == 'CO':
        #bad_mm = ['ultraviolet photometry','nonsensical instrument','blank']
        bad_mm = [1,10,23,25,30,99,0]
    
    elif species == 'ISOP':
        #bad_mm = ['nonsensical instrument','blank']
        bad_mm = [99,0]

    instru_samp_inv_count = 0
    instru_test = np.array([x in bad_mm for x in p_mm_grid])
    if True in instru_test:
        #p_mm_grid[instru_test] = 'na'
        p_mm_grid[instru_test] = -1

#-----------------------------------------------------------
#SAMPLING CHECK - REMOTE SAMPLING REMOVED
    #st_test1 = np.array([x == 'remote' for x in p_st_grid])
    st_test1 = np.array([x == 2 for x in p_st_grid])
    if True in st_test1:
        #p_st_grid[st_test1] = 'na'
        p_st_grid[st_test1] = -1
#-----------------------------------------------------------
#SAMPLING CHECK - FILTER SAMPLING REMOVED
    #st_test2 = np.array([x == 'filter' for x in p_st_grid])
    st_test2 = np.array([x == 4 for x in p_st_grid])
    if True in st_test2:
        #p_st_grid[st_test2] = 'na'
        p_st_grid[st_test2] = -1
    
    combined_test = instru_test+st_test1+st_test2
    instru_samp_inv_count += len(full_data[np.where(combined_test == True)])-np.sum(n_dup_array[np.where(combined_test == True)])
    #valid_samp = p_st_grid[np.where(p_st_grid != 'na')]
    #valid_instru = p_mm_grid[np.where(p_mm_grid != 'na')]
    valid_samp = p_st_grid[np.where(p_st_grid != -1)]
    valid_instru = p_mm_grid[np.where(p_mm_grid != -1)]
    
    full_data[combined_test] = -99999
    #p_st_grid[combined_test] = 'na'
    #p_mm_grid[combined_test] = 'na'
    p_st_grid[combined_test] = -1
    p_mm_grid[combined_test] = -1
    
    if (len(valid_samp) == 0) or (len(valid_instru) == 0):
        data_valid = False
        inv_badmeasurementmethod+=1
        exit_code = 'badmeasurementmethod'
        print 'Site Invalid, Bad Measurement Method or sampling type.'
        #for molybdenum no2 process only save valid data counts 
        if no2_type == 'MOLYBDENUM':
            n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod = 0,0,0,0,0,0,0,0,0,0,0,0,0
        #for photolytic no2 process if exiting due to molybdenum instrument solely then do not add counts
        elif no2_type == 'PHOTOLYTIC':
            if (len(valid_instru) == 0) & (len(valid_samp) > 0) & (4 in raw_mm_grid):
                print 'Altering n counts as run is NO2 PHOTOLYTIC AND CANNOT HAVE COUNTS ASSOCIATED WITH MOLYBDENUM INSTRUMENTATION- AND RUN EXIT IS DUE TO MOLYBDENUM INSTRUMENTATION'
                n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod = 0,0,0,0,0,0,0,0,0,0,0,0,0
                exit_code = 'molybdenum'
    
        return data_valid,full_data,valid_hours_dup,p_st_grid,p_mm_grid,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod,exit_code
    
    #get updated n obs valid after removal of bad sampling/instrument points, 
    #being careful to remove n obs that are duplicated because of D or M resolution
    updated_n_obs_valid2 = updated_n_obs_valid1 - instru_samp_inv_count
    n_obs_after_badmeasurementmethod += updated_n_obs_valid2
    
    #need to be careful as molybdenum no2 process should only contain counts for molybdenum
    #so if there is a mixture of molybdenum and something else in process need to backdate counts to only those for valid molybdenum
    #firstly backdate n counts equal to n_obs_after_badmeasurementmethod
    #secondly add counts removed due to flagsandlod only associated with molybdenum

    if (no2_type == 'MOLYBDENUM') & (len(list(set(raw_mm_grid))) > 2):
        print 'Altering n counts as run is NO2 MOLYBDENUM AND ONLY WANT COUNTS ASSOCIATED WITH MOLYBDENUM INSTRUMENTATION- AND THERE IS DATA REMOVED THAT IS NON-MOLYBDENUM'
        test1 = raw_mm_grid == 4
        #test2 = p_mm_grid_afterflagsandlod == 'chemiluminescence (conversion-molybdenum)'
        #test3 = p_mm_grid_afterresolution == 'chemiluminescence (conversion-molybdenum)'
        test2 = p_mm_grid_afterflagsandlod == 4
        test3 = p_mm_grid_afterduplicate == 4
        test4 = p_mm_grid_afterresolution == 4
        
        n_dup_orig_sum = np.sum(n_dup_array[test1])
        n_dup_afterflagsandlod_sum = np.sum(n_dup_array[test2])
        n_dup_afterduplicate_sum = np.sum(n_dup_array[test3])
        n_dup_afterresolution_sum = np.sum(n_dup_array[test4])
        
        nmol1 = len(raw_mm_grid[test1]) - n_dup_orig_sum
        nmol2 = len(p_mm_grid_afterflagsandlod[test2]) - n_dup_afterflagsandlod_sum
        nmol3 = len(p_mm_grid_afterduplicate[test3]) - n_dup_afterduplicate_sum
        nmol4 = len(p_mm_grid_afterresolution[test4]) - n_dup_afterresolution_sum

        print n_all,nmol1,nmol2,nmol3,nmol4
        
        n_all,n_after_nometa = nmol1,nmol1
        n_after_flagsandlod = nmol2
        n_after_duplicate,n_obs_after_anyvaliddata,n_obs_after_nokeymeta = nmol3,nmol3,nmol3
        n_obs_after_resolution = nmol4
    
    elif (no2_type == 'PHOTOLYTIC') & (4 in raw_mm_grid):
        print 'Altering n counts as run is NO2 PHOTOLYTIC AND DO NOT WANT COUNTS ASSOCIATED WITH MOLYBDENUM INSTRUMENTATION- AND THERE IS DATA REMOVED THAT IS MOLYBDENUM'
        test1 = (raw_mm_grid != 4) & (raw_mm_grid != -1)
        #test2 = (p_mm_grid_afterflagsandlod != 'chemiluminescence (conversion-molybdenum)') & (p_mm_grid_afterflagsandlod != 'na')
        #test3 = (p_mm_grid_afterresolution != 'chemiluminescence (conversion-molybdenum)') & (p_mm_grid_afterresolution != 'na')
        test2 = (p_mm_grid_afterflagsandlod != 4) & (p_mm_grid_afterflagsandlod != -1)
        test3 = (p_mm_grid_afterduplicate != 4) & (p_mm_grid_afterduplicate != -1)
        test4 = (p_mm_grid_afterresolution != 4) & (p_mm_grid_afterresolution != -1)
        
        n_dup_orig_sum = np.sum(n_dup_array[test1])
        n_dup_afterflagsandlod_sum = np.sum(n_dup_array[test2])
        n_dup_afterduplicate_sum = np.sum(n_dup_array[test3])
        n_dup_afterresolution_sum = np.sum(n_dup_array[test4])
        
        nmol1 = len(raw_mm_grid[test1]) - n_dup_orig_sum
        nmol2 = len(p_mm_grid_afterflagsandlod[test2]) - n_dup_afterflagsandlod_sum
        nmol3 = len(p_mm_grid_afterduplicate[test3]) - n_dup_afterduplicate_sum
        nmol4 = len(p_mm_grid_afterresolution[test4]) - n_dup_afterresolution_sum

        print n_all,nmol1,nmol2,nmol3,nmol4
        
        n_all,n_after_nometa = nmol1,nmol1
        n_after_flagsandlod = nmol2
        n_after_duplicate,n_obs_after_anyvaliddata,n_obs_after_nokeymeta = nmol3,nmol3,nmol3
        n_obs_after_resolution = nmol4
    
    #update valid hours duplicated count
    test = full_data >= 0
    valid_hours_dup = np.sum(n_dup_array[test])
        
    return data_valid,full_data,valid_hours_dup,p_st_grid,p_mm_grid,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod,'na'


#DO SECONDARY QUALITY CONTROL PROCESS CHECKS ON SURFACE OBSEVRATIONAL DATA - IN MERGE STEP
def secondary_quality_control(site_ref,species,n_dup_array,flask_flag,instru_data,meta,full_data,n_obs_valid,output_set,process_group,start_year,end_year,grid_dates,datetime_cut,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap):

    data_valid = True

    lat,lon,alt,raw_class,tz = meta

    #GET ANTHROME CLASS FOR NON-RESTRAINED
    if 'N' in output_set:
        anthfile = '/work/home/db876/plotting_tools/AtmosChem_Tools/obs_process/obs_surface_process/anthro2_a2000.nc'
        anthload = Dataset(anthfile)
        class_valid,anthrome_class_name = anthrome_classify(anthload,[lat],[lon])

    if 'N' not in output_set:
        #-----------------------------------------------------------
        #RAW CLASS CHECK - REMOVE SITE IF URBAN IN CLASS
        #FOR AIRBASE,AQS,CAPMON,CASTNET,EANET,SEARCH: IF 'URBAN','TRAFFIC' OR 'INDUSTRIAL' IN RAW CLASS, SITE INVALID
        #FOR NAPS: IF RAW CLASS EQUAL TO C OR I, SITE INVALID
        if (process_group == 'AirBase') or (process_group == 'EPA AQS') or (process_group == 'CAPMON') or (process_group == 'CASTNET') or (process_group == 'EANET') or (process_group == 'SEARCH'):
            if ('urban' in raw_class.lower()) or ('traffic' in raw_class.lower()) or ('industrial' in raw_class.lower()):
                data_valid = False  
                inv_rawclass_count+=1
                print 'Site Invalid, raw class is urban.'
                exit_r = 'rawclass'
                return 'na',data_valid,'na',instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap
        elif (process_group == 'NAPS'):
            if ('i' == raw_class.lower()) or ('c' == raw_class.lower()):
                data_valid = False  
                inv_rawclass_count+=1
                print 'Site Invalid, raw class is urban.'
                exit_r = 'rawclass'
                return 'na',data_valid,'na',instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap
    n_obs_after_rawclass += n_obs_valid
    
    if 'N' not in output_set:
        #-----------------------------------------------------------
        #ANTHROME CLASS CHECK - REMOVE SITE IF 5 or more of SITES IN 3X3 GRID CENTERED ON OBS SITE URBAN
        anthfile = '/work/home/db876/plotting_tools/AtmosChem_Tools/obs_process/obs_surface_process/anthro2_a2000.nc'
        anthload = Dataset(anthfile)
        class_valid,anthrome_class_name = anthrome_classify(anthload,[lat],[lon])
        if class_valid == 'invalid':
            data_valid = False
            inv_anthromeclass_count+=1
            print 'Site Invalid, site classed as urban by anthrome map.'
            exit_r = 'anthromeclass'
            return 'na',data_valid,'na',instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap
    
        #----------------------------------------------------------
        # remove sites that have too high an average by year
        #all manually screened
        inv_refs = ['abat2ka71','abat2kl17','abat2m226','abat2sv24','abat2vk26','abat30402','abat30406','abat30603','abat30901','abat30992','abat4s184','abat4s401','abat4s406','abat4s412','abat4s413','abat4s414','abat60103','abat60104','abat60106','abat60114','abat60115','abat60119','abat60139','abat60140','abat60141','abat60145','abat60160','abat60161','abat60162','abat60172','abat60197','abat72126','abat82707','abba0034a','abba0035a',
                    'abbg0019a','abbg0056a','abbg0057a','abbg0058a','abbg0074a','abch0023a','abch0036a','abch0040a','abch0041a','abcz0asan','abcz0ljnm','abcz0molj','abcz0molo','abcz0skls','abcz0tovk','abcz0tozr','abcz0uchm','abcz0udcm','abdebb016','abdebb037','abdeby009','abdeby067','abdeby069','abdehe002','abdehh054','abdeni026','abdenw066','abderp042','abdest006','abdeth012','abdeth017','abee0021a','abes1180a','abes1364a', 'abes1400a','abes1538a','abes1851a','abfr0037a','abfr02005','abfr02037','abfr03014','abfr03038','abfr03047','abfr03062','abfr03066','abfr05070','abfr06018','abfr07002','abfr07003','abfr07035','abfr08003','abfr08711','abfr10007','abfr10014','abfr10016','abfr10026','abfr10041','abfr1054a','abfr1056a','abfr12003','abfr1250a','abfr13007','abfr14003','abfr14005','abfr14006','abfr14007','abfr14012','abfr15002','abfr15007','abfr16054','abfr16056','abfr16057','abfr18053','abfr22004','abfr22016','abfr23058','abfr23065','abfr23093','abfr23174','abfr23179','abfr23188','abfr24030','abfr25028','abfr25048','abfr26003','abfr26006','abfr27001','abfr28001','abfr28003','abfr28011','abfr28131','abfr32014','abfr33305','abfr41001','abgb0036r','abgb0040r','abgb0051a','abgb0962a','abgr0033a','abgr0042a','abgr0110r','abie0132a','abit0440a','abit0441a','abit0447a','abit0558a','abit0659a','abit0663a','abit0684a',                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                    'abit0692a','abit0709a','abit0782a','abit0838a','abit0862a','abit0906a','abit0982a','abit1042a','abit1061a','abit1065a','abit1103a','abit1156a','abit1222a','abit1225a','abit1238a','abit1292a','abit1328a','abit1340a','abit1439a','abit1459a','abit1468a','abit1557a','abit1571a','abit1670a','abit1725a','abit1845a','abit1868a','abit1919a',                                                                                                                                             
                    'abit1925a','abit1944a','abit2032a','abit2051a','ablt00001','ablt00012','ablt0021a','abmk0030a','abmk0048a','abnl00107','abnl00442','abno0075a','abpl0011a','abpl0021a','abpl0033a','abpl0039a','abpl0043a','abpl0051a','abpl0126a','abpl0115a','abpl0119a','abpl0122a','abpl0126a','abpl0134a','abpl0136a','abpl0144a','abpl0148a','abpl0162a','abpl0198a','abpl0221a','abpl0256a','abpl0276a','abpl0321a','abpt03101','abpt03103','abro0074a','abro0075a','abro0077a','abro0079a','abro0084a','abro0191a','abro0208a','abrs0030a','abrs0031a','abrs0033a','abrs0034a','abrs0039a','abrs0040a','abrs0043a','abrs0044a','abrs0045a','abrs0046a','abrs0048a','abse0053a','absk0003a','absk0004a','absk0008a','absk0015a','absk0028a','absk0047a','absk0051a','absk0267a','aqs060290011','aqs060651016','aqs180830004','aqs181290001','aqs320030007','aqs400219002','aqs400719003','aqs420110001','canaps051802','canaps060807','canaps061004','canaps065001','canaps090402','canaps090601','canaps101101','canaps101701','canaps103205','canaps104003','canaps104003','canaps104301','bsc','bur','car','log','pld','plv','sof','vrn',
                    'abat60180','abbg0071a','abes1811a','abfr01014','abfr03032','abfr07038','abfr12052','abfr14022','abfr14042','abfr14051','abfr15018','abfr17016','abfr19010','abfr19015','abfr20061','abfr24033','abfr26019','abfr29439','abfr36019','abfr37003','abfr39007','abfr39008','abfr39009','abit0824a','aqs320030023','canaps040302','canaps060809','canaps092901','canaps101501','canaps102401','canaps102701','aqs291830008','ablv0rim1','abse0054a','abse0055a',
                    'abbg0026a','abbg0038a','abdebb010','abdest074','abdeub020','abfi00066','abfi00424','abfi00582','abfr02026','abfr05086','abfr06134','abfr08401','abfr12053','abfr15001','abfr18049','abfr22022','abfr23181','abfr23182','abfr34034','abfr34046','abfr41002','abgb0929a','abgb0995a','abgb0999a','abgb1005a','abgr0120a','abie0140a','abit1307a','abit1663a','abit1729a','abit2054a','abpl0151a','abpl0152a','abpl0187a','abpl0210a','abro0126a','abse0030a','abse0080a','absi0033a','absk0031a','aqs050350005','aqs350451233','canaps060610','canaps061104','canaps061201','canaps061802','canaps062001','canaps090502','canaps102102','canaps102801','canaps107100','canaps119004','canaps010401','abit1311a','ch0031r','canaps060602','canaps060806','snl','abis0015a','canaps030201','canaps094301','abpl0212a','abno0080a','abro0072a','aqs481210034','aqs481830001','aqs371590021','abpt03091','abch0033a','abch0045a','abch0003r','abdenw068','abdenw081','abat4s407','abat2sp10','aqs480290059','abch0002r','nl0011r','aqs482030002','pay','abro0153a','abnl00738','abch0002r','ablt00044','abdehe042','abee0022a','abro0086a']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
        if site_ref.lower() in inv_refs:
            data_valid = False
            inv_anthromeclass_count+=1
            print 'Site Invalid. Site has been manually checked using maps and found to be urban/road located'
            exit_r = 'anthromeclass'
            return 'na',data_valid,'na',instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap
    
    n_obs_after_anthromeclass += n_obs_valid
    
    if 'N' not in output_set:
        #-----------------------------------------------------------
        #Remove sites above 1500m from sea level
        if np.float64(alt) >= 1500:
            data_valid = False
            inv_altitude_count+=1
            print 'Site Invalid, Site is over 1500m from sea level.' 
            exit_r = 'altitude'
            return 'na',data_valid,'na',instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap
    n_obs_after_altitude += n_obs_valid
    
    #-----------------------------------------------------------
    if 'S' in output_set:
        #REMOVE NIGHTTIME NO - do not do check for non-restrained datasets
        if (species == 'NO'):
            day_inds,night_inds = get_daynight_inds(datetime_cut,full_data,lat,lon)
            full_data[night_inds] = -99999
            valid_data = full_data[np.where(full_data>=0)]
            if len(valid_data) == 0:
                data_valid = False
                inv_night_count+=1
                print 'Site Invalid, No valid data after removal of nighttime NO'
                exit_r = 'night'
                return 'na',data_valid,'na',instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap
    
    valid_hours_dup = np.sum(n_dup_array[full_data >= 0])
    n_obs_valid = int(len(full_data[full_data >= 0]) - valid_hours_dup)
    n_obs_after_night += n_obs_valid
    
    #-----------------------------------------------------------
    #CHECK DATA REPRESENTIVENESS OF DATASET
    #IF ONLY HAVE 6 OR HOURS LESS IN EACH DAY OF YEAR THEN SET ALL VALS IN DAY AS INVALID
    #MAKE SURE DATA FOR YEAR IS FOR MORE THAN 1 DAY
    #DO NOT INCLUDE CO FLASK DATA IN THIS TEST
    #if species == 'NO':
    #    min_vals = 3
    #else:
    #    min_vals = 6
    min_vals = 6
    if (species == 'CO') & (flask_flag == 'Yes'):
        pass
    else:
        start = 0
        end = 24
        for y in range(len(grid_dates)/24):
            day_vals = full_data[start:end]
            valid_day_inds = np.where(day_vals >= 0)[0]
            #if length of valid hours in day is > 1 check representativeness
            #if n vals between 1 and 5 cannot be representative of day on 4 hour interval so invalid day
            if (len(valid_day_inds) >= 1) & (len(valid_day_inds) < min_vals):
                full_data[start:end] = -99999
                instru_data[start:end] = -1
            #if n vals >= 6 then can be representative of day but check spacing between points to make sure
            #minimum spacing is not greater than 4
            #* n_vals changes to 3 for NO (due to removal of night NO)
            elif len(valid_day_inds) >= min_vals:
                min_day_diff = np.min(np.diff(valid_day_inds))
                #if min hours between valid measurements in day is > 4 then day is invalid
                #add to inv count
                if min_day_diff > 4:
                    full_data[start:end] = -99999
                    instru_data[start:end] = -1
            #else no n valid in day so pass over
            else:
                pass
        
            start+=24
            end+=24
        
        #CHECK THERE IS AT LEAST 1 VALID DATAPOINT IN DATASET AFTER DATA REMOVAL
        valid_data = full_data[np.where(full_data>=0)]
        if len(valid_data) == 0:
            data_valid = False
            inv_representativeness_count+=1
            print 'Site Invalid, Data Representativeness is not sufficient.'
            exit_r = 'representativeness'
            return 'na',data_valid,'na',instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap
    
    valid_hours_dup = np.sum(n_dup_array[full_data >= 0])
    n_obs_valid = int(len(full_data[full_data >= 0]) - valid_hours_dup)
    n_obs_after_representativeness += n_obs_valid

    #screen extreme data
    
    #remove bad sites for all data types
    if species == 'O3': 
        inv_refs = ['abpt04005','aqs720330008','mk0007r','ang','fun','aqs180030001','aqs180591001','aqs180850001','aqs180891016','aqs180892008']
    elif species == 'NO': 
        inv_refs = ['abes0365a','abit1272a','each539','abit1103a','canaps090901','abit1901a','aqs320030073','aqs400019009','aqs220470009','aqs390090004','aqs250154002','abes0016r','abes1827a','gb0002r','gb0036r','abdenw065','abat60156']
    elif species == 'NO2-MOLYBDENUM': 
        inv_refs = ['aqs300870761','aqs181470008','aqs260430901','aqs260430902','mk0007r','abit0857a','abit0898a','abit1491a','abit1504a','aqs180190003','aqs371190034','aqs500070003','aqs100032002']
    elif species == 'NO2-PHOTOLYTIC': 
        inv_refs = []
    elif species == 'CO': 
        inv_refs = ['abes1256a','abes1328a','abes1860a','abes1875a','aqs051190007','abmk0030a','abmk0039a','abmk0045a','abmk0046a','abes1033a','abes1034a','abes1130a','abes1137a','abes1181a','abes1356a','abes1365a','abes1418a','abes1419a','abes1421a','abes1601a','abes1610a','abes1688a','abrs0036a','aqs060190008','aqs150030010','aqs560391012','abbg0070a','kmw']
    elif species == 'ISOP': 
        inv_refs = []
      
    if site_ref.lower() in inv_refs:
        data_valid = False
        inv_extreme_count+=1
        print 'Site Invalid. Site data has been determined to be bad'
        exit_r = 'extreme'
        return 'na',data_valid,'na',instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap

    #remove extreme data based on manually screening
    full_data[full_data < 0] = np.NaN
    if species == 'O3':
        inv_refs_limit = {'abdeub006':150,'abit1865a':150,'aqs060831012':175,'aqs380530002':100,'aqs450190046':150,'aqs450210002':150,'aqs450790021':150,'aqs550270001':150,'easa209':150,'de0011r':150,'glh':150,'puy':150,'canaps091401':125,'abes1131a':125,'abfr0250a':150,'abmk0031a':150,'abmk0037a':125,'aqs150031004':100,'aqs380171004':100,'aqs420070003':200,'aqs420730015':150,'aqs460990008':125,'aqs470651011':150,'aqs530570020':100,'aqs720770001':100,'aqs800020010':250,'abbg0041a':150,'abbg0043a':150,'abbg0053r':110,'abes1117a':150,'abes1193a':125,'abit0991a':150,'abit1360a':150,'abmk0034a':150,'abmk0039a':100,'abro0063a':150,'canaps050103':125,'aqs040038001':100,'nmy':50}
        inv_refs_inds = {'abit1407a':[250000,265000],'abmk0043a':[315000,319200],'ablv0rke2':[300000,320000],'abpt04003':[369114,380000],'abpt04003':[369114,380000],'aqs470651011':[100000,145000],'beo':[320000,355000],'canaps102701':[290000,311000]}
    elif species == 'NO':
        inv_refs_limit = {'abdeby068':1000,'abdehb005':1000,'abdeni059':200,'abpt03027':1000,'abpt03083':1000,'aqs261630019':1000,'aqs320031001':1000,'aqs320031001':1000,'aqs471250009':500,'aqs550790026':500,'fr0013r':100,'abat30407':200,'abba0001g':100,'abes1518a':200,'abes1810a:150,abit1729a':150,'abit1902a':150,'aqs010710020':100,'aqs050350005':300,'aqs380530002':50,'aqs560090819':100,'abat31906':200,'abbg0026a':100,'abes1810a':100,'abes1923a':100,'abit1061a':350,'abit1288a':350,'abit1729a':150,'abit2014a':200,'ablu0103a':150,'canaps129203':100,'abit0782a':600,'aqs060290014':700,'aqs060670006':700,'aqs295100062':500,'aqs540390004':700,'aqs550790041':600,'canaps010102':400,'aqs181270015':200,'abit1193a':100,'aqs550790041':600,'abbg0074a':400,'abit1425a':150,'aqs060711234':150,'aqs160390002':150,'aqs181270015':200,'abdebb018':500,'abes1472a':700,'abfr21016':600,'abit0837a':500,'abit1088a':600,'abrs0036a':600,'aqs060290014':600,'aqs060670006':600,'aqs170314002':600,'aqs250270019':600,'aqs295100062':500,'aqs360290005':800,'aqs540390004':800,'aqs550790041':600,'canaps010102':400,'canaps050102':500}
        inv_refs_inds = {'abat4s108':[226000,237000],'abrs0036a':[371754,371980],'aqs320030078':[280000,280944],'aqs400310647':[207000,208000],'aqs550591001':[110000,119800],'abes0373a':[333000,336183],'aqs230090103':[346450,347700]}
    elif species == 'NO2-MOLYBDENUM':
        inv_refs_limit = {'abat0hbg1':100,'abcz0uchm':100,'abcz0uvse':100,'abdehe042':100,'abdehe043':100,'abdeni059':100,'abdeub017':100,'abes1831a':100,'abfr39007':100,'abpt04002':80,'abro0079a':300,'aqs021221004':60,'aqs160390002':60,'aqs160770011':60,'aqs181270015':100,'aqs181470002:100,aqs181631001':150,'aqs181631002':150,'aqs210150007':150,'aqs220470002':100,'aqs220470012':100,'aqs391291001':100,'aqs401010167':150,'aqs420990301':150,'aqs470430009':100,'aqs471131001':100,
                          'de0043g':100,'es0005r':100,'es0006r':100,'es0007r':100,'es0008r':100,'es0010r':100,'es0011r':100,'es0012r':100,'es0013r':100,'es0014r':100,'es0016r':100,'es0017r':100,'fr0013r':100,'gb0002r':100,'gb0013r':100,'gb0015r':100,'gb0031r':100,'gb0033r':100,'gb0037r':100,'gb0038r':100,'gb0043r':100,'gb0045r':100,'gb0050r':100,'gb0051r':100,'gb0053r':100,'gr0001r':100,'nl0007r':100,'nl0010r':100,'nl0011r':100,'nl0091r':100,'nl0644r':100,'canaps070203':200,'canaps090806':100,
                          'canaps090901':125,'canaps091401':100,'canaps091601':100,'abat4s406':150,'abat4s418':100,'abat60141':100,'abat60145':100,'abat72110':150,'abcz0arie':200,'abcz0avys':200,'abcz0bbnf':100,'abcz0uulk':100,'abcz0zsnv':150,'abdeby068':150,'abdehh002':100,'abdehh007':100,'abdehh008':150,'abdehh009':150,'abdehh010':100,'abdehh014':150,'abdehh015':150,'abdehh022':150,'abdenw116':150,'abdesn005':150,'abdk0051a':200,'abes0886a':100,'abes0915a':200,'abes1054a':200,'abes1139a':200,'abes1172a':200,'abes1365a':200,'abes1846a':150,
                          'abfr30034':300,'abit1381a':150,'abpt01021':200,'abpt03075':260,'abrs0028a':150,'aqs080051002':195,'aqs180891016':150,'aqs211451024':125,'aqs250251003':200,'aqs260990009':100,'aqs270530957':150,'aqs291890001':200,'aqs320031001':350,'aqs390614002':400,'aqs470110102':175,'aqs530570018':100,'aqs550590016':250,'canaps010102':175,'aqs181470002':200,'abch0017a':100,'abcz0ultt':100,'abdebw010':100,'abdehe003':150,'abdehe036':150,'abes0095a':80,'abit1044a':150,'abit1050a':100,'abpt03027':300,'aqs090031003':150,'aqs120864002':250,'aqs132470001':100,'aqs480290046':100,'aqs800020017':150,'abch0017a':100,'abcz0ppla':80,'abcz0ultt':100,'abdebe056':100,'abdebw008':125,'abdebw010':100,'abdebw072':100,'abdebw093':90,'abdeby078':100,'abdeby116':100,'abdehe003':150,'abdehe004':150,'abdehe036':150,'abdest002':150,'abdest025':100,
                          'abdest081':200,'abdk0014a':100,'abes0095a':100,'abes1130a':100,'abes1163a':150,'abes1596a':100,'abes1630a':120,'abes1643a':150,'abes1865a':150,'abfi00031':150,'abfr10024':150,'abfr11032':100,'abfr18089':100,'abfr19003':100,'abfr24017':80,'abhu0020a':125,'abie0028a':100,'abie0036a':125,'abie0135a':200,'abis0005a':150,'abit0714a':150,'abit0781a':150,'abit1269a':125,'abit1273a':60,'abmk0028a':150,'abmk0037a':100,'abno0058a':150,'abpt01022':150,'abpt03027':300,'abpt04001':100,'abpt05010':100,'abro0079a':150,'aqs061130004':100,'aqs120864002':250,'aqs132470001':100,'aqs160050015':100,'aqs171190017':150,'aqs181631001':150,'aqs211010013':150,'aqs245100040':200,'aqs250270020':200,'aqs261630019':250,'aqs291890014':100,'aqs300870700':100,'aqs320030016':200,'aqs360010012':100,'aqs360610010':300,'aqs390290016':200,'aqs390350043':200,'aqs480290046':100,'aqs550890005':100}
        inv_refs_inds = {'abbg0074a':[350000,352000],'abcz0zsnv':[210000,212000],'abbg0013a':[300000,322000],'dk0005r':[377000,390000],'dk0008r':[370000,390000],'dk0012r':[377000,390000],'abes0365a':[309100,325500],'abes0624a':[320000,328400],'abes0625a':[310000,319500],'abes1741a':[311400,315500],'abes1748a':[300000,315700],'abno0076a':[310000,324600],'aqs060290010':[271000,276100],'aqs060731002':[100000,150000],'aqs060850004':[200000,208200],'aqs370670022':[131800,230000]}   
    elif species == 'NO2-PHOTOLYTIC':
        inv_refs_limit = {'ablv0rim1':100,'searchctr':50,'abro0058a':175,'ch0001g':40,'searchpns':80,'abro0058a':120,'abro0061a':90}
        inv_refs_inds = {'abro0061a':[325000,332775]}
    elif species == 'CO':
        inv_refs_limit = {'abes1599a':2000,'abes1660a':1500,'abfr16057':2500,'abro0072a':7000,'abro0086a':2500,'aqs230090103':500,'aqs240190004':500,'aqs240230002':750,'aqs330115001':700,'aqs550270001':2000,'abbetr512':10000,'abbg0012a':12500,'abdesl010':5000,'abes0817a':30000,'abes0818a':20000,'abes0917a':30000,'abes1250a':30000,'abes1251a':20000,'abes1254a':20000,'abes1260a':20000,'abes1340a':20000,'abes1364a':20000,'abes1394a':20000,'abes1425a':20000,'abes1675a':1000,'abit1199a':4000,'abit1263a':15000,'abpt01031':5000,'abro0068a':15000,'abro0159a':5000,'aqs800020010':15000,'abnl00230':3000,'abcz0avyn':4000,'abmk0044a':15000,'abbetr701':6000,'abbetb004':6000,'abbetr701':6000,'abdebb032':3000,'abderp024':3000,'abee0020a':3000,'bes1050a':3000,'abes1353a': 4000,'abes1615a':8000,'abes1499a':4000,'abfr19007':6000,'abit1246a':8000,'abit1233a':4000,'abhu0041a':4000,'abro0069a':15000,'aqs090030017':6000,'aqs130890002':6000,'abro0071a':8000,'aqs490353006':6000,'abat60172':3000,'mnm':600,'abpt06004':5000,'aqs380171004':3000,'glh':1000}
        inv_refs_inds = {'abderp002':[100000,150000],'abbetn043':[315400,333000],'abbetwol1':[315500,325000],'abes1417a':[305000,320000],'abes1821a':[341000,344300],'abat2ka21':[236578,245734],'abat51000':[236578,245734],'abat90mba':[220000,245734],'abbetb005':[320000,333194],'abbetr001':[315400,325000],'abbetr002':[315400,325000],'abbetr223':[315400,333200],'abbetr501':[315400,333200],'abbetr512':[315400,325000],'abcz0topr':[298000,300272],'abdeth072':[350000,380000],'abes1616a':[357223,363310],'abit1339a':[260000,275000],'aqs330115001':[355000,367551],'aqs461270001':[381280,390000],'abes1633a':[294586,299806],'abbetr222':[310000,333100],'abbetr501':[260000,266800],'abbetr512':[260000,266800],'abbetr750':[315500,326000],'abdest015':[330000,351265],'abdest076':[330000,350000],'abdeth020':[359818,380000],'abee0018a':[300000,315494],'abee0019a':[300000,315494],'abes0124a':[356436,360000],'abes1320a':[300000,320000],'abes1615a':[280000,317200],'abes1624a':[309000,350000],'abfr29429':[250000,269417],'abpl0026a':[270000,290000],'abpt03010':[240000,270863],'abpt03071':[240000,270069],'abpt03072':[240000,267394],'abpt03075':[240000,266164],'aqs110010043':[384232,385416],'aqs220330009':[320000,327231],'aqs250250042':[330000,339335] ,'aqs330150018':[355000,364326],'aqs350010023':[355000,361000],'aqs360810124':[345430,350712],'aqs371830014':[340000,346153],'aqs410510080':[330000,339714],'aqs481130069':[330000,336000],'aqs482451035':[320000,332021],'ee0009r':[385000,395000]}
    elif species == 'ISOP':
        inv_refs_limit = {'ch0005r':80}
        inv_refs_inds = {}
        
    #remove extreme data - checked manually
    if site_ref.lower() in inv_refs_limit:
        full_data[full_data >= inv_refs_limit[site_ref.lower()]] = np.NaN
        
    #remove strange data period - checked manually
    if site_ref.lower() in inv_refs_inds:
        start_i = inv_refs_inds[site_ref.lower()][0]
        end_i = inv_refs_inds[site_ref.lower()][1]
        full_data[start_i:end_i] = np.NaN
    
    #do counts after extreme removal
    valid_hours_dup = np.sum(n_dup_array[~np.isnan(full_data)])
    n_obs_valid = int(len(full_data[~np.isnan(full_data)]) - valid_hours_dup)
    n_obs_after_extreme += n_obs_valid
    #make nans again -99999
    full_data[np.isnan(full_data)] = -99999

    #-----------------------------------------------------------
    if 'T' not in output_set:
        #CHECK PARTIAL YEARS, only don't check for trend set
        #check there is gaps > 2 months in each year
        #if is first year then see if next year is invalid before removing data

        first_year_flag = False
        inv_year_count = 0
        year_range = range(start_year,end_year)
        all_years = [i[:4] for i in grid_dates]
        all_years = np.array(all_years)
        for i in range(len(year_range)):
            gap_count = 0
            max_gap_count = 0
            year_test = all_years == str(year_range[i])
            year_vals = full_data[year_test]
        
            #test data 
            #if year is all invalid move to next year
            if len(year_vals[year_vals >= 0]) == 0:
                continue
            else:
                for j in year_vals:
                    if j <= 0:
                        gap_count+=1
                        if gap_count > max_gap_count:
                            max_gap_count = gap_count
                    else:
                        gap_count = 0
    
            #check if have gap of > 2 months in year of data
            #if first year then see if next year is invalid too before removing data
            if max_gap_count > 1392:
                if 'P' not in output_set:
                    if i == 0:
                        first_year_flag = True
                    elif i == 1:
                        if first_year_flag == True:
                            year_test = all_years == str(year_range[0])
                            full_data[year_test] = -99999
                        year_test = all_years == str(year_range[1])
                        full_data[year_test] = -99999
                    else:
                        full_data[year_test] = -99999
                else:
                    inv_year_count+=1
    
        #for periodic output if more than half of years have gaps then all data is invalid 
        if 'P' in output_set:
            if inv_year_count >= (len(year_range)/2):
                full_data[:] = -99999
            
        #CHECK THERE IS AT LEAST 1 VALID DATAPOINT IN DATASET AFTER PARTIAL YEAR REMOVAL
        valid_data = full_data[np.where(full_data>=0)]
        unique_vals = np.unique(full_data[full_data>=0])
        if (len(valid_data) == 0) or (len(unique_vals) < 3):
            data_valid = False
            inv_partialyear_count+=1
            print 'Site Invalid. Site has no years with valid data after remove partial data years'
            exit_r = 'partialyear'
            return 'na',data_valid,'na',instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap
    
    valid_hours_dup = np.sum(n_dup_array[full_data >= 0])
    n_obs_valid = int(len(full_data[full_data >= 0]) - valid_hours_dup)
    n_obs_after_partialyear += n_obs_valid
    
#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------
    #if output type contains P: 'PERIODIC', CHECK TIMEZONE IS EVEN
    #if all timezones recorded are not whole number, site is invalid
    if ('P' in output_set):
        if (float(tz) % 1) != 0:
            data_valid = False
            inv_timezone_count+=1
            print 'Site Invalid, all recorded timezones for site are not whole numbers.'
            exit_r = 'timezone'
            return 'na',data_valid,'na',instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap
    n_obs_after_timezone += n_obs_valid
    
#-----------------------------------------------------------
    #if output type contains P: 'PERIODIC' THEN TEST FOR DATA GAPS > 1 YEAR IN TOTAL DATA RECORD 
    if ('P' in output_set):
        max_count = 0
        count = 0
        for i in range(len(full_data)):
            if full_data[i] <= 0:
                count+=1
                if count > max_count:
                    max_count = count
            else:
                count = 0

        if max_count >= 8760:
            data_valid = False
            inv_bigdatagap_count+=1
            print 'Site Invalid. Site has data gap > 1 year'
            exit_r = 'bigdatagap'
            return 'na',data_valid,'na',instru_data,exit_r,inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap
    n_obs_after_bigdatagap += n_obs_valid

    return full_data,data_valid,anthrome_class_name,instru_data,'na',inv_rawclass_count,n_obs_after_rawclass,inv_anthromeclass_count,n_obs_after_anthromeclass,inv_altitude_count,n_obs_after_altitude,inv_night_count,n_obs_after_night,inv_representativeness_count,n_obs_after_representativeness,inv_extreme_count,n_obs_after_extreme,inv_partialyear_count,n_obs_after_partialyear,inv_timezone_count,n_obs_after_timezone,inv_bigdatagap_count,n_obs_after_bigdatagap

def remove_duplicate_points(ref,time,var,mm,st,dup,output_res):
    #Check if any times are duplicate or overlapping
    #if duplicate delete all but first
    #if overlap, then sort list and then check for duplicates
    
    #check for overlaps
    check_over = all(time[i] <= time[i+1] for i in xrange(len(time)-1))
    if check_over == False:
        print 'Times overlap'
        #print num2date(time[:200], units='days since 1970-01-01 00:00:00', calendar='gregorian')
        #for x in range(len(time)):
        #    for y in range(x+1,len(time)):
        #        if time[y] < time[x]:
        #            print num2date(time[x], units='days since 1970-01-01 00:00:00', calendar='gregorian'),num2date(time[y], units='days since 1970-01-01 00:00:00', calendar='gregorian'),x,y
                    #break
        
        if output_res == 'HDM':
            file = open("overlapdata.txt", "a")
            file.write("%s\n"%(ref))
            file.close()
        #sort all lists in ascending order
        var = np.array([v for (t,v) in sorted(zip(time,var))])
        mm = np.array([v for (t,v) in sorted(zip(time,mm))])
        st = np.array([v for (t,v) in sorted(zip(time,st))])
        if dup != 'blank':
            dup = np.array([v for (t,v) in sorted(zip(time,dup))])
        time = np.sort(time)
    
    #check for duplicates
    del_list_dup = defaultdict(list)
    for i,item in enumerate(time):
        del_list_dup[item].append(i)
    del_list_dup = {k:v for k,v in del_list_dup.items() if len(v)>1}
    #get all indices except the first matched as duplicate
    try:
        del_list_d = del_list_dup.values()[0][1:]
    except:
        del_list_d = []
    
    if len(del_list_dup) > 0:
        print 'Deleting duplicate timepoints' 
        print np.sort(list(del_list_dup.keys()))
        a = np.sort(list(del_list_dup.keys()))
        print num2date(a,units='days since 1970-01-01 00:00:00', calendar='gregorian')
        
        time = np.delete(time,del_list_d)
        var = np.delete(var,del_list_d)
        mm = np.delete(mm,del_list_d)
        st = np.delete(st,del_list_d)
        if dup != 'blank':
            dup = np.delete(dup,del_list_d)
        if output_res == 'HDM':
            file = open("duplicatedata.txt", "a")
            file.write("%s\n"%(ref))
            file.close()
        
    return time,var,mm,st,dup
        

def remove_duplicate_sites(refs,lats,lons,alts,all_data,native_res,all_mm):
    
    match_refs = []
    match_refs_i = []
    
    #round all to 2 dp
    round_lats = np.around(lats,decimals=2)
    round_lons = np.around(lons,decimals=2)
    round_alts = np.around(alts,decimals=2)
    
    #find refs (and inds) whose lat/lons are equivalent at 2dp and altitudes within 50m,
    #and have same measurement methods(as lists per site)
    for i in range(len(refs)):
        match = False
        match_refs_list = [refs[i]]
        match_refs_i_list = [i]
        
        current_lat = round_lats[i]
        current_lon = round_lons[i]
        current_alt = round_alts[i]
        current_mm = all_mm[i]
    
        for j in range(i+1,len(refs)):
            next_lat = round_lats[j]
            next_lon = round_lons[j]
            next_alt = round_alts[j]
            next_mm = all_mm[j]
            
            diff_alt = np.abs(next_alt - current_alt)
            
            if (current_lat == next_lat) & (current_lon == next_lon) & (diff_alt < 50) & (current_mm == next_mm):
                match = True
                match_refs_list.append(refs[j])
                match_refs_i_list.append(j)
                
        if match == True:
            match_refs.append(match_refs_list)
            match_refs_i.append(match_refs_i_list)
        
    #iterate through duplicate lists
    #preferentially keep ref with lowest native resolution
    #then test which ref has most valid data
    #if all same, then keep first ref in list
    del_inds = []
    for x in range(len(match_refs)):
        print match_refs[x]
        current_del_inds = []
        keep_inds = []
        native_list = []
        len_list = []
        start_hours = []
        end_hours = []
        for z in range(len(match_refs[x])):
            current_data = np.array(all_data[match_refs_i[x][z]])
            current_res = native_res[match_refs_i[x][z]]
            test = current_data != -99999
            valid_len = len(current_data[test])
            native_list.append(current_res)
            len_list.append(valid_len)
            #print match_refs[x][z]
            #print len(current_data)
            #print len(current_data[test])
            #print np.max(current_data)
            
            start_hours.append(next(x[0] for x in enumerate(current_data) if x[1] >= 0))
            end_hours.append(len(all_data[i]) - next(x[0] for x in enumerate(current_data[::-1]) if x[1] >= 0))
        
        current_del_inds = np.array(current_del_inds)
        #preferentially keep refs with lowest native resolution
        if ('H' in native_list) & ('D' in native_list):
            current_del_inds=np.append(current_del_inds,np.where(np.array(native_list) == 'D')[0])
        if ('H' in native_list) & ('M' in native_list):
            current_del_inds=np.append(current_del_inds,np.where(np.array(native_list) == 'M')[0])
        if ('D' in native_list) & ('M' in native_list):
            current_del_inds=np.append(current_del_inds,np.where(np.array(native_list) == 'M')[0])
        current_del_inds = current_del_inds.astype(int)
        current_del_inds = list(set(current_del_inds))
        
        #test if start/end times both are significant different (> 5 year (26280 hours)) and data overlap not greater than 50%
        #if so keep site
        for s in range(len(start_hours)):
            current_start = start_hours[s]
            current_end = end_hours[s]
            current_range = current_end-current_start
            
            it_list = range(len(start_hours))
            it_list.remove(s)
            
            flag = True
            for t in it_list:
                next_start = start_hours[t]
                next_end = end_hours[t]
                next_range = next_end-next_start
                #calculate overlap relative to current list range
                if (next_start >= current_start) & (next_start <= current_end):
                    overlap_range = current_end-next_start
                    if (100./current_range)*overlap_range > 50:
                        flag = False
                    
                elif (next_end >= current_start) & (next_end <= current_end):
                    overlap_range = next_end-current_start
                    if (100./current_range)*overlap_range > 50:
                        flag = False
                
                elif (next_end <= current_end) & (next_start >= current_start):
                    overlap_range = next_end-next_start
                    if (100./current_range)*overlap_range > 50:
                        flag = False
                    
                elif (current_end <= next_end) & (current_start >= next_start):
                    overlap_range = current_end-current_start
                    if (100./current_range)*overlap_range > 50:
                        flag = False
                
                #no overlap
                else:
                    pass
                    
            if flag == True:
                if s not in current_del_inds:
                    keep_inds.append(s)
            
        keep_inds = list(set(keep_inds))
        
        #next test which ref has most valid data
        if len(match_refs[x]) - (len(current_del_inds)) > 1:
            sort_valid_i = np.argsort(len_list)[::-1]
            got_ref = False
            for a in range(len(len_list)):
                if (sort_valid_i[a] not in current_del_inds) & (sort_valid_i[a] not in keep_inds):
                    if got_ref == True:
                        current_del_inds.append(sort_valid_i[a])
                    got_ref = True
            del_inds = np.append(del_inds,np.array(match_refs_i[x])[current_del_inds])
        else:
            del_inds = np.append(del_inds,np.array(match_refs_i[x])[current_del_inds])
               
    del_inds = np.sort(list(set(del_inds))).astype(int)

    return del_inds

        
def mad_based_outlier(points, thresh=5.):
    median = np.nanmedian(points)
    diff = points - median
    med_abs_deviation = np.nanmedian(np.abs(diff))
    modified_z_score = 1.4826 * med_abs_deviation
    #lower_limit = median-(modified_z_score*thresh)
    upper_limit = median+(modified_z_score*thresh)
    #print lower_limit
    #lower_test = points < lower_limit 
    upper_test = points > upper_limit
    return upper_test


#1.07
#----------

#TAKE AVERAGE OF HIGHER FREQUENCY OBS DATA TO GET DAILY OR MONTHLY OUTPUT

def take_average(full_data,obs_time_pd,output_res):

    #make sure all invalids set to np.nan
    invalid_test = full_data < 0
    full_data[invalid_test] = np.nan
    
    #take average of dataset
    obs_var_pd = pd.Series(full_data, index=obs_time_pd)

    if output_res == 'D':
        obs_key = lambda x: pd.Period(str(x.year)+'-'+str(x.month)+'-'+str(x.day))
    elif output_res == 'M':
        obs_key = lambda x: pd.Period(str(x.year)+'-'+str(x.month))
    elif output_res == 'Y':
        obs_key = lambda x: pd.Period(str(x.year))
        
    #group the data
    var_group=obs_var_pd.groupby(obs_key)
    raw_st_group = raw_st_pd.groupby(obs_key)

    #calculate the mean and write out
    obs_mean = var_group.mean()
    full_data = np.array(obs_mean.values.tolist())

    #move nan's back to -99999
    invalid_test = np.isnan(full_data)
    full_data[invalid_test] = -99999
        
    return full_data
    
#1.08
#----------

#WRITE OUT ALL PRIMARY OBS PROCESS DATA TO NETCDF
    
def write_out_data(site_ref,process_group,root_grp,species,full_data,p_st_grid,p_mm_grid,data_valid,meta,n_dup):
    
    if data_valid == True:
        #unpack meta
        lat,lon,alt,raw_class_name,file_res,unit,p_unit,data_tz,local_tz,site_name,country,contact = meta       

        #save out netcdf file
        ref = root_grp.createGroup('%s'%(site_ref.lower()))

        #set variables
        spec = ref.createVariable(species.lower(), 'f8', ('flex',))
        mm = ref.createVariable('measurement method', 'int', ('flex',))
        dup = ref.createVariable('n_duplicate_array', 'int', ('flex',))
        
        #set core meta attributes
        ref.site_name = site_name
        ref.country = country
        ref.data_contact = contact
        ref.latitude = lat
        ref.longitude = lon
        ref.altitude = alt
        ref.data_timezone = data_tz
        ref.local_timezone = local_tz
        ref.process_group = process_group
        ref.raw_site_class = raw_class_name
        ref.raw_units = unit
        ref.processed_units = p_unit
        ref.native_resolution = file_res
        ref.contact = contact
        if (3 in p_st_grid) & (len(list(set(p_st_grid)))==2):
            ref.flask_flag = 'Yes'
        elif (3 in p_st_grid) & (len(list(set(p_st_grid)))>2):
            1+'a'
        else:
            ref.flask_flag = 'No'
        
        spec[:] = full_data
        dup[:] = n_dup
        mm[:] = p_mm_grid
    
        print 'Site is Valid\n'
    
    else:
        print 'Site is Invalid\n'   
    
    return
    
def grid_data_badc(raw_species,species):
    rawdatafile = 'GLOBAL_SURFACE_%s_1970_2015_H_HDMS.nc'%(raw_species)

    #1x1 centered on degrees
    print '1x1a'
    lonarr = np.arange(-180.5,180.5,1)
    latarr = np.arange(-90.5,91.5,1)
    filename='FINAL_GRID/Surface%sObs_1x1a_Monthly_v1.0.nc'%(species.upper()) 
    filenameA='FINAL_GRID/Surface%sObs_1x1a_Annual_v1.0.nc'%(species.upper())  
    do_grid(raw_species,species, rawdatafile, filename, lonarr=lonarr, latarr=latarr, timeres = 'MS')
    do_grid(raw_species,species, rawdatafile, filenameA, lonarr=lonarr, latarr=latarr, timeres = 'AS')

    # 2x2.5 GC
    print '2x2.5GC'
    lonarr = np.arange(-181.25,178.75+2.5,2.5)
    latarr = np.arange(-91,92,2.0)
    latarr[0]=-90.
    latarr[-1]=90.
    filename='FINAL_GRID/Surface%sObs_2x25_GEOS-Chem_Monthly_v1.0.nc'%(species.upper())  
    filenameA='FINAL_GRID/Surface%sObs_2x25_GEOS-Chem_Annual_v1.0.nc'%(species.upper())  
    do_grid(raw_species,species, rawdatafile, filename, lonarr=lonarr, latarr=latarr, timeres = 'MS')
    do_grid(raw_species,species, rawdatafile, filenameA, lonarr=lonarr, latarr=latarr, timeres = 'AS')

    #4x5 GC
    print '4x4.5GC'
    lonarr = np.arange(-182.5,177.5+2.5,5)
    latarr = np.arange(-92,92+1,4.)
    latarr[0]=-90.
    latarr[-1]=90.
    filename='FINAL_GRID/Surface%sObs_4x5_GEOS-Chem_Monthly_v1.0.nc'%(species.upper())  
    filenameA='FINAL_GRID/Surface%sObs_4x5_GEOS-Chem_Annual_v1.0.nc'%(species.upper())   
    do_grid(raw_species,species, rawdatafile, filename, lonarr=lonarr, latarr=latarr, timeres = 'MS')
    do_grid(raw_species,species, rawdatafile, filenameA, lonarr=lonarr, latarr=latarr, timeres = 'AS')

    # 4x5 GCAP
    print '4x4.5GCAP'
    lonarr = np.arange(-182.5,177.5+2.5,5)
    latarr = np.arange(-90,92+1,4.)
    filename='FINAL_GRID/Surface%sObs_4x5_GCAP_Monthly_v1.0.nc'%(species.upper())  
    filenameA='FINAL_GRID/Surface%sObs_4x5_GCAP_Annual_v1.0.nc'%(species.upper())   
    do_grid(raw_species,species, rawdatafile, filename, lonarr=lonarr, latarr=latarr, timeres = 'MS')
    do_grid(raw_species,species, rawdatafile, filenameA, lonarr=lonarr, latarr=latarr, timeres = 'AS')

    # 2x2.5 GISS
    print '2x2.5GISS'
    lonarr = np.arange(-178.75,182,2.5)
    latarr = np.arange(-90,91,2)
    filename='FINAL_GRID/Surface%sObs_2x25_GISSE_Monthly_v1.0.nc'%(species.upper())  
    filenameA='FINAL_GRID/Surface%sObs_2x25_GISSE_Annual_v1.0.nc'%(species.upper()) 
    do_grid(raw_species,species, rawdatafile, filename, lonarr=lonarr, latarr=latarr, timeres = 'MS')
    do_grid(raw_species,species, rawdatafile, filenameA, lonarr=lonarr, latarr=latarr, timeres = 'AS')

    #UKCA 2.5x3.75
    print '2.5x3.75UKCA'
    lonarr = np.arange(-180,182,3.75)
    latarr = np.arange(-90,91,2.5)
    filename='FINAL_GRID/Surface%sObs_25x375_UKCA_Monthly_v1.0.nc'%(species.upper())  
    filenameA='FINAL_GRID/Surface%sObs_25x375_UKCA_Annual_v1.0.nc'%(species.upper())  
    do_grid(raw_species,species, rawdatafile, filename, lonarr=lonarr, latarr=latarr, timeres = 'MS')
    do_grid(raw_species,species, rawdatafile, filenameA, lonarr=lonarr, latarr=latarr, timeres = 'AS')

    #ACCMIP
    print '2x2ACCMIP'
    lonarr = np.arange(-180,182,2.)
    latarr = np.arange(-90,91,2.)
    filename='FINAL_GRID/Surface%sObs_2x2_ACCMIP_Monthly_v1.0.nc'%(species.upper()) 
    filenameA='FINAL_GRID/Surface%sObs_2x2_ACCMIP_Annual_v1.0.nc'%(species.upper()) 
    do_grid(raw_species,species, rawdatafile, filename, lonarr=lonarr, latarr=latarr, timeres = 'MS')
    do_grid(raw_species,species, rawdatafile, filenameA, lonarr=lonarr, latarr=latarr, timeres = 'AS')

    #1x1 centered on half-degrees
    print '1x1'
    lonarr = np.arange(-180,180+1,1)
    latarr = np.arange(-90,90+1,1)
    filename='FINAL_GRID/Surface%sObs_1x1_Monthly_v1.0.nc'%(species.upper())  
    filenameA='FINAL_GRID/Surface%sObs_1x1_Annual_v1.0.nc'%(species.upper()) 
    do_grid(raw_species,species, rawdatafile, filename, lonarr=lonarr, latarr=latarr, timeres = 'MS')
    do_grid(raw_species,species, rawdatafile, filenameA, lonarr=lonarr, latarr=latarr, timeres = 'AS')

def do_grid(raw_species,species, rawdatafile,filename,lonarr=None, latarr=None, timeres = None):
    
    #read in obs data netcdf
    refs,raw_time,ref_time,datetime_time,obs_data,obs_lats,obs_lons,obs_alt,obs_process_group,obs_raw_class,obs_anthrome_class,gap_inds = read_obs_all(rawdatafile,raw_species,1970,2015)

    #output extra metadata - not typically used
    root_grp = Dataset(rawdatafile)
    
    obs_contact = []
    obs_site_name = []
    obs_data_tz = []
    obs_local_tz = []
    obs_raw_units = []
    obs_processed_units = []
    obs_native_res = []
    obs_methods = []
    
    for ref in refs:
        ref_group = root_grp.groups[ref]
        obs_contact = np.append(obs_contact,ref_group.data_contact)
        obs_site_name = np.append(obs_site_name,ref_group.site_name)
        obs_data_tz = np.append(obs_data_tz,ref_group.data_timezone)
        obs_local_tz = np.append(obs_local_tz,ref_group.local_timezone)
        obs_raw_units = np.append(obs_raw_units,ref_group.raw_units)
        obs_processed_units = np.append(obs_processed_units,ref_group.processed_units)
        obs_native_res = np.append(obs_native_res,ref_group.native_resolution)
        obs_methods.append(ref_group.variables['measurement method'][:])

    #set metadata list
    meta = np.array([refs,obs_site_name,obs_lats,obs_lons,obs_alt,obs_raw_class,obs_anthrome_class,obs_process_group,obs_data_tz,obs_local_tz,obs_raw_units,obs_processed_units,obs_native_res,obs_contact])
    meta_str = np.array(['site_id','site_name','latitude','longitude','altitude','raw_class','secondary_class','data_network','data_timezone','local_timezone','raw_units','processed_units','native_resolution','data_contact'])
    #print meta.shape

    datetime_out = pd.date_range(start=datetime.datetime(1980,1,1,0,0),end=datetime.datetime(2015,1,1,0,0), freq = timeres)[:-1]
    #print datetime_out
    
    longrid,latgrid = np.meshgrid(lonarr,latarr)
    longrid = longrid.T
    latgrid = latgrid.T
    
    #set output lat,lon,time (in seconds from 1980)
    lon_c = lonarr[:-1] + (np.diff(lonarr)/2.)
    lat_c = latarr[:-1] + (np.diff(latarr)/2.)
    t_out = []
    for tt in range(len(datetime_out)):
        t_out = np.append(t_out,(datetime_out[tt]-datetime.datetime(1980,1,1,0,0)).total_seconds())
         
    #if grid longitudes less than -180 then convert longitude accordingly
    #test = obs_lons
    #if np.amax(lonarr)>180. and np.amin(lonarr)>=0.:
    #    datalon = (datalon+360.0)%360
    
    #-------------------------------------------------------------------------
    #create netcdf
    root_grp = Dataset(filename,'w')
    
    #set attributes
    import time as t
    root_grp.comment = 'This netCDF file holds gridded surface %s data from multiple databases (AirBase, CASTNET, EANET, EMEP, EPA AQS, NAPS, SEARCH, WMO GAW), plus metadata.'%(species)
    root_grp.title = 'Global Gridded Surface %s Dataset'%(species)
    root_grp.cdm_data_type = 'Grid'
    root_grp.creator_name = 'Dr. Dene R. Bowdalo, Prof. Mathew J. Evans'
    root_grp.creater_email = 'denebowdalo@googlemail.com, mat.evans@york.ac.uk'
    root_grp.institution = 'University of York'
    root_grp.conventions = 'CF-1.6'
    root_grp.date_modified = t.ctime()
    root_grp.history= t.ctime()
    root_grp.data_version= '1.0'
    root_grp.id = 'global_gridded_%s'%(species)
    root_grp.reference = 'Dene R. Bowdalo, Mathew J. Evans, and Eric D. Sofen., Gridded global surface ozone precursor metrics for global atmospheric chemistry model evaluation., Submitted to ESSD.'
    root_grp.acknowledgment = 'Funded by the NERC BACCHUS and CAST projects NE/L01291X/1, NE/J006165/1. Data provided by various providers specified in metadata.'
    root_grp.source = 'Gridded metrics derived from hourly surface observations.'
    root_grp.licence = 'NERC Data Use Policy.  Freely available for use with attribution.'
    
    #set dims
    root_grp.createDimension('stations', len(refs))
    root_grp.createDimension('metadata', len(meta_str))
    root_grp.createDimension('time', len(datetime_out))
    root_grp.createDimension('longitude', len(lon_c))
    root_grp.createDimension('latitude', len(lat_c))
    root_grp.createDimension('methods', None)
    
    #create variables
    t = root_grp.createVariable('time','int', ('time',))
    lon = root_grp.createVariable('longitude','f8', ('longitude',))
    lat = root_grp.createVariable('latitude','f8', ('latitude',))
    m1 = root_grp.createVariable('meta_fields',str, ('metadata',))
    m2 = root_grp.createVariable('meta_values',str, ('metadata','stations',))
    mtd1 = root_grp.createVariable('methods',str, ('methods','stations',))
    mtd2 = root_grp.createVariable('methods_times',str, ('methods','stations',))
    pri1 = root_grp.createVariable('Mean_Gridded','f8', ('time','longitude','latitude',))
    pri2 = root_grp.createVariable('Median_Gridded','f8', ('time','longitude','latitude',))
    pri3 = root_grp.createVariable('Skewness_Gridded','f8', ('time','longitude','latitude',))
    pri4 = root_grp.createVariable('Kurtosis_Gridded','f8', ('time','longitude','latitude',))
    pri5 = root_grp.createVariable('25_Percentile_Gridded','f8', ('time','longitude','latitude',))
    pri6 = root_grp.createVariable('75_Percentile_Gridded','f8', ('time','longitude','latitude',))
    pri7 = root_grp.createVariable('95_Percentile_Gridded','f8', ('time','longitude','latitude',))
    pri8 = root_grp.createVariable('99_Percentile_Gridded','f8', ('time','longitude','latitude',))
    aux1 = root_grp.createVariable('Std_Dev_Time','f8', ('time','longitude','latitude',))
    aux2 = root_grp.createVariable('Std_Dev_Sites','f8', ('time','longitude','latitude',))
    aux3 = root_grp.createVariable('Count_Nsites','f8', ('time','longitude','latitude',))
    aux4 = root_grp.createVariable('DataFrac','f8', ('time','longitude','latitude',))
    aux5 = root_grp.createVariable('Mean_Altitude','f8', ('time','longitude','latitude',))
    
    #variable attributes
    t.units = 'seconds since 1970-01-01 00:00:00.0'
    t.standard_name = 'time'
    t.long_name = 'time'
    t.calendar = 'standard'
    t.axis = 'T'
    t.tz = 'UTC'
    t.valid_min = '0.0'
    
    lon.units = 'degrees_east'
    lon.standard_name = 'longitude'
    lon.long_name = 'longitude'
    lon.axis = 'X'
    lon.valid_min = np.min(lon_c)
    lon.valid_max = np.max(lon_c)
    
    lat.units = 'degrees_north'
    lat.standard_name = 'latitude'
    lat.long_name = 'latitude'
    lat.axis = 'Y'
    lat.valid_min = np.min(lat_c)
    lat.valid_max = np.max(lat_c)
    
    m1.standard_name = 'meta_fields'
    m1.long_name = 'metadata fields'
    
    m2.standard_name = 'meta_values'
    m2.long_name = 'metadata parameter values'
    
    mtd1.standard_name = 'methods'
    mtd1.long_name = 'unique measurement methodologies applied at each respective site'
    
    mtd2.units = 'seconds since 1970-01-01 00:00:00.0'
    mtd2.standard_name = 'methods_times'
    mtd2.long_name = 'range when unique measurement methodologies are applied at each respective site'
    
    pri1.units = 'ppbv'
    pri1.standard_name = 'Mean_Gridded'
    pri1.long_name = 'Mean of gridded hourly ozone mixing ratios (ppbv)'
    pri1.valid_min = '0.0'
    
    pri2.units = 'ppbv'
    pri2.standard_name = 'Median_Gridded'
    pri2.long_name = 'Median (50th Percentile) of gridded hourly ozone mixing ratios (ppbv)'
    pri2.valid_min = '0.0'
    
    pri3.units = ''
    pri3.standard_name = 'Skewness_Gridded'
    pri3.long_name = 'Skewness of gridded hourly ozone mixing ratios (ppbv)'
    
    pri4.units = ''
    pri4.standard_name = 'Kurtosis_Gridded'
    pri4.long_name = 'Kurtosis of gridded hourly ozone mixing ratios; zero centered (ppbv)'
    
    pri5.units = 'ppbv'
    pri5.standard_name = '25_Percentile_Gridded'
    pri5.long_name = '25th Percentile of gridded hourly ozone mixing ratios (ppbv)'
    pri5.valid_min = '0.0'
    
    pri6.units = 'ppbv'
    pri6.standard_name = '75_Percentile_Gridded'
    pri6.long_name = '75th Percentile of gridded hourly ozone mixing ratios (ppbv)'
    pri6.valid_min = '0.0'
    
    pri7.units = 'ppbv'
    pri7.standard_name = '95_Percentile_Gridded'
    pri7.long_name = '95th Percentile of gridded hourly ozone mixing ratios (ppbv)'
    pri7.valid_min = '0.0'
    
    pri8.units = 'ppbv'
    pri8.standard_name = '99_Percentile_Gridded'
    pri8.long_name = '99th Percentile of gridded hourly ozone mixing ratios (ppbv)'
    pri8.valid_min = '0.0'
    
    aux1.units = 'ppbv'
    aux1.standard_name = 'Std_Dev_Time'
    aux1.long_name = 'long_name: Standard Deviation of Gridded Data with respect to time (hours) within an averaging time interval (See eq. 3, Bowdalo et al., ESSD, 2017)'
    aux1.valid_min = '0.0'
    
    aux2.units = 'ppbv'
    aux2.standard_name = 'Std_Dev_Sites'
    aux2.long_name = 'Standard Deviation between sites going into a grid box averaged over hours within a time interval (See eq. 4, Bowdalo et al., ESSD, 2017)'
    aux2.valid_min = '0.0'
    
    aux3.units = ''
    aux3.standard_name = 'Count_Nsites'
    aux3.long_name = 'number of sites with data in a gridbox for each time interval'
    aux3.valid_min = '0.0'
    
    aux4.units = ''
    aux4.standard_name = 'DataFrac'
    aux4.long_name = 'Data completeness fraction - number of hours and sites with data divided by total number of site x hours'
    aux4.valid_min = '0.0'

    aux5.units = 'm'
    aux5.standard_name = 'Mean_Altitude'
    aux5.long_name = 'Mean site altitude (m)'
    aux5.valid_min = '0.0'
    
    #fill variables
    t[:] = t_out
    lon[:] = lon_c 
    lat[:] = lat_c
    m1[:] = meta_str
    m2[:] = meta
    pri1[:] = np.NaN
    pri2[:] = np.NaN
    pri3[:] = np.NaN
    pri4[:] = np.NaN
    pri5[:] = np.NaN
    pri6[:] = np.NaN
    pri7[:] = np.NaN
    pri8[:] = np.NaN
    aux1[:] = np.NaN
    aux2[:] = np.NaN
    aux3[:] = np.NaN
    aux4[:] = np.NaN
    aux5[:] = np.NaN
    
    #-------------------------------------------------------------------------
    #get unique methods by site
    
    unique_method_dict = {1:'ultraviolet photometry',2:'vacuum ultraviolet fluorescence',3:'chemiluminescence',
                          4:'chemiluminescence (conversion-molybdenum)',5:'chemiluminescence (conversion-photolysis)',
                          6:'ethylene chemiluminescence',7:'potassium iodide',8:'electrochemical concentration cell',
                          9:'non-dispersive infrared spectroscopy',10:'differential optical absorption spectrosocopy',
                          11:'cavity ringdown spectroscopy',12:'cavity attenuated phase shift spectroscopy',
                          13:'gas chromatography flame ionisation detection',14:'gas chromatography mercury oxide reduction',
                          15:'gas chromatography mass spectrometry',17:'gas chromatography photo ionization detection',
                          18:'gas chromatography fourier transform infrared spectroscopy',23:'photometric flame photometry',
                          24:'proton transfer reaction mass spectrometry',25:'ion chromatography',26:'spectrometry',
                          27:'spectrophotometry',28:'colorimetric',29:'fluid injection analysis',30:'coulometry',
                          31:'diffusive sampler'}
                          
    unique_methods = []
    for l in obs_methods:
        l = l[(l != -1) & (l != 0) & (l != 99)]
        unique_methods.append(list(set(l)))
    
    for s in range(len(refs)):
        site_methods = obs_methods[s]
        unique_method_str = []
        method_t = []
        for u_m in unique_methods[s]:
            inds = np.where(site_methods==u_m)[0]
            #get start and end time 0f each method (in seconds from 1980)
            start_t = int((inds[0]*60.)*60.)
            end_t = int((inds[-1]*60.)*60.)
            #then put into string
            method_t=np.append(method_t,'%s-%s'%(start_t,end_t))
            #convert unique method number to string
            unique_method_str=np.append(unique_method_str,unique_method_dict[u_m])
            
        mtd1[:,s] = unique_method_str
        mtd2[:,s] = method_t
    
    #-------------------------------------------------------------------------
    
    count_sites_found = 0
    # loop over all gridboxes, finding data in that gridbox
    #Note that lon and lat dimensions should be one element shorter than lonarr, latarr and shift the grid indices by 1/2 grid length
    import itertools
    lola_ind = itertools.product(xrange(len(lonarr)-1),xrange(len(latarr)-1))
    for i,j in lola_ind:
        inlonrange = np.logical_and(obs_lons>=longrid[i,j], obs_lons<longrid[i+1,j])
        inlatrange = np.logical_and(obs_lats>=latgrid[i,j], obs_lats<latgrid[i,j+1])
        in_gridbox = np.logical_and(inlonrange,inlatrange)
        
        if in_gridbox.sum() > 0:
            print in_gridbox.sum(), i, j, longrid[i,j],longrid[i+1,j],latgrid[i,j],latgrid[i,j+1]
            for t in range(len(datetime_out)):
                print 
                s_ind = int((((datetime_out[t] - datetime.datetime(1970,1,1,0,0)).total_seconds())/60.)/60.)
                try:
                    e_ind = int((((datetime_out[t+1] - datetime.datetime(1970,1,1,0,0)).total_seconds())/60.)/60.)
                except:
                    e_ind = int((((datetime.datetime(2015,1,1,0,0) - datetime.datetime(1970,1,1,0,0)).total_seconds())/60.)/60.)
                
                #get data in time window, in gridbox
                cut_data = obs_data[in_gridbox,s_ind:e_ind]
                
                site_inds = np.arange(len(refs))[in_gridbox]
                
                print len(cut_data[cut_data >= 0]), s_ind,e_ind
                
                #convert -99999's to NaNs
                cut_data[cut_data == -99999] = np.NaN
                
                #loop through sites and determine n sites with valid data
                valid_site_inds = []
                for s,ss in enumerate(site_inds):
                    if np.isnan(np.nanmean(cut_data[s,:])) != True:
                        valid_site_inds.append(ss)
                valid_site_inds = np.array(valid_site_inds)
                
                print np.array(refs)[site_inds]
                
                #if have valid data for time period then proceed
                if len(valid_site_inds) > 0:
                    
                    #primary stats
                    average = np.nanmean(cut_data)
                    median = np.nanmedian(cut_data)
                    kurtosis = stats.kurtosis(cut_data[~np.isnan(cut_data)])
                    skewness = stats.skew(cut_data[~np.isnan(cut_data)])
                    percentile25 = np.nanpercentile(cut_data, 25)
                    percentile75 = np.nanpercentile(cut_data, 75)
                    percentile95 = np.nanpercentile(cut_data, 95)
                    percentile99 = np.nanpercentile(cut_data, 99)
                    
                    #auxiliary stats
                    count_sites = len(valid_site_inds)
                    ave_alt = np.nanmean(obs_alt[valid_site_inds])
                    
                    merge_cut_data = np.nanmean(cut_data,axis=0)
                    datafrac = (100./len(merge_cut_data))*len(merge_cut_data[~np.isnan(merge_cut_data)])
                    
                    #calculate temporal stdev
                    #calculates the standard deviation of the gridded data over the time (hours) within each timestep
                    std_dev_t = np.nanstd(cut_data)
                    
                    #calculate spatial stdev - only calculate if more than 1 site
                    #calculates the average standard deviation of individual sites within a gridbox 
                    #first calculate the standard deviation for each hour across the sites in a gridbox   
                    if len(valid_site_inds) > 1:
                        std_dev_hours = []
                        for x in range(len(cut_data[0,:])):
                            std_dev_hours = np.append(std_dev_hours,np.nanstd(cut_data[:,x]))
                        #then take std dev over that time interval
                        std_dev_sites = np.nanmean(std_dev_hours)
                    else:
                       std_dev_sites = np.NaN 
                        
                    #write data to netcdf variables
                    pri1[t,i,j] = average
                    pri2[t,i,j] = median
                    pri3[t,i,j] = skewness
                    pri4[t,i,j] = kurtosis
                    pri5[t,i,j] = percentile25
                    pri6[t,i,j] = percentile75
                    pri7[t,i,j] = percentile95
                    pri8[t,i,j] = percentile99
                    
                    aux1[t,i,j] = std_dev_t
                    aux2[t,i,j] = std_dev_sites
                    aux3[t,i,j] = count_sites
                    aux4[t,i,j] = datafrac
                    aux5[t,i,j] = ave_alt
                    
                    #print cut_data
                    
                    print average,median,skewness,kurtosis,percentile25,percentile75,percentile95,percentile99
                    print std_dev_t,std_dev_sites,count_sites,datafrac,ave_alt      
        
        count_sites_found+=in_gridbox.sum() # to keep track of whether all sites are used after doing global gridding
    #print count_sites_found
    
    root_grp.close()


#1.10
#----------
#SET DEFINED AREAS FOR AVERAGE SITE DATA
def def_areas():
    areas = ['ANT','S_O','OC','SA','AF','SW_NA','NW_NA','NE_NA','CE_NA','SE_NA','C_NA','S_EU','SW_EU','CW_EU','CS_EU','C_EU','E_EU','NW_EU','N_EU','S_AS','SE_AS','C_AS','NE_AS','N_O','AL','ARC']
    return areas

def model_areas():
    areas = ['ANT','CE_NA','SE_NA','S_NA','SW_NA','C_NA','NE_NA','NW_NA','AL','N_EU','C_EU','SW_EU','S_EU','E_EU','NW_EU','NE_AS','SW_AS','NW_AS','SE_AS','NE_AS','JAP','N_O','S_O','ARC']
    return areas    

#1.11
#----------
#SET DICTIONARIES FOR SPLITTING SITES INTO DIFFERENT AREAS

def area_dicts():  
    area_boundaries = {'CE_NA':[37,43.5,-90,-50],'SE_NA':[25,37,-90,-60],'S_NA':[25,41,-108,-90],'SW_NA':[25,40,-130,-108],'C_NA':[41,70,-108,-90],'NE_NA':[43.5,80,-90,-50],'NW_NA':[40,80,-141,-108],'AL':[50,80,-170,-141],'N_EU':[55,80,4,45],'C_EU':[46,55,4,12.5],'SW_EU':[35,46,-10,0],'S_EU':[30,46,0,40],'E_EU':[46,55,12.5,40],'NW_EU':[46,70,-15,4],'NE_AS':[30,60,90,155],'SE_AS':[-10,30,90,130],'C_AS':[30,60,60,90],'S_AS':[0,30,60,90],'N_O':[0,0,0,0],'S_O':[0,0,0,0],'OC':[0,0,0,0],'AF':[0,0,0,0],'SA':[0,0,0,0],'ANT':[0,0,0,0],'ARC':[0,0,0,0]}
    area_tag = {'CE_NA':['NA'],'SE_NA':['NA'],'S_NA':['NA'],'SW_NA':['NA'],'C_NA':['NA'],'NE_NA':['NA'],'NW_NA':['NA'],'AL':['NA'],'C_EU':['EU'],'N_EU':['EU'],'S_EU':['EU'],'SW_EU':['EU'],'E_EU':['EU'],'NW_EU':['EU'],'NE_AS':['EA'],'SE_AS':['EA'],'C_AS':['EA'],'S_AS':['EA'],'N_O':['O'],'S_O':['O'],'OC':['OC'],'AF':['AF'],'SA':['SA'],'ANT':['ANT'],'ARC':['ARC']}
    area_labels = {'CE_NA':'CE NA','SE_NA':'SE NA','S_NA':'S NA','SW_NA':'SW NA','C_NA':'C NA','NE_NA':'NE NA','NW_NA':'NW NA','AL':'Alaska','N_EU':'N EU','NW_EU':'NW EU','C_EU':'C EU','E_EU':'E EU','S_EU':'S EU','SW_EU':'SW EU','NE_AS':'NE Asia','SE_AS':'SE Asia','C_AS':'C Asia','S_AS':'S Asia','N_O':'NH Oceanic','S_O':'SH Oceanic','OC':'Oceania','AF':'Africa','SA':'S America','ANT':'Antarctica','ARC':'Arctic'}
    return area_boundaries,area_tag,area_labels

def model_area_dicts():
    area_boundaries = {'ANT':[-90,-60,-180,180],'CE_NA':[37,45,-90,-50],'SE_NA':[25,37,-90,-60],'S_NA':[25,37,-105,-90],'SW_NA':[25,40,-130,-105],'C_NA':[37,70,-115,-90],'NE_NA':[45,70,-90,-50],'NW_NA':[40,70,-141,-115],'AL':[50,70,-170,-141],'N_EU':[55,70,5,42],'C_EU':[44,55,4,16],'SW_EU':[35,46,-15,4],'S_EU':[30,44,4,37],'E_EU':[44,55,16,42],'NW_EU':[46,60,-15,4],'NE_AS':[30,60,90,155],'SW_AS':[-10,32,65,90],'NW_AS':[32,50,65,100],'SE_AS':[-10,32,90,123],'NE_AS':[32,50,100,123],'JAP':[30,50,123,150],'N_O':[0,90,-180,180],'S_O':[-90,0,-180,180],'ARC':[70,90,-180,180]}
    area_labels = {'ANT':'ANT','CE_NA':'CE NA','SE_NA':'SE NA','S_NA':'S NA','SW_NA':'SW NA','C_NA':'C NA','NE_NA':'NE NA','NW_NA':'NW NA','AL':'AL','N_EU':'N EU','C_EU':'C EU','SW_EU':'SW EU','S_EU':'S EU','E_EU':'E EU','NW_EU':'NW EU','NE_AS':'NE AS','SW_AS':'SW AS','NW_AS':'NW AS','SE_AS':'SE AS','NE_AS':'NE AS','JAP':'JAP','N_O':'N O','S_O':'S O','ARC':'ARC'}
    return area_boundaries,area_labels

#1.12
#----------
#GET AREA CUT FOR SITES USING LAT,LON AND TAG
 
def area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag):
    if (area == 'C_EU') or (area == 'E_EU') or (area == 'S_EU') or (area == 'NW_EU') or (area == 'N_EU') or (area == 'CW_EU') or (area == 'SW_EU') or (area == 'CS_EU'):
        cut_test = (obs_lats >= area_grid[0]) & (obs_lats < area_grid[1]) & (obs_lons >= area_grid[2]) & (obs_lons < area_grid[3]) & (tags == area_tag[0])
    elif (area == 'OC') or (area == 'AF') or (area == 'ANT') or (area == 'ARC') or (area == 'SA'):
        cut_test = tags == area_tag[0]
    elif (area == 'N_O'):
        cut_test = (tags == area_tag[0]) & (obs_lats >= 0)
    elif (area == 'S_O'):
        cut_test = (tags == area_tag[0]) & (obs_lats < 0)
    elif (area == 'CE_NA') or (area == 'SE_NA') or (area == 'S_NA') or (area == 'SW_NA') or (area == 'C_NA') or (area == 'NE_NA') or (area == 'NW_NA') or (area == 'AL'):              
        cut_test = (obs_lats >= area_grid[0]) & (obs_lats < area_grid[1]) & (obs_lons >= area_grid[2]) & (obs_lons < area_grid[3]) & (tags == area_tag[0])
    elif (area == 'NE_AS') or (area == 'SE_AS') or (area == 'S_AS') or (area == 'C_AS'):
        cut_test = (obs_lats >= area_grid[0]) & (obs_lats < area_grid[1]) & (obs_lons >= area_grid[2]) & (obs_lons < area_grid[3]) & (tags == area_tag[0])

    return cut_test
   
#1.13
#----------
#GET AREA FOR SITES USING LAT,LON AND TAG
 
def get_area(area_boundaries,tag,obs_lat,obs_lon,obs_ref):
    
    if tag == 'NA':
        if (obs_lat >= area_boundaries['CE_NA'][0]) & (obs_lat < area_boundaries['CE_NA'][1]) & (obs_lon >= area_boundaries['CE_NA'][2]) & (obs_lon < area_boundaries['CE_NA'][3]):
            area = 'CE_NA'
        elif (obs_lat >= area_boundaries['SE_NA'][0]) & (obs_lat < area_boundaries['SE_NA'][1]) & (obs_lon >= area_boundaries['SE_NA'][2]) & (obs_lon < area_boundaries['SE_NA'][3]):
            area = 'SE_NA'
        elif (obs_lat >= area_boundaries['S_NA'][0]) & (obs_lat < area_boundaries['S_NA'][1]) & (obs_lon >= area_boundaries['S_NA'][2]) & (obs_lon < area_boundaries['S_NA'][3]):
            area = 'S_NA'
        elif (obs_lat >= area_boundaries['SW_NA'][0]) & (obs_lat < area_boundaries['SW_NA'][1]) & (obs_lon >= area_boundaries['SW_NA'][2]) & (obs_lon < area_boundaries['SW_NA'][3]):
            area = 'SW_NA'
        elif (obs_lat >= area_boundaries['C_NA'][0]) & (obs_lat < area_boundaries['C_NA'][1]) & (obs_lon >= area_boundaries['C_NA'][2]) & (obs_lon < area_boundaries['C_NA'][3]):
            area = 'C_NA'
        elif (obs_lat >= area_boundaries['NE_NA'][0]) & (obs_lat < area_boundaries['NE_NA'][1]) & (obs_lon >= area_boundaries['NE_NA'][2]) & (obs_lon < area_boundaries['NE_NA'][3]):
            area = 'NE_NA'
        elif (obs_lat >= area_boundaries['NW_NA'][0]) & (obs_lat < area_boundaries['NW_NA'][1]) & (obs_lon >= area_boundaries['NW_NA'][2]) & (obs_lon < area_boundaries['NW_NA'][3]):
            area = 'NW_NA'
        elif (obs_lat >= area_boundaries['AL'][0]) & (obs_lat < area_boundaries['AL'][1]) & (obs_lon >= area_boundaries['AL'][2]) & (obs_lon < area_boundaries['AL'][3]):
            area = 'AL'
        else:
            area = 'NO GIVEN AREA'
            print tag,obs_lat,obs_lon,obs_ref
            
    elif tag == 'EU':
        if (obs_lat >= area_boundaries['N_EU'][0]) & (obs_lat < area_boundaries['N_EU'][1]) & (obs_lon >= area_boundaries['N_EU'][2]) & (obs_lon < area_boundaries['N_EU'][3]):
            area = 'N_EU'
        elif (obs_lat >= area_boundaries['C_EU'][0]) & (obs_lat < area_boundaries['C_EU'][1]) & (obs_lon >= area_boundaries['C_EU'][2]) & (obs_lon < area_boundaries['C_EU'][3]):
            area = 'C_EU'
        elif (obs_lat >= area_boundaries['CS_EU'][0]) & (obs_lat < area_boundaries['CS_EU'][1]) & (obs_lon >= area_boundaries['CS_EU'][2]) & (obs_lon < area_boundaries['CS_EU'][3]):
            area = 'CS_EU'  
        elif (obs_lat >= area_boundaries['CW_EU'][0]) & (obs_lat < area_boundaries['CW_EU'][1]) & (obs_lon >= area_boundaries['CW_EU'][2]) & (obs_lon < area_boundaries['CW_EU'][3]):
            area = 'CW_EU'
        elif (obs_lat >= area_boundaries['SW_EU'][0]) & (obs_lat < area_boundaries['SW_EU'][1]) & (obs_lon >= area_boundaries['SW_EU'][2]) & (obs_lon < area_boundaries['SW_EU'][3]):
            area = 'SW_EU'
        elif (obs_lat >= area_boundaries['S_EU'][0]) & (obs_lat < area_boundaries['S_EU'][1]) & (obs_lon >= area_boundaries['S_EU'][2]) & (obs_lon < area_boundaries['S_EU'][3]):
            area = 'S_EU'
        elif (obs_lat >= area_boundaries['NW_EU'][0]) & (obs_lat < area_boundaries['NW_EU'][1]) & (obs_lon >= area_boundaries['NW_EU'][2]) & (obs_lon < area_boundaries['NW_EU'][3]):
            area = 'NW_EU'
        elif (obs_lat >= area_boundaries['E_EU'][0]) & (obs_lat < area_boundaries['E_EU'][1]) & (obs_lon >= area_boundaries['E_EU'][2]) & (obs_lon < area_boundaries['E_EU'][3]):
            area = 'E_EU'
        else:
            area = 'NO GIVEN AREA'
            print tag,obs_lat,obs_lon,obs_ref
     
    elif (tag == 'MA') or (tag == 'EA'):
        if (obs_lat >= area_boundaries['NE_AS'][0]) & (obs_lat < area_boundaries['NE_AS'][1]) & (obs_lon >= area_boundaries['NE_AS'][2]) & (obs_lon < area_boundaries['NE_AS'][3]):
            area = 'NE_AS'
        elif (obs_lat >= area_boundaries['SE_AS'][0]) & (obs_lat < area_boundaries['SE_AS'][1]) & (obs_lon >= area_boundaries['SE_AS'][2]) & (obs_lon < area_boundaries['SE_AS'][3]):
            area = 'SE_AS'    
        elif (obs_lat >= area_boundaries['C_AS'][0]) & (obs_lat < area_boundaries['C_AS'][1]) & (obs_lon >= area_boundaries['C_AS'][2]) & (obs_lon < area_boundaries['C_AS'][3]):
            area = 'C_AS'
        elif (obs_lat >= area_boundaries['S_AS'][0]) & (obs_lat < area_boundaries['S_AS'][1]) & (obs_lon >= area_boundaries['S_AS'][2]) & (obs_lon < area_boundaries['S_AS'][3]):
            area = 'S_AS'
        else:
            area = 'NO GIVEN AREA'
            print tag,obs_lat,obs_lon,obs_ref
        
    elif tag == 'OC':
        area = 'OC'
        
    elif tag == 'ANT':
        area = 'ANT'
      
    elif tag == 'ARC':
        area = 'ARC'  
    
    elif tag == 'AF':
        area = 'AF'
    
    elif tag == 'SA':
        area = 'SA'
    
    elif (tag == 'O') & (obs_lat >= 0):
        area = 'N_O'
        
    elif (tag == 'O') & (obs_lat < 0):
        area = 'S_O'
                
    return area

#1.14
#----------
#GET COUNTRY FROM LAT,LON USING GOOGLE API

def get_country(lat, lon):
    data = json.load(urllib2.urlopen('http://maps.googleapis.com/maps/api/geocode/json?latlng=%s,%s&sensor=false' % (lat, lon)))
    if len(data['results']) >= 1:
        for result in data['results']:
            for component in result['address_components']:
                if 'country' in component['types']:
                    return component['long_name']        
    else:
        return 'Ocean'

#1.15
#----------
#GET ANTHROME CLASSIFICATION FROM LAT,LON AND ANTHROME MAP 

def anthrome_classify(anthload,obs_lat,obs_lon):
    
    anthromes = {
    11: 'Urban',
    12: 'Mixed settlements',
    21: 'Rice villages',
    22: 'Irrigated villages',
    23: 'Rainfed villages',
    24: 'Pastoral villages',
    31: 'Residential irrigated croplands',
    32: 'Residential rainfed croplands',
    33: 'Populated croplands',
    34: 'Remote croplands',
    41: 'Residential rangelands',
    42: 'Populated rangelands',
    43: 'Remote rangelands',
    51: 'Residential woodlands',
    52: 'Populated woodlands',
    53: 'Remote woodlands',
    54: 'Inhabited treeless and barren lands',
    61: 'Wild woodlands',
    62: 'Wild treeless and barren lands'
    }

    #Anthrome classes may be grouped for analysis into Anthrome levels:
    simple_anthromes = {
    11:	'Dense Settlements',
    12:	'Dense Settlements',
    21:	'Villages',
    22:	'Villages',
    23:	'Villages',
    24:	'Villages',
    31:	'Croplands',
    32:	'Croplands',
    33:	'Croplands',
    34:	'Croplands',
    41:	'Rangelands',
    42:	'Rangelands',
    43:	'Rangelands',
    51:	'Seminatural',
    52:	'Seminatural',
    53:	'Seminatural',
    54:	'Seminatural',
    61:	'Wildlands',
    62:	'Wildlands',
    }

    simple_groups = {
    1: 'Dense Settlement',
    2: 'Village',
    3: 'Cropland',
    4: 'Rangeland',
    5: 'Seminaturnal',
    6: 'Wildland',
    }
    
    anth_lon_c = anthload.variables['lon'][:]
    anth_lat_c = anthload.variables['lat'][:]
    anthdata = anthload.variables['anthro2_a2000'][:]
    
    anth_lat_c = anth_lat_c[::-1]
    
    anthdata = anthdata[::-1,:]
    anthdata = anthdata//10
    
    lon_diff = (anth_lon_c[1]-anth_lon_c[0])/2.
    lat_diff = (anth_lat_c[1] - anth_lat_c[0])/2.
    
    anth_lon_e = anth_lon_c-lon_diff
    anth_lon_e = np.append(anth_lon_e,180.00)
    anth_lat_e = anth_lat_c-lat_diff
    anth_lat_e = np.append(anth_lat_e,90.00)

    first_lon_i = 0
    first_lat_i = 0
    last_lon_i = len(anth_lon_c)-1
    last_lat_i = len(anth_lat_c)-1

    if len(obs_lat) > 1:
        validity = []
        class_name = []

        #check anthdata grid around gridbox site is in, if more than 50% of 9 gridboxes 
        #are dense settlement then site is invalid

        for o_lat,o_lon in zip(obs_lat,obs_lon): 
            anthdata_grids = []
            class_ns = []
            
            lat_i,lon_i = obs_model_gridbox(anth_lat_e,anth_lon_e,o_lat,o_lon) 
            
            lat_is = range(lat_i-1,lat_i+2)
            lon_is = range(lon_i-1,lon_i+2)
            
            #take care of wrapping if neccessary
            if lat_is[0] < first_lat_i:
                lat_is[0] = last_lat_i
            if lat_is[2] > last_lat_i:
                lat_is[2] = first_lat_i
            if lon_is[0] < first_lon_i:
                lon_is[0] = last_lon_i
            if lon_is[2] > last_lon_i:
                lon_is[2] = first_lon_i      
            
            for la in range(len(lat_is)):
                for lo in range(len(lon_is)):
                    anthdata_grid = anthdata[lat_is[la],lon_is[lo]]  
                    try:
                        class_ns=np.append(class_ns,simple_groups[anthdata_grid])
                    except:
                        anthdata_grid = 0
                        class_ns = np.append(class_ns,'Ocean')
                    anthdata_grids=np.append(anthdata_grids,anthdata_grid)
    
            test = anthdata_grids == 1
            n_urban = len(anthdata_grids[test])
            if n_urban >= 5:
                valid = 'invalid'
                class_n = 'Dense Settlement'
            else:
                valid = 'valid'
                test = class_ns != 'Dense Settlement'
                class_n = stats.mode(class_ns[test])[0][0]
                
            validity.append(valid)
            class_name.append(class_n)
        validity = np.array(validity)
        class_name = np.array(class_name)
    
    
    else:
        obs_lon = obs_lon[0]
        obs_lat = obs_lat[0]
        
        anthdata_grids = []
        class_ns = []
        
        lat_i,lon_i = obs_model_gridbox(anth_lat_e,anth_lon_e,obs_lat,obs_lon)

        lat_is = range(lat_i-1,lat_i+2)
        lon_is = range(lon_i-1,lon_i+2)
        
        #take care of wrapping if neccessary
        if lat_is[0] < first_lat_i:
            lat_is[0] = last_lat_i
        if lat_is[2] > last_lat_i:
            lat_is[2] = first_lat_i
        if lon_is[0] < first_lon_i:
            lon_is[0] = last_lon_i
        if lon_is[2] > last_lon_i:
            lon_is[2] = first_lon_i  

        for la in range(len(lat_is)):
            for lo in range(len(lon_is)):
                anthdata_grid = anthdata[lat_is[la],lon_is[lo]]  
                try:
                    class_ns=np.append(class_ns,simple_groups[anthdata_grid])
                except:
                    anthdata_grid = 0
                    class_ns = np.append(class_ns,'Ocean')
                anthdata_grids=np.append(anthdata_grids,anthdata_grid)
        
        test = anthdata_grids == 1
        n_urban = len(anthdata_grids[test])
        if n_urban >= 5:
            validity = 'invalid'
            class_name = 'Dense Settlement'
        else:
            validity = 'valid'
            test = class_ns != 'Dense Settlement'
            class_name = stats.mode(class_ns[test])[0][0]
            #print 'Site is %s'%(class_name)

        #print obs_lon,obs_lat
        #print anthdata_grid

        #fig, ax = plt.subplots(1,1)
        #m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
        #                llcrnrlon=-180,\
        #                 urcrnrlon=180,\
        #                resolution='h')
                        
        #m.drawcoastlines(linewidth=0.5)
        #m.drawcountries()
        #x,y=m(anth_lon_c,anth_lat_c)
        #cmap = mpl.cm.get_cmap('rainbow', 7)    # 7 discrete colors
        #norm = mpl.colors.BoundaryNorm(np.arange(8),7)
        #plt.pcolormesh(x,y,anthdata, cmap=cmap, norm=norm)
        #plt.colorbar()
        #plt.show()
    
    return validity,class_name

#1.16
#----------

#get continental tags from obs ref (manually compiled)

def get_tags(obs_refs):
                #GAW
    tag_dict = {'arh':'ANT','ask':'AF','bhd':'OC','bmw':'O','brt':'NA','brw':'NA','cgo':'OC','cmn':'EU','cpt':'AF','cvo':'O','dcc':'ANT','deu':'EU','dig':'EU','glh':'EU','hpb':'EU','irb':'EU','izo':'O','jfj':'EU', \
                'kmw':'EU','kos':'EU','kps':'EU','kvk':'EU','kvv':'EU','lau':'OC','lgb':'EU','mhd':'EU','mlo':'O','mnm':'O','ngl':'EU','nmy':'ANT','nwr':'NA','pal':'EU','pay':'EU','pyr':'MA','rcv':'EU', \
                'rig':'EU','rpb':'O','ryo':'EA','smo':'O','snb':'EU','snl':'SA','spo':'ANT','ssl':'EU','sum':'ARC','syo':'ANT','thd':'NA','tkb':'EA','ush':'SA','vdl':'EU','bkt':'OC', \
                'wes':'EU','yon':'EA','zgt':'EU','zpt':'EU','zrn':'EU','est':'NA','egb':'NA','ela':'NA','sat':'NA','kej':'NA','ang':'O','cas':'EU','prs':'EU', \
                'dak':'EU','dbl':'EU','dmv':'OC','lon':'NA','mbi':'ANT','mcm':'ANT','mvh':'EU','shp':'EU','sja':'SA','wkt':'NA','fra':'NA','mtm':'NA','beo':'EU','mkn':'AF', \
                'ams':'O','asc':'O','azr':'O','bal':'EU','bme':'O','bsc':'EU','cba':'NA','cfa':'OC','chr':'O','cmo':'NA','cri':'MA','crz':'O','cya':'ANT','eic':'O','gmi':'O','goz':'EU','hba':'ANT','hun':'EU','itn':'NA','kum':'O','kzd':'MA','lef':'NA','lmp':'EU','maa':'ANT', \
                'mbc':'ARC','mid':'O','mqa':'O','nmb':'AF','poc900n':'O','poc905n':'O','poc905s':'O','poc910n':'O','poc910s':'O','poc915n':'EU','poc915s':'EU','poc920n':'O','poc920s':'O','poc925n':'O','poc925s':'O', \
                'poc930n':'O','poc930s':'O','poc935n':'O','poc935s':'O','psa':'ANT','pta':'NA','scs903n':'O','scs906n':'O','scs909n':'O','scs912n':'O','scs915n':'O','scs918n':'O','scs921n':'O','sey':'O','sgp':'NA', \
                'shm':'O','sis':'O','stm':'O','tap':'EA','tdf':'SA','uum':'MA','wis':'EU','cha':'NA','ice':'EU','alt':'ARC','edm':'EU','zep':'EU', 'cdl':'NA','chm':'NA','esp':'NA','llb':'NA','bur':'EU','ivn':'EU','jcz':'EU','kam':'EU','leb':'EU','log':'EU','plm':'SA','plv':'EU','zsn':'EU', \
                'aht':'EU','oul':'EU','uto':'EU','vir':'EU','ktb':'EU','mhn':'EU','nia':'EU','roq':'EU','spm':'EU','puy':'MA','lis':'EU','isk':'MA','fun':'EU','hkg':'MA','tll':'SA','zug':'EU','lqo':'SA','tar':'OC','zsf':'EU','cai':'AF','don':'EU','fdt':'EU','pdi':'MA','pen':'EU','pil':'SA','wlg':'MA', \
                #CAPMON
                'capmonalg':'NA','capmonbra':'NA','capmonalt':'NA','capmongos':'NA','capmonmin':'NA','capmonpkl':'NA', \
                #CHILE
                'ch004':'SA','ch006':'SA','ch007':'SA', \
                #CASTNET
                'qak172':'NA','eve419':'NA','vin140':'NA','irl141':'NA','oxf122':'NA','prk134':'NA','lav410':'NA','dcp114':'NA','par107':'NA','thr422':'NA','cat175':'NA','gth161':'NA','bbe401':'NA','bel116':'NA','mac426':'NA','ckt136':'NA','bvl130':'NA','knz184':'NA', \
                'pin414':'NA','wsp144':'NA','ash135':'NA','cad150':'NA','cha467':'NA','hwf187':'NA','cdr119':'NA','voy413':'NA','cvl151':'NA','gas153':'NA','wst109':'NA','cnt169':'NA','sum156':'NA','egb181':'NA', \
                'bwr139':'NA','lrl117':'NA','mev405':'NA','grb411':'NA','cth110':'NA','yel408':'NA','cow137':'NA','pet427':'NA','pnd165':'NA','sal133':'NA','hox148':'NA','vpi120':'NA','can407':'NA','spd111':'NA', \
                'uvl124':'NA','mkg113':'NA','ped108':'NA','mor409':'NA','alc188':'NA','che185':'NA','kef112':'NA','cnd125':'NA','san189':'NA','psu106':'NA','aca416':'NA','grs420':'NA','shn418':'NA','glr468':'NA', \
                'rom206':'NA','mck131':'NA','cdz171':'NA','mck231':'NA','alh157':'NA','yos404':'NA','bft142':'NA','ana115':'NA','are128':'NA','jot403':'NA','esp127':'NA','abt147':'NA','snd152':'NA','sek430':'NA', \
                'rom406':'NA','stk138':'NA','den417':'NA','how132':'NA','pnf126':'NA','wnc429':'NA','grc474':'NA','con186':'NA','wel149':'NA','lyk123':'NA','are228':'NA','dcp214':'NA','lcw121':'NA','onl102':'NA', \
                'pbf129':'NA','sum256':'NA','wfm105':'NA','wpa103':'NA','wpb104':'NA','ash235':'NA','dev412':'NA','ncs415':'NA','vii423':'O','pal190':'NA','anl146':'NA','rtp101':'NA','rck163':'NA','sav164':'NA','uin162':'NA', \
                'sek402':'NA','oly421':'NA','hvt424':'NA','bas601':'NA','nec602':'NA','how191':'NA','lye145':'NA','she604':'NA','yos204':'NA',
                #EANET
                'eahe275':'EA','eaij231':'EA','eaog198':'EA','eaok242':'EA','eayu264':'EA','eari165':'EA','eaoc176':'EA','easa209':'EA','eaba253':'EA','eaha220':'EA','eata187':'EA','eaba495':'EA','each429':'EA','each539':'EA','eaim440':'EA','eaka418':'EA','eakh528':'EA','eamo451':'EA','easa506':'EA','eata308':'EA', \
                #EMEP
                'gb0014r':'EU','no0043r':'EU','no0056r':'EU','no0042g':'ARC','gb0037r':'EU', \
                'gb0002r':'EU','fi0009r':'EU','gb0006r':'EU','fi0096g':'EU','es0008r':'EU','gb0052r':'EU','gb0038r':'EU','ie0031r':'EU','gb0039r':'EU','gb0035r':'EU','gb0033r':'EU','gb0031r':'EU','fi0022r':'EU', \
                'gb0043r':'EU','gb0036r':'EU','no0039r':'EU','fi0037r':'EU','ie0001r':'EU','gr0001r':'EU','gb0049r':'EU','gb0051r':'EU','gb0013r':'EU','gb0015r':'EU','no0015r':'EU','gb0050r':'EU','fi0017r':'EU', \
                'es0012r':'EU','gb0048r':'EU','no0052r':'EU','gb0045r':'EU','min':'NA','fre':'NA','pkl':'NA','sna':'NA','gos':'NA','cps':'NA','rtr':'NA','bon':'NA','etl':'NA','fsd':'NA','wsa':'NA', \
                'cz0005r':'EU','fr0008r':'EU','fr0009r':'EU','fr0010r':'EU','fr0013r':'EU','fr0014r':'EU','fr0015r':'EU','fr0016r':'EU','fr0017r':'EU','mk0007r':'EU','nl0009r':'EU','ch0001g':'EU','cz0003r':'EU','nl0091r':'EU', \
                'at0004r':'EU','ca0420g':'ARC','de0004r':'EU','de0011r':'EU','de0012r':'EU','de0017r':'EU','de0026r':'EU','de0035r':'EU','es0001r':'EU','es0003r':'EU',  \
                'es0004r':'EU','es0005r':'EU','fi0004r':'EU','it0004r':'EU','no0030r':'EU','no0041r':'EU','no0044r':'EU','no0045r':'EU','no0048r':'EU','no0488r':'EU','se0002r':'EU', \
                'at0002r':'EU','at0005r':'EU','at0030r':'EU','at0032r':'EU','at0034g':'EU','at0037r':'EU','at0038r':'EU','at0040r':'EU','at0049r':'EU','fr0019r':'EU', \
                'at0041r':'EU','at0042r':'EU','at0043r':'EU','at0044r':'EU','at0045r':'EU','at0046r':'EU','at0047r':'EU','at0048r':'EU','be0001r':'EU','be0032r':'EU','be0035r':'EU','bg0053r':'EU', \
                'ch0002r':'EU','ch0003r':'EU','ch0004r':'EU','cy0002r':'EU','cz0001r':'EU','de0001r':'EU','de0002r':'EU','de0003r':'EU','de0007r':'EU','de0008r':'EU','de0009r':'EU','dk0005r':'EU', \
                'dk0031r':'EU','dk0041r':'EU','ee0009r':'EU','ee0011r':'EU','es0007r':'EU','es0009r':'EU','es0010r':'EU','es0011r':'EU','es0013r':'EU','es0014r':'EU','es0016r':'EU','hu0002r':'EU', \
                'it0001r':'EU','lt0015r':'EU','lv0010r':'EU','no0001r':'EU','no0055r':'EU','pl0002r':'EU','pl0003r':'EU','pl0004r':'EU','pl0005r':'EU','pt0004r':'EU','se0005r':'EU','se0011r':'EU', \
                'se0012r':'EU','se0013r':'EU','se0014r':'EU','se0032r':'EU','se0035r':'EU','se0039r':'EU','si0008r':'EU','si0031r':'EU','si0032r':'EU','sk0002r':'EU','sk0004r':'EU','sk0006r':'EU', \
                'sk0007r':'EU','at0003r':'EU','ch0031r':'EU','de0006r':'EU','de0013r':'EU','de0018r':'EU','de0038r':'EU','de0042r':'EU','de0046r':'EU', \
                'dk0010g':'ARC','fr0011r':'EU','fr0018r':'EU','fr0032r':'EU','gb0041r':'EU','gr0002r':'EU','no0762r':'EU','ru0018r':'EU','se0003r':'EU','se0033r':'EU','se0034r':'EU','se0094r':'EU', \
                'si0033r':'EU','nl0010r':'EU','de0031r':'EU','de0039r':'EU','gb0044r':'EU','mt0001r':'EU','de0045r':'EU','de0047r':'EU','nl0007r':'EU', \
                'be0011r':'EU','be0013r':'EU','de0043g':'EU','de0044r':'EU','dk0008r':'EU','dk0009r':'EU','dk0012r':'EU','es0006r':'EU','es0017r':'EU','gb0053r':'EU','lv0016r':'EU','md0013r':'EU','nl0002r':'EU','nl0011r':'EU','no0002r':'EU','no0008r':'EU', \
                'no0047r':'EU','no0057r':'EU','no0091r':'EU','no0099r':'EU','no0796r':'EU','no0797r':'EU','no1010r':'EU','ro0008r':'EU','rs0005r':'EU','ru0001r':'EU','ru0013r':'EU','ru0014r':'EU','ru0016r':'EU','se0008r':'EU','sk0005r':'EU', \
                'ch0005r':'EU','fr0012r':'EU','fr0030r':'EU','no0058g':'EU','no0489r':'EU','de0005r':'EU','de0014r':'EU','dk0032r':'EU','no0492r':'EU','am0001r':'EU','at0050r':'EU','cz0007r':'EU','es0002r':'EU','fr0031r':'EU','gr0003r':'EU',
                'nl0644r':'EU','no0059g':'EU','no0977r':'EU','no1007r':'EU','no1083r':'EU','no1200r':'EU','ro0003r':'EU', \
                #NAPS
                'canaps054703':'NA','canaps090602':'NA','canaps093801':'NA','canaps093901':'NA','canaps104803':'NA','canaps060706':'NA','canaps030118':'NA','canaps101002':'NA',\
                'canaps010601':'NA','canaps030501':'NA','canaps031001':'NA','canaps040203':'NA','canaps040207':'NA','canaps040302':'NA','canaps040401':'NA','canaps040501':'NA','canaps040601':'NA','canaps040801':'NA','canaps040901':'NA','canaps041201':'NA','canaps041302':'NA','canaps050103':'NA', \
                'canaps050113':'NA','canaps050116':'NA','canaps050119':'NA','canaps050121':'NA','canaps050126':'NA','canaps050129':'NA','canaps050204':'NA','canaps050308':'NA','canaps050310':'NA','canaps050311':'NA','canaps050404':'NA','canaps050504':'NA','canaps050604':'NA','canaps051501':'NA', \
                'canaps052001':'NA','canaps052201':'NA','canaps052301':'NA','canaps052401':'NA','canaps052601':'NA','canaps052801':'NA','canaps053201':'NA','canaps053301':'NA','canaps053401':'NA','canaps053501':'NA','canaps053601':'NA','canaps053701':'NA','canaps053901':'NA','canaps054102':'NA', \
                'canaps054201':'NA','canaps054801':'NA','canaps054901':'NA','canaps055001':'NA','canaps055101':'NA','canaps055301':'NA','canaps055501':'NA','canaps055601':'NA','canaps055701':'NA','canaps060211':'NA','canaps060410':'NA','canaps060428':'NA','canaps060429':'NA','canaps060513':'NA', \
                'canaps060809':'NA','canaps061004':'NA','canaps061104':'NA','canaps061201':'NA','canaps061603':'NA','canaps061802':'NA','canaps062001':'NA','canaps062501':'NA','canaps062601':'NA','canaps063001':'NA','canaps063301':'NA','canaps063701':'NA','canaps064001':'NA','canaps064101':'NA', \
                'canaps064401':'NA','canaps065001':'NA','canaps065101':'NA','canaps065201':'NA','canaps065401':'NA','canaps065901':'NA','canaps066101':'NA','canaps070118':'NA','canaps070203':'NA','canaps080901':'NA','canaps090120':'NA','canaps090222':'NA','canaps090502':'NA','canaps090701':'NA', \
                'canaps090702':'NA','canaps090801':'NA','canaps091001':'NA','canaps091301':'NA','canaps091401':'NA','canaps091501':'NA','canaps091801':'NA','canaps091901':'NA','canaps092001':'NA','canaps092601':'NA','canaps100110':'NA','canaps100118':'NA','canaps100119':'NA','canaps100125':'NA', \
                'canaps100127':'NA','canaps100128':'NA','canaps100132':'NA','canaps100134':'NA','canaps100135':'NA','canaps101003':'NA','canaps101101':'NA','canaps101202':'NA','canaps101301':'NA','canaps101401':'NA','canaps101501':'NA','canaps101701':'NA','canaps102001':'NA','canaps102102':'NA', \
                'canaps102401':'NA','canaps102801':'NA','canaps103302':'NA','canaps104003':'NA','canaps129102':'NA','canaps040103':'NA','canaps040701':'NA','canaps060104':'NA','canaps060204':'NA','canaps100111':'NA','canaps100202':'NA','canaps041101':'NA','canaps052701':'NA','canaps065301':'NA','canaps065701':'NA', \
                'canaps091601':'NA','canaps093101':'NA','canaps100314':'NA','canaps054401':'NA','canaps051901':'NA','canaps052101':'NA','canaps052901':'NA','canaps053001':'NA','canaps053101':'NA','canaps053801':'NA','canaps060302':'NA','canaps061005':'NA','canaps062201':'NA','canaps062401':'NA','canaps062701':'NA','canaps063601':'NA','canaps064201':'NA','canaps090901':'NA','canaps054301':'NA','canaps090605':'NA', \
                'canaps030801':'NA','canaps030701':'NA','canaps030901':'NA','canaps055201':'NA','canaps060303':'NA','canaps066201':'NA','canaps091101':'NA','canaps020201':'NA','canaps054001':'NA','canaps062101':'NA','canaps063401':'NA','canaps063901':'NA','canaps064301':'NA','canaps066001':'NA','canaps129401':'ARC','canaps129501':'NA', \
                'canaps080801':'NA','canaps090606':'NA','canaps090805':'NA','canaps090806':'NA','canaps092201':'NA','canaps092401':'NA','canaps092501':'NA','canaps092901':'NA','canaps093401':'NA','canaps093501':'NA','canaps093701':'NA','canaps094201':'NA','canaps094202':'NA','canaps094401':'NA','canaps100315':'NA','canaps101803':'NA','canaps104301':'NA','canaps105301':'NA', \
                'canaps105601':'NA','canaps105604':'NA','canaps106800':'NA','canaps129103':'NA','canaps129202':'NA','canaps030120':'NA','canaps090402':'NA','canaps090601':'NA','canaps094301':'NA','canaps101601':'NA','canaps102301':'NA','canaps102701':'NA','canaps103205':'NA','canaps105001':'NA', \
                'canaps010401':'NA','canaps031101':'NA','canaps050801':'NA','canaps100316':'NA','canaps102201':'NA','canaps100138':'NA','canaps060401':'NA','canaps060501':'NA','canaps061001':'NA','canaps070101':'NA','canaps100106':'NA','canaps060101':'NA','canaps060403':'NA','canaps060413':'NA','canaps060414':'NA','canaps060602':'NA','canaps060901':'NA','canaps061501':'NA','canaps061601':'NA','canaps100108':'NA','canaps100109':'NA', \
                'canaps030310':'NA','canaps040206':'NA','canaps050102':'NA','canaps050104':'NA','canaps050109':'NA','canaps050110':'NA','canaps050115':'NA','canaps050128':'NA','canaps054501':'NA','canaps060421':'NA','canaps060430':'NA','canaps060432':'NA',
                'canaps060433':'NA','canaps060512':'NA','canaps060609':'NA','canaps060709':'NA','canaps060903':'NA','canaps061302':'NA','canaps061402':'NA','canaps061502':'NA','canaps061702':'NA','canaps065601':'NA','canaps065801':'NA','canaps070119':'NA',
                'canaps080211':'NA','canaps080402':'NA','canaps090121':'NA','canaps090130':'NA','canaps090218':'NA','canaps090227':'NA','canaps090302':'NA','canaps100112':'NA','canaps100121':'NA','canaps100402':'NA','canaps100701':'NA','canaps129003':'NA','canaps060402':'NA','canaps060415':'NA','canaps080109':'NA',
                'canaps030601':'NA','canaps060417':'NA','canaps060418':'NA','canaps060419':'NA','canaps060514':'NA','canaps060515':'NA','canaps061602':'NA','canaps061701':'NA','canaps063101':'NA','canaps063201':'NA','canaps080110':'NA','canaps080209':'NA','canaps090122':'NA',
                'canaps100120':'NA','canaps100122':'NA','canaps100124':'NA','canaps101201':'NA','canaps100302':'NA','canaps061301':'NA','canaps090219':'NA','canaps060105':'NA','canaps060806':'NA','canaps030115':'NA','canaps030116':'NA','canaps100123':'NA','canaps101001':'NA',\
                'canaps030117':'NA','canaps060707':'NA','canaps060807':'NA','canaps060420':'NA','canaps060422':'NA','canaps060607':'NA','canaps062801':'NA','canaps062901':'NA','canaps100129':'NA','canaps060423':'NA','canaps060424':'NA','canaps100303':'NA','canaps062902':'NA','canaps064501':'NA',
                'canaps050123':'NA','canaps050306':'NA','canaps050307':'NA','canaps050309':'NA','canaps064502':'EU','canaps091201':'NA','canaps091701':'NA','canaps100304':'NA','canaps129002':'NA','canaps100307':'NA','canaps060106':'NA','canaps030302':'NA','canaps050134':'NA','canaps060434':'NA','canaps090228':'NA','canaps081001':'NA',
                'canaps030201':'NA','canaps100308':'NA','canaps103502':'NA','canaps107100':'NA','canaps010102':'NA','canaps010301':'NA','canaps010501':'NA','canaps010602':'NA','canaps010701':'NA','canaps010801':'NA','canaps011001':'NA','canaps020101':'NA','canaps020104':'NA','canaps030113':'NA','canaps030301':'NA',
                'canaps030311':'NA','canaps040202':'NA','canaps050107':'NA','canaps050108':'NA','canaps050112':'NA','canaps050120':'NA','canaps050130':'NA','canaps050135':'NA','canaps050203':'NA','canaps050304':'NA','canaps050802':'NA','canaps051401':'NA','canaps051601':'NA','canaps051701':'NA','canaps054101':'NA',
                'canaps060425':'NA','canaps060435':'NA','canaps060610':'NA','canaps060904':'NA','canaps061103':'NA','canaps064701':'NA','canaps070301':'NA','canaps070501':'NA','canaps080701':'NA','canaps090114':'NA','canaps092701':'NA','canaps093001':'NA','canaps093202':'NA','canaps093301':'NA','canaps093601':'NA',
                'canaps094101':'NA','canaps100105':'NA','canaps100130':'NA','canaps100131':'NA','canaps100140':'NA','canaps100313':'NA','canaps100401':'NA','canaps100801':'NA','canaps101004':'NA','canaps101602':'NA','canaps101702':'NA','canaps102706':'NA','canaps119003':'NA','canaps119004':'NA','canaps129203':'NA',
                #SEARCH
                'searchctr':'NA','searchoak':'NA','searchyrk':'NA','searchgfp':'NA','searchjst':'NA','searcholf':'NA','searchpns':'NA', \
                #AQS
                'aqs010510001':'NA','aqs010731009':'NA','aqs010735003':'NA','aqs010972005':'NA','aqs011170004':'NA','aqs011190002':'NA','aqs020680003':'NA', \
                'aqs040038001':'NA','aqs040058001':'NA','aqs040070010':'NA','aqs040128000':'NA','aqs040131010':'NA','aqs040132001':'NA','aqs040134008':'NA','aqs040134011':'NA','aqs040137003':'NA','aqs040137020':'NA','aqs040137021':'NA', \
                'aqs040137022':'NA','aqs040137024':'NA','aqs040139508':'NA','aqs040139702':'NA','aqs040170119':'NA','aqs040190021':'NA','aqs040191018':'NA','aqs040191020':'NA','aqs040213007':'NA','aqs040217001':'NA','aqs040218001':'NA', \
                'aqs050350005':'NA','aqs051010002':'NA','aqs051130003':'NA','aqs051191002':'NA','aqs060070007':'NA','aqs060090001':'NA','aqs060111002':'NA','aqs060131002':'NA','aqs060194001':'NA','aqs060254003':'NA','aqs060254004':'NA', \
                'aqs060270101':'NA','aqs060290007':'NA','aqs060290011':'NA','aqs060310500':'NA','aqs060390004':'NA','aqs060430003':'NA','aqs060470003':'NA','aqs060650008':'NA','aqs060651016':'NA','aqs060652002':'NA','aqs060670011':'NA', \
                'aqs060675003':'NA','aqs060690002':'NA','aqs060690003':'NA','aqs060710005':'NA','aqs060710012':'NA','aqs060711234':'NA','aqs060719002':'NA','aqs060732007':'NA','aqs060794002':'NA','aqs060798006':'NA','aqs060830008':'NA', \
                'aqs060831013':'NA','aqs060831014':'NA','aqs060831018':'NA','aqs060831021':'NA','aqs060831025':'NA','aqs060833001':'NA','aqs060834003':'NA','aqs060852006':'NA','aqs060893003':'NA','aqs060950005':'NA','aqs061070006':'NA', \
                'aqs061070009':'NA','aqs061110009':'NA','aqs061113001':'NA','aqs061130004':'NA','aqs080013001':'NA','aqs080350004':'NA','aqs080410013':'NA','aqs080590006':'NA','aqs080671004':'NA','aqs080677003':'NA','aqs080690007':'NA', \
                'aqs080690011':'NA','aqs080830101':'NA','aqs090050005':'NA','aqs100010002':'NA','aqs100031007':'NA','aqs100031010':'NA','aqs100051003':'NA','aqs120013011':'NA','aqs120030002':'NA','aqs120050006':'NA','aqs120310077':'NA', \
                'aqs120550003':'NA','aqs120570081':'NA','aqs120573002':'NA','aqs120590004':'NA','aqs120860029':'NA','aqs121010005':'NA','aqs121035002':'NA','aqs121290001':'NA','aqs130210012':'NA','aqs130550001':'NA','aqs130850001':'NA', \
                'aqs131510002':'NA','aqs132130003':'NA','aqs132230003':'NA','aqs132470001':'NA','aqs132611001':'NA','aqs170230001':'NA','aqs170491001':'NA','aqs170650002':'NA','aqs170971007':'NA','aqs171170002':'NA','aqs171570001':'NA', \
                'aqs171971011':'NA','aqs190170011':'NA','aqs191370002':'NA','aqs191630014':'NA','aqs191770006':'NA','aqs201070002':'NA','aqs201730001':'NA','aqs201910002':'NA','aqs201950001':'NA','aqs210150003':'NA','aqs210430500':'NA', \
                'aqs210470006':'NA','aqs210610501':'NA','aqs210910012':'NA','aqs211010014':'NA','aqs211390003':'NA','aqs211850004':'NA','aqs212130004':'NA','aqs212218001':'NA','aqs212270008':'NA','aqs220050004':'NA','aqs220170001':'NA', \
                'aqs220190008':'NA','aqs220190009':'NA','aqs220330013':'NA','aqs220470009':'NA','aqs220470012':'NA','aqs220630002':'NA','aqs220770001':'NA','aqs220890003':'NA','aqs220930002':'NA','aqs220950002':'NA','aqs230052003':'NA', \
                'aqs230090102':'NA','aqs230090103':'NA','aqs230130004':'NA','aqs230173001':'NA','aqs230194008':'NA','aqs230290032':'NA','aqs230310038':'NA','aqs240030014':'NA','aqs240090011':'NA','aqs240130001':'NA','aqs240150003':'NA', \
                'aqs240170010':'NA','aqs240230002':'NA','aqs240251001':'NA','aqs240290002':'NA','aqs240430009':'NA','aqs250010002':'NA','aqs250070001':'NA','aqs250154002':'NA','aqs261130001':'NA','aqs261530001':'NA','aqs270750005':'NA', \
                'aqs271370034':'NA','aqs280110001':'NA','aqs280470008':'NA','aqs290370003':'NA','aqs290390001':'NA','aqs290470003':'NA','aqs290470005':'NA','aqs290470006':'NA','aqs290490001':'NA','aqs291130003':'NA','aqs291370001':'NA', \
                'aqs291570001':'NA','aqs291831002':'NA','aqs291831004':'NA','aqs291860005':'NA','aqs291890005':'NA','aqs291890014':'NA','aqs300298001':'NA','aqs311090016':'NA','aqs320010002':'NA','aqs320030022':'NA','aqs320030023':'NA', \
                'aqs320031019':'NA','aqs320330101':'NA','aqs330074001':'NA','aqs330115001':'NA','aqs340071001':'NA','aqs340110007':'NA','aqs340150002':'NA','aqs340190001':'NA','aqs340230011':'NA','aqs340273001':'NA','aqs340290006':'NA', \
                'aqs340315001':'NA','aqs350010029':'NA','aqs350130008':'NA','aqs350130020':'NA','aqs350130022':'NA','aqs350151005':'NA','aqs350171003':'NA','aqs350431001':'NA','aqs350439004':'NA','aqs350450009':'NA','aqs350450018':'NA', \
                'aqs350451005':'NA','aqs360130011':'NA','aqs360270007':'NA','aqs360310002':'NA','aqs360310003':'NA','aqs360337003':'NA','aqs360410005':'NA','aqs360430005':'NA','aqs360450002':'NA','aqs360530006':'NA','aqs360631006':'NA', \
                'aqs360650004':'NA','aqs360750003':'NA','aqs360790005':'NA','aqs360830004':'NA','aqs360910004':'NA','aqs361010003':'NA','aqs361030004':'NA','aqs361111005':'NA','aqs361173001':'NA','aqs370330001':'NA','aqs370650099':'NA', \
                'aqs370670028':'NA','aqs370671008':'NA','aqs370750001':'NA','aqs370870035':'NA','aqs370870036':'NA','aqs371010002':'NA','aqs371090004':'NA','aqs371170001':'NA','aqs371191009':'NA','aqs371450003':'NA','aqs371570099':'NA', \
                'aqs371590021':'NA','aqs380070002':'NA','aqs380130004':'NA','aqs380250003':'NA','aqs380530002':'NA','aqs380570004':'NA','aqs380650002':'NA','aqs390230001':'NA','aqs390271002':'NA','aqs390410002':'NA','aqs390550004':'NA', \
                'aqs390610010':'NA','aqs390830002':'NA','aqs390950034':'NA','aqs390970007':'NA','aqs391090005':'NA','aqs391351001':'NA','aqs391550009':'NA','aqs391550011':'NA','aqs400019009':'NA','aqs400159008':'NA','aqs400219002':'NA', \
                'aqs400370144':'NA','aqs400430860':'NA','aqs400719010':'NA','aqs400871073':'NA','aqs400979014':'NA','aqs401159004':'NA','aqs401210415':'NA','aqs401430174':'NA','aqs420010002':'NA','aqs420070002':'NA','aqs420070005':'NA', \
                'aqs420270100':'NA','aqs420290100':'NA','aqs420334000':'NA','aqs420430401':'NA','aqs420550001':'NA','aqs420590002':'NA','aqs420630004':'NA','aqs420890002':'NA','aqs420990301':'NA','aqs421174000':'NA','aqs421255001':'NA', \
                'aqs450010001':'NA','aqs450150002':'NA','aqs450190046':'NA','aqs450210002':'NA','aqs450250001':'NA','aqs450290002':'NA','aqs450310003':'NA','aqs450370001':'NA','aqs450730001':'NA','aqs450770002':'NA','aqs450790021':'NA', \
                'aqs450791001':'NA','aqs450830009':'NA','aqs460330132':'NA','aqs470010101':'NA','aqs470090101':'NA','aqs470090102':'NA','aqs470370026':'NA','aqs470651011':'NA','aqs470654003':'NA','aqs470890002':'NA','aqs470930021':'NA', \
                'aqs471210104':'NA','aqs471490101':'NA','aqs471550101':'NA','aqs471571004':'NA','aqs471632002':'NA','aqs471650007':'NA','aqs471650101':'NA','aqs471870106':'NA','aqs471890103':'NA','aqs480290052':'NA','aqs480290059':'NA', \
                'aqs480430101':'NA','aqs481210034':'NA','aqs481830001':'NA','aqs482010029':'NA','aqs482030002':'NA','aqs482090614':'NA','aqs482450101':'NA','aqs483670081':'NA','aqs484230007':'NA','aqs484390075':'NA','aqs484530020':'NA', \
                'aqs490037001':'NA','aqs490370101':'NA','aqs490530130':'NA','aqs500030004':'NA','aqs500070007':'NA','aqs510330001':'NA','aqs510410004':'NA','aqs510610002':'NA','aqs511130003':'NA','aqs511390004':'NA','aqs511630003':'NA', \
                'aqs511970002':'NA','aqs530090013':'NA','aqs530530012':'NA','aqs530531010':'NA','aqs550090026':'NA','aqs550210015':'NA','aqs550290004':'NA','aqs550390006':'NA','aqs550410007':'NA','aqs550590019':'NA','aqs550610002':'NA', \
                'aqs550710007':'NA','aqs550730012':'NA','aqs550850004':'NA','aqs550890008':'NA','aqs550890009':'NA','aqs551110007':'NA','aqs551170006':'NA','aqs551250001':'NA','aqs551270005':'NA','aqs560050123':'NA','aqs560050456':'NA', \
                'aqs560350099':'NA','aqs560350100':'NA','aqs560370200':'NA','aqs560391011':'NA','aqs080677001':'NA','aqs090019003':'NA','aqs150030011':'O','aqs270370020':'NA','aqs270370423':'NA','aqs380570102':'NA','aqs380570124':'NA','aqs460710001':'NA', \
                'aqs230230006':'NA','aqs230290019':'NA','aqs230310040':'NA','aqs240338003':'NA','aqs270953051':'NA','aqs290190011':'NA','aqs370510008':'NA','aqs371290002':'NA','aqs371630003':'NA','aqs371630004':'NA','aqs120310086':'NA', \
                'aqs040133003':'NA','aqs040134003':'NA','aqs040139997':'NA','aqs040190002':'NA','aqs040191028':'NA','aqs040213003':'NA','aqs360912001':'NA','aqs400270044':'NA','aqs490353001':'NA', \
                'aqs051190007':'NA','aqs060050002':'NA','aqs060070002':'NA','aqs060130002':'NA','aqs060190007':'NA','aqs060195001':'NA','aqs060290008':'NA','aqs060290014':'NA','aqs060290232':'NA','aqs060295001':'NA','aqs060296001':'NA', \
                'aqs060333001':'NA','aqs060370002':'NA','aqs060370016':'NA','aqs060370113':'NA','aqs060371002':'NA','aqs060371103':'NA','aqs060371201':'NA','aqs060371701':'NA','aqs060372005':'NA','aqs060374002':'NA','aqs060375005':'NA', \
                'aqs060376012':'NA','aqs060379033':'NA','aqs060410001':'NA','aqs060450008':'NA','aqs060530002':'NA','aqs060531003':'NA','aqs060550003':'NA','aqs060570005':'NA','aqs060590007':'NA','aqs060591003':'NA','aqs060592022':'NA', \
                'aqs060595001':'NA','aqs060610002':'NA','aqs060610004':'NA','aqs060610006':'NA','aqs060650012':'NA','aqs060655001':'NA','aqs060656001':'NA','aqs060670002':'NA','aqs060670010':'NA','aqs060792006':'NA','aqs060793001':'NA', \
                'aqs060971003':'NA','aqs061111004':'NA','aqs061112002':'NA','aqs061131003':'NA','aqs080050002':'NA','aqs080410016':'NA','aqs080590005':'NA','aqs080691004':'NA','aqs100031013':'NA','aqs120118002':'NA','aqs120210004':'NA', \
                'aqs120330004':'NA','aqs120571035':'NA','aqs120571065':'NA','aqs120830003':'NA','aqs120830004':'NA','aqs120860027':'NA','aqs121012001':'NA','aqs121030004':'NA','aqs121111002':'NA','aqs121130015':'NA','aqs121151005':'NA', \
                'aqs121171002':'NA','aqs121272001':'NA','aqs121275002':'NA','aqs170310076':'NA','aqs170314007':'NA','aqs170314201':'NA','aqs170317002':'NA','aqs170436001':'NA','aqs171193007':'NA','aqs171431001':'NA','aqs171613002':'NA', \
                'aqs171630010':'NA','aqs220330009':'NA','aqs220331001':'NA','aqs220470007':'NA','aqs220511001':'NA','aqs220550007':'NA','aqs220570004':'NA','aqs220730004':'NA','aqs240053001':'NA','aqs250170009':'NA','aqs250250042':'NA', \
                'aqs320030020':'NA','aqs320030043':'NA','aqs320030071':'NA','aqs320030072':'NA','aqs320030073':'NA','aqs320030075':'NA','aqs320030538':'NA','aqs320030601':'NA','aqs320032002':'NA','aqs320310020':'NA','aqs320310025':'NA', \
                'aqs320311005':'NA','aqs320312009':'NA','aqs340170006':'NA','aqs340210005':'NA','aqs350010023':'NA','aqs350010024':'NA','aqs350010027':'NA','aqs350011012':'NA','aqs350011013':'NA','aqs361192004':'NA','aqs371190041':'NA', \
                'aqs371830014':'NA','aqs401091037':'NA','aqs401430137':'NA','aqs420031005':'NA','aqs450790007':'NA','aqs480610006':'NA','aqs480850005':'NA','aqs481130087':'NA','aqs482010026':'NA','aqs482010046':'NA','aqs482010051':'NA', \
                'aqs482010055':'NA','aqs482010062':'NA','aqs482010066':'NA','aqs482011035':'NA','aqs482011039':'NA','aqs482210001':'NA','aqs482311006':'NA','aqs482450009':'NA','aqs482450011':'NA','aqs482450022':'NA','aqs482450102':'NA', \
                'aqs482450628':'NA','aqs482510003':'NA','aqs483390078':'NA','aqs483550025':'NA','aqs483550026':'NA','aqs483611001':'NA','aqs483611100':'NA','aqs483970001':'NA','aqs484391002':'NA','aqs484392003':'NA','aqs484393009':'NA', \
                'aqs490353006':'NA','aqs510595001':'NA','aqs551091002':'NA','aqs800020001':'NA','aqs060250006':'NA','aqs060731006':'NA','aqs060870003':'NA','aqs110010025':'NA','aqs120972002':'NA','aqs121152002':'NA','aqs220190002':'NA', \
                'aqs350130021':'NA','aqs380171004':'NA','aqs400170101':'NA','aqs482570005':'NA','aqs550030010':'NA','aqs550270007':'NA','aqs550790026':'NA','aqs800020005':'NA','aqs800020018':'NA', \
                'aqs010010003':'NA','aqs010270001':'NA','aqs010499991':'NA','aqs010610001':'NA','aqs010790002':'NA','aqs010970025':'NA', \
                'aqs011011001':'NA','aqs040190013':'NA','aqs050199991':'NA','aqs060110002':'NA','aqs060150002':'NA','aqs060150003':'NA','aqs060172003':'NA','aqs060192009':'NA','aqs060379034':'NA','aqs060390500':'NA','aqs060410002':'NA','aqs060651004':'NA', \
                'aqs060651010':'NA','aqs060651999':'NA','aqs060731201':'NA','aqs061070005':'NA','aqs061070008':'NA','aqs090159991':'NA','aqs120310106':'NA','aqs120350004':'NA','aqs120730002':'NA','aqs120779991':'NA','aqs120990006':'NA','aqs131910002':'NA', \
                'aqs132319991':'NA','aqs170190007':'NA','aqs170191001':'NA','aqs170331001':'NA','aqs170859991':'NA','aqs171170001':'NA','aqs171199991':'NA','aqs171332001':'NA','aqs180030001':'NA','aqs180839991':'NA','aqs181630013':'NA','aqs181699991':'NA', \
                'aqs191210003':'NA','aqs191530024':'NA','aqs191770004':'NA','aqs201619991':'NA','aqs201730018':'NA','aqs210350004':'NA','aqs211759991':'NA','aqs211770005':'NA','aqs212210001':'NA','aqs212219991':'NA','aqs212299991':'NA','aqs220190006':'NA', \
                'aqs230039991':'NA','aqs230090003':'NA','aqs230194007':'NA','aqs230199991':'NA','aqs240190004':'NA','aqs240199991':'NA','aqs240210034':'NA','aqs240451004':'NA','aqs260492001':'NA','aqs260630905':'NA','aqs260770906':'NA','aqs261330901':'NA', \
                'aqs261579991':'NA','aqs261611001':'NA','aqs261619991':'NA','aqs261659991':'NA','aqs270031001':'NA','aqs270710101':'NA','aqs271376317':'NA','aqs271636015':'NA','aqs280370001':'NA','aqs280450719':'NA','aqs280890002':'NA','aqs281250001':'NA', \
                'aqs281619991':'NA','aqs310550032':'NA','aqs311079991':'NA','aqs320038000':'NA','aqs330099991':'NA','aqs340090003':'NA','aqs360291005':'NA','aqs360310005':'NA','aqs360319991':'NA','aqs361099991':'NA','aqs370230004':'NA','aqs370319991':'NA', \
                'aqs371139991':'NA','aqs371230099':'NA','aqs371239991':'NA','aqs380570101':'NA','aqs380650101':'NA','aqs390030002':'NA','aqs390479991':'NA','aqs391219991':'NA','aqs391532004':'NA','aqs400190296':'NA','aqs400190297':'NA','aqs400310649':'NA', \
                'aqs400330680':'NA','aqs400670671':'NA','aqs400770440':'NA','aqs400819005':'NA','aqs400819024':'NA','aqs400850300':'NA','aqs400892001':'NA','aqs401430177':'NA','aqs410290008':'NA','aqs420019991':'NA','aqs420279991':'NA','aqs420479991':'NA', \
                'aqs420859991':'NA','aqs421119991':'NA','aqs450030004':'NA','aqs450070005':'NA','aqs450230002':'NA','aqs450290001':'NA','aqs450451003':'NA','aqs450890001':'NA','aqs460110003':'NA','aqs461270003':'NA','aqs470250001':'NA','aqs470259991':'NA', \
                'aqs470419991':'NA','aqs470550001':'NA','aqs471190016':'NA','aqs480710900':'NA','aqs480710902':'NA','aqs480710903':'NA','aqs483739991':'NA','aqs484570101':'NA','aqs510150004':'NA','aqs510719991':'NA','aqs511479991':'NA','aqs511870002':'NA', \
                'aqs530530027':'NA','aqs540219991':'NA','aqs540939991':'NA','aqs550170001':'NA','aqs550210005':'NA','aqs550210008':'NA','aqs550210013':'NA','aqs550270001':'NA','aqs550330003':'NA','aqs550370001':'NA','aqs550430003':'NA','aqs551199991':'NA', \
                'aqs551230008':'NA','aqs560370898':'NA','aqs180730002':'NA','aqs180770001':'NA','aqs180830004':'NA','aqs181270015':'NA','aqs181270016':'NA','aqs181291002':'NA','aqs181530001':'NA','aqs181631001':'NA','aqs181631002':'NA','aqs181671012':'NA', \
                'aqs210150006':'NA','aqs210150007':'NA','aqs212230004':'NA','aqs261470904':'NA','aqs261470905':'NA','aqs550110004':'NA','aqs550110007':'NA','aqs550350010':'NA','aqs550591001':'NA','aqs551230003':'NA','aqs350290003':'NA','aqs481211032':'NA', \
                'aqs481390016':'NA','aqs050970001':'NA','aqs060530006':'NA','aqs060675002':'NA','aqs060770009':'NA','aqs060831012':'NA','aqs060831015':'NA','aqs060831016':'NA','aqs060831019':'NA','aqs061010002':'NA','aqs061110004':'NA','aqs061110005':'NA','aqs120813002':'NA', \
                'aqs120860030':'NA','aqs121091003':'NA','aqs191131015':'NA','aqs191632011':'NA','aqs210610500':'NA','aqs220110002':'NA','aqs220190007':'NA','aqs220470002':'NA','aqs230090101':'NA','aqs340010005':'NA','aqs450870001':'NA', \
                'aqs010710020':'NA','aqs010830004':'NA','aqs021221004':'NA','aqs040139993':'NA','aqs060012005':'NA','aqs060231005':'NA','aqs060610007':'NA','aqs060773002':'NA', \
                'aqs060791004':'NA','aqs060794001':'NA','aqs060831011':'NA','aqs060831017':'NA','aqs060831026':'NA','aqs060831027':'NA','aqs060831030':'NA','aqs060834004':'NA','aqs060835001':'NA','aqs060991005':'NA','aqs061110006':'NA','aqs090131001':'NA', \
                'aqs160210003':'NA','aqs171050002':'NA','aqs180050007':'NA','aqs180330002':'NA','aqs180510010':'NA','aqs180730003':'NA','aqs181090004':'NA','aqs181290001':'NA','aqs181290002':'NA','aqs181410012':'NA','aqs181470002':'NA','aqs181470006':'NA', \
                'aqs181470008':'NA','aqs181550001':'NA','aqs181830003':'NA','aqs211771004':'NA','aqs220330008':'NA','aqs220470006':'NA','aqs260050001':'NA','aqs260050002':'NA','aqs260210013':'NA','aqs260430901':'NA','aqs260430902':'NA','aqs261050005':'NA', \
                'aqs261110941':'NA','aqs261390006':'NA','aqs270177416':'NA','aqs270495302':'NA','aqs270757630':'NA','aqs270757631':'NA','aqs270757632':'NA','aqs270757634':'NA','aqs270757635':'NA','aqs271710007':'NA','aqs280030004':'NA','aqs280190001':'NA', \
                'aqs280450001':'NA','aqs280930001':'NA','aqs281070001':'NA','aqs291290001':'NA','aqs291830010':'NA','aqs292110001':'NA','aqs300710010':'NA','aqs300750001':'NA','aqs300830001':'NA','aqs300870001':'NA','aqs301110086':'NA','aqs330150013':'NA', \
                'aqs350015010':'NA','aqs350281002':'NA','aqs350451233':'NA','aqs370370004':'NA','aqs371050002':'NA','aqs380130001':'NA','aqs380130002':'NA','aqs380570103':'NA','aqs380570104':'NA','aqs381010114':'NA','aqs390010013':'NA','aqs391030004':'NA', \
                'aqs391291001':'NA','aqs391450014':'NA','aqs400719003':'NA','aqs400770441':'NA','aqs401010160':'NA','aqs401010167':'NA','aqs401110152':'NA','aqs401110153':'NA','aqs401359015':'NA','aqs410591003':'NA','aqs420010001':'NA','aqs420150011':'NA', \
                'aqs440030002':'NA','aqs450750003':'NA','aqs450791006':'NA','aqs461094003':'NA','aqs461270001':'NA','aqs461270002':'NA','aqs470110004':'NA','aqs470310004':'NA','aqs470430009':'NA','aqs470630003':'NA','aqs470750002':'NA','aqs470750003':'NA', \
                'aqs470850020':'NA','aqs471050106':'NA','aqs471190106':'NA','aqs471250009':'NA','aqs471251010':'NA','aqs471310004':'NA','aqs471410004':'NA','aqs471451020':'NA','aqs471572005':'NA','aqs481390017':'NA','aqs483491051':'NA','aqs490530007':'NA', \
                'aqs511650002':'NA','aqs530150014':'NA','aqs530330023':'NA','aqs540250001':'NA','aqs540990003':'NA','aqs550210016':'NA','aqs560090819':'NA','aqs380910001':'NA','aqs060773003':'NA','aqs060792004':'NA','aqs220430001':'NA', \
                'aqs150010006':'O','aqs120570110':'NA','aqs460711001':'NA','aqs530570013':'NA','aqs060832012':'NA','aqs450110001':'NA','aqs780200001':'O','aqs060798005':'NA','aqs481391044':'NA','aqs483091037':'NA', \
                'aqs040270006':'NA','aqs060231004':'NA','aqs060251003':'NA','aqs060311004':'NA','aqs060650009':'NA','aqs060710001':'NA','aqs060773005':'NA','aqs060792002':'NA','aqs060798001':'NA','aqs060832004':'NA','aqs060990006':'NA', \
                'aqs061010003':'NA','aqs150030010':'NA','aqs160550003':'NA','aqs180630003':'NA','aqs180890022':'NA','aqs210290006':'NA','aqs210590005':'NA','aqs211451024':'NA','aqs230031100':'NA','aqs280590006':'NA','aqs290770036':'NA','aqs350151004':'NA', \
                'aqs380150003':'NA','aqs390090004':'NA','aqs420130801':'NA','aqs420210011':'NA','aqs420490003':'NA','aqs420730015':'NA','aqs450030003':'NA','aqs460990007':'NA','aqs460990008':'NA','aqs470110102':'NA','aqs471070101':'NA','aqs480290622':'NA', \
                'aqs480391016':'NA','aqs481390015':'NA','aqs500210002':'NA','aqs510360002':'NA','aqs511530009':'NA','aqs511650003':'NA','aqs530570018':'NA','aqs800020004':'NA','aqs800020010':'NA','aqs800020017':'NA','aqscc0110002':'NA', \
                'aqs010030010':'NA','aqs010550011':'NA','aqs010690004':'NA','aqs010731005':'NA','aqs010731010':'NA','aqs011130002':'NA','aqs011210003':'NA','aqs011250010':'NA','aqs040135100':'NA','aqs040139704':'NA','aqs040191030':'NA','aqs040191034':'NA', \
                'aqs040213001':'NA','aqs060170010':'NA','aqs060210002':'NA','aqs060210003':'NA','aqs060310004':'NA','aqs060430033':'NA','aqs060450009':'NA','aqs060530005':'NA','aqs060530008':'NA','aqs060659003':'NA','aqs060710017':'NA','aqs060790005':'NA', \
                'aqs060792001':'NA','aqs060850002':'NA','aqs060890004':'NA','aqs060890007':'NA','aqs060890009':'NA','aqs060932001':'NA','aqs061030005':'NA','aqs061090005':'NA','aqs081230009':'NA','aqs120010025':'NA','aqs120230002':'NA','aqs120690002':'NA', \
                'aqs120730012':'NA','aqs120910002':'NA','aqs130210013':'NA','aqs132151003':'NA','aqs170010006':'NA','aqs170010007':'NA','aqs170190004':'NA','aqs170831001':'NA','aqs171132003':'NA','aqs171150013':'NA','aqs171670010':'NA','aqs200450004':'NA', \
                'aqs201030003':'NA','aqs210130002':'NA','aqs210670001':'NA','aqs210890007':'NA','aqs210930006':'NA','aqs211110051':'NA','aqs211130001':'NA','aqs211490001':'NA','aqs211830032':'NA','aqs211930003':'NA','aqs211950002':'NA','aqs211990003':'NA', \
                'aqs212210013':'NA','aqs221010003':'NA','aqs221030002':'NA','aqs230191100':'NA','aqs230290031':'NA','aqs280010004':'NA','aqs280330002':'NA','aqs280450003':'NA','aqs280470009':'NA','aqs280590007':'NA','aqs280750003':'NA','aqs280810005':'NA', \
                'aqs320030007':'NA','aqs330050007':'NA','aqs330090010':'NA','aqs350130023':'NA','aqs350250008':'NA','aqs360130006':'NA','aqs360150003':'NA','aqs371070004':'NA','aqs371470006':'NA','aqs400310647':'NA','aqs400310651':'NA','aqs400670670':'NA', \
                'aqs401090096':'NA','aqs450910006':'NA','aqs460930001':'NA','aqs461030020':'NA','aqs470990002':'NA','aqs471050108':'NA','aqs471050109':'NA','aqs471632003':'NA','aqs484690003':'NA','aqs484790016':'NA','aqs490050004':'NA','aqs550630012':'NA', \
                'aqs560130232':'NA','aqs560050892':'NA','aqs060295002':'NA','aqs490472003':'NA','aqs040191024':'NA','aqs060831020':'NA','aqs080030002':'NA','aqs080050003':'NA','aqs080410012':'NA','aqs080770600':'NA','aqs120095001':'NA','aqs120712001':'NA','aqs170310050':'NA','aqs490190101':'NA','aqs541071002':'NA','aqs551171002':'NA',\
                'aqs010730023':'NA','aqs010890014':'NA','aqs040130019':'NA','aqs040131004':'NA','aqs040132005':'NA','aqs040133002':'NA','aqs040191011':'NA','aqs040191032':'NA','aqs060010007':'NA','aqs060011001':'NA','aqs060131004':'NA','aqs060133001':'NA',
                'aqs060190008':'NA','aqs060190242':'NA','aqs060250005':'NA','aqs060290010':'NA','aqs060371301':'NA','aqs060658001':'NA','aqs060659001':'NA','aqs060670006':'NA','aqs060670012':'NA','aqs060670013':'NA','aqs060710306':'NA','aqs060711004':'NA',
                'aqs060712002':'NA','aqs060714001':'NA','aqs060714003':'NA','aqs060719004':'NA','aqs060730001':'NA','aqs060730003':'NA','aqs060730006':'NA','aqs060731001':'NA','aqs060731002':'NA','aqs060731008':'NA','aqs060731010':'NA','aqs060750005':'NA',
                'aqs060771002':'NA','aqs060811001':'NA','aqs060830011':'NA','aqs060831008':'NA','aqs060832011':'NA','aqs060850005':'NA','aqs060851001':'NA','aqs060870004':'NA','aqs060870006':'NA','aqs060870007':'NA','aqs060950004':'NA','aqs060953003':'NA',
                'aqs060970003':'NA','aqs060990005':'NA','aqs061072002':'NA','aqs061110007':'NA','aqs061112003':'NA','aqs080130011':'NA','aqs080310002':'NA','aqs080310014':'NA','aqs080590002':'NA','aqs080590011':'NA','aqs110010041':'NA','aqs110010043':'NA',
                'aqs120090007':'NA','aqs120094001':'NA','aqs120110031':'NA','aqs120112003':'NA','aqs120310100':'NA','aqs120330018':'NA','aqs120330024':'NA','aqs120712002':'NA','aqs120713002':'NA','aqs120730013':'NA','aqs120814012':'NA','aqs120814013':'NA',
                'aqs120950008':'NA','aqs120952002':'NA','aqs120990009':'NA','aqs120990020':'NA','aqs121030018':'NA','aqs121056005':'NA','aqs121056006':'NA','aqs121151006':'NA','aqs150031004':'O','aqs170310064':'NA','aqs170890005':'NA','aqs170971002':'NA',
                'aqs171110001':'NA','aqs171430024':'NA','aqs172010009':'NA','aqs172012001':'NA','aqs200910010':'NA','aqs201730010':'NA','aqs202090021':'NA','aqs220150008':'NA','aqs220330003':'NA','aqs221210001':'NA','aqs240330030':'NA','aqs250092006':'NA',
                'aqs250130008':'NA','aqs270031002':'NA','aqs320310016':'NA','aqs340070003':'NA','aqs350010019':'NA','aqs350011014':'NA','aqs350130017':'NA','aqs350431003':'NA','aqs360010012':'NA','aqs360050110':'NA','aqs360290002':'NA','aqs360551007':'NA',
                'aqs360671015':'NA','aqs360810124':'NA','aqs360930003':'NA','aqs361030002':'NA','aqs361030009':'NA','aqs400270049':'NA','aqs401090033':'NA','aqs401430178':'NA','aqs401431127':'NA','aqs420030008':'NA','aqs421010004':'NA','aqs421010014':'NA',
                'aqs421010024':'NA','aqs421010136':'NA','aqs480290032':'NA','aqs480391004':'NA','aqs481130069':'NA','aqs481130075':'NA','aqs481410029':'NA','aqs481410037':'NA','aqs481410044':'NA','aqs481410055':'NA','aqs481410057':'NA','aqs481410058':'NA',
                'aqs482010024':'NA','aqs482010047':'NA','aqs482010070':'NA','aqs482010075':'NA','aqs482011015':'NA','aqs482011034':'NA','aqs482011050':'NA','aqs482150042':'NA','aqs482150043':'NA','aqs484393011':'NA','aqs484530014':'NA','aqs560350098':'NA',
                'aqs720330008':'NA','aqs800020012':'NA','aqs800060004':'NA','aqs800060006':'NA','aqs800060007':'NA','aqs040130013':'NA','aqs040132004':'NA','aqs040190019':'NA','aqs060831010':'NA','aqs080131001':'NA','aqs080410004':'NA','aqs081230007':'NA',
                'aqs120110003':'NA','aqs120110004':'NA','aqs120571055':'NA','aqs120860021':'NA','aqs170310037':'NA','aqs170970001':'NA','aqs170973001':'NA','aqs171190008':'NA','aqs171191009':'NA','aqs171192007':'NA','aqs171610003':'NA','aqs191530058':'NA',
                'aqs220191003':'NA','aqs220330004':'NA','aqs220512001':'NA','aqs220550003':'NA','aqs220710012':'NA','aqs220730002':'NA','aqs220870002':'NA','aqs240050010':'NA','aqs240330002':'NA','aqs270371007':'NA','aqs310550035':'NA','aqs350011003':'NA',
                'aqs360271003':'NA','aqs360551004':'NA','aqs360610063':'NA','aqs360670014':'NA','aqs360810004':'NA','aqs390811012':'NA','aqs400870073':'NA','aqs420031001':'NA','aqs420070014':'NA','aqs420110001':'NA','aqs420110009':'NA','aqs420170012':'NA',
                'aqs420431100':'NA','aqs420450103':'NA','aqs420690100':'NA','aqs420692006':'NA','aqs420770004':'NA','aqs420791100':'NA','aqs420791101':'NA','aqs420810403':'NA','aqs420850100':'NA','aqs420910013':'NA','aqs420950100':'NA','aqs421250005':'NA',
                'aqs421250200':'NA','aqs421330008':'NA','aqs470650028':'NA','aqs490110001':'NA','aqs490350003':'NA','aqs510590018':'NA','aqs510591004':'NA','aqs511611004':'NA','aqs516000005':'NA','aqs516500004':'NA','aqs720540001':'NA','aqs040130015':'NA',
                'aqs040190018':'NA','aqs080051002':'NA','aqs170310003':'NA','aqs171190012':'NA','aqs171971007':'NA','aqs171971008':'NA','aqs171990001':'NA','aqs180591001':'NA','aqs180892008':'NA','aqs180970031':'NA','aqs180970042':'NA','aqs220710005':'NA',
                'aqs240030019':'NA','aqs261630016':'NA','aqs261632002':'NA','aqs311090011':'NA','aqs350010015':'NA','aqs350010021':'NA','aqs390250002':'NA','aqs390490015':'NA','aqs390490081':'NA','aqs390610006':'NA','aqs390610019':'NA','aqs550890005':'NA',
                'aqs051191003':'NA','aqs170310045':'NA','aqs170431002':'NA','aqs170910001':'NA','aqs171190006':'NA','aqs240430007':'NA','aqs390171004':'NA','aqs180891016':'NA','aqs180970057':'NA','aqs240010012':'NA','aqs261630001':'NA','aqs261630019':'NA',
                'aqs271230003':'NA','aqs390931003':'NA','aqs391032001':'NA','aqs391511009':'NA','aqs518100007':'NA','aqs040131003':'NA','aqs110010017':'NA','aqs170310053':'NA','aqs371190034':'NA','aqs420450002':'NA','aqs420490010':'NA','aqs420710007':'NA','aqs420950017':'NA',
                'aqs040010012':'NA','aqs060371902':'NA','aqs060430004':'NA','aqs121151002':'NA','aqs211110027':'NA','aqs211111021':'NA','aqs240313001':'NA','aqs245100018':'NA','aqs250251003':'NA','aqs420070003':'NA','aqs481090101':'NA','aqs540110006':'NA',
                'aqs540291004':'NA','aqs540690007':'NA','aqs060659002':'NA','aqs120310070':'NA','aqs120990007':'NA','aqs120992004':'NA','aqs191130020':'NA','aqs220570002':'NA','aqs450150042':'NA','aqs530090012':'NA','aqs560391010':'NA',
                'aqs060430005':'NA','aqs061070007':'NA','aqs100031003':'NA','aqs120713001':'NA','aqs120814010':'NA','aqs150010005':'NA','aqs270376018':'NA','aqs320050004':'NA','aqs350131012':'NA','aqs350611001':'NA',
                'aqs060832002':'NA','aqs061090004':'NA','aqs150090101':'NA','aqs450070003':'NA','aqs450791002':'NA','aqs480290036':'NA','aqs480391003':'NA','aqs481130045':'NA','aqs481410027':'NA','aqs481410028':'NA','aqs481671002':'NA','aqs481990002':'NA',
                'aqs482010059':'NA','aqs482011037':'NA','aqs484530003':'NA','aqs160230101':'NA','aqs181270020':'NA','aqs202090001':'NA','aqs260770905':'NA','aqs340230006':'NA','aqs401430127':'NA','aqs510590005':'NA','aqs010970003':'NA','aqs010970028':'NA',
                'aqs011011002':'NA','aqs060730005':'NA','aqs060731007':'NA','aqs061111003':'NA','aqs350431002':'NA','aqs481210033':'NA','aqs482011003':'NA','aqs060010003':'NA','aqs060010005':'NA','aqs060010006':'NA','aqs060130003':'NA','aqs060170011':'NA',
                'aqs060371601':'NA','aqs060375001':'NA','aqs060376002':'NA','aqs060379002':'NA','aqs060531002':'NA','aqs060571001':'NA','aqs060590001':'NA','aqs060592001':'NA','aqs060613001':'NA','aqs060650002':'NA','aqs060710014':'NA','aqs060710015':'NA',
                'aqs060830010':'NA','aqs060831007':'NA','aqs060850004':'NA','aqs060851002':'NA','aqs060852005':'NA','aqs060950002':'NA','aqs060974001':'NA','aqs061130005':'NA','aqs220550005':'NA','aqs250094004':'NA','aqs250270015':'NA','aqs320030016':'NA',
                'aqs481210054':'NA','aqs484230004':'NA','aqs060953002':'NA','aqs061090006':'NA','aqs250051005':'NA','aqs250130003':'NA','aqs350130019':'NA','aqs800060001':'NA','aqs360050083':'NA','aqs481670014':'NA','aqs800020003':'NA',
                'aqs060131003':'NA','aqs060631006':'NA','aqs325100004':'NA','aqs470370011':'EU','aqs800020013':'NA','aqs090070007':'NA','aqs120574004':'NA','aqs170310075':'NA','aqs320030021':'NA','aqs360610010':'NA','aqs360810097':'NA',
                'aqs360810098':'NA','aqs483150050':'NA','aqs484390057':'NA','aqs720330006':'NA','aqs800020014':'NA','aqs040250005':'NA','aqs170310001':'NA','aqs320050008':'NA','aqs320312002':'NA','aqs480850010':'NA','aqs060170012':'NA',
                'aqs121130014':'NA','aqs320030078':'NA','aqs340030005':'NA','aqs481133003':'NA','aqs060250004':'NA','aqs090031003':'NA','aqs051191005':'NA','aqs470931020':'NA','aqs230050027':'NA','aqs510591005':'NA','aqs051430005':'NA','aqs060371602':'NA',
                'aqs170311601':'NA','aqs201770013':'NA','aqs220870009':'NA','aqs320031021':'NA','aqs482010416':'NA','aqs482451035':'NA','aqs530330080':'NA','aqs550250041':'NA','aqs051191008':'NA','aqs060010009':'NA','aqs060012004':'NA','aqs340010006':'NA',
                'aqs340030006':'NA','aqs350490021':'NA','aqs360050133':'NA','aqs360610135':'NA','aqs360850067':'NA','aqs410510080':'NA','aqs560130099':'NA','aqs560410101':'NA','aqs720770001':'NA','aqs060371302':'NA','aqs060650004':'NA','aqs060658005':'NA',
                'aqs060670014':'NA','aqs080310025':'NA','aqs080450012':'NA','aqs080770020':'NA','aqs080830006':'NA','aqs090090027':'NA','aqs230050029':'NA','aqs360715001':'NA','aqs481671034':'NA','aqs482151048':'NA','aqs490490002':'NA','aqs490570002':'NA',
                'aqs080050006':'NA','aqs080450016':'NA','aqs080590013':'NA','aqs080690012':'NA','aqs080970002':'NA','aqs120110033':'NA','aqs130890002':'NA','aqs180970078':'NA','aqs260810020':'NA','aqs330150014':'NA','aqs340130003':'NA','aqs340250005':'NA',
                'aqs350010032':'NA','aqs350610008':'NA','aqs360870005':'NA','aqs390350060':'NA','aqs480271047':'NA','aqs490472002':'NA','aqs560070100':'NA','aqs560350101':'NA','aqs040134004':'NA','aqs040134005':'NA','aqs040134010':'NA','aqs040139706':'NA',
                'aqs060010011':'NA','aqs060392010':'NA','aqs060731016':'NA','aqs060852009':'NA','aqs061072010':'NA','aqs081030005':'NA','aqs081030006':'NA','aqs100051002':'NA','aqs120850007':'NA','aqs191630015':'NA','aqs211110067':'NA','aqs260650012':'NA',
                'aqs295100085':'NA','aqs390610040':'NA','aqs440071010':'NA','aqs490471002':'NA','aqs560030002':'NA','aqs560351002':'NA','aqs560370300':'NA','aqs020900034':'NA','aqs060610003':'NA','aqs060650016':'NA','aqs080519991':'NA','aqs080570003':'NA',
                'aqs080810002':'NA','aqs100032004':'NA','aqs120619991':'NA','aqs171670014':'NA','aqs220870004':'NA','aqs240339991':'NA','aqs300490004':'NA','aqs300630024':'NA','aqs310550019':'NA','aqs320030540':'NA','aqs330150018':'NA','aqs340219991':'NA',
                'aqs340410007':'NA','aqs370119991':'NA','aqs390179991':'NA','aqs420050001':'NA','aqs420110011':'NA','aqs420750100':'NA','aqs420810100':'NA','aqs420950025':'NA','aqs421290008':'NA','aqs471570075':'NA','aqs483819991':'NA','aqs490071003':'NA',
                'aqs490131001':'NA','aqs490170004':'NA','aqs530570020':'NA','aqs560019991':'NA','aqs560210100':'NA','aqs560350700':'NA','aqs560359991':'NA','aqs560370077':'NA','aqs560390008':'NA','aqs051430006':'NA','aqs060070008':'NA','aqs060190011':'NA',
                'aqs060292012':'NA','aqs060374006':'NA','aqs060612002':'NA','aqs060731011':'NA','aqs250051004':'NA','aqs250094005':'NA','aqs250095005':'NA','aqs250213003':'NA','aqs250270024':'NA','aqs300270006':'NA','aqs340070002':'NA','aqs401359021':'NA',
                'aqs420110006':'NA','aqs420690101':'NA','aqs420710012':'NA','aqs421255200':'NA','aqs480611023':'NA','aqs490110004':'NA','aqs490130002':'NA','aqs490471003':'NA','aqs510870014':'NA','aqs560090801':'NA','aqs010010002':'NA','aqs010330044':'NA',
                'aqs010331002':'NA','aqs010331003':'NA','aqs010430003':'NA','aqs010550001':'NA','aqs010731003':'NA','aqs010732006':'NA','aqs010735002':'NA','aqs010736002':'NA','aqs010890021':'NA','aqs010890022':'NA','aqs010990001':'NA','aqs011030007':'NA',
                'aqs011030009':'NA','aqs011030011':'NA','aqs011170005':'NA','aqs011190003':'NA','aqs011251001':'NA','aqs011270003':'NA','aqs020200018':'NA','aqs020200024':'NA','aqs020201004':'NA','aqs021700013':'NA','aqs021850041':'NA','aqs021850042':'NA',
                'aqs022610004':'NA','aqs040032003':'NA','aqs040051004':'NA','aqs040051008':'NA','aqs040079993':'NA','aqs040130016':'NA','aqs040130018':'NA','aqs040131006':'NA','aqs040133004':'NA','aqs040133006':'NA','aqs040133009':'NA','aqs040134007':'NA',
                'aqs040137002':'NA','aqs040139604':'NA','aqs040139701':'NA','aqs040139707':'NA','aqs040139805':'NA','aqs040139994':'NA','aqs040139995':'NA','aqs040139998':'NA','aqs040213009':'NA','aqs040213010':'NA','aqs040250002':'NA','aqs040258033':'NA',
                'aqs040270003':'NA','aqs040270005':'NA','aqs040278011':'NA','aqs050350006':'NA','aqs050350007':'NA','aqs050930005':'NA','aqs051190005':'NA','aqs051410001':'NA','aqs060010010':'NA','aqs060012001':'NA','aqs060130010':'NA','aqs060132007':'NA',
                'aqs060170013':'NA','aqs060170020':'NA','aqs060171002':'NA','aqs060172002':'NA','aqs060190010':'NA','aqs060190243':'NA','aqs060190244':'NA','aqs060210001':'NA','aqs060250007':'NA','aqs060290004':'NA','aqs060291005':'NA','aqs060299000':'NA',
                'aqs060370001':'NA','aqs060370030':'NA','aqs060370031':'NA','aqs060370206':'NA','aqs060371106':'NA','aqs060372401':'NA','aqs060379006':'NA','aqs060390002':'NA','aqs060390003':'NA','aqs060430006':'NA','aqs060431003':'NA','aqs060431004':'NA',
                'aqs060431005':'NA','aqs060510001':'NA','aqs060570007':'NA','aqs060593002':'NA','aqs060596001':'NA','aqs060610810':'NA','aqs060611003':'NA','aqs060611004':'NA','aqs060612001':'NA','aqs060650003':'NA','aqs060650006':'NA','aqs060651002':'NA',
                'aqs060670003':'NA','aqs060670005':'NA','aqs060671001':'NA','aqs060710006':'NA','aqs060710217':'NA','aqs060711001':'NA','aqs060719000':'NA','aqs060719003':'NA','aqs060719006':'NA','aqs060719007':'NA','aqs060719008':'NA','aqs060731009':'NA',
                'aqs060731014':'NA','aqs060750006':'NA','aqs060830002':'NA','aqs060852004':'NA','aqs060852007':'NA','aqs060870005':'NA','aqs060871001':'NA','aqs060871003':'NA','aqs060950006':'NA','aqs060970004':'NA','aqs060990010':'NA','aqs061010004':'NA',
                'aqs061030004':'NA','aqs061071001':'NA','aqs061118001':'NA','aqs080010600':'NA','aqs080017015':'NA','aqs080130007':'NA','aqs080130601':'NA','aqs080137001':'NA','aqs080137002':'NA','aqs080150001':'NA','aqs080190003':'NA','aqs080190004':'NA',
                'aqs080190005':'NA','aqs080310009':'NA','aqs080310026':'NA','aqs080350002':'NA','aqs080350603':'NA','aqs080450013':'NA','aqs080450014':'NA','aqs080450015':'NA','aqs080450017':'NA','aqs080510008':'NA','aqs080590004':'NA','aqs080590012':'NA',
                'aqs080590600':'NA','aqs080590601':'NA','aqs080677002':'NA','aqs080755002':'NA','aqs080770021':'NA','aqs080770022':'NA','aqs080771001':'NA','aqs080771004':'NA','aqs080875001':'NA','aqs080930001':'NA','aqs080970007':'NA','aqs081015004':'NA',
                'aqs081130007':'NA','aqs081130008':'NA','aqs081190003':'NA','aqs081230005':'NA','aqs081230011':'NA','aqs081231001':'NA','aqs081231007':'NA','aqs090010017':'NA','aqs090010113':'NA','aqs090010123':'NA','aqs090011123':'NA','aqs090012004':'NA',
                'aqs090013007':'NA','aqs090031010':'NA','aqs090031123':'NA','aqs090050006':'NA','aqs090052001':'NA','aqs090090123':'NA','aqs090091123':'NA','aqs090093002':'NA','aqs090099002':'NA','aqs090099005':'NA','aqs090110005':'NA','aqs090110008':'NA',
                'aqs090110124':'NA','aqs100010001':'NA','aqs100030018':'NA','aqs100033001':'NA','aqs100052001':'NA','aqs110010050':'NA','aqs120010010':'NA','aqs120010016':'NA','aqs120111001':'NA','aqs120311003':'NA','aqs120710005':'NA','aqs120730003':'NA',
                'aqs120731004':'NA','aqs120810002':'NA','aqs120860026':'NA','aqs120996001':'NA','aqs121110013':'NA','aqs121150010':'NA','aqs130510021':'NA','aqs130570001':'NA','aqs130590002':'NA','aqs130670003':'NA','aqs130730001':'NA','aqs130770002':'NA',
                'aqs130893001':'NA','aqs130970002':'NA','aqs130970004':'NA','aqs131110094':'NA','aqs131130001':'NA','aqs131210053':'NA','aqs131210055':'NA','aqs131270006':'NA','aqs131350002':'NA','aqs132150008':'NA','aqs132230002':'NA','aqs132450091':'NA',
                'aqs150031001':'NA','aqs160010010':'NA','aqs160010017':'NA','aqs160010019':'NA','aqs160010030':'NA','aqs160050015':'NA','aqs160270007':'NA','aqs160310001':'NA','aqs160390010':'NA','aqs160550008':'NA','aqs170010005':'NA','aqs170310009':'NA',
                'aqs170310025':'NA','aqs170310032':'NA','aqs170310038':'NA','aqs170310039':'NA','aqs170310042':'NA','aqs170310044':'NA','aqs170310062':'NA','aqs170310072':'NA','aqs170311002':'NA','aqs170311003':'NA','aqs170312301':'NA','aqs170313005':'NA',
                'aqs170313103':'NA','aqs170314002':'NA','aqs170314003':'NA','aqs170314006':'NA','aqs170315001':'NA','aqs170315002':'NA','aqs170316002':'NA','aqs170318003':'NA','aqs170430003':'NA','aqs170438002':'NA','aqs170650001':'NA','aqs170890003':'NA',
                'aqs170971003':'NA','aqs170990001':'NA','aqs171050001':'NA','aqs171070001':'NA','aqs171150006':'NA','aqs171192008':'NA','aqs171432001':'NA','aqs171670007':'NA','aqs171670013':'NA','aqs180030002':'NA','aqs180030004':'NA','aqs180030006':'NA',
                'aqs180050004':'NA','aqs180110001':'NA','aqs180150001':'NA','aqs180150002':'NA','aqs180190003':'NA','aqs180190008':'NA','aqs180350010':'NA','aqs180390002':'NA','aqs180390007':'NA','aqs180431004':'NA','aqs180510011':'NA','aqs180550001':'NA',
                'aqs180570004':'NA','aqs180570005':'NA','aqs180570006':'NA','aqs180571001':'NA','aqs180590001':'NA','aqs180590002':'NA','aqs180590003':'NA','aqs180590004':'NA','aqs180630004':'NA','aqs180690001':'NA','aqs180690002':'NA','aqs180710001':'NA',
                'aqs180810001':'NA','aqs180810002':'NA','aqs180830003':'NA','aqs180850001':'NA','aqs180850002':'NA','aqs180890024':'NA','aqs180890030':'NA','aqs180892001':'NA','aqs180892002':'NA','aqs180910005':'NA','aqs180910010':'NA','aqs180930003':'NA',
                'aqs180950009':'NA','aqs180950010':'NA','aqs180970030':'NA','aqs180970050':'NA','aqs180970070':'NA','aqs180970073':'NA','aqs180970087':'NA','aqs181090001':'NA','aqs181090003':'NA','aqs181090005':'NA','aqs181230008':'NA','aqs181230009':'NA',
                'aqs181270021':'NA','aqs181270024':'NA','aqs181270026':'NA','aqs181270903':'NA','aqs181271004':'NA','aqs181271005':'NA','aqs181290003':'NA','aqs181410008':'NA','aqs181410010':'NA','aqs181410015':'NA','aqs181410016':'NA','aqs181411007':'NA',
                'aqs181411008':'NA','aqs181450001':'NA','aqs181570006':'NA','aqs181571001':'NA','aqs181630008':'NA','aqs181630012':'NA','aqs181630021':'NA','aqs181670018':'NA','aqs181670024':'NA','aqs181698001':'NA','aqs181730002':'NA','aqs181730008':'NA',
                'aqs181730009':'NA','aqs181730011':'NA','aqs181750001':'NA','aqs190450021':'NA','aqs190850007':'NA','aqs190851101':'NA','aqs191031001':'NA','aqs191130028':'NA','aqs191130033':'NA','aqs191130040':'NA','aqs191471002':'NA','aqs191530030':'NA',
                'aqs191690011':'NA','aqs191770005':'NA','aqs191810022':'NA','aqs191931006':'NA','aqs191931007':'NA','aqs200870001':'NA','aqs200870002':'NA','aqs200910005':'NA','aqs201030004':'NA','aqs201210001':'NA','aqs201330003':'NA','aqs201450001':'NA',
                'aqs202090017':'NA','aqs210131002':'NA','aqs210190010':'NA','aqs210190015':'NA','aqs210190017':'NA','aqs210210004':'NA','aqs210290004':'NA','aqs210370003':'NA','aqs210371001':'NA','aqs210373002':'NA','aqs210490003':'NA','aqs210491002':'NA',
                'aqs210510005':'NA','aqs210670012':'NA','aqs210710002':'NA','aqs210751001':'NA','aqs210830003':'NA','aqs210930004':'NA','aqs210930005':'NA','aqs211010008':'NA','aqs211010013':'NA','aqs211030001':'NA','aqs211070001':'NA','aqs211110034':'NA',
                'aqs211110041':'NA','aqs211170007':'NA','aqs211231001':'NA','aqs211270003':'NA','aqs211390004':'NA','aqs211410003':'NA','aqs211450027':'NA','aqs211451020':'NA','aqs211670001':'NA','aqs211690001':'NA','aqs211770009':'NA','aqs211830031':'NA',
                'aqs211910002':'NA','aqs211930002':'NA','aqs211950001':'NA','aqs211990002':'NA','aqs212090001':'NA','aqs212270009':'NA','aqs220050003':'NA','aqs220111001':'NA','aqs220190001':'NA','aqs220191001':'NA','aqs220330006':'NA','aqs220510003':'NA',
                'aqs220550004':'NA','aqs220570001':'NA','aqs220770002':'NA','aqs220890100':'NA','aqs220930001':'NA','aqs221010002':'NA','aqs230010014':'NA','aqs230013001':'NA','aqs230031008':'NA','aqs230050014':'NA','aqs230052001':'NA','aqs230053002':'NA',
                'aqs230071001':'NA','aqs230071002':'NA','aqs230090001':'NA','aqs230090301':'NA','aqs230090401':'NA','aqs230092002':'NA','aqs230092003':'NA','aqs230092004':'NA','aqs230112001':'NA','aqs230112005':'NA','aqs230130002':'NA','aqs230130003':'NA',
                'aqs230150002':'NA','aqs230171006':'NA','aqs230194005':'NA','aqs230194006':'NA','aqs230210002':'NA','aqs230210003':'NA','aqs230230003':'NA','aqs230230004':'NA','aqs230230005':'NA','aqs230252003':'NA','aqs230270005':'NA','aqs230290016':'NA',
                'aqs230310037':'NA','aqs230310039':'NA','aqs230312002':'NA','aqs230312003':'NA','aqs230313002':'NA','aqs230314753':'NA','aqs230316001':'NA','aqs230318001':'NA','aqs239010001':'NA','aqs240010006':'NA','aqs240030001':'NA','aqs240031003':'NA',
                'aqs240032002':'NA','aqs240050003':'NA','aqs240051007':'NA','aqs240090010':'NA','aqs240150080':'NA','aqs240210037':'NA','aqs240250080':'NA','aqs240259001':'NA','aqs240270005':'NA','aqs240310005':'NA','aqs240311001':'NA','aqs240330004':'NA',
                'aqs240338001':'NA','aqs240430004':'NA','aqs240471001':'NA','aqs245100004':'NA','aqs245100011':'NA','aqs245100019':'NA','aqs245100036':'NA','aqs245100040':'NA','aqs245100050':'NA','aqs245100051':'NA','aqs245100053':'NA','aqs245100054':'NA',
                'aqs250030005':'NA','aqs250030007':'NA','aqs250034002':'NA','aqs250050002':'NA','aqs250050003':'NA','aqs250050004':'NA','aqs250051001':'NA','aqs250051002':'NA','aqs250051006':'NA','aqs250090005':'NA','aqs250091002':'NA','aqs250092001':'NA',
                'aqs250093102':'NA','aqs250094001':'NA','aqs250094003':'NA','aqs250114001':'NA','aqs250130002':'NA','aqs250130014':'NA','aqs250150002':'NA','aqs250150103':'NA','aqs250155001':'NA','aqs250171005':'NA','aqs250171102':'NA','aqs250171401':'NA',
                'aqs250171801':'NA','aqs250171901':'NA','aqs250172003':'NA','aqs250173003':'NA','aqs250174003':'NA','aqs250176001':'NA','aqs250211001':'NA','aqs250212002':'NA','aqs250230005':'NA','aqs250232001':'NA','aqs250250021':'NA','aqs250250041':'NA',
                'aqs250270012':'NA','aqs250270019':'NA','aqs250271001':'NA','aqs260050003':'NA','aqs260170006':'NA','aqs260190002':'NA','aqs260190003':'NA','aqs260210014':'NA','aqs260270001':'NA','aqs260270003':'NA','aqs260330901':'NA','aqs260370001':'NA',
                'aqs260410009':'NA','aqs260410912':'NA','aqs260490011':'NA','aqs260490021':'NA','aqs260550903':'NA','aqs260610101':'NA','aqs260630006':'NA','aqs260630007':'NA','aqs260630904':'NA','aqs260651002':'NA','aqs260770008':'NA','aqs260810022':'NA',
                'aqs260811002':'NA','aqs260812001':'NA','aqs260890001':'NA','aqs260910007':'NA','aqs260990009':'NA','aqs260991003':'NA','aqs261010922':'NA','aqs261030930':'NA','aqs261050006':'NA','aqs261050007':'NA','aqs261070901':'NA','aqs261150037':'NA',
                'aqs261170002':'NA','aqs261210006':'NA','aqs261210022':'NA','aqs261210038':'NA','aqs261210039':'NA','aqs261250001':'NA','aqs261250902':'NA','aqs261251002':'NA','aqs261270001':'NA','aqs261390005':'NA','aqs261390007':'NA','aqs261430901':'NA',
                'aqs261452001':'NA','aqs261470003':'NA','aqs261470005':'NA','aqs261470030':'NA','aqs261590901':'NA','aqs261610005':'NA','aqs261610007':'NA','aqs261610008':'NA','aqs261630015':'NA','aqs261630025':'NA','aqs261630062':'NA','aqs270030002':'NA',
                'aqs270052013':'NA','aqs270076301':'NA','aqs270130002':'NA','aqs270353203':'NA','aqs270353204':'NA','aqs270530962':'NA','aqs270750017':'NA','aqs270834210':'NA','aqs271090016':'NA','aqs271095008':'NA','aqs271230001':'NA','aqs271371011':'NA',
                'aqs271377550':'NA','aqs271390505':'NA','aqs271410008':'NA','aqs271453052':'NA','aqs271630027':'NA','aqs271636016':'NA','aqs271713201':'NA','aqs280010734':'NA','aqs280450002':'NA','aqs280490010':'NA','aqs280490019':'NA','aqs280490020':'NA',
                'aqs280590005':'NA','aqs280730002':'NA','aqs280730720':'NA','aqs280750002':'NA','aqs280750733':'NA','aqs280810003':'NA','aqs280810004':'NA','aqs280810738':'NA','aqs280870722':'NA','aqs280890001':'NA','aqs281490004':'NA','aqs281490721':'NA',
                'aqs281551001':'NA','aqs290030001':'NA','aqs290050001':'NA','aqs290270002':'NA','aqs290370002':'NA','aqs290470004':'NA','aqs290470025':'NA','aqs290770014':'NA','aqs290770026':'NA','aqs290770042':'NA','aqs290950036':'NA','aqs290970003':'NA',
                'aqs290970004':'NA','aqs290990012':'NA','aqs290990019':'NA','aqs291650023':'NA','aqs291830008':'NA','aqs291890001':'NA','aqs291890004':'NA','aqs291890006':'NA','aqs291893001':'NA','aqs291895001':'NA','aqs291897001':'NA','aqs291897002':'NA',
                'aqs291897003':'NA','aqs292130004':'NA','aqs295100007':'NA','aqs295100061':'NA','aqs295100062':'NA','aqs295100063':'NA','aqs295100064':'NA','aqs295100072':'NA','aqs295100080':'NA','aqs295100086':'NA','aqs300130302':'NA','aqs300230007':'NA',
                'aqs300351001':'NA','aqs300630036':'NA','aqs300930018':'NA','aqs301110059':'NA','aqs301110073':'NA','aqs310410001':'NA','aqs310430001':'NA','aqs310550028':'NA','aqs311570005':'NA','aqs311651001':'NA','aqs320031001':'NA','aqs320031007':'NA',
                'aqs320037771':'NA','aqs320037772':'NA','aqs320037773':'NA','aqs320037774':'NA','aqs320037775':'NA','aqs320037776':'NA','aqs320037777':'NA','aqs320037778':'NA','aqs320037779':'NA','aqs320050003':'NA','aqs320190006':'NA','aqs325100001':'NA',
                'aqs325100002':'NA','aqs325100020':'NA','aqs330012003':'NA','aqs330012004':'NA','aqs330031002':'NA','aqs330050006':'NA','aqs330070014':'NA','aqs330074002':'NA','aqs330074003':'NA','aqs330090008':'NA','aqs330110012':'NA','aqs330110016':'NA',
                'aqs330110019':'NA','aqs330110020':'NA','aqs330111005':'NA','aqs330111010':'NA','aqs330111011':'NA','aqs330130006':'NA','aqs330130007':'NA','aqs330131004':'NA','aqs330131007':'NA','aqs330150009':'NA','aqs330150012':'NA','aqs330150015':'NA',
                'aqs330150016':'NA','aqs330173002':'NA','aqs330190003':'NA','aqs340030001':'NA','aqs340031001':'NA','aqs340053001':'NA','aqs340112002':'NA','aqs340130008':'NA','aqs340130011':'NA','aqs340130016':'NA','aqs340131003':'NA','aqs340170003':'NA',
                'aqs340211002':'NA','aqs340232002':'NA','aqs340250004':'NA','aqs340351001':'NA','aqs340353001':'NA','aqs340390008':'NA','aqs340391002':'NA','aqs340395001':'NA','aqs350010008':'NA','aqs350010017':'NA','aqs350010028':'NA','aqs350130010':'NA',
                'aqs350153001':'NA','aqs350390026':'NA','aqs360050073':'NA','aqs360050080':'NA','aqs360070007':'NA','aqs360470007':'NA','aqs360470011':'NA','aqs360470076':'NA','aqs360550005':'NA','aqs360590005':'NA','aqs360590010':'NA','aqs360594002':'NA',
                'aqs360632006':'NA','aqs360651004':'NA','aqs360670015':'NA','aqs360671011':'NA','aqs360671014':'NA','aqs360810070':'NA','aqs360831001':'NA','aqs361037001':'NA','aqs361130003':'NA','aqs361191002':'NA','aqs361195003':'NA','aqs370030003':'NA',
                'aqs370030004':'NA','aqs370030005':'NA','aqs370110001':'NA','aqs370110002':'NA','aqs370210029':'NA','aqs370210030':'NA','aqs370270003':'NA','aqs370290099':'NA','aqs370370098':'NA','aqs370510001':'NA','aqs370511002':'NA','aqs370511003':'NA',
                'aqs370590002':'NA','aqs370590003':'NA','aqs370590099':'NA','aqs370610002':'NA','aqs370630013':'NA','aqs370630015':'NA','aqs370670004':'NA','aqs370670006':'NA','aqs370670007':'NA','aqs370670022':'NA','aqs370670027':'NA','aqs370670030':'NA',
                'aqs370690001':'NA','aqs370770001':'NA','aqs370810011':'NA','aqs370810013':'NA','aqs370870004':'NA','aqs370870008':'NA','aqs370990005':'NA','aqs371010099':'NA','aqs371090099':'NA','aqs371170099':'NA','aqs371191005':'NA','aqs371310002':'NA',
                'aqs371450099':'NA','aqs371470099':'NA','aqs371510004':'NA','aqs371550099':'NA','aqs371572001':'NA','aqs371590022':'NA','aqs371730002':'NA','aqs371730007':'NA','aqs371790003':'NA','aqs371830015':'NA','aqs371830016':'NA','aqs371830017':'NA',
                'aqs371832001':'NA','aqs371990003':'NA','aqs371990004':'NA','aqs380171002':'NA','aqs380171003':'NA','aqs380570001':'NA','aqs381050003':'NA','aqs390030001':'NA','aqs390030009':'NA','aqs390071001':'NA','aqs390170004':'NA','aqs390170018':'NA',
                'aqs390230003':'NA','aqs390250020':'NA','aqs390250022':'NA','aqs390271001':'NA','aqs390350002':'NA','aqs390350034':'NA','aqs390350064':'NA','aqs390354003':'NA','aqs390355002':'NA','aqs390490004':'NA','aqs390490028':'NA','aqs390490029':'NA',
                'aqs390490037':'NA','aqs390570006':'NA','aqs390610020':'NA','aqs390610035':'NA','aqs390610037':'NA','aqs390730002':'NA','aqs390771001':'NA','aqs390810016':'NA','aqs390810017':'NA','aqs390850001':'NA','aqs390850003':'NA','aqs390850007':'NA',
                'aqs390853002':'NA','aqs390870006':'NA','aqs390870011':'NA','aqs390870012':'NA','aqs390890005':'NA','aqs390911001':'NA','aqs390930017':'NA','aqs390930018':'NA','aqs390931002':'NA','aqs390950006':'NA','aqs390950022':'NA','aqs390950024':'NA',
                'aqs390950027':'NA','aqs390950081':'NA','aqs390970006':'NA','aqs390990007':'NA','aqs390990009':'NA','aqs390990013':'NA','aqs390991001':'NA','aqs391030003':'NA','aqs391130007':'NA','aqs391130019':'NA','aqs391130025':'NA','aqs391130033':'NA',
                'aqs391130037':'NA','aqs391331001':'NA','aqs391430017':'NA','aqs391510016':'NA','aqs391510019':'NA','aqs391510021':'NA','aqs391510022':'NA','aqs391510023':'NA','aqs391511001':'NA','aqs391514005':'NA','aqs391530020':'NA','aqs391550008':'NA',
                'aqs391570003':'NA','aqs391591001':'NA','aqs391650006':'NA','aqs391650007':'NA','aqs391651002':'NA','aqs391670004':'NA','aqs391730003':'NA','aqs400130380':'NA','aqs400690323':'NA','aqs400950311':'NA','aqs400950312':'NA','aqs400970185':'NA',
                'aqs401019019':'NA','aqs401139020':'NA','aqs401259023':'NA','aqs410050004':'NA','aqs410050102':'NA','aqs410051006':'NA','aqs410052001':'NA','aqs410052002':'NA','aqs410053001':'NA','aqs410054001':'NA','aqs410090004':'NA','aqs410170122':'NA',
                'aqs410290010':'NA','aqs410290201':'NA','aqs410390008':'NA','aqs410390060':'NA','aqs410391007':'NA','aqs410430103':'NA','aqs410470004':'NA','aqs410470006':'NA','aqs410597001':'NA','aqs410670005':'NA','aqs410671004':'NA','aqs420030010':'NA',
                'aqs420030067':'NA','aqs420030080':'NA','aqs420030081':'NA','aqs420030088':'NA','aqs420031008':'NA','aqs420070004':'NA','aqs420070501':'NA','aqs420110010':'NA','aqs420190501':'NA','aqs420274000':'NA','aqs420290050':'NA','aqs420450102':'NA',
                'aqs420691100':'NA','aqs420770003':'NA','aqs420792016':'NA','aqs420810402':'NA','aqs420814000':'NA','aqs420890001':'NA','aqs420910101':'NA','aqs420950701':'NA','aqs420958000':'NA','aqs421010048':'NA','aqs421011002':'NA','aqs421250501':'NA',
                'aqs421290006':'NA','aqs421290101':'NA','aqs421330011':'NA','aqs440070012':'NA','aqs440090007':'NA','aqs450130003':'NA','aqs450130090':'NA','aqs450190045':'NA','aqs450210003':'NA','aqs450310002':'NA','aqs450450007':'NA','aqs450450009':'NA',
                'aqs450450016':'NA','aqs450510003':'NA','aqs450511001':'NA','aqs450770003':'NA','aqs450791004':'NA','aqs450911004':'NA','aqs450918001':'NA','aqs450918002':'NA','aqs461030016':'NA','aqs461030018':'NA','aqs470090008':'NA','aqs470430007':'NA',
                'aqs470450003':'NA','aqs470450005':'NA','aqs470470103':'NA','aqs470550002':'NA','aqs470890001':'NA','aqs470931012':'NA','aqs470931013':'NA','aqs470931030':'NA','aqs471050003':'NA','aqs471070003':'NA','aqs471131001':'NA','aqs471131003':'NA',
                'aqs471170001':'NA','aqs471310003':'NA','aqs471450004':'NA','aqs471550102':'NA','aqs471570021':'NA','aqs471570024':'NA','aqs471570032':'NA','aqs471630007':'NA','aqs471630009':'NA','aqs471631001':'NA','aqs471870103':'NA','aqs471870105':'NA',
                'aqs480271045':'NA','aqs480290055':'NA','aqs480850004':'NA','aqs481130052':'NA','aqs481130055':'NA','aqs481210002':'NA','aqs481350002':'NA','aqs481390082':'NA','aqs481710001':'NA','aqs482010027':'NA','aqs482010028':'NA','aqs482011017':'NA',
                'aqs482011041':'NA','aqs482017001':'NA','aqs482330001':'NA','aqs482450010':'NA','aqs482731001':'NA','aqs483390089':'NA','aqs483750010':'NA','aqs483751001':'NA','aqs484530613':'NA','aqs484730001':'NA','aqs490030003':'NA','aqs490050001':'NA',
                'aqs490050002':'NA','aqs490090001':'NA','aqs490110002':'NA','aqs490137011':'NA','aqs490350004':'NA','aqs490350009':'NA','aqs490352004':'NA','aqs490353003':'NA','aqs490353007':'NA','aqs490353008':'NA','aqs490450003':'NA','aqs490470014':'NA',
                'aqs490471001':'NA','aqs490475632':'NA','aqs490477022':'NA','aqs490493001':'NA','aqs490495007':'NA','aqs490495008':'NA','aqs490495010':'NA','aqs490530004':'NA','aqs490530006':'NA','aqs490570001':'NA','aqs490570007':'NA','aqs490571001':'NA',
                'aqs490571002':'NA','aqs490571003':'NA','aqs500070003':'NA','aqs500070006':'NA','aqs500250002':'NA','aqs500272001':'NA','aqs510030001':'NA','aqs510090006':'NA','aqs510130020':'NA','aqs510590030':'NA','aqs510690010':'NA','aqs510850001':'NA',
                'aqs510850003':'NA','aqs510890006':'NA','aqs511071005':'NA','aqs511310001':'NA','aqs511730005':'NA','aqs511790001':'NA','aqs511910002':'NA','aqs515100009':'NA','aqs515100021':'NA','aqs515500008':'NA','aqs516500008':'NA','aqs517000013':'NA',
                'aqs517750001':'NA','aqs518000004':'NA','aqs518000005':'NA','aqs530090005':'NA','aqs530090016':'NA','aqs530091004':'NA','aqs530110009':'NA','aqs530110011':'NA','aqs530150013':'NA','aqs530330010':'NA','aqs530330017':'NA','aqs530330018':'NA',
                'aqs530330088':'NA','aqs530337001':'NA','aqs530337002':'NA','aqs530390003':'NA','aqs530410003':'NA','aqs530410007':'NA','aqs530450005':'NA','aqs530530028':'NA','aqs530531008':'NA','aqs530570011':'NA','aqs530612001':'NA','aqs530630001':'NA',
                'aqs530630021':'NA','aqs530630046':'NA','aqs530670002':'NA','aqs530670005':'NA','aqs530730005':'NA','aqs540030003':'NA','aqs540250002':'NA','aqs540250003':'NA','aqs540390004':'NA','aqs540390010':'NA','aqs540610003':'NA','aqs540690009':'NA',
                'aqs540690010':'NA','aqs550090012':'NA','aqs550090014':'NA','aqs550090024':'NA','aqs550170002':'NA','aqs550250007':'NA','aqs550250022':'NA','aqs550250026':'NA','aqs550250034':'NA','aqs550310003':'NA','aqs550310007':'NA','aqs550330001':'NA',
                'aqs550350014':'NA','aqs550390004':'NA','aqs550390005':'NA','aqs550430007':'NA','aqs550450001':'NA','aqs550550001':'NA','aqs550550002':'NA','aqs550550009':'NA','aqs550590002':'NA','aqs550590021':'NA','aqs550590022':'NA','aqs550590025':'NA',
                'aqs550610001':'NA','aqs550630003':'NA','aqs550630006':'NA','aqs550630008':'NA','aqs550631007':'NA','aqs550710002':'NA','aqs550710003':'NA','aqs550710004':'NA','aqs550710008':'NA','aqs550730007':'NA','aqs550730010':'NA','aqs550730991':'NA',
                'aqs550790010':'NA','aqs550790040':'NA','aqs550790041':'NA','aqs550790042':'NA','aqs550790044':'NA','aqs550790045':'NA','aqs550790048':'NA','aqs550790085':'NA','aqs550791025':'NA','aqs550830002':'NA','aqs550850002':'NA','aqs550870009':'NA',
                'aqs550870010':'NA','aqs550950001':'NA','aqs550990998':'NA','aqs550990999':'NA','aqs551010017':'NA','aqs551010031':'NA','aqs551050009':'NA','aqs551050017':'NA','aqs551050018':'NA','aqs551050020':'NA','aqs551050021':'NA','aqs551050024':'NA',
                'aqs551050030':'NA','aqs551051003':'NA','aqs551051005':'NA','aqs551090001':'NA','aqs551170001':'NA','aqs551170004':'NA','aqs551170005':'NA','aqs551170007':'NA','aqs551250002':'NA','aqs551270002':'NA','aqs551310003':'NA','aqs551310007':'NA',
                'aqs551310009':'NA','aqs551330017':'NA','aqs551330027':'NA','aqs551390007':'NA','aqs551390011':'NA','aqs551410998':'NA','aqs560050800':'NA','aqs560070099':'NA','aqs560071000':'NA','aqs560090008':'NA','aqs560111013':'NA','aqs560130900':'NA',
                'aqs560136001':'NA','aqs560250100':'NA','aqs560252601':'NA','aqs560330004':'NA','aqs560350097':'NA','aqs560370100':'NA','aqs560370870':'NA','aqs720210010':'NA','aqs721370002':'NA','aqs721370003':'NA','aqscc0040002':'NA','aqscc0087004':'NA',
                'aqs300870760':'NA','aqs300870762':'NA',
                #AIRBASE
                'abat0enk1':'EU','abat0ill1':'EU','abat0pil1':'EU','abat0zoe2':'EU','abat10002':'EU','abat2f202':'EU','abat2m121':'EU','abat2sp10':'EU','abat2sp20':'EU','abat2vk26':'EU','abat2wo35':'EU', \
                'abat30103':'EU','abat30202':'EU','abat30302':'EU','abat30403':'EU','abat30502':'EU','abat30603':'EU','abat30701':'EU','abat30801':'EU','abat31102':'EU','abat31204':'EU','abat31502':'EU', \
                'abat31701':'EU','abat31904':'EU','abat32101':'EU','abat4s108':'EU','abat4s165':'EU','abat4s418':'EU','abat4s420':'EU','abat52100':'EU','abat53055':'EU','abat56071':'EU','abat60137':'EU', \
                'abat60151':'EU','abat60157':'EU','abat60185':'EU','abat60190':'EU','abat72538':'EU','abat72705':'EU','abat82801':'EU','abbetn016':'EU','abbetn029':'EU','abbetn040':'EU','abbetn046':'EU', \
                'abbetn051':'EU','abbetn052':'EU','abbetn054':'EU','abbetn063':'EU','abbetn073':'EU','abbetn085':'EU','abbetn093':'EU','abbetn100':'EU','abbetn113':'EU','abbetn121':'EU','abbetn132':'EU', \
                'abbetr222':'EU','abbg0055a':'EU','abch0002r':'EU','abch0003r':'EU','abch0016a':'EU','abch0019a':'EU','abch0020a':'EU','abch0024a':'EU','abch0028a':'EU','abch0031a':'EU','abch0033a':'EU', \
                'abch0037a':'EU','abch0040a':'EU','abch0041a':'EU','abch0045a':'EU','abch0051a':'EU','abcy0002r':'EU','abcz0avyn':'EU','abcz0bkuc':'EU','abcz0bmis':'EU','abcz0chvo':'EU','abcz0ckoc':'EU', \
                'abcz0esvr':'EU','abcz0jkmy':'EU','abcz0jkos':'EU','abcz0kprb':'EU','abcz0lsou':'EU','abcz0mjes':'EU','abcz0ppla':'EU','abcz0pprm':'EU','abcz0sonr':'EU','abcz0tbkr':'EU','abcz0tcer':'EU', \
                'abcz0tstd':'EU','abcz0ulom':'EU','abcz0urvh':'EU','abcz0usnz':'EU','abcz0utus':'EU','abcz0uval':'EU','abcz0zsnv':'EU','abdebb053':'EU','abdebb065':'EU','abdebb066':'EU','abdebe027':'EU', \
                'abdebe032':'EU','abdebe056':'EU','abdebe062':'EU','abdebw004':'EU','abdebw030':'EU','abdebw031':'EU','abdebw087':'EU','abdebw103':'EU','abdeby013':'EU','abdeby047':'EU','abdeby049':'EU', \
                'abdeby072':'EU','abdeby081':'EU','abdeby109':'EU','abdehe023':'EU','abdehe024':'EU','abdehe026':'EU','abdehe028':'EU','abdehe042':'EU','abdehe043':'EU','abdehe046':'EU','abdehe051':'EU', \
                'abdehe052':'EU','abdehe060':'EU','abdemv004':'EU','abdemv012':'EU','abdemv017':'EU','abdeni019':'EU','abdeni031':'EU','abdeni051':'EU','abdeni058':'EU','abdeni059':'EU','abdeni060':'EU', \
                'abdeni063':'EU','abdenw064':'EU','abdenw065':'EU','abdenw068':'EU','abdenw074':'EU','abdenw081':'EU','abderp013':'EU','abderp014':'EU','abderp015':'EU','abderp016':'EU','abderp017':'EU', \
                'abderp028':'EU','abdesh001':'EU','abdesh008':'EU','abdesh011':'EU','abdesh017':'EU','abdesl019':'EU','abdesn049':'EU','abdesn051':'EU','abdesn052':'EU','abdesn074':'EU','abdesn076':'EU', \
                'abdesn079':'EU','abdesn080':'EU','abdest068':'EU','abdest069':'EU','abdest089':'EU','abdest098':'EU','abdeth026':'EU','abdeth027':'EU','abdeth040':'EU','abdeth042':'EU','abdeth061':'EU', \
                'abdeub001':'EU','abdeub005':'EU','abdeub028':'EU','abdeub029':'EU','abdeub030':'EU','abdk0031r':'EU','abdk0041a':'EU','abdk0054a':'EU','abee0009r':'EU','abee0011r':'EU','abee0016a':'EU', \
                'abee0020a':'EU','abes0008r':'EU','abes0010r':'EU','abes0011r':'EU','abes0012r':'EU','abes0013r':'EU','abes0014r':'EU','abes0016r':'EU','abes0094a':'EU','abes0316a':'EU','abes0324a':'EU', \
                'abes0330a':'EU','abes0774a':'EU','abes0813a':'EU','abes1173a':'EU','abes1201a':'EU','abes1217a':'EU','abes1222a':'EU','abes1248a':'EU','abes1285a':'EU','abes1311a':'EU','abes1347a':'EU', \
                'abes1359a':'EU','abes1378a':'EU','abes1379a':'EU','abes1400a':'EU','abes1488a':'EU','abes1489a':'EU','abes1491a':'EU','abes1517a':'EU','abes1518a':'EU','abes1531a':'EU','abes1542a':'EU', \
                'abes1543a':'EU','abes1588a':'EU','abes1599a':'EU','abes1616a':'EU','abes1648a':'EU','abes1654a':'EU','abes1660a':'EU','abes1661a':'EU','abes1670a':'EU','abes1671a':'EU','abes1688a':'EU', \
                'abes1689a':'EU','abes1690a':'EU','abes1746a':'EU','abes1753a':'EU','abes1754a':'EU','abes1779a':'EU','abes1793a':'EU','abes1794a':'EU','abfi00208':'EU','abfi00293':'EU','abfi00349':'EU', \
                'abfi00351':'EU','abfi00352':'EU','abfi00356':'EU','abfi00368':'EU','abfi00372':'EU','abfi00428':'EU','abfi00453':'EU','abfi00460':'EU','abfi00564':'EU','abfr02023':'EU','abfr02024':'EU', \
                'abfr03027':'EU','abfr03031':'EU','abfr04038':'EU','abfr04066':'EU','abfr04142':'EU','abfr04158':'EU','abfr04322':'EU','abfr05053':'EU','abfr07020':'EU','abfr07022':'EU','abfr08204':'EU', \
                'abfr08209':'EU','abfr09022':'EU','abfr12020':'EU','abfr12029':'EU','abfr12031':'EU','abfr13011':'EU','abfr14008':'EU','abfr15001':'EU','abfr15012':'EU','abfr16017':'EU','abfr16031':'EU', \
                'abfr18010':'EU','abfr18026':'EU','abfr18037':'EU','abfr18039':'EU','abfr18045':'EU','abfr19001':'EU','abfr20049':'EU','abfr22017':'EU','abfr23124':'EU','abfr24023':'EU','abfr25045':'EU', \
                'abfr26012':'EU','abfr30033':'EU','abfr31008':'EU','abfr32008':'EU','abfr34003':'EU','abfr34043':'EU','abfr34054':'EU','abfr35012':'EU','abfr36005':'EU','abgb0002r':'EU','abgb0006r':'EU', \
                'abgb0013r':'EU','abgb0014r':'EU','abgb0015r':'EU','abgb0031r':'EU','abgb0033r':'EU','abgb0035r':'EU','abgb0036r':'EU','abgb0037r':'EU','abgb0038r':'EU','abgb0039r':'EU','abgb0051a':'EU', \
                'abgr0110r':'EU','abhu0002r':'EU','abhu0040a':'EU','abie0001r':'EU','abie0031r':'EU','abie0090a':'EU','abie0091a':'EU','abie0102a':'EU','abie0111a':'EU','abis0007a':'EU','abit0741a':'EU', \
                'abit0824a':'EU','abit0842a':'EU','abit0952a':'EU','abit0988a':'EU','abit0989a':'EU','abit0992a':'EU','abit1121a':'EU','abit1174a':'EU','abit1179a':'EU','abit1188a':'EU','abit1233a':'EU', \
                'abit1236a':'EU','abit1371a':'EU','abit1373a':'EU','abit1451a':'EU','abit1474a':'EU','abit1519a':'EU','abit1665a':'EU','abit1678a':'EU','abit1685a':'EU','abit1695a':'EU','abli0002a':'EU', \
                'ablt00051':'EU','ablt00053':'EU','ablt00054':'EU','ablu0102a':'EU','ablu0103a':'EU','ablu0104a':'EU','ablu0105a':'EU','abnl00107':'EU','abnl00131':'EU','abnl00227':'EU','abnl00230':'EU', \
                'abnl00235':'EU','abnl00301':'EU','abnl00318':'EU','abnl00437':'EU','abnl00444':'EU','abnl00538':'EU','abnl00620':'EU','abnl00631':'EU','abnl00633':'EU','abnl00641':'EU','abnl00722':'EU', \
                'abnl00738':'EU','abnl00807':'EU','abnl00818':'EU','abnl00918':'EU','abnl00929':'EU','abnl00934':'EU','abno0001r':'EU','abno0015r':'EU','abno0039r':'EU','abno0042r':'EU','abno0043r':'EU', \
                'abno0052r':'EU','abno0055r':'EU','abno0056r':'EU','abpl0002r':'EU','abpl0004r':'EU','abpl0005r':'EU','abpl0014a':'EU','abpl0028a':'EU','abpl0077a':'EU','abpl0094a':'EU','abpl0105a':'EU', \
                'abpl0121a':'EU','abpl0128a':'EU','abpl0150a':'EU','abpl0182a':'EU','abpl0211a':'EU','abpl0243a':'EU','abpl0246a':'EU','abpl0247a':'EU','abpt01044':'EU','abpt02019':'EU','abpt02020':'EU', \
                'abpt03091':'EU','abpt03092':'EU','abpt03096':'EU','abpt04003':'EU','abpt04006':'EU','abro0072a':'EU','abse0005r':'EU','abse0011r':'EU','abse0012r':'EU','abse0013r':'EU','abse0014r':'EU', \
                'abse0026a':'EU','abse0035r':'EU','abse0039r':'EU','abse0054a':'EU','absi0008r':'EU','absi0033a':'EU','absk0004r':'EU','absk0006r':'EU','absk0007r':'EU','absk0041a':'EU','abfi00357':'EU', \
                'abat0ach1':'EU','abat11002':'EU','abat30102':'EU','abat32201':'EU','abat60147':'EU','abat82707':'EU','abdehe002':'EU','abdehe027':'EU','abdehe034':'EU','abdemv001':'EU','abdenw063':'EU', \
                'abdenw075':'EU','abdenw076':'EU','abdesl008':'EU','abdeub002':'EU','abdeub006':'EU','abdeub007':'EU','abdeub012':'EU','abdeub013':'EU','abdeub019':'EU','abdeub022':'EU','abdeub023':'EU', \
                'abdeub024':'EU','abdeub025':'EU','abdeub026':'EU','abdeub027':'EU','abdeub031':'EU','abdeub032':'EU','abdeub033':'EU','abdeub035':'EU','abfi00350':'EU','abnl00724':'EU','abnl00733':'EU', \
                'abnl00913':'EU','abnl00928':'EU','abdesh013':'EU','abdesh014':'EU','abdeub021':'EU','abnl00232':'EU','abat60184':'EU','abcz0uvse':'EU','abdesn057':'EU','abat30992':'EU','abat31496':'EU', \
                'abcz0mbup':'EU','abdehe039':'EU','abdest070':'EU','abgb0044r':'EU','abgb0617a':'EU','abdeub034':'EU','abes1200a':'EU','abes1436a':'EU','abgb0045r':'EU','abpl0024a':'EU','abdeub038':'EU', \
                'abdeub039':'EU','abgb0043r':'EU','ablv00010':'EU','abno0041r':'EU','abno0045r':'EU','abpl0026a':'EU','abdenw093':'EU','abdeub017':'EU','abdeub040':'EU','abes1216a':'EU','abfr04324':'EU', \
                'abfr04328':'EU','abfr16201':'EU','abfr19009':'EU','abfr22010':'EU','abit0774a':'EU','abpl0030a':'EU','abgb0745a':'EU','abes1605a':'EU','abes1607a':'EU','abes1614a':'EU','abgb0754a':'EU', \
                'abes1538a':'EU','absk0008a':'EU','absk0009a':'EU','abat30407':'EU','abit1522a':'EU','abro0060a':'EU','abbg0056a':'EU','abes1662a':'EU','abes1802a':'EU','abes1805a':'EU','abes1806a':'EU', \
                'abes1808a':'EU','abes1810a':'EU','abes1811a':'EU','abes1813a':'EU','abfr10041':'EU','abfr23175':'EU','abfr24030':'EU','abfr26019':'EU', \
                'abfr34034':'EU','abit1464a':'EU','abit1681a':'EU','abit1736a':'EU','ablt00052':'EU','abpt02021':'EU','abro0082a':'EU','abes1827a':'EU','abes1851a':'EU','abfr07038':'EU','abfr12046':'EU', \
                'abfr25049':'EU','abfr38011':'O','abfr38012':'O','abfr41007':'EU','abit0659a':'EU','abit1061a':'EU','abit1149a':'EU','abit1214a':'EU','abit1288a':'EU','abit1340a':'EU','abit1553a':'EU', \
                'abit1619a':'EU','abit1848a':'EU','abit1863a':'EU','abit1865a':'EU','abmt00007':'EU','abpt01047':'EU','abse0032r':'EU','abse0066a':'EU', \
                'abat0hbg1':'EU','abat2m226':'EU', 'abat2sv24':'EU','abat30402':'EU','abat30602':'EU','abat30901':'EU','abat31703':'EU','abat31902':'EU','abat31903':'EU','abat31905':'EU','abat31906':'EU','abat32604':'EU','abat4s199':'EU','abat4s407':'EU', \
                'abat60103':'EU','abat60104':'EU','abat60114':'EU','abat60115':'EU','abat60119':'EU','abat60144':'EU','abat60154':'EU','abat60180':'EU','abat60197':'EU','abba0001g':'EU','abbelhr01':'EU','abbelld02':'EU', \
                'abbetn015':'EU','abbetn027':'EU','abbetn060':'EU','abbg0019a':'EU','abbg0026a':'EU','abbg0038a':'EU','abbg0057a':'EU','abbg0058a':'EU','abbg0071a':'EU','abcz0htrm':'EU','abcz0jtre':'EU','abcz0kchm':'EU', \
                'abcz0knal':'EU','abcz0lbli':'EU','abcz0lclm':'EU','abcz0lfru':'EU','abcz0ljnm':'EU','abcz0molj':'EU','abcz0molo':'EU','abcz0ppls':'EU','abcz0sdub':'EU','abcz0skls':'EU','abcz0skry':'EU','abcz0tbom':'EU', \
                'abcz0thar':'EU','abcz0tora':'EU','abcz0torv':'EU','abcz0tovk':'EU','abcz0tver':'EU','abcz0uchm':'EU','abcz0ukru':'EU','abcz0umed':'EU','abcz0utpm':'EU','abdebb010':'EU','abdebb016':'EU','abdebb037':'EU', \
                'abdebb051':'EU','abdeby030':'EU','abdeby067':'EU','abdeby069':'EU','abdeby0691':'EU','abdeby092':'EU','abdeby124':'EU','abdehe048':'EU','abdehe050':'EU','abdehh063':'EU','abdemv024':'EU','abdeni014':'EU', \
                'abdeni030':'EU','abdeni0311':'EU','abdeni077':'EU','abdenw001':'EU','abdenw005':'EU','abdenw032':'EU','abdenw033':'EU','abdest104':'EU','abdeth017':'EU','abdeth0261':'EU','abdeth037':'EU','abdeub0011':'EU', \
                'abdeub0021':'EU','abdeub0051':'EU','abdeub0061':'EU','abdeub0071':'EU','abdeub008':'EU','abdeub0121':'EU','abdeub0131':'EU','abdeub014':'EU','abdeub0171':'EU','abdeub018':'EU','abdeub036':'EU','abdeub037':'EU', \
                'abdeub041':'EU','abdeub042':'EU','abdk0012r':'EU','abdk0031r1':'EU','abdk0048a':'EU','abee0009r1':'EU','abes0001r':'EU','abes0005r':'EU','abes0006r':'EU','abes0017r':'EU','abes0296a':'EU','abes0297a':'EU', \
                'abes1716a':'EU','abes1774a':'EU','abes1778a':'EU','abes1831a':'EU','abes1853a':'EU','abes1854a':'EU','abes1878a':'EU','abes1882a':'EU','abes1898a':'EU','abes1913a':'EU','abes1917a':'EU','abes1923a':'EU', \
                'abes1987a':'EU','abes1990a':'EU','abes1996a':'EU','abes2014a':'EU','abes2018a':'EU','abfi00424':'EU','abfi00582':'EU','abfr01004':'EU','abfr01014':'EU','abfr02026':'EU','abfr05087':'EU','abfr06018':'EU', \
                'abfr06134':'EU','abfr08401':'EU','abfr09302':'EU','abfr10007':'EU','abfr10014':'EU','abfr10029':'EU','abfr12028':'EU','abfr12052':'EU','abfr12053':'EU','abfr14042':'EU','abfr14051':'EU','abfr15002':'EU', \
                'abfr15018':'EU','abfr16054':'EU','abfr16056':'EU','abfr16058':'EU','abfr16302':'EU','abfr17016':'EU','abfr18049':'EU','abfr20039':'EU','abfr20044':'EU','abfr20068':'EU','abfr21004':'EU','abfr21031':'EU', \
                'abfr22004':'EU','abfr22016':'EU','abfr22022':'EU','abfr23096':'EU','abfr23097':'EU','abfr23156':'EU','abfr23177':'EU','abfr24033':'EU','abfr24039':'EU','abfr27001':'EU','abfr27101':'EU','abfr30029':'EU', \
                'abfr30030':'EU','abfr30037':'EU','abfr33302':'EU','abfr35001':'EU','abfr350011':'EU','abfr36019':'EU','abfr39008':'O','abfr41005':'EU','abfr41024':'EU','abgb0737a':'EU','abgb0838a':'EU','abgb0957a':'EU', \
                'abgb0998a':'EU','abgb0999a':'EU','abgb1017a':'EU','abgr0033a':'EU','abgr0120a':'EU','abie0113a':'EU','abie0115a':'EU','abie0120a':'EU','abie0122a':'EU','abie0124a':'EU','abie0130a':'EU','abie0134a':'EU', \
                'abie0138a':'EU','abie0139a':'EU','abie0144a':'EU','abie0146a':'EU','abie0147a':'EU','abie0148a':'EU','abis0016a':'EU','abis0020a':'EU','abit0267a':'EU','abit0440a':'EU','abit0558a':'EU','abit0692a':'EU', \
                'abit0709a':'EU','abit0838a':'EU','abit0864a':'EU','abit0870a':'EU','abit0982a':'EU','abit0988a1':'EU','abit1085a':'EU','abit1184a':'EU','abit1210a':'EU','abit1212a':'EU','abit1220a':'EU','abit1222a':'EU', \
                'abit1225a':'EU','abit1292a':'EU','abit1294a':'EU','abit1307a':'EU','abit1312a':'EU','abit1418a':'EU','abit1459a':'EU','abit1474a1':'EU','abit1542a':'EU','abit1543a':'EU','abit1548a':'EU','abit1582a':'EU', \
                'abit1595a':'EU','abit1596a':'EU','abit1646a':'EU','abit1663a':'EU','abit1666a':'EU','abit1725a':'EU','abit1729a':'EU','abit1773a':'EU','abit1806a':'EU','abit1823a':'EU','abit1831a':'EU','abit1832a':'EU', \
                'abit1861a':'EU','abit1868a':'EU','abit1870a':'EU','abit1898a':'EU','abit1902a':'EU','abit1904a':'EU','abit1911a':'EU','abit1912a':'EU','abit1914a':'EU','abit1915a':'EU','abit1917a':'EU','abit1919a':'EU', \
                'abit1921a':'EU','abit1924a':'EU','abit1925a':'EU','abit1942a':'EU','abit1944a':'EU','abit1948a':'EU','abit1960a':'EU','abit1998a':'EU','abit2011a':'EU','abit2013a':'EU','abit2014a':'EU','abit2023a':'EU', \
                'abit2027a':'EU','abit2032a':'EU','abit2054a':'EU','abit2060a':'EU','abit2063a':'EU','abit2069a':'EU','abit2071a':'EU','abit2072a':'EU','abit2074a':'EU','abit2076a':'EU','abli0001a':'EU','ablt00021':'EU', \
                'ablu0099a':'EU','ablv000o1':'EU','ablv00ng1':'EU','abnl001071':'EU','abnl001311':'EU','abnl002271':'EU','abnl002351':'EU','abnl00245':'EU','abnl00246':'EU','abnl00247':'EU','abnl003011':'EU','abnl003181':'EU', \
                'abnl004371':'EU','abnl005381':'EU','abnl006201':'EU','abnl006311':'EU','abnl006331':'EU','abnl00644':'EU','abnl00703':'EU','abnl007221':'EU','abnl007241':'EU','abnl008071':'EU','abnl009181':'EU','abnl009281':'EU', \
                'abnl009341':'EU','abno0063a':'EU','abpl0027a':'EU','abpl0043a':'EU','abpl0062a':'EU','abpl0068a':'EU','abpl0115a':'EU','abpl0157a':'EU','abpl0157a1':'EU','abpl0191a':'EU','abpl0198a':'EU','abpl0212a':'EU', \
                'abpl0222a':'EU','abpl0273a':'EU','abpl0276a':'EU','abpl0313a':'EU','abpl0314a':'EU','abpl0317a':'EU','abpl0321a':'EU','abpl0349a':'EU','abpl0437a':'EU','abpl0450a':'EU','abpl0468a':'EU','abpl0469a':'EU', \
                'abpl0470a':'EU','abpl0473a':'EU','abpl0503a':'EU','abpl0505a':'EU','abpl0506a':'EU','abpl0525a':'EU','abpl0552a':'EU','abpl0560a':'EU','abpl0561a':'EU','abpl0568a':'EU','abpl0573a':'EU','abpl0581a':'EU', \
                'abpt01051':'EU','abpt01054':'EU','abpt02022':'EU','abpt03099':'EU','abpt03101':'EU','abpt03102':'EU','abpt04002':'EU','abpt05012':'EU','abpt07001':'O','abro0008r':'EU','abro0077a':'EU','abro0087a':'EU', \
                'abro0126a':'EU','abro0138a':'EU','abro0153a':'EU','abro0191a':'EU','abro0198a':'EU','abro0210a':'EU','abrs0005r':'EU','abrs0031a':'EU','absi0046a':'EU','absk0006a':'EU','absk0007a':'EU','absk0010a':'EU', \
                'absk0011a':'EU','absk0013a':'EU','absk0028a':'EU','absk0030a':'EU','absk0031a':'EU','absk0051a':'EU','abdebw0311':'EU','abdest0701':'EU','abes1662a1':'EU','abit1019a':'EU','abpl0397a':'EU','abro0212a':'EU', \
                'abes1540a':'EU','abhu0017a':'EU','abhu0018a':'EU','abhu0019a':'EU','ablv000101':'EU','ablv00016':'EU','abno0002r':'EU','abno0008r':'EU','abpl0067a':'EU','abpl0081a':'EU','abpl0084a':'EU','abpl0088a':'EU', \
                'abpl0090a':'EU','abpl0103a':'EU','abpl0116a':'EU','abpl0131a':'EU','abpl0153a':'EU','abpl0154a':'EU','abpl0155a':'EU','abpl0158a':'EU','abpl0160a':'EU','abpl0161a':'EU','abpl0163a':'EU','abpl0185a':'EU','abpl0201a':'EU','abpl0208a':'EU','abpl0229a':'EU','abpl0230a':'EU', \
                'abpl0232a':'EU','abpl0233a':'EU','abpl0255a':'EU','abpl0259a':'EU','abpl0277a':'EU','abpl0281a':'EU','abpl0282a':'EU','abpl0287a':'EU','abpl0323a':'EU','abpl0324a':'EU','abpl0332a':'EU','abpl0338a':'EU', \
                'abpl0340a':'EU','abpl0341a':'EU','abpl0342a':'EU','abpl0344a':'EU','abpl0345a':'EU','abpl0347a':'EU','abpl0351a':'EU','abpl0363a':'EU','abpl0366a':'EU','abpl0370a':'EU','abpl0376a':'EU','abpl0378a':'EU', \
                'abpl0379a':'EU','abpl0380a':'EU','abpl0384a':'EU','abpl0389a':'EU','abpl0390a':'EU','abpl0391a':'EU','abpl0392a':'EU','abpl0393a':'EU','abpl0394a':'EU','abpl0418a':'EU','abpl0427a':'EU','abpl0428a':'EU', \
                'abpl0430a':'EU','abpl0433a':'EU','abpl0435a':'EU','abpl0439a':'EU','abpl0440a':'EU','abpl0442a':'EU','abpl0451a':'EU','abpl0452a':'EU','abpl0453a':'EU','abpl0454a':'EU','abpl0456a':'EU','abpl0459a':'EU', \
                'abpl0460a':'EU','abpl0462a':'EU','abpl0464a':'EU','abpl0471a':'EU','abpl0511a':'EU','abro0001a':'EU','abro0006a':'EU','abro0034a':'EU','abro0038a':'EU','abrs0010a':'EU','abrs0011a':'EU','abrs0013a':'EU', \
                'abrs0019a':'EU','abrs0034a':'EU','abrs0040a':'EU','abrs0046a':'EU','abse0002r':'EU','abse0008r':'EU','abse0010a':'EU','abse0018a':'EU','abse0061a':'EU', \
                'abrs0044a':'EU','abes1774a1':'EU','abdeby0690':'EU','abdeth0260':'EU','abdk0031r0':'EU','abee0009r0':'EU','abit1474a0':'EU','abpl0157a0':'EU','abdebw0310':'EU','abdest0700':'EU','abes1662a0':'EU','abes1774a0':'EU', \
                'abfr27006':'EU','ablv000100':'EU','abfr08025':'EU','abfr20204':'EU','abfr36010':'EU','abes1793a1':'EU','abat10007':'EU','abbetn041':'EU','abdeby026':'EU','abdeby093':'EU','abdeni0310':'EU','abdest017':'EU','abdest018':'EU','abdest106':'EU','abdeub0010':'EU','abdeub0020':'EU','abdeub0050':'EU', \
                'abdeub0060':'EU','abdeub0070':'EU','abdeub0120':'EU','abdeub0130':'EU','abdeub0170':'EU','abes2001a':'EU','abfi00354':'EU','abfi00426':'EU','abfr03085':'EU','abfr03088':'EU','abfr05088':'EU','abfr06133':'EU', \
                'abfr08618':'EU','abfr10004':'EU','abfr10132':'EU','abfr1227a':'EU','abfr13001':'EU','abfr19008':'EU','abfr21050':'EU','abfr23102':'EU','abfr23230':'EU','abfr24001':'EU','abfr29437':'EU','abfr34026':'EU', \
                'abfr34038':'EU','abgb0039r0':'EU','abgb0041r':'EU','abgb0048r':'EU','abgb0794a':'EU','abgb0881a':'EU','abie0092a':'EU','abie0109a':'EU','abit0988a0':'EU','abit1121a0':'EU','abit1726a':'EU','abit1758a':'EU', \
                'abit1759a':'EU','abit1982a':'EU','abit2022a':'EU','ablt000510':'EU','ablt000530':'EU','ablt000540':'EU','ablt0052a':'EU','abmt00001':'EU','abnl001070':'EU','abnl001310':'EU','abnl002270':'EU','abnl002350':'EU', \
                'abnl003010':'EU','abnl003180':'EU','abnl004370':'EU','abnl005380':'EU','abnl006200':'EU','abnl006310':'EU','abnl006330':'EU','abnl007220':'EU','abnl007240':'EU','abnl008070':'EU','abnl009180':'EU','abnl009280':'EU', \
                'abnl009340':'EU','abno0056r0':'EU','abpl0383a':'EU','abpt040020':'EU','absi0053a':'EU','absk0010r':'EU','abat0zil1':'EU','abat4s184':'EU','abat4s406':'EU','abat4s412':'EU','abat4s414':'EU','abat60105':'EU','abat60106':'EU', \
                'abat60127':'EU','abat60139':'EU','abat60141':'EU','abat60145':'EU','abat60160':'EU','abba0035a':'EU','abcz0tctn':'EU','abcz0tfmi':'EU','abcz0tozr':'EU','abcz0ttrk':'EU','abcz0udcm':'EU','abdeby009':'EU','abdehh054':'EU', \
                'abdenw066':'EU','abderp042':'EU','abdest006':'EU','abdest074':'EU','abee0021a':'EU','abee0022a':'EU','abes1829a':'EU','abfi00066':'EU','abfi00431':'EU','abfr02005':'EU','abfr02037':'EU','abfr03014':'EU','abfr03038':'EU', \
                'abfr03047':'EU','abfr03062':'EU','abfr05070':'EU','abfr05086':'EU','abfr07002':'EU','abfr07035':'EU','abfr07037':'EU','abfr08003':'EU','abfr10016':'EU','abfr10026':'EU','abfr12003':'EU','abfr13007':'EU','abfr14003':'EU', \
                'abfr14006':'EU','abfr14012':'EU','abfr14022':'EU','abfr15007':'EU','abfr15045':'EU','abfr16057':'EU','abfr19015':'EU','abfr20061':'EU','abfr23174':'EU','abfr23179':'EU','abfr23181':'EU','abfr23182':'EU','abfr23188':'EU', \
                'abfr25048':'EU','abfr28003':'EU','abfr28022':'EU','abfr32014':'EU','abfr33220':'EU','abfr33305':'EU','abfr34046':'EU','abfr37003':'EU','abfr39007':'EU','abfr39009':'EU','abfr41001':'EU','abfr41002':'EU','abfr41017':'EU', \
                'abgb0920a':'EU','abgb0929a':'EU','abgb0962a':'EU','abie0132a':'EU','abie0133a':'EU','abit0441a':'EU','abit0447a':'EU','abit0663a':'EU','abit0684a':'EU','abit0782a':'EU','abit0816a':'EU','abit0862a':'EU','abit0906a':'EU', \
                'abit1042a':'EU','abit1065a':'EU','abit1156a':'EU','abit1187a':'EU','abit1238a':'EU','abit1245a':'EU','abit1328a':'EU','abit1425a':'EU','abit1468a':'EU','abit1557a':'EU','abit1571a':'EU','abit1601a':'EU','abit1670a':'EU', \
                'abit1762a':'EU','abit1845a':'EU','abit1997a':'EU','ablt00001':'EU','ablt00012':'EU','ablt00044':'EU','abmk0030a':'EU','abno0075a':'EU','abno0080a':'EU','abpl0033a':'EU','abpl0039a':'EU','abpl0046a':'EU','abpl0051a':'EU', \
                'abpl0119a':'EU','abpl0122a':'EU','abpl0126a':'EU','abpl0134a':'EU','abpl0136a':'EU','abpl0144a':'EU','abpl0148a':'EU','abpl0151a':'EU','abpl0152a':'EU','abpl0162a':'EU','abpl0187a':'EU','abpl0210a':'EU','abpl0221a':'EU', \
                'abpl0238a':'EU','abpl0256a':'EU','abpl0306a':'EU','abpl0310a':'EU','abpl0311a':'EU','abpl0360a':'EU','abro0074a':'EU','abro0075a':'EU','abro0079a':'EU','abro0084a':'EU','abro0086a':'EU','abro0206a':'EU','abse0030a':'EU', \
                'abse0080a':'EU','absi0054a':'EU','absk0004a':'EU','absk0015a':'EU','absk0027a':'EU','absk0047a':'EU','absk0134a':'EU','abch0036a':'EU','abno0066a':'EU','abro0042a':'EU','abse0053a':'EU','abse0055a':'EU', \
                'abat0sto1':'EU','abat0vor1':'EU','abat2f101':'EU','abat52055':'EU','abba0034a':'EU','abcz0cchu':'EU','abcz0hkry':'EU','abcz0hser':'EU','abdest039':'EU','abes0007r':'EU','abes0009r':'EU','abes0015r':'EU','abes1310a':'EU', \
                'abes1348a':'EU','abes1441a':'EU','abes1606a':'EU','abes1608a':'EU','abes1883a':'EU','abfr07015':'EU','abfr07023':'EU','abfr07029':'EU','abfr07031':'EU','abfr11015':'EU','abfr15111':'EU','abfr16032':'EU','abfr16041':'EU', \
                'abfr24014':'EU','abfr30028':'EU','abfr31027':'EU','abfr38014':'EU','abit0979a':'EU','abit1672a':'EU','abit1812a':'EU','abit1842a':'EU','abpl0058a':'EU','abpl0186a':'EU','absk0040a':'EU','absk0042a':'EU','absk0050a':'EU', \
                'abat60156':'EU','abat80503':'EU','abbg0070a':'EU','abch0004r':'EU','abch0005r':'EU','abdeub004':'EU','abit1791a':'EU','abpt01048':'EU','abfi00586':'EU','abfr24032':'EU','abgb0995a':'EU','absk0263a':'EU','abfr02015':'EU','abmk0042a':'EU', \
                'abgb0040r':'EU','abat0sig1':'EU','abat10001':'EU','abat10003':'EU','abat2ka11':'EU','abat2ka41':'EU','abat2sp18':'EU','abat2vi12':'EU','abat2vl52':'EU','abat2wo15':'EU','abat30065':'EU','abat30101':'EU','abat30201':'EU','abat30301':'EU','abat30601':'EU', \
                'abat30902':'EU','abat31301':'EU','abat31401':'EU','abat31402':'EU','abat31501':'EU','abat32301':'EU','abat32401':'EU','abat32501':'EU','abat32603':'EU','abat32701':'EU','abat4s125':'EU','abat4s156':'EU','abat4s404':'EU','abat4s409':'EU', \
                'abat4s416':'EU','abat4s417':'EU','abat51066':'EU','abat51200':'EU','abat54057':'EU','abat55018':'EU','abat55032':'EU','abat60018':'EU','abat60107':'EU','abat60118':'EU','abat60138':'EU','abat60143':'EU','abat60150':'EU','abat60170':'EU', \
                'abat60181':'EU','abat60182':'EU','abat60188':'EU','abat60189':'EU','abat60194':'EU','abat72106':'EU','abat72113':'EU','abat72123':'EU','abat72218':'EU','abat72547':'EU','abat72807':'EU','abat72908':'EU','abat80706':'EU','abat82708':'EU', \
                'abat900za':'EU','abat90laa':'EU','abat90lob':'EU','abat9jaeg':'EU','abat9stef':'EU','abbetb006':'EU','abbetb011':'EU','abbetm705':'EU','abbetn012':'EU','abbetn035':'EU','abbetn043':'EU','abbetn045':'EU','abbetn070':'EU','abbetr001':'EU', \
                'abbetr012':'EU','abbetr201':'EU','abbetr240':'EU','abbetr502':'EU','abbetr701':'EU','abbetr710':'EU','abbetr740':'EU','abbetr801':'EU','abbetr811':'EU','abbetr831':'EU','abbetwol1':'EU','abbg0012a':'EU','abbg0013a':'EU','abbg0041a':'EU', \
                'abbg0043a':'EU','abbg0045a':'EU','abbg0050a':'EU','abbg0051a':'EU','abbg0052a':'EU','abbg0053r':'EU','abch0001a':'EU','abch0005a':'EU','abch0008a':'EU','abch0010a':'EU','abch0011a':'EU','abch0013a':'EU','abch0014a':'EU','abch0017a':'EU', \
                'abch0022a':'EU','abch0042a':'EU','abch0046a':'EU','abch0047a':'EU','abch0048a':'EU','abch0049a':'EU','abch0050a':'EU','abcy0003a':'EU','abcz0akob':'EU','abcz0alib':'EU','abcz0asmi':'EU','abcz0asto':'EU','abcz0asuc':'EU','abcz0avel':'EU', \
                'abcz0bbnd':'EU','abcz0bbny':'EU','abcz0ccbd':'EU','abcz0ctab':'EU','abcz0epao':'EU','abcz0epau':'EU','abcz0hhkb':'EU','abcz0hhko':'EU','abcz0jjih':'EU','abcz0ksom':'EU','abcz0llim':'EU','abcz0mprr':'EU','abcz0mpst':'EU','abcz0pplb':'EU', \
                'abcz0ppll':'EU','abcz0pplv':'EU','abcz0sklm':'EU','abcz0smbo':'EU','abcz0tkar':'EU','abcz0toff':'EU','abcz0ttro':'EU','abcz0ultt':'EU','abcz0umom':'EU','abcz0utem':'EU','abcz0uulk':'EU','abcz0uulm':'EU','abcz0zzln':'EU','abdebb021':'EU', \
                'abdebb026':'EU','abdebb029':'EU','abdebb031':'EU','abdebb032':'EU','abdebb042':'EU','abdebb048':'EU','abdebb050':'EU','abdebb055':'EU','abdebb063':'EU','abdebb064':'EU','abdebb067':'EU','abdebe034':'EU','abdebe051':'EU','abdebw010':'EU', \
                'abdebw019':'EU','abdebw024':'EU','abdebw027':'EU','abdebw029':'EU','abdebw032':'EU','abdebw034':'EU','abdebw037':'EU','abdebw039':'EU','abdebw042':'EU','abdebw046':'EU','abdebw052':'EU','abdebw056':'EU','abdebw059':'EU','abdebw073':'EU', \
                'abdebw076':'EU','abdebw081':'EU','abdebw107':'EU','abdebw110':'EU','abdebw111':'EU','abdebw112':'EU','abdeby001':'EU','abdeby002':'EU','abdeby004':'EU','abdeby005':'EU','abdeby020':'EU','abdeby031':'EU','abdeby032':'EU','abdeby037':'EU', \
                'abdeby039':'EU','abdeby052':'EU','abdeby053':'EU','abdeby062':'EU','abdeby063':'EU','abdeby068':'EU','abdeby075':'EU','abdeby077':'EU','abdeby079':'EU','abdeby088':'EU','abdeby089':'EU','abdeby099':'EU','abdeby113':'EU','abdeby118':'EU', \
                'abdehb001':'EU','abdehb002':'EU','abdehb003':'EU','abdehb005':'EU','abdehe001':'EU','abdehe005':'EU','abdehe008':'EU','abdehe011':'EU','abdehe014':'EU','abdehe018':'EU','abdehe020':'EU','abdehe022':'EU','abdehe030':'EU','abdehe032':'EU', \
                'abdehe044':'EU','abdehe045':'EU','abdehe058':'EU','abdehh008':'EU','abdehh021':'EU','abdehh033':'EU','abdehh047':'EU','abdehh049':'EU','abdehh050':'EU','abdemv002':'EU','abdemv003':'EU','abdemv005':'EU','abdemv006':'EU','abdemv007':'EU', \
                'abdemv018':'EU','abdemv019':'EU','abdeni011':'EU','abdeni016':'EU','abdeni020':'EU','abdeni028':'EU','abdeni029':'EU','abdeni038':'EU','abdeni041':'EU','abdeni042':'EU','abdeni043':'EU','abdeni052':'EU','abdeni053':'EU','abdeni054':'EU', \
                'abdeni062':'EU','abdenw006':'EU','abdenw008':'EU','abdenw015':'EU','abdenw021':'EU','abdenw029':'EU','abdenw030':'EU','abdenw038':'EU','abdenw042':'EU','abdenw050':'EU','abdenw058':'EU','abdenw067':'EU','abdenw078':'EU','abdenw079':'EU', \
                'abdenw080':'EU','abdenw094':'EU','abdenw095':'EU','abdenw096':'EU','abdenw114':'EU','abdenw179':'EU','abderp001':'EU','abderp007':'EU','abderp018':'EU','abderp021':'EU','abderp022':'EU','abderp023':'EU','abderp024':'EU','abderp025':'EU', \
                'abderp027':'EU','abderp034':'EU','abderp040':'EU','abdesh016':'EU','abdesh023':'EU','abdesl002':'EU','abdesl011':'EU','abdesl017':'EU','abdesl018':'EU','abdesn001':'EU','abdesn004':'EU','abdesn011':'EU','abdesn012':'EU','abdesn017':'EU', \
                'abdesn019':'EU','abdesn024':'EU','abdesn045':'EU','abdesn050':'EU','abdesn059':'EU','abdesn061':'EU','abdesn081':'EU','abdesn082':'EU','abdest002':'EU','abdest011':'EU','abdest015':'EU','abdest028':'EU','abdest030':'EU','abdest044':'EU', \
                'abdest050':'EU','abdest057':'EU','abdest066':'EU','abdest072':'EU','abdest077':'EU','abdest078':'EU','abdest090':'EU','abdest097':'EU','abdeth005':'EU','abdeth009':'EU','abdeth013':'EU','abdeth018':'EU','abdeth020':'EU','abdeth025':'EU', \
                'abdeth036':'EU','abdeth041':'EU','abdeth060':'EU','abdk0030a':'EU','abdk0034a':'EU','abdk0046a':'EU','abdk0052a':'EU','abdk0053a':'EU','abee0015a':'EU','abee0018a':'EU','abee0019a':'EU','abes0110a':'EU','abes0113a':'EU','abes0115a':'EU', \
                'abes0116a':'EU','abes0117a':'EU','abes0119a':'EU','abes0120a':'EU','abes0121a':'EU','abes0122a':'EU','abes0123a':'EU','abes0124a':'EU','abes0126a':'EU','abes0201a':'EU','abes0373a':'EU','abes0377a':'EU','abes0584a':'EU','abes0587a':'EU', \
                'abes0588a':'EU','abes0625a':'EU','abes0651a':'EU','abes0694a':'EU','abes0752a':'EU','abes0805a':'EU','abes0822a':'EU','abes0824a':'EU','abes0825a':'EU','abes0879a':'EU','abes0880a':'EU','abes0890a':'EU','abes0892a':'EU','abes0905a':'EU', \
                'abes0971a':'EU','abes1038a':'EU','abes1044a':'EU','abes1048a':'EU','abes1050a':'EU','abes1072a':'EU','abes1076a':'EU','abes1090a':'EU','abes1094a':'EU','abes1096a':'EU','abes1117a':'EU','abes1123a':'EU','abes1125a':'EU','abes1126a':'EU', \
                'abes1131a':'EU','abes1132a':'EU','abes1133a':'EU','abes1135a':'EU','abes1137a':'EU','abes1138a':'EU','abes1148a':'EU','abes1161a':'EU','abes1162a':'EU','abes1163a':'EU','abes1164a':'EU','abes1169a':'EU','abes1172a':'EU','abes1181a':'EU', \
                'abes1182a':'EU','abes1185a':'EU','abes1192a':'EU','abes1193a':'EU','abes1208a':'EU','abes1215a':'EU','abes1224a':'EU','abes1225a':'EU','abes1231a':'EU','abes1239a':'EU','abes1240a':'EU','abes1244a':'EU','abes1246a':'EU','abes1247a':'EU', \
                'abes1262a':'EU','abes1268a':'EU','abes1269a':'EU','abes1270a':'EU','abes1271a':'EU','abes1272a':'EU','abes1275a':'EU','abes1277a':'EU','abes1278a':'EU','abes1279a':'EU','abes1284a':'EU','abes1286a':'EU','abes1287a':'EU','abes1291a':'EU', \
                'abes1292a':'EU','abes1295a':'EU','abes1296a':'EU','abes1297a':'EU','abes1320a':'EU','abes1333a':'EU','abes1339a':'EU','abes1346a':'EU','abes1349a':'EU','abes1350a':'EU','abes1351a':'EU','abes1353a':'EU','abes1358a':'EU','abes1363a':'EU', \
                'abes1365a':'EU','abes1370a':'EU','abes1372a':'EU','abes1380a':'EU','abes1386a':'EU','abes1387a':'EU','abes1393a':'EU','abes1405a':'EU','abes1421a':'EU','abes1422a':'EU','abes1423a':'EU','abes1424a':'EU','abes1426a':'EU','abes1427a':'EU', \
                'abes1428a':'EU','abes1432a':'EU','abes1433a':'EU','abes1435a':'EU','abes1437a':'EU','abes1438a':'EU','abes1442a':'EU','abes1443a':'EU','abes1445a':'EU','abes1449a':'EU','abes1450a':'EU','abes1453a':'EU','abes1464a':'EU','abes1472a':'EU', \
                'abes1479a':'EU','abes1480a':'EU','abes1490a':'EU','abes1492a':'EU','abes1494a':'EU','abes1496a':'EU','abes1497a':'EU','abes1498a':'EU','abes1499a':'EU','abes1501a':'EU','abes1502a':'EU','abes1516a':'EU','abes1519a':'EU','abes1525a':'EU', \
                'abes1529a':'EU','abes1530a':'EU','abes1532a':'EU','abes1535a':'EU','abes1536a':'EU','abes1537a':'EU','abes1544a':'EU','abes1548a':'EU','abes1549a':'EU','abes1560a':'EU','abes1564a':'EU','abes1567a':'EU','abes1568a':'EU','abes1569a':'EU', \
                'abes1572a':'EU','abes1573a':'EU','abes1574a':'EU','abes1576a':'EU','abes1577a':'EU','abes1578a':'EU','abes1587a':'EU','abes1589a':'EU','abes1593a':'EU','abes1595a':'EU','abes1596a':'EU','abes1597a':'EU','abes1598a':'EU','abes1602a':'EU', \
                'abes1604a':'EU','abes1610a':'EU','abes1611a':'EU','abes1612a':'EU','abes1613a':'EU','abes1615a':'EU','abes1617a':'EU','abes1618a':'EU','abes1619a':'EU','abes1620a':'EU','abes1623a':'EU','abes1624a':'EU','abes1630a':'EU','abes1632a':'EU', \
                'abes1633a':'EU','abes1635a':'EU','abes1636a':'EU','abes1638a':'EU','abes1640a':'EU','abes1641a':'EU','abes1642a':'EU','abes1643a':'EU','abes1646a':'EU','abes1650a':'EU','abes1651a':'EU','abes1653a':'EU','abes1656a':'EU','abes1657a':'EU', \
                'abes1658a':'EU','abes1666a':'EU','abes1672a':'EU','abes1675a':'EU','abes1676a':'EU','abes1677a':'EU','abes1679a':'EU','abes1684a':'EU','abes1685a':'EU','abes1696a':'EU','abes1697a':'EU','abes1709a':'EU','abes1710a':'EU','abes1711a':'EU', \
                'abes1713a':'EU','abes1740a':'EU','abes1747a':'EU','abes1749a':'EU','abes1750a':'EU','abes1751a':'EU','abes1752a':'EU','abes1765a':'EU','abes1783a':'EU','abes1786a':'EU','abes1787a':'EU','abes1788a':'EU','abes1790a':'EU','abes1791a':'EU', \
                'abes1799a':'EU','abes1800a':'EU','abfi00226':'EU','abfi00363':'EU','abfi00425':'EU','abfr01001':'EU','abfr01006':'EU','abfr01008':'EU','abfr01011':'EU','abfr01012':'EU','abfr01015':'EU','abfr01017':'EU','abfr01018':'EU','abfr01019':'EU', \
                'abfr01020':'EU','abfr01021':'EU','abfr02001':'EU','abfr02004':'EU','abfr02012':'EU','abfr02013':'EU','abfr02016':'EU','abfr02017':'EU','abfr02019':'EU','abfr02020':'EU','abfr02021':'EU','abfr03019':'EU','abfr03029':'EU','abfr03036':'EU', \
                'abfr03037':'EU','abfr03043':'EU','abfr03048':'EU','abfr03063':'EU','abfr03064':'EU','abfr03067':'EU','abfr03069':'EU','abfr03080':'EU','abfr03083':'EU','abfr04002':'EU','abfr04004':'EU','abfr04008':'EU','abfr04017':'EU','abfr04018':'EU', \
                'abfr04023':'EU','abfr04034':'EU','abfr04037':'EU','abfr04048':'EU','abfr04049':'EU','abfr04055':'EU','abfr04063':'EU','abfr04069':'EU','abfr04098':'EU','abfr04100':'EU','abfr04101':'EU','abfr04145':'EU','abfr04149':'EU','abfr04160':'EU', \
                'abfr04299':'EU','abfr04319':'EU','abfr05010':'EU','abfr05040':'EU','abfr05064':'EU','abfr05074':'EU','abfr05079':'EU','abfr05082':'EU','abfr05084':'EU','abfr06001':'EU','abfr06003':'EU','abfr06004':'EU','abfr06007':'EU','abfr06008':'EU', \
                'abfr06009':'EU','abfr06011':'EU','abfr06012':'EU','abfr07001':'EU','abfr07004':'EU','abfr07006':'EU','abfr07008':'EU','abfr07009':'EU','abfr07012':'EU','abfr07013':'EU','abfr07014':'EU','abfr07016':'EU','abfr07017':'EU','abfr07018':'EU', \
                'abfr07030':'EU','abfr07032':'EU','abfr08005':'EU','abfr08016':'EU','abfr08017':'EU','abfr08018':'EU','abfr08022':'EU','abfr08023':'EU','abfr08614':'EU','abfr08617':'EU','abfr08712':'EU','abfr09002':'EU','abfr09003':'EU','abfr09008':'EU', \
                'abfr09009':'EU','abfr09010':'EU','abfr09013':'EU','abfr09014':'EU','abfr09015':'EU','abfr09016':'EU','abfr09017':'EU','abfr09019':'EU','abfr09020':'EU','abfr09103':'EU','abfr10009':'EU','abfr10015':'EU','abfr10025':'EU','abfr10032':'EU', \
                'abfr11002':'EU','abfr11016':'EU','abfr11022':'EU','abfr11026':'EU','abfr11027':'EU','abfr11029':'EU','abfr11030':'EU','abfr11033':'EU','abfr12001':'EU','abfr12004':'EU','abfr12017':'EU','abfr12021':'EU','abfr12024':'EU','abfr12025':'EU', \
                'abfr12026':'EU','abfr12027':'EU','abfr12030':'EU','abfr12034':'EU','abfr12041':'EU','abfr12042':'EU','abfr13002':'EU','abfr13008':'EU','abfr13009':'EU','abfr13012':'EU','abfr13014':'EU','abfr14002':'EU','abfr14004':'EU','abfr14009':'EU', \
                'abfr14010':'EU','abfr14021':'EU','abfr14031':'EU','abfr14032':'EU','abfr14033':'EU','abfr14049':'EU','abfr15013':'EU','abfr15017':'EU','abfr15031':'EU','abfr15038':'EU','abfr15043':'EU','abfr15044':'EU','abfr16001':'EU','abfr16029':'EU', \
                'abfr16036':'EU','abfr16038':'EU','abfr16053':'EU','abfr16060':'EU','abfr16064':'EU','abfr16065':'EU','abfr16066':'EU','abfr17004':'EU','abfr17009':'EU','abfr17013':'EU','abfr17014':'EU','abfr17018':'EU','abfr17023':'EU','abfr18008':'EU', \
                'abfr18019':'EU','abfr18025':'EU','abfr18032':'EU','abfr18035':'EU','abfr18036':'EU','abfr18038':'EU','abfr18040':'EU','abfr18042':'EU','abfr18043':'EU','abfr18044':'EU','abfr19003':'EU','abfr19004':'EU','abfr19006':'EU','abfr19011':'EU', \
                'abfr19012':'EU','abfr19021':'EU','abfr19031':'EU','abfr19032':'EU','abfr19051':'EU','abfr19061':'EU','abfr19071':'EU','abfr19081':'EU','abfr19091':'EU','abfr20004':'EU','abfr20017':'EU','abfr20036':'EU','abfr20037':'EU','abfr20045':'EU', \
                'abfr20046':'EU','abfr20047':'EU','abfr20048':'EU','abfr21001':'EU','abfr21017':'EU','abfr21019':'EU','abfr21021':'EU','abfr21030':'EU','abfr21035':'EU','abfr21040':'EU','abfr22012':'EU','abfr22014':'EU','abfr22015':'EU','abfr22018':'EU', \
                'abfr22019':'EU','abfr23076':'EU','abfr23078':'EU','abfr23110':'EU','abfr23123':'EU','abfr23128':'EU','abfr23152':'EU','abfr23157':'EU','abfr23163':'EU','abfr23172':'EU','abfr24007':'EU','abfr24009':'EU','abfr24011':'EU','abfr24013':'EU', \
                'abfr24015':'EU','abfr24017':'EU','abfr24018':'EU','abfr24019':'EU','abfr24020':'EU','abfr24024':'EU','abfr25032':'EU','abfr25036':'EU','abfr25038':'EU','abfr25039':'EU','abfr25040':'EU','abfr25043':'EU','abfr25046':'EU','abfr26001':'EU', \
                'abfr26002':'EU','abfr26005':'EU','abfr26007':'EU','abfr26010':'EU','abfr26011':'EU','abfr26015':'EU','abfr26016':'EU','abfr26017':'EU','abfr27002':'EU','abfr27003':'EU','abfr27004':'EU','abfr27005':'EU','abfr27007':'EU','abfr28002':'EU', \
                'abfr28006':'EU','abfr28010':'EU','abfr28014':'EU','abfr28018':'EU','abfr28019':'EU','abfr28020':'EU','abfr28028':'EU','abfr28030':'EU','abfr29421':'EU','abfr29423':'EU','abfr29424':'EU','abfr29425':'EU','abfr29426':'EU','abfr30016':'EU', \
                'abfr30017':'EU','abfr30018':'EU','abfr30019':'EU','abfr30020':'EU','abfr30021':'EU','abfr30022':'EU','abfr30023':'EU','abfr30024':'EU','abfr30026':'EU','abfr30034':'EU','abfr31001':'EU','abfr31002':'EU','abfr31004':'EU','abfr31007':'EU', \
                'abfr31013':'EU','abfr31014':'EU','abfr31016':'EU','abfr31018':'EU','abfr31021':'EU','abfr31030':'EU','abfr31031':'EU','abfr31032':'EU','abfr31033':'EU','abfr31034':'EU','abfr31035':'EU','abfr31036':'EU','abfr32001':'EU','abfr32002':'EU', \
                'abfr32004':'EU','abfr32005':'EU','abfr32006':'EU','abfr32007':'EU','abfr33101':'EU','abfr33102':'EU','abfr33103':'EU','abfr33111':'EU','abfr33120':'EU','abfr33121':'EU','abfr33201':'EU','abfr33202':'EU','abfr33211':'EU','abfr33212':'EU', \
                'abfr33260':'EU','abfr33301':'EU','abfr33303':'EU','abfr34011':'EU','abfr34012':'EU','abfr34014':'EU','abfr34017':'EU','abfr34018':'EU','abfr34021':'EU','abfr34024':'EU','abfr34025':'EU','abfr34031':'EU','abfr34032':'EU','abfr34033':'EU', \
                'abfr34041':'EU','abfr34042':'EU','abfr34044':'EU','abfr34051':'EU','abfr34052':'EU','abfr34061':'EU','abfr34062':'EU','abfr35002':'EU','abfr35003':'EU','abfr35004':'EU','abfr35005':'EU','abfr35006':'EU','abfr35007':'EU','abfr36001':'EU', \
                'abfr36002':'EU','abfr36004':'EU','abfr37001':'EU','abfr37002':'EU','abfr38002':'EU','abfr38008':'EU','abfr39001':'EU','abfr39002':'EU','abfr40001':'EU','abgb0032r':'EU','abgb0566a':'EU','abgb0567a':'EU','abgb0568a':'EU','abgb0569a':'EU', \
                'abgb0580a':'EU','abgb0583a':'EU','abgb0584a':'EU','abgb0586a':'EU','abgb0597a':'EU','abgb0598a':'EU','abgb0613a':'EU','abgb0615a':'EU','abgb0640a':'EU','abgb0641a':'EU','abgb0642a':'EU','abgb0643a':'EU','abgb0673a':'EU','abgb0682a':'EU', \
                'abgr0001a':'EU','abgr0002a':'EU','abgr0022a':'EU','abgr0027a':'EU','abgr0031a':'EU','abgr0035a':'EU','abgr0037a':'EU','abgr0038a':'EU','abgr0039a':'EU','abgr0043a':'EU','abhu0022a':'EU','abhu0023a':'EU','abhu0025a':'EU','abhu0026a':'EU', \
                'abhu0029a':'EU','abhu0032a':'EU','abhu0033a':'EU','abhu0034a':'EU','abhu0035a':'EU','abhu0036a':'EU','abhu0038a':'EU','abhu0039a':'EU','abhu0042a':'EU','abie0028a':'EU','abit0451a':'EU','abit0459a':'EU','abit0461a':'EU','abit0462a':'EU', \
                'abit0463a':'EU','abit0499a':'EU','abit0501a':'EU','abit0503a':'EU','abit0505a':'EU','abit0506a':'EU','abit0507a':'EU','abit0508a':'EU','abit0545a':'EU','abit0554a':'EU','abit0589a':'EU','abit0591a':'EU','abit0611a':'EU','abit0612a':'EU', \
                'abit0614a':'EU','abit0620a':'EU','abit0703a':'EU','abit0710a':'EU','abit0727a':'EU','abit0732a':'EU','abit0740a':'EU','abit0743a':'EU','abit0753a':'EU','abit0754a':'EU','abit0775a':'EU','abit0822a':'EU','abit0839a':'EU','abit0854a':'EU', \
                'abit0856a':'EU','abit0858a':'EU','abit0869a':'EU','abit0883a':'EU','abit0902a':'EU','abit0908a':'EU','abit0912a':'EU','abit0929a':'EU','abit0940a':'EU','abit0948a':'EU','abit0977a':'EU','abit0980a':'EU','abit0983a':'EU','abit0990a':'EU', \
                'abit0991a':'EU','abit0995a':'EU','abit1010a':'EU','abit1017a':'EU','abit1037a':'EU','abit1048a':'EU','abit1076a':'EU','abit1079a':'EU','abit1087a':'EU','abit1088a':'EU','abit1110a':'EU','abit1112a':'EU','abit1120a':'EU','abit1125a':'EU', \
                'abit1137a':'EU','abit1166a':'EU','abit1167a':'EU','abit1168a':'EU','abit1178a':'EU','abit1180a':'EU','abit1182a':'EU','abit1191a':'EU','abit1203a':'EU','abit1216a':'EU','abit1246a':'EU','abit1247a':'EU','abit1259a':'EU','abit1269a':'EU', \
                'abit1279a':'EU','abit1335a':'EU','abit1374a':'EU','abit1377a':'EU','abit1395a':'EU','abit1397a':'EU','abit1407a':'EU','abit1421a':'EU','abit1423a':'EU','abit1427a':'EU','abit1465a':'EU','abit1466a':'EU','abit1478a':'EU','abit1515a':'EU', \
                'abit1518a':'EU','abit1523a':'EU','abit1524a':'EU','abit1529a':'EU','abit1532a':'EU','abit1538a':'EU','abit1551a':'EU','abit1576a':'EU','abit1586a':'EU','abit1587a':'EU','abit1598a':'EU','abit1606a':'EU','abit1650a':'EU','abit1654a':'EU', \
                'abit1662a':'EU','abit1696a':'EU','abit1697a':'EU','ablt00002':'EU','ablt00003':'EU','ablt00041':'EU','ablt00042':'EU','ablt00043':'EU','ablu0100a':'EU','ablu0101a':'EU','ablv00rm4':'EU','abmk0028a':'EU','abmk0031a':'EU','abmk0034a':'EU', \
                'abmk0035a':'EU','abmk0036a':'EU','abmk0037a':'EU','abmk0038a':'EU','abmk0040a':'EU','abmk0041a':'EU','abmk0043a':'EU','abnl00137':'EU','abnl00241':'EU','abnl00411':'EU','abnl00441':'EU','abnl00447':'EU','abnl00537':'EU','abnl00544':'EU', \
                'abnl00742':'EU','abnl00938':'EU','abpl0003r':'EU','abpl0010a':'EU','abpl0038a':'EU','abpl0044a':'EU','abpl0048a':'EU','abpl0049a':'EU','abpl0052a':'EU','abpl0096a':'EU','abpl0106a':'EU','abpl0129a':'EU','abpl0138a':'EU','abpl0139a':'EU', \
                'abpl0149a':'EU','abpl0171a':'EU','abpl0175a':'EU','abpl0176a':'EU','abpl0189a':'EU','abpl0192a':'EU','abpl0193a':'EU','abpl0194a':'EU','abpl0209a':'EU','abpl0213a':'EU','abpl0218a':'EU','abpl0234a':'EU','abpl0237a':'EU','abpl0241a':'EU', \
                'abpl0245a':'EU','abpl0248a':'EU','abpt01022':'EU','abpt01023':'EU','abpt01024':'EU','abpt01025':'EU','abpt01028':'EU','abpt01031':'EU','abpt01042':'EU','abpt01045':'EU','abpt02004':'EU','abpt02016':'EU','abpt02018':'EU','abpt03063':'EU', \
                'abpt03070':'EU','abpt03071':'EU','abpt03082':'EU','abpt03083':'EU','abpt03084':'EU','abpt03085':'EU','abpt03087':'EU','abpt03089':'EU','abpt03093':'EU','abpt03095':'EU','abpt03097':'EU','abpt04001':'EU','abpt05006':'EU','abpt05007':'EU', \
                'abpt05008':'EU','abpt06005':'EU','abpt06006':'EU','abro0011a':'EU','abro0056a':'EU','abro0063a':'EU','abro0064a':'EU','abro0065a':'EU','abro0066a':'EU','abro0067a':'EU','abro0068a':'EU','abro0069a':'EU','abro0071a':'EU','abse0001a':'EU', \
                'abse0004a':'EU','abse0052a':'EU','absi0001a':'EU','absi0002a':'EU','absi0003a':'EU','absi0034a':'EU','absi0035a':'EU','absi0036a':'EU','absi0037a':'EU','absi0038a':'EU','absk0001a':'EU','absk0002r':'EU','absk0020a':'EU','absk0025a':'EU', \
                'absk0037a':'EU','absk0048a':'EU','abnl00404':'EU','abnl00418':'EU','abnl00433':'EU','abnl00520':'EU','abfr0250a':'EU','abfr0526a':'EU','abes0596a':'EU','abes0691a':'EU','abgr0032a':'EU','abnl00636':'EU','abnl00638':'EU','abnl00639':'EU', \
                'abnl00640':'EU','abdebe016':'EU','abdebe033':'EU','abdehh006':'EU','abdeni036':'EU','abderp020':'EU','abdesh006':'EU','abdeub009':'EU','abgb0034r':'EU','abgr0029a':'EU','abpt00314':'EU','abpt04005':'EU',
                'abdebw012':'EU','abdehe010':'EU','abnl00236':'EU','abbetn050':'EU','abdebe001':'EU','abdebe014':'EU','abdebw008':'EU','abdebw016':'EU','abdebw020':'EU','abdebw021':'EU','abdebw028':'EU','abdebw033':'EU','abdebw035':'EU','abdebw036':'EU',
                'abdebw040':'EU','abdebw041':'EU','abdebw044':'EU','abdebw045':'EU','abdebw047':'EU','abdebw048':'EU','abdebw049':'EU','abdebw050':'EU','abdebw053':'EU','abdebw054':'EU','abdebw057':'EU','abdebw060':'EU','abdebw065':'EU','abdeby008':'EU',
                'abdeby012':'EU','abdeby034':'EU','abdeby061':'EU','abdeby082':'EU','abdehe003':'EU','abdehe007':'EU','abdehe012':'EU','abdehe025':'EU','abdehe029':'EU','abdehe031':'EU','abdehe033':'EU','abdehh005':'EU','abdenw004':'EU','abdenw010':'EU',
                'abdenw013':'EU','abdenw018':'EU','abdenw036':'EU','abdenw047':'EU','abdenw062':'EU','abderp003':'EU','abdesl009':'EU','abdeub016':'EU','abdeub020':'EU','abfi00204':'EU','abgb0536a':'EU','abgr0030a':'EU','abnl00133':'EU','abnl00238':'EU','abnl00729':'EU',
                'abat51002':'EU','abdebb003':'EU','abdebe044':'EU','abdehe019':'EU','abderp026':'EU','abch0002a':'EU','abch0003a':'EU','abch0004a':'EU','abch0018a':'EU','abch0021a':'EU','abch0023a':'EU','abch0030a':'EU','abdebw072':'EU','abdeby073':'EU',
                'abdehe009':'EU','abdehe016':'EU','abdehe017':'EU','abdehe035':'EU','abdesn014':'EU','abdesn025':'EU','abdesn036':'EU','abdest005':'EU','abdeth031':'EU','abgb0581a':'EU','abgr0028a':'EU','abat2he19':'EU','abat2kl17':'EU','abat2vk16':'EU',
                'abat2vl12':'EU','abat4s410':'EU','abcz0avyc':'EU','abdebb001':'EU','abdebb006':'EU','abdebb024':'EU','abdeby017':'EU','abdehe004':'EU','abdehe015':'EU','abdehh022':'EU','abdest022':'EU','abdest025':'EU','abdest031':'EU','abdest032':'EU',
                'abdest033':'EU','abdest034':'EU','abgb0585a':'EU','abgb0594a':'EU','abgb0595a':'EU','abat55016':'EU','abcz0bbnf':'EU','abdebb020':'EU','abdebb028':'EU','abdebb030':'EU','abdebb038':'EU','abdesh015':'EU','abdesn028':'EU','abdeth024':'EU',
                'abgb0596a':'EU','abgb0608a':'EU','abgb0609a':'EU','abch0027a':'EU','abch0038a':'EU','abch0039a':'EU','abcz0hhks':'EU','abdeby064':'EU','abderp008':'EU','abdesn002':'EU','abdesn005':'EU','abdesn006':'EU','abdesn020':'EU','abdesn034':'EU',
                'abdesn047':'EU','abdesn048':'EU','abdesn060':'EU','abgb0614a':'EU','abat2fe13':'EU','abcz0arep':'EU','abdebb012':'EU','abdebb018':'EU','abdebb034':'EU','abdebb036':'EU','abdeni039':'EU','abdeni040':'EU','abdest036':'EU','abdest042':'EU',
                'abdest052':'EU','abdest061':'EU','abdest063':'EU','abdest071':'EU','abdest073':'EU','abdest076':'EU','abdeth016':'EU','abdeth021':'EU','abgb0616a':'EU','abgb0620a':'EU','abgb0621a':'EU','abgb0622a':'EU','abgb0638a':'EU','abgb0644a':'EU',
                'abgb0645a':'EU','abgb0646a':'EU','abgb0649a':'EU','abdenw051':'EU','abes0817a':'EU','abes0975a':'EU','abes1016a':'EU','abes1019a':'EU','abes1046a':'EU','abes1047a':'EU','abes1119a':'EU','abes1171a':'EU','abes1177a':'EU','abes1178a':'EU',
                'abes1179a':'EU','abes1183a':'EU','abes1184a':'EU','abes1195a':'EU','abes1238a':'EU','abes1250a':'EU','abes1251a':'EU','abes1254a':'EU','abes1280a':'EU','abes1281a':'EU','abes1329a':'EU','abes1357a':'EU','abes1364a':'EU','abes1371a':'EU',
                'abes1394a':'EU','abes1404a':'EU','abes1417a':'EU','abes1418a':'EU','abes1419a':'EU','abgb0650a':'EU','abgb0651a':'EU','abgb0652a':'EU','abgb0656a':'EU','abgb0658a':'EU','abgb0660a':'EU','abgb0672a':'EU','abgb0677a':'EU','abgb0679a':'EU',
                'abgb0683a':'EU','abgb0684a':'EU','abgb0687a':'EU','abgb0689a':'EU','ablt0013a':'EU','abpl0008a':'EU','abpl0013a':'EU','abpl0017a':'EU','abpl0022a':'EU','abpl0023a':'EU','abpt02002':'EU','abpt03072':'EU','abpt04004':'EU',
                'abdebb009':'EU','abdesn075':'EU','abes0125a':'EU','abes1051a':'EU','abes1260a':'EU','abgb0698a':'EU','abpl0011a':'EU','abpl0031a':'EU','abpt03027':'EU','abbg0044a':'EU','abdeth019':'EU','abdeth032':'EU','abee0013a':'EU','abes1018a':'EU',
                'abes1027a':'EU','abes1084a':'EU','abes1122a':'EU','abes1399a':'EU','abes1482a':'EU','abfr01003':'EU','abfr02008':'EU','abfr02022':'EU','abfr02025':'EU','abfr02031':'EU','abfr03022':'EU','abfr03024':'EU','abfr03028':'EU','abfr03035':'EU',
                'abfr03045':'EU','abfr03046':'EU','abfr03061':'EU','abfr03065':'EU','abfr03081':'EU','abfr03084':'EU','abfr04060':'EU','abfr05005':'EU','abfr06002':'EU','abfr06005':'EU','abfr07003':'EU','abfr07019':'EU','abfr08613':'EU','abfr08615':'EU',
                'abfr08713':'EU','abfr08714':'EU','abfr09011':'EU','abfr10002':'EU','abfr10023':'EU','abfr14001':'EU','abfr16027':'EU','abfr17002':'EU','abfr17006':'EU','abfr18006':'EU','abfr18034':'EU','abfr20003':'EU','abfr20020':'EU','abfr21007':'EU',
                'abfr21009':'EU','abfr23059':'EU','abfr23065':'EU','abfr23092':'EU','abfr23093':'EU','abfr23099':'EU','abfr23107':'EU','abfr23112':'EU','abfr23113':'EU','abfr23118':'EU','abfr23119':'EU','abfr23120':'EU','abfr24004':'EU','abfr24005':'EU',
                'abfr24010':'EU','abfr25030':'EU','abfr28001':'EU','abfr28011':'EU','abfr30031':'EU','abfr32003':'EU','abhu0020a':'EU','abie0093a':'EU','abit0466a':'EU','abit0476a':'EU','abit0594a':'EU','abit0705a':'EU','abit0771a':'EU','abit0777a':'EU',
                'abit0892a':'EU','abit0905a':'EU','abit1034a':'EU','abit1101a':'EU','abit1160a':'EU','abit1161a':'EU','abit1208a':'EU','abit1488a':'EU','ablv00vn1':'EU','ablv0rke2':'EU','abpl0019a':'EU','abpt01021':'EU','abpt02005':'EU','abpt03055':'EU',
                'abse0022a':'EU','abat4s425':'EU','abdebb043':'EU','abdebb052':'EU','abdebw094':'EU','abdebw102':'EU','abdest014':'EU','abfr23103':'EU','abfr29427':'EU','abgb0727a':'EU','abgb0728a':'EU','abgb0730a':'EU','abgb0731a':'EU','abit1033a':'EU',
                'abit1360a':'EU','abpl0047a':'EU','abpt01034':'EU','absk0005a':'EU','absk0012a':'EU','absk0016a':'EU','absk0018a':'EU','absk0022a':'EU','absk0038a':'EU','abbetb004':'EU','abdk0047a':'EU','abes1237a':'EU','abes1420a':'EU','abes1451a':'EU',
                'abes1481a':'EU','abes1566a':'EU','abes1600a':'EU','abes1609a':'EU','abfr06015':'EU','abfr16037':'EU','abfr23125':'EU','abfr23126':'EU','abfr25023':'EU','abfr27102':'EU','abfr34023':'EU','abgb0739a':'EU','abgb0743a':'EU','abgr0020a':'EU',
                'abgr0041a':'EU','abgr0042a':'EU','abgr0044a':'EU','abgr0045a':'EU','abgr0046a':'EU','abgr0047a':'EU','abis0004a':'EU','abit0448a':'EU','abit0963a':'EU','abit1192a':'EU','abit1343a':'EU','abdest029':'EU','abdest080':'EU','abes0114a':'EU',
                'abes1045a':'EU','abes1288a':'EU','abes1289a':'EU','abes1298a':'EU','abes1354a':'EU','abes1485a':'EU','abes1521a':'EU','abes1580a':'EU','abes1601a':'EU','abes1627a':'EU','abfr06017':'EU','abfr10005':'EU','abfr11032':'EU','abfr23131':'EU',
                'abfr35008':'EU','abgb0776a':'EU','abpl0029a':'EU','abpt01027':'EU','abpt01032':'EU','abbetr841':'EU','abcz0topr':'EU','abes0556a':'EU','abes0789a':'EU','abes1049a':'EU','abes1644a':'EU','abfr23137':'EU','abfr26018':'EU','abfr38001':'EU',
                'abgb0729a':'EU','abgb0733a':'EU','abgb0738a':'EU','abgb0741a':'EU','abgb0744a':'EU','abgb0777a':'EU','abgb0839a':'EU','abgb0840a':'EU','abie0107a':'EU','abit0804a':'EU','abit1215a':'EU','abit1218a':'EU','abit1232a':'EU','abit1270a':'EU',
                'abit1393a':'EU','abit1453a':'EU','abit1469a':'EU','ablt00022':'EU','ablt00031':'EU','ablv00rp5':'EU','abpt01026':'EU','abpt01033':'EU','abrs0007a':'EU','absk0046a':'EU','abba0031a':'EU','abes1674a':'EU','abfr23150':'EU','abgb0837a':'EU',
                'abgb0851a':'EU','abgb0860a':'EU','abgb0863a':'EU','abgb0864a':'EU','abit1199a':'EU','abmk0039a':'EU','abpl0141a':'EU','abro0059a':'EU','abro0061a':'EU','abro0070a':'EU','abat30104':'EU','abat60186':'EU','abat60195':'EU','abba0029a':'EU',
                'abdemv021':'EU','abes0365a':'EU','abes1095a':'EU','abes1397a':'EU','abes1579a':'EU','abes1649a':'EU','abes1771a':'EU','abes1801a':'EU','abes1803a':'EU','abes1804a':'EU','abes1807a':'EU','abes1809a':'EU','abes1812a':'EU','abes1814a':'EU',
                'abes1815a':'EU','abes1816a':'EU','abes1818a':'EU','abes1819a':'EU','abes1820a':'EU','abes1821a':'EU','abes1822a':'EU','abes1823a':'EU','abes1824a':'EU','abes1825a':'EU','abes1826a':'EU','abes1832a':'EU','abes1833a':'EU','abes1834a':'EU',
                'abes1836a':'EU','abes1837a':'EU','abfr06010':'EU','abhu0031a':'EU','abit0524a':'EU','abit0813a':'EU','abit0880a':'EU','abit1144a':'EU','abit1152a':'EU','abit1193a':'EU','abit1282a':'EU','abit1480a':'EU','abit1485a':'EU','abit1550a':'EU',
                'abit1573a':'EU','abit1578a':'EU','abit1590a':'EU','abit1594a':'EU','abit1602a':'EU','abit1618a':'EU','abit1644a':'EU','abit1674a':'EU','abit1734a':'EU','abit1735a':'EU','abit1743a':'EU','abit1744a':'EU','abit1745a':'EU','abit1771a':'EU',
                'abit1801a':'EU','ablt00033':'EU','ablv0rvl7':'EU','abmt00004':'EU','abmt00005':'EU','abno0081a':'EU','abpl0184a':'EU','abpl0252a':'EU','abpl0253a':'EU','abpt01040':'EU','abro0081a':'EU','abro0085a':'EU','abat31901':'EU','abat60198':'EU',
                'abch0043a':'EU','abdebb007':'EU','abdebb075':'EU','abdehe013':'EU','abdesh033':'EU','abes0622a':'EU','abes0624a':'EU','abes1033a':'EU','abes1034a':'EU','abes1130a':'EU','abes1625a':'EU','abes1639a':'EU','abes1691a':'EU','abes1781a':'EU',
                'abes1817a':'EU','abes1828a':'EU','abes1835a':'EU','abes1838a':'EU','abes1848a':'EU','abes1849a':'EU','abes1857a':'EU','abes1858a':'EU','abes1860a':'EU','abes1867a':'EU','abes1868a':'EU','abes1869a':'EU','abes1889a':'EU','abfr03086':'EU',
                'abfr20062':'EU','abfr41003':'EU','abgr0008a':'EU','abhu0037a':'EU','abit0944a':'EU','abit0950a':'EU','abit1027a':'EU','abit1030a':'EU','abit1067a':'EU','abit1089a':'EU','abit1239a':'EU','abit1342a':'EU','abit1357a':'EU','abit1424a':'EU',
                'abit1535a':'EU','abit1555a':'EU','abit1593a':'EU','abit1656a':'EU','abit1795a':'EU','abit1796a':'EU','abit1835a':'EU','abit1836a':'EU','abit1840a':'EU','abit1857a':'EU','abit1866a':'EU','abit1876a':'EU','abit1880a':'EU','abit1883a':'EU',
                'ablv000l1':'EU','abmk0044a':'EU','abno0015a':'EU','abpl0239a':'EU','abpl0240a':'EU','abpl0295a':'EU','abpt03098':'EU','abro0080a':'EU','abat4s173':'EU','abbg0040a':'EU','abch0044a':'EU','abdebb082':'EU','abdemv022':'EU','abdenw247':'EU',
                'abdesn092':'EU','abes0893a':'EU','abes0901a':'EU','abes1785a':'EU','abes1846a':'EU','abes1856a':'EU','abes1859a':'EU','abes1879a':'EU','abes1880a':'EU','abes1881a':'EU','abes1884a':'EU','abes1885a':'EU','abes1886a':'EU','abes1890a':'EU',
                'abes1905a':'EU','abes1907a':'EU','abes1919a':'EU','abes1920a':'EU','abfr09301':'EU','abfr13016':'EU','abfr15048':'EU','abfr41004':'EU','abgb0861a':'EU','abgb0884a':'EU','abgb0885a':'EU','abgb0906a':'EU','abhr0005a':'EU','abhr0009a':'EU',
                'abie0121a':'EU','abie0127a':'EU','abit0706a':'EU','abit0770a':'EU','abit0903a':'EU','abit0934a':'EU','abit0935a':'EU','abit0956a':'EU','abit1035a':'EU','abit1043a':'EU','abit1177a':'EU','abit1186a':'EU','abit1213a':'EU','abit1463a':'EU',
                'abit1486a':'EU','abit1648a':'EU','abit1692a':'EU','abit1790a':'EU','abit1798a':'EU','abit1799a':'EU','abit1826a':'EU','abit1827a':'EU','abit1851a':'EU','abit1852a':'EU','abit1864a':'EU','abit1871a':'EU','abit1872a':'EU','abit1901a':'EU',
                'abit1916a':'EU','abit1918a':'EU','abit1922a':'EU','abit1930a':'EU','abno0062a':'EU','abpl0242a':'EU','abro0096a':'EU','abro0104a':'EU','abro0105a':'EU','abro0107a':'EU','abro0112a':'EU','abro0140a':'EU','abro0151a':'EU','abro0160a':'EU',
                'abro0165a':'EU','abro0188a':'EU','abro0204a':'EU','abro0205a':'EU','abbg0073a':'EU','abcy0004a':'EU','abdebb083':'EU','abes1129a':'EU','abes1741a':'EU','abes1742a':'EU','abes1745a':'EU','abes1759a':'EU','abes1792a':'EU','abes1891a':'EU',
                'abes1892a':'EU','abes1901a':'EU','abes1912a':'EU','abes1914a':'EU','abes1915a':'EU','abes1922a':'EU','abes1953a':'EU','abes1954a':'EU','abes1955a':'EU','abes1957a':'EU','abfi00301':'EU','abfi00781':'EU','abfr03060':'EU','abgb0960a':'EU',
                'abie0118a':'EU','abit0885a':'EU','abit1614a':'EU','abit1624a':'EU','abit1658a':'EU','abit1667a':'EU','abit1668a':'EU','abit1688a':'EU','abit1782a':'EU','abit1817a':'EU','abit1818a':'EU','abit1822a':'EU','abit1906a':'EU','abit1928a':'EU',
                'abit1935a':'EU','abit1953a':'EU','abit1964a':'EU','abit1965a':'EU','abit1967a':'EU','abit1975a':'EU','abmk0045a':'EU','abmt00003':'EU','abpl0487a':'EU','abpt01050':'EU','abro0076a':'EU','abro0092a':'EU','abro0093a':'EU','abro0094a':'EU',
                'abro0095a':'EU','abro0121a':'EU','abro0127a':'EU','abro0150a':'EU','abro0177a':'EU','abro0178a':'EU','abro0195a':'EU','abro0196a':'EU','abse0058a':'EU','abat60179':'EU','abat72912':'EU','abcz0arie':'EU','abdeby122':'EU','abdehb013':'EU',
                'abdesl012':'EU','abdesl020':'EU','abdeth095':'EU','abdeth096':'EU','abes0354a':'EU','abes0360a':'EU','abes0761a':'EU','abes0823a':'EU','abes0873a':'EU','abes1253a':'EU','abes1591a':'EU','abes1645a':'EU','abes1766a':'EU','abes1767a':'EU',
                'abes1768a':'EU','abes1863a':'EU','abes1864a':'EU','abes1910a':'EU','abes1911a':'EU','abes1921a':'EU','abes1924a':'EU','abes1926a':'EU','abes1927a':'EU','abes1939a':'EU','abes1941a':'EU','abes1943a':'EU','abes1945a':'EU','abes1946a':'EU',
                'abes1947a':'EU','abes1952a':'EU','abes1956a':'EU','abes1962a':'EU','abes1963a':'EU','abes1967a':'EU','abes1974a':'EU','abes1975a':'EU','abes1978a':'EU','abes1979a':'EU','abes1980a':'EU','abes1986a':'EU','abes1989a':'EU','abes1991a':'EU',
                'abgb0882a':'EU','abgr0119a':'EU','abie0140a':'EU','abie0141a':'EU','abit1536a':'EU','abit1724a':'EU','abit1843a':'EU','abit1897a':'EU','abit1955a':'EU','abit1961a':'EU','abit1962a':'EU','abit1992a':'EU','abit1993a':'EU','abit1994a':'EU',
                'abit1996a':'EU','abit2003a':'EU','abit2004a':'EU','abit2006a':'EU','abit2018a':'EU','abit2021a':'EU','abmt00008':'EU','abnl00442':'EU','abpl0079a':'EU','abpl0143a':'EU','abpl0283a':'EU','abpl0398a':'EU','abpl0501a':'EU','abpl0504a':'EU',
                'abpl0507a':'EU','abpl0509a':'EU','abpl0520a':'EU','abpt01052':'EU','abro0106a':'EU','abro0108a':'EU','abro0128a':'EU','abro0157a':'EU','abro0163a':'EU','abro0174a':'EU','abro0181a':'EU','abro0194a':'EU','abrs0028a':'EU','abdebw084':'EU',
                'abdebw015':'EU','abdebw006':'EU','abdebw005':'EU','abdebw013':'EU','abdebw007':'EU','abdebw001':'EU','abdebw011':'EU','abdebw038':'EU','abdebe010':'EU','abdebw009':'EU','abdebw023':'EU','abdebw022':'EU','abdebw026':'EU','abdebw002':'EU',
                'abdebw101':'EU','abdebw025':'EU','abdebw014':'EU','abdebw083':'EU','abba0032a':'EU','abch0025a':'EU','abat2aut1':'EU','abat2wo25':'EU','abat60140':'EU','abdebb040':'EU','abdebw113':'EU','abdebb008':'EU','abat60036':'EU','abat2ka71':'EU',
                'abat80807':'EU','abdebw105':'EU','abdebw104':'EU','abat72710':'EU','abat72530':'EU','abdebw074':'EU','abdebw003':'EU','abdebw108':'EU','abbe0312a':'EU','abdebw067':'EU','abdebw051':'EU','abdebw091':'EU','abdebw092':'EU','abdebw058':'EU',
                'abdebw077':'EU','abdebw068':'EU','abcz0avex':'EU','abdebw099':'EU','abat31406':'EU','abdebb035':'EU','abbetr501':'EU','abal0202a':'EU','abal0201a':'EU','abdeby021':'EU','abdebw070':'EU','abdeby058':'EU','abch0029a':'EU','abdeby007':'EU',
                'abdebw075':'EU','abdebw122':'EU','abdebw090':'EU','abdebw086':'EU','abdeby035':'EU','abbetand3':'EU','abat72519':'EU','abdebw089':'EU','abcz0avys':'EU','abderp019':'EU','abdenw053':'EU','abdehb004':'EU','abdesl003':'EU','abdenw071':'EU',
                'abdenw034':'EU','abdenw059':'EU','abdesn053':'EU','abdenw028':'EU','abfr02010':'EU','abes0692a':'EU','abdesh005':'EU','abdeth011':'EU','abdesh007':'EU','abes1565a':'EU','abes1563a':'EU','abes0693a':'EU','abdehe021':'EU','abdesh010':'EU',
                'abdesh012':'EU','abdenw073':'EU','abdenw049':'EU','abfr08716':'EU','abdenw056':'EU','abfr24012':'EU','abdesh021':'EU','abes1398a':'EU','abdehe006':'EU','abdesn038':'EU','abfr22011':'EU','abes1054a':'EU','abdk0032r':'EU','abfr1054a':'EU',
                'abfr12037':'EU','abdenw133':'EU','abdeby076':'EU','abfr17001':'EU','abfr15008':'EU','abes1988a':'EU','abes1865a':'EU','abfr17007':'EU','abdest009':'EU','abdest007':'EU','abfr04132':'EU','abdest027':'EU','abdest001':'EU','abes1136a':'EU',
                'abes1139a':'EU','abes0886a':'EU','abfr19033':'EU','abdest023':'EU','abdesn085':'EU','abfr23129':'EU','abes0918a':'EU','abfr04019':'EU','abes1429a':'EU','abfr03032':'EU','abdesh009':'EU','abfr11035':'EU','abes1388a':'EU','abes1336a':'EU',
                'abfr22007':'EU','abdenw009':'EU','abfr19010':'EU','abfr18029':'EU','abfi00297':'EU','abes2000a':'EU','abes1995a':'EU','abes1970a':'EU','abes1968a':'EU','abes1755a':'EU','abdesh035':'EU','abdeni022':'EU','abdeni021':'EU','abdeni015':'EU',
                'abdeni008':'EU','abes1982a':'EU','abes1981a':'EU','abes1971a':'EU','abes1903a':'EU','abes1718a':'EU','abes1439a':'EU','abfr18053':'EU','abes1969a':'EU','abes1997a':'EU','abes1258a':'EU','abdeni001':'EU','abes1664a':'EU','abdest105':'EU',
                'abes1261a':'EU','abderp031':'EU','abes0915a':'EU','abes2008a':'EU','abes1772a':'EU','abes1764a':'EU','abes1756a':'EU','abdest082':'EU','abdenw143':'EU','abes1960a':'EU','abdemv023':'EU','abes1448a':'EU','abfr1185a':'EU','abes2013a':'EU',
                'abdeth055':'EU','abdehh030':'EU','abes2004a':'EU','abdk0045a':'EU','abfr08019':'EU','abes0193a':'EU','abdeth086':'EU','abes0078a':'EU','abfr03030':'EU','abes1483a':'EU','abes1906a':'EU','abes2003a':'EU','abes1340a':'EU','abes1897a':'EU',
                'abes2020a':'EU','abfr08021':'EU','abes1647a':'EU','abes1020a':'EU','abfr22054':'EU','abdest091':'EU','abfr10027':'EU','abfr14041':'EU','abes1086a':'EU','abdeby090':'EU','abdenw054':'EU','abdehh062':'EU','abdehe038':'EU','abfr04143':'EU',
                'abdeth063':'EU','abdeth064':'EU','abdeth062':'EU','abfi00867':'EU','abes2019a':'EU','abes1983a':'EU','abdest065':'EU','abdeby123':'EU','abes1992a':'EU','abfr04026':'EU','abdesn008':'EU','abes1930a':'EU','abes1929a':'EU','abes1032a':'EU',
                'abdenw174':'EU','abes2016a':'EU','abdeth067':'EU','abdeth065':'EU','abfr04024':'EU','abdeth066':'EU','abdeth073':'EU','abes1168a':'EU','abfr10031':'EU','abdenw180':'EU','abfr08711':'EU','abdeby091':'EU','abes1973a':'EU','abfr04011':'EU',
                'abes2022a':'EU','abes1796a':'EU','abes2025a':'EU','abes1584a':'EU','abfr23138':'EU','abdenw181':'EU','abdesn037':'EU','abes0197a':'EU','abfr15113':'EU','abfr10033':'EU','abes1682a':'EU','abdenw192':'EU','abfr20018':'EU','abfr05063':'EU',
                'abdeby188':'EU','abfr02014':'EU','abfr23053':'EU','abdenw193':'EU','abdenw113':'EU','abdenw182':'EU','abfr15112':'EU','abfr02007':'EU','abdeby187':'EU','abfr23031':'EU','abfr02035':'EU','abes1079a':'EU','abdenw206':'EU','abdenw208':'EU',
                'abdenw204':'EU','abdehe037':'EU','abfr15100':'EU','abfr10010':'EU','abdehh018':'EU','abdehh015':'EU','abdehh012':'EU','abdehh007':'EU','abdehh001':'EU','abdehh019':'EU','abes1294a':'EU','abdehh020':'EU','abdehh014':'EU','abdehh013':'EU',
                'abdehh009':'EU','abdehh004':'EU','abdehh002':'EU','abdehh010':'EU','abfr20030':'EU','abes1312a':'EU','abes0738a':'EU','abes1180a':'EU','abfr12023':'EU','abes1290a':'EU','abfr18051':'EU','abfr25013':'EU','abfr25028':'EU','abfr25033':'EU',
                'abfr28024':'EU','abfr29432':'EU','abfr29436':'EU','abfr29439':'EU','abfr31010':'EU','abfr33233':'EU','abfr33411':'EU','abfr33414':'EU','abfr36008':'EU','abfr36009':'EU','abfr36015':'EU','abfr36016':'EU','abfr36017':'EU','abfr36018':'EU',
                'abfr38004':'EU','abfr38009':'EU','abfr38016':'EU','abgb0204a':'EU','abgb0206a':'EU','abgb0654a':'EU','abgb0655a':'EU','abgb0680a':'EU','abgb0681a':'EU','abgb0736a':'EU','abgb0896a':'EU','abgb1013a':'EU','abgb1019a':'EU','abgb1020a':'EU',
                'abgb1023a':'EU','abgr0009a':'EU','abgr0010a':'EU','abgr0018a':'EU','abgr0019a':'EU','abgr0230a':'EU','abhr0007a':'EU','abhu0002a':'EU','abhu0021a':'EU','abhu0027a':'EU','abhu0028a':'EU','abis0005a':'EU','abis0006a':'EU','abit0186a':'EU',
                'abit0439a':'EU','abit0510a':'EU','abit0522a':'EU','abit0544a':'EU','abit0553a':'EU','abit0590a':'EU','abit0638a':'EU','abit0640a':'EU','abit0680a':'EU','abit0738a':'EU','abit0757a':'EU','abit0760a':'EU','abit0764a':'EU','abit0776a':'EU',
                'abit0781a':'EU','abit0800a':'EU','abit0821a':'EU','abit0825a':'EU','abit0826a':'EU','abit0828a':'EU','abit0829a':'EU','abit0863a':'EU','abit0867a':'EU','abit0873a':'EU','abit0881a':'EU','abit0884a':'EU','abit0888a':'EU','abit0897a':'EU',
                'abit0920a':'EU','abit0921a':'EU','abit0953a':'EU','abit0957a':'EU','abit0966a':'EU','abit0976a':'EU','abit1005a':'EU','abit1011a':'EU','abit1018a':'EU','abit1023a':'EU','abit1041a':'EU','abit1052a':'EU','abit1071a':'EU','abit1084a':'EU',
                'abit1091a':'EU','abit1099a':'EU','abit1103a':'EU','abit1131a':'EU','abit1138a':'EU','abit1141a':'EU','abit1143a':'EU','abit1145a':'EU','abit1147a':'EU','abit1153a':'EU','abit1165a':'EU','abit1170a':'EU','abit1176a':'EU','abit1190a':'EU',
                'abit1194a':'EU','abit1198a':'EU','abit1202a':'EU','abit1204a':'EU','abit1206a':'EU','abit1217a':'EU','abit1221a':'EU','abit1243a':'EU','abit1309a':'EU','abit1310a':'EU','abit1311a':'EU','abit1334a':'EU','abit1339a':'EU','abit1345a':'EU',
                'abit1347a':'EU','abit1355a':'EU','abit1356a':'EU','abit1364a':'EU','abit1365a':'EU','abit1372a':'EU','abit1380a':'EU','abit1385a':'EU','abit1388a':'EU','abit1411a':'EU','abit1414a':'EU','abit1419a':'EU','abit1420a':'EU','abit1431a':'EU',
                'abit1439a':'EU','abit1460a':'EU','abit1491a':'EU','abit1492a':'EU','abit1493a':'EU','abit1495a':'EU','abit1496a':'EU','abit1497a':'EU','abit1498a':'EU','abit1504a':'EU','abit1510a':'EU','abit1533a':'EU','abit1539a':'EU','abit1556a':'EU',
                'abit1564a':'EU','abit1580a':'EU','abit1608a':'EU','abit1611a':'EU','abit1635a':'EU','abit1636a':'EU','abit1637a':'EU','abit1638a':'EU','abit1641a':'EU','abit1645a':'EU','abit1673a':'EU','abit1675a':'EU','abit1679a':'EU','abit1680a':'EU',
                'abit1682a':'EU','abit1683a':'EU','abit1684a':'EU','abit1721a':'EU','abit1727a':'EU','abit1728a':'EU','abit1739a':'EU','abit1740a':'EU','abit1741a':'EU','abit1742a':'EU','abit1746a':'EU','abit1756a':'EU','abit1760a':'EU','abit1765a':'EU',
                'abit1766a':'EU','abit1783a':'EU','abit1803a':'EU','abit1805a':'EU','abit1829a':'EU','abit1830a':'EU','abit1846a':'EU','abit1847a':'EU','abit1856a':'EU','abit1873a':'EU','abit1875a':'EU','abit1878a':'EU','abit1881a':'EU','abit1888a':'EU',
                'abit1889a':'EU','abit1890a':'EU','abit1891a':'EU','abit1895a':'EU','abit1899a':'EU','abit1905a':'EU','abit1908a':'EU','abit1909a':'EU','abit1910a':'EU','abit1926a':'EU','abit1927a':'EU','abit1929a':'EU','abit1933a':'EU','abit1938a':'EU',
                'abit1940a':'EU','abit1952a':'EU','abit1963a':'EU','abit1988a':'EU','abit2002a':'EU','abit2005a':'EU','abit2007a':'EU','abit2012a':'EU','abit2020a':'EU','abit2026a':'EU','abit2028a':'EU','abit2029a':'EU','abit2031a':'EU','abit2033a':'EU',
                'abit2036a':'EU','abit2040a':'EU','abit2042a':'EU','abit2049a':'EU','abit2052a':'EU','abit2055a':'EU','abit2056a':'EU','abit2061a':'EU','abit2075a':'EU','abit2078a':'EU','abit2103a':'EU','ablt00011':'EU','ablt00023':'EU','ablv000d1':'EU',
                'ablv000j1':'EU','ablv000v1':'EU','ablv00rc3':'EU','ablv00rz1':'EU','ablv0rim1':'EU','abme0003a':'EU','abme0004a':'EU','abme0008a':'EU','abme0009a':'EU','abmk0047a':'EU','abmk0048a':'EU','abmt00002':'EU','abmt00006':'EU','abnl00003':'EU',
                'abnl00012':'EU','abnl00014':'EU','abnl00485':'EU','abnl00489':'EU','abnl00493':'EU','abnl00494':'EU','abnl00495':'EU','abnl00496':'EU','abnl00637':'EU','abnl00643':'EU','abnl00701':'EU','abnl00937':'EU','abno0088a':'EU','abpl0015a':'EU',
                'abpl0041a':'EU','abpl0053a':'EU','abpl0064a':'EU','abpl0066a':'EU','abpl0071a':'EU','abpl0072a':'EU','abpl0085a':'EU','abpl0097a':'EU','abpl0104a':'EU','abpl0190a':'EU','abpl0205a':'EU','abpl0224a':'EU','abpl0236a':'EU','abpl0296a':'EU',
                'abpl0312a':'EU','abpl0488a':'EU','abpl0502a':'EU','abpl0518a':'EU','abpl0538a':'EU','abpl0558a':'EU','abpl0559a':'EU','abpl0565a':'EU','abpl0575a':'EU','abpt00119':'EU','abpt01017':'EU','abpt01018':'EU','abpt01020':'EU','abpt01053':'EU',
                'abpt03103':'EU','abpt05010':'EU','abro0062a':'EU','abro0091a':'EU','abro0099a':'EU','abro0100a':'EU','abro0101a':'EU','abro0102a':'EU','abro0109a':'EU','abro0111a':'EU','abro0114a':'EU','abro0116a':'EU','abro0117a':'EU','abro0119a':'EU',
                'abro0120a':'EU','abro0122a':'EU','abro0123a':'EU','abro0132a':'EU','abro0133a':'EU','abro0135a':'EU','abro0136a':'EU','abro0137a':'EU','abro0139a':'EU','abro0142a':'EU','abro0143a':'EU','abro0144a':'EU','abro0145a':'EU','abro0154a':'EU',
                'abro0155a':'EU','abro0159a':'EU','abro0162a':'EU','abro0166a':'EU','abro0168a':'EU','abro0171a':'EU','abro0172a':'EU','abro0176a':'EU','abro0180a':'EU','abro0182a':'EU','abro0184a':'EU','abro0185a':'EU','abro0186a':'EU','abro0187a':'EU',
                'abro0189a':'EU','abro0192a':'EU','abro0201a':'EU','abro0203a':'EU','abrs0029a':'EU','abrs0032a':'EU','abrs0036a':'EU','abrs0037a':'EU','abse0021a':'EU','abse0028a':'EU','absk0023a':'EU','absk0039a':'EU','absk0043a':'EU','absk0049a':'EU',
                'abme0002a':'EU','abrs0042a':'EU','abrs0047a':'EU',
                #ozonesonde sites
                #------------------------
                #NDACC
                #SHADOZ
                #WOUDC
                'stn012':'AS','stn014':'AS','stn018':'ARC','stn021':'NA','stn024':'ARC','stn029':'O','stn043':'O','stn053':'EU','stn076':'NA','stn077':'NA','stn089':'ARC','stn099':'EU','stn101':'ANT','stn107':'NA','stn109':'O', \
                'stn156':'EU','stn175':'AF','stn191':'O','stn233':'ANT','stn308':'EU','stn315':'ARC','stn316':'EU','stn318':'EU','stn323':'ANT','stn330':'AS','stn338':'NA','stn348':'EU','stn394':'OC','stn435':'SA','stn436':'O', \
                'stn437':'OC','stn443':'OC','stn456':'NA','stn457':'NA','stn458':'NA','ascen':'O','costarica':'SA','hanoi':'AS','hilo':'O','java':'OC','kuala':'OC','nairobi':'AF','natal':'SA','paramaribo':'SA','reunion':'O', \
                'samoa':'O'}
                
    for i in range(len(obs_refs)):
        if '_' in obs_refs[i]:
            split_res = obs_refs[i].split('_')
            obs_refs[i] = split_res[0]
            
    inv = []
    for i in obs_refs:
        if i not in tag_dict:
            inv.append(i)

    print ("\n%s\n"%("':'EU','".join(i for i in inv)))
        
    tags = [tag_dict[key] for key in obs_refs]
    tags = np.array(tags)		

    return tags
    
def get_spec_specs(species):
    if (species == 'ISOP'):
        data_resolution = 0.21
        mol_mass = 68.11702
        aqs_code =  43243 
        airbase_code = '00451'

    elif (species == 'O3'):
        data_resolution = 1.1
        mol_mass = 47.9982
        aqs_code = 44201
        airbase_code = '00007'

    elif (species == 'NO'):
        data_resolution = 0.51
        mol_mass = 30.00614
        aqs_code = 42601
        airbase_code = '00038'
    
    elif (species == 'NO2'):
        data_resolution = 1.1
        mol_mass = 46.00554
        aqs_code = 42602
        airbase_code = '00008'
        
    elif species == 'CO':
        data_resolution = 10.1
        mol_mass = 28.0101
        aqs_code = 42101
        airbase_code = '00010'
        
    return data_resolution,mol_mass,aqs_code,airbase_code
  
def get_periodic_specific(obs_fname,model_fname):
    #set filenames as 'na' if don't want to read
    
    if obs_fname != 'na':
        #######
        #read in obs.
        obs_root = Dataset(obs_fname)

        obs_diurnal_mag = obs_root.variables['diurnal_amplitude'][:]
        obs_diurnal_ph = obs_root.variables['diurnal_phase'][:]
        obs_seasonal_mag = obs_root.variables['seasonal_amplitude'][:]
        obs_seasonal_ph = obs_root.variables['seasonal_phase'][:]
        obs_mean = obs_root.variables['average'][:]
        obs_p1 = obs_root.variables['percentile1'][:]
        obs_p5 = obs_root.variables['percentile5'][:]
        obs_p25 = obs_root.variables['percentile25'][:]
        obs_p50 = obs_root.variables['median'][:]
        obs_p75 = obs_root.variables['percentile75'][:]
        obs_p95 = obs_root.variables['percentile95'][:]
        obs_p99 = obs_root.variables['percentile99'][:]
        obs_diurnal_ave_waveform = obs_root.variables['diurnal_average_waveform'][:]
        #obs_diurnal_ave_waveform_extended = obs_root.variables['diurnal_average_waveform_extended'][:]
        #obs_diurnal_season_waveform_extended = obs_root.variables['diurnal_season_waveform_extended'][:]
        obs_seasonal_waveform = obs_root.variables['seasonal_waveform'][:]
        #obs_seasonal_waveform_extended = obs_root.variables['seasonal_waveform_extended'][:]
        obs_full_ave_waveform = obs_root.variables['full_average_waveform'][:]
        #obs_full_season_waveform = obs_root.variables['full_season_waveform'][:]
        obs_pc_var_daily = obs_root.variables['diurnal_periodic_variability'][:]
        obs_pc_var_seasonal = obs_root.variables['seasonal_periodic_variability'][:]
        obs_pc_var_full = obs_root.variables['full_periodic_variability'][:]
        obs_pc_var_noise = obs_root.variables['noise_periodic_variability'][:]
        obs_total_var = obs_root.variables['total_variability'][:]
        obs_diurnal_mag_spring = obs_root.variables['diurnal_amplitude_spring'][:]
        obs_diurnal_ph_spring = obs_root.variables['diurnal_phase_spring'][:]
        obs_mean_spring = obs_root.variables['average_spring'][:]
        obs_p1_spring = obs_root.variables['percentile1_spring'][:]
        obs_p5_spring = obs_root.variables['percentile5_spring'][:]
        obs_p25_spring = obs_root.variables['percentile25_spring'][:]
        obs_p50_spring = obs_root.variables['median_spring'][:]
        obs_p75_spring = obs_root.variables['percentile75_spring'][:]
        obs_p95_spring = obs_root.variables['percentile95_spring'][:]
        obs_p99_spring = obs_root.variables['percentile99_spring'][:]
        obs_diurnal_waveform_spring = obs_root.variables['diurnal_waveform_spring'][:]
        obs_diurnal_mag_summer = obs_root.variables['diurnal_amplitude_summer'][:]
        obs_diurnal_ph_summer = obs_root.variables['diurnal_phase_summer'][:]
        obs_mean_summer = obs_root.variables['average_summer'][:]
        obs_p1_summer = obs_root.variables['percentile1_summer'][:]
        obs_p5_summer = obs_root.variables['percentile5_summer'][:]
        obs_p25_summer = obs_root.variables['percentile25_summer'][:]
        obs_p50_summer = obs_root.variables['median_summer'][:]
        obs_p75_summer = obs_root.variables['percentile75_summer'][:]
        obs_p95_summer = obs_root.variables['percentile95_summer'][:]
        obs_p99_summer = obs_root.variables['percentile99_summer'][:]
        obs_diurnal_waveform_summer = obs_root.variables['diurnal_waveform_summer'][:]
        obs_diurnal_mag_autumn = obs_root.variables['diurnal_amplitude_autumn'][:]
        obs_diurnal_ph_autumn = obs_root.variables['diurnal_phase_autumn'][:]
        obs_mean_autumn = obs_root.variables['average_autumn'][:]
        obs_p1_autumn = obs_root.variables['percentile1_autumn'][:]
        obs_p5_autumn = obs_root.variables['percentile5_autumn'][:]
        obs_p25_autumn = obs_root.variables['percentile25_autumn'][:]
        obs_p50_autumn = obs_root.variables['median_autumn'][:]
        obs_p75_autumn = obs_root.variables['percentile75_autumn'][:]
        obs_p95_autumn = obs_root.variables['percentile95_autumn'][:]
        obs_p99_autumn = obs_root.variables['percentile99_autumn'][:]
        obs_diurnal_waveform_autumn = obs_root.variables['diurnal_waveform_autumn'][:]
        obs_diurnal_mag_winter = obs_root.variables['diurnal_amplitude_winter'][:]
        obs_diurnal_ph_winter = obs_root.variables['diurnal_phase_winter'][:]
        obs_mean_winter = obs_root.variables['average_winter'][:]
        obs_p1_winter = obs_root.variables['percentile1_winter'][:]
        obs_p5_winter = obs_root.variables['percentile5_winter'][:]
        obs_p25_winter = obs_root.variables['percentile25_winter'][:]
        obs_p50_winter = obs_root.variables['median_winter'][:]
        obs_p75_winter = obs_root.variables['percentile75_winter'][:]
        obs_p95_winter = obs_root.variables['percentile95_winter'][:]
        obs_p99_winter = obs_root.variables['percentile99_winter'][:]
        obs_diurnal_waveform_winter = obs_root.variables['diurnal_waveform_winter'][:]
        obs_seasonal_mag_day = obs_root.variables['seasonal_amplitude_day'][:]
        obs_seasonal_ph_day = obs_root.variables['seasonal_phase_day'][:]
        obs_mean_day = obs_root.variables['average_day'][:]
        obs_p1_day = obs_root.variables['percentile1_day'][:]
        obs_p5_day = obs_root.variables['percentile5_day'][:]
        obs_p25_day = obs_root.variables['percentile25_day'][:]
        obs_p50_day = obs_root.variables['median_day'][:]
        obs_p75_day = obs_root.variables['percentile75_day'][:]
        obs_p95_day = obs_root.variables['percentile95_day'][:]
        obs_p99_day = obs_root.variables['percentile99_day'][:]
        obs_seasonal_waveform_day = obs_root.variables['seasonal_waveform_day'][:]
        obs_seasonal_mag_night = obs_root.variables['seasonal_amplitude_night'][:]
        obs_seasonal_ph_night = obs_root.variables['seasonal_phase_night'][:]
        obs_mean_night = obs_root.variables['average_night'][:]
        obs_p1_night = obs_root.variables['percentile1_night'][:]
        obs_p5_night = obs_root.variables['percentile5_night'][:]
        obs_p25_night = obs_root.variables['percentile25_night'][:]
        obs_p50_night = obs_root.variables['median_night'][:]
        obs_p75_night = obs_root.variables['percentile75_night'][:]
        obs_p95_night = obs_root.variables['percentile95_night'][:]
        obs_p99_night = obs_root.variables['percentile99_night'][:]
        obs_seasonal_waveform_night = obs_root.variables['seasonal_waveform_night'][:]
        obs_daily_h11_mag = obs_root.variables['diurnal_amplitude_harmonic12'][:]
        obs_daily_h10_mag = obs_root.variables['diurnal_amplitude_harmonic11'][:]
        obs_daily_h9_mag = obs_root.variables['diurnal_amplitude_harmonic10'][:]
        obs_daily_h8_mag = obs_root.variables['diurnal_amplitude_harmonic9'][:]
        obs_daily_h7_mag = obs_root.variables['diurnal_amplitude_harmonic8'][:]
        obs_daily_h6_mag = obs_root.variables['diurnal_amplitude_harmonic7'][:]
        obs_daily_h5_mag = obs_root.variables['diurnal_amplitude_harmonic6'][:]
        obs_daily_h4_mag = obs_root.variables['diurnal_amplitude_harmonic5'][:]
        obs_daily_h3_mag = obs_root.variables['diurnal_amplitude_harmonic4'][:]
        obs_daily_h2_mag = obs_root.variables['diurnal_amplitude_harmonic3'][:]
        obs_daily_h1_mag = obs_root.variables['diurnal_amplitude_harmonic2'][:]
        obs_daily_mag = obs_root.variables['diurnal_amplitude_harmonic1'][:]
        #obs_annual_h11_mag = obs_root.variables['seasonal_amplitude_harmonic12'][:]
        #obs_annual_h10_mag = obs_root.variables['seasonal_amplitude_harmonic11'][:]
        #obs_annual_h9_mag = obs_root.variables['seasonal_amplitude_harmonic10'][:]
        #obs_annual_h8_mag = obs_root.variables['seasonal_amplitude_harmonic9'][:]
        #obs_annual_h7_mag = obs_root.variables['seasonal_amplitude_harmonic8'][:]
        #obs_annual_h6_mag = obs_root.variables['seasonal_amplitude_harmonic7'][:]
        #obs_annual_h5_mag = obs_root.variables['seasonal_amplitude_harmonic6'][:]
        #obs_annual_h4_mag = obs_root.variables['seasonal_amplitude_harmonic5'][:]
        obs_annual_h3_mag = obs_root.variables['seasonal_amplitude_harmonic4'][:]
        obs_annual_h2_mag = obs_root.variables['seasonal_amplitude_harmonic3'][:]
        obs_annual_h1_mag = obs_root.variables['seasonal_amplitude_harmonic2'][:]
        obs_annual_mag = obs_root.variables['seasonal_amplitude_harmonic1'][:]

    if model_fname != 'na':
        ######
        #read in model
        model_root = Dataset(model_fname)
    
        model_diurnal_mag = model_root.variables['diurnal_amplitude'][:]
        model_diurnal_ph = model_root.variables['diurnal_phase'][:]
        model_seasonal_mag = model_root.variables['seasonal_amplitude'][:]
        model_seasonal_ph = model_root.variables['seasonal_phase'][:]
        model_mean = model_root.variables['average'][:]
        model_p1 = model_root.variables['percentile1'][:]
        model_p5 = model_root.variables['percentile5'][:]
        model_p25 = model_root.variables['percentile25'][:]
        model_p50 = model_root.variables['median'][:]
        model_p75 = model_root.variables['percentile75'][:]
        model_p95 = model_root.variables['percentile95'][:]
        model_p99 = model_root.variables['percentile99'][:]
        model_diurnal_ave_waveform = model_root.variables['diurnal_average_waveform'][:]
        #model_diurnal_ave_waveform_extended = model_root.variables['diurnal_average_waveform_extended'][:]
        #model_diurnal_season_waveform_extended = model_root.variables['diurnal_season_waveform_extended'][:]
        model_seasonal_waveform = model_root.variables['seasonal_waveform'][:]
        #model_seasonal_waveform_extended = model_root.variables['seasonal_waveform_extended'][:]
        model_full_ave_waveform = model_root.variables['full_average_waveform'][:]
       # model_full_season_waveform = model_root.variables['full_season_waveform'][:]
        model_pc_var_daily = model_root.variables['diurnal_periodic_variability'][:]
        model_pc_var_seasonal = model_root.variables['seasonal_periodic_variability'][:]
        model_pc_var_full = model_root.variables['full_periodic_variability'][:]
        model_pc_var_noise = model_root.variables['noise_periodic_variability'][:]
        model_total_var = model_root.variables['total_variability'][:]
        model_diurnal_mag_spring = model_root.variables['diurnal_amplitude_spring'][:]
        model_diurnal_ph_spring = model_root.variables['diurnal_phase_spring'][:]
        model_mean_spring = model_root.variables['average_spring'][:]
        model_p1_spring = model_root.variables['percentile1_spring'][:]
        model_p5_spring = model_root.variables['percentile5_spring'][:]
        model_p25_spring = model_root.variables['percentile25_spring'][:]
        model_p50_spring = model_root.variables['median_spring'][:]
        model_p75_spring = model_root.variables['percentile75_spring'][:]
        model_p95_spring = model_root.variables['percentile95_spring'][:]
        model_p99_spring = model_root.variables['percentile99_spring'][:]
        model_diurnal_waveform_spring = model_root.variables['diurnal_waveform_spring'][:]
        model_diurnal_mag_summer = model_root.variables['diurnal_amplitude_summer'][:]
        model_diurnal_ph_summer = model_root.variables['diurnal_phase_summer'][:]
        model_mean_summer = model_root.variables['average_summer'][:]
        model_p1_summer = model_root.variables['percentile1_summer'][:]
        model_p5_summer = model_root.variables['percentile5_summer'][:]
        model_p25_summer = model_root.variables['percentile25_summer'][:]
        model_p50_summer = model_root.variables['median_summer'][:]
        model_p75_summer = model_root.variables['percentile75_summer'][:]
        model_p95_summer = model_root.variables['percentile95_summer'][:]
        model_p99_summer = model_root.variables['percentile99_summer'][:]
        model_diurnal_waveform_summer = model_root.variables['diurnal_waveform_summer'][:]
        model_diurnal_mag_autumn = model_root.variables['diurnal_amplitude_autumn'][:]
        model_diurnal_ph_autumn = model_root.variables['diurnal_phase_autumn'][:]
        model_mean_autumn = model_root.variables['average_autumn'][:]
        model_p1_autumn = model_root.variables['percentile1_autumn'][:]
        model_p5_autumn = model_root.variables['percentile5_autumn'][:]
        model_p25_autumn = model_root.variables['percentile25_autumn'][:]
        model_p50_autumn = model_root.variables['median_autumn'][:]
        model_p75_autumn = model_root.variables['percentile75_autumn'][:]
        model_p95_autumn = model_root.variables['percentile95_autumn'][:]
        model_p99_autumn = model_root.variables['percentile99_autumn'][:]
        model_diurnal_waveform_autumn = model_root.variables['diurnal_waveform_autumn'][:]
        model_diurnal_mag_winter = model_root.variables['diurnal_amplitude_winter'][:]
        model_diurnal_ph_winter = model_root.variables['diurnal_phase_winter'][:]
        model_mean_winter = model_root.variables['average_winter'][:]
        model_p1_winter = model_root.variables['percentile1_winter'][:]
        model_p5_winter = model_root.variables['percentile5_winter'][:]
        model_p25_winter = model_root.variables['percentile25_winter'][:]
        model_p50_winter = model_root.variables['median_winter'][:]
        model_p75_winter = model_root.variables['percentile75_winter'][:]
        model_p95_winter = model_root.variables['percentile95_winter'][:]
        model_p99_winter = model_root.variables['percentile99_winter'][:]
        model_diurnal_waveform_winter = model_root.variables['diurnal_waveform_winter'][:]
        model_seasonal_mag_day = model_root.variables['seasonal_amplitude_day'][:]
        model_seasonal_ph_day = model_root.variables['seasonal_phase_day'][:]
        model_mean_day = model_root.variables['average_day'][:]
        model_p1_day = model_root.variables['percentile1_day'][:]
        model_p5_day = model_root.variables['percentile5_day'][:]
        model_p25_day = model_root.variables['percentile25_day'][:]
        model_p50_day = model_root.variables['median_day'][:]
        model_p75_day = model_root.variables['percentile75_day'][:]
        model_p95_day = model_root.variables['percentile95_day'][:]
        model_p99_day = model_root.variables['percentile99_day'][:]
        model_seasonal_waveform_day = model_root.variables['seasonal_waveform_day'][:]
        model_seasonal_mag_night = model_root.variables['seasonal_amplitude_night'][:]
        model_seasonal_ph_night = model_root.variables['seasonal_phase_night'][:]
        model_mean_night = model_root.variables['average_night'][:]
        model_p1_night = model_root.variables['percentile1_night'][:]
        model_p5_night = model_root.variables['percentile5_night'][:]
        model_p25_night = model_root.variables['percentile25_night'][:]
        model_p50_night = model_root.variables['median_night'][:]
        model_p75_night = model_root.variables['percentile75_night'][:]
        model_p95_night = model_root.variables['percentile95_night'][:]
        model_p99_night = model_root.variables['percentile99_night'][:]
        model_seasonal_waveform_night = model_root.variables['seasonal_waveform_night'][:]
        
        
        model_daily_h11_mag = model_root.variables['diurnal_amplitude_harmonic12'][:]
        model_daily_h10_mag = model_root.variables['diurnal_amplitude_harmonic11'][:]
        model_daily_h9_mag = model_root.variables['diurnal_amplitude_harmonic10'][:]
        model_daily_h8_mag = model_root.variables['diurnal_amplitude_harmonic9'][:]
        model_daily_h7_mag = model_root.variables['diurnal_amplitude_harmonic8'][:]
        model_daily_h6_mag = model_root.variables['diurnal_amplitude_harmonic7'][:]
        model_daily_h5_mag = model_root.variables['diurnal_amplitude_harmonic6'][:]
        model_daily_h4_mag = model_root.variables['diurnal_amplitude_harmonic5'][:]
        model_daily_h3_mag = model_root.variables['diurnal_amplitude_harmonic4'][:]
        model_daily_h2_mag = model_root.variables['diurnal_amplitude_harmonic3'][:]
        model_daily_h1_mag = model_root.variables['diurnal_amplitude_harmonic2'][:]
        model_daily_mag = model_root.variables['diurnal_amplitude_harmonic1'][:]
        #model_annual_h11_mag = model_root.variables['seasonal_amplitude_harmonic12'][:]
        #model_annual_h10_mag = model_root.variables['seasonal_amplitude_harmonic11'][:]
        #model_annual_h9_mag = model_root.variables['seasonal_amplitude_harmonic10'][:]
        #model_annual_h8_mag = model_root.variables['seasonal_amplitude_harmonic9'][:]
        #model_annual_h7_mag = model_root.variables['seasonal_amplitude_harmonic8'][:]
        #model_annual_h6_mag = model_root.variables['seasonal_amplitude_harmonic7'][:]
        #model_annual_h5_mag = model_root.variables['seasonal_amplitude_harmonic6'][:]
        #model_annual_h4_mag = model_root.variables['seasonal_amplitude_harmonic5'][:]
        model_annual_h3_mag = model_root.variables['seasonal_amplitude_harmonic4'][:]
        model_annual_h2_mag = model_root.variables['seasonal_amplitude_harmonic3'][:]
        model_annual_h1_mag = model_root.variables['seasonal_amplitude_harmonic2'][:]
        model_annual_mag = model_root.variables['seasonal_amplitude_harmonic1'][:]

    #only output model files
    if obs_fname == 'na':
        return (#model_diurnal_mag,model_diurnal_ph,model_seasonal_mag,model_seasonal_ph,model_mean,model_p1,model_p5,model_p25,model_p50,model_p75,model_p95,model_p99,model_diurnal_ave_waveform,model_seasonal_waveform,model_full_ave_waveform,model_pc_var_daily,model_pc_var_seasonal,model_pc_var_full,model_pc_var_noise,model_total_var,
                #model_diurnal_mag_spring,model_diurnal_ph_spring,model_mean_spring,model_p1_spring,model_p5_spring,model_p25_spring,model_p50_spring,model_p75_spring,model_p95_spring,model_p99_spring,model_diurnal_waveform_spring,
                #model_diurnal_mag_summer,model_diurnal_ph_summer,model_mean_summer,model_p1_summer,model_p5_summer,model_p25_summer,model_p50_summer,model_p75_summer,model_p95_summer,model_p99_summer,model_diurnal_waveform_summer, 
                #model_diurnal_mag_autumn,model_diurnal_ph_autumn,model_mean_autumn,model_p1_autumn,model_p5_autumn,model_p25_autumn,model_p50_autumn,model_p75_autumn,model_p95_autumn,model_p99_autumn,model_diurnal_waveform_autumn, 
                #model_diurnal_mag_winter,model_diurnal_ph_winter,model_mean_winter,model_p1_winter,model_p5_winter,model_p25_winter,model_p50_winter,model_p75_winter,model_p95_winter,model_p99_winter,model_diurnal_waveform_winter,
                #model_seasonal_mag_day,model_seasonal_ph_day,model_mean_day,model_p1_day,model_p5_day,model_p25_day,model_p50_day,model_p75_day,model_p95_day,model_p99_day,model_seasonal_waveform_day,
                #model_seasonal_mag_night,model_seasonal_ph_night,model_mean_night,model_p1_night,model_p5_night,model_p25_night,model_p50_night,model_p75_night,model_p95_night,model_p99_night,model_seasonal_waveform_night,
                #model_daily_h11_mag,model_daily_h10_mag,model_daily_h9_mag,model_daily_h8_mag,model_daily_h7_mag,model_daily_h6_mag,model_daily_h5_mag,model_daily_h4_mag,model_daily_h3_mag,model_daily_h2_mag,model_daily_h1_mag,model_daily_mag,model_annual_h3_mag,model_annual_h2_mag,model_annual_h1_mag,model_annual_mag)
    			model_diurnal_mag,model_diurnal_ph,model_seasonal_mag,model_seasonal_ph,model_mean,model_p1,model_p5,model_p25,model_p50,model_p75,model_p95,model_p99,model_diurnal_ave_waveform,model_seasonal_waveform,model_full_ave_waveform,model_pc_var_daily,model_pc_var_seasonal,model_pc_var_full,model_pc_var_noise,model_total_var,
                model_diurnal_mag_spring,model_diurnal_ph_spring,model_mean_spring,model_p1_spring,model_p5_spring,model_p25_spring,model_p50_spring,model_p75_spring,model_p95_spring,model_p99_spring,model_diurnal_waveform_spring,
                model_diurnal_mag_summer,model_diurnal_ph_summer,model_mean_summer,model_p1_summer,model_p5_summer,model_p25_summer,model_p50_summer,model_p75_summer,model_p95_summer,model_p99_summer,model_diurnal_waveform_summer, 
                model_diurnal_mag_autumn,model_diurnal_ph_autumn,model_mean_autumn,model_p1_autumn,model_p5_autumn,model_p25_autumn,model_p50_autumn,model_p75_autumn,model_p95_autumn,model_p99_autumn,model_diurnal_waveform_autumn, 
                model_diurnal_mag_winter,model_diurnal_ph_winter,model_mean_winter,model_p1_winter,model_p5_winter,model_p25_winter,model_p50_winter,model_p75_winter,model_p95_winter,model_p99_winter,model_diurnal_waveform_winter,
                model_seasonal_mag_day,model_seasonal_ph_day,model_mean_day,model_p1_day,model_p5_day,model_p25_day,model_p50_day,model_p75_day,model_p95_day,model_p99_day,model_seasonal_waveform_day,
                model_seasonal_mag_night,model_seasonal_ph_night,model_mean_night,model_p1_night,model_p5_night,model_p25_night,model_p50_night,model_p75_night,model_p95_night,model_p99_night,model_seasonal_waveform_night,
                model_daily_h3_mag,model_daily_h2_mag,model_daily_h1_mag,model_daily_mag,model_annual_h3_mag,model_annual_h2_mag,model_annual_h1_mag,model_annual_mag)
    
    
    #only output obs files
    elif model_fname == 'na':
        return (obs_diurnal_mag,obs_diurnal_ph,obs_seasonal_mag,obs_seasonal_ph,obs_mean,obs_p1,obs_p5,obs_p25,obs_p50,obs_p75,obs_p95,obs_p99,obs_diurnal_ave_waveform,obs_seasonal_waveform,obs_full_ave_waveform,obs_pc_var_daily,obs_pc_var_seasonal,obs_pc_var_full,obs_pc_var_noise,obs_total_var,
                obs_diurnal_mag_spring,obs_diurnal_ph_spring,obs_mean_spring,obs_p1_spring,obs_p5_spring,obs_p25_spring,obs_p50_spring,obs_p75_spring,obs_p95_spring,obs_p99_spring,obs_diurnal_waveform_spring,
                obs_diurnal_mag_summer,obs_diurnal_ph_summer,obs_mean_summer,obs_p1_summer,obs_p5_summer,obs_p25_summer,obs_p50_summer,obs_p75_summer,obs_p95_summer,obs_p99_summer,obs_diurnal_waveform_summer,
                obs_diurnal_mag_autumn,obs_diurnal_ph_autumn,obs_mean_autumn,obs_p1_autumn,obs_p5_autumn,obs_p25_autumn,obs_p50_autumn,obs_p75_autumn,obs_p95_autumn,obs_p99_autumn,obs_diurnal_waveform_autumn,
                obs_diurnal_mag_winter,obs_diurnal_ph_winter,obs_mean_winter,obs_p1_winter,obs_p5_winter,obs_p25_winter,obs_p50_winter,obs_p75_winter,obs_p95_winter,obs_p99_winter,obs_diurnal_waveform_winter,
                obs_seasonal_mag_day,obs_seasonal_ph_day,obs_mean_day,obs_p1_day,obs_p5_day,obs_p25_day,obs_p50_day,obs_p75_day,obs_p95_day,obs_p99_day,obs_seasonal_waveform_day,
                obs_seasonal_mag_night,obs_seasonal_ph_night,obs_mean_night,obs_p1_night,obs_p5_night,obs_p25_night,obs_p50_night,obs_p75_night,obs_p95_night,obs_p99_night,obs_seasonal_waveform_night,
                obs_daily_h3_mag,obs_daily_h2_mag,obs_daily_h1_mag,obs_daily_mag,obs_annual_h3_mag,obs_annual_h2_mag,obs_annual_h1_mag,obs_annual_mag)
    
    #output both files
    else:
        return (obs_diurnal_mag,obs_diurnal_ph,obs_seasonal_mag,obs_seasonal_ph,obs_mean,obs_p1,obs_p5,obs_p25,obs_p50,obs_p75,obs_p95,obs_p99,obs_diurnal_ave_waveform,obs_seasonal_waveform,obs_full_ave_waveform,obs_pc_var_daily,obs_pc_var_seasonal,obs_pc_var_full,obs_pc_var_noise,obs_total_var,
                obs_diurnal_mag_spring,obs_diurnal_ph_spring,obs_mean_spring,obs_p1_spring,obs_p5_spring,obs_p25_spring,obs_p50_spring,obs_p75_spring,obs_p95_spring,obs_p99_spring,obs_diurnal_waveform_spring,
                obs_diurnal_mag_summer,obs_diurnal_ph_summer,obs_mean_summer,obs_p1_summer,obs_p5_summer,obs_p25_summer,obs_p50_summer,obs_p75_summer,obs_p95_summer,obs_p99_summer,obs_diurnal_waveform_summer,
                obs_diurnal_mag_autumn,obs_diurnal_ph_autumn,obs_mean_autumn,obs_p1_autumn,obs_p5_autumn,obs_p25_autumn,obs_p50_autumn,obs_p75_autumn,obs_p95_autumn,obs_p99_autumn,obs_diurnal_waveform_autumn,
                obs_diurnal_mag_winter,obs_diurnal_ph_winter,obs_mean_winter,obs_p1_winter,obs_p5_winter,obs_p25_winter,obs_p50_winter,obs_p75_winter,obs_p95_winter,obs_p99_winter,obs_diurnal_waveform_winter,
                obs_seasonal_mag_day,obs_seasonal_ph_day,obs_mean_day,obs_p1_day,obs_p5_day,obs_p25_day,obs_p50_day,obs_p75_day,obs_p95_day,obs_p99_day,obs_seasonal_waveform_day,
                obs_seasonal_mag_night,obs_seasonal_ph_night,obs_mean_night,obs_p1_night,obs_p5_night,obs_p25_night,obs_p50_night,obs_p75_night,obs_p95_night,obs_p99_night,obs_seasonal_waveform_night,
                obs_daily_h3_mag,obs_daily_h2_mag,obs_daily_h1_mag,obs_daily_mag,obs_annual_h3_mag,obs_annual_h2_mag,obs_annual_h1_mag,obs_annual_mag,
                model_diurnal_mag,model_diurnal_ph,model_seasonal_mag,model_seasonal_ph,model_mean,model_p1,model_p5,model_p25,model_p50,model_p75,model_p95,model_p99,model_diurnal_ave_waveform,model_seasonal_waveform,model_full_ave_waveform,model_pc_var_daily,model_pc_var_seasonal,model_pc_var_full,model_pc_var_noise,model_total_var,
                model_diurnal_mag_spring,model_diurnal_ph_spring,model_mean_spring,model_p1_spring,model_p5_spring,model_p25_spring,model_p50_spring,model_p75_spring,model_p95_spring,model_p99_spring,model_diurnal_waveform_spring,
                model_diurnal_mag_summer,model_diurnal_ph_summer,model_mean_summer,model_p1_summer,model_p5_summer,model_p25_summer,model_p50_summer,model_p75_summer,model_p95_summer,model_p99_summer,model_diurnal_waveform_summer, 
                model_diurnal_mag_autumn,model_diurnal_ph_autumn,model_mean_autumn,model_p1_autumn,model_p5_autumn,model_p25_autumn,model_p50_autumn,model_p75_autumn,model_p95_autumn,model_p99_autumn,model_diurnal_waveform_autumn, 
                model_diurnal_mag_winter,model_diurnal_ph_winter,model_mean_winter,model_p1_winter,model_p5_winter,model_p25_winter,model_p50_winter,model_p75_winter,model_p95_winter,model_p99_winter,model_diurnal_waveform_winter,
                model_seasonal_mag_day,model_seasonal_ph_day,model_mean_day,model_p1_day,model_p5_day,model_p25_day,model_p50_day,model_p75_day,model_p95_day,model_p99_day,model_seasonal_waveform_day,
                model_seasonal_mag_night,model_seasonal_ph_night,model_mean_night,model_p1_night,model_p5_night,model_p25_night,model_p50_night,model_p75_night,model_p95_night,model_p99_night,model_seasonal_waveform_night,
                model_daily_h3_mag,model_daily_h2_mag,model_daily_h1_mag,model_daily_mag,model_annual_h3_mag,model_annual_h2_mag,model_annual_h1_mag,model_annual_mag)
          
    
#get percent periodic variance of total time series variance 
#(accounts for gaps by including gaps from raw time series in periodic waveform)

def periodic_variance_percent(raw_ts,full_periodic_waveform,full_periodic_amps,cut_periodic_amps,valid):
    #calculate total periodic variance from all periodic amplitudes
    #time domain variance is related to amplitudes by **2/2.
    full_periodic_var = np.sum([i**2/2. for i in full_periodic_amps])
    
    #calculate total periodic variance from cut of periodic amplitudes
    #time domain variance is related to amplitudes by **2/2.
    cut_periodic_var = np.sum([i**2/2. for i in cut_periodic_amps])
    
    #put gaps into full periodic waveform
    full_periodic_waveform = full_periodic_waveform[valid]
    
    #subtract full periodic waveform from raw time series 
    #remianing in just noise of raw time series
    noise_waveform = raw_ts - full_periodic_waveform
    
    #get variance of noise time series
    noise_var = np.var(noise_waveform)
    
    #get total variance of time series - done this way as LSP does follow Parseval's theorem exactly.
    #but are confident of periodic cycle variance, as done work to ensure this accuracy
    #also if took total variance in standard way of time series, can fluctuate due to its length
    #this way aims to get all % to add as close to 100% as possible
    #add periodic variance to noise variance
    total_var = full_periodic_var+noise_var
    
    #get percentage of cut periodic variance in total variance
    periodic_var_pc = (100./total_var)*cut_periodic_var
    #get percentage of noise variance in total variance
    noise_var_pc = (100./total_var)*noise_var
 
    return periodic_var_pc,noise_var_pc,total_var
 
#collection of modules that naturally sort (human intuitive)
#natsorted is module to call to run sort

def try_int(s):
    "Convert to integer if possible."
    try: return int(s)
    except: return s

def natsort_key(s):
    "Used internally to get a tuple by which s is sorted."
    import re
    return map(try_int, re.findall(r'(\d+|\D+)', s))

def natcmp(a, b):
    "Natural string comparison, case sensitive."
    return cmp(natsort_key(a), natsort_key(b))

def natcasecmp(a, b):
    "Natural string comparison, ignores case."
    return natcmp(a.lower(), b.lower())

def natsort(seq, cmp=natcmp):
    "In-place natural string sort."
    seq.sort(cmp)
    
def natsorted(seq, cmp=natcmp):
    "Returns a copy of seq, sorted by natural string sort."
    import copy
    temp = copy.copy(seq)
    natsort(temp, cmp)
    return temp

#module that takes date/time ints and outputs a relative time from a start year in hours. 

def date_process(date,time,start_year):
	year=(date//10000)
	month=((date-year*10000)//100)
	day=(date-year*10000-month*100)

	hour=time//100
	min=(time-hour*100)

	doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(int(start_year),1,1,0,0,0) \
              for i in range(len(year))]

	processed_dates=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
	
	return np.array(processed_dates)

#-------------------------------------------------------------
#-------------------------------------------------------------
#SECTION 2 - MODEL PROCESSING

#2.01 - get model filename (from directory structure)
#2.02 - read model data from netcdf for 1 site
#2.03 - read model data from netcdf for all sites
#2.04 - quality checks for _ALL tagged observational data
#2.05 - quality checks for _NR tagged observational data
#2.06 - quality checks for _PERIODIC tagged observational data
#2.07 - average observational data into daily/monthly resolution 
#2.08 - write out _ALL/_NR tagged observational data to netcdf
#2.09 - write out _PERIODIC tagged observational data to netcdf
#2.10 - set defined areas for site averaging
#2.11 - set dicts for site area classification
#2.12 - get area cut for sites (using lat,lon and tag)
#2.13 - get area for sites (using lat,lon and tag)
#2.14 - get country (from lat,lon) using google api 
#2.15 - get anthrome classification
#2.16 - get observational continental tags (manually put together)


#2.01
#----------

#GET MODEL FILE INFO, NEED DIRECTORIES STRUCTURED IN FOLLOWING WAY:
#..../SPECIES/STARTYEAROFFILE_ENDYEAROFFILE/STARTYEARWANTED_ENDYEARWANTED/MODELNAME_VERTRES_VERSION_HORIZRES_MET_TIMERES_ADDITIONAL
#USE * IF NO DETAILS FOR CERTAIN ASPECT, i.e ACCMIP ONLY HAVE 1 VERSION,HORIZRES,MET 

def get_model_info(present_dir):
    
    paths = present_dir.split("/")
    species = paths[-4]
    file_years = paths[-3]
    file_start_year = file_years[:4]
    file_end_year = file_years[5:9]
    years = paths[-2]
    start_year = years[:4]
    end_year = years[5:9]
    model_path = paths[-1] 

    data_split = model_path.split('_') 
    model = data_split[0]
    vres = data_split[1]
    version = data_split[2]
    hres = data_split[3]
    met = data_split[4]
    timeres = data_split[5]
    additional = data_split[6]

    fname = '/work/home/db876/plotting_tools/model_files/%s_%s_%s_%s_%s_%s_%s_%s_%s.nc'%(model,vres,file_start_year,file_end_year,version,hres,met,timeres,additional)
    
    if species[:4] == 'GMAO':
        species = species[:4]+'_'+species[4:]
    if species == 'WINDDIRECTION':
        species = 'WIND_DIRECTION'
    if species == 'WINDSPEED':
        species = 'WIND_SPEED'
    
    return fname,species,start_year,end_year,vres,timeres

#2.02
#----------

#READ ONE MODEL GRIDBOX DATA FILE FROM NETCDF (USING OBS LAT/LON)

def read_model_one(model_file,species,obs_lat,obs_lon,start_year,end_year):

    if (species == 'NO2-PHOTOLYTIC') or (species == 'NO2-MOLYBDENUM'):
        species = 'NO2'
    elif (species == 'NOADJ'):
    	species = 'NO'

    root_grp = Dataset(model_file)
    std_var = root_grp.variables[species.lower()][:]
    raw_time = root_grp.variables['time'][:]
    lat_e = root_grp.variables['lat_edges'][:]
    lon_e = root_grp.variables['lon_edges'][:]
    lat_c = root_grp.variables['lat_centre'][:]
    lon_c = root_grp.variables['lon_centre'][:]
    grid_size = root_grp.variables['grid_size'][:]
    grid_size = grid_size[0]

    gridbox_count = len(lat_c)*len(lon_c)
    
    start_time = date2num(datetime.datetime(int(start_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    end_time = date2num(datetime.datetime(int(end_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    start_ind = np.searchsorted(raw_time,start_time)
    end_ind = np.searchsorted(raw_time,end_time)
    if (start_time < raw_time[0]) & (end_time > raw_time[-1]):
        raw_time = raw_time[:]
        std_var = std_var[:]
    elif start_time < raw_time[0]:
        raw_time = raw_time[:end_ind]
        std_var = std_var[:end_ind]
    elif end_time > raw_time[-1]:
        raw_time = raw_time[start_ind:]
        std_var = std_var[start_ind:]
    else:
        raw_time = raw_time[start_ind:end_ind]
        std_var = std_var[start_ind:end_ind]
    
    ref_time = raw_time - raw_time[0]
    datetime_time = num2date(raw_time,units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    
    #get raw and ref time in days
    raw_time = raw_time / 24.
    ref_time = ref_time / 24.
    
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    std_var = std_var[:,lat_n,lon_n]
    std_var_mask = np.ma.masked_where(std_var<=0,std_var)
    pd_var = pd.Series(std_var_mask, index=datetime_time)

    return raw_time,ref_time,datetime_time,std_var,pd_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count
     
#2.03
#----------

#READ ALL MODEL GRIDBOX DATA FROM NETCDF

def read_model_all(model_file,species,start_year,end_year):

    if (species == 'NO2-PHOTOLYTIC') or (species == 'NO2-MOLYBDENUM'):
        species = 'NO2'
    elif (species == 'NOADJ'):
    	species = 'NO'

    root_grp = Dataset(model_file)
    std_var = root_grp.variables[species.lower()][:]
    raw_time = root_grp.variables['time'][:]
    lat_e = root_grp.variables['lat_edges'][:]
    lon_e = root_grp.variables['lon_edges'][:]
    lat_c = root_grp.variables['lat_centre'][:]
    lon_c = root_grp.variables['lon_centre'][:]
    grid_size = root_grp.variables['grid_size'][:]
    grid_size = grid_size[0]

    gridbox_count = len(lat_c)*len(lon_c)
    
    start_time = date2num(datetime.datetime(int(start_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    end_time = date2num(datetime.datetime(int(end_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    start_ind = np.searchsorted(raw_time,start_time)
    end_ind = np.searchsorted(raw_time,end_time)
    if (start_time < raw_time[0]) & (end_time > raw_time[-1]):
        raw_time = raw_time[:]
        std_var = std_var[:]
    elif start_time < raw_time[0]:
        raw_time = raw_time[:end_ind]
        std_var = std_var[:end_ind]
    elif end_time > raw_time[-1]:
        raw_time = raw_time[start_ind:]
        std_var = std_var[start_ind:]
    else:
        raw_time = raw_time[start_ind:end_ind]
        std_var = std_var[start_ind:end_ind]
    
    ref_time = raw_time - raw_time[0]
    datetime_time = num2date(raw_time,units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    
    #get raw and ref time in days
    raw_time = raw_time / 24.
    ref_time = ref_time / 24.
    
    return raw_time,ref_time,datetime_time,std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count

#2.04
#----------

#GET MODEL LAT/LON INDEX FOR GRIDBOX THAT CONTAINS OBS LAT/LON POINT

def obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon): 
#convert obs_lon to same grid as model, if neccessary!
    if obs_lon > lon_e[-1]:
        diff = obs_lon - lon_e[-1]
        obs_lon = lon_e[0] + diff

#check which gridboxes each obs_lat & obs_lon convergance point lie in 
    lat_i = np.searchsorted(lat_e,obs_lat,side='left')
    lon_i = np.searchsorted(lon_e,obs_lon,side='left')
    lat_i = lat_i-1
    lon_i = lon_i-1
	
    return lat_i,lon_i

#2.05
#----------

#GET MODEL LAT/LON INDICES FOR GRIDBOXES THAT CONTAINS OBS LAT/LON POINTS AND GRIDBOX CENTRES

def obs_model_gridboxes(lat_e,lon_e,obs_lats,obs_lons): 
    #convert obs_lon to same grid as model, if neccessary!
    for i in range(len(obs_lons)):
        if obs_lons[i] > lon_e[-1]:
            diff = obs_lons[i] - lon_e[-1]
            obs_lons[i] = lon_e[0] + diff
            print obs_lons[i]
            
    #check which gridboxes each obs_lat & obs_lon convergance point lie in & then return central convergance point of lat and lon in that gridbox
    obs_lats_centre = []
    obs_lons_centre = []
    lat_indices = []
    lon_indices = []
    model_indices = []

    for a,b in zip(obs_lats, obs_lons):
        lat_i = np.searchsorted(lat_e,a,side='left')
        lon_i = np.searchsorted(lon_e,b,side='left')

        lower_lat = lat_e[lat_i-1]
        upper_lat = lat_e[lat_i]
        lower_lon = lon_e[lon_i-1]
        upper_lon = lon_e[lon_i]

        lat_diff = upper_lat - lower_lat
        lon_diff = upper_lon - lower_lon
        centre_lat = lower_lat + (lat_diff/2.)
        centre_lon = lower_lon + (lon_diff/2.)

        obs_lats_centre=np.append(obs_lats_centre,centre_lat)
        obs_lons_centre=np.append(obs_lons_centre,centre_lon)
        lat_indices = np.append(lat_indices,int(lat_i)-1)
        lon_indices = np.append(lon_indices,int(lon_i)-1)
        model_indices = np.append(model_indices,[int(lat_i)-1,int(lon_i)-1])
    model_indices = np.reshape(model_indices, (-1,2))
    lat_indices = lat_indices.astype(int)
    lon_indices = lon_indices.astype(int)
    model_indices = model_indices.astype(int)
    
	
    return obs_lats_centre, obs_lons_centre, model_indices,lat_indices,lon_indices

#2.06
#----------
  
#GET MODEL GRID FOR MODEL TYPE

def model_grids(type):
    if type == 'GEOSCHEM_4x5':
        #lat lon dims for 4x5 grid  
        lat_c = np.arange(-86.,87.,4)       
        lat_c = np.insert(lat_c,0,-89)
        lat_c = np.append(lat_c,89)

        lon_c = np.arange(-180,176,5)

        lat_e = np.arange(-88.,90,4)
        lat_e = np.insert(lat_e,0,-90.)
        lat_e = np.append(lat_e,90.)

        lon_e = np.arange(-182.5,178,5)

    elif type == 'GEOSCHEM_2x2.5':
        #lat lon dims for 2x2.5 grid

        lat_c = np.arange(-88.,89.,2)
        lat_c = np.insert(lat_c,0,-89.5)
        lat_c = np.append(lat_c,89.5)

        lon_c = np.arange(-180,178,2.5)

        lat_e = np.arange(-89.,90,2)
        lat_e = np.insert(lat_e,0,-90.)
        lat_e = np.append(lat_e,90.)

        lon_e = np.arange(-181.25,179,2.5)
    
    elif type == 'GFDLAM3':
        #lat lon dims for GFDL 2x2.5 grid]
        
        lat_e = np.arange(-90,91,2)                                                                                                                                                                                                               
        lon_e = np.arange(-180,181,2.5)
           
        lat_c = np.arange(-89,90,2)                                                                                                                                                                                                               
        lon_c = np.arange(-178.75,179,2.5)
        
    return lat_c,lat_e,lon_c,lon_e

#2.07
#----------
   
#REGRID MODEL INTO NEW GRID USING AVERAGING BETWEEN BOXES

def regrid_interp(var,orig_lat_e,orig_lat_c,orig_lon_e,orig_lon_c,min_lat,min_lon,data_type,max_ph):

    new_lat_e = np.arange(-90,91,min_lat)
    new_lat_c = np.arange((-90+(min_lat/2.)),90,min_lat)
    new_lon_e = np.arange(-180,181,min_lon)
    new_lon_c = np.arange((-180+(min_lon/2.)),180,min_lon)
    n_boxes = len(new_lat_c)*len(new_lon_c)

    z = np.empty((len(new_lat_c),len(new_lon_c)))

    #find edge of boxes for new lat and new lon edges to left of gridbox new lat and lon edges are in. 
    for lat_i in range(len(new_lat_c)):
        for lon_i in range(len(new_lon_c)): 
        
            exception_lat = False
            exception_lon = False
            
            #print 'new lon edges', new_lon_e[lon_i],new_lon_e[lon_i+1]  
            #print 'new lat edges', new_lat_e[lat_i],new_lat_e[lat_i+1] 
        
            if new_lat_e[lat_i] > orig_lat_e[-1]:
                diff = new_lat_e[lat_i] - orig_lat_e[-1]
                corrected_new_lat_e = orig_lat_e[0] + diff
                lat_box_left_i = np.where(orig_lat_e == np.max(orig_lat_e[orig_lat_e <= corrected_new_lat_e]))[0][0]
            else:
                lat_box_left_i = np.where(orig_lat_e == np.max(orig_lat_e[orig_lat_e <= new_lat_e[lat_i]]))[0][0]
            
            if new_lat_e[lat_i+1] > orig_lat_e[-1]:
                diff = new_lat_e[lat_i+1] - orig_lat_e[-1]
                corrected_new_lat_e = orig_lat_e[0] + diff
                lat_box_right_i = np.where(orig_lat_e == np.min(orig_lat_e[orig_lat_e >= corrected_new_lat_e]))[0][0]
            else:
                lat_box_right_i = np.where(orig_lat_e == np.min(orig_lat_e[orig_lat_e >= new_lat_e[lat_i+1]]))[0][0]
            if lat_box_right_i - lat_box_left_i < 0: 
                #print 'lat exception'
                exception_lat = True
            
            if new_lon_e[lon_i] > orig_lon_e[-1]:
                diff = new_lon_e[lon_i] - orig_lon_e[-1]
                corrected_new_lon_e = orig_lon_e[0] + diff
                lon_box_left_i = np.where(orig_lon_e == np.max(orig_lon_e[orig_lon_e <= corrected_new_lon_e]))[0][0]   
            else:
                lon_box_left_i = np.where(orig_lon_e == np.max(orig_lon_e[orig_lon_e <= new_lon_e[lon_i]]))[0][0]
            
            if new_lon_e[lon_i+1] > orig_lon_e[-1]:   
                diff = new_lon_e[lon_i+1] - orig_lon_e[-1]
                corrected_new_lon_e = orig_lon_e[0] + diff
                lon_box_right_i = np.where(orig_lon_e == np.min(orig_lon_e[orig_lon_e >= corrected_new_lon_e]))[0][0]
            else:
                lon_box_right_i = np.where(orig_lon_e == np.min(orig_lon_e[orig_lon_e >= new_lon_e[lon_i+1]]))[0][0]
            if lon_box_right_i - lon_box_left_i < 0: 
                #print 'lon exception'
                exception_lon = True
            
            #print 'orig lon edges',orig_lon_e[lon_box_left_i],orig_lon_e[lon_box_right_i]
            #print 'orig lat edges',orig_lat_e[lat_box_left_i],orig_lat_e[lat_box_right_i]
            
            orig_lat_left = orig_lat_e[lat_box_left_i]
            orig_lon_left = orig_lon_e[lon_box_left_i]
            
            new_lat_left = new_lat_e[lat_i]
            new_lat_right = new_lat_e[lat_i+1]
            new_lon_left = new_lon_e[lon_i]
            new_lon_right = new_lon_e[lon_i+1]
            
            if exception_lat == True:
                lat_diff_i = len(orig_lat_e) - lat_box_left_i
            else:
                orig_lat_right = orig_lat_e[lat_box_right_i]
                lat_diff_i = lat_box_right_i - lat_box_left_i
            
            if exception_lon == True:    
                lon_diff_i = len(orig_lon_e) - lon_box_left_i
            else: 
                orig_lon_right = orig_lon_e[lon_box_right_i]
                lon_diff_i = lon_box_right_i - lon_box_left_i
                
            #print lat_box_left_i, lat_box_right_i,lon_box_left_i, lon_box_right_i 
            #print lat_diff_i,lon_diff_i
            
            #how many boxes does new boxes overlap?
            #1 box
            if (lat_diff_i == 1) & (lon_diff_i == 1):
                #print '1 box'
            
                if exception_lat == True:
                    orig_box1_lat_diff = 180. -  np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_right_i]])[0])
                else:
                    orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_right_i]])[0])
                
                lat_box1_n = lat_box_left_i

                if exception_lon == True:
                    orig_box1_lon_diff = 360. - np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_right_i]])[0])
                else:
                    orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_right_i]])[0])
                
                if lon_box_left_i == len(orig_lon_c):
                    lon_box1_n = 0
                else:
                    lon_box1_n = lon_box_left_i
        
                z[lat_i,lon_i] = var[lat_box1_n,lon_box1_n]
                
                #print var[lat_box1_n,lon_box1_n]
                #print z[lat_i,lon_i]
    
            #2 boxes
            elif (lat_diff_i) + (lon_diff_i) == 3:
                #print '2 boxes'
                #2 lat boxes, 1 lon box
                if (lat_diff_i) == 2:
                    if exception_lat == True:
                        orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_left_i+1]])[0])
                        orig_box2_lat_diff = 180. - np.abs(np.diff([orig_lat_e[lat_box_left_i+1] , orig_lat_e[lat_box_right_i]])[0])
                    else:
                        orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_left_i+1]])[0])
                        orig_box2_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i+1] , orig_lat_e[lat_box_right_i]])[0])
                    
                    lat_box1_n = lat_box_left_i
                    lat_box2_n = lat_box_right_i-1
                
                    if exception_lon == True:
                        orig_box1_lon_diff = 360. - np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_right_i]])[0])
                    else:
                        orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_right_i]])[0])
                        

                    if lon_box_left_i == len(orig_lon_c):
                        lon_box1_n = 0
                    else:
                        lon_box1_n = lon_box_left_i
                
                    actual_box1_lat_diff = np.abs(np.diff([new_lat_e[lat_i],orig_lat_e[lat_box_left_i+1]])[0])
                    actual_box2_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i+1],new_lat_e[lat_i+1]])[0])
                    actual_box1_lon_diff = np.abs(np.diff([new_lon_e[lon_i],new_lon_e[lon_i+1]])[0])
                
                    box1_area = orig_box1_lat_diff * orig_box1_lon_diff
                    actual_box1_area = actual_box1_lat_diff * actual_box1_lon_diff
                    box1_fract = (1./box1_area)*actual_box1_area
                    #print 'actual box 1 fract = ',box1_fract
            
                    box2_area = orig_box2_lat_diff * orig_box1_lon_diff
                    actual_box2_area = actual_box2_lat_diff * actual_box1_lon_diff
                    box2_fract = (1./box2_area)*actual_box2_area
                    #print 'actual box 2 fract = ',box2_fract
                    
                    if data_type == 'ph':
                        data = [var[lat_box1_n,lon_box1_n],var[lat_box2_n,lon_box1_n]]
                        weights=[box1_fract,box2_fract]
                        z[lat_i,lon_i] = cylic_weighted_average(data,weights,max_ph)
                    else:
                        z[lat_i,lon_i] = np.average([var[lat_box1_n,lon_box1_n],
                                                     var[lat_box2_n,lon_box1_n]],weights=[box1_fract,box2_fract])
                       
                    #print var[lat_box1_n,lon_box1_n],var[lat_box2_n,lon_box1_n]                          
                    #print z[lat_i,lon_i]
        
                #2 lon boxes, 1 lat box
                elif (lon_diff_i) == 2:
                    if exception_lat == True:
                        orig_box1_lat_diff = 180. - np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_right_i]])[0])
                    else:
                        orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_right_i]])[0])
                    
                    lat_box1_n = lat_box_left_i
                
                    if exception_lon == True:
                        orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_left_i+1]])[0])
                        orig_box2_lon_diff = 360. - np.abs(np.diff([orig_lon_e[lon_box_left_i+1] , orig_lon_e[lon_box_right_i]])[0])
                    else:
                        orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_left_i+1]])[0])
                        orig_box2_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i+1] , orig_lon_e[lon_box_right_i]])[0])
                    
                    if lon_box_left_i == len(orig_lon_c):
                        lon_box1_n = 0
                        lon_box2_n = 1
                    elif lon_box_left_i+1 == len(orig_lon_c):
                        lon_box1_n = lon_box_left_i
                        lon_box2_n = 0
                    else:
                        lon_box1_n = lon_box_left_i
                        lon_box2_n = lon_box_left_i+1 
                    
                    actual_box1_lat_diff = np.abs(np.diff([new_lat_e[lat_i],new_lat_e[lat_i+1]])[0])
                    actual_box1_lon_diff = np.abs(np.diff([new_lon_e[lon_i],orig_lon_e[lon_box_left_i+1]])[0])
                    actual_box2_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i+1],new_lon_e[lon_i+1]])[0])
                
                    box1_area = orig_box1_lat_diff * orig_box1_lon_diff
                    actual_box1_area = actual_box1_lat_diff * actual_box1_lon_diff
                    box1_fract = (1./box1_area)*actual_box1_area
                    #print 'actual box 1 fract = ',box1_fract
            
                    box2_area = orig_box1_lat_diff * orig_box2_lon_diff
                    actual_box2_area = actual_box1_lat_diff * actual_box2_lon_diff
                    box2_fract = (1./box2_area)*actual_box2_area 
                    #print 'actual box 2 fract = ',box2_fract
            
                    if data_type == 'ph':
                        data = [var[lat_box1_n,lon_box1_n],var[lat_box1_n,lon_box2_n]]
                        weights=[box1_fract,box2_fract]
                        z[lat_i,lon_i] = cylic_weighted_average(data,weights,max_ph)
                    else:
                        z[lat_i,lon_i] = np.average([var[lat_box1_n,lon_box1_n],
                                                     var[lat_box1_n,lon_box2_n]],weights=[box1_fract,box2_fract])
                    
                    #print var[lat_box1_n,lon_box1_n],var[lat_box1_n,lon_box2_n]
                    #print z[lat_i,lon_i]
                    
    
            #4 boxes
            elif (lat_diff_i) + (lon_diff_i) == 4:    
            
                #print '4 boxes'
                if exception_lat == True:
                    orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_left_i+1]])[0])
                    orig_box2_lat_diff = 180. - np.abs(np.diff([orig_lat_e[lat_box_left_i+1] , orig_lat_e[lat_box_right_i]])[0])
                else:
                    orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_left_i+1]])[0])
                    orig_box2_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i+1] , orig_lat_e[lat_box_right_i]])[0])
                
                lat_box1_n = lat_box_left_i
                lat_box2_n = lat_box_right_i-1
            
                if exception_lon == True:
                    orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_left_i+1]])[0])
                    orig_box2_lon_diff = 360. - np.abs(np.diff([orig_lon_e[lon_box_left_i+1] , orig_lon_e[lon_box_right_i]])[0])
                else:
                    orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_left_i+1]])[0])
                    orig_box2_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i+1] , orig_lon_e[lon_box_right_i]])[0])
                #print orig_box2_lon_diff
                
                if lon_box_left_i == len(orig_lon_c):
                    lon_box1_n = 0
                    lon_box2_n = 1
                elif lon_box_left_i+1 == len(orig_lon_c):
                    lon_box1_n = lon_box_left_i
                    lon_box2_n = 0
                else:
                    lon_box1_n = lon_box_left_i
                    lon_box2_n = lon_box_left_i+1
                
                actual_box1_lat_diff = np.abs(np.diff([new_lat_e[lat_i],orig_lat_e[lat_box_left_i+1]])[0])
                actual_box2_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i+1],new_lat_e[lat_i+1]])[0])
                actual_box1_lon_diff = np.abs(np.diff([new_lon_e[lon_i],orig_lon_e[lon_box_left_i+1]])[0])
                actual_box2_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i+1],new_lon_e[lon_i+1]])[0])
                
                #lower left box of grid
                box1_area = orig_box1_lat_diff * orig_box1_lon_diff
                actual_box1_area = actual_box1_lat_diff * actual_box1_lon_diff
                box1_fract = (1./box1_area)*actual_box1_area
                #print 'actual box 1 fract = ',box1_fract
    
                #lower right box of grid
                box2_area = orig_box1_lat_diff * orig_box2_lon_diff
                actual_box2_area = actual_box1_lat_diff * actual_box2_lon_diff
                box2_fract = (1./box2_area)*actual_box2_area
                #print 'actual box 2 fract = ',box2_fract
                
                #upper left box of grid
                box3_area = orig_box2_lat_diff * orig_box1_lon_diff
                actual_box3_area = actual_box2_lat_diff * actual_box1_lon_diff
                box3_fract = (1./box3_area)*actual_box3_area
                #print 'actual box 3 fract = ',box3_fract
        
                #upper right box of grid
                box4_area = orig_box2_lat_diff * orig_box2_lon_diff
                actual_box4_area = actual_box2_lat_diff * actual_box2_lon_diff
                box4_fract = (1./box4_area)*actual_box4_area
                #print 'actual box 4 fract = ',box4_fract
                
                if data_type == 'ph':
                    data = [var[lat_box1_n,lon_box1_n],var[lat_box1_n,lon_box2_n],var[lat_box2_n,lon_box1_n],var[lat_box2_n,lon_box2_n]]
                    weights=[box1_fract,box2_fract,box3_fract,box4_fract]
                    z[lat_i,lon_i] = cylic_weighted_average(data,weights,max_ph)
                else:
                    z[lat_i,lon_i] = np.average([var[lat_box1_n,lon_box1_n],
                                                 var[lat_box1_n,lon_box2_n],
                                                 var[lat_box2_n,lon_box1_n],
                                                 var[lat_box2_n,lon_box2_n]],weights=[box1_fract,box2_fract,box3_fract,box4_fract])
                                  
                #print var[lat_box1_n,lon_box1_n],var[lat_box1_n,lon_box2_n],var[lat_box2_n,lon_box1_n],var[lat_box2_n,lon_box2_n]
                #print z[lat_i,lon_i]
        
            #print ''                                      
    
    
    return new_lat_e,new_lat_c,new_lon_e,new_lon_c,z
 
#2.08
#----------
   
#GET PRESSURE EGES FOR GEOS5 AND MERRA MET DATA

def get_pressure_edges(psfc):

    # !  GEOS-4, GEOS-5, GEOS-5.7, and MERRA (hybrid grids):
    # !  ----------------------------------------------------------------------------
    # !  For GEOS-4/GEOS-5/MERRA met data products, the pressure at the bottom edge 
    # !  of grid box (I,J,L) is defined as follows:
    # !                                                                             .
    # !     Pedge(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
    # !                                                                             .
    # !  where
    # !                                                                             .
    # !     Psurface(I,J) is  the "true" surface pressure at lon,lat (I,J)
    # !     Ap(L)         has the same units as surface pressure [hPa]
    # !     Bp(L)         is  a unitless constant given at level edges
    # !                                                                             .
    # !  Ap(L) and Bp(L) are given to us by GMAO.
    # 
    # 
    # !  The following are true for GCAP, GEOS-3, GEOS-4, GEOS-5, MERRA:
    # !  ----------------------------------------------------------------------------
    # !  (1) Bp(LLPAR+1) = 0.0          (L=LLPAR+1 is the atmosphere top)
    # !  (2) Bp(1)       = 1.0          (L=1       is the surface       )
    # !  (3) PTOP        = Ap(LLPAR+1)  (L=LLPAR+1 is the atmosphere top) 
    # 
    # 
    # !-----------------------------------------------------------------
    # ! GEOS-5/MERRA 47-level reduced vertical grid
    # !  
    # !  Bottom   Bottom    # levels
    # !  edge of  edge prs  lumped 
    # !  level    (hPa)     together
    # !
    # !   PTOP       0.010   
    # !    47        0.066     4
    # !    46        0.211     4
    # !    45        0.617     4
    # !    44        1.651     4
    # !    43        4.077     4
    # !    42        9.293     4
    # !    41       19.792     4
    # !    40       28.368     2
    # !    39       40.175     2
    # !    38       56.388     2
    # !    37       78.512     2
    # ! %%%% START LUMPING LEVELS ABOVE HERE %%%%%
    # !    36       92.366       
    # !    35      108.663
    # !    34      127.837
    # !    33      150.393
    # !    32      176.930
    # ! %%%% FIXED-PRESSURE LEVELS BEGIN HERE %%%%
    # !-----------------------------------------------------------------

    # Ap [hPa] for 47 levels (48 edges)
    AP = np.array([0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,
            1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,
            4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
            7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,
            1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,
            1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
            2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,
            2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,
            1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
            7.851231e+01, 5.638791e+01, 4.017541e+01, 2.836781e+01, 
            1.979160e+01, 9.292942e+00, 4.076571e+00, 1.650790e+00, 
            6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02])

    # Bp [unitless] for 47 levels (48 edges)
    BP = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,
            9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,
            8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
            7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,
            6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,
            4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
            2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,
            6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])

    pressure_edges = np.empty((psfc.shape[0],48,psfc.shape[1],psfc.shape[2]))
    lat_inds = range(psfc.shape[1])
    lon_inds = range(psfc.shape[2])
    for i in range(psfc.shape[0]):
        print i
        for lat_i in lat_inds:
            for lon_i in lon_inds:
                pressure_edges[i,:,lat_i,lon_i] = AP + (BP * psfc[i,lat_i,lon_i])
        print psfc[i,lat_i,lon_i]
        print pressure_edges[i,:,lat_i,lon_i]
    
    return pressure_edges
    
    
#-------------------------------------------------------------
#-------------------------------------------------------------
#SECTION 3 - SPECTRAL ANALYSIS
    
#TAKE FULL STANDARD FFT (using scipy package)
def fft_spectra(time,var,ofac,SAMP_R,w=False,kp=[]):
    var_mean = np.mean(var)
    var = var - var_mean

    if w == True:
        window = signal.hanning(len(var))
        var = var*window
        amp_corr = 1./(sum(window)/len(window))
    else:
        amp_corr = 1.

    #add zero padding, oversampling of n
    orig_time = np.copy(time)
    n_zeros = (len(time)*ofac) - (len(time))
    zeros_array = [0]*int(n_zeros)
    time = np.append(time,zeros_array)
    var = np.append(var,zeros_array)

    #set frequencies (with zero padded oversampling
    n_points = len(time)
    N = np.copy(n_points)
    T = SAMP_R*N
    df = 1./T
    fft_freqs = np.array([df*n if n<N/2 else df*(n-N) for n in range(N)])
    fft_array = fft(var,n_points)
    fft_mag = np.abs(fft_array)
    ratio = n_points/len(var)
    valid = fft_freqs > 0
    fft_periods = 1./fft_freqs
    fft_mag = fft_mag/len(orig_time)
    fft_fr = fft_array.real
    fft_fi = fft_array.imag 
    fft_periods,fft_mag,fft_fr,fft_fi = fft_periods[valid],fft_mag[valid],fft_fr[valid],fft_fi[valid]
    fft_mag = fft_mag*2
    fft_ph = np.arctan2(fft_fi,fft_fr)

    #if have key periods calculate amp/phase at specific frequencies and impose into full spectrum where closest point is
    if len(kp) > 0:
        s_periods,s_mag,s_ph,s_fr,s_fi = take_lomb_spec(time,var,False,kp)
        for count in range(len(kp)):
            closest_period = min(range(len(fft_periods)), key=lambda i: abs(fft_periods[i]-kp[count]))
        
            fft_periods[closest_period] = s_periods[count]
            fft_mag[closest_period] = s_mag[count]
            fft_ph[closest_period] = s_ph[count]
            fft_fr[closest_period] = s_fr[count]
            fft_fi[closest_period] = s_fi[count]

    if w == True:
        fft_mag = fft_mag * amp_corr   


    return fft_periods,fft_mag,fft_ph,fft_fr,fft_fi,fft_array,amp_corr


#TAKE FULL STANDARD LOMB-SCARGLE - NEED FORTRAN FILE COMPILED WITH F2PY (lomb_phase.so)
def lomb_spectra(time,var,OFAC,SAMP_R,w=False,kp=[]):
#Note - Sampling rate can be lower than average nyquist frequency, as when data is gapped nyquist breaks down. Thus if have something that is non-uniformly spaced, just choose a low samp_r.

    if w == True:
        window = signal.hanning(len(var))
        #window = signal.flattop(len(var))
        var_mean = np.mean(var)
        var = var - var_mean
        var = var*window
        amp_corr = 1./(sum(window)/len(window))
    else:
        amp_corr = 1.

    #set frequencies (equivalent to FFT freqs)
    N = (len(time)*OFAC)
    T = SAMP_R*N
    df = 1./T
    freqs = np.array([df*n if n<N/2 else df*(n-N) for n in range(N)])
    test = freqs > 0
    freqs = freqs[test]

    #take lomb
    fb, mag, ph, fr, fi = lomb_phase.lomb(time,var,freqs)

    periods = 1./freqs

    #CORRECT IMAGINARY COMPONENTS TO FFT EQUIVALENTS
    for i in range(len(fi)):
        if fi[i] < 0:
            fi[i] = fi[i]*-1
        elif fi[i] > 0:
            fi[i] = -fi[i]
    

    #if have key periods calculate amp/phase at specific frequencies and impose into full spectrum where closest point is
    if len(kp) > 0:
        s_periods,s_mag,s_ph,s_fr,s_fi = lomb_specific(time,var,False,kp)
        for count in range(len(kp)):
            closest_period = min(range(len(periods)), key=lambda i: abs(periods[i]-kp[count]))
        
            periods[closest_period] = s_periods[count]
            mag[closest_period] = s_mag[count]
            ph[closest_period] = s_ph[count]
            fr[closest_period] = s_fr[count]
            fi[closest_period] = s_fi[count]

    if w == True:
        mag = mag * amp_corr
    
    return periods,mag,ph,fr,fi,amp_corr
    
    
#3.02
#----------

#TAKE SPECIFIC LOMB-SCARGLE - NEED FORTRAN FILE COMPILED WITH F2PY (lomb_phase_spec.so)
def lomb_specific(time,var,w = False,key_periods=[1./4.,1./3.,1./2.,1.,365.25/4.,365.25/3.,365.25/2.,365.25]):
    if w == True:
        window = signal.hanning(len(var))
        var_mean = np.mean(var)
        var = var - var_mean
        var = var*window
        amp_corr = 1./(sum(window)/len(window))
	
    freqs = 1./np.array(key_periods)
    
	#take lomb
    fb, mag, ph, fr, fi = lomb_phase_spec.lomb(time,var,freqs)
    
    if w == True:
        mag = mag * amp_corr
        
    periods = 1./freqs
    
    #CORRECT IMAGINARY COMPONENTS TO FFT EQUIVALENTS
    for i in range(len(fi)):
        if fi[i] < 0:
            fi[i] = fi[i]*-1
        elif fi[i] > 0:
            fi[i] = -fi[i]
    	
    return periods,mag,ph,fr,fi

#3.03
#----------

#TAKE IFFT OF FULL STANDARD LOMB-SCARGLE
def lomb_spectra_ifft(key_periods,time,vals,periods,fr,fi,ofac,amp_corr,w=False):
    #note using windowing reduces length of output time series as need to cut data at start and end of time series. This loss is severe for LSP.
    
    #ONLY WORKS ACCURATELY FOR NON-GAPPED DATA 

    closest_periods = []
    
    for period in key_periods:
        #get indices of key periods, setting real and imag parts to zero that are not relevant.
        #Need to do this as LSP ifft of noise is set equivalent to FFT ifft of noise
        #reconstruction only works for significant periodicities
        closest_period = min(range(len(periods)), key=lambda i: abs(periods[i]-period))
        closest_periods.append(closest_period)

    F = [0]*((len(fr)*2)+2)
    
    #set first real value to average, imag to 0 
    F[0] = complex((np.average(vals)*len(vals)*ofac),0)

    #Get reverse real and imaginary values
    rev_fr=np.copy(fr[::-1])
    rev_fi=np.copy(fi[::-1])

    f_index = 1

    #Fill Fourier Spectrum real and imaginary values
    for i in range(len(fr)):
        F[f_index] = complex(fr[i],fi[i])
        f_index+=1

    F[f_index] = complex(0,0)
    f_index+=1
    
    for i in range(len(fr)):
        F[f_index] = complex(rev_fr[i],-rev_fi[i])
        f_index+=1

    F = np.array(F) 
    
    #Take ifft and just take real values
    ifft_ts = np.fft.ifft(F)
    ifft_ts = ifft_ts.astype('float64')

    #cut reconstructed time series if oversampled
    if ofac > 1:
        actual_len = len(time)
        ifft_ts = ifft_ts[:actual_len]
        if w == True:
            window = signal.hanning(len(vals))
            ifft_ts = ifft_ts/window
            #cut off 40 percent of ends if windowed. Yes loss is severe for LSP
            len_cut = int((len(ifft_ts)/100.)*40.)
            ifft_ts = ifft_ts[len_cut:-len_cut]
            time = time[len_cut:-len_cut]
    
    return time,ifft_ts

#3.04
#----------
    
#TAKE IFFT OF SPECIFIC LOMB-SCARGLE 
def lomb_specific_ifft(time,vals,periods,mag,ph,samp_r):
    #Reconstruction of time series based on known significant components

    pi2 = np.pi*2.
    
    #get max time length
    max_time = time[-1] - time[0]
    
    #get time array based on max time length and sampling frequency.
    t = []
    t = np.arange(time[0],time[-1]+(samp_r/2.),samp_r)
    
    waveform = np.array([0]*len(t))
        
    #put together waveform with significant component amplitudes and phases.
    for i in range(len(mag)):
        waveform = waveform + (mag[i]*(np.cos((pi2*t/periods[i])-(ph[i]))))
    
    waveform = waveform+np.average(vals) 
        
    return t,waveform

#-----------

#TAKE LOMB AT SPECIFIC A PRIORI KNOWN FREQUENCIES AND SUBTRACT SUPERPOSED
#WAVEFORM FROM RAW TIME SERIES
def lomb_remove_periodicity(time,var,w = False,key_periods=[1./4.,1./3.,1./2.,1.,365.25/4.,365.25/3.,365.25/2.,365.25]):
    valid = var > 0
    cut_time = time[valid]
    cut_var = var[valid]
    mean = np.average(cut_var)
    
    periods,mag,ph,fr,fi = take_lomb_spec(cut_time,cut_var,w,key_periods)
    
    print mag
    
    full_mag,full_min_ph,full_max_ph,full_waveform,full_ff = period_convolution(key_periods,time,mag,ph,mean)
    
    plt.plot(time,full_waveform)
    
    var[~valid] = np.NaN
    var = (var - full_waveform)+mean
    return time,var

#FFT
#-------------------------------------------------------------



#3.06
#----------

#TAKE IFFT (using scipy package)

def fft_ifft(time,raw_ts,fft_array,ofac,w=False):  
   
    #Take ifft and just take real values
    ifft_ts = np.fft.ifft(fft_array)
    ifft_ts = ifft_ts.astype('float64')

    #cut reconstructed time series if oversampled
    if ofac > 1:
        actual_len = len(time)
        ifft_ts = ifft_ts[:actual_len]
        
        if w == True:
            window = signal.hanning(len(raw_ts))
            ifft_ts = ifft_ts/window
            #cut off 2 percent of ends if windowed.
            len_cut = int((len(ifft_ts)/100.)*2.)
            ifft_ts = ifft_ts[len_cut:-len_cut]
            time = time[len_cut:-len_cut]
            ifft_ts = ifft_ts + np.average(raw_ts)
    
    return time,ifft_ts
    
#SPECTRAL ANALYSIS PROCESSING TOOLS
#-------------------------------------------------------------
#3.
#----------

#OUTPUT TIME SERIES STATS
def raw_stats(ts):
    mean = np.average(ts)
    p1 = np.percentile(ts, 1)
    p5 = np.percentile(ts, 5)
    p25 = np.percentile(ts, 25)
    p50 = np.percentile(ts, 50)
    p75 = np.percentile(ts, 75)
    p95 = np.percentile(ts, 95)
    p99 = np.percentile(ts, 99)
    
    return mean,p1,p5,p25,p50,p75,p95,p99


def roundTime(dt=None, roundTo=60):
   """Round a datetime object to any time laps in seconds
   dt : datetime.datetime object, default now.
   roundTo : Closest number of seconds to round to, default 1 minute.
   Author: Thierry Husson 2012 - Use it as you want but don't blame me.
   """
   if dt == None : dt = datetime.datetime.now()
   seconds = (dt - dt.min).seconds
   # // is a floor division, not a comment on following line:
   rounding = (seconds+roundTo/2) // roundTo * roundTo
   return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)

#cut full time series into 4 separate seasonal time series
        
def cut_season(orig_datetime,orig_time,orig_ts):
    #cut time series into seasonal chunks
    #Dec,Jan,Feb = Winter, Mar,Apr,May = Spring, Jun,Jul,Aug = Summer, Sep,Oct,Nov = Autumn

    #do correction so there are no yearly gaps in time series (for spectral anlalsysi purposes)
    
    leap_years = range(1852,2200,4)
    #If the year can be evenly divided by 100, it is NOT a leap year, unless the year is also evenly divisible by 400. Then it is a leap year
    #Thus remove follwing years
    leap_years.remove(1900)
    leap_years.remove(2100)
    
    years = []
    months = []
    
    for i in range(len(orig_datetime)):
        months.append(orig_datetime[i].month)
        years.append(orig_datetime[i].year)
    
    spring_ts = []
    summer_ts = []
    autumn_ts = []
    winter_ts = []
    
    spring_inds = []
    summer_inds = []
    autumn_inds = []
    winter_inds = []
    
    spring_time = []
    summer_time = []
    autumn_time = []
    winter_time = []
    
    winter_years = []
    
    correction = 0
    for i in range(len(months)):
        if (months[i] == 3) or (months[i] == 4) or (months[i] == 5):
            spring_time.append(orig_time[i]-correction)
            spring_ts.append(orig_ts[i])
            spring_inds.append(i)
        elif (months[i] == 6) or (months[i] == 7) or (months[i] == 8):
            summer_time.append(orig_time[i]-correction)
            summer_ts.append(orig_ts[i])
            summer_inds.append(i)
        elif (months[i] == 9) or (months[i] == 10) or (months[i] == 11):
            autumn_time.append(orig_time[i]-correction)
            autumn_ts.append(orig_ts[i])
            autumn_inds.append(i)
        elif (months[i] == 12) or (months[i] == 1) or (months[i] == 2):
            try:
                if (months[i-1] < 12) & (months[i-1] > 2):
                    if winter_years[-1] in leap_years:
                        correction+=366
                        winter_time.append(orig_time[i]-correction)
                        winter_ts.append(orig_ts[i])
                        winter_years.append(years[i])
                        winter_inds.append(i)
                    else:
                        correction+=365
                        winter_time.append(orig_time[i]-correction)
                        winter_ts.append(orig_ts[i])
                        winter_years.append(years[i])
                        winter_inds.append(i)
                else:
                    winter_time.append(orig_time[i]-correction)
                    winter_ts.append(orig_ts[i])
                    winter_years.append(years[i])
                    winter_inds.append(i)
            except:
                winter_time.append(orig_time[i]-correction)
                winter_ts.append(orig_ts[i])
                winter_years.append(years[i])
                winter_inds.append(i)
    
    spring_inds = np.array(spring_inds)
    summer_inds = np.array(summer_inds)
    autumn_inds = np.array(autumn_inds)
    winter_inds = np.array(winter_inds)

    return spring_time,spring_ts,spring_inds,summer_time,summer_ts,summer_inds,autumn_time,autumn_ts,autumn_inds,winter_time,winter_ts,winter_inds
    
def cut_daynight(orig_datetime,orig_time,orig_ts,site_lat,site_lon,timeres):

    #cut time series into day/night chunks
    #round localtime at each site to closest hour to sunrise/sunset
    #i.e if sunset is at minute 29, sunset is at current hour.
    #if sunset is at minute 31 sunset is at next hour
    
    if (timeres == 'D') or (timeres == 'M'):
        return orig_time,orig_ts,orig_time,orig_ts
    else:
        #set ephem location instance
        my_location = ep.Observer()
        my_location.lat = str(site_lat)
        my_location.lon = str(site_lon)
        my_location.pressure = 0
        my_location.horizon = '-0:34'  # per USNO definition
        
        #create grid time array for 1 reference year with leap day (year 2000)
        start = datetime.datetime(year = 2000, month = 1, day = 1, hour = 0, minute = 0)
        end = datetime.datetime(year = 2001, month = 1, day = 1, hour = 0, minute = 0)
        ref_date_dt = pd.date_range(start,end,freq='H')[:-1]
        ref_date_str = np.array([d.strftime('%m%d%H') for d in ref_date_dt]).astype('int')
        
        ref_state_array = []
        
        #for ref year determine for each time if day/night
        for i in range(len(ref_date_dt)):
        
            #update date
            current_dt = ref_date_dt[i]
            my_location.date = current_dt
            
            #get prev sunrise and sunset for current time
            try:
                next_sunrise = my_location.next_rising(ep.Sun()).datetime()
                next_sunset = my_location.next_setting(ep.Sun()).datetime()
            
                #determine if day or night
                if next_sunset < next_sunrise:
                    #if delta between current datetime and next sunrise/sunset is < 30 mins
                    #then set state is not representative of hour and need to correct
                    #(i.e time could be 7am and sunrise is at 7.01am, state would be set as night, but really there are 59 minutes of day)
                    if (next_sunset - current_dt) <  datetime.timedelta(minutes=30):
                        current_state = 'night'
                    else:
                        current_state = 'day'
                        
                else:
                    if (next_sunrise - current_dt) <  datetime.timedelta(minutes=30):
                        current_state = 'day'
                    else: 
                        current_state = 'night'
                    
                #append to ref state array
                ref_state_array.append(current_state)
        
            #exceptions for days when sun never rises or sets
            except ep.AlwaysUpError:
                current_state = 'day'
                ref_state_array.append(current_state)
            
            except ep.NeverUpError:
                current_state = 'night'
                ref_state_array.append(current_state)
    
            #print current_dt,current_state,next_sunrise,next_sunset
    
        ref_state_array = np.array(ref_state_array)
        
        #create str array from actual datetimes of month+day+hour
        orig_date_str = np.array([d.strftime('%m%d%H') for d in orig_datetime]).astype('int')
                
        #cut original timeseries now into day and night chunks
        #based on ref year determined states at each hour
        match_inds = np.searchsorted(ref_date_str,orig_date_str)
        ref_day_inds = np.where(ref_state_array == 'day')[0]
        ref_night_inds = np.where(ref_state_array == 'night')[0]
        
        day_inds = np.array([i for i, item in enumerate(match_inds) if item in ref_day_inds])
        night_inds = np.array([i for i, item in enumerate(match_inds) if item in ref_night_inds])
        
        day_time = orig_time[day_inds]
        day_ts = orig_ts[day_inds]
        night_time = orig_time[night_inds]
        night_ts = orig_ts[night_inds]
        
        return day_time,day_ts,night_time,night_ts
        
def get_daynight_inds(orig_datetime,orig_ts,site_lat,site_lon):

    #cut time series into day/night chunks
    #round localtime at each site to closest hour to sunrise/sunset
    #i.e if sunset is at minute 29, sunset is at current hour.
    #if sunset is at minute 31 sunset is at next hour
    
    #set ephem location instance
    my_location = ep.Observer()
    my_location.lat = str(site_lat)
    my_location.lon = str(site_lon)
    my_location.pressure = 0
    my_location.horizon = '-0:34'  # per USNO definition
    
    #create grid time array for 1 reference year with leap day (year 2000)
    start = datetime.datetime(year = 2000, month = 1, day = 1, hour = 0, minute = 0)
    end = datetime.datetime(year = 2001, month = 1, day = 1, hour = 0, minute = 0)
    ref_date_dt = pd.date_range(start,end,freq='H')[:-1]
    ref_date_str = np.array([d.strftime('%m%d%H') for d in ref_date_dt]).astype('int')
    
    ref_state_array = []
    
    #for ref year determine for each time if day/night
    for i in range(len(ref_date_dt)):
    
        #update date
        current_dt = ref_date_dt[i]
        my_location.date = current_dt
        
        #get prev sunrise and sunset for current time
        try:
            next_sunrise = my_location.next_rising(ep.Sun()).datetime()
            next_sunset = my_location.next_setting(ep.Sun()).datetime()
        
            #determine if day or night
            if next_sunset < next_sunrise:
                #if delta between current datetime and next sunrise/sunset is < 30 mins
                #then set state is not representative of hour and need to correct
                #(i.e time could be 7am and sunrise is at 7.01am, state would be set as night, but really there are 59 minutes of day)
                if (next_sunset - current_dt) <  datetime.timedelta(minutes=30):
                    current_state = 'night'
                else:
                    current_state = 'day'
                    
            else:
                if (next_sunrise - current_dt) <  datetime.timedelta(minutes=30):
                    current_state = 'day'
                else: 
                    current_state = 'night'
                
            #append to ref state array
            ref_state_array.append(current_state)
    
        #exceptions for days when sun never rises or sets
        except ep.AlwaysUpError:
            current_state = 'day'
            ref_state_array.append(current_state)
        
        except ep.NeverUpError:
            current_state = 'night'
            ref_state_array.append(current_state)

        #print current_dt,current_state,next_sunrise,next_sunset

    ref_state_array = np.array(ref_state_array)
    
    #create str array from actual datetimes of month+day+hour
    orig_date_str = np.array([d.strftime('%m%d%H') for d in orig_datetime]).astype('int')
            
    #cut original timeseries now into day and night chunks
    #based on ref year determined states at each hour
    match_inds = np.searchsorted(ref_date_str,orig_date_str)
    ref_day_inds = np.where(ref_state_array == 'day')[0]
    ref_night_inds = np.where(ref_state_array == 'night')[0]
    
    day_inds = np.array([i for i, item in enumerate(match_inds) if item in ref_day_inds])
    night_inds = np.array([i for i, item in enumerate(match_inds) if item in ref_night_inds])
    
    return day_inds,night_inds

#3.
#----------

#SUPERPOSE FUNDAMENTAL AND HARMONIC WAVEFORMS

#THE WAY THE LSP/FFT returns the phase is not directly associated with sampling times
#I.E IF MEASURE EVERY HOUR, PHASE RETURNED IN RADIANS DOES NOT HAVE TO LINE UP WITH AN HOUR
#IT COULD BE BETWEEN HOURS, THIS IS BECAUSE THE PHASE IS SIMPLY WHERE IT MINIMISES THE LEAST SQUARES
#IF ONLY SAMPLE EVERY HOUR THIS ONLY HAS A MEANIGFUL EFFECT FOR THE DIURNAL CYCLE AS PHASE CAN BE BETWEEN HOURS

#SOLUTION: OVERSAMPLE CALCULATED WAVEFORM ENOUGH TO DETERMINE ACTUAL PHASE (AND THUS AMPLITUDE ALSO)
#THEN TAKE SUBSAMPLE OF OVERSAMPLED WAVEFORM TO GET  BACK TO SAMPLED TIMES  

def period_convolution(key_periods,time,mag,ph,mean):
    pi2 = np.pi*2
    
    if np.max(key_periods) == 1.:
        max_time = 1.
    elif np.max(key_periods) == 365.25:
        max_time = 365.25
    else:
        max_time = np.max(key_periods)
        
    #oversample time array (100 times size of original time array)
    oversampled_time = np.linspace(0,max_time,len(time)*100.,endpoint=False)
    
    for i in range(len(key_periods)):
        try:
            waveform = waveform + (mag[i]*(np.cos((pi2*oversampled_time/key_periods[i])-(ph[i]))))
        except:
            waveform = mag[i]*(np.cos((pi2*oversampled_time/key_periods[i])-(ph[i])))
    
    waveform = mean+waveform
    waveform_min = np.min(waveform)
    waveform_max = np.max(waveform)
    waveform_min_ind = np.argmin(waveform)
    waveform_max_ind = np.argmax(waveform)
    mag = (waveform_max-waveform_min)/2.
    ph_max_time = oversampled_time[waveform_max_ind]
    ratio = (pi2)/max_time
    ph_max_rads = ratio*ph_max_time
    
    #get waveform on sampled times by undersampling oversampled array
    waveform = waveform[np.searchsorted(oversampled_time,time)]
    
    #get phase on sampled timescale
    #ph_max_time = time[min(range(len(time)), key=lambda i: abs(time[i]-ph_max_time))]
    #ph_max_rads = ratio*ph_max_time
    
    return mag,ph_max_rads,waveform

#3.
#----------
    
#CONVERT UTC TIME PHASE TO SOLAR TIME PHASE (IN RADIANS)
    
def solar_time_phase_corrector(current_phase,hours,time_diff):
    pi2 = np.pi*2.
    ph_ratio = pi2/hours
    ph_offset = time_diff * ph_ratio
    
    new_phase = current_phase + ph_offset
    remainder = np.mod(np.abs(new_phase),pi2)
    
    if new_phase >= pi2:
        phase = remainder
    elif new_phase < 0:
        phase = pi2 - remainder
    else:
        phase = new_phase
    
    return phase

#3.09
#----------
    
#CONVERT PHASE (IN TIME) TO PHASE (RADIANS)
#RADIANS GO FROM 0 TO 2PI

def time_to_radians(obs,period):
    pi2 = np.pi*2.
    ratio = pi2/period
    obs = obs*ratio
   
    if obs>= np.pi:
        obs = -np.pi+(obs-np.pi)
   
    return obs

#CONVERT PHASE (IN RADIANS) TO PHASE (TIME)

#3.10
#----------

def radians_to_time(obs,con_num):			
    if obs < 0:
        obs = np.pi+ (np.pi-np.abs(obs))
    elif obs > 0:
		obs = obs

    convert_factor = con_num/(2*np.pi)
    obs = obs*convert_factor
	
    return obs
 
#3.11
#----------
 
def radians_to_time_spectra(obs,periods):			
    for i in range(len(obs)):
        if obs[i] < 0:
            obs[i] = np.pi+ (np.pi-np.abs(obs[i]))
        elif obs[i] > 0:
            obs[i] = obs[i]

        con_num = periods[i]

        convert_factor = con_num/(2*np.pi)
        obs[i] = obs[i]*convert_factor
	
    return obs

#3.12
#----------
    
#make time start from zero as spectral analysis is more accurate when this is the case
    
def time_from0(time):
    time = time-time[0]
    return time

#3.13
#----------
    
#phase is relative to start point in time, i.e if time starts in July
#the phase would be from July, we need to correct for this relativity
#so that start of daily and seasonal cycles are 00:00 hours and Jan 1 00:00 
#respectively

#PHASE INPUT/OUTPUT IS IN RADIANS
    
def phase_start_time_relative_correct(periods,phase,valid_times): 
    full_cycle = np.pi*2
    
    for i in range(len(periods)):
        #correct phases for offset from start point
        start_time = valid_times[0]
        if start_time >= periods[i]:
            while start_time >= periods[i]:
                start_time = start_time - periods[i]
        fraction_off = start_time/periods[i]
        offset = full_cycle*fraction_off
    
        phase[i] = phase[i] + offset

        if phase[i] > np.pi:
            phase[i] = -np.pi + (phase[i] - np.pi)
    
    return phase

#3.14
#----------
    
#TOOL THAT INTERPOLATES TO GET ACCURATE PEAK AMPLITUDE AND PHASE (IF PERIOD OF VARIABILITY KNOWN A PRIORI)
#INTERPOLATES USING SHAPE OF PEAK, DETERMINED BY WINDOWING METHOD
#BY DEFAULT IS HANNING
   
def periodic_interp(fr,fi,zoomfact,periods,key_period,n,amp_corr,window='hanning', alpha=6.0):
            
    closest_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-key_period))

    peak_fr = fr[closest_period_index-5:closest_period_index+6]
    peak_fi = fi[closest_period_index-5:closest_period_index+6]
    
    peak_periods = periods[closest_period_index-5:closest_period_index+6]
    
    peak_real = peak_fr/n
    peak_real = peak_real*2
    peak_real = peak_real*amp_corr

    peak_imag = peak_fi/n
    peak_imag = peak_imag*2
    peak_imag = peak_imag*amp_corr
    
    reverse_peak_periods = peak_periods[::-1]
    data = peak_real[::-1]
    data2 = peak_imag[::-1]
    
    #interpolate real component
    zoomfact = int(zoomfact)
    if (zoomfact < 1):
        print "zoomfact must be >= 1."
        return 0.0
    elif zoomfact==1:
        return data
    newN = len(data)*zoomfact
    # Space out the data
    comb = np.zeros((zoomfact, len(data)), dtype='d')
    comb[0] += data
    comb = np.reshape(np.transpose(comb), (newN,))
    # Compute the offsets
    xs = np.zeros(newN, dtype='d')
    xs[:newN/2+1] = np.arange(newN/2+1, dtype='d')/zoomfact
    xs[-newN/2:]  = xs[::-1][newN/2-1:-1]
    # Calculate the sinc times window for the kernel
    if window.lower()=="kaiser":
        win = _window_function[window](xs, len(data)/2, alpha)
    else:
        win = hanning_window(xs, len(data)/2)
    kernel = win * sinc(xs)
    newreal = FFT.irfft(FFT.rfft(kernel) * FFT.rfft(comb))
    
    #interpolate imag component
    if (zoomfact < 1):
        print "zoomfact must be >= 1."
        return 0.0
    elif zoomfact==1:
        return data2
    newN = len(data2)*zoomfact
    # Space out the data
    comb = np.zeros((zoomfact, len(data2)), dtype='d')
    comb[0] += data2
    comb = np.reshape(np.transpose(comb), (newN,))
    # Compute the offsets
    xs = np.zeros(newN, dtype='d')
    xs[:newN/2+1] = np.arange(newN/2+1, dtype='d')/zoomfact
    xs[-newN/2:]  = xs[::-1][newN/2-1:-1]
    # Calculate the sinc times window for the kernel
    if window.lower()=="kaiser":
        win = _window_function[window](xs, len(data2)/2, alpha)
    else:
        win = hanning_window(xs, len(data2)/2)
    kernel = win * sinc(xs)
    newimag = FFT.irfft(FFT.rfft(kernel) * FFT.rfft(comb))
    
    #get new periods for interpolation
    new_x = []
    ns = 0
    ne = 1
    over_its = len(reverse_peak_periods)-1
    for i in range(over_its):
        new_x = np.append(new_x,np.linspace(reverse_peak_periods[ns],reverse_peak_periods[ne],zoomfact,endpoint=False))
        ns+=1
        ne+=1

    new_x = np.append(new_x,reverse_peak_periods[-1])

    #get max mag in selected range (and associated phase)
    newreal = newreal[:len(new_x)]
    newimag = newimag[:len(new_x)]
    newmag = []
    newphase = []
    for i in range(len(newreal)):
        newmag.append(np.abs(complex(newreal[i],newimag[i])))
        newphase.append(np.angle(complex(newreal[i],newimag[i])))
    peakind = np.argmax(newmag)

    amp = newmag[peakind]
    phase = newphase[peakind]
    
    return amp,phase

#3.15
#----------
    
#REMOVE KNOWN PERIODIC COMPONENTS FROM REAL AND IMAG COMPONENTS SPECTRA
    
def remove_periodic(bp_periods,bp_mag,ofac):
    
    rm_periods = [1./4,1./3,1./2,1.,365.25/4,365.25/3,365.25/2,365.25]
    rm_n_points = [200,150,100,50,4,3,1,1]
        
    rm_points = []
    
    for i in range(len(rm_periods)):
        rm_i = min(range(len(bp_periods)), key=lambda x: abs(bp_periods[x]-rm_periods[i]))
        n_points = int((rm_n_points[i]*ofac)+np.floor(ofac/2.))
        rm_points_solo = range(rm_i-n_points,(rm_i+n_points)+1)
        if rm_points_solo[0] < 0:
            rm_points_solo = range(0,(rm_i+n_points)+1)
        rm_points = np.append(rm_points,rm_points_solo)
    
    #can either interpolate spectra using points either side of cut, set points to be 0 or delete
    #if going to delete points need to delete equivalent periods too
    
    rm_points = rm_points.astype(int)
    
    #1.
    #bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    #2.
    bp_mag[rm_points] = 0
    #3.
    #bp_mag = np.delete(bp_mag,rm_points)
    #bp_periods = np.delete(bp_periods,rm_points)
    
    return bp_periods,bp_mag

#3.16
#----------
    
#PROGRAM THAT CALCULATES AND WRITES OUT DIFFERENCES IN PERIODIC COMPONENTS BETWEEN OBSERVATIONS AND MODEL
   
def bulge_calc(obs_refs,obs_waveforms,model_waveforms,max_ph):
    #save out data to netcdf
    root_grp_period = Dataset('bulge.nc', 'w')
    root_grp_period.description = 'Bulge info found between obs and model - Program written by Dene Bowdalo'
        
    for r in range(len(obs_refs)):
        ref = obs_refs[r]
        
        print '\n',ref
        
        ref_period = root_grp_period.createGroup('%s'%(ref.lower()))
    
        obs_waveform = obs_waveforms[r]
        model_waveform = model_waveforms[r]
        waveform_len = len(obs_waveform)
        
        #process date and time, make cuts of time for day and year respectively    
        #date_l = date.astype(int)
        #time_l = time.astype(int)
        #full_times = date_process(date_l,time_l,start_year)
        #full_times_year = full_times[:8766]
        #full_times_day = full_times[:24]
        
        #get spacing
        #spacing = full_times_day[1] - full_times_day[0]
        

        #calc residual
        residual = model_waveform - obs_waveform
        
        #make ref bulge arrays
        b_sign = []
        b_maxconcoffset = []
        b_maxconcoffsetph= []
        b_aveconcoffset = []
        b_totalconcoffset = []
        b_startphaseoffset = []
        b_endphaseoffset = []
        b_totalphaseoffset = []
        b_pcbulgearea = []
        b_pcbulgetime = []
        
        #calc summary stats for period
        #) N bulges
        #) N + bulges
        #) N - bulges
        #) Ave conc. offset (+ or -) over whole period
        #) Total conc. offset over whole period (abs)
        #) Ave + conc. offset over whole period
        #) Total + area bulge over period
        #) Ave - conc. offset over whole period
        #) Total - area bulge over period
        #) Percent + area bulge of all bulges
        #) Percent - area bulge of all bulges
        #) Percent + area bulge of time
        #) Percent - area bulge of time
        
        ph_ratio = max_ph/waveform_len
        
        a_aveconcoffset = np.average(residual)
        a_totalconcoffset = np.sum(np.abs(residual))
        
        test = residual > 0
        residual_plus = residual[test]
        residual_plus_len = len(residual_plus)
        a_aveplusconcoffset = np.average(residual_plus)
        a_totalplusconcoffset = np.sum(residual_plus)
        
        test = residual < 0
        residual_neg = residual[test]
        residual_neg_len = len(residual_neg)
        a_avenegconcoffset = np.average(residual_neg)
        a_totalnegconcoffset = np.sum(residual_neg)  
        
        a_pcplusbulgearea = (100./a_totalconcoffset)*a_totalplusconcoffset
        a_pcplusbulgetime = ph_ratio * residual_plus_len
        a_pcnegbulgearea = (100./a_totalconcoffset)*np.abs(a_totalnegconcoffset)
        a_pcnegbulgetime = ph_ratio *residual_neg_len
        
        
        #find first point in residual where there is a zero crossover and iterative inds from that
        zero_crossings = np.where(np.diff(np.signbit(residual)))[0]
        
        #zero crossings will be zero if there is only 1 bulge, thus test for this
        if len(zero_crossings) > 0:
            zero_crossing_first = zero_crossings[0]+1
        else:
            zero_crossing_first = 0
        
        iter_inds = range(zero_crossing_first,len(residual),1) + range(0,zero_crossing_first,1)
        
        #calc different bulges in turn and write out stats.
        bulge_out = False
        current_sign = ''
        nbulges = 0
        nplusbulges = 0
        nnegbulges = 0
        
        last_ind = iter_inds[-1]
        
        for i in iter_inds:
            
            #if ind is last ind then finish by writing out bulge on currently
            if i == last_ind:
                save_start_i = np.copy(start_i)
                end_i = i+1
                bulge_out = True
                nbulges+=1
                if current_sign == '+':
                    nplusbulges+=1
                    bulge_flag = '+'
                else:
                    nnegbulges+=1
                    bulge_flag = '-'
                current_sign = 'na'
                print nbulges
            
            if current_sign == '+':
                if residual[i] > 0:
                    pass
                else:
                    save_start_i = np.copy(start_i)
                    end_i = i
                    bulge_out = True
                    nbulges+=1
                    nplusbulges+=1
                    current_sign = ''
                    bulge_flag = '+'
                    print nbulges
                    
            elif current_sign == '-':
                if residual[i] < 0:
                    pass
                else:
                    save_start_i = np.copy(start_i)
                    end_i = i
                    bulge_out = True
                    nbulges+=1
                    nnegbulges+=1
                    current_sign = ''
                    bulge_flag = '-'
                    print nbulges
            
            if current_sign == '':
                #get sign for next i 
                try:
                    if residual[i] > 0:
                        current_sign = '+'
                        start_i = i
                    elif residual[i] < 0:
                        current_sign = '-'
                        start_i = i
                    else:
                        current_sign = ''
                    
                except:
                    pass
                    
            if bulge_out == True:
                print bulge_flag
            
                #calc specific bulge stats
                #) Max conc offset
                #) Point of max conc offset (phase)
                #) Average conc offset
                #) Total conc offset
                #) Start phase offset
                #) End phase offset
                #) Total phase offset
                #) Percent of total bulge taken up by current bulge
                #) Percent of time taken up by bulge
        
                print save_start_i,end_i
                #make array of inds to take from residual from start_i and end_i
                if end_i < save_start_i:
                    ind_array = range(save_start_i,len(residual),1) + range(0,end_i,1)
                else:
                    ind_array = range(save_start_i,end_i,1) 
        
        
                b_sign.append(bulge_flag)
                b_maxconcoffset.append(np.max(np.abs(residual[ind_array])))
                ind_max = ind_array[np.argmax(np.abs(residual[ind_array]))]
                #print ph_ratio, ind_max
                b_maxconcoffsetph.append(ph_ratio * ind_max) 
                b_aveconcoffset.append(np.average(np.abs(residual[ind_array])))
                b_totalconcoffset.append(np.sum(np.abs(residual[ind_array])))
                b_startphaseoffset.append(ph_ratio*save_start_i)
                b_endphaseoffset.append(ph_ratio*end_i)
                b_totalphaseoffset.append(b_endphaseoffset[-1] - b_startphaseoffset[-1])
                b_pcbulgearea.append((100./a_totalconcoffset) * b_totalconcoffset[-1])
                b_pcbulgetime.append((max_ph/100.) * b_totalphaseoffset[-1])
                
                bulge_out = False
    

        for nb in range(nbulges):
            #write out from biggest total bulge to smallest 
            
            #sort all bulge arrays based on total bulge (biggest to smallest)
            sorted_inds = sorted(range(len(b_totalconcoffset)),key=lambda x:b_totalconcoffset[x])
            sorted_inds = sorted_inds[::-1]
            
            b_sign = np.array(b_sign)
            b_maxconcoffset = np.array(b_maxconcoffset)
            b_maxconcoffsetph = np.array(b_maxconcoffsetph)
            b_aveconcoffset = np.array(b_aveconcoffset)
            b_totalconcoffset = np.array(b_totalconcoffset)
            b_startphaseoffset = np.array(b_startphaseoffset)
            b_endphaseoffset = np.array(b_endphaseoffset)
            b_totalphaseoffset = np.array(b_totalphaseoffset)
            b_pcbulgearea = np.array(b_pcbulgearea)
            b_pcbulgetime = np.array(b_pcbulgetime)
            
            b_sign = b_sign[sorted_inds]
            b_maxconcoffset = b_maxconcoffset[sorted_inds]                
            b_maxconcoffsetph = b_maxconcoffsetph[sorted_inds]
            b_aveconcoffset = b_aveconcoffset[sorted_inds]
            b_totalconcoffset = b_totalconcoffset[sorted_inds]
            b_startphaseoffset = b_startphaseoffset[sorted_inds]
            b_endphaseoffset = b_endphaseoffset[sorted_inds]
            b_totalphaseoffset = b_totalphaseoffset[sorted_inds]
            b_pcbulgearea = b_pcbulgearea[sorted_inds]
            b_pcbulgetime = b_pcbulgetime[sorted_inds]
            
            bulge_grp = ref_period.createGroup('bulge%s'%(nb+1))

            bulge_grp.b_sign = b_sign[nb]
            bulge_grp.b_maxconcoffset = b_maxconcoffset[nb]
            bulge_grp.b_maxconcoffsetph = b_maxconcoffsetph[nb]
            bulge_grp.b_aveconcoffset = b_aveconcoffset[nb]
            bulge_grp.b_totalconcoffset = b_totalconcoffset[nb]
            bulge_grp.b_startphaseoffset = b_startphaseoffset[nb]
            bulge_grp.b_endphaseoffset = b_endphaseoffset[nb]
            bulge_grp.b_totalphaseoffset = b_totalphaseoffset[nb]
            bulge_grp.b_pcbulgearea = b_pcbulgearea[nb]
            bulge_grp.b_pcbulgetime = b_pcbulgetime[nb]
            

        #write out summary stats
        ref_period.nbulges = nbulges
        ref_period.nplusbulges = nplusbulges
        ref_period.nnegbulges = nnegbulges

        ref_period.a_aveconcoffset = a_aveconcoffset
        ref_period.a_aveplusconcoffset = a_aveplusconcoffset
        ref_period.a_totalplusconcoffset = a_totalplusconcoffset
        ref_period.a_avenegconcoffset = a_avenegconcoffset
        ref_period.a_totalnegconcoffset = a_totalnegconcoffset
        
        ref_period.a_pcplusbulgearea = a_pcplusbulgearea
        ref_period.a_pcnegbulgearea = a_pcnegbulgearea
        ref_period.a_pcplusbulgetime = a_pcplusbulgetime
        ref_period.a_pcnegbulgetime = a_pcnegbulgetime
        
    root_grp_period.close()
        
    return

#3.17
#----------
    
def maxmin_calc(d_waveforms,s_waveforms):
    #save out data to netcdf
    root_grp_period = Dataset('maxmin.nc', 'w')
    root_grp_period.description = 'MaxMin info - Program written by Dene Bowdalo'
    
    root_grp_period.createDimension('waveform_len',len(d_waveforms))
    
    d_n_peaks = []
    d_peak1_conc = []
    d_peak2_conc = []
    d_peak3_conc = []
    d_peak4_conc = []
    d_peak1_ph = []
    d_peak2_ph = []
    d_peak3_ph = []
    d_peak4_ph = []
    s_n_peaks = []
    s_peak1_conc = []
    s_peak2_conc = []
    s_peak3_conc = []
    s_peak4_conc = []
    s_peak1_ph = []
    s_peak2_ph = []
    s_peak3_ph = []
    s_peak4_ph = []
    
    print 'daily processing'
    #calculate diurnal
    max_ph = 24.
    ph_ratio = max_ph/len(d_waveforms[0])

    for wf in d_waveforms:
        peak_conc = []
        peak_ph = []
        
        n_peaks = 1
    
        max_ind = np.argmax(wf)
        iter_inds_prev = range(max_ind,len(wf),1) + range(0,max_ind-1,1)
        iter_inds_next = range(max_ind+2,len(wf),1) + range(0,max_ind+1,1)
        iter_inds = range(max_ind+1,len(wf),1) + range(0,max_ind,1)
        peak_conc.append(wf[max_ind])
        peak_ph.append(ph_ratio*max_ind)
        for i in range(len(iter_inds)):
            if (wf[iter_inds_prev[i]] < wf[iter_inds[i]]) & (wf[iter_inds_next[i]] < wf[iter_inds[i]]):
                peak_conc.append(wf[iter_inds[i]])
                peak_ph.append(ph_ratio*iter_inds[i])
                n_peaks+=1
        
        if n_peaks == 1:
            d_n_peaks.append(1)
            d_peak1_ph.append(peak_ph[0])
            d_peak2_ph.append(-99999)
            d_peak3_ph.append(-99999)
            d_peak4_ph.append(-99999)
            d_peak1_conc.append(peak_conc[0])
            d_peak2_conc.append(-99999)
            d_peak3_conc.append(-99999)
            d_peak4_conc.append(-99999)
        
        elif n_peaks == 2:
            d_n_peaks.append(2)
            d_peak1_ph.append(peak_ph[0])
            d_peak2_ph.append(peak_ph[1])
            d_peak3_ph.append(-99999)
            d_peak4_ph.append(-99999)
            d_peak1_conc.append(peak_conc[0])
            d_peak2_conc.append(peak_conc[1])
            d_peak3_conc.append(-99999)
            d_peak4_conc.append(-99999)
        
        elif n_peaks == 3:
            d_n_peaks.append(3)
            d_peak1_ph.append(peak_ph[0])
            d_peak2_ph.append(peak_ph[1])
            d_peak3_ph.append(peak_ph[2])
            d_peak4_ph.append(-99999)
            d_peak1_conc.append(peak_conc[0])
            d_peak2_conc.append(peak_conc[1])
            d_peak3_conc.append(peak_conc[2])
            d_peak4_conc.append(-99999)
        
        elif n_peaks == 4:
            d_n_peaks.append(4)
            d_peak1_ph.append(peak_ph[0])
            d_peak2_ph.append(peak_ph[1])
            d_peak3_ph.append(peak_ph[2])
            d_peak4_ph.append(peak_ph[3])
            d_peak1_conc.append(peak_conc[0])
            d_peak2_conc.append(peak_conc[1])
            d_peak3_conc.append(peak_conc[2])
            d_peak4_conc.append(peak_conc[3])
        
    print 'seasonal processing'
    #calculate seasonal
    max_ph = 12.
    ph_ratio = max_ph/len(s_waveforms[0])
    
    for wf in s_waveforms:
        n_peaks = 1
        peak_conc = []
        peak_ph = []
    
        max_ind = np.argmax(wf)
        iter_inds_prev = range(max_ind,len(wf),1) + range(0,max_ind-1,1)
        iter_inds_next = range(max_ind+2,len(wf),1) + range(0,max_ind+1,1)
        iter_inds = range(max_ind+1,len(wf),1) + range(0,max_ind,1)
        peak_conc.append(wf[max_ind])
        peak_ph.append(ph_ratio*max_ind)
        for i in range(len(iter_inds)):
            if (wf[iter_inds_prev[i]] < wf[iter_inds[i]]) & (wf[iter_inds_next[i]] < wf[iter_inds[i]]):
                peak_conc.append(wf[iter_inds[i]])
                peak_ph.append(ph_ratio*iter_inds[i])
                n_peaks+=1
        
            
        if n_peaks == 1:
            s_n_peaks.append(1)
            s_peak1_ph.append(peak_ph[0])
            s_peak2_ph.append(-99999)
            s_peak3_ph.append(-99999)
            s_peak4_ph.append(-99999)
            s_peak1_conc.append(peak_conc[0])
            s_peak2_conc.append(-99999)
            s_peak3_conc.append(-99999)
            s_peak4_conc.append(-99999)
        
        elif n_peaks == 2:
            s_n_peaks.append(2)
            s_peak1_ph.append(peak_ph[0])
            s_peak2_ph.append(peak_ph[1])
            s_peak3_ph.append(-99999)
            s_peak4_ph.append(-99999)
            s_peak1_conc.append(peak_conc[0])
            s_peak2_conc.append(peak_conc[1])
            s_peak3_conc.append(-99999)
            s_peak4_conc.append(-99999)
        
        elif n_peaks == 3:
            s_n_peaks.append(3)
            s_peak1_ph.append(peak_ph[0])
            s_peak2_ph.append(peak_ph[1])
            s_peak3_ph.append(peak_ph[2])
            s_peak4_ph.append(-99999)
            s_peak1_conc.append(peak_conc[0])
            s_peak2_conc.append(peak_conc[1])
            s_peak3_conc.append(peak_conc[2])
            s_peak4_conc.append(-99999)
        
        elif n_peaks == 4:
            s_n_peaks.append(4)
            s_peak1_ph.append(peak_ph[0])
            s_peak2_ph.append(peak_ph[1])
            s_peak3_ph.append(peak_ph[2])
            s_peak4_ph.append(peak_ph[3])
            s_peak1_conc.append(peak_conc[0])
            s_peak2_conc.append(peak_conc[1])
            s_peak3_conc.append(peak_conc[2])
            s_peak4_conc.append(peak_conc[3])
            
            
    d_n_p = root_grp_period.createVariable('diurnal_n_peaks', 'f8', ('waveform_len',))
    d_p1_ph = root_grp_period.createVariable('diurnal_peak1_phase', 'f8', ('waveform_len',))
    d_p2_ph = root_grp_period.createVariable('diurnal_peak2_phase', 'f8', ('waveform_len',))
    d_p3_ph = root_grp_period.createVariable('diurnal_peak3_phase', 'f8', ('waveform_len',))
    d_p4_ph = root_grp_period.createVariable('diurnal_peak4_phase', 'f8', ('waveform_len',))
    d_p1_conc = root_grp_period.createVariable('diurnal_peak1_conc', 'f8', ('waveform_len',))
    d_p2_conc = root_grp_period.createVariable('diurnal_peak2_conc', 'f8', ('waveform_len',))
    d_p3_conc = root_grp_period.createVariable('diurnal_peak3_conc', 'f8', ('waveform_len',))
    d_p4_conc = root_grp_period.createVariable('diurnal_peak4_conc', 'f8', ('waveform_len',))           
    s_n_p = root_grp_period.createVariable('seasonal_n_peaks', 'f8', ('waveform_len',))
    s_p1_ph = root_grp_period.createVariable('seasonal_peak1_phase', 'f8', ('waveform_len',))
    s_p2_ph = root_grp_period.createVariable('seasonal_peak2_phase', 'f8', ('waveform_len',))
    s_p3_ph = root_grp_period.createVariable('seasonal_peak3_phase', 'f8', ('waveform_len',))
    s_p4_ph = root_grp_period.createVariable('seasonal_peak4_phase', 'f8', ('waveform_len',))
    s_p1_conc = root_grp_period.createVariable('seasonal_peak1_conc', 'f8', ('waveform_len',))
    s_p2_conc = root_grp_period.createVariable('seasonal_peak2_conc', 'f8', ('waveform_len',))
    s_p3_conc = root_grp_period.createVariable('seasonal_peak3_conc', 'f8', ('waveform_len',))
    s_p4_conc = root_grp_period.createVariable('seasonal_peak4_conc', 'f8', ('waveform_len',))
    
    d_n_p[:] = d_n_peaks
    d_p1_ph[:] = d_peak1_ph
    d_p2_ph[:] = d_peak2_ph
    d_p3_ph[:] = d_peak3_ph
    d_p4_ph[:] = d_peak4_ph
    d_p1_conc[:] = d_peak1_conc
    d_p2_conc[:] = d_peak2_conc
    d_p3_conc[:] = d_peak3_conc
    d_p4_conc[:] = d_peak4_conc
    s_n_p[:] = s_n_peaks
    s_p1_ph[:] = s_peak1_ph
    s_p2_ph[:] = s_peak2_ph
    s_p3_ph[:] = s_peak3_ph
    s_p4_ph[:] = s_peak4_ph
    s_p1_conc[:] = s_peak1_conc
    s_p2_conc[:] = s_peak2_conc
    s_p3_conc[:] = s_peak3_conc
    s_p4_conc[:] = s_peak4_conc

        
    root_grp_period.close()
        
    return
  
#3.18
#----------
    
#SPECTRA MODEL FITTING
#-------------------------------------------------------------

#FUNCTION FOR FITTING OF 2 FIXED LINEAR PIECEWISE FUNCTIONS TO SPECTRA
 
def spectra_f_2(x,grad1,grad2,y_inter):

    x1_r = x[0]
    x2_r = x[1]
   
    diff1 = x2_r[0] - x1_r[-1]
   
    x1_r = x1_r - x1_r[0]
    x2_r = x2_r - x2_r[0]
   
    p1 = grad1*x1_r+y_inter
    y_inter2 = grad1*(x1_r[-1]+diff1)+y_inter
    p2 = grad2*x2_r+y_inter2

    p = [p1]+[p2]
    p = [item for sublist in p for item in sublist]
    
    return p   

#3.19
#----------
    
#FUNCTION FOR FITTING OF 3 FIXED BREAKPOINT LINEAR PIECEWISE FUNCTIONS TO SPECTRA

def spectra_f_3(x,grad1,grad2,grad3,y_inter):

    x1_r = x[0]
    x2_r = x[1]
    x3_r = x[2]  
   
    diff1 = x2_r[0] - x1_r[-1]
    diff2 = x3_r[0] - x2_r[-1]
   
    x1_r = x1_r - x1_r[0]
    x2_r = x2_r - x2_r[0]
    x3_r = x3_r - x3_r[0]
   
    p1 = grad1*x1_r+y_inter
    y_inter2 = grad1*(x1_r[-1]+diff1)+y_inter
    p2 = grad2*x2_r+y_inter2
    y_inter3 = grad2*(x2_r[-1]+diff2)+y_inter2
    p3 = grad3*x3_r+y_inter3

    p = [p1]+[p2]+[p3]
    p = [item for sublist in p for item in sublist]
    
    return p

#3.20
#----------
    
#FUNCTION FOR FITTING OF 2 NON-FIXED LINEAR PIECEWISE LINEAR FUNCTIONS TO SPECTRA

def spectra_f_parts_2(x,grad1,grad2,y_inter):
    
    x1_r = x[0]
    x2_r = x[1]
   
    diff1 = x2_r[0] - x1_r[-1]
   
    x1_r = x1_r - x1_r[0]
    x2_r = x2_r - x2_r[0]
   
    p1 = grad1*x1_r+y_inter
    y_inter2 = grad1*(x1_r[-1]+diff1)+y_inter
    p2 = grad2*x2_r+y_inter2
  
    return p1,p2

#3.21
#----------
    
#FUNCTION FOR FITTING OF 3 NON-FIXED PIECEWISE LINEAR FUNCTIONS TO SPECTRA

def spectra_f_parts_3(x,grad1,grad2,grad3,y_inter):
    
    x1_r = x[0]
    x2_r = x[1]
    x3_r = x[2]  
   
    diff1 = x2_r[0] - x1_r[-1]
    diff2 = x3_r[0] - x2_r[-1]
   
    x1_r = x1_r - x1_r[0]
    x2_r = x2_r - x2_r[0]
    x3_r = x3_r - x3_r[0]
   
    p1 = grad1*x1_r+y_inter
    y_inter2 = grad1*(x1_r[-1]+diff1)+y_inter
    p2 = grad2*x2_r+y_inter2
    y_inter3 = grad2*(x2_r[-1]+diff2)+y_inter2
    p3 = grad3*x3_r+y_inter3
  
    return p1,p2,p3
    
#3.22
#----------

#FIT FIXED PIECEWISE MODEL TO SPECTRA (CURRENTLY USES 2 PIECES - REGIMES ARE WEATHER AND MACROWEATHER)
    
def spectra_fit_fixed_piecewise(periods,mag,ofac,lower_period_limit,upper_period_limit,breakpoint):

    start_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-upper_period_limit))
    end_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-lower_period_limit))
    
    if periods[end_cut_i] > lower_period_limit:
        end_cut_i = end_cut_i + 1.
    bp_periods = periods[start_cut_i:end_cut_i]
    bp_mag = mag[start_cut_i:end_cut_i]
    
    if (lower_period_limit < 1./10.) & (1./10. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./10.))
        n_points = int((500*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./9.) & (1./9. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./9.))
        n_points = int((450*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./8.) & (1./8. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./8.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./7.) & (1./7. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./7.))
        n_points = int((350*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./6.) & (1./6. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./6.))
        n_points = int((300*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./5.) & (1./5. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./5.))
        n_points = int((250*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./4.) & (1./4. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./4.))
        n_points = int((200*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./3.) & (1./3. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./3.))
        n_points = int((150*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./2.) & (1./2. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./2.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 1.) & (1. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/4.) & (365.25/4. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/4.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 365.25/3.) & (365.25/3. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/3.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/2.) & (365.25/2. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/2.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
     
    if (lower_period_limit < 365.25) & (365.25 < upper_period_limit):
        annual_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(annual_i-n_points,(annual_i+n_points)+1)
        if rm_points[0] < 0:
            rm_points = range(0,(annual_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    #make periods (and mag) go in ascending order.
    bp_periods = bp_periods[::-1]
    bp_mag = bp_mag[::-1]
    
    #log periods and mag, will keep log distributed nature after interpolation
    bp_periods = np.log(bp_periods)
    bp_mag = np.log(bp_mag)
    
    bp1 = np.log(breakpoint)
    
    test = bp_periods < bp1

    p1 = bp_periods[test]
    bp_mag1 = bp_mag[test]
    len1 = len(p1)
    test = bp_periods >= bp1
    p2 = bp_periods[test]
    bp_mag2 = bp_mag[test]
    len2 = len(p2)
    
    x_lens = [p1,p2]
    
    popt, pcov = curve_fit(spectra_f_2,x_lens,bp_mag)
    
    grad1 = popt[0]
    grad2 = popt[1]
    yinter = popt[2]

    full_periods = np.log(periods[::-1])
    test = full_periods < bp1
    p1 = full_periods[test]
    test = full_periods >= bp1
    p2 = full_periods[test]
    x_lens = [p1,p2]
    
    calc_1,calc_2 = spectra_f_parts_2(x_lens,grad1,grad2,yinter)
    
    calc_1 = calc_1.tolist()
    calc_2 = calc_2.tolist()
    
    calc_y = calc_1+calc_2 
    
    #calculate stats of fitted lines
    #-----------------------------------------
  
    ave1 = np.average(calc_1)
    ave2 = np.average(calc_2)
    med1 = np.median(calc_1)
    med2 = np.median(calc_2)
    sum1 = np.sum(calc_1)
    sum2 = np.sum(calc_2)
    
    y1_s = calc_1[0]
    y1_e = calc_1[-1]
    y2_s = calc_2[0]
    y2_e = calc_2[-1]
    
    plt.loglog(periods,mag)
    plt.loglog(np.exp(p1),np.exp(calc_1))
    plt.loglog(np.exp(p2),np.exp(calc_2))
    plt.show()
    
    return grad1,grad2,np.exp(bp1),np.exp(p1),np.exp(calc_1),np.exp(p2),np.exp(calc_2),np.exp(ave1),np.exp(ave2),np.exp(med1),np.exp(med2),np.exp(sum1),np.exp(sum2),np.exp(y1_s),np.exp(y1_e),np.exp(y2_s),np.exp(y2_e)

#3.23
#----------
    
#FIT FIXED SEPARATE LINEAR REGRESSION FUNCTIONS TO SPECTRA (NON-PIECEWISE) (CURRENTLY USES 2 PIECES - REGIMES ARE WEATHER AND MACROWEATHER)

def spectra_fit_fixed_linear(periods,mag,ofac,lower_period_limit,upper_period_limit,breakpoint):

    start_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-upper_period_limit))
    end_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-lower_period_limit))
    
    if periods[end_cut_i] > lower_period_limit:
        end_cut_i = end_cut_i + 1.
    bp_periods = periods[start_cut_i:end_cut_i]
    bp_mag = mag[start_cut_i:end_cut_i]
    
    if (lower_period_limit < 1./10.) & (1./10. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./10.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./9.) & (1./9. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./9.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./8.) & (1./8. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./8.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./7.) & (1./7. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./7.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./6.) & (1./6. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./6.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./5.) & (1./5. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./5.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./4.) & (1./4. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./4.))
        n_points = int((320*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./3.) & (1./3. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./3.))
        n_points = int((240*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./2.) & (1./2. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./2.))
        n_points = int((160*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit >= 1.) & (1. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit >= 365.25/4.) & (365.25/4. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/4.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit >= 365.25/3.) & (365.25/3. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/3.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit >= 365.25/2.) & (365.25/2. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/2.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
     
    if (lower_period_limit >= 365.25) & (365.25 < upper_period_limit):
        annual_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(annual_i-n_points,(annual_i+n_points)+1)
        if rm_points[0] < 0:
            rm_points = range(0,(annual_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    #make periods (and mag) go in ascending order.
    bp_periods = bp_periods[::-1]
    bp_mag = bp_mag[::-1]
    
    #log periods and mag, will keep log distributed nature after interpolation
    bp_periods = np.log(bp_periods)
    bp_mag = np.log(bp_mag)
    
    bp1 = np.log(breakpoint)

    test = bp_periods < bp1
    p1 = bp_periods[test]
    bp_mag1 = bp_mag[test]
    test = bp_periods >= bp1
    p2 = bp_periods[test]
    bp_mag2 = bp_mag[test]
    
    grad1, yinter1, r_value, p_value, std_err = stats.linregress(p1,bp_mag1)
    grad2, yinter2, r_value, p_value, std_err = stats.linregress(p2,bp_mag2)
    
    full_periods = np.log(periods[::-1])
    test = full_periods < bp1
    p1 = full_periods[test]
    test = full_periods >= bp1
    p2 = full_periods[test]   
    
    calc_1 = grad1*p1+yinter1
    calc_2 = grad2*p2+yinter2
    
    calc_1 = calc_1.tolist()
    calc_2 = calc_2.tolist()
    
    calc_y = calc_1+calc_2 
    
    orig_periods_log = np.log(periods)[::-1]
    
    #calculate stats of fitted lines
    #-----------------------------------------
  
    ave1 = np.average(calc_1)
    ave2 = np.average(calc_2)
    med1 = np.median(calc_1)
    med2 = np.median(calc_2)
    sum1 = np.sum(calc_1)
    sum2 = np.sum(calc_2)
    
    y1_s = calc_1[0]
    y1_e = calc_1[-1]
    y2_s = calc_2[0]
    y2_e = calc_2[-1]
    
    return grad1,grad2,np.exp(bp1),np.exp(p1),np.exp(calc_1),np.exp(p2),np.exp(calc_2),np.exp(ave1),np.exp(ave2),np.exp(med1),np.exp(med2),np.exp(sum1),np.exp(sum2),np.exp(y1_s),np.exp(y1_e),np.exp(y2_s),np.exp(y2_e)

#3.24
#----------
    
#FIT FREE PIECEWISE MODEL TO SPECTRA (CURRENTLY USES 2 PIECES - REGIMES ARE WEATHER AND MACROWEATHER)

def spectra_fit_free_piecewise(periods,mag,ofac,lower_period_limit,upper_period_limit,lower_break,upper_break):
    
    start_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-upper_period_limit))
    end_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-lower_period_limit))
    
    if periods[end_cut_i] > lower_period_limit:
        end_cut_i = end_cut_i + 1.
    bp_periods = periods[start_cut_i:end_cut_i]
    bp_mag = mag[start_cut_i:end_cut_i]
    
    if (lower_period_limit < 1./10.) & (1./10. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./10.))
        n_points = int((500*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./9.) & (1./9. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./9.))
        n_points = int((450*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./8.) & (1./8. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./8.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./7.) & (1./7. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./7.))
        n_points = int((350*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./6.) & (1./6. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./6.))
        n_points = int((300*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./5.) & (1./5. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./5.))
        n_points = int((250*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./4.) & (1./4. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./4.))
        n_points = int((200*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./3.) & (1./3. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./3.))
        n_points = int((150*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./2.) & (1./2. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./2.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 1.) & (1. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/4.) & (365.25/4. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/4.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 365.25/3.) & (365.25/3. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/3.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/2.) & (365.25/2. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/2.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
     
    if (lower_period_limit < 365.25) & (365.25 < upper_period_limit):
        annual_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(annual_i-n_points,(annual_i+n_points)+1)
        if rm_points[0] < 0:
            rm_points = range(0,(annual_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    #make periods (and mag) go in ascending order.
    bp_periods = bp_periods[::-1]
    bp_mag = bp_mag[::-1]
    
    #log periods and mag, will keep log distributed nature after interpolation
    bp_periods = np.log(bp_periods)
    bp_mag = np.log(bp_mag)
    
    #set limiting ranges for breakpoints
    breakpoint_test = (bp_periods >= np.log(lower_break)) & (bp_periods < np.log(upper_break))
    bp_periods_1 = bp_periods[breakpoint_test]
    bp1 = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-bp_periods_1[0]))

    all_error1 = []
    all_error2 = []
    all_bp1 = []
    all_grad1 = []
    all_grad2 = []
    all_yinter = []

    for i in range(len(bp_periods_1)):
        bp1 = bp_periods_1[i]
    
        test = bp_periods < bp1
        bp_periods1 = bp_periods[test]
        bp_mag1 = bp_mag[test]
        len1 = len(bp_periods1)
        test = bp_periods >= bp1
        bp_periods2 = bp_periods[test]
        bp_mag2 = bp_mag[test]
        len2 = len(bp_periods2)
    
        x_lens = [bp_periods1,bp_periods2]
    
        popt, pcov = curve_fit(spectra_f_2,x_lens,bp_mag)
    
        grad1 = popt[0]
        grad2 = popt[1]
        yinter = popt[2]
    
        bp_mag_est = spectra_f_2(x_lens,grad1,grad2,yinter)    
        
        line1,line2 = spectra_f_parts_2(x_lens,grad1,grad2,yinter)
        
        err1=np.sqrt(np.sum((np.abs(line1-bp_mag1))**2)/len(bp_mag1))
        err2=np.sqrt(np.sum((np.abs(line2-bp_mag2))**2)/len(bp_mag2))
    
        all_error1.append(err1)
        all_error2.append(err2)
        all_bp1.append(np.average((bp_periods1[-1],bp_periods2[0])))
        all_grad1.append(grad1)
        all_grad2.append(grad2)
        all_yinter.append(yinter)
        
    
    #get best fit based on best rank of errors for each line
    all_error1 = np.array(all_error1)
    all_error2 = np.array(all_error2)

    temp = all_error1.argsort()
    error1_ranks = np.empty(len(all_error1), int)
    error1_ranks[temp] = np.arange(len(all_error1))
    temp = all_error2.argsort()
    error2_ranks = np.empty(len(all_error2), int)
    error2_ranks[temp] = np.arange(len(all_error2))

    bp_joint_error_rank = error1_ranks+error2_ranks
    bp_ind_err = np.argmin(bp_joint_error_rank)

    bp1 = all_bp1[bp_ind_err]
    grad1 = all_grad1[bp_ind_err]
    grad2 = all_grad2[bp_ind_err]
    yinter = all_yinter[bp_ind_err]
    
    full_periods = np.log(periods[::-1])
    test = full_periods < bp1
    p1 = full_periods[test]
    test = full_periods >= bp1
    p2 = full_periods[test]    
    
    x_lens = [p1,p2]
    
    calc_1,calc_2 = spectra_f_parts_2(x_lens,grad1,grad2,yinter)

    calc_1 = calc_1.tolist()
    calc_2 = calc_2.tolist()
    
    calc_y = calc_1+calc_2 
    
    orig_periods_log = np.log(periods)[::-1]
    
    #calculate stats of fitted lines
    #-----------------------------------------
    ave1 = np.average(calc_1)
    ave2 = np.average(calc_2)
    med1 = np.median(calc_1)
    med2 = np.median(calc_2)
    sum1 = np.sum(calc_1)
    sum2 = np.sum(calc_2)
    
    y1_s = calc_1[0]
    y1_e = calc_1[-1]
    y2_s = calc_2[0]
    y2_e = calc_2[-1]

    return grad1,grad2,np.exp(bp1),np.exp(p1),np.exp(calc_1),np.exp(p2),np.exp(calc_2),np.exp(ave1),np.exp(ave2),np.exp(med1),np.exp(med2),np.exp(sum1),np.exp(sum2),np.exp(y1_s),np.exp(y1_e),np.exp(y2_s),np.exp(y2_e)

#3.25
#----------
    
#FIT FREE SEPARATE LINEAR REGRESSION FUNCTIONS TO SPECTRA (NON-PIECEWISE) (CURRENTLY USES 2 PIECES - REGIMES ARE WEATHER AND MACROWEATHER)

def spectra_fit_free_linear(periods,mag,ofac,lower_period_limit,upper_period_limit,lower_break,upper_break):
    
    start_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-upper_period_limit))
    end_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-lower_period_limit))
    
    if periods[end_cut_i] > lower_period_limit:
        end_cut_i = end_cut_i + 1.
    bp_periods = periods[start_cut_i:end_cut_i]
    bp_mag = mag[start_cut_i:end_cut_i]
    
    if (lower_period_limit < 1./12.) & (1./12. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./12.))
        n_points = int((600*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        if rm_points[-1] >= len(bp_periods):
            rm_points = range(rm_points[0],len(bp_periods)-1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./11.) & (1./11. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./11.))
        n_points = int((550*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./10.) & (1./10. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./10.))
        n_points = int((500*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./9.) & (1./9. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./9.))
        n_points = int((450*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./8.) & (1./8. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./8.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./7.) & (1./7. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./7.))
        n_points = int((350*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./6.) & (1./6. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./6.))
        n_points = int((300*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./5.) & (1./5. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./5.))
        n_points = int((250*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./4.) & (1./4. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./4.))
        n_points = int((200*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./3.) & (1./3. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./3.))
        n_points = int((150*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./2.) & (1./2. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./2.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 1.) & (1. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/4.) & (365.25/4. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/4.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 365.25/3.) & (365.25/3. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/3.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/2.) & (365.25/2. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/2.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
     
    if (lower_period_limit < 365.25) & (365.25 < upper_period_limit):
        annual_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(annual_i-n_points,(annual_i+n_points)+1)
        if rm_points[0] < 0:
            rm_points = range(0,(annual_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    #make periods (and mag) go in ascending order.
    bp_periods = bp_periods[::-1]
    bp_mag = bp_mag[::-1]
    
    #log periods and mag, will keep log distributed nature after interpolation
    bp_periods = np.log(bp_periods)
    bp_mag = np.log(bp_mag)
    
    #set limiting ranges for breakpoints
    breakpoint_test = (bp_periods >= np.log(lower_break)) & (bp_periods < np.log(upper_break))
    bp_periods_1 = bp_periods[breakpoint_test]
    bp1 = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-bp_periods_1[0]))

    all_error1 = []
    all_error2 = []
    all_bp1 = []
    all_grad1 = []
    all_grad2 = []
    all_yinter = []
    all_yinter2 = []

    for i in range(len(bp_periods_1)):
        bp1 = bp_periods_1[i]
    
        test = bp_periods < bp1
        bp_periods1 = bp_periods[test]
        bp_mag1 = bp_mag[test]
        len1 = len(bp_periods1)
        test = bp_periods >= bp1
        bp_periods2 = bp_periods[test]
        bp_mag2 = bp_mag[test]
        len2 = len(bp_periods2)
    
        grad1, yinter1, r_value, p_value, std_err = stats.linregress(bp_periods1,bp_mag1)
        grad2, yinter2, r_value, p_value, std_err = stats.linregress(bp_periods2,bp_mag2)
    
        line1 = grad1*bp_periods1+yinter1
        line2 = grad2*bp_periods2+yinter2
        
        err1=np.sqrt(np.sum((np.abs(line1-bp_mag1))**2)/len(bp_mag1))
        err2=np.sqrt(np.sum((np.abs(line2-bp_mag2))**2)/len(bp_mag2))
    
        all_error1.append(err1)
        all_error2.append(err2)
        all_bp1.append(np.average((bp_periods1[-1],bp_periods2[0])))
        all_grad1.append(grad1)
        all_grad2.append(grad2)
        all_yinter.append(yinter1)
        all_yinter2.append(yinter2)
        
    #get best fit based on best rank of errors for each line
    all_error1 = np.array(all_error1)
    all_error2 = np.array(all_error2)

    temp = all_error1.argsort()
    error1_ranks = np.empty(len(all_error1), int)
    error1_ranks[temp] = np.arange(len(all_error1))
    temp = all_error2.argsort()
    error2_ranks = np.empty(len(all_error2), int)
    error2_ranks[temp] = np.arange(len(all_error2))

    bp_joint_error_rank = error1_ranks+error2_ranks
    bp_ind_err = np.argmin(bp_joint_error_rank)

    bp1 = all_bp1[bp_ind_err]
    grad1 = all_grad1[bp_ind_err]
    grad2 = all_grad2[bp_ind_err]
    yinter = all_yinter[bp_ind_err]
    yinter2 = all_yinter2[bp_ind_err]
    
    full_periods = np.log(periods[::-1])
    test = full_periods < bp1
    p1 = full_periods[test]
    test = full_periods >= bp1
    p2 = full_periods[test]    
    
    calc_1 = grad1*p1+yinter1
    calc_2 = grad2*p2+yinter2
    
    calc_1 = calc_1.tolist()
    calc_2 = calc_2.tolist()
    
    calc_y = calc_1+calc_2 
    
    orig_periods_log = np.log(periods)[::-1]
    
    #calculate stats of fitted lines
    #-----------------------------------------
    ave1 = np.average(calc_1)
    ave2 = np.average(calc_2)
    med1 = np.median(calc_1)
    med2 = np.median(calc_2)
    sum1 = np.sum(calc_1)
    sum2 = np.sum(calc_2)
    
    y1_s = calc_1[0]
    y1_e = calc_1[-1]
    y2_s = calc_2[0]
    y2_e = calc_2[-1]

    return grad1,grad2,np.exp(bp1),np.exp(p1),np.exp(calc_1),np.exp(p2),np.exp(calc_2),np.exp(ave1),np.exp(ave2),np.exp(med1),np.exp(med2),np.exp(sum1),np.exp(sum2),np.exp(y1_s),np.exp(y1_e),np.exp(y2_s),np.exp(y2_e)

    
#-------------------------------------------------------------
#-------------------------------------------------------------
#SECTION 4 - STATISTICS

#4.01
#----------

#MOVING AVERAGE FUNCTION WHERE N IS SIZE OF WINDOW FOR AVERAGE
#TYPE CAN BE SIMPLE OR EXPONENTIAL

def moving_average(x, n, type='exponential'):
    x = np.asarray(x)
    if type=='simple':
        weights = np.ones(n)
    else:
        weights = np.exp(np.linspace(-1., 0., n))

    weights /= weights.sum()

    a =  np.convolve(x, weights, mode='full')[:len(x)]
    a[:n] = a[n]
    return a

#4.02
#----------
    
#EXPONENTIAL MOVING AVERAGE FUNCTION

def exp_smooth(s, periods, n):
	"""
	returns an n period exponential moving average for
	the time series s

	s is a list ordered from oldest (index 0) to most
	recent (index -1)
	n is an integer

	returns a numeric array of the exponential
	moving average
	"""
	s = np.array(s)
	ema = []
	
	cut_periods = []
	j = 1

    #get n sma first and calculate the next n period ema
	sma = sum(s[:n]) / n
	multiplier = 2 / float (1 + n)
	ema.append(sma)
	periods_mid = float(n/2)
	
	if np.mod(n,2) == 0:
		periods_mid = int(periods_mid)
		period_mid_val =  periods[periods_mid]
		first_period = periods_mid-1
		first_period_val =  periods[first_period]
		valid_period = np.average((first_period_val,period_mid_val))
		cut_periods.append(valid_period)
	else:
		valid_period = periods_mid-0.5
		cut_periods.append(periods[valid_period])

	#EMA(current) = ( (Price(current) - EMA(prev) ) x Multiplier) + EMA(prev)
	ema.append(( (s[n] - sma) * multiplier ) + sma)
	cut_periods.append(periods[n])

	adjust_f = (n-len(ema))+1

	#now calculate the rest of the values
	periods_counter = n+1
	for i in s[n+1:]:
		tmp = ( (i - ema[j]) * multiplier ) + ema[j]
		j = j + 1
		ema.append(tmp)
		cut_periods.append(periods[periods_counter])
		periods_counter +=1

	return cut_periods,ema

#4.03
#----------
    
#USE SAVITZY GOLAY FILTER TO REMOVE HIGH FREQUENCY NOISE - WORKS V WELL - A LITTLE BIT MAGIC...
	
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

#4.05
#----------
    
#TAKE LOG BINNED MEAN OF XY DATA (USUALLY SPECTRA) - SMOOTHS BEAUTIFULLY

def log_bin_mean(limit, n,x_array,y_array):
    #to get exact numbers areas of len n, must add 1 to n. 
    #eg 10 log spaced numbers returns 9 areas to smooth between. Must add 1 to get 10 areas
    n+=1
    
    #generate log spaced numbers
    result = [1]
    if n>1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result)<n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # problem! same integer. we need to find next_value by artificially incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining values will scale correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    log_spaced =  np.array(map(lambda x: round(x)-1, result),dtype=np.uint64)
    #amend last value to make it size of input array
    print log_spaced
    log_spaced = log_spaced.tolist()
    log_spaced[-1] = np.uint64(limit)
    
    #start smoothing of arrays
    smoothed_x_array = np.array([])
    smoothed_y_array = np.array([])
    
    x_array = x_array[::-1]
    y_array = y_array[::-1]
    
    for i in range(len(log_spaced)-1):
        try:
            start_i+=1
            end_i+=1
            print np.mean(x_array[log_spaced[start_i]:log_spaced[end_i]])
            smoothed_x_array = np.append(smoothed_x_array,np.mean(x_array[log_spaced[start_i]:log_spaced[end_i]]))
            smoothed_y_array = np.append(smoothed_y_array,np.mean(y_array[log_spaced[start_i]:log_spaced[end_i]]))
        except:
            start_i = 0
            end_i = 1   
            print np.mean(x_array[log_spaced[start_i]:log_spaced[end_i]])
            smoothed_x_array = np.append(smoothed_x_array,np.mean(x_array[log_spaced[start_i]:log_spaced[end_i]]))
            smoothed_y_array = np.append(smoothed_y_array,np.mean(y_array[log_spaced[start_i]:log_spaced[end_i]]))
            
    return smoothed_x_array,smoothed_y_array			

#4.06
#----------
    
#TAKE LOG BINNED MEDIAN OF XY DATA (USUALLY SPECTRA) - SMOOTHS BEAUTIFULLY

def log_bin_median(limit, n,x_array,y_array):
    #to get exact numbers areas of len n, must add 1 to n. 
    #eg 10 log spaced numbers returns 9 areas to smooth between. Must add 1 to get 10 areas
    n+=1
    
    #generate log spaced numbers
    result = [1]
    if n>1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result)<n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # problem! same integer. we need to find next_value by artificially incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining values will scale correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    log_spaced =  np.array(map(lambda x: round(x)-1, result),dtype=np.uint64)
    #amend last value to make it size of input array
    print log_spaced
    log_spaced = log_spaced.tolist()
    log_spaced[-1] = np.uint64(limit)
    
    #start smoothing of arrays
    smoothed_x_array = np.array([])
    smoothed_y_array = np.array([])
    
    for i in range(len(log_spaced)-1):
        try:
            start_i+=1
            end_i+=1
            smoothed_x_array = np.append(smoothed_x_array,np.median(x_array[log_spaced[start_i]:log_spaced[end_i]]))
            smoothed_y_array = np.append(smoothed_y_array,np.median(y_array[log_spaced[start_i]:log_spaced[end_i]]))
        except:
            start_i = 0
            end_i = 1    
            smoothed_x_array = np.append(smoothed_x_array,np.median(x_array[log_spaced[start_i]:log_spaced[end_i]]))
            smoothed_y_array = np.append(smoothed_y_array,np.median(y_array[log_spaced[start_i]:log_spaced[end_i]]))
            
    return smoothed_x_array,smoothed_y_array

#4.07
#----------
    
#REMOVE SIDELOBES ON SPECTRA OK KEY PEAKS
    
def sidelobe_peak_remove(fb,fr,fi,closest,crit_percent,periods):
    critical_val = (fb[closest]/100.) * crit_percent # certain % of top peak
    up_peak_val = fb[closest]
    down_peak_val = fb[closest]
    i=1
    while up_peak_val >= critical_val:
        if i != 1:
            try:
                up_inds.append(closest+(i-1))
            except:
                up_inds = [closest+(i-1)]
    
        up_peak_val = fb[closest+i]
        i+=1 
    if i == 2:
        up_inds = np.array([])
    i=1
    while down_peak_val >= critical_val:
        if i != 1:
            try:
                down_inds.append(closest-(i-1))
            except:
                down_inds = [closest-(i-1)]
        down_peak_val = fb[closest-i]
        i+=1

    if i == 2:
        down_inds = []
    all_inds = np.concatenate((down_inds,up_inds))
    all_inds = np.sort(all_inds)
    all_inds = [int(i) for i in all_inds]
    
    #set significant indices for mag, fr, and fi to 0
    fb[closest] = 0
    fr[closest] = 0
    fi[closest] = 0
    fb[all_inds] = 0
    fr[all_inds] = 0
    fi[all_inds] = 0
    
    print 'Altering periods: ', periods[closest], periods[all_inds] 
    
    return fb,fr,fi
 
 #4.08
#----------
 
#DEFINE SINC WINDOW FUNCTION
    
def sinc(xs):
    pxs = np.pi*xs
    return np.where(np.fabs(pxs)<1e-3, 1.0-pxs*pxs/6.0, np.sin(pxs)/pxs)

#4.09
#----------
    
#DEFINE HANNING WINDOW FUNCTION

def hanning_window(xs, halfwidth):
    win =  0.5 + 0.5*np.cos(np.pi*xs/halfwidth)
    return np.where(np.fabs(xs)<=halfwidth, win, 0.0)

#4.10
#----------
    
#CALCULATE GREAT CIRCLE DISTANCE BETWEEN 2 POINTS ON EARTH
   
def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2]) 
    # haversine formula 
    dlon = lon2 - lon1  
    dlat = lat2 - lat1  
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    return km

#4.11
#----------
    
#GET ORTHOGONAL REGRESSION RESIDUAL VARIANCE
   
def orthogonal_1_1_res_var(x_array,y_array): # x3,y3 is the point
    
    x = np.arange(0,1000,1)
    y = np.arange(0,1000,1)
    
    x1 = x[0] 
    y1 = y[0]
    x2 = x[-1]
    y2 = y[-1]
    
    all_dists = []
    
    for i in range(len(x_array)):
        x3 = x_array[i] 
        y3 = y_array[i]
    
        px = x2-x1
        py = y2-y1

        something = px*px + py*py

        u =  ((x3 - x1) * px + (y3 - y1) * py) / float(something)

        if u > 1:
            u = 1
        elif u < 0:
            u = 0

        x = x1 + u * px
        y = y1 + u * py

        dx = x - x3
        dy = y - y3

        # Note: If the actual distance does not matter,
        # if you only want to compare what this function
        # returns to other results of this function, you
        # can just return the squared distance instead
        # (i.e. remove the sqrt) to gain a little performance

        dist = math.sqrt(dx*dx + dy*dy)
        all_dists = np.append(all_dists,dist)
    
    dist_sum = np.sum(all_dists)
    dist_ave = np.average(all_dists)
        
    return np.around(dist_sum,2),np.around(dist_ave,2)
 
#4.12
#----------
 
#DETECT PEAKS ON SPECTRA BASED ON AMPLITUDES
   
def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):

    """Detect peaks in data based on their amplitude and other features.

    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`
    
    The function can handle NaN's 

    See this IPython Notebook [1]_.

    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    Examples
    --------
    >>> from detect_peaks import detect_peaks
    >>> x = np.random.randn(100)
    >>> x[60:81] = np.nan
    >>> # detect all peaks and plot data
    >>> ind = detect_peaks(x, show=True)
    >>> print(ind)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # set minimum peak height = 0 and minimum peak distance = 20
    >>> detect_peaks(x, mph=0, mpd=20, show=True)

    >>> x = [0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0]
    >>> # set minimum peak distance = 2
    >>> detect_peaks(x, mpd=2, show=True)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # detection of valleys instead of peaks
    >>> detect_peaks(x, mph=0, mpd=20, valley=True, show=True)

    >>> x = [0, 1, 1, 0, 1, 1, 0]
    >>> # detect both edges
    >>> detect_peaks(x, edge='both', show=True)

    >>> x = [-2, 1, -2, 2, 1, 1, 3, 0]
    >>> # set threshold = 2
    >>> detect_peaks(x, threshold = 2, show=True)
    """

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])

    if show:
        if indnan.size:
            x[indnan] = np.nan
        if valley:
            x = -x
        _plot(x, mph, mpd, threshold, edge, valley, ax, ind)

    return ind

#4.13
#----------
    
#GET DEGRESS OF FREEDOM OF WINDOW USED IN SPECTRAL ANALYSIS
   
def getdof(iwin):
    #c50 relates to type of window used, 0 is rectangular, 2 is hanning
    c50 = [0.500, 0.344, 0.167, 0.250, 0.096]

    c2 = 2.0 * c50[iwin] * c50[iwin]
    denom = 1.0 + c2 - c2/ 1.
    neff = 1. / denom
    dof = 2.0 * neff
    return dof

#4.14
#----------
    
#GET CHI SQUARED USING DEGREES OF FREEDOM

def getchi2(dof,alpha):
    tol = 1.0e-3
    itmax = 100
    ierr = 0
    
#use approximation for dof > 30 (Eq. 1.132 in Sachs (1984))
#----------------------------------------------------------
    if dof > 30.0:
        #za = -getz(alpha)   # NB: Eq. requires change of sign for percentile
        if ierr == 1:
            return
        x = 2.0 / 9.0 / dof
        chi2 = dof * (1.0 - x + za * np.sqrt(x))**3.0
    else:
        iter = 0
        lm = 0.0
        rm = 1000.0
        if alpha > 0.5:
            eps = (1.0 - alpha) * tol
        else:
            eps = alpha * tol
        while iter <= itmax:
            iter= iter + 1
            chi2 = 0.5 * (lm + rm)
            ac = 1.0 - scipy.special.gammainc(0.5*dof, 0.5*chi2)
            if np.abs(ac - alpha) <= eps:
                break
            if ac > alpha:
                lm = chi2
            else:
                rm = chi2

    if iter > itmax:
        print "Error in GETCHI2: Iter > ItMax"
        ierr = 1
        return
    
    getchi2 = chi2

    return getchi2
    
#4.15
#----------
    
#GET CHI SQUARED SIGNIFICANCE VALUES FOR SPECTRA AT VARIOUS DIFFERENT PERCENTILES OF CONFIDENCE

def chi_signif(periods,window_t_or_f):
    if window_t_or_f == 't':
        num = 2
    else:
        num = 0
    
    #get degrees of freedom
    dof = getdof(num) # 0 is rectangular,2 is hanning window

    #get scaling factors for red noise model
    fac80 = getchi2(dof, 0.20) / dof
    fac85 = getchi2(dof, 0.15) / dof
    fac90 = getchi2(dof, 0.10) / dof
    fac95 = getchi2(dof, 0.05) / dof
    fac99 = getchi2(dof, 0.01) / dof
    fac99_9 = getchi2(dof, 0.001) / dof
    fac99_99 = getchi2(dof, 0.0001) / dof
    
    # critical false alarm level after Thomson (1990)
    # -----------------------------------------------
    ave_seg_len = len(periods)
    alphacrit = 1.0 / ave_seg_len
    faccrit = getchi2(dof, alphacrit) / dof
    
    return fac80,fac85,fac90,fac95,fac99,fac99_9,fac99_99,faccrit

#4.16
#----------
    
#LIST MUST BE SORTED

def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx

#4.17
#----------
    
#GET AVERAGE AND STANDARD DEVIATION OF PHASE 

def yamartino_ave_std(data,max_ph):
    #remove nans
    data = data[~np.isnan(data)]
    
    if len(data) > 0:
        #convert phase into radians
        pi2 = np.pi*2.
        ratio = pi2/max_ph
        rads = ratio*np.array(data)
    
        sum_sin = 0
        sum_cos = 0
    
        ave = []
        
        for i in range(len(rads)):
            sum_sin = sum_sin + np.sin(rads[i])        
            sum_cos = sum_cos + np.cos(rads[i])
    
        ave = np.arctan2(sum_sin,sum_cos)
    
        #convert radians to 0 to 2pi (from -np.pi:pi)
        if ave < 0:
            ave = np.pi+ (np.pi-np.abs(ave))

        #convert radians back to phase
        ratio = max_ph/(np.pi*2.)
        mean = ratio*ave
    
        #calculate stdev
        sum_cos = sum_cos/len(data)
        sum_sin = sum_sin/len(data)
    
        e = np.sqrt(1. - sum_sin * sum_sin - sum_cos * sum_cos)
        std = np.arcsin(e) * (1 + 0.1547 * e * e * e)
        std = max_ph/(pi2/std)

        return mean,std
    else:
        return np.NaN,np.NaN
        
def difference_matrix(a):
    x = np.reshape(a, (len(a), 1))
    y = x - x.transpose()
    z = np.ravel([np.delete(y[i],i) for i in range(len(y))])
    return z

def flatten(container):
    for i in container:
        if isinstance(i, list) or isinstance(i, tuple):
            for j in flatten(i):
                yield j
        else:
            yield i

#-------------------------------------------------------------
#-------------------------------------------------------------
#SECTION 5 - INTERACTIVE PLOTTING

#5.01
#----------

#INTERACTIVE CLICKER FOR MAP MULTI OBS/MODEL PLOT
#PLOTS TIMESERIES AT SITE, DIURNAL AND SEASONAL CYCLES, AND JOINT PERIODIC CYCLE FOR OBS AND MODEL
def clicker_interactive_map_obsmodel(event,species,lat_e,lon_e,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_d_waveform,obs_s_waveform,obs_all_waveform,model_d_waveform,model_s_waveform,model_all_waveform,fig,all_m):
    
    #ind of point clicked near
    ind = event.ind
    print event
    print ind
    print obs_refs[ind]
    print len(obs_refs)
    print obs_refs[:20]
    ind = ind[0]
    ref = obs_refs[ind]
    
    
    #read in obs and model periodic parts
    obs_site_group = obs_period_grp.groups[ref]
    obs_d_mag = obs_site_group.daily_amplitude
    obs_d_ph = obs_site_group.daily_phase
    obs_s_mag = obs_site_group.seasonal_amplitude
    obs_s_ph = obs_site_group.seasonal_phase
    obs_ave = obs_site_group.average
    obs_d_waveform = obs_d_waveform[ind]
    obs_s_waveform = obs_s_waveform[ind]
    obs_all_waveform = obs_all_waveform[ind]
    
    model_site_group = model_period_grp.groups[ref]
    model_d_mag = model_site_group.daily_amplitude
    model_d_ph = model_site_group.daily_phase
    model_s_mag = model_site_group.seasonal_amplitude
    model_s_ph = model_site_group.seasonal_phase
    model_ave = model_site_group.average
    model_d_waveform = model_d_waveform[ind]
    model_s_waveform = model_s_waveform[ind]
    model_all_waveform = model_all_waveform[ind]
      
    #get obs timeseries for site clicked and obs details
    obs_site_group = obs_ts_grp.groups[ref]
    obs_var = obs_site_group.variables[species.lower()][:]
    group = obs_site_group.process_group
    lat = obs_site_group.latitude
    lon = obs_site_group.longitude
    lon = obs_site_group.longitude
    alt = obs_site_group.altitude
    complete = obs_site_group.data_completeness
    a_class = obs_site_group.anthrome_class
    r_class = obs_site_group.raw_class
    continent = loc_dict[tags[ind]]
    country = obs_site_group.country
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
    obs_var_pd = pd.Series(obs_var_mask, index=obs_time_pd)
    
    print lat,lon
    
    #highlight points on click
    pl = all_m[0].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20)
    pl2 = all_m[1].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20)
    pl3 = all_m[2].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20) 
    pl4 = all_m[3].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20) 
    
    #get model timeseries for site clicked
    model_var = model_ts_grp.variables[species.lower()][:]*1e9
    lat_n,lon_n = obs_model_gridbox(lat_e,lon_e,lat,lon)
    model_var_pick = model_var[:,lat_n,lon_n]
    model_var_mask = np.ma.masked_where(model_var_pick<=0,model_var_pick)
    model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
    model_var_pd = pd.Series(model_var_mask, index=model_time_pd)     
    
    #get cut times and datetimes
    d_datetime = obs_datetimes[:24]
    s_datetime = obs_datetimes[:8766]
    d_time_pd = pd.date_range(start = d_datetime[0],end = d_datetime[-1], freq = 'H')
    s_time_pd = pd.date_range(start = s_datetime[0],end = s_datetime[-1], freq = 'H')
    
    #make panda series
    obs_d_waveform_pd = pd.Series(obs_d_waveform, index=d_time_pd)
    obs_s_waveform_pd = pd.Series(obs_s_waveform, index=s_time_pd)
    obs_all_waveform_pd = pd.Series(obs_all_waveform, index=obs_time_pd)
    
    model_d_waveform_pd = pd.Series(model_d_waveform, index=d_time_pd)
    model_s_waveform_pd = pd.Series(model_s_waveform, index=s_time_pd)
    model_all_waveform_pd = pd.Series(model_all_waveform, index=model_time_pd)
    
    fig.canvas.draw()
        
    fig2 = plt.figure(figsize =(24,12))
    fig2.patch.set_facecolor('white')
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=1)
    ax3 = plt.subplot2grid((2,2), (1,1), rowspan=1)

    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_var_pd, color='black',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Raw Data')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='red',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Raw Data')
    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_all_waveform_pd, color='gray',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Significant Waveform')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_all_waveform_pd, color='lightcoral',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Significant Waveform')
    ax2.plot_date(d_time_pd.to_pydatetime(), obs_d_waveform_pd, color='black',marker = 'o', markersize = 6,markeredgecolor = None, mew = 0, linestyle = '--', label = 'Obs. Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), obs_s_waveform_pd, color='black',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Obs. Seasonal Waveform') 
    ax2.plot_date(d_time_pd.to_pydatetime(), model_d_waveform_pd, color='red',marker = 'o', markersize = 6, markeredgecolor = None, mew = 0,linestyle = '--', label = 'Model Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), model_s_waveform_pd, color='red',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Model Seasonal Waveform')   
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=model_d_waveform_pd>obs_d_waveform_pd, facecolor='yellow', interpolate=True)
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=obs_d_waveform_pd>model_d_waveform_pd, facecolor='blue', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=model_s_waveform_pd>obs_s_waveform_pd, facecolor='yellow', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=obs_s_waveform_pd>model_s_waveform_pd, facecolor='blue', interpolate=True)
    
    ax1.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s, Grid Index = %s,%s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,lat_n,lon_n))
    
    #plt.legend(loc = 'lower right') 
    plt.tight_layout()
    ax1.grid()
    ax2.grid()
    ax3.grid()
    
    pl.pop(0).remove()
    pl2.pop(0).remove()
    pl3.pop(0).remove()
    pl4.pop(0).remove()
    
    plt.show()
    
    return

#5.02
#----------
    
#INTERACTIVE CLICKER FOR MAP SINGLE MODEL PLOT
#PLOTS TIMESERIES AT SITE, DIURNAL AND SEASONAL CYCLES, AND JOINT PERIODIC CYCLE FOR OBS AND MODEL
def clicker_interactive_map_model(event,plot_type,lat_e,lon_e,linear_lats,linear_lons,datetimes,all_periods,var,d_waveform,s_waveform,all_waveform,fig,m):
    
    #ind of point clicked near
    ind = event.ind
    ind = ind[0]
    
    #highlight points on click
    pl = m.plot([linear_lons[ind]], [linear_lats[ind]], 's', ms=20, alpha=0.6, color='yellow',zorder=20)
    if plot_type == 'multi':
        pl2 = 1
        pl3 = 1
        pl4 = 1
    
    #get model timeseries for site clicked
    lat_n,lon_n = obs_model_gridbox(lat_e,lon_e,linear_lats[ind],linear_lons[ind])
    var_pick = var[:,lat_n,lon_n]   
    var_mask = np.ma.masked_where(var_pick<=0,var_pick)
    time_pd = pd.date_range(start = datetimes[0],end = datetimes[-1], freq = 'H')
    var_pd = pd.Series(var_mask, index=time_pd) 
    
    #read in periodic parts
    d_mag = all_periods[4,lat_n,lon_n]
    d_ph = all_periods[9,lat_n,lon_n]
    s_mag = all_periods[14,lat_n,lon_n]
    s_ph = all_periods[19,lat_n,lon_n]
    ave = all_periods[20,lat_n,lon_n]
    d_waveform = d_waveform[ind]
    s_waveform = s_waveform[ind]
    all_waveform = all_waveform[ind]
    
    #get cut times and datetimes
    d_datetime = datetimes[:24]
    s_datetime = datetimes[:8766]
    d_time_pd = pd.date_range(start = d_datetime[0],end = d_datetime[-1], freq = 'H')
    s_time_pd = pd.date_range(start = s_datetime[0],end = s_datetime[-1], freq = 'H')
    
    #make panda series
    d_waveform_pd = pd.Series(d_waveform, index=d_time_pd)
    s_waveform_pd = pd.Series(s_waveform, index=s_time_pd)
    all_waveform_pd = pd.Series(all_waveform, index=time_pd)
    
    fig.canvas.draw()
        
    fig2 = plt.figure(figsize =(24,12))
    fig2.patch.set_facecolor('white')
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=1)
    ax3 = plt.subplot2grid((2,2), (1,1), rowspan=1)

    ax1.plot_date(time_pd.to_pydatetime(), var_pd, color='red', markersize = 3,marker = 'o', markeredgecolor = None, mew = 0, label = 'Model Raw Data')
    ax1.plot_date(time_pd.to_pydatetime(), all_waveform_pd, color='lightcoral', markersize = 3, marker = 'o', markeredgecolor = None, mew = 0, label = 'Model Significant Waveform')
    ax2.plot_date(d_time_pd.to_pydatetime(), d_waveform_pd, color='red', markersize = 6, marker = 'o', markeredgecolor = None, mew = 0, label = 'Model Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), s_waveform_pd, color='red', markersize = 3, marker = 'o', markeredgecolor = None, mew = 0, label = 'Model Seasonal Waveform') 
    
    #plt.legend(loc = 'lower right') 
    plt.tight_layout()
    ax1.grid()
    ax2.grid()
    ax3.grid()
    
    pl.pop(0).remove()
    
    plt.show()
    
    return

#5.03
#----------
    
#INTERACTIVE CLICKER FOR XY SINGLE OBS/MODEL PLOT
#PLOTS TIMESERIES AT SITE, DIURNAL AND SEASONAL CYCLES, AND JOINT PERIODIC CYCLE FOR OBS AND MODEL

def clicker_interactive_xy_obsmodel_single(event,species,lat_e,lon_e,obs_lats,obs_lons,date,time,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_chosen,model_chosen,obs_d_waveform,obs_s_waveform,obs_all_waveform,model_d_waveform,model_s_waveform,model_all_waveform,fig,ax):
    
    #ind of point clicked near
    ind = event.ind
    ind = ind[0]
    ref = obs_refs[ind]
    
    #read in obs and model periodic parts
    obs_site_group = obs_period_grp.groups[ref]
    obs_d_mag = obs_site_group.daily_amplitude
    obs_d_ph = obs_site_group.daily_phase
    obs_s_mag = obs_site_group.seasonal_amplitude
    obs_s_ph = obs_site_group.seasonal_phase
    obs_ave = obs_site_group.average
    obs_d_waveform = obs_d_waveform[ind]
    obs_s_waveform = obs_s_waveform[ind]
    obs_all_waveform = obs_all_waveform[ind]
    obs_d_ff = obs_site_group.daily_ff
    obs_s_ff = obs_site_group.seasonal_ff
    
    model_site_group = model_period_grp.groups[ref]
    model_d_mag = model_site_group.daily_amplitude
    model_d_ph = model_site_group.daily_phase
    model_s_mag = model_site_group.seasonal_amplitude
    model_s_ph = model_site_group.seasonal_phase
    model_ave = model_site_group.average
    model_d_waveform = model_d_waveform[ind]
    model_s_waveform = model_s_waveform[ind]
    model_all_waveform = model_all_waveform[ind]
    model_d_ff = model_site_group.daily_ff
    model_s_ff = model_site_group.seasonal_ff
    
    #highlight points on click
    pl = ax.plot(obs_chosen[ind], model_chosen[ind], 'o', ms=20, alpha=0.6, color='yellow',zorder=20)        
    
    #get obs timeseries for site clicked and obs details
    obs_site_group = obs_ts_grp.groups[ref]
    obs_var = obs_site_group.variables[species.lower()][:]
    group = obs_site_group.process_group
    lat = obs_site_group.latitude
    lon = obs_site_group.longitude
    lon = obs_site_group.longitude
    alt = obs_site_group.altitude
    complete = obs_site_group.data_completeness
    a_class = obs_site_group.anthrome_class
    r_class = obs_site_group.raw_class
    continent = loc_dict[tags[ind]]
    country = obs_site_group.country
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
    obs_var_pd = pd.Series(obs_var_mask, index=obs_time_pd)
    
    #get model timeseries for site clicked
    model_var = model_ts_grp.variables[species.lower()][:]*1e9
    lat_n,lon_n = obs_model_gridbox(lat_e,lon_e,obs_lats[ind],obs_lons[ind])
    model_var_pick = model_var[:,lat_n,lon_n]
    model_var_mask = np.ma.masked_where(model_var_pick<=0,model_var_pick)
    model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
    model_var_pd = pd.Series(model_var_mask, index=model_time_pd)     
    
    #get cut times and datetimes
    d_datetime = obs_datetimes[:24]
    s_datetime = obs_datetimes[:8766]
    d_time_pd = pd.date_range(start = d_datetime[0],end = d_datetime[-1], freq = 'H')
    s_time_pd = pd.date_range(start = s_datetime[0],end = s_datetime[-1], freq = 'H')
    
    #make panda series
    obs_d_waveform_pd = pd.Series(obs_d_waveform, index=d_time_pd)
    obs_s_waveform_pd = pd.Series(obs_s_waveform, index=s_time_pd)
    obs_all_waveform_pd = pd.Series(obs_all_waveform, index=obs_time_pd)
    
    model_d_waveform_pd = pd.Series(model_d_waveform, index=d_time_pd)
    model_s_waveform_pd = pd.Series(model_s_waveform, index=s_time_pd)
    model_all_waveform_pd = pd.Series(model_all_waveform, index=model_time_pd)
    
    fig.canvas.draw()
        
    fig2 = plt.figure(figsize =(24,12))
    fig2.patch.set_facecolor('white')
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=1)
    ax3 = plt.subplot2grid((2,2), (1,1), rowspan=1)

    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_var_pd, color='black',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Raw Data')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='red',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Raw Data')
    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_all_waveform_pd, color='gray',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Significant Waveform')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_all_waveform_pd, color='lightcoral',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Significant Waveform')
    ax2.plot_date(d_time_pd.to_pydatetime(), obs_d_waveform_pd, color='black',marker = 'o', markersize = 6,markeredgecolor = None, mew = 0, linestyle = '--', label = 'Obs. Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), obs_s_waveform_pd, color='black',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Obs. Seasonal Waveform') 
    ax2.plot_date(d_time_pd.to_pydatetime(), model_d_waveform_pd, color='red',marker = 'o', markersize = 6, markeredgecolor = None, mew = 0,linestyle = '--', label = 'Model Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), model_s_waveform_pd, color='red',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Model Seasonal Waveform')   
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=model_d_waveform_pd>obs_d_waveform_pd, facecolor='yellow', interpolate=True)
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=obs_d_waveform_pd>model_d_waveform_pd, facecolor='blue', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=model_s_waveform_pd>obs_s_waveform_pd, facecolor='yellow', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=obs_s_waveform_pd>model_s_waveform_pd, facecolor='blue', interpolate=True)
    
    ax1.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s, Grid Index = %s,%s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,lat_n,lon_n))
    
    #plt.legend(loc = 'lower right') 
    plt.tight_layout()
    ax1.grid()
    ax2.grid()
    ax3.grid()
    
    pl.pop(0).remove()
    
    plt.show()
    
    return

#5.04
#----------
    
#INTERACTIVE CLICKER FOR MAP MULTI OBS/MODEL PLOT
#PLOTS TIMESERIES AT SITE, DIURNAL AND SEASONAL CYCLES, AND JOINT PERIODIC CYCLE FOR OBS AND MODEL

def clicker_interactive_xy_obsmodel_multi(event,species,lat_e,lon_e,obs_lats,obs_lons,date,time,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_d_waveform,obs_s_waveform,obs_all_waveform,model_d_waveform,model_s_waveform,model_all_waveform,fig,ax,ax2,ax3,ax4,ax5,ax7,ax8):
    
    #ind of point clicked near
    ind = event.ind
    ind = ind[0]
    ref = obs_refs[ind]
    
    #read in obs and model periodic parts
    obs_site_group = obs_period_grp.groups[ref]
    obs_d_mag = obs_site_group.daily_amplitude
    obs_d_ph = obs_site_group.daily_phase
    obs_s_mag = obs_site_group.seasonal_amplitude
    obs_s_ph = obs_site_group.seasonal_phase
    obs_ave = obs_site_group.average
    obs_d_ff = obs_site_group.daily_ff
    obs_s_ff = obs_site_group.seasonal_ff
    obs_d_waveform = obs_d_waveform[ind]
    obs_s_waveform = obs_s_waveform[ind]
    obs_all_waveform = obs_all_waveform[ind]
    
    model_site_group = model_period_grp.groups[ref]
    model_d_mag = model_site_group.daily_amplitude
    model_d_ph = model_site_group.daily_phase
    model_s_mag = model_site_group.seasonal_amplitude
    model_s_ph = model_site_group.seasonal_phase
    model_ave = model_site_group.average
    model_d_ff = model_site_group.daily_ff
    model_s_ff = model_site_group.seasonal_ff
    model_d_waveform = model_d_waveform[ind]
    model_s_waveform = model_s_waveform[ind]
    model_all_waveform = model_all_waveform[ind]
    
    
    #highlight points on click
    pl = ax.plot(obs_d_mag, model_d_mag, 'o', ms=12, alpha=0.9, color='black',zorder=20)        
    pl2 = ax2.plot(obs_s_mag, model_s_mag, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    pl3 = ax3.plot(obs_ave, model_ave, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    pl4 = ax4.plot(obs_d_ph, model_d_ph, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    pl5 = ax5.plot(obs_s_ph, model_s_ph, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    pl6 = ax7.plot(obs_d_ff, model_d_ff, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    pl7 = ax8.plot(obs_s_ff, model_s_ff, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    
    #get obs timeseries for site clicked and obs details
    obs_site_group = obs_ts_grp.groups[ref]
    obs_var = obs_site_group.variables[species.lower()][:]
    group = obs_site_group.process_group
    lat = obs_site_group.latitude
    lon = obs_site_group.longitude
    lon = obs_site_group.longitude
    alt = obs_site_group.altitude
    complete = obs_site_group.data_completeness
    a_class = obs_site_group.anthrome_class
    r_class = obs_site_group.raw_class
    continent = loc_dict[tags[ind]]
    country = obs_site_group.country
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
    obs_var_pd = pd.Series(obs_var_mask, index=obs_time_pd)
    
    #get model timeseries for site clicked
    model_var = model_ts_grp.variables[species.lower()][:]*1e9
    lat_n,lon_n = obs_model_gridbox(lat_e,lon_e,obs_lats[ind],obs_lons[ind])
    model_var_pick = model_var[:,lat_n,lon_n]
    model_var_mask = np.ma.masked_where(model_var_pick<=0,model_var_pick)
    model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
    model_var_pd = pd.Series(model_var_mask, index=model_time_pd) 
    
    #get cut times and datetimes
    d_datetime = obs_datetimes[:24]
    s_datetime = obs_datetimes[:8766]
    d_time_pd = pd.date_range(start = d_datetime[0],end = d_datetime[-1], freq = 'H')
    s_time_pd = pd.date_range(start = s_datetime[0],end = s_datetime[-1], freq = 'H')
    
    #make panda series
    obs_d_waveform_pd = pd.Series(obs_d_waveform, index=d_time_pd)
    obs_s_waveform_pd = pd.Series(obs_s_waveform, index=s_time_pd)
    obs_all_waveform_pd = pd.Series(obs_all_waveform, index=obs_time_pd)
    
    model_d_waveform_pd = pd.Series(model_d_waveform, index=d_time_pd)
    model_s_waveform_pd = pd.Series(model_s_waveform, index=s_time_pd)
    model_all_waveform_pd = pd.Series(model_all_waveform, index=model_time_pd)
    
    fig.canvas.draw()
        
    fig2 = plt.figure(figsize =(24,12))
    fig2.patch.set_facecolor('white')
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=1)
    ax3 = plt.subplot2grid((2,2), (1,1), rowspan=1)

    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_var_pd, color='black',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Raw Data')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='red',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Raw Data')
    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_all_waveform_pd, color='gray',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Significant Waveform')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_all_waveform_pd, color='lightcoral',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Significant Waveform')
    ax2.plot_date(d_time_pd.to_pydatetime(), obs_d_waveform_pd, color='black',marker = 'o', markersize = 6,markeredgecolor = None, mew = 0, label = 'Obs. Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), obs_s_waveform_pd, color='black',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Obs. Seasonal Waveform') 
    ax2.plot_date(d_time_pd.to_pydatetime(), model_d_waveform_pd, color='red',marker = 'o', markersize = 6, markeredgecolor = None, mew = 0, label = 'Model Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), model_s_waveform_pd, color='red',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Model Seasonal Waveform')   
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=model_d_waveform_pd>obs_d_waveform_pd, facecolor='yellow', interpolate=True)
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=obs_d_waveform_pd>model_d_waveform_pd, facecolor='blue', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=model_s_waveform_pd>obs_s_waveform_pd, facecolor='yellow', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=obs_s_waveform_pd>model_s_waveform_pd, facecolor='blue', interpolate=True)
    
    ax1.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s, Grid Index = %s,%s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,lat_n,lon_n))
    
    #plt.legend(loc = 'lower right') 
    plt.tight_layout()
    ax1.grid()
    ax2.grid()
    ax3.grid()
    
    pl.pop(0).remove()
    pl2.pop(0).remove()
    pl3.pop(0).remove()
    pl4.pop(0).remove()
    pl5.pop(0).remove()
    
    plt.show()
    
    return

#5.05
#----------
    
#PLOTS TIMESERIES AT SITE, DIURNAL AND SEASONAL CYCLES, AND JOINT PERIODIC CYCLE FOR OBS AND MODEL

def clicker_interactive_map_obsmodel_loc(event,species,lat_e,lon_e,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_d_waveform,obs_s_waveform,obs_all_waveform,model_d_waveform,model_s_waveform,model_all_waveform,fig,all_m):
    
    #ind of point clicked near
    ind = event.ind
    print event
    print ind
    print obs_refs[ind]
    print len(obs_refs)
    print obs_refs[:20]
    ind = ind[0]
    ref = obs_refs[ind]
    
    
    #read in obs and model periodic parts
    obs_site_group = obs_period_grp.groups[ref]
    obs_d_mag = obs_site_group.daily_amplitude
    obs_d_ph = obs_site_group.daily_phase
    obs_s_mag = obs_site_group.seasonal_amplitude
    obs_s_ph = obs_site_group.seasonal_phase
    obs_ave = obs_site_group.average
    obs_d_waveform = obs_d_waveform[ind]
    obs_s_waveform = obs_s_waveform[ind]
    obs_all_waveform = obs_all_waveform[ind]
    
    model_site_group = model_period_grp.groups[ref]
    model_d_mag = model_site_group.daily_amplitude
    model_d_ph = model_site_group.daily_phase
    model_s_mag = model_site_group.seasonal_amplitude
    model_s_ph = model_site_group.seasonal_phase
    model_ave = model_site_group.average
    model_d_waveform = model_d_waveform[ind]
    model_s_waveform = model_s_waveform[ind]
    model_all_waveform = model_all_waveform[ind]
      
    #get obs timeseries for site clicked and obs details
    obs_site_group = obs_ts_grp.groups[ref]
    obs_var = obs_site_group.variables[species.lower()][:]
    group = obs_site_group.process_group
    lat = obs_site_group.latitude
    lon = obs_site_group.longitude
    lon = obs_site_group.longitude
    alt = obs_site_group.altitude
    complete = obs_site_group.data_completeness
    a_class = obs_site_group.anthrome_class
    r_class = obs_site_group.raw_class
    continent = loc_dict[tags[ind]]
    country = obs_site_group.country
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
    obs_var_pd = pd.Series(obs_var_mask, index=obs_time_pd)
    
    print lat,lon
    
    #highlight points on click
    pl = all_m[0].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20)
    pl2 = all_m[1].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20)
    pl3 = all_m[2].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20) 
    pl4 = all_m[3].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20) 
    
    #get model timeseries for site clicked
    model_var = model_ts_grp.variables[species.lower()][:]*1e9
    lat_n,lon_n = obs_model_gridbox(lat_e,lon_e,lat,lon)
    model_var_pick = model_var[:,lat_n,lon_n]
    model_var_mask = np.ma.masked_where(model_var_pick<=0,model_var_pick)
    model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
    model_var_pd = pd.Series(model_var_mask, index=model_time_pd)     
    
    #get cut times and datetimes
    d_datetime = obs_datetimes[:24]
    s_datetime = obs_datetimes[:8766]
    d_time_pd = pd.date_range(start = d_datetime[0],end = d_datetime[-1], freq = 'H')
    s_time_pd = pd.date_range(start = s_datetime[0],end = s_datetime[-1], freq = 'H')
    
    #make panda series
    obs_d_waveform_pd = pd.Series(obs_d_waveform, index=d_time_pd)
    obs_s_waveform_pd = pd.Series(obs_s_waveform, index=s_time_pd)
    obs_all_waveform_pd = pd.Series(obs_all_waveform, index=obs_time_pd)
    
    model_d_waveform_pd = pd.Series(model_d_waveform, index=d_time_pd)
    model_s_waveform_pd = pd.Series(model_s_waveform, index=s_time_pd)
    model_all_waveform_pd = pd.Series(model_all_waveform, index=model_time_pd)
    
    fig.canvas.draw()
        
    fig2 = plt.figure(figsize =(24,12))
    fig2.patch.set_facecolor('white')
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=1)
    ax3 = plt.subplot2grid((2,2), (1,1), rowspan=1)

    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_var_pd, color='black',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Raw Data')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='red',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Raw Data')
    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_all_waveform_pd, color='gray',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Significant Waveform')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_all_waveform_pd, color='lightcoral',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Significant Waveform')
    ax2.plot_date(d_time_pd.to_pydatetime(), obs_d_waveform_pd, color='black',marker = 'o', markersize = 6,markeredgecolor = None, mew = 0, linestyle = '--', label = 'Obs. Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), obs_s_waveform_pd, color='black',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Obs. Seasonal Waveform') 
    ax2.plot_date(d_time_pd.to_pydatetime(), model_d_waveform_pd, color='red',marker = 'o', markersize = 6, markeredgecolor = None, mew = 0,linestyle = '--', label = 'Model Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), model_s_waveform_pd, color='red',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Model Seasonal Waveform')   
   
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=obs_d_waveform_pd>model_d_waveform_pd, facecolor='blue', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=model_s_waveform_pd>obs_s_waveform_pd, facecolor='yellow', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=obs_s_waveform_pd>model_s_waveform_pd, facecolor='blue', interpolate=True)
    
    ax1.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s, Grid Index = %s,%s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,lat_n,lon_n))
    
    #plt.legend(loc = 'lower right') 
    plt.tight_layout()
    ax1.grid()
    ax2.grid()
    ax3.grid()
    
    pl.pop(0).remove()
    pl2.pop(0).remove()
    pl3.pop(0).remove()
    pl4.pop(0).remove()
    
    plt.show()
    
    return