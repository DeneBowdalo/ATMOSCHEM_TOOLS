import csv
from netCDF4 import Dataset
import glob
import modules
import datetime
import numpy as np
import pandas as pd

start_year = 2005
end_year = 2011

year_range = range(start_year,end_year+1)

s_dt = datetime.datetime(start_year,1,1,0,0)
e_dt = datetime.datetime(end_year+1,1,1,0,0)

dt_range = pd.date_range(start = s_dt,end = e_dt, freq = 'H')[:-1]
n_hours = len(dt_range)

def date_process(date,time):
    year=(date//10000)
    month=((date-year*10000)//100)
    day=(date-year*10000-month*100)

    hour=time//100
    min=(time-hour*100)

    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
        np.int(hour[i]),np.int(min[i]),0)- \
        datetime.datetime(2005,1,1,0,0,0) \
        for i in range(len(year))]
    
    processed_dates=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
    return processed_dates

obs_lat = 16.848
obs_lon = -24.871

model_root_grp = Dataset('/work/home/db876/plotting_tools/model_files/GEOS_CHEM_SURFACE_O3_2005_2010_v90103_2x2.5_GEOS5.nc')
lat_e = model_root_grp.variables['lat_edges'][:]
lon_e = model_root_grp.variables['lon_edges'][:]
lat_c = model_root_grp.variables['lat_centre'][:]
lon_c = model_root_grp.variables['lon_centre'][:]
n_lats = len(lat_c)
n_lons = len(lon_c)
grid_size = n_lats * n_lons 
lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
cvo_index =  (lon_n*n_lons)+lat_n

start_dt = datetime.datetime(start_year,1,1,0,0)                                                                                                                                                                                          
end_dt = datetime.datetime(end_year,12,31,23,0)
dt_range = pd.date_range(start = start_dt,end = end_dt, freq = 'H')
yyyymmdd = []
hhmm = []
for i in dt_range:
    yyyymmdd.append(i.strftime("%Y%m%d"))
    hhmm.append(i.strftime("%H%M"))
yyyymmdd_int = np.array(yyyymmdd).astype('int')
hhmm_int = np.array(hhmm).astype('int')
all_time = date_process(yyyymmdd_int,hhmm_int)

all_files = glob.glob('/work/home/db876/modelling/GEOS_CHEM/v90103/2x2.5/2x2.5_GEOS5_fullgrid/run/planeflight/logs/2*')
all_files = sorted(all_files)

o3 = np.empty(n_hours)
co = np.empty(n_hours)
no = np.empty(n_hours)
no2 = np.empty(n_hours)
c2h6 = np.empty(n_hours)
c3h8 = np.empty(n_hours)
isoprene = np.empty(n_hours)
oh = np.empty(n_hours)
ro2 = np.empty(n_hours)

start = 0 
for num in range(len(all_files)):  
    print num                                                                                                                                                                                                   
    data = np.load(all_files[num])
    current_o3 = data['O3'][cvo_index::grid_size]*1e9
    try:
        start+=len(current_o3)
        end+=len(current_o3)
    except:
        start = 0
        end = len(current_o3)
    o3[start:end] = current_o3
    co[start:end] = data['CO'][cvo_index::grid_size]*1e9
    no[start:end] = data['NO'][cvo_index::grid_size]*1e9
    no2[start:end] = data['NO2'][cvo_index::grid_size]*1e9
    c2h6[start:end] = data['C2H6'][cvo_index::grid_size]*1e9
    c3h8[start:end] = data['C3H8'][cvo_index::grid_size]*1e9
    isoprene[start:end] = data['ISOP'][cvo_index::grid_size]*1e9
    oh[start:end] = data['OH'][cvo_index::grid_size]*1e9
    ro2[start:end] = data['RO2'][cvo_index::grid_size]*1e9
    
all_data = np.vstack((yyyymmdd,hhmm,all_time,o3,co,c2h6,c3h8,isoprene,no,no2,oh,ro2))
all_data = np.transpose(all_data)

speciesheader = ['YYYYMMDD','HHMM','DAYS_SINCE_2006','O3(ppbv)','CO(ppbv)','C2H6(ppbv)','C3H8(ppbv)','ISOPRENE(ppbv)','NO(ppbv)','NO2(ppbv)','OH(ppbv)','RO2(ppbv)']
b = open('GEOSCHEM_v90103_2005_2012_ALLSPECIES_HOURLY.txt', 'w')
a = csv.writer(b)
a.writerow(speciesheader)
a.writerows(all_data)
b.close()



