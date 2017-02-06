import numpy as np
#import modules
import numpy.fft
import matplotlib.pyplot as plt
import scipy.signal
import glob
import csv
from scipy import stats
import datetime
import pandas as pd

def date_process(date,time):
    year=(date//10000)
    month=((date-year*10000)//100)
    day=(date-year*10000-month*100)

    hour=time//100
    min=(time-hour*100)

    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
        np.int(hour[i]),np.int(min[i]),0)- \
        datetime.datetime(2006,1,1,0,0,0) \
        for i in range(len(year))]
    
    processed_dates=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
    return processed_dates


#all time in days since 2006

start_datetime = datetime.datetime(2006,1,1,0,0)

#create hourly array of times and then average all species between this.
start_datetime = datetime.datetime(2006,1,1,0,0)
end_datetime = datetime.datetime(2015,1,1,0,0)

all_datetimes = pd.date_range(start_datetime,end_datetime,freq='H')
all_datetimes = all_datetimes[:-1]
all_yyyymmdd =np.array([int(d.strftime('%Y%m%d')) for d in all_datetimes])
all_hhmm =np.array([int(d.strftime('%H%M')) for d in all_datetimes])
all_time = date_process(all_yyyymmdd,all_hhmm)

print len(all_datetimes)

start_check_datetime = datetime.datetime(2005,12,31,23,30)
end_check_datetime = datetime.datetime(2014,12,31,23,30)
check_datetimes = pd.date_range(start_check_datetime,end_check_datetime,freq='H')
check_yyyymmdd =np.array([int(d.strftime('%Y%m%d')) for d in check_datetimes])
check_hhmm =np.array([int(d.strftime('%H%M')) for d in check_datetimes])
check_time = date_process(check_yyyymmdd,check_hhmm)

#load in hydrocarbon data
#converting from pptv to ppbv

data = np.loadtxt('hydrocarbon_cvo.txt',skiprows=41)


hydro_time = data[:,0]
ethane = data[:,1]/1e3
ethane_flags = data[:,2]
propane = data[:,3]/1e3
propane_flags = data[:,4]
isobutane = data[:,5]/1e3
isobutane_flags = data[:,6]
nbutane = data[:,7]/1e3
nbutane_flags = data[:,8]
acetylene = data[:,9]/1e3
acetylene_flags = data[:,10]
isopentane = data[:,11]/1e3
isopentane_flags = data[:,12]
npentane = data[:,13]/1e3
npentane_flags = data[:,14]

test = ethane_flags != 0
ethane[test] = -99999
test = propane_flags != 0
propane[test] = -99999
test = isobutane_flags != 0
isobutane[test] = -99999
test = nbutane_flags != 0
nbutane[test] = -99999
test = acetylene_flags != 0
acetylene[test] = -99999
test = isopentane_flags != 0
isopentane[test] = -99999
test = npentane_flags != 0
npentane[test] = -99999

hydro_datetime = []
hydro_yyyymmdd = []
hydro_hhmm = []
for t in hydro_time:
    hydro_datetime.append(start_datetime + datetime.timedelta(days=t))
hydro_datetime = np.array(hydro_datetime)

#load in co data
co_files = glob.glob('CO_CVO/*')
co_files = sorted(co_files)

co_time = []
co = []
co_flags = []

print 'CO'
for i in co_files:
    print i
    skip_n = -1
    with open(i, 'rU') as f:
        reader = csv.reader(f,delimiter=' ')
        counter = 0
        for row in reader:
            row = [var for var in row if var]
            
            if counter == 0:
                skip_n = int(row[0])+1
            if counter >= skip_n:
                current_time = float(row[0])
                current_var = float(row[1])
                current_var_flag = float(row[2])
                
                co_time.append(current_time)
                co.append(current_var)
                co_flags.append(current_var_flag)
            counter+=1

co_time = np.array(co_time)
co = np.array(co)
co_flags = np.array(co_flags)

co_datetime = []
co_yyyymmdd = []
co_hhmm = []
for t in co_time:
    co_datetime.append(start_datetime + datetime.timedelta(days=t))
co_datetime = np.array(co_datetime)

test = co_flags != 0
co[test] = -99999


#load in o3 data
o3_files = glob.glob('O3_CVO/*')
o3_files = sorted(o3_files)

o3_time = []
o3 = []
o3_flags = []

print 'O3'
for i in o3_files:
    print i
    skip_n = -1
    with open(i, 'rU') as f:
        reader = csv.reader(f,delimiter=' ')
        counter = 0
        for row in reader:
            row = [var for var in row if var]
            
            if counter == 0:
                skip_n = int(row[0])+1
            if counter >= skip_n:
                current_time = float(row[0])
                current_var = float(row[1])
                try:
                    current_var_flag = float(row[2])
                except:
                    current_var_flag = 0.
                o3_time.append(current_time)
                o3.append(current_var)
                o3_flags.append(current_var_flag)
            counter+=1

o3_time = np.array(o3_time)
o3 = np.array(o3)
o3_flags = np.array(o3_flags)

test = o3_flags != 0
o3[test] = -99999


o3_datetime = []
o3_yyyymmdd = []
o3_hhmm = []
for t in o3_time:
    o3_datetime.append(start_datetime + datetime.timedelta(days=t))
o3_datetime = np.array(o3_datetime)

#load in noxy data
#converting from pptv to ppbv
nox_files = glob.glob('NOX_CVO/*')
nox_files = sorted(nox_files)

nox_time = []
no = []
no_flags = []
no2 = []
no2_flags = []
noy = []
noy_flags = []

print 'NOX'
for i in nox_files:
    print i
    skip_n = -1
    with open(i, 'rU') as f:
        reader = csv.reader(f,delimiter=' ')
        counter = 0
        for row in reader:
            row = [var for var in row if var]
            
            if counter == 0:
                skip_n = int(row[0])+1
            if counter >= skip_n:
                current_time = float(row[0])
                current_var_1 = float(row[1])
                current_var_flag_1 = float(row[2])
                current_var_2 = float(row[3])
                current_var_flag_2 = float(row[4])
                current_var_3 = float(row[5])
                current_var_flag_3 = float(row[6])
                
                
                nox_time.append(current_time)
                no.append(current_var_1)
                no_flags.append(current_var_flag_1)
                no2.append(current_var_2)
                no2_flags.append(current_var_flag_2)
                noy.append(current_var_3)
                noy_flags.append(current_var_flag_3)
            counter+=1


nox_time = np.array(nox_time)
no = np.array(no)/1e3
no_flags = np.array(no_flags)
no2 = np.array(no2)/1e3
no2_flags = np.array(no2_flags)
noy = np.array(noy)/1e3
noy_flags = np.array(noy_flags)

test = no_flags != 0
no[test] = -99999
test = no2_flags != 0
no2[test] = -99999
test = noy_flags != 0
noy[test] = -99999


nox_datetime = []
nox_yyyymmdd = []
nox_hhmm = []
for t in nox_time:
    nox_datetime.append(start_datetime + datetime.timedelta(days=t))
nox_datetime = np.array(nox_datetime)

#average all data into hours

yyyymmdd = [d.strftime('%Y%m%d') for d in all_datetimes]
hhmm = [d.strftime('%H%M') for d in all_datetimes]

ave_ethane = np.empty(len(all_datetimes))
ave_propane = np.empty(len(all_datetimes))
ave_isobutane = np.empty(len(all_datetimes))
ave_nbutane = np.empty(len(all_datetimes))
ave_acetylene = np.empty(len(all_datetimes))
ave_isopentane = np.empty(len(all_datetimes))
ave_npentane = np.empty(len(all_datetimes))
ave_co = np.empty(len(all_datetimes))
ave_o3 = np.empty(len(all_datetimes))
ave_no = np.empty(len(all_datetimes))
ave_no2 = np.empty(len(all_datetimes))
ave_noy = np.empty(len(all_datetimes))

hydro_indices = np.searchsorted(check_time,hydro_time)
o3_indices = np.searchsorted(check_time,o3_time)
co_indices = np.searchsorted(check_time,co_time)
nox_indices = np.searchsorted(check_time,nox_time)

for i in range(len(all_datetimes)):
    print i

    pull_hydro = hydro_indices == i
    pull_o3 = o3_indices == i   
    pull_co = co_indices == i  
    pull_nox = nox_indices == i 
    
    ethane_pulled = ethane[pull_hydro]
    propane_pulled = propane[pull_hydro]
    isobutane_pulled = isobutane[pull_hydro]
    nbutane_pulled = nbutane[pull_hydro]
    acetylene_pulled = acetylene[pull_hydro]
    isopentane_pulled = isopentane[pull_hydro]
    npentane_pulled = npentane[pull_hydro]
    co_pulled = co[pull_co]
    o3_pulled = o3[pull_o3]
    no_pulled = no[pull_nox]
    no2_pulled = no2[pull_nox]
    noy_pulled = noy[pull_nox]
    
    if len(ethane_pulled) > 0:
        ave_ethane[i] = np.average(ethane_pulled)
    else:
        ave_ethane[i] = -99999
        
    if len(propane_pulled) > 0:
        ave_propane[i] = np.average(propane_pulled)
    else:
        ave_propane[i] = -99999
        
    if len(isobutane_pulled) > 0:
        ave_isobutane[i] = np.average(isobutane_pulled)
    else:
        ave_isobutane[i] = -99999
                
    if len(nbutane_pulled) > 0:
        ave_nbutane[i] = np.average(nbutane_pulled)
    else:
        ave_nbutane[i] = -99999

    if len(acetylene_pulled) > 0:
        ave_acetylene[i] = np.average(acetylene_pulled)
    else:
        ave_acetylene[i] = -99999

    if len(isopentane_pulled) > 0:
        ave_isopentane[i] = np.average(isopentane_pulled)
    else:
        ave_isopentane[i] = -99999
        
    if len(npentane_pulled) > 0:
        ave_npentane[i] = np.average(npentane_pulled)
    else:
        ave_npentane[i] = -99999
        
    if len(co_pulled) > 0:
        ave_co[i] = np.average(co_pulled)
    else:
        ave_co[i] = -99999
        
    if len(o3_pulled) > 0:
        ave_o3[i] = np.average(o3_pulled)
    else:
        ave_o3[i] = -99999
        
    if len(no_pulled) > 0:
        ave_no[i] = np.average(no_pulled)
    else:
        ave_no[i] = -99999
        
    if len(no2_pulled) > 0:
        ave_no2[i] = np.average(no2_pulled)
    else:
        ave_no2[i] = -99999
        
    if len(noy_pulled) > 0:
        ave_noy[i] = np.average(noy_pulled)
    else:
        ave_noy[i] = -99999

speciesheader = ['YYYYMMDD','HHMM','DAYS_SINCE_2006','ETHANE(ppbv)','PROPANE(ppbv)','ISOBUTANE(ppbv)','NBUTANE(ppbv)','ACETYLENE(ppbv)','ISOPENTANE(ppbv)','NPENTANE(ppbv)','CO(ppbv)','O3(ppbv)','NO(ppbv)','NO2(ppbv)','NOY(ppbv)']

all_data = np.vstack((yyyymmdd,hhmm,all_time,ave_ethane,ave_propane,ave_isobutane,ave_nbutane,ave_acetylene,ave_isopentane,ave_npentane,ave_co,ave_o3,ave_no,ave_no2,ave_noy))
all_data = np.transpose(all_data)

np.save('CVO_ALLSPECIES_HOURLY_2006_2014',all_data)

b = open('CVO_ALLSPECIES_HOURLY_2006_2014.txt', 'w')
a = csv.writer(b)
a.writerow(speciesheader)
a.writerows(all_data)
b.close()









