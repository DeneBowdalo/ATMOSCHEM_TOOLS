import numpy as np
import scipy.stats
import scipy.signal
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
import numpy.fft
import glob
import sys
import multiprocessing
import operator
import collections
import re
import matplotlib.pyplot as plt
import datetime
import csv
from datetime import datetime

#files = glob.glob('original_files/L3_tropo_ozone_vmr_may05.txt')
files = glob.glob('original_files/L3*')
amended_files = []

print files

remove_chars = {'jan':'01','feb':'02','mar':'03','apr':'04','may':'05','jun':'06','jul':'07','aug':'08','sep':'09','oct':'10','nov':'11','dec':'12','L3_tropo_ozone_vmr_':'','.txt':'','original_files/':''}
for filename in files:
    for char in remove_chars:
        try:
            stripped = stripped.replace(char,remove_chars[char])
        except:
            stripped = filename.replace(char,remove_chars[char])

    amended_files.append(stripped)
    del stripped

print amended_files
amended_files = [datetime(int(i[2:4]),int(i[0:2]),1) for i in amended_files]
amended_inds = np.argsort(amended_files)
files = np.array(files)
files = files[amended_inds]

counter_2 = 0

#create big array#
all_data = np.empty((len(files),120,288))

print files
for file_i in range(len(files)):
    print files[file_i]
    reader=csv.reader(open(files[file_i],'rb'), delimiter=',')
    counter_2=0
    lat_band=[]
    lat_counter = 0
    for cut in reader:
         
        #skip header info
        if counter_2 < 3:
            counter_2+=1
            continue
        counter_2+=1
        
        #read each lat band 
        cut =  str(cut)

        print cut
        cut = cut[3:-2]
        split_j = [cut[i:i+3] for i in range(0, len(cut), 3)]

        #determine end of lat band, and remove text
        #then put lat band into big array
        if len(split_j) == 18:
            split_j = split_j[0:13]
            #print split_j
            lat_band = np.append(lat_band,split_j)
            #print lat_band
            print lat_band
            print len(lat_band)
            #print file_i,lat_counter
            #check for invalid data points and convert values to floats
            lat_band = np.float64(lat_band)
            inv_test = lat_band == 999
            valid_test = lat_band != 999
            lat_band[inv_test] = -99999 
            #divide by 5 to get to ppbv
            lat_band[valid_test] = lat_band[valid_test]/5. 
            all_data[file_i,lat_counter,:] = lat_band
            lat_band=[]
            lat_counter+=1
            continue
        
        lat_band = np.append(lat_band,split_j)   
        

np.save('OMI_O3_trop_ave',all_data) 

   
