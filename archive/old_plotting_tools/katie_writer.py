import csv
import struct
import numpy as np
import glob


filename = glob.glob("logs/plane.log.20*")

count=0
text_file = open("GEOS_hourly_O3.txt", "w")

for files in filename:
    if count == 0:
        lines = csv.reader(open(files,'rb'), delimiter=' ',skipinitialspace=True)
        for row in lines:
            if count == 0:
                test_names = row 
                names = test_names[:-1]
                count+=1            
              
valid_name = 'O3'

if 

for files in filename:
    reader =np.genfromtxt(files,dtype=str, skip_header=1)
    
    try:
        data = np.concatenate((data,reader))      
    except:
        data = names
        
        print len(data)
        print reader.shape
        data = np.vstack((data,reader)) 
 
np.save('GEOS_v90103_4x5_CV_logs', data)




