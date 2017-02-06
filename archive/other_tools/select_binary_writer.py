import csv
import struct
import numpy as np
import glob
from operator import itemgetter

filename = glob.glob("logs/plane.log.20*")

count=0

for files in filename:
    if count == 0:
        lines = csv.reader(open(files,'rb'), delimiter=' ',skipinitialspace=True)
        for row in lines:
            if count == 0:
                test_names = row
                names = test_names[:-1]
                count+=1

o3 = names.index('O3')
temp = names.index('GMAO_TEMP')
pres = names.index('GMAO_PSFC')

for files in filename:
    print files
    reader =np.genfromtxt(files,dtype=str, skip_header=1, usecols=(1,2,3,o3))
    try:
        data = np.concatenate((data,reader))      
    except:
        data = itemgetter(1,2,3,o3)(names)
        data = np.vstack((data,reader))
np.save('gaw_logs_O3', data)




