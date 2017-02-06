import numpy as np
import csv
def NOAA_data_reader_mace_head(filename):
    reader=csv.reader(open(filename,'rb'), delimiter=',', skipinitialspace=True)
    for row in reader:
        new = row[:]
        try:
            data.append(new)
        except:
            data=[new]

    data=np.array(data)
    date = data[:,0]
    time = data[:,1]
    vals_p = data[:,2]
    vals=[]
    for i in vals_p:
        if i == '':
            i = '-999'
        vals.append(i)
    vals = np.array(vals)
    vals = vals.astype(str)
    return date,time,vals

def readfile(filename, location):
    read = np.load(filename)
    names = read[0,2:]
    locs = np.where(read[:,1] == location)
    big = read[locs]
    valid_data = big[:,2:]
    names = names.tolist()
    valid_data = np.float64(valid_data)

    return valid_data, names

model , names = readfile("GEOS_v90103_4x5_CV_logs.npy","112") #112 represents Mace Head
date, time, vals =  NOAA_data_reader_mace_head('o3_hourly_maceh_2006-2013.txt')

k=names.index('GMAO_PSFC')
i=names.index('GMAO_TEMP')
pressure = model[:,k]
temperature = model[:,i]

counter = 0
#conversion

molecular_weight = 48 # for o3  
new_vals = []

print len(vals)
print len(pressure)

for i in vals:
    pressure = pressure[counter]/1013.25   #convert into atm
    #temperature = temperature[counter]+273.15 #convert into K

    calc = 0.082*(i/molecular_weight)*temperature/pressure
    converted =  calc/1000
    new_vals.append(converted)

    counter+=1

#Write to new file 
text_file = open("O3_mace_head_ppbV.txt", "w")
counter2=0
for i in new_vals:
    text_file.write(date[counter2]+','+time[counter2]+','+i+'\n')
    counter2+=1
text_file.close()
