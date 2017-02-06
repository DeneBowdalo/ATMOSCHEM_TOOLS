import glob
from collections import Counter
import numpy as np
import csv
import datetime as datetime
from itertools import izip

#find valid files, over half of years present
read_file = 'macehead_03041987-29102013.txt' 


def date_process(date,time):
    print date
    year=(date//10000)
    print year
    month=((date-year*10000)//100)
    day=(date-year*10000-month*100)

    hour=time//100
    print hour
    print np.min(hour)
    min=(time-hour*100)

    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(1988,1,1,0,0,0) \
              for i in range(len(year))]

    processed_dates=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
    return processed_dates

name = []
lat = []
lon = []
units = []

#for each location remove headers and process each year and concatenate 
site_counter = 1
read = np.loadtxt(read_file,delimiter=',',dtype="S10,S5,S4,S7",usecols=(0,1,2,3),unpack =True) 	
#read = np.loadtxt(read_file,unpack =True)
big_list = np.array(read)
# for big_list, extract site info and do checks to test if conc needs converting
	
#remove blank data
good_data = big_list[2] != ''
big_list = big_list[:,good_data]

#remove suspect data
good_data = big_list[3] != 'S'
big_list = big_list[:,good_data]		

# convert from microg/m3	
print 'converting units'
# assume 1 atm and 25 degreesC
o3_counter = 0
for conc in big_list[2]:
	conc = np.float64(conc)
	convert = 24.45*conc/48
	big_list[2,o3_counter] = convert
	o3_counter+=1
			
#rearrange date
dates = big_list[0]
year = []
month = []
day = []
for i in dates:
	print i
	year = np.append(year,i[6:10])
	month = np.append(month, i[3:5]) 
	day = np.append(day, i[0:2])

for i in range(len(year)):
	big_list[0,i] = year[i]+month[i]+day[i]

#process dates from date, time to days since start: 03/04/1987
date_con = [s.replace('-', '') for s in big_list[0]]			
date_con = np.int32(date_con)
time_con = [s.replace(':', '') for s in big_list[1]]
time_con = np.int32(time_con)

#print date_con
#print time_con

tester = False
#some times go from 0100 to 2400, as opposed to 0000 to 2300, need to test if this is case, and if so correct
time_count = 0
for i in time_con:
	if i == 2400:
		tester = True
if tester == True:
	for i in time_con:
		if i == 0:
			break
		time_con[time_count] = i-100
		time_count +=1
	converted_time = date_process(date_con,time_con)
	converted_time = np.round(converted_time,decimals=5)

print converted_time
print len(converted_time)
print len(big_list[2,:])
big_array = np.zeros((2,len(converted_time)))
big_array[0,:] = converted_time
big_array[1,:] = big_list[2,:]

#cut data between 1988 and 2013
test =  big_array[0,:] >= 0
big_array = big_array[:,test]

test = big_array[0,:] < 9132
big_array = big_array[:,test]

#create running 5 year chunks
#1988 is leap year
start_year = 1988
end_year = 1993

start_day = -365
end_day = 1461

#have checks to add 1 or 2 days due to differet configs of leap years
leap_year_count = -1
for i in range(21):
	leap_year_count+=1

	if leap_year_count == 0:
		start_day+=365
		end_day+=366
	if leap_year_count == 1:
		start_day+=366
		end_day+=365
	if leap_year_count == 2:
		start_day+=365
		end_day+=365
	if leap_year_count == 3:
		leap_year_count = -1
		start_day+=365
		end_day+=365

	print start_day
	print end_day
	#limit data to 5 year chunk
	test = (big_array[0,:] >= start_day) & (big_array[0,:] < end_day) 
	chunk = big_array[:,test]

	np.save('chunked_binary_files/MH_O3_SFC_%s_%s'%(start_year,end_year),chunk)
	start_year+=1
	end_year+=1
