from netCDF4 import Dataset
import numpy as np
import datetime
import modules
from scipy import interpolate
import matplotlib.pyplot as plt 
from scipy import signal
import cmath

start_year = 1971
end_year = 2008

d0 = datetime.date(start_year, 1, 1)
d1 = datetime.date(end_year+1, 1, 1)
delta = d1 - d0
n_days = delta.days

all_hours = np.arange(0,n_days,1./24.)

group = Dataset('GAW_SURFACE_O3_1971_2009.nc')
site_group = group.groups['hpb']
vals = site_group.variables['o3'][:]
date = site_group.variables['date'][:]
time = site_group.variables['time'][:]

current_time = modules.date_process(date,time,start_year) 

print current_time
print vals

valid = vals > 0
vals = vals[valid]
current_time = current_time[valid]


print current_time

all_hours = np.arange(np.min(current_time),np.max(current_time)+1./48.,1./24.)

f = interpolate.interp1d(current_time, vals)
vals =  f(all_hours)

ofac=4
periods,mag,ph,fr,fi,fft_array,amp_corr = modules.take_fft(all_hours,vals,1./24.,ofac,w=False)

print len(fft_array)
print len(mag)
print len(fr)
print len(fi)
print fr

plt.loglog(periods,mag)

#remove points n around fifth diurnal peak                                                                                                                                      
diurnal_i = min(range(len(periods)), key=lambda i: abs(periods[i]-0.2))
n_points = int((170*ofac)+np.floor(ofac/2.))
rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
mag[rm_points] = 0
fr[rm_points] = 0
fi[rm_points] = 0

#remove points n around quarter diurnal peak
diurnal_i = min(range(len(periods)), key=lambda i: abs(periods[i]-0.25))
n_points = int((190*ofac)+np.floor(ofac/2.))
rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
mag[rm_points] = 0
fr[rm_points] = 0
fi[rm_points] = 0
 
#remove points n around third diurnal peak
diurnal_i = min(range(len(periods)), key=lambda i: abs(periods[i]-0.33333333))
n_points = int((210*ofac)+np.floor(ofac/2.))
rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
mag[rm_points] = 0
fr[rm_points] = 0
fi[rm_points] = 0
 
#remove points n around half diurnal peak
diurnal_i = min(range(len(periods)), key=lambda i: abs(periods[i]-0.5))
n_points = int((230*ofac)+np.floor(ofac/2.))
rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
mag[rm_points] = 0
fr[rm_points] = 0
fi[rm_points] = 0 

#remove points n around diurnal peak
diurnal_i = min(range(len(periods)), key=lambda i: abs(periods[i]-1.))
n_points = int((250*ofac)+np.floor(ofac/2.))
rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
mag[rm_points] = 0
fr[rm_points] = 0
fi[rm_points] = 0 

#remove points n around quarter-annual peak
ha_i = min(range(len(periods)), key=lambda i: abs(periods[i]-91.325))
n_points = int((4*ofac)+np.floor(ofac/2.))
rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
mag[rm_points] = 0
fr[rm_points] = 0
fi[rm_points] = 0

#remove points n around tri-annual peak
ha_i = min(range(len(periods)), key=lambda i: abs(periods[i]-121.75))
n_points = int((4*ofac)+np.floor(ofac/2.))
rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
mag[rm_points] = 0
fr[rm_points] = 0
fi[rm_points] = 0 

#remove points n around half-annual peak
ha_i = min(range(len(periods)), key=lambda i: abs(periods[i]-182.625))
n_points = int((4*ofac)+np.floor(ofac/2.))
rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
mag[rm_points] = 0
fr[rm_points] = 0
fi[rm_points] = 0 

#remove points n around annual peak
annual_i = min(range(len(periods)), key=lambda i: abs(periods[i]-365.25))
n_points = int((4*ofac)+np.floor(ofac/2.))
rm_points = range(annual_i-n_points,(annual_i+n_points)+1)
if rm_points[0] < 0:
    rm_points = range(0,(annual_i+n_points)+1)
mag[rm_points] = 0
fr[rm_points] = 0
fi[rm_points] = 0

#fr = []
#fi = []

#for i in range(len(mag)):
#    m = mag[i]
#    p = ph[i]
#    fr.append(cmath.rect(m,p).real)
#    fi.append(cmath.rect(m,p).imag)

#fr = np.array(fr)
#fi = np.array(fi)

plt.loglog(periods,mag)
plt.show()

F = [0]*((len(fr)*2)+2)
    
#set first real value to average, imag to 0 
F[0] = complex((np.average(vals)*len(vals)),0)

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

print fft_array
print F
print len(fft_array)
print len(F)
#Take ifft and just take real values
ifft_ts = np.fft.ifft(F)
ifft_ts = ifft_ts.astype('float64')

actual_len = len(vals)
ifft_ts = ifft_ts[:actual_len]
ifft_ts = np.average(vals)+ifft_ts

#plt.loglog(periods,mag)
plt.plot(all_hours,vals)
#plt.show()

data = np.vstack((all_hours,ifft_ts))

print data.shape

np.save('hpb_deperiodised_1971_2009',data)

plt.plot(all_hours,ifft_ts)

plt.show()

#periods,mag,ph,fr,fi,fft_array,amp_corr = modules.take_fft(all_hours,ifft_ts,1./24.,ofac,w=False)
#plt.loglog(periods,mag)
#plt.show()

