import numpy as np
import csv
import glob
import datetime as datetime 
import matplotlib.pyplot as plt
import logging as log
import nappy
import scipy.signal as sig
import scipy.stats as stats
import lomb_edit
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft

#reads in the model data
def readfile(filename, location):

    for files in filename:
        print files
        reader=csv.reader(open(files,'rb'), delimiter=' ', skipinitialspace = True)
        for row in reader:
            if row[1] == location: 
                new=row[2:]
              
                try:    
                    big.append(new)
                    
                except:
                    big=[new]
            
            if row[0] == 'POINT':
                names = row[2:]
   
    big=np.array(big)
    big=np.float64(big)

  
                 
    return big, names

def NOAA_data_reader(filename):
    for files in filename:
        print files
        reader=csv.reader(open(files,'rb'), delimiter=' ', skipinitialspace=True)
        for row in reader:
            new = row[:]
            try:
                data.append(new)
            except:
                data=[new]

    data=np.array(data)
    year = data[:,1]  
    month= data[:,2] 
    day  = data[:,3]
    almost_date = [a+b for a,b in zip(year, month)]
    date = [a+b for a,b in zip(almost_date, day)] 
    date = np.array(date)
    date = date.astype(int)
    hour = data[:,4] 
    time = [i+'00' for i in hour]
    time = np.array(time)
    time = time.astype(int)
    vals = data[:,5]
    vals = np.float64(vals)
    return date, time, vals

#do I need to read everything in 
try: 
    names
except NameError:
# Readin the model output
    model , names = readfile(glob.glob("../plane.log.20*"),'010')
# Processes the date for Model 
    year=(model[:,0]//10000)
    month=((model[:,0]-year*10000)//100)
    day=(model[:,0]-year*10000-month*100)
    
    hour=model[:,1]//100
    min=(model[:,1]-hour*100)
    
    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2006,1,1,0,0,0) \
              for i in range(len(year))]
    
    since2006=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]


#now read in the observations
    date, time, vals =  NOAA_data_reader(glob.glob("o3*"))
    valid = vals >= 0
    vals = vals[valid]
    date = date[valid]
    time = time[valid]
# Process NOAA obs time
    year=(date//10000)
    month=((date-year*10000)//100)
    day=(date-year*10000-month*100)

    hour=time//100
    min=(time-hour*100)

    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2006,1,1,0,0,0) \
              for i in range(len(year))]

    since2006_2=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]

#Convert model time array into numpy array

since2006=np.array(since2006)
since2006_2=np.array(since2006_2)

#Plot them all up. 
#Plot them all up. 
fig=plt.figure(figsize=(14,12))

#Plot up standard conc. v time plot
k=names.index('O3')
ax1= fig.add_subplot(3, 1, 1)
fig.subplots_adjust(hspace=0.3)
plt.plot(since2006_2,vals, color='black', label='Arrival Heights Obs.')
plt.plot(since2006, model[:,k]*1e9, color='green', label='GEOS v9.01.03 4x5 ')
plt.grid(True)
leg=plt.legend(loc=3)
leg.get_frame().set_alpha(0.4)
plt.xlabel('Time (Days since 2006)')
plt.ylabel('Concentration, ppbV' )
plt.title('O3 Conc. V Time')

#Plot up Fourier plot

ax2= fig.add_subplot(3, 1, 2)
fig.subplots_adjust(hspace=0.3)

#Plot axis period lines and labels
annotate_line_y=np.arange(1e-1,1e3,1)
freq_max= (float(1./24))*2
freq_max_2 = [freq_max]*len(annotate_line_y)
freq_day = [1]*len(annotate_line_y)
freq_year = [345]*len(annotate_line_y)
plt.plot(freq_max_2, annotate_line_y,'r--', alpha=0.4)
plt.plot(freq_day, annotate_line_y,'r--', alpha=0.4)
plt.plot(freq_year, annotate_line_y,'r--',alpha=0.4)
plt.text(freq_max, 0.1, '2 Hours', fontweight='bold')
plt.text(1, 0.1, '1 Day', fontweight='bold')
plt.text(345, 0.1, '1 Year', fontweight='bold')

fft_obs=np.abs(fft(vals))**2
samp_spacing=float(1./24.)
frequencies = fftfreq(vals.size, d = samp_spacing)
keep_freq = frequencies>0
frequencies, fft_obs = frequencies[keep_freq], fft_obs[keep_freq]
periods = 1 / frequencies
plt.loglog(periods,fft_obs,'kx', markersize=2, label='Arrival Heights Obs.')
model_cut=model[:,k]*1e9
fft_model=np.abs(fft(model_cut))**2
frequencies = fftfreq(model_cut.size, d = samp_spacing)
keep_freq = frequencies>0
frequencies, fft_model = frequencies[keep_freq], fft_model[keep_freq]
periods = 1 / frequencies
plt.loglog(periods,fft_model,'gx', alpha =0.75,  markersize=2, label='GEOS v9.01.03 4x5 ')
plt.grid(True)
leg=plt.legend(loc=7)
leg.get_frame().set_alpha(0.4)
plt.xlabel('Frequency (days-1)')
plt.ylabel('Power')
plt.title('Fourier O3 Power V Freq.')


#Plot up lomb-scargle plot
ax3= fig.add_subplot(3, 1, 3)

#Plot axis period lines and labels
annotate_line_y=np.arange(1e-10,1e4,1)
freq_max= (float(1./24))*2
freq_max_2 = [freq_max]*len(annotate_line_y)
freq_day = [1]*len(annotate_line_y)
freq_year = [345]*len(annotate_line_y)
plt.plot(freq_max_2, annotate_line_y,'r--', alpha=0.4)
plt.plot(freq_day, annotate_line_y,'r--', alpha=0.4)
plt.plot(freq_year, annotate_line_y,'r--',alpha=0.4)
plt.text(freq_max, 1e-10, '2 Hours', fontweight='bold')
plt.text(1, 1e-10, '1 Day', fontweight='bold')
plt.text(345, 1e-10, '1 Year', fontweight='bold')

fa, fb, nout, jmax, prob = lomb_edit.fasper(since2006_2, vals, 2.0001, 6.)
plt.loglog(1./fa, fb,'kx',markersize=2, label='Arrival Heights Obs.')
fx, fy, nout, jmax, prob2 = lomb_edit.fasper(since2006, model[:,k]*1e9, 2.0001, 6.)
plt.loglog(1./fx, fy, 'gx', alpha = 0.75,markersize=2, label='GEOS v9.01.03 4x5')
plt.grid(True)
leg=plt.legend(loc=7)
leg.get_frame().set_alpha(0.4)
#plt.ylim(1e-10,1e4)
plt.xlabel('Frequency (days-1)')
plt.ylabel('Power')
plt.title('Lomb-Scargle O3 Power V Freq.')

#plt.savefig('O3_capeverde_comparison_plots.ps', dpi = 200)
plt.show()
