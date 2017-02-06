import numpy as np
import matplotlib.pyplot as plt
import modules

pi2 = np.pi*2.
mag = 10
ph = 2.3
period = 365.25
time = np.arange(0,365.25*5,1./24)
waveforma = (mag*(np.cos((pi2*time/period)-(ph))))

mag = 6
ph = 1.
period = 365.25/2.
waveformb = (mag*(np.cos((pi2*time/period)-(ph))))

mag = 2
ph = 3.
period = 365.25/3.
waveformc = (mag*(np.cos((pi2*time/period)-(ph))))

mag = 5
ph = 0.5
period=1.
waveform2a = (mag*(np.cos((pi2*time/period)-(ph))))

mag = 2
ph = 2.
period=1./2
waveform2b = (mag*(np.cos((pi2*time/period)-(ph)))) 

mag = 1
ph = 0.5
period=1./3
waveform2c = (mag*(np.cos((pi2*time/period)-(ph)))) 

#mag = 4
#ph = 2.5
#period = 20
#waveform3 = (mag*(np.cos((pi2*time/period)-(ph))))

#mag = 2
#ph = 2.4
#period = 10
#waveform4 = (mag*(np.cos((pi2*time/period)-(ph))))

noise = np.random.normal(0,20,len(waveforma))
noise_var = np.nanvar(noise)
#waveform = waveform+noise

#noise[0:500] = np.NaN

#waveforma =60+ (waveform+waveform2+waveform3)
#waveforma =(waveform2)

#waveform_sum_beforegap = waveforma+waveformb+waveformc+waveform2a+waveform2b+waveform2c+noise

#print np.nanvar(waveform_sum_beforegap)

#waveforma[0:600] = np.NaN
#waveformb[0:600] = np.NaN
#waveformc[0:600] = np.NaN
#waveform2a[0:600] = np.NaN
#waveform2b[0:600] = np.NaN
#waveform2c[0:600] = np.NaN
#waveforma[2000:6200] = np.NaN
#waveformb[2000:6200] = np.NaN 
#waveformc[2000:6200] = np.NaN 
#waveform2a[2000:6200] = np.NaN 
#waveform2b[2000:6200] = np.NaN 
#waveform2c[2000:6200] = np.NaN 
#time[0:600] = np.NaN
#time[2000:6200] = np.NaN
#noise[0:600] = np.NaN
#noise[2000:6200] = np.NaN

waveform_sum = waveforma+waveformb+waveformc+waveform2a+waveform2b+waveform2c+noise
#waveform_sum[0:600] = np.NaN
#waveform_sum[2000:5200] = np.NaN
#time[0:600] = np.NaN
#time[2000:5200] = np.NaN

waveform_sum[0:600] = -99999
waveform_sum[2000:20200] = -99999
time[0:600] = -99999
time[2000:20200] = -99999
new_time = time[time > -99999]
waveform_sum = waveform_sum[waveform_sum > -99999]

#print np.nanvar(waveform_sum)

#waveform_sum[np.isnan(waveform_sum)] = waveform_sum_beforegap[np.isnan(waveform_sum)]

#print np.nanmean(waveform_sum)
#print (100./(365.25*24))*(4200+600)

#print (10**2+5**2+4**2+2**2)/2. , (10**2+5**2)/2. , (4**2+2**2)/2.
#print (5**2+2.5**2+3**2+1**2)/2.

#print np.nanmean([np.abs(np.diff(i,0))**2 for i in waveform_sum])
#print np.nanvar(waveform_sum),np.nanvar(waveforma+waveformb+waveformc),np.nanvar(waveform2a+waveform2b+waveform2c),np.nanvar(waveforma+waveformb+waveformc)+np.nanvar(waveform2a+waveform2b+waveform2c)

est_noise_waveform =  waveform_sum - (waveforma[time > -99999]+waveformb[time > -99999]+waveformc[time > -99999]+waveform2a[time > -99999]+waveform2b[time > -99999]+waveform2c[time > -99999])
#est_noise_waveform =  waveform_sum - (waveforma+waveformb+waveformc+waveform2a+waveform2b+waveform2c)

periodic_var = (10**2/2.)+(6**2/2.)+(2**2/2.)+(5**2/2.)+(2**2/2.)+(1**1/2.)
actual_var = noise_var+periodic_var
est_noise_var = np.nanvar(est_noise_waveform)
print 'actual var = ', actual_var
print 'periodic var = ', periodic_var
print 'noise var = ', noise_var
print 'est noise var = ', est_noise_var

plt.plot(time,est_noise_waveform)
plt.show()

#np.nanvar(waveforma)+np.nanvar(waveformb)+np.nanvar(waveformc)+np.nanvar(waveform2a)+np.nanvar(waveform2b)+np.nanvar(waveform2c)#+np.nanvar(noise)

#time = time[~np.isnan(time)]
#waveform_sum = waveform_sum[~np.isnan(waveform_sum)] 

#periods,mag,ph,fr,fi,amp_corr = modules.take_lomb(time,waveform_sum,1,1./24.)


#print np.sum(mag**2)/2.

#periods,key_mag,ph,fr,fi = modules.take_lomb_spec(time,waveform_sum,w=True,key_periods=[10.,20.,30.,50.])

#print np.sum(key_mag**2)/2.
