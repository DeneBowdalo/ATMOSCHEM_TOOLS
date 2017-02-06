import matplotlib.pyplot as plt 
from scipy.fftpack import fftfreq,fft
from scipy import signal
import numpy as np
from cmath import *
import operator
import lomb_phase
import lomb
import old_lomb
import modules
from multiprocessing import Pool
import time as now
from operator import itemgetter
import test_lomb
import fast_lomb_phase
import lomb_townsend
import matplotlib.ticker
from matplotlib.ticker import FuncFormatter

#Time array in hourly intervals



hourly_step = 1./24

first_year = np.arange(0,365.25*4,hourly_step)
last_year = np.arange(365.25*2,365.25*4,hourly_step)

#a = np.arange(0,365.25*1,hourly_step)
a2 = np.arange(365.25*4,365.25*8,hourly_step)

a = np.arange(0,24,0.1)

print len(a)

#a = np.delete(a,-1)

print len(a)

b = np.arange(365.25*2,365.25*4,hourly_step)


combined_time = np.append(a,b)

pi2 = np.pi*2

daily_amplitude = 2
half_annual_amplitude = 5
annual_amplitude = 10

cos_waveform1 = daily_amplitude*(np.cos((pi2*a/1)-(0)))
cos_waveform2 = half_annual_amplitude*(np.cos((pi2*a/(365.25/2))-0))
cos_waveform3 = annual_amplitude*(np.cos((pi2*a/365.25)+(0)))

cs_daily = annual_amplitude*(np.cos((pi2*a/24)+(0)))

sin_waveform3 = annual_amplitude*(np.sin((pi2*a/365.25)+(0)))

cos_waveform1_a = daily_amplitude*(np.sin((pi2*a2/1)-(0)))
cos_waveform2_a = half_annual_amplitude*(np.sin((pi2*a2/(365.25/2))-0))
cos_waveform3_a = annual_amplitude*(np.sin((pi2*a2/365.25)+(0)))


#cos_waveform3=np.append(cos_waveform3,shifted_waveform)

waveform = 50 +  (cs_daily)
#waveform2 = 50 +(sin_waveform3)

#waveform = 50 + (cos_waveform1+cos_waveform3)
#waveform = 50 + (cos_waveform1 + cos_waveform2 + cos_waveform3)
#waveform2 = 50 + (cos_waveform1_a +cos_waveform2_a+cos_waveform3_a)

#waveform = np.append(waveform,waveform2)

#waveform = 50 + (cos_waveform1)
#waveform2 = 50 + (cos_waveform4)


#
remove = np.arange(100,2100,1)

remove2 = np.arange(2500,4500,1)

remove3 = np.arange(6200,8200,1)

remove4 = np.arange(12000,14000,1)

remove5 = np.arange(19000,21000,1)

remove6 = np.arange(23000,25000,1)

remove7 = np.arange(30000,32000,1)

#for i in remove:
#	waveform[i] = -99

#for i in remove2:
#    waveform[i] = -99

#for i in remove3:
#   waveform[i] = -99

#for i in remove4:
#    waveform[i] = -99

#for i in remove5:
#    waveform[i] = -99

#for i in remove6:
#    waveform[i] = -99

#for i in remove7:
#    waveform[i] = -99

#test = waveform >=0

#a = a[test]
#waveform = waveform[test]

#first year/last year waveform

#year_waveform  = annual_amplitude*(np.cos((pi2*first_year/365.25)-(np.pi/2)))
#last_year_waveform = annual_amplitude*(np.cos((pi2*last_year/365.25)-(0)))

#do FFT
#window = signal.flattop(len(waveform), sym=False)
#window = np.hanning(len(waveform))
#waveform = waveform*window

#window= np.hamming(len(waveform))
#waveform_mean = np.mean(waveform)
#waveform  = waveform - waveform_mean

#waveform = waveform*window

val = ((32768*2)*2)*2

#diff = val - len(waveform)

#append_vals = [0]*diff

#waveform = np.append(waveform,append_vals)

print waveform

fft_array = fft(waveform)
fft_cos_amp = fft_array.real
fft_sin_amp = fft_array.imag
fft_mag = np.abs(fft_array)
fft_phase = [phase(i) for i in fft_array]
fft_phase = np.array(fft_phase)
fft_freqs = fftfreq(len(waveform),hourly_step) 

valid = fft_freqs > 0
fft_mag, fft_phase, fft_freqs, fft_cos_amp, fft_sin_amp = fft_mag[valid], fft_phase[valid], fft_freqs[valid], fft_cos_amp[valid], fft_sin_amp[valid]

fft_periods = 1./fft_freqs
fft_cos_amp = fft_cos_amp/(len(a)/2)
fft_sin_amp = fft_sin_amp/(len(a)/2)
fft_mag = fft_mag/len(a)
fft_mag = fft_mag*2
#fft_mag = fft_mag* np.sqrt(8./3)
#amp_corr = 1./(sum(window)/len(window))
#fft_mag = fft_mag * amp_corr

daily_index = np.argmin(np.abs(fft_periods - 1))
half_annual_index = np.argmin(np.abs(fft_periods - 365.25/2))
annual_index = np.argmin(np.abs(fft_periods - 365.25))

#for num in range(len(fft_phase)):
#   if fft_phase[num] < 0:
#       diff = np.pi- np.abs(fft_phase[num]) 
#       fft_phase[num] = np.pi - diff
#   elif fft_phase[num] > 0:
#       diff = np.pi- np.abs(fft_phase[num])
#       fft_phase[num] = np.pi + diff


print 'Daily phase = ', fft_phase[daily_index]
print 'Half annual phase = ', fft_phase[half_annual_index]
print 'Lomb Annual phase = ', fft_phase[annual_index]



#fft_mag = fft_mag*5
#fft_mag = fft_mag/(len(a)/10)
#fft_periods = 1./fft_freqs
#for num in range(len(fft_phase)):
#	if fft_phase[num] < 0:
#		diff = np.pi- np.abs(fft_phase[num]) 
#		fft_phase[num] = np.pi - diff
#	elif fft_phase[num] > 0:
#		diff = np.pi- np.abs(fft_phase[num])
#		fft_phase[num] = np.pi + diff

#print 'mag sum', np.sum(fft_mag)

	
#wk1, wk2, a, ph = lomb_phase.lomb(a,waveform)  

#window = np.kaiser(len(waveform),5)
#window = signal.flattop(len(waveform), sym=False)
#window= np.hamming(len(waveform))
#waveform_mean = np.mean(waveform)
#waveform  = waveform - waveform_mean

#waveform = waveform*window

wk1,wk2, amp, l_phase = lomb.fasper(a,waveform)
lomb_periods = 1./wk1

#amp_corr = 1./(sum(window)/len(window))
#amp = amp * amp_corr

#window = signal.flattop(len(waveform), sym=False)

#test = waveform >= 0
#a = a[test]
#waveform = waveform[test]

#test = waveform >= 0
#a = a[test]
#waveform = waveform[test]

#window = np.kaiser(len(waveform),5)
#window = signal.flattop(len(waveform), sym=False)
#window= np.hamming(len(waveform))
#waveform_mean = np.mean(waveform)
#waveform  = waveform - waveform_mean

#waveform = waveform*window

#wk1, wk2, amp, ph, fft_phase = lomb_phase.lomb(a,waveform)
#lomb_periods = 1./wk1
#amp_corr = 1./(sum(window)/len(window))
#amp = amp * amp_corr

#daily_index = np.argmin(np.abs(lomb_periods - 1))
#half_annual_index = np.argmin(np.abs(lomb_periods - 365.25/2))
#annual_index = np.argmin(np.abs(lomb_periods - 365.25))

#print 'Lomb Daily phase = ', ph[daily_index]
#print 'Lomb Half annual phase = ', ph[half_annual_index]
#print 'Lomb Annual phase = ', ph[annual_index]

#print 'Lomb Daily phase = ', fft_phase[daily_index]
#print 'Lomb Half annual phase = ', fft_phase[half_annual_index]
#print 'Lomb Annual phase = ', fft_phase[annual_index]
#plot up

fig = plt.figure(figsize=(14, 10 ))
fig.patch.set_facecolor('white')
ax = plt.subplot(111)
#plt.plot(a,cos_waveform1)
#plt.plot(a,cos_waveform2)
#plt.plot(a,waveform)

print waveform
plt.plot(a,waveform, color='red')
#plt.plot(a,waveform2,color='blue',label='Sine')
#plt.semilogx(fft_periods, fft_mag, color='purple')
#plt.semilogx(fft_periods, fft_cos_amp, color = 'purple')
#plt.semilogx(fft_periods,fft_sin_amp, color = 'purple')
#plt.semilogx(lomb_periods, amp, color='purple')
#print np.min(fft_phase)
#print np.max(fft_phase)
#plt.plot(len_s,save)
#plt.semilogx(lomb_periods_2,amp_2)
#plt.plot(a,waveform)

#ax1.semilogx(lomb_periods_phase, wk2_phase)

#ax2 = plt.subplot(212)
#ax2.semilogx(fft_periods, fft_mag)
#ax2.set_xlim([0.9,5])
#ax2.set_ylim([0, 10])
#plt.grid()
#plt.xlabel('Period (Days)', fontsize = 16)
#plt.ylabel('Magnitude (ppbV)', fontsize = 16)
#plt.title('Lomb-Scargle derived FFT', fontsize = 16)
#plt.plot(a, waveform)
#plt.plot(data_points,data)
#plt.plot(a, waveform2)
#plt.plot(combined_time, waveform)
#plt.semilogx(fft_periods, fft_sin_amp)
#plt.semilogx(full_freqs, full_real)
#plt.semilogx(fft_periods,fft_mag)

plt.legend(fontsize=22)
#ax2.semilogx(lomb_periods, amp)
#plt.ylim(0,15)

#plt.xlim(0.99,1.01)
#plt.ylim(0,2.5)
#plt.ylim(0,1)

#plt.ylim(0,12)

def form2(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.2f' % x

xformatter = FuncFormatter(form2)
plt.tick_params(axis='x', which='major', labelsize=14)
plt.tick_params(axis='y', which='major', labelsize=14)
plt.tick_params(direction='out', pad=0.00005)
#ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
xlabel = plt.xlabel('Period (Days)', fontsize = 16, labelpad = 8)
#ylabel = plt.ylabel('Amplitude (ppbV)', fontsize = 16, labelpad = 8)
#ylabel = plt.ylabel('Magnitude (ppbV)', fontsize = 16, labelpad = 8)
ylabel = plt.ylabel('Concentration (ppbV)', fontsize = 16, labelpad = 8)
title = plt.title('Phase', fontsize = 22)
#title = plt.title(r'Synthetic Time Series of $O_3$', fontsize = 22)
#title = plt.title(r'FFT of Synthetic Diurnal + Annual $O_3$ Waveform', fontsize = 22)
#title = plt.title(r'FFT of Synthetic Annual $O_3$ Waveform - Real Amplitudes', fontsize = 22)
#title = plt.title(r'FFT of Synthetic Daily $O_3$ Waveform - Imaginary Amplitudes', fontsize = 22)
#title = plt.title('Lomb-Scargle derived FFT', fontsize = 22)
title.set_y(1.01) 

plt.grid()
#plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
plt.axvline(0, color = 'orange', linestyle = '--')
plt.axvline(365.25/2, color = 'orange', linestyle = '--')
plt.axvline(365.25, color = 'orange', linestyle = '--')
plt.axvline(365.25*1.5, color = 'orange', linestyle = '--')
plt.axvline(365.25*2, color = 'orange', linestyle = '--')

#plt.axvline(1, color = 'red', linestyle = '--')
#plt.axvline(365.25/2, color = 'red', linestyle = '--')
#plt.axvline(365.25, color = 'red', linestyle = '--')
plt.xlim(0,24)

#plt.text(0.6,14,"Daily Period")
#plt.text(50,13,"Half Annual Period")
#plt.text(230,14,"Annual Period")

plt.show()
