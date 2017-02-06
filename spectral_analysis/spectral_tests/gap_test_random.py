import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from cmath import *
from math import *
import csv
import datetime
import lomb_phase
from scipy import signal
import multiprocessing
import datetime
import time
import modules
from scipy import interpolate
import random

hourly_step = 1./24

b = np.arange(0,214.5,hourly_step)

#b = b[:-1]

print len(b)

pi2 = np.pi*2


p1 = 1.
p2 = 50.
amp1 = 3.
amp2 = 10.
phase_1 = 0.
phase_2 = 0.

a1_ratio = 100./amp1 
a2_ratio = 100./amp2

cos_waveform1 = amp1*(np.cos((pi2*b/p1)-(phase_1)))
cos_waveform2 = amp2*(np.cos((pi2*b/p2)-(phase_2)))

vals = 50 + (cos_waveform1+cos_waveform2)

print '\nActual Amplitude 1 = %s'%(amp1)
print 'Actual Phase 1 = %s'%(phase_1)
print 'Actual Amplitude 2 = %s'%(amp2)
print 'Actual Phase 2 = %s'%(phase_2)
print '------------------------------------------\n'

noise = np.random.normal(0,2,len(vals))
vals = vals+noise

b = np.array(b)
vals = np.array(vals)
orig_b = np.copy(b)
orig_vals = np.copy(vals)

#plt.plot(b,vals)
#plt.show()

#set arrays to put data into
lomb1_a1 = []
lomb1_a2 = []
lomb1_p1 = []
lomb1_p2 = []

lomb2_a1 = []
lomb2_a2 = []
lomb2_p1 = []
lomb2_p2 = []

lomb3_a1 = []
lomb3_a2 = []
lomb3_p1 = []
lomb3_p2 = []

lomb4_a1 = []
lomb4_a2 = []
lomb4_p1 = []
lomb4_p2 = []

lomb5_a1 = []
lomb5_a2 = []
lomb5_p1 = []
lomb5_p2 = []

fft1_a1 = []
fft1_a2 = []
fft1_p1 = []
fft1_p2 = []

fft2_a1 = []
fft2_a2 = []
fft2_p1 = []
fft2_p2 = []

fft3_a1 = []
fft3_a2 = []
fft3_p1 = []
fft3_p2 = []

fft4_a1 = []
fft4_a2 = []
fft4_p1 = []
fft4_p2 = []

#-------------------------------
#Put in gaps
all_frac_points = []
all_lomb_diff = []
all_fft_diff = []

for n in np.arange(0,100,1):
    print n
    all_frac_points.append(n)
    
    n_points = (len(orig_vals)/100.)*n
    n_points = int(n_points)
    
    if n == 0:
        b = orig_b[:]
        vals = orig_vals[:]
    else:
        #make sure range has to not leave in first value and last value, to make sure interpolation can be done to leave full time series.
        frac_points = random.sample(xrange(1,len(orig_vals)-1),n_points)
        frac_points = np.sort(frac_points)

        b = np.delete(orig_b,frac_points)
        vals = np.delete(orig_vals,frac_points)
    
    #put spectral output for each percent gaps in monte carlo loop
    for i in range(1):
        #define arrays to store data in
        lsp1_diff_a1 = []
        lsp1_diff_a2 = []
        lsp1_diff_p1 = []
        lsp1_diff_p2 = []
        
        lsp2_diff_a1 = []
        lsp2_diff_a2 = []
        lsp2_diff_p1 = []
        lsp2_diff_p2 = []
        
        lsp3_diff_a1 = []
        lsp3_diff_a2 = []
        lsp3_diff_p1 = []
        lsp3_diff_p2 = []
        
        lsp4_diff_a1 = []
        lsp4_diff_a2 = []
        lsp4_diff_p1 = []
        lsp4_diff_p2 = []
        
        lsp5_diff_a1 = []
        lsp5_diff_a2 = []
        lsp5_diff_p1 = []
        lsp5_diff_p2 = []
        
        fft1_diff_a1 = []
        fft1_diff_a2 = []
        fft1_diff_p1 = []
        fft1_diff_p2 = []
        
        fft2_diff_a1 = []
        fft2_diff_a2 = []
        fft2_diff_p1 = []
        fft2_diff_p2 = []
        
        fft3_diff_a1 = []
        fft3_diff_a2 = []
        fft3_diff_p1 = []
        fft3_diff_p2 = []
        
        fft4_diff_a1 = []
        fft4_diff_a2 = []
        fft4_diff_p1 = []
        fft4_diff_p2 = []
    
        #Lomb OUT OF BOX
        #-----------------------------------------------
        ofac = 1
        periods,mag,ph,fr,fi,amp_corr = modules.lomb_spectra(b,vals,ofac,hourly_step,w=False)
        
        closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
        closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
        mag1 = mag[closest_period_1]
        phase1 = ph[closest_period_1]
        mag2 = mag[closest_period_2]
        phase2 = ph[closest_period_2]
        
        diff_a1 = (mag1/amp1)*100.
        diff_a2 = (mag2/amp2)*100. 
        diff_p1 = phase_1 - phase1
        if diff_p1 > np.pi:
            diff_p1 = np.abs(-np.pi + (diff_p1 - np.pi))
        elif diff_p1 < -np.pi:
            diff_p1 = np.abs(np.pi - (np.abs(diff_p1) - np.pi)) 
        diff_p2 = phase_2 - phase2
        if diff_p2 > np.pi:
            diff_p2 = np.abs(-np.pi + (diff_p2 - np.pi))
        elif diff_p2 < -np.pi:
            diff_p2 = np.abs(np.pi - (np.abs(diff_p2) - np.pi))
            
        lsp1_diff_a1.append(diff_a1)
        lsp1_diff_a2.append(diff_a2)
        lsp1_diff_p1.append(diff_p1)
        lsp1_diff_p2.append(diff_p2)
        

        #print '\nLOMB OUT OF BOX'
        #print 'Est. Amplitude 1 = ',mag1
        #print 'Est. Phase 1 = ',phase1
        #print 'Est. Amplitude 2 = ',mag2
        #print 'Est. Phase 2 = ',phase2
        #print '------------------------------------------\n'
           
        #Lomb WINDOWED
        #-----------------------------------------------
        ofac = 1
        periods,mag,ph,fr,fi,amp_corr = modules.lomb_spectra(b,vals,ofac,hourly_step,w=True)
        
        closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
        closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
        mag1 = mag[closest_period_1]
        phase1 = ph[closest_period_1]
        mag2 = mag[closest_period_2]
        phase2 = ph[closest_period_2]
        
        diff_a1 = (mag1/amp1)*100.
        diff_a2 = (mag2/amp2)*100.
        diff_p1 = phase_1 - phase1
        if diff_p1 > np.pi:
            diff_p1 = np.abs(-np.pi + (diff_p1 - np.pi))
        elif diff_p1 < -np.pi:
            diff_p1 = np.abs(np.pi - (np.abs(diff_p1) - np.pi)) 
        diff_p2 = phase_2 - phase2
        if diff_p2 > np.pi:
            diff_p2 = np.abs(-np.pi + (diff_p2 - np.pi))
        elif diff_p2 < -np.pi:
            diff_p2 = np.abs(np.pi - (np.abs(diff_p2) - np.pi))
            
        lsp2_diff_a1.append(diff_a1)
        lsp2_diff_a2.append(diff_a2)
        lsp2_diff_p1.append(diff_p1)
        lsp2_diff_p2.append(diff_p2)

        #print '\nLOMB WINDOWED'
        #print 'Est. Amplitude 1 = ',mag1
        #print 'Est. Phase 1 = ',phase1
        #print 'Est. Amplitude 2 = ',mag2
        #print 'Est. Phase 2 = ',phase2
        #print '------------------------------------------\n'
        
        #Lomb WINDOWED, SPECIFIC FREQUENCIES
        #-----------------------------------------------
        
        periods,mag,ph,fr,fi = modules.lomb_specific(b,vals,w=True,key_periods=[p1,p2])
        
        closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
        closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))

        mag1 = mag[closest_period_1]
        phase1 = ph[closest_period_1]
        mag2 = mag[closest_period_2]
        phase2 = ph[closest_period_2]
        
        diff_a1 = (mag1/amp1)*100.
        diff_a2 = (mag2/amp2)*100.
        diff_p1 = phase_1 - phase1
        if diff_p1 > np.pi:
            diff_p1 = np.abs(-np.pi + (diff_p1 - np.pi))
        elif diff_p1 < -np.pi:
            diff_p1 = np.abs(np.pi - (np.abs(diff_p1) - np.pi)) 
        diff_p2 = phase_2 - phase2
        if diff_p2 > np.pi:
            diff_p2 = np.abs(-np.pi + (diff_p2 - np.pi))
        elif diff_p2 < -np.pi:
            diff_p2 = np.abs(np.pi - (np.abs(diff_p2) - np.pi))
            
        lsp3_diff_a1.append(diff_a1)
        lsp3_diff_a2.append(diff_a2)
        lsp3_diff_p1.append(diff_p1)
        lsp3_diff_p2.append(diff_p2)

        #print '\nLOMB WINDOWED, SPECIFIC FREQUENCIES'

        #print 'Est. Amplitude 1 = ',mag1
        #print 'Est. Phase 1 = ',phase1
        #print 'Est. Amplitude 2 = ',mag2
        #print 'Est. Phase 2 = ',phase2
        #print '------------------------------------------\n'
                                                                                                
        #Do FFT interpolation if needed
        if n > 0:
            f = interpolate.interp1d(b, vals)
            b = np.arange(np.min(b),np.max(b),hourly_step) 
            vals =  f(b)
                                                                                                    
        #FFT - OUT OF BOX
        #----------------------------------------------------------------
        ofac=1
        periods,mag,ph,fr,fi,fft_array,amp_corr = modules.fft_spectra(b,vals,ofac,hourly_step,w=False)
        
        closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
        closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
        mag1 = mag[closest_period_1]
        phase1 = ph[closest_period_1]
        mag2 = mag[closest_period_2]
        phase2 = ph[closest_period_2] 
        
        diff_a1 = (mag1/amp1)*100.
        diff_a2 = (mag2/amp2)*100.
        diff_p1 = phase_1 - phase1
        if diff_p1 > np.pi:
            diff_p1 = np.abs(-np.pi + (diff_p1 - np.pi))
        elif diff_p1 < -np.pi:
            diff_p1 = np.abs(np.pi - (np.abs(diff_p1) - np.pi)) 
        diff_p2 = phase_2 - phase2
        if diff_p2 > np.pi:
            diff_p2 = np.abs(-np.pi + (diff_p2 - np.pi))
        elif diff_p2 < -np.pi:
            diff_p2 = np.abs(np.pi - (np.abs(diff_p2) - np.pi))
            
        fft1_diff_a1.append(diff_a1)
        fft1_diff_a2.append(diff_a2)
        fft1_diff_p1.append(diff_p1)
        fft1_diff_p2.append(diff_p2)
        
        #print '\nFFT OUT OF BOX'
        #print 'Est. Amplitude 1 = ',mag1
        #print 'Est. Phase 1 = ',phase1
        #print 'Est. Amplitude 2 = ',mag2
        #print 'Est. Phase 2 = ',phase2
        #print '------------------------------------------\n'
        
        #FFT - WINDOWED
        #----------------------------------------------------------------
        ofac=1
        periods,mag,ph,fr,fi,fft_array,amp_corr = modules.fft_spectra(b,vals,ofac,hourly_step,w=True)
        
        closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
        closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
        mag1 = mag[closest_period_1]
        phase1 = ph[closest_period_1]
        mag2 = mag[closest_period_2]
        phase2 = ph[closest_period_2] 
        
        diff_a1 = (mag1/amp1)*100.
        diff_a2 = (mag2/amp2)*100.
        diff_p1 = phase_1 - phase1
        if diff_p1 > np.pi:
            diff_p1 = np.abs(-np.pi + (diff_p1 - np.pi))
        elif diff_p1 < -np.pi:
            diff_p1 = np.abs(np.pi - (np.abs(diff_p1) - np.pi)) 
        diff_p2 = phase_2 - phase2
        if diff_p2 > np.pi:
            diff_p2 = np.abs(-np.pi + (diff_p2 - np.pi))
        elif diff_p2 < -np.pi:
            diff_p2 = np.abs(np.pi - (np.abs(diff_p2) - np.pi))
            
        fft2_diff_a1.append(diff_a1)
        fft2_diff_a2.append(diff_a2)
        fft2_diff_p1.append(diff_p1)
        fft2_diff_p2.append(diff_p2)
        
        #print '\nFFT WINDOWED'
        #print 'Est. Amplitude 1 = ',mag1
        #print 'Est. Phase 1 = ',phase1
        #print 'Est. Amplitude 2 = ',mag2
        #print 'Est. Phase 2 = ',phase2
        #print '------------------------------------------\n'
        
        #FFT - WINDOWED,OVERSAMPLED
        #----------------------------------------------------------------
        ofac=4
        periods,mag,ph,fr,fi,fft_array,amp_corr = modules.fft_spectra(b,vals,ofac,hourly_step,w=True)
        
        closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
        closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
        mag1 = mag[closest_period_1]
        phase1 = ph[closest_period_1]
        mag2 = mag[closest_period_2]
        phase2 = ph[closest_period_2] 
        
        diff_a1 = (mag1/amp1)*100.
        diff_a2 = (mag2/amp2)*100.
        diff_p1 = phase_1 - phase1
        if diff_p1 > np.pi:
            diff_p1 = np.abs(-np.pi + (diff_p1 - np.pi))
        elif diff_p1 < -np.pi:
            diff_p1 = np.abs(np.pi - (np.abs(diff_p1) - np.pi)) 
        diff_p2 = phase_2 - phase2
        if diff_p2 > np.pi:
            diff_p2 = np.abs(-np.pi + (diff_p2 - np.pi))
        elif diff_p2 < -np.pi:
            diff_p2 = np.abs(np.pi - (np.abs(diff_p2) - np.pi))
            
        fft3_diff_a1.append(diff_a1)
        fft3_diff_a2.append(diff_a2)
        fft3_diff_p1.append(diff_p1)
        fft3_diff_p2.append(diff_p2)
        
        #print '\nFFT WINDOWED, OVERSAMPLED'
        #print 'Est. Amplitude 1 = ',mag1
        #print 'Est. Phase 1 = ',phase1
        #print 'Est. Amplitude 2 = ',mag2
        #print 'Est. Phase 2 = ',phase2
        #print '------------------------------------------\n'

    #append average of each monte carlo processed output
    lomb1_a1.append(np.average(lsp1_diff_a1))
    lomb1_a2.append(np.average(lsp1_diff_a2))
    lomb1_p1.append(np.average(lsp1_diff_p1))
    lomb1_p2.append(np.average(lsp1_diff_p2))

    lomb2_a1.append(np.average(lsp2_diff_a1))
    lomb2_a2.append(np.average(lsp2_diff_a2))
    lomb2_p1.append(np.average(lsp2_diff_p1))
    lomb2_p2.append(np.average(lsp2_diff_p2))

    lomb3_a1.append(np.average(lsp3_diff_a1))
    lomb3_a2.append(np.average(lsp3_diff_a2))
    lomb3_p1.append(np.average(lsp3_diff_p1))
    lomb3_p2.append(np.average(lsp3_diff_p2))

    #lomb4_a1.append(np.average(lsp4_diff_a1))
    #lomb4_a2.append(np.average(lsp4_diff_a2))
    #lomb4_p1.append(np.average(lsp4_diff_p1))
    #lomb4_p2.append(np.average(lsp4_diff_p2))

    #lomb5_a1.append(np.average(lsp5_diff_a1))
    #lomb5_a2.append(np.average(lsp5_diff_a2))
    #lomb5_p1.append(np.average(lsp5_diff_p1))
    #lomb5_p2.append(np.average(lsp5_diff_p2))

    fft1_a1.append(np.average(fft1_diff_a1))
    fft1_a2.append(np.average(fft1_diff_a2))
    fft1_p1.append(np.average(fft1_diff_p1))
    fft1_p2.append(np.average(fft1_diff_p2))

    fft2_a1.append(np.average(fft2_diff_a1))
    fft2_a2.append(np.average(fft2_diff_a2))
    fft2_p1.append(np.average(fft2_diff_p1))
    fft2_p2.append(np.average(fft2_diff_p2))

    fft3_a1.append(np.average(fft3_diff_a1))
    fft3_a2.append(np.average(fft3_diff_a2))
    fft3_p1.append(np.average(fft3_diff_p1))
    fft3_p2.append(np.average(fft3_diff_p2))

    #fft4_a1.append(np.average(fft4_diff_a1))
    #fft4_a2.append(np.average(fft4_diff_a2))
    #fft4_p1.append(np.average(fft4_diff_p1))
    #fft4_p2.append(np.average(fft4_diff_p2))

print '-----------------------------'
print np.average(lomb1_a1)
print np.average(lomb2_a1)
print np.average(lomb3_a1)   
#print np.average(lomb4_a1)
#print np.average(lomb5_a1)
print np.average(fft1_a1)
print np.average(fft2_a1)
print np.average(fft3_a1)
#print np.average(fft4_a1)
print '-----------------------------'

plt.plot(all_frac_points,lomb1_a1,label='LOMB STD',color='black',marker='o')
plt.plot(all_frac_points,lomb2_a1,label='LOMB WIN',color='blue',marker='o')
plt.plot(all_frac_points,lomb3_a1,label='LOMB WIN,SPEC',color='green',marker='o')
#plt.plot(all_frac_points,lomb4_a1,label='LOMB WIN,OVER,INTERP',color='brown',marker='<')
#plt.plot(all_frac_points,lomb5_a1,label='LOMB WIN,SPEC',color='purple',marker='^')
plt.plot(all_frac_points,fft1_a1,label='FFT STD',color='red',marker='x')
plt.plot(all_frac_points,fft2_a1,label='FFT WIN',color='yellow',marker='x')
plt.plot(all_frac_points,fft3_a1,label='FFT WIN,OVER',color='orange',marker='x')
#plt.plot(all_frac_points,fft4_a1,label='FFT WIN,OVER,INTERP',color='pink',marker='<')
plt.axhline(linewidth=1, color='r',linestyle='--')
plt.title('Amplitude diffs, period = 1')
plt.xlabel('% gaps')
plt.ylabel('% Difference')
plt.legend(loc=2)
plt.show()

print '-----------------------------'
print np.average(lomb1_a2)
print np.average(lomb2_a2)
print np.average(lomb3_a2)   
#print np.average(lomb4_a2)
#print np.average(lomb5_a2)
print np.average(fft1_a2)
print np.average(fft2_a2)
print np.average(fft3_a2)
#print np.average(fft4_a2)
print '-----------------------------'

plt.plot(all_frac_points,lomb1_a2,label='LOMB STD',color='black',marker='x')
plt.plot(all_frac_points,lomb2_a2,label='LOMB WIN',color='blue',marker='x')
plt.plot(all_frac_points,lomb3_a2,label='LOMB WIN,SPEC',color='green',marker='x')
#plt.plot(all_frac_points,lomb4_a2,label='LOMB WIN,OVER,INTERP',color='brown',marker='<')
#plt.plot(all_frac_points,lomb5_a2,label='LOMB WIN,SPEC',color='purple',marker='^')
plt.plot(all_frac_points,fft1_a2,label='FFT STD',color='red',marker='o')
plt.plot(all_frac_points,fft2_a2,label='FFT WIN',color='yellow',marker='o')
plt.plot(all_frac_points,fft3_a2,label='FFT WIN,OVER',color='orange',marker='o')
#plt.plot(all_frac_points,fft4_a2,label='FFT WIN,OVER,INTERP',color='pink',marker='<')
plt.axhline(linewidth=1, color='r',linestyle='--')
plt.title('Amplitude diffs, period = 50')
plt.xlabel('% gaps')
plt.ylabel('% Difference')
plt.legend(loc=2)  
plt.show()

print '-----------------------------'
print np.average(lomb1_p1)
print np.average(lomb2_p1)
print np.average(lomb3_p1)   
#print np.average(lomb4_p1)
#print np.average(lomb5_p1)
print np.average(fft1_p1)
print np.average(fft2_p1)
print np.average(fft3_p1)
#print np.average(fft4_p1)
print '-----------------------------'

plt.plot(all_frac_points,lomb1_p1,label='LOMB STD',color='black',marker='o')
plt.plot(all_frac_points,lomb2_p1,label='LOMB WIN',color='blue',marker='o')
plt.plot(all_frac_points,lomb3_p1,label='LOMB WIN,SPEC',color='green',marker='o')
#plt.plot(all_frac_points,lomb4_p1,label='LOMB WIN,OVER,INTERP',color='brown',marker='<')
#plt.plot(all_frac_points,lomb5_p1,label='LOMB WIN,SPEC',color='purple',marker='^')
plt.plot(all_frac_points,fft1_p1,label='FFT STD',color='red',marker='x')
plt.plot(all_frac_points,fft2_p1,label='FFT WIN',color='yellow',marker='x')
plt.plot(all_frac_points,fft3_p1,label='FFT WIN,OVER',color='orange',marker='x')
#plt.plot(all_frac_points,fft4_p1,label='FFT WIN,OVER,INTERP',color='pink',marker='<')
plt.axhline(linewidth=1, color='r',linestyle='--')
plt.title('Phase diffs, period = 1')
plt.xlabel('% gaps')
plt.ylabel('Radian diff')
plt.legend(loc=2)  
plt.show()

print '-----------------------------'
print np.average(lomb1_p2)
print np.average(lomb2_p2)
print np.average(lomb3_p2)   
#print np.average(lomb4_p2)
#print np.average(lomb5_p2)
print np.average(fft1_p2)
print np.average(fft2_p2)
print np.average(fft3_p2)
#print np.average(fft4_p2)
print '-----------------------------'

plt.plot(all_frac_points,lomb1_p2,label='LOMB STD',color='black',marker='o')
plt.plot(all_frac_points,lomb2_p2,label='LOMB WIN',color='blue',marker='o')
plt.plot(all_frac_points,lomb3_p2,label='LOMB WIN,SPEC',color='green',marker='o')
#plt.plot(all_frac_points,lomb4_p2,label='LOMB WIN,OVER,INTERP',color='brown',marker='<')
#plt.plot(all_frac_points,lomb5_p2,label='LOMB WIN,SPEC',color='purple',marker='^')
plt.plot(all_frac_points,fft1_p2,label='FFT STD',color='red',marker='x')
plt.plot(all_frac_points,fft2_p2,label='FFT WIN',color='yellow',marker='x')
plt.plot(all_frac_points,fft3_p2,label='FFT WIN,OVER',color='orange',marker='x')
#plt.plot(all_frac_points,fft4_p2,label='FFT WIN,OVER,INTERP',color='pink',marker='<')
plt.axhline(linewidth=1, color='r',linestyle='--')
plt.title('Phase diffs, period = 50')
plt.xlabel('% gaps')
plt.ylabel('Radian diff')                                                                                                                                                                                                                     
plt.legend(loc=2)
plt.show()


