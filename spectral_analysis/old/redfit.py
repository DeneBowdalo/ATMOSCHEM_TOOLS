import numpy as np
import math
import random
import scipy.special
import scipy.signal
import modules
import matplotlib.pyplot as plt
from pandas import *
import math

#----------------------------------------------------------------------
#  Manfred Mudelsee's code for tau estimation
#----------------------------------------------------------------------
# TAUEST: Routine for persistence estimation for unevenly spaced time series
#----------------------------------------------------------------------
#       Main variables
#
#       t       :       time
#       x       :       time series value
#      npo      :       number of points
#      dt       :       average spacing
#   scalt       :       scaling factor (time)
#     rho       :       in the case of equidistance, rho = autocorr. coeff.
#      ls       :       LS function
#   brent       :       Brent's search, minimum LS value
#    mult       :       flag (multiple solution)
#    amin       :       estimated value of a = exp(-scalt/tau)
#

def tauest(t, x, npo,ave,var,ierr):
    calc = 'normal'
    
    #----------------------------------------------------------------------
    #Correct time direction; assume that ages are input
    #--------------------------------------------------
    tscal = -t[::-1]
    xscal = x[::-1]
    
    #Scaling of x
    #------------
    fac = np.sqrt(var)
    xscal = xscal / fac

    # Scaling of t (=> start value of a = 1/e)
    #---------------------------------------
    dt = (tscal[npo-1]-tscal[0]) / (npo-1.)
  
    rho = rhoest(npo, xscal)
    print 'rho = ', rho
    
    if rho <= 0.0:
        rho = 0.05
        print 'Warning: rho estimation: < 0'
        ierr = 2
    elif rho > 1.0:
        rho = 0.95
        print 'Warning: rho estimation: > 1'
        ierr = 2
    scalt = -np.log(rho)/dt
    tscal = tscal * scalt
    
    #Estimation
    #----------
    amin, mult, ierr = minls(npo, tscal, xscal, ierr)
    if (ierr == 1):
        print ' Error in MINLS'
        calc = 'invalid' 
    if (mult == 1):
        print ' Estimation problem: LS function has > 1 minima'
        calc = 'invalid'
    if (amin <= 0.0):
        print ' Estimation problem: a_min =< 0'
        calc = 'invalid'
    elif (amin >= 1.0):
        print ' Estimation problem: a_min >= 1'
        calc = 'invalid'

    if calc != 'invalid':
        #determine tau
        #-------------
        tau = -1.0 /(scalt*np.log(amin))

        # determine rho, corresponding to tau
        # -----------------------------------
        rhoavg = np.exp(-dt / tau)

    else:
        tau = 'invalid'
        rhoavg = 'invalid'

    return tau, rhoavg, ierr
    
#----------------------------------------------------------------------
# Autocorrelation coefficient estimation (equidistant data).
#----------------------------------------------------------------------
def rhoest(n,x):
    sum1=0.0
    sum2=0.0
    for i in np.arange(1,n):
        sum1=sum1+x[i]*x[i-1]
        sum2=sum2+x[i]**2.0
    rho=sum1/sum2
    return rho

#----------------------------------------------------------------------
# Least-squares function
#----------------------------------------------------------------------
def ls(a,t,x,n):  
    ls=0.0
    for i in np.arange(1,n):
        ls=ls+(x[i]-x[i-1]*np.sign(a)*np.abs(a)**(t[i]-t[i-1]))**2.0
    return ls
  
#----------------------------------------------------------------------
# Minimization of least-squares function ls.
#----------------------------------------------------------------------
def minls(n, t, x, ierr):
    a_ar1 = 1./math.exp(1) # 1/e
    tol = 3.0e-8          # Brent's search, precision
    tol2 = 1.0e-6         # multiple solutions, precision
  
    nmu_=0
    dum1, a_ar11, ierr = brent(-2.0, a_ar1, +2.0, tol, t, x, n, ierr)
    dum2, a_ar12, ierr = brent( a_ar1, 0.5*(a_ar1+1.0), +2.0, tol, t, x, n, ierr)
    dum3, a_ar13, ierr = brent(-2.0, 0.5*(a_ar1-1.0),  a_ar1, tol, t, x, n, ierr)
    if (ierr == 1):
        print ' Error in MINLS (call to brent)'
        return
    if ((np.abs(a_ar12-a_ar11) > tol2) & (np.abs(a_ar12-a_ar1) > tol2)) or ((np.abs(a_ar13-a_ar11) > tol2) & (np.abs(a_ar13-a_ar1) > tol2)):
        nmu_=1
    dum4=np.min((dum1,dum2,dum3))
    if dum4 == dum2:
        amin=a_ar12
    elif dum4 == dum3:
        amin=a_ar13
    else:
        amin=a_ar11
    return amin, nmu_, ierr

#----------------------------------------------------------------------
# Numerical Recipes (modified): Brent's search in one direction.
#----------------------------------------------------------------------
def brent(ax,bx,cx,tol,xfunc,yfunc,nfunc,ierr):
    itmax=100
    cgold=.3819660
    zeps=1.e-18

    a=np.min((ax,cx))
    b=np.max((ax,cx))
    v=bx
    w=v
    x=v
    e=0.
    fx=ls(x,xfunc,yfunc,nfunc)
    fv=fx
    fw=fx
    for i in range(itmax):
        xm=0.5*(a+b)
        #print 'x=', np.abs(x)
        tol1=tol*np.abs(x)+zeps
        tol2=2.*tol1
        
        #print 'a= ',a
        #print 'b= ',b
        #print np.abs(x-xm)
        #print tol2-0.5*(b-a)
        if np.abs(x-xm) <= tol2-0.5*(b-a): 
            brent, xmin = brent_3(x,fx)
            print 'exit brent'
            return brent, xmin, ierr
        if np.abs(e) > tol1:
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.*(q-r)
            if q > 0.: 
                p=-p
            q=np.abs(q)
            etemp=e
            e=d
            if(np.abs(p) >= np.abs(0.5*q*etemp)) or (p <= q*(a-x)) or (p >= q*(b-x)):
                #print '1'
                d,e = brent_1(x,xm,e,a,b,cgold)
                u, fu, a, b, v, fv, w, fw, x, fx = brent_2(x,d,tol1,xfunc,yfunc,nfunc,fx,w,fw,v,fv,a,b)
                continue 
            d=p/q
            u=x+d
            if(u-a < tol2) or (b-u < tol2): 
                #print '2'
                d = np.abs(tol1)*np.sign(xm-x)
                u, fu, a, b, v, fv, w, fw, x, fx = brent_2(x,d,tol1,xfunc,yfunc,nfunc,fx,w,fw,v,fv,a,b) 
                continue
        #print '3'
        d,e = brent_1(x,xm,e,a,b,cgold)
        u, fu, a, b, v, fv, w, fw, x, fx = brent_2(x,d,tol1,xfunc,yfunc,nfunc,fx,w,fw,v,fv,a,b)
            
    ierr = 1
    print ' brent: exceed maximum iterations'
    brent, xmin = brent_3(x,fx)
    return brent, xmin, ierr
    
def brent_1(x,xm,e,a,b,cgold):
    if x >= xm:
        e=a-x
    else:
        e=b-x
    d=cgold*e
    return d,e
    
def brent_2(x,d,tol1,xfunc,yfunc,nfunc,fx,w,fw,v,fv,a,b):
    if np.abs(d) >= tol1:
        u=x+d
    else:
        u=x+(np.abs(tol1)*np.sign(d))
    fu=ls(u,xfunc,yfunc,nfunc)
    if fu <= fx:
        if u >= x: 
            a=x
        else:
            b=x
        v=w
        fv=fw
        w=x
        fw=fx
        x=u
        fx=fu
    else:
        if u <x:
            a=u
        else:
            b=u
        
        if (fu <= fw) or (w == x):
            v=w
            fv=fw
            w=u
            fw=fu
        elif(fu <= fv) or (v == x) or (v == w):
            v=u
            fv=fu
    
    return u, fu, a, b, v, fv, w, fw, x, fx
    
def brent_3(x,fx):
    xmin=x
    brent=fx
    return brent, xmin
    
#--------------------------------------------------------------------------
#  subroutine makear1(t, npo, tau, idum, red)
#--------------------------------------------------------------------------

def makear1(t,npo,tau):
    # set up AR(1) time series
    red = np.zeros(npo)
    red[0] = random.gauss(0, 1)
    for i in np.arange(1,npo):
        dt = t[i] - t[i-1]
        sigma = 1.0 - np.exp(-2.0 * dt / tau)
        sigma = np.sqrt(sigma)
        red[i] = np.exp(-dt/tau) * red[i-1] + sigma * random.gauss(0, 1)
    return red

def getdof(iwin):
    #c50 relates to type of window used, 0 is rectangular, 2 is hanning
    c50 = [0.500, 0.344, 0.167, 0.250, 0.096]

    c2 = 2.0 * c50[iwin] * c50[iwin]
    denom = 1.0 + c2 - c2/ 1.
    neff = 1. / denom
    dof = 2.0 * neff
    return dof
    
    
def getchi2(dof,alpha):
    tol = 1.0e-3
    itmax = 100
    ierr = 0
    
#use approximation for dof > 30 (Eq. 1.132 in Sachs (1984))
#----------------------------------------------------------
    if dof > 30.0:
        #za = -getz(alpha)   # NB: Eq. requires change of sign for percentile
        if ierr == 1:
            return
        x = 2.0 / 9.0 / dof
        chi2 = dof * (1.0 - x + za * np.sqrt(x))**3.0
    else:
        iter = 0
        lm = 0.0
        rm = 1000.0
        if alpha > 0.5:
            eps = (1.0 - alpha) * tol
        else:
            eps = alpha * tol
        while iter <= itmax:
            iter= iter + 1
            chi2 = 0.5 * (lm + rm)
            ac = 1.0 - scipy.special.gammainc(0.5*dof, 0.5*chi2)
            if np.abs(ac - alpha) <= eps:
                break
            if ac > alpha:
                lm = chi2
            else:
                rm = chi2

    if iter > itmax:
        print "Error in GETCHI2: Iter > ItMax"
        ierr = 1
        return
    
    getchi2 = chi2

    return getchi2

def rmtrend(x, y):
# determine linear trend by means of least squares and subtract
# this trend from a data set
#
# parameters:  x, y   : real arrays for data, on output y is replaced
#                       by the detrended values
#              n      : number of data pairs
#
# ref.: after numerical recipes 14.2
#
# written:       23/07/94
# modifications: 29/01/99 - allow n <> array dimension

    sx = 0.0
    sy = 0.0
    st2 = 0.0
    b = 0.0
    sx = sx + x
    sy = sy + y
    z = x - sx
    st2 = st2 + z*z
    b = b + z * y
    b = b / st2
    a = (sy - sx * b) 
    y = y - (a + b * x)
    return y

def welch_lomb(t,x,ofac,n50,SAMP_R,detrend=True):        
    ave_seg_len = 2. * len(t) / (n50 + 1)
    if ave_seg_len%1 != 0:
        i = 1
        while ave_seg_len%1 != 0:
            ave_seg_len = 2. * (len(t)-i) / (n50 + 1)
            i+=1
        
    half_ave_seg_len = ave_seg_len/2
    
    start = 0
    end = ave_seg_len

    for i in range(n50):
    
    #split data into 50% overlapping segments
        print 'WOSA segment start =', start
        print 'WOSA segment end = ',end
        twk = t[start:end]
        xwk = x[start:end]
        start = start+half_ave_seg_len
        end = end+half_ave_seg_len

        #if detrend == True:
        #    xwk = scipy.signal.detrend(xwk, type='linear')

        #take lomb
        model_periods,model_mag,model_ph,fr,fi,amp_corr = modules.take_lomb(twk,xwk,ofac,SAMP_R)
        full_n = len(model_periods)

        #sum raw spectra
        try:
            gxx = gxx + model_mag
            fr_ave = fr_ave + fr
            fi_ave = fi_ave + fi
        except:
            gxx = model_mag   
            fr_ave = fr
            fi_ave = fi

    model_freq = 1./model_periods 
    df = model_freq[1] - model_freq[0]
    gxx  = gxx/n50
    
    #trend_test = model_periods < 10000
    #corr_periods = model_periods[trend_test]
    #corr_freq = model_freq[trend_test]
    #corr_gxx = gxx[trend_test]
    #corr_fr = fr_ave[trend_test]
    #corr_fi = fi_ave[trend_test]
    
    return model_periods,gxx,full_n,fr_ave,fi_ave

def gettau(t,x,n50,avgdt):

    rhosum = 0.0

    #ave_seg_len = 2. * len(t) / (n50 + 1)
    #if ave_seg_len%1 != 0:
    #    i = 1
    #    while ave_seg_len%1 != 0:
    #        ave_seg_len = 2. * (len(t)-i) / (n50 + 1)
    #        i+=1
        
    #half_ave_seg_len = ave_seg_len/2
    
    #start = 0
    #end = ave_seg_len

    #for i in range(n50):
    ##split data into 50% overlapping segments
    #print 'WOSA Segment %i'%(i) 
    #print 'WOSA segment start =', start
    #print 'WOSA segment end = ',end
    #twk = t[start:end]
    #xwk = x[start:end]
    #start = start+half_ave_seg_len
    #end = end+half_ave_seg_len

    #xwk = scipy.signal.detrend(x, type='linear')

    tau, rho, ierr = tauest(t, x, len(x),np.average(x),np.var(x),0)

    if tau == 'invalid':
        return tau, rho, ierr
            
    #bias correction for rho (Kendall & Stuart, 1967; Vol. 3))
    rho = (rho * (len(x) - 1.0) + 1.0) / (len(x) - 4.0)
    rhosum = rhosum + rho
  
    #average rho
    rho = rhosum / float(n50)
    
    #average tau
    tau = -avgdt / np.log(rho)
    
    return tau, rho, ierr

def significance_tests(x,red_periods,red_mag,nsim,n50,corr,mctest):
    ci80 = np.zeros(len(red_periods))
    ci90 = np.zeros(len(red_periods))
    ci95 = np.zeros(len(red_periods))
    ci99 = np.zeros(len(red_periods))

    #red-noise false-alarm levels from percentiles of MC simulation
    #--------------------------------------------------------------
    if mctest == True: 
        for i in range(len(red_periods)):
            red_mag[:nsim,i] = np.sort(red_mag[:nsim, i])
    
        #set percentil indices
        #---------------------
        idx80 = int(0.80 * nsim)
        idx90 = int(0.90 * nsim)
        idx95 = int(0.95 * nsim)
        idx99 = int(0.99 * nsim)
    
        #find frequency-dependent percentil and apply bias correction
        #------------------------------------------------------------
        for i in range(len(red_periods)):
            ci80[i] = red_mag[idx80, i] / corr[i]
            ci90[i] = red_mag[idx90, i] / corr[i]
            ci95[i] = red_mag[idx95, i] / corr[i]
            ci99[i] = red_mag[idx99, i] / corr[i]
     

    #get degrees of freedom
    dof = getdof(0) # 0 is rectangular,2 is hanning window

    #set significance grid
    sig_levels = np.arange(95,99.9,0.1)
    sig_levels = np.append(sig_levels,np.arange(99.9,99.99,0.01))
    sig_levels = np.append(sig_levels,np.arange(99.99,99.999,0.001))
    sig_levels = np.append(sig_levels,np.arange(99.999,99.9999,0.0001))
    sig_levels = np.append(sig_levels,np.arange(99.9999,99.99999,0.00001))
    
    grid_range = 1.-(sig_levels/100.)
    fac_grid = np.empty(len(sig_levels))

    #get scaling factors for red noise model
    fac80 = getchi2(dof, 0.20) / dof
    fac85 = getchi2(dof, 0.15) / dof
    fac90 = getchi2(dof, 0.10) / dof
    fac95 = getchi2(dof, 0.05) / dof
    fac99 = getchi2(dof, 0.01) / dof
    fac99_9 = getchi2(dof, 0.001) / dof
    fac99_99 = getchi2(dof, 0.0001) / dof
    fac99_999 = getchi2(dof, 0.00001) / dof

    for i in range(len(grid_range)):
        fac = getchi2(dof, grid_range[i]) / dof
        fac_grid[i] = fac

    # critical false alarm level after Thomson (1990)
    # -----------------------------------------------
    #split = len(x)/n50
    #half_split = split/2.
    #ave_seg_len = split+half_split
    ave_seg_len = len(x)
    alphacrit = 1.0 / ave_seg_len
    faccrit = getchi2(dof, alphacrit) / dof

    return fac95,fac99,fac99_9,faccrit,fac_grid,sig_levels,ci99


def red_background(nsim,mctest,t,x,ofac,all_t,all_x ):
    
    #average dt of entire time series
    diffs = [t[i+1]-t[i] for i in range(len(t)-1)]  
    avgdt = np.average(diffs)

    ave = np.mean(x)
    #subtract mean from data
    x = x - ave

    #GET TAU
    tau, rhoavg, ierr = gettau(t,x,1,avgdt)

    nout = int(0.5*int(ofac)*1*len(t))

    if tau == 'invalid':
        model_periods,model_mag,corr_model_mag,model_fr,model_fi,red_periods,red_mag_avg,gredth,fac95,fac99,fac99_9,faccrit,fac_grid,sig_levels,tau,corr = np.zeros(nout),np.zeros(nout),np.zeros(nout),np.zeros(nout),np.zeros(nout),np.zeros(nout),np.zeros(nout),np.zeros(nout),0,0,0,0,np.zeros(nout),np.zeros(nout),0,1
        return model_periods,model_mag,corr_model_mag,model_fr,model_fi,red_periods,red_mag_avg,gredth,fac95,fac99,fac99_9,faccrit,fac_grid,sig_levels,tau,corr

    #make sure that tau is non-negative
    if (tau < 0.0):
        print 'Negative tau is forced to zero.'
        tau = 0.0
        
    x = x + ave
    
    #determine lag-1 autocorrelation coefficient
    rho = np.exp(-avgdt / tau)    # avg. autocorrelation coefficient
    rhosq = rho * rho
    
    #t = np.copy(all_t)
    #x = np.copy(all_x)
    
    xdif = np.max(t)-np.min(t)
    
    #Calculate model spectrum
    
    
    model_periods,model_mag,model_ph,model_fr,model_fi = modules.take_lomb_unwindowed(t,x,ofac,avgdt)
    fft_periods_all,fft_mag_all,fft_fr,fft_fi,fft_array = modules.take_fft_unwindowed(t,x,avgdt)
    fft_freq_all = 1./fft_periods_all
    
    
    full_n = len(model_periods)
    model_freq = 1./model_periods

    # estimate data variance from data spectrum
    # ----------------------------------------
    varx = (model_freq[0]) * np.sum(model_mag)  # NB: freq[1] = df
    
    varx_fft = (fft_freq_all[0]) * np.sum(fft_mag_all)
    
    red_mag = np.zeros((nsim,len(model_periods)))
    red_mag_sum = np.zeros(len(model_periods))
    fft_mag = np.zeros((nsim,len(fft_periods_all)))
    fft_mag_sum = np.zeros(len(fft_periods_all))

    #create AR1 spectrum nsim times
    for i in range(nsim):
        print 'Nsim = ', i+1
        red = makear1(t,len(x),tau)
        if mctest == True:
            #red_periods,red_mag[i,:],red_ph,red_fr,red_fi = modules.take_lomb_unwindowed(t,red,ofac,avgdt)
            #red_freq = 1./red_periods
            
            fft_periods,fft_mag[i,:],fft_fr,fft_fi,fft_array = modules.take_fft_unwindowed(t,red,avgdt)
            fft_freq = 1./fft_periods
            
        else:
            #red_periods,red_mag[0,:],red_ph,red_fr,red_fi = modules.take_lomb_unwindowed(t,red,ofac,avgdt)
            #red_freq = 1./red_periods
            
            fft_periods,fft_mag[0,:],fft_fr,fft_fi,fft_array = modules.take_fft_unwindowed(t,red,avgdt)
            fft_freq = 1./fft_periods
            
            #plt.loglog(red_periods,red_mag[0,:],color='black')
            #plt.loglog(fft_periods,fft_mag[0,:],color='red')
            
            #plt.show()

        #red_periods = red_periods[cp:]
        #red_mag = red_mag[0,cp:]

    #scale and sum red-noise spectra
    #-------------------------------
        if mctest == True:    
            #varr = (red_freq[0]) * np.sum(red_mag[i,:])  # NB: freq[1] = df
            #fac = varx / varr
            #red_mag[i,:] = fac * red_mag[i,:]
            #red_mag_sum = red_mag_sum + red_mag[i,:]
            
            varr = (fft_freq[0]) * np.sum(fft_mag[i,:])  # NB: freq[1] = df
            fac = varx_fft / varr
            fft_mag[i,:] = fac * fft_mag[i,:]
            fft_mag_sum = fft_mag_sum + fft_mag[i,:]
            
            
        else:
            #varr = (red_freq[0]) * np.sum(red_mag[0,:])  # NB: freq[1] = df
            #fac = varx / varr
            #red_mag_sum = red_mag_sum + fac * red_mag[0,:]
            
            varr = (fft_freq[0]) * np.sum(fft_mag[0,:])  # NB: freq[1] = df
            fac = varx_fft / varr
            fft_mag_sum = fft_mag_sum + fac * fft_mag[0,:]
        
    #determine average red-noise spectrum; scale average again to
    #make sure that roundoff errors do not affect the scaling
    #------------------------------------------------------------
    #red_mag_avg = red_mag_sum / float(nsim)
    #varr = (red_freq[0]) * np.sum(red_mag_avg)
    #fac = varx / varr
    #red_mag_avg = fac * red_mag_avg

    fft_mag_avg = fft_mag_sum / float(nsim)
    varr = (fft_freq[0]) * np.sum(fft_mag_avg)
    fac = varx_fft / varr
    fft_mag_avg = fac * fft_mag_avg

    # set theoretical spectrum (e.g., Mann and Lees, 1996, Eq. 4)
    # make area equal to that of the input time series
    # -----------------------------------------------------------
    #print red_freq[-1]
    #fnyq = red_freq[-1]     
    #gredth = (1.0-rhosq) / (1.0-2.0*rho*np.cos(np.pi*red_freq/fnyq)+rhosq)
    #varr = red_freq[0] * np.sum(gredth)
    #fac = varx / varr
    #gredth = fac * gredth
    
    
    print fft_freq[-1] 
    fnyq = fft_freq[-1]     
    gredth_fft = (1.0-rhosq) / (1.0-2.0*rho*np.cos(np.pi*fft_freq/fnyq)+rhosq)
    varr = fft_freq[0] * np.sum(gredth_fft)
    fac = varx_fft / varr
    gredth_fft = fac * gredth_fft
    
    #ratio = float(len(act_periods))/float(len(model_periods))
    #print 'ratio = ', ratio
    
    #model_mag = model_mag/ratio
    #red_mag_avg = red_mag_avg/ratio
    #gredth = gredth/ratio

    
    #determine correction factor
    #---------------------------
    #corr = red_mag_avg / gredth

    #correct for bias in autospectrum
    #--------------------------------
    #corr_model_mag = model_mag / corr
    #corr_model_mag = model_mag
    #gredth = gredth*corr
    
    
    #print 'model freq 0 = ',  model_freq[1]
    print 'varx = ', varx
    print 'ofac = ', ofac
    print 'avgdt = ', avgdt
    print 'Avg. autocorr. coeff., rho = ', rho 
    print 'Avg. tau = ', tau
 
    
    npoints = len(t)                           # of data points
    nseg = int(2 * npoints / (1 + 1))         # points per segment
    #avgdt = (t[-1] - t[0]) / (npoints-1.0)    # average sampling interval
    tp = avgdt * nseg                      # average period of a segment
    #df = 1.0 / (ofac * tp)                 # freq. spacing
    df = model_freq[0]
    wz = 2.0 * np.pi * df                     # omega = 2*pi*f
    fnyq = 1.0 * 1.0 / (2.0 * avgdt)     # average Nyquist freq.
    nfreq = fnyq / df + 1                  # f(1) = f0; f(nfreq) = fNyq
    lfreq = nfreq * 2
    nout = nfreq
    
    scal = 2.0 / (1.0 * nseg * df * ofac)
    
    print 'df = ', df
    print 'Nseg = ', nseg
    print 'Scal = ', scal
    
    print model_freq[0:20]
    print model_freq[-10:]
    
    #fft_mag_all = fft_mag_all*scal
    #fft_mag_all = fft_mag_all/4.
    
    #winbw = (model_freq[1]-model_freq[0]) * ofac * 1.21
    
    #model_mag = model_mag*len(t)
    
    #model_mag = 20*np.log10(model_mag)
    #gredth = 20*np.log10(gredth)

    #plt.loglog(red_periods,red_mag_avg)
    #plt.loglog(fft_periods,fft_mag_avg)
    #plt.loglog(fft_periods_all,fft_mag_avg)
    #plt.loglog(fft_periods,gredth_fft)
    #plt.plot(model_freq,model_mag)
    #plt.plot(red_freq,red_mag_avg)
    plt.loglog(fft_periods_all,fft_mag_all)
    plt.loglog(fft_periods,fft_mag_avg)
    #plt.plot(red_freq,gredth)
    #plt.loglog(fft_periods,gredth_fft)
    #plt.plot(model_freq,corr_model_mag)
    
    #plt.show()
    
    corr = [1.]*len(fft_freq)
    fac95,fac99,fac99_9,faccrit,fac_grid,sig_levels,ci99 = significance_tests(x,fft_periods,fft_mag,nsim,1,corr,mctest)

    plt.loglog(fft_periods,fft_mag_avg*fac95)
    plt.loglog(fft_periods,fft_mag_avg*fac99)
    plt.loglog(fft_periods,fft_mag_avg*faccrit)
    plt.show()
    
    return model_periods,model_mag,corr_model_mag,model_fr,model_fi,red_periods,red_mag_avg,gredth,fac95,fac99,fac99_9,faccrit,fac_grid,sig_levels,tau,corr

def red_ar1_fit(nsim,t,x,ofac):
    
    #average dt of entire time series
    diffs = [t[i+1]-t[i] for i in range(len(t)-1)]  
    avgdt = np.average(diffs)

    ave = np.mean(x)
    #subtract mean from data
    x = x - ave

    #GET TAU
    tau, rhoavg, ierr = gettau(t,x,1,avgdt)

    nout = int(0.5*int(ofac)*1*len(t))

    if tau == 'invalid':
        model_periods,model_mag,corr_model_mag,model_fr,model_fi,red_periods,red_mag_avg,gredth,fac95,fac99,fac99_9,faccrit,fac_grid,sig_levels,tau,corr = np.zeros(nout),np.zeros(nout),np.zeros(nout),np.zeros(nout),np.zeros(nout),np.zeros(nout),np.zeros(nout),np.zeros(nout),0,0,0,0,np.zeros(nout),np.zeros(nout),0,1
        return model_periods,model_mag,corr_model_mag,model_fr,model_fi,red_periods,red_mag_avg,gredth,fac95,fac99,fac99_9,faccrit,fac_grid,sig_levels,tau,corr

    #make sure that tau is non-negative
    if (tau < 0.0):
        print 'Negative tau is forced to zero.'
        tau = 0.0
        
    x = x + ave
    
    #determine lag-1 autocorrelation coefficient
    rho = np.exp(-avgdt / tau)    # avg. autocorrelation coefficient
    rhosq = rho * rho
    
    xdif = np.max(t)-np.min(t)
    
    #Calculate model spectrum
    #model_periods,model_mag,model_ph,model_fr,model_fi = modules.take_lomb_unwindowed(t,x,ofac)
    fft_periods_all,fft_mag_all,fft_fr,fft_fi,fft_array = modules.take_fft_unwindowed(t,x,avgdt)
    fft_freq_all = 1./fft_periods_all
    
    #full_n = len(model_periods)
    #model_freq = 1./model_periods

    # estimate data variance from data spectrum
    # ----------------------------------------
    #varx = (model_freq[0]) * np.sum(model_mag)  # NB: freq[1] = df
    
    varx_fft = (fft_freq_all[0]) * np.sum(fft_mag_all)
    
    #red_mag = np.zeros((nsim,len(model_periods)))
    #red_mag_sum = np.zeros(len(model_periods))
    fft_mag = np.zeros((nsim,len(fft_periods_all)))
    fft_mag_sum = np.zeros(len(fft_periods_all))

    #create AR1 spectrum nsim times
    for i in range(nsim):
        print 'Nsim = ', i+1
        red = makear1(t,len(x),tau)
        #red_periods,red_mag[i,:],red_ph,red_fr,red_fi = modules.take_lomb_unwindowed(t,red,ofac,avgdt)
        #red_freq = 1./red_periods
            
        fft_periods,fft_mag[i,:],fft_fr,fft_fi,fft_array = modules.take_fft_unwindowed(t,red,avgdt)
        fft_freq = 1./fft_periods
            

    #scale and sum red-noise spectra
    #-------------------------------    
        #varr = (red_freq[0]) * np.sum(red_mag[i,:])  # NB: freq[1] = df
        #fac = varx / varr
        #red_mag[i,:] = fac * red_mag[i,:]
        #red_mag_sum = red_mag_sum + red_mag[i,:]
            
        varr = (fft_freq[0]) * np.sum(fft_mag[i,:])  # NB: freq[1] = df
        fac = varx_fft / varr
        fft_mag[i,:] = fac * fft_mag[i,:]
        fft_mag_sum = fft_mag_sum + fft_mag[i,:]
            
    #determine average red-noise spectrum; scale average again to
    #make sure that roundoff errors do not affect the scaling
    #------------------------------------------------------------
    #red_mag_avg = red_mag_sum / float(nsim)
    #varr = (red_freq[0]) * np.sum(red_mag_avg)
    #fac = varx / varr
    #red_mag_avg = fac * red_mag_avg

    fft_mag_avg = fft_mag_sum / float(nsim)
    varr = (fft_freq[0]) * np.sum(fft_mag_avg)
    fac = varx_fft / varr
    fft_mag_avg = fac * fft_mag_avg


    #plt.plot(model_freq,model_mag)
    #plt.plot(red_freq,red_mag_avg)
    #plt.loglog(fft_periods_all,fft_mag_all)
    #plt.loglog(fft_periods,fft_mag_avg)
    
    corr = [1.]*len(fft_freq)
    fac95,fac99,fac99_9,faccrit,fac_grid,sig_levels,ci99 = significance_tests(x,fft_periods,fft_mag,nsim,1,corr,True)

    #plt.loglog(fft_periods,fft_mag_avg*fac95)
    #plt.loglog(fft_periods,fft_mag_avg*fac99)
    #plt.loglog(fft_periods,fft_mag_avg*faccrit)
    #plt.show()
    
    model_periods = np.copy(fft_periods_all)
    model_mag = np.copy(fft_mag_all)
    red_periods = np.copy(fft_periods)
    red_mag_avg = np.copy(fft_mag_avg)
    
    return model_periods,model_mag,red_periods,red_mag_avg,fac95,fac99,fac99_9,faccrit,fac_grid,sig_levels,tau,corr

def sidelobe_peak_remove(fb,fr,fi,closest_inds,crit_percent,periods):    
    peak_alter_inds = []
    
    for closest in closest_inds:
        print 'period =',periods[closest]
        print 'ind = ',closest
    
        critical_val = (fb[closest]/100.) * crit_percent # certain % of top peak
        
        try:
            up_peak_val = fb[closest+1]
        except:
            up_peak_val = 'na'
        try:
            next_up_peak_val = fb[closest+2]
        except:
            next_up_peak_val = 'na'
        try:
            down_peak_val = fb[closest-1]
        except:
            down_peak_val = 'na'
        try:
            next_down_peak_val = fb[closest-2]
        except:
            next_down_peak_val = 'na'
        
        down_inds = [closest]
        up_inds = [closest]
        
        i=1
        while (up_peak_val >= critical_val) & (next_up_peak_val < up_peak_val) & (up_peak_val != 'na') & (next_up_peak_val != 'na') :
            up_inds.append(closest+i)    
            try:
                i+=1
                up_peak_val = fb[closest+i]
                next_up_peak_val = fb[closest+(i+1)]
            except:
                next_up_peak_val = 'na'
                if up_peak_val >= critical_val:
                    up_inds.append(closest+i) 
        
        i=1
        while (down_peak_val >= critical_val) & (next_down_peak_val < down_peak_val) & (down_peak_val != 'na') & (next_down_peak_val != 'na') :
            down_inds.append(closest-i)
            try:
                i+=1
                down_peak_val = fb[closest-i]
                next_down_peak_val = fb[closest-(i+1)]
            except:
                next_down_peak_val = 'na'
                if down_peak_val >= critical_val:
                    down_inds.append(closest-i)
                
            
        all_inds = np.concatenate((down_inds,up_inds))
        all_inds = np.sort(all_inds)
        #all_inds = [int(i) for i in all_inds]
        all_inds = set(all_inds)
        all_inds = [i for i in all_inds]
        
        print all_inds
        
        
        fb[all_inds] = np.nan
        fr[all_inds] = np.nan
        fi[all_inds] = np.nan
            
        #peak_alter_inds = np.append(peak_alter_inds,closest)
    
    
    not_nan = np.logical_not(np.isnan(fb))
    indices = np.arange(len(fb))
    
    fb = np.interp(indices, indices[not_nan], fb[not_nan])
    fr = np.interp(indices, indices[not_nan], fr[not_nan])
    fi = np.interp(indices, indices[not_nan], fi[not_nan])
    
    #fb = Series(fb).interpolate().values
    #fr = Series(fr).interpolate().values
    #fi = Series(fi).interpolate().values
          
    #peak_alter_inds = peak_alter_inds.astype(int)
    
    return fb,fr,fi
    
def sidelobe_percent_remove(fb,fr,fi,closest_inds,percent_rm,periods):    

    array_inds = range(len(fb))
    del_inds = []
    
    for closest in closest_inds:
        print 'period =',periods[closest]
        print 'ind = ',closest
    
        crit_val = (periods[closest]/100.)*percent_rm
        lower_period = periods[closest]-crit_val
        upper_period = periods[closest]+crit_val
        
        closest_lower_index = min(range(len(periods)), key=lambda i: abs(periods[i]-lower_period))
        closest_upper_index = min(range(len(periods)), key=lambda i: abs(periods[i]-upper_period))      
        
        all_inds = range(closest_upper_index,closest_lower_index+1)
        
        all_inds = [i for i in all_inds if i in array_inds]
        
        del_inds.append(all_inds)
        
        fb[all_inds] = np.nan
        fr[all_inds] = np.nan
        fi[all_inds] = np.nan
    
    not_nan = np.logical_not(np.isnan(fb))
    indices = np.arange(len(fb))
    
    fb = np.interp(indices, indices[not_nan], fb[not_nan])
    fr = np.interp(indices, indices[not_nan], fr[not_nan])
    fi = np.interp(indices, indices[not_nan], fi[not_nan])
    
    del_inds = [item for sublist in del_inds for item in sublist]
    
    return fb,fr,fi,del_inds
    
def sidelobe_n_remove(fb,fr,fi,closest_inds,n_remove,periods):    

    array_inds = range(len(fb))
    del_inds = []
    
    for closest in closest_inds:
        print 'period =',periods[closest]
        print 'ind = ',closest

        #if periods[closest] >= 50:
        #    n_remove = 1
        #if periods[closest] < 50:
        #    n_remove = 0

        closest_lower_index = closest-n_remove
        closest_upper_index = closest+n_remove
        
        
        #closest_lower_index = min(range(len(periods)), key=lambda i: abs(periods[i]-lower_period))
        #closest_upper_index = min(range(len(periods)), key=lambda i: abs(periods[i]-upper_period))      
        
        print closest_lower_index
        print closest_upper_index
        
        all_inds = range(closest_lower_index,closest_upper_index+1)
        
        all_inds = [i for i in all_inds if i in array_inds]
        del_inds.append(all_inds)
        
        fb[all_inds] = np.nan
        fr[all_inds] = np.nan
        fi[all_inds] = np.nan
    
    not_nan = np.logical_not(np.isnan(fb))
    indices = np.arange(len(fb))
    
    fb = np.interp(indices, indices[not_nan], fb[not_nan])
    fr = np.interp(indices, indices[not_nan], fr[not_nan])
    fi = np.interp(indices, indices[not_nan], fi[not_nan])
    
    del_inds = [item for sublist in del_inds for item in sublist]
    
    print del_inds
    
    return fb,fr,fi,del_inds
    
