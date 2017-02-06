from numpy import *
from numpy.fft import *

def spread(y, yy, n, x, m):
#    """ 
#    Given an array yy(0:n-1), extirpolate (spread) a value y into 
#    m actual array elements that best approximate the "fictional" 
#    (i.e., possible noninteger) array element number x. The weights 
#   used are coefficients of the Lagrange interpolating polynomial 
#    Arguments: 
#    y :  
#    yy :  
#    n :  
#    x :  
#    m :  
#    Returns: 
     
#    """
    nfac=[0,1,1,2,6,24,120,720,5040,40320,362880]
    if m > 10. :
        print 'factorial table too small in spread'
        return

    ix=long(x)
    if x == float(ix):
       yy[ix]=yy[ix]+y
    else:
       ilo = long(x-0.5*float(m)+1.0)
       ilo = max(ilo,1)
       ilo = min(ilo, n-m+1 )
       #print 'x = ', x
       #print 'ilo = ', ilo
       if ilo+m-1 == n:
           ilo = ilo-1
       ihi = ilo+m-1
       nden = nfac[m]
       fac=x-ilo
       #print 'fac = ', fac
       for j in range(ilo+1,ihi+1): fac = fac*(x-j)
       #print 'fac = ', fac
       yy[ihi] = yy[ihi] + y*fac/(nden*(x-ihi))
       #print yy[ihi] + y*fac/(nden*(x-ihi))
       #print 'ihi = ', ihi
       for j in range(ihi-1,ilo-1,-1):
          nden=(nden/(j+1-ilo))*(j-ihi)
          #print yy[j] + y*fac/(nden*(x-j))
          yy[j] = yy[j] + y*fac/(nden*(x-j))
          #print 'j = ', j


def fasper(x, y, ofac=4, hifac=1 , MACC=4):

    """Calc Lomb periodogram of the data in 'x' and 'y'.
    The values in 'x' are the positions of the corresponding amplitude values
    in 'y'. For best results, normalize 'x' values to fit into [0.0..1.0) and
    (since the data should be periodic anyway) copy the value at 0.0 to 1.0.
    """
    
    #if len(x) % 2 == 1:
    #    x = x[:-1]
    #    y = y[:-1]

    n = len(x)
    print 'N In = ', n
    if n != len(y):
        print 'Incompatible arrays.'
        return


  #Set min and max periods 
    pmin = x[2]-x[0]
    pmax = x[-1]-x[0]


  #Max Frequency set
    fmax = float64(1.0/pmin)
  #Min Frequency set
    fmin = float64(1.0/pmax)

    print 'Minimum Period= ', pmin
    print 'Maximum Period= ', pmax
    print 'Minimum Frequency= ', fmin
    print 'Maximum Frequency= ', fmax


    #if np.mod(n,2) == 0: 
    #    nout  = int(0.5*ofac*hifac*n)-1
    #else:
    #    nout  = int(0.5*ofac*hifac*n)
    nout  = int(0.5*ofac*hifac*n)
    print 'Nout = ', nout
    nfreqt = int(ofac*hifac*n*MACC)   #Size the FFT as next power  
    nfreq = 64             # of 2 above nfreqt.  

    while nfreq < nfreqt:
        nfreq = 2*nfreq
    nwk = nfreq

  #Compute the mean, variance  
    ave = y.mean()
  ##sample variance because the divisor is N-1  
    var = ((y-y.mean())**2).sum()/(len(y)-1)
  # and range of the data.  
    #xmin = x.min()
    xmin = min(x)
    xmax = max(x)
    xdif = xmax-xmin
    fac  = nwk/(xdif*ofac)
    fndim = nwk

    xave=0.5*(xmax+xmin)

    print 'xdif=', xdif
  #extirpolate the data into the workspaces  
    wk1 = zeros(nwk, dtype='complex')
    wk2 = zeros(nwk, dtype='complex')

    for j in xrange(n):
        ck = (x[j] - xmin) * fac % nwk
        ckk = 2.0 * ck % nwk

        spread(y[j] - ave, wk1, nwk, ck, MACC)
        spread(1.0, wk2,nwk, ckk, MACC)

    #print len(wk2)
    #Take the Fast Fourier Transforms  
    #Take the FFT's 
    print wk1[0:100]
    print len(wk1)
    W = fft(wk1);
    # % only positive frequencies, without f=0
    print W[1:11]
    print W[nout+2:nout+12]
    wk1 = W[1:nout+1]
    rwk1 = wk1.real
    iwk1 = wk1.imag

    W = fft(wk2)
    wk2 = W[1:nout+1]
    rwk2 = wk2.real
    iwk2 = wk2.imag 

    #print 'first wk1 10',wk1[0:10]
    #save = wk1
    #print ifft(wk1)[0:10]
    #wk1 = ifft( wk1 )*len(wk1)
    #wk2 = ifft( wk2 )*len(wk1)
    #wk1 = wk1[1: +1]
    #wk2 = wk2[1: +1]
    #wk1 = wk1[1:nout+1]
    #print 'wk1 cut ', wk1[0:10]
    #wk2 = wk2[1:nout+1]
    #rwk1 = wk1.real
    #iwk1 = wk1.imag
    #rwk2 = wk2.real
    #iwk2 = wk2.imag


    #print 'iwk1', iwk1[720:740]
    #print 'rwk1', rwk1[0:10]
    
    #print 'iwk1 len =', len(iwk1)
    #sumsh = np.add.reduce(rwk1*iwk1)
    #sumc = np.add.reduce((rwk1-iwk1)*(rwk1+iwk1))

    #wtau = 0.5*np.arctan2(2*sumsh,sumc)
    #print wtau	

    df  = 1.0/(xdif*ofac)
    
	
  #Compute the Lomb value for each frequency  
    hypo2 = 2.0 *  np.abs(wk2) 
    hc2wt = rwk2/hypo2
    hs2wt = iwk2/hypo2

    print 'hc2wt =', hc2wt
    print 'hs2wt =', hs2wt

    wtau = []
    #wtau2 = []
    cwt  = np.sqrt(0.5+hc2wt)
    swt  = np.sign(hs2wt)*(np.sqrt(0.5-hc2wt))
    print 'cwt', cwt[0:10]
    print 'swt', swt[0:10]
	#for i in swt:
    #for i in cwt:
    #    wtau.append(math.acos(i))
    for i in swt:
        wtau.append(math.asin(i))

    den  = 0.5*n+hc2wt*rwk2+hs2wt*iwk2
    cterm = (cwt*rwk1+swt*iwk1)**2./den
    sterm = (cwt*iwk1-swt*rwk1)**2./(n-den)
    print 'sterm ', sterm[0:10]
    print 'cterm ', cterm[0:10]

    print 'alt ss', iwk1[0:10]
    print 'ss ',(cwt[0:10]*iwk1[0:10]-swt[0:10]*rwk1[0:10])	
    samp_spacing = 1./24
    
    previous = 0
    wk1 = df*(np.arange(nout, dtype='float')+1.)
    #wk1 = np.arange(1,nout+1,1,dtype='float')/(samp_spacing*n)
	#min_freq = 2*(1./24) 
    #wk1 = np.linspace(1./xdif,12,nout)
    wk2 = (cterm+sterm)/(2*var)
    pmax = wk2.max()
    jmax = wk2.argmax()

    #Limit the frequencies between min and max frequency
    fminpos = min(range(len(wk1)), key=lambda i: abs(wk1[i]-fmin))
    fmaxpos = min(range(len(wk1)), key=lambda i: abs(wk1[i]-fmax))
    wk1=wk1[fminpos:fmaxpos]
    wk2=wk2[fminpos:fmaxpos]
    wtau = wtau[fminpos:fmaxpos]
    sterm = sterm[fminpos:fmaxpos]
    cterm = cterm[fminpos:fmaxpos]

    #calculate amplitude
    dim  = 2*(len(x)+1)
    #dim=2*(nout+1)
    fac=np.sqrt(var*dim/2)
    amp=fac*np.sqrt(wk2)
    amp = amp/len(x)
    amp = amp*2  

    #twopi = 2*np.pi
    #calculate phase
    #phLS = np.arctan2(sterm**0.5,cterm**0.5)
    #arg0 = twopi*(xave+xmin)*wk1+wtau 
    #print 'iy', np.sqrt(sterm[0:10])
    #print 'ry', np.sqrt(cterm[0:10])
    #print 'iy', iy[0:10]
    #print 'ry', ry[0:10]
    #print 'Fx ',wk1[0:10]
    #print 'Wtau ',wtau[0:10]
    #print 'Arg0 ',arg0[0:10]
	#arg0 = twopi*ave*wk1+wtau
    #ph = np.mod(phLS+arg0, twopi)
    #print 'phLS ',phLS[0:10]
    #print 'ph ',ph[0:10]
    ph=np.mod(ph +5*twopi, twopi)

    #print np.min(ph)
    #print np.max(ph)
   #print 'ph2 ',ph[0:10]
   #print 'py ',wk2[0:10]

    return wk1, wk2, amp
