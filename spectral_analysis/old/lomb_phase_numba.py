# Lomb-Scargle periodogram (work version for reconstruction)  
# procedure is based on the Numerical Recipe's programs 
# period.f (please see there for comments and explanations; 
# Numerical Recipes, 2nd edition, Chapter 13.8, by Press et al., Cambridge, 1992)
# Here, the program code is adopted from the related IDL program lnp.pro
# and is translated to Matlab. New features are the 
# phase determination (Hocke, Ann. Geophys. 16, 356-358,1998) 
# and the output of a complex Fourier spectrum F. 
# This spectrum can be used for inverse FFT and reconstruction of an evenly 
# spaced time series (Scargle, 1989).
#    
# ATTENTION: 
# -> Because of the long story of program development and some open problems  
# -> of phase definition and construction of the FFT spectrum, 
# -> the program must be regarded as a working and discussion version 
# -> without any warranty!  
# -> Particularly the phase determination with the Lomb-Scargle
# -> periodogram has been done in a heuristic manner. Switching between 
# -> the phase reference systems may introduce errors which I am not aware yet.   
# -> Scargle (1989) gives more informations on the details  of the problem. 
#    (K. Hocke, Nov. 2007). 
#
#  input:
#  x: e.g. time vector 
#  y: observational data y(x)
#  ofac: optional (oversampling factor , integer, default is 4)
#
#  output:
#  wk1: frequency axis ( a vector of increasing linear frequencies)
#  wk2: Lomb normalized power as function of wk1-vector
#  ph:  phase vector as function of wk1-vector. The phase 
#       is in radians and is defined to be the argument 
#       of the cosine wave at the time x=0 ! 
#  vari: sigma^2,  variance of y, necessary to derive 
#       the amplitude from the normalized power wk2 
#  F:   complex Pseudo-Fourier spectrum 

#  please check the phases and their signs before interpreting the phases! 
#
#  keywords:
#  ofac: oversampling factor , integer  
#        The default is 4
#  hifac: integer, 1 for frequencies up to the Nyquist frequency 
#         (2 for 2*Nyquist frequency)
# 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

import numpy as np
#from numba import autojit, jit, float_

#@autojit
def lomb(X,Y,NOUT,OFAC,SAMP_R):
      PX = np.zeros(NOUT)
      PY = np.zeros(NOUT)
      PH = np.zeros(NOUT)
      PH1 = np.zeros(NOUT)
     
      N = len(X)
     
      TWOPI = np.pi*2
     
      XSTART = X[0]
 
      XMAX=np.max(X) 
      XMIN=np.min(X)
      XDIF=XMAX-XMIN
      
      #Create paramaters to limit periodogram, between max_frequency(min_period) = 1./2 x sampling rate, due to nyquist theorem
      #And Min_frequency(max_period) = 1./range of x input array
      P_MIN = 2.*SAMP_R
      P_MAX = XDIF
      F_MAX = 1./P_MIN
      F_MIN = 1./P_MAX

      XAVE=0.5*(XMAX+XMIN)
      PYMAX=0
      
      START_SWITCH = 'Y'
      END_SWITCH = 'Y'
      START_I = 0
      END_I = NOUT

      PNOW=1./(XDIF*OFAC)
      ARG=TWOPI*((X-XAVE)*PNOW)
      WPR=-2.*np.sin(0.5*ARG)**2
      WPI=np.sin(ARG)
      WR=np.cos(ARG)
      WI=np.sin(ARG)
      AVE=np.sum(Y)/len(Y)
      YY=(Y-AVE)
      VARI=np.var(Y)
    
      for NUM in range(NOUT):
          PX[NUM]=PNOW
          if (PNOW >= F_MIN) & (START_SWITCH == 'Y'):
              START_I = NUM
              START_SWITCH = 'N'

          if ((PNOW >= F_MAX) & (END_SWITCH == 'Y')):
              END_I = NUM-1
              END_SWITCH = 'N'
              
          SUMSH = np.add.reduce(WR*WI)
          SUMC = np.add.reduce((WR-WI)*(WR+WI))
   
          WTAU=0.5*np.arctan2(2.*SUMSH,SUMC)
          SWTAU=np.sin(WTAU)
          CWTAU=np.cos(WTAU)

          SS=WI*CWTAU-WR*SWTAU
          CC=WR*CWTAU+WI*SWTAU

          SUMS = np.add.reduce((SS**2))
          SUMC = np.add.reduce((CC**2))
          SUMSY = np.add.reduce((YY*SS))
          SUMCY = np.add.reduce((YY*CC))
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
          IY=SUMSY/(SUMS**.5)  
          RY=SUMCY/(SUMC**.5) 
          PY[NUM]=0.5*(RY**2+IY**2)/VARI

          PHLS=np.arctan2(IY,RY)            
          ARG0=TWOPI*(XAVE+XSTART)*PNOW +WTAU  
          ARG1=TWOPI*XAVE*PNOW +WTAU   
          PH[NUM]=np.mod((PHLS+ARG0), TWOPI)  
          PH1[NUM]=np.mod((PHLS+ARG1), TWOPI)   
          PNOW=PNOW+1./(OFAC*XDIF)   
      
      DIMI=2.*(N+1.)      
      FAC=np.sqrt(VARI*DIMI/2.)
      A=FAC*np.sqrt(PY) 
      CORR_A = A/N
      CORR_A = CORR_A*2.
      PH=np.mod(PH+5.*TWOPI, TWOPI)    # for value range 0,..., 2 pi
      PH1=np.mod(PH1+5.*TWOPI, TWOPI)
      FX=A*np.cos(PH1) 
      FY=A*np.sin(PH1)

      return PX,PY,CORR_A,PH1,START_I,END_I,FX,FY
