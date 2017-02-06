* Lomb-Scargle periodogram (work version for reconstruction)  
* procedure is based on the Numerical Recipe's programs 
* period.f (please see there for comments and explanations; 
* Numerical Recipes, 2nd edition, Chapter 13.8, by Press et al., Cambridge, 1992)
* Here, the program code is adopted from the related IDL program lnp.pro
* and is translated to Matlab. New features are the 
* phase determination (Hocke, Ann. Geophys. 16, 356-358,1998) 
* and the output of a complex Fourier spectrum F. 
* This spectrum can be used for inverse FFT and reconstruction of an evenly 
* spaced time series (Scargle, 1989).
*    
* ATTENTION: 
* -> Because of the long story of program development and some open problems  
* -> of phase definition and construction of the FFT spectrum, 
* -> the program must be regarded as a working and discussion version 
* -> without any warranty!  
* -> Particularly the phase determination with the Lomb-Scargle
* -> periodogram has been done in a heuristic manner. Switching between 
* -> the phase reference systems may introduce errors which I am not aware yet.   
* -> Scargle (1989) gives more informations on the details  of the problem. 
*    (K. Hocke, Nov. 2007). 
*
*  input:
*  x: e.g. time vector 
*  y: observational data y(x)
*  ofac: optional (oversampling factor , integer, default is 4)
*
*  output:
*  wk1: frequency axis ( a vector of increasing linear frequencies)
*  wk2: Lomb normalized power as function of wk1-vector
*  ph:  phase vector as function of wk1-vector. The phase 
*       is in radians and is defined to be the argument 
*       of the cosine wave at the time x=0 ! 
*  vari: sigma^2,  variance of y, necessary to derive 
*       the amplitude from the normalized power wk2 
*  F:   complex Pseudo-Fourier spectrum 

*  please check the phases and their signs before interpreting the phases! 
*
*  keywords:
*  ofac: oversampling factor , integer  
*        The default is 4
*  hifac: integer, 1 for frequencies up to the Nyquist frequency 
*         (2 for 2*Nyquist frequency)
* 
* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      subroutine lomb(X,Y,PX,PY,CORR_A,PH1,FX,FY,NIN,NOUT)
      IMPLICIT NONE
      
Cf2py intent(in) X(NIN)
Cf2py intent(in) Y(NIN)
Cf2py intent(in) PX(NOUT)
Cf2py intent(out) PY(NOUT)
Cf2py intent(out) CORR_A(NOUT)
Cf2py intent(out) PH1(NOUT)
Cf2py intent(out) FX(NOUT)
Cf2py intent(out) FY(NOUT)
Cf2py intent(hide) NIN
Cf2py intent(hide) NOUT

*  Input arguments
*    
      INTEGER NIN, NOUT, OFAC
      DOUBLE PRECISION X(NIN), Y(NIN), PX(NOUT)
      DOUBLE PRECISION PY(NOUT), CORR_A(NOUT), PH1(NOUT)
      DOUBLE PRECISION FX(NOUT),FY(NOUT)
      
*  Local variables
* 
      INTEGER HIFAC, PYMAX, DIMI, I, NUM
      DOUBLE PRECISION WI(NIN), WPI(NIN), WPR(NIN), WR(NIN), WTEMP(NIN)
      DOUBLE PRECISION SS(NIN), CC(NIN), YY(NIN), ARG(NIN)
      DOUBLE PRECISION PH(NOUT),PH2(NOUT),DIFF(NOUT),A(NOUT)
      DOUBLE PRECISION FFT_PHASE(NOUT),EXPY(NOUT)
      DOUBLE PRECISION SUMCY, SUMSY, IY, RY, PHLS, AVE, VARI
      DOUBLE PRECISION NMAX, XMAX, XMIN, XDIF
      DOUBLE PRECISION XAVE, PNOW, SUMSH, FAC, ARG0, ARG1, WTAU
      DOUBLE PRECISION CWTAU, SWTAU, SUMS, SUMC
      DOUBLE PRECISION XSTART, P
      REAL*8, PARAMETER :: TWOPI = 3.141592653589793238462643*2

      PRINT *,'N = ', NIN
      PRINT *,'NOUT =', NOUT
     
      XSTART = X(1)
      XMAX=MAXVAL(X) 
      XMIN=MINVAL(X)
      XDIF=XMAX-XMIN
      
      PRINT *, 'XDIF = ', XDIF
      XAVE=0.5*(XMAX+XMIN)
      PYMAX=0

      !PNOW=1./(XDIF*OFAC)
      PNOW = PX(1)
      ARG=TWOPI*((X-XAVE)*PNOW)
      WPR=-2.*SIN(0.5*ARG)**2
      WPI=SIN(ARG)
      WR=COS(ARG)
      WI=SIN(ARG)
      AVE=SUM(Y)/SIZE(Y)
      YY=(Y-AVE)
      VARI=SUM((Y-AVE)**2)/SIZE(Y)
    
      DO NUM = 1, NOUT
          PNOW=PX(NUM)

          SUMSH = SUM(WR*WI)
          SUMC = SUM((WR-WI)*(WR+WI))
   
          WTAU=0.5*ATAN2(2.*SUMSH,SUMC)
          SWTAU=SIN(WTAU)
          CWTAU=COS(WTAU)

          SS=WI*CWTAU-WR*SWTAU
          CC=WR*CWTAU+WI*SWTAU

          SUMS = SUM((SS**2))
          SUMC = SUM((CC**2))
          SUMSY = SUM((YY*SS))
          SUMCY = SUM((YY*CC))
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
          IY=SUMSY/(SUMS**.5)  
          RY=SUMCY/(SUMC**.5) 
          PY(NUM)=0.5*(RY**2+IY**2)/VARI

          PHLS=ATAN2(IY,RY)            
          ARG0=TWOPI*(XAVE+XSTART)*PNOW +WTAU  
          ARG1=TWOPI*XAVE*PNOW +WTAU   
          PH(NUM)=MOD((PHLS+ARG0), TWOPI)  
          PH1(NUM)=MOD((PHLS+ARG1), TWOPI)   
          !PNOW=PNOW+1./(OFAC*XDIF)   
      END DO
      
      DIMI=2.*(NIN+1.)      
      FAC=SQRT(VARI*DIMI/2.)
      A=FAC*SQRT(PY) 
      CORR_A = A/NIN
      CORR_A = CORR_A*2.
      PH=MOD(PH+5.*TWOPI, TWOPI)    
      PH1=MOD(PH1+5.*TWOPI, TWOPI)
      FX=A*COS(PH1) 
      FY=A*SIN(PH1)
      
      END SUBROUTINE lomb
