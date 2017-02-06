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

      subroutine lomb(X,Y,N,N2,PX,PY,A,PH)
      IMPLICIT NONE

Cf2py intent(in) X
Cf2py intent(in) Y

Cf2py intent(out) PX
Cf2py intent(out) PY
Cf2py intent(out) A
Cf2py intent(out) PH
Cf2py depend(N2) PX
Cf2py depend(N2) PY
Cf2py depend(N2) A
Cf2py depend(N2) PH


*  Input arguments
* 
      INTEGER N, N2
      REAL*4 :: X(N), Y(N)
      REAL*4 :: PX(N2), PY(N2), A(N2), PH(N2)

*  Local variables
* 
      INTEGER HIFAC, OFAC, PYMAX, NOUT, DIMI, NUM
      REAL*4 :: WI(N), WPI(N), WPR(N), WR(N), WTEMP(N)
      REAL*4 :: PH1(N), SUMCY, SUMSY, IY
      REAL*4 :: RY, SS(N), CC(N), YY(N)   
      REAL*4 :: PHLS, ARG(N)
      REAL*4 :: AVE, VARI, NMAX, XMAX, XMIN, XDIF
      REAL*4 :: XAVE, PNOW, SUMSH, FAC, ARG0, ARG1
      REAL*4 :: WTAU, CWTAU, SWTAU, SUMS, SUMC, XSTART
      REAL, PARAMETER :: TWOPI = 3.1415927*2

      print *,'N = ', N
      PRINT *,'N2 =', N2
      HIFAC=1
      OFAC=1
      NMAX=0.5*OFAC*HIFAC*N
     
      XSTART = X(1)
 
      XMAX=MAXVAL(X)
      !print *, XMAX 
      XMIN=MINVAL(X)
      XDIF=XMAX-XMIN
      XAVE=0.5*(XMAX+XMIN)
      PYMAX=0
      PNOW=1/(XDIF*OFAC)
      NOUT = INT(NMAX)
      PRINT *, 'NOUT = ', NOUT
      ARG=TWOPI*((X-XAVE)*PNOW)
      !print *, ARG
      WPR=-2*SIN(0.5*ARG)**2
      WPI=SIN(ARG)
      WR=COS(ARG)
      WI=SIN(ARG)
      AVE=SUM(Y)/SIZE(Y)
      YY=(Y-AVE)
      VARI=SUM((Y-AVE)**2)/SIZE(Y)

      DO NUM = 1, NOUT
          !print *, NUM
          !PRINT *, 'PERIODS =', PNOW
          PX(NUM)=PNOW
          SUMSH = SUM(WR*WI)
          !print *, SUMSH
          SUMC = SUM((WR-WI)*(WR+WI))
   
          WTAU=0.5*ATAN2(2*SUMSH,SUMC)
          !print *, wtau
          SWTAU=SIN(WTAU)
          CWTAU=COS(WTAU)

          SS=WI*CWTAU-WR*SWTAU
          CC=WR*CWTAU+WI*SWTAU

          SUMS = SUM((SS**2))
          SUMC = SUM((CC**2))
          SUMSY = SUM((YY*SS))
          SUMCY = SUM((YY*CC))
          !print *, SUMS
          !PRINT *, SUMSY
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
          IY=SUMSY/(SUMS**.5)  
          RY=SUMCY/(SUMC**.5) 
          PY(NUM)=0.5*(RY**2+IY**2)/VARI
          !print *, iy
          !print *, PY(NUM)

          PHLS=ATAN2(IY,RY)            
          ARG0=TWOPI*(XAVE+XSTART)*PNOW +WTAU  
          ARG1=TWOPI*XAVE*PNOW +WTAU   
          PH(NUM)=MOD((PHLS+ARG0), TWOPI)  
          PH1(NUM)=MOD((PHLS+ARG1), TWOPI)   
          PNOW=PNOW+1./(OFAC*XDIF)   
          !print *, PH(NUM)
      END DO
      !print *, size(px)
      !print *, size(py)
      DIMI=2*(N+1)      
      FAC=SQRT(VARI*DIMI/2)
      A=FAC*SQRT(PY) 
      A = A/N
      A = A*2
      !print *, size(A)
      !print *, size(PH)
      !print *, A
*      FX=A*COS(PH1) 
*      FY=S*SIN(PH1)  
      !PH2=MOD(PH+5*TWOPI, TWOPI)
      !print *, PH     
      !print *, 'Phase =', PH
      END SUBROUTINE lomb
