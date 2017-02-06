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

      subroutine lomb(X,Y,N,N2,OFAC,SAMP_R,PX,PY,A,FFT_PHASE,PROB, 
     +                EFFM,START_I,END_I)
      
      IMPLICIT NONE

      !integer, parameter :: QR_K = selected_real_kind (32)

*  Input arguments
* 
      !INTEGER N2, OFAC
      INTEGER*8 :: N, N2, OFAC
      REAL*8 :: X(N), Y(N), FX(N2), SAMP_R
      REAL*8 :: PX(N2), PY(N2), A(N2), FFT_PHASE(N2)
      REAL *8 :: PROB(N2),EFFM
      

*  Local variables
* 
      INTEGER*8 :: HIFAC, PYMAX, DIMI, I, START_I, END_I
      REAL*8 :: WI(N), WPI(N), WPR(N), WR(N), WTEMP(N)
      REAL*8 :: PH1(N2), SUMCY, SUMSY, IY, FY(N2)
      REAL*8 :: RY, SS(N), CC(N), YY(N), PH(N2), PH2(N2)   
      REAL*8 :: PHLS, ARG(N)
      REAL*8 :: AVE, VARI, NMAX, XMAX, XMIN, XDIF
      REAL*8 :: P_MIN, P_MAX, F_MIN, F_MAX
      REAL*8 :: XAVE, PNOW, SUMSH, FAC, ARG0, ARG1
      REAL*8 :: WTAU, CWTAU, SWTAU, SUMS, SUMC, XSTART
      REAL*8:: EXPY(N2), P, GAP
      INTEGER*8:: NOUT, NUM
      CHARACTER(1) :: START_SWITCH, END_SWITCH
      REAL*8, PARAMETER :: TWOPI = 3.141592653589793238462643*2
      

Cf2py intent(in) X
Cf2py intent(in) Y
Cf2py intent(in) N
Cf2py intent(in) N2
Cf2py intent(in) OFAC
Cf2py intent(in) SAMP_R

Cf2py intent(out) PX
Cf2py intent(out) PY
Cf2py intent(out) A
Cf2py intent(out) FFT_PHASE
Cf2py intent(out) PROB
Cf2py intent(out) EFFM
Cf2py intent(out) START_I
Cf2py intent(out) END_I
Cf2py depend(N2) PX
Cf2py depend(N2) PY
Cf2py depend(N2) A
Cf2py depend(N2) FFT_PHASE
Cf2py depend(N2) PROB


      print *,'N = ', N
      PRINT *,'NOUT =', N2
     
      XSTART = X(1)
 
      XMAX=MAXVAL(X) 
      XMIN=MINVAL(X)
      XDIF=XMAX-XMIN
      
      PRINT *, 'XDIF = ', XDIF
      XAVE=0.5*(XMAX+XMIN)
      PYMAX=0
      NOUT=N2
      
      PX(1) = 12.
      PX(NOUT) = 0.001
      
      GAP = PX(1) - PX(NOUT)
      GAP = GAP/(NOUT-2)
      
      DO NUM = 2, NOUT-1
          PX(NUM) = PX(NUM-1)-GAP
      END DO 
                
      PNOW= PX(1)
      ARG=TWOPI*((X-XAVE)*PNOW)
      WPR=-2*SIN(0.5*ARG)**2
      WPI=SIN(ARG)
      WR=COS(ARG)
      PRINT *, ARG
      WI=SIN(ARG)
      AVE=SUM(Y)/SIZE(Y)
      YY=(Y-AVE)
      VARI=SUM((Y-AVE)**2)/SIZE(Y)
    
      DO NUM = 1, NOUT
          PNOW = PX(NUM)
          print *, PNOW
          SUMSH = SUM(WR*WI)
          SUMC = SUM((WR-WI)*(WR+WI))
   
          WTAU=0.5*ATAN2(2*SUMSH,SUMC)
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
          PRINT *, PY(NUM)

          PHLS=ATAN2(IY,RY)            
          ARG0=TWOPI*(XAVE+XSTART)*PNOW +WTAU  
          ARG1=TWOPI*XAVE*PNOW +WTAU   
          PH(NUM)=MOD((PHLS+ARG0), TWOPI)  
          PH1(NUM)=MOD((PHLS+ARG1), TWOPI)   
          !PNOW=PNOW+1./(OFAC*XDIF)   
      END DO
      DIMI=2*(N+1)      
      FAC=SQRT(VARI*DIMI/2)
      A=FAC*SQRT(PY) 
      
      
      A = A/N
      A = A*2
      FX=A*COS(PH1) 
      FY=A*SIN(PH1)  
      PH2=MOD((PH+5*TWOPI), TWOPI)
      FFT_PHASE = atan2(FY,FX)


!return the peak false alarm probabilities - the lower is the probability the more significant is the peak
      EXPY = EXP(-PY)  
      !EFFM=2.*NOUT/OFAC

      EFFM = (-6.362D0)+(1.193D0*N)+(0.00098D0*(N**2))
      EFFM = DBLE(INT(EFFM))

      PROB = EFFM*EXPY
      
      !loop designed to stop arithmetic errors
      DO I = 1, NOUT
          P = PROB(I)
          IF (P.gt.0.01d0) then
              PROB(I) = 1.d0-(1.d0-EXPY(I))**EFFM 
          END IF
      END DO

      !convert FAP to percentage significance
      PROB = (1.d0 - PROB)*100.d0
      
      END SUBROUTINE lomb
