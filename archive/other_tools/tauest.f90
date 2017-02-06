!----------------------------------------------------------------------
!  Manfred Mudelsee's code for tau estimation
!----------------------------------------------------------------------
! TAUEST: Routine for persistence estimation for unevenly spaced time series
!----------------------------------------------------------------------
!       Main variables
!
!       t       :       time
!       x       :       time series value
!       np      :       number of points
!      dt       :       average spacing
!   scalt       :       scaling factor (time)
!     rho       :       in the case of equidistance, rho = autocorr. coeff.
!      ls       :       LS function
!   brent       :       Brent's search, minimum LS value
!    mult       :       flag (multiple solution)
!    amin       :       estimated value of a = exp(-scalt/tau)
!
!----------------------------------------------------------------------
  subroutine tauest(t, x, np,ave,var,tau,rhoavg)
!

  use const
  use error
!
  implicit none
!
  integer, intent(in) :: np
  real, dimension(np), intent(in) :: t, x
  real, intent (out)  :: tau
  real, dimension(np) :: tscal, xscal
  real :: fac, avg, var, dt, rho, scalt, amin, rhoavg
  double precision :: damin
  integer :: i, mult
!
  external brent, ls, minls, rhoest
!

!f2py intent(in) t
!f2py intent(in) x
!f2py intent(in) np
!f2py intent(in) ave
!f2py intent(in) var

!f2py intent(out) tau
!f2py intent(out) rhoavg


!
! Correct time direction; assume that ages are input
! --------------------------------------------------
  do i = 1, np
      tscal(i) = -t(np+1-i)
      xscal(i) = x(np+1-i)
  end do
!
! Scaling of x
! ------------
  fac = sqrt(var)
  xscal(1:np) = xscal(1:np) / fac
!
! Scaling of t (=> start value of a = 1/e)
! ---------------------------------------
  dt = (tscal(np)-tscal(1)) / real((np-1))
  call rhoest(np, xscal(1:np), rho)
  if (rho .le. 0.0) then
      rho = 0.05
      write(errio,*) 'Warning: rho estimation: < 0'
      ierr = 2
  else if (rho .gt. 1.0) then
      rho = 0.95
      write(errio,*) 'Warning: rho estimation: > 1'
      ierr = 2
  end if
  scalt = -log(rho)/dt
  tscal(1:np) = tscal(1:np) * scalt
!
! Estimation
! ----------
  call minls(np, dble(tscal(1:np)), dble(xscal(1:np)), damin, mult)
  if (ierr .eq. 1) then
     write(errio,*) ' Error in MNILS'
     return
  end if
  amin = sngl(damin)
  if (mult .eq. 1) then
     write(errio,*) ' Estimation problem: LS function has > 1 minima'
     return
  end if
  if (amin .le. 0.0) then
     write(errio,*) ' Estimation problem: a_min =< 0'
     return
  else if (amin .ge. 1.0) then
     write(errio,*) ' Estimation problem: a_min >= 1'
     return
  end if
!
! determine tau
! -------------
  tau = -1.0 /(scalt*log(amin))
!
! determine rho, corresponding to tau
! -----------------------------------
  rhoavg = exp(-dt / tau)
!
  end subroutine tauest
!
!
!----------------------------------------------------------------------
! Numerical Recipes (modified): Brent's search in one direction.
!----------------------------------------------------------------------
  function brent(ax,bx,cx,f,tol,xmin,xfunc,yfunc,nfunc)
!
  use error
  implicit none
!
  integer nfunc
  double precision xfunc(1:nfunc),yfunc(1:nfunc)
  integer itmax
  double precision brent,ax,bx,cx,tol,xmin,f,cgold,zeps
  external f
  parameter (itmax=100,cgold=.3819660d0,zeps=1.d-18)
  integer iter
  double precision a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
!
  a=min(ax,cx)
  b=max(ax,cx)
  v=bx
  w=v
  x=v
  e=0.d0
  fx=f(x,xfunc,yfunc,nfunc)
  fv=fx
  fw=fx
  do iter=1,itmax
    xm=0.5d0*(a+b)
    tol1=tol*abs(x)+zeps
    tol2=2.d0*tol1
    if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
    if(abs(e).gt.tol1) then
      r=(x-w)*(fx-fv)
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=2.d0*(q-r)
      if(q.gt.0.d0) p=-p
      q=abs(q)
      etemp=e
      e=d
      if(abs(p).ge.abs(.5d0*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x))goto 1
      d=p/q
      u=x+d
      if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
      goto 2
    endif
1   if(x.ge.xm) then
      e=a-x
    else
      e=b-x
    endif
    d=cgold*e
2   if(abs(d).ge.tol1) then
      u=x+d
    else
      u=x+sign(tol1,d)
    endif
    fu=f(u,xfunc,yfunc,nfunc)
    if(fu.le.fx) then
      if(u.ge.x) then
        a=x
      else
        b=x
      endif
      v=w
      fv=fw
      w=x
      fw=fx
      x=u
      fx=fu
    else
      if(u.lt.x) then
        a=u
      else
        b=u
      endif
      if(fu.le.fw .or. w.eq.x) then
        v=w
        fv=fw
        w=u
        fw=fu
      else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
        v=u
        fv=fu
      endif
    endif
  end do
  ierr = 1
  write(errio,*) ' brent: exceed maximum iterations'
3 xmin=x
  brent=fx
  return
  end
!
!
!----------------------------------------------------------------------
! Least-squares function
!----------------------------------------------------------------------
  double precision function ls(a,t,x,n)
  implicit none
  integer n
  double precision t(1:n),x(1:n)
  double precision a
  integer i
  ls=0.0d0
  do i=2,n
     ls=ls+(x(i)-x(i-1)*dsign(1.0d0,a)* dabs(a)**(t(i)-t(i-1)))**2.0d0
  end do
  return
  end
!
!
!----------------------------------------------------------------------
! Minimization of least-squares function ls.
!----------------------------------------------------------------------
  subroutine minls(n, t, x, amin, nmu_)
!
  use error
!
  implicit none
!
  double precision, parameter :: a_ar1 = 0.367879441d0 ! 1/e
  double precision, parameter :: tol = 3.0d-8          ! Brent's search, precision
  double precision, parameter :: tol2 = 1.0d-6         ! multiple solutions, precision
  integer n
  double precision t(1:n),x(1:n)
  double precision amin
  integer nmu_
  double precision dum1,dum2,dum3,dum4,a_ar11,a_ar12,a_ar13
  double precision ls,brent
  external ls,brent
!
  nmu_=0
  dum1=brent(-2.0d0, a_ar1, +2.0d0, ls, tol, a_ar11, t, x, n)
  dum2=brent( a_ar1, 0.5d0*(a_ar1+1.0d0), +2.0d0, ls, tol, a_ar12, t, x, n)
  dum3=brent(-2.0d0, 0.5d0*(a_ar1-1.0d0),  a_ar1, ls, tol, a_ar13, t, x, n)
  if (ierr .eq. 1) then
     write(errio, *) ' Error in MINLS (call to brent)'
     return
  end if
  if  ((dabs(a_ar12-a_ar11).gt.tol2.and.dabs(a_ar12-a_ar1).gt.tol2) &
  .or.(dabs(a_ar13-a_ar11).gt.tol2.and.dabs(a_ar13-a_ar1).gt.tol2)) &
  nmu_=1
  dum4=dmin1(dum1,dum2,dum3)
  if (dum4.eq.dum2) then
     amin=a_ar12
  else if (dum4.eq.dum3) then
     amin=a_ar13
  else
     amin=a_ar11
  end if
  return
  end
!
!
!----------------------------------------------------------------------
! Autocorrelation coefficient estimation (equidistant data).
!----------------------------------------------------------------------
  subroutine rhoest(n,x,rho)
!
  implicit none
!
  integer n
  real x(1:n)
  real rho
  integer i
  real sum1,sum2
!
  sum1=0.0
  sum2=0.0
  do i=2,n
     sum1=sum1+x(i)*x(i-1)
     sum2=sum2+x(i)**2.0
  end do
  rho=sum1/sum2
  return
  end
!
!
!----------------------------------------------------------------------
!     Numerical Recipes code
!----------------------------------------------------------------------
  include 'c:\msdata\src\f90\recipes\recipes\avevar.f90'
  include 'c:\msdata\src\f77\recipes\recipes\gasdev.f'
  include 'c:\msdata\src\f90\recipes\recipes\ran.f90'
  include 'c:\msdata\src\f90\recipes\recipes\gammp.f90'
  include 'c:\msdata\src\f90\recipes\recipes\gser.f90'
  include 'c:\msdata\src\f90\recipes\recipes\gcf.f90'
  include 'c:\msdata\src\f90\recipes\recipes\gammln.f90'
  include 'c:\msdata\src\f90\recipes\recipes\erfcc.f90'
  include 'c:\msdata\src\f90\recipes\recipes\sort.f90'
