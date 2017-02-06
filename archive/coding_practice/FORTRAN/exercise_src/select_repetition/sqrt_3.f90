PROGRAM sqrt_3
IMPLICIT None

wp = SELECTED_REAL_KIND(12,70)

np = SELECTED_REAL_KIND(6,35)

integer :: precis, xnew, x
real(wp) :: sixdp_precis
real(np) :: twelvedp_precis

write(*,*) 'Please write "6" for 6dp precision' &
           'or "12" for 12dp precison'
read(*,*) precis

write(*,*) 'What Number would to like to find the square root of?'
read(*,*) sixdp_precis & twelvedp_precis
 
if(precis == 6) then
    xnew=0.5(3/sixdp_precis + sixdp_precis)

else if(precis == 12) then
   xnew=0.5(3/twelvedp_precis + twelvedp_precis 

end if

END PROGRAM 
