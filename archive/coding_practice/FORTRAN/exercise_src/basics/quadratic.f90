PROGRAM quadratic
IMPLICIT None

integer, parameter :: wp = selected_real_kind(12,70)
real(wp) :: a, b, c, q, root_1, root_2, factor 


print*, 'Enter a, b and c'
read(*,*) a, b, c
 
factor = sqrt(b**2.0_wp-4.0_wp*a*c) 
q = -0.5_wp*((b+sign(1.0_wp,b))*factor)

root_1 = q/a
root_2 = c/q
print*,'The roots are:', root_1,'And', root_2 

END PROGRAM quadratic


