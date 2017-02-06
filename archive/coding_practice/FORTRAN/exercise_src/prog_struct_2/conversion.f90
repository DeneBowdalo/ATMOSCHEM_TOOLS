MODULE conversion

USE numbers
Implicit None

contains
subroutine deg_to_rad(value,rads)

Real(wp), Intent(In) :: value
Real(wp), Intent(Out) :: rads

!degrees = radians x (180/3.14)
rads = value*(180_wp/3.14_wp)
print*,value,'degrees equals', rads,'radians'

end subroutine deg_to_rad

subroutine rad_to_deg(value,degs)

Real(wp), Intent(In) :: value
Real(wp), intent(Out) :: degs
!radians = degrees * (3.14/180)
degs = value*(3.14_wp/180_wp)
print*,value,'radians equals', degs,'degrees'

end subroutine rad_to_deg

subroutine convert_time(value,con_hours,con_minutes,con_seconds) 
Real(wp), Intent(In) :: value
Real(wp), Intent(Out) :: con_hours, con_minutes, con_seconds
con_minutes = value / 60.0_wp
con_hours = value / 60.0_wp / 60.0_wp
con_seconds = mod (value, 60.0_wp)
WRITE(*,*) 'The time is', con_hours, 'hours', con_minutes, 'minutes and', con_seconds, 'seconds'

end subroutine convert_time

end module conversion
