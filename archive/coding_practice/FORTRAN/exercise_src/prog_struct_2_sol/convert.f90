MODULE Convert
  USE numbers_module
  IMPLICIT NONE
  Public :: convert_time, deg_to_rad, rad_to_deg
  PRIVATE

CONTAINS

  SUBROUTINE Convert_time(timein,hours,mins,secs)
    IMPLICIT NONE
    REAL(WP), INTENT(IN) :: timein 
    REAL(WP), INTENT(OUT):: hours, mins, secs
    REAL(WP) :: temp
!
! if the time is positive
!
    IF(timein.GE.0)THEN
       hours = FLOOR(timein/3600)
       mins = FLOOR((timein-(hours*3600))/60)
       secs = timein - (hours*3600) - (mins*60)
!
! if the time is negative
!
    ELSE 
       temp = ABS(timein)
       hours = FLOOR(temp/3600)
       mins = FLOOR((temp-(hours*3600))/60)
       secs = temp - (hours*3600) - (mins*60)
    END IF
  END SUBROUTINE Convert_time

  SUBROUTINE Deg_to_Rad(angle)
    IMPLICIT NONE
    REAL(WP), INTENT(INOUT) :: angle
    REAL(WP) :: pi
!
! calculate pi
!
    pi = 4.0_wp*ATAN(1.0_wp)
!
! convert the angle from degrees to radians
!
    angle = angle*pi/180.0_wp

  END SUBROUTINE Deg_to_Rad

  SUBROUTINE Rad_to_Deg(angle)
    IMPLICIT NONE
    REAL(WP), INTENT(INOUT) :: angle
    REAL(WP) :: pi
!
! calculate pi
!
    pi = 4.0_wp*ATAN(1.0_wp)

!
! convert the angle from radians to degrees
!    
    angle = angle*180.0_wp/pi
    
  END SUBROUTINE RAD_TO_DEG


END MODULE CONVERT
