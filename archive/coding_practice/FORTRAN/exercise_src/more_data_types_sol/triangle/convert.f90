MODULE Convert
  IMPLICIT NONE
  PRIVATE
  PUBLIC Convert_time, Deg_to_rad, Rad_to_deg
CONTAINS

  SUBROUTINE Convert_time(timein,hours,mins,secs)
    IMPLICIT NONE
    REAL, INTENT(IN) :: timein 
    REAL, INTENT(OUT):: hours, mins, secs
    REAL :: temp
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

    REAL, INTENT(INOUT) :: angle
    REAL :: pi
!
! calculate pi
!
    pi = 4.0*ATAN(1.0)
!
! convert the angle from degrees to radians
!
    angle = angle*pi/180.0

  END SUBROUTINE deg_to_rad

  SUBROUTINE rad_to_deg(angle)  
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: angle
    REAL :: pi
!
! calculate pi
!
    pi = 4.0*ATAN(1.0)

!
! convert the angle from radians to degrees
!    
    angle = angle*180.0/pi
    
  END SUBROUTINE rad_to_deg


END MODULE Convert
