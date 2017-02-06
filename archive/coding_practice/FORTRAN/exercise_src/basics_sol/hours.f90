PROGRAM Convert_time
  IMPLICIT NONE
  REAL :: timein, hours, mins, secs         
  WRITE(*,*) "TYPE THE TIME IN SECONDS:"
  READ (*,*) timein
!
!if the time is positive
!
  IF(timein.GE.0)THEN
     hours = FLOOR(timein/3600)
     mins = FLOOR((timein-(hours*3600))/60)
     secs = timein - (hours*3600) - (mins*60)
     WRITE(*,*) INT(hours),":",INT(mins),":",INT(secs)
!
!if the time is negative
!
  ELSE 
     timein = ABS(timein)
     hours = FLOOR(timein/3600)
     mins = FLOOR((timein-(hours*3600))/60)
     secs = timein - (hours*3600) - (mins*60)
     WRITE(*,*) -INT(hours),":",-INT(mins),":",-INT(secs)
  END IF
END PROGRAM Convert_time
