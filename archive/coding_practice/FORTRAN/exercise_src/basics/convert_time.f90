PROGRAM convert_time
IMPLICIT None

INTEGER :: seconds, con_minutes, con_hours, con_seconds

READ (*,*) seconds

con_minutes = seconds / 60
con_hours = seconds / 60 / 60 
con_seconds = mod (seconds, 60) 
WRITE(*,*) 'The time is', con_hours, 'hours', con_minutes, 'minutes and', con_seconds, 'seconds' 

END PROGRAM convert_time
