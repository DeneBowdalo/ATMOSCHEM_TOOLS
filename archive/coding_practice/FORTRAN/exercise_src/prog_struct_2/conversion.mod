VE 0 0 3 0 0 0
MODULE CONVERSION,0 0
FILE 0,conversion.f90
USE NUMBERS 2
PROC CONVERT_TIME,4,8,0,17: 8,0,0,0,0,40000,1,0,0,0
VAR VALUE,3,,: 2,2,5,0,0,103,0,0,0,1
VAR CON_HOURS,3,,: 2,2,5,0,0,183,0,0,0,2
VAR CON_MINUTES,3,,: 2,2,5,0,0,183,0,0,0,2
VAR CON_SECONDS,3,,: 2,2,5,0,0,183,0,0,0,2
ENDPROC
PROC RAD_TO_DEG,2,8,0,17: 8,0,0,0,0,40000,1,0,0,0
VAR VALUE,3,,: 2,2,5,0,0,103,0,0,0,1
VAR DEGS,3,,: 2,2,5,0,0,183,0,0,0,2
ENDPROC
PROC DEG_TO_RAD,2,8,0,17: 8,0,0,0,0,40000,1,0,0,0
VAR VALUE,3,,: 2,2,5,0,0,103,0,0,0,1
VAR RADS,3,,: 2,2,5,0,0,183,0,0,0,2
ENDPROC
END
