PROGRAM Periodic_table
   IMPLICIT NONE
   INTEGER, PARAMETER   :: elements = 10
   INTEGER              :: i_err
!
!format specifier
!
   CHARACTER(len=44), PARAMETER  :: style ='(10(1x,a12,1x,a2,1x,i3,1x,f9.5,1x,f8.5,/))'
!
!create a derived type to store the atomic data
!
   TYPE atomic_data 
      CHARACTER(len=12) :: name
      CHARACTER(len=2)  :: symbol
      INTEGER           :: number
      REAL              :: weight
      REAL              :: radius
   END TYPE atomic_data
   TYPE(atomic_data)  periodic(elements)
!
!open the file with the atomic data in it
!use status and action since we only want to read input from the file
!
   OPEN(unit=1,file='atomic_data1.dat',status='old',action='read',iostat=i_err)
   IF(i_err.ne.0)THEN
      WRITE(*,*)'Error opening file',i_err
      STOP
   END IF
   READ(1,*)periodic
   WRITE(*,style)periodic
!
!close the file once finished with it
!
   CLOSE(unit=1,iostat=i_err)
   IF(i_err.ne.0)THEN
       WRITE(*,*)'An error occurred closing the file',i_err
       STOP
   END IF
END PROGRAM Periodic_table
