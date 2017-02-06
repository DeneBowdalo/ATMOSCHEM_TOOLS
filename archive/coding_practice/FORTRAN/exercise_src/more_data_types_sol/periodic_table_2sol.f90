PROGRAM Periodic_table_2
   IMPLICIT NONE
   INTEGER, PARAMETER   :: elements = 10
   INTEGER              :: i_err, i, alloc_err, j, dealloc_err
!format specifier
!
   CHARACTER(len=38), PARAMETER  :: style ='(1x,a12,1x,a2,1x,i3,1x,f8.5,1x,i2,/)'
!
!create a derived type to store the atomic data
!
   TYPE atomic_data
      CHARACTER(len=12) :: name
      CHARACTER(len=2)  :: symbol
      INTEGER           :: number
      REAL              :: radius
      INTEGER           :: isotopes
      REAL, DIMENSION(:), POINTER :: weight
   END TYPE atomic_data
   TYPE(atomic_data)  periodic(elements)
!
!open the file with the atomic data in it
!use status and action since we only want to read input from the file
!
   OPEN(unit=1,file='atomic_data2.dat',status='old',action='read',iostat=i_err)
   IF(i_err.ne.0)THEN
      WRITE(*,*)'Error opening file',i_err
      STOP
   END IF
   DO i = 1, elements
      READ(1,*)periodic(i)%name, periodic(i)%symbol, periodic(i)%number, &
               periodic(i)%radius, periodic(i)%isotopes
!
!allocate the pointer
!
      ALLOCATE(periodic(i)%weight(periodic(i)%isotopes),stat=alloc_err)
      IF(alloc_err.ne.0)THEN
          WRITE(*,*)'Error in allocation'
          STOP
      END IF
      READ(1,*)periodic(i)%weight
   END DO
   DO i = 1, elements
      WRITE(*,style)periodic(i)%name, periodic(i)%symbol, periodic(i)%number, &
               periodic(i)%radius, periodic(i)%isotopes
      DO j = 1, periodic(i)%isotopes
         WRITE(*,*)periodic(i)%weight(j)
      END DO
   END DO
!
!close the file once finished with it
!
   CLOSE(unit=1,iostat=i_err)
   IF(i_err.ne.0)THEN
       WRITE(*,*)'An error occurred closing the file'
       STOP
   END IF
   DO i=1,elements
      DEALLOCATE(periodic(i)%weight ,stat=dealloc_err)
   END DO
    IF(dealloc_err.ne.0)THEN
          WRITE(*,*)'Error in deallocation'
          STOP
      END IF
END PROGRAM Periodic_table_2
