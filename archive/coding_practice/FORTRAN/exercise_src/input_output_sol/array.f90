PROGRAM array
  IMPLICIT NONE
  INTEGER :: i,j, n, stat
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: a
  CHARACTER(LEN = 10) :: edit_str
  
  WRITE(*,*)'Enter n'
  READ(*,*) n
  
  ALLOCATE(a(n,n), STAT=stat)

  IF(stat /= 0) STOP 'Failure to allocate'
  
  DO j = 1, n
     DO i = 1, n
        a(i,j) = i*10+j
     END DO
  END DO
!  
! Write out the array with non-advancing IO
!
  DO i = 1,n
     DO j = 1, n
        WRITE(UNIT=*, FMT='(I4)', ADVANCE='NO') a(i,j)
     END DO
     WRITE(*,*)
  END DO
  
  WRITE(*,*)
!  
! Write out the array with a dynamic edit descriptor
!
  edit_str = gen_edit_str(n)
  DO i = 1,n
     WRITE(UNIT=*, FMT=edit_str) a(i,:)
  END DO
  
CONTAINS
!  
!  Function to generate edit descriptor for n I4 integers!
!
  FUNCTION gen_edit_str(n)
    INTEGER :: n
    CHARACTER(LEN=8) :: gen_edit_str
    CHARACTER(LEN=20) :: temp
!
! Write the value of n to the string temp
!    
    WRITE(temp,*) n
!
! use concatenate to build up the edit descriptor
!    
    gen_edit_str = '(' // TRIM(ADJUSTL(temp)) // 'I4)'
    
  END FUNCTION gen_edit_str
  
END PROGRAM array
