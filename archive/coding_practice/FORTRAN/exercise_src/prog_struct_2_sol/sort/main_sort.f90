PROGRAM Sort_overload
  use numbers_module
  USE sort_it
  IMPLICIT NONE
  INTEGER :: n, alloc_err, i
  INTEGER, DIMENSION(:), ALLOCATABLE :: x
  REAL( wp ), DIMENSION(:), ALLOCATABLE :: y
!
! input the number of elements in the array
!
   DO
      WRITE(*,*)'Input number of elements in array'
      READ(*,*)n
      IF(n.LE.1)THEN
          WRITE(*,*)'The array should have more than 1 element'
      ELSE
          EXIT
      END IF
   END DO
!
! Allocate memory for the arrays and input
! store the input as both reals and ints to test our
! overloaded sort routine
!
   ALLOCATE(x(n),stat=alloc_err)
   IF(alloc_err.NE.0)THEN
       WRITE(*,*)'Error allocating array'
      STOP
   END IF
   ALLOCATE(y(n),stat=alloc_err)
   IF(alloc_err.NE.0)THEN
       WRITE(*,*)'Error allocating array'
      STOP
   END IF
   DO i = 1, n
      WRITE(*,*)'Input element ' ,i,' of array'
      READ(*,*)x(i)
      y(i)=REAL(x(i),wp) 
   END DO
!
! Output the unsorted arrays
!
   WRITE(*,*)'The unsorted integer array is;'
   WRITE(*,*)x
   WRITE(*,*)'The unsorted real array is;'
   WRITE(*,*)y

   CALL sort(x)
   CALL sort(y)

!
! and output the sorted arrays
!
   WRITE(*,*)'The sorted integer array is;'
   WRITE(*,*)x
   WRITE(*,*)'The sorted real array is;'
   WRITE(*,*)y

!
! and tidy up by deallocating memory
!
   DEALLOCATE(x)
   DEALLOCATE(y)

END PROGRAM Sort_overload
