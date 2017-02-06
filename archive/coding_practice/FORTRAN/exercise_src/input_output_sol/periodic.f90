PROGRAM Periodic_table
  IMPLICIT NONE
  CHARACTER(len=2),DIMENSION(36) :: atoms
  INTEGER                        :: i_err
!
! open the file containing the list of elements
! 
  OPEN(unit=1,file='atoms.dat',status='old',action='read',iostat=i_err)
  IF(i_err.NE.0)THEN
     WRITE(*,*)'Error opening file',i_err
     STOP
  END IF
!
! read the elements into an array
!
  READ(1,*)atoms
!
! write them out formatted like a standard periodic table
!
  WRITE(*,'(a2,49x,a2)') atoms( 1:2 )
  WRITE(*,'(2(a2,1x),30x,6(a2,1x))') atoms( 3:10 )
  WRITE(*,'(2(a2,1x),30x,6(a2,1x))') atoms( 11:18 )
  WRITE(*,'(17(a2,1x),a2)') atoms( 19:36 )
!
! close the file
!
  CLOSE(unit=1,iostat=i_err)
  IF(i_err.NE.0)THEN
     WRITE(*,*)'Error closing file',i_err
     STOP
  END IF
END PROGRAM Periodic_table

