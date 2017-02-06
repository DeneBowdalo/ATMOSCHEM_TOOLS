PROGRAM Leed_input
   IMPLICIT NONE
   REAL, DIMENSION(:,:,:), ALLOCATABLE :: phase_shift
   INTEGER :: no_energy, lmax, no_elements, i, j, alloc_err
   INTEGER :: i_err
   REAL, DIMENSION(:), ALLOCATABLE :: energy
!
!open the file containing the data then check if it was successful
!
   OPEN(UNIT=1,FILE='phase.dat',STATUS='OLD',ACTION='READ',IOSTAT = i_err)
   IF(i_err.NE.0)THEN
      WRITE(*,*)'AN ERROR OCCURRED OPENING THE FILE',i_err
      STOP
   END IF
   READ(1,*) no_energy,lmax,no_elements
!
!allocate the rank 3 array using the info from phase.dat 
!
   ALLOCATE(phase_shift(no_energy,(lmax+1),no_elements), &
      STAT = alloc_err)
   IF(alloc_err.ne.0)then
      WRITE(*,*)'AN ERROR OCCURRED DURING ALLOCATION OF phase_shift',alloc_err
      STOP
   END IF
!
!allocate the rank 1 array
!
   ALLOCATE(energy(no_energy), stat = alloc_err)
   IF(alloc_err.ne.0)then
      WRITE(*,*)'AN ERROR OCCURRED DURING ALLOCATION OF energy', alloc_err
      STOP
   END IF
   DO i = 1, no_energy
      READ(1,*)energy(i)
      DO j = 1, no_elements
         READ (1,'(7f7.4)') phase_shift(i,1:lmax + 1,j)
         WRITE(*,*) phase_shift(i,1:lmax + 1,j)
      END DO
   END DO
!
!close the file once it is finished with and deallocate the arrays
!
   CLOSE(UNIT=1,IOSTAT = i_err)
   IF(i_err.NE.0)THEN
       WRITE(*,*)'ERROR CLOSING FILE',i_err
   END IF
   DEALLOCATE(energy)
   DEALLOCATE(phase_shift)

END PROGRAM Leed_input
