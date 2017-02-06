PROGRAM linked_list

  IMPLICIT NONE
  TYPE node
    INTEGER :: value                    ! data filed
    TYPE (node), POINTER :: next        ! pointer field
  END TYPE node
  INTEGER :: num, status
  TYPE (node), POINTER :: head, tail, current

  ! build up the list

  NULLIFY(head, tail)                   ! initially nullify the list (empty)

  WRITE(*,*) 'Type-in an integer to build a linked list (0 to terminate)'

  READ(*,*) num                         ! read num from keyboard
  IF (num /= 0) then                    ! if 0 is entered, do nothing
    ALLOCATE(head, STAT = status)       ! create the head of the list

    IF (status > 0) STOP 'Fail to allocate a new node'

    head%value = num                    ! give the value
    NULLIFY(head%next)                  ! point to null
    tail => head                        ! update tail of list

    DO                                  ! create rest of list
      READ *, num                       ! read num from keyboard
      IF (num == 0) EXIT                ! until 0 is entered
      ALLOCATE(current, STAT = status)  ! create new node

      IF (status > 0) STOP 'Fail to allocate a new node'

      current%value = num               ! giving the value
      NULLIFY(current%next)             ! point to null (end of list)
      tail%next => current              ! link to tail of list
      tail => current                   ! update tail of list
    END DO
  END IF

  ! transverse the list and print the values

  WRITE(*,*) 'Transverse the list built up and print the values'

  current => head                       ! make current an alias of list
  DO
    IF (.NOT. ASSOCIATED(current)) EXIT ! exit if null pointer
    PRINT *, current%value              ! print the value
    current => current%next             ! make current alias of next node
  END DO

END PROGRAM linked_list
