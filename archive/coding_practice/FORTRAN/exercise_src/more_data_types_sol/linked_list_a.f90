PROGRAM simple_linked_list

  IMPLICIT NONE
  TYPE node
    INTEGER :: value                    ! data filed
    TYPE (node), POINTER :: next        ! pointer field
  END TYPE node
  INTEGER :: num, status
  TYPE (node), POINTER :: head, current

  ! build up the list

  NULLIFY(head)                         ! initially nullify list (empty)

  WRITE(*,*) 'Type-in an integer to build a linked list (0 to terminate)'

  DO
    READ(*,*) num                       ! read num from keyboard
    IF (num == 0) EXIT                  ! until 0 is entered
    ALLOCATE(current, STAT = status)    ! create new node

    IF (status > 0) STOP 'Fail to allocate a new node'

    current%value = num                 ! give the value
    current%next => head                ! point to previous one
    head => current                     ! update head of list
  END DO

  ! transverse the list and print the values

  WRITE(*,*) 'Transverse the list built up and print the values'

  current => head                       ! make current an alias of list
  DO
    IF (.NOT. ASSOCIATED(current)) EXIT ! exit if null pointer
    PRINT *, current%value              ! print the value
    current => current%next             ! make current alias of next node
  END DO

END PROGRAM simple_linked_list
