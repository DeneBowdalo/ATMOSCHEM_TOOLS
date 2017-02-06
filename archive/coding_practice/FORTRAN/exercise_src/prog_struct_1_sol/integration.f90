
program integrate

  ! Program to integrate a function FUNC by the trapezium
  ! rule

  Implicit None

  ! Definition of the working precision

  integer, parameter :: dr = selected_real_kind(p=12,r=70) 

  ! Bounds on definite integral - integrate from  0.0 -> 1.0

  real(kind=dr), parameter :: lower = 0.0_dr
  real(kind=dr), parameter :: upper = 1.0_dr

  real(kind=dr) :: value

  integer :: steps

  ! Choose how many steps to take

  call read_input( steps )

  call integrate_trapezium( lower, upper, steps, value )
 
  Call write_result( value )

contains

  subroutine read_input( steps )

    ! Read the input

    integer, intent( out ) :: steps

    write( *, * ) 'How many steps should I use ?'
    read ( *, * ) steps

  end subroutine read_input


  subroutine integrate_trapezium( lower, upper, steps, value )

    ! Integrate a function from LOWER to UPPER with steps STEPS
    ! Use the trapezium rule

    real(kind=dr), Intent( In    ) :: lower
    real(kind=dr), Intent( In    ) :: upper
    integer      , Intent( In    ) :: steps
    real(kind=dr), Intent(   Out ) :: value

    real(kind=dr) :: func1, func2
    real(kind=dr) :: x, dx

    integer :: i

    dx = real( upper - lower, kind=dr ) / real( steps, kind=dr )

    value = 0.0_dr

    do i = 1, steps - 1

       x = lower + real( i - 1, kind=dr ) * dx

       func1 = func( x      )
       func2 = func( x + dx )

       value = value + dx * 0.5_dr * ( func1 + func2 )

    end do

  end subroutine integrate_trapezium


  subroutine write_result( value )

    real(kind=dr), intent( In ) :: value

    write( *, * ) 'The integral evaluates to ', value

  end subroutine write_result

  real(kind=dr) function func( x )

    real(kind=dr), intent( In ) :: x

    ! We will integrate this function

    if ( x == 0.0_dr ) then 

      func = 1.0_dr

    else

      func = sin(x)/x 
    
    end if

  end function func

end program integrate
