MODULE constants

	implicit none

	integer, parameter :: mp = SELECTED_REAL_KIND(15)
	real(mp), parameter :: pi = 4.0_mp*ATAN(1.0_mp)
    logical :: print_simulation_detail = .false.
    integer :: Ng
    real(mp) :: Time, A0, A1
    character(len=256) :: dir

end MODULE
