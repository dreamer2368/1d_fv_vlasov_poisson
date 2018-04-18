MODULE constants

	implicit none

	integer, parameter :: mp = SELECTED_REAL_KIND(15)
	real(mp), parameter :: pi = 4.0_mp*ATAN(1.0_mp)
    logical :: print_simulation_detail = .false.
    integer :: Ng_
    real(mp) :: Time_, A0_, A1_
    character(len=256) :: dir_

end MODULE
