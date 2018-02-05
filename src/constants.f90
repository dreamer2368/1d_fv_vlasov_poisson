MODULE constants

	implicit none

	integer, parameter :: mp = SELECTED_REAL_KIND(15)
	real(mp), parameter :: pi = 4.0_mp*ATAN(1.0_mp)

end MODULE
