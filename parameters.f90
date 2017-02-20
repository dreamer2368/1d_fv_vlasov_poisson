module parameters

	use machinePrecision
	use modPlasma
	use modRecord
!these modules don't have to be added in the dependency of parameters.o, unless parameters.o uses their function or subroutines

	implicit none

	! parameters setup
	integer, parameter :: N = 10000
	real(mp), parameter :: pi = 4.0_mp*ATAN(1.0_mp)
	real(mp), parameter :: L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
	real(mp), parameter :: Lv = 1.0_mp
	real(mp), parameter :: v0 = 0.2_mp
	real(mp), parameter :: eps0 = 1.0_mp
	real(mp), parameter :: wp = 1.0_mp
	real(mp), parameter :: qe = -wp*wp/(N/L)
	real(mp), parameter :: me = -qe
	real(mp), parameter :: rho_back = -qe*N/L

	real(mp) :: Ti = 40.0_mp
	real(mp) :: Tp = 2.0_mp*pi/wp
!	real(mp) :: T = 40.0_mp
	real(mp) :: T = 40.0_mp + 10.0_mp*(2.0_mp*pi/wp)

	! time step parameter
	real(mp) :: dt, CFL=0.25_mp				!smaller initial perturbation A0 would require smaller CFL.
	integer :: Nt, Ni

	!grid and operators setup
	integer, parameter :: Nx = 256, Nv = 127
	real(mp) :: dx = L/Nx, dv = Lv/Nv
	real(mp) :: xg(Nx), vg(2*Nv+1)
	real(mp) :: rhs(Nx-1)

	!initial spatial distribution
	real(mp) :: f0(Nx,2*Nv+1)
	real(mp), parameter :: A0 = 0.001_mp, B0 = 1.0_mp, sigma = 0.01_mp
	real(mp) :: A = A0
	real(mp) :: B = B0
	integer, parameter :: mode = 1

	!time-stepping particle & field variable are composed of type fluid
	type(plasma) :: twostream
	type(history) :: plasmaRecord

	!error convergence variable
!	real(mp) :: fDA(30)
!	real(mp) :: ek(30)
!	real(mp) :: e(8,30)

!	character(32) :: filename ! You can make this longer than 32 characters, if needed

end module