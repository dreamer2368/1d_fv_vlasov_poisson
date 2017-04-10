module modQoI

	use modPlasma
	use modCircuit

	implicit none

	abstract interface
		subroutine QoI(p,c,k,j)
			use modPlasma
			use modCircuit
			type(plasma), intent(in) :: p
			type(circuit), intent(in) :: c
			integer, intent(in) :: k
			real(mp), intent(inout) :: j
		end subroutine
	end interface

contains

	subroutine DO_NOTHING(p,c,k,j)
		type(plasma), intent(in) :: p
		type(circuit), intent(in) :: c
		integer, intent(in) :: k
		real(mp), intent(inout) :: j
	end subroutine

	subroutine Screening_distance(p,c,k,j)
		type(plasma), intent(in) :: p
		type(circuit), intent(in) :: c
		integer, intent(in) :: k
		real(mp), intent(inout) :: j
		real(mp), dimension(p%nx) :: n
		real(mp) :: xc
		xc=0.5_mp*p%Lx

		n = integrate_dv(p%f,p%dv)
		j = SUM( (p%xg-xc)**2*n )*p%dx
	end subroutine

end module
