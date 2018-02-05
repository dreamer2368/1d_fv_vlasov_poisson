module modQoI

	use modPlasma
	use modCircuit

	implicit none

	abstract interface
		subroutine QoI(p,c,k,j)
			use modPlasma
			use modCircuit
			type(plasma), intent(in) :: p(:)
			type(circuit), intent(in) :: c
			integer, intent(in) :: k
			real(mp), intent(inout) :: j
		end subroutine
	end interface

contains

	subroutine DO_NOTHING(p,c,k,j)
		type(plasma), intent(in) :: p(:)
		type(circuit), intent(in) :: c
		integer, intent(in) :: k
		real(mp), intent(inout) :: j
	end subroutine

	subroutine Screening_distance(p,c,k,j)
		type(plasma), intent(in) :: p(:)
		type(circuit), intent(in) :: c
		integer, intent(in) :: k
		real(mp), intent(inout) :: j
		real(mp), dimension(c%nx) :: n
		real(mp) :: xc
		xc=0.5_mp*c%Lx

		n = integrate_dv(p(1)%f,p(1)%dv)
		j = SUM( (c%xg-xc)**2*n )*c%dx
	end subroutine

end module
