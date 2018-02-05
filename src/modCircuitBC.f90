module modCircuitBC

	use modCircuit

	implicit none

contains

	subroutine Electrode(this,p,dt)
		class(circuit), intent(inout) :: this
		type(plasma), intent(inout) :: p(:)
		real(mp), intent(in) :: dt
		real(mp) :: NetQdot
		NetQdot = p(1)%qs*p(1)%A(2)+p(2)%qs*p(2)%A(2)

		this%rho_back(this%nx) = this%rho_back(this%nx) + NetQdot*dt
	end subroutine

end module
