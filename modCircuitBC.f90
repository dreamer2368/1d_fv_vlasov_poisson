module modCircuitBC

	use modCircuit

	implicit none

contains

	subroutine Electrode(this,A,dt)
		class(circuit), intent(inout) :: this
		real(mp), intent(inout) :: A(:,:)
		real(mp), intent(in) :: dt

		this%rho_back(this%nx) = this%rho_back(this%nx) + SUM(A(2,:))*dt
	end subroutine

end module
