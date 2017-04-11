module modCircuitBC

	use modCircuit

	implicit none

contains

	subroutine Electrode(this,A)
		class(circuit), intent(inout) :: this
		real(mp), intent(inout) :: A(:,:)

		this%rho_back(this%nx) = this%rho_back(this%nx) + SUM(A(2,:))
	end subroutine

end module
