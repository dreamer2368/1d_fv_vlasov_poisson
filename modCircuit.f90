module modCircuit

	use MatrixVector
!	use constants
!	use modBC

	implicit none

	type circuit
		real(mp) :: Lx				![0, Lx]*[-Lv, Lv]
		integer :: nx
		real(mp) :: dx
		real(mp) :: eps0=1.0_mp
		real(mp) :: A
		real(mp), allocatable :: E(:)
		real(mp), allocatable :: rho(:), rho_back(:)
		real(mp), allocatable :: phi(:)
		real(mp), allocatable :: xg(:)
!		procedure(BC), nopass, pointer :: PtrBC=>Periodic
	end type

contains

	subroutine buildCircuit(this,Lx,nx,eps0)
		type(circuit), intent(out) :: this
		real(mp), intent(in) :: Lx
		integer, intent(in) :: nx
		real(mp), intent(in), optional :: eps0
		integer :: i
		if( present(eps0) ) this%eps0 = eps0
		this%Lx = Lx
		this%nx = nx
		this%dx = Lx/nx

		allocate(this%xg(this%nx))
		this%xg = (/ ( (i-0.5_mp)/nx*Lx, i=1,nx ) /)

		allocate(this%E(this%nx))
		allocate(this%phi(this%nx))
		allocate(this%rho(this%nx))
		allocate(this%rho_back(this%nx))

		this%E = 0.0_mp
		this%phi = 0.0_mp
		this%rho = 0.0_mp
		this%rho_back = 0.0_mp
	end subroutine

	subroutine setBackGround(this,rho_back)
		type(circuit), intent(inout) :: this
		real(mp), intent(in) :: rho_back(this%nx)

		this%rho_back = rho_back
	end subroutine

	subroutine destroyCircuit(this)
		type(circuit), intent(inout) :: this
		deallocate(this%E)
		deallocate(this%rho)
		deallocate(this%rho_back)
		deallocate(this%phi)
		deallocate(this%xg)
	end subroutine

	subroutine Efield(this)
		type(circuit), intent(inout) :: this
		real(mp) :: dx
		real(mp) :: rhs(this%nx-1)
		real(mp) :: phi1(this%nx-1)
		dx = this%dx

		rhs = -1.0_mp/this%eps0*( this%rho(2:this%nx) + this%rho_back(2:this%nx) )
		call CG_K(phi1,rhs,dx)
		this%phi(2:this%nx) = phi1
		this%phi(1) = 0.0_mp
		this%E = - multiplyD(this%phi,dx)
	end subroutine

end module
