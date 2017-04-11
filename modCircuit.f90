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
		procedure(MeshSolver), pass(this), pointer :: Efield=>PeriodicCircuit
		procedure(updateCircuit), pass(this), pointer :: updateCircuit=>Null_Circuit
	end type

	abstract interface
		subroutine MeshSolver(this)
			import circuit
			class(circuit), intent(inout) :: this
		end subroutine
	end interface

	abstract interface
		subroutine updateCircuit(this,A,dt)
			use constants
			import circuit
			class(circuit), intent(inout) :: this
			real(mp), intent(inout) :: A(:,:)
			real(mp), intent(in) :: dt
		end subroutine
	end interface

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

!=============  MeshSolver subroutines  =============================================

	subroutine PeriodicCircuit(this)
		class(circuit), intent(inout) :: this
		real(mp) :: dx
		real(mp), dimension(this%nx-1) :: rhs, phi1, co1, co2, co3
		dx = this%dx
		co1 = 1.0_mp/dx/dx
		co2 = -2.0_mp/dx/dx
		co3 = 1.0_mp/dx/dx

		rhs = -1.0_mp/this%eps0*( this%rho(2:this%nx) + this%rho_back(2:this%nx) )
!		call CG_K(phi1,rhs,dx)
		call solve_tridiag(co1,co2,co3,rhs,phi1,this%nx-1)
		this%phi(2:this%nx) = phi1
		this%phi(1) = 0.0_mp
		this%E = - multiplyD(this%phi,dx)
	end subroutine

	subroutine DirichletNeumann(this)						!D(i=1), N(i=N)
		class(circuit), intent(inout) :: this
		real(mp) :: dx, eps
		integer :: ng
		real(mp), dimension(this%nx) :: rhs, phi1, co1, co2, co3
		dx = this%dx
		ng = this%nx
		eps = this%eps0
		co1 = 1.0_mp/dx/dx
		co2 = -2.0_mp/dx/dx
		co3 = 1.0_mp/dx/dx
		co2(ng) = -1.0_mp/dx/dx

		!rho_back(ng): surface charge
		rhs = -( this%rho + this%rho_back/dx )/eps
		call solve_tridiag(co1,co2,co3,rhs,phi1,ng)
		this%phi = phi1

		!Efield
		this%E(ng) = 0.5_mp*( -(phi1(ng)-phi1(ng-1))/dx-this%rho_back(ng)/eps )
		this%E(2:ng-1) = -( phi1(3:ng)-phi1(1:ng-2) )/2.0_mp/dx
		this%E(1) = -( phi1(2) )/2.0_mp/dx
	end subroutine

!===============    BASE CIRCUIT UPDATE SUBROUTINE   ===============================

	subroutine	Null_Circuit(this,A,dt)
		class(circuit), intent(inout) :: this
		real(mp), intent(inout) :: A(:,:)
		real(mp), intent(in) :: dt
	end subroutine

end module
