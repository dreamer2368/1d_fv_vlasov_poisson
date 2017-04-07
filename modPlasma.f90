module modPlasma

	use modPlasmaBC

	implicit none

	type plasma
		real(mp) :: Lx, Lv				![0, Lx]*[-Lv, Lv]
		integer :: nx, nv
		real(mp) :: dx, dv
		real(mp) :: qs=-1.0_mp, ms=1.0_mp, eps0=1.0_mp
		real(mp) :: A
		real(mp), allocatable :: f(:,:)
		real(mp), allocatable :: E(:)
		real(mp), allocatable :: rho(:), rho_back(:)
		real(mp), allocatable :: phi(:)
		real(mp), allocatable :: xg(:)
		real(mp), allocatable :: vg(:)
		procedure(PlasmaBC), nopass, pointer :: PtrBC=>Periodic
	end type

contains

	subroutine buildPlasma(this,Lx,Lv,nx,nv,qs,ms,eps0)
		type(plasma), intent(out) :: this
		real(mp), intent(in) :: Lx, Lv
		integer, intent(in ) :: nx, nv
		real(mp), intent(in), optional :: qs,ms,eps0
		integer :: i
		if( present(qs) ) this%qs = qs
		if( present(ms) ) this%ms = ms
		if( present(eps0) ) this%eps0 = eps0
		this%Lx = Lx
		this%Lv = Lv
		this%nx = nx
		this%nv = nv
		this%dx = Lx/nx
		this%dv = Lv/nv

		allocate(this%xg(this%nx))
		allocate(this%vg(2*this%nv+1))
		this%xg = (/ ( (i-0.5_mp)/nx*Lx, i=1,nx ) /)
		this%vg = (/ ( Lv*(i-nv-1)/nv, i=1,2*nv+1 ) /)

		allocate(this%f(this%nx,2*this%nv+1))
		allocate(this%E(this%nx))
		allocate(this%phi(this%nx))
		allocate(this%rho(this%nx))
		allocate(this%rho_back(this%nx))

		this%f = 0.0_mp
		this%E = 0.0_mp
		this%phi = 0.0_mp
		this%rho = 0.0_mp
		this%rho_back = 0.0_mp
	end subroutine

	subroutine setPlasma(this,f0)
		type(plasma), intent(inout) :: this
		real(mp), intent(in) :: f0(this%nx,2*this%nv+1)

		this%f = f0
	end subroutine

	subroutine setBackGround(this,rho_back)
		type(plasma), intent(inout) :: this
		real(mp), intent(in) :: rho_back(this%nx)

		this%rho_back = rho_back
	end subroutine

	subroutine destroyPlasma(this)
		type(plasma), intent(inout) :: this
		deallocate(this%f)
		deallocate(this%E)
		deallocate(this%rho)
		deallocate(this%phi)
		deallocate(this%xg)
		deallocate(this%vg)
	end subroutine

end module
