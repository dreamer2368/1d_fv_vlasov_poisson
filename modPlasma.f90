module modPlasma

	use machinePrecision

	implicit none

	type plasma
		integer :: nx
		integer :: nv
		real(mp) :: time
		real(mp), allocatable :: f(:,:)
		real(mp), allocatable :: E(:)
		real(mp), allocatable :: rho(:)
		real(mp), allocatable :: phi(:)
		real(mp), allocatable :: xg(:)
		real(mp), allocatable :: vg(:)
	end type

contains

	subroutine buildPlasma(this,xg,vg,f0)

		type(plasma), intent(out) :: this
		real(mp), intent(in) :: xg(:)
		real(mp), intent(in) :: vg(:)
		real(mp), intent(in) :: f0(:,:)

		this%nx = size(xg)
		this%nv = ( size(vg)-1 )/2

		allocate(this%f(this%nx,2*this%nv+1))
		allocate(this%E(this%nx))
		allocate(this%phi(this%nx))
		allocate(this%rho(this%nx))

		this%f = f0
		this%time = 0.0_mp
		this%xg = xg
		this%vg = vg

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