module modPlasma

	use MatrixVector
	use modPlasmaBC

	implicit none

	type plasma
		real(mp) :: Lx, Lv				![0, Lx]*[-Lv, Lv]
		integer :: nx, nv
		real(mp) :: dx, dv
		real(mp) :: qs=-1.0_mp, ms=1.0_mp
		real(mp), allocatable :: A(:)
		real(mp), allocatable :: f(:,:), n(:)
		real(mp), allocatable :: xg(:)
		real(mp), allocatable :: vg(:)
		procedure(PlasmaBC), nopass, pointer :: PtrBC=>Periodic
	end type

contains

	subroutine buildPlasma(this,Lx,Lv,nx,nv,qs,ms,A0)
		type(plasma), intent(out) :: this
		real(mp), intent(in) :: Lx, Lv
		integer, intent(in ) :: nx, nv
		real(mp), intent(in), optional :: qs,ms,A0(:)
		integer :: i
		if( present(qs) ) this%qs = qs
		if( present(ms) ) this%ms = ms
		if( present(A0) ) then
			allocate(this%A(SIZE(A0)))
			this%A = A0
		else
			allocate(this%A(1))
			this%A = 0.0_mp
		end if
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
		allocate(this%n(this%nx))

		this%f = 0.0_mp
		this%n = 0.0_mp
	end subroutine

	subroutine setPlasma(this,f0)
		type(plasma), intent(inout) :: this
		real(mp), intent(in) :: f0(this%nx,2*this%nv+1)

		this%f = f0
	end subroutine

	subroutine destroyPlasma(this)
		type(plasma), intent(inout) :: this
		deallocate(this%f)
		deallocate(this%n)
		deallocate(this%xg)
		deallocate(this%vg)

		deallocate(this%A)
	end subroutine

	subroutine NumberDensity(f,dv,n)
		real(mp), intent(in) :: f(:,:), dv
		real(mp), intent(out) :: n(size(f,1))
		n = integrate_dv(f,dv)
	end subroutine

end module
