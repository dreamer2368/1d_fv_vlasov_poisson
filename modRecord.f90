module modRecord

	use machinePrecision
	use modPlasma
	use MatrixVector

	implicit none

	type history
		integer :: nt
!		integer :: ni
		integer :: nx
		integer :: nv
		integer :: nr = 1600						!number of taking record

		real(mp), allocatable :: xg(:)
		real(mp), allocatable :: vg(:)

		real(mp) :: dx,dv,dt

		real(mp) :: KE
		real(mp) :: PE
		real(mp) :: TE
	end type

contains

	subroutine buildRecord(this,xg,vg,nt,dx,dv,dt)
		type(history), intent(out) :: this
		real(mp), intent(in) :: xg(:)
		real(mp), intent(in) :: vg(:)
		real(mp), intent(in) :: dx,dv,dt
		integer, intent(in) :: nt

		this%nx = size(xg)
		this%nv = ( size(vg)-1 )/2
		this%nt = nt
		this%nr = nt/(nt/this%nr)

		this%xg = xg
		this%vg = vg

		this%dx = dx
		this%dv = dv
		this%dt = dt
	end subroutine

	subroutine destroyRecord(this)
		type(history), intent(inout) :: this

	end subroutine

	subroutine initRecord(this)
		type(history), intent(inout) :: this

		open(unit=1,file='data/record_readme.out',status='replace')
		write(1,*) this%nt
		write(1,*) this%nr
		write(1,*) this%nx
		write(1,*) this%nv
		write(1,*) this%dx
		write(1,*) this%dv
		write(1,*) this%dt
		close(1)

		open(unit=1,file='data/xg.bin',status='replace',form='unformatted',access='stream')
		write(1) this%xg
		close(1)

		open(unit=1,file='data/vg.bin',status='replace',form='unformatted',access='stream')
		write(1) this%vg
		close(1)

		open(unit=1,file='data/T.bin',status='replace',form='unformatted',access='stream')
		open(unit=2,file='data/f.bin',status='replace',form='unformatted',access='stream')
		open(unit=3,file='data/E.bin',status='replace',form='unformatted',access='stream')
		open(unit=4,file='data/rho.bin',status='replace',form='unformatted',access='stream')
		open(unit=5,file='data/phi.bin',status='replace',form='unformatted',access='stream')
		open(unit=6,file='data/KE.bin',status='replace',form='unformatted',access='stream')
		open(unit=7,file='data/PE.bin',status='replace',form='unformatted',access='stream')
		open(unit=8,file='data/TE.bin',status='replace',form='unformatted',access='stream')
	end subroutine

	subroutine closeRecord()
		close(1)
		close(2)
		close(3)
		close(4)
		close(5)
		close(6)
		close(7)
		close(8)
	end subroutine

	subroutine recordPlasma(this,p,k,me,eps0)
		type(plasma), intent(in) :: p
		type(history), intent(inout) :: this
		integer, intent(in) :: k
		real(mp), intent(in) :: me,eps0
		integer :: i
		real(mp) :: w(this%nx)

		if( MOD(k,this%nt/this%nr)==0 ) then
			if( k>this%nt ) then
				print *, 'ERROR: attempt to save beyond the ending time T'
			end if

			w = 0.0_mp
			do i=1,size(this%vg)-1
				w = w + 0.5_mp*( p%f(:,i)+p%f(:,i+1) )*( 0.5_mp*(this%vg(i)+this%vg(i+1)) )**2
			end do
			this%KE = 0.5_mp*me*SUM(w)*this%dx*this%dv
			this%PE = 0.5_mp*eps0*SUM(p%E*p%E)*this%dx
			this%TE = this%KE + this%PE

			write(1) p%time
			write(2) p%f
			write(3) p%E
			write(4) p%rho
			write(5) p%phi
			write(6) this%KE
			write(7) this%PE
			write(8) this%TE
		end if
	end subroutine

end module