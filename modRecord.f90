module modRecord

	use modPlasma
	use MatrixVector

	implicit none

	type history
		integer :: nt, nmod = 20
		real(mp) :: CFL, dt, T

		character(len=:), allocatable :: dir

		real(mp) :: KE
		real(mp) :: PE
		real(mp) :: TE
	end type

contains

	subroutine buildRecord(this,p,CFL,T,input_dir,nmod)
		type(history), intent(out) :: this
		type(plasma), intent(inout) :: p
		real(mp), intent(in) :: CFL, T
		character(len=*), intent(in), optional :: input_dir
		integer, intent(in), optional :: nmod
!		integer :: nr
		if( present(nmod) ) then
			this%nmod = nmod
		end if
		if( present(input_dir) ) then
			allocate(character(len=len(input_dir)) :: this%dir)
			this%dir = input_dir
		else
			allocate(character(len=0) :: this%dir)
			this%dir = ''
		end if

		!set timestep size according to CFL criterion with initial condition
		this%CFL = CFL
		call Efield_record(p)
		this%dt = CFL*MIN( p%dx/p%Lv, p%dv/MAXVAL(ABS(p%E))/ABS(p%qs)*p%ms )
		this%T = T
		this%nt = CEILING(T/this%dt)

		print *, 'T=',this%T,', CFL=',this%CFL,', Nt=',this%nt,', dt=',this%dt
		call system('mkdir -p data/'//this%dir//'/f')
		call system('rm data/'//this%dir//'/f/*.*')

!		nr = nt/this%nmod+1
	end subroutine

	subroutine destroyRecord(this)
		type(history), intent(inout) :: this

		deallocate(this%dir)
	end subroutine

	subroutine initRecord(this,p)
		type(history), intent(inout) :: this
		type(plasma), intent(in) :: p

		open(unit=301,file='data/'//this%dir//'/record_readme.out',status='replace')
		write(301,*) this%nt
		write(301,*) this%nmod
		write(301,*) this%dt
		write(301,*) p%nx
		write(301,*) p%nv
		write(301,*) p%Lx
		write(301,*) p%Lv
		close(301)

		open(unit=301,file='data/'//this%dir//'/xg.bin',status='replace',form='unformatted',access='stream')
		write(301) p%xg
		close(301)

		open(unit=301,file='data/'//this%dir//'/vg.bin',status='replace',form='unformatted',access='stream')
		write(301) p%vg
		close(301)

		open(unit=303,file='data/'//this%dir//'/E.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/'//this%dir//'/rho.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/'//this%dir//'/phi.bin',status='replace',form='unformatted',access='stream')
		open(unit=306,file='data/'//this%dir//'/KE.bin',status='replace',form='unformatted',access='stream')
		open(unit=307,file='data/'//this%dir//'/PE.bin',status='replace',form='unformatted',access='stream')
		open(unit=308,file='data/'//this%dir//'/TE.bin',status='replace',form='unformatted',access='stream')
	end subroutine

	subroutine closeRecord()
		close(303)
		close(304)
		close(305)
		close(306)
		close(307)
		close(308)
	end subroutine

	subroutine recordPlasma(this,p,k)
		type(plasma), intent(in) :: p
		type(history), intent(inout) :: this
		integer, intent(in) :: k
		integer :: j, kr
		character(len=100) :: kstr
		integer :: i
		real(mp) :: w(p%nx)

		if( (this%nmod.eq.1) .or. (mod(k,this%nmod).eq.0) ) then
			kr = merge(k,k/this%nmod,this%nmod.eq.1)
			write(kstr,*) kr

			w = 0.0_mp
			do i=1,2*p%nv
				w = w + 0.5_mp*( p%f(:,i)+p%f(:,i+1) )*( 0.5_mp*(p%vg(i)+p%vg(i+1)) )**2
			end do
			this%KE = 0.5_mp*p%ms*SUM(w)*p%dx*p%dv
			this%PE = 0.5_mp*p%eps0*SUM(p%E*p%E)*p%dx
			this%TE = this%KE + this%PE

			open(unit=302,file='data/'//this%dir//'/f/'//trim(adjustl(kstr))//		&
					'.bin',status='replace',form='unformatted',access='stream')
			write(302) p%f
			close(302)

			write(303) p%E
			write(304) p%rho
			write(305) p%phi
			write(306) this%KE
			write(307) this%PE
			write(308) this%TE
			print *, 'Time: ', k*this%dt, ', KE: ',this%KE,', PE: ',this%PE,', TE: ',this%TE
		end if
	end subroutine

	subroutine Efield_record(this)
		type(plasma), intent(inout) :: this
		real(mp) :: dx, dv
		real(mp) :: rhs(this%nx-1)
		real(mp) :: phi1(this%nx-1)
		dx = this%dx
		dv = this%dv

		this%rho = this%qs*integrate_dv(this%f,dv) + this%rho_back
		rhs = -1.0_mp/this%eps0*this%rho(2:this%nx)
		call CG_K(phi1,rhs,dx)
		this%phi(2:this%nx) = phi1
		this%phi(1) = 0.0_mp
		this%E = - multiplyD(this%phi,dx)
	end subroutine

end module