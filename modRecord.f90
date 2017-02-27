module modRecord

	use modPlasma
	use MatrixVector

	implicit none

	type history
		integer :: nt, nmod = 20
		real(mp) :: CFL, dt, T

		character(len=:), allocatable :: dir

		real(mp), allocatable :: KE(:), PE(:), TE(:)

		real(mp), allocatable :: rho(:,:), phi(:,:), E(:,:)
	end type

contains

	subroutine buildRecord(this,p,T,CFL,dt,input_dir,nmod)
		type(history), intent(out) :: this
		type(plasma), intent(inout) :: p
		real(mp), intent(in) :: T
		real(mp), intent(in), optional :: CFL, dt
		character(len=*), intent(in), optional :: input_dir
		integer, intent(in), optional :: nmod
		integer :: nr
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

		this%T = T
		!set timestep size according to CFL criterion with initial condition
		if( PRESENT(CFL) ) then
			this%CFL = CFL
			call Efield_record(p)
			if( MAXVAL(ABS(p%E)).ne.0.0_mp ) then
				print *, 'dx/v=',p%dx/p%Lv,', dv/acc=',p%dv/MAXVAL(ABS(p%E))/ABS(p%qs)*p%ms
				this%dt = CFL*MIN( p%dx/p%Lv, p%dv/MAXVAL(ABS(p%E))/ABS(p%qs)*p%ms )
			else
				this%dt = CFL*p%dx/p%Lv
			end if
			this%nt = CEILING(T/this%dt)
		!measure CFL according to timestep size
		elseif( PRESENT(dt) ) then
			this%nt = CEILING(T/dt)
			this%dt = T/this%nt
			call Efield_record(p)
			if( MAXVAL(ABS(p%E)).ne.0.0_mp ) then
				this%CFL = MAX( p%Lv*this%dt/p%dx, MAXVAL(ABS(p%E))*ABS(p%qs)/p%ms*this%dt/p%dv )
			else
				this%CFL = p%Lv*this%dt/p%dx
			end if
		else
			print *, 'ERROR: required to input either CFL or dt'
			stop
		end if

		print *, 'T=',this%T,', CFL=',this%CFL,', Nt=',this%nt,', dt=',this%dt
		call system('mkdir -p data/'//this%dir//'/f')
		call system('rm data/'//this%dir//'/f/*.*')

		nr = this%nt/this%nmod+1
		allocate(this%rho(p%nx,nr))
		allocate(this%phi(p%nx,nr))
		allocate(this%E(p%nx,nr))
		allocate(this%KE(nr))
		allocate(this%PE(nr))
		allocate(this%TE(nr))

		!initRecord: save grid information
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
	end subroutine

	subroutine destroyRecord(this)
		type(history), intent(inout) :: this

		deallocate(this%dir)
		deallocate(this%KE)
		deallocate(this%PE)
		deallocate(this%TE)
		deallocate(this%rho)
		deallocate(this%phi)
		deallocate(this%E)
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
			this%KE(kr+1) = 0.5_mp*p%ms*SUM(w)*p%dx*p%dv
			this%PE(kr+1) = 0.5_mp*p%eps0*SUM(p%E*p%E)*p%dx
			this%TE(kr+1) = this%KE(kr+1) + this%PE(kr+1)

			open(unit=302,file='data/'//this%dir//'/f/'//trim(adjustl(kstr))//		&
					'.bin',status='replace',form='unformatted',access='stream')
			write(302) p%f
			close(302)

			this%rho(:,kr+1)=p%rho
			this%phi(:,kr+1)=p%phi
			this%E(:,kr+1)=p%E
			print *, 'Time: ', k*this%dt, ', KE: ',this%KE(kr+1),', PE: ',this%PE(kr+1),', TE: ',this%TE(kr+1)
		end if
	end subroutine

	subroutine printPlasma(this)
		type(history), intent(in) :: this

		open(unit=303,file='data/'//this%dir//'/E.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/'//this%dir//'/rho.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/'//this%dir//'/phi.bin',status='replace',form='unformatted',access='stream')
		open(unit=306,file='data/'//this%dir//'/KE.bin',status='replace',form='unformatted',access='stream')
		open(unit=307,file='data/'//this%dir//'/PE.bin',status='replace',form='unformatted',access='stream')
		open(unit=308,file='data/'//this%dir//'/TE.bin',status='replace',form='unformatted',access='stream')
		write(303) this%E
		write(304) this%rho
		write(305) this%phi
		write(306) this%KE
		write(307) this%PE
		write(308) this%TE
		close(303)
		close(304)
		close(305)
		close(306)
		close(307)
		close(308)
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
