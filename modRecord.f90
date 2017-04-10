module modRecord

	use modPlasma
	use modCircuit

	implicit none

	type history
		integer :: nt, nmod = 20
		real(mp) :: CFL, dt, T

		character(len=:), allocatable :: dir

		real(mp), allocatable :: KE(:), PE(:), TE(:)

		real(mp), allocatable :: rho(:,:), phi(:,:), E(:,:)

		real(mp), allocatable :: j(:)

		real(mp), allocatable :: cpt_time(:,:)
		real(mp) :: cpt_temp(6)
	end type

contains

	subroutine buildRecord(this,p,c,T,CFL,dt,input_dir,nmod)
		type(history), intent(out) :: this
		type(plasma), intent(inout) :: p
		type(circuit), intent(inout) :: c
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
			call Efield_record(p,c)
			if( MAXVAL(ABS(c%E)).ne.0.0_mp ) then
				print *, 'dx/v=',p%dx/p%Lv,', dv/acc=',p%dv/MAXVAL(ABS(c%E))/ABS(p%qs)*p%ms
				this%dt = CFL*MIN( p%dx/p%Lv, p%dv/MAXVAL(ABS(c%E))/ABS(p%qs)*p%ms )
			else
				this%dt = CFL*p%dx/p%Lv
			end if
			this%nt = CEILING(T/this%dt)
		!measure CFL according to timestep size
		elseif( PRESENT(dt) ) then
			this%nt = CEILING(T/dt)
			this%dt = T/this%nt
			call Efield_record(p,c)
			if( MAXVAL(ABS(c%E)).ne.0.0_mp ) then
				this%CFL = MAX( p%Lv*this%dt/p%dx, MAXVAL(ABS(c%E))*ABS(p%qs)/p%ms*this%dt/p%dv )
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

		!Quantity of Interest: save every timestep
		allocate(this%j(this%nt))
		this%j=0.0_mp

		!Computation time measurement
		allocate(this%cpt_time(6,nr))
		this%cpt_time = 0.0_mp
		this%cpt_temp = 0.0_mp
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
		deallocate(this%j)
	end subroutine

	subroutine recordPlasma(this,p,c,k)
		type(plasma), intent(in) :: p
		type(circuit), intent(in) :: c
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
			this%PE(kr+1) = 0.5_mp*c%eps0*SUM(c%E*c%E)*c%dx
			this%TE(kr+1) = this%KE(kr+1) + this%PE(kr+1)

			open(unit=302,file='data/'//this%dir//'/f/'//trim(adjustl(kstr))//		&
					'.bin',status='replace',form='unformatted',access='stream')
			write(302) p%f
			close(302)

			this%rho(:,kr+1)=c%rho
			this%phi(:,kr+1)=c%phi
			this%E(:,kr+1)=c%E
			print *, 'Time: ', k*this%dt, ', KE: ',this%KE(kr+1),', PE: ',this%PE(kr+1),', TE: ',this%TE(kr+1)
			print *, 'MEAN(j): ', SUM(this%j(1:k))/k

			!Computation time measurement
			this%cpt_time(:,kr+1) = this%cpt_temp
			this%cpt_temp = 0.0_mp
		end if
	end subroutine

	subroutine printPlasma(this)
		type(history), intent(in) :: this
		real(mp) :: total,mean,pct

		open(unit=303,file='data/'//this%dir//'/E.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/'//this%dir//'/rho.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/'//this%dir//'/phi.bin',status='replace',form='unformatted',access='stream')
		open(unit=306,file='data/'//this%dir//'/KE.bin',status='replace',form='unformatted',access='stream')
		open(unit=307,file='data/'//this%dir//'/PE.bin',status='replace',form='unformatted',access='stream')
		open(unit=308,file='data/'//this%dir//'/TE.bin',status='replace',form='unformatted',access='stream')
		open(unit=309,file='data/'//this%dir//'/j.bin',status='replace',form='unformatted',access='stream')
		open(unit=310,file='data/'//this%dir//'/cpt_time.bin',status='replace',form='unformatted',access='stream')
		write(303) this%E
		write(304) this%rho
		write(305) this%phi
		write(306) this%KE
		write(307) this%PE
		write(308) this%TE
      write(309) this%j
		write(310) this%cpt_time
		close(303)
		close(304)
		close(305)
		close(306)
		close(307)
		close(308)
      close(309)
		close(310)

701	FORMAT	(A, F10.3,'	',F10.3,'	', F10.2,'%')
		if( SUM(this%cpt_time(5,:)).eq.0.0_mp ) then
			open(unit=301,file='data/'//this%dir//'/original_cpt_summary.dat',status='replace')
			write(301,*) 'Step	Total	Mean	Percentage'
			print *, "================ Computation Time Summary ==================================="
			print *, "Original simulation	   	     Total            Mean	 Percentage	"
			total = SUM(this%cpt_time(1,:))*this%nmod
			mean = total/this%nt
			pct = total/this%nmod/SUM(this%cpt_time)*100.0_mp
			print 701, "TransportSpace			", total, mean, pct
			write(301,701) 'Transport-Space	', total, mean, pct
			total = SUM(this%cpt_time(2,:))*this%nmod
			mean = total/this%nt
			pct = total/this%nmod/SUM(this%cpt_time)*100.0_mp
			print 701, "Efield				", total, mean, pct
			write(301,701) 'Efield	', total, mean, pct
			total = SUM(this%cpt_time(3,:))*this%nmod
			mean = total/this%nt
			pct = total/this%nmod/SUM(this%cpt_time)*100.0_mp
			print 701, "TransportVelocity		", total, mean, pct
			write(301,701) 'Transport-Velocity	', total, mean, pct
         total = SUM(this%cpt_time(6,:))*this%nmod
         mean = total/this%nt
         pct = total/this%nmod/SUM(this%cpt_time)*100.0_mp
         print 701, "Record			            ", total, mean, pct
			write(301,701) 'Record	', total, mean, pct
			print *, "============================================================================="
			close(301)
		else
			open(unit=301,file='data/'//this%dir//'/sensitivity_cpt_summary.dat',status='replace')
			write(301,*) 'Step	Total	Mean	Percentage'
			print *, "================ Computation Time Summary ==================================="
			print *, "Sensitivity simulation	  	     Total            Mean   	 Percentage	"
			total = SUM(this%cpt_time(1,:))*this%nmod
			mean = total/this%nt
			pct = total/this%nmod/SUM(this%cpt_time)*100.0_mp
			print 701, "TransportSpace1/2		", total, mean, pct
			write(301,701) 'Transport-Space1/2	', total, mean, pct
			total = SUM(this%cpt_time(2,:))*this%nmod
			mean = total/this%nt
			pct = total/this%nmod/SUM(this%cpt_time)*100.0_mp
			print 701, "Efield				", total, mean, pct
			write(301,701) 'Efield	', total, mean, pct
			total = SUM(this%cpt_time(3,:))*this%nmod
			mean = total/this%nt
			pct = total/this%nmod/SUM(this%cpt_time)*100.0_mp
			print 701, "TransportVelocity		", total, mean, pct
			write(301,701) 'Transport-Velocity	', total, mean, pct
			total = SUM(this%cpt_time(4,:))*this%nmod
			mean = total/this%nt
			pct = total/this%nmod/SUM(this%cpt_time)*100.0_mp
			print 701, "TransportSpace1/2		", total, mean, pct
			write(301,701) 'Transport-Space1/2	', total, mean, pct
			total = SUM(this%cpt_time(5,:))*this%nmod
			mean = total/this%nt
			pct = total/this%nmod/SUM(this%cpt_time)*100.0_mp
			print 701, "TransportSource			", total, mean, pct
			write(301,701) 'TransportSource	', total, mean, pct
         total = SUM(this%cpt_time(6,:))*this%nmod
         mean = total/this%nt
         pct = total/this%nmod/SUM(this%cpt_time)*100.0_mp
         print 701, "Record  	  		        ", total, mean, pct
			write(301,701) 'Record	', total, mean, pct
			print *, "============================================================================="
			close(301)
		end if
	end subroutine

	subroutine Efield_record(p,c)
		type(plasma), intent(inout) :: p
		type(circuit), intent(inout) :: c
		real(mp) :: dx, dv
		real(mp) :: rhs(p%nx-1)
		real(mp) :: phi1(p%nx-1)
		dx = p%dx
		dv = p%dv

		c%rho = p%qs*integrate_dv(p%f,dv) + c%rho_back
		rhs = -1.0_mp/c%eps0*c%rho(2:p%nx)
		call CG_K(phi1,rhs,dx)
		c%phi(2:p%nx) = phi1
		c%phi(1) = 0.0_mp
		c%E = - multiplyD(c%phi,dx)
	end subroutine

end module
