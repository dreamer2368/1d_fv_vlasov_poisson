module modRecord

	use modPlasma
	use modCircuit

	implicit none

	type history
		integer :: nt, nr, nmod = 20
		real(mp) :: CFL, dt, T

		integer :: ns

		character(len=:), allocatable :: dir

		real(mp), allocatable :: KE(:,:), PE(:), TE(:)

		real(mp), allocatable :: rho(:,:), phi(:,:), E(:,:)

		real(mp), allocatable :: j(:)

		real(mp), allocatable :: cpt_time(:,:)
		real(mp) :: cpt_temp(6)
	end type

contains

	subroutine buildRecord(this,p,c,T,Emax,CFL,dt,input_dir,nmod)
		type(history), intent(out) :: this
		type(plasma), intent(inout) :: p(:)
		type(circuit), intent(inout) :: c
		real(mp), intent(in) :: T
		real(mp), intent(in), optional :: Emax, CFL, dt
		character(len=*), intent(in), optional :: input_dir
		integer, intent(in), optional :: nmod
		integer :: ns, nr, Nv(SIZE(p))
		real(mp), dimension(SIZE(p)) :: qs, ms, Lv, dv
		real(mp), allocatable :: A(:,:)
		real(mp) :: nu_x, nu_v
		integer :: i
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

		ns = SIZE(p)
		this%ns = ns
		allocate(A(SIZE(p(1)%A),ns))
		do i=1,ns
			qs(i) = p(i)%qs
			ms(i) = p(i)%ms
			Lv(i) = p(i)%Lv
			Nv(i) = p(i)%nv
			dv(i) = p(i)%dv
			A(:,i) = p(i)%A
		end do

		this%T = T
		!set timestep size according to CFL criterion with initial condition
		if( PRESENT(CFL) ) then
			this%CFL = CFL
			if( PRESENT(Emax) ) then
				nu_x = c%dx/MAXVAL(Lv)
				nu_v = MAXVAL(dv/ABS(qs)*ms)/Emax
				print *, 'dx/v=',nu_x,', dv/acc=',nu_v
				this%dt = CFL*MIN( nu_x,nu_v )
			else
				nu_x = c%dx/MAXVAL(Lv)
				this%dt = CFL*nu_x
			end if
			this%nt = CEILING(T/this%dt)
		!measure CFL according to timestep size
		elseif( PRESENT(dt) ) then
			this%nt = CEILING(T/dt)
			this%dt = T/this%nt
			if( PRESENT(Emax) ) then
				nu_x = c%dx/MAXVAL(Lv)
				nu_v = MAXVAL(dv/ABS(qs)*ms)/Emax
				this%CFL = MAX( this%dt/nu_x, this%dt/nu_v )
			else
				nu_x = c%dx/MAXVAL(Lv)
				this%CFL = this%dt/nu_x
			end if
		else
			print *, 'ERROR: required to input either CFL or dt'
			stop
		end if

		print *, 'T=',this%T,', CFL=',this%CFL,', Nt=',this%nt,', dt=',this%dt
		call system('mkdir -p data/'//this%dir//'/f')
		call system('rm data/'//this%dir//'/f/*.*')

		nr = this%nt/this%nmod+1
		this%nr = nr
		allocate(this%rho(c%nx,nr))
		allocate(this%phi(c%nx,nr))
		allocate(this%E(c%nx,nr))
		allocate(this%KE(ns,nr))
		allocate(this%PE(nr))
		allocate(this%TE(nr))

		!initRecord: save grid information
		open(unit=301,file='data/'//this%dir//'/record_readme.out',status='replace')
		write(301,*) this%nt
		write(301,*) this%nmod
		write(301,*) this%dt
		write(301,*) this%nr
		write(301,*) this%ns
		write(301,*) c%nx
		write(301,*) c%Lx
		close(301)

		open(unit=301,file='data/'//this%dir//'/xg.bin',status='replace',form='unformatted',access='stream')
		write(301) c%xg
		close(301)

		open(unit=301,file='data/'//this%dir//'/Lv.bin',status='replace',form='unformatted',access='stream')
		write(301) Lv
		close(301)

		open(unit=301,file='data/'//this%dir//'/Nv.bin',status='replace',form='unformatted',access='stream')
		write(301) Nv
		close(301)

		open(unit=301,file='data/'//this%dir//'/vg.bin',status='replace',form='unformatted',access='stream')
		do i=1,ns
			write(301) p(i)%vg
		end do
		close(301)

		!Quantity of Interest: save every timestep
		allocate(this%j(this%nt))
		this%j=0.0_mp

		!Computation time measurement
		allocate(this%cpt_time(6,nr))
		this%cpt_time = 0.0_mp
		this%cpt_temp = 0.0_mp

		deallocate(A)
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
		type(plasma), intent(in) :: p(:)
		type(circuit), intent(in) :: c
		type(history), intent(inout) :: this
		integer, intent(in) :: k
		integer :: j, kr
		character(len=100) :: kstr, jstr
		integer :: i
		real(mp) :: w(c%nx)

		if( (this%nmod.eq.1) .or. (mod(k,this%nmod).eq.0) ) then
			kr = merge(k,k/this%nmod,this%nmod.eq.1)
			write(kstr,*) kr

			do j=1,this%ns
				write(jstr,*) j
				w = 0.0_mp
				do i=1,2*p(j)%nv
					w = w + 0.5_mp*( p(j)%f(:,i)+p(j)%f(:,i+1) )*( 0.5_mp*(p(j)%vg(i)+p(j)%vg(i+1)) )**2
				end do
				this%KE(j,kr+1) = 0.5_mp*p(j)%ms*SUM(w)*p(j)%dx*p(j)%dv

				open(unit=302,file='data/'//this%dir//'/f/'	&
						//trim(adjustl(jstr))//'_'//trim(adjustl(kstr))//'.bin',	&
						status='replace',form='unformatted',access='stream')
				write(302) p(j)%f
				close(302)
			end do
			this%PE(kr+1) = 0.5_mp*c%eps0*SUM(c%E*c%E)*c%dx
			this%TE(kr+1) = SUM( this%KE(:,kr+1) ) + this%PE(kr+1)

			this%rho(:,kr+1)=c%rho
			this%phi(:,kr+1)=c%phi
			this%E(:,kr+1)=c%E
			print *, 'Time: ', k*this%dt, ', KE: ',SUM( this%KE(:,kr+1) ),', PE: ',this%PE(kr+1),', TE: ',this%TE(kr+1)
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
         total = SUM(this%cpt_time(4,:))*this%nmod
         mean = total/this%nt
         pct = total/this%nmod/SUM(this%cpt_time)*100.0_mp
         print 701, "Source			        ", total, mean, pct
			write(301,701) 'Source	', total, mean, pct
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

end module
