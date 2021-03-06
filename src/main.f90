program main

	use testmodules
    use modInputHelper

	implicit none

    character(len=STRING_LENGTH), parameter :: PROJECT_NAME='PASS'
    character(len=STRING_LENGTH) :: filename
    real(mp) :: start, finish, dummy(2)

    ! initiate MPI
    call mpih%buildMPIHandler

	! print to screen
    if( mpih%my_rank .eq. 0 )           &
    	print *, 'calling program main'

    ! Parse options from the input file.
    filename = trim(PROJECT_NAME) // ".inp"
    call parseInputFile(filename)
    print_simulation_detail = getOption('print_simulation_detail',.false.)
    Ng_ = getOption('number_of_grids',1024)
    Time_ = getOption('simulation_time',30.0_mp)
    dir_ = getOption('base_directory','Debye')
    A0_ = getOption('parameters_of_interest/001',1.5_mp)
    A1_ = getOption('parameters_of_interest/002',0.0_mp)

	call cpu_time(start)
!	call sheath
	call debye(.false.,A0_+A1_,Time_,Ng_,dummy,trim(dir_))
!	call twostream
!	call manufactured_solution
!	call debye_sensitivity
!	call BoundaryTest
!	call DNsolverTest
!    call QoI_curve(debye)
!    call QoI_convergence(debye)
	call cpu_time(finish)


	! print to screen
    if( mpih%my_rank .eq. 0) then
    	print *, 'Elapsed time = ',finish-start
	    print *, 'program main...done.'
    end if

    ! finish MPI
    call mpih%destroyMPIHandler

contains

	! You can add custom subroutines/functions here later, if you want

	subroutine sheath
		type(plasma) :: p(2)
		type(circuit) :: c
		type(history) :: r
		integer, parameter :: Nx=256, Nv=128
		real(mp), parameter :: CFL = 0.5_mp
		real(mp), parameter :: Kb = 1.38065E-23, EV_TO_K = 11604.52_mp, eps = 8.85418782E-12
		real(mp), parameter :: Te = 50.0_mp*EV_TO_K, tau = 100.0_mp
		real(mp), parameter :: me = 9.10938215E-31, qe = 1.602176565E-19, mu = 1836
		real(mp), parameter :: n0 = 2.0E14
		real(mp) :: mi, Ti, wp0, lambda0, dt, dx, L
		real(mp) :: ve0, vi0, Lv_e, Lv_i, Time_f
		real(mp) :: A
		integer :: i

		mi = mu*me
		Ti = Te/tau
		wp0 = sqrt(n0*qe*qe/me/eps)
		lambda0 = sqrt(eps*Kb*Te/n0/qe/qe)
		L = 20.0_mp*lambda0

		print *, 'L = ',L,', lambda0 = ',lambda0,' e = lambda/L = ',lambda0/L

		ve0 = sqrt(Kb*Te/me)
		Lv_e = 4.0_mp*ve0
		vi0 = sqrt(Kb*Ti/mi)
		Lv_i = 35.0_mp*vi0
		Time_f = 1.0_mp*L/vi0

		call buildPlasma(p(1),L,Lv_e,Nx,Nv,-qe,me,A0=(/ve0,0.0_mp,0.0_mp/))
		call buildPlasma(p(2),L,Lv_i,Nx,Nv,qe,mi,A0=(/vi0,0.0_mp,0.0_mp/))
		call buildCircuit(c,L,Nx,eps)
		call initial_sheath(p,c,n0,ve0,vi0)

		call buildRecord(r,p,c,Time_f,CFL=CFL,input_dir='sheath',nmod=500)
		call forward_sweep(p,c,r,inputSource=ConstantIon)
		call printPlasma(r)

		call destroyRecord(r)
		call destroyPlasma(p(1))
		call destroyPlasma(p(2))
		call destroyCircuit(c)
	end subroutine

	subroutine debye_sensitivity
		type(plasma) :: p(1),dp(1)
		type(circuit) :: c,dc
		type(history) :: r,dr
		real(mp), parameter :: L = 20.0_mp, Lv = 9.0_mp
		real(mp), parameter :: vT = 1.5_mp, Q = 2.0_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = -1.0_mp, me = 1.0_mp
		integer, parameter :: Nx = 128, Nv = 64
		real(mp), parameter :: T = 150.0_mp, CFL = 0.5_mp

		call buildPlasma(p(1),L,Lv,Nx,Nv,qe,me)
		call buildCircuit(c,L,Nx,eps0)
		call buildPlasma(dp(1),L,Lv,Nx,Nv,qe,me)
		call buildCircuit(dc,L,Nx,eps0)

		call initial_debye(p(1),c,vT,Q)
		call initial_debye_sensitivity(dp(1),dc,vT,'vT')

		call buildRecord(r,p,c,T,CFL=CFL,input_dir='debye',nmod=500)
		call buildRecord(dr,dp,dc,T,dt=r%dt,input_dir='debye/f_A',nmod=500)

		call forward_sensitivity(p,c,r,dp,dc,dr,Screening_distance)

		call printPlasma(r)
		call printPlasma(dr)

		call destroyRecord(r)
		call destroyRecord(dr)
		call destroyPlasma(p(1))
		call destroyPlasma(dp(1))
		call destroyCircuit(c)
		call destroyCircuit(dc)
	end subroutine

	subroutine debye(base,fk,time,Ng,output,dir)
        !<< argument >>
        logical, intent(in) :: base
        real(mp), intent(in) :: fk
        real(mp), intent(in) :: time
        integer, intent(in) :: Ng
        real(mp), intent(out) :: output(2)
        character(len=*), intent(in), optional :: dir

        !<< local variable >>
		type(plasma) :: p(1), dp(1)
		type(circuit) :: c, dc
		type(history) :: r, dr
        real(mp) :: vT, Lv, J, dJ
        character(len=STRING_LENGTH) :: dir_
		real(mp), parameter :: L = 20.0_mp
		real(mp), parameter :: Q = 2.0_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = -1.0_mp, me = 1.0_mp
		real(mp), parameter :: CFL = 0.5_mp
		integer :: Nx, Nv
        vT = fk
        Lv = MAX( 6.0_mp*vT, 9.0_mp )
        if( present(dir) ) then
            dir_ = trim(dir)
        else
            dir_ = 'debyeNg1024'
        end if
        Nx = Ng
        Nv = Ng/2

		call buildPlasma(p(1),L,Lv,Nx,Nv,qe,me)
		call buildCircuit(c,L,Nx,eps0)
		call initial_debye(p(1),c,vT,Q)

        if( base ) then
    		call buildPlasma(dp(1),L,Lv,Nx,Nv,qe,me)
	    	call buildCircuit(dc,L,Nx,eps0)
		    call initial_debye_sensitivity(dp(1),dc,vT,'vT')
        end if

		call buildRecord(r,p,c,time,CFL=CFL,input_dir=dir_,nmod=500)
        if( base ) call buildRecord(dr,dp,dc,time,dt=r%dt,input_dir=trim(dir_)//'/f_A',nmod=500)

        if( .not. base ) then
    		call forward_sweep(p,c,r,Screening_Distance)
        else
		    call forward_sensitivity(p,c,r,dp,dc,dr,Screening_distance)
        end if

        if( print_simulation_detail ) then
            call printPlasma(r)
            if( base ) call printPlasma(dr)
        end if

        J = sum( r%j )*r%dt/time
        dJ = sum( dr%j )*dr%dt/time
        if( base ) then
            output = (/ J, dJ /)
        else
            output = (/ vT, J /)
        end if

		call destroyRecord(r)
		call destroyPlasma(p(1))
		call destroyCircuit(c)
        if( base ) then
		    call destroyRecord(dr)
    		call destroyPlasma(dp(1))
    		call destroyCircuit(dc)
        end if
	end subroutine

	subroutine twostream
		type(plasma) :: p(1)
		type(circuit) :: c
		type(history) :: r
		real(mp), parameter :: L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
		real(mp), parameter :: Lv = 1.0_mp
		real(mp), parameter :: a=0.0001_mp, v0 = 0.2_mp, vT = 0.01_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = -1.0_mp, me = 1.0_mp
		integer, parameter :: Nx = 128, Nv = 64
		real(mp), parameter :: T=60.0_mp, CFL = 0.5_mp

		call buildPlasma(p(1),L,Lv,Nx,Nv,qe,me)
		call buildCircuit(c,L,Nx,eps0)
		call initial_twostream(p(1),c,a,v0,vT)
		call buildRecord(r,p,c,T,CFL=CFL,input_dir='twostream',nmod=20)
		call forward_sweep(p,c,r)
		call printPlasma(r)
		call destroyRecord(r)
		call destroyPlasma(p(1))
		call destroyCircuit(c)
	end subroutine

	subroutine QoI_curve(problem)
		real(mp) ::  vT_min, vT_max
        real(mp), allocatable :: vT(:)
		integer :: Nsample, Ng
		real(mp) :: Time
		integer :: i, thefile, idx, input
		character(len=100):: dir, filename
		interface
			subroutine problem(base,fk,time,Ng,output,dir)
				use modPlasma
				use modCircuit
				use modRecord
                logical, intent(in) :: base
				real(mp), intent(in) :: fk, time
                integer, intent(in) :: Ng
				real(mp), intent(out) :: output(2)
				character(len=*), optional, intent(in) :: dir
				type(plasma) :: p
                type(circuit) :: c
				type(history) :: r
			end subroutine
		end interface
        Time = getOption('QoI_curve/time',150.0_mp)
        dir = getOption('QoI_curve/directory','Debye_curve')
        filename = getOption('QoI_curve/filename','J.bin')
        vT_min = getOption('QoI_curve/min_parameter_value',1.49_mp)
        vT_max = getOption('QoI_curve/max_parameter_value',1.51_mp)
        Nsample = getOption('QoI_curve/number_of_sample',1001)
        Ng = getOption('QoI_curve/number_of_grids',1024)

        allocate(vT(Nsample))
		vT = (/ ((vT_max-vT_min)*(i-1)/(Nsample-1)+vT_min,i=1,Nsample) /)

		call allocateBuffer(Nsample,2,mpih)
        thefile = MPIWriteSetup(mpih,'data/'//trim(dir),filename)

		do i=1,mpih%sendcnt
            call problem( .false., vT(mpih%displc(mpih%my_rank)+i),         &
                            Time, Ng, mpih%writebuf,                        &
                            trim(dir)//'/'//trim(adjustl(mpih%rank_str)) )

            call MPI_FILE_WRITE(thefile, mpih%writebuf, 2, MPI_DOUBLE, & 
                                MPI_STATUS_IGNORE, mpih%ierr)
            call MPI_FILE_SYNC(thefile,mpih%ierr)

            print ('(A,I5,A,I5,A,F8.3,A,F8.3)'), 'Rank-',mpih%my_rank,      &
                                                 ' Sample-',i,              &
                                                 ', vT=',mpih%writebuf(1),  &
                                                 ', J=',mpih%writebuf(2)
		end do

        deallocate(vT)

        call MPI_FILE_CLOSE(thefile, mpih%ierr)            
    end subroutine

	subroutine QoI_convergence(problem)
		real(mp) ::  vT0, f0
        real(mp), allocatable :: vT(:)
		integer :: Nsample, Ng
		real(mp) :: Time
		integer :: i, thefile, idx, input
		character(len=100):: dir, filename
        logical :: base
		interface
			subroutine problem(base,fk,time,Ng,output,dir)
				use modPlasma
				use modCircuit
				use modRecord
                logical, intent(in) :: base
				real(mp), intent(in) :: fk, time
                integer, intent(in) :: Ng
				real(mp), intent(out) :: output(2)
				character(len=*), optional, intent(in) :: dir
				type(plasma) :: p
                type(circuit) :: c
				type(history) :: r
			end subroutine
		end interface
        Time = getOption('QoI_convergence/time',150.0_mp)
        dir = getOption('QoI_convergence/directory','Debye_convergence')
        filename = getOption('QoI_convergence/filename','J.bin')
        vT0 = getOption('QoI_convergence/evaluation_point',1.5_mp)
        f0 = getOption('QoI_convergence/initial_perturbation',1.0_mp)
        Nsample = getOption('QoI_convergence/number_of_sample',70)
        Ng = getOption('QoI_convergence/number_of_grids',1024)

        allocate(vT(Nsample))
		vT = (/ (f0 - 0.25_mp*(i-2),i=1,Nsample) /)
        vT = 10.0_mp**vT
        vT(1) = 0.0_mp
        vT = vT0 + vT

		call allocateBuffer(Nsample,2,mpih)
        thefile = MPIWriteSetup(mpih,'data/'//trim(dir),filename)

		do i=1,mpih%sendcnt
            idx = mpih%displc(mpih%my_rank)+i
            base = (idx.eq.1)
            call problem( base, vT(idx),                                    &
                          Time, Ng, mpih%writebuf,                          &
                          trim(dir)//'/'//trim(adjustl(mpih%rank_str)) )

            call MPI_FILE_WRITE(thefile, mpih%writebuf, 2, MPI_DOUBLE, & 
                                MPI_STATUS_IGNORE, mpih%ierr)
            call MPI_FILE_SYNC(thefile,mpih%ierr)

            if( base ) then
                print ('(A,I5,A,I5,3(A,F8.3))'), 'Rank-',mpih%my_rank,      &
                                                 ' Sample-',i,              &
                                                 ', vT=',vT(1),             &
                                                 ', J=',mpih%writebuf(1),   &
                                                 ', dJ=',mpih%writebuf(2)
            else
                print ('(A,I5,A,I5,2(A,F8.3))'), 'Rank-',mpih%my_rank,      &
                                                 ' Sample-',i,              &
                                                 ', vT=',mpih%writebuf(1),  &
                                                 ', J=',mpih%writebuf(2)
            end if
		end do

        deallocate(vT)

        call MPI_FILE_CLOSE(thefile, mpih%ierr)            
    end subroutine

end program
