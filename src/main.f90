program main

	use testmodules

	implicit none

	real(mp) :: start, finish

	! print to screen
	print *, 'calling program main'
	call cpu_time(start)
!	call sheath
!	call debye
!	call twostream
!	call manufactured_solution
	call debye_sensitivity
!	call BoundaryTest
!	call DNsolverTest
	call cpu_time(finish)
	print *, 'Elapsed time = ',finish-start

	! print to screen
	print *, 'program main...done.'

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

	subroutine debye
		type(plasma) :: p(1)
		type(circuit) :: c
		type(history) :: r
		real(mp), parameter :: L = 20.0_mp, Lv = 9.0_mp
		real(mp), parameter :: vT = 1.5_mp, Q = 2.0_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = -1.0_mp, me = 1.0_mp
		integer, parameter :: Nx = 512, Nv = 256
		real(mp), parameter :: T = 150.0_mp, CFL = 0.5_mp

		call buildPlasma(p(1),L,Lv,Nx,Nv,qe,me)
		call buildCircuit(c,L,Nx,eps0)
		call initial_debye(p(1),c,vT,Q)
		call buildRecord(r,p,c,T,CFL=CFL,input_dir='debye',nmod=300)
		call forward_sweep(p,c,r)
		call printPlasma(r)
		call destroyRecord(r)
		call destroyPlasma(p(1))
		call destroyCircuit(c)
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

end program
