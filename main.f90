program main

	use testmodules

	implicit none

	real(mp) :: start, finish

	! print to screen
	print *, 'calling program main'
	call cpu_time(start)
!	call debye
!	call twostream
!	call manufactured_solution
!	call debye_sensitivity
	call BoundaryTest
	call cpu_time(finish)
	print *, 'Elapsed time = ',finish-start

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want
	subroutine debye_sensitivity
		type(plasma) :: p,dp
		type(circuit) :: c,dc
		type(history) :: r,dr
		real(mp), parameter :: L = 20.0_mp, Lv = 9.0_mp
		real(mp), parameter :: vT = 1.5_mp, Q = 2.0_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = -1.0_mp, me = 1.0_mp
		integer, parameter :: Nx = 512, Nv = 256
		real(mp), parameter :: T = 150.0_mp, CFL = 0.5_mp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me)
		call buildCircuit(c,L,Nx,eps0)
		call buildPlasma(dp,L,Lv,Nx,Nv,qe,me)
		call buildCircuit(dc,L,Nx,eps0)

		call initial_debye(p,c,vT,Q)
		call initial_debye_sensitivity(dp,dc,vT,'qp')

		call buildRecord(r,p,c,T,CFL=CFL,input_dir='debye',nmod=500)
		call buildRecord(dr,dp,dc,T,dt=r%dt,input_dir='debye/f_A',nmod=500)

		call forward_sensitivity(p,c,r,dp,dc,dr,Screening_distance)

		call printPlasma(r)
		call printPlasma(dr)

		call destroyRecord(r)
		call destroyRecord(dr)
		call destroyPlasma(p)
		call destroyPlasma(dp)
		call destroyCircuit(c)
		call destroyCircuit(dc)
	end subroutine

	subroutine debye
		type(plasma) :: p
		type(circuit) :: c
		type(history) :: r
		real(mp), parameter :: L = 20.0_mp, Lv = 9.0_mp
		real(mp), parameter :: vT = 1.5_mp, Q = 2.0_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = -1.0_mp, me = 1.0_mp
		integer, parameter :: Nx = 512, Nv = 256
		real(mp), parameter :: T = 150.0_mp, CFL = 0.5_mp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me)
		call buildCircuit(c,L,Nx,eps0)
		call initial_debye(p,c,vT,Q)
		call buildRecord(r,p,c,T,CFL=CFL,input_dir='debye',nmod=300)
		call forward_sweep(p,c,r)
		call printPlasma(r)
		call destroyRecord(r)
		call destroyPlasma(p)
		call destroyCircuit(c)
	end subroutine

	subroutine twostream
		type(plasma) :: p
		type(circuit) :: c
		type(history) :: r
		real(mp), parameter :: L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
		real(mp), parameter :: Lv = 1.0_mp
		real(mp), parameter :: a=0.0001_mp, v0 = 0.2_mp, vT = 0.01_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = -1.0_mp, me = 1.0_mp
		integer, parameter :: Nx = 128, Nv = 64
		real(mp), parameter :: T=60.0_mp, CFL = 0.5_mp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me)
		call buildCircuit(c,L,Nx,eps0)
		call initial_twostream(p,c,a,v0,vT)
		call buildRecord(r,p,c,T,CFL=CFL,input_dir='twostream',nmod=20)
		call forward_sweep(p,c,r)
		call printPlasma(r)
		call destroyRecord(r)
		call destroyPlasma(p)
		call destroyCircuit(c)
	end subroutine

end program
