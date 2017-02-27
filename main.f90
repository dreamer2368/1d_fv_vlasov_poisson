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
	call debye_sensitivity
	call cpu_time(finish)
	print *, 'Elapsed time = ',finish-start

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want
	subroutine debye_sensitivity
		type(plasma) :: p,dp
		type(history) :: r,dr
		real(mp), parameter :: L = 20.0_mp, Lv = 9.0_mp
		real(mp), parameter :: vT = 1.5_mp, Q = 2.0_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = -1.0_mp, me = 1.0_mp
		integer, parameter :: Nx = 512, Nv = 256
		real(mp), parameter :: T = 0.05_mp, CFL = 0.5_mp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me,eps0)
		call buildPlasma(dp,L,Lv,Nx,Nv,qe,me,eps0)

		call initial_debye(p,vT,Q)
		call initial_debye_sensitivity(dp,vT)

		call buildRecord(r,p,T,CFL=CFL,input_dir='debye',nmod=1)
		call buildRecord(dr,dp,T,dt=r%dt,input_dir='debye/f_A',nmod=1)

		call forward_sensitivity(p,r,dp,dr)

		call printPlasma(r)
		call printPlasma(dr)

		call destroyRecord(r)
		call destroyRecord(dr)
		call destroyPlasma(p)
		call destroyPlasma(dp)
	end subroutine

	subroutine debye
		type(plasma) :: p
		type(history) :: r
		real(mp), parameter :: L = 20.0_mp, Lv = 9.0_mp
		real(mp), parameter :: vT = 1.5_mp, Q = 2.0_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = -1.0_mp, me = 1.0_mp
		integer, parameter :: Nx = 512, Nv = 256
		real(mp), parameter :: T = 150.0_mp, CFL = 0.5_mp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me,eps0)
		call initial_debye(p,vT,Q)
		call buildRecord(r,p,T,CFL=CFL,input_dir='debye',nmod=300)
		call forward_sweep(p,r)
		call printPlasma(r)
		call destroyRecord(r)
		call destroyPlasma(p)
	end subroutine

	subroutine twostream
		type(plasma) :: p
		type(history) :: r
		real(mp), parameter :: L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
		real(mp), parameter :: Lv = 1.0_mp
		real(mp), parameter :: a=0.0001_mp, v0 = 0.2_mp, vT = 0.01_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = -1.0_mp, me = 1.0_mp
		integer, parameter :: Nx = 128, Nv = 64
		real(mp), parameter :: T=80.0_mp, CFL = 0.5_mp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me,eps0)
		call initial_twostream(p,a,v0,vT)
		call buildRecord(r,p,T,CFL=CFL,input_dir='twostream',nmod=20)
		call forward_sweep(p,r)
		call printPlasma(r)
		call destroyRecord(r)
		call destroyPlasma(p)
	end subroutine

end program
