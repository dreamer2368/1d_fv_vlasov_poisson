program main

	use timeStep
	use init
	use modRecord

	implicit none

	real(mp) :: start, finish

	! print to screen
	print *, 'calling program main'

	call cpu_time(start)
!	call setup
!	call initial_condition(f0,A,B,sigma,xg,vg)
!	call buildPlasma(twostream,xg,vg,f0)
!	call forward_sweep(twostream,plasmaRecord)
!	call destroyPlasma(twostream)
!	call destroyRecord(plasmaRecord)
	call twostream
	call cpu_time(finish)
	print *, 'Elapsed time = ',finish-start

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

	subroutine twostream
		type(plasma) :: p
		type(history) :: r
		real(mp), parameter :: L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
		real(mp), parameter :: Lv = 1.0_mp
		real(mp), parameter :: a=0.001_mp, v0 = 0.2_mp, vT = 0.01_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = -1.0_mp, me = 1.0_mp
		integer, parameter :: Nx = 128, Nv = 64
		real(mp), parameter :: T=80.0_mp, CFL = 0.5_mp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me,eps0)
		call initial_twostream(p,a,v0,vT)
		call buildRecord(r,p,CFL,T,'twostream',20)
		call forward_sweep(p,r)
		call destroyRecord(r)
		call destroyPlasma(p)
	end subroutine

end program
