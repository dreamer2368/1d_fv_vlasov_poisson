program main

	use timeStep
	use init
	use modRecord

	implicit none

	real(mp) :: start, finish

	! print to screen
	print *, 'calling program main'

	call cpu_time(start)
	call setup
	call initial_condition(f0,A,B,sigma,xg,vg)
	call buildPlasma(twostream,xg,vg,f0)
	call forward_sweep(twostream,plasmaRecord)
	call destroyPlasma(twostream)
	call destroyRecord(plasmaRecord)
	call cpu_time(finish)
	print *, 'Elapsed time = ',finish-start

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

end program
