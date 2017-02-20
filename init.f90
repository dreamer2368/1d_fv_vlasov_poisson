module init

	use parameters
	use MatrixVector
	use modRecord

	implicit none

contains

	subroutine setup()
		!iteration variable
		integer :: i
		!Efield Evaluation
		real(mp) :: rhs(Nx), phi1(Nx-1), phi(Nx), E(Nx)

		!grid setup
		xg = (/ ( i*dx, i=0,Nx-1 ) /)
		vg = (/ ( i*dv, i=-Nv,Nv ) /)

		!reference perturbation amplitude
		A = A0
		B = B0

		!E field evaluation to determine CFL condition
		call initial_condition(f0,A,B,sigma,xg,vg)

		rhs = -1.0_mp/eps0*( qe*integrate_dv(f0,dv) + rho_back )
		call CG_K(phi1,rhs(2:Nx),dx)
		phi(2:Nx) = phi1
		phi(1) = 0.0_mp
		E = - multiplyD(phi,dx)

		!time step parameter
		dt = CFL*MIN( dx/Lv, dv/MAXVAL(ABS(E))/ABS(qe)*me )
		Nt = CEILING(T/dt)
		dt = T/Nt
		Ni = FLOOR(Ti/dt) + 1
		print *, 'T=',T,'Ni=',Ni,', Nt=',Nt,', dt=',dt
		call buildRecord(plasmaRecord,xg,vg,Nt,dx,dv,dt)

	end subroutine

	subroutine initial_condition(f,a,b,r,xg,vg)
		real(mp), intent(in) :: a, b, r				!perturbation amplitude
		real(mp), intent(in) :: xg(:), vg(:)				!mesh grid
		real(mp), intent(out) :: f(size(xg),size(vg))
		integer :: i,j
		real(mp) :: fx(size(xg)), fv(size(vg))

		fx = ( 1.0_mp + a*COS( 2.0_mp*pi*xg/L ) )
		fv = 0.5_mp/sqrt(2.0_mp*pi)/r*( exp( -(vg-v0)**2/2.0_mp/r/r ) + exp( -(vg+v0)**2/2.0_mp/r/r ) )
		do i=1,size(xg)
			do j=1, size(vg)
				f(i,j) = fx(i)*fv(j)
			end do
		end do
		f = f/( SUM(integrate_dv(f,dv))*dx )*N
	end subroutine

end module