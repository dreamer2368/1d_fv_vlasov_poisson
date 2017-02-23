module timeStep

	use MatrixVector
	use modPlasma
	use modRecord
	use Limiter

	implicit none

contains

	subroutine forward_sweep(p,r)
		type(plasma), intent(inout) :: p
		type(history), intent(inout) :: r
		integer :: i

		call initRecord(r,p)
		do i=1,r%nt
			call updatePlasma(p,r%dt)
			call recordPlasma(r,p,i)
!			if( MOD(i,1000) == 0 ) then
!				print *, i
!			end if
		end do
		call closeRecord
	end subroutine

	subroutine updatePlasma(this,dt)					!each time step
		type(plasma), intent(inout) :: this
		real(mp), intent(in) :: dt

		call transportSpace(this,dt)
!		call transportSpace(this,0.5_mp*dt,this%dx,this%dv)
		call Efield(this)
		call transportVelocity(this,dt)
!		call transportSpace(this,0.5_mp*dt,this%dx,this%dv)
	end subroutine

	subroutine transportSpace(this,h)
		type(plasma), intent(inout) :: this
		real(mp), intent(in) :: h					!time step
		real(mp) :: dx, dv
		integer :: i,j,nx
		real(mp) :: nu, theta_p, theta_m
		real(mp) :: fp1,f0,fm1,f2
		real(mp), dimension(this%nx,2*this%nv+1) :: newf
		dx = this%dx
		dv = this%dv
		newf = 0.0_mp

		nx = this%nx

		do i=1,2*this%nv+1
			nu = this%vg(i)*h/dx
			if( this%vg(i)>=0 )	then
				do j=1,nx
					fp1 = this%f( MODULO(j,nx)+1, i)
					f0 = this%f(j,i)
					fm1 = this%f( MODULO(j-2,nx)+1, i)
					f2 = this%f( MODULO(j-3,nx)+1, i)
					theta_m = (fm1-f2)/(f0-fm1)
					theta_p = (f0-fm1)/(fp1-f0)
					newf(j,i) = f0 - nu*(f0-fm1) - 0.5_mp*nu*(1.0_mp-nu)*( FluxLimiter(theta_p,'SB')*(fp1-f0)	&
																									- FluxLimiter(theta_m,'SB')*(f0-fm1) )
				end do
			else
				do j=1,nx
					f2 = this%f( MODULO(j+1,nx)+1, i)
					fp1 = this%f( MODULO(j,nx)+1, i)
					f0 = this%f(j,i)
					fm1 = this%f( MODULO(j-2,nx)+1, i)
					theta_m = (fp1-f0)/(f0-fm1)
					theta_p = (f2-fp1)/(fp1-f0)
					newf(j,i) = f0 - nu*(fp1-f0) + 0.5_mp*nu*(1.0_mp+nu)*( FluxLimiter(theta_p,'SB')*(fp1-f0)	&
																									- FluxLimiter(theta_m,'SB')*(f0-fm1) )
				end do
			end if
		end do
		this%f = newf
	end subroutine

	subroutine transportVelocity(this,h)
		type(plasma), intent(inout) :: this
		real(mp), intent(in) :: h					!time step
		real(mp) :: dx, dv
		integer :: i,j
		real(mp) ::  acc, nu, theta_p, theta_m
		real(mp) :: fp1,f0,fm1,f2
		real(mp), dimension(this%nx,2*this%nv+1) :: newf
		integer :: NV
		dx = this%dx
		dv = this%dv
		NV = 2*this%nv+1
		newf = 0.0_mp

		acc = 0.0_mp
		do i=1,this%nx
			acc = this%qs*this%E(i)/this%ms
			nu = acc*h/dv
			if( acc>=0 )	then
				do j=1,NV
					fp1 = MERGE( this%f(i,j+1), 0.0_mp, j<NV )
					f0 = this%f(i,j)
					fm1 = MERGE( this%f(i,j-1), 0.0_mp, j>1 )
					f2 = MERGE( this%f(i,j-2), 0.0_mp, j>2 )
					theta_m = (fm1-f2)/(f0-fm1)
					theta_p = (f0-fm1)/(fp1-f0)
					newf(i,j) = f0 - nu*(f0-fm1) - 0.5_mp*nu*(1.0_mp-nu)*( FluxLimiter(theta_p,'SB')*(fp1-f0)	&
																									- FluxLimiter(theta_m,'SB')*(f0-fm1) )
				end do
			else
				do j=1,NV
					f2 = MERGE( this%f(i,j+2), 0.0_mp, j<NV-1 )
					fp1 = MERGE( this%f(i,j+1), 0.0_mp, j<NV )
					f0 = this%f(i,j)
					fm1 = MERGE( this%f(i,j-1), 0.0_mp, j>1 )
					theta_m = (fp1-f0)/(f0-fm1)
					theta_p = (f2-fp1)/(fp1-f0)
					newf(i,j) = f0 - nu*(fp1-f0) + 0.5_mp*nu*(1.0_mp+nu)*( FluxLimiter(theta_p,'SB')*(fp1-f0)	&
																									- FluxLimiter(theta_m,'SB')*(f0-fm1) )
				end do
			end if
		end do
		this%f = newf
	end subroutine

	subroutine Efield(this)
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
