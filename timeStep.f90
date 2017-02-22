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
		integer :: i,j,pnx
		real(mp) :: FxL(this%nx), FxR(this%nx)
		dx = this%dx
		dv = this%dv

		pnx = this%nx

		do i=1,2*this%nv+1
			FxL = 0.0_mp
			FxR = 0.0_mp
			if( this%vg(i)>=0 )	then
				do j=2,this%nx
					FxL(j) = this%vg(i)*this%f(j-1,i) + 0.5_mp*this%vg(i)*( 1.0_mp - this%vg(i)*h/dx )*spaceLimiter(this,j-1,i,'SB')
					FxR(j) = this%vg(i)*this%f(j,i) + 0.5_mp*this%vg(i)*( 1.0_mp - this%vg(i)*h/dx )*spaceLimiter(this,j,i,'SB')
				end do
				FxL(1) = this%vg(i)*this%f(this%nx,i) + 0.5_mp*this%vg(i)*( 1.0_mp - this%vg(i)*h/dx )*spaceLimiter(this,this%nx,i,'SB')
				FxR(1) = this%vg(i)*this%f(1,i) + 0.5_mp*this%vg(i)*( 1.0_mp - this%vg(i)*h/dx )*spaceLimiter(this,1,i,'SB')
			else
				do j=1,this%nx-1
					FxL(j) = this%vg(i)*this%f(j,i) - 0.5_mp*this%vg(i)*( 1.0_mp + this%vg(i)*h/dx )*spaceLimiter(this,j,i,'SB')
					FxR(j) = this%vg(i)*this%f(j+1,i) - 0.5_mp*this%vg(i)*( 1.0_mp + this%vg(i)*h/dx )*spaceLimiter(this,j+1,i,'SB')
				end do
				FxL(this%nx) = this%vg(i)*this%f(this%nx,i) - 0.5_mp*this%vg(i)*( 1.0_mp + this%vg(i)*h/dx )*spaceLimiter(this,pnx,i,'SB')
				FxR(this%nx) = this%vg(i)*this%f(1,i) - 0.5_mp*this%vg(i)*( 1.0_mp + this%vg(i)*h/dx )*spaceLimiter(this,1,i,'SB')
			end if
			this%f(:,i) = this%f(:,i) - h/dx*( FxR - FxL )
		end do
	end subroutine

	subroutine transportVelocity(this,h)
		type(plasma), intent(inout) :: this
		real(mp), intent(in) :: h					!time step
		real(mp) :: dx, dv
		integer :: i,j
		real(mp) :: FvL(2*this%nv+1), FvR(2*this%nv+1), acc
		dx = this%dx
		dv = this%dv

		acc = 0.0_mp
		do i=1,this%nx
			acc = this%qs*this%E(i)/this%ms

			FvL = 0.0_mp
			FvR = 0.0_mp
			if( acc>=0 )	then
				do j=2,2*this%nv+1
					FvL(j) = acc*this%f(i,j-1) + 0.5_mp*acc*( 1.0_mp - acc*h/dv )*velLimiter(this,i,j-1,'SB')
					FvR(j) = acc*this%f(i,j) + 0.5_mp*acc*( 1.0_mp - acc*h/dv )*velLimiter(this,i,j,'SB')
				end do
!				FvL(1) = 0.0_mp
				FvR(1) = acc*this%f(i,1) + 0.5_mp*acc*( 1.0_mp - acc*h/dv )*velLimiter(this,i,1,'SB')
			else
				do j=1,2*this%nv
					FvL(j) = acc*this%f(i,j) - 0.5_mp*acc*( 1.0_mp + acc*h/dv )*velLimiter(this,i,j,'SB')
					FvR(j) = acc*this%f(i,j+1) - 0.5_mp*acc*( 1.0_mp + acc*h/dv )*velLimiter(this,i,j+1,'SB')
				end do
				FvL(2*this%nv+1) = acc*this%f(i,2*this%nv+1) - 0.5_mp*acc*( 1.0_mp + acc*h/dv )*velLimiter(this,i,2*this%nv+1,'SB')
!				FvR(2*this%nv+1) = 0.0_mp
			end if
			this%f(i,:) = this%f(i,:) - h/dv*( FvR - FvL )
		end do
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
