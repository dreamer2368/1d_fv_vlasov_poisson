module timeStep

	use modQoI
	use modRecord
	use Limiter

	implicit none

contains

	subroutine forward_sweep(p,c,r,inputQoI)
		use modQoI
		type(plasma), intent(inout) :: p(:)
		type(circuit), intent(inout) :: c
		type(history), intent(inout) :: r
		procedure(QoI), optional :: inputQoI
		procedure(QoI), pointer :: targetQoI=>DO_NOTHING
		integer :: i

		if( PRESENT(inputQoI) ) then
			targetQoI=>inputQoI
		end if

		do i=1,r%nt
			call updatePlasma(p,c,r)
			call targetQoI(p,c,i,r%j(i))
			call recordPlasma(r,p,c,i)
!			if( MOD(i,1000) == 0 ) then
!				print *, i
!			end if
		end do
	end subroutine

	subroutine updatePlasma(p,c,r)											!each time step
		type(plasma), intent(inout) :: p(:)
		type(circuit), intent(inout) :: c
		type(history), intent(inout) :: r
		integer :: i
		real(mp) :: dt
		real(mp) :: time1,time2
		dt=r%dt

		call CPU_TIME(time1)
!		call transportSpace(this,0.5_mp*dt,this%dx,this%dv)
		do i=1,SIZE(p)
			call transportSpace(p(i),dt)
		end do
		call CPU_TIME(time2)
		r%cpt_temp(1) = r%cpt_temp(1) + (time2-time1)/r%nmod

		c%rho = 0.0_mp
		do i=1,SIZE(p)
			call NumberDensity(p(i)%f,p(i)%dv,p(i)%n)
			c%rho = c%rho + p(i)%qs*p(i)%n
		end do
		call c%Efield
		call CPU_TIME(time1)
		r%cpt_temp(2) = r%cpt_temp(2) + (time1-time2)/r%nmod

		do i=1,SIZE(p)
			call transportVelocity(p(i),c%E,dt)
		end do
		call CPU_TIME(time2)
		r%cpt_temp(3) = r%cpt_temp(3) + (time2-time1)/r%nmod

!		call transportSpace(this,0.5_mp*dt,this%dx,this%dv)
	end subroutine

	subroutine transportSpace(this,h)
		type(plasma), intent(inout) :: this
		real(mp), intent(in) :: h					!time step
		real(mp) :: dx, dv
		integer :: i,j,nx
		real(mp) :: nu
		real(mp), dimension(this%nx) :: theta_p, theta_m, Df0, Df1, Df2
		real(mp), dimension(-1:this%nx+2,2*this%nv+1) :: tempf
		procedure(FluxLimiter), pointer :: PtrFluxLimiter=>MC
		dx = this%dx
		dv = this%dv
		nx = this%nx

		call this%PtrBC(this%f,tempf,this%vg,this%dx,h,this%A)

		!v<0
		do i=1,this%nv
			nu = this%vg(i)*h/dx
			Df0 = tempf(1:nx,i)-tempf(0:nx-1,i)
			Df1 = tempf(2:nx+1,i)-tempf(1:nx,i)
			Df2 = tempf(3:nx+2,i)-tempf(2:nx+1,i)
			theta_m = Df1/Df0
			theta_p = Df2/Df1
			this%f(:,i) = tempf(1:nx,i) - nu*Df1 + 0.5_mp*nu*(1.0_mp+nu)*( PtrFluxLimiter(theta_p)*Df1	&
																								- PtrFluxLimiter(theta_m)*Df0 )
		end do
		!v>=0
		do i=this%nv+1,2*this%nv+1
			nu = this%vg(i)*h/dx
			Df0 = tempf(2:nx+1,i)-tempf(1:nx,i)
			Df1 = tempf(1:nx,i)-tempf(0:nx-1,i)
			Df2 = tempf(0:nx-1,i)-tempf(-1:nx-2,i)
			theta_m = Df2/Df1
			theta_p = Df1/Df0
			this%f(:,i) = tempf(1:nx,i) - nu*Df1 - 0.5_mp*nu*(1.0_mp-nu)*( PtrFluxLimiter(theta_p)*Df0	&
																									- PtrFluxLimiter(theta_m)*Df1 )
		end do
	end subroutine

	subroutine transportVelocity(this,E,h)
		type(plasma), intent(inout) :: this
		real(mp), dimension(this%nx), intent(in) :: E
		real(mp), intent(in) :: h					!time step
		real(mp) :: dx, dv
		integer :: i,j
		real(mp) ::  acc, nu
		real(mp), dimension(this%nx,-this%nv-2:this%nv+2) :: tempf
		real(mp), dimension(2*this%nv+1) :: theta_p, theta_m, Df0, Df1, Df2
		integer :: NV
		procedure(FluxLimiter), pointer :: PtrFluxLimiter=>MC
		dx = this%dx
		dv = this%dv
		NV = this%nv

		tempf(:,-NV:NV) = this%f
		tempf(:,-NV-2:-NV-1) = 0.0_mp
		tempf(:,NV+1:NV+2) = 0.0_mp

		acc = 0.0_mp
		do i=1,this%nx
			acc = this%qs*E(i)/this%ms
			nu = acc*h/dv
			if( acc.ge.0.0_mp )	then
				Df0 = tempf(i,-NV+1:NV+1)-tempf(i,-NV:NV)
				Df1 = tempf(i,-NV:NV)-tempf(i,-NV-1:NV-1)
				Df2 = tempf(i,-NV-1:NV-1)-tempf(i,-NV-2:NV-2)
				theta_m = Df2/Df1
				theta_p = Df1/Df0
				this%f(i,:) = tempf(i,-NV:NV) - nu*Df1 - 0.5_mp*nu*(1.0_mp-nu)*( PtrFluxLimiter(theta_p)*Df0	&
																										- PtrFluxLimiter(theta_m)*Df1 )
			else
				Df0 = tempf(i,-NV:NV)-tempf(i,-NV-1:NV-1)
				Df1 = tempf(i,-NV+1:NV+1)-tempf(i,-NV:NV)
				Df2 = tempf(i,-NV+2:NV+2)-tempf(i,-NV+1:NV+1)
				theta_m = Df1/Df0
				theta_p = Df2/Df1
				this%f(i,:) = tempf(i,-NV:NV) - nu*Df1 + 0.5_mp*nu*(1.0_mp+nu)*( PtrFluxLimiter(theta_p)*Df1	&
																										- PtrFluxLimiter(theta_m)*Df0 )
			end if
		end do
	end subroutine



!==============  Continuous-Forward Sensitivity  ========================================

	subroutine forward_sensitivity(p,c,r,dp,dc,dr,inputQoI)
		type(plasma), intent(inout) :: p(:),dp(:)
		type(circuit), intent(inout) :: c,dc
		type(history), intent(inout) :: r,dr
		procedure(QoI) :: inputQoI
		integer :: i
		real(mp) :: time1,time2

		do i=1,r%nt
			call updateSensitivity(p,c,r,dp,dc,dr)
			call inputQoI(p,c,i,r%j(i))
			call inputQoI(dp,dc,i,dr%j(i))

         call CPU_TIME(time1)
			call recordPlasma(r,p,c,i)
         call CPU_TIME(time2)
         r%cpt_temp(6) = r%cpt_temp(6) + (time2-time1)/r%nmod

			call recordPlasma(dr,dp,dc,i)
         call CPU_TIME(time1)
         dr%cpt_temp(6) = dr%cpt_temp(6) + (time1-time2)/r%nmod
		end do
	end subroutine

	subroutine updateSensitivity(p,c,r,dp,dc,dr)
		type(plasma), intent(inout) :: p(:),dp(:)
		type(circuit), intent(inout) :: c,dc
		type(history), intent(inout) :: r,dr
		real(mp) :: time1,time2
		integer :: i

		!Original simulation
		call CPU_TIME(time1)
		do i=1,SIZE(p)
			call transportSpace(p(i),r%dt)
		end do
		call CPU_TIME(time2)
		r%cpt_temp(1) = r%cpt_temp(1) + (time2-time1)/r%nmod

		c%rho=0.0_mp
		do i=1,SIZE(p)
			call NumberDensity(p(i)%f,p(i)%dv,p(i)%n)
			c%rho = c%rho + p(i)%qs*p(i)%n
		end do
		call c%Efield
		call CPU_TIME(time1)
		r%cpt_temp(2) = r%cpt_temp(2) + (time1-time2)/r%nmod

		do i=1,SIZE(p)
			call transportVelocity(p(i),c%E,r%dt)
		end do
		call CPU_TIME(time2)
		r%cpt_temp(3) = r%cpt_temp(3) + (time2-time1)/r%nmod

		!Sensitivity simulation
		do i=1,SIZE(dp)
			call transportSpace(dp(i),0.5_mp*r%dt)
		end do
		call CPU_TIME(time1)
		dr%cpt_temp(1) = dr%cpt_temp(1) + (time1-time2)/dr%nmod

		dc%rho = 0.0_mp
		do i=1,SIZE(dp)
			call NumberDensity(dp(i)%f,dp(i)%dv,dp(i)%n)
    		dc%rho = dp(i)%qs*dp(i)%n + p(i)%n
		end do
		dc%rho = dc%rho + dc%rho_back
		call dc%Efield
		call CPU_TIME(time2)
		dr%cpt_temp(2) = dr%cpt_temp(2) + (time2-time1)/dr%nmod

		do i=1,SIZE(dp)
			call transportVelocity(dp(i),c%E,r%dt)
		end do
		call CPU_TIME(time1)
		dr%cpt_temp(3) = dr%cpt_temp(3) + (time1-time2)/dr%nmod

		do i=1,SIZE(dp)
			call transportSpace(dp(i),0.5_mp*r%dt)
		end do
		call CPU_TIME(time2)
		dr%cpt_temp(4) = dr%cpt_temp(4) + (time2-time1)/dr%nmod

		do i=1,SIZE(dp)
			call sourceSensitivity(dp(i),p(i),dc%E,c%E,r%dt)
		end do
		call CPU_TIME(time1)
		dr%cpt_temp(5) = dr%cpt_temp(5) + (time1-time2)/dr%nmod
	end subroutine

	subroutine sourceSensitivity(dp,p,dE,E,dt)
		type(plasma), intent(inout) :: dp
		type(plasma), intent(in) :: p
		real(mp), dimension(dp%nx), intent(in) :: dE, E
		real(mp), intent(in) :: dt
		real(mp) :: dx, dv
		integer :: i,j
		real(mp) ::  acc, nu, theta_p, theta_m
		real(mp) :: fp1,f0,fm1,f2
		real(mp), dimension(dp%nx,2*dp%nv+1) :: newf
		integer :: NV
		procedure(FluxLimiter_temp), pointer :: PtrFluxLimiter=>MC_temp
		dx = dp%dx
		dv = dp%dv
		NV = 2*dp%nv+1
		newf = 0.0_mp

		acc = 0.0_mp
		do i=1,dp%nx
			!dJdqp
			acc = (dp%qs*dE(i)+E(i))/dp%ms
			nu = acc*dt/dv
			if( acc.ge.0.0_mp )	then
				do j=1,NV
					fp1 = MERGE( p%f(i,j+1), 0.0_mp, j<NV )
					f0 = p%f(i,j)
					fm1 = MERGE( p%f(i,j-1), 0.0_mp, j>1 )
					f2 = MERGE( p%f(i,j-2), 0.0_mp, j>2 )
					theta_m = (fm1-f2)/(f0-fm1)
					theta_p = (f0-fm1)/(fp1-f0)
					newf(i,j) = dp%f(i,j) - nu*(f0-fm1) - 0.5_mp*nu*(1.0_mp-nu)*( PtrFluxLimiter(theta_p)*(fp1-f0)	&
																									- PtrFluxLimiter(theta_m)*(f0-fm1) )
				end do
			else
				do j=1,NV
					f2 = MERGE( p%f(i,j+2), 0.0_mp, j<NV-1 )
					fp1 = MERGE( p%f(i,j+1), 0.0_mp, j<NV )
					f0 = p%f(i,j)
					fm1 = MERGE( p%f(i,j-1), 0.0_mp, j>1 )
					theta_m = (fp1-f0)/(f0-fm1)
					theta_p = (f2-fp1)/(fp1-f0)
					newf(i,j) = dp%f(i,j) - nu*(fp1-f0) + 0.5_mp*nu*(1.0_mp+nu)*( PtrFluxLimiter(theta_p)*(fp1-f0)	&
																									- PtrFluxLimiter(theta_m)*(f0-fm1) )
				end do
			end if
		end do
		dp%f = newf
	end subroutine

end module
