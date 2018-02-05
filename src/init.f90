module init

	use modPlasma
	use modCircuitBC	
	use modRecord

	implicit none

contains

	subroutine initial_sheath(p,c,n0,vT_e,vT_i)
		type(plasma), intent(inout) :: p(2)
		type(circuit), intent(inout) :: c
		real(mp), intent(in) :: n0, vT_e, vT_i
		integer :: i

		do i=1,c%nx
			p(1)%f(i,:) = n0/SQRT(2.0_mp*pi)/vT_e*EXP( -p(1)%vg**2/2.0_mp/vT_e/vT_e )
			p(2)%f(i,:) = n0/SQRT(2.0_mp*pi)/vT_i*EXP( -p(2)%vg**2/2.0_mp/vT_i/vT_i )
		end do
		c%rho_back = 0.0_mp

		p(1)%PtrBC=>Refluxing_Absorbing
		p(2)%PtrBC=>Refluxing_Absorbing
		c%Efield=>DirichletNeumann
		c%updateCircuit=>Electrode
	end subroutine

	subroutine initial_debye_sensitivity(dp,dc,vT,input_str)
		type(plasma), intent(inout) :: dp
		type(circuit), intent(inout) :: dc
		real(mp), intent(in) :: vT
      character(len=*), intent(in) :: input_str
		integer :: i,j
		real(mp) :: fv(2*dp%nv+1), f(dp%nx,2*dp%nv+1)
		real(mp) :: w, rho_back(dp%nx)
		w = dp%Lx/10.0_mp

      SELECT CASE(input_str)
         CASE('vT')
		      fv = (dp%vg**2/vT/vT-1.0_mp)/sqrt(2.0_mp*pi)/vT/vT*exp( -dp%vg**2/2.0_mp/vT/vT )
		      do i=1,size(dp%xg)
		      	f(i,:) = fv
	      	end do
            rho_back = 0.0_mp
         CASE('qp')
            f = 0.0_mp
            rho_back = -1.0_mp
      END SELECT

		call setPlasma(dp,f)
		call setBackGround(dc,rho_back)
	end subroutine

	subroutine initial_debye(p,c,vT,Q)
		type(plasma), intent(inout) :: p
		type(circuit), intent(inout) :: c
		real(mp), intent(in) :: vT,Q
		integer :: i,j
		real(mp) :: fv(2*p%nv+1), f(p%nx,2*p%nv+1)
		real(mp) :: w, rho_back(p%nx)
		w = p%Lx/10.0_mp

		fv = 1.0_mp/sqrt(2.0_mp*pi)/vT*exp( -p%vg**2/2.0_mp/vT/vT )
		do i=1,size(p%xg)
			f(i,:) = fv
		end do
		f = f/( SUM(integrate_dv(f,p%dv))*p%dx )*p%Lx

		rho_back = 1.0_mp - Q/p%Lx + Q/SQRT(2.0_mp*pi)/w*EXP( -(p%xg-0.5_mp*p%Lx)**2/2.0_mp/w/w )

		call setPlasma(p,f)
		call setBackGround(c,rho_back)
	end subroutine

	subroutine initial_twostream(p,c,a,v0,vT)
		type(plasma), intent(inout) :: p
		type(circuit), intent(inout) :: c
		real(mp), intent(in) :: a, v0, vT				!perturbation amplitude
		integer :: i,j
		real(mp) :: fx(p%nx), fv(2*p%nv+1), f(p%nx,2*p%nv+1)
		real(mp) :: rho_back(p%nx)

		fx = ( 1.0_mp + a*COS( 2.0_mp*pi*p%xg/p%Lx ) )
		fv = 0.5_mp/sqrt(2.0_mp*pi)/vT*( exp( -(p%vg-v0)**2/2.0_mp/vT/vT )	&
											+ exp( -(p%vg+v0)**2/2.0_mp/vT/vT ) )
		do i=1,size(p%xg)
			do j=1, size(p%vg)
				f(i,j) = fx(i)*fv(j)
			end do
		end do
		f = f/( SUM(integrate_dv(f,p%dv))*p%dx )*p%Lx

		rho_back = 1.0_mp

		call setPlasma(p,f)
		call setBackGround(c,rho_back)
	end subroutine

end module
