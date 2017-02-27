module init

	use modPlasma
	use modRecord

	implicit none

contains

	subroutine initial_debye(p,vT,Q)
		type(plasma), intent(inout) :: p
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
		call setBackGround(p,rho_back)
	end subroutine

	subroutine initial_twostream(p,a,v0,vT)
		type(plasma), intent(inout) :: p
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
		call setBackGround(p,rho_back)
	end subroutine

end module