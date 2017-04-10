module modPlasmaBC

	use Limiter

	implicit none

	abstract interface
		subroutine PlasmaBC(f,f_gc,vg,dx,dt,A)
			use constants
			real(mp), intent(in) :: f(:,:)
			real(mp), intent(out) :: f_gc(-1:size(f,1)+2,size(f,2))
			real(mp), intent(in) :: vg(size(f,2)), dx, dt
			real(mp), intent(in), optional :: A
		end subroutine
	end interface

contains

	subroutine Periodic(f,f_gc,vg,dx,dt,A)
		real(mp), intent(in) :: f(:,:)
		real(mp), intent(out) :: f_gc(-1:size(f,1)+2,size(f,2))
		real(mp), intent(in) :: vg(size(f,2)), dx, dt
		real(mp), intent(in), optional :: A
		integer :: nx
		nx = size(f,1)
		
		f_gc(1:nx,:) = f
		f_gc(-1:0,:) = f(nx-1:nx,:)
		f_gc(nx+1:nx+2,:) = f(1:2,:)
	end subroutine

	subroutine testBC(f,f_gc,vg,dx,dt,t)
		real(mp), intent(in) :: f(:,:)
		real(mp), intent(out) :: f_gc(-1:size(f,1)+2,size(f,2))
		real(mp), intent(in) :: vg(size(f,2)), dx, dt
		real(mp), intent(in), optional :: t
		real(mp), dimension((size(f,2)-1)/2) :: vp, Tf, w
		integer :: nx,nv
		nx = size(f,1)
		nv = (size(f,2)-1)/2
		
		f_gc(1:nx,:) = f

		vp = -vg(1:nv)
		Tf = 4.0_mp/vp
		w = 3.0_mp/Tf*2.0_mp*pi
		f_gc(nx+1,1:nv) = vp/dx*0.5_mp*( dx/vp - 1.0_mp/w*( COS(w*(t+dx/vp)-pi/2.0_mp) - COS(w*t-pi/2.0_mp) ) )
		f_gc(nx+2,1:nv) = vp/dx*0.5_mp*( dx/vp - 1.0_mp/w*( COS(w*(t+2.0_mp*dx/vp)-pi/2.0_mp)	&
																								- COS(w*(t+dx/vp)-pi/2.0_mp) ) )
		f_gc(0,1:nv+1) = f(1,nv+1:2*nv+1)

		vp = vg(nv+2:2*nv+1)
		Tf = 4.0_mp/vp
		w = 3.0_mp/Tf*2.0_mp*pi
		f_gc(0,nv+2:2*nv+1) = vp/dx*0.5_mp*( dx/vp - 1.0_mp/w*( COS(w*(t+dx/vp)-pi/2.0_mp) - COS(w*t-pi/2.0_mp) ) )
		f_gc(-1,nv+2:2*nv+1) = -vp/dx*0.5_mp*( dx/vp - 1.0_mp/w*( COS(w*(t+2.0_mp*dx/vp)-pi/2.0_mp)	&
																								- COS(w*(t+dx/vp)-pi/2.0_mp) ) )
		f_gc(nx+1,nv+1:2*nv+1) = f(nx,nv+1:2*nv+1)
	end subroutine

	subroutine testRefluxing(f,f_gc,vg,dx,dt,vT)
		real(mp), intent(in) :: f(:,:)
		real(mp), intent(out) :: f_gc(-1:size(f,1)+2,size(f,2))
		real(mp), intent(in) :: vg(size(f,2)), dx, dt
		real(mp), intent(in), optional :: vT
		real(mp) :: OutFlux
		procedure(FluxLimiter), pointer :: PtrFluxLimiter=>MC
		integer :: nx,nv
		real(mp), dimension((size(f,2)-1)/2) :: nu, theta, vp
		real(mp) :: dv, vc, vc_app, w
		nx = size(f,1)
		nv = (size(f,2)-1)/2
		dv = vg(2)-vg(1)

		f_gc(1:nx,:) = f
		!x=0, Outflux v<=0
		f_gc(0,1:nv+1) = f(1,1:nv+1)
		nu = vg(1:nv)*dt/dx
		theta = ( f_gc(2,1:nv)-f_gc(1,1:nv) )/( f_gc(1,1:nv)-f_gc(0,1:nv) )
		OutFlux = dx/dt*dv*SUM( -nu*f_gc(1,1:nv) - 0.5_mp*nu*(1.0_mp+nu)*PtrFluxLimiter(theta)*( f_gc(1,1:nv)-f_gc(0,1:nv) ) )
print *, 'x=0, outflux: ',OutFlux

		!x=0, Influx v>0
		vp = vg(nv+2:2*nv+1)
		vc = dv*SUM( vp*EXP( -vp**2/2.0_mp/vT/vT ) )
		f_gc(0,nv+2:2*nv+1) = OutFlux/vc*EXP( -vp**2/2.0_mp/vT/vT )
		f_gc(-1,nv+2:2*nv+1) = f_gc(0,nv+2:2*nv+1)
		nu = vp*dt/dx
		theta = ( f_gc(1,nv+2:2*nv+1)-f_gc(0,nv+2:2*nv+1) )/( f_gc(0,nv+2:2*nv+1)-f_gc(-1,nv+2:2*nv+1) )
		OutFlux = dx/dt*dv*SUM( nu*f_gc(0,nv+2:2*nv+1) + 0.5_mp*nu*(1.0_mp-nu)*PtrFluxLimiter(theta)	&
															*( f_gc(0,nv+2:2*nv+1)-f_gc(-1,nv+2:2*nv+1) ) )
print *, 'x=0, influx: ',OutFlux

		!x=L, Outflux v>=0
		f_gc(nx+1,nv+1:2*nv+1) = f(nx,nv+1:2*nv+1)
		nu = vg(nv+2:2*nv+1)*dt/dx
		theta = ( f_gc(nx,nv+2:2*nv+1)-f_gc(nx-1,nv+2:2*nv+1) )/( f_gc(nx+1,nv+2:2*nv+1)-f_gc(nx,nv+2:2*nv+1) )
		OutFlux = dx/dt*dv*SUM( nu*f_gc(nx,nv+2:2*nv+1) + 0.5_mp*nu*(1.0_mp-nu)*PtrFluxLimiter(theta)	&
																		*( f_gc(nx+1,nv+2:2*nv+1)-f_gc(nx,nv+2:2*nv+1) ) )
print *, 'x=L, outflux: ',OutFlux

		!x=L, Influx v<0
		vp = vg(1:nv)
		w = 2.0_mp*dv
		vc = -2.0_mp*vT/SQRT(2.0_mp*pi)
		vc_app = -dv*SUM( vp*EXP(-(vp-vc)**2/2.0_mp/w/w) )
		f_gc(nx+1,1:nv) = OutFlux/vc_app*EXP( -(vp-vc)**2/2.0_mp/w/w )
		f_gc(nx+2,1:nv) = f_gc(nx+1,1:nv)
		nu = vp*dt/dx
		theta = ( f_gc(nx+2,1:nv)-f_gc(nx+1,1:nv) )/( f_gc(nx+1,1:nv)-f_gc(nx,1:nv) )
		OutFlux = dx/dt*dv*SUM( -nu*f_gc(nx+1,1:nv) - 0.5_mp*nu*(1.0_mp+nu)*PtrFluxLimiter(theta)*( f_gc(nx+1,1:nv)-f_gc(nx,1:nv) ) )
print *, 'x=L, influx: ',OutFlux
	end subroutine

end module
