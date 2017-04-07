module modPlasmaBC

	use constants

	implicit none

	abstract interface
		subroutine PlasmaBC(f,f_gc,vg,dx,A)
			use constants
			real(mp), intent(in) :: f(:,:)
			real(mp), intent(out) :: f_gc(-1:size(f,1)+2,size(f,2))
			real(mp), intent(in) :: vg(size(f,2)), dx
			real(mp), intent(in), optional :: A
		end subroutine
	end interface

contains

	subroutine Periodic(f,f_gc,vg,dx,A)
		real(mp), intent(in) :: f(:,:)
		real(mp), intent(out) :: f_gc(-1:size(f,1)+2,size(f,2))
		real(mp), intent(in) :: vg(size(f,2)), dx
		real(mp), intent(in), optional :: A
		integer :: nx
		nx = size(f,1)
		
		f_gc(1:nx,:) = f
		f_gc(-1:0,:) = f(nx-1:nx,:)
		f_gc(nx+1:nx+2,:) = f(1:2,:)
	end subroutine

	subroutine testBC(f,f_gc,vg,dx,t)
		real(mp), intent(in) :: f(:,:)
		real(mp), intent(out) :: f_gc(-1:size(f,1)+2,size(f,2))
		real(mp), intent(in) :: vg(size(f,2)), dx
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

end module
