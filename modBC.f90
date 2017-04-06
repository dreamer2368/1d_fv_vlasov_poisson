module modBC

	use constants

	implicit none

	abstract interface
		subroutine BC(f,f_gc,A)
			use constants
			real(mp), intent(in) :: f(:,:)
			real(mp), intent(out) :: f_gc(-1:size(f,1)+2,size(f,2))
			real(mp), intent(in), optional :: A
		end subroutine
	end interface

contains

	subroutine Periodic(f,f_gc,A)
		real(mp), intent(in) :: f(:,:)
		real(mp), intent(out) :: f_gc(-1:size(f,1)+2,size(f,2))
		real(mp), intent(in), optional :: A
		integer :: nx
		nx = size(f,1)
		
		f_gc(1:nx,:) = f
		f_gc(-1:0,:) = f(nx-1:nx,:)
		f_gc(nx+1:nx+2,:) = f(1:2,:)
	end subroutine

end module
