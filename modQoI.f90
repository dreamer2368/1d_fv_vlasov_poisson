module modQoI

	use modPlasma
	use MatrixVector

	implicit none

contains

	subroutine Screening_distance(p,k,j)
		type(plasma), intent(in) :: p
		integer, intent(in) :: k
		real(mp), intent(inout) :: j
		real(mp), dimension(p%nx) :: n
		real(mp) :: xc
		xc=0.5_mp*p%Lx

		n = integrate_dv(p%f,p%dv)
		j = SUM( (p%xg-xc)**2*n )*p%dx
print *, xc, MAXVAL(n), j
	end subroutine

end module
