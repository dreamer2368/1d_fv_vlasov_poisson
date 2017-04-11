module MatrixVector

	use constants

	implicit none

contains

	function multiplyD(x,dx) result(y)						!Derivative with periodic BC
		real(mp), intent(in) :: x(:)
		real(mp), intent(in) :: dx
		real(mp) :: y(size(x))
		integer :: i

		y=0.0_mp
		do i=2,size(x)-1
			y(i) = 0.5_mp/dx*( x(i+1) - x(i-1) )
		end do
		y(1) = 0.5_mp/dx*( x(2) - x(size(x)) )
		y(size(x)) = 0.5_mp/dx*( x(1) - x(size(x)-1) )
	end function

	function multiplyK(x,dx) result(y)
		real(mp), intent(in) :: x(:)
		real(mp), intent(in) :: dx
		real(mp) :: y(size(x))
		integer :: i

		y=0.0_mp
		do i=2,size(x)-1
			y(i) = 1.0_mp/dx/dx*( x(i+1) - 2.0_mp*x(i) + x(i-1) )
		end do
		y(1) = 1.0_mp/dx/dx*( x(2) - 2.0_mp*x(1) )
		y(size(x)) = 1.0_mp/dx/dx*( - 2.0_mp*x(size(x)) + x(size(x)-1) )		!periodic BC + Dirichlet BC
	end function

	subroutine CG_K(x,b,dx)								!Kx = b
		real(mp), intent(in) :: b(:)
		real(mp), intent(in) :: dx
		real(mp), intent(out) :: x(:)
		real(mp) :: r(size(x)), p(size(x)), r1(size(x))
		real(mp) :: alpha, beta
		real(mp) :: tol
		integer :: iter = 0

		select case (mp)
			case (SELECTED_REAL_KIND(4))
				tol = 1.0e-16
			case (SELECTED_REAL_KIND(15))
				tol = 1.0e-32
			case (SELECTED_REAL_KIND(33))
				tol = 10.0_mp**(-64.0_mp)
			case default
				tol = 1.0e-32
		end select

		if( size(b)/=size(x) ) then
			print *, '===================================================='
			print *, '====================  FAULT  ======================='
			print *, '=========  x AND b HAVE NOT EQUAL SIZES  ==========='
			print *, '====================  Ax=b  ========================'
			print *, '===================================================='
		end if

		x = 0.0_mp
		r = b - multiplyK(x,dx)
		p = r

		iter = 0
		do while( DOT_PRODUCT(r,r)>tol )
			alpha = DOT_PRODUCT(r,r)/DOT_PRODUCT(p,multiplyK(p,dx))
			x = x + alpha*p
			r1 = r
			r = r - alpha*multiplyK(p,dx)

			beta = DOT_PRODUCT(r,r)/DOT_PRODUCT(r1,r1)
			p = r + beta*p
			iter = iter+1
			if( iter > 1e8 ) then
				print *, '====================================='
				print *, '=========  CG METHOD FAILS   ========'
				print *, '====================================='
				exit
			end if
		end do
	end subroutine

	subroutine solve_tridiag(a,b,c,d,x,n)
		!	 a - sub-diagonal (means it is the diagonal below the main diagonal)
		!	 b - the main diagonal
		!	 c - sup-diagonal (means it is the diagonal above the main diagonal)
		!	 d - right part
		!	 x - the answer
		!	 n - number of equations
		integer,intent(in) :: n
		real(mp),dimension(n),intent(in) :: a,b,c,d
		real(mp),dimension(n),intent(out) :: x
		real(mp),dimension(n) :: cp,dp
		real(mp) :: m
		integer i

		! initialize c-prime and d-prime
		cp(1) = c(1)/b(1)
		dp(1) = d(1)/b(1)
		! solve for vectors c-prime and d-prime
		do i = 2,n
			m = b(i)-cp(i-1)*a(i)
			cp(i) = c(i)/m
			dp(i) = (d(i)-dp(i-1)*a(i))/m
		enddo
		! initialize x
		x(n) = dp(n)
		! solve for x from the vectors c-prime and d-prime
		do i = n-1, 1, -1
			x(i) = dp(i)-cp(i)*x(i+1)
		end do
	end subroutine

	function integrate_dv(f,dv) result(rho)
		!recommended to interpolate into G-L gridpoint for accurate integration
		real(mp), intent(in) :: f(:,:)
		real(mp), intent(in) :: dv
		real(mp) :: rho(size(f,1))
		rho = dv*SUM( f(:, 2:(size(f,2)-1) ), 2 )
		rho = rho + 0.5_mp*dv*( f(:,1) + f(:,size(f,2)) )
	end function

end module
