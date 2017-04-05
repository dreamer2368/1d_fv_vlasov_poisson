module Limiter

	use constants
!	use modPlasma

	implicit none

	abstract interface
		function FluxLimiter(theta) result(rho)
			use constants
			real(mp), intent(in) :: theta
			real(mp) :: rho
		end function
	end interface

contains

	function minmod(theta) result(rho)
		real(mp), intent(in) :: theta
		real(mp) :: rho
		if( theta<0.0_mp ) then
			rho = 0.0_mp
		else
			rho = MIN( 1.0_mp, theta )
		end if
	end function

	function MC(theta) result(rho)
		real(mp), intent(in) :: theta
		real(mp) :: rho
		real(mp) :: a
		a = MINVAL( (/(1.0_mp+theta)/2.0_mp,2.0_mp,2.0_mp*theta/) )
		rho = MAX(0.0_mp,a)
	end function

	function SB(theta) result(rho)
		real(mp), intent(in) :: theta
		real(mp) :: rho
		real(mp) :: a,b
		a = MIN(1.0_mp,2.0_mp*theta)
		b = MIN(2.0_mp,theta)
		rho = MAXVAL( (/0.0_mp,a,b/) )
	end function

!	function FluxLimiter(theta,kind) result(rho)
!		real(mp), intent(in) :: theta
!		character(len=*), intent(in) :: kind
!		real(mp) :: rho
!		real(mp) :: a,b
!
!		select case (kind)
!			case ('minmod')
!				if( theta<0.0_mp ) then
!					rho = 0.0_mp
!				else
!					rho = MIN(1.0_mp,theta)
!				end if
!			case ('MC')
!				a = MINVAL( (/(1.0_mp+theta)/2.0_mp,2.0_mp,2.0_mp*theta/) )
!				rho = MAX(0.0_mp,a)
!			case ('SB')	!superbee
!				a = MIN(1.0_mp,2.0_mp*theta)
!				b = MIN(2.0_mp,theta)
!				rho = MAXVAL( (/0.0_mp,a,b/) )
!		end select
!	end function

!
!	function spaceLimiter(p,i,j,kind) result(rho)
!		type(plasma), intent(in) :: p
!		integer, intent(in) :: i,j				!i:x grid point, j:v grid point
!		character(len=*), intent(in) :: kind
!		real(mp) :: rho
!
!		select case (kind)
!			case ('minmod')
!				rho = minModSpace(p,i,j)
!			case ('MC')
!				rho = MCSpace(p,i,j)
!			case ('SB')	!superbee
!				rho = SBSpace(p,i,j)
!		end select
!	end function
!
!	function velLimiter(p,i,j,kind) result(rho)
!		type(plasma), intent(in) :: p
!		integer, intent(in) :: i,j				!i:x grid point, j:v grid point
!		character(len=*), intent(in) :: kind
!		real(mp) :: rho
!
!		select case (kind)
!			case ('minmod')
!				rho = minModVel(p,i,j)
!			case ('MC')
!				rho = MCVel(p,i,j)
!			case ('SB') !superbee
!				rho = SBVel(p,i,j)
!		end select
!	end function
!
!!minMod Limiter
!
!	function minModSpace(p,i,j) result(rho)
!		type(plasma), intent(in) :: p
!		integer, intent(in) :: i,j				!i:x grid point, j:v grid point
!		real(mp) :: rho, a, b					!a:upwind, b:downwind
!
!		if( i==1 ) then
!			a = ( p%f(i+1,j) - p%f(i,j) )
!			b = ( p%f(1,j) - p%f(p%nx,j) )
!		elseif( i==p%nx ) then
!			a = ( p%f(1,j) - p%f(p%nx,j) )
!			b = ( p%f(i,j) - p%f(i-1,j) )
!		else
!			a = ( p%f(i+1,j) - p%f(i,j) )
!			b = ( p%f(i,j) - p%f(i-1,j) )
!		end if
!
!		if( a*b<0 ) then
!			rho = 0.0_mp
!		elseif( ABS(a)>=ABS(b) ) then
!			rho = b
!		else
!			rho = a
!		end if
!	end function
!
!	function minModVel(p,i,j) result(rho)
!		type(plasma), intent(in) :: p
!		integer, intent(in) :: i,j				!i:x grid point, j:v grid point
!		real(mp) :: rho, a, b
!
!		if( j==1 ) then
!			a = ( p%f(i,j+1) - p%f(i,j) )
!			b = p%f(i,1)
!		elseif( j==2*p%nv+1 ) then
!			a = -p%f(i,2*p%nv+1)
!			b = ( p%f(i,j) - p%f(i,j-1) )
!		else
!			a = ( p%f(i,j+1) - p%f(i,j) )
!			b = ( p%f(i,j) - p%f(i,j-1) )
!		end if
!
!		if( a*b<0 ) then
!			rho = 0.0_mp
!		elseif( ABS(a)>=ABS(b) ) then
!			rho = b
!		else
!			rho = a
!		end if
!	end function
!
!!MC Limiter
!
!	function MCSpace(p,i,j) result(rho)
!		type(plasma), intent(in) :: p
!		integer, intent(in) :: i,j				!i:x grid point, j:v grid point
!		real(mp) :: rho, a, b, c					!a:upwind*2, b:downwind*2, c:centered
!
!		if( i==1 ) then
!			a = ( p%f(i+1,j) - p%f(i,j) )*2.0_mp
!			b = ( p%f(1,j) - p%f(p%nx,j) )*2.0_mp
!			c = ( p%f(i+1,j) - p%f(p%nx,j) )/2.0_mp
!		elseif( i==p%nx ) then
!			a = ( p%f(1,j) - p%f(p%nx,j) )*2.0_mp
!			b = ( p%f(i,j) - p%f(i-1,j) )*2.0_mp
!			c = ( p%f(1,j) - p%f(i-1,j) )/2.0_mp
!		else
!			a = ( p%f(i+1,j) - p%f(i,j) )*2.0_mp
!			b = ( p%f(i,j) - p%f(i-1,j) )*2.0_mp
!			c = ( p%f(1+1,j) - p%f(i-1,j) )/2.0_mp
!		end if
!
!		if( a*b<0 ) then
!			rho = 0.0_mp
!		else
!			rho = a
!			if( ABS(b)<ABS(a) ) then
!				rho = b
!				if( ABS(c)<ABS(b) ) then
!					rho = c
!				end if
!			elseif( ABS(c)<ABS(a) ) then
!				rho = c
!			end if
!		end if
!	end function
!
!	function MCVel(p,i,j) result(rho)
!		type(plasma), intent(in) :: p
!		integer, intent(in) :: i,j				!i:x grid point, j:v grid point
!		real(mp) :: rho, a, b, c				!a:upwind*2, b:downwind*2, c:centered
!
!		if( j==1 ) then
!			a = ( p%f(i,j+1) - p%f(i,j) )*2.0_mp
!			b = p%f(i,1)*2.0_mp
!			c = p%f(i,j+1)/2.0_mp
!		elseif( j==2*p%nv+1 ) then
!			a = -p%f(i,2*p%nv+1)*2.0_mp
!			b = ( p%f(i,j) - p%f(i,j-1) )*2.0_mp
!			c = -p%f(i,j-1)/2.0_mp
!		else
!			a = ( p%f(i,j+1) - p%f(i,j) )*2.0_mp
!			b = ( p%f(i,j) - p%f(i,j-1) )*2.0_mp
!			c = ( p%f(i,j+1) - p%f(i,j-1) )/2.0_mp
!		end if
!
!		if( a*b<0 ) then
!			rho = 0.0_mp
!		else
!			rho = a
!			if( ABS(b)<ABS(a) ) then
!				rho = b
!				if( ABS(c)<ABS(b) ) then
!					rho = c
!				end if
!			elseif( ABS(c)<ABS(a) ) then
!				rho = c
!			end if
!		end if
!	end function
!
!!SuperBee Limiter
!
!	function SBSpace(p,i,j) result(rho)
!		type(plasma), intent(in) :: p
!		integer, intent(in) :: i,j				!i:x grid point, j:v grid point
!		real(mp) :: rho, rho1, rho2, a, b					!a:upwind, b:downwind
!
!		if( i==1 ) then
!			a = ( p%f(i+1,j) - p%f(i,j) )
!			b = ( p%f(1,j) - p%f(p%nx,j) )
!		elseif( i==p%nx ) then
!			a = ( p%f(1,j) - p%f(p%nx,j) )
!			b = ( p%f(i,j) - p%f(i-1,j) )
!		else
!			a = ( p%f(i+1,j) - p%f(i,j) )
!			b = ( p%f(i,j) - p%f(i-1,j) )
!		end if
!
!		if( a*b<0 ) then
!			rho1 = 0.0_mp
!		elseif( ABS(a)>=ABS(2.0_mp*b) ) then
!			rho1 = 2.0_mp*b
!		else
!			rho1 = a
!		end if
!
!		if( a*b<0 ) then
!			rho2 = 0.0_mp
!		elseif( ABS(2.0_mp*a)>=ABS(b) ) then
!			rho2 = b
!		else
!			rho2 = 2.0_mp*a
!		end if
!
!		if( rho1*rho2<0 ) then
!			rho = 0.0_mp
!		elseif( ABS(rho1)>=ABS(rho2) ) then
!			rho = rho1
!		else
!			rho = rho2
!		end if
!	end function
!
!	function SBVel(p,i,j) result(rho)
!		type(plasma), intent(in) :: p
!		integer, intent(in) :: i,j				!i:x grid point, j:v grid point
!		real(mp) :: rho, rho1, rho2, a, b					!a:upwind, b:downwind
!
!		if( j==1 ) then
!			a = ( p%f(i,j+1) - p%f(i,j) )
!			b = p%f(i,1)
!		elseif( j==2*p%nv+1 ) then
!			a = -p%f(i,2*p%nv+1)
!			b = ( p%f(i,j) - p%f(i,j-1) )
!		else
!			a = ( p%f(i,j+1) - p%f(i,j) )
!			b = ( p%f(i,j) - p%f(i,j-1) )
!		end if
!
!		if( a*b<0 ) then
!			rho1 = 0.0_mp
!		elseif( ABS(a)>=ABS(2.0_mp*b) ) then
!			rho1 = 2.0_mp*b
!		else
!			rho1 = a
!		end if
!
!		if( a*b<0 ) then
!			rho2 = 0.0_mp
!		elseif( ABS(2.0_mp*a)>=ABS(b) ) then
!			rho2 = b
!		else
!			rho2 = 2.0_mp*a
!		end if
!
!		if( rho1*rho2<0 ) then
!			rho = 0.0_mp
!		elseif( ABS(rho1)>=ABS(rho2) ) then
!			rho = rho1
!		else
!			rho = rho2
!		end if
!	end function

end module
