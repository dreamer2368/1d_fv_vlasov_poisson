module testmodules

	use timeStep
	use init
	use modRecord

	implicit none

contains

	! You can add custom subroutines/functions here later, if you want
	subroutine manufactured_solution
		type(plasma) :: p
		type(history) :: r
		real(mp), parameter :: L=270_mp, Lv=5.2_mp
		integer, parameter :: Nx=512, Nv=64
		real(mp), parameter :: vT = 1.0_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = 1.0_mp, me = 1.0_mp
		real(mp), parameter :: T=0.5_mp, CFL = 0.5_mp
		real(mp), dimension(Nx,2*Nv+1) :: f0, src
		real(mp) :: temp, error
		integer :: i,j

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me,eps0)
		do i=1,Nx
			do j=1,2*Nv+1
				f0(i,j) = SIN(2.0_mp*pi*p%xg(i)/p%Lx)/SQRT(2.0_mp*pi)/vT*EXP( -p%vg(j)**2/2.0_mp/vT/vT )
			end do
		end do
		call setPlasma(p,f0)
		p%rho_back = 0.0_mp
		call buildRecord(r,p,T,CFL=CFL,input_dir='test',nmod=20)
		do i=1,Nx
			do j=1,2*Nv+1
				src(i,j) = 0.5_mp*r%dt*( 2.0_mp*pi/L + L/2.0_mp/pi/vT/vT*SIN(2.0_mp*pi*p%xg(i)/L) )	&
								*COS(2.0_mp*pi*p%xg(i)/L)*p%vg(j)/SQRT(2.0_mp*pi)/vT*EXP( -p%vg(j)**2/2.0_mp/vT/vT )
			end do
		end do
		open(unit=302,file='data/'//r%dir//'/f0.bin',status='replace',form='unformatted',access='stream')
		write(302) p%f
		close(302)
		open(unit=302,file='data/'//r%dir//'/src.bin',status='replace',form='unformatted',access='stream')
		write(302) src
		close(302)

		error = 0.0_mp
		do i=1,r%nt
			p%f = p%f + src
			call transportSpace(p,0.5_mp*r%dt)
			call Efield(p)
			call transportVelocity(p,p%E,r%dt)
			call transportSpace(p,0.5_mp*r%dt)
			p%f = p%f + src
			call recordPlasma(r,p,i)

			temp = SQRT(SUM( (p%f-f0)**2 )*p%dx*p%dv)
			if( error<temp ) error = temp
		end do

		call printPlasma(r)

		print *, 'Nx, Nv: ',Nx, Nv
		print *, 'Error: ', error

		call destroyRecord(r)
		call destroyPlasma(p)
	end subroutine

	subroutine manufactured_solution_x
		type(plasma) :: p
		type(history) :: r
		real(mp), parameter :: L=2.7_mp, Lv=5.2_mp
		integer, parameter :: Nx=2048, Nv=32
		real(mp), parameter :: w = 0.1_mp*L
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = 1.0_mp, me = 1.0_mp
		real(mp), parameter :: Tf=0.5_mp, CFL = 0.5_mp
		real(mp) :: t, xc, xw(Nx)
		integer :: cgrid(Nx)
		real(mp), dimension(Nx,2*Nv+1) :: f0
		integer :: i,j,k
		real(mp) :: error=0.0_mp, temp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me,eps0)
		do j=1,2*Nv+1
			xc = L/2.0_mp
			xc = xc - FLOOR(xc/L)*L
			p%f(:,j) = 1.0_mp/SQRT(2.0_mp*pi)/w*EXP( -(p%xg-xc)**2/2.0_mp/w/w )
		end do
		p%rho_back = 0.0_mp
		call buildRecord(r,p,Tf,CFL=CFL,input_dir='testx',nmod=100)
		open(unit=302,file='data/'//r%dir//'/f0.bin',status='replace',form='unformatted',access='stream')
		write(302) p%f
		close(302)

		t=0.0_mp
		do k=1,r%nt
			do j=1,2*Nv+1
				xc = L/2.0_mp+p%vg(j)*(t+r%dt)
				xc = xc - FLOOR(xc/L)*L
				xw = p%xg-xc
				cgrid = FLOOR(xw/0.5_mp/L)
				xw = (xw-cgrid*0.5_mp*L)*(-1.0_mp)**cgrid + 0.25_mp*L*( 1.0_mp+(-1.0_mp)**(cgrid+1) )
				f0(:,j) = 1.0_mp/SQRT(2.0_mp*pi)/w*EXP( -xw**2/2.0_mp/w/w )
			end do
			call transportSpace(p,r%dt)
			call recordPlasma(r,p,k)
!			if( MOD(i,1000) == 0 ) then
!				print *, i
!			end if
			t=t+r%dt

			temp = SQRT(SUM( (p%f-f0)**2 )*p%dx*p%dv)
			if( error<temp ) error = temp
		end do
		call printPlasma(r)

		print *, 'Nx,Nv: ', Nx, Nv
		print *, 'error: ',error

		call destroyRecord(r)
		call destroyPlasma(p)
	end subroutine

	subroutine manufactured_solution_v
		type(plasma) :: p
		type(history) :: r
		real(mp), parameter :: L=2.7_mp, Lv=5.2_mp
		integer, parameter :: Nx=32, Nv=1024
		real(mp), parameter :: w = 0.1_mp*Lv
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = 1.0_mp, me = 1.0_mp
		real(mp), parameter :: Tf=1.0_mp, CFL = 0.5_mp
		real(mp) :: t, vc, vw(2*Nv+1)
		real(mp), dimension(Nx,2*Nv+1) :: f0
		integer :: i,j,k
		real(mp) :: error=0.0_mp, temp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me,eps0)
		do j=1,Nx
			p%f(j,:) = 1.0_mp/SQRT(2.0_mp*pi)/w*EXP( -p%vg**2/2.0_mp/w/w )
		end do
		p%E = 2.7_mp*SIN(2.0_mp*pi*p%xg/L)
		p%rho_back = 0.0_mp
		call buildRecord(r,p,Tf,CFL=CFL,input_dir='testv',nmod=20)
		open(unit=302,file='data/'//r%dir//'/f0.bin',status='replace',form='unformatted',access='stream')
		write(302) p%f
		close(302)

		t=0.0_mp
		do k=1,r%nt
			do j=1,Nx
				vc = p%E(j)*(t+r%dt)
				vw = p%vg-vc
				f0(j,:) = 1.0_mp/SQRT(2.0_mp*pi)/w*EXP( -vw**2/2.0_mp/w/w )
			end do
			call transportVelocity(p,p%E,r%dt)
			call recordPlasma(r,p,k)
!			if( MOD(i,1000) == 0 ) then
!				print *, i
!			end if
			t=t+r%dt

			temp = SQRT(SUM( (p%f-f0)**2 )*p%dx*p%dv)
			if( error<temp ) error = temp
		end do
		call printPlasma(r)

		print *, 'Nx,Nv: ', Nx, Nv
		print *, 'error: ',error

		call destroyRecord(r)
		call destroyPlasma(p)
	end subroutine

end module
