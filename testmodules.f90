module testmodules

	use timeStep
	use init
	use modRecord

	implicit none

contains

	! You can add custom subroutines/functions here later, if you want
	subroutine BoundaryTest
		type(plasma) :: p
		type(circuit) :: c
		type(history) :: r
		real(mp), parameter :: L=3.7_mp, Lv=1.3_mp
		integer, parameter :: Nx=128, Nv=128
		real(mp), parameter :: vT = 0.3_mp*Lv, w = Lv/Nv*2.0_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = 1.0_mp, me = 1.0_mp
		real(mp), parameter :: Tf=50.0_mp, CFL = 0.5_mp
		real(mp) :: t
		real(mp), dimension(Nx,2*Nv+1) :: f0
		integer :: i,j,k
		real(mp) :: error=0.0_mp, temp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me)
		call buildCircuit(c,L,Nx,eps0)
		p%A = vT
		do i=1,Nx
			do j=1,Nv
				f0(i,j) = 0.5_mp/SQRT(2.0_mp*pi)/w*EXP( -(p%vg(j)+2.0_mp*vT/SQRT(2.0_mp*pi))**2/2.0_mp/w/w )
			end do
			do j=Nv+1,2*Nv+1
				f0(i,j) = 1.0_mp/SQRT(2.0_mp*pi)/vT*EXP( -p%vg(j)**2/2.0_mp/vT/vT )
			end do
		end do
		p%f = f0
		c%rho_back = 0.0_mp
		p%PtrBC=>testRefluxing
		call buildRecord(r,p,c,Tf,CFL=CFL,input_dir='testBC',nmod=50)

		t=0.0_mp
		do k=1,r%nt
			call transportSpace(p,r%dt)
			call recordPlasma(r,p,c,k)
!			if( MOD(i,1000) == 0 ) then
!				print *, i
!			end if
			t=t+r%dt

!			temp = SQRT(SUM( (p%f-f0)**2 )*p%dx*p%dv)
!			if( error<temp ) error = temp
		end do
		call printPlasma(r)

		print *, 'Nx,Nv: ', Nx, Nv
!		print *, 'error: ',error

		call destroyRecord(r)
		call destroyPlasma(p)
		call destroyCircuit(c)
	end subroutine

	subroutine manufactured_solution
		type(plasma) :: p
		type(circuit) :: c
		type(history) :: r
		real(mp), parameter :: L=270_mp, Lv=5.2_mp
		integer, parameter :: Nx=1024, Nv=64
		real(mp), parameter :: vT = 1.0_mp
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = 1.0_mp, me = 1.0_mp
		real(mp), parameter :: T=0.5_mp, CFL = 0.5_mp
		real(mp), dimension(Nx,2*Nv+1) :: f0, src
		real(mp) :: temp, error
		integer :: i,j

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me)
		call buildCircuit(c,L,Nx,eps0)
		do i=1,Nx
			do j=1,2*Nv+1
				f0(i,j) = SIN(2.0_mp*pi*p%xg(i)/p%Lx)/SQRT(2.0_mp*pi)/vT*EXP( -p%vg(j)**2/2.0_mp/vT/vT )
			end do
		end do
		call setPlasma(p,f0)
		c%rho_back = 0.0_mp
		call buildRecord(r,p,c,T,CFL=CFL,input_dir='test',nmod=20)
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
			call NumberDensity(p%f,p%dv,p%n)
			c%rho = 0.0_mp
			c%rho = c%rho + p%qs*p%n
			call Efield(c)
			call transportVelocity(p,c%E,r%dt)
			call transportSpace(p,0.5_mp*r%dt)
			p%f = p%f + src
			call recordPlasma(r,p,c,i)

			temp = SQRT(SUM( (p%f-f0)**2 )*p%dx*p%dv)
			if( error<temp ) error = temp
		end do

		call printPlasma(r)

		print *, 'Nx, Nv: ',Nx, Nv
		print *, 'Error: ', error

		call destroyRecord(r)
		call destroyPlasma(p)
		call destroyCircuit(c)
	end subroutine

	subroutine manufactured_solution_x
		type(plasma) :: p
		type(circuit) :: c
		type(history) :: r
		real(mp), parameter :: L=2.7_mp, Lv=5.2_mp
		integer, parameter :: Nx=1024, Nv=32
		real(mp), parameter :: w = 0.1_mp*L
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = 1.0_mp, me = 1.0_mp
		real(mp), parameter :: Tf=0.5_mp, CFL = 0.5_mp
		real(mp) :: t, xc, xw(Nx)
		integer :: cgrid(Nx)
		real(mp), dimension(Nx,2*Nv+1) :: f0
		integer :: i,j,k
		real(mp) :: error=0.0_mp, temp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me)
		call buildCircuit(c,L,Nx,eps0)
		do j=1,2*Nv+1
			xc = L/2.0_mp
			xc = xc - FLOOR(xc/L)*L
			p%f(:,j) = 1.0_mp/SQRT(2.0_mp*pi)/w*EXP( -(p%xg-xc)**2/2.0_mp/w/w )
		end do
		c%rho_back = 0.0_mp
		call buildRecord(r,p,c,Tf,CFL=CFL,input_dir='testx',nmod=100)
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
			call recordPlasma(r,p,c,k)
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
		call destroyCircuit(c)
	end subroutine

	subroutine manufactured_solution_v
		type(plasma) :: p
		type(circuit) :: c
		type(history) :: r
		real(mp), parameter :: L=2.7_mp, Lv=5.2_mp
		integer, parameter :: Nx=32, Nv=512
		real(mp), parameter :: w = 0.1_mp*Lv
		real(mp), parameter :: eps0 = 1.0_mp, wp = 1.0_mp
		real(mp), parameter :: qe = 1.0_mp, me = 1.0_mp
		real(mp), parameter :: Tf=1.0_mp, CFL = 0.5_mp
		real(mp) :: t, vc, vw(2*Nv+1)
		real(mp), dimension(Nx,2*Nv+1) :: f0
		integer :: i,j,k,kr
		character(len=100) :: kstr
		real(mp) :: error=0.0_mp, temp

		call buildPlasma(p,L,Lv,Nx,Nv,qe,me)
		call buildCircuit(c,L,Nx,eps0)
		do j=1,Nx
			p%f(j,:) = 1.0_mp/SQRT(2.0_mp*pi)/w*EXP( -p%vg**2/2.0_mp/w/w )
		end do
		c%E = 0.7_mp*SIN(2.0_mp*pi*p%xg/L)

		c%rho_back = 0.0_mp
		call buildRecord(r,p,c,Tf,CFL=CFL,input_dir='testv',nmod=20)
		c%E = 0.7_mp*SIN(2.0_mp*pi*p%xg/L)
		open(unit=302,file='data/'//r%dir//'/f0.bin',status='replace',form='unformatted',access='stream')
		write(302) p%f
		close(302)

		t=0.0_mp
		do k=1,r%nt
			do j=1,Nx
				vc = c%E(j)*(t+r%dt)
				vw = p%vg-vc
				f0(j,:) = 1.0_mp/SQRT(2.0_mp*pi)/w*EXP( -vw**2/2.0_mp/w/w )
			end do
			call transportVelocity(p,c%E,r%dt)
			call recordPlasma(r,p,c,k)
			if( (r%nmod.eq.1) .or. (mod(k,r%nmod).eq.0) ) then
				kr = merge(k,k/r%nmod,r%nmod.eq.1)
				write(kstr,*) kr
				open(unit=402,file='data/'//r%dir//'/f0_'//trim(adjustl(kstr))//'.bin',status='replace',form='unformatted',access='stream')
				write(402) f0
				close(402)
			end if
			t=t+r%dt

			temp = SQRT(SUM( (p%f-f0)**2 )*p%dx*p%dv)
			if( error<temp ) error = temp
		end do
		call printPlasma(r)

		print *, 'Nx,Nv: ', Nx, Nv
		print *, 'error: ',error

		call destroyRecord(r)
		call destroyPlasma(p)
		call destroyCircuit(c)
	end subroutine

end module
