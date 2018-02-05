module modSource

	use modPlasma
	use modCircuit

	implicit none

	abstract interface
		subroutine Source(p,c,dt)
			import plasma, circuit, mp
			type(plasma), intent(inout) :: p(:)
			type(circuit), intent(inout) :: c
			real(mp), intent(in) :: dt
		end subroutine
	end interface

contains

	subroutine NullSource(p,c,dt)
		type(plasma), intent(inout) :: p(:)
		type(circuit), intent(inout) :: c
		real(mp), intent(in) :: dt
	end subroutine

	subroutine ConstantIon(p,c,dt)
		type(plasma), intent(inout) :: p(:)
		type(circuit), intent(inout) :: c
		real(mp), intent(in) :: dt
		real(mp) :: Ls, vT_e, vT_i, N
		integer :: NGs, k
		real(mp), allocatable :: f_src(:,:)
		real(mp) :: W
		vT_e = p(1)%A(1)
		vT_i = p(2)%A(1)
		N = p(2)%A(3)
		p(1)%A(3) = 0.0_mp
		p(2)%A(3) = 0.0_mp
		Ls = 0.2_mp*c%Lx
		NGs = CEILING(Ls/c%dx)
		Ls = NGs*c%dx

		allocate(f_src(NGs,2*p(1)%nv+1))
		W = p(1)%dv/SQRT(2.0_mp*pi)/vT_e*SUM( EXP(-p(1)%vg**2/2.0_mp/vT_e/vT_e) )
		f_src(1,:) = N/Ls/W/SQRT(2.0_mp*pi)/vT_e*EXP( -p(1)%vg**2/2.0_mp/vT_e/vT_e )
		do k=2,NGs
			f_src(k,:) = f_src(1,:)
		end do
		p(1)%f(1:NGs,:) = p(1)%f(1:NGs,:) + f_src

		deallocate(f_src)
		allocate(f_src(NGs,2*p(2)%nv+1))
		W = p(2)%dv/SQRT(2.0_mp*pi)/vT_i*SUM( EXP(-p(2)%vg**2/2.0_mp/vT_i/vT_i) )
		f_src(1,:) = N/Ls/W/SQRT(2.0_mp*pi)/vT_i*EXP( -p(2)%vg**2/2.0_mp/vT_i/vT_i )
		do k=2,NGs
			f_src(k,:) = f_src(1,:)
		end do
		p(2)%f(1:NGs,:) = p(2)%f(1:NGs,:) + f_src
	end subroutine

end module
