! TODO: how to produce reciprocal lattice

INCLUDE "lat_par.f90"
PROGRAM pro_lat
	use lat_par
	IMPLICIT NONE

	call init_lat
	call allo_ary
	call produce_site
	call produce_bond
	call produce_neighbors_of_site
	call output_lat
	!call output_lat_py
	!call pl_lat

	STOP

CONTAINS
	INCLUDE "chain.f90"              !  1
	INCLUDE "square.f90"             !  2
	INCLUDE "triangular.f90"         !  3
	INCLUDE "honeycomb.f90"          !  4
	INCLUDE "kagome.f90"             !  5
	INCLUDE "checkboard.f90"         !  6
	INCLUDE "cubic.f90"              !  7
	INCLUDE "octahedra.f90"          !  8
	INCLUDE "pyrochlore.f90"         !  9
	INCLUDE "diamond.f90"            ! 10
	INCLUDE "packingAB.f90"          ! 11
	INCLUDE "packingABC.f90"         ! 12
	INCLUDE "produce_bond.f90"
	INCLUDE "calculation.f90"
	INCLUDE "output_lat.f90"
	SUBROUTINE init_lat
		IMPLICIT NONE

		!print *, "Please choose one type of lattices:"
		!print *, " 1. chain"
		!print *, " 2. square"
		!print *, " 3. triangular"
		!print *, " 4. honeycomb"
		!print *, " 5. kagome"
		!print *, " 6. checkboard"
		!print *, " 7. cubic"
		!print *, " 8. octahedra "
		!print *, " 9. pyrochlore "
		!print *, "10. diamond "
		!print *, "your choice -> "
		!read  *, laty
		!print *, "Please enter the size of lattice, Lx Ly Lz"
		!print *, "Lx, Ly, Lz"
		!read  *, Lx, Ly, Lz
		open(1,file="lat_inp",action="read")
		read(1,*) laty
		read(1,*) Lx, Ly, Lz
		close(1)

		LatVec    = 0.d0
		SubLatVec = 0.d0
		Lxyz      = (/Lx,Ly,Lz/)
		NumNeig = 0

		select case (laty)
		case ( 1); call chain
		case ( 2); call square
		case ( 3); call triangular
		case ( 4); call honeycomb
		case ( 5); call kagome
		case ( 6); call checkboard; stop
		case ( 7); call cubic
		case ( 8); call octahedra
		case ( 9); call pyrochlore
		case (10); call diamond
		case (11); call packingAB
		case (12); call packingABC
		case default
			print *, "your enter is wrong"
			stop
		end select
		write(*,'(A16,A16)') "Your chose is: ", trim(NameLat)

	END SUBROUTINE init_lat

	SUBROUTINE allo_ary
		implicit none
		allocate(Site(SubLat,NSite/SubLat))

		allocate(Bond(NNr,NSite*maxval(NumNeig)/2,2))
		allocate(NNst(NNr,NSite,maxval(NumNeig)))
		allocate(NNbk(NNr,NSite,maxval(NumNeig)))
		allocate(Relt(NNr,NSite*maxval(NumNeig)/2,NNr))
	END SUBROUTINE allo_ary

	SUBROUTINE produce_site
		IMPLICIT NONE

		integer(4)              :: ix, iy, iz
		integer(4)              :: subi
		integer(4)              :: num
		integer(4)              :: bx=1, by=2, bz=4 ! 1:x  2:y  4:z
		integer(4)              :: i

		if( Dimen==1 ) then
			bx=1; by=0; bz=0
		else if( Dimen==2 ) then
			bx=1; by=2; bz=0
		else if( Dimen==3 ) then
			bx=1; by=2; bz=4
		end if

		call cal_reclatvec

		do iz=1, Lz
		do iy=1, Ly
		do ix=1, Lx
			num = ix+(iy-1)*Lx+(iz-1)*Lx*Ly
			do subi=1, SubLat
				Site(subi,num)%coord(:)       = 0.d0
				Site(subi,num)%coord(1:Dimen) = &
					(ix-1)*LatVec(1,1:Dimen)  + &
					(iy-1)*LatVec(2,1:Dimen)  + &
					(iz-1)*LatVec(3,1:Dimen)  + &
					SubLatVec(subi,1:Dimen)
				Site(subi,num)%icoord(1)      = ix
				Site(subi,num)%icoord(2)      = iy
				Site(subi,num)%icoord(3)      = iz
				Site(subi,num)%sn     = (num-1)*SubLat +subi
				!---- boundary sites ----!
				!-- 1D: the boundary is a site
				!-- 2D: the boundary is a line
				!-- 3D: the boundary is a surface
				!-- usage: measure the winding number in undirect worm --!
				if( ix==1 ) then    ! lie in the boundary of axis x
					if( abs(dot_product(Site(subi,num)%coord,RecLatVec(1,1:3)))<1.d-5 ) then
						Site(subi,num)%tb = Site(subi,num)%tb+bx
					end if
				end if
				if( iy==1 ) then    ! lie in the boundary of axis y
					if( abs(dot_product(Site(subi,num)%coord,RecLatVec(2,1:3)))<1.d-5 ) then
						Site(subi,num)%tb = Site(subi,num)%tb+by
					end if
				end if
				if( iz==1 ) then    ! lie in the boundary of axis z
					if( abs(dot_product(Site(subi,num)%coord,RecLatVec(3,1:3)))<1.d-5 ) then
						Site(subi,num)%tb = Site(subi,num)%tb+bz
					end if
				end if
				!------------------------!
			end do
		end do
		end do
		end do

	END SUBROUTINE produce_site

	SUBROUTINE produce_neighbors_of_site
		implicit none
		integer(4), allocatable :: ns(:,:)
		integer(4)              :: n1, n2
		integer(4)              :: i, j

		! 1st: site number; 2: the n_th neighbors
		allocate(ns(NNr,NSite))
		ns = 0

		!do j=1, 1
		do j=1, NNr
		do i=1, NSite*NumNeig(j)/2
			n1 = Bond(j,i,1)
			n2 = Bond(j,i,2)
			ns(j,n1) = ns(j,n1) + 1
			ns(j,n2) = ns(j,n2) + 1
			if( ns(j,n1)>NumNeig(j) .or. ns(j,n2)>NumNeig(j) ) stop "Err: ns"
			NNst(j,n1,ns(j,n1)) = n2
			NNst(j,n2,ns(j,n2)) = n1
			NNbk(j,n1,ns(j,n1)) = ns(j,n2)
			NNbk(j,n2,ns(j,n2)) = ns(j,n1)
		end do
		end do

		deallocate(ns)

	END SUBROUTINE produce_neighbors_of_site


END PROGRAM pro_lat
