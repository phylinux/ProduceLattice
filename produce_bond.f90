SUBROUTINE produce_bond
	implicit none

	select case (laty)
	case ( 1); call produce_bond_general              ! chain
	case ( 2); call produce_bond_general              ! square
	case ( 3); call produce_bond_general              ! triangular
	case ( 4); call produce_bond_general              ! honeycomb
	case ( 5); call produce_bond_general              ! kagome
	case ( 6); call produce_bond_general              ! checkboard
	case ( 7); call produce_bond_general              ! cubic
	case ( 8); call produce_bond_general              ! octahedra
	case ( 9); call produce_bond_general              ! pyrochlore
	case (10); call produce_bond_general              ! diamond
	case (11); call produce_bond_general              ! packingAB
	end select
END SUBROUTINE produce_bond

SUBROUTINE produce_bond_general
	IMPLICIT NONE

	! 1st: the nth neighbors
	! 2nd: the mth of nth neighbor
	! 3rd: 1-3: integer coordinate, 4-5: sublattice
	integer(4), allocatable :: neig(:,:,:)
	type(la_si)             :: tempsite(SubLat)
	integer(4)              :: nbincell, nbbecell, nbsum

	integer(4)              :: ix, iy, iz
	integer(4)              :: nx, ny, nz
	integer(4)              :: ncoord(3)
	integer(4)              :: subi
	integer(4)              :: num, num_n             ! 
	integer(4)              :: nbcount                ! bond number
	real(8)                 :: dist

	integer(4)              :: nr
	integer(4)              :: i,j,k
	integer(4)              :: si, sj

	allocate(neig(NNr,maxval(NumNeig)*SubLat/2,5))
	neig = -1e6
	!C0:do nr = 1, 1
	C0:do nr = 1, NNr
		if( NumNeig(nr)==0 ) cycle C0
	!if( SubLat==1 ) then
	!	distofsite = cal_distance(LatVec(1,1:3))
	!else
	!	distofsite = cal_distance(SubLatVec(2,1:3))
	!end if
	nbsum    = NumNeig(nr)*SubLat/2 ! the number of bonds linked to a cell not a site

	!---- determine the neighbor in the same cell ----!
	nbcount  = 0
	nbincell = 0
	do i=1, SubLat
		do j=i+1, SubLat
			dist = cal_distance( Site(j,1)%coord-Site(i,1)%coord )
			if( abs(dist-NNdt(nr))<1.d-4 ) then
				nbincell = nbincell+1
				nbcount  = nbcount+1
				neig(nr,nbcount,1:3) = (/0,0,0/)
				neig(nr,nbcount,4:5) = (/i,j/)
			end if
		end do
	end do
	!write(*,*) nbincell
	!nbbecell = nbsum - nbincell
	!---- determine neighbors between different cells ----!
	num = NSite/SubLat  !-- This number is arbitrary --!
	C1: do iz= 0,nr
		C2: do iy=-nr,nr
			C3: do ix=-nr,nr
				if( ix==0 .and. iy==0 .and. iz==0 ) cycle C3
				! select nbbecell bonds in which any couple are not in the same line
				select case (Dimen)
				case(1); if( iy/=0 .or. iz/=0 ) then; cycle C3; end if
				case(2); if(            iz/=0 ) then; cycle C3; end if
				case(3); if(          .false. ) then; cycle C3; end if
				case default
					write(*,*) "Err: Dimen /= 1,2,3"; stop
				end select
				do i=nbincell+1, nbcount
					select case (Dimen)
					case(1); if( neig(nr,i,1)==(-ix)                                                     ) then; cycle C3; end if
					case(2); if( neig(nr,i,1)==(-ix) .and. neig(nr,i,2)==(-iy)                           ) then; cycle C3; end if
					case(3); if( neig(nr,i,1)==(-ix) .and. neig(nr,i,2)==(-iy) .and. neig(nr,i,3)==(-iz) ) then; cycle C3; end if
					end select
				end do
				if( Dimen==1 ) then
					tempsite(1)%icoord = (/ix,0,0/)
				else if( Dimen==2 ) then
					tempsite(1)%icoord = (/ix,iy,0/)
				else
					tempsite(1)%icoord = (/ix,iy,iz/)
				end if
				do i=1, SubLat
					tempsite(i)%coord = Site(i,num)%coord &
						& + tempsite(1)%icoord(1)*LatVec(1,:) &
						& + tempsite(1)%icoord(2)*LatVec(2,:) &
						& + tempsite(1)%icoord(3)*LatVec(3,:)
				end do
				!-- bonds between two cells --!
				do i=1, SubLat
					do j=1, SubLat
						dist = cal_distance( tempsite(j)%coord-Site(i,num)%coord )
						if( abs(dist-NNdt(nr))<1.d-4 ) then
							!write(*,*) tempsite(1)%icoord(1:3)
							nbcount  = nbcount+1
							neig(nr,nbcount,1:3) = tempsite(1)%icoord(1:3)
							neig(nr,nbcount,4:5) = (/i,j/)
							!if( nbcount==nbsum ) exit C1
							!write(*,*) nr, nbcount, neig(nr,nbcount,1:3)
						end if
					end do
				end do
			end do C3
		end do C2
	end do C1
	if( nbcount/=nbsum ) then
		write(*,*) nr, nbcount, nbsum
		write(*,*) "Err: nbcount /= nbsum"
		stop
	end if
	!-----------------------------------------------------!

	!-- produce bonds in same and betweendifferent cell --!
	nbcount = 0
	do num=1, NSite/SubLat
		do i=1, nbsum
			ncoord = Site(1,num)%icoord(:)+neig(nr,i,1:3)
			nx = mod(ncoord(1)-1+Lxyz(1), Lxyz(1))+1
			ny = mod(ncoord(2)-1+Lxyz(2), Lxyz(2))+1
			nz = mod(ncoord(3)-1+Lxyz(3), Lxyz(3))+1
			num_n = nx+(ny-1)*Lx+(nz-1)*Lx*Ly
			si = neig(nr,i,4)
			sj = neig(nr,i,5)
			nbcount = nbcount+1
			Bond(nr,nbcount,1) = Site(si,num)%sn
			Bond(nr,nbcount,2) = Site(sj,num_n)%sn
			Relt(nr,nbcount,1:3) = neig(nr,i,1:3)
		end do
	end do
	if( nbcount /= (NSite*NumNeig(nr)/2) ) stop "produce bond wrong"
	!-------------------------------------------------------!
	end do C0

END SUBROUTINE produce_bond_general

FUNCTION cal_distance(vec)
	IMPLICIT NONE
	real(8)          :: cal_distance
	real(8)          :: vec(3)

	cal_distance = vec(1)**2.d0+vec(2)**2.d0+vec(3)**2.d0
	cal_distance = cal_distance**0.5d0

END FUNCTION cal_distance
