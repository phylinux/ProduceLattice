SUBROUTINE produce_bond(laty)
	implicit none
	integer(4), intent(in)      :: laty

	call produce_bond_general
END SUBROUTINE produce_bond

SUBROUTINE produce_bond_general
	IMPLICIT NONE

	real(8)                 :: distofsite    ! distance of neighbor sites
	integer(4), allocatable :: neig(:,:,:)   ! 1-Dimen: axis, Dimen+1: which sublattice
	type(la_si)             :: tempsite(SubLat)
	integer(4)              :: nbincell, nbbecell, nbsum

	integer(4)              :: ix, iy, iz
	integer(4)              :: nx, ny, nz
	integer(4)              :: ncoord(3)
	integer(4)              :: subi
	integer(4)              :: num, num_n
	integer(4)              :: nbcount       ! bond count
	real(8)                 :: dist

	integer(4)              :: nr
	integer(4)              :: i,j,k
	integer(4)              :: si, sj

	allocate(neig(NNr,maxval(NumNeig)*,5))
	nr = 1
	!if( SubLat==1 ) then
	!	distofsite = cal_distance(LatVec(1,1:3))
	!else
	!	distofsite = cal_distance(SubLatVec(2,1:3))
	!end if
	nbsum    = NumNeig(nr)*SubLat/2
	nbincell = (SubLat-1)*SubLat/2
	nbbecell = nbsum - nbincell

	!---- determine the neighbor in same cell ----!
	nbcount = 0
	do i=1, SubLat
		do j=i+1, SubLat
			dist = cal_distance( Site(j,1)%coord-Site(i,1)%coord )
			if( abs(dist-NNdt(nr))<1.d-4 ) then
				nbcount  = nbcount+1
				neig(nr,nbcount,1:3) = (/0,0,0/)
				neig(nbcount,4:5) = (/i,j/)
			end if
		end do
	end do
	!---- determine the neighbor between different cell ----!
	num = NSite/SubLat  !-- This number is arbitrary --!
	C1: do iz= 0,NNr
		C2: do iy=-NNr,NNr
			C3: do ix=-NNr,NNr
				if( ix==0 .and. iy==0 .and. iz==0 ) cycle C3
				! select nbbecell bonds in which any couple are not in same line
				do i=nbincell+1, nbcount
					if( neig(i,1)==(-ix) .and. neig(i,2)==(-iy) .and. neig(i,3)==(-iz) ) cycle C3
				end do
				if( Dimen==1 ) then
					tempsite(1)%icoord = (/ix,0,0/)
				else if( Dimen==2 ) then
					tempsite(1)%icoord = (/ix,iy,0/)
				else
					tempsite(1)%icoord = (/ix,iy,iz/)
				end if
				do i=1, SubLat
					tempsite(i)%coord = Site(i,num)%coord+ix*LatVec(1,:)+iy*LatVec(2,:)+iz*LatVec(3,:)
				end do
				!-- bonds between two cells --!
				do i=1, SubLat
					do j=1, SubLat
						dist = cal_distance( tempsite(j)%coord-Site(i,num)%coord )
						if( abs(dist-distofsite)<1.d-4 ) then
							nbcount  = nbcount+1
							neig(nbcount,1:3) = tempsite(1)%icoord(1:3)
							neig(nbcount,4:5) = (/i,j/)
							if( nbcount==nbsum ) exit C1
						end if
					end do
				end do
			end do C3
		end do C2
	end do C1
	!-----------------------------------------------------!

	!-- produce bonds in same and betweendifferent cell --!
	nbcount = 0
	do num=1, NSite/SubLat
		do i=1, nbsum
			ncoord = Site(1,num)%icoord(:)+neig(i,1:3)
			nx = mod(ncoord(1)-1+Lxyz(1), Lxyz(1))+1
			ny = mod(ncoord(2)-1+Lxyz(2), Lxyz(2))+1
			nz = mod(ncoord(3)-1+Lxyz(3), Lxyz(3))+1
			num_n = nx+(ny-1)*Lx+(nz-1)*Lx*Ly
			si = neig(i,4)
			sj = neig(i,5)
			nbcount = nbcount+1
			Bond(nbcount,1) = Site(si,num)%sn
			Bond(nbcount,2) = Site(sj,num_n)%sn
			Relt(nbcount,1:3) = neig(i,1:3)
		end do
	end do
	if( nbcount /= (NSite*NumNeig/2) ) stop "produce bond wrong"
	!-------------------------------------------------------!

END SUBROUTINE produce_bond_general
