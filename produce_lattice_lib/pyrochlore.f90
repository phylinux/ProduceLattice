
!============ Pyrochlore lattice ============!
SUBROUTINE pyrochlore
	use lat_par
	IMPLICIT NONE

	NameLat = 'pyrochlore lattice'
	Dimen   = 3
	SubLat  = 4
	NumNeig = 6
	NumNeig2= 12
	NumNeig3= 6
	NSite   = Lx*Ly*Lz*SubLat

	!**** For produce sites ****!
	LatVec(1,1) = 1.d0
	LatVec(1,2) = 1.d0
	LatVec(1,3) = 0.d0
	LatVec(2,1) = 0.d0
	LatVec(2,2) = 1.d0
	LatVec(2,3) = 1.d0
	LatVec(3,1) = 1.d0
	LatVec(3,2) = 0.d0
	LatVec(3,3) = 1.d0
	SubLatVec(1,1) = 0.d0
	SubLatVec(1,2) = 0.d0
	SubLatVec(1,3) = 0.d0
	SubLatVec(2,1) = 0.5d0
	SubLatVec(2,2) = 0.5d0
	SubLatVec(2,3) = 0.d0
	SubLatVec(3,1) = 0.d0
	SubLatVec(3,2) = 0.5d0
	SubLatVec(3,3) = 0.5d0
	SubLatVec(4,1) = 0.5d0
	SubLatVec(4,2) = 0.d0
	SubLatVec(4,3) = 0.5d0
	allocate(Site(SubLat,NSite/SubLat))

	!**** For produce bonds ****!
	allocate(Bond(NSite*NumNeig/2,2))
	allocate(Relt(NSite*NumNeig/2,3))

	RETURN
END SUBROUTINE pyrochlore

SUBROUTINE pyrochlore_second_neighbor
	use lat_par
	IMPLICIT NONE

	real, external       :: cal_distance

	real                 :: distofsite    ! distance of neighbor sites
	integer              :: neig(100,5)     ! 1-Dimen: axis, Dimen+1: which sublattice
	type(la_si)          :: tempsite(4)

	integer              :: ix, iy, iz
	integer              :: nx, ny, nz
	integer              :: ncoord(3)
	integer              :: subi
	integer              :: num, num_n
	integer              :: nbcount       ! bond count
	integer              :: nbsum         ! bond count
	real                 :: dist

	integer              :: i,j,k
	integer              :: si, sj

	if( allocated(Bond) ) then
		deallocate(Bond)
	end if
	if( allocated(Relt) ) then
		deallocate(Relt)
	end if
	allocate(Bond(NSite*NumNeig2/2,2))
	allocate(Relt(NSite*NumNeig2/2,3))

	distofsite = cal_distance(LatVec(1,1:3))*sqrt(3.d0)/2.d0
	nbcount = 0

	!---- determine the neighbor between different cell ----!
	num = NSite/SubLat  ! arbitrary
	C1: do iz=0,Lz/2
		C2: do iy=-Ly/2,Ly/2
			C3: do ix=-Lx/2,Lx/2
				if( ix==0 .and. iy==0 .and. iz==0 ) cycle C3
				! select bonds between sites in which any couple are not in same line
				do i=1, nbcount
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
				do i=1, SubLat
					do j=1, SubLat
						dist = cal_distance( tempsite(j)%coord-Site(i,num)%coord )
						if( abs(dist-distofsite)<1.d-4 ) then
							nbcount  = nbcount+1
							neig(nbcount,1:3) = tempsite(1)%icoord(1:3)
							neig(nbcount,4:5) = (/i,j/)
						end if
					end do
				end do
			end do C3
		end do C2
	end do C1
	!-------------------------------------------------------!

	!open(1, file='temp', action='write')
	!do i=1, nbcount
	!	write(1,'(5I4)') neig(i,:)
	!end do
	!close(1)

	!**** produce bonds ***************************************!
	nbsum = 0
	do num=1, NSite/SubLat
		do i=1, nbcount
			ncoord = Site(1,num)%icoord(:)+neig(i,1:3)
			nx = mod(ncoord(1)-1+Lxyz(1), Lxyz(1))+1
			ny = mod(ncoord(2)-1+Lxyz(2), Lxyz(2))+1
			nz = mod(ncoord(3)-1+Lxyz(3), Lxyz(3))+1
			num_n = nx+(ny-1)*Lx+(nz-1)*Lx*Ly
			si = neig(i,4)
			sj = neig(i,5)
			nbsum = nbsum+1
			Bond(nbsum,1) = Site(si,num)%sn
			Bond(nbsum,2) = Site(sj,num_n)%sn
			Relt(nbsum,1:3) = neig(i,1:3)
		end do
	end do
	if( nbsum /= (NSite*NumNeig2/2) ) stop "Second neighbor wrong"
	!*********************************************************!

	RETURN
END SUBROUTINE pyrochlore_second_neighbor

SUBROUTINE pyrochlore_third_neighbor
	use lat_par
	IMPLICIT NONE

	interface
		function cross_pro(vec1, vec2)
			real           :: cross_pro(3)
			real           :: vec1(3), vec2(3)
		end function
	end interface

	real, external       :: cal_distance

	real                 :: distofsite      ! distance of neighbor sites
	integer              :: neig(100,5)     ! 1-Dimen: axis, Dimen+1: which sublattice
	type(la_si)          :: tempsite(4)

	integer              :: ix, iy, iz
	integer              :: nx, ny, nz
	integer              :: ncoord(3)
	integer              :: subi
	integer              :: num, num_n
	integer              :: nbcount       ! bond count
	integer              :: nbsum         ! bond count
	real                 :: dist
	real                 :: relvec(4,4,3)

	integer              :: i,j,k,l
	integer              :: si, sj
	real                 :: vec(3)
	real                 :: crvec(3)

	if( allocated(Bond) ) then
		deallocate(Bond)
	end if
	if( allocated(Relt) ) then
		deallocate(Relt)
	end if
	allocate(Bond(NSite*NumNeig3/2,2))
	allocate(Relt(NSite*NumNeig3/2,3))

	do i=1, 4
		do j=1, 4
			relvec(i,j,1:3) = Site(i,1)%coord - Site(j,1)%coord
		end do
	end do

	!!!!!! NOT DONE !!!!!!
	distofsite = cal_distance(LatVec(1,1:3))
	nbcount = 0

	!---- determine the neighbor between different cell ----!
	num = NSite/SubLat  ! arbitrary
	C1: do iz=0,Lz/2
		C2: do iy=-Ly/2,Ly/2
			C3: do ix=-Lx/2,Lx/2
				if( ix==0 .and. iy==0 .and. iz==0 ) cycle C3
				! select bonds between sites in which any couple are not in same line
				do i=1, nbcount
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
				do i=1, SubLat
					do j=1, SubLat
						dist = cal_distance( tempsite(j)%coord-Site(i,num)%coord )
						vec  = Site(i,num)%coord - tempsite(j)%coord
						k    = 1
						crvec = 1.d0
						do while ( k < 4 )
							do l=1, 4
								if ( l /= i ) then
									crvec(k) = cal_distance(cross_pro(vec,relvec(i,l,1:3)))
									k = k+1
								end if
							end do
						end do
						if( abs(dist-distofsite)<1.d-4 ) then
							if ( crvec(1)<1e-4 .or. crvec(2)<1e-4 .or. crvec(3)<1e-4 ) then
								nbcount  = nbcount+1
								neig(nbcount,1:3) = tempsite(1)%icoord(1:3)
								neig(nbcount,4:5) = (/i,j/)
							end if
						end if
					end do
				end do
			end do C3
		end do C2
	end do C1
	!-------------------------------------------------------!

	!open(1, file='temp', action='write')
	!do i=1, nbcount
	!	write(1,'(5I4)') neig(i,:)
	!end do
	!close(1)

	!**** produce bonds ***************************************!
	nbsum = 0
	do num=1, NSite/SubLat
		do i=1, nbcount
			ncoord = Site(1,num)%icoord(:)+neig(i,1:3)
			nx = mod(ncoord(1)-1+Lxyz(1), Lxyz(1))+1
			ny = mod(ncoord(2)-1+Lxyz(2), Lxyz(2))+1
			nz = mod(ncoord(3)-1+Lxyz(3), Lxyz(3))+1
			num_n = nx+(ny-1)*Lx+(nz-1)*Lx*Ly
			si = neig(i,4)
			sj = neig(i,5)
			nbsum = nbsum+1
			Bond(nbsum,1) = Site(si,num)%sn
			Bond(nbsum,2) = Site(sj,num_n)%sn
			Relt(nbsum,1:3) = neig(i,1:3)
		end do
	end do
	if( nbsum /= (NSite*NumNeig3/2) ) stop "Third neighbor wrong"
	!*********************************************************!

	RETURN
END SUBROUTINE pyrochlore_third_neighbor
