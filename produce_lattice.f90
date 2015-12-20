! TODO: how to produce reciprocal lattice

PROGRAM pro_lat
	use lat_par
	IMPLICIT NONE

	call init_lat
	call pro_site
	call pro_bond
	call output_lat
	!call output_lat_py
	!call pl_lat

	STOP

CONTAINS
	INCLUDE "chain.f90"              ! 1
	INCLUDE "square.f90"             ! 2
	INCLUDE "triangular.f90"         ! 3
	INCLUDE "honeycomb.f90"          ! 4
	INCLUDE "kagome.f90"             ! 5
	INCLUDE "cubic.f90"              ! 8
	INCLUDE "pyrochlore.f90"         ! 9
	SUBROUTINE init_lat
		IMPLICIT NONE

		integer       :: laty

		print *, "Please choose one type of lattices:"
		print *, " 1. chain"
		print *, " 2. square"
		print *, " 3. triangular"
		print *, " 4. honeycomb"
		print *, " 5. kagome"
		print *, " 6. "
		print *, " 7. "
		print *, " 8. cubic "
		print *, " 9. pyrochlore "
		print *, "your choice -> "
		read  *, laty
		print *, "Please enter the size of lattice, evev*even*1 or even*even*even"
		print *, "Lx, Ly, Lz"
		read  *, Lx, Ly, Lz

		!if( mod(Lx,2)/=0 .or. mod(Ly,2)/=0 .or. (mod(Lz,2)/=0 .and. Lz/=1) ) then
		!	print *, "Lx, Ly must be even. For 2D, Lz=1; for 3D, Lz must be even"
		!	stop
		!end if

		LatVec    = 0.d0
		SubLatVec = 0.d0
		Lxyz      = (/Lx,Ly,Lz/)

		select case (laty)
		case (1)
			call chain
			write(*,*) "chain"
		case (2)
			call square
			write(*,*) "square"
		case (3)
			call triangular
			write(*,*) "triangular"
		case (4)
			call honeycomb
			write(*,*) "honeycomb"
		case (5)
			call kagome 
			write(*,*) "kagome"
			stop
		case (6)
			!call 
			write(*,*) " "
			stop
		case (7)
			!call 
			write(*,*) " "
			stop
		case (8)
			call cubic
			write(*,*) "cubic"
		case (9)
			call pyrochlore
			write(*,*) "pyrochlore"
		case default
			print *, "your enter is wrong"
			stop
		end select

		RETURN
	END SUBROUTINE init_lat

	SUBROUTINE pro_site
		IMPLICIT NONE

		integer         :: ix, iy, iz
		integer         :: subi
		integer         :: num
		integer         :: bx=1, by=2, bz=4 ! 1:x  2:y  4:z
		integer         :: i

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
						Site(subi,num)%sn     = num*SubLat-SubLat+subi
						!---- boundary sites ----!
						if( ix==1 ) then    ! lie in x boundary
							if( abs(dot_product(Site(subi,num)%coord,RecLatVec(1,1:3)))<1.d-5 ) then
								Site(subi,num)%tb = Site(subi,num)%tb+bx
							end if
						end if
						if( iy==1 ) then    ! lie in y boundary
							if( abs(dot_product(Site(subi,num)%coord,RecLatVec(2,1:3)))<1.d-5 ) then
								Site(subi,num)%tb = Site(subi,num)%tb+by
							end if
						end if
						if( iz==1 ) then    ! lie in z boundary
							if( abs(dot_product(Site(subi,num)%coord,RecLatVec(3,1:3)))<1.d-5 ) then
								Site(subi,num)%tb = Site(subi,num)%tb+bz
							end if
						end if
						!------------------------!
					end do
				end do
			end do
		end do

		RETURN
	END SUBROUTINE pro_site

	SUBROUTINE pro_bond
		IMPLICIT NONE

		real                 :: distofsite    ! distance of neighbor sites
		integer, allocatable :: neig(:,:)     ! 1-Dimen: axis, Dimen+1: which sublattice
		type(la_si)          :: tempsite(4)
		integer              :: nbincell, nbbecell, nbsum

		integer              :: ix, iy, iz
		integer              :: nx, ny, nz
		integer              :: ncoord(3)
		integer              :: subi
		integer              :: num, num_n
		integer              :: nbcount       ! bond count
		real                 :: dist

		integer              :: i,j,k
		integer              :: si, sj

		if( SubLat==1 ) then
			distofsite = cal_distance(LatVec(1,1:3))
		else
			distofsite = cal_distance(SubLatVec(2,1:3))
		end if
		nbsum    = NumNeig*SubLat/2
		nbincell = (SubLat-1)*SubLat/2
		nbbecell = nbsum - nbincell
		allocate(neig(nbsum,5))

		!---- determine the neighbor in same cell ----!
		nbcount = 0
		do i=1, SubLat
			do j=i+1, SubLat
				nbcount = nbcount+1
				neig(nbcount,1:3) = (/0,0,0/)
				neig(nbcount,4:5) = (/i,j/)
			end do
		end do
		if( nbcount/=nbincell ) stop "bond in cell wrong"
		!---- determine the neighbor between different cell ----!
		num = NSite/SubLat
		C1: do iz=0,1
			C2: do iy=-1,1
				C3: do ix=-1,1
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
		!-------------------------------------------------------!

		!**** produce bonds in same and betweendifferent cell ****!
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
		!*********************************************************!

		RETURN
	END SUBROUTINE pro_bond

	FUNCTION cal_distance(vec)
		IMPLICIT NONE
		real          :: cal_distance
		real          :: vec(3)

		cal_distance = vec(1)**2.d0+vec(2)**2.d0+vec(3)**2.d0
		cal_distance = cal_distance**0.5d0

		RETURN
	END FUNCTION cal_distance

	SUBROUTINE cal_reclatvec
		IMPLICIT NONE

		integer        :: i
		real           :: rec_vol

		rec_vol = dot_product( LatVec(1,1:3), cross_pro(LatVec(2,1:3),LatVec(3,1:3)) )
		RecLatVec(1,1:3) = 2.d0*PI*cross_pro(LatVec(2,1:3),LatVec(3,1:3))/rec_vol
		RecLatVec(2,1:3) = 2.d0*PI*cross_pro(LatVec(3,1:3),LatVec(1,1:3))/rec_vol
		RecLatVec(3,1:3) = 2.d0*PI*cross_pro(LatVec(1,1:3),LatVec(2,1:3))/rec_vol
		!write(*,*) rec_vol
		!write(*,*) RecLatVec(1,1:3), dot_product(RecLatVec(1,1:3),LatVec(2,1:3)), dot_product(RecLatVec(1,1:3),LatVec(3,1:3)), dot_product(RecLatVec(1,1:3),LatVec(1,1:3))
		!write(*,*) RecLatVec(2,1:3), dot_product(RecLatVec(2,1:3),LatVec(3,1:3)), dot_product(RecLatVec(2,1:3),LatVec(1,1:3)), dot_product(RecLatVec(2,1:3),LatVec(2,1:3))
		!write(*,*) RecLatVec(3,1:3), dot_product(RecLatVec(3,1:3),LatVec(1,1:3)), dot_product(RecLatVec(3,1:3),LatVec(2,1:3)), dot_product(RecLatVec(3,1:3),LatVec(3,1:3))
		!stop

		RETURN
	END SUBROUTINE cal_reclatvec

	FUNCTION cross_pro( vec1, vec2 )
		IMPLICIT NONE

		real           :: cross_pro(3)
		real           :: vec1(3), vec2(3)

		cross_pro(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
		cross_pro(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
		cross_pro(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)

		RETURN
	END FUNCTION cross_pro

	SUBROUTINE output_lat
		IMPLICIT NONE

		integer              :: i, subi
		integer, allocatable :: coubond(:) 
		call system('rm -f lattice')
		call system('rm -f site.list')
		call system('rm -f bond.list')
		open(1, file='lattice',access='append')
		write(1,'(1X,6I4,A20)') Dimen, Lx, Ly, Lz, SubLat, NumNeig, trim(NameLat)
		close(1)
		open(1,file='site.list',access='append')
		open(2,file='bond.list',access='append')
		write(1,'(1X, I5)') NSite
		!write(2,'(1X,A5,2X,I5)') 'Bond=', NSite*NumNeig/2
		write(2,'(1X, I5)') NSite*NumNeig/2
		!**** output sites ****!
		do i=1, NSite/SubLat
			do subi=1, SubLat
				write(1,'(1X,I5,1X,I5,1X,I2,1X,I2,1X,3I4,1X,3F10.4)') &
					Site(subi,i)%sn, i, subi, Site(subi,i)%tb, Site(subi,i)%icoord(:), Site(subi,i)%coord(:)
			end do
		end do
		!**** output bonds ****!
		allocate(coubond(NSite))
		coubond=0
		do i=1, NSite*NumNeig/2
			write(2,'(1X,I5,2X,I5,1X,3I3)') Bond(i,1), Bond(i,2), Relt(i,1:3)
			coubond(Bond(i,1))=coubond(Bond(i,1))+1
			coubond(Bond(i,2))=coubond(Bond(i,2))+1
		end do
		close(1)
		close(2)

		!*** check neighbor number ****!
		do i=1, NSite
			if( coubond(NSite) /= NumNeig ) then
				write(*,*) 'neighbor wrong'
			end if
		end do

		RETURN
	END SUBROUTINE output_lat

	SUBROUTINE output_lat_py
		IMPLICIT NONE

		integer              :: i, subi
		integer              :: s1,s2,sb1,sb2,si1,si2
		integer, allocatable :: coubond(:) 
		integer              :: nbcount
		real                 :: distofsite

		if( SubLat==1 ) then
			distofsite = cal_distance(LatVec(1,1:3))
		else
			distofsite = cal_distance(SubLatVec(2,1:3))
		end if

		call system('rm -f coordinate.list')
		open(1,file='coordinate.list',access='append')

		write(1,*) "{'Interaction': [],"
		write(1,*) "'Lines': ["
		nbcount = NSite*NumNeig/2
		do i=1, nbcount
			s1  = Bond(i,1)
			s2  = Bond(i,2)
			si1 = (s1-1)/SubLat+1
			si2 = (s2-1)/SubLat+1
			!Site(subi,num)%sn = num*SubLat-SubLat+subi
			sb1 = mod(s1-1,SubLat)+1
			sb2 = mod(s2-1,SubLat)+1
			if( Site(sb1,si1)%sn/=s1 ) stop "s1 wrong"
			if( Site(sb2,si2)%sn/=s2 ) stop "s2 wrong"
			if( abs(cal_distance(Site(sb1,si1)%coord-Site(sb2,si2)%coord)-distofsite)<1.d-4 ) then
				write(1,'("[(",I4,",",I4,"),",I3,"],")') Bond(i,1)-1, Bond(i,2)-1,0
			end if
		end do
		close(1)


		open(1,file='coordinate.list',access='append')
		write(1,*) "'Points': ["
		!**** output sites ****!
		do i=1, NSite/SubLat
			do subi=1, SubLat
				write(1,'(1X,"[(",F6.3,",",F6.3,",",F6.3,"),", &
					"(",I3,",",I3,",",I3,"),", I2,"],")') &
					Site(subi,i)%coord(1),Site(subi,i)%coord(2),Site(subi,i)%coord(3), &
					Site(subi,i)%icoord(1)-1,Site(subi,i)%icoord(2)-1,Site(subi,i)%icoord(3)-1, &
					subi
			end do
		end do
		close(1)

		RETURN
	END SUBROUTINE output_lat_py

END PROGRAM pro_lat
