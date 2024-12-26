! TODO: how to produce reciprocal lattice

INCLUDE "lat_par.f90"
PROGRAM produce_lattice
	use lat_par
	implicit none

	call read_parameter
	call allo_ary
	call cal_reclatvec
	call produce_site
	call search_neighbors
	call calculate_bond
	call calculate_neighbor
	call write_lattice

	STOP
CONTAINS

SUBROUTINE read_parameter
	implicit none
	integer(4)                  :: i

	!open(1,file="lat_inp", action="read")
	!read(1,*) NameLat
	!read(1,*) Dimen
	!read(1,*) Lx, Ly, Lz
	!read(1,*) SubLat
	!read(1,*) LatVec(1,1:3)
	!read(1,*) LatVec(2,1:3)
	!read(1,*) LatVec(3,1:3)
	!allocate(SubLatVec(SubLat,3))
	!do i=1, SubLat
	!	read(1,*) SubLatVec(i,1:3)
	!end do
	!close(1)

	open(1,file="lat_inp",action="read")
	read(1,*) laty
	read(1,*) Lx, Ly, Lz
	read(1,*) MaxNr
	close(1)

	Lxyz  = (/Lx,Ly,Lz/)

	select case (laty)
	case ( 1); call chain
	case ( 2); call square
	case ( 3); call triangular
	case ( 4); call honeycomb
	case ( 5); call kagome
	case ( 6); call checkboard; stop
	case ( 7); call ruby
	case ( 8); call cubic
	case ( 9); call octahedra; stop
	case (10); call pyrochlore
	case (11); call diamond; stop
	case (12); call packingAB; stop
	case (13); call packingABC; stop
	case default
		print *, "your enter is wrong"
		stop
	end select
	write(*,'(A16,A16)') "Your chose is: ", trim(NameLat)

END SUBROUTINE read_parameter

SUBROUTINE allo_ary
	implicit none
	allocate(Site(NSite))
END SUBROUTINE allo_ary

SUBROUTINE produce_site
	IMPLICIT NONE

	integer(4)              :: ix, iy, iz
	integer(4)              :: sb
	integer(4)              :: icell
	integer(4)              :: isite
	integer(4)              :: bx, by, bz  ! boundary types
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
		icell = ix+(iy-1)*Lx+(iz-1)*Lx*Ly
		do sb=1, SubLat
			isite = (icell-1)*SubLat + sb
			Site(isite)%cell = icell
			Site(isite)%coord(:)       = 0.d0
			Site(isite)%coord(1:Dimen) = &
				(ix-1)*LatVec(1,1:Dimen)  + &
				(iy-1)*LatVec(2,1:Dimen)  + &
				(iz-1)*LatVec(3,1:Dimen)  + &
				SubLatVec(sb,1:Dimen)
			!---- boundary sites ----!
			!-- 1D: the boundary is a site
			!-- 2D: the boundary is a line
			!-- 3D: the boundary is a surface
			!-- usage: measure the winding number in undirect worm --!
			if( ix==1 ) then    ! lie in the boundary of axis x
				if( abs(dot_product(Site(isite)%coord,RecLatVec(1,1:3)))<1.d-5 ) then
					Site(isite)%tb = Site(isite)%tb+bx
				end if
			end if
			if( iy==1 ) then    ! lie in the boundary of axis y
				if( abs(dot_product(Site(isite)%coord,RecLatVec(2,1:3)))<1.d-5 ) then
					Site(isite)%tb = Site(isite)%tb+by
				end if
			end if
			if( iz==1 ) then    ! lie in the boundary of axis z
				if( abs(dot_product(Site(isite)%coord,RecLatVec(3,1:3)))<1.d-5 ) then
					Site(isite)%tb = Site(isite)%tb+bz
				end if
			end if
			!------------------------!
		end do
	end do
	end do
	end do

END SUBROUTINE produce_site

SUBROUTINE search_neighbors
	implicit none
	integer(4)                  :: dx, dy, dz
	integer(4)                  :: dx_s, dx_e, dy_s, dy_e, dz_s, dz_e
	integer(4)                  :: ix, iy, iz
	integer(4)                  :: sb1, sb2
	integer(4)                  :: ir
	real(8)                     :: dr, dvec(3)
	real(8)                     :: smallnumber
	integer(4)                  :: i

	allocate(Dist(SubLat, MaxNr))
	allocate(NumNeig(SubLat, MaxNr))
	allocate(NeigiCoord(SubLat, MaxNr, 20, 0:3))

	smallnumber = cal_distance(LatVec(1,1:3))/1.d5

	select case(Dimen)
	case(1)
		dx_s = -MaxNr; dx_e = MaxNr
		dy_s = 0; dy_e = 0
		dz_s = 0; dz_e = 0
	case(2)
		dx_s = -MaxNr; dx_e = MaxNr
		dy_s = -MaxNr; dy_e = MaxNr
		dz_s = 0; dz_e = 0
	case(3)
		dx_s = -MaxNr; dx_e = MaxNr
		dy_s = -MaxNr; dy_e = MaxNr
		dz_s = -MaxNr; dz_e = MaxNr
	case default
		write(*,*) "Err: Dim/=1, 2 or 3"
		stop
	end select

	Dist = 1.d10
	NumNeig = 0
	do sb1 = 1, SubLat
		do dz = dz_s, dz_e
		do dy = dy_s, dy_e
		do dx = dx_s, dx_e
		do sb2 = 1, SubLat
			if( dz==0 .and. dy==0 .and. dx==0 .and. sb1==sb2 ) cycle
			dvec = 0.d0
			dvec(1:Dimen) = SubLatVec(sb1,1:Dimen) &
				& - dx*LatVec(1,1:Dimen) &
				& - dy*LatVec(2,1:Dimen) &
				& - dz*LatVec(3,1:Dimen) &
				& - SubLatVec(sb2,1:Dimen)
			dr = cal_distance( dvec(1:3) )

			do ir=1, MaxNr
				if( dr<Dist(sb1,ir) .and. abs(dr-Dist(sb1,ir))>smallnumber ) then
					if( ir<MaxNr ) then
						Dist(sb1,ir+1:MaxNr)    = Dist(sb1,ir:MaxNr-1)
						NumNeig(sb1,ir+1:MaxNr) = NumNeig(sb1,ir:MaxNr-1)
						NeigiCoord(sb1,ir+1:MaxNr,:,:) = NeigiCoord(sb1,ir:MaxNr-1,:,:)
					end if
					Dist(sb1,ir) = dr
					NumNeig(sb1,ir) = 1
					NeigiCoord(sb1,ir,NumNeig(sb1,ir),0:3) = (/sb2,dx,dy,dz/)
					exit
				else if( abs(dr-Dist(sb1,ir))<smallnumber ) then
					NumNeig(sb1,ir) = NumNeig(sb1,ir) + 1
					if( NumNeig(sb1,ir)>20 ) then
						write(*,*) "Err: NumNeig(sb1,ir)>20"
						stop
					end if
					NeigiCoord(sb1,ir,NumNeig(sb1,ir),0:3) = (/sb2,dx,dy,dz/)
					exit
				end if
			end do
		end do
		end do
		end do
		end do
	end do

	open(1,file="Dist.txt",action="write")
	do sb1=1, SubLat
	do ir=1, MaxNr
		write(1,'(1X,I2,1X,I4,F12.6,1X,I4)') sb1, ir, Dist(sb1,ir), NumNeig(sb1,ir)
	end do
	end do
	close(1)

	open(1,file="NeigiCoord.txt",action="write")
	do sb1=1, SubLat
	do ir=1, MaxNr
		do i=1, NumNeig(sb1,ir)
			write(1,'(1X,I2,1X,I4,1X,I4,4I4)') sb1, ir, i, NeigiCoord(sb1, ir, i, 0:3)
		end do
	end do
	end do
	close(1)

	do ir=1, MaxNr
		do sb1=2, SubLat
			if( NumNeig(1,ir)/=NumNeig(sb1,ir) ) then
				write(*,*) "Err: NumNeig", ir
				write(*,*) NumNeig(1:SubLat,ir)
				call write_site
				stop
			end if
		end do
	end do

END SUBROUTINE search_neighbors

SUBROUTINE calculate_bond
	implicit none
	integer(4)                  :: sb1, icell
	integer(4)                  :: ix1, iy1, iz1
	integer(4)                  :: ix2, iy2, iz2
	integer(4)                  :: s1, s2
	integer(4)                  :: ir, i
	integer(4)                  :: nb(MaxNr)

	allocate(Bond(MaxNr,NSite*maxval(NumNeig)/2,2))

	nb = 0
	do ir=1, MaxNr
	do sb1=1, SubLat
		do icell=1, Nsite/SubLat
			ix1 = mod(icell-1, Lx)
			iy1 = mod((icell-1)/Lx, Ly)
			iz1 = (icell-1)/(Lx*Ly)
			s1 = (icell-1)*SubLat + sb1
			!write(*,'(A4,3I8)') "----", ix1, iy1, iz1
			do i=1, NumNeig(sb1,ir)
				ix2 = mod(ix1+NeigiCoord(sb1,ir,i,1)+Lx, Lx)
				iy2 = mod(iy1+NeigiCoord(sb1,ir,i,2)+Ly, Ly)
				iz2 = mod(iz1+NeigiCoord(sb1,ir,i,3)+Lz, Lz)
				s2 = (ix2+iy2*Lx+iz2*Lx*Ly)*SubLat + NeigiCoord(sb1,ir,i,0)
				if( s1<s2 ) then
					nb(ir) = nb(ir)+1
					Bond(ir,nb(ir),1:2) = (/s1,s2/)
					!write(*,'(1X,I4,2I8,A4,3I8,3I4)') ir, s1, s2, "::::", ix2,iy2,iz2, NeigiCoord(sb1,ir,i,1:3)
				end if
			end do
		end do
	end do
	end do
END SUBROUTINE calculate_bond

SUBROUTINE calculate_neighbor
	implicit none
	integer(4)              :: ns(MaxNr,Nsite)
	integer(4)              :: s1, s2
	integer(4)              :: ir, ib
	integer(4)              :: i, j

	allocate(NNst(MaxNr,NSite,maxval(NumNeig)))
	allocate(NNbk(MaxNr,NSite,maxval(NumNeig)))
	! ns: 1st: site number; 2: the n_th neighbors
	ns = 0

	do ir=1, MaxNr
	do ib=1, NSite*NumNeig(1,ir)/2
		s1 = Bond(ir,ib,1)
		s2 = Bond(ir,ib,2)
		ns(ir,s1) = ns(ir,s1) + 1
		ns(ir,s2) = ns(ir,s2) + 1
		if( ns(ir,s1)>NumNeig(1,ir) .or. ns(ir,s2)>NumNeig(1,ir) ) then
			write(*,*) "Err: ns"
			write(*,*) NumNeig(1,1:MaxNr)
			write(*,*) s1, s2
			call write_bond
			stop
		end if
		NNst(ir,s1,ns(ir,s1)) = s2
		NNst(ir,s2,ns(ir,s2)) = s1
		NNbk(ir,s1,ns(ir,s1)) = ns(ir,s2)
		NNbk(ir,s2,ns(ir,s2)) = ns(ir,s1)
	end do
	end do
END SUBROUTINE calculate_neighbor

SUBROUTINE write_lattice
	IMPLICIT NONE

	call write_site
	call write_bond
	call write_nnsite
END SUBROUTINE write_lattice

SUBROUTINE write_site
	IMPLICIT NONE

	integer(4)           :: sb
	integer              :: ir
	character(50)        :: charsite
	integer(4)           :: i, j

	call system('rm -f lattice')
	call system('rm -f latvec.dat')
	call system('rm -f reclatvec.dat')
	call system('rm -f site.dat')

	open(1, file='lattice',action='write')
	write(1,'(1X,5I4,A20)') Dimen, Lx, Ly, Lz, SubLat, trim(NameLat)
	do ir = 1, MaxNr
		write(1,'(1X,I4)') NumNeig(1,ir)
	end do
	close(1)

	open(1, file='latvec.dat',action='write')
	write(1,'(3(1X,F16.8))') LatVec(1,1:3)
	write(1,'(3(1X,F16.8))') LatVec(2,1:3)
	write(1,'(3(1X,F16.8))') LatVec(3,1:3)
	close(1)
	open(1, file='sublatvec.dat',action='write')
	do i=1, SubLat
		write(1,'(3(1X,F16.8))') SubLatVec(i,1:3)
	end do
	close(1)
	open(1, file='reclatvec.dat',action='write')
	write(1,'(3(1X,F16.8))') RecLatVec(1,1:3)
	write(1,'(3(1X,F16.8))') RecLatVec(2,1:3)
	write(1,'(3(1X,F16.8))') RecLatVec(3,1:3)
	close(1)

	!**** write sites ****!
	open(1,file="site.dat",action='write')
	write(1,'(1X, I8)') NSite
	do i=1, NSite
		sb = mod(i-1, SubLat)+1
		write(1,'(1X,I8,1X,I8,1X,I2,1X,I2,1X,3F16.8)') &
			& i, Site(i)%cell, sb, Site(i)%tb, Site(i)%coord(:)
	end do
	close(1)
END SUBROUTINE write_site

SUBROUTINE write_bond
	IMPLICIT NONE

	integer(4)           :: sb
	integer, allocatable :: coubond(:) 
	integer              :: ir
	character(50)        :: charbond
	integer(4)           :: i, j

	call system('rm -f bond*.dat')

	allocate(coubond(NSite))

	!**** write bonds ****!
	do ir = 1, MaxNr
		if( ir<10 ) then
			write(charbond,'(A4,I1,A4)') "bond",ir,".dat"
		else if( ir<100 ) then
			write(charbond,'(A4,I2,A4)') "bond",ir,".dat"
		end if

		open(1,file=trim(charbond),action='write')
		!write(2,'(1X,A5,2X,I5)') 'Bond=', NSite*NumNeig/2
		write(1,'(1X, I8)') NSite*NumNeig(1,ir)/2
		coubond=0

		do i=1, NSite*NumNeig(1,ir)/2
			write(1,'(1X,I8,1X,I8)') Bond(ir,i,1), Bond(ir,i,2)
			coubond(Bond(ir,i,1))=coubond(Bond(ir,i,1))+1
			coubond(Bond(ir,i,2))=coubond(Bond(ir,i,2))+1
		end do
		close(1)

		!*** check neighbor number ****!
		do i=1, NSite
			if( coubond(i) /= NumNeig(1,ir) ) then
				write(*,*) 'neighbor wrong'
			end if
		end do
	end do

	deallocate(coubond)
END SUBROUTINE write_bond

SUBROUTINE write_nnsite
	IMPLICIT NONE

	integer(4)           :: sb
	integer              :: ir
	character(50)        :: charnnsite, charbackdir
	integer(4)           :: i, j

	call system('rm -f nnsite*.dat')

	!**** write neighbor list ****!
	do ir = 1, MaxNr
		if( ir<10 ) then
			write(charnnsite,'(A6,I1,A4)') "nnsite",ir,".dat"
			write(charbackdir,'(A7,I1,A4)') "backdir",ir,".dat"
		else if( ir<100 ) then
			write(charnnsite,'(A6,I2,A4)') "nnsite",ir,".dat"
			write(charbackdir,'(A7,I2,A4)') "backdir",ir,".dat"
		end if
		open(1,file=trim(charnnsite),action='write')
		open(2,file=trim(charbackdir),action='write')

		do j=1, NSite
			write(1,'(I8,$)') j
			write(2,'(I8,$)') j
			do i=1, NumNeig(1,ir)
				write(1,'(I8,$)') NNst(ir,j,i)
				write(2,'(I8,$)') NNbk(ir,j,i)
			end do
			write(1,'()')
			write(2,'()')
		end do
		close(1)
		close(2)
	end do
END SUBROUTINE write_nnsite

INCLUDE "calculation.f90"

INCLUDE "lattice_def/chain.f90"              !  1
INCLUDE "lattice_def/square.f90"             !  2
INCLUDE "lattice_def/triangular.f90"         !  3
INCLUDE "lattice_def/honeycomb.f90"          !  4
INCLUDE "lattice_def/kagome.f90"             !  5
INCLUDE "lattice_def/checkboard.f90"         !  6
INCLUDE "lattice_def/ruby.f90"               !  7
INCLUDE "lattice_def/cubic.f90"              !  8
INCLUDE "lattice_def/octahedra.f90"          !  9
INCLUDE "lattice_def/pyrochlore.f90"         ! 10
INCLUDE "lattice_def/diamond.f90"            ! 11
INCLUDE "lattice_def/packingAB.f90"          ! 12
INCLUDE "lattice_def/packingABC.f90"         ! 13

END PROGRAM produce_lattice
