
SUBROUTINE pl_lat
	use lat_par
	IMPLICIT NONE

	!**** setting ****!
	real        :: wi=600.0, hi=600.0
	integer     :: BLACK=0, RED=1, GREEN=2, YELLOW=3, BLUE=4
	real        :: xfact, yfact
	real        :: xmin, xmax, ymin, ymax, zmin, zmax

	integer     :: i, subi
	real        :: ix, iy, iz
	real, allocatable :: sl(:,:)
	real        :: len
	character(10)  :: num

	allocate(sl(NSite,3))
	sl=0.d0

	xmin=-1.0
	ymin=-2.0
	xmax= 1.0*Lx*2
	ymax= 1.0*Ly*2
	zmin=-1.0
	zmax= 1.d0*Lz*2
	call prefsize(int(wi+36), int(hi+36))
	call vinit('X11')
	!call getfactors(xfact,yfact)
	call viewport(-1.0,1.0,-1.0,1.0)
	call ortho(xmin, xmax, ymin, ymax, zmin, zmax)
	call swapbuffers()
	call clear

	!if( Dimen==2 ) then
	!	call lookat(0.0,0.0,0.0, 0.0,0.0,1.0, 0.0)
	!else if( Dimen==3 ) then
	!	call lookat(0.0,0.0,0.0, 0.2,0.4,0.2, 0.0)
	!end if

	call color(RED)
	call polyfill(.TRUE.)
	do i=1, NSite/SubLat
		do subi=1, SubLat
			ix = Site(subi,i)%coord(1)
			iy = Site(subi,i)%coord(2)
			iz = Site(subi,i)%coord(3)
			write(num,'(I3)') Site(subi,i)%sn
			call move(ix,iy,iz)
			call drawstr(num)
			call circle(ix,iy,0.1)
			sl(Site(subi,i)%sn,:) = Site(subi,i)%coord(:)
		end do
	end do

	open(1,file='distant.dat',access='append')
	call color(GREEN)
	do i=1, NSite*NumNeig/2
		call move(sl(Bond(i,1),1),sl(Bond(i,1),2),sl(Bond(i,1),3))
		call draw(sl(Bond(i,2),1),sl(Bond(i,2),2),sl(Bond(i,2),3))
		len = sum((sl(Bond(i,1),:)-sl(Bond(i,2),:))**2)
		write(1,'(I3,2X,I3,2X,F12.4)') Bond(i,1), Bond(i,2), len
	end do
	close(1)

	call swapbuffers()
	call getkey()
	call vexit()
	RETURN
END SUBROUTINE pl_lat
