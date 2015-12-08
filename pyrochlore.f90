
!============ Pyrochlore lattice ============!
SUBROUTINE pyrochlore
	IMPLICIT NONE

	NameLat = 'pyrochlore lattice'
	Dimen   = 3
	SubLat  = 4
	NumNeig = 6
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
