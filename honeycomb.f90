
!============ Honeycomb ============!
SUBROUTINE honeycomb
	IMPLICIT NONE

	NameLat = 'honeycomb lattice'
	Dimen   = 2
	SubLat  = 2
	NumNeig = 3
	NSite   = Lx*Ly*Lz*SubLat

	!**** For produce sites ****!
	LatVec(1,1) = SQRT3/2.d0
	LatVec(1,2) = -0.5d0
	LatVec(1,3) = 0.d0
	LatVec(2,1) = 0.d0
	LatVec(2,2) = 1.d0
	LatVec(2,3) = 0.d0
	SubLatVec(1,1) = 0.d0
	SubLatVec(1,2) = 0.d0
	SubLatVec(1,3) = 0.d0
	SubLatVec(2,1) = 1.d0/SQRT3/2.d0
	SubLatVec(2,2) = 0.5d0
	SubLatVec(2,3) = 0.d0
	allocate(Site(SubLat,NSite/SubLat))

	!**** For produce bonds ****!
	allocate(Bond(NSite*NumNeig/2,2))
	allocate(Relt(NSite*NumNeig/2,3))

	RETURN
END SUBROUTINE honeycomb
