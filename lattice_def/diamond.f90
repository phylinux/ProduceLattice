
!============ Diamond lattice ============!
SUBROUTINE diamond
	IMPLICIT NONE

	NameLat = 'diamond'
	Dimen   = 3
	SubLat  = 2
	NSite   = Lx*Ly*Lz*SubLat

	allocate(SubLatVec(SubLat,3))
	SubLatVec = 0.d0

	LatVec = 0.d0
	!**** For produce sites ****!
	LatVec(1,1) = 0.d0
	LatVec(1,2) = 1.d0
	LatVec(1,3) = 1.d0
	LatVec(2,1) = 1.d0
	LatVec(2,2) = 0.d0
	LatVec(2,3) = 1.d0
	LatVec(3,1) = 1.d0
	LatVec(3,2) = 1.d0
	LatVec(3,3) = 0.d0
	SubLatVec(1,1) = 0.d0
	SubLatVec(1,2) = 0.d0
	SubLatVec(1,3) = 0.d0
	SubLatVec(2,1) = 0.5d0
	SubLatVec(2,2) = 0.5d0
	SubLatVec(2,3) = 0.5d0
END SUBROUTINE diamond
