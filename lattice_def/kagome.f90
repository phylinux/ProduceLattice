
!============ Kagome lattice ============!
SUBROUTINE kagome
	IMPLICIT NONE

	NameLat = 'kagome'
	Dimen   = 2
	SubLat  = 3
	NSite   = Lx*Ly*Lz*SubLat

	allocate(SubLatVec(SubLat,3))
	SubLatVec = 0.d0

	LatVec = 0.d0
	!-- base vectors and sublattice's vectors --!
	LatVec(1,1) = 1.d0
	LatVec(1,2) = 0.d0
	LatVec(1,3) = 0.d0
	LatVec(2,1) = 0.5d0
	LatVec(2,2) = SQRT3/2
	LatVec(2,3) = 0.d0
	LatVec(3,1) = 0.d0
	LatVec(3,2) = 0.d0
	LatVec(3,3) = 1.d0
	SubLatVec(1,1) = 0.d0
	SubLatVec(1,2) = 0.d0
	SubLatVec(1,3) = 0.d0
	SubLatVec(2,1:3) = LatVec(1,1:3)/2.d0
	SubLatVec(3,1:3) = LatVec(2,1:3)/2.d0
END SUBROUTINE kagome
