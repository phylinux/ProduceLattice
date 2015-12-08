
!============ Triangular lattice ============!
SUBROUTINE triangular
	IMPLICIT NONE

	NameLat = 'triangular lattice'
	Dimen   = 2
	SubLat  = 1
	NumNeig = 6
	NSite   = Lx*Ly*Lz*SubLat

	!**** For produce sites ****!
	LatVec(1,1) = 1.d0
	LatVec(1,2) = 0.d0
	LatVec(2,1) = 0.5d0
	LatVec(2,2) = SQRT3/2
	SubLatVec(1,1) = 0.d0
	SubLatVec(1,2) = 0.d0
	allocate(Site(SubLat,NSite/SubLat))

	!**** For produce bonds ****!
	allocate(Bond(NSite*NumNeig/2,2))
	allocate(Relt(NSite*NumNeig/2,3))

	RETURN
END SUBROUTINE triangular