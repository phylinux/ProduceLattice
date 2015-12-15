
!============ Chain  ============!
SUBROUTINE chain
	IMPLICIT NONE

	NameLat = 'onechain'
	Dimen   = 1
	SubLat  = 1
	NumNeig = 2
	NSite   = Lx*Ly*Lz*SubLat

	!**** For produce sites ****!
	LatVec(1,1) = 1.d0
	LatVec(1,2) = 0.d0
	SubLatVec(1,1) = 0.d0
	SubLatVec(1,2) = 0.d0
	allocate(Site(SubLat,NSite/SubLat))

	!**** For produce bonds ****!
	allocate(Bond(NSite*NumNeig/2,2))
	allocate(Relt(NSite*NumNeig/2,3))

	RETURN
END SUBROUTINE chain
