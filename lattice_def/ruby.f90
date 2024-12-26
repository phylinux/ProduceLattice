
!============ Ruby lattice ============!
SUBROUTINE ruby
	IMPLICIT NONE

	NameLat = 'ruby'
	Dimen   = 2
	SubLat  = 6
	NSite   = Lx*Ly*Lz*SubLat

	allocate(SubLatVec(SubLat,3))
	SubLatVec = 0.d0

	LatVec = 0.d0
	!-- base vectors and sublattice's vectors --!
	LatVec(1,1:3) = (/1.d0, 0.d0, 0.d0/)
	LatVec(2,1:3) = (/0.5d0, SQRT3/2, 0.d0/)
	LatVec(3,1:3) = (/0.d0, 0.d0, 1.d0/)
	SubLatVec(1,1:3) = (/0.d0, 1.d0, 0.d0/)
	SubLatVec(2,1:3) = (/-SQRT3/2.d0,  0.5d0, 0.d0/)
	SubLatVec(3,1:3) = (/-SQRT3/2.d0, -0.5d0, 0.d0/)
	SubLatVec(4,1:3) = -SubLatVec(1,1:3)
	SubLatVec(5,1:3) = -SubLatVec(2,1:3)
	SubLatVec(6,1:3) = -SubLatVec(3,1:3)
	SubLatVec = SubLatVec * SQRT3/4.d0
END SUBROUTINE ruby
