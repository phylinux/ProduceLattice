
!============ Kagome lattice ============!
SUBROUTINE kagome
	IMPLICIT NONE
	real(8)                     :: nvector(3)

	NameLat = 'kagome'
	Dimen   = 2
	SubLat  = 3
	Lz      = 1
	NSite   = Lx*Ly*Lz*SubLat

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

	!-- the n_th neighbor's distance --!
	!NumNeig = (/4,4,2/)
	NumNeig = (/4,4,6/)
	nvector = SubLatVec(2,1:3)
	NNdt(1) = sqrt(dot_product(nvector,nvector))
	nvector = SubLatVec(2,1:3)*SQRT3
	NNdt(2) = sqrt(dot_product(nvector,nvector))
	nvector = SubLatVec(2,1:3)*2.d0
	NNdt(3) = sqrt(dot_product(nvector,nvector))

END SUBROUTINE kagome
