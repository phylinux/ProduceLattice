
!============ Chain  ============!
SUBROUTINE chain
	IMPLICIT NONE

	NameLat = 'chain'
	Dimen   = 1
	SubLat  = 1
	Ly = 1; Lz = 1
	NSite   = Lx*Ly*Lz*SubLat

	!-- base vectors and sublattice's vectors --!
	LatVec(1,1) = 1.d0
	LatVec(1,2) = 0.d0
	LatVec(1,3) = 0.d0
	LatVec(2,1) = 0.d0
	LatVec(2,2) = 1.d0
	LatVec(2,3) = 0.d0
	LatVec(3,1) = 0.d0
	LatVec(3,2) = 0.d0
	LatVec(3,3) = 1.d0
	SubLatVec(1,1) = 0.d0
	SubLatVec(1,2) = 0.d0

	!-- the n_th neighbor's distance --!
	NumNeig = (/2,2,2/)
	NNdt(1) = 1.d0*dsqrt(dot_product(LatVec(1,1:3),LatVec(1,1:3)))
	NNdt(2) = 2.d0*dsqrt(dot_product(LatVec(1,1:3),LatVec(1,1:3)))
	NNdt(3) = 3.d0*dsqrt(dot_product(LatVec(1,1:3),LatVec(1,1:3)))

END SUBROUTINE chain
