
!============ Octahedra lattice ============!
SUBROUTINE octahedra
	IMPLICIT NONE

	NameLat = 'octahedra'
	Dimen   = 3
	SubLat  = 3
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
	SubLatVec(1,1) = 0.5d0
	SubLatVec(1,2) = 0.d0
	SubLatVec(1,3) = 0.d0
	SubLatVec(2,1) = 0.d0
	SubLatVec(2,2) = 0.5d0
	SubLatVec(2,3) = 0.d0
	SubLatVec(3,1) = 0.d0
	SubLatVec(3,2) = 0.d0
	SubLatVec(3,3) = 0.5d0

	!-- the n_th neighbor's distance --!
	NumNeig = (/8,0,0/)
	NNdt(1) = SQRT2/2.d0*dsqrt(dot_product(LatVec(1,1:3),LatVec(1,1:3)))
	NNdt(2) = 2.d0*dsqrt(dot_product(LatVec(1,1:3),LatVec(1,1:3)))
	NNdt(3) = 2.d0*SQRT2*dsqrt(dot_product(LatVec(1,1:3),LatVec(1,1:3)))

END SUBROUTINE octahedra

