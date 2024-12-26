
!============ Honeycomb ============!
SUBROUTINE honeycomb
	IMPLICIT NONE
	real(8)                     :: nvector(3)

	NameLat = 'honeycomb'
	Dimen   = 2
	SubLat  = 2
	Lz      = 1
	NSite   = Lx*Ly*Lz*SubLat

	!-- base vectors and sublattice's vectors --!
	LatVec(1,1) = 1.d0
	LatVec(1,2) = 0.d0
	LatVec(1,3) = 0.d0
	LatVec(2,1) = 0.5d0
	LatVec(2,2) = SQRT3/2.d0
	LatVec(2,3) = 0.d0
	LatVec(3,1) = 0.d0
	LatVec(3,2) = 0.d0
	LatVec(3,3) = 1.d0
	SubLatVec(1,1) = 0.d0
	SubLatVec(1,2) = 0.d0
	SubLatVec(1,3) = 0.d0
	SubLatVec(2,1) = 0.d0
	SubLatVec(2,2) = SQRT3/3.d0
	SubLatVec(2,3) = 0.d0

	!-- the n_th neighbor's distance --!
	NumNeig = (/3,6,3/)
	nvector = SubLatVec(2,1:3)
	NNdt(1) = sqrt(dot_product(nvector,nvector))
	nvector = LatVec(1,1:3)
	NNdt(2) = sqrt(dot_product(nvector,nvector))
	nvector = 2.d0*SubLatVec(2,1:3)
	NNdt(3) = sqrt(dot_product(nvector,nvector))

	!!!TODO: how to search neighbors !!!

END SUBROUTINE honeycomb
