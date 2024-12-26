
!============ Diamond lattice ============!
SUBROUTINE diamond
	IMPLICIT NONE
	real(8)                     :: neighborvec(3)

	NameLat = 'diamond'
	Dimen   = 3
	SubLat  = 2
	NSite   = Lx*Ly*Lz*SubLat

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

	!-- the n_th neighbor's distance --!
	NumNeig = (/4,12,0/)
	neighborvec(1:3) = SubLatVec(2,1:3)
	NNdt(1) = 1.d0*dsqrt(dot_product(neighborvec(1:3),neighborvec(1:3)))
	neighborvec(1:3) = LatVec(1,1:3)
	NNdt(2) = 1.d0*dsqrt(dot_product(neighborvec(1:3),neighborvec(1:3)))

END SUBROUTINE diamond
