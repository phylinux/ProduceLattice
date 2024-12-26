
!============ Pyrochlore lattice ============!
SUBROUTINE pyrochlore
	IMPLICIT NONE

	NameLat = 'pyrochlore'
	Dimen   = 3
	SubLat  = 4
	NSite   = Lx*Ly*Lz*SubLat

	!**** For produce sites ****!
	LatVec(1,1) = 1.d0
	LatVec(1,2) = 1.d0
	LatVec(1,3) = 0.d0
	LatVec(2,1) = 0.d0
	LatVec(2,2) = 1.d0
	LatVec(2,3) = 1.d0
	LatVec(3,1) = 1.d0
	LatVec(3,2) = 0.d0
	LatVec(3,3) = 1.d0
	SubLatVec(1,1) = 0.d0
	SubLatVec(1,2) = 0.d0
	SubLatVec(1,3) = 0.d0
	SubLatVec(2,1) = 0.5d0
	SubLatVec(2,2) = 0.5d0
	SubLatVec(2,3) = 0.d0
	SubLatVec(3,1) = 0.d0
	SubLatVec(3,2) = 0.5d0
	SubLatVec(3,3) = 0.5d0
	SubLatVec(4,1) = 0.5d0
	SubLatVec(4,2) = 0.d0
	SubLatVec(4,3) = 0.5d0

	!-- the n_th neighbor's distance --!
	!NumNeig = (/6,12,6/)
	NumNeig = (/6,12,12/)
	NNdt(1) = dsqrt(dot_product(SubLatVec(2,1:3),SubLatVec(2,1:3)))
	NNdt(2) = SQRT3*dsqrt(dot_product(SubLatVec(2,1:3),SubLatVec(2,1:3)))
	NNdt(3) = 2.d0*dsqrt(dot_product(SubLatVec(2,1:3),SubLatVec(2,1:3)))

END SUBROUTINE pyrochlore
