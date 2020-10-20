
!============ Sphere packing AB ============!
SUBROUTINE packingAB
	IMPLICIT NONE
	real(8)                     :: neighborvec(3)
	real(8)                     :: height

	NameLat = 'packingAB'
	Dimen   = 3
	SubLat  = 2
	NSite   = Lx*Ly*Lz*SubLat

	!A- - - - -
	!|    |
	!B    h
	!|    |
	!A- - - - -
	
	! 1<(h/2)^2+(1/sqrt(3))^2<(sqrt(3))^2
	! sqrt(8/3)=1.633 < h < sqrt(32/3)=3.266
	! A-A: nearest neighbor;  A-B: next nearest neighbor
	! height=2.d0

	! (h/2)^2+(1/sqrt(3))^2<1
	! h < sqrt(8/3)=1.633  && 1.0 < h
	! A-A: nearest neighbor;  A-B: next nearest neighbor
	  height=1.5d0

	!-- base vectors and sublattice's vectors --!
	LatVec(1,1) = 1.d0
	LatVec(1,2) = 0.d0
	LatVec(1,3) = 0.d0
	LatVec(2,1) = 0.5d0
	LatVec(2,2) = SQRT3/2
	LatVec(2,3) = 0.d0
	LatVec(3,1) = 0.d0
	LatVec(3,2) = 0.d0
	LatVec(3,3) = height
	SubLatVec(1,1) = 0.d0
	SubLatVec(1,2) = 0.d0
	SubLatVec(1,3) = 0.d0
	SubLatVec(2,1) = 0.5d0
	SubLatVec(2,2) = 0.5d0/SQRT3
	SubLatVec(2,3) = height/2.d0

	!-- the n_th neighbor's distance --!
	NumNeig = (/6,6,0/)
	neighborvec(1:3) = SubLatVec(2,1:3)
	NNdt(1) = 1.d0*dsqrt(dot_product(neighborvec(1:3),neighborvec(1:3)))
	neighborvec(1:3) = LatVec(2,1:3)
	NNdt(2) = 1.d0*dsqrt(dot_product(neighborvec(1:3),neighborvec(1:3)))

END SUBROUTINE packingAB
