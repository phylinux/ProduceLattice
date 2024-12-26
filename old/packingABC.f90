
!============ Sphere packing ABC ============!
SUBROUTINE packingABC
	IMPLICIT NONE
	real(8)                     :: neighborvec(3)
	real(8)                     :: height

	NameLat = 'packingABC'
	Dimen   = 3
	SubLat  = 3
	NSite   = Lx*Ly*Lz*SubLat

	!A- - - - -
	!|    |
	!C    !
	!|   3h
	!B    |
	!|    |
	!A- - - - -
	
	!       h = sqrt(2/3)=0.816496
	! A-B-C: nearest neighbor
	! height=(2.d0/3.d0)**0.5*3.d0

	! 0.5 < h < sqrt(2/3)=0.816496
	! A-B: nearest neighbor;  A-A: next nearest neighbor
	! height=0.70*3.d0

	! sqrt(2/3)=0.816496 < h < sqrt(5/3)=1.29099
	! A-A: nearest neighbor;  A-B: next nearest neighbor
	  height=0.9d0*3.d0

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
	SubLatVec(2,3) = height/3.d0
	SubLatVec(3,:) = SubLatVec(2,:)*2

	!-- the n_th neighbor's distance --!
	NumNeig = (/6,6,0/)
	if( height/3.d0 < sqrt(2.d0/3.d0) ) then
		neighborvec(1:3) = SubLatVec(2,1:3)
		NNdt(1) = 1.d0*dsqrt(dot_product(neighborvec(1:3),neighborvec(1:3)))
		neighborvec(1:3) = LatVec(2,1:3)
		NNdt(2) = 1.d0*dsqrt(dot_product(neighborvec(1:3),neighborvec(1:3)))
	else
		neighborvec(1:3) = LatVec(2,1:3)
		NNdt(1) = 1.d0*dsqrt(dot_product(neighborvec(1:3),neighborvec(1:3)))
		neighborvec(1:3) = SubLatVec(2,1:3)
		NNdt(2) = 1.d0*dsqrt(dot_product(neighborvec(1:3),neighborvec(1:3)))
	end if

END SUBROUTINE packingABC
