
SUBROUTINE cal_reclatvec
	IMPLICIT NONE

	integer        :: i
	real(8)        :: rec_vol

	rec_vol = dot_product( LatVec(1,1:3), cross_pro(LatVec(2,1:3),LatVec(3,1:3)) )
	RecLatVec(1,1:3) = 2.d0*PI*cross_pro(LatVec(2,1:3),LatVec(3,1:3))/rec_vol
	RecLatVec(2,1:3) = 2.d0*PI*cross_pro(LatVec(3,1:3),LatVec(1,1:3))/rec_vol
	RecLatVec(3,1:3) = 2.d0*PI*cross_pro(LatVec(1,1:3),LatVec(2,1:3))/rec_vol
	!write(*,*) rec_vol
	!write(*,*) RecLatVec(1,1:3), dot_product(RecLatVec(1,1:3),LatVec(2,1:3)), dot_product(RecLatVec(1,1:3),LatVec(3,1:3)), dot_product(RecLatVec(1,1:3),LatVec(1,1:3))
	!write(*,*) RecLatVec(2,1:3), dot_product(RecLatVec(2,1:3),LatVec(3,1:3)), dot_product(RecLatVec(2,1:3),LatVec(1,1:3)), dot_product(RecLatVec(2,1:3),LatVec(2,1:3))
	!write(*,*) RecLatVec(3,1:3), dot_product(RecLatVec(3,1:3),LatVec(1,1:3)), dot_product(RecLatVec(3,1:3),LatVec(2,1:3)), dot_product(RecLatVec(3,1:3),LatVec(3,1:3))
	!stop

END SUBROUTINE cal_reclatvec

FUNCTION cross_pro( vec1, vec2 )
	IMPLICIT NONE

	real(8)        :: cross_pro(3)
	real(8)        :: vec1(3), vec2(3)

	cross_pro(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
	cross_pro(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
	cross_pro(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)

END FUNCTION cross_pro

FUNCTION cal_distance(vec)
	IMPLICIT NONE
	real(8)          :: cal_distance
	real(8)          :: vec(3)

	cal_distance = sqrt(dot_product(vec(1:3),vec(1:3)))
END FUNCTION cal_distance
