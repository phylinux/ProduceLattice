
MODULE lat_par
	IMPLICIT NONE
	character(20)        :: NameLat       ! the name of lattice
	integer              :: Lx, Ly, Lz    ! the size of sublattice
	integer              :: Lxyz(3)
	integer              :: NSite
	integer              :: Dimen         ! the dimension of lattice
	integer              :: SubLat        ! the number of sublatties
	real                 :: LatVec(3,3)   ! basis vector
	real                 :: SubLatVec(6,3)! basis vector of sublattice

	type la_si
		integer          :: sn            ! site number
		real             :: coord(3)      ! site location
		integer          :: icoord(3)     ! site location
		integer          :: tb=0          ! in which boundary, 0: not in boundary
		                                  ! 1:x  2:y  4:z
	end type

	type(la_si), allocatable   :: Site(:,:)     ! note the sites' coordinates
	integer, allocatable :: Bond(:,:)     ! note the bonds
	integer, allocatable :: Relt(:,:)     ! bond between i and j, Relt means j is in direction of i Relt(,1:3)
	integer              :: NumNeig       ! the number of neighbors

	!==== some constant ====!
	real, parameter   :: PI    = 3.141592653589d0
	real, parameter   :: SQRT2 = 2.d0**0.5d0
	real, parameter   :: SQRT3 = 3.d0**0.5d0
	real, parameter   :: SQRT6 = 6.d0**0.5d0
END MODULE lat_par
