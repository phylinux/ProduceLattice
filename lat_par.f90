
MODULE lat_par
	IMPLICIT NONE
	character(50)               :: NameLat         ! the name of lattice
	integer(4)                  :: laty
	integer(4)                  :: Lx=1, Ly=1, Lz=1! the size of sublattice
	integer(4)                  :: Lxyz(3)
	integer(4)                  :: NSite
	integer(4)                  :: Dimen           ! the dimension of lattice
	integer(4)                  :: SubLat          ! the number of sublatties
	real(8)                     :: LatVec(3,3)     ! basis vector
	real(8), allocatable        :: SubLatVec(:,:)  ! basis vector of sublattice
	real(8)                     :: RecLatVec(3,3)  ! Reciprocal lattice vector

	type lat_info
		integer(4)       :: cell          ! cell integer location
		real(8)          :: coord(3)      ! site real location
		integer(4)       :: tb=0          ! in which boundary, 0: not in boundary
		                                  ! 1:x  2:y  4:z
	end type

	type(lat_info), allocatable :: Site(:)       ! note the sites' coordinates
	integer(4)                  :: MaxNr         ! the furthest neighbors
	real(8), allocatable        :: Dist(:,:)     ! note distances of different neighbors
	integer(4), allocatable     :: Bond(:,:,:)   ! note the bonds
	integer(4), allocatable     :: NNst(:,:,:)   ! neighbors' list
	integer(4), allocatable     :: NNbk(:,:,:)   ! back direction
	integer(4), allocatable     :: Relt(:,:,:)   ! bond between i and j, Relt means j is in direction of i Relt(,1:3)
	integer(4), allocatable     :: NumNeig(:,:)  ! the number of ith neighbors
	integer(4), allocatable     :: NeigiCoord(:,:,:,:)     ! the relative integer coordination

	!==== some constant ====!
	real(8), parameter          :: PI    = 3.141592653589d0
	real(8), parameter          :: SQRT2 = 2.d0**0.5d0
	real(8), parameter          :: SQRT3 = 3.d0**0.5d0
	real(8), parameter          :: SQRT6 = 6.d0**0.5d0
END MODULE lat_par
