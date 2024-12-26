
SUBROUTINE output_lat
	IMPLICIT NONE

	integer(4)           :: subi
	integer, allocatable :: coubond(:) 
	integer              :: nr
	character(50)        :: charsite, charbond, charnnsite, charbackdir
	integer(4)           :: i, j

	call system('rm -f lattice')
	call system('rm -f reclatvec.dat')
	call system('rm -f latvec.dat')
	call system('rm -f site.dat')
	call system('rm -f bond*.dat')
	call system('rm -f nnsite*.dat')

	allocate(coubond(NSite))

	open(1, file='lattice',action='write')
	write(1,'(1X,5I4,A20)') Dimen, Lx, Ly, Lz, SubLat, trim(NameLat)
	do nr = 1, NNr
		write(1,'(1X,I4)') NumNeig(nr)
	end do
	close(1)

	open(1, file='latvec.dat',action='write')
	write(1,'(3(1X,F16.8))') LatVec(1,1:3)
	write(1,'(3(1X,F16.8))') LatVec(2,1:3)
	write(1,'(3(1X,F16.8))') LatVec(3,1:3)
	close(1)
	open(1, file='sublatvec.dat',action='write')
	do i=1, SubLat
		write(1,'(3(1X,F16.8))') SubLatVec(i,1:3)
	end do
	close(1)
	open(1, file='reclatvec.dat',action='write')
	write(1,'(3(1X,F16.8))') RecLatVec(1,1:3)
	write(1,'(3(1X,F16.8))') RecLatVec(2,1:3)
	write(1,'(3(1X,F16.8))') RecLatVec(3,1:3)
	close(1)

	!**** output sites ****!
	open(1,file="site.dat",action='write')
	write(1,'(1X, I8)') NSite
	do i=1, NSite/SubLat
		do subi=1, SubLat
			write(1,'(1X,I8,1X,I8,1X,I2,1X,I2,1X,3I4,1X,3F16.8)') &
				Site(subi,i)%sn, i, subi, Site(subi,i)%tb, Site(subi,i)%icoord(:), Site(subi,i)%coord(:)
		end do
	end do
	close(1)

	!**** output bonds ****!
	!do nr = 1, 1
	do nr = 1, NNr
		if( NumNeig(nr)==0 ) cycle
		write(charbond,'(A4,I1,A4)') "bond",nr,".dat"

		open(1,file=trim(charbond),action='write')
		!write(2,'(1X,A5,2X,I5)') 'Bond=', NSite*NumNeig/2
		write(1,'(1X, I8)') NSite*NumNeig(nr)/2
		coubond=0

		do i=1, NSite*NumNeig(nr)/2
			write(1,'(1X,I8,2X,I8,1X,3I3)') Bond(nr,i,1), Bond(nr,i,2), Relt(nr,i,1:3)
			coubond(Bond(nr,i,1))=coubond(Bond(nr,i,1))+1
			coubond(Bond(nr,i,2))=coubond(Bond(nr,i,2))+1
		end do
		close(1)

		!*** check neighbor number ****!
		do i=1, NSite
			if( coubond(i) /= NumNeig(nr) ) then
				write(*,*) 'neighbor wrong'
			end if
		end do
	end do

	!**** output neighbor list ****!
	!do nr = 1, 1
	do nr = 1, NNr
		if( NumNeig(nr)==0 ) cycle
		write(charnnsite,'(A6,I1,A4)') "nnsite",nr,".dat"
		write(charbackdir,'(A7,I1,A4)') "backdir",nr,".dat"
		open(1,file=trim(charnnsite),action='write')
		open(2,file=trim(charbackdir),action='write')

		do j=1, NSite
			write(1,'(I8,$)') j
			write(2,'(I8,$)') j
			do i=1, NumNeig(nr)-1
				write(1,'(I8,$)') NNst(nr,j,i)
				write(2,'(I8,$)') NNbk(nr,j,i)
			end do
			i=NumNeig(nr)
			write(1,'(I8)') NNst(nr,j,i)
			write(2,'(I8)') NNbk(nr,j,i)
		end do
		close(1)
		close(2)
	end do

	deallocate(coubond)

END SUBROUTINE output_lat

!SUBROUTINE output_lat_py
!	IMPLICIT NONE
!
!	integer              :: i, subi
!	integer              :: s1,s2,sb1,sb2,si1,si2
!	integer, allocatable :: coubond(:) 
!	integer              :: nbcount
!	real                 :: distofsite
!
!	if( SubLat==1 ) then
!		distofsite = cal_distance(LatVec(1,1:3))
!	else
!		distofsite = cal_distance(SubLatVec(2,1:3))
!	end if
!
!	call system('rm -f coordinate.dat')
!	open(1,file='coordinate.dat',access='append')
!
!	write(1,*) "{'Interaction': [],"
!	write(1,*) "'Lines': ["
!	nbcount = NSite*NumNeig/2
!	do i=1, nbcount
!		s1  = Bond(i,1)
!		s2  = Bond(i,2)
!		si1 = (s1-1)/SubLat+1
!		si2 = (s2-1)/SubLat+1
!		!Site(subi,num)%sn = num*SubLat-SubLat+subi
!		sb1 = mod(s1-1,SubLat)+1
!		sb2 = mod(s2-1,SubLat)+1
!		if( Site(sb1,si1)%sn/=s1 ) stop "s1 wrong"
!		if( Site(sb2,si2)%sn/=s2 ) stop "s2 wrong"
!		if( abs(cal_distance(Site(sb1,si1)%coord-Site(sb2,si2)%coord)-distofsite)<1.d-4 ) then
!			write(1,'("[(",I4,",",I4,"),",I3,"],")') Bond(i,1)-1, Bond(i,2)-1,0
!		end if
!	end do
!	close(1)
!
!
!	open(1,file='coordinate.dat',access='append')
!	write(1,*) "'Points': ["
!	!**** output sites ****!
!	do i=1, NSite/SubLat
!		do subi=1, SubLat
!			write(1,'(1X,"[(",F6.3,",",F6.3,",",F6.3,"),", &
!				"(",I3,",",I3,",",I3,"),", I2,"],")') &
!				Site(subi,i)%coord(1),Site(subi,i)%coord(2),Site(subi,i)%coord(3), &
!				Site(subi,i)%icoord(1)-1,Site(subi,i)%icoord(2)-1,Site(subi,i)%icoord(3)-1, &
!				subi
!		end do
!	end do
!	close(1)
!
!	RETURN
!END SUBROUTINE output_lat_py
