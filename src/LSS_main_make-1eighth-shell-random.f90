
program main 

use mpi
use LSS_cosmo_funs

	implicit none
	! HR3 file, data
        integer :: i, iran, nran, ranseed
        real(dl) :: rmin, rmax, r,x,y,z
        character(len=char_len) :: randatafile, tmpstr

	print *, '## Generating randoms uniformly distributed within 1/8 shell'
        i = iargc()
        if(i.ne.5) then
        	print *, 'Useage:'
        	print *, '  EXE randatafile ranseed nran rmin rmax'
        	stop
        endif
        
        call getarg(1,randatafile)
        call getarg(2,tmpstr)
        read(tmpstr,*) ranseed
        randatafile = trim(adjustl(randatafile))//trim(adjustl(tmpstr))
        call getarg(3,tmpstr)
        read(tmpstr,*) nran
        call getarg(4,tmpstr)
        read(tmpstr,*) rmin
        call getarg(5,tmpstr)
        read(tmpstr,*) rmax
        
        write(*,'(i9,A,f10.2,A,f10.2,A,A)') nran, ' randoms in 1/8 shell, r ', rmin, ' to', rmax, ', file: ', trim(adjustl(randatafile))
        
        call random_seed(put=(/0,ranseed/))

        iran = 0
        
        open(unit=1,file=randatafile)
!        open(unit=2,file=trim(adjustl(randatafile))//'.1percent')
        do while(iran.lt.nran)
	        call random_number(x)
	        call random_number(y)
	        call random_number(z)
	        x = x*rmax
	        y = y*rmax
	        z = z*rmax
	        r = sqrt(x*x+y*y+z*z)
	        if(r<rmax .and. r>rmin) then
	        	write(1,'(3(e15.7,1x))') x,y,z
	        	iran = iran+1
!	                if(mod(iran,100).eq.1) then
!	        	        write(2,'(3(e15.7,1x))') x,y,z
 !       		endif
	        endif
	enddo
	close(1); 
        !close(2)
end program main
