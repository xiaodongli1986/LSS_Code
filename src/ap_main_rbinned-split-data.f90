
! July 23: This is an old code. We can:
!	# test it;
!	# update it (distance rather than redshift split; no consideration of ra/dec split and just redshift split; ...)

!####################################
! scan chisqs
!####################################
module functions

use ap_tools
implicit none
contains

	integer function find_ibin_rbinned(binned_z,  nbin, z)
		! Dummy
		integer :: nbin
		real(dl) :: binned_z(nbin), z
		! Local
		integer :: i1
		i1 = 1
		do while(i1.le.nbin)
                        if(binned_z(i1).le.z) then
	        		i1 = i1+1
                        else
                                exit
                        endif
		enddo
		find_ibin_rbinned= i1-1
	end function find_ibin_rbinned

end module functions


program main
use functions
use ap_tools
implicit none

      character(len=char_len) :: datafile, outputinfo
      
      real(dl) :: x,y,z, ra,dec,r

      !###################################################
      ! Settings: ra, dec, red bins, number & edges
      integer :: numrbin
      real(dl) :: rcutlist(1000)
      logical :: routputlist(1000)
      !###################################################

      integer :: fileunits(1000), binnedr_num(1000), &
              ndat, numnotinbin, i,j, rbin, fileunit
      character(len=char_len) :: filenames(1000), tmpstr1, tmpstr2, tmpstr3, tmpstr

      ! first check of #-arg
      i = iargc()
      if(i.lt.2) then
              print *, 'Usage:             EXE datafile numrbin edge1 edge2 ... '
              print *, 'fmt of datafile:   x,y,z,***'
              print *, 'Example: ./Rbinned-Split-Data   haharandom10  2  0.0 0.500 1.000 T T ## 2bins; fist 0-0.5; second 0.5-1; do output for both first and second '
              stop
      endif
      
      ! read in: datafile, number of bins
      call getarg(1,datafile)
      write(*,'(A,A)') 'input datafile: ', trim(adjustl(datafile))
      print *, '		(must in fmt: x,y,z, ***) '
      call getarg(2,tmpstr)
      read(tmpstr,*) numrbin

      outputinfo = trim(adjustl(datafile))//'.rbinned-split-data.info'
      open(unit=100,file=outputinfo)

      ! check number of argc
      i = iargc()
      if(i.ne.2+2*numrbin+1) then
              print *, 'ERROR: wrong number of iarg: ', i, 2+numrbin+1
              print *, 'Usage:             ./split_data datafile numrbin edge1 edge2 ... '
              print *, 'fmt of datafile:   x,y,z,***'
	      print *, 'Example: ./Rbinned-Split-Data   haharandom10  2  0.0 0.500 1.000 T T ## 2bins; fist 0-0.5; second 0.5-1; do output for both first and second '
              stop
      endif

      do i = 1, numrbin+1
      	call getarg(2+i,tmpstr)
      	read(tmpstr,*) rcutlist(i)
      enddo
      do i = 1, numrbin
         call getarg(2+numrbin+1+i,tmpstr)
         read(tmpstr,*) routputlist(i)
      enddo

      print *, 'Information also output to file: ', trim(adjustl(outputinfo))
      write(*,'(A,i3,A,(<numrbin+1>(f14.5,1x)))')  '   Edges of ',numrbin,' r   bins: ',   rcutlist(1:numrbin+1)
      write(*,'(A,(<numrbin>(L3)))')  '   Whether output them: ',   routputlist(1:numrbin)
      write(100,'(A,i3,A,(<numrbin+1>(f14.5,1x)))')  '   Edges of ',numrbin,' r   bins: ',   rcutlist(1:numrbin+1)
      write(100,'(A,(<numrbin>(L3)))')  '   Whether output them: ',   routputlist(1:numrbin)
      
      print *, 'Creating files saving the results...'
      write(100,'(A)') 'Creating files saving the results...'
      fileunit = 1099
      do rbin = 1, numrbin
		fileunits(rbin) = fileunit 
		write(tmpstr1,'(f10.2)') rcutlist(rbin)
		write(tmpstr2,'(f10.2)') rcutlist(rbin+1)
                filenames(rbin) =  trim(adjustl(datafile))//'.r'//trim(adjustl(tmpstr1))//'to'//trim(adjustl(tmpstr2))
                if(routputlist(rbin)) then
	                open(unit=fileunits(rbin),file=filenames(rbin))
        	        write(*,'(A,A)') '   ', trim(adjustl(filenames(rbin)))
        	        write(100,'(A,A)') '   ', trim(adjustl(filenames(rbin)))
        	endif
                fileunit = fileunit + 1
      enddo

      binnedr_num = 0
      ndat = 0
      numnotinbin = 0
      open(unit=1,file=datafile)
      do while(.true.)
            read(1,'(A)',end=100) tmpstr
            read(tmpstr,*) x,y,z
!            call xyz_to_radecr_dlpre(x,y,z,ra,dec,r)
	    r = sqrt(x*x+y*y+z*z)
            ndat = ndat+1
            ! ra/dec bin            ! redshift bin
            if(r>rcutlist(1).and.r<rcutlist(numrbin+1)) then
                    rbin = find_ibin_rbinned(rcutlist, numrbin+1, r)
                    if(r>rcutlist(rbin+1).or.r<rcutlist(rbin)) then
                    	print *, 'ERROR: r not in bin!: ', rbin, r, rcutlist(rbin), rcutlist(rbin+1)
                    	stop
                    endif
                    if(routputlist(rbin)) then
                    	write(fileunits(rbin),'(A)') trim(adjustl(tmpstr))
                    endif
                    binnedr_num(rbin) = binnedr_num(rbin)+1
            else
            	    numnotinbin = numnotinbin+1
            endif
            cycle
100         exit            
      enddo
      close(1)
      print *, ' tot #-gal: ', ndat, '; # not in any bin: ', numnotinbin
      write(100,*) ' tot #-gal: ', ndat, '; # not in any bin: ', numnotinbin
      
      do rbin = 1, numrbin
                     print *, '#-gal in bin ', rbin, ':', binnedr_num(rbin)
                     write(100,*), '#-gal in bin ', rbin, ':', binnedr_num(rbin)
      enddo
      close(100)
end program main
