
program main

use LSS_tools

implicit none

	real(dl) :: rmin, rmax, massmin, massmax, x,y,z,mass,r, tmp(1000), rmin2, rmax2, massmin2, massmax2, mass1,mass2,masscut, &
		xmin,ymin,zmin,xmax,ymax,zmax
	real(dl), allocatable :: redges(:), massedges(:)
	integer :: numarg, numrbin,nummassbin, xcol,ycol,zcol,masscol,maxcol, i,j,nowrbin,nowmassbin
	integer(8) :: nlines, nlines_degrade, numoutrange, numbigmass, numdegrade, num,num1,num2
	integer(8), allocatable :: binnednum(:,:)
	logical :: logmass, dodegrade
	character(len=char_len) :: inputfilename, outputfilename, outputfilename1, outputfilename2, printstr, tmpstr1,tmpstr2
	
	print *
	printstr = "Usage: EXE -input intpufilename -rmin rmin "//&
		'-rmax rmax -massmin massmin -massmax massmax '//&
		'-logmass logmass -numrbin numrbin -nummassbin nummassbin '//&
		'-xcol xcol -ycol ycol -zcol zcol -masscol masscol '//&
		'-dodegrade dodegrade -numdegrade numdegrade '//&
		'### dodegrade will choose a suitable masscut and '//&
		'degrade the data to a subsample with number numdegrade'

	! Default values
	rmin = 0.0; rmax = 100000.0d0; 
	massmin = 1.0d10; massmax = 1.0d16; logmass = .true.
	xcol=1; ycol=2; zcol=3; masscol=4
	numrbin = 1; nummassbin = 100
	
	numarg = iargc()
	if(numarg.le.1) then
		write(*,'(A)') trim(adjustl(printstr))
		stop
	endif
	
	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfilename
		elseif(trim(adjustl(tmpstr1)).eq.'-rmin') then
			read(tmpstr2,*) rmin
		elseif(trim(adjustl(tmpstr1)).eq.'-rmax') then
			read(tmpstr2,*) rmax
		elseif(trim(adjustl(tmpstr1)).eq.'-numrbin') then
			read(tmpstr2,*) numrbin
		elseif(trim(adjustl(tmpstr1)).eq.'-massmin') then
			read(tmpstr2,*) massmin
		elseif(trim(adjustl(tmpstr1)).eq.'-massmax') then
			read(tmpstr2,*) massmax
		elseif(trim(adjustl(tmpstr1)).eq.'-logmass') then
			read(tmpstr2,*) logmass
		elseif(trim(adjustl(tmpstr1)).eq.'-nummassbin') then
			read(tmpstr2,*) nummassbin
		elseif(trim(adjustl(tmpstr1)).eq.'-xcol') then
			read(tmpstr2,*) xcol
		elseif(trim(adjustl(tmpstr1)).eq.'-ycol') then
			read(tmpstr2,*) ycol
		elseif(trim(adjustl(tmpstr1)).eq.'-zcol') then
			read(tmpstr2,*) zcol
		elseif(trim(adjustl(tmpstr1)).eq.'-masscol') then
			read(tmpstr2,*) masscol
		elseif(trim(adjustl(tmpstr1)).eq.'-dodegrade') then
			read(tmpstr2,*) dodegrade
		elseif(trim(adjustl(tmpstr1)).eq.'-numdegrade') then
			read(tmpstr2,*) numdegrade
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	maxcol = max(xcol,ycol,zcol,masscol)
	if(maxcol>size(tmp)) then
		print *, 'Overflow: increase size of tmp!'
		print *, '  maxcol, size(tmp) = ', maxcol, size(tmp)
		stop
	endif
	outputfilename = trim(adjustl(inputfilename))//'.rmassinfo'
	
	print *, '#####################################'
	print *, ' Settings'
	print *, '  r: range, numrbin = ', real(rmin), real(rmax), numrbin
	print *, '  mass: range, numbin, logmass = ', real(massmin), real(massmax), nummassbin, logmass
	print *, '  cols of x,y,z,mass: ', xcol,ycol,zcol,masscol
	write(*,'(A,A)') '   inputfilename:  ', trim(adjustl(inputfilename))
	write(*,'(A,A)') '   outputfilename: ', trim(adjustl(outputfilename))
	if(dodegrade) then
		print *, '  dodegrade with #-goal: ', numdegrade
	endif
	print *, '#####################################'


	allocate(redges(numrbin+1),massedges(nummassbin+1),binnednum(numrbin,nummassbin))
	
	do i = 1, numrbin+1
		redges(i) = rmin + (rmax-rmin)/(dble(numrbin))*(i-1)
	enddo
	
	do i = 1, nummassbin+1
		if(logmass) then
			massedges(i) = exp(log(massmin) + (log(massmax)-log(massmin))/(dble(nummassbin))*(i-1))
		else
			massedges(i) = massmin + (massmax-massmin)/(dble(nummassbin))*(i-1)
		endif
	enddo
	
	write(*,*)     '  Edges of r:    ', redges(1:numrbin+1)
	write(*,*) '  Edges of mass: ', massedges(1:nummassbin+1)
	
	binnednum = 0
	
	nlines =0 
	numoutrange =0
	numbigmass =0
	rmin2=1.0e30;rmax2=-rmin2;
	xmin = 1.0e30; ymin = 1.0e30; zmin = 1.0e30
	xmax = -xmin; ymax = -ymin; zmax = -zmin
	massmin2=rmin2;massmax2=-massmin2
	open(unit=1,file=inputfilename,action='read')
	do while(.true.)
		read(1,*,end=101) tmp(1:maxcol)
		nlines = nlines+1
		x=tmp(xcol); y=tmp(ycol); z=tmp(zcol); mass=tmp(masscol)
		r = sqrt(x*x+y*y+z*z)
		rmin2=min(r,rmin2); rmax2=max(r,rmax2)
		xmin = min(x,xmin); xmax = max(x,xmax)
		ymin = min(y,ymin); ymax = max(y,ymax)
		zmin = min(z,zmin); zmax = max(z,zmax)
		massmin2=min(mass,massmin2); massmax2=max(mass,massmax2)
		
		if(numrbin<100) then
			nowrbin = find_ibin(r,redges,numrbin+1)
		else
			nowrbin = find_ibin_2split(r,redges,numrbin+1)
		endif
		if(nummassbin<100) then
			nowmassbin = find_ibin(mass,massedges,nummassbin+1)
		else
			nowmassbin = find_ibin_2split(mass,massedges,nummassbin+1)
		endif

		if(nowrbin.ge.1.and.nowrbin.le.numrbin.and.nowmassbin.ge.1.and.nowmassbin.le.nummassbin) then
			binnednum(nowrbin,nowmassbin) = binnednum(nowrbin,nowmassbin)+1
		else
			numoutrange = numoutrange+1
		endif
		
		if(mass>massedges(nummassbin+1)) numbigmass=numbigmass+1
			
		cycle
101		exit
	enddo		
	close(1)
	
	write(*,'(A,i10,A,i10)') ' Finishing processing ', nlines, ' lines; numoutrange = ', numoutrange
	write(*,'(A,3("(",2f15.7,"); "), A,2f15.7,A,2e15.7,A,i10,A,i10)') &
		'Done. range of x,y,z: ', xmin,xmax,ymin,ymax,zmin,zmax, '; min/max of r:', real(rmin2),real(rmax2), &
		';  min/max of mass:', real(massmin2),real(massmax2), &
		';  numoutrange = ', numoutrange, '; numbigmass = ', numbigmass
	open(unit=2,file=outputfilename)
	print *, 'Result: rmin, rmax, massmin, massmax, num'
	write(2,'(A,3("(",2f15.7,"); "), A,2f15.7,A,2e15.7,A,i10,A,i10)') &
		'#Result: rmin, rmax, massmin, massmax, num. p.s. range of x,y,z: ',&
		xmin,xmax,ymin,ymax,zmin,zmax, 'min/max of r:', real(rmin2),real(rmax2), &
		';  min/max of mass:', real(massmin2),real(massmax2), &
		';  numoutrange = ', numoutrange, '; numbigmass = ', numbigmass
	do i = 1, numrbin
		do j = 1, nummassbin
			write(*,'(2f15.7,2e15.7, i10)') redges(i), redges(i+1), massedges(j), massedges(j+1), binnednum(i,j)
			write(2,'(2f15.7,2e15.7, i10)') redges(i), redges(i+1), massedges(j), massedges(j+1), binnednum(i,j)
		enddo
	enddo
	close(1)
	
	if(.not.dodegrade) stop

	outputfilename1 = trim(adjustl(inputfilename))//'.degraded'
	outputfilename2 = trim(adjustl(inputfilename))//'.degraded.info'	
	if(numdegrade.gt.nlines) then
		print *, 'numdegrade larger than nline: no need to do degrade! Programme will just copy file!'
		open(unit=3,file=outputfilename2)
		write(3,'(A,A,A,i10,A,i10,A,A)') '# This is for file ', trim(adjustl(outputfilename)), '. This file has ',nlines,&
			' lines. We want to degrade it to ',numdegrade,' halos. New sample written to: ', trim(adjustl(outputfilename1))
		write(3,'(A)') '# numdegrade larger than nline: no need to do degrade! Programme will just copy file!'
		close(3)
		call system('cp '//trim(adjustl(inputfilename))//' '//trim(adjustl(outputfilename1)))
		
		stop
	endif
	
	do j = nummassbin, 1, -1
		num = sum(binnednum(1:numrbin,j:nummassbin)) + numbigmass
		if(num>numdegrade) exit
	enddo
	
	mass2 = massedges(j); num2 = sum(binnednum(1:numrbin,j:nummassbin)) + numbigmass
	
	do i = j+1, nummassbin
		num = sum(binnednum(1:numrbin,i:nummassbin)) + numbigmass
		if(num.ne.num2) exit
	enddo
	mass1 = massedges(i); num1 = sum(binnednum(1:numrbin,i:nummassbin)) + numbigmass
	masscut = mass1 + dble(numdegrade-num1)/dble(num2-num1)*(mass2-mass1)
	
	write(*,'(1x,A,e15.7)') 'Do degrade: Applying masscut ', masscut
	
	! TBD...
	nlines_degrade =0
	open(unit=1,file=inputfilename,action='read')
	open(unit=2,file=outputfilename1)
	open(unit=3,file=outputfilename2)
	do while(.true.)
		read(1,'(A)',end=102) tmpstr1
		read(tmpstr1,*) tmp(1:maxcol)
		x=tmp(xcol); y=tmp(ycol); z=tmp(zcol); mass=tmp(masscol)
		if(mass > masscut) then
			write(2,'(A)') trim(adjustl(tmpstr1))
			nlines_degrade = nlines_degrade+1
		endif
		cycle
102		exit
	enddo
	close(2);

	write(*,'(A,A,A,i10,A,i10,A,A)') '# This is for file ', trim(adjustl(outputfilename)), '. This file has ',nlines,&
		' lines. We want to degrade it to ',numdegrade,' halos. New sample written to: ', trim(adjustl(outputfilename1))
	write(3,'(A,A,A,i10,A,i10,A,A)') '# This is for file ', trim(adjustl(outputfilename)), '. This file has ',nlines,&
		' lines. We want to degrade it to ',numdegrade,' halos. New sample written to: ', trim(adjustl(outputfilename1))

	write(*,'(A,e15.7,A,i10,A)') '# When selecting halos with mass>masscut=', mass1, ', there are ', num1, ' halos'
	write(3,'(A,e15.7,A,i10,A)') '# When selecting halos with mass>masscut=', mass1, ', there are ', num1, ' halos'

	write(*,'(A,e15.7,A,i10,A)') '# When selecting halos with mass>masscut=', mass2, ', there are ', num2, ' halos'
	write(3,'(A,e15.7,A,i10,A)') '# When selecting halos with mass>masscut=', mass2, ', there are ', num2, ' halos'

	write(*,'(A,e15.7,A,i10,A,f15.7)') '# We interploate and adopt  mass>masscut=', masscut, ', we obtain ', nlines_degrade, &
		 ' halos. Rato (num-obtaned/num-goal): ', dble(nlines_degrade)/dble(numdegrade)
	write(3,'(A,e15.7,A,i10,A,f15.7)') '# We interploate and adopt  mass>masscut=', masscut, ', we obtain ', nlines_degrade, &
		 ' halos. Rato (num-obtaned/num-goal): ', dble(nlines_degrade)/dble(numdegrade)
	close(3)
end program main
