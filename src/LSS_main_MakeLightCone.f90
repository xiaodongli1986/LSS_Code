
!! Tomorrow Task
!!	1. Check this code; (everything right?)
!!	2. Running this code on baekdu

program main

use LSS_settings_init

implicit none

	real(dl) :: x,y,z,r, xyzmin,xyzmax, omegam,w, tmp(1000),constantmasscut, dxyz,maxdist, nowdiff,mindiff
	real(dl), allocatable :: obs_pos(:,:), slicereds(:), slicedists(:), masscuts(:), masscutreds(:), &
		tmpinfos(:,:), tmpdists(:), tmpxyzs(:,:)
	integer :: numobs, numslice, nummasscut, numskiplowmass, numshift, xcol,ycol,zcol,vxcol,vycol,vzcol,masscol,skiprow,maxcol,numarg, &
		i,j,ndat, ishiftx,ishifty,ishiftz,islice, numLCdat
	integer, allocatable :: outputfileunits(:)
	character(len=char_len) :: inputfilename, outputfilename, outputfilenameinfo, sliceinfofile, obsinfofile, evolvmasscutfile, &
		suffix, tmpstr1, tmpstr2, printstr1,printstr2,printstr3,printstr4,printstr5,printstr6,printstr7,printstr8
	character(len=char_len), allocatable :: outputfiles(:)
	logical :: doevolvmasscut, rsdshiftinclu, endoffile, withinLC
!	const_Mpctokm
	printstr1 = 'Usage: EXE -input inputfilename -omegam omegam -w w -xyzmin xyzmin -xyzmax xyzmax -maxdist maxdist -xcol xcol -ycol ycol -zcol zcol -vxcol vxcol -vycol vycol -vzcol vzcol -masscol masscol -skiprow skiprow -suffix suffix -sliceinfo sliceinfofile -obsinfo obsinfofile -constantmasscut constantmasscut -doevolvmasscut doevolvmasscut  -evolvmasscut evolvmasscutfile -rsdshiftinclu rsdshiftinclu -numshift numshift'
	printstr2 = '### # sliceinfo is information of slice; fmt: redshift'
	printstr3 = '# obsinfo is information of observer; fmt: xloc,yloc,zloc (in percentage, e.g. 0.5 0.5 0.5 means an observer at the center)'
	printstr4 = '# maxdist is the maximal distance that the LC reaches. ' 
	printstr5 = '# A positive mass cut will be applied to avoid zero mass halo (so far doevolvmasscut not supported...)'
	printstr6 = '# numshift: shift the periodical box to reduplicate the data. if numshift = 1, will making (2*1+1)^3 = 9 copies'
	printstr7 = '### Example: ./MakeLightCone -input testdats/HR4LC/SubSample00.dat -xyzmin 0.0 -xyzmax 3150.0 -obsinfo HR4LC/obsinfo.txt -sliceinfo HR4LC/sliceinfo.txt -omegam 0.26 -w -1.0 -maxdist 1800.0 -numshift 1'
! E.g.  obspos.txt: 
!	 0.00 0.00 0.00
!	 0.00 0.00 0.50
!	 0.00 0.50 0.00
!	 0.50 0.00 0.00
!	 0.50 0.00 0.50
!	 0.50 0.50 0.00
!	 0.00 0.50 0.50
!	 0.50 0.50 0.50

!	'If doevolvmasscut, applying evolving masscut depending on redshift. fmt of evolvmasscutfile: redshift masscut. if rsdshiftinclu is .true., then also considering rsdshift when applying evolving masscut (rsd shift r to r_shift, then store this halo if either mass > mass_cut(r) .or. mass > mass_cut(r_shift))'

	! Default values
	xyzmin = 0.0; xyzmax = 0.0; maxdist=0.0 ! must set them!!!
	omegam = 0.26; w = -1.0; om_dft = omegam; w_dft = w;
	xcol=1; ycol=2; zcol=3; vxcol=4; vycol=5; vzcol=6; masscol=7 
	suffix = '.LC'
	numshift = 1
	skiprow=0
	doevolvmasscut = .false.; rsdshiftinclu = .false.; 
	constantmasscut = 1.0d0 ! avoid zero mass !

	numarg = iargc()
	if(numarg.le.1) then
		write(*,'(5(/A/))') trim(adjustl(printstr1)),trim(adjustl(printstr2)),trim(adjustl(printstr3)),&
			trim(adjustl(printstr4)),trim(adjustl(printstr5)),trim(adjustl(printstr6)), trim(adjustl(printstr7))
		stop
	endif
	
	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfilename
		elseif(trim(adjustl(tmpstr1)).eq.'-omegam') then
			read(tmpstr2,*) omegam
		elseif(trim(adjustl(tmpstr1)).eq.'-w') then
			read(tmpstr2,*) w
		elseif(trim(adjustl(tmpstr1)).eq.'-xyzmin') then
			read(tmpstr2,*) xyzmin
		elseif(trim(adjustl(tmpstr1)).eq.'-xyzmax') then
			read(tmpstr2,*) xyzmax
		elseif(trim(adjustl(tmpstr1)).eq.'-maxdist') then
			read(tmpstr2,*) maxdist
		elseif(trim(adjustl(tmpstr1)).eq.'-xcol') then
			read(tmpstr2,*) xcol
		elseif(trim(adjustl(tmpstr1)).eq.'-ycol') then
			read(tmpstr2,*) ycol
		elseif(trim(adjustl(tmpstr1)).eq.'-zcol') then
			read(tmpstr2,*) zcol
		elseif(trim(adjustl(tmpstr1)).eq.'-vxcol') then
			read(tmpstr2,*) vxcol
		elseif(trim(adjustl(tmpstr1)).eq.'-vycol') then
			read(tmpstr2,*) vycol
		elseif(trim(adjustl(tmpstr1)).eq.'-vzcol') then
			read(tmpstr2,*) vzcol
		elseif(trim(adjustl(tmpstr1)).eq.'-masscol') then
			read(tmpstr2,*) masscol
		elseif(trim(adjustl(tmpstr1)).eq.'-skiprow') then
			read(tmpstr2,*) skiprow
		elseif(trim(adjustl(tmpstr1)).eq.'-suffix') then
			read(tmpstr2,'(A)') suffix
		elseif(trim(adjustl(tmpstr1)).eq.'-numshift') then
			read(tmpstr2,*) numshift
		elseif(trim(adjustl(tmpstr1)).eq.'-sliceinfo') then
			read(tmpstr2,'(A)') sliceinfofile
		elseif(trim(adjustl(tmpstr1)).eq.'-obsinfo') then
			read(tmpstr2,'(A)') obsinfofile
		elseif(trim(adjustl(tmpstr1)).eq.'-constantmasscut') then
			read(tmpstr2,*) constantmasscut
		elseif(trim(adjustl(tmpstr1)).eq.'-doevolvmasscut') then
			read(tmpstr2,*) doevolvmasscut
		elseif(trim(adjustl(tmpstr1)).eq.'-evolvmasscut') then
			read(tmpstr2,'(A)') evolvmasscutfile
		elseif(trim(adjustl(tmpstr1)).eq.'-rsdshiftinclu') then
			read(tmpstr2,*) rsdshiftinclu
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
		write(*,'(5(/A/))') trim(adjustl(printstr1)),trim(adjustl(printstr2)),trim(adjustl(printstr3)),&
			trim(adjustl(printstr4)),trim(adjustl(printstr5)),trim(adjustl(printstr6)), trim(adjustl(printstr7))
			stop
		endif
	enddo

	maxcol = max(xcol,ycol,zcol,vxcol,vycol,vzcol,masscol)
	if(maxcol>size(tmp)) then
		print *, 'Overflow: increase size of tmp!'
		print *, '  maxcol, size(tmp) = ', maxcol, size(tmp)
		stop
	endif
	
	dxyz=xyzmax-xyzmin
	
	call cosmo_funs_init(.true.)
	gb_omegam = omegam; gb_w = w; gb_h = 0.73_dl;
	call de_calc_comovr()
!	do i = 1, 10
!		x = exp((i+0.0)/5.0) + 0.0_dl
!		print *, x, t_age(x) !* const_Mpctokm / const_gyrtos
!	enddo

	! infor of observer
	call read_in_revfmt(obsinfofile, 3, numobs, obs_pos)
	allocate(outputfiles(numobs),outputfileunits(numobs))
	do i = 1, numobs
		write(tmpstr1,*) i-1
		if(i.le.9) then
			outputfiles(i) = trim(adjustl(inputfilename))//trim(adjustl(suffix))//'0'//trim(adjustl(tmpstr1))
		else
			outputfiles(i) = trim(adjustl(inputfilename))//trim(adjustl(suffix))//trim(adjustl(tmpstr1))
		endif
		outputfileunits(i) = 9870+i
		open(unit=outputfileunits(i),file=outputfiles(i))
	enddo
	
	! info of slice...
	call count_line_number(sliceinfofile, numslice)
	allocate(slicereds(numslice),slicedists(numslice))
	open(unit=1987,file=sliceinfofile)
	do i = 1, numslice
		read(1987,*) slicereds(i)
		slicedists(i) = de_get_comovr(slicereds(i))
	enddo
	close(1987)
	
	! information
	outputfilenameinfo = trim(adjustl(inputfilename))//trim(adjustl(suffix))//'.info'
	open(unit=100,file=outputfilenameinfo)
		
	print *, '#####################################'
	write(*,'(A)')   	'  Settings:'
	write(*,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(*,'(A,2f16.7)')   '   bounary of box = ', xyzmin, xyzmax
	write(*,'(A,2f16.7)') 	'   omegam, w = ', real(omegam), real(w)
	write(*,'(A,7i3)')  	'   cols of x,y,z,vx,vy,vz,mass: ', xcol,ycol,zcol,vxcol,vycol,vzcol,masscol
	write(*,'(A,i5)')	'   skip rows ', skiprow
	write(*,'(2(A,2x))')	'   suffix = ', trim(adjustl(suffix))
	write(*,'(A,i3,i5)')	'   numshift = ', numshift, (2*numshift+1)**3
	write(*,'(3(A,2x))')	'   obsinfo, sliceinfo = ', trim(adjustl(obsinfofile)), trim(adjustl(sliceinfofile))
	write(*,'(A,i3,A)')	'   Positions of ',numobs, ' observers (rat, loc):'
	do i = 1, numobs
		write(*,'(8x,i3,3(f10.5),5x,3(f10.3,3x),A,A)') i,obs_pos(1:3,i), xyzmin+dxyz*obs_pos(1:3,i), &
			'   file: ', trim(adjustl(outputfiles(i)))
	enddo
	write(*,'(A,i3,A)') '   Redshift/Distances of ', numslice, ' slices:'
	do i = 1, numslice
		write(*,'(12x,f12.7,3x,f12.3)') slicereds(i), slicedists(i)
	enddo
	write(*,'(A,e15.7)') 	'   applying constant mass cut = ', constantmasscut
	if(doevolvmasscut) then
		write(*,'(A,A,L2,A)')	'   evolvmasscutfile, rsdshiftinclu = ', trim(adjustl(evolvmasscutfile)), &
			rsdshiftinclu, '. Info (reds, masscuts):'
		do i = 1, nummasscut
			write(*,'(8x,f12.7,e15.7)') masscutreds(i), masscuts(i)
		enddo
	endif

	write(100,'(A)')   	'  Settings:'
	write(100,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(100,'(A,2f16.7)')   '   bounary of box = ', xyzmin, xyzmax
	write(100,'(A,2f16.7)') 	'   omegam, w = ', real(omegam), real(w)
	write(100,'(A,7i3)')  	'   cols of x,y,z,vx,vy,vz,mass: ', xcol,ycol,zcol,vxcol,vycol,vzcol,masscol
	write(100,'(A,i5)')	'   skip rows ', skiprow
	write(100,'(2(A,2x))')	'   suffix = ', trim(adjustl(suffix))
	write(100,'(A,i3,i5)')	'   numshift = ', numshift, (2*numshift+1)**3
	write(100,'(3(A,2x))')	'   obsinfo, sliceinfo = ', trim(adjustl(obsinfofile)), trim(adjustl(sliceinfofile))
	write(100,'(A,i3,A)')	'   Positions of ',numobs, ' observers (rat, loc):'
	do i = 1, numobs
		write(100,'(8x,i3,3(f10.5),5x,3(f10.3,3x),A,A)') i,obs_pos(1:3,i), xyzmin+dxyz*obs_pos(1:3,i), &
			'   file: ', trim(adjustl(outputfiles(i)))
		obs_pos(1:3,i) = xyzmin+dxyz*obs_pos(1:3,i)
		print *, obs_pos(1:3,i)
	enddo
	write(100,'(A,i3,A)') '   Redshift/Distances of ', numslice, ' slices:'
	do i = 1, numslice
		write(100,'(12x,f12.7,3x,f12.3)') slicereds(i), slicedists(i)
	enddo
	write(100,'(A,e15.7)') 	'   applying constant mass cut = ', constantmasscut
	if(doevolvmasscut.or.rsdshiftinclu) then
		print *, 'ERROR (MakeLightCone)!!! doevolvmasscut, rsdshiftinclu not supported yet!!!', &
			doevolvmasscut, rsdshiftinclu
		stop
		write(100,'(A,A,L2,A)')	'   evolvmasscutfile, rsdshiftinclu = ', trim(adjustl(evolvmasscutfile)), &
			rsdshiftinclu, '. Info (reds, masscuts):'
		do i = 1, nummasscut
			write(100,'(8x,f12.7,e15.7)') masscutreds(i), masscuts(i)
		enddo
	endif
	print *, '#####################################'


	! main loop
	allocate(tmpinfos(maxcol,numslice),tmpdists(numslice),tmpxyzs(3,numslice))
	ndat = 0
	numLCdat = 0
	numskiplowmass = 0
	endoffile = .false.
	open(unit=1,file=inputfilename,action='read')
	do while(.true.)
		do j = 1, numslice
			read(1,*,end=100) tmpinfos(1:maxcol,j)
			cycle
100			endoffile = .true.
			exit			
		enddo
		
		if(endoffile) then
			if(j.ne.1) then
				print *, 'Warning! Unexpected break with j.ne.1: ', j
				write(100,*) 'Warning! Unexpected break with j.ne.1: ', j
			else
				write(100,*) 'End of file. As expected j = ', j
			endif
			exit
		endif
		
!		if(ndat.ge.10000) exit
		
		ndat = ndat+1
		do i = 1, numobs
			do ishiftx = -numshift,numshift
			do ishifty = -numshift,numshift
			do ishiftz = -numshift,numshift
				withinLC = .false.
				do j = 1, numslice
					tmpxyzs(1,j) = tmpinfos(xcol,j) + ishiftx * dxyz
					tmpxyzs(2,j) = tmpinfos(ycol,j) + ishifty * dxyz
					tmpxyzs(3,j) = tmpinfos(zcol,j) + ishiftz * dxyz
					tmpdists(j) = sqrt((tmpxyzs(1,j)-obs_pos(1,i))**2.0 + &
							(tmpxyzs(2,j)-obs_pos(2,i))**2.0 + &
							(tmpxyzs(3,j)-obs_pos(3,i))**2.0)
					if(tmpdists(j) < maxdist) withinLC = .true. ! if any of them close enough to the center...
				enddo
				if(.not.withinLC) cycle
				islice = 1
				mindiff = abs(tmpdists(1) - slicedists(1))
				do j = 2, numslice
					nowdiff = abs(tmpdists(j) - slicedists(j))
					if(nowdiff < mindiff) then
						islice = j
						mindiff = nowdiff
					endif
				enddo
				if(tmpinfos(masscol,islice) .ge. constantmasscut) then
					write(outputfileunits(i),'(7e15.7,i10,i3,5x,3i3,5x,2e15.7,3x,f10.3)') &
						tmpxyzs(1:3,islice)-obs_pos(1:3,i), &
						tmpinfos(vxcol,islice),tmpinfos(vycol,islice),tmpinfos(vzcol,islice),&
						tmpinfos(masscol,islice), ndat, islice, ishiftx,ishifty,ishiftz, &
						tmpdists(islice), slicedists(islice), mindiff
					numLCdat = numLCdat+1
				else
					numskiplowmass = numskiplowmass+1
				endif
			enddo
			enddo
			enddo
		enddo
	enddo		
	
	write(*,'(A,3i10)') ' # of total data / light cone data / skip low mass:  ', ndat, numLCdat, numskiplowmass
	write(100,'(A,3i10)') ' # of total data / light cone data / skip low mass:  ', ndat, numLCdat, numskiplowmass
	
	close(100)
	do i = 1, numobs
		close(unit=outputfileunits(i))
	enddo

end program main

