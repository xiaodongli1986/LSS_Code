
program main

use LSS_cosmo_funs

implicit none

	integer :: xcol,ycol,zcol,wcol, maxcol, skiprow,i
	integer, parameter :: bignum = 50000000
	real(dl) :: omin,win,omout,wout, x,y,z,weight,r,redshift, randreds(2,bignum),rat, tmp(1000),&
		xmin,xmax,ymin,ymax,zmin,zmax,rmin,rmax,redmin,redmax
	logical :: hasweight
	integer(8) :: numarg, nlines
	character(len=char_len) :: suffix, inputfilename, outputfilename, outputinfo, tmpstr1,tmpstr2, printstr
	
	printstr = 'Usage: EXE -input intpufilename -omin -omin -win -win '//&
		'-omout -omout -wout -wout -xcol xcol -ycol ycol -zcol zcol '//&
		'-hasweight hasweight -wcol wcol -skiprow skiprow -suffix suffix  '//&
		'### fmt of output: x, y, z, weight ### Converting the input file into a different cosmology'

	! Default values
	omin = 0.26; win = -1.0; 
	omout = 0.30; wout = -1.0;
	hasweight = .false.
	xcol=1; ycol=2; zcol=3; wcol=4; 
	skiprow=0;
	suffix = '.cosmo-converted'
	
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
		elseif(trim(adjustl(tmpstr1)).eq.'-omin') then
			read(tmpstr2,*) omin
		elseif(trim(adjustl(tmpstr1)).eq.'-win') then
			read(tmpstr2,*) win
		elseif(trim(adjustl(tmpstr1)).eq.'-omout') then
			read(tmpstr2,*) omout
		elseif(trim(adjustl(tmpstr1)).eq.'-wout') then
			read(tmpstr2,*) wout
		elseif(trim(adjustl(tmpstr1)).eq.'-xcol') then
			read(tmpstr2,*) xcol
		elseif(trim(adjustl(tmpstr1)).eq.'-ycol') then
			read(tmpstr2,*) ycol
		elseif(trim(adjustl(tmpstr1)).eq.'-zcol') then
			read(tmpstr2,*) zcol
		elseif(trim(adjustl(tmpstr1)).eq.'-wcol') then
			read(tmpstr2,*) wcol
		elseif(trim(adjustl(tmpstr1)).eq.'-skiprow') then
			read(tmpstr2,*) skiprow
		elseif(trim(adjustl(tmpstr1)).eq.'-hasweight') then
			read(tmpstr2,*) hasweight
		elseif(trim(adjustl(tmpstr1)).eq.'-suffix') then
			read(tmpstr2,'(A)') suffix
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	if(.not.hasweight) then
		maxcol = max(xcol,ycol,zcol)
	else
		maxcol = max(xcol,ycol,zcol,wcol)
		if(wcol.eq.xcol .or. wcol.eq.ycol .or. wcol.eq.zcol) then
			print *, 'WARNING!!! Found same wcol and xyzcol; Check carefully!!!', xcol, ycol, zcol, wcol
		endif
	endif
	
	if(maxcol>size(tmp)) then
		print *, 'Overflow: increase size of tmp!'
		print *, '  maxcol, size(tmp) = ', maxcol, size(tmp)
		stop
	endif

	! intput cosmology
	call cosmo_funs_init(.true.)
	gb_omegam = omin; gb_w = win; gb_h = 0.73_dl;
	call de_calc_comovr()
	
	outputfilename = trim(adjustl(inputfilename))//trim(adjustl(suffix))//'.'//trim(adjustl(get_omwstr(omout,wout)))
	outputinfo = trim(adjustl(outputfilename))//'.info'

	open(unit=100,file=outputinfo)

	print *, '#####################################'
	write(*,'(A)')   	'  Settings:'
	write(*,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(*,'(A,2f16.7)') 	'   omin, win = ', real(omin), real(win)
	write(*,'(A,2f16.7)') 	'   omout, wout = ', real(omout), real(wout)
	write(*,'(A,L3)')	'   hasweight = ', hasweight
	write(*,'(A,4i3)')  	'   cols of x,y,z,weight: ', xcol,ycol,zcol,wcol
	write(*,'(A,i5)')	'   skip rows ', skiprow
	write(*,'(A,A)')	'   suffix = ', trim(adjustl(suffix))
	write(*,'(A,i3)')	'   maxcol = ', maxcol

	write(100,'(A)')   	'  Settings:'
	write(100,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(100,'(A,2f16.7)') 	'   omin, win = ', real(omin), real(win)
	write(100,'(A,2f16.7)') 	'   omout, wout = ', real(omout), real(wout)
	write(100,'(A,L3)')	'   hasweight = ', hasweight
	write(100,'(A,4i3)')  	'   cols of x,y,z,weight: ', xcol,ycol,zcol,wcol
	write(100,'(A,i5)')	'   skip rows ', skiprow
	write(100,'(A,A)')	'   suffix = ', trim(adjustl(suffix))
	write(100,'(A,i3)')	'   maxcol = ', maxcol
	print *, '#####################################'

	nlines =0; 
	xmin = 1.0e30; ymin = 1.0e30; zmin = 1.0e30; rmin = 1.0e30; redmin = 1.0e30
	xmax = -xmin; ymax = -ymin; zmax = -zmin; rmax = -rmin; redmax = -redmin
	
	open(unit=1,file=inputfilename,action='read')
	do i = 1, skiprow
		read(1,*) tmpstr1
	enddo
	do while(.true.)
		read(1,*,end=101) tmp(1:maxcol)
		nlines = nlines+1
		if(nlines>bignum) then
			print *, 'ERROR of overflow: increase bignum!: ', nlines, bignum
			stop
		endif
		x=tmp(xcol); y=tmp(ycol); z=tmp(zcol); r = sqrt(x*x+y*y+z*z)
		xmin = min(x,xmin); xmax = max(x,xmax)
		ymin = min(y,ymin); ymax = max(y,ymax)
		zmin = min(z,zmin); zmax = max(z,zmax)
		rmin = min(r,rmin); rmax = max(r,rmax)
	        redshift = de_zfromintpl(r)
		redmin = min(redshift,redmin); redmax = max(redshift,redmax)
	        randreds(1,nlines) = r
	        randreds(2,nlines) = redshift
		cycle
101		exit
	enddo
	close(1)
	write(*,'(A,i10,A)') ' Finishing processing ', nlines, ' lines.'
	write(100,'(A,i10,A)') ' Finishing processing ', nlines, ' lines.'

	write(*,'(A,5("(",2f15.7,"); "))') 'Range of x,y,z,r,red: ', xmin,xmax,ymin,ymax,zmin,zmax, rmin,rmax, redmin,redmax
	write(100,'(A,5("(",2f15.7,"); "))') 'Range of x,y,z,r,red: ', xmin,xmax,ymin,ymax,zmin,zmax, rmin,rmax, redmin,redmax

	close(100)
	
	print *, 'Outputing information to file ', trim(adjustl(outputfilename))
	
	! output cosmology
	call cosmo_funs_init(.true.)
	gb_omegam = omout; gb_w = wout; gb_h = 0.73_dl;
	call de_calc_comovr()

	! output to file
	open(unit=1,file=inputfilename,action='read')
	open(unit=2,file=outputfilename)
	do i = 1, skiprow
		read(1,*) tmpstr1
	enddo
	i =0
	do while(.true.)
		read(1,*,end=102) tmp(1:maxcol)
		i = i+1
		if(i>nlines) then
			print *, 'ERROR of #-line overflow! ', i, nlines
			stop
		endif
		x=tmp(xcol); y=tmp(ycol); z=tmp(zcol);
		if(hasweight) then
			weight=tmp(wcol)
		endif
		r = randreds(1,i)
		rat = de_get_comovr(randreds(2,i))/r
		if(.not.hasweight) then
			write(2,'(3(1x,f13.6))')  x*rat, y*rat, z*rat
		else
			write(2,'(3(1x,f13.6),1x,e15.7)')  x*rat, y*rat, z*rat, weight
		endif
		cycle
102		exit
	enddo
	close(1)
	close(2)
	
	if(i.ne.nlines) then
		print *, 'ERROR of #-line mismatches! ', i, nlines
		stop
	else
		print *, 'Expected #-line matches: ', i, nlines
	endif

end program main
