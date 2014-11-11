
program main

use LSS_cosmo_funs

implicit none

	real(dl) :: omegam,w, tmp(1000), x,y,z,vx,vy,vz,mass, r,robs,rat, vlos,redshift,redobs, &
		xmin,xmax,ymin,ymax,zmin,zmax,rmin,rmax
	integer :: xcol,ycol,zcol,vxcol,vycol,vzcol,masscol,maxcol,skiprow, i
	integer(8) :: numarg, nlines
	character(len=char_len) :: printstr, suffix, inputfilename, outputfilename, outputfilenameinfo, tmpstr1,tmpstr2
	
	printstr = 'Usage: EXE -input intpufilename -omegam omegam -w w -xcol xcol -ycol ycol -zcol zcol -vxcol vxcol -vycol vycol -vzcol vzcol -masscol masscol -skiprow skiprow -suffix suffix  ### fmt of output: x, y, z, vx, vy, vz, mass, redshift, redobs = redshift + vlos * (1.0+redshift) / const_c, rat = robs / r  ### Generating "compact" sample, matching the format of creat-mock'

	! Default values
	omegam = 0.26; w = -1.0; 
	xcol=1; ycol=2; zcol=3; vxcol=4; vycol=5; vzcol=6; masscol=7; skiprow=0
	suffix = '.compact'
	
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
		elseif(trim(adjustl(tmpstr1)).eq.'-omegam') then
			read(tmpstr2,*) omegam
		elseif(trim(adjustl(tmpstr1)).eq.'-w') then
			read(tmpstr2,*) w
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
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	maxcol = max(xcol,ycol,zcol,vxcol,vycol,vzcol,masscol)
	if(maxcol>size(tmp)) then
		print *, 'Overflow: increase size of tmp!'
		print *, '  maxcol, size(tmp) = ', maxcol, size(tmp)
		stop
	endif

	call cosmo_funs_init(.true.)
	gb_omegam = omegam; gb_w = w; gb_h = 0.73_dl;
	call de_calc_comovr()
	
	outputfilename = trim(adjustl(inputfilename))//trim(adjustl(suffix))
	outputfilenameinfo = trim(adjustl(outputfilename))//'.info'

	open(unit=100,file=outputfilenameinfo)

	print *, '#####################################'
	write(*,'(A)')   	'  Settings:'
	write(*,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(*,'(A,2f16.7)') 	'   omegam, w = ', real(omegam), real(w)
	write(*,'(A,7i3)')  	'   cols of x,y,z,vx,vy,vz,mass: ', xcol,ycol,zcol,vxcol,vycol,vzcol,masscol
	write(*,'(A,i5)')	'   skip rows ', skiprow
	write(*,'(A,A)')	'   suffix = ', trim(adjustl(suffix))

	write(100,'(A)')   	'  Settings:'
	write(100,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(100,'(A,2f16.7)') '   omegam, w = ', real(omegam), real(w)
	write(100,'(A,7i3)')  	'   cols of x,y,z,vx,vy,vz,mass: ', xcol,ycol,zcol,vxcol,vycol,vzcol,masscol
	write(100,'(A,i5)')	'   skip rows ', skiprow
	write(100,'(A,A)')	'   suffix = ', trim(adjustl(suffix))
	print *, '#####################################'

	nlines =0; 
	xmin = 1.0e30; ymin = 1.0e30; zmin = 1.0e30; rmin = 1.0e30
	xmax = -xmin; ymax = -ymin; zmax = -zmin; rmax = -rmin
	open(unit=1,file=inputfilename,action='read')
	open(unit=2,file=outputfilename)
	do i = 1, skiprow
		read(1,*) tmpstr1
	enddo
	do while(.true.)
		read(1,*,end=101) tmp(1:maxcol)
		nlines = nlines+1
		x=tmp(xcol); y=tmp(ycol); z=tmp(zcol); r = sqrt(x*x+y*y+z*z)
		vx=tmp(vxcol); vy=tmp(vycol); vz=tmp(vzcol); mass = tmp(masscol)
		xmin = min(x,xmin); xmax = max(x,xmax)
		ymin = min(y,ymin); ymax = max(y,ymax)
		zmin = min(z,zmin); zmax = max(z,zmax)
		rmin = min(r,rmin); rmax = max(r,rmax)
	        redshift = de_zfromintpl(r)
	        vlos = (vx*x + vy*y + vz*z) / r
	        redobs = redshift + vlos * (1.0+redshift) / const_c
	        robs = de_get_comovr(redobs)
	        rat = robs / r
		write(2,'(6(1x,f13.6), e14.7, 3(1x,f10.7) )')  x, y, z, vx, &
		      vy, vz, mass, redshift, redobs, rat
		cycle
101		exit
	enddo
	
	close(2)
	write(*,'(A,i10,A)') ' Finishing processing ', nlines, ' lines.'
	write(100,'(A,i10,A)') ' Finishing processing ', nlines, ' lines.'

	write(*,'(A,4("(",2f15.7,"); "))') 'Range of x,y,z,r: ', xmin,xmax,ymin,ymax,zmin,zmax, rmin,rmax
	write(100,'(A,4("(",2f15.7,"); "))') 'Range of x,y,z,r: ', xmin,xmax,ymin,ymax,zmin,zmax, rmin,rmax

	close(100)

end program main
