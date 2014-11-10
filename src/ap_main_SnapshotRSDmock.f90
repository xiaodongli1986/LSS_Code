! Aug 29: We consider rescale of volume due to change of the cosmology!!! (dovoleff)


program main

use ap_cosmo_funs

implicit none

	real(dl) :: xyzmin,xyzmax,dxyz, redshift,omegam,w,dist, x,y,z,vx,vy,vz, xmin,xmax,ymin,ymax,zmin,zmax, tmp(1000), &
		zobs, shiftdist, shiftedcord, volscale, voldft, volnow, rescalefac
	integer :: xcol,ycol,zcol,vxcol,vycol,vzcol,maxcol, i,j, numshiftback, numshiftbackorig, skiprow
	integer(8) :: numarg, nlines
	character(len=char_len) :: printstr, inputfilename, outputname, outputfilenamex,outputfilenamey,outputfilenamez,&
		outputfilenameinfo, tmpstr1,tmpstr2, suffix
	logical :: shiftx,shifty,shiftz, dovoleff
	
	printstr = 'Usage: EXE -input intpufilename -xyzmin xyzmin -xyzmax xyzmax -xcol xcol -ycol ycol -zcol zcol -vxcol vxcol -vycol vycol -vzcol vzcol -redshift redshift -omegam omegam -w w -shiftx shiftx -shifty shifty -shiftz shiftz -skiprow skiprow -suffix suffix -dovoleff dovoleff ### xyzmin/xyzmax used to do periodical boundary condition (shift inside to box if shifted outside by RSD); please tell cosmology of the simulation, redshift of the snapshot; please choose which direction to shift (x,y,z; three possibilities) ### dovoleff only respect to Om=0.26/w=-1.0 cosmology!!!'

	! Default values
	xyzmin = -1.0e30; xyzmax = 1.0e30;
	xcol=1; ycol=2; zcol=3; vxcol=4; vycol=5; vzcol=6; skiprow=0
	shiftx = .true.; shifty=.false.; shiftz=.false.
	omegam = 0.26; w = -1.0; redshift = 0.0
	suffix = ''
	dovoleff = .false.
	
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
		elseif(trim(adjustl(tmpstr1)).eq.'-xyzmin') then
			read(tmpstr2,*) xyzmin
		elseif(trim(adjustl(tmpstr1)).eq.'-xyzmax') then
			read(tmpstr2,*) xyzmax
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
		elseif(trim(adjustl(tmpstr1)).eq.'-redshift') then
			read(tmpstr2,*) redshift
		elseif(trim(adjustl(tmpstr1)).eq.'-omegam') then
			read(tmpstr2,*) omegam
		elseif(trim(adjustl(tmpstr1)).eq.'-w') then
			read(tmpstr2,*) w
		elseif(trim(adjustl(tmpstr1)).eq.'-shiftx') then
			read(tmpstr2,*) shiftx
		elseif(trim(adjustl(tmpstr1)).eq.'-shifty') then
			read(tmpstr2,*) shifty
		elseif(trim(adjustl(tmpstr1)).eq.'-shiftz') then
			read(tmpstr2,*) shiftz
		elseif(trim(adjustl(tmpstr1)).eq.'-skiprow') then
			read(tmpstr2,*) skiprow
		elseif(trim(adjustl(tmpstr1)).eq.'-suffix') then
			read(tmpstr2,*) suffix
		elseif(trim(adjustl(tmpstr1)).eq.'-dovoleff') then
			read(tmpstr2,*) dovoleff
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	if(.not.shiftx.and..not.shifty.and..not.shiftz) then
		print *, 'No action! All shifts are .false.: ', shiftx,shifty,shiftz
		stop
	endif

	maxcol = max(xcol,ycol,zcol,vxcol,vycol,vzcol)
	if(maxcol>size(tmp)) then
		print *, 'Overflow: increase size of tmp!'
		print *, '  maxcol, size(tmp) = ', maxcol, size(tmp)
		stop
	endif

	dxyz = xyzmax-xyzmin
	
	gb_h = 0.73_dl;
	if(dovoleff) then
		gb_omegam = 0.26; gb_w = -1.0; 
		call de_calc_comovr()
		dist = comov_r(dble(redshift))
		voldft = dist**2.0/Hz(dble(redshift))
	endif

	call cosmo_funs_init(.true.)
	gb_omegam = omegam; gb_w = w
	call de_calc_comovr()
	dist = comov_r(dble(redshift))
	if(dovoleff) then
		volnow = dist**2.0/Hz(dble(redshift))
		volscale = volnow / voldft
		rescalefac = volscale**(1.0d0/3.0d0)
		write(*,'(4x,A,2e14.7,2f10.5)') &
			'Doing vol rescale w.r.t. Om=0.26/-1.0 cosmology: voldft, volnow, volscale, rescalefac = ', &
			voldft, volnow, volscale, rescalefac
	else
		rescalefac = 1.0
	endif
	
	if(trim(adjustl(suffix)).ne.'') then
		outputname = trim(adjustl(inputfilename))//trim(adjustl(suffix))
	else
		outputname = inputfilename
	endif
	
	outputfilenamex = trim(adjustl(outputname))//'.shiftx'
	outputfilenamey = trim(adjustl(outputname))//'.shifty'
	outputfilenamez = trim(adjustl(outputname))//'.shiftz'
	outputfilenameinfo = trim(adjustl(outputname))//'.shift.info'

	open(unit=100,file=outputfilenameinfo)

	print *, '#####################################'
	write(*,'(A)')   	'  Settings:'
	write(*,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(*,'(A,2f16.7)')   '   bounary of box = ', xyzmin, xyzmax
	write(*,'(A,3L2)') 	'   RSD shift at x,y,z = ', shiftx, shifty, shiftz
	write(*,'(A,4f16.7)') 	'   omegam, w, redshift, comov_r(redshift) = ', real(omegam), real(w), real(redshift), real(dist)
	write(*,'(A,6i3)')  	'   cols of x,y,z,vx,vy,vz: ', xcol,ycol,zcol,vxcol,vycol,vzcol
	write(*,'(A,i5)')	'   skip rows ', skiprow
	write(*,'(A,f10.5)')	'   rescalefac ', rescalefac

	write(100,'(A)')   	'  Settings:'
	write(100,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(100,'(A,2f16.7)')	'   bounary of box = ', xyzmin, xyzmax
	write(100,'(A,3L2)') 	'   RSD shift at x,y,z = ', shiftx, shifty, shiftz
	write(100,'(A,4f16.7)') '   omegam, w, redshift, comov_r(redshift) = ', real(omegam), real(w), real(redshift), real(dist)
	write(100,'(A,6i3)')  	'   cols of x,y,z,vx,vy,vz: ', xcol,ycol,zcol,vxcol,vycol,vzcol
	write(100,'(A,f10.5)')	'   rescalefac ', rescalefac
	print *, '#####################################'

	if(shiftx) open(unit=10,file=outputfilenamex)
	if(shifty) open(unit=20,file=outputfilenamey)
	if(shiftz) open(unit=30,file=outputfilenamez)

	nlines =0; numshiftback=0; numshiftbackorig=0
	xmin = 1.0e30; ymin = 1.0e30; zmin = 1.0e30
	xmax = -xmin; ymax = -ymin; zmax = -zmin
	open(unit=1,file=inputfilename,action='read')
	do i = 1, skiprow
		read(1,*) tmpstr1
	enddo
	do while(.true.)
		read(1,*,end=101) tmp(1:maxcol)
		nlines = nlines+1
		x=tmp(xcol); y=tmp(ycol); z=tmp(zcol); 
		vx=tmp(vxcol); vy=tmp(vycol); vz=tmp(vzcol); 
		xmin = min(x,xmin); xmax = max(x,xmax)
		ymin = min(y,ymin); ymax = max(y,ymax)
		zmin = min(z,zmin); zmax = max(z,zmax)
		if(x < xyzmin) then
			x = x + dxyz
			numshiftbackorig = numshiftbackorig+1
		elseif(x > xyzmax) then
			x = x - dxyz
			numshiftbackorig = numshiftbackorig+1
		endif
		if(y < xyzmin) then
			y = y + dxyz
			numshiftbackorig = numshiftbackorig+1
		elseif(y > xyzmax) then
			y = y - dxyz
			numshiftbackorig = numshiftbackorig+1
		endif
		if(z < xyzmin) then
			z = z + dxyz
			numshiftbackorig = numshiftbackorig+1
		elseif(z > xyzmax) then
			z = z - dxyz
			numshiftbackorig = numshiftbackorig+1
		endif
		! rescale x,y,z: take into consideration the rescale of volume due to different cosmology!!! fix Om=0.26/w=-1.0 as right
		x = x/rescalefac
		y = y/rescalefac
		z = z/rescalefac
		if(shiftx) then
			zobs = redshift + vx * (1.0_dl+redshift) / const_c 
			shiftdist = de_get_comovr(zobs) - dist
			shiftedcord = x + shiftdist
			if(shiftedcord < xyzmin/rescalefac) then
				shiftedcord = shiftedcord + dxyz
				numshiftback = numshiftback+1
			elseif(shiftedcord > xyzmax/rescalefac) then
				shiftedcord = shiftedcord - dxyz
				numshiftback = numshiftback+1
			endif
			write(10,'(3e15.7)') shiftedcord, y, z
		endif
		if(shifty) then
			zobs = redshift + vy * (1.0_dl+redshift) / const_c 
			shiftdist = de_get_comovr(zobs) - dist
			shiftedcord = y + shiftdist
			if(shiftedcord < xyzmin/rescalefac) then
				shiftedcord = shiftedcord + dxyz
				numshiftback = numshiftback+1
			elseif(shiftedcord > xyzmax/rescalefac) then
				shiftedcord = shiftedcord - dxyz
				numshiftback = numshiftback+1
			endif
			write(20,'(3e15.7)') x, shiftedcord, z
		endif
		if(shiftz) then
			zobs = redshift + vz * (1.0_dl+redshift) / const_c 
			shiftdist = de_get_comovr(zobs) - dist
			shiftedcord = z + shiftdist
			if(shiftedcord < xyzmin/rescalefac) then
				shiftedcord = shiftedcord + dxyz
				numshiftback = numshiftback+1
			elseif(shiftedcord > xyzmax/rescalefac) then
				shiftedcord = shiftedcord - dxyz
				numshiftback = numshiftback+1
			endif
			write(30,'(3e15.7)') x, y, shiftedcord
		endif
		cycle
101		exit
	enddo
	
	if(shiftx) close(10)
	if(shifty) close(20)
	if(shiftz) close(30)

	write(*,'(A,i10,A,2i10)') ' Finishing processing ', nlines, ' lines; num shifted back (orig, due to RSD) = ', &
		numshiftbackorig, numshiftback
	write(100,'(A,i10,A,2i10)') ' Finishing processing ', nlines, ' lines; num shifted back (orig, due to RSD) = ', &
		numshiftbackorig, numshiftback

	write(*,'(A,3("(",2f15.7,"); "))') 'Range of x,y,z: ', xmin,xmax,ymin,ymax,zmin,zmax
	write(100,'(A,3("(",2f15.7,"); "))') 'Range of x,y,z: ', xmin,xmax,ymin,ymax,zmin,zmax

	close(100)

end program main
