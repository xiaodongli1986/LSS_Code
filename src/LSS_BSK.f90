
module LSS_BSK
use LSS_chisq
implicit none


	! Type beta skeleton connection
	type BSK_Connect
		integer :: i1, i2 ! index of the two points 
		real(dl) :: r, length, mu ! distance (middle point), length, mu of the connection
	end type

	! Global collection of BSK information
	type(BSK_Connect), allocatable :: gb_BSKs(:)
	
	! number of BSK connection
	integer :: gb_numBSKConnect

contains

!!! Beta Skeleton
  !------------------------------------------
  ! beta skeleton: check wether the two points can be connected
  !   Focus on the circle centered at xA,yA,zA; check 
  !------------------------------------------
  	subroutine BSK_Connected(allparindex1, allparindex2,  beta, Connected, erflaglog)
		! Dummy
  		integer, intent(in) :: allparindex1, allparindex2
  		real(dl), intent(in) :: beta
  		logical, intent(out) :: Connected, erflaglog
  		! local
  		integer :: ixl, ixr, iyl, iyr, izl, izr, ix,iy,iz, dix,diy,diz, l, allparindex
  		real(dl) :: x1,y1,z1,x2,y2,z2, xA,yA,zA, xB,yB,zB, pA(3),pB(3),nowP(3) ! centers of the two spheres
  		real(dl) :: betaover2, radiusq, radius
  		real(dl) :: xl,xr,yl,yr,zl,zr, nowrAsq,nowrBsq, x,y,z
  		
		erflaglog = .false.
  		
  		if(beta<1.0_dl) then
  			print *, 'ERROR (BSK_Connected)! Beta < 1 code not finished yet!: beta = ', beta
  			stop
  		endif
  		
		x1 = gb_xyz_list(1,allparindex1); y1 = gb_xyz_list(2,allparindex1); z1 = gb_xyz_list(3,allparindex1)
		x2 = gb_xyz_list(1,allparindex2); y2 = gb_xyz_list(2,allparindex2); z2 = gb_xyz_list(3,allparindex2)

  		betaover2 = beta * 0.5_dl
  		! centers, radius of the two circles
  		xA = (1.0_dl-betaover2)*x1 + betaover2*x2
  		yA = (1.0_dl-betaover2)*y1 + betaover2*y2
  		zA = (1.0_dl-betaover2)*z1 + betaover2*z2
  		xB = betaover2 * x1 + (1.0_dl -betaover2)*x2
  		yB = betaover2 * y1 + (1.0_dl -betaover2)*y2
  		zB = betaover2 * z1 + (1.0_dl -betaover2)*z2
  		
  		radius  = sqrt((x2-x1)**2.0 + (y2-y1)**2.0 + (z2-z1)**2.0)*betaover2
  		radiusq = radius**2.0
  		
!  		open(unit=1,file='BSK_Spheres.txt')
!  		write(1,'(3(e14.7,1x))') x1, y1, z1
!  		write(1,'(3(e14.7,1x))') x2, y2, z2
!  		write(1,'(3(e14.7,1x))') xA, yA, zA 
!  		write(1,'(3(e14.7,1x))') xB, yB, zB
!  		write(1,'(3(e14.7,1x))') radius, radius, radius
!  		close(1)
  		
  		pA(1)=xA; pA(2)=yA; pA(3)=zA
  		pB(1)=xB; pB(2)=yB; pB(3)=zB
  		  		
  		! xl,xr,yl,yr,zl,zr is a sphere which cover the sphere near pA
  		! xl,yl,zl are the boundary which is relatively more close to point B
  		! we will start scanning points at that side
		if(xB>xA) then
			xl=xA+radius; xr=xA-radius; dix=-1
		else
			xl=xA-radius; xr=xA+radius; dix=1
		endif
		if(yB>yA) then
			yl=yA+radius; yr=yA-radius; diy=-1
		else
			yl=yA-radius; yr=yA+radius; diy=1
		endif
		if(zB>zA) then
			zl=zA+radius; zr=zA-radius; diz=-1
		else
			zl=zA-radius; zr=zA+radius; diz=1
		endif

		! range of ix,iy,iz determined by xyz l/r
		ixl = int((xl-gbgridxmin)/gbdeltax+1.0_dl);ixr = int((xr-gbgridxmin)/gbdeltax+1.0_dl)
		iyl = int((yl-gbgridymin)/gbdeltay+1.0_dl);iyr = int((yr-gbgridymin)/gbdeltay+1.0_dl)		
		izl = int((zl-gbgridzmin)/gbdeltaz+1.0_dl);izr = int((zr-gbgridzmin)/gbdeltaz+1.0_dl)

		ixl = min(gb_n_cellx,ixl); ixl = max(1,ixl); iyl = min(gb_n_celly,iyl); iyl = max(1,iyl)
		izl = min(gb_n_cellz,izl); izl = max(1,izl); ixr = min(gb_n_cellx,ixr); ixr = max(1,ixr)
		iyr = min(gb_n_celly,iyr); iyr = max(1,iyr); izr = min(gb_n_cellz,izr); izr = max(1,izr)

		do ix = ixl, ixr, dix
		do iy = iyl, iyr, diy
		do iz = izl, izr, diz
			do l = 1, gb_cell_mat(ix,iy,iz)%numdata
				allparindex = gb_cell_mat(ix,iy,iz)%idatalist(l)
				if(allparindex.eq.allparindex1.or.allparindex.eq.allparindex2) cycle
				x = gb_xyz_list(1,allparindex)
				y = gb_xyz_list(2,allparindex)
				z = gb_xyz_list(3,allparindex)
				nowrAsq = (x-xA)*(x-xA) + (y-yA)*(y-yA) + (z-zA)*(z-zA)
				if(nowrAsq>radiusq) cycle
				nowrBsq = (x-xB)*(x-xB) + (y-yB)*(y-yB) + (z-zB)*(z-zB)
				if(nowrBsq>radiusq) cycle
				Connected = .false.
				return
			enddo				
		enddo
		enddo
		enddo 
		Connected = .true.		
  	end subroutine BSK_Connected

  !------------------------------------------
  ! beta skeleton: do beta skeleton for all data points...
  !------------------------------------------
	subroutine BSK(inputfilename, omegam, w, beta, outputfilename, printinfo,numNNB, minimalr_cut, maximalr_cut, &
		do_nglcrosscheck, ngldir, outputBSKindex, outputBSKinfo, outputxyz, do_init)
		! Dummy 
		real(dl), intent(in) :: omegam, w, beta
		character(len=char_len), intent(in) :: inputfilename, outputfilename, ngldir
		real(dl), intent(in), optional :: minimalr_cut, maximalr_cut
		logical, intent(in) :: printinfo, do_nglcrosscheck, outputBSKindex, outputBSKinfo, outputxyz, do_init
		integer, intent(in) :: numNNB
		! Local
		integer :: i, j, i1, ii2,i2, idiv, di1, tmp_numBSK
		logical :: Connected,erflag
		integer, allocatable :: selected_list(:)
		character(len=char_len) :: omwstr,tmpstr1,tmpstr2, file_index, file_ngl, file_xyz, file_BSKinfo, cmdstr
		real(dl) :: time1, time2, x,y,z, dx,dy,dz
		type(BSK_Connect), allocatable :: tmp_BSKs(:),tmp_BSKs2(:) ! temporary collection of BSKS

		! minimal, maximal cut
		if(present(minimalr_cut)) then
			gb_minimalr_cut = minimalr_cut
		else
			gb_minimalr_cut = -1.0e30
		endif
		if(present(maximalr_cut)) then
			gb_maximalr_cut = maximalr_cut
		else
			gb_maximalr_cut = 1.0e30
		endif
		
		! initialization
		print *, '  (BSK) Begins.'
		gb_datafile = inputfilename
		if(do_init) then
			call cosmo_funs_init(printinfo)
			call readin_dataran(printinfo)
			call init_cosmo(omegam,w,h_dft,printinfo)
			call init_mult_lists(printinfo)
			call do_cell_init((real(gb_numdata)**0.33_dl), printinfo)		
		endif

		! make directory
		cmdstr = 'mkdir -p '//trim(adjustl(outputfilename))
		call system(cmdstr)

		! names of files
		write(tmpstr1,'(f8.4)') omegam
		write(tmpstr2,'(f8.4)') w
		omwstr = 'om'//trim(adjustl(tmpstr1))//'_w'//trim(adjustl(tmpstr2))
		file_index = trim(adjustl(outputfilename))//'/'//trim(adjustl(omwstr))
		file_ngl = trim(adjustl(file_index))//'.nglBSK'
		file_index = trim(adjustl(file_index))//'.BSKIndex'

		! Calculation: preparation
		tmp_numBSK = gb_numdata * 1.5
		allocate(tmp_BSKs(tmp_numBSK))
		call cpu_time(time1)
		print *, '  (BSK) Calculating BSK...'
		if(outputBSKindex) then
			write(*,'(A,A)'), '  *** Result wrote to: ', trim(adjustl(file_index))
			open(unit=93408,file=file_index)
		endif
		idiv = 100
		j = 0
		di1 = floor(gb_numdata/dble(idiv))
		
		! Calculation: loop
		do i1 = 1, gb_numdata
			! avoid redudant calculation
			if(mod(i1,di1).eq.1) then
				write(*,'(i2,A\)') (i1/di1)*100/idiv, '%==>'
			endif
!			call NNBSearch_data(gb_xyz_list(1,i1),gb_xyz_list(2,i1),gb_xyz_list(3,i1),&
!				numNNB,1.0_dl,selected_list) 
			if(allocated(selected_list)) deallocate(selected_list); allocate(selected_list(numNNB))
			call NNBExactSearch_data(gb_xyz_list(1,i1),gb_xyz_list(2,i1),gb_xyz_list(3,i1),&
				numNNB,selected_list,erflag) 
			if(erflag) then
				print *, 'Unexpected Error Encountered (i1 = )', i1,'; switch to Non-Exact NNBSearch_data'
				if(allocated(selected_list)) deallocate(selected_list)
				call NNBSearch_data(gb_xyz_list(1,i1),gb_xyz_list(2,i1),gb_xyz_list(3,i1),&
					numNNB*4,1.0_dl,selected_list) 
			endif
			do ii2 = 1, size(selected_list)
				i2 = selected_list(ii2)
				if (i1.ge.i2) cycle
				call BSK_Connected(i1,i2,beta,Connected,erflag)
				if(Connected) then
					if(outputBSKindex) then
						write(tmpstr1,*) i1-1
						write(tmpstr2,*) i2-1
						write(93408,'(A,A,A)') trim(adjustl(tmpstr1)), ' ', trim(adjustl(tmpstr2))
					endif
					j = j+1
					if(j > tmp_numBSK ) then
						print *
						write(*,*) 'Large numer of connections; re-allocate array: ', &
							tmp_numBSK,int(tmp_numBSK*1.5)
						allocate(tmp_BSKs2(tmp_numBSK)) 
						tmp_BSKs2 = tmp_BSKs
						tmp_numBSK = tmp_numBSK*1.5; deallocate(tmp_BSKs); allocate(tmp_BSKs(tmp_numBSK))
						tmp_BSKs(1:j-1) = tmp_BSKs2(1:j-1); deallocate(tmp_BSKs2)
					endif
					! index
					tmp_BSKs(j)%i1 = i1; tmp_BSKs(j)%i2 = i2
					! distance
					x = (gb_xyz_list(1,i1)+gb_xyz_list(1,i2))*0.5
					y = (gb_xyz_list(2,i1)+gb_xyz_list(2,i2))*0.5
					z = (gb_xyz_list(3,i1)+gb_xyz_list(3,i2))*0.5
					tmp_BSKs(j)%r  = sqrt(x*x+y*y+z*z)
					! length
					dx = (gb_xyz_list(1,i1)-gb_xyz_list(1,i2))*0.5
					dy = (gb_xyz_list(2,i1)-gb_xyz_list(2,i2))*0.5
					dz = (gb_xyz_list(3,i1)-gb_xyz_list(3,i2))*0.5
					tmp_BSKs(j)%length  = sqrt(dx*dx+dy*dy+dz*dz)
					! mu
					tmp_BSKs(j)%mu = (x*dx+y*dy+z*dz)/tmp_BSKs(j)%r/tmp_BSKs(j)%length
				endif
			enddo
		enddo
		if(outputBSKindex) close(93408)
		call cpu_time(time2)
		print *, '       Time used in BSK calculation: ', time2-time1

		! store the data into gb_BSKs
		gb_numBSKConnect = j
		if(allocated(gb_BSKs)) deallocate(gb_BSKs)
		allocate(gb_BSKs(gb_numBSKConnect))
		gb_BSKs(1:gb_numBSKConnect) = tmp_BSKs(1:gb_numBSKConnect) 
		deallocate(tmp_BSKs)

		! information of BSK Connections		
		print *, '       # of Connections: ', gb_numBSKConnect
		file_BSKinfo =  trim(adjustl(outputfilename))//'/'//trim(adjustl(omwstr))//'.BSKinfo'
		if(outputBSKinfo) then
			write(*,'(A,A)'), '  *** information of BSK (fmt: distance, redshift, length, mu) wrote to: ', &
				trim(adjustl(file_BSKinfo))
			open(unit=2198,file=file_BSKinfo)
			do i = 1, gb_numBSKConnect
				write(2198,'(4e15.7)') gb_BSKs(i)%r, de_zfromintpl(gb_BSKs(i)%r), gb_BSKs(i)%length, gb_BSKs(i)%mu
			enddo
			close(2198)
		endif
		
		! output xyz (cosmology may be different from the input cosmology, so we write xyz calcuated in the new cosmology)
		file_xyz = trim(adjustl(outputfilename))//'/'//trim(adjustl(omwstr))//'_xyz.txt'
		if(outputxyz.or.do_nglcrosscheck) then
			print *, ' *** xyz wrote to: ', trim(adjustl(file_xyz))
			write(*,'(7x,A,A)'), '  (cosmology may be different from the input cosmology, ',&
				'so we write xyz calcuated in the new cosmology)'
			open(unit=198,file=file_xyz)
			do i =1, gb_numdata
				write(198,*) gb_xyz_list(1:3,i)
			enddo
			close(198)
		endif

		!!! Corss-Check using ngl
		! do the same calculation using ngl_beta ! for check of result
		! comment it in real calculation
!			subroutine gd_mldprho_chi2s(cs, changenuminx, chisqlist, dfchisqlist, multdfchisqlist, numchisq)
		if(do_nglcrosscheck) then
			print *, '       Now doing similar job using ngl!!! (declare your ngl directory in LSS_smooth.f90/BSK)'
			print *, ' *** Result of ngl wrote to: ', trim(adjustl(file_ngl))
			cmdstr = 'rm '//trim(adjustl(file_ngl))
			call system(cmdstr)
			write(tmpstr1,'(f5.1)') beta		
			cmdstr = trim(adjustl(ngldir))//' -i '//trim(adjustl(file_xyz))//' -d 3 -b '//trim(adjustl(tmpstr1))//' >> '//trim(adjustl(file_ngl))
			write(*,'(7x, A)') trim(adjustl(cmdstr))
			call system(cmdstr)
		endif
		print *, '  (BSK) Done.'
	end subroutine BSK
	
  !------------------------------------------
  ! beta skeleton: statistical analysis
  !------------------------------------------	
	subroutine BSK_stat(omegam, w, beta, outputfilename, numdrop, dropstep, printinfo)
		! Dummy
		real(dl), intent(in) :: omegam, w, beta, dropstep
		character(len=char_len) :: outputfilename
		integer, intent(in) :: numdrop
		logical, intent(in) :: printinfo
		! Local
		type(chisq_settings) :: cs
		real(dl) :: chisqlist(numdrop), dfchisqlist(numdrop), multdfchisqlist(nbinchisq,numdrop), x,y,z,dx,dy,dz
		integer :: i,j, i1,i2
		character(len=char_len) :: tmpstr1, tmpstr2

		cs%print_info = printinfo
		call init_cosmo(omegam,w,h_dft,printinfo)
		call init_mult_lists(printinfo)
		call do_cell_init((real(gb_numdata)**0.33_dl), printinfo)	

		! gb_rho/drho lists... 
		! pos will be position; rho will be length of connection; drho will be dx,dy,dz... Oh my god.
		if(allocated(gb_pos_list)) deallocate(gb_pos_list)
		if(allocated(gb_rho_list)) deallocate(gb_rho_list)
		if(allocated(gb_drho_list)) deallocate(gb_drho_list)
		allocate(gb_pos_list(3,gb_numBSKConnect), gb_rho_list(gb_numBSKConnect), gb_drho_list(3,gb_numBSKConnect))
		do i = 1, gb_numBSKConnect
			i1 = gb_BSKs(i)%i1  
			i2 = gb_BSKs(i)%i2
			gb_pos_list(1,i) = (gb_xyz_list(1,i1)+gb_xyz_list(1,i2))*0.5
			gb_pos_list(2,i) = (gb_xyz_list(2,i1)+gb_xyz_list(2,i2))*0.5
			gb_pos_list(3,i) = (gb_xyz_list(3,i1)+gb_xyz_list(3,i2))*0.5
			! length
			gb_drho_list(1,i) = (gb_xyz_list(1,i1)-gb_xyz_list(1,i2))*0.5
			gb_drho_list(2,i) = (gb_xyz_list(2,i1)-gb_xyz_list(2,i2))*0.5
			gb_drho_list(3,i) = (gb_xyz_list(3,i1)-gb_xyz_list(3,i2))*0.5
		enddo
		
		! outputname
		write(tmpstr1,'(f8.4)') omegam
		write(tmpstr2,'(f8.4)') w
		gb_omwstr = 'om'//trim(adjustl(tmpstr1))//'_w'//trim(adjustl(tmpstr2))
		gb_suboutput_name = trim(adjustl(outputfilename))//'/'//trim(adjustl(gb_omwstr))

		! initilize dropping		
		cs%numdrop = numdrop
		allocate(cs%dropval(cs%numdrop),  cs%lowdropvalratio(cs%numdrop),  cs%highdropvalratio(cs%numdrop), &
			cs%dropdval(cs%numdrop), cs%lowdropdvalratio(cs%numdrop), cs%highdropdvalratio(cs%numdrop))
		do i = 1, cs%numdrop
			cs%dropval(i) = .true.;	cs%lowdropvalratio(i) = min((i-1)*dropstep,0.95_dl) 
			cs%dropdval(i) = .false.; cs%lowdropdvalratio(i) = 0.0
			cs%highdropvalratio(i) = 0.0
			cs%highdropdvalratio(i) = 0.0
		enddo

		write(*,'(A,i3,A)') '  Droppings (',cs%numdrop,'):'
		do i = 1, cs%numdrop
			write(*,'(A,L2,f6.3,f6.3)') '    Length of edge: ', cs%dropval(i), real(cs%lowdropvalratio(i)),  &
				real(cs%highdropvalratio(i))
		enddo

		call gd_mldprho_chi2s(cs, 0.0_dl, chisqlist, dfchisqlist, multdfchisqlist, numdrop, BSKmode=.true.)
		
		write(*,*) '  (BSK_stat) results of chisqs:'
		write(*,'(36x,<cs%numdrop>(f8.3,1x))')	chisqlist(1:cs%numdrop)
		write(*,'(A,12x,<cs%numdrop>(f8.3,1x))') '    2 bins (old method):', dfchisqlist(1:cs%numdrop)
		do j = 1, nbinchisq	
			write(*,'(4x,i4,A,24x,<cs%numdrop>(f8.3,1x))') nbinchisqlist(j),' bins:', multdfchisqlist(j,1:cs%numdrop)
		enddo
		
	end subroutine BSK_stat

end module LSS_BSK
