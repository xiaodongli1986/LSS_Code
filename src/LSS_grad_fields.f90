
!!! 

!!! TBF (June 17)
	! Code to do weight check in the grid_rho_drho_list: Done!
	! This has been graphically tested: roughly we find the code shall be correct.
	! Next we shall improve the code to do weight calculation & check in 2D: Done!
	! Check the normalization: properties of nbar before/after normalized; check the boundary carefully!!!:

! Please check dft_ra_ratio in LSS_smooth.f90.
!
! Quick find is successful...
!  Now, we can make quick rho_delta...
!   For each pixel, we save its xyzrlist & indexlist, for further usage ...
!   Adjust # of cells to make nb_list more efficient...

!####################################
!This module does smooth
!####################################
module LSS_grad_fields
use LSS_smooth
use LSS_mu_tools
	implicit none

!!! Fixed Settings

	! List of nbins (separate the whole sample into nbins and calculate \bar\mu in each bin...)
	integer, parameter :: nbinchisq = 19
	integer, parameter :: nbinchisqlist(nbinchisq) = (/2, 3, 4, 5, 6, 7, 8, &
                9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 /)!, 21, 22, 23, 24, 25, 26, 27, 28, &
!                29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40/)!, 41, 42, 43, 44, 45, 46, 47, 48, &
!                49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, &
 !               69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, &
  !              89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100 /)

	! Fixed maximal distance smoothing (Not used now!)
	logical,  public, parameter :: gb_do_seg_cut = .false.
	real(dl), public, parameter :: gb_seg_cut_dist = 5.0_dl

!!! Adjustable & To be fixed Settings	
! All variables in the "Adjustable & To be fixed Settings" section must be assigned a value in the main**.f90!!!

	! smoothing method: fixed radius of smoothing sphere (maximal distance we can go away from the position)
	logical, public :: gb_use_fixmd
	real(dl), public :: gb_fixmd

	! Keep gb_bddist (distance between smoothing-sphere surface and bd) / gb_bufferdist (smoothing sphere center and boundary) 
	!  away from boundary (x,y,z,r,ra,dec) to avoid bd effect;
	!  for the boundary r, keep away with distance **_rextra
	real(dl), public :: gb_bddist, gb_bddist_rextra, gb_bufferdist, gb_bufferdist_rextra

!!! Important Variables

	! global collection of postions, rhos, drhos
	real(dl), allocatable :: gb_pos_list(:,:), gb_rho_list(:), gb_drho_list(:,:)
	
contains	


  !------------------------------------------
  ! check whether a point has boundary effect
  !  for BOSSDR12
  !------------------------------------------
	logical function has_boundary_effect(datatype, x,y,z,max_dist,printinfo)
		! Dummy
		real(dl), intent(in) :: x,y,z,max_dist !max_dist is the radius of smoothing sphere
		logical, intent(in), optional :: printinfo
		integer, intent(in) :: datatype
		!local 
		real(dl) :: ra,dec,r,dr, rdr ! An extra distance away from the boundary
                logical :: checkxyz = .true.
                real(dl) :: ddl,ddr,ddu,ddd,minimal_dradec
                integer :: i_check

		if(present(printinfo)) print *, ' (has_boundary_effect) Start. Checking max_dist: ', max_dist

		if(max_dist .ge. 1.0e5) then !smoothing sphere cannot be ridiculously large; ...
!			print *, 'Encountered max_dist ', max_dist
			has_boundary_effect = .true.; return
		endif

		call xyz_to_radecr(x,y,z,ra,dec,r)
		if(present(printinfo)) print *, ' (has_boundary_effect) Converting xyz to radecr: ', ra,dec,r
		! ra dec within the range of data ra/dec
                if(ra>gbramindata.and.ra<gbramaxdata.and.dec>gbdecmindata.and.dec<gbdecmaxdata) then
			goto 861
		endif
		has_boundary_effect = .true.; return ! this will reject the Nan of ra,dec,r. AMTB...
861		continue		

		! Special check of boundary for SGC (there the ra range is 0-360; in fact there is large empty region between that...)
		!/home/lixiaodong/SparseFilaments/data/input/DR12/cmass-dr12v1-S-Anderson.dat_t1.txt
		!45.0657669266 316.361423817
		!/home/lixiaodong/SparseFilaments/data/input/DR12/lowz-dr12v1-S-Anderson.dat_t1.txt
		!45.0813082213 316.381216501
		!/home/lixiaodong/SparseFilaments/data/input/DR12/cmass-dr12v4-S-Reid.dat_t1.txt
		!45.4518156117 316.361423817
		!/home/lixiaodong/SparseFilaments/data/input/DR12/lowz-dr12v4-S-Reid.dat_t1.txt
		!45.441787537 316.381216501
		if(.true.) then
			if(trim(adjustl(gb_catalogue_name)).eq.'DR12v1-CMASS-S' &
				.or.trim(adjustl(gb_catalogue_name)).eq.'DR12v1-LOWZ-S' &
				.or.trim(adjustl(gb_catalogue_name)).eq.'DR12v4-CMASS-S' &
				.or.trim(adjustl(gb_catalogue_name)).eq.'DR12v4-LOWZ-S') then
				if(datatype.eq.gb_dt_nbarz_NoRSDM .or. datatype.eq.gb_dt_nbarz_RSDM .or. &
					 datatype.eq.gb_dt_realdata) then
					if(ra>46.0 .and.ra<316.0) then
						has_boundary_effect = .true.; return
					endif
				endif	
			endif
		endif


		do i_check = 1, 2
			if(i_check.eq. 1) then
				dr = max_dist + gb_bddist! distance for the center of sphere to keep away from boundary
				rdr = dr + gb_bddist_rextra ! for r boundary
			elseif(i_check.eq.2) then
				dr = gb_bufferdist
				rdr = dr + gb_bufferdist_rextra
			endif

			! check r boundary
			if (abs(r-gbrmindata)<rdr .or.abs(r-gbrmaxdata)<rdr .or.abs(r-gbgridrmax)<rdr) then
				has_boundary_effect = .true.; return
			endif
			! check x,y,z boundary
		        if(checkxyz) then
				if(abs(x-gbxmindata)<dr.or.abs(y-gbymindata)<dr.or.abs(z-gbzmindata)<dr &
		       			.or.abs(x-gbxmaxdata)<dr.or.abs(y-gbymaxdata)<dr.or.abs(z-gbzmaxdata)<dr &
					.or.abs(x-gbgridxmin)<dr.or.abs(y-gbgridymin)<dr.or.abs(z-gbgridzmin)<dr &
					.or.abs(x-gbgridxmax)<dr.or.abs(y-gbgridymax)<dr.or.abs(z-gbgridzmax)<dr) then
						has_boundary_effect = .true.; return
				endif
		        endif

			if(present(printinfo)) print *, ' (has_boundary_effect) Before checking ra,dec buffer'
			! check ra,dec buffer
			if(allocated(gb_radecbdmat_raleft)) then
				minimal_dradec = dr / r * 180.0d0/const_pi
				call dist_to_radecbd(ra,dec,ddl,ddr,ddu,ddd)
				if(.not.(ddl>minimal_dradec .and. ddr>minimal_dradec &
					.and. ddu>minimal_dradec .and. ddd>minimal_dradec)) then
					has_boundary_effect = .true.; return
				endif
			endif
			
			if(present(printinfo)) print *, ' (has_boundary_effect) After checking ra,dec buffer'
		enddo

		if(present(printinfo)) print *, ' (has_boundary_effect) Done '
		has_boundary_effect=.false.
	end function has_boundary_effect

  !------------------------------------------
  ! estimating rho/drho at grid points
  !------------------------------------------
	subroutine grid_rho_drho_list(smnum, printinfo, pos_list, rho_list, drho_list)
		!Dummy
		integer, intent(in) :: smnum 
		logical, intent(in) :: printinfo
		real(dl), allocatable, intent(out) :: pos_list(:,:),  rho_list(:),  drho_list(:,:)
		! Local
		real(dl) :: x,y,z,ra,dec,r,deltaradec, rho,drhox,drhoy,drhoz,max_dist,weirat, &
			numhasbd, numsmallsmnum, numavg, numsqavg, numranweirej, insphere_numran_avg, &
			numranweifail1, numranweifail2, sumwei_acpt, sumwei_all, drho,mu
		integer, parameter :: n_ransmnum = 6, ransmnum(n_ransmnum) = (/10, 20, 50, 100, 200, 300/)
		integer :: i,j, ix, iy, iz, n, iscan, n0, n1, n2, fixmdsmnum, erflagint, insphere_numran, percent, &
			insphere_numran_min, insphere_numran_max, numransmnum(n_ransmnum), num_scan, check_iscan
		logical ::  erflaglog, tmplogical
		character(len=char_len) :: tmpstr1, tmpstr2, tmpstr3
		
		if(printinfo) then
			print *, '  (grid_rho_drho_list begins) Estimating rho_list and its gradient...'
			write(*,'(26x,A,L2,L2)') 'use_fixmd / gb_do_seg_cut = ', gb_use_fixmd, gb_do_seg_cut
			if(gb_use_fixmd) then
				write(*,'(26x,A,f10.3)') 'Fixed smoothing length gb_fixmd = ', real(gb_fixmd)
			else
				write(*,'(26x,A,i4)') 'Using #-NNB = ', smnum
				write(*,'(26x,A,f6.3,A,L3)') 'NNB Over-Searching Rat = ', gb_NNBosrat, &
					'; Do NNB-ReSearch = ', gb_do_NNBReSearch
			endif
			if(gb_do_seg_cut)  then
				write(*,'(26x,A,f10.3)') 'Applying minimal-r cut: ', real(gb_seg_cut_dist)
				write(*,*) 'ERROR (grid_rho_drho_list): gb_do_seg_cut not supported now!'; stop
			endif
			write(*,'(26x,A)') 'Settings of Boundary Correction: '
			write(*,'(30x,A,4f10.3)') 'bd, bd_rextra, buffer, buffer_rextra = ', &
				gb_bddist, gb_bddist_rextra, gb_bufferdist, gb_bufferdist_rextra
		endif
		
		n = gb_n_cellx*gb_n_celly*gb_n_cellz
		write(*,'(25x,A,$)') ' '
		! n2 counts how many pixels in shell and has no boundary effect (pass has_boundary_effect)
		numransmnum = 0; num_scan= gb_n_cellx*gb_n_celly*gb_n_cellz
		iscan=0; n0=0; n1=0; n2=0; gb_num_NNBSearch=0; gb_num_NNBReSearch = 0;
		numavg = 0.0; numsqavg = 0.0; 
		insphere_numran_avg = 0.0; insphere_numran_min = max_insphere_num; insphere_numran_max = -max_insphere_num
		numsmallsmnum = 0.0; numhasbd = 0.0; numranweirej=0.0; numranweifail1=0.0; numranweifail2=0.0
		
		check_iscan = -1
		
!		open(unit=98,file='checkiscan.txt') !!! Useful sentence: check where the programe break
		do ix = 1, gb_n_cellx
		do iy = 1, gb_n_celly
		do iz = 1, gb_n_cellz
			iscan = iscan+1
!			write(98,*) iscan
			if(iscan.eq.check_iscan) print *, ' (grid_rho_drho_list) Checking iscan: ', iscan, 'Start'
			if(mod(iscan,(num_scan/10)).eq.0.and.iscan>num_scan/10-10) then
                                percent = int(iscan/dble(num_scan/10.0)+0.5)*10
                                write(tmpstr1,*) percent
                                write(tmpstr2,'(f6.3)') gb_omegam
                                write(tmpstr3,'(f6.3)') gb_w
                                write(*,'(A,$)') trim(adjustl(tmpstr1))//'%('//trim(adjustl(tmpstr2))&
					//','//trim(adjustl(tmpstr3))//')==>'
                        endif
                        if(iscan .eq. num_scan) then
                                write(*,'(A)') 'Finish'
                        endif
			gb_cell_mat(ix,iy,iz)%rhodrhoindex = -1
			gb_cell_mat(ix,iy,iz)%weirat = -1.0
			call cell_pos(ix,iy,iz,x,y,z)
			r = sqrt(x*x+y*y+z*z)
			if(r<gbrmindata.or.r>gbrmaxdata) then
				cycle
			endif
			n0=n0+1

			if(iscan.eq.check_iscan) print *, ' (grid_rho_drho_list) Before Compute rho,drho...'
			! Computing rho / drhos ...
			if(gb_use_fixmd) then
				print *, 'Check use_fixmd bd correction carefully! especially check combination 2d&3d!'; stop
				max_dist = gb_fixmd
				if(has_boundary_effect(gb_i_datatype, x,y,z,max_dist)) then
					numhasbd = numhasbd + 1.0
					cycle
				else
					call nb_fixmd_list(x,y,z,fixmdsmnum,rho,drhox,drhoy,drhoz,gb_fixmd,erflaglog)
				endif
			else
				if(iscan.eq.check_iscan) print *, iscan, 'Before nb_list0: x,y,z = ', x,y,z
				if(iscan.eq.check_iscan) call nb_list0(x,y,z,smnum,rho,drhox,drhoy,drhoz,max_dist,erflaglog,.true.)
				call nb_list0(x,y,z,smnum,rho,drhox,drhoy,drhoz,max_dist,erflaglog)
				if(iscan.eq.check_iscan) then
					print *, ' (grid_rho_drho_list) has_boundary_effect...'
					tmplogical = has_boundary_effect(gb_i_datatype, x,y,z,max_dist,.true.)
				endif
				if(has_boundary_effect(gb_i_datatype, x,y,z,max_dist) .or. erflaglog) then
					numhasbd = numhasbd + 1.0; cycle
				endif
				if(iscan.eq.check_iscan) print *, ' (grid_rho_drho_list) After has_boundary_effect'
			endif
			if(iscan.eq.check_iscan) print *, iscan, ' (grid_rho_drho_list) After Compute rho,drho'

			if(gb_i_rantype.eq.gb_3dran) then
				call nb_fixmd_ranwei(x,y,z,max_dist,insphere_numran,weirat,erflagint,sumwei_acpt, sumwei_all)
			elseif(gb_i_rantype.eq.gb_radecran) then
				call xyz_to_radecr(x,y,z,ra,dec,r)
				deltaradec = max_dist/r * 180.0d0/const_pi
				call nb_fixmd_radecranwei(ra,dec,deltaradec,insphere_numran,weirat,erflagint,sumwei_acpt,sumwei_all)
			elseif(gb_i_rantype.eq.gb_radec3dran) then ! combination of 3d/radec random
				call xyz_to_radecr(x,y,z,ra,dec,r)
				if(r<gb_3dranmaxr) then ! use 3dran for small radius
					call nb_fixmd_ranwei(x,y,z,max_dist,insphere_numran,weirat,erflagint,sumwei_acpt, sumwei_all)
				else
					deltaradec = max_dist/r * 180.0d0/const_pi
					call nb_fixmd_radecranwei(ra,dec,deltaradec,insphere_numran,weirat,&
						erflagint,sumwei_acpt,sumwei_all)
				endif
			elseif(gb_i_rantype.eq.gb_noran) then
				erflagint = 0
				weirat = 1.0
			else
				print *, 'ERROR (grid_rho_drho_list)!!: Unknown value of gb_i_rantype: ', gb_i_rantype
				stop
			endif
			
			if(erflagint.eq.1) then
				numranweifail1 = numranweifail1+1; cycle
			elseif(erflagint.eq.2) then
				numranweifail2 = numranweifail2+1; cycle
			elseif(erflagint.eq.0) then
				gb_cell_mat(ix,iy,iz)%weirat = weirat
				if(weirat < gb_ranwei_tol) then
					numranweirej = numranweirej+1;
					!!!! cycle ### CAUTION! Comment this to output test information!!
				endif
			else
				print *, 'WARNING (grid_rho_drho_list): rare/unexpected erflag from nb_fixmd_ranwei: ' , erflagint
!				stop
			endif

			n1 = n1+1
			gb_cell_mat(ix,iy,iz)%rhodrhos(1) = rho
			gb_cell_mat(ix,iy,iz)%rhodrhos(2) = drhox
			gb_cell_mat(ix,iy,iz)%rhodrhos(3) = drhoy
			gb_cell_mat(ix,iy,iz)%rhodrhos(4) = drhoz
			gb_cell_mat(ix,iy,iz)%maxdist = max_dist
			gb_cell_mat(ix,iy,iz)%insphere_numran = insphere_numran
			gb_cell_mat(ix,iy,iz)%sumwei_acpt = sumwei_acpt
			gb_cell_mat(ix,iy,iz)%sumwei_all = sumwei_all
			insphere_numran_avg = insphere_numran_avg+insphere_numran
			insphere_numran_min = min(insphere_numran_min,insphere_numran)
			insphere_numran_max = max(insphere_numran_max,insphere_numran)
			do i = 1, n_ransmnum
				if(insphere_numran < ransmnum(i)) numransmnum(i) = numransmnum(i)+1
			enddo
			if(weirat >= gb_ranwei_tol) then
				n2=n2+1	
				gb_cell_mat(ix,iy,iz)%rhodrhoindex = n2
			endif
			if(iscan.eq.check_iscan) print *, ' (grid_rho_drho_list)  Checking iscan ', iscan, 'End'
		enddo
		enddo
		enddo
		insphere_numran_avg = insphere_numran_avg/dble(n2)

		if(gb_outputinfo_gfcell) then
			tmpstr1 = '_gfcellrlt.txt'
		
			if(printinfo) then
				write(*,'(26x,A,A)') 'output information: ', trim(adjustl(tmpstr1))
			endif

			open(file=trim(adjustl(gb_suboutput_name))//trim(adjustl(tmpstr1)),unit=5)
		
			write(5,'(A)') '# fmt: 1 x,2 y,3 z,4 ra,5 dec,'//&
				'6 r,7 gb_cell_mat(ix,iy,iz)%maxdist, '//&
				'8 deltaradec, 9 gb_cell_mat(ix,iy,iz)%insphere_numran,'//&
				' 10 gb_cell_mat(ix,iy,iz)%weirat, '//&
				'11 gb_cell_mat(ix,iy,iz)%sumwei_acpt,'//&
				' 12 gb_cell_mat(ix,iy,iz)%sumwei_all, '//&
				'13 rho, 14 drhodx, 15 drhody, 16 drhodz,'//&
				' 17 drho, 18 mu, 19 redshift, 20 drho_los'

			i = 0
			do ix = 1, gb_n_cellx
			do iy = 1, gb_n_celly
			do iz = 1, gb_n_cellz
				i = i+1
				call cell_pos(ix,iy,iz,x,y,z)
				call xyz_to_radecr(x,y,z,ra,dec,r)
				deltaradec = gb_cell_mat(ix,iy,iz)%maxdist/r * 180.0d0/const_pi
				if(gb_cell_mat(ix,iy,iz)%weirat .ge. gb_ranwei_tol) then
					drhox = gb_cell_mat(ix,iy,iz)%rhodrhos(2)
					drhoy = gb_cell_mat(ix,iy,iz)%rhodrhos(3)
					drhoz = gb_cell_mat(ix,iy,iz)%rhodrhos(4)
					drho = sqrt(drhox**2.0+drhoy**2.0+drhoz**2.0)
					mu = get_mu_of_gradient(x,y,z,drhox,drhoy,drhoz)
					write(5,'(8e15.7,i10,11e15.7)') x,y,z,ra,dec,r, gb_cell_mat(ix,iy,iz)%maxdist, deltaradec, &
						gb_cell_mat(ix,iy,iz)%insphere_numran, gb_cell_mat(ix,iy,iz)%weirat, &
						gb_cell_mat(ix,iy,iz)%sumwei_acpt, gb_cell_mat(ix,iy,iz)%sumwei_all, &
						gb_cell_mat(ix,iy,iz)%rhodrhos(1:4), drho, mu, de_zfromintpl(r), drho*mu

!					write(5,'(i10,3i4,3x,i8,i6,f10.5,4e15.7)') &
!						i, ix,iy,iz, gb_cell_mat(ix,iy,iz)%rhodrhoindex, &
!							gb_cell_mat(ix,iy,iz)%insphere_numran, gb_cell_mat(ix,iy,iz)%weirat, &
!							gb_cell_mat(ix,iy,iz)%rhodrhos(1:4)
				endif
			enddo
			enddo
			enddo
			close(5)
			close(6)
		endif
		
		if(printinfo.or.gbtp) then
			if(gb_do_NNBReSearch) then
				write(*,'(26x,i7,A,f5.2,A,i7,A,f5.2,A)') int(gb_num_NNBReSearch+0.1), ' (',&
					gb_num_NNBReSearch/dble(gb_num_NNBSearch)*100.0,&
					'%) NNB ReSearch conducted; Large osrat Error: ', &
					gb_num_NNBReSearchEr, '(', gb_num_NNBReSearchEr/dble(gb_num_NNBSearch)*100.0, '%)'
			endif
		
			write(*,'(26x,i7,A,f5.2,A)') int(numhasbd+0.1), ' (',numhasbd/dble(n0)*100.0,&
					'%) rejected - has_boundary_effect()'
			write(*,'(26x,i7,A,f5.2,A)') int(numranweifail2+0.1), ' (',numranweifail2/dble(n0)*100.0,&
					'%) rejected - random weight estimation failure 2'
			write(*,'(26x,i7,A,f5.2,A)') int(numranweifail1+0.1), ' (',numranweifail1/dble(n0)*100.0,&
					'%) rejected - random weight estimation failure 1'
			write(*,'(26x,i7,A,f5.2,A,f7.2,i3,i6)') int(numranweirej+0.1), ' (',numranweirej/dble(n0)*100.0,&
					'%) rejected - weirat < weitol; mean/min/max-# of insphere random  = ', &
						insphere_numran_avg, int(insphere_numran_min+0.5), int(insphere_numran_max+0.5)
			do i = 1, n_ransmnum
				write(*,'(30x,i7,A,f5.2,A,i5)') numransmnum(i), '(', numransmnum(i)/dble(n1)*100.0, &
					'%) pixels has #-ran smaller than ', ransmnum(i)
			enddo
		endif

		allocate(pos_list(3,n2),rho_list(n2),drho_list(3,n2))
		i = 0
		do ix = 1, gb_n_cellx
		do iy = 1, gb_n_celly
		do iz = 1, gb_n_cellz
			if(gb_cell_mat(ix,iy,iz)%rhodrhoindex .ge. 1) then
				i =i+1
				call cell_pos(ix,iy,iz,x,y,z)
				pos_list(1,i)=x; pos_list(2,i)=y; pos_list(3,i)=z; 
				rho_list(i) = gb_cell_mat(ix,iy,iz)%rhodrhos(1)
				drho_list(1,i) = gb_cell_mat(ix,iy,iz)%rhodrhos(2)
				drho_list(2,i) = gb_cell_mat(ix,iy,iz)%rhodrhos(3)
				drho_list(3,i) = gb_cell_mat(ix,iy,iz)%rhodrhos(4)
			endif
		enddo
		enddo
		enddo

		if(printinfo) then
			write(*,'(24x,A,1x,i9)')    'Total # of pixels: ', n2
			write(*,'(26x,A,f8.4)')     'ratio of boundary skipped = ', dble(n1-n2)/(n1+0.0)
			write(*,'(26x,A,f8.4)')     'ratio-of-used in cell     = ', dble(n2)/dble(gb_n_cellx*gb_n_celly*gb_n_cellz)
			print *, '  (grid_rho_drho_list done)'
		endif
	end subroutine grid_rho_drho_list

  !------------------------------------------
  ! Mark the pixels that shall be dropped
  !------------------------------------------
	subroutine mark_drop_pixels(reflist, markdrop, n, i1, i2)
		!dummy
		real(dl) :: reflist(n)
		integer :: markdrop(n), n, i1, i2
		!local
		integer :: i, indexarray(n)
		if(i1.le.0 .and. i2.ge.n) then
!			print *, 'Warning (mark_drop_pixels)! No need to drop: ', i1, i2
			return
		endif
		do i = 1, n
			indexarray(i) = i
		enddo
		call Qsort2(reflist, indexarray, n)
		if(i1 .ge. 1) then
			do i = 1, i1
				markdrop(indexarray(i)) = 1
			enddo
		endif
		if(i2 .le. n) then
			do i = i2, n
				markdrop(indexarray(i)) = 1
			enddo
		endif
	end subroutine mark_drop_pixels
	

  !------------------------------------------
  ! dropping pixels in lists
  !------------------------------------------
  	subroutine drop_pixels(val_list, pos_list, dval_list, markdrop)
  		real(dl), allocatable :: val_list(:), distance_list(:), pos_list(:,:), dval_list(:,:), reflist(:), tmp(:,:)
  		integer :: markdrop(:)
  		integer :: i, j,  n, m
		n = size(val_list)
		if(size(pos_list,2).ne.n .or. size(dval_list,2).ne. n.or.size(markdrop).ne.n) then
			print *, 'ERROR (drop_pixels)! Check length of arrays: ', n,size(pos_list,2),size(dval_list,2),size(markdrop)
			stop
		endif
		m = 0
		do i = 1, n
			if(markdrop(i) .eq. 0) m = m+1
		enddo
		if(m.eq.n) return
		allocate(tmp(7,n))
		do i = 1, n
			tmp(1,i)=val_list(i); 
			tmp(2:4,i)=pos_list(1:3,i); tmp(5:7,i)=dval_list(1:3,i)	
		enddo
		deallocate(val_list,dval_list,pos_list)
		allocate(val_list(m),dval_list(3,m),pos_list(3,m))
		j = 0
		do i = 1,n
			if(markdrop(i) .eq. 0) then
				j = j+1
				val_list(j) = tmp(1,i)
				pos_list(1:3,j) = tmp(2:4,i)
				dval_list(1:3,j) = tmp(5:7,i)
			endif
		enddo
		deallocate(tmp)
	end subroutine drop_pixels
	

  !------------------------------------------
  ! dropping pixels in lists
  !------------------------------------------
  	subroutine drop_pixels2(pos_list, val_list, dval_list, markdrop, pos_list2, val_list2, dval_list2)
		! Dummy
  		real(dl), intent(in) :: pos_list(:,:), val_list(:), dval_list(:,:)
  		integer, intent(in) :: markdrop(:)
  		real(dl), intent(out), allocatable :: pos_list2(:,:), val_list2(:), dval_list2(:,:)
		! local
  		integer :: i, j,  n, m
		n = size(val_list)
		if( size(pos_list,2).ne.n .or. size(dval_list,2).ne. n.or.size(markdrop).ne.n) then
			print *, 'ERROR (drop_pixels2)! Check length of arrays: ', n,size(pos_list,2),size(dval_list,2),size(markdrop)
			stop
		endif
		m = 0
		do i = 1, n
			if(markdrop(i) .eq. 0) m = m+1
		enddo
		allocate(pos_list2(3,m),val_list2(m),dval_list2(3,m))
		j = 0
		do i = 1,n
			if(markdrop(i) .eq. 0) then
				j = j+1
				pos_list2(1:3,j) = pos_list(1:3,i)	
				val_list2(j) = val_list(i)
				dval_list2(1:3,j) = dval_list(1:3,i)	
			endif
		enddo
	end subroutine drop_pixels2


  !------------------------------------------
  ! est recommended # of points
  !------------------------------------------
  	real(dl) function est_num_in_x(smnum)
  		integer :: smnum, i
  		real(dl) :: rmin, rmax, x1,x2,y1,y2,z1,z2
  		rmin = gb_datalist(1)%r; rmax = rmin
  		x1 = gb_datalist(1)%x; x2 = x1
  		y1 = gb_datalist(1)%y; y2 = y1
  		z1 = gb_datalist(1)%z; z2 = z1
  		do i = 2, gb_numallpar
  			rmin = min(rmin, gb_datalist(i)%r)
  			rmax = max(rmax, gb_datalist(i)%r)
  			x1 = min(x1, gb_datalist(i)%x)
  			x2 = max(x2, gb_datalist(i)%x)
  			y1 = min(y1, gb_datalist(i)%y)
  			y2 = max(y2, gb_datalist(i)%y)
  			z1 = min(z1, gb_datalist(i)%z)
  			z2 = max(z2, gb_datalist(i)%z)
  		enddo
		
!		print *, x1,x2,y1,y2,z1,z2, (dble(gb_numallpar)/dble(smnum)), (x2-x1)*(y2-y1)*(z2-z2), vol_fun(rmin,rmax)
		
		est_num_in_x = ( (dble(gb_numallpar)/dble(smnum)) * (x2-x1)*(y2-y1)*(z2-z1) / vol_fun(rmin,rmax) )**(1.0_dl/3.0_dl)
  	
  	end function est_num_in_x
end module LSS_grad_fields	
