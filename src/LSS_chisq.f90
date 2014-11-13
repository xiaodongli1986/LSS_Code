
! July 1: remove gfmethod sampling...

module LSS_chisq

use LSS_grad_fields

implicit none

!!! Fixed Settings

	! access type when outputing info of mu/rho/absrho
	character(len=char_len), parameter :: gb_muinfo_access = 'append'

	! jackknife n
	integer, parameter :: gb_numjkf = 1
	integer, parameter :: gb_jkf_njkfs(gb_numjkf) = (/ 30 /)

!!! Adjustable & To be fixed Settings		
! All variables in the "Adjustable & To be fixed Settings" section must be assigned a value in the main**.f90!!!

	integer, public :: gb_num_changenuminx
	real(dl), public :: gb_amp_changenuminx

!!! Important Variables
	logical, public :: gb_chisq_initied = .false.
	integer, public :: gb_iom, gb_iw
	
contains

  !------------------------------------------
  ! Jackknife result of muav/muer
  !------------------------------------------
	subroutine muaver_jackknife(poslist, drholist, ndat, njkf, nbin, ravs, zavs, muavs, muers, printinfo)
		! Dummy
		integer, intent(in) :: ndat, njkf, nbin
		real(dl), intent(in) :: poslist(3,ndat), drholist(3,ndat)
		real(dl), intent(out) :: ravs(nbin), zavs(nbin), muavs(nbin), muers(nbin)
		logical, intent(in) :: printinfo
		! Local
		integer :: i,j,ijkf, i1,i2,i3,i4, nowndat, jkforder(ndat)
		real(dl) :: rlist(ndat), mulist(ndat), jkfquan(ndat),  nowrlist(ndat), nowmulist(ndat), &
			muavjkfs(nbin,njkf), ravjkfs(nbin,njkf), zavjkfs(nbin,njkf), tmpmuers(nbin), tmpmuvars(nbin), rers(nbin), &
			ra,dec,r, rmin,rmax
!		muaver_jackknife
		! jkfquan: this defines an order following which we split the sample
		do i = 1, ndat
			jkforder(i) = i
			call xyz_to_radecr(poslist(1,i),poslist(2,i),poslist(3,i),ra,dec,r)
			jkfquan(i) = ra
		enddo
		
		! sort rlist, mulist according the order
		call QSort2(jkfquan,jkforder,ndat)
		do i = 1, ndat
			j = jkforder(i)
			rlist(i) = rms(poslist(1:3,j),3)
			mulist(i) = get_mu_of_gradient(poslist(1,j),poslist(2,j),poslist(3,j),&
				drholist(1,j),drholist(2,j),drholist(3,j))
			mulist(i) = abs(mulist(i))
		enddo
		
		! jackknife estimation of mean/variance
		do ijkf = 1, njkf
			! partly selected data to do estimation
			i1 = floor((ijkf-1)/dble(njkf)*ndat+0.5)
			i2 = ceiling((ijkf)/dble(njkf)*ndat+0.5)
			i2 = min(i2,ndat); i1 = max(i1,0); 
			nowndat = i1 + (ndat-i2)+1
			j = 1;
			do i = 1,i1
				nowmulist(j) = mulist(i);
				nowrlist(j) = rlist(i);
				j = j+1
			enddo
			do i = i2, ndat
				nowmulist(j) = mulist(i);
				nowrlist(j) = rlist(i);
				j = j+1
			enddo
			nowndat = j-1;
			! calculate mu for this part of data;
			call find_min_max(nowrlist(1:nowndat),nowndat,rmin,rmax)
			call eqvl_binned_quan(nowmulist(1:nowndat), nowrlist(1:nowndat), rmin, rmax, nowndat, nbin, &	
				muavjkfs(1:nbin,ijkf), tmpmuers, ravjkfs(1:nbin,ijkf), tmpmuvars)
!			write(*,'(19x,A,<nbin>(f6.3,1x))') 'muavjkfs = ', muavjkfs(1:nbin,ijkf)
!			if(printinfo) then
!				print *, '(muaver_jackknife) ijkf, nowndat, muavs = ', ijkf, nowndat, muavjkfs(1:nbin,ijkf)
!			endif
		enddo
		
		! calculate mean/variance of mu/r using jackknife method
		muavs = 0.0; ravs = 0.0; muers = 0.0; rers = 0.0;
		do i = 1, nbin		
			do ijkf = 1, njkf
				muavs(i) = muavs(i)+muavjkfs(i,ijkf)
				ravs(i) = ravs(i)+ravjkfs(i,ijkf)
				muers(i) = muers(i)+muavjkfs(i,ijkf)**2.0
				rers(i) = rers(i)+ravjkfs(i,ijkf)**2.0
			enddo
			muavs(i) = muavs(i) / dble(njkf)
			muers(i) = sqrt((muers(i) - muavs(i)**2.0*dble(njkf)) * dble(njkf-1.0)/dble(njkf))
			ravs(i) = ravs(i) / dble(njkf)
			rers(i) = sqrt((rers(i) - ravs(i)**2.0*dble(njkf)) * dble(njkf-1.0)/dble(njkf))
			zavs(i) = de_zfromintpl(ravs(i))
		enddo
		if(printinfo) then
			write(*,*) '(muaver_jackknife) muavs = ', muavs(1:nbin)
			write(*,*) 'muers = ', muers(1:nbin)
			write(*,*) 'ravs  = ', ravs(1:nbin)
			write(*,*) 'rers  = ', rers(1:nbin)
		endif
	end subroutine muaver_jackknife
		
		
  !------------------------------------------
  ! estimating gradient of fields of 
  !  rho/delta/normed delta at grid points
  !------------------------------------------
	subroutine gd_mldprho_chi2s(cs, changenuminx, chisqlist, dfchisqlist, multdfchisqlist, numchisq, BSKmode)
		! Dummy
		integer, intent(in) :: numchisq
		type(chisq_settings) :: cs
		real(dl), intent(in) :: changenuminx
		real(dl), intent(out) :: chisqlist(numchisq), dfchisqlist(numchisq)
		real(dl), intent(out) :: multdfchisqlist(nbinchisq,numchisq)
		logical, optional, intent(in) :: BSKmode
		! Local
		real(dl) :: rmin,rmax
		integer :: i,j,k,i_split,nownbin, idrop, n1,n2, num_bins, num, num_points
		real(dl), allocatable :: & !gb_pos_list(:,:), gb_rho_list(:), gb_drho_list(:,:), 
			drho_mu_data(:), drho_mu_data_orig(:), &
			pos_list(:,:), rho_list(:), drho_list(:,:),  absdrholist(:), reflist(:)
		integer, allocatable :: markdrop(:)
		logical :: print_time = .false.
		real(dl) :: tmp, time0, time1, time2, time3, time4, time5, time6, time7
		!variables used to estimate mu at different r bins
		real(dl), allocatable :: r_list(:), muavlist(:), muerlist(:), ravlist(:), muvarlist(:), &
			muavlist_test(:), muerlist_test(:), ravlist_test(:), muvarlist_test(:), zlist_test(:), &
			absdrhoavlist(:), absdrhoerlist(:), absdrhovarlist(:), &
			absdrhoavlist_test(:), absdrhoerlist_test(:),  absdrhovarlist_test(:), &
			absrhoavlist_test(:),  absrhoerlist_test(:),   absrhovarlist_test(:)
		!tmp variables used for outputing info of gradient field
		real(dl) :: xyz1,xyz2,r1,r2,tsbedropratio,nowx,nowy,nowz,nowrsq,rmaxsq,invrsq,rhomean, rhovar, tmpvar, chisq_div_frac
		character(len=char_len) :: str, str1, tmpstr1, tsbefile, file1, file2, mufile1, mufile2, mufile3, &
			absdrhofile1, absdrhofile2, absdrhofile3, kwstr, muinfo_access
	
		if(numchisq.ne.cs%numdrop) then
			print *, 'ERROR (gd_mldprho_chi2s)! size of chisqlist mismathces than ', cs%numdrop
			stop
		endif

		call cpu_time(time0)
		
		if(cs%print_info) then
			write(*,'(A)') '   (gd_mldprho_chi2s) Estimating chisq: method = \"grid\" '
			write(*,'(A,f9.4,A,i3)') '   (gd_mldprho_chi2s) changenuminx = ', real(changenuminx), &
				'; # of drop = ', cs%numdrop
			write(*,'(22x,A,i2,A)') 'gradient field bin-#  (', nbinchisq, ' options): '
			write(*,'(22x,$)') 
			write(*,*) nbinchisqlist
			write(*,'(22x,A,i2,A)') '    jackknife  bin-#  (', gb_numjkf, ' options): '
			write(*,'(22x,$)') 
			write(*,*)  gb_jkf_njkfs
		endif

		muinfo_access = gb_muinfo_access
		if(present(BSKmode)) then
			if(BSKmode) then 
				muinfo_access = 'sequential'
				goto 100
			endif
		endif

		call init_mult_lists(cs%print_info)
		call do_cell_init(real(cs%num_in_x)*(1.0_dl+changenuminx), cs%print_info)		
		
		call grid_rho_drho_list(cs%smnum, cs%print_info, gb_pos_list, gb_rho_list, gb_drho_list)
100		continue	

		chisq_div_frac = 1.0 ! do not consider chisq division for grid_rho_drho_list
		call cpu_time(time1)	

		call get_mean_var(gb_rho_list, rhomean, rhovar)
		n1 = size(gb_rho_list)
!		do i = 1, n1
!			gb_rho_list(i) = gb_rho_list(i) - 1
!		enddo
		if(cs%print_info) then
			write(*,'(A,i8,2(1x,e14.7,1x))') '   (gd_mldprho_chi2s) # / MEAN / VAR of deltas = ', &
				 n1, real(rhomean), real(rhovar)
		endif
		
!		print *, '(gd_mldprho_chi2s) Test C...'
		allocate(absdrholist(n1),markdrop(n1),reflist(n1))
		markdrop = 0
		do i = 1, n1
			absdrholist(i)=sqrt(gb_drho_list(1,i)**2.0+gb_drho_list(2,i)**2.0+gb_drho_list(3,i)**2.0)
		enddo

		num_bins = 2
		allocate(muavlist(num_bins), muerlist(num_bins), ravlist(num_bins), muvarlist(num_bins))

		if(cs%print_info) write(*,'(A)') '   (gd_mldprho_chi2s) Estimating multi-dropped chisqs...'
		do idrop = 1, cs%numdrop
			! cycle if no drop
			markdrop = 0
			if(cs%dropval(idrop)) then
				reflist = gb_rho_list
				call mark_drop_pixels(reflist, markdrop, n1,&
					int(cs%lowdropvalratio(idrop)*n1)-1, int(n1*( 1.0_dl-cs%highdropvalratio(idrop) )) + 1)
			endif
			if(cs%dropdval(idrop)) then
				reflist = absdrholist
				call mark_drop_pixels(reflist, markdrop, n1,&
					int(cs%lowdropdvalratio(idrop)*n1)-1, int((1-cs%highdropdvalratio(idrop))*n1)+1)
			endif
			
			n2 = 0
			do i = 1, n1
				if(markdrop(i).eq.0) n2 = n2+1
			enddo
			if(cs%print_info) then
				write(*,'(20x,A,f6.2,1x,f6.2,A,f6.2,1x,f6.2,A,i9,f6.2)') &
					'  Drop rho (low/high)  = ', &
					cs%lowdropvalratio(idrop), cs%highdropvalratio(idrop), &
					'; Drop drho = ', &
					cs%lowdropdvalratio(idrop), cs%highdropdvalratio(idrop), &
					'; Left_# / drop_ratio = ', n2, (n1-n2)/dble(n1)
			endif

			call drop_pixels2(gb_pos_list, gb_rho_list, gb_drho_list, markdrop, pos_list, rho_list, drho_list)
			call get_mu_from_gradient_list(pos_list, drho_list, drho_mu_data_orig)
			num_points = size(rho_list)
			allocate(r_list(num_points))
			do i = 1, num_points
				r_list(i) = sqrt(pos_list(1,i)**2.0+pos_list(2,i)**2.0+pos_list(3,i)**2.0)
			enddo
			call find_min_max(r_list, size(r_list), rmin, rmax)
			
			! chisq of mu
			chisqlist(idrop) = chisq_of_mu_data2(drho_mu_data_orig) / chisq_div_frac
			! abs mu, drho
			if(allocated(drho_mu_data)) deallocate(drho_mu_data)
			allocate(drho_mu_data(num_points))
			do i = 1, num_points
				drho_mu_data(i) = abs(drho_mu_data_orig(i))
				if(rho_list(i)<0.0) then
					print *, 'WARNING (gd_mldprho_chi2s)!!! Encountered negative rho: ', i, rho_list(i)
				endif
			enddo

			!############################
			! chisq from redshift dependent mu...
		  	call eqvl_binned_quan(drho_mu_data, r_list, rmin, rmax, num_points, num_bins, &	
				muavlist, muerlist, ravlist, muvarlist)
!			if(cs%print_info) & write(*,'(22x,A,2e14.7)') '  <mu> at two vlumes:            ', &
!						muavlist(1), muavlist(2)
			dfchisqlist(idrop) = abs(muavlist(2) - muavlist(1))**2.0 &
					/ ((muerlist(1)**2.0 + (muerlist(2)**2.0)))  / chisq_div_frac
			call dfchisqs(nbinchisq, nbinchisqlist, drho_mu_data, r_list, n2, multdfchisqlist(1:nbinchisq,idrop))
			
			do i = 1, nbinchisq
				multdfchisqlist(i,idrop) = multdfchisqlist(i,idrop) / chisq_div_frac
			enddo 

			if(gb_outputinfo_mu) then ! only output mu for idrop 1
!				if(cs%print_info.and.idrop.eq.1) then
!					write(*,'(21x,A)') 'output information: _muinfo_idrop?.txt; _muinfo_idrop?.noabs.txt; '
!				endif
	                        write(str,*) idrop
	                        
	                        kwstr = '_muinfo_idrop'
	                        if(present(BSKmode)) then
	                        	if(BSKmode) kwstr = '_BSKmuinfo_idrop'
	                        endif
	                        
                                mufile1 = trim(adjustl(gb_suboutput_name))//trim(adjustl(kwstr))//trim(adjustl(str))//'.txt'
				mufile2 = trim(adjustl(gb_suboutput_name))//trim(adjustl(kwstr))//trim(adjustl(str))//'.noabs.txt'

				open(unit=898,file=mufile1,access=muinfo_access)
				open(unit=899,file=mufile2,access=muinfo_access)

				if(cs%print_info.and.idrop.eq.1) then
					write(*,'(A)') '  output info of mu: '//trim(adjustl(mufile1))
					write(*,'(A)') '                     '//trim(adjustl(mufile2))
				endif
				write(898,'(2i4,2(e14.7,1x),A)') gb_iom, gb_iw, gb_omegam, gb_w, ' # i_om, i_w, omegam, w'
				write(899,'(2i4,2(e14.7,1x),A)') gb_iom, gb_iw, gb_omegam, gb_w, ' # i_om, i_w, omegam, w'

				! output jkf mus...
				do i = 1, gb_numjkf
 	                            write(str,*) idrop
				    write(tmpstr1,*) gb_jkf_njkfs(i)
				    mufile3 = trim(adjustl(gb_suboutput_name))//trim(adjustl(kwstr))//trim(adjustl(str))//'.jkf' &
						//trim(adjustl(tmpstr1))//'.txt'
				    open(unit=900,file=mufile3,access=muinfo_access)
				    if(cs%print_info.and.idrop.eq.1) then
					write(*,'(A)') '                     '//trim(adjustl(mufile3))
				    endif
				    write(900,'(2i4,2(e14.7,1x),A)') gb_iom, gb_iw, gb_omegam, gb_w, ' # i_om, i_w, omegam, w'
				    do i_split = 1, nbinchisq
					nownbin = nbinchisqlist(i_split)
					write(900,*) nownbin, ' # num-of-bin'
					allocate(muavlist_test(nownbin), muerlist_test(nownbin),zlist_test(nownbin),&
						 ravlist_test(nownbin), muvarlist_test(nownbin)) ! for tests
					! mu calculated from jackknife
					call muaver_jackknife(pos_list, drho_list, num_points, gb_jkf_njkfs(i), nownbin, &
							ravlist_test, zlist_test, muavlist_test, muerlist_test, .false.)
					write(900,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,zlist_test(1:nownbin)))), ' # redshift'
					write(900,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,ravlist_test(1:nownbin)))), ' # r'
					write(900,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,muavlist_test(1:nownbin)))), ' # \bar mu'
					write(900,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,muerlist_test(1:nownbin)))), ' # er of \bar mu'
					deallocate(muavlist_test, muerlist_test,zlist_test, ravlist_test, muvarlist_test)
				    enddo
				    close(900)
				enddo
				
				do i_split = 1, nbinchisq
					nownbin = nbinchisqlist(i_split)
					write(898,*) nownbin, ' # num-of-bin'
					write(899,*) nownbin, ' # num-of-bin'

					allocate(muavlist_test(nownbin), muerlist_test(nownbin),zlist_test(nownbin),&
						 ravlist_test(nownbin), muvarlist_test(nownbin)) ! for tests
					! mu usually used to do chisq calculating (after taking tha absolute)
					call eqvl_binned_quan(drho_mu_data, r_list, rmin, rmax, size(drho_mu_data), &
						nownbin, muavlist_test, muerlist_test, ravlist_test, muvarlist_test)
					do k = 1, nownbin
						zlist_test(k) = de_zfromintpl(ravlist_test(k))
					enddo
					write(898,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,zlist_test(1:nownbin)))), ' # redshift'
					write(898,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,ravlist_test(1:nownbin)))), ' # r'
					write(898,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,muavlist_test(1:nownbin)))), ' # \bar mu'
					write(898,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,&
					  muerlist_test(1:nownbin)* sqrt(chisq_div_frac)))), ' # er of \bar mu'

					! "most original" mu (without taking tha absolute)
					call eqvl_binned_quan(drho_mu_data_orig, r_list, rmin, rmax, size(drho_mu_data), &
						nownbin, muavlist_test, muerlist_test, ravlist_test, muvarlist_test)
					do k = 1, nownbin
						zlist_test(k) = de_zfromintpl(ravlist_test(k))
					enddo
					write(899,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,zlist_test(1:nownbin)))), ' # redshift'
					write(899,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,ravlist_test(1:nownbin)))), ' # r'
					write(899,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,muavlist_test(1:nownbin)))), ' # \bar mu'
					write(899,'(A,A)') &
					 trim(adjustl(WriteFmtEfloat(nownbin,&
					  muerlist_test(1:nownbin)* sqrt(chisq_div_frac)))), ' # er of \bar mu'
					deallocate(muavlist_test, muerlist_test,zlist_test, ravlist_test, muvarlist_test)
				enddo
				close(898); close(899);
			endif

200			continue
			deallocate(r_list,pos_list,rho_list,drho_list,drho_mu_data,drho_mu_data_orig)
		enddo
		deallocate(absdrholist,markdrop,reflist)
		deallocate(muavlist,muerlist,ravlist,muvarlist)
		call cpu_time(time2)		

		if(print_time .or. cs%print_info) then
			write(*,'(A,f12.5)') '   (gd_mldprho_chi2s) Time used in grid_rho_drho_list: ', real(time1-time0)
			write(*,'(A,f12.5)') '   (gd_mldprho_chi2s) Time used in drop: ', real(time2-time1)
		endif
		if(cs%print_info) write(*,'(A)') '   (gd_mldprho_chi2s) Done.'
	end subroutine gd_mldprho_chi2s


  !------------------------------------------
  ! calculating chisq
  !------------------------------------------
	  subroutine gf_mldprho_chi2s(omegam, w, h, cs, chisqlist, dfchisqlist, multdfchisqlist, numchisq, calc_comvr)
		! dummy
		real(dl), intent(in) :: omegam, w, h
		integer, intent(in) :: numchisq
		type(chisq_settings), intent(inout) :: cs
		real(dl), intent(out) :: chisqlist(numchisq), dfchisqlist(numchisq), multdfchisqlist(nbinchisq,numchisq)
		logical, intent(in) :: calc_comvr 
		! local
		integer :: i, j, k, numchanges, num
		real(dl), allocatable :: changenuminxlist(:)
		real(dl) :: tmpchisqlist(numchisq), tmpdfchisqlist(numchisq), tmpmultdfchisqlist(nbinchisq,numchisq)
		! variables used in r-coorection
		real(dl) :: numinx_float, x,y,z
		logical :: touchbdflag,  outputinfo = .false. !!CHECK
		character(len=char_len) :: str1, file1

		if(cs%print_info) then
			write(*,'(A)') '   (gf_mldprho_chi2s) Averaged-chisq estimated from multiple adjustments of grid: '
			print *, '                      # of chisqs   = ',  2*gb_num_changenuminx+1
			if(gb_num_changenuminx.eq.0) then
				print *, '                      maximal amp change = ', 0
			else
				print *, '                      maximal amp change = ', real(gb_amp_changenuminx)
			endif
		endif

		if (numchisq.ne.cs%numdrop) then
			print *, 'ERROR (gf_mldprho_chi2s)! size of chisqlist mismatches with ', cs%numdrop
			stop
		endif

		if (.not. gb_chisq_initied) then
			call cosmo_funs_init(cs%print_info)
			call readin_dataran(cs%print_info)
			gb_chisq_initied = .true.
		endif

		! take the first cosmology as the input cosmology
		call init_cosmo(omegam,w,h,cs%print_info)

		!get the chisqs
		numchanges = 2*gb_num_changenuminx+1
		allocate(changenuminxlist(numchanges))
		if(cs%print_info) write(*,'(A)') '   (gf_mldprho_chi2s) Estimating multiple chisqs from different adjustments...'
		changenuminxlist(1) = 0.0_dl
		do i = 1, gb_num_changenuminx
			changenuminxlist(i+1) = gb_amp_changenuminx*i/(gb_num_changenuminx+0.0)
			changenuminxlist(2*gb_num_changenuminx-i+2) = -1.0_dl *  gb_amp_changenuminx*i/(gb_num_changenuminx+0.0)
		enddo
		chisqlist = 0.0_dl
		dfchisqlist = 0.0_dl
		multdfchisqlist = 0.0_dl
		do i = 1, numchanges
			if(cs%print_info) write(*,'(A,f7.3,A)') '   (gf_mldprho_chi2s) chisqs with adjustment ratio ', &
				changenuminxlist(i),':'
			tmpchisqlist = 0.0
			tmpdfchisqlist = 0.0
			tmpmultdfchisqlist = 0.0
			call gd_mldprho_chi2s(cs, changenuminxlist(i),tmpchisqlist,tmpdfchisqlist,tmpmultdfchisqlist,cs%numdrop)
			!print *, 'Calling gd_mldprho_chi2s done.'
			do j = 1, cs%numdrop
				chisqlist(j) = chisqlist(j) + tmpchisqlist(j)/dble(numchanges)
				dfchisqlist(j) = dfchisqlist(j) + tmpdfchisqlist(j)/dble(numchanges)
				do k = 1, nbinchisq
					multdfchisqlist(k,j) = multdfchisqlist(k,j) + tmpmultdfchisqlist(k,j)/dble(numchanges)
				enddo
			enddo
		enddo
		deallocate(changenuminxlist)
		
		! clean up ...
                if(allocated(gb_pos_list)) deallocate(gb_pos_list)
                if(allocated(gb_rho_list)) deallocate(gb_rho_list)
                if(allocated(gb_drho_list)) deallocate(gb_drho_list)
                
		if(cs%print_info) then
			write(*,'(A)') '   (gf_mldprho_chi2s) Done.'
		endif
	end subroutine gf_mldprho_chi2s
end module LSS_chisq
