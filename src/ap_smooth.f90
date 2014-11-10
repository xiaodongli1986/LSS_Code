
! Now we are confident that our NNB and NNRandom search subroutines shall be right (we compare the results with direct search!!!)

!####################################
!This module does smooth
!####################################
module ap_smooth
use ap_cell

	implicit none

!!! Fixed Settings
	! Over-searching ratio for nearest neighbors.
	real(dl), parameter :: gb_NNBosrat = 2.5d0 ! Searching for *gb_NNBosrat numbers of nearest neighbors lying within a nearest box 
						   !  Oversearch to ensure we will find enough number in nearest sphere!
	! Skip graident vectors estimated with nearest # <= gb_min_smnum: only works for fixed_sm
	integer, parameter :: gb_min_smnum = 0 
	logical, parameter :: gb_keepzerorho = .true. ! do not skip pixels with gb_num_allpar <= gb_min_smnum but keep them (set rho as 0)

	! maximal # of neary points (related with array length)
	integer, parameter :: max_insphere_num = 200000
	
	! check the validity of NNB (sometimes there is error due to the small gb_NNBosrat)
	logical, parameter :: gb_do_NNBReSearch = .true.
	real(dl), parameter :: gb_maxNNBosrat = 100.0_dl ! maximal Osrat allowd (avoid endless cycling...)
	
!!! Important Variables
	! Number count of how many NNB re-search happens
	integer :: gb_num_NNBSearch, gb_num_NNBReSearch, gb_num_NNBReSearchEr
	
contains

  !------------------------------------------
  ! N nearest neighbours search
  ! Select out data in the cells nearby the
  !  given position. 
  ! Set gb_NNBosrat >>1 to make sure the number is 
  !  far larger than required
  !------------------------------------------
  	subroutine NNBSearch_data(x,y,z,num,given_NNBosrat,selected_list,touchbdflag,bddist) 
		! Dummy
  		real(dl), intent(in) :: x,y,z
  		real(dl), optional, intent(in) :: given_NNBosrat
  		integer, intent(in) :: num
  		integer, allocatable, intent(out) :: selected_list(:)
  		real(dl), optional, intent(out) :: bddist
  		logical, optional, intent(out) :: touchbdflag
  		! Local
  		integer :: i,j,k,l1,l2,ix,iy,iz,requirednum
  		integer :: di, i1,i2, j1,j2, k1,k2, now_num
  		real(dl) :: dev, NNBosrat, ratio, xl,xr,yl,yr,zl,zr
  		
  		ix = max(min(int((x-gbgridxmin)/gbdeltax)+1,gb_n_cellx),1)
  		iy = max(min(int((y-gbgridymin)/gbdeltay)+1,gb_n_celly),1)
  		iz = max(min(int((z-gbgridzmin)/gbdeltaz)+1,gb_n_cellz),1)
  		dev = max(abs(abs(gbgridxmin + (ix-1)*gbdeltax - x)/gbdeltax - 0.5), &
  			abs(abs(gbgridymin + (iy-1)*gbdeltay - y)/gbdeltay - 0.5), &
  			abs(abs(gbgridzmin + (iz-1)*gbdeltaz - z)/gbdeltaz - 0.5) )
  		if(present(given_NNBosrat)) then
  			NNBosrat = given_NNBosrat
  		else
  			NNBosrat = gb_NNBosrat
  		endif
  		ratio = NNBosrat * (1.0_dl)! + dev * 1.5)
  		requirednum = ceiling(num*ratio) ! at least 40 data!!!
		di = 0
		do while(1.eq.1)
			call ilist(ix,di,1,gb_n_cellx,i1,i2)
			call ilist(iy,di,1,gb_n_celly,j1,j2)
			call ilist(iz,di,1,gb_n_cellz,k1,k2)
			now_num = 0
			do i = i1, i2
			do j = j1, j2
			do k = k1, k2
				now_num = now_num + gb_cell_mat(i,j,k)%numdata
			enddo
			enddo
			enddo
			if(now_num .ge. requirednum) exit
			di = di + 1
		enddo
		!i1,i2,j1,j2,k1,k2 are the index range of cells surround the point and containing enough points...
		allocate(selected_list(now_num))
		l1 = 1
		do i = i1, i2
		do j = j1, j2
		do k = k1, k2
			do l2 = 1, gb_cell_mat(i,j,k)%numdata
				selected_list(l1) = gb_cell_mat(i,j,k)%idatalist(l2)
				l1 = l1 + 1
			enddo
		enddo
		enddo
		enddo
		if(present(touchbdflag)) then
			if(i1.eq.1.or.i2.eq.gb_n_cellx.or.j1.eq.1.or.j2.eq.gb_n_celly.or.k1.eq.1.or.k2.eq.gb_n_cellz) then
				touchbdflag = .true.
			else
				touchbdflag = .false.
			endif
		endif
		if(present(bddist)) then
			xl = gbgridxmin + (i1-1.0_dl)*gbdeltax
			xr = gbgridxmin + (i2+1.0_dl)*gbdeltax
			yl = gbgridymin + (j1-1.0_dl)*gbdeltay
			yr = gbgridymin + (j2+1.0_dl)*gbdeltay
			zl = gbgridzmin + (k1-1.0_dl)*gbdeltaz
			zr = gbgridzmin + (k2+1.0_dl)*gbdeltaz
			bddist = min(abs(x-xl),abs(x-xr),abs(y-yl),abs(y-yr),abs(z-zl),abs(z-zr))
		endif
	end subroutine NNBSearch_data

  !------------------------------------------
  ! N nearest neighbours search (exact)
  !------------------------------------------  
	subroutine NNBExactSearch_data(x,y,z,insphere_num,selected_list,erflag) 	
!  	subroutine nb_list0(x,y,z, insphere_num, rho,drhodx,drhody,drhodz, max_dist, erflag)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z ! position where you want to calculate rho and drho
  		integer, intent(in) :: insphere_num ! # of nearest neighbors you want to use
  		logical, intent(out) :: erflag ! whether there is error in the result
  		integer, intent(out) :: selected_list(insphere_num) ! array of index of nearby data points selected 
  		! LOCAL VARIABLES
  		logical :: touchbdflag
  		real(dl) :: bddist, ratrat, h, r0(3), r, max_dist
  		real(dl), allocatable :: distance_array(:), tmp(:), smda(:)
  		integer :: i,n,allparindex
  		integer, allocatable :: tmpindex(:), smlablist(:), index_array(:)

		ratrat = 1.0_dl ! given_NNBosrat=gb_NNBosrat*ratrat; increase if touch boundary happens frequently
		do while(.true.)
			if(allocated(index_array)) deallocate(index_array)
			if(allocated(tmpindex)) deallocate(tmpindex)
			if(allocated(distance_array)) deallocate(distance_array)
			if(allocated(smlablist)) deallocate(smlablist)
			if(allocated(smda)) deallocate(smda)

			! Near x,y,z, we require insphere_num of NNB, then get a rough list of n (much larger then insphere_num) NBs
	  		call NNBSearch_data(x,y,z,insphere_num,gb_NNBosrat*ratrat, index_array, touchbdflag, bddist)
	  		n=size(index_array)

			! Calculate the distances to these NBs;
			r0(1)=x; r0(2)=y; r0(3)=z;
	  		allocate(tmpindex(n), distance_array(n))
			do i = 1, n
				distance_array(i) = distance(gb_xyz_list(1:3,index_array(i)),r0,3)
				tmpindex(i) = index_array(i)
			enddo
			! From these NBs select out the insphere_num NNBs; 
			allocate(smlablist(insphere_num),smda(insphere_num))
			call ltlablist2(distance_array,n,insphere_num,smda,smlablist)
	  		
	  		max_dist = maxval(smda)
	  		if(.not.touchbdflag.and.(max_dist > bddist)) then
!	  			write(*,'(A,f10.3,i7,f10.3,f10.1)') 'TESTING (nb_list0): Found max_dist > bddist:'//&
!	  				' max_dist, insphere_num, bddist, ratrat = ', max_dist, insphere_num, bddist, ratrat
	  			ratrat = ratrat*1.2+1.0_dl
	  			if(ratrat>gb_maxNNBosrat) then
!	  				print *, 'WARNING (nb_list0): Very large NNBosrat encountered:', ratrat, '; skip with error.'
	  				erflag = .true.
	  				return
	  			endif
	  		else
	  			exit
	  		endif
		enddo
		do i = 1, insphere_num
			selected_list(i) = tmpindex(smlablist(i))
		enddo
	end subroutine NNBExactSearch_data
		

  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel
  !------------------------------------------
  	subroutine nb_list0(x,y,z, insphere_num, rho,drhodx,drhody,drhodz, max_dist, erflag, printinfo)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z ! position where you want to calculate rho and drho
  		integer, intent(in) :: insphere_num ! # of nearest neighbors you want to use
  		real(dl), intent(out) :: rho,drhodx,drhody,drhodz, max_dist ! results; max_dist is the radius of the sphere
  		logical, intent(out) :: erflag ! whether there is error in the result
  		logical, optional, intent(in) :: printinfo
  		! LOCAL VARIABLES
  		integer, allocatable :: index_array(:) ! array of index of nearby data points selected 
  		real(dl), allocatable :: distance_array(:), xyz_mass_array(:,:), tmp(:), smda(:)
  		integer, allocatable :: tmpindex(:), smlablist(:)
  		integer :: i,n,allparindex
  		real(dl) :: ratrat,bddist, h, r0(3), r, mass, dweight
  		logical :: touchbdflag

		erflag = .false.
		if(present(printinfo)) then
			print *, ' (nb_list0) Start'
		endif

		ratrat = 1.0_dl ! given_NNBosrat=gb_NNBosrat*ratrat; increase if touch boundary happens frequently
		do while(.true.)
			if(allocated(index_array)) deallocate(index_array)
			if(allocated(tmpindex)) deallocate(tmpindex)
			if(allocated(distance_array)) deallocate(distance_array)
			if(allocated(smlablist)) deallocate(smlablist)
			if(allocated(smda)) deallocate(smda)
			if(allocated(xyz_mass_array)) deallocate(xyz_mass_array)

			if(present(printinfo)) then
				print *, ' (nb_list0) Before NNBSearch_data'
			endif

			! Near x,y,z, we require insphere_num of NNB, then get a rough list of n (much larger then insphere_num) NBs
	  		call NNBSearch_data(x,y,z,insphere_num,gb_NNBosrat*ratrat, index_array, touchbdflag, bddist)
	  		n=size(index_array)

			if(present(printinfo)) then
				print *, ' (nb_list0) n = ', n
			endif

			! Calculate the distances to these NBs;
			r0(1)=x; r0(2)=y; r0(3)=z;
	  		allocate(tmpindex(n), distance_array(n))
			do i = 1, n
				distance_array(i) = distance(gb_xyz_list(1:3,index_array(i)),r0,3)
				tmpindex(i) = index_array(i)
			enddo

			if(present(printinfo)) then
				print *, ' (nb_list0) Before ltlablist2'
			endif

			! From these NBs select out the insphere_num NNBs; 
			allocate(smlablist(insphere_num),smda(insphere_num))
			call ltlablist2(distance_array,n,insphere_num,smda,smlablist)
			if(present(printinfo)) then
				print *, ' (nb_list0) After ltlablist2'
			endif
			! Get the x,y,z,mass of these NNBs
			allocate(xyz_mass_array(4,insphere_num))
			do i = 1, insphere_num
				allparindex = tmpindex(smlablist(i)) ! this headache code get the index of the point in the gb_***_list...
				xyz_mass_array(1:3,i) = gb_xyz_list(1:3,allparindex)
				xyz_mass_array(4,i) = gb_mass_list(allparindex)
			enddo
	  		max_dist = maxval(smda)
	  		gb_num_NNBSearch = gb_num_NNBSearch+1
	  		if(gb_do_NNBReSearch.and.(.not.touchbdflag.and.(max_dist > bddist))) then
!	  			write(*,'(A,f10.3,i7,f10.3,f10.1)') 'TESTING (nb_list0): Found max_dist > bddist:'//&
!	  				' max_dist, insphere_num, bddist, ratrat = ', max_dist, insphere_num, bddist, ratrat
	  			ratrat = ratrat*1.2+1.0_dl
	  			gb_num_NNBReSearch = gb_num_NNBReSearch+1
	  			if(ratrat>gb_maxNNBosrat) then
		  			gb_num_NNBReSearchEr = gb_num_NNBReSearchEr+1
!	  				print *, 'WARNING (nb_list0): Very large NNBosrat encountered:', ratrat, '; skip with error.'
	  				erflag = .true.
	  				return
	  			endif
	  		else
	  			exit
	  		endif
		enddo

		if(present(printinfo)) then
			print *, ' (nb_list0) Before final computation of rho, drhos'
		endif

  		h = max_dist / 2.0
		rho = 0; drhodx=0;drhody=0;drhodz=0;
		do i = 1, insphere_num
			r = distance(xyz_mass_array(1:3,i),r0,3) !distance_array(smlablist(i))
			mass = xyz_mass_array(4,i)
			rho = rho + mass*w_kernel(r, h)
			dweight = der_w_kernel(r,h)
			drhodx = drhodx + mass*(x-xyz_mass_array(1,i)) / r * dweight
			drhody = drhody + mass*(y-xyz_mass_array(2,i)) / r * dweight
			drhodz = drhodz + mass*(z-xyz_mass_array(3,i)) / r * dweight
		enddo
		if(present(printinfo)) then
			print *, ' (nb_list0) End'
		endif
	end subroutine nb_list0			


  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel; fixed radius
  ! TBCheck... !!! Not used...
  !------------------------------------------
  	subroutine nb_fixmd_list(x,y,z,num,rho,drhodx,drhody,drhodz,fixmd, erflag)!Testing
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z, fixmd
  		integer, intent(out) :: num
  		real(dl), intent(out) :: rho,drhodx,drhody,drhodz
  		logical, intent(out) :: erflag ! whether there is error in the result
  		! LOCAL VARIABLES
  		integer, parameter :: max_num = 100000 !maximal # of neary halos 
  		real(dl) :: distance_array(max_num), xyz_mass_array(4,max_num)
  		integer :: i,j,k,l, imin,imax,jmin,jmax,kmin,kmax, nownum,allparindex
  		real(dl) :: r0(3), h, nowr, mass, dweight

  		imin = int((x-fixmd-gbgridxmin)/gbdeltax +1.0_dl)
  		imax = int((x+fixmd-gbgridxmin)/gbdeltax +1.0_dl)
  		jmin = int((y-fixmd-gbgridymin)/gbdeltay +1.0_dl)
  		jmax = int((y+fixmd-gbgridymin)/gbdeltay +1.0_dl)
  		kmin = int((z-fixmd-gbgridzmin)/gbdeltaz +1.0_dl)
  		kmax = int((z+fixmd-gbgridzmin)/gbdeltaz +1.0_dl)
  		
		num = 0
		r0(1)=x; r0(2)=y; r0(3)=z;

		if(imin<1.or.imin>gb_n_cellx.or.jmin<1.or.jmax>gb_n_celly.or.kmin<1.or.kmax>gb_n_cellz) then
			erflag = .true. 
			return
		endif

		do i = max(1,imin), min(gb_n_cellx,imax)
		do j = max(1,jmin), min(gb_n_celly,jmax)
		do k = max(1,kmin), min(gb_n_cellz,kmax)
			do l = 1, gb_cell_mat(i,j,k)%numdata
				allparindex = gb_cell_mat(i,j,k)%idatalist(l)
				nowr = distance(gb_xyz_list(1:3,allparindex),r0,3)
				if(nowr < fixmd) then
					num = num+1
					distance_array(num) = nowr
					xyz_mass_array(1:3, num) = gb_xyz_list(1:3,allparindex)
					xyz_mass_array(4, num) = gb_mass_list(allparindex)
				endif
			enddo
		enddo
		enddo
		enddo

		if(num .eq. 0) then
			erflag = .true.
			return
		endif

		if(num > max_num) then
			print *, 'ERROR (nb_fixmd_list): # of data overflow: ', num, max_num
			erflag = .true. 
			return
		endif

		rho = 0; drhodx=0;drhody=0;drhodz=0;		
		if(num .le. gb_min_smnum) then
			erflag = .true. 
			return
		endif			
			
  		h = fixmd / 2.0
		do i = 1, num
			nowr = distance_array(i)
			mass = xyz_mass_array(4,i)
			rho = rho + mass*w_kernel(nowr, h)
			dweight = der_w_kernel(nowr,h)
			drhodx = drhodx + mass*(x-xyz_mass_array(1,i)) / nowr * dweight
			drhody = drhody + mass*(y-xyz_mass_array(2,i)) / nowr * dweight
			drhodz = drhodz + mass*(z-xyz_mass_array(3,i)) / nowr * dweight
		enddo
		erflag = .false.
	end subroutine nb_fixmd_list

  !------------------------------------------
  ! calculating weirat based on randoms
  !------------------------------------------
  	subroutine nb_fixmd_ranwei(x,y,z,fixmd,insphere_numran,weirat,erflag,opt_sumwei_acpt,opt_sumwei_all)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z, fixmd 
  		integer, intent(out) :: insphere_numran
  		real(dl), intent(out) :: weirat ! this mark the boundary effect: ratio of wei, randoms-inside-boundary .vs. all-randoms
  		integer, intent(out):: erflag ! 0: noerror; 1: out of grid; 2: no random found
  		real(dl), optional, intent(out) :: opt_sumwei_acpt,opt_sumwei_all
  		! LOCAL VARIABLES
  		real(dl) :: distance_array(max_insphere_num), xyz_mass_array(4,max_insphere_num)
  		logical :: acptarray(max_insphere_num)
  		integer :: i,j,k,l, imin,imax,jmin,jmax,kmin,kmax, nowinsphere_num, allparindex,ranindex
  		real(dl) :: r0(3), h, nowr, mass, wei,sumwei_acpt, sumwei_all

		weirat = 0.0
  		imin = int((x-fixmd-gbgridxmin)/gbdeltax +1.0_dl)
  		imax = int((x+fixmd-gbgridxmin)/gbdeltax +1.0_dl)
  		jmin = int((y-fixmd-gbgridymin)/gbdeltay +1.0_dl)
  		jmax = int((y+fixmd-gbgridymin)/gbdeltay +1.0_dl)
  		kmin = int((z-fixmd-gbgridzmin)/gbdeltaz +1.0_dl)
  		kmax = int((z+fixmd-gbgridzmin)/gbdeltaz +1.0_dl)
  		
  		erflag = 0
		if(imin<1.or.imin>gb_n_cellx.or.jmin<1.or.jmax>gb_n_celly.or.kmin<1.or.kmax>gb_n_cellz) then
			erflag = 1 ! 1 means out of grid range
		else
			erflag = .false.
		endif

		insphere_numran = 0
		r0(1)=x; r0(2)=y; r0(3)=z;

		do i = max(1,imin), min(gb_n_cellx,imax)
		do j = max(1,jmin), min(gb_n_celly,jmax)
		do k = max(1,kmin), min(gb_n_cellz,kmax)
			do l = 1, gb_cell_mat(i,j,k)%numran
				allparindex = gb_cell_mat(i,j,k)%iranlist(l)
				nowr = distance(gb_xyz_list(1:3,allparindex),r0,3)
				if(nowr < fixmd) then
					insphere_numran = insphere_numran+1
					distance_array(insphere_numran) = nowr
					xyz_mass_array(1:3,insphere_numran) = gb_xyz_list(1:3,allparindex)
					xyz_mass_array(4,insphere_numran) = gb_mass_list(allparindex)
					ranindex = allparindex - gb_numdata
					acptarray(insphere_numran) = gb_ranlist(ranindex)%acpt
				endif
			enddo
		enddo
		enddo
		enddo

		if(insphere_numran > max_insphere_num) then
			print *, 'WARNING (nb_fixmd_ranwei): # of random overflow: ', insphere_numran, max_insphere_num
			erflag = 3
			return
		endif

		if(insphere_numran .eq. 0) then
			erflag = 2 ! 2 means no random found!
			return
		endif

		sumwei_acpt = 0.0; sumwei_all = 0.0
  		h = fixmd / 2.0
		do i = 1, insphere_numran
			nowr = distance_array(i)
			mass = xyz_mass_array(4,i)
			wei = mass*w_kernel(nowr, h)
			sumwei_all = sumwei_all + wei
			if(acptarray(i)) then
				sumwei_acpt = sumwei_acpt + wei
			endif
		enddo
		weirat = sumwei_acpt / sumwei_all
		if(present(opt_sumwei_acpt)) then
			opt_sumwei_acpt = sumwei_acpt
		endif
		if(present(opt_sumwei_all)) then
			opt_sumwei_all  = sumwei_all
		endif
	end subroutine nb_fixmd_ranwei

  !------------------------------------------
  ! calculating weirat based on ra/dec randoms
  !------------------------------------------
  	subroutine nb_fixmd_radecranwei(ra,dec,fixmd,insphere_numradecran,weirat,erflag,opt_sumwei_acpt,opt_sumwei_all)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: ra,dec, fixmd 
  		integer, intent(out) :: insphere_numradecran
  		real(dl), intent(out) :: weirat ! this mark the boundary effect: ratio of wei, randoms-inside-boundary .vs. all-randoms
  		integer, intent(out):: erflag ! 0: noerror; 1: out of grid; 2: no random found
  		real(dl), optional, intent(out) :: opt_sumwei_acpt,opt_sumwei_all
  		! LOCAL VARIABLES
  		real(dl) :: RDSR,raspan,RDSRsq, distance_array(max_insphere_num), radec_mass_array(3,max_insphere_num)
  		logical :: acptarray(max_insphere_num), randomlyselect
  		integer :: i,j,l, imin,imax,jmin,jmax, nowinsphere_num, radecranindex, est_insphere_numradecran
  		real(dl) :: h, nowra,nowdec,nowr, mass, wei,sumwei_acpt, sumwei_all, rat, tmpx

		RDSR = 1.0_dl! get_RDSR(dec)
		RDSRsq = RDSR*RDSR
		raspan = RDSR * fixmd

		weirat = 0.0
  		imin = int((ra-raspan-gbgridramin)/gbdeltara +1.0_dl)
  		imax = int((ra+raspan-gbgridramin)/gbdeltara +1.0_dl)
  		jmin = int((dec-fixmd-gbgriddecmin)/gbdeltadec +1.0_dl)
  		jmax = int((dec+fixmd-gbgriddecmin)/gbdeltadec +1.0_dl)
  		
  		erflag = 0
		if(imin<1.or.imin>gb_n_cellra.or.jmin<1.or.jmax>gb_n_celldec) then
			erflag = 1 ! 1 means out of grid range
		else
			erflag = .false.
		endif

		est_insphere_numradecran = 0
		do i = max(1,imin), min(gb_n_cellra,imax)
		do j = max(1,jmin), min(gb_n_celldec,jmax)
			est_insphere_numradecran = est_insphere_numradecran + gb_radecrancell_mat(i,j)%numradecran
		enddo
		enddo
		
		! We donot need so many randoms; 500 must be enough; randomly select a subset of randoms if number > 500
		if(est_insphere_numradecran > 500) then 
			rat = 500.0 / dble(est_insphere_numradecran); randomlyselect = .true.
		else
			rat = 1.0_dl; randomlyselect = .false.
		endif

		insphere_numradecran = 0
		do i = max(1,imin), min(gb_n_cellra,imax)
		do j = max(1,jmin), min(gb_n_celldec,jmax)
			do l = 1, gb_radecrancell_mat(i,j)%numradecran
				if(randomlyselect) then
					call random_number(tmpx)
					if(tmpx > rat) cycle
				endif
				radecranindex = gb_radecrancell_mat(i,j)%iradecranlist(l)
				nowra = gb_radecranlist(radecranindex)%ra
				nowdec= gb_radecranlist(radecranindex)%dec
				nowr = sqrt((nowra-ra)**2.0/RDSRsq+(nowdec-dec)**2.0)
				if(nowr < fixmd) then
					insphere_numradecran = insphere_numradecran+1
					distance_array(insphere_numradecran) = nowr
					radec_mass_array(1,insphere_numradecran) = nowra
					radec_mass_array(2,insphere_numradecran) = nowdec
					radec_mass_array(3,insphere_numradecran) = 1.0_dl !!!Do not consider masses for ra/decran!!! TBU
					acptarray(insphere_numradecran) = gb_radecranlist(radecranindex)%acpt
				endif
			enddo
		enddo
		enddo

		if(insphere_numradecran .eq. 0) then
			erflag = 2 ! 2 means no random found!
			return
		endif

		if(insphere_numradecran > max_insphere_num) then
			print *, 'WARNING (nb_fixmd_radecranwei): # of random overflow: ', &
				insphere_numradecran, max_insphere_num
			erflag = 3
			return
		endif

		sumwei_acpt = 0.0; sumwei_all = 0.0
  		h = fixmd / 2.0
		do i = 1, insphere_numradecran
			nowr = distance_array(i)
			mass = radec_mass_array(3,i)
			wei = mass*w_2dkernel(nowr, h)
			sumwei_all = sumwei_all + wei
			if(acptarray(i)) then
				sumwei_acpt = sumwei_acpt + wei
			endif
!			write(*,'(i3,3f10.3,L2,2f10.3)') i, nowr, mass, wei, acptarray(i), sumwei_acpt, sumwei_all
		enddo
		weirat = sumwei_acpt / sumwei_all
!		print *, weirat
!		call system('sleep 3')
		if(present(opt_sumwei_acpt)) then
			opt_sumwei_acpt = sumwei_acpt
		endif
		if(present(opt_sumwei_all)) then
			opt_sumwei_all  = sumwei_all
		endif
	end subroutine nb_fixmd_radecranwei
 

!!! Test subroutines
  !------------------------------------------
  ! Input a series of points, output the result...
  !------------------------------------------
	subroutine nb_test_smooth()
		integer, parameter :: numpar = 5
		real(dl) :: xyz_mass_array(4,numpar), dxy = 0.3, r0(3), &
			x,y,z,r, mass, rho, dweight, drhodx, drhody, drhodz, h, mu, summu
		integer :: i,j
	
		do j = 1, 4
			r0 = 0.0_dl; x=r0(1); y=r0(2); z=r0(3)
			h = 1.0_dl
			if(j.eq.1) then
				xyz_mass_array(1:4,1) = (/0.0_dl+dxy,  -1.0_dl+dxy, 0.0_dl, 4.0_dl/)
			elseif(j.eq.2) then
				xyz_mass_array(1:4,1) = (/0.0_dl+dxy,  -1.0_dl-dxy, 0.0_dl, 4.0_dl/)
			elseif(j.eq.3) then
				xyz_mass_array(1:4,1) = (/0.0_dl-dxy,  -1.0_dl+dxy, 0.0_dl, 4.0_dl/)
			elseif(j.eq.4) then
				xyz_mass_array(1:4,1) = (/0.0_dl-dxy,  -1.0_dl-dxy, 0.0_dl, 4.0_dl/)
			endif
			xyz_mass_array(1:4,2) = (/0.0_dl+dxy,  1.0_dl+dxy, 0.0_dl, 1.0_dl/)
			xyz_mass_array(1:4,3) = (/0.0_dl+dxy,  1.0_dl-dxy, 0.0_dl, 1.0_dl/)
			xyz_mass_array(1:4,4) = (/0.0_dl-dxy,  1.0_dl+dxy, 0.0_dl, 1.0_dl/)
			xyz_mass_array(1:4,5) = (/0.0_dl-dxy,  1.0_dl-dxy, 0.0_dl, 1.0_dl/)
			rho=0; drhodx=0; drhody=0; drhodz = 0
			do i = 1, numpar
				if(i.eq.1) then
					write(*,'("info",i1,"=[[",3(f10.4,","),i4,"],\")') &
						j,xyz_mass_array(1:3,i), int(xyz_mass_array(4,i))
				elseif(i.lt.numpar) then
					write(*,'(8x,"["3(f10.4,","),i4,"],\")') xyz_mass_array(1:3,i), int(xyz_mass_array(4,i))
				else
					write(*,'(8x,"["3(f10.4,","),i4,"]]")') xyz_mass_array(1:3,i), int(xyz_mass_array(4,i))
				endif
				r = distance(xyz_mass_array(1:3,i),r0,3)
				mass = xyz_mass_array(4,i)
				rho = rho + mass*w_kernel(r, h)
				dweight = der_w_kernel(r,h)
				drhodx = drhodx + mass*(x-xyz_mass_array(1,i)) / r * dweight
				drhody = drhody + mass*(y-xyz_mass_array(2,i)) / r * dweight
				drhodz = drhodz + mass*(z-xyz_mass_array(3,i)) / r * dweight
			enddo
			mu = drhody*10.0/sqrt(drhodx**2.0+drhody**2.0+drhodz**2.0)/10.0
			write(*,'("rlt",i1,"=[",3(f10.4,","),f10.4,"]")') j, drhodx, drhody, drhodz, mu
		enddo
	end subroutine 
	
  !------------------------------------------
  ! Test the NNB search: compare with direct search 
  !------------------------------------------
	subroutine NNBTest()
		integer :: idata, insphere_num, insphere_numran, erflag
		real(dl) :: x,y,z,delta,maxdist,weirat
		print *, 'Check the NNB Search!!!'
		do while(.true.)
			print *, '################################'
			print *, 'input the index of data near which you want to search.'
			print *, '     Range: 1 to ', gb_numdata
			read(*,*) idata
			print *, '  xyz: ', real(gb_xyz_list(1:3,idata))
			print *, 'Will search at x+delta,y+delta,z+delta.'
			print *, 'input delta,#-of-NNBs...'
			read(*,*) delta,insphere_num
			x = gb_xyz_list(1,idata)+delta
			y = gb_xyz_list(2,idata)+delta
			z = gb_xyz_list(3,idata)+delta
			write(*,'(A,3f15.3)') 'Searching for nearest NNBs around ', x,y,z, '...'
			call NNBTest_sub(x,y,z,insphere_num,maxdist)
			print *, 'PRESS ENTER TO DO RNADOM TEST'
			read(*,*) 
			call NNBRanWeiTest_sub(x,y,z,maxdist,insphere_numran,weirat)
			print *, 'Random wei rat (test routine): ', weirat
			call nb_fixmd_ranwei(x,y,z,maxdist,insphere_numran,weirat,erflag)
			if(erflag .eq. 0) then
				print *, 'Random #/weirat from real routine:', insphere_numran, real(weirat)
			endif
		enddo
	end subroutine NNBTest
	
  	subroutine NNBTest_sub(x,y,z,insphere_num,maxdist)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z ! position where you want to calculate rho and drho
  		integer, intent(in) :: insphere_num ! # of nearest neighbors you want to use
  		real(dl), intent(out) :: maxdist
  		! LOCAL VARIABLES
  		integer, allocatable :: index_array(:) ! array of index of nearby data points selected 
  		real(dl), allocatable :: distance_array(:), xyz_mass_array(:,:), tmp(:), smda(:)
  		integer, allocatable :: tmpindex(:), smlablist(:)
  		integer :: i,j,n,allparindex
  		real(dl) :: r0(3), r, mass, dweight, touchbdmass, maxdist1, maxdist2
  		real(dl), allocatable :: all_distance_array(:)

		! Near x,y,z, we require insphere_num of NNB, then get a rough list of n (usually much larger then insphere_num) NBs
  		call NNBSearch_data(x,y,z,insphere_num,selected_list=index_array)
  		n=size(index_array)

		! Calculate the distances to these NBs;
		r0(1)=x; r0(2)=y; r0(3)=z;
  		allocate(tmpindex(n), distance_array(n))
		do i = 1, n
			distance_array(i) = distance(gb_xyz_list(1:3,index_array(i)),r0,3)
			tmpindex(i) = index_array(i)
		enddo

		! From these NBs select out the insphere_num NNBs; 
		allocate(smlablist(insphere_num),smda(insphere_num))
		call ltlablist2(distance_array,n,insphere_num,smda,smlablist)

		! Get the x,y,z,mass of these NNBs
		print *, ' List of NNBs:'
		maxdist1 = 0.0
		allocate(xyz_mass_array(4,insphere_num))
		do i = 1, insphere_num
			allparindex = tmpindex(smlablist(i)) ! this headache sentence get the index of the point in the gb_***_list...
			xyz_mass_array(1:3,i) = gb_xyz_list(1:3,allparindex)
			xyz_mass_array(4,i) = gb_mass_list(allparindex)
			write(*,'(A,i10,4f15.3)') '  Index, x/y/z/distance = ', &
				allparindex, xyz_mass_array(1:3,i), distance(xyz_mass_array(1:3,i),r0,3)
			maxdist1 = max(maxdist1,distance(xyz_mass_array(1:3,i),r0,3))
		enddo

		! Get the NNBs from direct search...		
		print *, ' List of NNBs (direct search):'
		maxdist2 = 0.0
		allocate(all_distance_array(gb_numdata))
		do i = 1, gb_numdata
			all_distance_array(i) = distance(gb_xyz_list(1:3,i),r0,3)
		enddo
		call ltlablist2(all_distance_array,gb_numdata,insphere_num,smda,smlablist)
		do i = 1, insphere_num
			allparindex = smlablist(i)	
			write(*,'(A,i10,4f15.3)') '  Index, x/y/z/distance = ', &
				allparindex, gb_xyz_list(1:3,allparindex), distance(gb_xyz_list(1:3,allparindex),r0,3)
			maxdist2 = max(maxdist2,distance(xyz_mass_array(1:3,i),r0,3))
		enddo
		print *, 'Maximal distance of the two searches: ', real(maxdist1), real(maxdist2)
		maxdist = maxdist1
	end subroutine NNBTest_sub
	
	
  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel; fixed maximal distance
  !------------------------------------------
  	subroutine NNBRanWeiTest_sub(x,y,z,fixmd,insphere_num,weirat)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z, fixmd 
  		integer, intent(out) :: insphere_num
  		real(dl), intent(out) :: weirat ! this mark the boundary effect: ratio of wei, randoms-inside-boundary .vs. all-randoms
  		! LOCAL VARIABLES
  		real(dl) :: distance_array(max_insphere_num), xyz_mass_array(4,max_insphere_num)
  		logical :: acptarray(max_insphere_num)
  		integer :: i,j,k,l, imin,imax,jmin,jmax,kmin,kmax, nowinsphere_num, allparindex,ranindex
  		real(dl) :: r0(3), h, nowr, mass, wei,sumwei_acpt, sumwei_all, maxdist1, maxdist2
		real(dl) :: touchbdmass

  		imin = int((x-fixmd-gbgridxmin)/gbdeltax +1.0)
  		imax = int((x+fixmd-gbgridxmin)/gbdeltax +1.0)
  		jmin = int((y-fixmd-gbgridymin)/gbdeltay +1.0)
  		jmax = int((y+fixmd-gbgridymin)/gbdeltay +1.0)
  		kmin = int((z-fixmd-gbgridzmin)/gbdeltaz +1.0)
  		kmax = int((z+fixmd-gbgridzmin)/gbdeltaz +1.0)
  		
		if(imin<1.or.imin>gb_n_cellx.or.jmin<1.or.jmax>gb_n_celly.or.kmin<1.or.kmax>gb_n_cellz) then
			print *, 'EXIT with error!'; return
		endif

		insphere_num = 0
		r0(1)=x; r0(2)=y; r0(3)=z; maxdist1=0.0

		print *, ' List of NNB Randoms:'
		do i = max(1,imin), min(gb_n_cellx,imax)
		do j = max(1,jmin), min(gb_n_celly,jmax)
		do k = max(1,kmin), min(gb_n_cellz,kmax)
			do l = 1, gb_cell_mat(i,j,k)%numran
				allparindex = gb_cell_mat(i,j,k)%iranlist(l)
				nowr = distance(gb_xyz_list(1:3,allparindex),r0,3)
				if(nowr < fixmd) then
					insphere_num = insphere_num+1
					distance_array(insphere_num) = nowr
					xyz_mass_array(1:3,insphere_num) = gb_xyz_list(1:3,allparindex)
					xyz_mass_array(4,insphere_num) = gb_mass_list(allparindex)
					ranindex = allparindex - gb_numdata
					acptarray(insphere_num) = gb_ranlist(ranindex)%acpt
					write(*,'(A,i5,i10,f15.5)') ' Count, Index, distance = ', insphere_num, allparindex, nowr
					maxdist1 = max(maxdist1, nowr)
				endif
			enddo
		enddo
		enddo
		enddo

		if(insphere_num > max_insphere_num) then
			print *, 'WARNING (NNBRanWeiTest_sub): # of random overflow: ', insphere_num, max_insphere_num
			return
		endif

		if(insphere_num .eq. 0) then
			print *, 'NO RANDOM IN CELL!'
			return
		endif

		sumwei_acpt = 0.0; sumwei_all = 0.0
  		h = fixmd / 2.0

		do i = 1, insphere_num
			nowr = distance_array(i)
			mass = xyz_mass_array(4,i)
			wei = mass*w_kernel(nowr, h)
			sumwei_all = sumwei_all + wei
			if(acptarray(i)) then
				sumwei_acpt = sumwei_acpt + wei
			endif
		enddo
		weirat = sumwei_acpt / sumwei_all
		
		print *, ' List of NNB Randoms (from direct search):'
		j = 0; maxdist2 = 0.0
		do i = gb_numdata, gb_numallpar
			nowr = distance(gb_xyz_list(1:3,i),r0,3)
			if(nowr < fixmd) then
				j = j+1
				write(*,'(A,i5,i10,f15.5)') ' Count, Index, distance = ', j, i, nowr
				maxdist2 = max(maxdist2, nowr)
			endif
		enddo
		print *, '#-ran selected from the two results: ', insphere_num, j
		print *, 'Maxdist: ', maxdist1,maxdist2
		print *, 'sumwei_acpt, sumwei_all, rat = ', sumwei_acpt, sumwei_all, weirat
	end subroutine NNBRanWeiTest_sub
end module ap_smooth	
	

