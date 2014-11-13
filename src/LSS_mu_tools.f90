
!####################################
!This module are tools for mu statistics
!####################################
module LSS_mu_tools
use LSS_tools
implicit none

contains

  !------------------------------------------
  ! mu of a gradient vector at given position
  !------------------------------------------
  	real(dl) function get_mu_of_gradient(x,y,z,dx,dy,dz)
  		real(dl) :: r1, r2, x,y,z,dx,dy,dz
  		r1 = sqrt(x*x+y*y+z*z)
  		r2 = sqrt(dx*dx+dy*dy+dz*dz)
  		get_mu_of_gradient = (x*dx+y*dy+z*dz)/r1/r2
  	end function get_mu_of_gradient
  	
  !------------------------------------------
  ! list of mu from list of pos/gradient
  !------------------------------------------
  	subroutine get_mu_from_gradient_list(pos_list, gradient_list, mu_list)
  		real(dl), intent(in) :: pos_list(:,:), gradient_list(:,:)
  		real(dl), allocatable, intent(out):: mu_list(:)
  		integer :: i, n
  		
  		n = size(pos_list,2)
  		if(size(gradient_list,2).ne.n) then
  		        print*, 'ERROR (get mu from gradient list)! Length of pos_list and gradient_list not matched:', &
  		        	n, size(gradient_list,2)
  		        stop
		endif
		allocate(mu_list(n))
		do i = 1, n
			mu_list(i) = get_mu_of_gradient(pos_list(1,i),pos_list(2,i),pos_list(3,i),gradient_list(1,i),&
				gradient_list(2,i),gradient_list(3,i))
		enddo
	end subroutine get_mu_from_gradient_list
	
  !------------------------------------------
  ! list of mu from list of pos/gradient
  !------------------------------------------
	subroutine histogram(x_data, num_bins, gv_x1, gv_x2, hist_list, hist_er_list, x_av, x_av_er)
		real(dl), intent(in) :: x_data(:)
		integer, intent(in) :: num_bins
		real(dl), optional, intent(in) :: gv_x1, gv_x2
		integer, allocatable, optional, intent(out) :: hist_list(:)
		real(dl), allocatable, optional, intent(out) :: hist_er_list(:)
		real(dl), optional, intent(out) :: x_av, x_av_er
		integer :: i, j, n, m
		real(dl) :: dx, x1, x2, x_mean, x_var, tmp
		
		n = size(x_data)
		if(present(gv_x1)) then
			x1 = gv_x1
			else
			x1 = minval(x_data) 
		endif
		if(present(gv_x2)) then
			x2 = gv_x2
			else
			x2 = maxval(x_data) 
		endif
		
		dx = (x2-x1)/(num_bins + 0.0_dl)
		
		allocate(hist_list(num_bins))
		hist_list = 0
		do i = 1, n
			tmp = x_data(i)
			if(tmp < x1 .or. tmp > x2) &
				cycle
			j = min(int((tmp-x1)/dx)+1,num_bins)
			hist_list(j) = hist_list(j) + 1
		enddo
		
		if(present(hist_er_list)) then
			if(.not.present(hist_list)) then
				print*, 'ERROR! MUST has hist if want error of hist!'
				stop
			endif
			m = sum(hist_list)
			allocate(hist_er_list(num_bins))
			do i = 1, num_bins
				!hist_er_list(i) = hist_list(i) / sqrt(hist_list(i) - 1.0_dl)
				hist_er_list(i) = sqrt(hist_list(i)*(1.0_dl-hist_list(i)/(m+0.0)) )
!				hist_er_list(i) = sqrt(n/(num_bins+0.0))*(1.0_dl-(1/(num_bins+0.0)) )
!				hist_er_list(i) = sqrt(hist_list(i)*(1.0_dl-hist_list(i)/(m+0.0)) )

			enddo
		endif
			
		if(present(x_av) .or. present(x_av_er)) then
			x_mean =  0.0_dl
			do i = 1, n
				x_mean = x_mean + x_data(i)
			enddo
			x_mean = x_mean / (n+0.0_dl)
		endif

		if(present(x_av)) x_av = x_mean
		
		if(present(x_av_er)) then
			x_var = 0.0_dl
			do i = 1, n
				x_var = x_var + (x_data(i)-x_mean)**2.0
			enddo
			x_var = x_var / (n-1.0_dl)
			x_av_er = sqrt(x_var / (n-1.0_dl))
		endif
	end subroutine histogram

  !------------------------------------------
  ! histogram of mu from list of pos/gradient
  !------------------------------------------	
	subroutine mu_histogram(mu_data, num_bins, abs_mu, mu_hist, mu_hist_er, mu_av, mu_av_er, chisq)
		real(dl), intent(inout) :: mu_data(:)
		integer, intent(in) :: num_bins
		logical, optional, intent(in) :: abs_mu
		integer, allocatable, optional, intent(out) :: mu_hist(:)
		real(dl), allocatable, optional, intent(out) :: mu_hist_er(:)
		real(dl), intent(out) :: mu_av, mu_av_er, chisq
		real(dl) :: mu1, mu2, tmp
		integer :: i, n
		
		n = size(mu_data)
		if(present(abs_mu).and.abs_mu) then
			do i = 1, n
				mu_data(i) = abs(mu_data(i))
			enddo
			mu1 = 0.0_dl; mu2 = 1.0_dl
		else
			mu1 = -1.0_dl; mu2 = 1.0_dl
		endif
		
		call histogram(mu_data, num_bins, mu1, mu2, mu_hist, mu_hist_er, mu_av, mu_av_er)
		
		tmp = n / (num_bins + 0.0_dl)
		chisq = 0.0_dl
		do i = 1, num_bins
			chisq = chisq + ((mu_hist(i)-tmp)/mu_hist_er(i))**2.0
		enddo
	end subroutine mu_histogram
	
  !------------------------------------------
  ! the same as above, but only calculate chisq
  !------------------------------------------	
	real(dl) function chisq_of_mu_data(mu_data, num_bins, abs_mu)
		real(dl), intent(inout) :: mu_data(:)
		integer, intent(in) :: num_bins
		logical, optional, intent(in) :: abs_mu
		integer, allocatable :: mu_hist(:)
		real(dl), allocatable :: mu_hist_er(:)
		real(dl) :: mu_av, mu_av_er, chisq
		real(dl) :: mu1, mu2, tmp
		integer :: i, n
		
		n = size(mu_data)
		if(present(abs_mu).and.abs_mu) then
			do i = 1, n
				mu_data(i) = abs(mu_data(i))
			enddo
			mu1 = 0.0_dl; mu2 = 1.0_dl
		else
			mu1 = -1.0_dl; mu2 = 1.0_dl
		endif
		
		call histogram(mu_data, num_bins, mu1, mu2, mu_hist, mu_hist_er, mu_av, mu_av_er)
		
		tmp = n / (num_bins + 0.0_dl)
		chisq = 0.0_dl
		do i = 1, num_bins
			chisq = chisq + ((mu_hist(i)-tmp)/mu_hist_er(i))**2.0
		enddo
		chisq_of_mu_data = chisq
	end function chisq_of_mu_data


  !------------------------------------------
  ! another definition of chisq
  !------------------------------------------	
	real(dl) function chisq_of_mu_data2(mu_data)
		real(dl), intent(in) :: mu_data(:)
		real(dl), allocatable :: tmp(:)
		real(dl) :: mumean, muvar, muer
		integer :: i, n
		n = size(mu_data)
		allocate(tmp(n))
		do i = 1, n
			tmp(i)= abs(mu_data(i))
		enddo
		call get_mean_var(tmp, mumean, muvar)
		muer = sqrt(muvar) / sqrt(n-1.0)
!		print *, '(chisq_of_mu_data2): mumean, muer = ', mumean, muer
		chisq_of_mu_data2 = (mumean-0.5_dl)**2.0/muer**2.0
	end function chisq_of_mu_data2
	
  !------------------------------------------
  ! df chisqs: comparing different volumes
  !------------------------------------------		
  	subroutine dfchisqs(numchisqs, nbinlist, absmudata, rlist, n, dfchisqlist)
  		! Dummy
  		integer, intent(in) :: numchisqs, n, nbinlist(numchisqs)
  		real(dl), intent(in) :: absmudata(n), rlist(n)
  		real(dl), intent(out) :: dfchisqlist(numchisqs)
  		! Local
  		integer :: num_bins, i, j
  		real(dl) :: rmin, rmax, nowmuav, nowchisq, summu, sumwei
  		real(dl), allocatable :: muavlist(:), muerlist(:), ravlist(:), muvarlist(:)
  		
  		call find_min_max(rlist, n, rmin, rmax)
		do i = 1, numchisqs
			num_bins = nbinlist(i)
			if(num_bins .lt. 2) then
				print *, 'ERROR (dfchisqs)!!! num_bins must >= 2: ', i, nbinlist(i)
				stop
			endif
			allocate(muavlist(num_bins),muerlist(num_bins),ravlist(num_bins),muvarlist(num_bins))
			call eqvl_binned_quan(absmudata, rlist, rmin, rmax, size(absmudata), num_bins, &	
  					muavlist, muerlist, ravlist, muvarlist)
  			summu = 0; sumwei = 0
  			do j = 1, num_bins
  				summu = summu + muavlist(j)/(muerlist(j))**2.0
  				sumwei = sumwei + 1.0/(muerlist(j))**2.0
  			enddo
  			nowmuav = summu / sumwei
!  			nowmuav = sum(muavlist)/dble(num_bins)
			nowchisq = 0.0_dl
			do j = 1, num_bins
				nowchisq = nowchisq + ((muavlist(j)-nowmuav)/muerlist(j))**2.0
			enddo
			dfchisqlist(i) = nowchisq
			deallocate(muavlist,muerlist,ravlist,muvarlist)
		enddo
	end subroutine dfchisqs
	
  !------------------------------------------
  ! another definition of chisq
  !------------------------------------------	
	real(dl) function weighted_chisq(mu_data, weight, n)
		! Dummy
		real(dl), intent(in) :: mu_data(n), weight(n)
		integer :: n
		! Output
		real(dl) :: sumw, sumwsq, summu, sumwmu, mumean, wmumean, wmean
		real(dl) :: muvar, wvar, wmuvar 
		real(dl) :: er
		integer :: i
		
		! weighted mean of mu
		sumw = 0.0_dl	
		sumwmu = 0.0_dl
		do i = 1, n
			sumw = sumw + weight(i)
			sumwmu = sumwmu + weight(i)*abs(mu_data(i))
		enddo
		wmumean = sumwmu / sumw

		! weighted variance of mu
		wmuvar = 0.0_dl
		do i = 1, n
			wmuvar = wmuvar + weight(i)*(abs(mu_data(i))-wmumean)**2.0 
			! when calculating weighted variance, also weight each sample
		enddo
		wmuvar = wmuvar / sumw
		
		! Error Estimator
		er = sqrt(wmuvar) / sqrt(n-1.0_dl)
		weighted_chisq = ((wmumean - 0.5_dl)/er)**2.0
	end function weighted_chisq	


  !------------------------------------------
  ! shift the mu data
  !------------------------------------------	
	subroutine shift_mu_data(mu_data, mumin, mumax, mushift)
		real(dl), intent(inout) :: mu_data(:)
		real(dl), intent(in) :: mumin, mumax, mushift
		real(dl) :: x
		integer :: i, n
!		logical :: changemudata
		
		n = size(mu_data)
		do i = 1, n
			x = mu_data(i) + mushift
			do while(x .gt. mumax)
				x = mumin + (x-mumax)
			enddo
			do while(x .lt. mumin)
				x = mumax - (mumin - x)
			enddo
			mu_data(i) = x
		enddo
		
	end subroutine shift_mu_data

  !------------------------------------------
  ! another definition of chisq
  !------------------------------------------	
	real(dl) function chisq_of_mu_data_shift(mu_data, num_bins, num_shift, gv_printinfo)
		! DUMMY
		real(dl), intent(in) :: mu_data(:)
		integer, intent(in) :: num_bins, num_shift
		logical, intent(in), optional :: gv_printinfo
		! LOCAL
		logical :: printinfo
		real(dl), allocatable :: tmp(:)
		real(dl) :: chisqlist(num_shift+1), dshift, mean, var
		integer :: i, n

		if(present(gv_printinfo)) then
			printinfo = gv_printinfo
		else
			printinfo = .false.
		endif
		
		dshift = 2.0_dl / (num_bins+0.0) / (num_shift+1.0_dl)
		
		n = size(mu_data)
		allocate(tmp(n))
		do i = 1, n
			tmp(i) = mu_data(i)
		enddo
		
		chisqlist(1) = chisq_of_mu_data(tmp, num_bins)
		
		do i = 1,num_shift
			call shift_mu_data(tmp, -1.0_dl, 1.0_dl, dshift)
			chisqlist(i+1) = chisq_of_mu_data(tmp, num_bins)
		enddo
		
		chisq_of_mu_data_shift = sum(chisqlist) / (num_shift + 1.0)
		if(printinfo) then
			write(*,*) ' list of chisqs:', real(chisqlist)
			call get_mean_var(chisqlist, mean, var)
			write(*,*) 'Mean, var, sqrt(var) = ', real(mean), real(var), real(sqrt(var))
		endif

	end function chisq_of_mu_data_shift

	
end module LSS_mu_tools
