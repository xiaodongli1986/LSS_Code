!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module LSS_constants_types

	implicit none

	!###################################################
	! constants
	integer, parameter :: char_len = 1000
	integer, parameter :: dl = KIND(1.0d0)
	integer, parameter :: sp = KIND(1.0)
	real(dl), parameter :: const_pi = 3.141592653589793d0
	real(dl), parameter :: const_c = 299792.458d0 !unit: km/s
	real(dl), parameter :: logzero = 1.0e30
	real(dl), parameter :: const_Mpctokm = 3.08567758d19 ! Mpc / km
	real(dl), parameter :: const_gyrtos = 31557600000000000.0d0 ! Gyr to seconds

	! types
	!! collection of settings 
	type :: chisq_settings
		integer :: smnum, num_in_x ! smooth-number; grid size in x-direction
		logical :: print_info = .false.
		integer :: numdrop 
		logical, allocatable :: dropval(:), dropdval(:)
		real(dl), allocatable :: lowdropvalratio(:), highdropvalratio(:), lowdropdvalratio(:), highdropdvalratio(:)
	end type
	
	!! allocatable array with given length
	type :: dy_array
		integer :: n
  		real(dl), allocatable :: array(:)
  	end type
	
	!! group structure (used in quick sort)
	type group
		integer :: order    ! original order of unsorted data
		real(dl) :: value       ! values to be sorted by
	end type group
end module LSS_constants_types
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module LSS_variables
use LSS_constants_types
	implicit none

	!###################################################
	! useful variables	

	!! a variable used to control global print for test; used for debug
	logical, parameter :: gbtp = .true. 

	! array used to get ra-dec-span-ratio as a function of dec
	!  Shape of elliptic when projecting a perfect sphere to ra/dec plan:
	!   	We will get an elliptic rather than a perfect circle.
	!	Main axis in DEC direction is simply radius-of-sphere/distance-of-sphere.
	!	Main axis in RA direction is longer than axis in DEC direction; the ratio is related with the DEC coordinate.
	!	DEC_vs_RAspan-to-DECspan.txt: relation between DEC and ratio. 
	!	For BOSS catalogues, DEC ranges from -12 to 69.73; in that range all errors controlled to <~0.5%.
	integer, public :: gb_RDSR_num ! RDSR means Ra Dec Span Rat!!!
	real(dl), allocatable :: gb_RDSR_decs(:), gb_RDSR_rats(:)
	real(dl) :: gb_RDSR_decmin, gb_RDSR_decmax, gb_RDSR_deltadec

end module LSS_variables
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
module LSS_tools
use LSS_variables

implicit none
			
contains

  !------------------------------------------
  ! how many lines are there in a file
  !------------------------------------------
	subroutine count_line_number (file_name, line_number)
		! Dummy
		integer :: line_number
		character(len=char_len) :: file_name
		! Local
		character(len=char_len) :: inline
    
		open(unit=1456,FILE=file_name,ERR=2)
    
		line_number = 0
		do While(1 .eq. 1)
			read(1456,*,ERR=3,end=100) inline
			line_number = line_number + 1
		enddo

2		write(*,*) "Error occurs when opening the file ", file_name
		close(1456)
		stop

3		write(*,*) "Error occurs when couting the lines of the file ", file_name
		close(1456)
		stop

100		close(1456)
	end subroutine count_line_number
  !------------------------------------------
  ! read in a file into allocatable array
  !------------------------------------------	
	subroutine read_in (file_name, n_column, n_lines, total_data)
		character(len=char_len), intent(in) :: file_name
		integer, intent(in)  :: n_column
		integer, intent(out) :: n_lines
		real(dl), allocatable :: total_data(:,:)
		integer :: i
		call count_line_number(file_name, n_lines)
		allocate(total_data(n_lines, n_column))
		open(unit=4321,file=file_name)
		do i = 1, n_lines
			read(4321, *) total_data(i,1:n_column)
		enddo
		close(4321)
	end subroutine read_in 
  !------------------------------------------
  ! read in a file into allocatable array,
  !  but with the opposite format of the previous function
  !------------------------------------------	
	subroutine read_in_revfmt (file_name, n_column, n_lines, total_data)
		character(len=char_len), intent(in) :: file_name
		integer, intent(in)  :: n_column
		integer, intent(out) :: n_lines
		real(dl), allocatable :: total_data(:,:)
		integer :: i
		call count_line_number(file_name, n_lines)
		allocate(total_data(n_column, n_lines))
		open(unit=4321,file=file_name)
		do i = 1, n_lines
			read(4321, *) total_data(1:n_column,i)
		enddo
		close(4321)
	end subroutine read_in_revfmt

  !-----------------------------------------------------------
  ! Fmt output
  !-----------------------------------------------------------
        character(len=char_len) function WriteFmtEFloat(num_float, floatvalues)
                ! Dummy
                integer, intent(in) :: num_float
                real(dl), intent(in) :: floatvalues(num_float)
                ! Local
                character(len=char_len) :: tmpstr, tmpstr0
                integer :: i
                WriteFmtEFloat = ''
                do i = 1, num_float
                        write(tmpstr0, '(e14.7)') floatvalues(i)
                        WriteFmtEFloat = trim(adjustl(WriteFmtEFloat))//' '//trim(adjustl(tmpstr0))
                enddo

        end function WriteFmtEFloat
  !-----------------------------------------------------------
  ! Output a 1-dimensional float array to a file.
  !-----------------------------------------------------------
	subroutine output_1d(file_name, output_data, n)
		! Dummy
		character(LEN=char_len), intent(in) :: file_name
		real(dl), intent(in) :: output_data(n)
		integer, intent(in) :: n
		! Local
		integer  :: i

		open(unit=60857,FILE=file_name,err=22)

		do i=1, n
			write(60857,'(e18.11)',err=33) output_data(i)
		enddo

		close(60857)
		return

22		write(*,*) "(output_1d) Error occurs when opening the file ", file_name
		close(60857)
		stop

33		write(*,*) "(output_1d) Error occurs when writing into the file ", file_name
		close(60857)
		stop
	end subroutine output_1d	
  !-----------------------------------------------------------
  ! Output a 2d float array to a file.
  !-----------------------------------------------------------
        subroutine output_2d(file_name, output_data)
                ! Dummy
                character(LEN=char_len) :: file_name
                real(dl) :: output_data(:,:)
                ! Local
                integer :: d1, d2, i, j

                d1 = size(output_data, 1)
                d2 = size(output_data, 2)

                open(unit=9873,FILE=file_name,err=22)

                do i=1, d1
			do j = 1, d2-1
                        	write(9873,'(e18.11,"  ",$)',err=33) output_data(i,j)
			enddo
                       	write(9873,'(e18.11)',err=33) output_data(i,d2)
                enddo

                close(9873)
                return

22              write(*,*) "Error occurs when opening the file ", file_name
                close(9873)
                stop

33              write(*,*) "Error occurs when writing into the file ", file_name
                close(9873)
                stop
	end subroutine output_2d
  !-----------------------------------------------------------
  ! Output a 2d float array to a file; 
  !  but the format is opposite to the previous function.
  !-----------------------------------------------------------
	subroutine output_2d_revfmt(file_name, output_data)
		character(LEN=char_len) :: file_name
		real(dl) :: output_data(:,:)
		integer :: d1, d2 
		integer  :: i,j

		d1 = size(output_data, 1)
		d2 = size(output_data, 2)

		open(unit=9873,FILE=file_name,err=22)

		do i=1, d2
			do j = 1, d1-1
				write(9873,'(e18.11,"  ")',advance='no',err=33) output_data(j,i)
			enddo
			write(9873,'(e18.11)',err=33) output_data(d1,i)
		enddo

		close(9873)
		return

22		write(*,*) "Error occurs when opening the file ", file_name
		close(9873)
		stop

33		write(*,*) "Error occurs when writing into the file ", file_name
		close(9873)
		stop
	end subroutine output_2d_revfmt	
	
	
  !------------------------------------------
  ! 2th order interpolating function.
  !------------------------------------------
	real(dl) function intpl_vl(x, x1, f1, x2, f2, x3, f3)
		real(dl) :: x, x1, f1, x2, f2, x3, f3
		real(dl) :: d1, d2, d3, d12, d13, d23
		d1  = x - x1
		d2  = x - x2
		d3  = x - x3
		d12 = x1 - x2
		d13 = x1 - x3
		d23 = x2 - x3
		intpl_vl = f1*d2*d3/(d12*d13) - f2*d1*d3/(d12*d23) + f3*d1*d2/(d13*d23)
	end function intpl_vl		
	! the same as before, but with real(dl)
	real(dl) function intpl_vl_dp(x, x1, f1, x2, f2, x3, f3)
		real(dl) :: x, x1, f1, x2, f2, x3, f3
		real(dl) :: d1, d2, d3, d12, d13, d23
		d1  = x - x1
		d2  = x - x2
		d3  = x - x3
		d12 = x1 - x2
		d13 = x1 - x3
		d23 = x2 - x3
		intpl_vl_dp = f1*d2*d3/(d12*d13) - f2*d1*d3/(d12*d23) + f3*d1*d2/(d13*d23)
	end function intpl_vl_dp
	

  !------------------------------------------
  ! get the theta, phi according to x,y. 2D
  !------------------------------------------
  	real(dl) function getangle(x,y)
  		real(dl) :: x, y
  		getangle = asin(y/sqrt(x*x+y*y))
  		if(x<0) then
  			getangle = const_pi - getangle
  		elseif(x>0 .and. y<0) then
  			getangle = 2*const_pi + getangle
  		endif
  	end function getangle
  !------------------------------------------
  ! get r, theta, phi according to x,y,z. 3D
  !------------------------------------------
  	subroutine getSC(x,y,z,r,theta,phi)
  		real(dl) :: x, y,z, r,theta,phi
		r=sqrt(x*x+y*y+z*z)
		theta=getangle(z,sqrt(x*x+y*y))
		phi=getangle(x,y)
  	end subroutine getSC
  !------------------------------------------
  ! get r, theta, phi according to x,y,z. 3D
  !------------------------------------------  	
	subroutine radecr_to_xyz(ra,dec,r, x,y,z)
		! DUMMY
		real(dl), intent(in) :: ra,dec,r
		real(dl), intent(out) :: x,y,z
		! LOCAL
		real(dl) :: ra2, dec2
		real(dl), parameter :: degreefac = const_pi/180.0
		ra2 = ra * degreefac
		dec2 = dec * degreefac
		x = r*cos(dec2) * cos(ra2)
		y = r*cos(dec2) * sin(ra2)
		z = r*sin(dec2)
	end subroutine radecr_to_xyz
  !------------------------------------------
  ! get r, ra, dec according to x,y,z
  !------------------------------------------
  	subroutine xyz_to_radecr(x,y,z,ra,dec,r)
		! DUMMY
		real(dl), intent(in) :: x,y,z
		real(dl), intent(out) :: ra,dec,r
		! LOCAL
		real(dl), parameter :: degreefac = const_pi/180.0
		r = sqrt(x*x+y*y+z*z)
		dec = (const_pi/2.0 - getangle(z,sqrt(x*x+y*y))) / degreefac
		ra = getangle(x,y) / degreefac
  	end subroutine xyz_to_radecr
  !------------------------------------------
  ! find root: f(x) == f0
  !------------------------------------------  	
  	real(dl) function findroot(f, f0, xl, xr, gv_tol, gv_maxstep)
		real(dl), external :: f
  		real(dl) :: f0, xl, xr, tol, xm, fl, fr, fm
  		real(dl), optional :: gv_tol, gv_maxstep
  		integer :: maxstep, step
  		
  		if(present(gv_tol)) then
  			tol = gv_tol
  			else
  			tol = 1.0e-5
		endif
		if(present(gv_maxstep)) then
			maxstep = gv_maxstep
			else
			maxstep = 500
		endif
				
		fl = f(xl); fr = f(xr)
		
		if(abs(fl-f0)<tol) then
			findroot = xl; return
		elseif(abs(fr-f0)<tol) then
			findroot = xr; return
		elseif( ((fl>f0).and.(fr>f0)) .or. (fl<f0).and.(fr<f0)) then
			print *, 'ERROR! fl, fr both larger/smaller than f0! fl, fr, f0 = ', fl, fr, f0
			findroot = -1; return
		endif
		xm = xl+(xr-xl)*(f0-fl)/(fr-fl); fm = f(xm)
		step = 0
		do while(abs(fm-f0)>tol)
			if((fm > f0 .and. fl > f0) .or. (fm < f0 .and. fl < f0)) then
				xl = xm; fl = fm;
				else
				xr = xm; fr = fm;
			endif
			xm = xl+(xr-xl)*(f0-fl)/(fr-fl); fm = f(xm)
			step = step + 1
			if(step > maxstep) then
				findroot = xm
				print *, 'exit at maxstep = ', maxstep, '. x, f, f0, diff = ', xm, fm, f0, abs(fm-f0)
			endif
		enddo
		findroot = xm
	end function findroot
		
  !------------------------------------------
  ! square summing
  !------------------------------------------	
	real(dl) function square_sum(p,n)
		! Dummy
		integer :: n
		real(dl) :: p(n)
		! Local
		integer :: i
		square_sum = 0
		do i = 1, n
			square_sum = square_sum + p(i)**2.0
		enddo
	end function square_sum
	
  !------------------------------------------
  ! root mean square
  !------------------------------------------	
	real(dl) function rms(p,n)
		! Dummy
		integer :: n
		real(dl) :: p(n)
		rms = sqrt(square_sum(p,n))
	end function rms

  !------------------------------------------
  ! distance of two vectors
  !------------------------------------------	
	real(dl) function distance(p0, p1,n)
		! Dummy
		integer :: n
		real(dl) :: p0(n), p1(n)
		!local
		integer :: i		
		distance = rms(p0(1:n)-p1(1:n),n)
	end function distance

  !------------------------------------------
  ! cubic spline kernel 
  !------------------------------------------		
	real(dl) function w_kernel_hfree(q)
		real(dl) :: q, sigma = 1.0/ const_pi
		if(q<0.0_dl) then
			w_kernel_hfree =0.0_dl; return
		endif
		if(q<1.0_dl) then
			w_kernel_hfree = sigma*(0.25_dl*(2.0_dl-q)**3.0_dl - (1.0_dl-q)**3.0_dl)
			return
		endif
		if(q<2.0_dl) then
			w_kernel_hfree = sigma*0.25_dl*(2.0_dl-q)**3.0_dl
			return
		endif
		w_kernel_hfree = 0.0
	end function w_kernel_hfree
	!kernel with h
	real(dl) function w_kernel(r, h)
		real(dl) r, h, q
		q = r/h
		w_kernel = w_kernel_hfree(q) / h**3.0
	end function w_kernel
	!kernel with h, 2d
	real(dl) function w_2dkernel(r, h)
		real(dl) :: r, h, q, coef = 10.0_dl/7.0_dl ! coef = sigma2d/sigma3d = (10.0/(7.0*const_pi)) / (1.0/ const_pi) = 10.0/7.0
		q = r/h
		w_2dkernel = coef*w_kernel_hfree(q) / h**2.0
	end function w_2dkernel
	
	
  !------------------------------------------
  ! derivative of cubic spline kernel
  !------------------------------------------		
	real(dl) function der_w_kernel_hfree(q)
		real(dl) :: q, sigma = 1.0/ const_pi
		if(q<0.0_dl) then
			der_w_kernel_hfree =0.0_dl; return
		endif
		if(q<1.0_dl) then
			der_w_kernel_hfree = sigma * (-0.75_dl*(2.0_dl-q)**2.0_dl + 3.0_dl*(1.0_dl-q)**2.0_dl)
			return
		endif
		if(q<2.0_dl) then
			der_w_kernel_hfree = sigma * (-0.75_dl)*(2.0_dl-q)**2.0_dl
			return
		endif
		der_w_kernel_hfree = 0.0_dl
	end function der_w_kernel_hfree
	!derivative of kernel with h
	real(dl) function der_w_kernel(r, h)
		real(dl) r, h, q
		q = r/h
		der_w_kernel = der_w_kernel_hfree(q) / h**4.0
	end function der_w_kernel

  !------------------------------------------
  !# return a list, with center i, number of elements 2*di+1, 
  !# but confined in the range (imin, imax)
  !------------------------------------------  
	subroutine ilist(i, di, imin, imax, i1, i2)
		integer :: i, di, imin, imax, i1, i2
	    	i1 = i - di 
		i2 = i + di 
		if(i1 < imin .and. i2 > imax) then
        		print *, 'ERROR! overflow: imin, imax, 2*di+1 = ', imin, imax, 2*di+1
        		stop
        	endif
        	if(i1 < imin) then
        		i1 = imin; i2 = imin + 2*di; return
        	endif 
		if(i2 > imax) then
        		i1 = imax - 2*di; i2 = imax; return
        	endif
        	!print *, i, di, i1, i2
	end subroutine ilist
	
  !---------------------------------------------------------------
  ! Simpson Integration
  !---------------------------------------------------------------
	real(dl) function Simpson(fun,xleft,xright,N)   
		real(dl), external :: fun
		real(dl) :: x1,x2,xleft,xright,BC
		real(dl) :: f1,f2
		integer :: i, N
		BC=(xright-xleft)/DBLE(N)
		x1=xleft;x2=x1+BC;
		f1=fun(x1);f2=fun(x2);Simpson=(f1+fun((x1+x2)*0.5_dl)*4.0_dl+f2)*BC/6.0_dl;
	
		do i = 2,N
			x1=x2;f1=f2;x2=x2+BC;
			f2=fun(x2);Simpson=Simpson+(f1+fun((x1+x2)*0.5_dl)*4.0_dl+f2)*BC/6.0_dl; 
		enddo
	end function Simpson
	
  !---------------------------------------------------------------
  ! 4rd order Runge-Kutta 
  !---------------------------------------------------------------
	real(dl) function RK(dfundz, zleft, funleft, zright, N)
		! Dummy
		real(dl), external :: dfundz
		real(dl), intent(in) :: zleft, funleft, zright
		integer, intent(in) :: N
		! Local
		real(dl) :: BC, K1, K2, K3, K4, nowz, nowfun
		integer :: i
		BC = (zright-zleft) / N
		nowz = zleft
		nowfun = funleft
		do i = 1, N
			K1 = dfundz(nowz, nowfun)
			K2 = dfundz(nowz + 0.5d0*BC, nowfun + 0.5d0*BC*K1)
			K3 = dfundz(nowz + 0.5d0*BC, nowfun + 0.5d0*BC*K2)
			K4 = dfundz(nowz + BC, nowfun + BC*K3)
			nowfun = nowfun + BC * (K1+2.0d0*K2+2.0d0*K3+K4) / 6.0d0
			nowz = nowz + BC
		enddo
		RK = nowfun
	end function RK


  !---------------------------------------------------------------
  ! Binning quantities according to their r.
  !---------------------------------------------------------------
  	subroutine binned_quan(quan_val_list, quan_r_list, rmin, rmax, ndat, num_bins, &
  		quan_val_av_list, quan_val_er_list, r_av_list, quan_val_var_list, quan_medval_list)
  		! Dummy
  		real(dl), intent(in) :: quan_val_list(ndat), quan_r_list(ndat) 
  		real(dl), intent(in) :: rmin, rmax
  		integer, intent(in) :: ndat, num_bins
  		real(dl), intent(out) :: quan_val_av_list(num_bins), quan_val_er_list(num_bins), r_av_list(num_bins)
  		real(dl), optional, intent(out) :: quan_val_var_list(num_bins), quan_medval_list(num_bins)
  		real(dl) :: deltar, r, quan_val_av, quan_val_var, quanval1, quanval2, quanvalmid
  		integer :: i, n, m, j, binned_i, numsmaller
  		integer :: binned_num_list(num_bins)
  		type(dy_array) :: binned_quan_val_list(num_bins), binned_r_list(num_bins)
  		
  		n = size(quan_val_list)
  		if(size(quan_r_list) .ne. n) then
  			print *, 'ERROR! # of quan and r mismatches: ', n, size(quan_r_list)
  			stop
  		endif
  		
  		deltar = (rmax - rmin) / (num_bins + 0.0_dl)
  				
		! how many elements in each bin
		binned_num_list = 0
  		do i = 1, n
  			r = quan_r_list(i)
  			binned_i = max(min(int((r-rmin)/deltar)+1,num_bins),1)
  			binned_num_list(binned_i) = binned_num_list(binned_i) + 1
  		enddo
  			
		! which elements in which bin
  		do binned_i = 1, num_bins
	  		m = binned_num_list(binned_i) 
  			binned_quan_val_list(binned_i)%n = m
  			binned_r_list(binned_i)%n = m
  			allocate(binned_quan_val_list(binned_i)%array(m))
  			allocate(binned_r_list(binned_i)%array(m))
  		enddo
		binned_num_list = 0  			
  		do i = 1, n
  			r = quan_r_list(i)
  			binned_i = max(min(int((r-rmin)/deltar)+1,num_bins),1)
  			binned_num_list(binned_i) = binned_num_list(binned_i) + 1
			j = binned_num_list(binned_i)
			binned_quan_val_list(binned_i)%array(j) = quan_val_list(i)
			binned_r_list(binned_i)%array(j) = quan_r_list(i)
  		enddo
  		
  		do binned_i = 1, num_bins
	  		m = binned_num_list(binned_i) 
	  		quan_val_av = sum(binned_quan_val_list(binned_i)%array) / (m+0.0)
	  		quan_val_var = 0.0
	  		do j = 1, m
	  			quan_val_var = quan_val_var + &
	  			(binned_quan_val_list(binned_i)%array(j)-quan_val_av)**2.0
	  		enddo
	  		quan_val_var = quan_val_var / (m-1.0) 
  			quan_val_av_list(binned_i) = quan_val_av
			if(present(quan_val_var_list)) &
				quan_val_var_list(binned_i) = quan_val_var
  			quan_val_er_list(binned_i) = sqrt(quan_val_var) / sqrt(m-1.0)
  			r_av_list(binned_i) = sum(binned_r_list(binned_i)%array) / (m+0.0)
  			if(present(quan_medval_list)) then
  				call find_min_max(binned_quan_val_list(binned_i)%array, m, quanval1, quanval2)
  				quanvalmid = ( quanval1 + quanval2 ) / 2.0
  				call count_num_smaller(binned_quan_val_list(binned_i)%array, m, quanvalmid,numsmaller)
  				do while(abs(numsmaller - m/2) > 2)
  					if(numsmaller < m/2) then
  						quanval1 = quanvalmid
  						quanvalmid = (quanval2+quanvalmid) / 2.0
  					else
  						quanval2 = quanvalmid
  						quanvalmid = (quanval1+quanvalmid) / 2.0
  					endif
	  				call count_num_smaller(binned_quan_val_list(binned_i)%array,m,quanvalmid,numsmaller)
	  			enddo
	  			quan_medval_list(binned_i) = quanvalmid
	  		endif
  		enddo  		
  	end subroutine binned_quan


  !---------------------------------------------------------------
  ! The same as binned_quan, but keep equal volume rather than distance
  !  Frequently used to estimate mu data...
  !---------------------------------------------------------------
	subroutine eqvl_binned_quan(mulist, rlist, rmin, rmax, n, numbins, &	
  					muavlist, muerlist, ravlist, muvarlist)
		! Dummy
		real(dl), intent(in) :: mulist(n), rlist(n), rmin, rmax
		integer, intent(in) :: n, numbins
		real(dl), intent(out) :: muavlist(numbins), muerlist(numbins), muvarlist(numbins), ravlist(numbins)
		! Local
		real(dl) :: mu, r, rmincube, deltavol
		integer :: i, i_bin, numinbin(numbins)
		
		rmincube = rmin**3.0
		deltavol = (rmax**3.0 - rmin**3.0) / dble(numbins)
		
		numinbin(1:numbins) = 0
		muavlist(1:numbins) = 0.0
		muvarlist(1:numbins) = 0.0
		muerlist(1:numbins) = 0.0
		ravlist(1:numbins) = 0.0
		
		do i = 1, n
			mu = mulist(i)
			r = rlist(i)
			i_bin = int((r**3.0-rmincube)/deltavol)+1
			i_bin = min(i_bin, numbins)
			i_bin = max(i_bin, 1)
			numinbin(i_bin) = numinbin(i_bin) + 1
			muavlist(i_bin) = muavlist(i_bin) + mu
			muvarlist(i_bin) = muvarlist(i_bin) + mu**2.0
			ravlist(i_bin) = ravlist(i_bin) + r
		enddo
		
		do i_bin = 1, numbins
			muavlist(i_bin) = muavlist(i_bin) / dble(numinbin(i_bin))
			muvarlist(i_bin) = muvarlist(i_bin) / dble(numinbin(i_bin)) - muavlist(i_bin)**2.0
			ravlist(i_bin) = ravlist(i_bin) / dble(numinbin(i_bin))
			muerlist(i_bin) = sqrt(muvarlist(i_bin)) / sqrt(dble(numinbin(i_bin)) - 1.0)
		enddo
	end subroutine eqvl_binned_quan
		

  !---------------------------------------------------------------
  ! Estimating mean and variance of an array
  !---------------------------------------------------------------
  	subroutine get_mean_var(val_list, mean, var)
  		real(dl), intent(in) :: val_list(:)
  		real(dl), intent(out) :: mean
  		real(dl), optional, intent(out) :: var
  		real(dl) :: x
  		integer :: i, n
  		
  		n = size(val_list)
  		x = 0.0_dl
		do i = 1, n
			x = x + val_list(i)
		enddo
		mean = x / (n + 0.0_dl)
		if(present(var)) then
			var = 0.0_dl
			do i = 1, n
				var = var + (val_list(i) - mean)**2.0
			enddo
			var = var / (n-1.0_dl)
		endif
  	end subroutine get_mean_var
  
  !---------------------------------------------------------------
  ! Random number in range
  !---------------------------------------------------------------  
  	real(dl) function randnum(xmin,xmax)
		! Dummy
  		real(dl), intent(in) :: xmin, xmax
  		! Local
  		real(dl) :: x
  		call random_number(x)
  		randnum = xmin + (xmax-xmin)*x
  	end function randnum
  !---------------------------------------------------------------
  ! Random integer in range
  !---------------------------------------------------------------  
  	integer function randint(imin,imax)
		! Dummy
  		integer, intent(in) :: imin, imax
  		! Local
  		real(dl) :: x
  		call random_number(x)
  		randint = imin + int((imax-imin+1)*x)
  	end function randint
	
  !! Subroutines do Sort
  !------------------------------------------
  ! 1d Bubble
  !------------------------------------------
	!Buuble sorting: 1d real(dl) array
	subroutine Bubble_1d(A,Aindex,max_num)
		real(dl) :: A(:), tmp
		integer, allocatable, optional, intent(out) :: Aindex(:)
		integer, optional, intent(in) :: max_num
		integer :: n, m, i, j, tmpindex
		n = size(A,1)
		if(present(Aindex)) then
			allocate(Aindex(n))
			do i = 1, n
				Aindex(i) = i
			enddo
		endif
		if(present(max_num)) then
			m = max_num
		else
			m = n
		endif
		do i = 1, m
		  do j = n, i+1, -1
		      if(A(j-1) > A(j)) then
		      	tmp = A(j-1)
		      	A(j-1) = A(j)
		      	A(j) = tmp
		      	if(present(Aindex)) then
		      	  tmpindex = Aindex(j-1)
		      	  Aindex(j-1) = Aindex(j)
		      	  Aindex(j) = tmpindex
		      	endif
		      endif
		   enddo
		 enddo
	end subroutine Bubble_1d

  !------------------------------------------
  ! 1d Bubble B
  !------------------------------------------
	!Buuble sorting: A is sorted. orderA is sorted along with A
	subroutine Bubble_1dB(A,orderA,nA,max_num)
		real(dl) :: A(nA)
		integer :: orderA(nA), nA
		integer, optional, intent(in) :: max_num
		real(dl) :: tmp
		integer :: n, m, i, j, tmpindex
		if(present(max_num)) then
			if(max_num > nA) then
				print *, 'ERROR (Bubble_1dB)! max_num shall not be larger than nA: ', max_num, nA
				stop
			endif
			m = max_num
		else
			m = nA
		endif
		do i = 1, m
		  do j = nA, i+1, -1
		      if(A(j-1) > A(j)) then
		      	tmp = A(j-1)
		      	A(j-1) = A(j)
		      	A(j) = tmp
		      	tmpindex = orderA(j-1)
		      	orderA(j-1) = orderA(j)
		      	orderA(j) = tmpindex
		      endif
		   enddo
		 enddo
	end subroutine Bubble_1dB


  !------------------------------------------
  ! 2d Bubble
  !------------------------------------------
	!Buuble sorting: 2d real(dl) array. B is the reference array.
	subroutine Bubble_2d(A,B,max_num)
		real(dl) :: A(:,:), B(:), tmpB
		real(dl), allocatable :: tmpA(:,:) 
		integer :: n, m, i, j
		integer, allocatable :: index_array(:)
		integer, optional, intent(in) :: max_num		
		! if max_num is given, only the first max_num elements are sorted
		n = size(A,2)
		m = size(A,1)
		if(size(B) .ne. n) then
			print*, 'ERROR! size mismatches. size of A, B = ', n, size(B)
			stop
		endif
		if(present(max_num)) then
			call Bubble_1d(B, index_array, max_num)
		else
			call Bubble_1d(B, index_array)
		endif
		if(present(max_num)) then
			allocate(tmpA(m,max_num))
			do i = 1, max_num
				tmpA(1:m,i) = A(1:m,index_array(i))
			enddo
			do i = 1, max_num
				A(1:m,i) = tmpA(1:m,i)
			enddo
			deallocate(tmpA)
			return
		else
			allocate(tmpA(m,n))
			do i = 1, n
				tmpA(1:m,i) = A(1:m,index_array(i))
			enddo
			A = tmpA
			deallocate(tmpA)
			return
		endif
	end subroutine Bubble_2d	
  	
  !---------------------------------------------------------------
  ! Quick Sorting 
  !---------------------------------------------------------------
	recursive subroutine QSort(A,nA)
		! DUMMY ARGUMENTS
		integer, intent(in) :: nA
		type (group), dimension(nA), intent(in out) :: A
 
		! LOCAL VARIABLES
		integer :: left, right
		real(dl) :: random
		real(dl) :: pivot
		type (group) :: temp
		integer :: marker
 
		if (nA > 1) then
			call random_number(random)
		        pivot = A(int(random*real(nA-1))+1)%value   ! random pivor (not best performance, but avoids worst-case)
        		left = 0
        		right = nA + 1
 
        		do while (left < right)
        		    right = right - 1
        		    do while (A(right)%value > pivot)
        		        right = right - 1
        		    end do
        		    left = left + 1
        		    do while (A(left)%value < pivot)
        		        left = left + 1
        		    end do
        		    if (left < right) then
        		        temp = A(left)
        		        A(left) = A(right)
        		        A(right) = temp
        		    end if
        		end do
 	
        		if (left == right) then
        		    marker = left + 1
        		else
        		    marker = left
        		end if
 	
        		call QSort(A(:marker-1),marker-1)
        		call QSort(A(marker:),nA-marker+1)
 	
		    end if
 
	end subroutine QSort 		
  		
  !---------------------------------------------------------------
  ! Quick Sorting 
  !---------------------------------------------------------------
	recursive subroutine QSort2(A,orderA,nA)
		! DUMMY ARGUMENTS
		real(dl) :: A(nA)
 		integer :: orderA(nA), nA
		! LOCAL VARIABLES
		integer :: left, right
		real(dl) :: random
		real(dl) :: pivot, tmpA, tmpi
		integer :: marker
 
! 		nA = size(A)
 !		if(size(orderA) .ne. nA) then
 !			print *, 'ERROR (Qsort2)! Size of order != size of A:', nA, size(orderA)
 !			stop
 !		endif
 
		if (nA > 1) then
			call random_number(random)
		        pivot = A(int(random*real(nA-1))+1)   ! random pivor (not best performance, but avoids worst-case)
        		left = 0
        		right = nA + 1
 
        		do while (left < right)
        		    right = right - 1
        		    do while (A(right) > pivot)
        		        right = right - 1
        		    end do
        		    left = left + 1
        		    do while (A(left) < pivot)
        		        left = left + 1
        		    end do
        		    if (left < right) then
        		        tmpA = A(left)
        		        A(left) = A(right)
        		        A(right) = tmpA
        		        tmpi = orderA(left)
        		        orderA(left) = orderA(right)
        		        orderA(right) = tmpi
        		    end if
        		end do
 	
        		if (left == right) then
        		    marker = left + 1
        		else
        		    marker = left
        		end if
 	
        		call QSort2(A(:marker-1),orderA(:marker-1),marker-1)
        		call QSort2(A(marker:),orderA(marker:),nA-marker+1)
 	
		    end if
 
	end subroutine QSort2 	
	
  !---------------------------------------------------------------
  ! Quick Sorting 3 (no index array)
  !---------------------------------------------------------------
	recursive subroutine QSort3(A,nA)
		! DUMMY ARGUMENTS
		real(dl) :: A(nA)
 		integer :: nA
		! LOCAL VARIABLES
		integer :: left, right
		real(dl) :: random
		real(dl) :: pivot, tmpA
		integer :: marker
 
		if (nA > 1) then
			call random_number(random)
		        pivot = A(int(random*real(nA-1))+1)   ! random pivor (not best performance, but avoids worst-case)
        		left = 0
        		right = nA + 1
 
        		do while (left < right)
        		    right = right - 1
        		    do while (A(right) > pivot)
        		        right = right - 1
        		    end do
        		    left = left + 1
        		    do while (A(left) < pivot)
        		        left = left + 1
        		    end do
        		    if (left < right) then
        		        tmpA = A(left)
        		        A(left) = A(right)
        		        A(right) = tmpA
        		    end if
        		end do
 	
        		if (left == right) then
        		    marker = left + 1
        		else
        		    marker = left
        		end if
 	
        		call QSort3(A(:marker-1),marker-1)
        		call QSort3(A(marker:),nA-marker+1)
 	
		    end if
	end subroutine QSort3 		
  		
  !! search within array  		
  !---------------------------------------------------------------
  ! Find out min/max in an array
  !---------------------------------------------------------------  		
  	subroutine find_min_max(A,nA,Amin,Amax)
  		real(dl) :: A(nA), Amin, Amax, x
  		integer :: nA, i
  		Amin = A(1)
  		Amax = A(1)
  		do i = 2, nA
  			x = A(i)
  			Amin = min(Amin, x)
  			Amax = max(Amax, x)
		enddo
	end subroutine find_min_max

  !---------------------------------------------------------------
  ! Count out how many elements are smaller than x
  !---------------------------------------------------------------  		
  	subroutine count_num_smaller(A,nA,x,num)
  		real(dl) :: A(nA), x
  		integer :: nA, i, num
		num = 0
  		do i = 1, nA
  			if(A(i) .le. x) num = num+1
		enddo
	end subroutine count_num_smaller
  !---------------------------------------------------------------
  ! Count out how many elements are larger than x
  !---------------------------------------------------------------  		
  	subroutine count_num_larger(A,nA,x,num)
  		real(dl) :: A(nA), x
  		integer :: nA, i, num
		num = 0
  		do i = 1, nA
  			if(A(i) .ge. x) num = num+1
		enddo
	end subroutine count_num_larger


  !---------------------------------------------------------------
  ! recursive subroutine find location with given ratio.
  !  not efficient...
  !---------------------------------------------------------------
	recursive subroutine findratloc(A,nA,numsm,numtol,xtry,ndepth)
		! DUMMy
		real(dl) :: A(nA), xtry
		integer :: nA, numsm, numtol
		! LOCAL
		real(dl) ::  Amin, Amax
		real(dl), allocatable :: B(:), C(:)
		integer :: i, j, k, num, ndepth

		if(nA .le. 5) then
!			print *, 'Bubble sort...'
			call bubble_1d(A,max_num=numsm)
			xtry = A(numsm)
			return
		endif
		
		call find_min_max(A, nA, Amin, Amax)
		
		xtry = Amin + ( real(numsm)/real(nA) ) * (Amax-Amin)
		num = 0
		j = 0
		k = 0
		allocate(B(nA),C(nA))
		do i = 1, nA
			if(A(i) .le. xtry) then
				num = num+1
				j = j + 1
				B(j) = A(i)
			else
				k = k + 1
				C(k) = A(i)
			endif
		enddo
		
!		print *, 'Amin, Amax, nA, numsm, xtry, num='
!s		write(*,'(5x,2e14.7,1x,i7,1x,i7,e14.7,1x,i7)'),  Amin, Amax, nA, numsm, xtry, num
		if(abs(num-numsm) .le. numtol) then
			deallocate(B,C)
			return
		elseif(num > numsm) then
			deallocate(C)
			call findratloc(B(1:j),j,numsm,numtol,xtry,ndepth+1)
		else
			deallocate(B)
			call findratloc(C(1:k),k,numsm-num,numtol,xtry,ndepth+1)
		endif
	end subroutine findratloc

	
  !---------------------------------------------------------------
  ! Count out how many elements are smaller than x
  !---------------------------------------------------------------  		
  	subroutine ltquan2(A,nA,numsm,dividecri,gvtol)
		! DUMMY
  		real(dl), intent(in) :: A(nA)
  		integer, intent(in) :: nA, numsm
  		real(dl), intent(out) :: dividecri
		integer, intent(in), optional :: gvtol
  		! LOCAL
  		real(dl) :: x, tol, x1, x2, xtry
  		real(dl), allocatable :: B(:)
  		integer, allocatable :: orderA(:), tmpint(:)
  		integer :: i, j, k, nstep, numtol, mark, num1, num2, numtry
  		logical :: print_test
  		
  		! consistency check
  		if(numsm .le. 0 .or. numsm .gt. nA) then
  			print *, 'ERROR (ltquan2)! numsm must be within 0, nA: ', numsm, nA
  			stop
  		endif
  		
  		! Setting of tolerance
  		if(present(gvtol)) then
  			if(gvtol < 0 .or. gvtol > nA) then
  				print *, 'ERROR (ltquan2)! tolerance # must be within ', 0, nA
  				stop
  			endif
  			numtol = gvtol
  		else
  			numtol = 1
  		endif
  		
  		call find_min_max(A,nA,x1,x2)

  		num1 = 1
  		num2 = nA
  		nstep = 1
  		do while(.true.)
  			if(nstep .gt. 100) then 
  			!Sometimes quickfind runs endless...(I donot know reason), and we have to switch to Qsort
  			!Code below are very hard understanding... Anyway, it works.
!  				print *, 'Warning (ltquan2)! Very large nstep! Switch to Qsort...'
  				allocate(B(nA),orderA(nA),tmpint(nA))
  				do i = 1, nA
  					orderA(i) = i
  					B(i) = A(i)
  				enddo
  				call Qsort2(B,orderA,nA)
 				dividecri = B(numsm)
  				return
  			endif
!  			
			if(nstep < 10) then
	  			xtry = x1 + ( dble(numsm-num1)/dble(num2-num1) ) * (x2-x1)
	  		else
	  			xtry = (x1 + x2)/2.0
	  		endif
	  		
  			call count_num_smaller(A,nA,xtry,numtry)
  			if(abs(numsm-numtry).le.numtol) then
  				dividecri = xtry
  				return  			
  			elseif(numtry .lt. numsm) then
  				x1 = xtry; num1 = numtry
  			else
  				x2 = xtry; num2 = numtry
  			endif
  			nstep = nstep+1
  		enddo
  		
  		print *, 'ERROR (ltlablist2)! failed: ', numsm, xtry, numtry, x1, num1, x2, num2
  		stop
	end subroutine ltquan2


  !---------------------------------------------------------------
  ! Count out how many elements are smaller than x
  !---------------------------------------------------------------  		
  	subroutine ltlablist1(A,nA,numsm,numtolin,lablist,splist,smlablist,Aminout,Amaxout,markin)
		! DUMMY
  		real(dl), intent(in) :: A(nA)
  		integer, intent(in) :: nA, numsm
  		integer, optional, intent(in) :: numtolin, markin
  		integer, optional, intent(out) :: lablist(nA) 
  		real(dl), optional, intent(out) :: Aminout, Amaxout
  		real(dl), optional, intent(out) :: splist(numsm)
  		integer, optional, intent(out) :: smlablist(numsm)
  		! LOCAL
  		real(dl) :: x, tol, x1, x2, xtry
  		real(dl), allocatable :: tmp(:)
  		integer :: i, j, k, nstep, numtol, mark, num1, num2, numtry
  		logical :: print_test
  		
  		if(present(numtolin)) then
  			numtol = numtolin
  		else
  			numtol = int(nA*0.0001)
  		endif	
  		
  		if(present(splist) .or. present(smlablist)) numtol = 0 ! must be exact if exist splist
  		
  		! elements satsifying the condition are marked by mark
  		if(present(markin)) then
  			mark = markin
  		else 
  			mark = 1 
  		endif 
  		
  		if(numsm .le. 0 .or. numsm .gt. nA) then
  			print *, 'ERROR (ltlablist1)! numsm must be within 0, nA: ', numsm, nA
  			stop
  		endif
  		
  		call find_min_max(A,nA,x1,x2)
 		if(present(Aminout)) Aminout = x1  		
 		if(present(Amaxout)) Amaxout = x2

  		num1 = 1
  		num2 = nA
  		nstep = 1
  		do while(num1 .ne. num2)
!  			
!			if(num2 - num1 .le. nA*0.01) then
!				j = 0
!				allocate(tmp(num2-num1+2))
!				do i = 1, nA
!					x = A(i)
!					if(x .ge. x1 .and. x .le. x2) then
!						j = j+1
!						if(j .gt. num2-num1) cycle 
!						tmp(j) = x
!					endif
!				enddo
!				call ltlablist1(tmp(1:j),j,numsm)... ! Recursive part Not Finished Yet
!			endif

			if(nstep < 10) then
	  			xtry = x1 + ( dble(numsm-num1)/dble(num2-num1) ) * (x2-x1)
	  		else
	  			xtry = (x1 + x2)/2.0
	  		endif
	  		
  			call count_num_smaller(A,nA,xtry,numtry)
!  			write(*,'(A,i7,i7,e14.7,i7,i7,i4)'),  'num1, num2, xtry, numtry, diff, numtol=', &
!  				num1, num2, xtry, numtry, abs(num-numtry), numtol
  			if(abs(numsm-numtry) .le. numtol) goto 100
  			
  			if(numtry .lt. numsm) then
  				x1 = xtry; num1 = numtry
  			else
  				x2 = xtry; num2 = numtry
  			endif
  			nstep = nstep+1
  			if(nstep .gt. 100) then
  				print *, 'Warning! Very large nstep!'
  				stop
  			endif
  		enddo
  		
  		print *, 'ERROR (ltlablist1)! failed: ', numsm, xtry, numtry, x1, num1, x2, num2
  		stop
  		
100		continue

		if(.not.present(lablist).and..not.present(splist)) return

		j = 0
		do i = 1, nA
			if(A(i) .le. xtry) then
				if(present(lablist)) lablist(i) = mark
				j = j + 1
				if(present(splist)) then
					splist(j) = A(i)
!					print *, j, splist(j)
				endif
				if(present(smlablist)) smlablist(j) = i
			endif
		enddo
	end subroutine ltlablist1

  !---------------------------------------------------------------
  ! Count out how many elements are smaller than x
  !---------------------------------------------------------------  		
  	subroutine ltlablist2(A,nA,numsm,splist,smlablist)
		! DUMMY
  		real(dl), intent(in) :: A(nA)
  		integer, intent(in) :: nA, numsm
  		real(dl), intent(out) :: splist(numsm)
  		integer, intent(out) :: smlablist(numsm)
  		! LOCAL
  		real(dl) :: x, tol, x1, x2, xtry
  		real(dl), allocatable :: B(:)
  		integer, allocatable :: orderA(:), tmpint(:)
  		integer :: i, j, k, nstep, numtol, mark, num1, num2, numtry
  		logical :: print_test
  		
  		if(numsm .le. 0 .or. numsm .gt. nA) then
  			print *, 'ERROR (Warning2)! numsm must be within 0, nA: ', numsm, nA
  			stop
  		endif
  		
  		call find_min_max(A,nA,x1,x2)

  		num1 = 1
  		num2 = nA
  		nstep = 1
  		do while(.true.)
  			if(nstep .gt. 100) then 
  			!Sometimes quickfind runs endless...(I donot know reason), and we have to switch to Qsort
  			!Code below are very hard understanding... Anyway, it works.
!  				print *, 'Warning (ltlablist2)! Very large nstep! Switch to Qsort...'
  				allocate(B(nA),orderA(nA),tmpint(nA))
  				do i = 1, nA
  					orderA(i) = i
  					B(i) = A(i)
  				enddo
  				call Qsort2(B,orderA,nA)
!  				smlablist(1:numsm) = orderA(1:numsm)
  				tmpint = 0
  				do i = 1, numsm
  					tmpint(orderA(i)) = 1
  				enddo
				j = 0
				do i = 1, nA
					if(tmpint(i) .eq. 1) then
						j = j + 1
						splist(j) = A(i)
						smlablist(j) = i
					endif
				enddo
  				return
  			endif
!  			
!			if(num2 - num1 .le. nA*0.01) then
!				j = 0
!				allocate(tmp(num2-num1+2))
!				do i = 1, nA
!					x = A(i)
!					if(x .ge. x1 .and. x .le. x2) then
!						j = j+1
!						if(j .gt. num2-num1) cycle 
!						tmp(j) = x
!					endif
!				enddo
!				call ltlablist1(tmp(1:j),j,numsm)... ! Recursive part Not Finished Yet
!			endif

			if(nstep < 10) then
	  			xtry = x1 + ( dble(numsm-num1)/dble(num2-num1) ) * (x2-x1)
	  		else
	  			xtry = (x1 + x2)/2.0
	  		endif
	  		
  			call count_num_smaller(A,nA,xtry,numtry)
!  			write(*,'(A,i7,i7,e14.7,i7,i7,i4)'),  'num1, num2, xtry, numtry, diff, numtol=', &
!  				num1, num2, xtry, numtry, abs(num-numtry), numtol
  			if(numsm .eq. numtry) goto 100
  			
  			if(numtry .lt. numsm) then
  				x1 = xtry; num1 = numtry
  			else
  				x2 = xtry; num2 = numtry
  			endif
  			nstep = nstep+1
  		enddo
  		
  		print *, 'ERROR (ltlablist2)! failed: ', numsm, xtry, numtry, x1, num1, x2, num2
  		stop
  		
100		continue

		j = 0
		do i = 1, nA
			if(A(i) .le. xtry) then
				j = j + 1
				splist(j) = A(i)
				smlablist(j) = i
			endif
		enddo
		if(j .ne. numsm) then
			print *, 'ERROR (ltlablist2)!!! j, numsm ne: ', j, numsm
		endif
	end subroutine ltlablist2

  !---------------------------------------------------------------
  ! Count out how many elements are smaller than x
  !---------------------------------------------------------------  		
  	subroutine ltlablist3(A,nA,numsm,smlablist)
		! DUMMY
  		real(dl), intent(in) :: A(nA)
  		integer, intent(in) :: nA, numsm
  		integer, intent(out) :: smlablist(numsm)
  		! LOCAL
  		real(dl) :: x, tol, x1, x2, xtry
  		real(dl), allocatable :: B(:)
  		integer, allocatable :: orderA(:), tmpint(:)
  		integer :: i, j, k, nstep, numtol, mark, num1, num2, numtry
  		logical :: print_test
  		
  		if(numsm .le. 0 .or. numsm .gt. nA) then
  			print *, 'ERROR (ltlablist3)! numsm must be within 0, nA: ', numsm, nA
  			stop
  		endif
  		
  		call find_min_max(A,nA,x1,x2)

  		num1 = 1
  		num2 = nA
  		nstep = 1
  		do while(.true.)
  			if(nstep .gt. 100) then 
  			!Sometimes quickfind runs endless...(I donot know reason), and we have to switch to Qsort
  			!Code below are very hard understanding... Anyway, it works.
  				print *, 'Warning (ltlablist3)! Very large nstep! Switch to Qsort...'
  				allocate(B(nA),orderA(nA),tmpint(nA))
  				do i = 1, nA
  					orderA(i) = i
  					B(i) = A(i)
  				enddo
  				call Qsort2(B,orderA,nA)
!  				smlablist(1:numsm) = orderA(1:numsm)
  				tmpint = 0
  				do i = 1, numsm
  					tmpint(orderA(i)) = 1
  				enddo
				j = 0
				do i = 1, nA
					if(tmpint(i) .eq. 1) then
						j = j + 1
						smlablist(j) = i
					endif
				enddo
  				return
  			endif

			if(nstep < 10) then
	  			xtry = x1 + ( dble(numsm-num1)/dble(num2-num1) ) * (x2-x1)
	  		else
	  			xtry = (x1 + x2)/2.0
	  		endif
	  		
  			call count_num_smaller(A,nA,xtry,numtry)
!  			write(*,'(A,i7,i7,e14.7,i7,i7,i4)'),  'num1, num2, xtry, numtry, diff, numtol=', &
!  				num1, num2, xtry, numtry, abs(num-numtry), numtol
  			if(numsm .eq. numtry) goto 100
  			
  			if(numtry .lt. numsm) then
  				x1 = xtry; num1 = numtry
  			else
  				x2 = xtry; num2 = numtry
  			endif
  			nstep = nstep+1
  		enddo
  		
  		print *, 'ERROR (ltlablist3)! failed: ', numsm, xtry, numtry, x1, num1, x2, num2
  		stop
  		
100		continue

		j = 0
		do i = 1, nA
			if(A(i) .le. xtry) then
				j = j + 1
				smlablist(j) = i
			endif
		enddo
	end subroutine ltlablist3

  !---------------------------------------------------------------
  ! Count out how many elements are smaller than x
  !---------------------------------------------------------------  		
  	subroutine gelablist2(A,nA,numge,gplist,gelablist)
		! DUMMY
  		real(dl), intent(in) :: A(nA)
  		real(dl), intent(out):: gplist(numge)
  		integer, intent(in) :: nA, numge
  		integer, intent(out) :: gelablist(numge)
  		! LOCAL
  		real(dl) :: x, tol, x1, x2, xtry
  		real(dl), allocatable :: B(:)
  		integer, allocatable :: orderA(:), tmpint(:)
  		integer :: i, j, k, nstep, numtol, mark, num1, num2, numtry
  		logical :: print_test
  		
  		if(numge .le. 0 .or. numge .gt. nA) then
  			print *, 'ERROR (gelablist2)! numge must be within 0, nA: ', numge, nA
  			stop
  		endif
  		
  		call find_min_max(A,nA,x2,x1)

  		num1 = 1
  		num2 = nA
  		nstep = 1
  		do while(.true.)
  			if(nstep .gt. 100) then 
  			!Sometimes quickfind runs endless...(I donot know reason), and we have to switch to Qsort
  			!Code below are very hard understanding... Anyway, it works.
  				print *, 'Warning (gelablist2)! Very large nstep! Switch to Qsort...'
  				allocate(B(nA),orderA(nA),tmpint(nA))
  				do i = 1, nA
  					orderA(i) = i
  					B(i) = A(i)
  				enddo
  				call Qsort2(B,orderA,nA)
!  				smlablist(1:numge) = orderA(1:numge)
  				tmpint = 0
  				do i = nA-numge+1,nA
  					tmpint(orderA(i)) = 1
  				enddo
				j = 0
				do i = 1, nA
					if(tmpint(i) .eq. 1) then
						j = j + 1
						gelablist(j) = i
						gplist(j) = A(i)
					endif
				enddo
  				return
  			endif

			if(nstep < 10) then
	  			xtry = x1 + ( dble(numge-num1)/dble(num2-num1) ) * (x2-x1)
	  		else
	  			xtry = (x1 + x2)/2.0
	  		endif
	  		
  			call count_num_larger(A,nA,xtry,numtry)
!  			write(*,'(A,i7,i7,e14.7,i7,i7,i4)'),  'num1, num2, xtry, numtry, diff, numtol=', &
!  				num1, num2, xtry, numtry, abs(num-numtry), numtol
  			if(numge .eq. numtry) goto 100
  			
  			if(numtry .gt. numge) then
  				x1 = xtry; num1 = numtry
  			else
  				x2 = xtry; num2 = numtry
  			endif
  			nstep = nstep+1
  		enddo
  		
  		print *, 'ERROR (gelablist2)! failed: ', numge, xtry, numtry, x1, num1, x2, num2
  		stop
  		
100		continue

		j = 0
		do i = 1, nA
			if(A(i) .ge. xtry) then
				j = j + 1
				gelablist(j) = i
				gplist(j) = A(i)
			endif
		enddo
	end subroutine gelablist2	
	
  !---------------------------------------------------------------
  ! Count out how many elements are smaller than x
  !---------------------------------------------------------------  		
  	subroutine gelablist3(A,nA,numge,gelablist)
		! DUMMY
  		real(dl), intent(in) :: A(nA)
  		integer, intent(in) :: nA, numge
  		integer, intent(out) :: gelablist(numge)
  		! LOCAL
  		real(dl) :: x, tol, x1, x2, xtry
  		real(dl), allocatable :: B(:)
  		integer, allocatable :: orderA(:), tmpint(:)
  		integer :: i, j, k, nstep, numtol, mark, num1, num2, numtry
  		logical :: print_test
  		
  		if(numge .le. 0 .or. numge .gt. nA) then
  			print *, 'ERROR (gelablist3)! numge must be within 0, nA: ', numge, nA
  			stop
  		endif
  		
  		call find_min_max(A,nA,x2,x1)

  		num1 = 1
  		num2 = nA
  		nstep = 1
  		do while(.true.)
  			if(nstep .gt. 100) then 
  			!Sometimes quickfind runs endless...(I donot know reason), and we have to switch to Qsort
  			!Code below are very hard understanding... Anyway, it works.
  				print *, 'Warning (gelablist3)! Very large nstep! Switch to Qsort...'
  				allocate(B(nA),orderA(nA),tmpint(nA))
  				do i = 1, nA
  					orderA(i) = i
  					B(i) = A(i)
  				enddo
  				call Qsort2(B,orderA,nA)
!  				smlablist(1:numge) = orderA(1:numge)
  				tmpint = 0
  				do i = nA-numge+1,nA
  					tmpint(orderA(i)) = 1
  				enddo
				j = 0
				do i = 1, nA
					if(tmpint(i) .eq. 1) then
						j = j + 1
						gelablist(j) = i
					endif
				enddo
  				return
  			endif

			if(nstep < 10) then
	  			xtry = x1 + ( dble(numge-num1)/dble(num2-num1) ) * (x2-x1)
	  		else
	  			xtry = (x1 + x2)/2.0
	  		endif
	  		
  			call count_num_larger(A,nA,xtry,numtry)
!  			write(*,'(A,i7,i7,e14.7,i7,i7,i4)'),  'num1, num2, xtry, numtry, diff, numtol=', &
!  				num1, num2, xtry, numtry, abs(num-numtry), numtol
  			if(numge .eq. numtry) goto 100
  			
  			if(numtry .gt. numge) then
  				x1 = xtry; num1 = numtry
  			else
  				x2 = xtry; num2 = numtry
  			endif
  			nstep = nstep+1
  		enddo
  		
  		print *, 'ERROR (gelablist3)! failed: ', numge, xtry, numtry, x1, num1, x2, num2
  		stop
  		
100		continue

		j = 0
		do i = 1, nA
			if(A(i) .ge. xtry) then
				j = j + 1
				gelablist(j) = i
			endif
		enddo
	end subroutine gelablist3

  !---------------------------------------------------------------
  ! inversion of a matrix
  !--------------------------------------------------------------- 	
	subroutine nizhen(aa,b,n)
		! Dummy
		real(dl), intent(in) :: aa(n,n)
		integer, intent(in) :: n
		real(dl), intent(out) :: b(n,n)
		! Local
	 	integer :: i,j,k
	  	real(dl) :: a(n,n)
	  	a=aa
		b=0.0d0
	  	do i=1,n
	   		b(i,i)=1
	  	enddo
	  	do i=1,n
	   		b(i,:)=b(i,:)/a(i,i)
	   		a(i,i:n)=a(i,i:n)/a(i,i) 
	   		do j=i+1,n     
	     			do k=1,n   
	   				b(j,k)=b(j,k)-b(i,k)*a(j,i)
	 			enddo
			  	a(j,i:n)=a(j,i:n)-a(i,i:n)*a(j,i)
			enddo
		enddo
	 	do i=n,1,-1
	  	do j=i-1,1,-1
	   	do k=1,n
	    		b(j,k)=b(j,k)-b(i,k)*a(j,i)
	    	enddo
	 	enddo
	  	enddo
 	end subroutine nizhen
 	
 	  	
  !------------------------------------------
  ! Polynomial Regression:
  !  polynomial fitting to data points
  ! Y(i) = A(1) + A(2)*X(1) + A(3)*X(2)^2 + ... 
  !        + A(n+1)*X(n)^n
  !------------------------------------------  	
	subroutine poly_fit(X,Y,A,ndat,n)
		! Dummy
		real(dl), intent(in) :: X(ndat),Y(ndat)
		integer, intent(in) :: ndat,n
		real(dl), intent(out) :: A(n+1)
		! Local
		real(dl) :: CapX(ndat,n+1), CapMatA(n+1,n+1), CapMatB(n+1,n+1), CapMatC(n+1,ndat)
		integer :: i,j,k
		! CapX(i,j) = X(i) ^ (j-1)
		do i = 1, ndat
		do j = 1, n+1
			CapX(i,j) = X(i)**(j-1)
		enddo
		enddo
		! CapMatA = ( CapX^T CapX)
		do i = 1,n+1
		do j = 1,n+1
			CapMatA(i,j) = 0.0_dl
			do k = 1,ndat
				CapMatA(i,j) = CapMatA(i,j) + CapX(k,i)*CapX(k,j)
			enddo
		enddo
		enddo
		! CapMatB = ( CapX^T CapX )^(-1)
		call nizhen(CapMatA,CapMatB,n+1)
		! CapMatC = ( CapX^T CapX )^(-1) CapX^T
		do i = 1, n+1
		do j = 1, ndat
			CapMatC(i,j) = 0.0_dl
			do k = 1, n+1
				CapMatC(i,j) = CapMatC(i,j) + CapMatB(i,k)*CapX(j,k)
			enddo
		enddo
		enddo
		! A = ( CapX^T CapX )^(-1) CapX^T Y
		do i = 1, n+1
			A(i) = 0.0_dl
			do j = 1, ndat
				A(i) = A(i) + CapMatC(i,j)*Y(j)
			enddo
		enddo
  	end subroutine poly_fit

  !------------------------------------------
  ! Value of a n-th polynomial 
  !  at some value of x
  !------------------------------------------   	
  	real(dl) function poly(x,A,n)
  		! Dummy
  		real(dl), intent(in) :: x, A(n+1)
  		integer, intent(in) :: n
  		! local
  		integer :: i
  		poly = A(1)
  		do i = 1, n
  			poly = poly + A(i+1)*x**(i)
  		enddo
  	end function poly
  	
  !---------------------------------------------------------------
  ! The same as binned_quan, but keep equal volume rather than distance ! not reliable
  !---------------------------------------------------------------
  	subroutine eqvl_binned_quan_old(quan_val_list, quan_r_list, rmin, rmax, ndat, num_bins, &
  		quan_val_av_list, quan_val_er_list, r_av_list, quan_val_var_list)
  		! Dummy
  		real(dl), intent(in) :: quan_val_list(ndat), quan_r_list(ndat) 
  		real(dl), intent(in) :: rmin, rmax
  		integer, intent(in) :: ndat, num_bins
  		real(dl), intent(out) :: quan_val_av_list(num_bins), &
			quan_val_er_list(num_bins), r_av_list(num_bins),  quan_val_var_list(num_bins)
  		! Local
  		real(dl) :: r, quan_val_av,quan_val_var, quanval1,quanval2,quanvalmid, rmincube,totvol,deltavol
  		integer :: i, n, m, j, binned_i, numsmaller
  		integer :: binned_num_list(num_bins)
  		type(dy_array) :: binned_quan_val_list(num_bins), binned_r_list(num_bins)
  		
!  		n = size(quan_val_list)
!  		if(size(quan_r_list) .ne. n) then
!  			print *, 'ERROR! # of quan and r mismatches: ', n, size(quan_r_list)
!  			stop
!  		endif
  		
  		rmincube = rmin**3.0
  		totvol = rmax**3.0 - rmincube
  		deltavol = totvol / dble(num_bins)
  		
  		print *, 'A'
  		
		! how many elements in each bin
		binned_num_list = 0
  		do i = 1, n
  			r = quan_r_list(i)
  			binned_i = max(min(int((r**3.0-rmincube)/deltavol)+1,num_bins),1)
  			binned_num_list(binned_i) = binned_num_list(binned_i) + 1
  		enddo
  		
  		print *, 'B'
  		
		! which elements in which bin
  		do binned_i = 1, num_bins
	  		m = binned_num_list(binned_i) 
  			binned_quan_val_list(binned_i)%n = m
  			binned_r_list(binned_i)%n = m
  			allocate(binned_quan_val_list(binned_i)%array(m))
  			allocate(binned_r_list(binned_i)%array(m))
  		enddo
  		print *, 'C'
		binned_num_list = 0  			
  		do i = 1, n
  			r = quan_r_list(i)
  			binned_i = max(min(int((r**3.0-rmincube)/deltavol)+1,num_bins),1)
  			binned_num_list(binned_i) = binned_num_list(binned_i) + 1
			j = binned_num_list(binned_i)
			binned_quan_val_list(binned_i)%array(j) = quan_val_list(i)
			binned_r_list(binned_i)%array(j) = quan_r_list(i)
  		enddo
  		print *, 'D'
  		do binned_i = 1, num_bins
	  		m = binned_num_list(binned_i) 
	  		quan_val_av = sum(binned_quan_val_list(binned_i)%array) / (m+0.0)
	  		quan_val_var = 0.0
	  		do j = 1, m
	  			quan_val_var = quan_val_var + &
	  			(binned_quan_val_list(binned_i)%array(j)-quan_val_av)**2.0
	  		enddo
	  		quan_val_var = quan_val_var / (m-1.0) 
  			quan_val_av_list(binned_i) = quan_val_av
			quan_val_var_list(binned_i) = quan_val_var
  			quan_val_er_list(binned_i) = sqrt(quan_val_var) / sqrt(m-1.0)
  			r_av_list(binned_i) = sum(binned_r_list(binned_i)%array) / (m+0.0)
  		enddo 
  		print *, 'E'
  	end subroutine eqvl_binned_quan_old
  	
  !---------------------------------------------------------------
  ! Calculate Ra-Dec-Span-Rat (shape of elliptic when projecting a perfect sphere to ra/dec plan)
  !---------------------------------------------------------------  	
	subroutine init_RDSR()
		character(len=char_len) :: readinfile = &
			"/home/xiaodongli/software/LSSCode/data/DEC_vs_RAspan-to-DECspan.txt"
		integer :: i
		call count_line_number(readinfile, gb_RDSR_num)
		gb_RDSR_num = gb_RDSR_num-1
		if(allocated(gb_RDSR_decs)) deallocate(gb_RDSR_decs)
		if(allocated(gb_RDSR_rats)) deallocate(gb_RDSR_rats)
		allocate(gb_RDSR_decs(gb_RDSR_num),gb_RDSR_rats(gb_RDSR_num))
		open(unit=1298403,file=readinfile)
		read(1298403,*)
		do i = 1, gb_RDSR_num
			read(1298403,*) gb_RDSR_decs(i), gb_RDSR_rats(i)
		enddo
		close(1298403)
		gb_RDSR_decmin = minval(gb_RDSR_decs); gb_RDSR_decmax = maxval(gb_RDSR_decs)
		gb_RDSR_deltadec = (gb_RDSR_decmax-gb_RDSR_decmin) / (gb_RDSR_num-1.0)
!		print *, 'gb_RDSR_decmin, gb_RDSR_decmax, gb_RDSR_num = ', gb_RDSR_decmin, gb_RDSR_decmax, gb_RDSR_num
!		print *, 'gb_RDSR_deltadec = ', gb_RDSR_deltadec, gb_RDSR_decs(2)-gb_RDSR_decs(1)
!		do i = 1, gb_RDSR_num
!			print *, i, gb_RDSR_decs(i),gb_RDSR_rats(i)
!		enddo
	end subroutine init_RDSR
	
	real(dl) function get_RDSR(dec)
		! Dummy
		real(dl), intent(in) :: dec
		! Local
		integer :: idec
		real(dl) :: x1, f1, x2, f2, x3, f3
		if(dec<gb_RDSR_decmin.or.dec>gb_RDSR_decmax) then
			print *, "WARNING (get_RDSR):: dec out of range: ", dec, gb_RDSR_decmin, gb_RDSR_decmax
!			get_RDSR = 1.0; return
		endif
		idec = int((dec-gb_RDSR_decmin) / gb_RDSR_deltadec + 0.5)
		idec = min(idec,gb_RDSR_num-1); idec = max(idec,2)
		x1 = gb_RDSR_decs(idec-1); f1 = gb_RDSR_rats(idec-1);
		x2 = gb_RDSR_decs(idec);   f2 = gb_RDSR_rats(idec);
		x3 = gb_RDSR_decs(idec+1); f3 = gb_RDSR_rats(idec+1);
		get_RDSR =  intpl_vl(dec, x1, f1, x2, f2, x3, f3)
	end function get_RDSR
	
        !! functons used to do rotation
        subroutine RotateMatrixX(theta, Matrix)        
                real(dl) :: theta, Matrix(3,3)
                Matrix = 0.0_dl
                Matrix(1,1) = 1.0
                Matrix(2,2) = dcos(dble(theta))
                Matrix(2,3) = -dsin(dble(theta))
                Matrix(3,2) = dsin(dble(theta))
                Matrix(3,3) = dcos(dble(theta))
        end subroutine RotateMatrixX
        subroutine RotateMatrixY(theta, Matrix)        
                real(dl) :: theta, Matrix(3,3)
                Matrix = 0.0_dl
                Matrix(2,2) = 1.0
                Matrix(1,1) = dcos(dble(theta))
                Matrix(1,3) = -dsin(dble(theta))
                Matrix(3,1) = dsin(dble(theta))
                Matrix(3,3) = dcos(dble(theta))
         end subroutine RotateMatrixY
        subroutine RotateMatrixZ(theta, Matrix)        
                real(dl) :: theta, Matrix(3,3)
                Matrix = 0.0_dl
                Matrix(3,3) = 1.0
                Matrix(1,1) = dcos(dble(theta))
                Matrix(1,2) = -dsin(dble(theta))
                Matrix(2,1) = dsin(dble(theta))
                Matrix(2,2) = dcos(dble(theta))
        end subroutine RotateMatrixZ
        subroutine RotateMatrix(thetax,thetay,thetaz,Matrix)
                ! Dummy
                real(dl) :: thetax, thetay, thetaz, Matrix(3,3)
                ! Local
                real(dl) :: tmp1(3,3), tmp2(3,3), tmp3(3,3), Matrix0(3,3)
                integer :: i
                logical :: print_info = .false.
                call RotateMatrixX(thetax,tmp1)
                call RotateMatrixY(thetay,tmp2)
                call RotateMatrixZ(thetaz,tmp3)
!               Matrix = matmul(tmp2,tmp1)
                call matrix_product(tmp1,tmp2,Matrix0,3,3,3)
                !Matrix = matmul(Matrix,tmp3)
                call matrix_product(Matrix0,tmp3,Matrix,3,3,3)
                if(print_info) then
                        call print_matrix(tmp1,3,3)
                        call print_matrix(tmp2,3,3)
                        call print_matrix(tmp3,3,3)
                        print *
                        call print_matrix(Matrix0,3,3)
                        print *
                        call print_matrix(Matrix,3,3)
                endif
        end subroutine RotateMatrix
        !!! To be checked: subroutine rotates a point ...
        subroutine RotateM(Point, Matrix, NewPoint)
                ! Dummy
                real(dl) :: Point(3), Matrix(3,3), NewPoint(3)
                ! Local
                integer :: i,j
                do i = 1,3
                        NewPoint(i) = 0.0
                        do j = 1, 3
                                NewPoint(i) = NewPoint(i) + Point(j)*Matrix(j,i)
                        enddo
                enddo
        end subroutine RotateM
        subroutine print_matrix(tmp,n1,n2)
                integer :: n1, n2, i,j
                real(dl) :: tmp(n1,n2)
                do i = 1, n1
                        print *, real(tmp(i,1:n2))
                enddo
        end subroutine print_matrix
        subroutine matrix_product(A,B,AB,n1,n2,n3)
                ! Dummy
                integer :: n1, n2, n3
                real(dl) :: A(n1,n2), B(n2,n3), AB(n1,n3)
                ! Local
                integer :: i,j,k
                do i = 1, n1
                        do j = 1, n3
                                AB(i,j) = 0
                                do k = 1, n2
                                        AB(i,j) = AB(i,j) + A(i,k)*B(k,j)
                                enddo
                        enddo
                enddo
        end subroutine matrix_product

  !---------------------------------------------------------------        
  ! x a variable; xedges an array; find out ibin so that xedges(ibin)<x<xedges(ibin+1)
  !---------------------------------------------------------------
	integer function find_ibin(x,xedges,nedge)
        	! Dummy
        	integer, intent(in) :: nedge
        	real(dl), intent(in) :: x, xedges(nedge)
        	! Local
        	integer :: i
        	if(x<xedges(1)) then
        		find_ibin = 0; return
        	elseif(x>xedges(nedge)) then
        		find_ibin = nedge; return
        	else
        		do i = 2, nedge
        			if(x<xedges(i)) then
        				find_ibin = i-1; return
        			endif
        		enddo
        	endif
        end function find_ibin
	! the same as before, use fast 2 split method 
	! can be significant faster than find_ibin if nedge > 1000
	integer function find_ibin_2split(x,xedges,nedge)
        	! Dummy
        	integer, intent(in) :: nedge
        	real(dl), intent(in) :: x, xedges(nedge)
        	! Local
        	integer :: i1, i2, midi
        	if(x<xedges(1)) then
        		find_ibin_2split = 0; return
        	elseif(x>xedges(nedge)) then
        		find_ibin_2split = nedge; return
        	else
        		i1 = 1; i2 = nedge; 
        		do while(.true.)
        			midi = (i1+i2)/2
!				midi = i1 + (i2-i1) * (x-xedges(i1)) / (xedges(i2)-xedges(i1))
!				midi = max(min(midi,i2-1),i1+1)
        			if(abs(i2-i1).le.15) exit
        			if(x<xedges(midi)) then
        				i2 = midi; cycle
        			else
        				i1 = midi; cycle
        			endif
        		enddo
        		find_ibin_2split = i1-1+find_ibin(x,xedges(i1:i2),i2-i1+1)
        	endif
        end function find_ibin_2split
	integer function find_ibin_revorder(x,xedges,nedge) ! same as find_ibin, but array is in decending rather than arising order
        	! Dummy
        	integer, intent(in) :: nedge
        	real(dl), intent(in) :: x, xedges(nedge)
        	! Local
        	integer :: i
        	if(x>xedges(1)) then
        		find_ibin_revorder = 0; return
        	elseif(x<xedges(nedge)) then
        		find_ibin_revorder = nedge; return
        	else
        		do i = 2, nedge
        			if(x>xedges(i)) then
        				find_ibin_revorder = i-1; return
        			endif
        		enddo
        	endif
        end function find_ibin_revorder
        
        character(len=char_len) function get_omwstr(omegam,w)
        	real(dl) :: omegam, w
        	character(len=char_len) :: tmpstr1,tmpstr2
        	write(tmpstr1,'(f12.4)') omegam
        	write(tmpstr2,'(f12.4)') w
        	get_omwstr = 'om'//trim(adjustl(tmpstr1))//'_w'//trim(adjustl(tmpstr2))
        end function get_omwstr
end module LSS_tools
