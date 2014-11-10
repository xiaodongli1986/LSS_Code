!############################3333333333333333333333
!
! definitions, functions related with cosmology
! All functions in this subroutine are in double precision
!
!############################3333333333333333333333

module ap_cosmo_funs
use ap_tools
implicit none

!!! Fixed settings

	! number of interpolation
	integer,  parameter	:: de_num_intpl = 2500

	! step size of redshift: used to get r from r(z) through interpolation
	double precision,  parameter	:: de_basenumber = 1.0_dl + 1.0_dl/512.0_dl
	double precision,  parameter	:: de_logbasenumber = log(de_basenumber)

	! step size of r: used to get redshift from r through interpolation
	double precision, parameter :: de_minr = 0.0_dl, de_maxr = 6000.0_dl, de_deltar = (de_maxr-de_minr)/dble(de_num_intpl-1)
	
!!! Important variables

	! current cosmological parameter
	real(dl):: gb_omegam, gb_w, gb_h

	! maximal reshift the interpolation can reach
	double precision :: de_maxintplz

	!Common used array, saving the interpolating data.
	! Frequently used in many dark energy models, so put them here as common variables.
	! Array for redshift.
	double precision :: de_zdata(de_num_intpl)
	! Array for ez, defined as H(z)/H(0)
	double precision :: de_comovr_data(de_num_intpl)
	double precision :: de_getz_data(de_num_intpl)
	double precision :: de_gfz_data(de_num_intpl)

	logical :: cosmo_funs_inited = .false.


contains

  !------------------------------------------
  ! Initialization for cosmo_funs
  !------------------------------------------
	subroutine cosmo_funs_init(printinfo)
		! Dummy
		logical, intent(in) :: printinfo
		! Local
		integer :: i
		
		if(printinfo) then
			write(*,*) "  (cosmo_funs_init) Initializing cosmo_funs..."
			write(*,*) "                      Redshift interpolating numbers: ", de_num_intpl
		endif

		! calculating Ra-Dec-Span-Rate
		call init_RDSR()

		!Initialize the redshift data
		do i = 1, de_num_intpl
			de_zdata(i) = de_zi(i)
		enddo
		
		!Range of the interpolation (maximal redshift)
		de_maxintplz = de_zdata(de_num_intpl)
		if(printinfo) then
			write(*,'(A,f9.3)') "                       Maximal redshift in interpolating: ", de_maxintplz
		endif
				
		cosmo_funs_inited = .true.
	end subroutine cosmo_funs_init

  !------------------------------------------
  ! Constructing the comovr table
  !------------------------------------------
  	subroutine de_calc_comovr()
  		integer :: i, maxi,num_bin, iofz30
  		double precision :: zleft, gfleft, zright
		! comvr(z)
  		de_comovr_data(1) = 0
  		do i = 2, de_num_intpl
  			de_comovr_data(i) = de_comovr_data(i-1) + &
		  		simpson(inv_Hz,de_zdata(i-1),de_zdata(i),1)*const_c*gb_h
		enddo
		! z(comvr)
		de_getz_data(1) = 0.0
		do i = 2, de_num_intpl
			de_getz_data(i) = get_z(de_ri(i))
		enddo
		! gf(z)
		!!  for the highest redshift point
		iofz30 = de_iz(30.0d0)
		if(iofz30 < de_num_intpl) then
			maxi = iofz30
			do i = maxi,de_num_intpl
				de_gfz_data(i) = 1.0_dl
			enddo
		else
			maxi = de_num_intpl
			zleft   =  30.0d0
			gfleft  =  1.0d0
			zright  =  de_zi(de_num_intpl)
			num_bin =  int((zright-zleft)*8.0)
			de_gfz_data(maxi) = RK(dgfdz, zleft, gfleft, zright, num_bin)
		endif
		do i = maxi-1,1,-1
			zleft  = de_zi(i+1)
			gfleft = de_gfz_data(i+1)
			zright = de_zi(i)
			de_gfz_data(i) = RK(dgfdz, zleft, gfleft, zright, 1)
		enddo
	end subroutine de_calc_comovr
	
  !------------------------------------------
  ! get z from interploation
  !------------------------------------------
  	double precision function de_zfromintpl(r)
  		integer :: i, i1, i2
  		double precision :: r
  		i = de_ir(r)
  		call ilist(i,1,1,de_num_intpl,i1,i2)
  		de_zfromintpl = intpl_vl(r, de_ri(i1), de_getz_data(i1), &
  			de_ri(i1+1), de_getz_data(i1+1),&
  			de_ri(i2), de_getz_data(i2))
	end function de_zfromintpl
  !------------------------------------------
  ! get comv_r from interploation
  !------------------------------------------
  	double precision function de_get_comovr(z)
  		integer :: i, i1, i2
  		double precision :: z
  		i = de_iz(z)
  		call ilist(i,1,1,de_num_intpl,i1,i2)
  		de_get_comovr = intpl_vl(z, de_zdata(i1), de_comovr_data(i1), &
  			de_zdata(i1+1), de_comovr_data(i1+1),&
  			de_zdata(i2), de_comovr_data(i2))
	end function de_get_comovr
  !------------------------------------------
  ! get comv_r from interploation
  !------------------------------------------
  	double precision function de_gfz_intpl(z)
  		integer :: i, i1, i2
  		double precision :: z
  		i = de_iz(z)
  		call ilist(i,1,1,de_num_intpl,i1,i2)
  		de_gfz_intpl = intpl_vl(z, de_zdata(i1), de_gfz_data(i1), &
  			de_zdata(i1+1), de_gfz_data(i1+1),&
  			de_zdata(i2), de_gfz_data(i2))
	end function de_gfz_intpl
	

  !------------------------------------------
  ! Get z from index
  !------------------------------------------
	double precision function de_zi(i)
		integer :: i
		de_zi = de_basenumber**DBLE(i-1) - 1.0_dl
	end function de_zi

  !------------------------------------------
  ! Get index from z
  !------------------------------------------
	integer function de_iz(z)
		double precision :: z
		de_iz = Ceiling(log(1.0_dl+z)/de_logbasenumber + 1.0_dl)
		if(de_iz < 2) de_iz = 2
	end function de_iz
	
  !------------------------------------------
  ! Get r from index
  !------------------------------------------
	double precision function de_ri(i)
		integer :: i
		de_ri = de_minr + de_deltar*(i-1)
	end function de_ri

  !------------------------------------------
  ! Get index from r
  !------------------------------------------
	integer function de_ir(r)
		double precision :: r
		de_ir = int((r-de_minr)/de_deltar+0.5)
		if(de_ir < 2) de_ir = 2
		if(de_ir > de_num_intpl) de_ir = de_num_intpl - 1
	end function de_ir

  !------------------------------------------
  ! Hubble parameter in unit of km/s/Mpc
  !------------------------------------------	
  	double precision function Hz(z)
  		double precision :: z
  		Hz = 100.0*gb_h*&
  		sqrt(gb_omegam*(1.0+z)**3.0 + (1.0-gb_omegam)*(1.0+z)**(3.0*(1.0+gb_w)))
  	end function Hz
  	double precision function inv_Hz(z)
  		double precision :: z
  		inv_Hz = 100.0*gb_h*sqrt(gb_omegam*(1.0+z)**3.0 + (1.0-gb_omegam)*(1.0+z)**(3.0*(1.0+gb_w)))
  		inv_Hz = 1.0 / inv_Hz
  	end function inv_Hz
  	double precision function inv_Hz1pz(z) ! 1/[H(z)*(1+z)]: function to be integrated by t_age(z)
  		double precision :: z
  		inv_Hz1pz = 100.0*gb_h*sqrt(gb_omegam*(1.0+z)**3.0 + (1.0-gb_omegam)*(1.0+z)**(3.0*(1.0+gb_w)))
  		inv_Hz1pz = 1.0 / (inv_Hz1pz*(1.0+z))
  	end function inv_Hz1pz

  !------------------------------------------
  ! ratio of omega matter
  !------------------------------------------
	double precision function omegamz(z)
		! Dummy
		double precision :: z
		! Local
		double precision :: ezsq
		ezsq = gb_omegam*(1.0+z)**3.0 + (1.0-gb_omegam)*(1.0+z)**(3.0*(1.0+gb_w))
		omegamz = gb_omegam*(1.0+z)**3.0 / ezsq
	end function omegamz
  	
  !------------------------------------------
  ! derivative of growth factor df(z) / dz 
  !------------------------------------------
	double precision function dgfdz(z, gf)
		! Dummy
		double precision, intent(in) :: z,gf
		! Local
		double precision ::  Omegamz, dlnHdz, esq

		esq     =  gb_omegam*(1.0d0+z)**3.0d0 + (1.0d0-gb_omegam)*(1.0d0+z)**(3.0d0*(1.0d0+gb_w))
		dlnHdz  =  gb_omegam*(1.0d0+z)**3.0d0 + (1.0d0+gb_w)*(1.0d0-gb_omegam)*(1.0d0+z)**(3.0d0*(1.0d0+gb_w))
		dlnHdz  =  1.5d0*dlnHdz/(1.0d0+z)/esq
		Omegamz =  gb_omegam * (1.0d0+z)**3.0d0 / esq

		dgfdz   =  (1.0d0/(1.0d0+z)) * (gf**2.0d0 - ((1.0d0+z)*dlnHdz-2.0d0)*gf - 1.5d0*Omegamz)
	end function dgfdz
  	
  !------------------------------------------
  ! comoving r in unit of Mpc/h
  !------------------------------------------	
	double precision function comov_r(z)
  		double precision :: z
  		integer :: n
  		n = max(4*ceiling(z / 0.125),4)
  		comov_r = simpson(inv_Hz,0.0d0,z,n)*const_c*gb_h
  	end function comov_r

  !------------------------------------------
  ! t_age in unit of s * Mpc/km
  !------------------------------------------	
	double precision function t_age(z)
  		double precision :: z
  		integer :: n
  		n = max(4*ceiling(z / 0.125),4)
  		t_age = simpson(inv_Hz1pz,0.0d0,z,n)
  	end function t_age

  !------------------------------------------
  !   	 get z according to comoving_r 
  !------------------------------------------	
	double precision function get_z(gv_comov_r, gv_zl, gv_zr)
		double precision :: gv_comov_r, zl, zr
		double precision, optional :: gv_zl, gv_zr
		if(present(gv_zl)) then
			zl = gv_zl
			else
			zl = 0.0
		endif
		if(present(gv_zr)) then
			zr = gv_zr
			else
			zr = 10.0
		endif
		get_z = findroot(comov_r, gv_comov_r, zl, zr)
	end function get_z
end module ap_cosmo_funs

