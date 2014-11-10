!July 1: remove wfkp norm, masscut norm; remove gb_treatbd_in_norm, gb_normn_pford, gb_treatbd_in_norm
!June 25: adding gb_radecbd_mat to mark the ra/dec boundary

!####################################
!This module does smooth
!####################################
module ap_cell
use ap_settings_init

	implicit none

!!! Fixed Settings

	! Normalize particle mass by the local density
	logical :: gb_dodensitynorm = .true.
	integer, parameter :: gb_normn_nbin = 20, gb_normn_nnorm = 2
	integer, parameter :: gb_radecmat_size = 800 ! size of the radecmat size

	! Fixed grid range with given values
	logical,  public, parameter :: use_gbfixgridrange = .false. ! use the range from some given values
	real(dl), public, parameter :: gbfixgridxmin=logzero, gbfixgridxmax=-logzero, &
				       gbfixgridymin=logzero, gbfixgridymax=-logzero, &
				       gbfixgridzmin=logzero, gbfixgridzmax=-logzero
	
	! Maximal length of array to save index of halos in one cell
	integer, parameter :: max_incell_point_num = 20000  

!!! Adjustable & To be fixed Settings	
! All variables in the "Adjustable & To be fixed Settings" section must be assigned a value in the main**.f90!!!

	! Number density or mass/weight density
	logical, public :: gb_usenumdensity

!!! Important Variables
	! Lists of x,y,z,mass/weight
	integer :: gb_num_xyz_mass
	real(dl), allocatable :: gb_xyz_list(:,:), gb_mass_list(:), gb_r_list(:)
	
	! 3D Grid
	integer :: gb_n_cellx, gb_n_celly, gb_n_cellz
	real(dl) ::  gb_cellvol, gb_cellwidth, gbdeltax, gbdeltay, gbdeltaz, gb_cell_vol

	type :: cell
		integer :: numdata = 0, numran = 0
		integer, allocatable :: idatalist(:), iranlist(:)
		integer :: rhodrhoindex=-1 ! index of the rho/drho of this cell saved in the gb_rhodrholist
		real :: rhodrhos(4) = -1.0, weirat = -1.0, maxdist = -1.0, sumwei_acpt=-1.0, sumwei_all=-1.0
		integer :: insphere_numran
!		real :: saddlerho 
		! rho at the nearby saddle point.
		! Coordinate of the nearby saddle point for cell labeled (ix,iy,iz) &
		!   will be (gbgridxmin + ix*gbdeltax, gbgridymin + iy*gbdeltay, gbgridzmin + iz*gbdeltaz)
	end type
	
	! array containing which halo in which cell
	type(cell), allocatable :: gb_cell_mat(:,:,:)
	
	! 2D grid capturing information of ra/dec randoms...
	integer :: gb_n_cellra, gb_n_celldec
	real(dl) :: gbgridramin, gbgridramax, gbgriddecmin, gbgriddecmax, &
		gb_radecrancell_area, gbdeltara, gbdeltadec
	type :: radecrancell
		integer :: numradecran = 0
		integer, allocatable :: iradecranlist(:)
		real(dl) :: acpt_rat = 0.0
	end type
	type(radecrancell), allocatable :: gb_radecrancell_mat(:,:)
	
	! mark the ra/dec boundary 
	real(dl), allocatable :: gb_radecbdmat_raleft(:),gb_radecbdmat_raright(:), &
		gb_radecbdmat_decleft(:), gb_radecbdmat_decright(:)
	
contains


  !------------------------------------------
  ! clean up allocatable arrays
  !------------------------------------------
	subroutine smooth_clean_up()
		if(allocated(gb_xyz_list)) &
			deallocate(gb_xyz_list)
		if(allocated(gb_mass_list)) &
			deallocate(gb_mass_list)
		if(allocated(gb_cell_mat)) &
			deallocate(gb_cell_mat)
		if(allocated(gb_radecrancell_mat)) &
			deallocate(gb_radecrancell_mat)
		if(allocated(gb_r_list)) &
			deallocate(gb_r_list)	
	end subroutine

  !------------------------------------------
  ! initialize the xyz_mass array
  !------------------------------------------
	subroutine init_mult_lists(printinfo)
		! Dummy
		logical, intent(in) :: printinfo
		! Local
		integer :: i, j
		real(dl) :: x, y, z, ratio, noRSDratio, r, theta, phi

		if(printinfo) write(*,'(A,i4,i4)'), '   (init_mult_list) Init lists of xyz/r/mass...'
		call smooth_clean_up()
		
		allocate(gb_xyz_list(3,gb_numallpar), gb_mass_list(gb_numallpar), gb_r_list(gb_numallpar))
		gb_num_xyz_mass = gb_numallpar

		! initializing data points (positive mass halos)
		do i = 1, gb_numdata
			! Mass
			if(gb_usenumdensity) then
				gb_mass_list(i) = 1.0
			else
				gb_mass_list(i) = gb_datalist(i)%mass
			endif
			ratio = gb_datalist(i)%rat
			! Position		
			x = gb_datalist(i)%x*ratio
			y = gb_datalist(i)%y*ratio
			z = gb_datalist(i)%z*ratio
			r = gb_datalist(i)%r*ratio
			gb_r_list(i) = r 
			gb_xyz_list(1,i) = x
			gb_xyz_list(2,i) = y
			gb_xyz_list(3,i) = z
		enddo
		
		! initializing mask rans
		do j = 1, gb_numran
			i = gb_numdata + j
			! Mass: will not consider number density for mask rans
			gb_mass_list(i) = 1.0  ! negative mass for mask rans!!!
			! Position		
			gb_xyz_list(1,i) = gb_ranlist(j)%x
			gb_xyz_list(2,i) = gb_ranlist(j)%y
			gb_xyz_list(3,i) = gb_ranlist(j)%z
!			gb_randata(4:6,j)
			call getSC(gb_xyz_list(1,i),gb_xyz_list(2,i),gb_xyz_list(3,i),r,theta,phi)
			gb_r_list(i) = r
			if(abs(r-gb_ranlist(j)%r) > 0.1)  then
				print *, 'ERROR (init_mult_lists)!: Mismatch r:', r, gb_ranlist(j)%r; stop
			endif
		enddo

		if(gb_dodensitynorm) then
			if(printinfo) write(*,'(29x,A)'), 'Normalizing mass according to <rho>(r)...'
			! This converts M_i into M_i / rhobar(r_i):
			!  the masses of all halos are divided by the rhobar(r), 
			!  to remove the global radial change of density
			call normn_mlist(printinfo,gb_normn_nbin,gb_normn_nnorm)
		endif
		
	end subroutine init_mult_lists


  !------------------------------------------
  ! get the list of x,y,z at centers of cells
  !------------------------------------------
	subroutine get_cell_pos_list(pos_list)
		real(dl), allocatable :: pos_list(:,:)
		real(dl) :: x,y,z
		integer :: ix, iy, iz, i
		allocate(pos_list(3,gb_n_cellx*gb_n_celly*gb_n_cellz))
		i = 1
		do ix = 1, gb_n_cellx
		do iy = 1, gb_n_celly
		do iz = 1, gb_n_cellz
			call cell_pos(ix,iy,iz,x,y,z)
			pos_list(1:3,i) = (/x,y,z/)
			i = i + 1
		enddo
		enddo
		enddo
	end subroutine get_cell_pos_list
	
  !------------------------------------------
  ! initialize the 3d cell (point+random)
  !------------------------------------------
	subroutine do_cell_init(rl_num_in_x, printinfo, do_not_init_cellmat)
		! Dummy
		real(dl), intent(in) :: rl_num_in_x
		logical, optional, intent(in) :: printinfo, do_not_init_cellmat
		! Local
		integer :: numset, i, j, k, ix, iy, iz, num_in_x, numskipdata,numskipran
		integer(2), allocatable :: mark_gb_numdata(:,:,:), mark_gb_numran(:,:,:)
		real(dl) :: x,y,z
		integer :: nownum
		integer, parameter :: re_alo_num = 5
		real(dl) :: numpixel, numhasdata, maxdatanum, avgdatanum, numhasran, maxrannum, avgrannum, extradx=0.02!1.001
		
		if(printinfo) print *, '  (do_cell_init begin) initializing Grid of Cells.'
		if(use_gbfixgridrange) then
			if(printinfo) print *, '  (do_cell_init) Use fixed range of grids:'
			gbgridxmin = gbfixgridxmin; gbgridxmax = gbfixgridxmax
			gbgridymin = gbfixgridymin; gbgridymax = gbfixgridymax
			gbgridzmin = gbfixgridzmin; gbgridzmax = gbfixgridzmax
		else
			if(printinfo) &
				write(*,'(3x,A)') '(do_cell_init) Set range of grids from xyz range of data (under current cosmology):'
			gbgridxmin = gbxmindata-extradx; gbgridxmax = gbxmaxdata+extradx
			gbgridymin = gbymindata-extradx; gbgridymax = gbymaxdata+extradx
			gbgridzmin = gbzmindata-extradx; gbgridzmax = gbzmaxdata+extradx
		endif

		gb_cellwidth = (gbgridxmax - gbgridxmin) / rl_num_in_x
		if(printinfo) write(*,'(20x,A,f12.5)')  'Width of cell: ', gb_cellwidth

		gb_n_cellx = ceiling( (gbgridxmax-gbgridxmin)/gb_cellwidth)
		gb_n_celly = ceiling( (gbgridymax-gbgridymin)/gb_cellwidth)
		gb_n_cellz = ceiling( (gbgridzmax-gbgridzmin)/gb_cellwidth)

		gbdeltax = gb_cellwidth; gbdeltay = gb_cellwidth; gbdeltaz = gb_cellwidth

		gbgridxmax = gbgridxmin + gbdeltax*dble(gb_n_cellx)
		gbgridymax = gbgridymin + gbdeltay*dble(gb_n_celly)
		gbgridzmax = gbgridzmin + gbdeltaz*dble(gb_n_cellz)
		
		gbgridrmax = sqrt(gbgridxmax**2.0+gbgridymax**2.0+gbgridzmax**2.0)

		if(printinfo) then
			write(*,'(20x,A,f12.5,2x,f12.5,A)')  'Grid Range of x: (', gbgridxmin, gbgridxmax, ')'
			write(*,'(20x,A,f12.5,2x,f12.5,A)')  '              y: (', gbgridymin, gbgridymax, ')'
			write(*,'(20x,A,f12.5,2x,f12.5,A)')  '              z: (', gbgridzmin, gbgridzmax, ')'
		endif
		
		gb_cell_vol = gbdeltax * gbdeltay * gbdeltaz
		if(printinfo) then
			write(*,'(20x,A,f12.5)')  'Effective # of cell in x direction: ', real(rl_num_in_x)
			write(*,'(20x,A,i3,A,i3,A,i3,A,i11)')  'Cell-#        = ', &
				gb_n_cellx,'*',gb_n_celly,'*',gb_n_cellz,' = ',gb_n_cellx*gb_n_celly*gb_n_cellz
			write(*,'(20x,A,3(f10.5,A),f15.5)') 'Cell-size = ', real(gbdeltax),' *',real(gbdeltay), &
				' *',real(gbdeltaz), ' =', real(gb_cell_vol)
		endif

		if(present(do_not_init_cellmat)) then
			if(do_not_init_cellmat .eq. .true.) return
		endif
		
		!allocating gb_cell_mat 
		if(allocated(gb_cell_mat)) deallocate(gb_cell_mat)
		allocate(gb_cell_mat(gb_n_cellx,gb_n_celly,gb_n_cellz),mark_gb_numdata(gb_n_cellx,gb_n_celly,gb_n_cellz),&
			mark_gb_numran(gb_n_cellx,gb_n_celly,gb_n_cellz))	
		
		! Count how many data/random points in each cell	
		numskipdata=0;numskipran=0
		do i = 1, gb_num_xyz_mass
			x= gb_xyz_list(1,i); y= gb_xyz_list(2,i); z= gb_xyz_list(3,i)
			if(x.ge.gbgridxmax .or. y.ge.gbgridymax .or. z.ge.gbgridzmax &
			   .or. x.le.gbgridxmin .or. y.le.gbgridymin .or. z.le.gbgridzmin) then
			   	if(i.le.gb_numdata) then
			   		numskipdata = numskipdata+1
			   	else
			   		numskipran = numskipran+1
			   	endif
				cycle
			endif
			ix=max(min(int((x-gbgridxmin)/gbdeltax)+1,gb_n_cellx),1)
			iy=max(min(int((y-gbgridymin)/gbdeltay)+1,gb_n_celly),1)
			iz=max(min(int((z-gbgridzmin)/gbdeltaz)+1,gb_n_cellz),1)
			if(i.le.gb_numdata) then
				gb_cell_mat(ix,iy,iz)%numdata = gb_cell_mat(ix,iy,iz)%numdata + 1
			else
				gb_cell_mat(ix,iy,iz)%numran = gb_cell_mat(ix,iy,iz)%numran + 1
			endif
		enddo
		
		if(printinfo) then
			write(*,'(20x,A,2(i10,A,f6.3,A))')  'Data/Random points lying outside the grid: ', &
				numskipdata,' (',numskipdata/dble(gb_numdata)*100.0,'%), ', &
				numskipran ,' (',numskipran/dble(gb_numran)*100.0,'%)'
		endif
				
		! allocate cell in pixels having halos
		numpixel = gb_n_cellx*gb_n_celly*gb_n_cellz
		numhasdata = 0.0; maxdatanum = 0.0; avgdatanum = 0.0
		numhasran = 0.0;  maxrannum = 0.0; avgrannum = 0.0
		
		do ix=1,gb_n_cellx
		do iy=1,gb_n_celly
		do iz=1,gb_n_cellz
			if(gb_cell_mat(ix,iy,iz)%numdata > 0) then
				if(gb_cell_mat(ix,iy,iz)%numdata > max_incell_point_num) then
					print *, 'ERROR! incell data num overflows!'
					print *, 'ix,iy,iz,num = ', ix,iy,iz,gb_cell_mat(ix,iy,iz)%numdata
					stop
				endif
				allocate(gb_cell_mat(ix,iy,iz)%idatalist(gb_cell_mat(ix,iy,iz)%numdata))
				numhasdata = numhasdata + 1.0_dl
				avgdatanum = avgdatanum + gb_cell_mat(ix,iy,iz)%numdata
				maxdatanum = max(maxdatanum,dble(gb_cell_mat(ix,iy,iz)%numdata))
			endif
			if(gb_cell_mat(ix,iy,iz)%numran > 0) then
				if(gb_cell_mat(ix,iy,iz)%numran > max_incell_point_num) then
					print *, 'ERROR! incell random num overflows!'
					print *, 'ix,iy,iz,num = ', ix,iy,iz,gb_cell_mat(ix,iy,iz)%numran
					stop
				endif
				allocate(gb_cell_mat(ix,iy,iz)%iranlist(gb_cell_mat(ix,iy,iz)%numran))
				numhasran = numhasran + 1.0_dl
				avgrannum = avgrannum + gb_cell_mat(ix,iy,iz)%numran
				maxrannum = max(maxrannum,dble(gb_cell_mat(ix,iy,iz)%numran))
			endif
		enddo
		enddo
		enddo
		avgdatanum = avgdatanum / dble(numhasdata)
		avgrannum = avgrannum / dble(numhasran)
		
		mark_gb_numdata = 0
		mark_gb_numran = 0
		! saving data into data lists
		do i = 1, gb_num_xyz_mass
			x= gb_xyz_list(1,i); y= gb_xyz_list(2,i); z= gb_xyz_list(3,i)
			if(x.ge.gbgridxmax .or. y.ge.gbgridymax .or. z.ge.gbgridzmax &
			   .or. x.le.gbgridxmin .or. y.le.gbgridymin .or. z.le.gbgridzmin) then
				cycle
			endif
			ix=max(min(int((x-gbgridxmin)/gbdeltax)+1,gb_n_cellx),1)
			iy=max(min(int((y-gbgridymin)/gbdeltay)+1,gb_n_celly),1)
			iz=max(min(int((z-gbgridzmin)/gbdeltaz)+1,gb_n_cellz),1)

			if(i.le.gb_numdata) then
				if(gb_cell_mat(ix,iy,iz)%numdata > 0) then
					mark_gb_numdata(ix,iy,iz) = mark_gb_numdata(ix,iy,iz)+1
					gb_cell_mat(ix,iy,iz)%idatalist(mark_gb_numdata(ix,iy,iz)) = i
				else
					print *, 'ERROR (do_cell_init)! Cell expected to have data inside!:',&
						ix,iy,iz,gb_cell_mat(ix,iy,iz)%numdata 
					stop
				endif
			else
				if(gb_cell_mat(ix,iy,iz)%numran > 0) then
					mark_gb_numran(ix,iy,iz) = mark_gb_numran(ix,iy,iz)+1
					gb_cell_mat(ix,iy,iz)%iranlist(mark_gb_numran(ix,iy,iz)) = i
				else
					print *, 'ERROR (do_cell_init)! Cell expected to have random inside!:',&
						ix,iy,iz,gb_cell_mat(ix,iy,iz)%numran
					stop
				endif
			endif
		enddo
		
		! Check matching of halo #
		do ix=1,gb_n_cellx
		do iy=1,gb_n_celly
		do iz=1,gb_n_cellz
			if(mark_gb_numdata(ix,iy,iz) .ne. gb_cell_mat(ix,iy,iz)%numdata .or. &
				mark_gb_numran(ix,iy,iz) .ne. gb_cell_mat(ix,iy,iz)%numran) then
				print *, 'ERROR (do_cell_init)!!! # mismatch: '
				print *, 'ix,iy,iz, #1, #2 = ',ix,iy,iz,mark_gb_numdata(ix,iy,iz),gb_cell_mat(ix,iy,iz)%numdata,&
					mark_gb_numran(ix,iy,iz),gb_cell_mat(ix,iy,iz)%numran
				stop
			endif
		enddo
		enddo
		enddo
		deallocate(mark_gb_numdata,mark_gb_numran)
		
		if(printinfo) then
			write(*,'(18x,i7,A,f6.2,A,i6,A,f10.3)') int(numhasdata+0.5),'(',&
				real(numhasdata/numpixel)*100.0,'%) cells have    data    inside; maximal # is ',&
				 int(maxdatanum+0.5), ', mean num is ', real(avgdatanum)
			write(*,'(18x,i7,A,f6.2,A,i6,A,f10.3)') int(numhasran+0.5),'(',&
				real(numhasran/numpixel)*100.0,'%) cells have   random   inside; maximal # is ',&
				 int(maxrannum+0.5), ', mean num is ', real(avgrannum)
			print *, '  (do_cell_init done)'
		endif
		
		if(gb_i_datatype .eq. gb_radecran .or. gb_i_datatype .eq. gb_radec3dran) then
			!call do_radecrancell_init(int(2*rl_num_in_x+0.5),printinfo)
			call do_radecrancell_init(gb_radecmat_size,printinfo)
		endif
	end subroutine do_cell_init	


  !------------------------------------------
  ! initialize the 2d cell for ra/dec random
  !------------------------------------------
	subroutine do_radecrancell_init(num_in_x, printinfo)
		! Dummy
		integer, intent(in) :: num_in_x
		logical, intent(in) :: printinfo
		! Local
		integer :: i,j,k,ira,idec, numskipradecran
		integer, allocatable :: mark_gb_numradecran(:,:)
		real(dl) :: ra,dec, numhasradecran, maxradecrannum, avgradecrannum, numpixel, ddl,ddr,ddu,ddd,minimal_dradec
		character(len=char_len) :: tmpstr1
		character, allocatable :: radecprintmat(:,:)
		
		if(printinfo) print *, '  (do_radecrancell_init begin) initializing 2d Grid of Cells for ra/dec random.'

		gbdeltara = (gbramaxradecran - gbraminradecran) / dble(num_in_x);  ! width in ra direction
		gbdeltadec = gbdeltara ! width in dec direction; keep the same value as ra direction 
		
		! how many cells in ra/dec direction
		gb_n_cellra  = ceiling( (gbramaxradecran -gbraminradecran )/gbdeltara) 
		gb_n_celldec = ceiling( (gbdecmaxradecran-gbdecminradecran)/gbdeltadec)
		
		! Grid range, in ra/dec directions
		gbgridramin  = gbraminradecran
		gbgriddecmin = gbdecminradecran
		gbgridramax = gbgridramin + gbdeltara*dble(gb_n_cellra)
		gbgriddecmax = gbgriddecmin + gbdeltara*dble(gb_n_celldec)

		if(printinfo) then
			write(*,'(20x,A,f12.5,2x,f12.5,A)')  'Grid Range of ra:  (', gbgridramin, gbgridramax, ')'
			write(*,'(20x,A,f12.5,2x,f12.5,A)')  '              dec: (', gbgriddecmin, gbgriddecmax, ')'
		endif
		
		! print info
		gb_radecrancell_area = gbdeltara * gbdeltadec
		if(printinfo) then
			write(*,'(20x,A,i4,A,i4,A,i11)')  'Cell-#        = ', &
				gb_n_cellra,'*',gb_n_celldec,' = ',gb_n_cellra*gb_n_celldec
			write(*,'(20x,A,2(f10.5,A),f15.5)') 'Cell-size = ', real(gbdeltara),' *',real(gbdeltadec), &
				' =', real(gb_radecrancell_area)
		endif

		!allocating matrix
		if(allocated(gb_radecrancell_mat)) deallocate(gb_radecrancell_mat)
		allocate(gb_radecrancell_mat(gb_n_cellra,gb_n_celldec),mark_gb_numradecran(gb_n_cellra,gb_n_celldec))	
		
		! Count how many ra/dec ran points in each cell	
		numskipradecran=0
		do i = 1, gb_numradecran
			ra = gb_radecranlist(i)%ra; dec = gb_radecranlist(i)%dec; 
			if(ra.ge.gbgridramax .or. dec.ge.gbgriddecmax &
			   .or. ra.le.gbgridramin .or. dec.le.gbgriddecmin) then
		   		numskipradecran = numskipradecran+1
				cycle
			endif
			ira =max(min(int((ra -gbgridramin)/gbdeltara)+1,gb_n_cellra),1)
			idec=max(min(int((dec-gbgriddecmin)/gbdeltadec)+1,gb_n_celldec),1)
			gb_radecrancell_mat(ira,idec)%numradecran = gb_radecrancell_mat(ira,idec)%numradecran + 1
			if(gb_radecranlist(i)%acpt) then
				 gb_radecrancell_mat(ira,idec)%acpt_rat = gb_radecrancell_mat(ira,idec)%acpt_rat + 1.0
			endif
		enddo
		
		if(printinfo) then
			write(*,'(20x,A,1(i10,A,f6.3,A))')  'ra/dec random points lying outside the grid: ', &
				numskipradecran,' (',numskipradecran/dble(gb_numradecran)*100.0,'%) '
		endif
				
		! allocate cell in pixels having ra/dec ran
		numpixel = gb_n_cellra * gb_n_celldec
		numhasradecran = 0.0;  maxradecrannum = 0.0; avgradecrannum = 0.0
		
		do ira =1,gb_n_cellra
		do idec=1,gb_n_celldec
			if(gb_radecrancell_mat(ira,idec)%numradecran > 0) then
				gb_radecrancell_mat(ira,idec)%acpt_rat = &
					gb_radecrancell_mat(ira,idec)%acpt_rat / dble(gb_radecrancell_mat(ira,idec)%numradecran)
				if(gb_radecrancell_mat(ira,idec)%numradecran > max_incell_point_num) then
					print *, 'ERROR! incell data num overflows!'
					print *, 'ira,idec,num = ', ira,idec,gb_radecrancell_mat(ira,idec)%numradecran
					stop
				endif
				allocate(gb_radecrancell_mat(ira,idec)%iradecranlist(gb_radecrancell_mat(ira,idec)%numradecran))
				numhasradecran = numhasradecran + 1.0_dl
				avgradecrannum = avgradecrannum + gb_radecrancell_mat(ira,idec)%numradecran
				maxradecrannum = max(maxradecrannum,dble(gb_radecrancell_mat(ira,idec)%numradecran))
			endif
		enddo
		enddo
		avgradecrannum = avgradecrannum / dble(numhasradecran)
		
		mark_gb_numradecran = 0
		! saving ra/dec rans into radecran lists
		do i = 1, gb_numradecran
			ra = gb_radecranlist(i)%ra; dec = gb_radecranlist(i)%dec; 
			if(ra.ge.gbgridramax .or. dec.ge.gbgriddecmax &
			   .or. ra.le.gbgridramin .or. dec.le.gbgriddecmin) then
				cycle
			endif
			ira =max(min(int((ra -gbgridramin)/gbdeltara)+1,gb_n_cellra),1)
			idec=max(min(int((dec-gbgriddecmin)/gbdeltadec)+1,gb_n_celldec),1)

			if(gb_radecrancell_mat(ira,idec)%numradecran > 0) then
				mark_gb_numradecran(ira,idec) = mark_gb_numradecran(ira,idec)+1
				gb_radecrancell_mat(ira,idec)%iradecranlist(mark_gb_numradecran(ira,idec)) = i
			else
				print *, 'ERROR (do_radecrancell_init)! Cell expected to have data inside!:',&
					ira,idec,gb_radecrancell_mat(ira,idec)%numradecran 
				stop
			endif
		enddo

		! gb_radecbdmats
		
		if(allocated(gb_radecbdmat_raleft)) deallocate(gb_radecbdmat_raleft)
		if(allocated(gb_radecbdmat_raright)) deallocate(gb_radecbdmat_raright)
		if(allocated(gb_radecbdmat_decleft)) deallocate(gb_radecbdmat_decleft)
		if(allocated(gb_radecbdmat_decright)) deallocate(gb_radecbdmat_decright)
		allocate(gb_radecbdmat_raleft(gb_n_celldec),gb_radecbdmat_raright(gb_n_celldec), &
			gb_radecbdmat_decleft(gb_n_cellra), gb_radecbdmat_decright(gb_n_cellra))  
		
		! searching for the left bd: start from the very left and ends until acpt_rat != 0.0 (I use > 0.01)
		do idec = 1, gb_n_celldec
			! left ra boundary for this idec
			do ira = 1, gb_n_cellra
				if(gb_radecrancell_mat(ira,idec)%acpt_rat .gt. 0.001_dl) exit
			enddo
			ira=max(min(ira,gb_n_cellra),1)
			call radecrancell_radec(ira,idec,ra,dec)
			gb_radecbdmat_raleft(idec) = (ra-0.5*gbdeltara) + (1.0-gb_radecrancell_mat(ira,idec)%acpt_rat)*gbdeltara
			! right ra boundary for this idec
			do ira = gb_n_cellra, 1, -1
				if(gb_radecrancell_mat(ira,idec)%acpt_rat .gt. 0.001_dl) exit
			enddo
			ira=max(min(ira,gb_n_cellra),1)
			call radecrancell_radec(ira,idec,ra,dec)
			gb_radecbdmat_raright(idec) = (ra+0.5*gbdeltara) - (1.0-gb_radecrancell_mat(ira,idec)%acpt_rat)*gbdeltara
		enddo

		do ira = 1, gb_n_cellra
			! left dec boundary for this ira
			do idec = 1, gb_n_celldec
				if(gb_radecrancell_mat(ira,idec)%acpt_rat .gt. 0.001_dl) exit
			enddo
			idec=max(min(idec,gb_n_celldec),1)
			call radecrancell_radec(ira,idec,ra,dec)
			gb_radecbdmat_decleft(ira) = (dec-0.5*gbdeltadec) + (1.0-gb_radecrancell_mat(ira,idec)%acpt_rat)*gbdeltadec
			! right dec boundary for this ira
			do idec = gb_n_celldec, 1, -1
				if(gb_radecrancell_mat(ira,idec)%acpt_rat .gt. 0.001_dl) exit
			enddo
			idec=max(min(idec,gb_n_celldec),1)
			call radecrancell_radec(ira,idec,ra,dec)
			gb_radecbdmat_decright(ira) = (dec+0.5*gbdeltadec) - (1.0-gb_radecrancell_mat(ira,idec)%acpt_rat)*gbdeltadec
		enddo

		! consistency check
		do ira =1,gb_n_cellra
		do idec=1,gb_n_celldec
			if(gb_radecrancell_mat(ira,idec)%numradecran .ne. mark_gb_numradecran(ira,idec)) then
				print *, 'ERROR (do_radecrancell_init)!!! # mismatch: '
				print *, 'ira,idec, #1, #2 = ',&
					ira,idec,mark_gb_numradecran(ira,idec),gb_radecrancell_mat(ira,idec)%numradecran
				stop
			endif
			if(gb_radecrancell_mat(ira,idec)%numradecran > 0) then
				do i = 1, gb_radecrancell_mat(ira,idec)%numradecran
					if(gb_radecrancell_mat(ira,idec)%iradecranlist(i)<= 0 .or. &
						gb_radecrancell_mat(ira,idec)%iradecranlist(i) > gb_numradecran) then
						print *, 'ERROR (do_radecrancell_init)!!! index out of range: ', i, 0, gb_numradecran
						stop
					endif
				enddo
			endif
		enddo
		enddo
		
		deallocate(mark_gb_numradecran)
		
		! output info if required
		if(gb_outputinfo_radeccell) then
			if(printinfo) write(*,'(20x,A)')  'output information: _radeccell.txt; _rabd.txt; _decbd.txt; _insideps_bd4/8/12'
			open(unit=198,file=trim(adjustl(gb_suboutput_name))//'_radeccell.txt')
			write(198,*) '## fmt: ira, idec, ra, dec, #-ran, rat-of-accpeted-ran'
			do ira=1,gb_n_cellra
				do idec=1,gb_n_celldec
					call radecrancell_radec(ira,idec,ra,dec)
					write(198,'(2i4,2f10.5,i7,f10.7)') ira, idec, ra, dec, &
						gb_radecrancell_mat(ira,idec)%numradecran, &
						gb_radecrancell_mat(ira,idec)%acpt_rat
				enddo
			enddo
			close(198)
			! information of ra/dec boundary: given ra/dec, where is the boundary of dec/ra?
			open(unit=199,file=trim(adjustl(gb_suboutput_name))//'_rabd.txt')
			open(unit=200,file=trim(adjustl(gb_suboutput_name))//'_decbd.txt')
			write(199,*) '## fmt: dec, left ra bd, right ra bd'
			write(200,*) '## fmt: ra, left dec bd, right dec bd'
			do idec=1,gb_n_celldec
				call radecrancell_radec(1,idec,ra,dec)
				write(199,'(3e16.7)') dec, gb_radecbdmat_raleft(idec), gb_radecbdmat_raright(idec)
			enddo
			do ira=1,gb_n_cellra
				call radecrancell_radec(ira,1,ra,dec)
				write(200,'(3e16.7)') ra, gb_radecbdmat_decleft(ira), gb_radecbdmat_decright(ira)
			enddo
			close(199); close(200)
			
			do k = 4, 12, 4
				minimal_dradec = k
				write(tmpstr1, *) int(minimal_dradec+0.5)
				open(unit=199,file=trim(adjustl(gb_suboutput_name))//'_insideps_bd'//trim(adjustl(tmpstr1))//'.txt')
				do i = 1, 100
					do j = 1, 100
						ra = gbgridramin + (gbgridramax-gbgridramin)/99.0*(i-1)
						dec = gbgriddecmin + (gbgriddecmax-gbgriddecmin)/99.0*(j-1)
						call dist_to_radecbd(ra,dec,ddl,ddr,ddu,ddd)
						if(ddl>minimal_dradec .and. ddr>minimal_dradec .and. &
						   ddu>minimal_dradec .and. ddd>minimal_dradec) then
							write(199,*) ra, dec
						endif
					enddo
				enddo
				close(199)
			enddo
		endif

				
		if(printinfo) then
			write(*,'(18x,i7,A,f6.2,A,i6,A,f10.3)') int(numhasradecran+0.5),'(',&
				real(numhasradecran/numpixel)*100.0,'%) cells have ra/dec ran inside; maximal # is ',&
				 int(maxradecrannum+0.5), ', mean num is ', real(avgradecrannum)
				 
			allocate(radecprintmat(gb_n_cellra,gb_n_celldec))
			do ira = 1, gb_n_cellra
				do idec = 1, gb_n_celldec
					if(gb_radecrancell_mat(ira,idec)%acpt_rat<0.5) then
						radecprintmat(ira,idec) = '#'
					else
						radecprintmat(ira,idec) = ' '
					endif
				enddo
			enddo
!			do idec = 1, gb_n_celldec
!				write(*,'(<gb_n_cellra>A)'), radecprintmat(1:gb_n_cellra,idec)
!			enddo

			print *, '  (do_radecrancell_init done)'
		endif
	end subroutine do_radecrancell_init

  !------------------------------------------
  ! central position of the cell
  !------------------------------------------	
  	subroutine cell_pos(ix,iy,iz,x,y,z)
		integer, intent(in) :: ix, iy, iz
  		real(dl), intent(out) :: x,y,z
  		x = gbgridxmin + (ix-0.5)*gbdeltax
  		y = gbgridymin + (iy-0.5)*gbdeltay
  		z = gbgridzmin + (iz-0.5)*gbdeltaz
	end subroutine cell_pos
  !------------------------------------------
  ! central position of the cell
  !------------------------------------------	
  	subroutine radecrancell_radec(ira,idec,ra,dec)
		integer, intent(in) :: ira,idec
  		real(dl), intent(out) :: ra,dec
 ! 		if(ira>gb_n_cellra.or.ira<=0.or.idec>gb_n_celldec.or.idec<=0) then
! ! 			print *, 'WARNING (radecrancell_radec): outflow ira,idec = ', ira, idec
!  		endif
  		ra  = gbgridramin  + (ira-0.5)*gbdeltara
  		dec = gbgriddecmin + (idec-0.5)*gbdeltadec
	end subroutine radecrancell_radec
	
  !------------------------------------------
  ! cell index 
  !------------------------------------------	
  	subroutine cell_index(ix,iy,iz,x,y,z)
		integer, intent(out) :: ix, iy, iz
  		real(dl), intent(in) :: x,y,z
  		ix = int((x-gbgridxmin) / gbdeltax +1.0)
  		iy = int((y-gbgridymin) / gbdeltay +1.0)
  		iz = int((z-gbgridzmin) / gbdeltaz +1.0)
	end subroutine cell_index
  !------------------------------------------
  ! cell index 
  !------------------------------------------	
  	subroutine radecrancell_index(ira,idec,ra,dec)
		integer, intent(out) :: ira,idec
  		real(dl), intent(in) :: ra,dec
  		ira = int((ra-gbgridramin) / gbdeltara +1.0)
  		idec = int((dec-gbgriddecmin) / gbdeltadec +1.0)
	end subroutine radecrancell_index
  !------------------------------------------
  ! initialize the 2d cell for ra/dec random
  !------------------------------------------
	subroutine dist_to_radecbd(ra,dec,ddl,ddr,ddu,ddd)
		! local
		real(dl), intent(in) :: ra,dec
		real(dl), intent(out) :: ddl,ddr,ddu,ddd
		! dummy
		integer :: ira,idec, ral, rar, decl, decr
		
		call radecrancell_index(ira,idec,ra,dec)
		ral = gb_radecbdmat_raleft(idec)
		rar = gb_radecbdmat_raright(idec)
		decl = gb_radecbdmat_decleft(ira)
		decr = gb_radecbdmat_decright(ira)
		ddl = max(0.0_dl,ra-ral)
		ddr = max(0.0_dl,rar-ra)
		ddu = max(0.0_dl,dec-decl)
		ddd = max(0.0_dl,decr-dec)
	end subroutine dist_to_radecbd

  !------------------------------------------
  ! Find out the index of point with nearest redshift
  !  used for w_fkp !
  !------------------------------------------
	integer function find_red_index(redshift, zwarray, numz)
		! Dummy
		integer :: numz
		real(dl) :: redshift, zwarray(2,numz)
		! Local
		integer :: i1, i2, i
		real(dl) :: z1, z2, z
		
		i1 = 1; i2 = numz; 
		z1 = zwarray(1,i1); z2 = zwarray(1,i2)
		if(redshift<z1.or.redshift>z2) then
			print *, 'ERROR(find_redindex)!!!: redshift out of range: ', redshift, z1,z2; stop
		endif
		do while(abs(i1-i2)>3)
			i = (i1+i2)/2; z = zwarray(1,i)
			if(z>redshift) then
				i2 = i; z2 = z
			else
				i1 = i; z1 = z
			endif
		enddo
		find_red_index = i1
	end function find_red_index
	

  !------------------------------------------
  ! Reset the mass_list: normalized by mass density;
  !------------------------------------------
  			! This function converts M_i into M_i / <rho>(r_i):
  			!  the masses of all halos are divided by the rhobar(r), 
			!  to remove the global radial change of density
			! Will ignore everthing with mass smaller than 0...
	subroutine normn_mlist(printinfo, nbins, nsm)
		!Dummy
		integer, intent(in) :: nbins, nsm
		logical, intent(in) :: printinfo
		!Local
		real(dl) :: quan_av_list(nbins),quan_er_list(nbins),binned_r_list(nbins),quan_var_list(nbins),&
			quan_medval_list(nbins),numinbin(nbins),redgelist(nbins+1)
		real(dl), allocatable :: mass_list(:), r_list(:)
		integer, parameter :: polyorder = 3
		real(dl) :: r, rmin, rmax, deltar, polyrho, polycoef(polyorder+1), pow = 1.0, binnedrho, redshift, &
			red1, red2, w1, w2, weight, deltaw, masscutnorm_redmins(1000), masscutnorm_redmaxs(1000), &
			masscutnorm_info(10,1000), masscutnorm_weis(1000)
		!, r1,r2,rho1,rho2
		integer :: i_sm, i, j, ibin, test_i, num_wfkp, i1,i2, masscutnorm_nbin
		character(len=char_len) :: tmpstr1, tmpstr2, tmpstr3, &
			w_fkpfile = '/home/lixiaodong/SparseFilaments/data/input/DR12/nbar/DR12v1-CMASS-N.z-wfkp.txt'
		real(dl), allocatable :: z_wfkp(:,:)
		
		if(printinfo) then
			write(*,'(A,i4,A,i3,A,2f10.3,A)'), '   (normn_mlist) Estimating medium rho for ', nbins, &
				' bins. Polyorder = ', polyorder ,'. Distance range: ', real(gbrmindata), real(gbrmaxdata), '...'
		endif
		
		! Construct list of halos with positive mass, and only use them for the normalization
  		do i_sm = 1, nsm
			if(.true.) then
				if(allocated(mass_list)) deallocate(mass_list)
				if(allocated(r_list)) deallocate(r_list)
				allocate(mass_list(gb_numdata), r_list(gb_numdata))
				j = 1
				do i = 1, gb_numdata
					mass_list(j) = gb_mass_list(i)
					r_list(j) = gb_r_list(i)**pow
					j = j+1
				enddo
			endif

			call find_min_max(r_list, gb_numdata, rmin, rmax)
			deltar = (rmax - rmin)  / dble(nbins)
			do ibin = 1, nbins+1
				redgelist(ibin) = rmin + deltar * (ibin-1)
			enddo

  			numinbin(1:nbins) = 0
  			quan_av_list(1:nbins) = 0.0
  			binned_r_list(1:nbins) = 0.0
  			do i = 1, gb_numdata
  				r = r_list(i)
  				ibin = int((r-rmin) / deltar) + 1
  				ibin = min(ibin, nbins)
				ibin = max(ibin, 1)
  				numinbin(ibin) = numinbin(ibin)+1
  				quan_av_list(ibin) = quan_av_list(ibin)+mass_list(i)
  				binned_r_list(ibin) = binned_r_list(ibin)+r_list(i)**(1.0/pow)
  			enddo
  			do ibin = 1, nbins
  				quan_av_list(ibin) = quan_av_list(ibin) / dble(numinbin(ibin))
  				binned_r_list(ibin) = binned_r_list(ibin) / dble(numinbin(ibin))
  			enddo
!  			call eqvl_binned_quan(mass_list, r_list, rmin, rmax, gb_numdata, nbins, &
!		  		quan_av_list, quan_er_list, binned_r_list, numinbinlist= numinbin, redgelist = redgelist)
!		  	call poly_fit(binned_r_list,quan_medval_list,polycoef,nbins,polyorder)
!			write(*,'(A)'), '   (normn_mlist) ATTENTION: Using mean value rather than middle value to normalize the mass!'
			! <rho>(r) = MASS / Vol = <M_i>*# / Vol
			do i = 1, nbins
				quan_av_list(i) = quan_av_list(i) * numinbin(i) / vol_fun(redgelist(i),redgelist(i+1))
			enddo
		  	call poly_fit(binned_r_list,quan_av_list,polycoef,nbins,polyorder)
			if(printinfo.and..false.) then
				print *
				write(*,'(A,<nbins>(e14.7,","))')  'binned_r_list: ', binned_r_list
				write(*,'(A,<nbins>(e14.7,","))')  'quan_av_list:  ', quan_av_list
				print *
				write(*,'(19x,A,6x,A,10x,A,12x,A,6x,A,3x,A)') ' bin     #-in-bin     distance','redshift',&
					'   range of r ','density','density(fit)','rat(fit/orig)'
				do i = 1, nbins
					write(*,'(20x,i3,2x,i10,4(2x,f12.3),4x,e14.7,4x,e14.7,4x,f8.5)') &
						i,int(numinbin(i)),binned_r_list(i),&
						de_zfromintpl(dble(binned_r_list(i))), redgelist(i), redgelist(i+1), &
						quan_av_list(i), poly(binned_r_list(i),polycoef,polyorder), &
						poly(binned_r_list(i),polycoef,polyorder) / quan_av_list(i)
				enddo
			endif
			
			if(gb_outputinfo_nbar) then
			  	write(tmpstr1,*) i_sm
			  	tmpstr1 = '_nbar_fitted_nsm'//trim(adjustl(tmpstr1))//'.txt'
			  	write(*,'(19x,A,A)') 'output information:', trim(adjustl(tmpstr1))
				open(unit=123948,file=trim(adjustl(gb_suboutput_name))//trim(adjustl(tmpstr1)))
				write(123948,*) '# fmt: r, redshift, nbar, nbar(fitted)'
				do i = 1, nbins
				  	write(123948,'(4e15.7)') binned_r_list(i), de_zfromintpl(dble(binned_r_list(i))), &
					  	quan_av_list(i), poly(binned_r_list(i),polycoef,polyorder)
				enddo
				close(123948)
			endif

	  		do i = 1, gb_numdata
  				r = gb_r_list(i)
				polyrho = poly(r,polycoef,polyorder)
  				ibin = int((r-rmin) / deltar) + 1
  				ibin = min(ibin,nbins); ibin = max(ibin,1)
  				if(ibin.ge.nbins-1.or.ibin.le.2) then
  					binnedrho = quan_av_list(ibin)
  				else
					binnedrho = intpl_vl(dble(r), dble(binned_r_list(ibin-1)), dble(quan_av_list(ibin-1)), &
						dble(binned_r_list(ibin)), dble(quan_av_list(ibin)), &
						dble(binned_r_list(ibin+1)), dble(quan_av_list(ibin+1)))
				endif
!  					r1 = binned_r_list(ibin);r2 = binned_r_list(ibin+1)
!  					rho1 = quan_av_list(ibin);rho2 = quan_av_list(ibin+1)
!  					binnedrho = rho1+ (r-r1)/(r2-r1)*(rho2-rho1)
!					binnedrho = quan_av_list(ibin)
  				gb_mass_list(i) = gb_mass_list(i) / binnedrho !polyrho
  			enddo
  		enddo
  	end subroutine normn_mlist	

!!! Test subroutines

  !------------------------------------------
  ! used to test the readin data
  !------------------------------------------	
	subroutine test_grid()
		integer :: ix,iy,iz,i,idata,iran
		real(dl) :: x,y,z
		do while(.true.)
			print *, 'input ix,iy,iz. Range: 1 to ', gb_n_cellx,gb_n_celly,gb_n_cellz
			read(*,*) ix,iy,iz
			call cell_pos(ix,iy,iz,x,y,z)
			print *, 'center position: ', x,y,z
			print *, 'Range of x/y/z:'
			print *, '   ',x-gbdeltax/2.0,x+gbdeltax/2.0
			print *, '   ',y-gbdeltay/2.0,y+gbdeltay/2.0
			print *, '   ',z-gbdeltaz/2.0,z+gbdeltaz/2.0
			print *, 'incell data:'
			print *, '   #:   ',gb_cell_mat(ix,iy,iz)%numdata
			print *, '   list:',gb_cell_mat(ix,iy,iz)%idatalist(1:gb_cell_mat(ix,iy,iz)%numdata)
			do i = 1, gb_cell_mat(ix,iy,iz)%numdata
				idata = gb_cell_mat(ix,iy,iz)%idatalist(i)
				print *, '      idata ',idata, ':',real(gb_xyz_list(1:3,idata))
			enddo
			print *, 'incell random:'
			print *, '   #:   ',gb_cell_mat(ix,iy,iz)%numran
			print *, '   list:',gb_cell_mat(ix,iy,iz)%iranlist(1:gb_cell_mat(ix,iy,iz)%numran)
			do i = 1, gb_cell_mat(ix,iy,iz)%numran
				iran = gb_cell_mat(ix,iy,iz)%iranlist(i)
				print *, '      iran  ',iran,iran-gb_numdata, ':',real(gb_xyz_list(1:3,iran))
			enddo
		enddo
	end subroutine test_grid
		
end module ap_cell	
	

