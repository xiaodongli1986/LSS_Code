
! What treatment we have adopted now:

	!-----------------------------------------------------------
	! We add special rejection of ra/dec range for SGC in has_boundary_effect()
	!  But in fact it does not save too much time

	!-----------------------------------------------------------
	! We find the grid size of SGCs are too large making the calculation too slow!
	!  Anyway since it's already fixed let's do not change it this time
	!  Next time we shall use 83/78 for CMASS-S/LOWZ-S

	!-----------------------------------------------------------
	! We start to test HR3 mbp mocks; we set gb_mockdir as 
	! 	'~/SparseFilaments/data/input/HR3/mbp_DR12_mock/'...

	!-----------------------------------------------------------
	! We add jackknife, 30 sub-samples
	! 

	!-----------------------------------------------------------
	! We reject a square near the hole (230,40); (check_reject_radecran)
	!   apply minimal/maximal r to mock/data; apply minimal/maxiaml redshift to data; (check_reject)
	!

	!
	! Use 3d random only for r < gb_3dranmaxr; read in 3d with r < gb_3dranmaxr+rext; we set rext=100.0_dl

	!-----------------------------------------------------------
	! For radially varing number density mock data we use wcp (rather than halo mass) as mass:
	!			gb_datalist(i)%mass=wcp!min(1.0,real(wcp));!wcp; !mass  
	! For constant number density mock (gb_i_datatype.eq.0) we use number density
	! For data we use: 	gb_datalist(i)%mass= (WCP+WNOZ-1.0)/Completeness*WSTAR*WSEE
			
	!-----------------------------------------------------------
	! outputinfo: for gb_cell_mat only output those pass ranwei_tol:
	! 		...
	!			if(gb_cell_mat(ix,iy,iz)%weirat .ge. gb_ranwei_tol) then
	
! July 1: remove numerous files, arrays, mpi-communications; remove 3d random; remove gfsampleing!!!
! June 30: zhengli chisq subroutine;... write code to calc drho chisq/ test; ...
! June 11: delete nuisance files like getcor/structure_count/grid;

program LSS_main

use mpi
use LSS_chisq

	implicit none

	! collection of settings
	type(chisq_settings) :: cs
	
	! loop indecies; numbers
	integer :: i, j, k, i_mask, i_mockname, i_mock, i_patch, i_bddist, i_vlos, i_wcp, i_sm, i_gridsize, icalc, &
		count_mockname, count_imockpatch

	! names, strings
	character(len=char_len) :: specnamestr, outputdirstr, &
		suffixstr, tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8
	
	! chi^2s
	real(dl), allocatable :: rhochi(:), rhodfchi(:), multdfchi(:,:)
	
	! cosmological parameters
	real(dl), allocatable :: om_w_list(:,:)
	real(dl) :: ommin, ommax, wmin, wmax, dom, dw, om, w
	integer :: num_omw, num_om, num_w
	integer, allocatable :: iomlist(:), iwlist(:)

	! mpi variables
	integer :: ierr, nproc, myid!, i1,i2,  mpistatus(MPI_STATUS_SIZE)
	
	! others
	real(dl) :: time1, time2, tmpx1,tmpx2,tmpx3, tmpx4,tmpx5,tmpx6, ra,dec!,r
	
	! for beta-skeleton test
	real(dl) :: tmpxs(6)
	logical :: Connected, erflaglog

	! Initialize MPI
	call mpi_init(ierr)
	call mpi_comm_size(mpi_comm_world,nproc,ierr)
	call mpi_comm_rank(mpi_comm_world,myid,ierr)

	!open(unit=1000,file='gb_outputnames.txt',access='append')
	! usual wcp
	!gb_wcpw = gb_wcpw_usual
	!gb_wcpw = gb_wcpw_le1
	!gb_wcpw = gb_wcpw_eq1
	icalc = -1
     do i_gridsize = 0,2!!!!TBC
     do gb_i_catalogue = 8,8
     do gb_i_datatype = 3,3!1,3
     	!-1: x,y,z,weight data
     	! 0: constant nbar, no RSD mock; 
     	! 1: radially varying nbar, no RSD mock; 
	! 2: radially varying nbar, with RSD mock; 
	! 3: real observational data

     do i_sm = 20,20 ! 2, 18,2
     do i_wcp = 1,1
     	if(i_wcp.eq.1) then
	     	gb_wcpw = gb_wcpw_usual
	elseif(i_wcp.eq.2) then
		gb_wcpw = gb_wcpw_eq1
	elseif(i_wcp.eq.3) then
		gb_wcpw = gb_wcpw_le1
	endif

     do i_mask = 2,2!,2
     	! 0: no mask; 1: Anderson mask applied; 2: Anderson+Veto masks applied
        count_mockname = 0
     do i_mockname = 2, 2
     	! 1: HR3 PSB; 2: HR3 PSB, mbp velocity; 3: HR4 J08Gal Snapshot; 4: HR4 J08
        count_mockname = count_mockname + 1
        count_imockpatch = 0
     do i_mock = 0,0
     do i_patch = 1,1
     	count_imockpatch = count_imockpatch + 1
     do i_vlos = 0, 0
!     do i_bddist = 0,0

	print *
 	print *, '###########################################################################'  
 	write(*,'(5(A,i3,";"))'), ' gb_i_datatype = ',gb_i_datatype, ' i_mask=',i_mask, &
 		'i_mock=',i_mock, ' i_patch=',i_patch, 'i_vlos=',i_vlos
 	print *, '###########################################################################'
 	print * 
 	
	if(gb_i_datatype.eq.0) then ! do not consider mask for type0
		if(i_mask.ne.0) cycle 
	elseif(gb_i_datatype.eq.3) then ! only consider one mock&patch for real data
		if(count_imockpatch.ne.1) cycle
		if(count_mockname.ne.1) cycle
	endif
	
	if(i_vlos.ne.0) then ! only consider vlos for RSD mock
		 if(gb_i_datatype.ne.2) cycle
	endif
	
	if(i_wcp.ne.1) then ! only consider wcp for RSD mock
		if(gb_i_datatype.ne.2) cycle 
	endif
	
	if(i_mockname.eq.1) then
		gb_mockname = gb_mockname_HR3PSB
	elseif(i_mockname.eq.2) then
!		gb_mockname = gb_mockname_HR4PSB
		gb_mockname = gb_mockname_HR3mbp
	elseif(i_mockname.eq.3) then
		gb_mockname = gb_mockname_HR4J08Gal
	elseif(i_mockname.eq.4) then
		gb_mockname = gb_mockname_HR4J08
	else
		stop
	endif

!##########################################################
!######  Sec A. Settings ##################################
!##########################################################

  !-------------------------------------------------------
  ! A0. Major directory of the output
  !-------------------------------------------------------
	outputdirstr = 'Aug18_FOGMocks/TestTime/'
	suffixstr = '' !'_info'

  !-------------------------------------------------------
  ! A1. Method: SPH length cut & # of smooth
  !-------------------------------------------------------
	
	! #-NNB used to compute rho/drho
	cs%smnum = i_sm
	gb_use_fixmd = .false.; gb_fixmd = 15.0*2.0 

  !-------------------------------------------------------
  ! A2. Number of grids in the x direction
  !	In other directions number will automatically determined by using the 
  !      same size of cell as the one used in x direction
  !-------------------------------------------------------

	! size of grid in x-direction; for boss data 120 shall be enough
	cs%num_in_x =  int( int(gb_suggested_numinx(gb_i_catalogue)+0.5) * (1.0 + 0.1 * i_gridsize))

	! slightly changing num_in_x to get multiple sets of results (#=2*gb_num_changenuminx+1); -> suppressing statsitical fluct
	gb_num_changenuminx = 0; gb_amp_changenuminx = 0.1_dl	

  !-------------------------------------------------------
  ! A3. cosmological parameters
  !-------------------------------------------------------

!	num_om = 1; num_w = 1; ommin = 0.26; ommax = 1.00; wmin = -1.0; wmax = -2.5;
	num_om = 1; num_w = 1; ommin = 0.26; ommax = 1.00; wmin = -1.0; wmax = -2.5;
	num_omw = num_om*num_w; allocate(om_w_list(2,num_omw),iomlist(num_omw), iwlist(num_omw))

!	num_omw = 6; allocate(om_w_list(2,num_omw),iomlist(num_omw), iwlist(num_omw))
!	om_w_list(1,1:6) = (/0.26_dl, 0.26_dl, 0.26_dl, 0.00_dl, 1.00_dl/)
!	om_w_list(2,1:6) = (/-1.0_dl, -0.5_dl, -3.0_dl, -1.0_dl, -1.0_dl/)
	
        do i = 1, num_om
		do j = 1, num_w
                        if(i.eq.1) then
			        om_w_list(1,(i-1)*num_w+j) = ommin 
                        else
	        		om_w_list(1,(i-1)*num_w+j) = ommin + dble(ommax-ommin)/dble(num_om-1)*dble(i-1)
                        endif
                        if(j.eq.1) then
        			om_w_list(2,(i-1)*num_w+j) = wmin
                        else
        			om_w_list(2,(i-1)*num_w+j) = wmin + dble(wmax-wmin)/dble(num_w-1)*dble(j-1)
                        endif
        	        iomlist((i-1)*num_w+j) = i-1; iwlist((i-1)*num_w+j) = j-1
		enddo
	enddo
	
  !-------------------------------------------------------
  ! A4. Data, masks, randoms
  !	
  !	 Here just set how many random files
  !-------------------------------------------------------

      !!!!----------------------------------------------
      !  how many random files?
      !!!!----------------------------------------------

	gb_numradecranfile = 2
	gb_numranfile = 10

      !!!!----------------------------------------------
      !  file names: data (real data, mock)
      !!!!----------------------------------------------
	! initialize the information of the catalogue...
	if(gb_i_catalogue.le.0.or.gb_i_catalogue.gt.gb_numcatalogue) then
		print *, 'ERROR (readin_dataran)!! wrong index of catalogue: ', gb_i_catalogue, gb_numcatalogue; stop
	endif
	gb_catalogue_name = gb_catalogue_names(gb_i_catalogue); gb_Seff = gb_Seffs(gb_i_catalogue)
	gb_dataredcut_min = gb_dataredcut_mins(gb_i_catalogue); gb_minimalr_cut = gb_minimalr_cuts(gb_i_catalogue)
	gb_dataredcut_max = gb_dataredcut_maxs(gb_i_catalogue); gb_maximalr_cut = gb_maximalr_cuts(gb_i_catalogue)
		
	! directory
	if(gb_mockname .eq. gb_mockname_HR3PSB) then
		gb_mockdir = '~/SparseFilaments/data/input/HR3/DR12_mock/'//trim(adjustl(gb_catalogue_name))//'/'  
		write(tmpstr1,*) i_mock
		if(i_mock .le. 9) then
		        gb_dataname = 'Xiao-dong.0000'//trim(adjustl(tmpstr1))//'.dat.compact.patch' ! Settings !!!
		else                
		        gb_dataname = 'Xiao-dong.000'//trim(adjustl(tmpstr1))//'.dat.compact.patch' ! Settings !!!
		endif
	elseif(gb_mockname .eq. gb_mockname_HR4PSB) then
		gb_mockdir = '~/SparseFilaments/data/input/HR4/PSB_DR12_mock/'//trim(adjustl(gb_catalogue_name))//'/' 
		write(tmpstr1,*) i_mock
		if(i_mock .le. 9) then
		        gb_dataname = 'n0000'//trim(adjustl(tmpstr1))//'halomass.dat.compact.B.C.D.patch' ! Settings !!!
		else                
		        gb_dataname = 'n000'//trim(adjustl(tmpstr1))//'halomass.dat.compact.B.C.D.patch' ! Settings !!!
		endif
	elseif(gb_mockname .eq. gb_mockname_HR3mbp) then
		gb_mockdir = '~/SparseFilaments/data/input/HR3/mbp_DR12_mock/'//trim(adjustl(gb_catalogue_name))//'/'
		write(tmpstr1,*) i_mock
		if(i_mock .le. 9) then
		        gb_dataname = 'Xiao-dong.pv.0000'//trim(adjustl(tmpstr1))//'.dat.compact.patch' ! Settings !!!
		else                
		        gb_dataname = 'Xiao-dong.pv.000'//trim(adjustl(tmpstr1))//'.dat.compact.patch' ! Settings !!!
		endif

	elseif(gb_mockname .eq. gb_mockname_HR4J08Gal) then
		gb_mockdir = '~/SparseFilaments/data/input/HR4/J08Gal/Red0Snapshot/'//trim(adjustl(gb_catalogue_name))//'/' 
		write(tmpstr1,*) i_mock+38
		if((i_mock+38) .le. 9) then
		        gb_dataname = 'RedshiftSlice0'//trim(adjustl(tmpstr1))//'.dat.compact.B.patch' ! Settings !!!
!		        gb_dataname = 'RedshiftSlice0'//trim(adjustl(tmpstr1))//'.dat.CompactObserved.patch'
		else                
		        gb_dataname = 'RedshiftSlice'//trim(adjustl(tmpstr1))//'.dat.compact.B.patch'
!      		        gb_dataname = 'RedshiftSlice'//trim(adjustl(tmpstr1))//'.dat.CompactObserved.patch'
		endif
	elseif(gb_mockname .eq. gb_mockname_HR4J08) then
		gb_mockdir = '~/SparseFilaments/data/input/HR4/'//trim(adjustl(gb_catalogue_name))//'/' 
		write(tmpstr1,*) i_mock
		if(i_mock .le. 9) then
		        gb_dataname = 'Xiao-dong.0000'//trim(adjustl(tmpstr1))//'.J08.dat.compact.patch' ! Settings !!!
		else                
		        gb_dataname = 'Xiao-dong.000'//trim(adjustl(tmpstr1))//'.J08.dat.compact.patch' ! Settings !!!
		endif
	else
		stop
	endif
	
	gb_bossdatadir = '~/SparseFilaments/data/input/DR12/OurRandoms/'//trim(adjustl(gb_catalogue_name))//'/Om0.260w-1.000/'
	gb_randomdir = '~/SparseFilaments/data/input/DR12/OurRandoms/'//trim(adjustl(gb_catalogue_name))
	
	! file names
	if(gb_i_datatype.ne.3.and.gb_i_datatype.ne.-1) then
		! which patch
       	        write(tmpstr2,*) i_patch
       	        gb_dataname = trim(adjustl(gb_dataname))//trim(adjustl(tmpstr2))
       	        ! suffcies        	        
		if(gb_i_datatype.eq.1) then
			gb_dataname = trim(adjustl(gb_dataname))//'.noRSD-radial-selected'
		elseif(gb_i_datatype.eq.2) then
			if(i_vlos.eq.0) then
				gb_dataname = trim(adjustl(gb_dataname))//'.RSD-radial-selected'
			else
				write(tmpstr1,*) i_vlos
				gb_dataname = trim(adjustl(gb_dataname))//'.RSD-radial-selected.vlos'//trim(adjustl(tmpstr1))
			endif
		endif
		! location of the file
       		gb_datafile = trim(adjustl(gb_mockdir))//trim(adjustl(gb_dataname))
	elseif(gb_i_datatype.eq.3) then
		gb_dataname = 	'Data-xyz-red-wei.txt' ! name of the data file
		gb_datafile = trim(adjustl(gb_bossdatadir))//trim(adjustl(gb_dataname))
	elseif(gb_i_datatype.eq.-1) then
		om_dft = 0.274_dl
		gb_catalogue_name = 'DR12v1-CMASS-N/xyzweight'
		gb_dataname = 'xyz_rec_3rec_it0.txt'
		gb_datafile = '~/SparseFilaments/data/input/RecMocks/QPM-PosRec/DR12-CMASS-N/'//trim(adjustl(gb_dataname))
	else
		print *, 'Wrong gb_i_datatype!', gb_i_datatype; continue
	endif

      !!!!----------------------------------------------
      !  file names: random
      !!!!----------------------------------------------
	! file names for 2d ra/dec randoms
	do i = 1, gb_numradecranfile
		write(tmpstr1,*) i-1
		gb_radecranfilelist(i) = &
			trim(adjustl(gb_randomdir))//'/RanRADEC'//trim(adjustl(tmpstr1))//'.txt '
	enddo
	! file names for 3d randoms 
	k = 0
	do i = 0, gb_numranfile/5+2
		do j = 0, 5
			k = k+1; if(k>gb_numranfile) exit
			write(tmpstr1,*) i
			write(tmpstr2,*) j
			gb_ranfilelist(k) = trim(adjustl(gb_randomdir))//'/Om0.260w-1.000/RanPoint/RADEC'&
				//trim(adjustl(tmpstr1))//'_Shuffle'//trim(adjustl(tmpstr2))//'.txt'
		enddo
	enddo

      !!!!----------------------------------------------
      !  number density or weighted density?
      !!!!----------------------------------------------
        ! if constant number density mock, then use number density; use weight for other types of data
        if(gb_i_datatype.eq.0) then
        	gb_usenumdensity = .true.
        else
        	gb_usenumdensity = .false.
        endif

      !!!!----------------------------------------------
      !  masks: 
      !    nothing; Anderson; Anderson+Veto. 
      !    They will affect the randoms, HR3 mocks
      !!!!----------------------------------------------
	if(i_mask.eq.0) then
		gb_AdsnRej = .false.; gb_VetoRej = .false.
	elseif(i_mask.eq.1) then
		gb_AdsnRej = .true.; gb_VetoRej = .false.
	elseif(i_mask.eq.2) then
		gb_AdsnRej = .true.; gb_VetoRej = .true.
	else
		print *, 'Wrong i_mask!', i_mask; continue
	endif

      !!!!---------------------------------------------- 
      !  randoms: 
      !    will use ra/dec random, much faster, more stable than 3D random; 
      !    two files shall be enough
      !!!!----------------------------------------------

	! use 2D rather than 3D; we may consider 3D for lowz...
	gb_i_rantype = gb_radecran

	! tolerance of ratio of weights of in-sphere randoms who pass the masks
	gb_ranwei_tol = 0.8 

    !********************************************************
    ! These are for tests
    !********************************************************
	if(.false.) then
		call readin_dataran(.true.)
	!	call test_readindata()
		gb_chisq_initied = .true.
		call init_mult_lists(.true.)
		call do_cell_init(50.0_dl, .true.)		
!		call test_grid()
!		call NNBTest()
		stop
	endif

  !-------------------------------------------------------
  ! A5. boundary correction
  !-------------------------------------------------------

	gb_bddist = 0
	gb_bddist_rextra = 10
	gb_bufferdist = 0
	gb_bufferdist_rextra = 0

  !-------------------------------------------------------
  ! A6. screen print of information
  !-------------------------------------------------------
	if(myid.eq.0) then
		cs%print_info = .true. 
	else
		cs%print_info = .false.
	endif

  !-------------------------------------------------------
  ! A7. Dropping high rho/drho
  !-------------------------------------------------------
	cs%numdrop = 1 			!!CHECK
	allocate(cs%dropval(cs%numdrop),  cs%lowdropvalratio(cs%numdrop),  cs%highdropvalratio(cs%numdrop), &
		cs%dropdval(cs%numdrop), cs%lowdropdvalratio(cs%numdrop), cs%highdropdvalratio(cs%numdrop))
	! Our conventional setting of 11 dropping ratios : 0,10,20,30,...,80,90,95.
	do i = 1, cs%numdrop
		cs%dropval(i) = .true.;	cs%lowdropvalratio(i) = 0.0
		cs%dropdval(i) = .true.; cs%lowdropdvalratio(i) = 0.0
		cs%highdropvalratio(i) = min((i-1)*0.20_dl,0.95_dl) !!CHECK
		cs%highdropdvalratio(i) = min((i-1)*0.20_dl,0.95_dl)!!CHECK
	enddo

!##########################################################
!######  Sec B. Main body #################################
!##########################################################
  !-------------------------------------------------------
  ! B0. Create a name for the run
  !-------------------------------------------------------

	!!CHECK !!MUST BE CHECKED ON CLUSTER
     !## basic name: output directory + name of catalogue + name of data file
	gb_output_name     = trim(adjustl(outputdirstr))//trim(adjustl(gb_catalogue_name))//'/'//trim(adjustl(gb_dataname))
	tmpstr1 = 'mkdir -p '//trim(adjustl(gb_output_name)); call system(tmpstr1); !call system('sleep 2')

     !## grid size
        write(tmpstr1,*) cs%num_in_x
        gb_output_name      = trim(adjustl(gb_output_name))//'/'//trim(adjustl(tmpstr1))//'cube_' ! redshift cut 0.36

     !## tolerance threshold of random weight ratio
        write(tmpstr1,'(f5.2)') gb_ranwei_tol
        gb_output_name      = trim(adjustl(gb_output_name))//'weitol'//trim(adjustl(tmpstr1)) ! redshift cut 0.36
     	
     !## mask
        if(gb_AdsnRej.and.gb_VetoRej)  then
        	gb_output_name      = trim(adjustl(gb_output_name))//'_AdsnVetoRej'
        elseif(gb_AdsnRej) then
	        gb_output_name      = trim(adjustl(gb_output_name))//'_AdsnRej'
	elseif(gb_VetoRej) then
		gb_output_name      = trim(adjustl(gb_output_name))//'_VetoRej'
	endif
	
     !## wcp settings
	if(gb_wcpw.eq.gb_wcpw_usual) then
		continue
	elseif(gb_wcpw.eq.gb_wcpw_eq1) then
		gb_output_name      = trim(adjustl(gb_output_name))//'_wcpeq1'
	elseif(gb_wcpw.eq.gb_wcpw_le1) then
		gb_output_name      = trim(adjustl(gb_output_name))//'_wcple1'
	endif

     !## smoothing: #-NNB, smoothing sphere radius        
	if(gb_use_fixmd) then
		write(tmpstr1,'(f4.1)') gb_fixmd/2.0
	        gb_output_name      = trim(adjustl(gb_output_name))//'_smh'//trim(adjustl(tmpstr1)) ! redshift cut 0.36
	else
	        write(tmpstr1,*) cs%smnum
		gb_output_name     = trim(adjustl(gb_output_name))//'_smnum'//trim(adjustl(tmpstr1))
	endif

     !## multi grids
	if(gb_num_changenuminx.ne.0) then
		write(tmpstr1,*) gb_num_changenuminx
		write(tmpstr2,'(f6.2)') gb_amp_changenuminx
		gb_output_name = trim(adjustl(gb_output_name))//'_multrlts-num'//trim(adjustl(tmpstr1))//'-amp'//trim(adjustl(tmpstr2))
	endif
	
     !## normalizatin of nbar
	if(.not. gb_dodensitynorm) then
		gb_output_name = trim(adjustl(gb_output_name))//'_nonorm'
	else
		write(tmpstr1,*) gb_normn_nnorm
		write(tmpstr2,*) gb_normn_nbin
		gb_output_name = trim(adjustl(gb_output_name))//'_norm'//trim(adjustl(tmpstr1))//'nbin'//trim(adjustl(tmpstr2))
	endif

     !## range of r
     	write(tmpstr1,'(f15.3)') gb_minimalr_cut
     	write(tmpstr2,'(f15.3)') gb_maximalr_cut
     	gb_output_name = trim(adjustl(gb_output_name))//'_r'//trim(adjustl(tmpstr1))//'to'//trim(adjustl(tmpstr2))

     !## boundary correction	
	write(tmpstr1,*) int(gb_bddist+0.5)
	write(tmpstr2,*) int(gb_bddist_rextra+0.5)
	write(tmpstr3,*) int(gb_bufferdist+0.5)
	write(tmpstr4,*) int(gb_bufferdist_rextra+0.5)
	write(tmpstr5,*) gb_normn_nnorm	
	gb_output_name = trim(adjustl(gb_output_name))//'_bddist'//trim(adjustl(tmpstr1))//'-rex'//trim(adjustl(tmpstr2))// &
		'_buffer'//trim(adjustl(tmpstr3))//'-rex'//trim(adjustl(tmpstr4))

     !## random
     	if(gb_i_rantype.eq.gb_radec3dran) then
     		gb_output_name = trim(adjustl(gb_output_name))//'_3dradecran'
     	elseif(gb_i_rantype.eq.gb_radecran) then
     		gb_output_name = trim(adjustl(gb_output_name))//'_radecran'
     	elseif(gb_i_rantype.eq.gb_3dran) then
     		gb_output_name = trim(adjustl(gb_output_name))//'_3dran'
     	endif
     	
     !## extra suffix
     	if(trim(adjustl(suffixstr)).ne.'') then
		gb_output_name = trim(adjustl(gb_output_name))//trim(adjustl(suffixstr))
	endif

     !## creat directory
        call system('mkdir -p '//trim(adjustl(gb_output_name)))

  !-------------------------------------------------------
  ! B1. Selected print of information
  !-------------------------------------------------------
     
	if(myid .eq. 0) then
		write(*,*) '######################################'
		write(*,*) ' Settings of the programm:'
		print * 
	        write(*,'(A)') '   Data file::: '//trim(adjustl(gb_datafile))
                if(gb_i_datatype.eq.-1) then
                	write(*,'(A)')  '   data type: xyzw'
                elseif(gb_i_datatype.eq.0) then
                	write(*,'(A)')  '   data type: constant number density'
                elseif(gb_i_datatype.eq.1) then
                	write(*,'(A)')  '   data type: varying radial density, no RSD'
                elseif(gb_i_datatype.eq.2) then
                	write(*,'(A)')  '   data type: varying radial density, with RSD'
                elseif(gb_i_datatype.eq.3) then
                	write(*,'(A)')  '   data type: real data'
                endif
                print *
		write(*,'(A,3f11.4)')   '   Om_dft  / w_dft / h_dft = ', om_dft, w_dft, h_dft
		write(*,'(A,e15.7,2f6.3,2f10.3)')   &
			'   gb_Seff, gb_dataredcut_min, gb_dataredcut_max, gb_minimalr_cut, gb_maximalr_cut = ', &
                	gb_Seff, gb_dataredcut_min, gb_dataredcut_max, gb_minimalr_cut, gb_maximalr_cut
		write(*,'(A,2L2)')   '   gb_AdsnRej, gb_VetoRej = ', gb_AdsnRej, gb_VetoRej
		if(gb_wcpw.eq.gb_wcpw_usual) then
			write(*,'(A)')   '   wcp: usual wcp'
		elseif(gb_wcpw.eq.gb_wcpw_eq1) then
			write(*,'(A)')   '   wcp: equal 1'
		elseif(gb_wcpw.eq.gb_wcpw_le1) then
			write(*,'(A)')   '   wcp: lessequal 1'
		endif

		print *
                write(*,'(A)') '   Outputname::: '//trim(adjustl(gb_output_name))
                print *
		write(*,'(A,i4)') '   Cube (x-direc)    = ', cs%num_in_x
                write(*,'(A,i3,f7.3)')  '   gb_num_changenuminx, gb_amp_changenuminx = ', gb_num_changenuminx, gb_amp_changenuminx
		write(*,'(A,l2)') '   use_num_density   = ', gb_usenumdensity 
                write(*,'(A,l2,2i3,i7)')   '   gb_dodensitynorm, gb_normn_nbin, gb_normn_nnorm, gb_radecmat_size = ', &
                	gb_dodensitynorm, gb_normn_nbin, gb_normn_nnorm, gb_radecmat_size
	!++++++++++++++++++++++++++++++++++                	
		print *
		write(*,'(3x,A,L2,L2)') 'Use_fixmd / gb_do_seg_cut = ', gb_use_fixmd, gb_do_seg_cut
		if(gb_use_fixmd) then
			write(*,'(3x,A,f10.3)') '  fixed smoothing length gb_fixmd = ', real(gb_fixmd)
		else
			write(*,'(3x,A,i4)') '  using #-NNB = ', cs%smnum
			write(*,'(3x,A,f6.3,A,L3)') '  NNB Over-Searching Rat = ', gb_NNBosrat, &
				'; do NNB-ReSearch = ', gb_do_NNBReSearch
		endif
		if(gb_do_seg_cut)  then
			write(*,'(3x,A,f10.3)') 'Applying minimal-r cut: ', real(gb_seg_cut_dist)
			write(*,*) 'ERROR (grid_rho_drho_list): gb_do_seg_cut not supported now!'; stop
		endif
		write(*,'(3x,A)') 'Settings of Boundary Correction: '
		write(*,'(6x,A,4f10.3)') 'bd, bd_rextra, buffer, buffer_rextra = ', &
			gb_bddist, gb_bddist_rextra, gb_bufferdist, gb_bufferdist_rextra
			
		print *
		
		write(*,'(A,i3,A)') '   Num of drop       = ', cs%numdrop, '; droppings:'
		do i = 1, cs%numdrop
			write(*,'(10x,A,L2,f6.3,f6.3,A,L2,f6.3,f6.3)') 'density: ', cs%dropval(i), real(cs%lowdropvalratio(i)),  &
				real(cs%highdropvalratio(i)), ';   density Gradient: ', &
				cs%dropdval(i), real(cs%lowdropdvalratio(i)), real(cs%highdropdvalratio(i))
		enddo
		print *
		write(*,'(A,i5)') '   Num of omegam/w       = ', num_omw	
                write(*,'(10x,A,$)') 'Lists of omegam/w: '
                do i = 1, num_omw
                        write(*,'(f6.3,"/",f6.3,$)') om_w_list(1:2,i)
                enddo
                print *
                print *

	endif

  !-------------------------------------------------------
  ! B2. Do initializations for the calculation
  !-------------------------------------------------------

	allocate(rhochi(cs%numdrop),rhodfchi(cs%numdrop), multdfchi(nbinchisq,cs%numdrop) )
	
!	write(1000,'(A,A,A,i3,A)') '[ "', trim(adjustl(gb_output_name(50:200))), ', "$', cs%num_in_x, '^3$ grid"],'
!	goto 1000
	
!	call mpi_barrier(mpi_comm_world,ierr)
	
	if(myid .eq. 0) then
		write(*,*) '######################################'
        	print *, 'Calculating chisqs...'	
        endif
	
	call cpu_time(time1)

  !-------------------------------------------------------
  ! B3. Calculation
  !-------------------------------------------------------
	do i =  1, num_omw
		icalc =icalc + 1
		if(mod(icalc,nproc).ne.myid) then
			cycle
		else
			write(*,'(A,i6,i4,2f10.4)'), '   @@@ Do calculation: icalc, myid, omw = ', icalc, myid, real(om_w_list(1:2,i))
		endif
				
		om = om_w_list(1,i)
		w  = om_w_list(2,i)
		
		gb_iom = iomlist(i); gb_iw = iwlist(i)

		write(tmpstr1,'(f8.4)') om
		write(tmpstr2,'(f8.4)') w
		gb_omwstr = 'om'//trim(adjustl(tmpstr1))//'_w'//trim(adjustl(tmpstr2))
		gb_suboutput_name = trim(adjustl(gb_output_name))//'/'//trim(adjustl(gb_omwstr))
		
		call cosmo_funs_init(cs%print_info)
		call readin_dataran(cs%print_info)
		gb_chisq_initied = .true.

		call gf_mldprho_chi2s(om, w, h_dft, cs, rhochi, rhodfchi, multdfchi, cs%numdrop, calc_comvr = .true.)
		write(*,'(A,i3,A,2(f8.3,1x))') '   Step ',i,'. Chisqs without RSD:', om, w
		write(*,'(36x,A)') trim(adjustl(WriteFmtEFloat(cs%numdrop, rhochi(1:cs%numdrop))))
		write(*,'(A,12x,A)') '    2 bins (old method):', &
		  trim(adjustl(WriteFmtEFloat(cs%numdrop,rhodfchi(1:cs%numdrop))))
		do j = 1, nbinchisq	
			write(*,'(4x,i4,A,24x,A)') nbinchisqlist(j),' bins:', &
			 trim(adjustl(WriteFmtEFloat(cs%numdrop,multdfchi(j,1:cs%numdrop))))
		enddo
	
		call cpu_time(time2)
		write(*,'(A,f10.4)') '  Time used in this step: ', time2-time1
		time1 = time2
		print *
		cs%print_info = .false.
	enddo

1000	continue
	deallocate(om_w_list,iomlist, iwlist, rhochi,rhodfchi,multdfchi)
	deallocate(cs%dropval,  cs%lowdropvalratio, cs%highdropvalratio, &
			cs%dropdval, cs%lowdropdvalratio, cs%highdropdvalratio)
1001	continue			
     enddo
     enddo
     enddo
     enddo
     enddo
     enddo
     enddo
     enddo
     enddo
     enddo
     	close(1000)
	call mpi_barrier(mpi_comm_world,ierr)
	call mpi_finalize(ierr)
end program LSS_main
