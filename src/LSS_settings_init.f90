! Some information... (we remove RSD shift code, but keep some note here for quick reference...)
!   Observed redshift of an object:
!      	  (1 + z_cosmology) * (1 + z_doppler)
! 	=  1 + z_cosmology + (1 + z_cosmology) * v_local / c 	!!! (v_local is the local pecular velocity)   
! 	=  1 + z_real + v_comov / c	!!! (v_comov = v_local / a = v_local * (1+z_cosmology))

!####################################
!This module does smooth
!####################################
module LSS_settings_init
use LSS_cosmo_funs

	implicit none

!!! Fixed Settings

	! Collection of information of catalogues
	integer, parameter :: gb_numcatalogue = 8

	! catalogue names
	character(len=char_len), parameter :: 	gb_catalogue_name1 = 'DR12v1-CMASS-N', &
						gb_catalogue_name2 = 'DR12v1-LOWZ-N', &
						gb_catalogue_name3 = 'DR12v1-CMASS-S', &
						gb_catalogue_name4 = 'DR12v1-LOWZ-S', &
						gb_catalogue_name5 = 'DR12v4-CMASS-N', &
						gb_catalogue_name6 = 'DR12v4-LOWZ-N', &
						gb_catalogue_name7 = 'DR12v4-CMASS-S', &
						gb_catalogue_name8 = 'DR12v4-LOWZ-S'
						
	! mock name: now we have HR3PSB, HR4PSB, HR4J08GalSnapshot...						
	character(len=char_len), parameter ::   gb_mockname_HR3PSB = 'HR3PSB', &
						gb_mockname_HR4PSB = 'HR4PSB', &
						gb_mockname_HR3mbp = 'HR3mbpPSB', &
						gb_mockname_HR4J08Gal = 'HR4J08GalSnapshot', &
						gb_mockname_HR4J08 = 'HR4J08'
						
	character(len=char_len), parameter :: gb_catalogue_names(gb_numcatalogue) = &
		(/ gb_catalogue_name1, gb_catalogue_name2, gb_catalogue_name3, gb_catalogue_name4, &
		   gb_catalogue_name5, gb_catalogue_name6, gb_catalogue_name7, gb_catalogue_name8 /)
	
	! information: Seff, minimal/maximal redshift and radius
	real(dl), parameter :: gb_Seffs(gb_numcatalogue) = (/ 6.918435E3_dl,  5.880969E+03_dl, 2.628270E+03_dl, 2.605330E+03_dl, &
							      6.851419E3_dl,  2.524669E+03_dl, 5.836207E+03_dl, 2.501263E+03_dl /)
	real(dl), parameter :: gb_dataredcut_mins(gb_numcatalogue) = (/  0.43_dl, 0.15_dl, 0.43_dl, 0.15_dl, & 
									 0.43_dl, 0.15_dl, 0.43_dl, 0.15_dl /)
	real(dl), parameter :: gb_dataredcut_maxs(gb_numcatalogue) = (/  0.70_dl, 0.43_dl, 0.70_dl, 0.43_dl, &
									 0.70_dl, 0.43_dl, 0.70_dl, 0.43_dl  /)
	real(dl), parameter :: gb_minimalr_cuts(gb_numcatalogue) = (/ 1172.7394_dl,  436.058_dl, 1172.7394_dl,  436.058_dl, &
	 							      1172.7394_dl,  436.058_dl, 1172.7394_dl,  436.058_dl/)
	real(dl), parameter :: gb_maximalr_cuts(gb_numcatalogue) = (/ 1787.3546_dl, 1172.8_dl, 1787.3546_dl, 1172.8_dl, &
								      1787.3546_dl, 1172.8_dl, 1787.3546_dl, 1172.8_dl /)
	
	real(dl), parameter :: gb_suggested_numinx(gb_numcatalogue) = (/ 120.0_dl, 107.0_dl, 83.00_dl, 78.0_dl, &
									 120.0_dl, 107.0_dl, 83.00_dl, 78.0_dl/)

!	real(dl), parameter :: gb_suggested_numinx(gb_numcatalogue) = (/ 120.0_dl, 107.0_dl, 100.00_dl, 101.0_dl, &
!									 120.0_dl, 107.0_dl, 100.00_dl, 101.0_dl/)
									 
				! They are good for smnum=20. If using smaller smnum, increasing the size:
					!GridSize2= GridSize1 * [(#-NNB1)/(#-NNB2)]^(1/3)
				! grid sizes and corresponded to the following cell widths and grid sizes: 
					! CMASS-N 13.94242 (120*230*127) 
					! LOWZ-N 13.93271 (83*150* 84) 
					! CMASS-S 13.83703 (68*179*101) 
					! LOWZ-S 13.99295 (78*115* 65) 

	! Type of data; -1: xyzw; 0: constant number density; 1: varying radial density, no RSD; 2: varying nbar, with RSD; 3: real data
	integer, parameter :: gb_dt_xyz = -2
	integer, parameter :: gb_dt_xyzw = -1
	integer, parameter :: gb_dt_ConstnbarM = 0
	integer, parameter :: gb_dt_nbarz_NoRSDM = 1
	integer, parameter :: gb_dt_nbarz_RSDM = 2
	integer, parameter :: gb_dt_realdata = 3
	
	! Types of randoms: 2d? 3d? Both? Neither? 
	integer, parameter :: gb_noran = 1, gb_radecran = 2, gb_3dran = 3, gb_radec3dran = 4

	real(dl), parameter :: gb_3dranmaxr = 500.0_dl 	! if chosed to use both raded/3drandom (gb_radec3dran), 
							!  then use 2d random if distance > gb_3dranmaxr
	
	!Global controls of outputing information
	logical, parameter :: gb_allow_output = .false.
	logical, parameter :: gb_outputinfo_nbar = .false..and.gb_allow_output, &
			gb_outputinfo_gfcell = .false..and.gb_allow_output, &
			gb_outputinfo_radeccell = .false..and.gb_allow_output, &
			gb_outputinfo_mu = .true., &
			gb_outputinfo_absdrho = .false.

!!! Adjustable & To be fixed Settings	
! All variables in the "Adjustable & To be fixed Settings" section must be assigned a value in the main**.f90!!!
	! Fiducial cosmology
	real(dl) :: om_dft = 0.26_dl,  w_dft = -1.0_dl, h_dft = 0.7_dl
	
	! type of data, type of ran
	integer :: gb_i_datatype = -10000
	integer :: gb_i_rantype = -1 !!! MUST SET it in main programe!

	! index of catalogue: must be a number from 1 to gb_numcatalogue
	integer :: gb_i_catalogue = -1

	! wcp weight: 1 means equal to 1; 2 means wcp; 0 means max(wcp,1) 
	integer, parameter :: gb_wcpw_le1 = 0, gb_wcpw_eq1 = 1, gb_wcpw_usual = 2
	integer :: gb_wcpw = gb_wcpw_usual

	!=========================================
	! fmt
		! -1: x,y,z, weight
		! 0 constant number density mock:  x,y,z,vx,vy,vz,mass, z_cosmology, z_obs, r(z_obs)/r(z_cosmology), ID (11 cols)
		! 1/2 varying radial density noRSD/withRSD mock:  x,y,z,wfkp, mass, ID, Andersonmask, vetomask, wcp (9 cols)
		! 3 data: x,y,z,redshift,ID(start from 1), w_tot =wsee*wstar*(wcp+wnoz-1) (without fkp), wfkp, completeness, wfkp2 (re-calculated in HR3 cosmology), wcp, wnoz, wstar, wsee (13 cols)

	! Tolerance of boundary enffect.
	real(dl) :: gb_ranwei_tol = 0.8 ! If (weight-of-acpt)/(weight-of-all) < gb_ranweitol, then reject this position.

	! Name of the project (used to creat outputs), files (data and random); information of catalgues, 
	character(len=char_len) :: gb_catalogue_name, gb_mockname, gb_mockdir, gb_bossdatadir, gb_randomdir, &
		gb_dataname, gb_datafile, gb_ranfilelist(1000), gb_radecranfilelist(1000) ! List of random files
	integer :: gb_numranfile, gb_numradecranfile
	real(dl) :: gb_Seff, gb_dataredcut_min, gb_dataredcut_max, gb_minimalr_cut=-1.0e30, gb_maximalr_cut=1.0e30
				
	! Format of data. Will determine what part of the data are loaded.
	logical :: gb_AdsnRej, gb_VetoRej
		
!!! Important Variables
	! Global strs used everywhere when outputing files
	character(len=char_len) :: gb_output_name, gb_suboutput_name, gb_omwstr

	! Number of data, random, all
	integer, public :: gb_numdata, gb_numran, gb_numallpar, gb_numradecran ! particle means data + 3dran; radecran are not particles!

	! Structure and array to restore data and random
	type :: datapoint
		real(dl) :: x,y,z, r,red, rat=1.0, mass
	end type
	type(datapoint), allocatable :: gb_datalist(:)
	
	type :: ranpoint
		real(dl) :: x,y,z, r,red, rat=1.0
		logical :: acpt
	end type
	type(ranpoint), allocatable :: gb_ranlist(:)
	
	type :: radecranpoint
		real(dl) :: ra,dec
		logical :: acpt
	end type
	type(radecranpoint), allocatable :: gb_radecranlist(:)
	
	! Range of xyzr, data, random, and grid.
	real(dl) :: gbxmindata, gbxmaxdata, gbymindata, gbymaxdata, gbzmindata, gbzmaxdata, gbrmindata, gbrmaxdata, &
		    gbxminran,  gbxmaxran,  gbyminran,  gbymaxran,  gbzminran,  gbzmaxran,  gbrminran,  gbrmaxran, &
		    gbgridxmin, gbgridxmax, gbgridymin, gbgridymax, gbgridzmin, gbgridzmax, gbgridrmin, gbgridrmax, &
		    gbramindata, gbramaxdata, gbdecmindata, gbdecmaxdata, &
		    gbraminran,  gbramaxran,  gbdecminran,  gbdecmaxran, &
    		    gbraminradecran,  gbramaxradecran,  gbdecminradecran,  gbdecmaxradecran
    	! Notice that gbgridra, gbgriddec are for the ra/dec grid, not for the 3d grid!!!!!!
		    
!	real(dl) :: gb_ramin_data, gb_decmin_data, gb_ramax_data, gb_decmax_data, &
!		    gb_ramin_ran,  gb_decmin_ran,  gb_ramax_ran,  gb_decmax_ran

contains
	
  !------------------------------------------
  ! volume between r1 and r2
  !------------------------------------------
	real(dl) function vol_fun(r1, r2)
		real(dl) :: r1, r2
		vol_fun = (4.0*const_pi/3.0) * (r2**3.0-r1**3.0) * gb_Seff / (41252.96124941928d0)
	end function vol_fun

  !------------------------------------------
  ! check reject: whether reject the point (data or 3d random)
  !------------------------------------------
	logical function check_reject(x,y,z)
		! Dummy
		real(dl), intent(in) :: x,y,z
		! Local
		real(dl) :: ra,dec,r,rsq

		check_reject = .false.
		r = sqrt(x**2.0+y**2.0+z**2.0)
		if( r.le. gb_minimalr_cut .or. r.ge.gb_maximalr_cut) then
			check_reject = .true.; return
		else			
			check_reject = .false.
		endif
	end function check_reject
  !------------------------------------------
  ! check reject: whether reject (2d ra/dec random)
  !------------------------------------------
	logical function check_reject_radecran(ra,dec)
		! Dummy
		real(dl), intent(in) :: ra,dec

		check_reject_radecran =  .false.
		if(ra>228.and.ra<232.and.dec>38.5.and.dec<41.5) then
			check_reject_radecran = .true.; return
		endif
	end function check_reject_radecran

	
  !------------------------------------------
  ! read in data: fmt x,y,z,weight
  !------------------------------------------
	subroutine readin_xyzwdata(printinfo,hasweight)
		! Dummy
		logical, intent(in) :: printinfo, hasweight
		! Local
		integer :: i,j, totnumdata
		real(dl), allocatable :: tmp(:,:)
		
		if(hasweight) then
			call read_in (gb_datafile, 4, totnumdata, tmp)
		else
			call read_in (gb_datafile, 3, totnumdata, tmp)
		endif
		
		gb_numdata = 0.0
		do i = 1, totnumdata
			if(check_reject(tmp(i,1),tmp(i,2),tmp(i,3))) cycle
			if(hasweight) then
				if(tmp(i,4).eq.0.0) cycle
			endif
			gb_numdata = gb_numdata + 1
		enddo
		
		if(allocated(gb_datalist)) deallocate(gb_datalist)

		if(printinfo) then
			write(*,'(A,A)') '   (readin_xyzwdata) read in from file:', trim(adjustl(gb_datafile))
			write(*,'(25x,A,2i8)')     '#-data (infile/readin) = ', totnumdata, gb_numdata
		endif

		j = 1
		if(allocated(gb_datalist)) deallocate(gb_datalist)
		allocate(gb_datalist(gb_numdata))
		do i = 1, totnumdata
			if(check_reject(tmp(i,1),tmp(i,2),tmp(i,3))) cycle
			if(hasweight) then
				if(tmp(i,4).eq.0.0) cycle
			endif
			gb_datalist(j)%x=tmp(i,1); gb_datalist(j)%y=tmp(i,2); gb_datalist(j)%z=tmp(i,3)
			gb_datalist(j)%r=rms(tmp(i,1:3),3) !these x,y,z,r are values in the fiducial cosmology!
			gb_datalist(j)%red=de_zfromintpl(gb_datalist(j)%r)
			if(hasweight) then
				gb_datalist(j)%mass=tmp(i,4)
			else
				gb_datalist(j)%mass=1.0
			endif
			j = j+1
		enddo
		if(j-1.ne.gb_numdata) then !consistency check
			print *, 'ERROR (readin_xyzwdata)! Mismatch of #: ', gb_numdata, j-1; stop
		endif
	end subroutine readin_xyzwdata
	
  !------------------------------------------
  ! read in the HR3 mock data
  !------------------------------------------
	subroutine readin_mockdata(printinfo)
		! Dummy
		logical, intent(in) :: printinfo
		! Local
		integer :: i, tot_numdata, numdata, ID, Andersonwei,vetowei,wei
		real(dl) :: tmp(10),wfkp,mass,wcp

		tot_numdata = 0; numdata = 0

		if(allocated(gb_datalist)) deallocate(gb_datalist)
		open(unit=4401,file=gb_datafile)
		do while(.true.)
       			read(4401,*,end=100,err=101) tmp(1:3),wfkp,mass,ID,Andersonwei,vetowei,wcp !!! fmt:x,y,z,mass,ID,two criteria
       			tot_numdata = tot_numdata+1
			wei=1; if(gb_AdsnRej) wei=wei*Andersonwei; if(gb_VetoRej) wei=wei*vetowei !wei=0 if any checked criteria fail
			if(gb_wcpw.eq.gb_wcpw_usual) then
				wcp = wcp
			elseif(gb_wcpw.eq.gb_wcpw_eq1) then
				wcp = 1.0_dl
			elseif(gb_wcpw.eq.gb_wcpw_le1) then
				wcp = max(wcp,1.0_dl)
			endif
			wei = wei * int(wcp+0.5)
			if(wei>0.5 .and. (.not.check_reject(tmp(1),tmp(2),tmp(3)))) numdata = numdata+1 !acpt ony if wei=1 & no other rej
			cycle
100			exit
101			print *, 'ERROR (read_in_HR3mockdata)! Error happes when read in ', trim(adjustl(gb_datafile)); stop
		enddo
		close(4401)
		
		gb_numdata = numdata
		if(printinfo) then
			write(*,'(A,A,A,2i8)') '   (read_in_HR3mockdata) '
			write(*,'(8x,A)')         trim(adjustl(gb_datafile))
			write(*,'(25x,A,2i8)')     '#-data (infile/readin) = ', tot_numdata, gb_numdata
		endif

		i = 1
		allocate(gb_datalist(gb_numdata))
		open(unit=4401,file=gb_datafile)
		do while(.true.)		
       			read(4401,*,end=103,err=104) tmp(1:3),wfkp,mass,ID,Andersonwei,vetowei,wcp !!! fmt:x,y,z,mass,ID,two criteria
       			tot_numdata=tot_numdata+1
			wei=1; if(gb_AdsnRej)wei=wei*Andersonwei; if(gb_VetoRej)wei=wei*vetowei !wei=0 if any checked criteria fail
			if(gb_wcpw.eq.gb_wcpw_usual) then
				wcp = wcp
			elseif(gb_wcpw.eq.gb_wcpw_eq1) then
				wcp = 1.0_dl
			elseif(gb_wcpw.eq.gb_wcpw_le1) then
				wcp = max(wcp,1.0_dl)
			endif
			wei = wei * int(wcp+0.5)
			if(wei>0.5 .and. (.not.check_reject(tmp(1),tmp(2),tmp(3))))  then !acpt ony if wei=1 & no other rej
				gb_datalist(i)%x=tmp(1); gb_datalist(i)%y=tmp(2); gb_datalist(i)%z=tmp(3)
				gb_datalist(i)%r=rms(tmp(1:3),3) !these x,y,z,r are values in the fiducial cosmology!
				wei = wei * int(wcp+0.5)
				gb_datalist(i)%mass=wcp!min(1.0,real(wcp));!wcp; !mass  
				gb_datalist(i)%red=de_zfromintpl(gb_datalist(i)%r) !calculate redshift
				i = i+1
			endif
			cycle
103			exit
104			print *, 'ERROR (read_in_HR3mockdata)! Error happes when read in ', trim(adjustl(gb_datafile)); stop
		enddo
		close(4401)
		if(i-1.ne.gb_numdata) then !consistency check
			print *, 'ERROR (read_in_HR3mockdata)! Mismatch of #: ', gb_numdata, i-1; stop
		endif
	end subroutine readin_mockdata


! Let's write a subrutin to read in BOSS data...
  !------------------------------------------
  ! read in the HR3 mock data
  !------------------------------------------
	subroutine readin_BOSSdata(printinfo)
		! Dummy
		logical, intent(in) :: printinfo
		! Local
		integer :: i, tot_numdata, numdata, ID
		real(dl) :: tmp(3), redshift,NOFKPwei,FKPwei,Completeness, FKPweiNew,WCP,WNOZ,WSTAR,WSEE

		tot_numdata = 0; numdata = 0

		if(allocated(gb_datalist)) deallocate(gb_datalist)
		open(unit=4401,file=gb_datafile)
		do while(.true.)		
       			read(4401,*,end=200,err=201) tmp(1:3), redshift
       			tot_numdata = tot_numdata+1
       			! We apply a redshift range cut 0.43-0.7
			if(redshift>gb_dataredcut_min .and. redshift<gb_dataredcut_max .and. &
				(.not.check_reject(tmp(1),tmp(2),tmp(3)))) numdata = numdata+1 
			cycle
200			exit
201			print *, 'ERROR (read_in_BOSSdata)! Error happes when read in ', trim(adjustl(gb_datafile)); stop
		enddo
		close(4401)
		
		gb_numdata = numdata
		if(printinfo) then
			write(*,'(A,A,A,2i8)') '   (read_in_BOSSdata) '
			write(*,'(8x,A)')         trim(adjustl(gb_datafile))
			write(*,'(25x,A,2i8)')     '#-data (infile/readin) = ', tot_numdata, gb_numdata
		endif

		i = 1
		allocate(gb_datalist(gb_numdata))
		open(unit=4401,file=gb_datafile)
		do while(.true.)		
       			read(4401,*,end=203,err=204) tmp(1:3), redshift,ID, NOFKPwei,FKPwei,Completeness,FKPweiNew,WCP,WNOZ,WSTAR,WSEE
       			! We apply a redshift range cut 0.43-0.7
			if(redshift>gb_dataredcut_min .and. redshift<gb_dataredcut_max .and. &
				(.not.check_reject(tmp(1),tmp(2),tmp(3)))) then
				gb_datalist(i)%x=tmp(1); gb_datalist(i)%y=tmp(2); gb_datalist(i)%z=tmp(3)
				gb_datalist(i)%r=rms(tmp(1:3),3) !these x,y,z,r are values in the fiducial cosmology!
				if(gb_wcpw.eq.gb_wcpw_usual) then
					WCP = WCP
				elseif(gb_wcpw.eq.gb_wcpw_eq1) then
					WCP = 1.0_dl
				elseif(gb_wcpw.eq.gb_wcpw_le1) then
					WCP = max(WCP,1.0_dl)
				endif
				gb_datalist(i)%mass= (WCP+WNOZ-1.0)/Completeness*WSTAR*WSEE
				gb_datalist(i)%red=de_zfromintpl(gb_datalist(i)%r) !calculate redshift
				if(abs(gb_datalist(i)%red-redshift).ge.1.0e-3)then
					print *, 'WARNING: Large diff in redshift: ', i, gb_datalist(i)%red, redshift
				endif
				i = i+1
			endif
			cycle
203			exit
204			print *, 'ERROR (read_in_HR3BOSSdata)! Error happes when read in ', trim(adjustl(gb_datafile)); stop
		enddo
		close(4401)
		if(i-1.ne.gb_numdata) then !consistency check
			print *, 'ERROR (read_in_BOSSdata)! Mismatch of #: ', gb_numdata, i-1; stop
		endif
	end subroutine readin_BOSSdata

  !------------------------------------------
  ! read in the random particles
  !------------------------------------------
	subroutine readin_HR3ran(printinfo)
		! Dummy
		logical, intent(in) :: printinfo
		! Local
		integer :: linenumber, i,i_ran,num_ranacpt, iAnder,iVeto,iRadecrange,iRedshiftrange
		real(dl) :: x,y,z,r,rext=100.0_dl

		! count how many randoms to read in (only read in randoms with r< gb_3dranmaxr+100.0)
		if(printinfo) write(*,'(A)',advance='no') '   (readin_HR3ran) Counting line-# of files...'
		gb_numran = 0
		do i_ran = 1, gb_numranfile
			open(unit=4401,file=gb_ranfilelist(i_ran))
			do while(.true.)
				read(4401,*,end=100) x,y,z
				r = sqrt(x**2.0+y**2.0+z**2.0)
				if(r>gb_3dranmaxr+rext) cycle
				gb_numran = gb_numran + 1
				cycle
100				exit
			enddo
			close(4401)
			write(*,'(1x,i8,$)') gb_numran
		enddo
		print *

		! read in the randoms
		if(allocated(gb_ranlist)) deallocate(gb_ranlist)
		allocate(gb_ranlist(gb_numran))
		i = 1; num_ranacpt=0
		do i_ran = 1, gb_numranfile
			if(printinfo.and.i_ran.eq.1) then
				write(*,'(A,A,A,2i8)') '                Loading: ran file = ',trim(adjustl(gb_ranfilelist(i_ran)))
			elseif(printinfo.and.i_ran.ne.1) then
				write(*,'(25x,A,A,A,2i8)') 'ran file = ',trim(adjustl(gb_ranfilelist(i_ran)))
			endif

			open(unit=4401,file=gb_ranfilelist(i_ran))
			do while(.true.)
				read(4401,*,end=105) x,y,z,iAnder,iVeto,iRadecrange,iRedshiftrange !fmt:x,y,z,criterias
				r = sqrt(x**2.0+y**2.0+z**2.0)
				if(r>gb_3dranmaxr+rext) cycle ! 
				gb_ranlist(i)%x=x; gb_ranlist(i)%y=y; gb_ranlist(i)%z=z; gb_ranlist(i)%r = r 
				! These x,y,z,r are values in the fiducial cosmology
				gb_ranlist(i)%red = de_zfromintpl(gb_ranlist(i)%r)
				!drop when redshift/radec out of range; * confirm dropping radec is correct *
				if(iRedshiftrange.le.0.5 .or. iRadecrange.le.0.5 & 
				   .or. (gb_VetoRej.and.iVeto.le.0.5) .or. (gb_AdsnRej.and.iAnder.le.0.5) & !drop veto/anderson
				   .or. check_reject(x,y,z)) then !drop due to some other reasons
					gb_ranlist(i)%acpt = .false.
				else
					gb_ranlist(i)%acpt = .true.
					num_ranacpt = num_ranacpt+1 
				endif
				i=i+1
				cycle
105				exit								
			enddo
			close(4401)
		enddo
		if(printinfo) then
			write(*,'(25x,A,2i10)')     '#-ran (readin/acpt) = ', gb_numran, num_ranacpt
		endif

		if(i-1.ne.gb_numran) then
			print *, 'ERROR (read_in_HR3ran)! Mismatch line #: ', i-1, gb_numran; stop
		endif
	end subroutine readin_HR3ran

	
  !------------------------------------------
  ! read in the radecrandom particles
  !------------------------------------------
	subroutine readin_HR3radecran(printinfo)
		! Dummy
		logical, intent(in) :: printinfo
		! Local
		integer :: linenumber, i,i_ran,num_ranacpt, iAnder,iVeto,iRadecrange
		real(dl) :: ra,dec, x,y,z

		! count how many radecrandoms to read in
		if(printinfo) write(*,'(A,A,A,2i8)') '   (readin_HR3radecran) Counting line-# of files...'
		gb_numradecran = 0
		do i_ran = 1, gb_numradecranfile
			call count_line_number(gb_radecranfilelist(i_ran), linenumber)
			gb_numradecran = gb_numradecran + linenumber
		enddo

		! read in the radecrandoms
		if(allocated(gb_radecranlist)) deallocate(gb_radecranlist)
		allocate(gb_radecranlist(gb_numradecran))
		i = 1; num_ranacpt=0
		do i_ran = 1, gb_numradecranfile
			if(printinfo.and.i_ran.eq.1) then
				write(*,'(A,A,A,2i8)') '                Loading: radecran file = ',trim(adjustl(gb_radecranfilelist(i_ran)))
			elseif(printinfo.and.i_ran.ne.1) then
				write(*,'(25x,A,A,A,2i8)') 'radecran file = ',trim(adjustl(gb_radecranfilelist(i_ran)))
			endif

			open(unit=4401,file=gb_radecranfilelist(i_ran))
			do while(.true.)
				read(4401,*,end=105) ra,dec,iAnder,iVeto,iRadecrange !fmt:x,y,z,criterias
				gb_radecranlist(i)%ra=ra; gb_radecranlist(i)%dec=dec
				call radecr_to_xyz(ra,dec,1.0_dl,x,y,z)
				! * confirm NOT reject iRadecrange .le. 0.5 *
				if(check_reject_radecran(ra,dec) &!.or. iRadecrange.le.0.5  & !drop when radec out of range or rejected
				   .or. (gb_VetoRej.and.iVeto.le.0.5) .or. (gb_AdsnRej.and.iAnder.le.0.5)) then !drop veto/anderson
					gb_radecranlist(i)%acpt = .false.
				else
					gb_radecranlist(i)%acpt = .true.
					num_ranacpt = num_ranacpt+1 
				endif
				i=i+1
				cycle
105				exit								
			enddo
			close(4401)
		enddo
		if(printinfo) then
			write(*,'(25x,A,2i10)')     '#-radecran (readin/acpt) = ', gb_numradecran, num_ranacpt
		endif

		if(i-1.ne.gb_numradecran) then
			print *, 'ERROR (read_in_HR3radecran)! Mismatch line #: ', i-1, gb_numradecran; stop
		endif
	end subroutine readin_HR3radecran

  !------------------------------------------
  ! get the min/max region of data/random
  !  not just fiducial values; depends on the cosmology
  !------------------------------------------
	subroutine get_minmaxs(printinfo)
		! Dummy
		logical, intent(in) :: printinfo
		! Local 
		integer :: i
		real(dl) :: ra, dec, r
		gbxmindata=logzero; gbxmaxdata=-logzero; gbymindata=logzero; gbymaxdata=-logzero; gbzmindata=logzero; gbzmaxdata=-logzero;
		gbxminran=logzero;  gbxmaxran=-logzero;  gbyminran=logzero;  gbymaxran=-logzero;  gbzminran=logzero;  gbzmaxran=-logzero;  
		gbrmindata=logzero; gbrmaxdata=-logzero; gbrminran=logzero;  gbrmaxran=-logzero;
	    	gbramindata=logzero; gbramaxdata=-logzero; gbdecmindata=logzero; gbdecmaxdata=-logzero; 
		gbraminran=logzero;  gbramaxran=-logzero;  gbdecminran=logzero;  gbdecmaxran=-logzero;
		gbraminradecran=logzero;  gbramaxradecran=-logzero;  gbdecminradecran=logzero;  gbdecmaxradecran=-logzero;

		do i = 1, gb_numdata
			gbxmindata = min(gbxmindata,gb_datalist(i)%x*gb_datalist(i)%rat)
			gbxmaxdata = max(gbxmaxdata,gb_datalist(i)%x*gb_datalist(i)%rat)
			gbymindata = min(gbymindata,gb_datalist(i)%y*gb_datalist(i)%rat)
			gbymaxdata = max(gbymaxdata,gb_datalist(i)%y*gb_datalist(i)%rat)
			gbzmindata = min(gbzmindata,gb_datalist(i)%z*gb_datalist(i)%rat)
			gbzmaxdata = max(gbzmaxdata,gb_datalist(i)%z*gb_datalist(i)%rat)
			gbrmindata = min(gbrmindata,gb_datalist(i)%r*gb_datalist(i)%rat)
			gbrmaxdata = max(gbrmaxdata,gb_datalist(i)%r*gb_datalist(i)%rat)
			call xyz_to_radecr(gb_datalist(i)%x,gb_datalist(i)%y,gb_datalist(i)%z,ra,dec,r)
			gbramindata = min(gbramindata,ra)
			gbramaxdata = max(gbramaxdata,ra)
			gbdecmindata = min(gbdecmindata,dec)
			gbdecmaxdata = max(gbdecmaxdata,dec)
		enddo
		do i = 1, gb_numran
			gbxminran = min(gbxminran,gb_ranlist(i)%x*gb_ranlist(i)%rat)
			gbxmaxran = max(gbxmaxran,gb_ranlist(i)%x*gb_ranlist(i)%rat)
			gbyminran = min(gbyminran,gb_ranlist(i)%y*gb_ranlist(i)%rat)
			gbymaxran = max(gbymaxran,gb_ranlist(i)%y*gb_ranlist(i)%rat)
			gbzminran = min(gbzminran,gb_ranlist(i)%z*gb_ranlist(i)%rat)
			gbzmaxran = max(gbzmaxran,gb_ranlist(i)%z*gb_ranlist(i)%rat)
			gbrminran = min(gbrminran,gb_ranlist(i)%r*gb_ranlist(i)%rat)
			gbrmaxran = max(gbrmaxran,gb_ranlist(i)%r*gb_ranlist(i)%rat)
			call xyz_to_radecr(gb_ranlist(i)%x,gb_ranlist(i)%y,gb_ranlist(i)%z,ra,dec,r)
			gbraminran = min(gbraminran,ra)
			gbramaxran = max(gbramaxran,ra)
			gbdecminran = min(gbdecminran,dec)
			gbdecmaxran = max(gbdecmaxran,dec)
		enddo
		do i = 1, gb_numradecran
			ra = gb_radecranlist(i)%ra; dec = gb_radecranlist(i)%dec
			gbraminradecran = min(gbraminradecran,ra)
			gbramaxradecran = max(gbramaxradecran,ra)
			gbdecminradecran = min(gbdecminradecran,dec)
			gbdecmaxradecran = max(gbdecmaxradecran,dec)
		enddo
		if(printinfo) then
			print *, '  (get_minmaxs)    Ranges. fmt:  (data), (random), (ra/dec random)'
			write(*,'(23x,A,2f10.3,A,2f10.3,A)') 'x:    (', real(gbxmindata), real(gbxmaxdata), &
				'), (', real(gbxminran), real(gbxmaxran), ')'
			write(*,'(23x,A,2f10.3,A,2f10.3,A)') 'y:   (', real(gbymindata), real(gbymaxdata), & 
				'), (', real(gbyminran), real(gbymaxran), ')'
			write(*,'(23x,A,2f10.3,A,2f10.3,A)') 'z:   (', real(gbzmindata), real(gbzmaxdata), &
				'), (', real(gbzminran), real(gbzmaxran), ')'
			write(*,'(23x,A,2f10.3,A,2f10.3,A)') 'r:   (', real(gbrmindata), real(gbrmaxdata), &
				'), (', real(gbrminran), real(gbrmaxran), ')'
			write(*,'(23x,A,2f10.3,A,2f10.3,A,2f10.3,A)') 'ra:  (', real(gbramindata), real(gbramaxdata), &
				'), (', real(gbraminran), real(gbramaxran), &
				'), (', real(gbraminradecran), real(gbramaxradecran), ')'
			write(*,'(23x,A,2f10.3,A,2f10.3,A,2f10.3,A)') 'dec: (', real(gbdecmindata), real(gbdecmaxdata), &
				'), (', real(gbdecminran), real(gbdecmaxran) , &
				'), (', real(gbdecminradecran), real(gbdecmaxradecran), ')'
		endif
	end subroutine get_minmaxs
	
  !------------------------------------------
  ! initializing the halo_info array
  !------------------------------------------
	subroutine readin_dataran(printinfo)
		! Dummy
		logical, intent(in) :: printinfo
		! local
		integer :: i,j,k
		character(len=char_len) :: tmpstr1, tmpstr2

		! initialize cosmology using fiducial values
		gb_omegam = om_dft; gb_w = w_dft; gb_h = h_dft
		if(printinfo) then
			print *, '  (readin_dataran) Loading and initializing data and randoms ...'
			print *, '                      datatype / Anderson-Rej / Veto-Rej = ', gb_i_datatype, gb_AdsnRej, gb_VetoRej
			write(*,'(23x,A,f8.3,A,f8.3,A,f8.3)') 'Fiducial Cosmology: Omegam/w/h = ',&
				 gb_omegam,' /',gb_w,' /',gb_h
			write(*,'(23x,A,A,A)') 'Catalogue = "', trim(adjustl(gb_catalogue_name)), '"'
			write(*,'(23x,A,1f12.3)') '   Effective-Area: ', gb_Seff
			write(*,'(23x,A,2f12.3,A)') '   Redshift range: (',gb_dataredcut_min,gb_dataredcut_max, ')'
			write(*,'(23x,A,2f12.3,A)') '   Distance range: (', gb_minimalr_cut, gb_maximalr_cut, ')'
		endif
		call cosmo_funs_init(printinfo)
		call de_calc_comovr()
		
		if(gb_i_datatype.eq.0) then
			if(gb_AdsnRej .or. gb_VetoRej) then
				print *, 'ERROR (readin_dataran)! Do not use Adsn-Rej or Veto-Rej when idatatype=0 (constant nbar, no mask data)!'; stop
			endif
		endif

		! read in data, random...
		gb_numdata = 0
		gb_numran = 0
		gb_numradecran = 0
		gb_numallpar = 0


		if(gb_i_datatype.eq.gb_dt_realdata) then
			call readin_BOSSdata(printinfo)
		elseif(gb_i_datatype.eq.gb_dt_xyzw) then
			call readin_xyzwdata(printinfo,hasweight=.true.)
		elseif(gb_i_datatype.eq.gb_dt_xyz) then
			call readin_xyzwdata(printinfo,hasweight=.false.)
		else
			call readin_mockdata(printinfo)
		endif
		
		if(gb_i_rantype.eq.gb_3dran) then
			call readin_HR3ran(printinfo)
		elseif(gb_i_rantype.eq.gb_radecran) then
			call readin_HR3radecran(printinfo)
		elseif(gb_i_rantype.eq.gb_radec3dran) then
			call readin_HR3ran(printinfo)
			call readin_HR3radecran(printinfo)
		elseif(gb_i_rantype.eq.gb_noran) then
			continue
		else
			print *, 'ERROR (readin_dataran)!!! gb_consider_ran must be one of ', &
				gb_3dran, gb_radecran, gb_noran, gb_radec3dran
			stop
		endif		

		! total number of data+random
		gb_numallpar = gb_numdata + gb_numran
				
		! get ranges of x/y/z
		call get_minmaxs(printinfo)
	end subroutine readin_dataran

  !------------------------------------------
  ! use given parameters to init cosmology
  !------------------------------------------		
	subroutine init_cosmo(omegam,w,h,printinfo)
		! Dummy
		real(dl), intent(in) :: omegam, w, h
		logical, intent(in) :: printinfo
		! Local
		integer :: i, j
		real(dl) :: r, ra, dec
		
		if(printinfo) &
			write(*,'(A,f7.3,f7.3)') '   (init_LSS_cosmo) Initializing using omegam/w = ', real(omegam), real(w)
		gb_omegam 	= omegam
		gb_w 		= w
		gb_h 		= h
		call de_calc_comovr()

		! Calculate rat between current-radius/original-radius for data, random
		do i = 1, gb_numdata
			gb_datalist(i)%rat = de_get_comovr(gb_datalist(i)%red) / gb_datalist(i)%r
		enddo		
		do i = 1, gb_numran
			gb_ranlist(i)%rat = de_get_comovr(gb_ranlist(i)%red) / gb_ranlist(i)%r
		enddo
		
		! Update the range of x,y,z,r
		call get_minmaxs(printinfo)
	end subroutine init_cosmo		
	
!!! Test subroutines

  !------------------------------------------
  ! used to test the readin data
  !------------------------------------------	
	subroutine test_readindata()
		integer :: idata
		idata = 1000
		print *, 'info of ', idata,'th data:'
		print *, gb_datalist(idata)
		print *, 'info of ', idata,'th random:'
		print *, gb_ranlist(idata)
		do while(.true.)
			print *, 'input your omegam, w, idata...'
			read(*,*) gb_omegam, gb_w, idata
			call init_cosmo(gb_omegam,gb_w,gb_h,.true.)
!			print *, 'info of ', idata,'th data:'
!			print *, gb_datalist(idata)
			print *, 'info of ', idata,'th random:'
			print *, gb_ranlist(idata)
		enddo
	end subroutine test_readindata
end module LSS_settings_init
