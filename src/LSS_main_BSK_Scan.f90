

program LSS_main

use mpi
use LSS_BSK
implicit none

	! Settings
	character(len=char_len) :: ngldir = '~/SparseFilaments/code/ngl-beta/bin/getNeighborGraph'
	integer :: numNNB = 100, numdrop=1	
	real(dl) :: dropstep=0.1
	logical :: printinfo = .true., do_nglcrosscheck = .false.

	! Useful variables
	real(dl) :: beta, omegam, w, ommin, ommax, wmin,wmax
	real(dl), allocatable :: om_w_list(:,:)
	integer :: numarg, nproc, ierr, myid, i,j, imock,idatatype,ibeta,iran,numran,nummock, num_om, num_w, num_omw
	integer, allocatable :: iomlist(:), iwlist(:)
	character(len=char_len) :: tmpstr1, tmpstr2, suffixstr
	
	! Initialize MPI
	call mpi_init(ierr)
	call mpi_comm_size(mpi_comm_world,nproc,ierr)
	call mpi_comm_rank(mpi_comm_world,myid,ierr)

	!------------------------
	! suffix string of the run
		suffixstr = '' 
	!------------------------
	
    do ibeta = 5,5
    do idatatype =  0,1
        ! idatatype = 0: noRSD data
    	! idatatype = 1: RSD data
    	! idatatype = 2: random
	nummock = 0
    do imock = 1, 1
    	nummock = nummock + 1
    	numran = 0
    do iran = 0, 0
    	numran = numran + 1

    	if(idatatype <=1 .and. numran .ne.1) cycle ! for one mock, will not redo calculation if numran > 1
    	if(idatatype >=2 .and. nummock .ne.1) cycle ! only consider random for the first mock

	print *
 	print *, '###########################################################################'  
 	write(*,'(3(A,i3,";"))'), ' imock = ',imock, ' idatatype=',idatatype, &
 		'iran=',iran
 	print *, '###########################################################################'
 	print * 

!##########################################################
!######  Sec A. Settings ##################################
!##########################################################

  !-------------------------------------------------------
  ! A0. Core Varialbes
  !-------------------------------------------------------

	!------------------------
	! values of beta
		beta = ibeta
	!------------------------
	
	!------------------------    	
    	! numdrop; dropstep 
    	!!! ignoring short connections in the connectin; this tells how many ignoring / step size of ignoring
    		numdrop = 2
    		dropstep = 0.5_dl
    	!------------------------	
    	
  !-------------------------------------------------------
  ! A1. Important Varialbes
  !-------------------------------------------------------
	!------------------------    	
    	! minimal/maximal cut applied to the data
    		gb_minimalr_cut = 100.0_dl
   	 	gb_maximalr_cut = 3000.0_dl
	!------------------------
	
	!------------------------    	
    	! # of NNB: 
    	!   beta = 10.0: 10-15 shall be enough
	!    beta = 5.0: 15-20 shall be enough; 
	!    beta = 3.0: 25 shall be enough
	!    beta = 1.0: 100 shall be enough;
	    	if(beta< 2.0) then
	    		numNNB = 100
	    	elseif(beta < 3.0) then
	    		numNNB = 25
	    	elseif(beta < 5.0) then
	    		numNNB = 20
	    	else
	    		numNNB = 15
	    	endif
 	!------------------------

	!------------------------	
    	! print information
		printinfo = .true. 
	! cross check using ngl
		do_nglcrosscheck = .false.
		ngldir = '~/SparseFilaments/code/ngl-beta/bin/getNeighborGraph'
  	!------------------------
  	  	
  !-------------------------------------------------------
  ! A2. Data, Random
  !	Names.
  !-------------------------------------------------------
	gb_i_datatype = -1 ! -1 means x,y,z,weight data
    	write(tmpstr1,*) imock
    	if(idatatype.eq.0) then
	    	gb_datafile = 'HR3WS_0to3100_'//trim(adjustl(tmpstr1))//'.txt.noRSD'
	elseif(idatatype.eq.1) then
	    	gb_datafile = 'HR3WS_0to3100_'//trim(adjustl(tmpstr1))//'.txt.RSD'
	else
		write(tmpstr2, *) iran
		gb_datafile = 'HR3WS_0to3100_'//trim(adjustl(tmpstr1))//'.txt.1ran'//trim(adjustl(tmpstr2))
	endif
  !-------------------------------------------------------
  ! A3. cosmological parameters
  !-------------------------------------------------------
	
!	num_om = 2; num_w = 1; ommin = 0.08_dl; ommax = 1.00_dl; wmin = -1.0_dl; wmax = -2.5_dl;
!	num_omw = num_om*num_w; allocate(om_w_list(2,num_omw),iomlist(num_omw), iwlist(num_omw))

!	num_omw = 6; allocate(om_w_list(2,num_omw),iomlist(num_omw), iwlist(num_omw))
!	om_w_list(1,1:6) = (/0.26_dl, 0.26_dl, 0.26_dl, 0.00_dl, 1.00_dl/)
!	om_w_list(2,1:6) = (/-1.0_dl, -0.5_dl, -3.0_dl, -1.0_dl, -1.0_dl/)

	num_omw = 1; allocate(om_w_list(2,num_omw),iomlist(num_omw), iwlist(num_omw))
	om_w_list(1,1:1) = (/ 0.00_dl /)!, 1.00_dl /) !, 0.26_dl, 0.00_dl, 1.00_dl/)
	om_w_list(2,1:1) = (/-1.0_dl /)!, -1.0_dl /) !, -3.0_dl, -1.0_dl, -1.0_dl/)

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

!##########################################################
!######  Sec B. Main body #################################
!##########################################################
  !-------------------------------------------------------
  ! B0. Create a name for the run
  !-------------------------------------------------------
     !## beta
     	write(tmpstr1,'(f10.1)') beta
	gb_output_name = trim(adjustl(gb_datafile))//'_BSKrlt/Beta'//trim(adjustl(tmpstr1))

      !## range of r 
	write(tmpstr1,'(f10.1)') gb_minimalr_cut
	write(tmpstr2,'(f10.1)') gb_maximalr_cut
	gb_output_name = trim(adjustl(gb_output_name))//'_r'//trim(adjustl(tmpstr1))//'to'//trim(adjustl(tmpstr2))
	
     !## numNNB
     	write(tmpstr1,*) numNNB
	gb_output_name = trim(adjustl(gb_output_name))//'_numNNB'//trim(adjustl(tmpstr1))

     !## dropping
     	write(tmpstr1,*) numdrop
     	write(tmpstr2,'(f5.1)') dropstep
	gb_output_name = trim(adjustl(gb_output_name))//'_'//trim(adjustl(tmpstr1))//'dropStep'//trim(adjustl(tmpstr2))
	
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
		print *, 'Default cosmology: omegam, w = ', om_dft, w_dft
		print *, '#####################################'
		print *, ' Settings of BSK'
		print *, '  gb_datafile:  ', trim(adjustl(gb_datafile))
		print *, '  gb_output_name: ', trim(adjustl(gb_output_name))
		print *, '  beta = ', beta
		print *, '  printinfo = ', printinfo
		print *, '  #-NNB searching for BSK: ', numNNB
		print *, '  gb_minimalr_cut: ', gb_minimalr_cut
		print *, '  gb_maximalr_cut: ', gb_maximalr_cut
		print *, '  do_nglcrosscheck: ', do_nglcrosscheck
		print *, '  ngldir: ', trim(adjustl(ngldir))
		print *, '  #-drop: ', numdrop
		print *, '  dropstep: ', dropstep
     		write(*,'(A,A)'), '   output directory: ', trim(adjustl(gb_output_name))
		print *, '#####################################'

		print *
		write(*,'(10x,A,$)') 'Lists of omegam/w: '
		do i = 1, num_omw
			write(*,'(f6.3,"/",f6.3,$)') om_w_list(1:2,i)
		enddo
		print *
		print *
	endif
	
  !-------------------------------------------------------
  ! B2. Do initializations calculation
  !-------------------------------------------------------

	gb_omegam = om_dft; gb_w = w_dft  
  	call cosmo_funs_init(printinfo)
	call readin_dataran(printinfo)
	
	do i =  myid+1, num_omw, nproc
		omegam = om_w_list(1,i)
		w  = om_w_list(2,i)
		gb_iom = iomlist(i); gb_iw = iwlist(i)

		write(tmpstr1,'(f8.4)') omegam
		write(tmpstr2,'(f8.4)') w
		gb_omwstr = 'om'//trim(adjustl(tmpstr1))//'_w'//trim(adjustl(tmpstr2))
		gb_suboutput_name = trim(adjustl(gb_output_name))//'/'//trim(adjustl(gb_omwstr))
		
		call init_cosmo(omegam,w,h_dft,printinfo)
		call init_mult_lists(printinfo)
		call do_cell_init((real(gb_numdata)**0.33_dl), printinfo)	

		call BSK(gb_datafile, omegam, w, beta, gb_output_name, printinfo, numNNB, gb_minimalr_cut, gb_maximalr_cut, &
			do_nglcrosscheck, ngldir, outputBSKindex=.true., outputBSKinfo=.false., outputxyz=.false.,do_init=.false.)

		call BSK_stat(omegam, w, beta, gb_output_name, numdrop, dropstep, printinfo)
	enddo
	deallocate(om_w_list,iomlist, iwlist)
     enddo
     enddo
     enddo
     enddo
	call mpi_barrier(mpi_comm_world,ierr)
	call mpi_finalize(ierr)
end program LSS_main
