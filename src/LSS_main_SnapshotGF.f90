
!!! Not finished;... Unexpected stack error in grid_rho_drho;...!!!! Need more test...

program main

use LSS_chisq

implicit none

	real(dl) ::omAP,wAP, muav,muer
	real(dl), allocatable :: drho_mu_data(:), pos_list(:,:)
	integer :: numarg, i, n, ixyz
	character(len=char_len) :: printstr, inputfilename, outputfilename, tmpstr1, tmpstr2
	logical :: hasweight, printinfo
	type(chisq_settings) :: cs
	
	printstr = 'Usage: EXE -input intpufilename -hasweight hasweight '//&
		'-bddist bddist -omdft omdft -wdft wdft -omAP omAP -wAP wAP '//&
		'-printinfo printinfo -usefixmd usefixmd -smnum smnum '//&
		'-fixmd fixmd -numinx numinx '//&
		'### Must be fmt of x,y,z, weight! By default only x,y,z '//&
		'## Keep a distance bddist from boundary, to avoid boundary effect'

	numarg = iargc()
	if(numarg.le.1) then
		write(*,'(A)') trim(adjustl(printstr))
		stop
	endif

	! Default values
	hasweight = .false.
	printinfo = .true.
	gb_bddist = 30.0_dl
	cs%smnum = 20
	gb_use_fixmd = .false.; gb_fixmd = 15.0*2.0 
	cs%num_in_x = 0
	omAP =0.26; wAP = -1.0
	
	! Important settings
	gb_minimalr_cut = -1.0e30
	gb_maximalr_cut = 1.0e30
	gb_i_datatype = gb_dt_xyzw ! for xyzweight fmt

	gb_i_rantype = gb_noran ! no random
	gb_ranwei_tol = -1.0 ! no tolerance; all pass
	gb_dodensitynorm = .false. ! no density normalization ! very important!!!
	gb_num_changenuminx = 0 ! no re-computing using different sizes of grid
	gb_bddist_rextra = 0; gb_bufferdist = 0; gb_bufferdist_rextra = 0 ! no extra bddist for r boundary; no buffer

	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-input') then
			read(tmpstr2,'(A)') inputfilename
		elseif(trim(adjustl(tmpstr1)).eq.'-hasweight') then
			read(tmpstr2,*) hasweight
		elseif(trim(adjustl(tmpstr1)).eq.'-omdft') then
			read(tmpstr2,*) om_dft
		elseif(trim(adjustl(tmpstr1)).eq.'-wdft') then
			read(tmpstr2,*) w_dft
		elseif(trim(adjustl(tmpstr1)).eq.'-omAP') then
			read(tmpstr2,*) omAP
		elseif(trim(adjustl(tmpstr1)).eq.'-wAP') then
			read(tmpstr2,*) wAP
		elseif(trim(adjustl(tmpstr1)).eq.'-bddist') then
			read(tmpstr2,*) gb_bddist
		elseif(trim(adjustl(tmpstr1)).eq.'-usefixmd') then
			read(tmpstr2,*) gb_use_fixmd
		elseif(trim(adjustl(tmpstr1)).eq.'-fixmd') then
			read(tmpstr2,*) gb_fixmd
		elseif(trim(adjustl(tmpstr1)).eq.'-smnum') then
			read(tmpstr2,*) cs%smnum
		elseif(trim(adjustl(tmpstr1)).eq.'-numinx') then
			read(tmpstr2,*) cs%num_in_x
		elseif(trim(adjustl(tmpstr1)).eq.'-printinfo') then
			read(tmpstr2,*) printinfo
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	if(hasweight) then 
		gb_i_datatype = gb_dt_xyzw
	else
		gb_i_datatype = gb_dt_xyz
	endif
	
	if(gb_use_fixmd) then
		write(tmpstr1,'(f4.1)') gb_fixmd/2.0
	        outputfilename      = trim(adjustl(inputfilename))//'_smh'//trim(adjustl(tmpstr1)) ! redshift cut 0.36
	else
	        write(tmpstr1,*) cs%smnum
		outputfilename     = trim(adjustl(inputfilename))//'_smnum'//trim(adjustl(tmpstr1))
	endif
	
	write(tmpstr1,'(f8.4)') omAP
	write(tmpstr2,'(f8.4)') wAP
	gb_omwstr = 'om'//trim(adjustl(tmpstr1))//'_w'//trim(adjustl(tmpstr2))
	write(tmpstr1,*) cs%smnum
	outputfilename = trim(adjustl(outputfilename))//'_'//trim(adjustl(gb_omwstr))//'.muinfo'

	print *, '#####################################'
	write(*,'(A)')   	'  Settings:'
	write(*,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(*,'(A,3L)') 	'   hasweigt, printinfo, usefixmd = ', hasweight, printinfo, gb_use_fixmd
	write(*,'(A,2i4,2f16.7)')'   gridsize, smnum, fixmd, bddist = ', cs%num_in_x, cs%smnum, gb_fixmd, gb_bddist
	write(*,'(A,2f16.7)') 	'   omegam, w (dft) = ', real(om_dft), real(w_dft)
	write(*,'(A,2f16.7)') 	'   omegam, w (AP)  = ', real(omAP), real(wAP)
	write(*,'(A,A)') 	'   outputfilename  = ', trim(adjustl(outputfilename))
	print *, '#####################################'
	
	! initialization
	print *, '  (SanpshotGF) Begins.'
	gb_datafile = inputfilename
	call cosmo_funs_init(printinfo)
	call readin_dataran(printinfo)
	call init_cosmo(omAP,wAP,h_dft,printinfo)
	call init_mult_lists(printinfo)
	if(cs%num_in_x.eq.0) cs%num_in_x = dble(gb_numdata)**0.33
	call do_cell_init(real(cs%num_in_x)+0.0_dl, printinfo)		
	
	call grid_rho_drho_list(cs%smnum, printinfo, gb_pos_list, gb_rho_list, gb_drho_list)
	n = size(gb_rho_list); deallocate(gb_rho_list)
	allocate(pos_list(3,n))
	open(unit=100,file=outputfilename)
	
	do ixyz = 1, 3
		pos_list = 0
		do i = 1, n
			pos_list(ixyz,i) = 1.0_dl
		enddo
		call get_mu_from_gradient_list(pos_list, gb_drho_list, drho_mu_data)
		call get_mean_var(drho_mu_data, muav, muer)
		write(*,'(A,i3,A)') 'ixyz = ', ixyz, '  (1 means observed from x direction; 2 means y; 3 means z):'
		write(100,'(A,i3,A)') 'ixyz = ', ixyz, '  (1 means observed from x direction; 2 means y; 3 means z):'
		write(*,'(3e15.7,A)') muav, muer, ((muav)/muer)**2.0, '  #  Without abs. muav, muer, chisq. '
		write(100,'(3e15.7,A)') muav, muer, ((muav)/muer)**2.0, '  #  Without abs. muav, muer, chisq. '
		do i = 1, n
			drho_mu_data(i) = abs(drho_mu_data(i))
		enddo
		call get_mean_var(drho_mu_data, muav, muer)
		muer = muer / sqrt(dble(n))
		write(*,'(3e15.7,A)') muav, muer, ((muav)/muer)**2.0, '  #  With abs. muav, muer, chisq. '
		write(100,'(3e15.7,A)') muav, muer, ((muav)/muer)**2.0, '  #  With abs. muav, muer, chisq. '
		deallocate(drho_mu_data)
	enddo
	deallocate(pos_list)
	
	write(100,'(A)') '#####################################'
	write(100,'(A)')   	'  Settings:'
	write(100,'(A,A)') 	'   inputfilename:  ', trim(adjustl(inputfilename))
	write(100,'(A,3L)') 	'   hasweigt, printinfo, usefixmd = ', hasweight, printinfo, gb_use_fixmd
	write(100,'(A,2i4,2f16.7)')'   gridsize, smnum, fixmd, bddist = ', cs%num_in_x, cs%smnum, gb_fixmd, gb_bddist
	write(100,'(A,2f16.7)') 	'   omegam, w (dft) = ', real(om_dft), real(w_dft)
	write(100,'(A,2f16.7)') 	'   omegam, w (AP)  = ', real(omAP), real(wAP)
	write(100,'(A,A)') 	'   outputfilename  = ', trim(adjustl(outputfilename))
	write(100,'(A,i10)') 	'   #-data  = ', gb_numdata
	write(100,'(A,i5)') 	'   size-of-grid  = ', cs%num_in_x
	write(100,'(A)') '#####################################'
	print *, ' (SnapshotGF) Done.'
end program main
