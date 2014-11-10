program ap_main_test

use ap_settings_init
implicit none

	character(len=char_len) :: dataname

	gb_HR3format = .true.
	gb_AdsnRej = .true.; gb_VetoRej = .true.

        dataname = 'Xiao-dong.00000.dat.compact.patch1.noRSD-radial-selected'
	gb_datafile = '~/SparseFilaments/data/input/HR3/DR12_mock/DR12v1-CMASS-N/'//trim(adjustl(dataname))

	gb_numranfile = 0
	call init_dataran(.true.)

end program ap_main_test
