module test_funs

use LSS_tools
implicit none
contains

	! I tried to format output float, but, still fails...
	character(len=char_len) function WriteFmtEFloat(num_float, floatvalues, E1_opt, E2_opt, num_space_opt)
		! Dummy
		integer, intent(in) :: num_float
		real(dl), intent(in) :: floatvalues(num_float)
		integer, intent(in), optional :: E1_opt, E2_opt, num_space_opt
		! Local
		integer :: E1, E2, num_space
		character(len=char_len) :: fmtstr, tmpstr
		
		E1 = 14; E2 = 7; num_space = 2
		if(present(E1_opt)) E1=E1_opt
		if(present(E2_opt)) E2=E2_opt
		if(present(num_space_opt)) num_space = num_space_opt

		write(fmtstr, '(A,i3,A,i3,A,i3,A,i3,A)') '''', num_float, '(e', E1, '.', E2, ',', num_space, 'x)'''
		print *, fmtstr

		write(WriteFmtEFloat, fmtstr) floatvalues

	end function WriteFmtEFloat

end module test_funs

program LSS_main_test

use test_funs
implicit none

	character(len=char_len) :: tmpstr

	real(dl) :: A(3)

	call random_number(A)

	tmpstr = WriteFmtEFloat(3, A)

end program LSS_main_test


