module test_funs

use LSS_tools
implicit none


end module test_funs

program LSS_main_test

use test_funs
implicit none

	character(len=char_len) :: tmpstr1, tmpstr2, tmpstr

	real(dl) :: A(100,3)

	integer :: i

	call random_number(A)

	tmpstr1 = 'old.txt'
	call output_2d(tmpstr1, A)

	open(unit=1,file='new2.txt')
	do i = 1, 100
		tmpstr = WriteFmtEFloat(3, A(i,1:3))
		write(1,'(A)') trim(adjustl(tmpstr))
		print *, trim(adjustl(tmpstr))
	enddo
	close(1)

end program LSS_main_test


