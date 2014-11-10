
! What treatment we have adopted now:


program ap_main_mpi_sh

use mpi
use ap_tools

	implicit none
	character(len=char_len) :: cmdfile, printstr, tmpstr1, tmpstr2, cmdstr
        character(len=10000), allocatable :: cmds(:)
	integer :: numarg, num_cmd, i

	! mpi variables
	integer :: ierr, nproc, myid

	printstr = '### Running shell commands using mpi ### Usage: mpirun -np ??? EXE -cmdfile cmdfile'

	call mpi_init(ierr)
	call mpi_comm_size(mpi_comm_world,nproc,ierr)
	call mpi_comm_rank(mpi_comm_world,myid,ierr)

	numarg = iargc()
	if(numarg.le.1) then
		write(*,'(A)') trim(adjustl(printstr))
		stop
	endif

	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-cmdfile') then
			read(tmpstr2,'(A)') cmdfile
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			write(*,'(A)') trim(adjustl(printstr))
			stop
		endif
	enddo

	call mpi_barrier(mpi_comm_world,ierr)
        call count_line_number (cmdfile, num_cmd)
        allocate(cmds(num_cmd))
        open(unit=100,file=cmdfile)
        do i = 1, num_cmd
                read(100,'(A)') cmds(i)
        enddo
        close(100)
               
	call mpi_barrier(mpi_comm_world,ierr)
        do i = 1, num_cmd
		cmdstr = cmds(i)
		if(mod(i,nproc).eq.myid) then
			write(*,'(A,i5,A,i4,A,A)') 'Running ', i, 'th command on process ', myid, ':   ', trim(adjustl(cmdstr))
			call system(trim(adjustl(cmdstr)))
		endif
	enddo

	call mpi_barrier(mpi_comm_world,ierr)
	call mpi_finalize(ierr)

end program ap_main_mpi_sh
