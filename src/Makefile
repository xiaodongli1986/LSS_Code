EXE0 =  ../bin/LSS_test
EXE1 =  ../bin/LSS_calc_grad_chisqs 
EXE2 =  ../bin/LSS_BSK_calc
EXE3 =  ../bin/LSS_BSK_Scan
EXE4 =  ../bin/LSS_rmass_stat
EXE5 =  ../bin/LSS_SnapshotRSDmock
EXE6 =  ../bin/LSS_SnapshotGF
EXE7 =  ../bin/LSS_MakeLightCone
EXE8 =  ../bin/LSS_MakeCompactSample
EXE9 =  ../bin/LSS_CosmoConvert
EXE10 = ../bin/LSS_Make-1eighth-Shell-Random
EXE11 = ../bin/LSS_Rbinned-Split-Data
EXE12 = ../bin/LSS_mpi_sh
OBJS = LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_settings_init.o LSS_cell.o LSS_smooth.o LSS_grad_fields.o LSS_chisq.o LSS_BSK.o
F90C = mpif90 -heap-arrays #gfortran #mpif90 -heap-arrays
#F90C = gfortran
#F90C = /home/xiaodongli/software/anaconda/bin/mpif90
F90FLAGS = 

default: $(EXE2)
nompi: $(EXE0)         $(EXE2)         $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11)
all:   $(EXE0) $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) 
test: $(EXE0)
bsk: $(EXE2)
bskscan:  $(EXE3)
rmass: $(EXE4)
SnapshotRSDmock: $(EXE5)
SnapshotGF: $(EXE6)
MakeLightCone: $(EXE7)
MakeCompactSample: $(EXE8)
CosmoConvert: $(EXE9)
MakeRandom: $(EXE10)
Split-Data: $(EXE11)
mpi_sh: $(EXE12)
$(EXE0): $(OBJS) LSS_main_test.o
	$(F90C) -o $(EXE0) $(OBJS) LSS_main_test.o
$(EXE1): $(OBJS) LSS_main_calc_grad_chisqs.o
	$(F90C) -o $(EXE1) $(OBJS) LSS_main_calc_grad_chisqs.o
$(EXE2): $(OBJS) LSS_main_BSK_calc.o
	$(F90C) -o $(EXE2) $(OBJS) LSS_main_BSK_calc.o
$(EXE3): $(OBJS) LSS_main_BSK_Scan.o
	$(F90C) -o $(EXE3) $(OBJS) LSS_main_BSK_Scan.o
$(EXE4): $(OBJS) LSS_main_r_mass_stat.o
	$(F90C) -o $(EXE4) $(OBJS) LSS_main_r_mass_stat.o
$(EXE5): $(OBJS) LSS_main_SnapshotRSDmock.o
	$(F90C) -o $(EXE5) $(OBJS) LSS_main_SnapshotRSDmock.o
$(EXE6): $(OBJS) LSS_main_SnapshotGF.o
	$(F90C) -o $(EXE6) $(OBJS) LSS_main_SnapshotGF.o
$(EXE7): $(OBJS) LSS_main_MakeLightCone.o
	$(F90C) -o $(EXE7) $(OBJS) LSS_main_MakeLightCone.o
$(EXE8): $(OBJS) LSS_main_MakeCompactSample.o
	$(F90C) -o $(EXE8) $(OBJS) LSS_main_MakeCompactSample.o
$(EXE9): LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_main_CosmoConvert.o
	$(F90C) -o $(EXE9) LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_main_CosmoConvert.o
$(EXE10): LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_main_make-1eighth-shell-random.o
	$(F90C) -o $(EXE10) LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o  LSS_main_make-1eighth-shell-random.o
$(EXE11): LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_main_rbinned-split-data.o
	$(F90C) -o $(EXE11) LSS_tools.o  LSS_mu_tools.o LSS_cosmo_funs.o LSS_main_rbinned-split-data.o
$(EXE12): LSS_tools.o LSS_main_mpi_sh.o
	$(F90C) -o $(EXE12) LSS_tools.o  LSS_main_mpi_sh.o

LSS_cosmo_funs.o: LSS_tools.o
LSS_mu_tools.o: LSS_tools.o
LSS_settings_init.o: LSS_cosmo_funs.o 
LSS_cell.o: LSS_settings_init.o
LSS_smooth.o: LSS_cell.o
LSS_grad_fields.o: LSS_smooth.o LSS_mu_tools.o
LSS_chisq.o: LSS_grad_fields.o 
LSS_main_test.o:  LSS_chisq.o 
LSS_main_calc_grad_chisqs.o:  LSS_chisq.o 
LSS_BSK.o : LSS_chisq.o
LSS_main_BSK_calc.o: LSS_BSK.o 
LSS_main_test.o: LSS_settings_init.o
LSS_main_BSK_Scan.o: LSS_BSK.o
LSS_main_r_mass_stat.o: LSS_tools.o
LSS_main_SnapshotRSDmock.o: LSS_cosmo_funs.o
LSS_main_SnapshotGF.o: LSS_chisq.o
LSS_main_MakeLightCone.o: LSS_settings_init.o
LSS_main_CosmoConvert.o: LSS_cosmo_funs.o
LSS_main_make-1eighth-shell-random.o: LSS_cosmo_funs.o
LSS_main_rbinned-split-data.o: LSS_cosmo_funs.o
LSS_main_mpi_sh.o: LSS_tools.o
#LSS_main.o: LSS_settings_init.o
#LSS_main.o: LSS_tools.o
#LSS_main.o: LSS_cosmo_funs.o

%.o: %.f90
	$(F90C) -c $*.f90 $(F90FLAGS)

clean :
	rm -rf *.o *.mod $(EXE0) $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(OBJS)
