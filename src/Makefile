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
OBJS = ap_tools.o  ap_mu_tools.o ap_cosmo_funs.o ap_settings_init.o ap_cell.o ap_smooth.o ap_grad_fields.o ap_chisq.o ap_BSK.o
F90C = mpif90 -heap-arrays
F90FLAGS = 

default: $(EXE1)
all: $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12)
bsk: $(EXE2) $(EXE3)
rmass: $(EXE4)
SnapshotRSDmock: $(EXE5)
SnapshotGF: $(EXE6)
MakeLightCone: $(EXE7)
MakeCompactSample: $(EXE8)
CosmoConvert: $(EXE9)
MakeRandom: $(EXE10)
Split-Data: $(EXE11)
mpi_sh: $(EXE12)
$(EXE1): $(OBJS) ap_main_calc_grad_chisqs.o
	$(F90C) -o $(EXE1) $(OBJS) ap_main_calc_grad_chisqs.o
$(EXE2): $(OBJS) ap_main_BSK_calc.o
	$(F90C) -o $(EXE2) $(OBJS) ap_main_BSK_calc.o
$(EXE3): $(OBJS) ap_main_BSK_Scan.o
	$(F90C) -o $(EXE3) $(OBJS) ap_main_BSK_Scan.o
$(EXE4): $(OBJS) ap_main_r_mass_stat.o
	$(F90C) -o $(EXE4) $(OBJS) ap_main_r_mass_stat.o
$(EXE5): $(OBJS) ap_main_SnapshotRSDmock.o
	$(F90C) -o $(EXE5) $(OBJS) ap_main_SnapshotRSDmock.o
$(EXE6): $(OBJS) ap_main_SnapshotGF.o
	$(F90C) -o $(EXE6) $(OBJS) ap_main_SnapshotGF.o
$(EXE7): $(OBJS) ap_main_MakeLightCone.o
	$(F90C) -o $(EXE7) $(OBJS) ap_main_MakeLightCone.o
$(EXE8): $(OBJS) ap_main_MakeCompactSample.o
	$(F90C) -o $(EXE8) $(OBJS) ap_main_MakeCompactSample.o
$(EXE9): ap_tools.o  ap_mu_tools.o ap_cosmo_funs.o ap_main_CosmoConvert.o
	$(F90C) -o $(EXE9) ap_tools.o  ap_mu_tools.o ap_cosmo_funs.o ap_main_CosmoConvert.o
$(EXE10): ap_tools.o  ap_mu_tools.o ap_cosmo_funs.o ap_main_make-1eighth-shell-random.o
	$(F90C) -o $(EXE10) ap_tools.o  ap_mu_tools.o ap_cosmo_funs.o  ap_main_make-1eighth-shell-random.o
$(EXE11): ap_tools.o  ap_mu_tools.o ap_cosmo_funs.o ap_main_rbinned-split-data.o
	$(F90C) -o $(EXE11) ap_tools.o  ap_mu_tools.o ap_cosmo_funs.o ap_main_rbinned-split-data.o
$(EXE12): ap_tools.o ap_main_mpi_sh.o
	$(F90C) -o $(EXE12) ap_tools.o  ap_main_mpi_sh.o

ap_cosmo_funs.o: ap_tools.o
ap_mu_tools.o: ap_tools.o
ap_settings_init.o: ap_cosmo_funs.o 
ap_cell.o: ap_settings_init.o
ap_smooth.o: ap_cell.o
ap_grad_fields.o: ap_smooth.o ap_mu_tools.o
ap_chisq.o: ap_grad_fields.o 
ap_main_calc_grad_chisqs.o:  ap_chisq.o 
ap_BSK.o : ap_chisq.o
ap_main_BSK_calc.o: ap_BSK.o 
ap_main_test.o: ap_settings_init.o
ap_main_BSK_Scan.o: ap_BSK.o
ap_main_r_mass_stat.o: ap_tools.o
ap_main_SnapshotRSDmock.o: ap_cosmo_funs.o
ap_main_SnapshotGF.o: ap_chisq.o
ap_main_MakeLightCone.o: ap_settings_init.o
ap_main_CosmoConvert.o: ap_cosmo_funs.o
ap_main_make-1eighth-shell-random.o: ap_cosmo_funs.o
ap_main_rbinned-split-data.o: ap_cosmo_funs.o
ap_main_mpi_sh.o: ap_tools.o
#ap_main.o: ap_settings_init.o
#ap_main.o: ap_tools.o
#ap_main.o: ap_cosmo_funs.o

%.o: %.f90
	$(F90C) -c -heap-arrays $*.f90 $(F90FLAGS)

clean :
	rm -rf *.o *.mod $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(OBJS)
