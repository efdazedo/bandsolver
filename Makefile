include make.inc
FSRC=\
	banded_mod.F90  \
	bandfactor.F90  \
	bandfactor_batched.F90  \
	bandsolve.F90  \
	bandsolve_batched.F90  \
	gen_banded.F90  \
	gen_banded_batched.F90  \
	test_band_batched.F90

CSRC=\
	ztrmv_sm.hpp \

banded_mod.o: $(FSRC)
	$(FC) $(FFLAGS) -c banded_mod.F90


test1: banded_mod.o test1.F90
	$(FC) $(FFLAGS) -o test1 test1.F90  banded_mod.o $(LIBS)
