include make.inc
CSRC=\
     trmv_sm.cpp
FSRC=\
	banded_mod.F90  \
	ztrmv_smf.F90 \
	bandfactor.F90  \
	bandfactor_batched.F90  \
	bandsolve.F90  \
	bandsolve_batched.F90  \
	gen_banded.F90  \
	gen_banded_batched.F90  \
	test_band_batched.F90

CSRC=\
     	bandsolve_sm.hpp \
	bandsolve_batched_sm.hpp \
	trmv_sm.hpp 

COBJ = \
       dmalloc.o \
       bandsolve_batched_sm.o  


test1: $(COBJ) banded_mod.o test1.F90  
	$(FC) $(FFLAGS) -o test1 $(COBJ) banded_mod.o test1.F90   $(LIBS)

banded_mod.o: $(FSRC)
	$(FC) $(FFLAGS) -c banded_mod.F90

dmalloc.o: $(CSRC) dmalloc.cpp
	$(CXX) $(CXXFLAGS) -c dmalloc.cpp

trmv_sm.o: $(CSRC) trmv_sm.cpp
	$(CXX) $(CXXFLAGS) -c trmv_sm.cpp

bandsolve_batched_sm.o: $(CSRC) bandsolve_batched_sm.cpp
	$(CXX) $(CXXFLAGS) -c bandsolve_batched_sm.cpp
clean:
	touch test1 banded_mod.o banded_mod.mod
	rm test1 *.o *.mod
