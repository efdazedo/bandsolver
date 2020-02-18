FC=gfortran
CXX=nvcc -x cu
FFLAGS=-O3 -fopenmp -DUSE_GPU -DUSE_DMALLOC

CUDADIR=/usr/local/cuda-10.1

CXXFLAGS=\
        -DUSE_GPU \
        -O3 -Xcompiler -fopenmp   -std=c++11 



# -O2 -Munroll=c:1 -Mlre -Mautoinline -Mvect=simd -Mcache_align   

# CXXFLAGS=-g     -mp   -std=c++11 -Minfo=all
LIBS= -llapack -lblas  -L$(CUDADIR)/lib64 -lcuda  -lcudart -lstdc++
