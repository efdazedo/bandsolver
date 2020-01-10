FC=gfortran
CXX=g++

FFLAGS_g=-g  -fbounds-check
CXXFLAGS_g=-g  -std=c++11


FFLAGS_O=-O3   -fopenmp
CXXFLAGS_O=-O3 -fopenmp  -std=c++11

CXXFLAGS=$(CXXFLAGS_O)
FFLAGS=$(FFLAGS_O)
LIBS= -lstdc++ -llapack -lblas
