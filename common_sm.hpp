#ifndef COMMON_SM_HPP
#define COMMON_SM_HPP 1



#include <complex>
#include <cassert>
#include <cstdlib>
#include <cmath>

#ifdef USE_GPU
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuComplex.h>

#define GLOBAL_FUNCTION  __global__ 
#define SYNC_THREADS __syncthreads()
#define SHARED_MEMORY __shared__
#define DEVICE_FUNCTION __device__
#define HOST_FUNCTION __host__

typedef cuDoubleComplex zcomplex;
typedef cuFloatComplex  ccomplex;

#else

#define GLOBAL_FUNCTION
#define SYNC_THREADS 
#define SHARED_MEMORY 
#define DEVICE_FUNCTION
#define HOST_FUNCTION

typedef std::complex<double> zcomplex;
typedef std::complex<float>  ccomplex;

#endif





static inline
HOST_FUNCTION DEVICE_FUNCTION
int indx2f( int const i, 
            int const j, 
            int const ld )
{
        return( (((i)-1) + ((j)-1)*(ld)) );
}


static inline
HOST_FUNCTION DEVICE_FUNCTION
int indx3f( int const i1, 
            int const i2, 
            int const i3,

            int const n1,
            int const n2 )
{
   return(indx2f(i1,i2,n1) + 
          ((i3)-1)*((n1)*(n2)) );
}

static inline
HOST_FUNCTION DEVICE_FUNCTION
int iceil(int const n, int const nb) {
	return(  (n + (nb-1))/nb );
}

static inline
HOST_FUNCTION DEVICE_FUNCTION
int max( int const i, int const j) {
	return(  (i >= j) ? i : j );
}

static inline
HOST_FUNCTION DEVICE_FUNCTION
int min( int const i, int const j) {
	return( (i <= j) ? i : j );
}

static inline
HOST_FUNCTION DEVICE_FUNCTION
double conj( double x ) {
	return(x);
}

static inline
HOST_FUNCTION DEVICE_FUNCTION
float conj( float x ) {
	return(x);
}


static inline
HOST_FUNCTION DEVICE_FUNCTION
zcomplex conj( zcomplex x ) {
        return( std::conj(x) );
}

static inline
HOST_FUNCTION DEVICE_FUNCTION
ccomplex conj( ccomplex x ) {
        return( std::conj(x)  );
}


#endif
