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
#ifdef USE_GPU
        return( cuConj(x) );
#else
        return( std::conj( x ) );
#endif
}


static inline
HOST_FUNCTION DEVICE_FUNCTION
ccomplex conj( ccomplex x ) {
#ifdef USE_GPU
        return(  cuConjf(x) );
#else
        return( std::conj(x) );
#endif
}


//  -----------------
//  code for make_val
//  -----------------

template<typename T>
inline
HOST_FUNCTION DEVICE_FUNCTION
T make_val( double x ) {
	return(x);
}


#ifdef USE_GPU

template<>
HOST_FUNCTION DEVICE_FUNCTION
zcomplex make_val<zcomplex>( double x_re ) {
        return( make_cuDoubleComplex(x_re, 0.0) );
}

template<>
HOST_FUNCTION DEVICE_FUNCTION
ccomplex make_val<ccomplex>( double x_re ) {
        return(  make_cuFloatComplex(x_re,0.0) );
}

#endif


// ------------
// code for fma
// ------------

template<typename T>
inline
HOST_FUNCTION DEVICE_FUNCTION
void fma( T& c, T const & a, T const & b ) {
	c += a*b;
}

#ifdef USE_GPU

template<>
void fma<ccomplex>( ccomplex& c, 
                    ccomplex const & a, ccomplex const & b )  {
        c = cuCaddf( c, cuCmulf(a,b) );
}

template<>
void fma<zcomplex>( zcomplex& c, 
                    zcomplex const & a, zcomplex const & b )  {
        c = cuCadd( c, cuCmul(a,b) );
}

#endif





#endif
