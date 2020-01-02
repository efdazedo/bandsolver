#ifndef TRMV_SM_HPP
#define TRMV_SM_HPP 1

#include <complex>
#include <cassert>

#ifdef USE_GPU
#define SYNC_THREADS sync_threads()
#define SHARED _shared_
#else
#define SYNC_THREADS
#define SHARED
#endif

typedef std::complex<double> zcomplex;
typedef std::complex<float>  ccomplex;

static inline
int iceil(int const n, int const nb) {
	return(  (n + (nb-1))/nb );
}

static inline
int max( int const i, int const j) {
	return(  (i >= j) ? i : j );
}

static inline
int min( int const i, int const j) {
	return( (i <= j) ? i : j );
}

static inline
double conj( double x ) {
	return(x);
}


static inline
zcomplex conj( zcomplex x ) {
        return( std::conj(x) );
}

static inline
ccomplex conj( ccomplex x ) {
        return( std::conj(x)  );
}

static inline
int indx2f( int const i, int const j, int const ld) {
	return(  (i-1) + ((j-1)*ld )  );
}

//  -----------------------------------------------
//  Perform   v = op(A) * x,  op(A) is m by n but 
//  may be a trapezoidal part
//  -----------------------------------------------
template<typename T>
void trmv_sm( char const uplo, 
	      char const trans, 
              int const m,
	      int const n, 
	      T const A_[], 
	      int const lda, 
	      T x_[], 
	      T v_[])
{
	bool const isupper = (uplo == 'U') || (uplo == 'u');
        bool const islower = (uplo == 'L') || (uplo == 'l');
	bool const istranspose = (trans == 'T') || (trans == 't');
	bool const isconj = (trans == 'C') || (trans == 'c');
	bool const isnotrans = (!istranspose) && (!isconj);


        assert( isupper || islower );

        int const nrow = (isnotrans) ? m : n;
        int const ncol = (isnotrans) ? n : m;

	auto A = [&] (int const ia, int const ja) -> T const & {

                assert( (1 <= ia) && (ia <= m) );
                assert( (1 <= ja) && (ja <= n) );
		return( A_[indx2f( ia, ja, lda ) ] );
	};

	auto x = [&] (int const ia) -> T& {
                assert( (1 <= ia) && (ia <= ncol) );
		return( x_[ ia-1 ] );
	};

	auto v = [&] (int const i) -> T& {
		assert( (1 <= i) && (i <= nrow) );
		return( v_[ i-1 ] );
	};

	

        SYNC_THREADS;

        for(int ia=1; ia <= nrow; ia++) {
           T vi = 0;
           for(int ja=1; ja <= ncol; ja++) {
                int ii = isnotrans ?  ia : ja;
                int jj = isnotrans ?  ja : ia;

                T const xj = x(jj);
                T aij = (islower && (ii >= jj)) ||
                        (isupper && (ii <= jj)) ? A(ii,jj) : 0;

                aij = (isconj) ? conj(aij) : aij;

                vi += aij * xj;
           };
           v(ia) = vi;
        };

        SYNC_THREADS;

}
#endif
