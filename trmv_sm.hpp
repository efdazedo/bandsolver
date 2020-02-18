#ifndef TRMV_SM_HPP
#define TRMV_SM_HPP 1

#include "common_sm.hpp"


//  -----------------------------------------------
//  Perform   v = op(A) * x,  op(A) is m by n but 
//  may be a trapezoidal part
//  -----------------------------------------------
template<typename T>
DEVICE_FUNCTION
void trmv_sm( char const uplo, 
	      char const trans, 
              char const diag,
              int const m,
	      int const n, 
	      T const A_[], 
	      int const lda, 
	      T x_[], 
	      T v_[])
{
#ifdef USE_GPU
        int const ix_start = 1 + threadIdx.x + 
                             threadIdx.y*blockDim.x + 
                             threadIdx.z*(blockDim.x*blockDim.y); 

        int const ix_size = blockDim.x*blockDim.y*blockDim.z;
#else
        int const ix_start = 1;
        int const ix_size = 1;
#endif
	bool const isupper = (uplo == 'U') || (uplo == 'u');
        bool const islower = (uplo == 'L') || (uplo == 'l');
	bool const istranspose = (trans == 'T') || (trans == 't');
	bool const isconj = (trans == 'C') || (trans == 'c');
	bool const isnotrans = (!istranspose) && (!isconj);
        bool const isunitdiag = (diag == 'U') || (diag == 'u');


        assert( isupper || islower );

        int const nrow = (isnotrans) ? m : n;
        int const ncol = (isnotrans) ? n : m;

	 auto A = [&] (int const ia, int const ja) -> T const & {

		return( A_[indx2f( ia, ja, lda ) ] );
	};

	 auto x = [&] (int const ia) -> T& {
		return( x_[ ia-1 ] );
	};

	 auto v = [&] (int const i) -> T& {
		return( v_[ i-1 ] );
	};


         T const one_T = make_val<T>(1.0);
         T const zero_T = make_val<T>(0.0);
	

         SYNC_THREADS;

        // ------------------------------------------------
        // check for special cases for further optimization
        // ------------------------------------------------
        bool const  is_triangular = (m == n);
        if (is_triangular && isnotrans) {
          for(int ia=ix_start; ia <= nrow; ia += ix_size) {
             T vi = zero_T;
             {
             // diagonal entry
             int const ja = ia;
             T const aij = isunitdiag ? one_T : A(ia,ja);
             T const xj = x(ja);
             // vi += aij * xj;
             fma<T>(vi,aij,xj);
             }


             int const jastart = (islower) ? 1      : (ia+1);
             int const jaend   = (islower) ? (ia-1) : ncol;

             T const *Ap = &( A(ia,jastart) );
             for(int ja=jastart; ja <= jaend; ja++) {
                     T const aij = *Ap; Ap += lda;


                      T const xj = x(ja);
                      // vi += aij * xj;
                      fma<T>(vi,aij,xj);
                 };
             v(ia) = vi;

          }; // end for ia


        }
        else {

           // -------------------
           // handle general case
           // -------------------

          for(int ia=ix_start; ia <= nrow; ia += ix_size) {
             T vi = zero_T;
             for(int ja=1; ja <= ncol; ja++) {
                  int ii = isnotrans ?  ia : ja;
                  int jj = isnotrans ?  ja : ia;
  
                  T const xj = x(ja);
                  T aij = (islower && (ii >= jj)) ||
                          (isupper && (ii <= jj)) ? A(ii,jj) : zero_T;
  
                  aij = (isunitdiag && (ii == jj)) ? one_T : aij;
                  aij = (isconj) ? conj(aij) : aij;
  
                  // vi += aij * xj;
                  fma<T>(vi,aij,xj);
             };
             v(ia) = vi;
          };

        };



        SYNC_THREADS;

}
#endif
