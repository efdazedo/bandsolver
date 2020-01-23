
#ifndef BANDSOLVE_SM_HPP
#define BANDSOLVE_SM_HPP 1

#include "common_sm.hpp"

#include "trmv_sm.hpp"

//  -----------------------------------------------
//  solve banded linear system
//  -----------------------------------------------
template<typename T>
DEVICE_FUNCTION
void bandsolve_sm( int const n,
                int const kl,
                int const ku_in,
                T const A_[],
                int const ldA,
                int const old2new_[],
                T const b_[],
                T       x_[],
                T       v_[],
                bool const is_full )
               
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
        int const ldAB = ldA - 1;

        int const ku = (is_full) ? ku_in : ku_in + kl;


        auto x = [&] (int const i) -> T& {
                return( x_[ (i)-1 ] );
        };

        auto AB =[&] (int const i, int const j) -> T const & {
                bool const is_in_band = (max(1,j-ku) <= i) && (i <= min(n,j+kl));
                assert( is_in_band );

                return(   A_[ indx2f( kl+ku_in+1+i-j,j, ldA)] );
                                        
        };

        auto A = [&] (int const ia, int const ja) -> T const & {
                return( A_[indx2f(ia,ja,ldA)] );
        };

        auto v = [&] (int const i) -> T& {
                return( v_[ (i)-1 ] );
        };

        auto old2new = [&] (int const i) -> int const & {
                return( old2new_[ (i)-1 ] );
        };

        auto b = [&] (int const i) -> T const & {
                return( b_[ (i)-1 ] );
        };
 
/*
! % --------------------
! % solve  A*x = b,
! % P*A = L * U
! % P*A*x = P*b
! % L*U*x = (P*b)
! % L*y = (P*b), y = U*x
! % --------------------
! 
! % -------------------
! % perform permutation
! % -------------------
       x(1:n) = b( old2new(1:n) )
*/

        SYNC_THREADS;

        for(int i=ix_start; i <= n; i += ix_size) {
                x(i) = b( old2new(i) );
        };

        SYNC_THREADS;

/*
% --------------
% solve L*y = P*b
% --------------
*/
  for(int istart=1; istart <= n; istart += kl) {
        int const iend = min(n, istart+kl-1);
        int const isize = iend - istart + 1;

        /*
          ---------------------------------------------------------------------
          x(istart:iend) =  tril(L( istart:iend, istart:iend) )* x(istart:iend);
          ---------------------------------------------------------------------
        */
        SYNC_THREADS;

        {

         int const nn = isize;
         int const mm = nn;

         char const uplo = 'L';
         char const trans = 'N';
         char const diag = 'U';
         T const * const pA = (is_full) ? &(A(istart,istart)) : 
                                          &(AB(istart,istart));
         int const ld1 = (is_full) ? ldA  : ldAB;

         trmv_sm<T>(uplo,trans,diag,mm,nn,
                    pA, ld1, &(x(istart)),&(v(1)) );
        };


        SYNC_THREADS;

         for(int i=ix_start; i <= isize; i += ix_size) {
                 x( (istart-1)+i) = v(i);
         };




        int const i1 = iend + 1;
        int const i2 = min(n,i1+kl-1);
        int const mm = i2-i1+1;
        int const nn = isize;

/*
!       ---------------------------------------------------------------------
!       x( i1:i2) = x(i1:i2) - triu(L( i1:i2, istart:iend)) * x(istart:iend);
!       computed as 
!        v(:) = triu(L( i1:i2, istart:iend)) * x(istart:iend);
!        x( i1:i2) = x(i1:i2) - v(:)
!       ---------------------------------------------------------------------
*/
        SYNC_THREADS;

        if ((mm >= 1) && (nn >= 1)) {

            char const uplo = 'U';
            char const trans = 'N';
            char const diag = 'N';

            T const * const pA = (is_full) ? &(A(i1,istart)) :
                                             &(AB(i1,istart));
            int const ld1 = (is_full) ? ldA : ldAB;

            trmv_sm<T>(uplo,trans,diag,mm,nn,
                       pA,ld1, &(x(istart)), &(v(1)) );

        };

           SYNC_THREADS;

           T const neg_one_T = make_val<T>(-1.0);

           for(int i=ix_start; i <= mm; i += ix_size) {
                // x( (i1-1)+i) = x( (i1-1)+i) - v(i);
                fma<T>(x( (i1-1)+i), neg_one_T, v(i) );
           };

           SYNC_THREADS;
  }; //  end for istart



/*
! % -------------
! % solve y = U*x
! % -------------
*/
  int const max_k = iceil(n,ku);

  for(int k=max_k; k >= 1; k--) {
          int const istart = 1 + (k-1)*ku;
          int const iend = min(n, istart+ku-1);
          int const isize = iend - istart + 1;

/*
!     ------------------------------------------------------------------------
!     x( istart:iend ) = triu(U( istart:iend, istart:iend) ) * x( istart:iend);
!     ------------------------------------------------------------------------
*/
        SYNC_THREADS;
        {


          int const nn = isize;
          int const mm = nn;
          char const uplo = 'U';
          char const trans = 'N';
          char const diag = 'N';

          T const * const pA = (is_full) ? &(A(istart,istart)) :
                                           &(AB(istart,istart));
          int const ld1 = (is_full) ? ldA : ldAB;

          trmv_sm<T>(uplo,trans,diag,mm,nn,
                     pA, ld1, &(x(istart)), &(v(1)) );
         };


        SYNC_THREADS;

          for(int i=ix_start; i <= isize; i += ix_size) {
                  x( (istart-1)+i) = v(i);
          };

        SYNC_THREADS;

/*
!        ------------------------------------------------------------------
!        x(i1:i2) = x(i1:i2) - tril(U( i1:i2, istart:iend))*x(istart:iend);
!
!        let   v(:) = tril( U(i1:i2, istart:iend)) * x(istart:iend)
!              x(i1:i2) = x(i1:i2) - v(:)
!        ------------------------------------------------------------------
*/
         int const i2 = istart-1;
         int const i1 = max(1, i2-ku+1);
         int const mm = i2-i1+1;
         int const nn = isize;

         if ((mm >= 1) && (nn >= 1)) {
          char const uplo = 'L';
          char const trans = 'N';
          char const diag = 'N';

          int const ld1 = (is_full)? ldA : ldAB;
          T const * const pA = is_full ? &(A(i1,istart)) :
                                         &(AB(i1,istart));

          trmv_sm<T>(uplo,trans,diag,mm,nn,
                     pA,ld1,&(x(istart)), &(v(1)) );
         }

        SYNC_THREADS;

/*
!              --------------------------
!              x(i1:i2) = x(i1:i2) - v(:)
!              --------------------------
*/
        T const neg_one_T = make_val<T>(-1.0);

        for(int i=ix_start; i <= (i2-i1+1); i += ix_size) {
                // x( (i1-1)+i) = x( (i1-1)+i) - v(i);
                fma<T>( x( (i1-1)+i), neg_one_T, v(i) );
        };

        SYNC_THREADS;
     }; // end for istart

        SYNC_THREADS;

        return;
}
#endif
