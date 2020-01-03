#ifndef BANDSOLVER_BATCHED_SM_HPP
#define BANDSOLVER_BATCHED_SM_HPP 1

#include "bandsolve_sm.hpp"

template<typename T>
GLOBAL_FUNCTION
void bandsolve_batched_sm( int const n,
                           int const kl_array_[],
                           int const ku_array_[],
                           T const A_[],
                           int const ldA,
                           int const old2new_[],
                           T const b_[],
                           int const ldB,
                           T       x_[],
                           int const ldX,
                           T       v_[],
                           int const ldV,
                           int const batchCount) 
{

        /*
         * -----------------------------------------------------
         * use 1-based indexing to match Fortran and matlab code
         * -----------------------------------------------------
         */
        auto v = [&] (int const i, int const ibatch) -> T& {
                return( v_[indx2f(i,ibatch,ldV)] );
        };

        auto x = [&] (int const i, int const ibatch) -> T& {
                return( x_[indx2f(i,ibatch,ldX)] );
        };

        auto old2new = [&] (int const i, int const ibatch) -> int const & {
                return( old2new_[indx2f(i,ibatch,n)] );
        };

        auto b = [&] (int const i, int const ibatch) -> T const & {
                return( b_[indx2f(i,ibatch,ldB)] );
        };

        auto A = [&] (int const i, int const j, int const ibatch) -> T const & {
                return( A_[ indx3f(i,j,ibatch,ldA,n) ] );
        };



        auto kl_array = [&] (int const ibatch) -> int const & {
                return( kl_array_[ ibatch-1 ] );
        };

        auto ku_array = [&] (int const ibatch) -> int const & {
                return( ku_array_[ ibatch-1 ] );
        };


#ifdef USE_GPU
        int const iz_start = blockIdx.x + 1;
        int const iz_size =  gridDim.x;

        assert( gridDim.y == 1);
        assert( gridDim.z == 1);
#else
        int const iz_start = 1;
        int const iz_size = 1;
#endif

        for(int ibatch=iz_start; ibatch <= batchCount; ibatch += iz_size) {
                int const kl = kl_array(ibatch);
                int const ku = ku_array(ibatch);
                bandsolve_sm<T>( n, kl, ku, 
                                 &(A(1,1,ibatch)), ldA, 
                                 &(old2new(1,ibatch)), 
                                 &(b(1,ibatch)),
                                 &(x(1,ibatch)),
                                 &(v(1,ibatch)) );
        };





}

#endif
