#include <iostream>
#include "dmalloc.hpp"
#include "bandsolve_batched_sm.hpp"


extern "C" {


void bandsolve_batched_sm( int const n,
                           int const kl_array_[],
                           int const ku_array_[],
                           zcomplex const A_[],
                           int const ldA,
                           int const old2new_[],
                           zcomplex const b_[],
                           int const ldB,
                           zcomplex       x_[],
                           int const ldX,
                           zcomplex       v_[],
                           int const ldV,
                           bool const is_full,
                           int const batchCount) 
{

#ifdef USE_GPU

   int constexpr idebug = 0;

   int constexpr warpsize = 32;
   int constexpr max_nthreads = 1024;
   int const max_klku = ldV;
   int const nwarps = iceil( max_klku, warpsize);
   int const nthreads = max(1, min(max_nthreads,nwarps * warpsize));
   // int const nthreads = 64;

   if (idebug >= 1) {
      std::cout << " max_klku =  " 
                <<   max_klku
                << " nthreads = "
                <<   nthreads
                << "\n";
   };

   dsync();
   bandsolve_batched_sm<zcomplex><<<batchCount,nthreads>>>(
                           n, kl_array_, ku_array_,
                           A_, ldA,
                           old2new_,
                           b_, ldB,
                           x_, ldX,
                           v_, ldV,
                           is_full,
                           batchCount );
   dsync();
#else

   bandsolve_batched_sm<zcomplex>(
                           n, kl_array_, ku_array_,
                           A_, ldA,
                           old2new_,
                           b_, ldB,
                           x_, ldX,
                           v_, ldV,
                           is_full,
                           batchCount );


#endif
};






}
