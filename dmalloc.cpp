#ifdef USE_GPU
#include <cuda.h>
#endif

#include <cassert>
#include <string.h>

extern "C" {
        void *dmalloc( size_t nbytes ) {
                void *ptr = 0;
#ifdef USE_GPU
                cudaError_t istat = cudaMallocManaged( &ptr, nbytes );
                assert( istat == cudaSuccess ); 

#else
                ptr = (void *) malloc( nbytes );
#endif
                assert(ptr != 0);
                return(ptr);
        }


        void dfree( void *ptr) {
                if (ptr != 0) {
#ifdef USE_GPU
                 cudaError_t istat = cudaFree( ptr );
                 assert( istat == cudaSuccess );
#else
                   free(ptr);
#endif
                }
        }


        void host2acc( void *dest, void *src, size_t nbytes ) {
#ifdef USE_GPU
                cudaError_t istat = cudaMemcpy( dest, src, nbytes, cudaMemcpyHostToDevice );
                assert( istat == cudaSuccess );
#else
                memcpy( dest, src, nbytes );
#endif
        }


        void acc2host( void *dest, void *src, size_t nbytes ) {
#ifdef USE_GPU
                cudaError_t istat = cudaMemcpy( dest, src, nbytes, cudaMemcpyDeviceToHost );
                assert( istat == cudaSuccess );

#else
                memcpy( dest, src, nbytes );
#endif
        }


}

