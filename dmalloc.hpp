#ifndef DMALLOC_HPP
#define DMALLOC_HPP 1

#include <cstddef>

extern "C" {
        void dsync();
        void *dmalloc( size_t const nbytes );
        void dfree( void *ptr );
        void host2acc( void *dest, void const * const src, size_t const nbytes );
        void acc2host( void *dest, void const * const src, size_t const nbytes );
}





#endif

