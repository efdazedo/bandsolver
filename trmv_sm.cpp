#include <complex>
#include "trmv_sm.hpp"




extern "C" {

void ztrmv_sm( char const uplo, 
	      char const trans, 
              int const m,
	      int const n, 
	      zcomplex const A_[], 
	      int const lda, 
	      zcomplex x_[], 
	      zcomplex v_[])
{
        trmv_sm<zcomplex>( uplo, trans, 
                           m,n,
                           A_, lda, x_, v_ );
}


void ctrmv_sm( char const uplo, 
	      char const trans, 
              int const m,
	      int const n, 
	      ccomplex const A_[], 
	      int const lda, 
	      ccomplex x_[], 
	      ccomplex v_[] )
{
        trmv_sm<ccomplex>( uplo, trans, 
                           m,n,
                           A_, lda, x_, v_ );

}


void dtrmv_sm( char const uplo, 
	      char const trans, 
              int const m,
	      int const n, 
	      double const A_[], 
	      int const lda, 
	      double x_[], 
	      double v_[] )
{
        trmv_sm<double>( uplo, trans, 
                           m,n,
                           A_, lda, x_, v_ );
}

void strmv_sm( char const uplo, 
	      char const trans, 
              int const m,
	      int const n, 
	      float const A_[], 
	      int const lda, 
	      float x_[], 
	      float v_[]) 
{
        trmv_sm<float>( uplo, trans, 
                           m,n,
                           A_, lda, x_, v_ );
}

}
