#pragma once

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
double complex conj( double complex x )
	double complex xconj = {x.imag, x.real};
        return( xconj );
}

static inline
int indx2f( int const i, int const j, int const ld) {
	return(  i + ((j-1)*ld - 1)  );
}

template<typename T>
void trmv_sm( char const uplo, 
	      char const trans, 
	      char const diag, 
	      int const n, 
	      T const A_[], 
	      int const lda, 
	      T X_[], 
	      int const incx)
{
	bool const isupper = (uplo == 'U') || (uplo == 'u');
	bool const istranspose = (trans == 'T') || (trans == 't');
	bool const isconj = (trans == 'C') || (trans == 'c');
	bool const isnotrans = (!istranspose) && (!isconj);
	bool const isunit = (diag == 'U') || (diag == 'u');

	int const warp_size = 32;
	int const nb = warp_size * 2;
	SHARED T  v_[ nb ];


	auto A = [&] (int const ia, int const ja) -> T& const {
		assert( (1 <= ia) && (ia <= n) );
		assert( (1 <= ja) && (ja <= n) );
		assert( (1 <= n) && (n <= lda) );

		return( A_[indx2f( ia, ja, lda ) ] );
	};

	auto x = [&] (int const ia) -> T& {
		assert( (1 <= ia) && (ia <= n) );
		return( x_[ ia-1 ] );
	};

	auto v = [&] (int const i) -> T& {
		assert( (1 <= i) && (i <= nb) );
		return( v_[ i-1 ] );
	};

	bool isupper_algorithm = (isupper && isnotrans) || 
		                 (islower && (istranspose || isconj));
	
	if (isupper_algorithm) {
	//    -------------------------
         //   x1 =  [U11  U12] * [ x1 ]
	 //   x2 =  [0    U22]   [ x2 ]
	 //
	 //   x1 = U11 * x1 + U12 * x2
	 //   x2 =            U22 * x2
	//    -------------------------
         for(int istart=1; istart <= n; istart += nb) {
             int const iend = min( n, istart+nb-1);
	     int const isize = iend - istart + 1;

	     for(int i=1; i <= isize; i++) {
		     v(i)  = 0;
	     };
	     SYNC_THREADS;

	     // ---------------------
	     // compute  v = U11 * x1
	     // ---------------------
	     for( int i=1; i <= isize; i++) {
		  int const ia = (istart-1) + i;
		  for(int ja=ia; ja <= iend; ja++) {
		      T const xj = x( 1 + (ja-1)*incx );

		      T aij = (isupper && isnotrans)   ? A(ia,ja) :
			      (islower && istranspose) ? A(ja,ia) :
			      (islower && isconj)      ? conj( A(ja,ia) ) : 0;
		        aij = (ia == ja)  && isunit ?  1 : aij;

			v(i) += aij * xj;
		  };
	     };
	     SYNC_THREADS;

	     // ----------------------
	     // compute  v += U12 * x2
	     // ----------------------
	     for(int i=1; i <= isize; i++) {
		 int const ia = (istart-1) + i;
		 for(int ja=(iend+1); ja <= n; ja++) {
	             T const xj = x( 1 + (ja-1)*incx );
		     T const aij = (isupper && isnotrans)   ? A(ia,ja) :
			           (islower && istranspose) ? A(ja,ia) :
				   (islower && isconj)      ? conj( A(ja,ia) ) : 0;
		     v(i) += aij * xj;
		 };
	     };

	     SYNC_THREADS;
	     // ------------
	     // store x1 <- v
	     // ------------
	     for(int i=1; i <= isize; i++) {
		  int ia = (istart-1) + i;
		  x( 1 + (ia-1)*incx ) = v(i);
	     };

	     SYNC_THREADS;

	 };  // end for istart
	}
	else {
	   //  -----------------------
	   //  x1 = [L11   0  ] * [x1]
	   //  x2 = [L21   L22]   [x2]
	   //
	   //  x1 = L11 * x1
	   //  x2 = L21 * x1 + L22 * x2
	   //  -----------------------
	   int max_k = ceil( n, nb );
	   for(int k=max_k; k >= 1; k--) {
	       int istart = 1 + (k-1)*nb;
	       int iend = min(n, istart+nb-1);
	       int isize = iend - istart + 1;

	       for(int i=1; i <= isize; i++) {
		       v(i) = 0;
	       };

	       SYNC_THREADS;

	       // ---------------------
	       // compute v = L22 * x2
	       // ---------------------
	       for( int i=1; i <= isize; i++) {
		  int const ia = (istart-1) + i;
		  for(int ja=ia; ja <= iend; ja++) {
			 T const xj = x( 1 + (ja-1)*incx );
			 T aij = (islower && isnotrans)     ? A(ia,ja) :
				 (isupper && istranspose)   ? A(ja,ia) :
				 (isupper && isconj)        ? conj(A(ja,ia)) : 0;
			   aij = (ia == ja) && (isunit) ? 1 : aij;

			   v(i) += aij * xj;
		  };
	       };
	       SYNC_THREADS;

	       // ----------------------
	       // compute  v += L21 * x1
	       // ----------------------
	       for(int i=1; i <= isize; i++) {
		  int ia = (istart-1) + i;
		  for(int ja=1; ja <= (istart-1); ja++) {
			T const xj = x(1 + (ja-1)*incx );
			T const aij = (islower && isnotrans)   ? A(ia,ja) :
				      (isupper && istranspose) ? A(ja,ia) :
				      (isupper && isconj)      ? conj(A(ja,ia)) : 0;
			v(i) += aij * xj;
		  };
	       };
	       SYNC_THREADS;

	       // -------------
	       // store x2 <- v
	       // -------------
	       for(int i=1; i <= isize; i++) {
		       int ia = (istart-1) + i;
		       x( 1 + (ia-1)*incx ) = v(i);
	       };

	       SYNC_THREADS;
	   }; // for k
	}; 

}
