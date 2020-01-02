      subroutine bandsolve( n, kl, ku, A, ldA, old2new, b)
!     --------------------------------------------------------------
!     Note  L, U are stored in  matrix A
!     diagonal blocks in L and U contain the explicit matrix inverse
!     --------------------------------------------------------------
      use iso_c_binding
      implicit none
      integer, intent(in) :: n, kl, ku,ldA
      integer, intent(in) :: old2new(*)
      complex(kind=wp), target, intent(in) :: A(lda,*)
      complex(kind=wp), intent(inout) :: b(*)



      integer, parameter :: idebug = 1
      integer :: istat
      logical :: is_square
      complex(kind=wp), target  :: x(n), v(max(kl,ku))
      integer :: inc, incx, max_k
      integer :: i1,i2,ib,i,j,k,istart,iend,isize,nn
      character :: uplo, trans, diag

!     ---------------
!     inline function
!     ---------------
      integer :: iceil,mm,kk
      iceil(mm,kk) = int( (mm + (kk-1))/kk )
! 
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

      do i=1,n
       j = old2new(i)
       x(j) = b(i)
      enddo


! % --------------
! % solve L*y = P*b
! % --------------
      do istart=1,n,kl
	 iend = min(n, istart+kl-1)
	 isize = iend - istart + 1

!
!     % -------------------------------------------------------
!     % Compute using  DTRMV  triangular matrix-vector multiply
!     % -------------------------------------------------------
!     ----------------------------------------------------------------------
!     x(istart:iend) =  tril(L( istart:iend, istart:iend) )* x(istart:iend);
!     ----------------------------------------------------------------------

	uplo = 'Lower'
	trans = 'NoTrans'
	diag = 'Unit'
	nn = isize
	incx = 1
        call Ztrmv(uplo,trans,diag,nn,A(istart,istart),ldA,              &
     &             x(istart),incx)

!       % -------------------------------------------------------
!       % Compute using  DTRMV  triangular matrix-vector multiply
!       % -------------------------------------------------------
!       ---------------------------------------------------------------------
!       x( i1:i2) = x(i1:i2) - triu(L( i1:i2, istart:iend)) * x(istart:iend);
!       ---------------------------------------------------------------------
        i1 = iend + 1;
        i2 = min(n, i1 + kl-1);
        mm = i2-i1 + 1
        nn = iend-istart+1

        if ((mm >= 1).and.(nn >= 1)) then

	  is_square  = (nn .eq. mm)
          if (is_square) then
	    uplo = 'Upper'
	    trans = 'NoTrans'
	    diag = 'NonUnit'
            do i=1,(iend-istart+1)
              v(i) = x( (istart-1)+i)
            enddo
            call Ztrmv( uplo, trans, diag, nn, A(i1,istart), ldA,        &
     &                v, incx)
          else
            call Ztrmv_sm(uplo,trans,mm,nn,                         &
!     &            c_loc(A(i1,istart)),ldA, c_loc(x(istart)), c_loc(v))
     &            A(i1,istart),ldA, x(istart), v)
          endif

          do i=1,(i2-i1+1)
               x( (i1-1)+i ) = x( (i1-1)+i) - v(i)
          enddo
	endif
       enddo
!
! % -------------
! % solve y = U*x
! % -------------
!
      max_k = iceil(n,ku)
      do k=max_k,1,-1
      istart = 1 + (k-1)*ku;
      iend = min(n, istart+ku-1);
      isize = iend - istart + 1

!     % -------------------------------------------------------
!     % Compute using  DTRMV  triangular matrix-vector multiply
!     % -------------------------------------------------------
!     ------------------------------------------------------------------------
!     x( istart:iend ) = triu(U( istart:iend, istart:iend) ) * x( istart:iend);
!     ------------------------------------------------------------------------
       uplo = 'Upper'
       trans = 'NoTrans'
       diag = 'NonUnit'
       nn = isize
       incx = 1
       call Ztrmv(uplo,trans,diag,nn,A(istart,istart),ldA,               &
     &            x(istart),incx)
!
!        % -------------------------------------------------------
!        % Compute using  DTRMV  triangular matrix-vector multiply
!        % -------------------------------------------------------
!        ------------------------------------------------------------------
!        x(i1:i2) = x(i1:i2) - tril(U( i1:i2, istart:iend))*x(istart:iend);
!        ------------------------------------------------------------------
       i2 = istart -1 ;
       i1 = max(1, i2 - ku + 1);
       isize = i2 - i1 + 1
       mm = i2-i1 +1
       nn = iend - istart + 1
       if ((mm >= 1).and.(nn >= 1)) then
         is_square = (mm .eq. nn)
         uplo = 'Lower'
	 trans = 'NoTrans'
	 diag = 'NonUnit'
         if (is_square) then
           do i=1,(iend-istart+1)
              v(i) = x( (istart-1) + i)
           enddo
           call Ztrmv( uplo,trans,diag,nn,                               &
     &                 A(i1,istart),ldA,v,incx)
         else
           call Ztrmv_sm(uplo,trans,mm,nn,                               &
!     &             c_loc(A(i1,istart)),ldA,c_loc(x(istart)),c_loc(v))
     &             A(i1,istart),ldA,x(i1),v)
         endif

         do i=1,(i2-i1+1)
              x( (i1-1)+i) = x( (i1-1)+i) - v(i)
         enddo
	endif
       enddo

!      -------------
!      copy solution 
!      -------------
       do i=1,n
          b(i) = x(i)
       enddo


       return
       end subroutine bandsolve
