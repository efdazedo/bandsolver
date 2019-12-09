      subroutine bandsolve( n, kl, ku, A, ldA, old2new, b,inc  )
!     --------------------------------------------------------------
!     Note  L, U are stored in  matrix A
!     diagonal blocks in L and U contain the explicit matrix inverse
!     --------------------------------------------------------------
      implicit none
      integer, intent(in) :: n, kl, ku,ldA,inc
      integer, intent(in) :: old2new(*)
      complex*16, intent(in) :: A(lda,*)
      complex*16, intent(inout) :: b(*)


      interface
	subroutine Ztrmv( uplo, trans, diag, n, A, ldA, X, incx )
	implicit none
	character uplo, trans, diag
	integer n, ldA, incx
	complex*16 A(ldA,*), X(*)
	end subroutine Ztrmv
      end interface

      integer, parameter :: idebug = 1
      logical :: isok
      complex*16 :: X(n)
      integer :: incx, max_k
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
! x(1:n) = b(old2new(1:n));

      do i=1,n
	j = old2new(i)
	ib = 1 + (j-1)*inc
	x(i) = b(ib)
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
        call Ztrmv(uplo,trans,diag,nn,A(istart,istart),ldA,                                   &
     &             x(istart),incx)

        i1 = iend + 1;
        i2 = min(n, i1 + kl-1);
        if (isize >= 1) then
!       % -------------------------------------------------------
!       % Compute using  DTRMV  triangular matrix-vector multiply
!       % -------------------------------------------------------
!       ---------------------------------------------------------------------
!       x( i1:i2) = x(i1:i2) - triu(L( i1:i2, istart:iend)) * x(istart:iend);
!       ---------------------------------------------------------------------

          nn = i2-i1 + 1
	  isok  = (iend-istart+1) .eq. (i2-i1+1)
	  if (.not.isok) then
           print*,'bandsolve:non-triangular matrix: i1,i2,istart,iend ',                    &
     &               i1,i2,istart,iend
            stop '** assertion error ** '
          endif
	  uplo = 'Upper'
	  trans = 'NoTrans'
	  diag = 'NonUnit'
          call Ztrmv( uplo, trans, diag, nn, A(i1,istart), ldA,                               &
     &                x(istart), incx)
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
       call Ztrmv(uplo,trans,diag,nn,A(istart,istart),ldA,                                    &
     &            x(istart),incx)
!
       i2 = istart -1 ;
       i1 = max(1, i2 - ku + 1);
       isize = i2 - i1 + 1
       if (isize >= 1) then
!        % -------------------------------------------------------
!        % Compute using  DTRMV  triangular matrix-vector multiply
!        % -------------------------------------------------------
!        ------------------------------------------------------------------
!        x(i1:i2) = x(i1:i2) - tril(U( i1:i2, istart:iend))*x(istart:iend);
!        ------------------------------------------------------------------
         uplo = 'Lower'
	 trans = 'NoTrans'
	 diag = 'NonUnit'
	 nn = i2-i1 +1
         call Ztrmv( uplo,trans,diag,nn,A(i1,istart),ldA,x(istart),incx)
	endif
       enddo

!      -------------
!      copy solution 
!      -------------
       do i=1,n
	 ib = 1 + (i-1)*inc
	 b(ib) = x(i)
       enddo

       return
       end subroutine bandsolve
