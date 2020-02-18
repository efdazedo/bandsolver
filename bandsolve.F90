      subroutine bandsolve( n, kl, ku, A, ldA, old2new, b,x,v)
!     --------------------------------------------------------------
!     Note  L, U are stored in  matrix A
!     diagonal blocks in L and U contain the explicit matrix inverse
!     --------------------------------------------------------------
      use iso_c_binding
      implicit none
      integer, intent(in) :: n, kl, ku,ldA
      integer, intent(in) :: old2new(*)
      complex(kind=wp), target, intent(in) :: A(lda,*)
      complex(kind=wp), intent(in) :: b(*)
      complex(kind=wp), intent(inout) :: x(*)
      complex(kind=wp), intent(inout) :: v(*)



      logical :: use_trmv_sm = .true.
      integer, parameter :: idebug = 1
      integer :: istat
      logical :: is_square
      integer :: inc, incx, max_k
      integer :: i1,i2,ib,i,j,k,istart,iend,isize,nn
      character :: uplo, trans, diag

      integer :: icount(1:n)
      logical :: isok

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

       if (idebug >= 1) then
!       -------------
!       check old2new
!       -------------
        icount(1:n) = 0
        do i=1,n
          j = old2new(i)
          icount(j) = icount(j) + 1
        enddo
        isok = all( icount(1:n) .eq. 1 )
        if (.not.isok) then
             print*,'bandsolve: error in old2new'
             stop '** error in bandsolve ** '
        endif
        
       endif
       x(1:n) = b( old2new(1:n) )

      if (idebug >= 2) then
         do j=1,n
           print*,'Pb(',j,') = ', x(j)
         enddo
      endif

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
        mm = nn
        incx = 1

        use_trmv_sm = .true.
        if (use_trmv_sm) then
           call Ztrmv_smf(uplo,trans,diag,mm,nn,A(istart,istart),ldA,         &
     &             x(istart), v )
           do i=1,mm
              x( (istart-1) + i) = v(i)
           enddo
        else
          call Ztrmv(uplo,trans,diag,nn,A(istart,istart),ldA,              &
     &             x(istart),incx)
        endif

!       % -------------------------------------------------------
!       % Compute using  DTRMV  triangular matrix-vector multiply
!       % -------------------------------------------------------
!       ---------------------------------------------------------------------
!       x( i1:i2) = x(i1:i2) - triu(L( i1:i2, istart:iend)) * x(istart:iend);
!       computed as 
!        v(:) = triu(L( i1:i2, istart:iend)) * x(istart:iend);
!        x( i1:i2) = x(i1:i2) - v(:)
!       ---------------------------------------------------------------------
        i1 = iend + 1;
        i2 = min(n, i1 + kl-1);
        mm = i2-i1 + 1
        nn = iend-istart+1

        if ((mm >= 1).and.(nn >= 1)) then
            uplo = 'Upper'
            trans = 'NoTrans'
            diag = 'NonUnit'

          is_square  = (nn .eq. mm)
          is_square = .false.
          if (is_square) then
!        -----------------------------------------------------
!        v(:) = triu(L( i1:i2, istart:iend)) * x(istart:iend);
!        -----------------------------------------------------
            do i=1,(iend-istart+1)
              v(i) = x( (istart-1)+i)
            enddo
            call Ztrmv( uplo, trans, diag, nn, A(i1,istart), ldA,        &
     &                v, incx)
          else
!        -----------------------------------------------------
!        v(:) = triu(L( i1:i2, istart:iend)) * x(istart:iend);
!        -----------------------------------------------------
            call Ztrmv_smf(uplo,trans,diag,mm,nn,                          &
     &            A(i1,istart),ldA, x(istart), v)

            if (idebug >= 3) then
             print*,'v(:) = triu(L(i1:i2,istart:iend))' //               &
     &               ' *x(istart:iend)'
             print*,'i1,i2,istart,iend ',i1,i2,istart,iend
             do j=istart,iend
               print*,'x(',j,') = ', x(j)
             enddo
             do j=1,(i2-i1+1)
               print*,'v(',j,') = ',v(j)
             enddo

             do j=istart,iend
             do i=i1,i2
               print*,'L(',i,',',j,') = ', A(i,j)
             enddo
             enddo
            endif

          endif

!        ---------------------------
!        x( i1:i2) = x(i1:i2) - v(:)
!        ---------------------------
          do i=1,(i2-i1+1)
               x( (i1-1)+i ) = x( (i1-1)+i) - v(i)
          enddo
        endif
       enddo

       if (idebug >= 2) then
          print*,'after solving L*x = P*b'
          do j=1,n
            print*,'x(',j,') = ', x(j)
          enddo
       endif
!
! % -------------
! % solve y = U*x
! % -------------
!
      max_k = iceil(n,ku)
      if (idebug >= 1) then
              print*,'bandsolve:n,ku,max_k', n,ku,max_k
      endif

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
       mm = nn
       incx = 1

       use_trmv_sm = .true.
       if (use_trmv_sm) then
         call Ztrmv_smf(uplo,trans,diag,mm,nn,A(istart,istart),lda,         &
     &            x(istart),v)
         do i=1,mm
           x((istart-1)+i) = v(i)
         enddo
       else
         call Ztrmv(uplo,trans,diag,nn,A(istart,istart),ldA,               &
     &            x(istart),incx)
       endif
!
!        % -------------------------------------------------------
!        % Compute using  DTRMV  triangular matrix-vector multiply
!        % -------------------------------------------------------
!        ------------------------------------------------------------------
!        x(i1:i2) = x(i1:i2) - tril(U( i1:i2, istart:iend))*x(istart:iend);
!
!        let   v(:) = tril( U(i1:i2, istart:iend)) * x(istart:iend)
!              x(i1:i2) = x(i1:i2) - v(:)
!        ------------------------------------------------------------------
       i2 = istart -1 ;
       i1 = max(1, i2 - ku + 1);
       isize = i2 - i1 + 1
       mm = i2-i1 +1
       nn = iend - istart + 1
       if ((mm >= 1).and.(nn >= 1)) then
         is_square = (mm .eq. nn)
         is_square = .false.
         uplo = 'Lower'
         trans = 'NoTrans'
         diag = 'NonUnit'
         if (is_square) then
!        -----------------------------------------------------------
!        let   v(:) = tril( U(i1:i2, istart:iend)) * x(istart:iend)
!        -----------------------------------------------------------
           do i=1,(iend-istart+1)
              v(i) = x( (istart-1) + i)
           enddo
           call Ztrmv( uplo,trans,diag,nn,                               &
     &                 A(i1,istart),ldA,v,incx)
         else
!        -----------------------------------------------------------
!        let   v(:) = tril( U(i1:i2, istart:iend)) * x(istart:iend)
!        -----------------------------------------------------------
           call Ztrmv_smf(uplo,trans,diag,mm,nn,                          &
     &             A(i1,istart),ldA,x(istart),v)

         endif

!              --------------------------
!              x(i1:i2) = x(i1:i2) - v(:)
!              --------------------------
         do i=1,(i2-i1+1)
              x( (i1-1)+i) = x( (i1-1)+i) - v(i)
         enddo
        endif
       enddo

       if (idebug >= 2) then
          print*,'after solving U*y = x '
          do j=1,n
           print*,'y(',j,') = ',x(j)
          enddo
       endif


       return
       end subroutine bandsolve
