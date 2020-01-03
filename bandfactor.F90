      subroutine bandfactor( n, A, lda, old2new, kl, ku, info )
      implicit none

      integer, parameter :: idebug = 0

      integer, intent(in) :: n, lda
      integer, intent(inout) :: old2new(*)
      integer, intent(inout)  :: kl, ku, info
      complex(kind=wp), intent(inout) :: A(lda,*)


      character :: side, uplo, trans, diag
      integer :: istart,iend,isize
      logical :: isupper,islower,isok
      integer :: mm, nn
      integer :: ipiv(n)
      integer :: i,j, itemp, ia,ja
      integer :: ldL, ldU
      complex(kind=wp), allocatable :: Linv(:,:), Uinv(:,:)
      complex(kind=wp) :: alpha
      complex(kind=wp) :: Uij, Lij


      mm = n
      nn = n


      call Zgetrf( mm,nn, A,lda, ipiv, info)
      if (info.ne.0) then
          print*,'bandfactor: Zgetrf return info=',info
          return
      endif

!     ---------------------------------------
!     compute lower (kl) and upper (ku) bands
!     ---------------------------------------
!
      kl = 0;
      ku = 0;
      do j=1,n
      do i=1,n
         islower = (i-j) >= 0
         if (islower) then
           Lij = A(i,j)
           if (Lij .ne. 0) then
             kl = max(kl, i-j)
           endif
         endif
         isupper = (j-i) >= 0
         if (isupper) then
             Uij = A(i,j)
             if (Uij .ne. 0) then
                ku = max(ku, j-i)
             endif
         endif
       enddo
       enddo

       if (idebug >= 1) then
            print*,'bandfactor:n,kl,ku', n,kl,ku
       endif
!
!% ------------------------------
!% precompute the explicit inverse for
!% triangular blocks in L and U
!% ------------------------------
!
!% ------------------------
!% note L has unit diagonal
!% ------------------------
      ldL = kl
      ldU = ku
      allocate( Linv(ldL,kl), Uinv(ldU,ku) )
      do istart=1,n,kl
       iend = min( n, istart+kl-1);
       isize = iend - istart + 1;
!
!    % -----------------------------------------------------------------
!    % Compute  using  DTRSM
!    % L( istart:iend, istart:iend) = inv( L(istart:iend,istart:iend) );
!    % -----------------------------------------------------------------
!    --------------------------------------------------------------------------------
!    L( istart:iend, istart:iend) =  L(istart:iend,istart:iend) \ eye( isize, isize);
!    --------------------------------------------------------------------------------
       do j=1,isize
       do i=1,isize
          Linv(i,j) = merge( 1, 0, (i.eq.j) )
       enddo
       enddo

       side = 'Left'
       uplo = 'Lower'
       trans = 'NoTrans'
       diag = 'Unit'
       mm = isize
       nn = isize
       alpha = 1
       call Ztrsm( side,uplo,trans,diag, mm,nn,alpha,                                           &
     &             A(istart,istart),ldA,Linv,ldL)

!       -------------------------------------------------------
!      copy lower triangular Linv back to A, note Linv also is unit diagonal
!       -------------------------------------------------------
        do j=1,nn
        do i=(j+1),mm
          ia = (istart-1) + i
          ja = (istart-1) + j
          A( ia,ja) = Linv(i,j)
        enddo
        enddo

       enddo
!
!
      do istart=1,n,ku
       iend = min(n, istart+ku-1)
       isize = iend - istart + 1

!    % -----------------------------------------------------------------
!    % Compute using DTRSM
!    % U( istart:iend, istart:iend ) = inv( U(istart:iend, istart:iend) );
!    % -----------------------------------------------------------------
!    -------------------------------------------------------------------------------
!    U( istart:iend, istart:iend ) = U(istart:iend, istart:iend) \ eye(isize,isize);
!    -------------------------------------------------------------------------------
        do j=1,isize
        do i=1,isize
          Uinv(i,j) = merge( 1, 0, i.eq.j)
        enddo
        enddo


        side = 'Left'
        uplo = 'Uppper'
        trans = 'NoTrans'
        diag = 'NonUnit'
        mm = isize
        nn = isize
        alpha = 1
        call Ztrsm( side, uplo, trans, diag, mm,nn,alpha,                                   &
     &              A(istart,istart),ldA,Uinv,ldU)

!       ---------------
!       copy upper triangular Uinv back to A
!       ---------------
        do j=1,nn
        do i=1,j
          ia = (istart-1) + i
          ja = (istart-1) + j
          A(ia,ja) = Uinv(i,j)
        enddo
        enddo

       enddo

       deallocate( Linv, Uinv )
!
!
!% ---------------------------
!% generate permutation vector
!% ---------------------------
! ------------------------
! ip = reshape( 1:n, n,1);
! old2new(1:n) = P*ip(1:n);
! ------------------------
       do i=1,n
         old2new(i) = i
       enddo
       do i=1,n-1
         j = ipiv(i)
         if (j.ne.i) then
               ! ------------
               ! perform swap
               ! ------------
               itemp = old2new(i)
               old2new(i) = old2new(j)
               old2new(j) = itemp
         endif
       enddo

       if (idebug >= 2) then
          do j=1,n
          do i=(j+1),n
           print *,'L(',i,',',j,') = ', A(i,j)
          enddo
          enddo

          do j=1,n
          do i=1,j
           print *,'U(',i,',',j,') = ',A(i,j)
          enddo
          enddo

          do j=1,n
            print*,'old2new(',j,') = ',old2new(j)
          enddo
       endif

       return
       end subroutine bandfactor
