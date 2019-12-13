      subroutine bandfactor( n, A, lda, ipiv, kl, ku, info )
      implicit none


      integer, intent(in) :: n, lda
      integer, intent(inout) :: ipiv(*)
      integer, intent(inout)  :: kl, ku, info
      complex(kind=wp), intent(inout) :: A(lda,*)


      character :: side, uplo, trans, diag
      integer :: istart,iend,isize
      logical :: isupper,islower,isok
      integer :: mm, nn
      integer :: old2new(n)
      integer :: i,j, itemp
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
       A(istart:iend,istart:iend) = Linv(1:isize,1:isize)

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

	A(istart:iend,istart:iend) = Uinv( 1:isize, 1:isize )
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

       return
       end subroutine bandfactor
