#define AB(i,j) A(kl_inout+ku_inout+1+(i)-(j),(j))

      subroutine bandfactor( n, A, lda, old2new,kl_inout,ku_inout,info,  &
     &                       is_full )
      use iso_c_binding
      implicit none

      integer, parameter :: idebug = 0

      integer, intent(in) :: n, lda
      integer, intent(inout) :: old2new(*)
      integer, intent(inout)  :: kl_inout, ku_inout, info
      complex(kind=wp), intent(inout) :: A(lda,*)
      logical(kind=c_bool), intent(in) :: is_full



      integer :: kl, ku

      character :: side, uplo, trans, diag
      integer :: istart,iend,isize
      logical :: isupper,islower,isok
      integer :: mm, nn
      integer :: ipiv(n)
      integer :: i,j, itemp, ia,ja
      integer :: ldL, ldU
      complex(kind=wp), allocatable :: Linv(:,:), Uinv(:,:)
      complex(kind=wp), allocatable :: Lmat(:,:), Umat(:,:)
      complex(kind=wp) :: alpha
      complex(kind=wp) :: Aij, Uij, Lij

      kl = kl_inout
      ku = ku_inout

      mm = n
      nn = n

      if (is_full) then
!     -------------------------------
!     matrix is in full dense storage
!     -------------------------------
        call Zgetrf( mm,nn, A,lda, ipiv, info)
        if (info.ne.0) then
            print*,'bandfactor: Zgetrf return info=',info
            return
        endif

      else
!     -----------------------------------------
!     matrix is in lapack banded storage format
!     note
!     A(kl+ku+1+i-j,j) = AB(i,j),   for max(1,j-ku) <= i <= min(m,j+kl)
!     -----------------------------------------
        isok = (ldA >= (2*kl+ku+1))
        if (.not.isok) then
            print*,'bandfactor: invalid ldA. kl,ku,ldA ',                &
     &                                       kl,ku,ldA
            stop '** error in bandfactor ** '
        endif


        call Zgbtrf_nopivot( mm,nn,kl,ku,A,ldA,ipiv,info)
        if (info.ne.0) then
                print*,'bandfactor: Zgbtrf return info=',info
                return
        endif
      endif

!     ----------
!     check ipiv
!     ----------
      do i=1,n-1
        j = ipiv(i)
        isok = (j >= i)
        if (.not.isok) then
           print*,'bandfactor: i,ipiv(i) ',i,ipiv(i)
           stop '** error in bandfactor '
        endif
      enddo



!     ---------------------------------------
!     compute lower (kl) and upper (ku) bands
!     ---------------------------------------
!
      if (is_full) then
        kl = 0;
        ku = 0;
        do j=1,n
        do i=1,n
           Aij = A(i,j)
           islower = (i-j) >= 0
           if (islower) then
             Lij = Aij
             if (Lij .ne. 0) then
               kl = max(kl, i-j)
             endif
           endif
           isupper = (j-i) >= 0
           if (isupper) then
               Uij = Aij
               if (Uij .ne. 0) then
                  ku = max(ku, j-i)
               endif
           endif
         enddo
         enddo
  
         if (idebug >= 1) then
              print*,'bandfactor:n,kl,ku', n,kl,ku
         endif

       else
!      ----------------------
!      note ku is increased
!      ----------------------
       ku = ku_inout + kl_inout 
       
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
      ldL = kl+1
      ldU = ku+1
      allocate( Lmat(ldL,kl), Linv(ldL,kl)) 
      Lmat(:,:) = 0
      Linv(:,:) = 0

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
          Linv(i,j) = merge(1,0,i.eq.j)
       enddo
       enddo

       do j=1,isize
         do i=1,(j-1)
         Lmat(i,j) = 0
         enddo

!      --------------------------
!      L matrix has unit diagonal
!      --------------------------
         i = j
         Lmat(i,j) = 1


        do i=(j+1),isize
          ia = (istart-1) + i
          ja = (istart-1) + j
          if (is_full) then
             Lmat(i,j) = A(ia,ja)
          else
             Lmat(i,j) = AB(ia,ja)
          endif
        enddo
       enddo

       if (idebug >= 2) then
          print*,'=== istart = ',istart,' ==== '
          do j=1,isize
          do i=1,isize
             Lij = Lmat(i,j)
             if (Lij .ne. 0) then
                print *, 'Lmat(',i,',',j,') = ', Lij
             endif
          enddo
          enddo
       endif


       side = 'Left'
       uplo = 'Lower'
       trans = 'NoTrans'
       diag = 'Unit'
       mm = isize
       nn = isize
       alpha = 1
       call Ztrsm( side,uplo,trans,diag, mm,nn,alpha,                                           &
     &             Lmat,ldL,Linv,ldL)

!       -------------------------------------------------------
!      copy lower triangular Linv back to A, note Linv also is unit diagonal
!       -------------------------------------------------------
        do j=1,nn
        do i=(j+1),mm
          ia = (istart-1) + i
          ja = (istart-1) + j
          if (is_full) then
            A( ia,ja) = Linv(i,j)
          else
            AB(ia,ja) = Linv(i,j)
          endif
        enddo
        enddo

       enddo
!
!
      deallocate( Lmat, Linv )

      allocate( Umat(ldU,ku), Uinv(ldU,ku) )
      Umat(:,:) = 0
      Uinv(:,:) = 0

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
          Uinv(i,j) = merge(1,0,i.eq.j)
        enddo
        enddo


        do j=1,isize
         do i=1,j
           ia = (istart-1) + i
           ja = (istart-1) + j
           if (is_full) then
              Umat(i,j) = A(ia,ja)
           else
              Umat(i,j) = AB(ia,ja)
           endif
          enddo

          do i=(j+1),isize
           Umat(i,j) = 0
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
     &              Umat,ldU,Uinv,ldU)

!       ---------------
!       copy upper triangular Uinv back to A
!       ---------------
        do j=1,nn
        do i=1,j
          ia = (istart-1) + i
          ja = (istart-1) + j
          if (is_full) then
            A(ia,ja) = Uinv(i,j)
          else
            AB(ia,ja) = Uinv(i,j)
          endif
        enddo
        enddo

       enddo

       deallocate( Umat, Uinv )
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
          do i=(j+1),min(n,kl+j)
           if (is_full) then
               Lij = A(i,j)
           else
               Lij = AB(i,j)
           endif
           if (Lij .ne. 0) then
             print *,'L(',i,',',j,') = ', Lij
           endif
          enddo
          enddo

          do j=1,n
          do i=max(1,j-ku),j
           if (is_full) then
               Uij = A(i,j)
           else
               Uij = AB(i,j)
           endif
           if (Uij .ne. 0) then
             print *,'U(',i,',',j,') = ',Uij
           endif
          enddo
          enddo

          do j=1,n
            print*,'ipiv(',j,') = ',ipiv(j)
          enddo

          do j=1,n
            print*,'old2new(',j,') = ',old2new(j)
          enddo
       endif

       if (is_full) then
         kl_inout = kl
         ku_inout = ku
       endif

       return
       end subroutine bandfactor
#undef AB
