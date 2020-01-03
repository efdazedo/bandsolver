      subroutine Ztrmv_smf(uplo,trans,diag,m,n,A_in,lda,x,v)
      implicit none
      character, intent(in) :: uplo, trans, diag
      integer, intent(in) :: m,n,lda
      complex*16, intent(in) :: A_in(lda,n)
      complex*16, intent(in) :: x(*)
      complex*16, intent(in) :: v(*)

!     ---------------------------------------
!     compute v(:) = op( A(1:m, 1:n) ) * x(:)
!     ---------------------------------------
      complex*16 :: A(m,n) 
      integer :: i,j
      logical :: islower,isupper,isok,isunitdiag
      logical :: istrans,isconj,isnotrans
      integer :: mm,nn,kk,ld1,ld2,ld3
      complex*16 :: alpha,beta
      character :: transA, transB

      integer, parameter :: idebug = 0

      islower = (uplo .eq. 'L').or.(uplo .eq. 'l')
      isupper = (uplo .eq. 'U').or.(uplo .eq. 'u')
      isok = islower .or. isupper
      if (.not.isok) then
          print*,'Ztrmv_smf: invalid uplo=',uplo
          stop '** error ** '
      endif

      istrans = (trans .eq. 'T').or.(trans .eq. 't')
      isconj  = (trans .eq. 'C').or.(trans .eq. 'c')
      isnotrans = (.not. istrans).and.(.not. isconj)
      isunitdiag = (diag .eq. 'U') .or. (diag .eq. 'u')

!     -------------------------------
!     copy A_in(1:m,1:n)  into A(:,:)
!     -------------------------------
      A(:,:) = 0
      if (islower) then
         do j=1,n
         do i=j,m
           A(i,j) = A_in(i,j)
           if (isunitdiag .and. (i.eq.j)) then
               A(i,j) = 1
           endif
         enddo
         enddo
      else
         do j=1,n
         do i=1,min(j,m)
           A(i,j) = A_in(i,j)
           if (isunitdiag .and. (i.eq.j)) then
               A(i,j) = 1
           endif
         enddo
         enddo
      endif

      mm = merge(m, n, isnotrans)
      nn = 1
      kk = merge(n, m, isnotrans)
      transA = trans
      transB = 'N'
      ld1 = size(A,1)
      ld2 = kk
      ld3 = mm
      alpha = 1
      beta = 0


      call Zgemm( transA, transB, mm,nn,kk,                              &
     &      alpha, A,ld1, x, ld2,                                        &
     &      beta,  v,ld3 )

      if (idebug >= 1) then
         print*,'ztrmv_smf: uplo,trans,m,n ',                            &
     &                      uplo,trans,m,n

         do j=1,n
         do i=1,m
           print*,'A(',i,',',j,') = ',A(i,j)
         enddo
         enddo
      endif

      return
      end subroutine Ztrmv_smf
