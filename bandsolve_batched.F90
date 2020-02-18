      subroutine bandsolve_batched(n, kl_array,ku_array,A,ldA,           &
     &                  old2new,B,ldB,X,ldX,v,ldV,batchCount)
!     -------------------------------------------------------
!     B(1:n,1:batchCount) is the rhs vector
!     X(1:n,1:batchCount) is the solution vector
!     v(ldV,1:batchCount) is temporary work space
!     ldV = max( maxval( kl_array(:)), maxval( ku_array(:)) )
!     -------------------------------------------------------
      implicit none
      integer, intent(in) :: n,ldA,batchCount
      integer, intent(in) :: kl_array(batchCount)
      integer, intent(in) :: ku_array(batchCount)
      complex(kind=wp), intent(in) :: A(ldA,n,batchCount)
      integer, intent(in) :: old2new(n,batchCount)
      integer, intent(in) :: ldB,ldX,ldV
      complex(kind=wp), intent(in) :: B(ldB,batchCount)
      complex(kind=wp), intent(inout) :: X(ldX,batchCount)
      complex(kind=wp), intent(inout) :: v(ldV,batchCount)
      logical :: isok_kl, isok_ku, isok

      integer :: ibatch,kl,ku

!$omp parallel do private(ibatch,kl,ku)
      do ibatch=1,batchCount
         kl = kl_array(ibatch)
         ku = ku_array(ibatch)

         isok_kl = (1 <= kl) .and. (kl <= n) 
         isok_ku = (1 <= ku) .and. (ku <= n)
         isok = isok_kl .and. isok_ku
         if (.not. isok) then
                 print*,'ibatch,kl,ku',ibatch,kl,ku
                 stop '** error in bandsolve_batched ** '
         endif

         call bandsolve(n,kl,ku,A(1,1,ibatch),ldA, old2new(1,ibatch),    &
     &                  B(1,ibatch), X(1,ibatch), v(1,ibatch) )
      enddo
      return
      end subroutine bandsolve_batched
