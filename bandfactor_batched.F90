      subroutine bandfactor_batched(n,A,lda,old2new,kl_array,ku_array,                        &
     &              info_array, batchCount)
      implicit none
      integer, intent(in) :: n, lda, batchCount
      integer, intent(out) :: kl_array(batchCount)
      integer, intent(out) :: ku_array(batchCount)
      complex(kind=wp), intent(inout) :: A(lda,n,batchCount)
      integer, intent(out) :: old2new(n,batchCount)
      integer, intent(out) :: info_array(batchCount)

      integer :: ibatch,kl,ku,info

      info = 0
      kl = 0
      ku = 0
!$omp parallel do private(ibatch,kl,ku,info)
      do ibatch=1,batchCount
         call bandfactor(n,A(1,1,ibatch),lda,old2new(1,ibatch),kl,ku,info)
         kl_array(ibatch) = kl
         ku_array(ibatch) = ku
         info_array(ibatch) = info
      enddo

      return
      end subroutine bandfactor_batched
