      subroutine bandfactor_batched(n,A,lda,old2new,kl_array,ku_array,   &
     &              info_array, is_full, batchCount )
      use iso_c_binding
      implicit none
      integer, intent(in) :: n, lda, batchCount
      integer, intent(out) :: kl_array(batchCount)
      integer, intent(out) :: ku_array(batchCount)
      complex(kind=wp), intent(inout) :: A(lda,n,batchCount)
      integer, intent(out) :: old2new(n,batchCount)
      integer, intent(out) :: info_array(batchCount)
      logical, intent(in) :: is_full

      logical(kind=c_bool) :: is_full_c

      integer :: ibatch,kl,ku,info

      info = 0
      kl = 0
      ku = 0
      is_full_c = is_full
!$omp parallel do private(ibatch,kl,ku,info)
      do ibatch=1,batchCount
         kl = kl_array(ibatch)
         ku = ku_array(ibatch)
         call bandfactor(n,A(1,1,ibatch),lda,                            &
     &                   old2new(1,ibatch),kl,ku,info,is_full_c)
         kl_array(ibatch) = kl
         ku_array(ibatch) = ku
         info_array(ibatch) = info
      enddo

      return
      end subroutine bandfactor_batched
