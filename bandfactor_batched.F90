      subroutine bandfactor_batched(n,A,lda,ipiv,kl_array,ku_array,                        &
     &              info_array, batchCount)
      implicit none
      integer, intent(in) :: n, lda, batchCount
      integer, intent(out) :: kl_array(batchCount)
      integer, intent(out) :: ku_array(batchCount)
      complex(kind=wp), intent(inout) :: A(lda,n,batchCount)
      integer, intent(out) :: ipiv(n,batchCount)
      integer, intent(out) :: info_array(batchCount)

      integer :: ibatch,kl,ku,info

!$omp parallel do private(ibatch,kl,ku,info)
      do ibatch=1,batchCount
	 call bandfactor(n,A(1,1,ibatch),lda,ipiv(1,ibatch),kl,ku,info)
	 kl_array(ibatch) = kl
	 ku_array(ibatch) = ku
	 info_array(ibatch) = info
      enddo

      return
      end subroutine bandfactor_batched
