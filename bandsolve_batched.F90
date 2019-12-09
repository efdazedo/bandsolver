      subroutine bandsolve_batched(n, kl_array,ku_array,A,ldA,                             &
     &                  old2new,b,inc,batchCount)
      implicit none
      integer, intent(in) :: n,ldA,batchCount
      integer, intent(in) :: kl_array(batchCount)
      integer, intent(in) :: ku_array(batchCount)
      complex(kind=wp), intent(in) :: A(ldA,n,batchCount)
      integer, intent(in) :: old2new(n,batchCount)
      integer, intent(in) :: inc
      complex(kind=wp), intent(in) :: b(inc*n,batchCount)

      integer :: ibatch,kl,ku

!$omp parallel do private(ibatch,kl,ku)
      do ibatch=1,batchCount
	 kl = kl_array(ibatch)
	 ku = ku_array(ibatch)
	 call bandsolve(n,kl,ku,A(1,1,ibatch),ldA, old2new(1,ibatch),                       &
     &                  B(1,ibatch), inc )
      enddo
      return
      end subroutine bandsolve_batched
