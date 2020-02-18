       subroutine gen_banded_batched(n,kl,ku,A,ldA,                      &
     &               is_full,is_diag_dominant,batchCount )
       implicit none
       integer, intent(in) :: n,kl,ku,ldA,batchCount
       logical, intent(in) :: is_full
       logical, intent(in) :: is_diag_dominant
       complex(kind=wp) :: A(ldA,n,batchCount)

       integer :: ibatch

!$omp  parallel do private(ibatch)
       do ibatch=1,batchCount
         call gen_banded(n,kl,ku,A(1,1,ibatch),ldA,                      &
     &                is_full, is_diag_dominant )
       enddo

       return
       end subroutine gen_banded_batched

