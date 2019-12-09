       subroutine gen_banded_batched(n,kl,ku,A,ldA, batchCount )
       implicit none
       integer, intent(in) :: n,kl,ku,ldA,batchCount
       complex(kind=wp) :: A(ldA,n,batchCount)

       integer :: ibatch

!$omp  parallel do private(ibatch)
       do ibatch=1,batchCount
	  call gen_banded(n,kl,ku,A(1,1,ibatch),ldA )
       enddo

       return
       end subroutine gen_banded_batched

