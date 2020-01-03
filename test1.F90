
      program test1
      use banded_mod
      implicit none
      integer :: n, kl, ku, batchCount
      real(kind=wp) :: max_err, max_res

      n = 100
      kl = 11
      ku = 13
      batchCount = 10
      print*,'n, kl, ku, batchCount ', n,kl,ku,batchCount

      call test_band_batched(n,kl,ku,                                    &
     &       max_err, max_res, batchCount )

      print*,'max_err = ',max_err
      print*,'max_res = ',max_res
      stop
      end program test1

