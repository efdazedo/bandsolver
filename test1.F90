
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

!     ----------------
!     performance test
!     ----------------
      n = 512
      batchCount = 256
      kl = 16
      ku = 17

      print*,'n, kl, ku, batchCount ', n,kl,ku,batchCount

      call test_band_batched(n,kl,ku,                                    &
     &       max_err, max_res, batchCount )

      print*,'max_err = ',max_err

      stop
      end program test1

