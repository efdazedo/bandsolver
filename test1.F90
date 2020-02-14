
      program test1
      use banded_mod
      implicit none
      integer :: n, kl, ku, batchCount
      integer :: icase
      logical :: is_full
      real(kind=wp) :: max_err, max_res
      integer, parameter :: idebug = 1

      if (idebug >= 1) then
           n = 100
           kl = 11
           ku = 13
           batchCount = 10
      else if (idebug >= 1) then
           n = 4
           kl = 1
           ku = 1
           batchCount = 1
      endif

      if (idebug >= 1) then
      print*,'n, kl, ku, batchCount ', n,kl,ku,batchCount

      do icase=0,1
        is_full = (icase.ne.0)
        call test_band_batched(n,kl,ku,                                    &
     &       max_err, max_res, is_full, batchCount )

        print*,'max_err = ',max_err
        print*,'max_res = ',max_res
      enddo
      endif


!     ----------------
!     performance test
!     ----------------
      n = 1024
      batchCount = 256
      kl = 32
      ku = 32

      print*,'n, kl, ku, batchCount ', n,kl,ku,batchCount

      do icase=0,1
        is_full = (icase.ne.0)
        print*,'is_full = ', is_full
        call test_band_batched(n,kl,ku,                                   &
     &       max_err, max_res, is_full, batchCount )

        print*,'max_err = ',max_err
      enddo

      stop
      end program test1

