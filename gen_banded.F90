      subroutine gen_banded(n,kl,ku, A, ldA )
      implicit none
! % ----------------------
! % generate banded matrix of size n by n
! % with lower bandwidth kl and upper bandwidth ku
! % ----------------------
      integer, intent(in) :: n, kl, ku, ldA
      complex(kind=wp), intent(inout) :: A(ldA,n)

      integer :: i,j
      real(kind=wp) :: x_re(n), x_im(n)
      complex(kind=wp) :: aij
      integer, parameter :: idebug = 1

      do j=1,n
	 call random_number( x_re(1:n) )
	 call random_number( x_im(1:n) )
	 x_re(1:n) = 2*x_re(1:n) - 1
	 x_im(1:n) = 2*x_im(1:n) - 1
	 do i=1,n
	   aij = cmplx( x_re(i), x_im(i), kind=wp)
           if (i - j > kl) then
	       aij = 0
	   endif
	   if (j - i > ku) then
	       aij = 0
	   endif
	   A(i,j) = aij
	 enddo
      enddo

      if (idebug >= 1) then
!      -----------------------------
!      make diagonal to be very large
!      to avoid pivoting
!      -----------------------------
       do j=1,n
         i = j
         A(i,j) = 2*n*n
       enddo
     endif

      return
      end subroutine gen_banded
