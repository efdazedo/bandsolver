#define AB(i,j) A(kl+ku+1+(i)-(j),(j))
      subroutine gen_banded(n,kl,ku, A, ldA, is_full )
      implicit none
! % ----------------------
! % generate banded matrix of size n by n
! % with lower bandwidth kl and upper bandwidth ku
! % ----------------------
      integer, intent(in) :: n, kl, ku, ldA
      complex(kind=wp), intent(inout) :: A(ldA,n)
      logical, intent(in) :: is_full

      integer :: i,j
      real(kind=wp) :: x_re(n), x_im(n)
      complex(kind=wp) :: aij
      integer, parameter :: idebug = 0
      logical :: isok
      integer :: istart,iend

      if (is_full) then
!  -----------------------------------
!  matrix in fully dense matrix format
!  -----------------------------------
        do j=1,n
         call random_number( x_re(1:n) )
         call random_number( x_im(1:n) )
         x_re(1:n) = 2*x_re(1:n) - 1
         x_im(1:n) = 2*x_im(1:n) - 1
         do i=1,n
           aij = cmplx( x_re(i), x_im(i), kind=wp)
           if (idebug >= 2) then
               aij = -(i + (j-1)*n)
           endif



           if (i - j > kl) then
               aij = 0
           endif
           if (j - i > ku) then
               aij = 0
           endif
           A(i,j) = aij
         enddo
        enddo
       else
         isok = (ldA >= 2*kl+ku + 1)
         if (.not.isok) then
              print*,'gen_banded: invalid ldA. kl,ku,ldA ',              &
     &                                         kl,ku,ldA
              return
         endif

!  ---------------------------------------------------------------
!  note: lapack band storage format
!  AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku) <= i <= min( m, j+kl)
!  ---------------------------------------------------------------
         do j=1,n
           A(:,j) = 0
           istart=max(1,j-ku)
           iend = min(n,j+kl)
           call random_number( x_re(istart:iend))
           call random_number( x_im(istart:iend))
           x_re(istart:iend) = 2*x_re(istart:iend) - 1
           x_im(istart:iend) = 2*x_im(istart:iend) - 1
           do i=istart,iend
             aij = cmplx( x_re(i), x_im(i), kind=wp)
             AB(i,j) = aij 
           enddo
          enddo


       endif

      if (idebug >= 1) then
!      -----------------------------
!      make diagonal to be very large
!      to avoid pivoting
!      -----------------------------
       do j=1,n
         i = j
         aij = max(1.0d4,2.0d0*n*n)
         if (is_full) then
           A(i,j) = aij
         else
           AB(i,j) = aij
         endif
       enddo
      endif

      return
      end subroutine gen_banded
#undef AB
