      subroutine zgbmv(m,n,kl,ku,AB,ldAB, x, incx, y, incy )
!     ------------------------------------------------------------
!     matrix is  m by n  in lapack banded format
!     AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!     ------------------------------------------------------------
      implicit none
      integer, intent(in) :: m,n,kl,ku,ldAB,incx,incy
      complex(kind=wp), intent(in) :: AB(ldAB,n)
      complex(kind=wp), intent(in) :: x(1+(n-1)*incx)
      complex(kind=wp), intent(inout) :: y(1+(m-1)*incy)

      integer :: i,j,ix,iy
      complex(kind=wp) :: aij, xj

      do i=1,m
         iy = 1 + (i-1)*incy
         y(iy) = 0
      enddo

      do j=1,n
         ix = 1 + (j-1)*incx
         xj = x(ix)

         do i=max(1,j-ku),min(m,j+kl)
           aij = AB(kl+ku+1+i-j,j)
           iy = 1 + (i-1)*incy
           y(iy) = y(iy) + aij * xj
         enddo
      enddo

      return
      end subroutine zgbmv
