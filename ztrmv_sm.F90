      subroutine Ztrmv_sm(uplo,trans,diag,m,n,A,lda,x,incx)
      use iso_c_binding
      implicit none
      character, intent(in) :: uplo, trans, diag
      integer, intent(in) :: m,n,lda,incx
      complex(kind=wp), intent(in) :: A(lda,*)
      complex(kind=wp), intent(inout) :: x(*)

      call ztrmv_sm( 
