      module banded_mod
      use iso_c_binding
      implicit none
      integer, parameter :: dp = selected_real_kind(15,100)
      integer, parameter :: sp = selected_real_kind(15,100)
      integer, parameter :: wp = dp

      interface
	subroutine Ztrsm(side,uplo,transA,diag,m,n,alpha,A,ldA,B,ldB)
	implicit none
	character side, uplo, transA, diag
	integer m,n,ldA,ldB
	complex*16 alpha, A(ldA,*), B(ldB,*)
	end subroutine Ztrsm

	subroutine Zgetrf(m,n,A,ldA,ipiv,info)
	implicit none
	integer m,n,ldA,info
	integer ipiv(*)
	complex*16 A(ldA,*)
	end subroutine

	subroutine Ztrmv( uplo, trans, diag, n, A, ldA, X, incx )
	implicit none
	character uplo, trans, diag
	integer n, ldA, incx
	complex*16 A(ldA,*), X(*)
	end subroutine Ztrmv

      end interface

      interface
         subroutine Ztrmv_sm( uplo,transA,                             &
     &      m,n,A,ldA,x,v) bind(C,name='ztrmv_sm')
         use iso_c_binding
         implicit none
         character(kind=c_char), value :: uplo,transA
         integer(kind=c_int), value :: m,n,ldA
         complex(kind=c_double_complex) :: A(*), x(*),v(*)
!         type(c_ptr), value :: A, x, v
         end subroutine Ztrmv_sm
      end interface

      contains
#include "gen_banded.F90"
#include "gen_banded_batched.F90"

#include "bandfactor.F90"
#include "bandsolve.F90"

#include "bandfactor_batched.F90"
#include "bandsolve_batched.F90"

#include "test_band_batched.F90"

      end module banded_mod
