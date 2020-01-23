      module banded_mod
      use iso_c_binding
      implicit none
      integer, parameter :: dp = selected_real_kind(15,100)
      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: wp = dp

      interface
        subroutine dsync() bind(c,name='dsync')
        use iso_c_binding
        end subroutine dsync

        function dmalloc( nbytes )                                       &
     &      bind(c,name='dmalloc') result(ptr)
        use iso_c_binding
        implicit none
        integer(kind=c_size_t), value :: nbytes
        type(c_ptr) :: ptr
        end function dmalloc

        subroutine dfree( ptr ) bind(c,name='dfree')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: ptr
        end subroutine dfree

        subroutine host2acc( dest, src, nbytes )                         &
     &       bind(c,name='host2acc')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: dest, src 
        integer(kind=c_size_t), value :: nbytes
        end subroutine host2acc


        subroutine acc2host( dest, src, nbytes )                         &
     &       bind(c,name='acc2host')
        use iso_c_binding
        implicit none
        type(c_ptr), value :: dest, src 
        integer(kind=c_size_t), value :: nbytes
        end subroutine acc2host


      end interface

      interface
       subroutine Zgemm(transA,transB,m,n,k,                             &
     &   alpha, A,lda, B,ldb, beta, C, ldc )
       implicit none
       character transA, transB
       integer m,n,k,lda,ldb,ldc
       complex*16 alpha,beta
       complex*16 A(lda,*)
       complex*16 B(ldb,*)
       complex*16 C(ldc,*)
       end subroutine Zgemm

       subroutine Cgemm(transA,transB,m,n,k,                             &
     &   alpha, A,lda, B,ldb, beta, C, ldc )
       implicit none
       character transA, transB
       integer m,n,k,lda,ldb,ldc
       complex*8 alpha,beta
       complex*8 A(lda,*)
       complex*8 B(ldb,*)
       complex*8 C(ldc,*)
       end subroutine Cgemm
      end interface

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
        end subroutine Zgetrf

        subroutine Ztrmv( uplo, trans, diag, n, A, ldA, X, incx )
        implicit none
        character uplo, trans, diag
        integer n, ldA, incx
        complex*16 A(ldA,*), X(*)
        end subroutine Ztrmv
      end interface



      interface
        subroutine Ctrsm(side,uplo,transA,diag,m,n,alpha,A,ldA,B,ldB)
        implicit none
        character side, uplo, transA, diag
        integer m,n,ldA,ldB
        complex*8 alpha, A(ldA,*), B(ldB,*)
        end subroutine Ctrsm

        subroutine Cgetrf(m,n,A,ldA,ipiv,info)
        implicit none
        integer m,n,ldA,info
        integer ipiv(*)
        complex*8 A(ldA,*)
        end subroutine Cgetrf

        subroutine Ctrmv( uplo, trans, diag, n, A, ldA, X, incx )
        implicit none
        character uplo, trans, diag
        integer n, ldA, incx
        complex*8 A(ldA,*), X(*)
        end subroutine Ctrmv
      end interface


      interface
        subroutine Dtrsm(side,uplo,transA,diag,m,n,alpha,A,ldA,B,ldB)
        implicit none
        character side, uplo, transA, diag
        integer m,n,ldA,ldB
        real*8 alpha, A(ldA,*), B(ldB,*)
        end subroutine Dtrsm

        subroutine Dgetrf(m,n,A,ldA,ipiv,info)
        implicit none
        integer m,n,ldA,info
        integer ipiv(*)
        real*8 A(ldA,*)
        end subroutine Dgetrf

        subroutine Dtrmv( uplo, trans, diag, n, A, ldA, X, incx )
        implicit none
        character uplo, trans, diag
        integer n, ldA, incx
        real*8 A(ldA,*), X(*)
        end subroutine Dtrmv
      end interface


      interface
        subroutine Strsm(side,uplo,transA,diag,m,n,alpha,A,ldA,B,ldB)
        implicit none
        character side, uplo, transA, diag
        integer m,n,ldA,ldB
        real*4 alpha, A(ldA,*), B(ldB,*)
        end subroutine Strsm

        subroutine Sgetrf(m,n,A,ldA,ipiv,info)
        implicit none
        integer m,n,ldA,info
        integer ipiv(*)
        real*4 A(ldA,*)
        end subroutine Sgetrf

        subroutine Strmv( uplo, trans, diag, n, A, ldA, X, incx )
        implicit none
        character uplo, trans, diag
        integer n, ldA, incx
        real*4 A(ldA,*), X(*)
        end subroutine Strmv
      end interface

#ifndef USE_GPU
      interface
         subroutine Ztrmv_sm( uplo,transA,diag,                          &
     &      m,n,A,ldA,x,v) bind(C,name='ztrmv_sm')
         use iso_c_binding
         implicit none
         character(kind=c_char), value :: uplo,transA,diag
         integer(kind=c_int), value :: m,n,ldA
         complex(kind=c_double_complex) :: A(*), x(*),v(*)
!         type(c_ptr), value :: A, x, v
         end subroutine Ztrmv_sm
      end interface


      interface
         subroutine Ctrmv_sm( uplo,transA,diag,                               &
     &      m,n,A,ldA,x,v) bind(C,name='ctrmv_sm')
         use iso_c_binding
         implicit none
         character(kind=c_char), value :: uplo,transA,diag
         integer(kind=c_int), value :: m,n,ldA
         complex(kind=c_float_complex) :: A(*), x(*),v(*)
         end subroutine Ctrmv_sm
      end interface
#endif

      interface

        subroutine bandsolve_batched_sm( n, kl_array, ku_array,          &
     &     A, ldA, old2new, b, ldB, x, ldX, v, ldV, is_full,batchCount)  &
     &     bind(C,name='bandsolve_batched_sm')
        use iso_c_binding
        implicit none
        integer(kind=c_int),value:: n, ldA, ldB, ldX, ldV, batchCount
        complex(kind=c_double_complex) :: A(*), b(*), x(*), v(*)
        integer(kind=c_int) :: kl_array(*), ku_array(*)
        integer(kind=c_int) :: old2new(*)
        logical(kind=c_bool) , value :: is_full
        !type(c_ptr), value :: A, b, x, v, kl_array, ku_array,old2new
        end subroutine bandsolve_batched_sm

      end interface

      contains
#include "gen_banded.F90"
#include "gen_banded_batched.F90"

#include "bandfactor.F90"
#include "bandfactor_batched.F90"


#ifndef USE_GPU
#include "ztrmv_smf.F90"
#include "bandsolve.F90"
#include "bandsolve_batched.F90"
#endif

#include "test_band_batched.F90"

      end module banded_mod
