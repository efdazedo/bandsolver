      subroutine test_band_batched(n,kl,ku,                                     &
     &                 max_err, max_res,batchCount )
! % ---------------------------------------------
! % simple test for bandfactor() and bandsolver()
! % ---------------------------------------------
      use iso_c_binding
      implicit none
      integer, intent(in) :: n, kl, ku, batchCount
      real(kind=wp) :: max_err, max_res

       integer, parameter :: idebug = 1
       complex(kind=wp), pointer :: A(:,:,:)
       complex(kind=wp), pointer :: Aorg(:,:,:)
       complex(kind=wp), pointer :: x(:,:), xnew(:,:),xdiff(:,:) 
       complex(kind=wp), pointer :: b(:,:),res(:,:),v(:,:)
       integer, pointer :: old2new(:,:)

       integer :: ldA,ldB,ldX,ldV
       integer :: ibatch, i, info
       logical :: isok

       real(kind=wp) :: err 
       real(kind=wp) :: x_re(n), x_im(n)
       integer, pointer :: kl_array(:)
       integer, pointer :: ku_array(:)
       integer :: info_array(batchCount)

       integer :: t1,t2,count_rate
       real(kind=wp) :: elapsed_time

       real(kind=wp) :: huge = 1.0d9

       complex(kind=wp),parameter:: one_z = cmplx(1.0d0,0.0d0,kind=wp)
       real(kind=wp),   parameter:: one_d = real(1.0d0,kind=wp)
       integer,         parameter:: one_i = 1

       integer(kind=c_size_t),parameter:: sizeof_real = c_sizeof(one_d)
       integer(kind=c_size_t),parameter:: sizeof_cmplx = c_sizeof(one_z)
       integer(kind=c_size_t),parameter:: sizeof_int = c_sizeof(one_i)

       integer(kind=c_size_t) :: nbytes
       type(c_ptr) :: d_v
       type(c_ptr) :: d_old2new, d_kl_array, d_ku_array
       type(c_ptr) :: d_A, d_Aorg, d_b, d_res, d_x, d_xdiff, d_xnew
       type(c_ptr) :: dest, dsrc

       if (idebug >= 1) then
               print*,'sizeof_cmplx,sizeof_real,sizeof_int',             &
     &                 sizeof_cmplx,sizeof_real,sizeof_int
       endif

       max_err = huge
       max_res = huge

      ldA = n

#ifdef USE_DMALLOC
      nbytes = sizeof_int * batchCount
      d_kl_array = dmalloc( nbytes )
      d_ku_array = dmalloc( nbytes )
      call c_f_pointer(d_kl_array, kl_array, (/ batchCount /) )
      call c_f_pointer(d_ku_array, ku_array, (/ batchCount /) )

      nbytes = sizeof_int * (n*batchCount)
      d_old2new = dmalloc( nbytes )
      call c_f_pointer( d_old2new, old2new, (/ n, batchCount /) )

      nbytes = (sizeof_cmplx * ldA * n) * batchCount
      d_A = dmalloc( nbytes )
      d_Aorg = dmalloc( nbytes )
      call c_f_pointer( d_A, A, (/ ldA, n, batchCount /) )
      call c_f_pointer( d_Aorg, Aorg, (/ ldA, n, batchCount /) )
#else
      allocate( kl_array(batchCount), ku_array(batchCount) )
      allocate( old2new(n,batchCount) )
      allocate( A(ldA, n, batchCount), Aorg(ldA,n,batchCount)  )

      d_kl_array = c_loc(kl_array)
      d_ku_array = c_loc(ku_array)
      d_A = c_loc(A)
      d_Aorg = c_loc(Aorg)
#endif


      call gen_banded_batched( n, kl, ku, A, lda, batchCount)

!$omp parallel do private(ibatch)
      do ibatch=1,batchCount
	 Aorg(:,:,ibatch) = A(:,:,ibatch)
      enddo

      ldB = n
      ldX = n

#ifdef USE_DMALLOC
      nbytes = sizeof_cmplx * ldB * batchCount
      d_b = dmalloc( nbytes )
      d_res = dmalloc( nbytes )
      call c_f_pointer( d_b, b, (/ ldB, batchCount /) )
      call c_f_pointer( d_res, res, (/ ldB, batchCount /) )

      nbytes = sizeof_cmplx * ldX * batchCount
      d_x = dmalloc( nbytes )
      d_xnew = dmalloc( nbytes )
      d_xdiff = dmalloc( nbytes )
      call c_f_pointer(d_x, x, (/ ldX, batchCount /) )
      call c_f_pointer(d_xnew, xnew, (/ ldX, batchCount /) )
      call c_f_pointer(d_xdiff, xdiff, (/ ldX, batchCount /) )

#else
      allocate( b(ldB,batchCount), x(ldX,batchCount),                    &
     &          xnew(ldX,batchCount) )
      allocate( xdiff(ldX,batchCount), res(ldB,batchCount))
      d_b = c_loc(b)
      d_res = c_loc(res)

      d_x = c_loc(x)
      d_xnew = c_loc(xnew)
      d_xdiff = c_loc(xdiff)
#endif

!$omp parallel do private(ibatch,x_re,x_im,i)
      do ibatch=1,batchCount
	 call random_number(x_re(1:n))
	 call random_number(x_im(1:n))
	 x_re(1:n) = 2*x_re(1:n) - 1
	 x_im(1:n) = 2*x_im(1:n) - 1
	 do i=1,n
	   x(i,ibatch) = cmplx( x_re(i), x_im(i), kind=wp)
	 enddo
         if (idebug >= 1) then
           do i=1,n
             x(i,ibatch) = i
           enddo
         endif
       enddo

! % -------------------------
! % generate solution and rhs
! % -------------------------
!$omp  parallel do private(ibatch)
       do ibatch=1,batchCount
	  b(1:n,ibatch) = matmul( Aorg(1:n,1:n,ibatch), x(1:n,ibatch))
       enddo
! 
! % ---------------------
! % perform factorization
! %
! % note new bandwidth may be larger due to pivoting
! % kl2 ~ 2*(kl+ku), ku2 ~ 2*ku
! % ---------------------

       call system_clock(t1,count_rate)
       call bandfactor_batched(n,A,lda,old2new,kl_array,ku_array,           &
     &              info_array, batchCount)
       call system_clock(t2,count_rate)
       elapsed_time = dble(t2-t1)/dble(count_rate)
       print*,'bandfactor_batched took ',elapsed_time,'sec'

       isok = all( info_array(1:batchCount).eq.0 )
       if (.not.isok) then
	  do ibatch=1,batchCount
            info = info_array(ibatch)
	    if (info.ne.0) then
	      print*,'ibatch,info ', ibatch,info
	    endif
	  enddo
	  return
	endif

        ldV = max( maxval(kl_array(1:batchCount)),                        &
     &             maxval(ku_array(1:batchCount)) )
#ifdef USE_DMALLOC
        nbytes = (sizeof_cmplx * ldV) * batchCount
        d_v = dmalloc( nbytes )
        call c_f_pointer( d_v, v, (/ ldV, batchCount /) )
#else
        allocate( v(ldV,batchCount) )
        d_v = c_loc(d_v)
#endif

       ldB = size(B,1)
       ldX = size(X,1)
       ldV = size(v,1)



!      --------------------------------------------------
!      1st call to warm up cache or perform data transfer
!      --------------------------------------------------
       call bandsolve_batched_sm(n, kl_array,ku_array,A,ldA,                &
     &                  old2new,b,ldB,xnew,ldX,v,ldV,batchCount)

       call system_clock(t1,count_rate)
       call bandsolve_batched_sm(n, kl_array,ku_array,A,ldA,                &
     &                  old2new,b,ldB,xnew,ldX,v,ldV,batchCount)
       call system_clock(t2,count_rate)
       elapsed_time = dble(t2-t1)/dble(count_rate)
       print*,'bandsolve_batched_sm took ', elapsed_time,'sec'





       max_err = 0
       max_res = 0

!      ----------------------------------
!      compute  xdiff(:) = xnew(:) - x(:)
!      ----------------------------------
       do ibatch=1,batchCount
	  do i=1,n
            xdiff(i,ibatch) = xnew(i,ibatch) - x(i,ibatch)
            err = abs( xdiff(i,ibatch) )
	    max_err = max( max_err, err )
	  enddo
       enddo

!      -------------------------------------------
!      compute residual vector, res = Aorg * xdiff
!      -------------------------------------------
       do ibatch=1,batchCount
	res(1:n,ibatch) = matmul( Aorg(1:n,1:n,ibatch), xdiff(1:n,ibatch))
       enddo

       do ibatch=1,batchCount
       do i=1,n
	  max_res = max( max_res, abs(res(i,ibatch)) )
       enddo
       enddo

!      ---------------------------------
!      deallocate storage 
!      ---------------------------------

#ifdef USE_DMALLOC

       call dfree( d_v )
       call dfree( d_old2new )
       call dfree( d_kl_array )
       call dfree( d_ku_array )

       call dfree( d_A )
       call dfree( d_Aorg )

       call dfree( d_b )
       call dfree( d_res )

       call dfree( d_xnew )
       call dfree( d_x )
       call dfree( d_xdiff )
#else
       deallocate( v )
       deallocate( old2new )
       deallocate( kl_array )
       deallocate( ku_array )

       deallocate( A )
       deallocate( Aorg )

       deallocate( b )
       deallocate( res )

       deallocate( xnew )
       deallocate( x )
       deallocate( xdiff )
#endif

       return
       end subroutine test_band_batched
