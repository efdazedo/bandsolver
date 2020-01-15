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
       integer, target :: kl_array(batchCount)
       integer, target :: ku_array(batchCount)
       integer :: info_array(batchCount)

       integer :: t1,t2,count_rate
       real(kind=wp) :: elapsed_time

       real(kind=wp) :: huge = 1.0d9

       complex(kind=wp), parameter :: one_z = dcmplx(1.0d0,0.0d0)
       real(kind=wp), parameter :: one_d = real(1.0d0,kind=wp)
       integer, parameter :: one_i = 1

       integer, parameter :: sizeof_real = c_sizeof(one_d)
       integer, parameter :: sizeof_complex = c_sizeof(one_z)
       integer, parameter :: sizeof_int = c_sizeof(one_i)

       integer(kind=c_size_t) :: nbytes
       type(c_ptr) :: d_v
       type(c_ptr) :: d_old2new, d_kl_array, d_ku_array
       type(c_ptr) :: d_A, d_B, d_xnew

       max_err = huge
       max_res = huge

      ldA = n
      allocate( old2new(n,batchCount) )
      allocate( A(ldA, n, batchCount), Aorg(ldA,n,batchCount)  )
      call gen_banded_batched( n, kl, ku, A, lda, batchCount)

!$omp parallel do private(ibatch)
      do ibatch=1,batchCount
	 Aorg(:,:,ibatch) = A(:,:,ibatch)
      enddo

      ldB = n
      ldX = n
      allocate( b(ldB,batchCount), x(ldX,batchCount), xnew(ldX,batchCount) )
      allocate( xdiff(ldX,batchCount), res(ldB,batchCount))

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
                   maxval(ku_array(1:batchCount)) )
        allocate( v(ldV,batchCount) )

       ldB = size(B,1)
       ldX = size(X,1)
       ldV = size(v,1)
!      -------------------------------
!      allocate storage on accelerator
!      -------------------------------
       nbytes = sizeof_complex
       nbytes = nbytes * size(v)
       d_v = dmalloc( nbytes )


       nbytes = sizeof_complex
       nbytes = nbytes * size(B)
       d_B = dmalloc( nbytes )
       call host2acc( d_B, c_loc(B), nbytes )

       nbytes = sizeof_complex
       nbytes = nbytes * size(A)
       d_A = dmalloc( nbytes )
       call host2acc( d_A, c_loc(A), nbytes )

       nbytes = sizeof_complex
       nbytes = nbytes * size(xnew)
       d_xnew = dmalloc( nbytes )
       call host2acc( d_xnew, c_loc(xnew), nbytes )

       nbytes = sizeof_int
       nbytes = nbytes * size(old2new)
       d_old2new = dmalloc( nbytes )
       call host2acc( d_old2new, c_loc(old2new), nbytes )

       nbytes = sizeof_int
       nbytes = nbytes * sizeof(kl_array)
       d_kl_array = dmalloc( nbytes )
       call host2acc( d_kl_array, c_loc(kl_array), nbytes )

       nbytes = sizeof_int
       nbytes = nbytes * sizeof(ku_array)
       d_ku_array = dmalloc( nbytes )
       call host2acc( d_ku_array, c_loc(ku_array), nbytes )




       call system_clock(t1,count_rate)
       call bandsolve_batched_sm(n, d_kl_array,d_ku_array,d_A,ldA,          &
     &                  d_old2new,d_B,ldB,d_xnew,ldX,d_v,ldV,batchCount)
!      call bandsolve_batched_sm(n, kl_array,ku_array,A,ldA,                &
!    &                  old2new,b,ldB,xnew,ldX,v,ldV,batchCount)
       call system_clock(t2,count_rate)
       elapsed_time = dble(t2-t1)/dble(count_rate)
       print*,'bandsolve_batched_sm took ', elapsed_time,'sec'


       nbytes = sizeof_complex
       nbytes = nbytes * size(xnew)
       call acc2host( c_loc(xnew), d_xnew, nbytes )

       deallocate( v )

!      ---------------------------------
!      deallocate storage on accelerator
!      ---------------------------------


       call dfree( d_v )
       call dfree( d_old2new )
       call dfree( d_kl_array )
       call dfree( d_ku_array )
       call dfree( d_A )
       call dfree( d_B )
       call dfree( d_xnew )

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

!      --------
!      clean up
!      --------
       deallocate( res, xdiff )
       deallocate( b, x, xnew )
       deallocate( A, Aorg )
       deallocate( old2new )

       return
       end subroutine test_band_batched
