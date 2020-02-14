      subroutine test_band_batched(n,kl,ku,                                     &
     &                 max_err, max_res, is_full,batchCount )
! % ---------------------------------------------
! % simple test for bandfactor() and bandsolver()
! % ---------------------------------------------
      use iso_c_binding
      implicit none
      integer, intent(in) :: n, kl, ku, batchCount
      logical, intent(in)  :: is_full 
      real(kind=wp), intent(out) :: max_err, max_res

       integer, parameter :: idebug = 0
       complex(kind=wp), pointer :: A(:,:,:)
       complex(kind=wp), pointer :: Aorg(:,:,:)
       complex(kind=wp), pointer :: x(:,:), xnew(:,:),xdiff(:,:) 
       complex(kind=wp), pointer :: b(:,:),res(:,:),v(:,:)
       integer, pointer :: old2new(:,:)

       integer :: ldA,ldB,ldX,ldV
       integer :: ibatch, i, j, info
       integer :: min_kl, max_kl, avg_kl
       integer :: min_ku, max_ku, avg_ku
       logical :: isok

       real(kind=wp) :: err 
       real(kind=wp) :: x_re(n), x_im(n)
       integer :: icount(n)
       integer, pointer :: kl_array(:)
       integer, pointer :: ku_array(:)
       integer :: info_array(batchCount)
       integer :: mm,nn,ld1, inc1, inc2
       integer :: istat

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

       logical, parameter :: is_diag_dominant = .true.
       logical(kind=c_bool) :: is_full_c  


       if (idebug >= 1) then
               print*,'sizeof_cmplx,sizeof_real,sizeof_int',             &
     &                 sizeof_cmplx,sizeof_real,sizeof_int
       endif

       max_err = huge
       max_res = huge

       if (is_full) then
          ldA = n
       else
!         ------------------------------
!         use lapack band storage format
!         ------------------------------
          ldA  = 2*kl + ku + 1
       endif

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
      if (.not.c_associated(d_A)) then
         print*,'dmalloc return null in d_A'
         stop '** out of memory **'
      endif
      d_Aorg = dmalloc( nbytes )
      if (.not.c_associated(d_Aorg)) then
         print*,'dmalloc return null in d_Aorg'
         stop '** out of memory **'
      endif

      call c_f_pointer( d_A, A, (/ ldA, n, batchCount /) )
      call c_f_pointer( d_Aorg, Aorg, (/ ldA, n, batchCount /) )
#else
      allocate( kl_array(batchCount), ku_array(batchCount) )
      allocate( old2new(n,batchCount) )
      allocate(A(ldA,n,batchCount),Aorg(ldA,n,batchCount),stat=istat)
      if (istat.ne.0) then
         print*,'allocate(A,Aorg) return istat'
         stop '** out of memory **'
      endif

      d_kl_array = c_loc(kl_array)
      d_ku_array = c_loc(ku_array)
      d_A = c_loc(A)
      d_Aorg = c_loc(Aorg)
#endif
      kl_array(:) = kl
      ku_array(:) = ku


      call gen_banded_batched( n, kl, ku, A, lda,                        &
     &              is_full, is_diag_dominant, batchCount)

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
       if (is_full) then
!$omp  parallel do private(ibatch)
         do ibatch=1,batchCount
	  b(1:n,ibatch) = matmul( Aorg(1:n,1:n,ibatch), x(1:n,ibatch))
         enddo
       else
         ld1 = size(Aorg,1)
         inc1 = 1
         inc2 = 1
         mm = n
         nn = n
!$omp    parallel do private(ibatch)
         do ibatch=1,batchCount
           call Zgbmv(mm,nn,kl,ku,Aorg(1:ld1,1:n,ibatch),ld1,                &
     &                x(1:n,ibatch),inc1,b(1:n,ibatch),inc2 )
         enddo
        endif
! 
! % ---------------------
! % perform factorization
! %
! % note new bandwidth may be larger due to pivoting
! % kl2 ~ 2*(kl+ku), ku2 ~ 2*ku
! % ---------------------

       call system_clock(t1,count_rate)
       call bandfactor_batched(n,A,lda,old2new,kl_array,ku_array,           &
     &              info_array, is_full, batchCount)
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

!       -------------
!       check old2new
!       -------------
        if (idebug >= 1) then
          do ibatch=1,batchCount
             icount(1:n) = 0
             do i=1,n
               j = old2new(i,ibatch)
               isok = (1 <= j) .and. (j <= n)
               if (.not.isok) then
                 print*,'test_band_batched:invalid old2new'
                 print*,'ibatch,i,old2new ',ibatch,i,old2new(i,ibatch)
                 stop '** error in test_band_batched '
               endif
               icount(j) = icount(j) + 1
             enddo

             do i=1,n
                isok = (icount(i).eq.1)
                if (.not.isok) then
                 print*,'test_band_batched: invalid icount'
                 print*,'i,icount(i) ',   i, icount(i)
                 stop '** error in test_band_batched '
                endif
             enddo
           enddo
          endif


        min_kl = minval( kl_array(1:batchCount) )
        max_kl = maxval( kl_array(1:batchCount) )
        avg_kl = sum( kl_array(1:batchCount) )/batchCount

        min_ku = minval( ku_array(1:batchCount) )
        max_ku = maxval( ku_array(1:batchCount) )
        avg_ku = sum( ku_array(1:batchCount) )/batchCount

        if (is_full) then
          ldV = max(1,max(max_kl,max_ku))
        else
          ldV = max(1,kl + ku)
        endif

        if (idebug >= 1) then
            min_kl = minval( kl_array(1:batchCount) )
            max_kl = maxval( kl_array(1:batchCount) )
            avg_kl = sum( kl_array(1:batchCount) )/batchCount

            min_ku = minval( ku_array(1:batchCount) )
            max_ku = maxval( ku_array(1:batchCount) )
            avg_ku = sum( ku_array(1:batchCount) )/batchCount


            print*,'min_kl, avg_kl, max_kl', min_kl, avg_kl, max_kl
            print*,'min_ku, avg_ku, max_ku', min_ku, avg_ku, max_ku
        endif

#ifdef USE_DMALLOC
        nbytes = (sizeof_cmplx * ldV) * batchCount
        d_v = dmalloc( nbytes )
        call c_f_pointer( d_v, v, (/ ldV, batchCount /) )
#else
        allocate( v(ldV,batchCount) )
        d_v = c_loc(v)
#endif

       ldB = size(B,1)
       ldX = size(X,1)
       ldV = size(v,1)



!      --------------------------------------------------
!      1st call to warm up cache or perform data transfer
!      --------------------------------------------------
       is_full_c = is_full
       call bandsolve_batched_sm(n, kl_array,ku_array,A,ldA,                &
     &             old2new,b,ldB,xnew,ldX,v,ldV,is_full_c,batchCount)

       call system_clock(t1,count_rate)
       call bandsolve_batched_sm(n, kl_array,ku_array,A,ldA,                &
     &             old2new,b,ldB,xnew,ldX,v,ldV,is_full_c,batchCount)
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
       if (is_full) then
!$omp  parallel do private(ibatch)
         do ibatch=1,batchCount
	  res(1:n,ibatch) = matmul( Aorg(1:n,1:n,ibatch), xdiff(1:n,ibatch))
         enddo
       else
         ld1 = size(Aorg,1)
         inc1 = 1
         inc2 = 1
         mm = n
         nn = n
!$omp    parallel do private(ibatch)
         do ibatch=1,batchCount
           call Zgbmv(mm,nn,kl,ku,Aorg(1:ld1,1:n,ibatch),ld1,             &
     &                xdiff(1:n,ibatch),inc1,res(1:n,ibatch),inc2 )
         enddo
        endif
            

       max_res = 0
!$omp  parallel do private(i,ibatch) reduction(max:max_res)
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
