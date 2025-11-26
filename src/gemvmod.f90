module gemvmod
   use iso_fortran_env, only: real64
   use omp_lib
   implicit none
   private

   public :: gemv_openmp_n, gemv_openmp_m, gemv_openmp_tiled, &
      gemv_openmp_tiled_improved, gemv_openmp_blocked

contains

   subroutine gemv_openmp_n(m, n, alpha, A, lda, x, beta, y)
      use iso_fortran_env, only: real64
      use omp_lib
      implicit none
      integer, intent(in) :: m, n, lda
      real(real64), intent(in) :: alpha, beta
      real(real64), intent(in) :: A(lda, *)
      real(real64), intent(in) :: x(*)
      real(real64), intent(inout) :: y(*)
      real(real64) :: ddot
      integer :: i
      real(real64) :: temp
      !$omp parallel do private(i,temp) shared(m, n, A, x, y, alpha, beta)
      do i = 1, m
         temp = ddot(n, A(i,1:n), 1, x, 1)
         y(i) = alpha * temp + beta * y(i)
      end do
      !$omp end parallel do
   end subroutine gemv_openmp_n

   subroutine gemv_openmp_m(m, n, alpha, A, lda, x, beta, y)
      use iso_fortran_env, only: real64
      use omp_lib
      implicit none
      integer, intent(in) :: m, n, lda
      real(real64), intent(in) :: alpha, beta
      real(real64), intent(in) :: A(lda, *)
      real(real64), intent(in) :: x(n)
      real(real64), intent(inout) :: y(m)
      integer :: i

      y = beta * y ! Update y with beta * y
      !$omp parallel do private(i) shared(A, x, alpha) reduction(+:y)
      do i = 1, n
         call daxpy(m, alpha*x(i), A(1:m,i), 1, y, 1)
      end do
      !$omp end parallel do
   end subroutine gemv_openmp_m

   subroutine gemv_openmp_tiled(m, n, alpha, A, lda, x, beta, y, n_x, n_y)
      use iso_fortran_env, only: real64
      use omp_lib
      implicit none

      integer, intent(in) :: m, n, lda
      real(real64), intent(in) :: alpha, beta
      real(real64), intent(in) :: A(lda, *)
      real(real64), intent(in) :: x(n)
      real(real64), intent(inout) :: y(m)
      integer, intent(in), optional :: n_x, n_y
      real(real64), allocatable :: yloc(:)


      integer :: n_x_, n_y_
      integer :: i, j, ti, mb, nb

      ! set tile sizes or defaults
      if (.not. present(n_x)) then
         n_x_ = 32
      else
         n_x_ = n_x
      end if
      if (.not. present(n_y)) then
         n_y_ = 32
      else
         n_y_ = n_y
      end if

      ! scale y by beta
      y = beta * y

      !$omp parallel default(none) &
      !$omp    shared(A, x, y, m, n, lda, alpha, n_x_, n_y_) &
      !$omp    private(i,j,ti,mb,nb,yloc)

      allocate(yloc(m))
      yloc = 0.0_real64

      ! Tile the i–j loops; collapse for better load balance
      !$omp do collapse(2) schedule(static)
      do i = 1, m, n_x_
         do j = 1, n, n_y_
            ! handle edge tiles
            mb = min(n_x_, m - i + 1)
            nb = min(n_y_, n - j + 1)
            ! perform the small GEMV into the thread‐local yloc
            call dgemv('N', mb, nb, alpha, &
               A(i, j), lda, &
               x(j), 1, &
               1.0_real64, yloc(i), 1)
         end do
      end do
      !$omp end do

      ! Safely accumulate thread‐local yloc into global y
      do ti = 1, m
         !$omp atomic
         y(ti) = y(ti) + yloc(ti)
      end do

      deallocate(yloc)
      !$omp end parallel

   end subroutine gemv_openmp_tiled

   subroutine gemv_openmp_tiled_improved(m, n, alpha, A, lda, x, beta, y, bm, bn)
      use iso_fortran_env
      implicit none
      integer,intent(in) :: m,n,lda,bm,bn
      real(real64),intent(in) :: alpha,beta
      real(real64),intent(in) :: A(lda,*), x(n)
      real(real64),intent(inout) :: y(m)
      integer :: i,j,mb,nb,t
      real(real64), allocatable :: ypriv(:,:)
      integer :: nth

      y = beta * y

!$omp parallel default(none) shared(m,n,lda,alpha,A,x,y,bm,bn,ypriv,nth)
      !$omp single
      nth = omp_get_num_threads()
      allocate(ypriv(m,nth))
      ypriv = 0.0_real64
      !$omp end single
      !$omp barrier

!$omp do schedule(static) private(i,j,mb,nb)
      do j = 1, n, bn
         nb = min(bn, n - j + 1)
         do i = 1, m, bm
            mb = min(bm, m - i + 1)
            ! Outer-product micro-kernel over tile (i:i+mb-1, j:j+nb-1)
            call dgemv('N', mb, nb, alpha, A(i,j), lda, x(j), 1, 1.0_real64, &
               ypriv(i, omp_get_thread_num()+1), 1)
         end do
      end do
!$omp end do

! Blocked reduction: each thread reduces a slice to avoid atomics
!$omp do schedule(static) private(t)
      do t = 1, nth
         y(:) = y(:) + ypriv(:,t)
      end do
!$omp end do

!$omp single
      deallocate(ypriv)
!$omp end single
!$omp end parallel
   end subroutine gemv_openmp_tiled_improved

   subroutine gemv_openmp_blocked(m, n, alpha, A, lda, x, beta, y)
      use iso_fortran_env, only: real64
      use omp_lib
      implicit none
      integer, intent(in) :: m, n, lda
      real(real64), intent(in) :: alpha, beta
      real(real64), intent(in) :: A(lda, *), x(*)
      real(real64), intent(inout) :: y(*)
      integer :: i, j, tid, nth, istart, iend
      integer :: base
      real(real64) :: xj

      ! First scale y by beta in parallel
      !$omp parallel default(none) shared(m,y,beta) private(i)
      !$omp do schedule(static)
      do i = 1, m
         y(i) = beta * y(i)
      end do
      !$omp end do
      !$omp end parallel

      !$omp parallel default(none) shared(m,n,A,lda,x,y,alpha) private(tid,nth,istart,iend,j,i,xj,base)
      tid = omp_get_thread_num()
      nth = omp_get_num_threads()

      ! partition rows [1..m] among threads: contiguous blocks
      base = (tid * m) / nth
      istart = base + 1
      iend  = ((tid + 1) * m) / nth

      ! iterate columns; A(:,j) is contiguous so A(istart:iend,j) is contiguous
      do j = 1, n
         xj = x(j)
         ! hint vectorization on inner loop
         !$omp simd
         do i = istart, iend
            y(i) = y(i) + alpha * A(i, j) * xj
         end do
      end do

      !$omp end parallel

   end subroutine gemv_openmp_blocked

end module gemvmod
