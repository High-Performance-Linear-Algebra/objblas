module gemmmod
   use iso_fortran_env, only: real64
   use omp_lib
   implicit none

   private
   public :: dgemm_openmp, matmul_jli, dgemm_tiled, dgemm_tiled_openmp &
      , matmul_ijl, matmul_jil

contains

   !> dgemm_openmp: OpenMP-parallelized general matrix-matrix multiply
   !> Computes C = alpha * A * B + beta * C
   !>
   !> Uses outer-product form (j-l-i loop order) for better cache locality
   !> and vectorization. The innermost loop is SIMD-vectorized and the outer
   !> loops are parallelized with OpenMP.
   !>
   !> @param m      Number of rows of A and C
   !> @param n      Number of columns of B and C
   !> @param k      Number of columns of A and rows of B
   !> @param alpha  Scalar multiplier for A*B
   !> @param A      Input matrix A (m x k)
   !> @param lda    Leading dimension of A
   !> @param B      Input matrix B (k x n)
   !> @param ldb    Leading dimension of B
   !> @param beta   Scalar multiplier for C
   !> @param C      Input/output matrix C (m x n)
   !> @param ldc    Leading dimension of C
   subroutine dgemm_openmp(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
      use iso_fortran_env, only: real64
      implicit none
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      real(real64), intent(in) :: alpha, beta
      real(real64), intent(in) :: A(lda, *)
      real(real64), intent(in) :: B(ldb, *)
      real(real64), intent(inout) :: C(ldc, *)

      ! Local variables
      integer :: i, j, l
      real(real64) :: blj

      ! Parallel scaling of C
      !$omp parallel default(none) shared(C,beta,m,n,ldc,A,B,alpha,k,lda,ldb) private(i,j,l,blj)
      !$omp do schedule(static)
      do j = 1, n
         !$omp simd
         do i = 1, m
            C(i,j) = beta * C(i,j)
         end do
      end do
      !$omp end do

      ! Parallelized and vectorized matrix multiply (outer-product form)
      !$omp do collapse(2) schedule(static)
      do j = 1, n
         do l = 1, k
            blj = alpha * B(l,j)
            !$omp simd
            do i = 1, m
               C(i,j) = C(i,j) + A(i,l) * blj
            end do
         end do
      end do
      !$omp end do
      !$omp end parallel

   end subroutine dgemm_openmp

   !> dgemm_tiled: Cache-blocked (tiled) matrix-matrix multiply
   !> Computes C = alpha * A * B + beta * C
   !>
   !> Uses cache blocking (tiling) to improve data reuse in cache hierarchy.
   !> The matrix is divided into smaller tiles that fit in cache, reducing
   !> cache misses and improving performance for large matrices.
   !>
   !> @param m       Number of rows of A and C
   !> @param n       Number of columns of B and C
   !> @param k       Number of columns of A and rows of B
   !> @param alpha   Scalar multiplier for A*B
   !> @param A       Input matrix A (m x k)
   !> @param lda     Leading dimension of A
   !> @param B       Input matrix B (k x n)
   !> @param ldb     Leading dimension of B
   !> @param beta    Scalar multiplier for C
   !> @param C       Input/output matrix C (m x n)
   !> @param ldc     Leading dimension of C
   !> @param tile_m  Optional tile size for m dimension (default: 64)
   !> @param tile_n  Optional tile size for n dimension (default: 64)
   !> @param tile_k  Optional tile size for k dimension (default: 64)
   subroutine dgemm_tiled(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, tile_m, tile_n, tile_k)
      use iso_fortran_env, only: real64
      implicit none
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      integer, intent(in), optional :: tile_m, tile_n, tile_k
      real(real64), intent(in) :: alpha, beta
      real(real64), intent(in) :: A(lda, *)
      real(real64), intent(in) :: B(ldb, *)
      real(real64), intent(inout) :: C(ldc, *)

      ! Local variables
      integer :: i, j, l, ii, jj, ll
      integer :: ts_m, ts_n, ts_k
      integer :: i_end, j_end, l_end
      real(real64) :: temp

      ! Set tile sizes (default 64)
      ts_m = 64
      ts_n = 64
      ts_k = 64
      if (present(tile_m)) ts_m = tile_m
      if (present(tile_n)) ts_n = tile_n
      if (present(tile_k)) ts_k = tile_k

      ! Scale C by beta
      do j = 1, n
         do i = 1, m
            C(i,j) = beta * C(i,j)
         end do
      end do

      ! Tiled matrix multiplication with non-square tiles
      do jj = 1, n, ts_n
         j_end = min(jj + ts_n - 1, n)
         do ll = 1, k, ts_k
            l_end = min(ll + ts_k - 1, k)
            do ii = 1, m, ts_m
               i_end = min(ii + ts_m - 1, m)

               ! Multiply tile
               do j = jj, j_end
                  do l = ll, l_end
                     temp = alpha * B(l,j)
                     do i = ii, i_end
                        C(i,j) = C(i,j) + A(i,l) * temp
                     end do
                  end do
               end do

            end do
         end do
      end do

   end subroutine dgemm_tiled

   !> dgemm_tiled_openmp: OpenMP-parallelized cache-blocked matrix-matrix multiply
   !> Computes C = alpha * A * B + beta * C
   !>
   !> Combines cache blocking with OpenMP parallelization for optimal performance.
   !> Uses a thread-private buffer (Cbuf) to accumulate results for each tile,
   !> minimizing cache conflicts and false sharing between threads.
   !>
   !> Key optimizations:
   !> - Cache blocking to fit tiles in L1/L2 cache
   !> - Thread-private tile buffers eliminate write conflicts
   !> - SIMD vectorization of innermost loop
   !> - Static scheduling for load balancing
   !>
   !> @param m       Number of rows of A and C
   !> @param n       Number of columns of B and C
   !> @param k       Number of columns of A and rows of B
   !> @param alpha   Scalar multiplier for A*B
   !> @param A       Input matrix A (m x k)
   !> @param lda     Leading dimension of A
   !> @param B       Input matrix B (k x n)
   !> @param ldb     Leading dimension of B
   !> @param beta    Scalar multiplier for C
   !> @param C       Input/output matrix C (m x n)
   !> @param ldc     Leading dimension of C
   !> @param tile_m  Optional tile size for m dimension (default: 64, max: 128)
   !> @param tile_n  Optional tile size for n dimension (default: 64, max: 128)
   !> @param tile_k  Optional tile size for k dimension (default: 64)
   subroutine dgemm_tiled_openmp(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, &
      tile_m, tile_n, tile_k)
      use iso_fortran_env, only: real64
      implicit none

      integer, intent(in) :: m, n, k, lda, ldb, ldc
      integer, intent(in), optional :: tile_m, tile_n, tile_k
      real(real64), intent(in) :: alpha, beta
      real(real64), intent(in) :: A(lda, *)
      real(real64), intent(in) :: B(ldb, *)
      real(real64), intent(inout) :: C(ldc, *)

      integer :: ts_m, ts_n, ts_k
      integer :: ii, jj, ll
      integer :: i, j, l
      integer :: i_end, j_end, l_end
      integer :: ib, jb
      real(real64) :: tmp

      ! ---- MAXIMUM tile sizes (adjust safely for your CPU cache) ----
      integer, parameter :: MAX_TS_M = 128
      integer, parameter :: MAX_TS_N = 128

      ! Local tile buffer, fixed size (thread-private due to OpenMP)
      real(real64) :: Cbuf(MAX_TS_M, MAX_TS_N)

      ! Default tile sizes
      ts_m = 64
      ts_n = 64
      ts_k = 64
      if (present(tile_m)) ts_m = min(tile_m, MAX_TS_M)
      if (present(tile_n)) ts_n = min(tile_n, MAX_TS_N)
      if (present(tile_k)) ts_k = tile_k

      !$omp parallel default(none) &
      !$omp shared(m,n,k,ts_m,ts_n,ts_k,A,B,C,alpha,beta,lda,ldb,ldc) &
      !$omp private(ii,jj,ll,i,j,l,i_end,j_end,l_end,ib,jb,Cbuf,tmp)
      !$omp do collapse(2) schedule(static)
      do jj = 1, n, ts_n
         do ii = 1, m, ts_m

            ! Work tile bounds
            i_end = min(ii + ts_m - 1, m)
            j_end = min(jj + ts_n - 1, n)

            ib = i_end - ii + 1   ! actual tile height
            jb = j_end - jj + 1   ! actual tile width

            ! -------------------------------
            ! Load and scale C tile: Cbuf = beta * C
            ! -------------------------------
            do j = 1, jb
               do i = 1, ib
                  Cbuf(i, j) = beta * C(ii + i - 1, jj + j - 1)
               end do
            end do

            ! -------------------------------
            ! Accumulate over all K tiles
            ! -------------------------------
            do ll = 1, k, ts_k
               l_end = min(ll + ts_k - 1, k)

               do l = ll, l_end
                  do j = 1, jb
                     ! scalar needed for whole column
                     tmp = alpha * B(l, jj + j - 1)

                     !$omp simd
                     do i = 1, ib
                        Cbuf(i, j) = Cbuf(i, j) + A(ii + i - 1, l) * tmp
                     end do
                  end do
               end do
            end do

            ! -------------------------------
            ! Write tile back to C
            ! -------------------------------
            do j = 1, jb
               do i = 1, ib
                  C(ii + i - 1, jj + j - 1) = Cbuf(i, j)
               end do
            end do

         end do
      end do
      !$omp end do
      !$omp end parallel

   end subroutine dgemm_tiled_openmp


   !> matmul_ijl: Simple matrix-matrix multiply with i-j-l loop order
   !> Computes C = alpha * A * B + beta * C
   !>
   !> This is a naive implementation with i-j-l loop ordering (row-wise access).
   !> The innermost loop over l has poor cache locality due to stride-k access
   !> patterns in both A and B matrices. This results in many cache misses.
   !>
   !> Performance: POOR - Use for pedagogical purposes or very small matrices only
   !>
   !> @param n      Number of rows of A and C
   !> @param m      Number of columns of B and C
   !> @param k      Number of columns of A and rows of B
   !> @param alpha  Scalar multiplier for A*B
   !> @param A      Input matrix A (n x k)
   !> @param B      Input matrix B (k x m)
   !> @param beta   Scalar multiplier for C
   !> @param C      Input/output matrix C (n x m)
   subroutine matmul_ijl(n,m,k,alpha,A,B,beta,C)
      use iso_fortran_env, only: real64
      implicit none
      integer, intent(in) :: n, m, k
      real(real64), intent(in) :: alpha
      real(real64), intent(in) :: A(n,k)
      real(real64), intent(in) :: B(k,m)
      real(real64), intent(in) :: beta
      real(real64), intent(inout) :: C(n,m)
      ! Local variables
      integer :: i, j, l
      real(real64) :: sum
      ! Matrix multiplication
      ! C = alpha * A * B + beta * C
      do i = 1, m
         do j = 1, n
            C(i,j) = beta * C(i,j)
            do l = 1, k
               C(i,j) = C(i,j) + alpha * A(i,l) * B(l,j)
            end do
         end do
      end do
   end subroutine matmul_ijl

   !> matmul_jil: Simple matrix-matrix multiply with j-i-l loop order
   !> Computes C = alpha * A * B + beta * C
   !>
   !> This implementation uses j-i-l loop ordering (column-wise access for C).
   !> Better cache locality than i-j-l since C is accessed column-wise (Fortran
   !> column-major order), but still suffers from poor locality in the innermost
   !> l loop which has stride-k access patterns.
   !>
   !> Performance: POOR to MODERATE - Better than i-j-l but not optimal
   !>
   !> @param n      Number of rows of A and C
   !> @param m      Number of columns of B and C
   !> @param k      Number of columns of A and rows of B
   !> @param alpha  Scalar multiplier for A*B
   !> @param A      Input matrix A (n x k)
   !> @param B      Input matrix B (k x m)
   !> @param beta   Scalar multiplier for C
   !> @param C      Input/output matrix C (n x m)
   subroutine matmul_jil(n,m,k,alpha,A,B,beta,C)
      use iso_fortran_env, only: real64
      implicit none
      integer, intent(in) :: n, m, k
      real(real64), intent(in) :: alpha
      real(real64), intent(in) :: A(n,k)
      real(real64), intent(in) :: B(k,m)
      real(real64), intent(in) :: beta
      real(real64), intent(inout) :: C(n,m)
      ! Local variables
      integer :: i, j, l
      real(real64) :: sum
      ! Matrix multiplication
      ! C = alpha * A * B + beta * C
      do j = 1, m
         do i = 1, n
            C(i,j) = beta * C(i,j)
            do l = 1, k
               C(i,j) = C(i,j) + alpha * A(i,l) * B(l,j)
            end do
         end do
      end do
   end subroutine matmul_jil

   !> matmul_jli: Matrix-matrix multiply with j-l-i loop order (outer-product form)
   !> Computes C = alpha * A * B + beta * C
   !>
   !> This implementation uses j-l-i loop ordering, which is the outer-product form.
   !> The innermost i loop has stride-1 access to both A(:,l) and C(:,j), providing
   !> excellent cache locality and enabling efficient SIMD vectorization.
   !>
   !> This is one of the best loop orderings for matrix multiplication in Fortran
   !> (column-major layout) and forms the basis for high-performance GEMM implementations.
   !>
   !> Performance: GOOD - Optimal cache locality and vectorization potential
   !>
   !> @param n      Number of rows of A and C
   !> @param m      Number of columns of B and C
   !> @param k      Number of columns of A and rows of B
   !> @param alpha  Scalar multiplier for A*B
   !> @param A      Input matrix A (n x k)
   !> @param B      Input matrix B (k x m)
   !> @param beta   Scalar multiplier for C
   !> @param C      Input/output matrix C (n x m)
   subroutine matmul_jli(n,m,k,alpha,A,B,beta,C)
      use iso_fortran_env, only: real64
      implicit none
      integer, intent(in) :: n, m, k
      real(real64), intent(in) :: alpha, A(n,k), B(k,m), beta
      real(real64), intent(inout) :: C(n,m)
      integer :: i, j, l
      C = beta * C
      do j = 1, n
         do l = 1, k
            do i = 1, m
               C(i,j) = C(i,j) + alpha * A(i,l) * B(l,j)
            end do
         end do
      end do
   end subroutine matmul_jli



end module gemmmod
