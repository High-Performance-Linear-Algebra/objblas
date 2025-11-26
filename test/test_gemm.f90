program test_gemm
   use iso_fortran_env, only: real64, output_unit, error_unit
   use omp_lib
   use gemmmod
   implicit none

   integer, parameter :: repetitions = 20
   real(real64), allocatable :: A(:,:), B(:,:), C(:,:)
   real(real64), allocatable :: C_openmp(:,:), C_naive(:,:), C_tiled(:,:), C_tiled_openmp(:,:)
   integer :: n,m,k
   integer :: i
   integer :: stat
   real(real64) :: tic, toc
   real(real64) :: time_blas, time_openmp, time_naive, time_tiled, time_tiled_openmp
   character(len=200) :: n_str, m_str, k_str

   ! Read matrix dimensions from command line
   if (command_argument_count() /= 3) then
      write(error_unit, '(A)') "Error: Expected 3 arguments (n, m, k)"
      write(error_unit, '(A)') "Usage: test_gemm n m k"
      stop 1
   end if

   call get_command_argument(1, n_str)
   call get_command_argument(2, m_str)
   call get_command_argument(3, k_str)

   read(n_str, *) n
   read(m_str, *) m
   read(k_str, *) k

   write(output_unit, '(A,I0,A,I0,A,I0)') "Matrix dimensions: n=", n, ", m=", m, ", k=", k

   allocate(A(n,k), B(k,m), C(n,m), C_openmp(n,m), C_naive(n,m), &
   & C_tiled(n,m), C_tiled_openmp(n,m), stat=stat)
   if (stat /= 0) then
      write(error_unit, '(A)') "Error: Allocation failed"
      stop 1
   end if

   ! Initialize matrices A and B
   call initialize_matrix(A, n, k)
   call initialize_matrix(B, k, m)
   C = 0.0_real64

   ! Perform matrix multiplication C = A * B + C using BLAS DGEMM
   time_blas = 0.0_real64
   do i = 1, repetitions
      C = 0.0_real64
      tic = omp_get_wtime()
      call dgemm('N', 'N', n, m, k, 1.0_real64, A, n, B, k, 1.0_real64, C, n)
      toc = omp_get_wtime()
      time_blas = time_blas + (toc - tic)
   end do
   time_blas = time_blas / repetitions
   write(output_unit, '(A,F8.4,A)') "BLAS DGEMM time: ", time_blas, " seconds"

   ! Perform matrix multiplication C_naive = A * B + C_naive using naive triple loop
   time_naive = 0.0_real64
   do i = 1, repetitions
      C_naive = 0.0_real64
      tic = omp_get_wtime()
      call matmul_jli(n, m, k, 1.0_real64, A, B, 1.0_real64, C_naive)
      toc = omp_get_wtime()
      time_naive = time_naive + (toc - tic)
   end do
   time_naive = time_naive / repetitions
   write(output_unit, '(A,F8.4,A)') "Naive DGEMM time: ", time_naive, " seconds"

   ! Verify correctness
   if (maxval(abs(C - C_naive)) < 1.0e-10_real64) then
      write(output_unit, '(A)') "Verification successful: Results match."
   else
      write(error_unit, '(A)') "Verification failed: Results do not match."
   end if

   ! Perform matrix multiplication C_openmp = A * B + C_openmp using custom OpenMP DGEMM
   time_openmp = 0.0_real64
   do i = 1, repetitions
      C_openmp = 0.0_real64
      tic = omp_get_wtime()
      call dgemm_openmp(n, m, k, 1.0_real64, A, n, B, k, 1.0_real64, C_openmp, n)
      toc = omp_get_wtime()
      time_openmp = time_openmp + (toc - tic)
   end do
   time_openmp = time_openmp / repetitions
   write(output_unit, '(A,F8.4,A)') "OpenMP DGEMM time: ", time_openmp, " seconds"

   ! Verify correctness
   if (maxval(abs(C - C_openmp)) < 1.0e-10_real64) then
      write(output_unit, '(A)') "Verification successful: Results match."
   else
      write(error_unit, '(A)') "Verification failed: Results do not match."
   end if

   ! Perform matrix multiplication C_tiled = A * B + C_tiled using custom tiled DGEMM
   time_tiled = 0.0_real64
   do i = 1, repetitions
      C_tiled = 0.0_real64
      tic = omp_get_wtime()
      call dgemm_tiled(n, m, k, 1.0_real64, A, n, B, k, 1.0_real64, C_tiled, n)
      toc = omp_get_wtime()
      time_tiled = time_tiled + (toc - tic)
   end do
   time_tiled = time_tiled / repetitions
   write(output_unit, '(A,F8.4,A)') "Tiled DGEMM time: ", time_tiled, " seconds"

   ! Verify correctness
   if (maxval(abs(C - C_tiled)) < 1.0e-10_real64) then
      write(output_unit, '(A)') "Verification successful: Results match."
   else
      write(error_unit, '(A)') "Verification failed: Results do not match."
   end if

   ! Perform matrix multiplication C_tiled_openmp = A * B + C_tiled_openmp using custom tiled OpenMP DGEMM
   time_tiled_openmp = 0.0_real64
   do i = 1, repetitions
      C_tiled_openmp = 0.0_real64
      tic = omp_get_wtime()
      call dgemm_tiled_openmp(n, m, k, 1.0_real64, A, n, B, k, 1.0_real64, C_tiled_openmp, n)
      toc = omp_get_wtime()
      time_tiled_openmp = time_tiled_openmp + (toc - tic)
   end do
   time_tiled_openmp = time_tiled_openmp / repetitions
   write(output_unit, '(A,F8.4,A)') "Tiled OpenMP DGEMM time: ", time_tiled_openmp, " seconds"

   ! Verify correctness
   if (maxval(abs(C - C_tiled_openmp)) < 1.0e-10_real64) then
      write(output_unit, '(A)') "Verification successful: Results match."
   else
      write(error_unit, '(A)') "Verification failed: Results do not match."
   end if

      ! Deallocate matrices
   deallocate(A, B, C, C_openmp, C_naive, C_tiled, C_tiled_openmp, stat=stat)
   if (stat /= 0) then
      write(error_unit, '(A)') "Error: Deallocation failed"
      stop 1
   end if

   return

contains

   subroutine initialize_matrix(A, rows, cols)
      real(real64), intent(out) :: A(rows, cols)
      integer, intent(in) :: rows, cols
      integer :: i, j

      do j = 1, cols
         do i = 1, rows
            A(i, j) = real(i + j - 1, real64)
         end do
      end do
   end subroutine initialize_matrix

end program test_gemm
