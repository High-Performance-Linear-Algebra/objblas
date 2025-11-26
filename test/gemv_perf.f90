program gemv_roofline
   use iso_fortran_env, only: dp => real64, output_unit, input_unit
   use omp_lib
   use blasomp
   use blas
   use gemvmod
   implicit none

   integer :: reps, n, m, k, nx, ny
   real(dp) :: bw_GBs, peak_GFLOPS

   integer :: nthreads
   logical :: file_exists
   real(dp), allocatable :: x(:), y(:), A(:,:), y0(:)

   real(dp) :: alpha, beta, dt
   integer :: i, rep, stat, ierr
   real(dp) :: t0, t1

   real(dp) :: total_flops, bytes_moved, bytes_per_elem, real_small
   real(dp) :: GBs_meas_blas, GBs_meas_row, GBs_meas_col, GBs_meas_tiled
   real(dp) :: gflops_dgemv_blas, gflops_openmp_col, gflops_openmp_row, gflops_openmp_tiled, best_gflops
   real(dp) :: AI, roof_bw, roofline, err_rel
   real(dp) :: pct_blas, pct_omp, pct_omp_simd, pct_omp_schedule

   real(dp) :: time_blas, time_gemv_openmp_n, time_gemv_openmp_m, time_gemv_openmp_tiled

   !$omp parallel
   !$omp single
   nthreads = omp_get_max_threads()
   write(output_unit, *) 'Using ', nthreads, ' OpenMP threads'
   !$omp end single
   !$omp end parallel

   write(output_unit, *) 'Enter memory bandwidth in GB/s: '
   read(input_unit, *) bw_GBs
   write(output_unit, *) 'Enter peak performance in GFLOPS: '
   read(input_unit, *) peak_GFLOPS
   write(output_unit, *) 'Enter number of rows (m): '
   read(input_unit, *) m
   write(output_unit, *) 'Enter number of columns (n): '
   read(input_unit, *) n
   write(output_unit, *) 'Enter size of the x-tile (nx): '
   read(input_unit, *) nx
   write(output_unit, *) 'Enter size of the y-tile (ny): '
   read(input_unit, *) ny
   write(output_unit, *) 'Enter number of repetitions: '
   read(input_unit, *) reps

   allocate(A(m,n), x(n), y(m), y0(m))
   call random_seed()
   alpha = 1.0_dp
   beta  = 1.0_dp
   do i=1,m
      do k=1,n
         call random_number(A(i,k))
      end do
   end do
   do i=1,n
      call random_number(x(i))
   end do
   do i=1,m
      call random_number(y0(i))
   end do

   inquire(file='results_gemv.csv', exist=file_exists)
   if (file_exists) then
      open(unit=10, file='results_gemv.csv', status='old', action='write', position='append', iostat=ierr)
      if (ierr /= 0) write(output_unit,*) 'Error opening existing results_gemv.csv, IOSTAT=', ierr
   else
      open(unit=10, file='results_gemv.csv', status='new', action='write', iostat=ierr)
      if (ierr /= 0) then
         write(output_unit,*) 'Error creating results_gemv.csv, IOSTAT=', ierr
      else
         write(10,'(A)') 'OI,GFLOPS_DGEMV_BLAS,GFLOPS_OPENMP_ROW,GFLOPS_OPENMP_COL,GFLOPS_OPENMP_TILED,ROOFLINE'
      end if
   end if

   write(output_unit, *) 'Testing GEMV with m = ', m, ' and n = ', n

   y(:) = y0(:)
   t0 = omp_get_wtime()
   do rep = 1, reps
      call dgemv('N', m, n, alpha, A, m, x, 1, beta, y, 1)
   end do
   t1 = omp_get_wtime()
   time_blas = t1 - t0

   ! gemv_openmp_m is parallel over columns
   y(:) = y0(:)
   t0 = omp_get_wtime()
   do rep = 1, reps
      call gemv_openmp_m(m, n, alpha, A, m, x, beta, y)   ! column-wise
   end do
   t1 = omp_get_wtime()
   time_gemv_openmp_m = t1 - t0

   ! gemv_openmp_n is parallel over rows
   y(:) = y0(:)
   t0 = omp_get_wtime()
   do rep = 1, reps
      call gemv_openmp_n(m, n, alpha, A, m, x, beta, y)   ! row-wise
   end do
   t1 = omp_get_wtime()
   time_gemv_openmp_n = t1 - t0

   y(:) = y0(:)
   t0 = omp_get_wtime()
   do rep = 1, reps
      call gemv_openmp_blocked(m, n, alpha, A, m, x, beta, y)
   end do
   t1 = omp_get_wtime()
   time_gemv_openmp_tiled = t1 - t0

   total_flops = 2.0_dp * real(m, dp) * real(n, dp) * real(reps, dp)
   bytes_per_elem = real(storage_size(A(1,1)), dp) / 8.0_dp
   bytes_moved = ( real(m,dp)*real(n,dp) + real(n,dp) + 2.0_dp*real(m,dp) ) &
      * bytes_per_elem * real(reps, dp)

   if (bytes_moved > 0.0_dp) then
      AI = total_flops / bytes_moved
   else
      AI = 0.0_dp
   end if

   real_small = 1.0e-15_dp
   if (time_blas < real_small) time_blas = real_small
   if (time_gemv_openmp_n < real_small) time_gemv_openmp_n = real_small
   if (time_gemv_openmp_m < real_small) time_gemv_openmp_m = real_small
   if (time_gemv_openmp_tiled < real_small) time_gemv_openmp_tiled = real_small

   ! Correct labeling: gemv_openmp_n => row-wise, gemv_openmp_m => column-wise
   gflops_dgemv_blas   = total_flops / time_blas / 1.0e9_dp
   gflops_openmp_row   = total_flops / time_gemv_openmp_n / 1.0e9_dp
   gflops_openmp_col   = total_flops / time_gemv_openmp_m / 1.0e9_dp
   gflops_openmp_tiled = total_flops / time_gemv_openmp_tiled / 1.0e9_dp

   GBs_meas_blas  = bytes_moved / time_blas / 1.0e9_dp
   GBs_meas_row   = bytes_moved / time_gemv_openmp_n / 1.0e9_dp
   GBs_meas_col   = bytes_moved / time_gemv_openmp_m / 1.0e9_dp
   GBs_meas_tiled = bytes_moved / time_gemv_openmp_tiled / 1.0e9_dp

   roof_bw  = bw_GBs * AI
   roofline = min( roof_bw, real(peak_GFLOPS, dp) )

   if (roofline > 0.0_dp) then
      pct_blas         = 100.0_dp * gflops_dgemv_blas / roofline
      pct_omp          = 100.0_dp * gflops_openmp_row / roofline
      pct_omp_simd     = 100.0_dp * gflops_openmp_col / roofline
      pct_omp_schedule = 100.0_dp * gflops_openmp_tiled / roofline

      best_gflops = max(gflops_dgemv_blas, max(gflops_openmp_row, max(gflops_openmp_col, gflops_openmp_tiled)))
      err_rel = abs(roofline - best_gflops) / roofline
   else
      pct_blas = 0.0_dp; pct_omp = 0.0_dp; pct_omp_simd = 0.0_dp; pct_omp_schedule = 0.0_dp
      err_rel = 0.0_dp
   end if

   write(output_unit,'(A,F10.6)') 'Arithmetic intensity (F/byte): ', AI
   write(output_unit,'(A,F10.3)') 'Roofline (GFLOPS): ', roofline
   write(output_unit,'(A,F10.3)') 'BLAS DGEMV GFLOPS: ', gflops_dgemv_blas
   write(output_unit,'(A,F10.3)') 'OpenMP row-wise GFLOPS: ', gflops_openmp_row
   write(output_unit,'(A,F10.3)') 'OpenMP column-wise GFLOPS: ', gflops_openmp_col
   write(output_unit,'(A,F10.3)') 'OpenMP tiled GFLOPS: ', gflops_openmp_tiled
   write(output_unit,'(A,F8.3)')  'BLAS % of roofline: ', pct_blas
   write(output_unit,'(A,F8.3)')  'OpenMP row-wise % of roofline: ', pct_omp
   write(output_unit,'(A,F8.3)')  'OpenMP column-wise % of roofline: ', pct_omp_simd
   write(output_unit,'(A,F8.3)')  'OpenMP tiled % of roofline: ', pct_omp_schedule
   write(output_unit,'(A,F8.6)')  'Relative error best vs roofline: ', err_rel

   write(10,'(F12.6,1X,F12.6,1X,F12.6,1X,F12.6,1X,F12.6,1X,F12.6)') AI, gflops_dgemv_blas, &
      gflops_openmp_row, gflops_openmp_col, gflops_openmp_tiled, roofline

   close(10)
   deallocate(A, x, y, y0)
end program gemv_roofline
