program axpy_roofline
   use iso_fortran_env, only: dp => real64, output_unit, input_unit
   use omp_lib
   use blasomp
   use blas
   implicit none

   ! Inputs
   integer :: reps, n, nlocal, k
   real(dp) :: bw_GBs,  peak_GFLOPS

   ! Data and timing
   real(dp), allocatable :: x(:), y(:)
   real(dp) :: alpha, dt
   integer :: i, rep, stat
   real(dp) :: t0, t1

   ! Metrics
   real(dp) :: total_flops, bytes_moved, bytes_per_elem
   real(dp) :: GBs_meas_blas, GBs_meas_omp, GBs_meas_omp_simd, GBs_meas_omp_schedule
   real(dp) :: gflops_blas, gflops_omp, gflops_omp_simd, gflops_omp_schedule
   real(dp) :: AI, roof_bw, roofline, err_rel
   real(dp) :: pct_blas, pct_omp, pct_omp_simd, pct_omp_schedule

   ! Additional variables for AXPY implementations
   real(dp) :: time_blas, time_omp, time_omp_simd, time_omp_schedule

   ! Read system parameters
   write(output_unit, *) 'Enter memory bandwidth in GB/s: '
   read(input_unit, *) bw_GBs
   write(output_unit, *) 'Enter peak performance in GFLOPS: '
   read(input_unit, *) peak_GFLOPS
   ! Read sizes and number of repetitions from user
   write(output_unit, *) 'Enter problem size (n): '
   read(input_unit, *) n
   write(output_unit, *) 'Enter number of repetitions: '
   read(input_unit, *) reps

   open(unit=10, file='results.csv', status='unknown', action='write')
   write(10, '(A)') 'OI,GFLOPS_BLAS,GFLOPS_OMP,GFLOPS_OMP_SIMD,GFLOPS_OMP_SCHEDULE'

   ! Test AXPY with the specified problem size
   write(output_unit, *) 'Testing AXPY with n = ', n
   allocate(x(n), y(n))
   x = [(real(i, dp), i = 1, n)]
   y = 0.0_dp

   alpha = 2.0_dp

   time_blas = 0.0_dp
   time_omp = 0.0_dp
   time_omp_simd = 0.0_dp
   time_omp_schedule = 0.0_dp

   ! Run AXPY implementations
   do rep = 1, reps
      ! Reset y for each repetition
      y = 0.0_dp
      
      ! BLAS implementation
      t0 = omp_get_wtime()
      call daxpy(n, alpha, x, 1, y, 1)
      t1 = omp_get_wtime()
      time_blas = time_blas + (t1 - t0)

      ! Reset y for next test
      y = 0.0_dp
      
      ! OpenMP implementation
      t0 = omp_get_wtime()
      call axpy_omp(n, alpha, x, 1, y, 1)
      t1 = omp_get_wtime()
      time_omp = time_omp + (t1 - t0)

      ! Reset y for next test
      y = 0.0_dp
      
      ! OpenMP + SIMD implementation
      t0 = omp_get_wtime()
      call axpy_omp_simd(n, alpha, x, 1, y, 1)
      t1 = omp_get_wtime()
      time_omp_simd = time_omp_simd + (t1 - t0)

      ! Reset y for next test
      y = 0.0_dp
      
      ! OpenMP with schedule implementation
      t0 = omp_get_wtime()
      call axpy_omp_schedule(n, alpha, x, 1, y, 1)
      t1 = omp_get_wtime()
      time_omp_schedule = time_omp_schedule + (t1 - t0)
   end do

   ! Average times over repetitions
   if (reps > 0) then
      time_blas = time_blas / real(reps, dp)
      time_omp = time_omp / real(reps, dp)
      time_omp_simd = time_omp_simd / real(reps, dp)
      time_omp_schedule = time_omp_schedule / real(reps, dp)
   end if

   ! Calculate metrics based on slide formulas
   bytes_per_elem = real(storage_size(x(1)) / 8, dp)
   
   ! AXPY: 2n FLOPs (1 multiplication + 1 addition per element)
   total_flops = 2.0_dp * real(n, dp)
   
   ! Data movement: 3n elements (read x, read y, write y) 
   bytes_moved = 3.0_dp * real(n, dp) * bytes_per_elem
   
   ! Operational intensity: 2n / (3n * 8 bytes) = 2 / 24 = 1/12 FLOPS/Byte
   ! This is independent of n for AXPY
   AI = 2.0_dp / (3.0_dp * bytes_per_elem)

   ! Maximum achievable performance (roofline model)
   ! Performance_max = min(Peak_GFLOPS, OI * Memory_BW)
   roof_bw = AI * bw_GBs  ! Memory-bound ceiling
   roofline = min(peak_GFLOPS, roof_bw)

   ! Performance (GFLOP/s) = (2n / t) / 10^9
   if (time_blas > 0.0_dp) then
      gflops_blas = (total_flops / time_blas) / 1.0e9_dp
      GBs_meas_blas = bytes_moved / (time_blas * 1.0e9_dp)
   else
      gflops_blas = 0.0_dp
      GBs_meas_blas = 0.0_dp
   end if

   if (time_omp > 0.0_dp) then
      gflops_omp = (total_flops / time_omp) / 1.0e9_dp
      GBs_meas_omp = bytes_moved / (time_omp * 1.0e9_dp)
   else
      gflops_omp = 0.0_dp
      GBs_meas_omp = 0.0_dp
   end if

   if (time_omp_simd > 0.0_dp) then
      gflops_omp_simd = (total_flops / time_omp_simd) / 1.0e9_dp
      GBs_meas_omp_simd = bytes_moved / (time_omp_simd * 1.0e9_dp)
   else
      gflops_omp_simd = 0.0_dp
      GBs_meas_omp_simd = 0.0_dp
   end if

   if (time_omp_schedule > 0.0_dp) then
      gflops_omp_schedule = (total_flops / time_omp_schedule) / 1.0e9_dp
      GBs_meas_omp_schedule = bytes_moved / (time_omp_schedule * 1.0e9_dp)
   else
      gflops_omp_schedule = 0.0_dp
      GBs_meas_omp_schedule = 0.0_dp
   end if

   ! Calculate percentage of maximum achievable performance
   if (roofline > 0.0_dp) then
      pct_blas = (gflops_blas / roofline) * 100.0_dp
      pct_omp = (gflops_omp / roofline) * 100.0_dp
      pct_omp_simd = (gflops_omp_simd / roofline) * 100.0_dp
      pct_omp_schedule = (gflops_omp_schedule / roofline) * 100.0_dp
   else
      pct_blas = 0.0_dp
      pct_omp = 0.0_dp
      pct_omp_simd = 0.0_dp
      pct_omp_schedule = 0.0_dp
   end if

   deallocate(x, y)

   ! Output results
   write(output_unit, *) '=================================================='
   write(output_unit, *) 'Roofline Model Analysis'
   write(output_unit, *) '=================================================='
   write(output_unit, *) 'System Parameters:'
   write(output_unit, *) '  Peak Performance: ', peak_GFLOPS, ' GFLOPS'
   write(output_unit, *) '  Memory Bandwidth: ', bw_GBs, ' GB/s'
   write(output_unit, *) ''
   write(output_unit, *) 'AXPY Characteristics:'
   write(output_unit, *) '  Operational Intensity: ', AI, ' FLOPS/Byte'
   write(output_unit, *) '  Memory-bound ceiling: ', roof_bw, ' GFLOPS'
   write(output_unit, *) '  Maximum achievable: ', roofline, ' GFLOPS'
   write(output_unit, *) ''
   write(output_unit, *) 'Measured Performance:'
   write(output_unit, '(A,F12.6,A,F12.2,A,F6.2,A)') '  BLAS:            ', gflops_blas, ' GFLOPS, ', &
                                                      GBs_meas_blas, ' GB/s (', pct_blas, '% of max)'
   write(output_unit, '(A,F12.6,A,F12.2,A,F6.2,A)') '  OpenMP:          ', gflops_omp, ' GFLOPS, ', &
                                                      GBs_meas_omp, ' GB/s (', pct_omp, '% of max)'
   write(output_unit, '(A,F12.6,A,F12.2,A,F6.2,A)') '  OpenMP SIMD:     ', gflops_omp_simd, ' GFLOPS, ', &
                                                      GBs_meas_omp_simd, ' GB/s (', pct_omp_simd, '% of max)'
   write(output_unit, '(A,F12.6,A,F12.2,A,F6.2,A)') '  OpenMP Schedule: ', gflops_omp_schedule, ' GFLOPS, ', &
                                                      GBs_meas_omp_schedule, ' GB/s (', pct_omp_schedule, '% of max)'
   write(output_unit, *) '=================================================='

   ! Write results to CSV file: OI, Perf Axpy1, Perf Axpy2, Perf Axpy3, Perf Axpy4
   write(10, '(F12.6,4(",",F12.6))') AI, gflops_blas, gflops_omp, gflops_omp_simd, gflops_omp_schedule

   close(10)

end program axpy_roofline
