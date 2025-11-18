program axpy_scaling
    use iso_fortran_env, only: dp => real64, output_unit, error_unit
    use omp_lib
    use blasomp
    use blas
    implicit none
    ! Inputs:
    integer :: reps, n
    character(len=20) :: n_str, reps_str
    character(len=100) :: filename
    ! Variables for timers
    real(dp) :: time_blas, time_omp, t0, t1
    ! Data arrays
    real(dp), allocatable :: x(:), y(:)
    real(dp) :: alpha
    ! Thread info
    integer :: nthreads
    ! Loop variable
    integer :: rep, i
    ! Info
    integer :: stat
    
    ! Read sizes and number of repetitions from command line arguments
    if (command_argument_count() < 2) then
        write(error_unit, *) 'Usage: axpy_scaling <problem_size> <repetitions>'
        stop
    end if
    call get_command_argument(1, n_str)
    call get_command_argument(2, reps_str)
    read(n_str, *) n
    read(reps_str, *) reps
    write(output_unit, *) 'Testing AXPY with n = ', n, ' for ', reps, ' repetitions'
    ! Get number of threads
    !$omp parallel
    !$omp single
    nthreads = omp_get_max_threads()
    write(output_unit, *) 'Using ', nthreads, ' OpenMP threads'
    !$omp end single
    !$omp end parallel

    ! Allocate and initialize data
    allocate(x(n), y(n),stat=stat)
    if (stat /= 0) then
        write(error_unit, *) 'Error allocating arrays of size ', n
        stop
    end if

    x = [(real(i, dp), i = 1, n)]
    y = 0.0_dp
    alpha = 2.0_dp  

    ! Benchmark BLAS AXPY
    time_blas = 0.0_dp
    do rep = 1, reps
        y = 0.0_dp
        t0 = omp_get_wtime()
        call daxpy(n, alpha, x, 1, y, 1)
        t1 = omp_get_wtime()
        time_blas = time_blas + (t1 - t0)
    end do

    ! Benchmark OpenMP AXPY
    time_omp = 0.0_dp
    do rep = 1, reps
        y = 0.0_dp
        t0 = omp_get_wtime()
        call axpy_omp(n, alpha, x, 1, y, 1)
        t1 = omp_get_wtime()
        time_omp = time_omp + (t1 - t0)
    end do

    ! Average times
    if (reps > 0) then
        time_blas = time_blas / real(reps, dp)
        time_omp = time_omp / real(reps, dp)
    end if

    ! Output results
    write(output_unit, *) 'Average time BLAS AXPY: ', time_blas, ' seconds'
    write(output_unit, *) 'Average time OpenMP AXPY: ', time_omp, ' seconds'

    ! Append to a csv file for the given problem size
    write(filename, '(A,I0,A)') 'axpy_scaling_results_', n, '.csv'
    open(unit=10, file=filename, status='unknown', action='write', position='append')
    write(10, '(I10,I10,1X,F15.6,F15.6)') n, nthreads, time_blas, time_omp
    close(10)

    ! Deallocate arrays
    deallocate(x, y, stat=stat)
    if (stat /= 0) then
        write(error_unit, *) 'Error deallocating arrays'
    end if


end program axpy_scaling