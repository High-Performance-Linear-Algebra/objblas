program roofline_demo
   use iso_fortran_env, only: wp => real64, output_unit
   implicit none
   integer :: i, j, k, N, nb, ntests, variant, method
   integer :: seed(4)
   real(wp), allocatable :: A(:,:), B(:,:), C(:,:)
   real(wp) :: t0, t1, elapsed
   real(wp) :: peak_gflops_default, bw_gb_s_default
   real(wp) :: peak_gflops, bw_gb_s
   real(wp) :: gflops_achieved, gflops_roofline
   real(wp) :: flops, bytes_min, intensity
   real(wp) :: knee, step
   integer :: repeat, r, warm_up_iters
   real(wp) :: repeat_real
   character(len=200) :: line
   integer, parameter :: max_tests = 24
   integer, allocatable :: sizes(:)
   logical :: use_blocking
   character(len=200) :: csv_filename
   integer :: csv_unit

   ! defaults (user should replace with their machine's measured peak floats and bandwidth)
   peak_gflops_default = 500.0_wp     ! GFLOPS (example)
   bw_gb_s_default      = 50.0_wp     ! memory bandwidth in GB/s (example)

   write(output_unit,*) "Roofline demo for dense matrix-matrix multiply (DGEMM)"
   write(output_unit,*) "-----------------------------------------------------"
   write(output_unit,'(A)',advance='no') "Enter peak FLOPS [GFLOPS] (default "
   write(output_unit,'(F0.1)',advance='no') peak_gflops_default
   write(output_unit,'(A)',advance='no') "): "
   read(*,'(A)') line
   if (len_trim(line) == 0) then
       peak_gflops = peak_gflops_default
   else
       read(line,*) peak_gflops
   end if

   write(output_unit,'(A)',advance='no') "Enter memory bandwidth [GB/s] (default "
   write(output_unit,'(F0.1)',advance='no') bw_gb_s_default
   write(output_unit,'(A)',advance='no') "): "
   read(*,'(A)') line
   if (len_trim(line) == 0) then
       bw_gb_s = bw_gb_s_default
   else
       read(line,*) bw_gb_s
   end if

   write(output_unit,*)
   write(output_unit,*) "Select matrix multiply method:"
   write(output_unit,*) "1 = naive triple-loop"
   write(output_unit,*) "2 = BLAS DGEMM"
   write(output_unit,*) "3 = blocked matmul"
   write(output_unit,'(A)',advance='no') "Enter choice (default 1): "
   read(*,'(A)') line
   if (len_trim(line) == 0) then
      method = 1
   else
      read(line,*) method
   end if

   if (method == 3) then
      write(output_unit,*)
      write(output_unit,'(A)',advance='no') "Try blocked algorithm? (y/N): "
      read(*,'(A)') line
      if (index('yY', line(1:1)) /= 0) then
          use_blocking = .true.
          write(output_unit,'(A)',advance='no') "Enter block size (typical 32..256, default 64): "
          read(*,'(A)') line
          if (len_trim(line) == 0) then
               nb = 64
          else
               read(line,*) nb
          end if
      else
          use_blocking = .false.
          nb = 0
      end if
   end if

    ! Calculate the knee point for reference
    knee = 48.0_wp * peak_gflops / bw_gb_s
    if (knee < 64.0_wp) knee = 64.0_wp

    ! Use N and variant as temporary variables for min and max sizes
    N = int(0.5*knee)        ! Default min size
    variant = int(2.0*knee)  ! Default max size

    write(output_unit,'(A,I0,A)',advance='no') "Enter minimum matrix size (default ", N, "): "
    read(*,'(A)') line
    if (len_trim(line) /= 0) then
        read(line,*) N
        ! Ensure minimum is at least 64
        if (N < 64) N = 64
    end if

    write(output_unit,'(A,I0,A)',advance='no') "Enter maximum matrix size (default ", variant, "): "
    read(*,'(A)') line
    if (len_trim(line) /= 0) then
        read(line,*) variant
        ! Ensure maximum is at least minimum
        if (variant < N) variant = N
    end if

   ! Size stepping
    write(output_unit,*)
    write(output_unit,'(A)',advance='no') "Enter number of matrix sizes to test (1..24, default 6): "
    read(*,'(A)') line
    if (len_trim(line) == 0) then
        ntests = 6
    else
        read(line,*) ntests
        if (ntests < 1) ntests = 1
        if (ntests > max_tests) ntests = max_tests
    end if

    allocate(sizes(ntests))
    sizes = 0

    ! Generate evenly spaced sizes from min to max
    if (ntests == 1) then
        sizes(1) = (N / 64) * 64  ! Round down to multiple of 64
    else
        step = real(variant - N) / real(ntests - 1, wp)
        do i = 1, ntests
            sizes(i) = int(N + step*real(i-1,wp))
            ! Round sizes to multiple of 64 for block-friendliness
            sizes(i) = (sizes(i) / 64) * 64  ! Round down to multiple of 64
        end do
    end if

   write(output_unit,*)
   write(output_unit,*) "Auto-selected matrix sizes based on peak & bandwidth:"
   do i = 1, ntests
       write(output_unit,'(A,I6)') "  N = ", sizes(i)
   end do
   
   write(output_unit,*)
   write(output_unit,*) "Using peak GFLOPS = ", peak_gflops, " GFLOPS"
   write(output_unit,*) "Using memory bandwidth = ", bw_gb_s, " GB/s"
   if (method == 3 .and. use_blocking) then
       write(output_unit,*) "Using blocked matmul with block size =", nb
   end if
   ! File for the output
   write(output_unit,'(A)',advance='no') "Enter CSV output filename (default: roofline.csv): "
   read(*,'(A)') csv_filename
   if (len_trim(csv_filename) == 0) csv_filename = "roofline.csv"
   ! Open it
   ! Assign a file unit (any unused integer)
   csv_unit = 10
   open(unit=csv_unit, file=csv_filename, status='replace', action='write', form='formatted')
   ! Optional: write CSV header
   write(csv_unit,'(A)') "N,elapsed_s,GFLOPS_achieved,intensity,roofline_bound_GFLOPS"

   write(output_unit,*) "-----------------------------------------------"
   write(output_unit,*) "Columns: N, time(s), achieved(GFLOPS), intensity(flops/byte), roofline_bound(GFLOPS)"
   write(output_unit,*) "-----------------------------------------------"

   ! Number of warm-up iterations
   warm_up_iters = 5

   do i = 1, ntests
       N = sizes(i)
       if (N <= 0) cycle

       ! allocate matrices
       allocate(A(N,N), B(N,N), C(N,N))
       call init_matrices(A,B,C,N)

       ! Compute repeat dynamically per N
       repeat_real = 1.0E9_wp / (real(N, wp)**3)  ! Increased factor for more accurate timing
       repeat = max(1, int(repeat_real))
       if (repeat > 20) repeat = 20  ! Cap maximum repeats

       ! warm-up (multiple iterations to get caches ready and CPU at full speed)
       do r = 1, warm_up_iters
            select case(method)
            case(1)
                call matmul_naive(A,B,C,N)
            case(2)
                call matmul_blas(N, A, B, C, 1.0d0, 0.0d0)
            case(3)
                call matmul_blocked(A,B,C,N,nb)
            end select
       end do

       ! Reset C matrix before timed run
       C = 0.0_wp

       ! Force synchronization before timing
       call flush()
       
       ! Use wall-clock time instead of CPU time for more accurate measurement
       call system_clock(count=k)  ! k is used as a temporary variable
       t0 = real(k, wp)
       
       do r = 1, repeat
            select case(method)
            case(1)
                call matmul_naive(A,B,C,N)
            case(2)
                call matmul_blas(N, A, B, C, 1.0d0, 0.0d0)
            case(3)
                call matmul_blocked(A,B,C,N,nb)
            end select
            ! Reset C for next iteration to avoid cache effects
            if (r < repeat) C = 0.0_wp
       end do
       
       call system_clock(count=k, count_rate=j)  ! j stores count_rate
       t1 = real(k, wp)
       elapsed = max(1.0e-12_wp, (t1 - t0)/real(j, wp)/repeat)

       ! Compute GFLOPS
       flops = 2.0_wp * real(N, wp)**3
       bytes_min = 3.0_wp * real(N, wp)**2 * 8.0_wp
       intensity = flops / bytes_min
       gflops_achieved = flops / elapsed / 1.0e9_wp
       gflops_roofline = min(peak_gflops, bw_gb_s * intensity)

       write(*,'(I6,2X,F10.4,2X,F10.3,2X,F10.4,2X,F10.3)') N, elapsed, gflops_achieved, intensity, gflops_roofline
       write(csv_unit,'(I0,A,F0.6,A,F0.6,A,F0.6,A,F0.6)') N, ",", elapsed, ",", gflops_achieved, ",", & 
         & intensity, ",", gflops_roofline

       deallocate(A,B,C)
   end do

   close(csv_unit)

contains

   subroutine init_matrices(A,B,C,N)
      integer, intent(in) :: N
      real(wp), intent(out) :: A(N,N), B(N,N), C(N,N)
      integer :: i,j
      ! Initialize with values that prevent compiler optimization
      do j = 1, N
          do i = 1, N
               A(i,j) = real(mod(i + j, 100) + 1, wp) * 1.0d-3
               B(i,j) = real(mod(i*j + 1, 100) + 1, wp) * 1.0d-3
               C(i,j) = 0.0_wp
          end do
      end do
   end subroutine init_matrices

   subroutine matmul_naive(A,B,C,N)
      integer, intent(in) :: N
      real(wp), intent(in) :: A(N,N), B(N,N)
      real(wp), intent(out) :: C(N,N)
      integer :: i,j,k
      real(wp) :: sum
      ! Optimized loop order for better cache utilization
      do j = 1, N
          do i = 1, N
               sum = 0.0_wp
               do k = 1, N
                   sum = sum + A(i,k) * B(k,j)
               end do
               C(i,j) = sum
          end do
      end do
   end subroutine matmul_naive

   subroutine matmul_blocked(A,B,C,N,nb)
      integer, intent(in) :: N, nb
      real(wp), intent(in) :: A(N,N), B(N,N)
      real(wp), intent(out) :: C(N,N)
      integer :: i,j,k,ii,jj,kk
      integer :: ib, jb, kb, iend, jend, kend
      real(wp) :: tmp
      
      ! Initialize C to zero
      C = 0.0_wp
      
      ! Improved blocked matmul with better cache locality
      do jb = 1, N, nb
          jend = min(N, jb + nb - 1)
          do kb = 1, N, nb
               kend = min(N, kb + nb - 1)
               do ib = 1, N, nb
                   iend = min(N, ib + nb - 1)
                   ! Block multiply with optimized inner loops
                   do j = jb, jend
                        do i = ib, iend
                            tmp = C(i,j)
                            do k = kb, kend
                                 tmp = tmp + A(i,k) * B(k,j)
                            end do
                            C(i,j) = tmp
                        end do
                   end do
               end do
          end do
      end do
   end subroutine matmul_blocked

   subroutine matmul_blas(N, A, B, C, alpha, beta)
      implicit none
      integer, intent(in) :: N
      real(8), intent(in) :: A(N,N), B(N,N)
      real(8), intent(out) :: C(N,N)
      real(8), intent(in) :: alpha, beta

      ! BLAS DGEMM interface
      external dgemm
      character(1) :: transa, transb
      integer :: lda, ldb, ldc

      transa = 'N'   ! No transpose
      transb = 'N'   ! No transpose
      lda = N
      ldb = N
      ldc = N

      call dgemm(transa, transb, N, N, N, alpha, A, lda, B, ldb, beta, C, ldc)
   end subroutine matmul_blas

end program roofline_demo
