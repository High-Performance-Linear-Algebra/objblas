module blasomp
  use iso_fortran_env, only: int32, real64
  use omp_lib
  implicit none

  private
  public :: axpy_omp, axpy_omp_simd, axpy_omp_schedule

contains

    subroutine axpy_omp(n, a, x, incx, y, incy)
        integer(int32), intent(in) :: n, incx, incy
        real(real64), intent(in) :: a
        real(real64), intent(in) :: x(*)
        real(real64), intent(inout) :: y(*)
        integer(int32) :: i

        !$omp parallel do private(i) shared(n, a, x, incx, y, incy)
        do i = 1, n
            y(1 + (i - 1) * incy) = y(1 + (i - 1) * incy) + a * x(1 + (i - 1) * incx)
        end do
        !$omp end parallel do

    end subroutine axpy_omp

    subroutine axpy_omp_simd(n, a, x, incx, y, incy)
        integer(int32), intent(in) :: n, incx, incy
        real(real64), intent(in) :: a
        real(real64), intent(in) :: x(*)
        real(real64), intent(inout) :: y(*)
        integer(int32) :: i

        !$omp parallel do simd private(i) shared(n, a, x, incx, y, incy)
        do i = 1, n
            y(1 + (i - 1) * incy) = y(1 + (i - 1) * incy) + a * x(1 + (i - 1) * incx)
        end do
        !$omp end parallel do simd

    end subroutine axpy_omp_simd

    subroutine axpy_omp_schedule(n, a, x, incx, y, incy)
        integer(int32), intent(in) :: n, incx, incy
        real(real64), intent(in) :: a
        real(real64), intent(in) :: x(*)
        real(real64), intent(inout) :: y(*)
        integer(int32) :: i

        !$omp parallel do private(i) shared(n, a, x, incx, y, incy) schedule(static)
        do i = 1, n
            y(1 + (i - 1) * incy) = y(1 + (i - 1) * incy) + a * x(1 + (i - 1) * incx)
        end do
        !$omp end parallel do

    end subroutine axpy_omp_schedule

end module blasomp