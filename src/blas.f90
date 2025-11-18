module blas
    use iso_fortran_env, only: real32, real64
    implicit none

    private

    interface axpy
        module procedure daxpy_blas, saxpy_blas
    end interface axpy

    public :: axpy

    contains

subroutine daxpy_blas(alpha, x, y, incx, incy)
    use iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in) :: alpha
    real(real64), intent(in) :: x(:)
    real(real64), intent(inout) :: y(:)
    integer, intent(in), optional :: incx, incy
    ! Local variables
    integer :: incx_, incy_
    incx_ = 1
    incy_ = 1
    if (present(incx)) incx_ = incx
    if (present(incy)) incy_ = incy
    call daxpy(size(x),alpha,x,incx_,y,incy_)
end subroutine daxpy_blas

subroutine saxpy_blas(alpha, x, y, incx, incy)
    use iso_fortran_env, only: real32
    implicit none
    real(real32), intent(in) :: alpha
    real(real32), intent(in) :: x(:)
    real(real32), intent(inout) :: y(:)
    integer, intent(in), optional :: incx, incy
    ! Local variables
    integer :: incx_, incy_
    incx_ = 1
    incy_ = 1
    if (present(incx)) incx_ = incx
    if (present(incy)) incy_ = incy
    call saxpy(size(x),alpha,x,incx_,y,incy_)
end subroutine saxpy_blas

end module blas