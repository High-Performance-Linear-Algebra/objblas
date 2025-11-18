program test_axpy
    use iso_fortran_env, only: real32, real64, output_unit
    use blas

    implicit none
    integer, parameter :: n = 10
    real(real32) :: x32(n),y32(n)
    real(real64) :: x64(n),y64(n)
    real(real32) :: alpha32
    real(real64) :: alpha64
    

    x32 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y32 = 0.0_real32
    alpha32 = 2.0_real32

    call axpy(alpha32, x32, y32)
    write(output_unit,*) "Single Precision AXPY Result:"
    write(output_unit,*) y32

    x64 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y64 = 0.0_real64
    alpha64 = 2.0_real64

    call axpy(alpha64, x64, y64)
    write(output_unit,*) "Double Precision AXPY Result:"
    write(output_unit,*) y64

end program test_axpy