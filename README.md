# objblas: Roofline Performance Analysis of AXPY Implementations

This repository contains Fortran implementations and performance analysis of the AXPY operation (`y ← αx + y`) using different approaches:
- Standard BLAS
- OpenMP parallel
- OpenMP with SIMD
- OpenMP with custom scheduling

## Project Structure

- `src/` — Fortran source files for BLAS and OpenMP implementations
- `test/axpy_perf.f90` — Main performance test program
- `test/axpy_perf.f90` — Measures and compares the performance of all AXPY implementations, writing results to `results.csv`
- `results.csv` — Output file with measured performance data

## How It Works

The test program:
1. Prompts for system memory bandwidth (GB/s), peak GFLOPS, problem size `n`, and number of repetitions
2. Allocates and initializes vectors `x` and `y`
3. Runs each AXPY implementation, measuring execution time
4. Computes:
   - Operational Intensity (OI): $\frac{1}{12}$ FLOPS/Byte (for double precision)
   - Measured GFLOPS for each implementation
   - Theoretical roofline (max achievable GFLOPS): $\min$(Peak GFLOPS, OI × Memory Bandwidth)
   - Percentage of roofline achieved by each implementation
5. Outputs results to the terminal and to `results.csv`

## Roofline Model
- **AXPY FLOPs:** $2n$ (1 multiply + 1 add per element)
- **Data movement:** $3n$ (read x, read y, write y)
- **Operational Intensity:** $\frac{2n}{3n \times 8}$ = $\frac{1}{12}$ FLOPS/Byte
- **Memory-bound ceiling:** OI × Memory Bandwidth
- **Maximum achievable performance:** $\min$(Peak GFLOPS, OI × Memory Bandwidth)

## Example Output
```
 Enter memory bandwidth in GB/s: 
89.6
 Enter peak performance in GFLOPS: 
844
 Enter problem size (n): 
100000
 Enter number of repetitions: 
50
 Testing AXPY with n =       100000
 ==================================================
 Roofline Model Analysis
 ==================================================
 System Parameters:
   Peak Performance:    844.00000000000000       GFLOPS
   Memory Bandwidth:    89.599999999999994       GB/s
 
 AXPY Characteristics:
   Operational Intensity:    8.3333333333333329E-002  FLOPS/Byte
   Memory-bound ceiling:    7.4666666666666659       GFLOPS
   Maximum achievable:    7.4666666666666659       GFLOPS
 
 Measured Performance:
  BLAS:                7.133314 GFLOPS,        85.60 GB/s ( 95.54% of max)
  OpenMP:              7.303652 GFLOPS,        87.64 GB/s ( 97.82% of max)
  OpenMP SIMD:         7.221572 GFLOPS,        86.66 GB/s ( 96.72% of max)
  OpenMP Schedule:     7.406552 GFLOPS,        88.88 GB/s ( 99.19% of max)
 ==================================================
```

## Usage
1. Build the project using CMake or your preferred Fortran build system
2. Run the `axpy_perf` executable
3. Enter the required parameters when prompted
4. Analyze the output in the terminal and in `results.csv`

## Requirements
- Fortran compiler with OpenMP support
- BLAS library
- CMake (for building)

## Exercise

You can play around with the implementation in the `blasOMP.f90` module to see
if you can obtain better measurements for your axpy.