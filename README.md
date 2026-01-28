# EFT Dot Fortran

A Fortran implementation of Error-Free Transformation (EFT) based dot products for high-precision computation.

## Features

- Support for both single precision (`real32`) and double precision (`real64`)
- Optional FMA (Fused Multiply-Add) support via compile flag

## Files

- `src/eftdot.f90`: Main EFT dot product module
- `tests/test.f90`: Test program for accuracy benchmarking
- `tests/gen_dot_data.py`: Generates ill-conditioned test vectors using `accupy` (accupy: https://github.com/sigma-py/accupy)
- `tests/plot_res.py`: Plots accuracy comparison results

## Usage

### Build Library
```bash
# Compile library and module
$ make

# Enable FMA optimization (uncomment in Makefile)
# FFLAGS += -D_FMA -mfma
```

### Run Tests
```bash
# Run full accuracy benchmark (generates data, runs tests, plots results)
$ make tests

# Clean build artifacts
$ make clean
```

## Non-installation Usage

Codes in `eftdot.f90` can be directly copied and used anywhere.

### Example: Using `dot2`
```Fortran
use eftdot
integer, parameter :: n = 3
real(sp) :: x(n), y(n), result

x = [1.0_sp, 2.0_sp, 3.0_sp]
y = [4.0_sp, 5.0_sp, 6.0_sp]

! Compensated dot product
call dot2(x, y, n, result)
```

### Example: Using `dotk` for higher accuracy
```Fortran
use eftdot
integer, parameter :: n = 100
real(sp) :: x(n), y(n), result

! k=3 provides approximately double precision accuracy for SP inputs
call dotk(x, y, n, 3, result)
```

## API Reference

### `dot2(x, y, n, res)`
Compensated dot product using TwoProduct and TwoSum error-free transformations.
- `x, y`: Input vectors (single or double precision)
- `n`: Vector length
- `res`: Result with improved accuracy

### `dotk(x, y, n, k, res)`
K-fold compensated dot product for higher accuracy.
- `x, y`: Input vectors
- `n`: Vector length
- `k`: Number of compensation passes (k=2 is equivalent to dot2)
- `res`: High-accuracy result

## Reference

Ogita, T., Rump, S. M., & Oishi, S. I. (2005). *Accurate sum and dot product.* 
SIAM Journal on Scientific Computing, 26(6), 1955-1988.
