! eftdot.f90
! Error-free transformation for dot product
!
! Compile Options:
! -D_FMA for enabling FMA
!
! Reference: 
! Ogita, T., Rump, S. M., & Oishi, S. I. (2005). Accurate sum and dot product. 
! SIAM Journal on Scientific Computing, 26(6), 1955-1988.
!
module eftdot
    use, intrinsic :: iso_fortran_env, only: real32, real64
#ifdef _FMA
    use, intrinsic :: ieee_arithmetic, only: ieee_fma
#endif
    implicit none

    private

    integer, parameter :: sp = real32
    integer, parameter :: dp = real64
    
    public :: sp, dp
    public :: dot2, dotk

    !> Generic interface for dot2, single and double precision
    interface dot2
        module procedure dot2_s
        module procedure dot2_d
    end interface

    !> Generic interface for dotk, single and double precision
    interface dotk
        module procedure dotk_s
        module procedure dotk_d
    end interface

contains

    !> Dot2 for single precision
    pure subroutine dot2_s(x, y, n, res)
        integer, intent(in) :: n
        real(sp), intent(in) :: x(:), y(:)
        real(sp), intent(out) :: res

        real(sp) :: p, s, h, r, q
        real(sp) :: tmp_p, sum_val, z
        integer :: i

#ifndef _FMA
        ! The factor for Dekker's split (only needed if no FMA)
        real(sp) :: c, a1, a2, b1, b2
        ! 2^12 + 1 = 4097.0 for single precision
        real(sp), parameter :: factor_s = 4097.0_sp
#endif

        p = 0.0_sp
        s = 0.0_sp
        res = 0.0_sp

        ! Check if n is valid
        if (n <= 0) return

        do i = 1, n
            ! TwoProduct (h, r) = x * y
#ifdef _FMA
            ! FMA version
            h = x(i) * y(i)
            r = ieee_fma(x(i), y(i), -h)
#else
            ! No FMA: with DEKKER SPLIT
            h = x(i) * y(i)
            ! Split x(i)
            c = factor_s * x(i)
            a1 = c - (c - x(i))
            a2 = x(i) - a1
            ! Split y(i)
            c = factor_s * y(i)
            b1 = c - (c - y(i))
            b2 = y(i) - b1
            ! Error term
            r = a2 * b2 - (((h - a1 * b1) - a2 * b1) - a1 * b2)
#endif
            ! TwoSum (p, q) = p + h
            if (i == 1) then
                p = h
                s = r
            else
                tmp_p = p
                sum_val = tmp_p + h
                z = sum_val - tmp_p
                q = (tmp_p - (sum_val - z)) + (h - z)
                p = sum_val
                
                ! Accumulate errors
                s = s + (q + r)
            end if
        end do

        res = p + s
    end subroutine dot2_s

    !> Dot2 for double precision
    pure subroutine dot2_d(x, y, n, res)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(:), y(:)
        real(dp), intent(out) :: res
        real(dp) :: p, s, h, r, q
        real(dp) :: tmp_p, sum_val, z

#ifndef _FMA
        real(dp) :: c, a1, a2, b1, b2
        ! 2^27 + 1 = 134217729.0 for double precision
        real(dp), parameter :: factor_d = 134217729.0_dp
#endif
        integer :: i

        p = 0.0_dp
        s = 0.0_dp
        res = 0.0_dp

        if (n <= 0) return

        do i = 1, n
            ! TwoProduct (h, r) = x * y
#ifdef _FMA
            h = x(i) * y(i)
            r = ieee_fma(x(i), y(i), -h)
#else
            ! No FMA: with DEKKER SPLIT
            h = x(i) * y(i)
            ! Split x(i)
            c = factor_d * x(i)
            a1 = c - (c - x(i))
            a2 = x(i) - a1
            ! Split y(i)
            c = factor_d * y(i)
            b1 = c - (c - y(i))
            b2 = y(i) - b1
            ! Error term
            r = a2 * b2 - (((h - a1 * b1) - a2 * b1) - a1 * b2)
#endif
            ! TwoSum (p, q) = p + h
            if (i == 1) then
                p = h
                s = r
            else
                tmp_p = p
                sum_val = tmp_p + h
                z = sum_val - tmp_p
                q = (tmp_p - (sum_val - z)) + (h - z)
                p = sum_val
                ! Accumulate errors
                s = s + (q + r)
            end if
        end do

        res = p + s
    end subroutine dot2_d

    !> DotK for single precision. Use k=3 to simulate DP accuracy.
    subroutine dotk_s(x, y, n, k, res)
        integer, intent(in) :: n, k
        real(sp), intent(in) :: x(:), y(:)
        real(sp), intent(out) :: res
        real(sp), allocatable :: z(:)
        real(sp) :: h, r, tmp_h
        integer :: i, j

        if (n <= 0) then
            res = 0.0_sp; return
        end if

        allocate(z(2*n))

        ! 1. Generate 2n terms (h_i, r_i)
        do i = 1, n
            call two_product_s(x(i), y(i), z(i), z(i+n))
        end do

        ! 2. K-1 passes of VecSum (Algorithm 4.7)
        do j = 1, k - 1
            do i = 2, 2*n
                ! TwoSum(z_i, z_{i-1})
                tmp_h = z(i) + z(i-1)
                r = tmp_h - z(i)
                z(i-1) = (z(i) - (tmp_h - r)) + (z(i-1) - r)
                z(i) = tmp_h
            end do
        end do

        ! 3. Final summation
        res = z(2*n)
        do i = 1, 2*n - 1
            res = res + z(i)
        end do
        deallocate(z)
    end subroutine dotk_s

    !> DotK for double precision.
    subroutine dotk_d(x, y, n, k, res)
        integer, intent(in) :: n, k
        real(dp), intent(in) :: x(:), y(:)
        real(dp), intent(out) :: res
        real(dp), allocatable :: z(:)
        real(dp) :: h, r, tmp_h
        integer :: i, j

        if (n <= 0) then
            res = 0.0_dp; return
        end if

        allocate(z(2*n))

        do i = 1, n
            call two_product_d(x(i), y(i), z(i), z(i+n))
        end do

        do j = 1, k - 1
            do i = 2, 2*n
                tmp_h = z(i) + z(i-1)
                r = tmp_h - z(i)
                z(i-1) = (z(i) - (tmp_h - r)) + (z(i-1) - r)
                z(i) = tmp_h
            end do
        end do

        res = z(2*n)
        do i = 1, 2*n - 1
            res = res + z(i)
        end do
        deallocate(z)
    end subroutine dotk_d

    ! --- Private EFT Helpers ---

    pure subroutine two_product_s(a, b, x, y)
        real(sp), intent(in) :: a, b
        real(sp), intent(out) :: x, y
#ifndef _FMA
        real(sp) :: c, a1, a2, b1, b2
        real(sp), parameter :: factor_s = 4097.0_sp
#endif
        x = a * b
#ifdef _FMA
        y = ieee_fma(a, b, -x)
#else
        c = factor_s * a; a1 = c - (c - a); a2 = a - a1
        c = factor_s * b; b1 = c - (c - b); b2 = b - b1
        y = a2 * b2 - (((x - a1 * b1) - a2 * b1) - a1 * b2)
#endif
    end subroutine two_product_s

    pure subroutine two_product_d(a, b, x, y)
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: x, y
#ifndef _FMA
        real(dp) :: c, a1, a2, b1, b2
        real(dp), parameter :: factor_d = 134217729.0_dp
#endif
        x = a * b
#ifdef _FMA
        y = ieee_fma(a, b, -x)
#else
        c = factor_d * a; a1 = c - (c - a); a2 = a - a1
        c = factor_d * b; b1 = c - (c - b); b2 = b - b1
        y = a2 * b2 - (((x - a1 * b1) - a2 * b1) - a1 * b2)
#endif
    end subroutine two_product_d

end module eftdot