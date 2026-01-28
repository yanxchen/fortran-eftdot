program test_eftdot
    use eftdot
    implicit none

    integer, parameter :: N = 100
    integer, parameter :: SAMPLES = 1000
    
    ! Buffers
    real(sp) :: x_s(N), y_s(N), res_dot3_s, res_int_s, res_loop_s
    real(dp) :: x_d(N), y_d(N), res_dot3_d, res_int_d, res_loop_d
    real(dp) :: res_dot3_s_d
    real(dp) :: accurate_ref, cond_val
    integer :: i, s, ios

    open(unit=10, file='input_data.bin', access='stream', form='unformatted', status='old')
    open(unit=20, file='results.csv', status='replace')
    
    write(20, '(A)') "cond,accurate,dot3_s,int_s,loop_s,dot3_d,int_d,loop_d"

    do s = 1, SAMPLES
        ! Read metadata (DP), then SP vectors, then DP vectors in sequence
        read(10, iostat=ios) cond_val, accurate_ref, x_s, y_s, x_d, y_d
        if (ios /= 0) exit

        ! --- Single Precision Logic ---
        call dotk(x_s, y_s, N, 3, res_dot3_s)
        res_int_s = dot_product(x_s, y_s)
        res_loop_s = 0.0_sp
        do i = 1, N
            res_loop_s = res_loop_s + (x_s(i) * y_s(i))
        end do

        ! --- Double Precision Logic ---
        call dotk(x_d, y_d, N, 3, res_dot3_d)
        res_int_d = dot_product(x_d, y_d)
        res_loop_d = 0.0_dp
        do i = 1, N
            res_loop_d = res_loop_d + (x_d(i) * y_d(i))
        end do

        ! Export using high-precision format without trailing comma
        write(20, '(7(E24.16, ","), E24.16)') &
            cond_val, accurate_ref, &
            real(res_dot3_s, dp), real(res_int_s, dp), real(res_loop_s, dp), &
            real(res_dot3_d, dp), real(res_int_d, dp), real(res_loop_d, dp)
    end do

    close(10); close(20)
end program