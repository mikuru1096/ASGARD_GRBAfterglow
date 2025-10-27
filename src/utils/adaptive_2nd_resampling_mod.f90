module adaptive_resampling_mod
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    
contains
    
    subroutine adaptive_resampling_log(x, f1, n, m, g, indices, n_resampled, info)

        integer, intent(in) :: n, m, g
        real(dp), intent(in) :: x(n), f1(n)
        
        integer, intent(out) :: indices(m), n_resampled, info
        
        real(dp) :: log_x(n), dlogx(n-1)
        real(dp) :: df_dlogx(n), d2f_dlogx2(n),f(n)
        real(dp) :: d2_abs(n), dx_adj(n), d2_weighted(n)
        real(dp) :: d2_smooth(n), w(n), weighted_w(n), c(n)
        real(dp) :: s, w_min, t(m), c_end
        integer :: i, k, idx, temp_indices(m)
        
        info = 0
        n_resampled = m
        
        ! check inputs
        if (n <= 0 .or. m <= 0 .or. g <= 0) then
            info = -1
            return
        end if
        
        if (m > n) then
            info = -2
            return
        end if
        
        log_x = log10(x)
        f = log10(f1)
        
        do i = 1, n-1
            dlogx(i) = log_x(i+1) - log_x(i)
        end do
        
        ! 1st-order derivative
        df_dlogx(1) = (f(2) - f(1)) / dlogx(1)  ! upwind
        df_dlogx(n) = (f(n) - f(n-1)) / dlogx(n-1)  ! backwind
        
        do i = 2, n-1
            df_dlogx(i) = (f(i+1) - f(i-1)) / (dlogx(i-1) + dlogx(i))
        end do
        
        ! 2nd-order derivative
        d2f_dlogx2(1) = (df_dlogx(2) - df_dlogx(1)) / dlogx(1)  
        d2f_dlogx2(n) = (df_dlogx(n) - df_dlogx(n-1)) / dlogx(n-1)
        
        do i = 2, n-1
            d2f_dlogx2(i) = (df_dlogx(i+1) - df_dlogx(i-1)) / (dlogx(i-1) + dlogx(i))
        end do
        
        d2_abs = abs(d2f_dlogx2)
        
        ! Consider adjusting the weight of the grid spacing
        dx_adj(1) = 1.0_dp / dlogx(1)
        dx_adj(n) = 1.0_dp / dlogx(n-1)
        
        do i = 2, n-1
            dx_adj(i) = 1.0_dp / ((dlogx(i-1) + dlogx(i)) / 2.0_dp)
        end do
        
        d2_weighted = d2_abs / dx_adj
        
        ! Smooth weight (moving)
        call moving_average(d2_weighted, n, 5, d2_smooth)
        
        s = sum(d2_smooth)
        
        ! Minimum weight
        w_min = s / (m * g)
        
        ! Adjust the weights
        w = max(d2_smooth, w_min)
        
        ! Consider the weighted weights of the grid spacing
        weighted_w(1) = w(1) * dlogx(1)
        weighted_w(n) = w(n) * dlogx(n-1)
        
        do i = 2, n-1
            weighted_w(i) = w(i) * (dlogx(i-1) + dlogx(i)) / 2.0_dp
        end do
        
        ! Cumulative weight
        c(1) = weighted_w(1)
        do i = 2, n
            c(i) = c(i-1) + weighted_w(i)
        end do
        
        c_end = c(n)
        
        ! Generate equally spaced target values
        do k = 1, m
            t(k) = (k-1) * c_end / (m-1)
        end do
        
        ! find indices
        indices = 0
        do k = 1, m
            idx = 1
            do while (idx <= n .and. c(idx) < t(k))
                idx = idx + 1
            end do
            
            if (idx > n) idx = n
            indices(k) = idx
        end do
        
        ! Ensure index is unique and ascending order
        call unique_sorted(indices, m, temp_indices, n_resampled)
        
        ! If index insufficient points due to duplicate indexes, add more points
        if (n_resampled < m) then
            call supplement_indices(temp_indices, n_resampled, n, m, indices)
            n_resampled = m
        else
            indices(1:m) = temp_indices(1:m)
        end if
        
    contains
    
        subroutine moving_average(input, n, window, output)
            integer, intent(in) :: n, window
            real(dp), intent(in) :: input(n)
            real(dp), intent(out) :: output(n)
            integer :: i, j, w_start, w_end, count
            real(dp) :: sum
            
            do i = 1, n
                w_start = max(1, i - window/2)
                w_end = min(n, i + window/2)
                count = w_end - w_start + 1
                sum = 0.0_dp
                
                do j = w_start, w_end
                    sum = sum + input(j)
                end do
                
                output(i) = sum / count
            end do
        end subroutine moving_average
        
        subroutine unique_sorted(input, n_input, output, n_output)
            integer, intent(in) :: n_input
            integer, intent(in) :: input(n_input)
            integer, intent(out) :: n_output
            integer, intent(out) :: output(n_input)
            integer :: i, j
            logical :: duplicate
            
            n_output = 0
            output = 0
            
            do i = 1, n_input
                if (input(i) == 0) cycle
                
                duplicate = .false.
                do j = 1, n_output
                    if (output(j) == input(i)) then
                        duplicate = .true.
                        exit
                    end if
                end do
                
                if (.not. duplicate) then
                    n_output = n_output + 1
                    output(n_output) = input(i)
                end if
            end do
            
            call sort_integers(output, n_output)
        end subroutine unique_sorted
        
        ! bubble sort
        subroutine sort_integers(arr, n)
            integer, intent(in) :: n
            integer, intent(inout) :: arr(n)
            integer :: i, j, temp
            
            do i = 1, n-1
                do j = i+1, n
                    if (arr(i) > arr(j)) then
                        temp = arr(i)
                        arr(i) = arr(j)
                        arr(j) = temp
                    end if
                end do
            end do
        end subroutine sort_integers
        
        subroutine supplement_indices(current_indices, n_current, n_total, m_target, new_indices)
            integer, intent(in) :: n_current, n_total, m_target
            integer, intent(in) :: current_indices(n_current)
            integer, intent(out) :: new_indices(m_target)
            integer :: i, j, k, missing_count
            integer, allocatable :: missing_indices(:)
            logical, allocatable :: present(:)
            
            ! 使用逻辑数组标记哪些索引已经存在
            allocate(present(n_total))
            present = .false.
            
            do j = 1, n_current
                if (current_indices(j) >= 1 .and. current_indices(j) <= n_total) then
                    present(current_indices(j)) = .true.
                end if
            end do
            
            ! 计算缺失的索引数量	
            missing_count = 0
            do i = 1, n_total
                if (.not. present(i)) then
                    missing_count = missing_count + 1
                end if
            end do
            
            ! 分配内存并收集缺失的索引
            if (missing_count > 0) then
                allocate(missing_indices(missing_count))
                missing_count = 0
                do i = 1, n_total
                    if (.not. present(i)) then
                        missing_count = missing_count + 1
                        missing_indices(missing_count) = i
                    end if
                end do
            else
                ! 如果没有缺失索引，创建一个空数组
                allocate(missing_indices(1))
                missing_count = 0
            end if
            
            ! 补充缺失的索引
            new_indices(1:n_current) = current_indices(1:n_current)
            
            if (missing_count > 0) then
                do k = 1, m_target - n_current
                    i = nint((k-1) * real(missing_count, dp) / max(1, m_target - n_current - 1)) + 1
                    if (i <= missing_count) then
                        new_indices(n_current + k) = missing_indices(i)
                    else
                        new_indices(n_current + k) = n_total
                    end if
                end do
            else
                ! 如果没有缺失索引，均匀补充
                do k = 1, m_target - n_current
                    i = nint((k-1) * real(n_total, dp) / max(1, m_target - n_current - 1)) + 1
                    new_indices(n_current + k) = min(i, n_total)
                end do
            end if
            
            ! 排序新索引
            call sort_integers(new_indices, m_target)
            
            ! 释放内存
            if (allocated(present)) deallocate(present)
            if (allocated(missing_indices)) deallocate(missing_indices)
        end subroutine supplement_indices
        
    end subroutine adaptive_resampling_log

end module adaptive_resampling_mod
