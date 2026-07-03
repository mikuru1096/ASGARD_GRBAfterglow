!f2py: skip
module adaptive_resampling_mod
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    private

    public :: adaptive_resampling_log

contains

    ! 对数空间自适应重采样：基于二阶导数曲率权重，将 n 点网格重采样为 m 点。
    ! Adaptive resampling in log space: use second-derivative curvature weights
    ! to reduce an n-point grid to m selected points.
    subroutine adaptive_resampling_log(x, f1, n, m, g, indices, n_resampled, info)

        integer, intent(in) :: n, m, g
        real(dp), intent(in), dimension(n) :: x,f1

        integer, intent(out), dimension(m) :: indices
        integer, intent(out) :: n_resampled,info

        real(dp), dimension(n) :: log_x
        real(dp), dimension(n-1) :: dlogx
        real(dp), dimension(n) :: df_dlogx,d2f_dlogx2,f
        real(dp), dimension(n) :: d2_abs,dx_adj,d2_weighted
        real(dp), dimension(n) :: d2_smooth,w,weighted_w,c
        real(dp), dimension(m) :: t
        real(dp) :: s,w_min,c_end
        integer, dimension(m) :: temp_indices
        integer :: i,k,idx

        info = 0
        n_resampled = m

        ! 检查 public 输入尺寸；这是系统边界，不是内部数值兜底。
        ! Check public input sizes; this is a system boundary, not an internal numerical fallback.
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

        ! 一阶导数：端点用单边差分，内部点用跨两侧的中心差分。
        ! First derivative: 1-sided differences at edges, centered across both neighbors inside.
        df_dlogx(1) = (f(2) - f(1)) / dlogx(1)  ! upwind
        df_dlogx(n) = (f(n) - f(n-1)) / dlogx(n-1)  ! backwind

        do i = 2, n-1
            df_dlogx(i) = (f(i+1) - f(i-1)) / (dlogx(i-1) + dlogx(i))
        end do

        ! 二阶导数：曲率越大，后续采样权重越高。
        ! Second derivative: larger curvature receives larger later sampling weight.
        d2f_dlogx2(1) = (df_dlogx(2) - df_dlogx(1)) / dlogx(1)
        d2f_dlogx2(n) = (df_dlogx(n) - df_dlogx(n-1)) / dlogx(n-1)

        do i = 2, n-1
            d2f_dlogx2(i) = (df_dlogx(i+1) - df_dlogx(i-1)) / (dlogx(i-1) + dlogx(i))
        end do

        d2_abs = abs(d2f_dlogx2)

        ! 用局部网格宽度修正曲率权重，避免宽网格被低估。
        ! Adjust curvature weights by local grid spacing so wide cells are not underestimated.
        dx_adj(1) = 1.0_dp / dlogx(1)
        dx_adj(n) = 1.0_dp / dlogx(n-1)

        do i = 2, n-1
            dx_adj(i) = 1.0_dp / ((dlogx(i-1) + dlogx(i)) / 2.0_dp)
        end do

        d2_weighted = d2_abs / dx_adj

        ! 平滑曲率权重，避免单点噪声独占采样预算。
        ! Smooth curvature weights so single-cell noise does not dominate the sampling budget.
        call moving_average(d2_weighted, n, 5, d2_smooth)

        s = sum(d2_smooth)

        ! 最小权重让低曲率区域仍然保留有限覆盖。
        ! Minimum weight keeps finite coverage in low-curvature regions.
        w_min = s / (m * g)

        ! 应用最小权重下限。
        ! Apply the minimum-weight floor.
        w = max(d2_smooth, w_min)

        ! 将点权重转成积分权重，匹配对数网格单元宽度。
        ! Convert point weights to integral weights matched to log-grid cell widths.
        weighted_w(1) = w(1) * dlogx(1)
        weighted_w(n) = w(n) * dlogx(n-1)

        do i = 2, n-1
            weighted_w(i) = w(i) * (dlogx(i-1) + dlogx(i)) / 2.0_dp
        end do

        ! 累计权重是反演采样坐标。
        ! Cumulative weight is the coordinate used for inverse sampling.
        c(1) = weighted_w(1)
        do i = 2, n
            c(i) = c(i-1) + weighted_w(i)
        end do

        c_end = c(n)

        ! 在累计权重坐标上均匀放置目标采样点。
        ! Place target samples uniformly in cumulative-weight coordinates.
        do k = 1, m
            t(k) = (k-1) * c_end / (m-1)
        end do

        ! 把目标累计权重映射回原始网格索引。
        ! Map target cumulative weights back to source-grid indices.
        indices = 0
        idx = 1
        do k = 1, m
            do while (idx <= n .and. c(idx) < t(k))
                idx = idx + 1
            end do

            if (idx > n) idx = n
            indices(k) = idx
        end do

        ! 去重并保持索引递增。
        ! Deduplicate indices and keep them ascending.
        call unique_sorted(indices, m, temp_indices, n_resampled)

        ! 如果曲率集中导致重复索引，就从缺失索引中均匀补足。
        ! If concentrated curvature creates duplicates, fill from missing indices uniformly.
        if (n_resampled < m) then
            call supplement_indices(temp_indices, n_resampled, n, m, indices)
            n_resampled = m
        else
            indices(1:m) = temp_indices(1:m)
        end if

    contains

        ! 滑动平均平滑：窗口半径为 window/2。
        ! Moving-average smoothing with radius window/2.
        subroutine moving_average(input, n, window, output)
            integer, intent(in) :: n, window
            real(dp), intent(in), dimension(n) :: input
            real(dp), intent(out), dimension(n) :: output
            integer :: i, j, w_start, w_end, count
            real(dp) :: window_sum

            do i = 1, n
                w_start = max(1, i - window/2)
                w_end = min(n, i + window/2)
                count = w_end - w_start + 1
                window_sum = 0.0_dp

                do j = w_start, w_end
                    window_sum = window_sum + input(j)
                end do

                output(i) = window_sum / count
            end do
        end subroutine moving_average

        ! 提取有序数组中的唯一值。
        ! Extract unique values from an already ordered array.
        subroutine unique_sorted(input, n_input, output, n_output)
            integer, intent(in) :: n_input
            integer, intent(in), dimension(n_input) :: input
            integer, intent(out) :: n_output
            integer, intent(out), dimension(n_input) :: output
            integer :: i

            n_output = 0
            output = 0

            do i = 1, n_input
                if (input(i) == 0) cycle

                if (n_output == 0 .or. output(n_output) /= input(i)) then
                    n_output = n_output + 1
                    output(n_output) = input(i)
                end if
            end do
        end subroutine unique_sorted

        ! 对很短的索引数组做原地整数排序。
        ! In-place integer sort for short index arrays.
        subroutine sort_integers(arr, n)
            integer, intent(in) :: n
            integer, intent(inout), dimension(n) :: arr
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

        ! 补充不足的索引点：从缺失索引中均匀选取，填充至目标数量 m_target。
        ! Supplement missing samples by selecting uniformly from absent indices until m_target is reached.
        subroutine supplement_indices(current_indices, n_current, n_total, m_target, new_indices)
            integer, intent(in) :: n_current, n_total, m_target
            integer, intent(in), dimension(n_current) :: current_indices
            integer, intent(out), dimension(m_target) :: new_indices
            integer :: i, j, k, missing_count, missing_needed, target_missing
            integer, dimension(n_total) :: missing_indices
            logical, dimension(n_total) :: present

            present = .false.
            do j = 1, n_current
                present(current_indices(j)) = .true.
            end do

            missing_count = 0
            do i = 1, n_total
                if (.not. present(i)) then
                    missing_count = missing_count + 1
                    missing_indices(missing_count) = i
                end if
            end do

            new_indices(1:n_current) = current_indices(1:n_current)
            missing_needed = m_target - n_current
            do k = 1, missing_needed
                target_missing = nint((k-1) * real(missing_count - 1, dp) / max(1, missing_needed - 1)) + 1
                new_indices(n_current + k) = missing_indices(target_missing)
            end do

            call sort_integers(new_indices, m_target)
        end subroutine supplement_indices

    end subroutine adaptive_resampling_log

end module adaptive_resampling_mod
