!f2py: skip
module interpolation_common
    use constants
    implicit none
    private

    public :: accum_logsed, accum_shifted

contains

! 在对数-线性空间中插值并累加 SED 到观测网格。
! Interpolate in log-linear space and accumulate the SED on the observer grid.
subroutine accum_logsed(src_x, src_y, num_src, dst_x, num_dst, accum)
    implicit none
    integer, intent(in) :: num_src, num_dst
    real(8), intent(in), dimension(num_src) :: src_x,src_y
    real(8), intent(in), dimension(num_dst) :: dst_x
    real(8), intent(inout), dimension(num_dst) :: accum
    integer :: i_dst, i_src
    real(8) :: ratio, y_lo, y_hi

    i_src = 1
    do i_dst = 1, num_dst
        if (dst_x(i_dst) < src_x(1)) cycle
        if (dst_x(i_dst) > src_x(num_src)) exit
        do while (i_src < num_src - 1 .and. dst_x(i_dst) > src_x(i_src + 1))
            i_src = i_src + 1
        end do
        if (dst_x(i_dst) >= src_x(i_src) .and. dst_x(i_dst) <= src_x(i_src + 1)) then
            ratio = (dst_x(i_dst) - src_x(i_src)) / (src_x(i_src + 1) - src_x(i_src))
            if (src_y(i_src) > -huge(1d0)/2d0) then
                y_lo = exp(src_y(i_src))
            else
                y_lo = 0d0
            end if
            if (src_y(i_src + 1) > -huge(1d0)/2d0) then
                y_hi = exp(src_y(i_src + 1))
            else
                y_hi = 0d0
            end if
            if (y_lo > 0d0 .and. y_hi > 0d0) then
                accum(i_dst) = accum(i_dst) + exp(src_y(i_src) + ratio * (src_y(i_src + 1) - src_y(i_src)))
            else
                accum(i_dst) = accum(i_dst) + (1d0-ratio)*y_lo + ratio*y_hi
            end if
        end if
    end do
end subroutine accum_logsed

! 源频率整体平移时，从线性 SED 对目标频率插值并累加。
! For a global source-frequency shift, interpolate linear SED values and accumulate.
subroutine accum_shifted(src_x, src_y, num_src, dst_x, num_dst, log_shift, log_weight, accum)
    implicit none
    integer, intent(in) :: num_src, num_dst
    real(8), intent(in), dimension(num_src) :: src_x,src_y
    real(8), intent(in), dimension(num_dst) :: dst_x
    real(8), intent(in) :: log_shift,log_weight
    real(8), intent(inout), dimension(num_dst) :: accum
    integer :: i_dst, i_src
    real(8) :: target_x, ratio, y_lo, y_hi, y_interp

    i_src = 1
    do i_dst = 1, num_dst
        target_x = dst_x(i_dst) + log_shift
        if (target_x < src_x(1)) cycle
        if (target_x > src_x(num_src)) exit
        do while (i_src < num_src - 1 .and. target_x > src_x(i_src + 1))
            i_src = i_src + 1
        end do
        if (target_x >= src_x(i_src) .and. target_x <= src_x(i_src + 1)) then
            ratio = (target_x - src_x(i_src)) / (src_x(i_src + 1) - src_x(i_src))
            y_lo = src_y(i_src)
            y_hi = src_y(i_src + 1)
            if (y_lo > 0d0 .and. y_hi > 0d0) then
                y_interp = dexp(dlog(y_lo) + ratio * (dlog(y_hi) - dlog(y_lo)))
            else
                y_interp = (1d0-ratio)*y_lo + ratio*y_hi
            end if
            accum(i_dst) = accum(i_dst) + y_interp*dexp(log_weight)
        end if
    end do
end subroutine accum_shifted

end module interpolation_common
