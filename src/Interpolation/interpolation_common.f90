!f2py: skip
module interpolation_common
    use constants
    implicit none

contains

! 在对数-线性空间中插值并累加SED到观测网格：源为(log x, log y)，目标为log x，输出累加线性y。
subroutine interpolation_accumulate_log_sed(src_x, src_y, num_src, dst_x, num_dst, accum)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: num_src, num_dst
    real(8), intent(in) :: src_x(num_src), src_y(num_src), dst_x(num_dst)
    real(8), intent(inout) :: accum(num_dst)
    integer :: i_dst, i_src
    real(8) :: ratio, y_lo, y_hi

    i_src = 1
    do i_dst = 1, num_dst
        if (dst_x(i_dst) <= src_x(1)) cycle
        if (dst_x(i_dst) > src_x(num_src)) exit
        do while (i_src < num_src - 1 .and. dst_x(i_dst) > src_x(i_src + 1))
            i_src = i_src + 1
        end do
        if (dst_x(i_dst) > src_x(i_src) .and. dst_x(i_dst) <= src_x(i_src + 1)) then
            ratio = (dst_x(i_dst) - src_x(i_src)) / (src_x(i_src + 1) - src_x(i_src))
            if (src_y(i_src) > -huge(one)/two) then
                y_lo = exp(src_y(i_src))
            else
                y_lo = zero
            end if
            if (src_y(i_src + 1) > -huge(one)/two) then
                y_hi = exp(src_y(i_src + 1))
            else
                y_hi = zero
            end if
            if (y_lo > zero .and. y_hi > zero) then
                accum(i_dst) = accum(i_dst) + exp(src_y(i_src) + ratio * (src_y(i_src + 1) - src_y(i_src)))
            else
                accum(i_dst) = accum(i_dst) + (one-ratio)*y_lo + ratio*y_hi
            end if
        end if
    end do
end subroutine interpolation_accumulate_log_sed

end module interpolation_common
