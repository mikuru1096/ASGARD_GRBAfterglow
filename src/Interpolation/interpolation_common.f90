module interpolation_common
    use constants
    implicit none

contains

subroutine interpolation_phi_setup(Num_Phi, dPhi, phi_scale)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_Phi
    real(8), intent(out) :: dPhi, phi_scale

    dPhi = pi / Num_Phi
    phi_scale = two
    if (Num_Phi == 1) then
        dPhi = pi / 1440d0
        phi_scale = two * 1440d0
    end if
end subroutine interpolation_phi_setup

subroutine interpolation_accumulate_log_sed(src_x, src_y, num_src, dst_x, num_dst, accum)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: num_src, num_dst
    real(8), intent(in) :: src_x(num_src), src_y(num_src), dst_x(num_dst)
    real(8), intent(inout) :: accum(num_dst)
    integer :: i_dst, i_src
    real(8) :: ratio

    i_src = 1
    do i_dst = 1, num_dst
        if (dst_x(i_dst) <= src_x(1)) cycle
        if (dst_x(i_dst) > src_x(num_src)) exit
        do while (i_src < num_src - 1 .and. dst_x(i_dst) > src_x(i_src + 1))
            i_src = i_src + 1
        end do
        if (dst_x(i_dst) > src_x(i_src) .and. dst_x(i_dst) <= src_x(i_src + 1)) then
            ratio = (dst_x(i_dst) - src_x(i_src)) / (src_x(i_src + 1) - src_x(i_src))
            accum(i_dst) = accum(i_dst) + exp(src_y(i_src) + ratio * (src_y(i_src + 1) - src_y(i_src)))
        end if
    end do
end subroutine interpolation_accumulate_log_sed

end module interpolation_common
