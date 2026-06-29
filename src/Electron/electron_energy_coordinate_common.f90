!f2py: skip
module electron_energy_coordinate_common
    use constants
    implicit none
    private

    integer, parameter, public :: electron_coord_log_gamma = 1
    integer, parameter, public :: electron_coord_log_four_velocity_sq = 2
    real(8), parameter, public :: electron_four_velocity_grid_gamma_scale = 2d0

    public :: electron_coord_from_xgamma, electron_xgamma_from_coord
    public :: electron_gamma_from_coord, electron_dxgamma_dcoord
    public :: electron_build_four_velocity_grid

contains

pure real(8) function electron_coord_from_xgamma(coord_kind, four_velocity_sq_scale, x_gamma) result(coord)
    integer, intent(in) :: coord_kind
    real(8), intent(in) :: four_velocity_sq_scale, x_gamma
    real(8) :: gamma, four_velocity_sq

    select case (coord_kind)
    case (electron_coord_log_gamma)
        coord = x_gamma
    case (electron_coord_log_four_velocity_sq)
        gamma = ten**x_gamma
        four_velocity_sq = gamma*gamma - one
        coord = dlog10(one + four_velocity_sq/four_velocity_sq_scale)
    end select
end function electron_coord_from_xgamma

pure real(8) function electron_xgamma_from_coord(coord_kind, four_velocity_sq_scale, coord) result(x_gamma)
    integer, intent(in) :: coord_kind
    real(8), intent(in) :: four_velocity_sq_scale, coord
    real(8) :: gamma

    select case (coord_kind)
    case (electron_coord_log_gamma)
        x_gamma = coord
    case (electron_coord_log_four_velocity_sq)
        gamma = dsqrt(one + four_velocity_sq_scale*(ten**coord - one))
        x_gamma = dlog10(gamma)
    end select
end function electron_xgamma_from_coord

pure real(8) function electron_gamma_from_coord(coord_kind, four_velocity_sq_scale, coord) result(gamma)
    integer, intent(in) :: coord_kind
    real(8), intent(in) :: four_velocity_sq_scale, coord

    select case (coord_kind)
    case (electron_coord_log_gamma)
        gamma = ten**coord
    case (electron_coord_log_four_velocity_sq)
        gamma = dsqrt(one + four_velocity_sq_scale*(ten**coord - one))
    end select
end function electron_gamma_from_coord

pure real(8) function electron_dxgamma_dcoord(coord_kind, four_velocity_sq_scale, coord) result(dxdy)
    integer, intent(in) :: coord_kind
    real(8), intent(in) :: four_velocity_sq_scale, coord
    real(8) :: gamma

    select case (coord_kind)
    case (electron_coord_log_gamma)
        dxdy = one
    case (electron_coord_log_four_velocity_sq)
        gamma = electron_gamma_from_coord(coord_kind, four_velocity_sq_scale, coord)
        dxdy = four_velocity_sq_scale*ten**coord/(two*gamma*gamma)
    end select
end function electron_dxgamma_dcoord

subroutine electron_build_four_velocity_grid(num_gamma, gamma_min, gamma_max, coord_gamma_scale, gamma_grid, coord_edge, x_edge)
    integer, intent(in) :: num_gamma
    integer :: i
    real(8), intent(in) :: gamma_min, gamma_max, coord_gamma_scale
    real(8), intent(out) :: gamma_grid(num_gamma), coord_edge(num_gamma+1), x_edge(num_gamma+1)
    real(8) :: coord_min, coord_max, coord_mid, four_velocity_sq_scale

    if (coord_gamma_scale <= one) error stop 'electron_build_four_velocity_grid requires coord_gamma_scale > 1.'
    if (gamma_max <= gamma_min) error stop 'electron_build_four_velocity_grid requires gamma_max > gamma_min.'
    four_velocity_sq_scale = coord_gamma_scale*coord_gamma_scale - one
    coord_min = electron_coord_from_xgamma(electron_coord_log_four_velocity_sq, four_velocity_sq_scale, dlog10(gamma_min))
    coord_max = electron_coord_from_xgamma(electron_coord_log_four_velocity_sq, four_velocity_sq_scale, dlog10(gamma_max))
    do i = 1, num_gamma + 1
        coord_edge(i) = coord_min + (coord_max - coord_min)*dble(i - 1)/dble(num_gamma)
        x_edge(i) = electron_xgamma_from_coord(electron_coord_log_four_velocity_sq, four_velocity_sq_scale, coord_edge(i))
    enddo
    do i = 1, num_gamma
        coord_mid = 0.5d0*(coord_edge(i) + coord_edge(i + 1))
        gamma_grid(i) = electron_gamma_from_coord(electron_coord_log_four_velocity_sq, four_velocity_sq_scale, coord_mid)
    enddo
end subroutine electron_build_four_velocity_grid

end module electron_energy_coordinate_common
