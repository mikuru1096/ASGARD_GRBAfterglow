module dynamics_density_profile
    use constants
    implicit none
    private
    public :: jump_max, jump_count, profile_count
    public :: jump_radius, jump_factor, jump_width
    public :: density_profile, set_density_profile

    integer, parameter :: jump_max = 8, profile_max = 96
    integer, parameter :: jump_slot = 28
    integer, parameter :: profile_slot = jump_slot+1+3*jump_max
    integer :: jump_count = 0, profile_count = 0
    real(8), dimension(jump_max) :: jump_radius= 0d0
    real(8), dimension(jump_max) :: jump_factor= 1d0
    real(8), dimension(jump_max) :: jump_width= 1d0
    real(8), dimension(profile_max) :: profile_radius= 0d0
    real(8), dimension(profile_max) :: profile_density= 0d0

contains

subroutine density_profile(A_star, dNe_ISM, RR, R0, apply_jump, R_tr, f_jump, f_wide, dNe)
    implicit none
    integer, intent(in) :: apply_jump
    integer :: i
    real(8), intent(in) :: A_star, dNe_ISM, RR, R0, R_tr, f_jump, f_wide
    real(8), intent(out) :: dNe
    real(8) :: dNe_base, dNe_wind, enhancement, width_cm

    select case (profile_count)
    case (0)
        if (A_star > 0d0) then
            dNe_wind = A_star*3.0d35/RR**2
            if (dNe_wind <= dNe_ISM/4d0) then
                dNe = dNe_ISM
            else
                dNe = dNe_wind
            end if
        else
            dNe = dNe_ISM
        end if
    case default
        call tab_density(RR, dNe)
    end select

    if (apply_jump /= 0 .and. profile_count == 0) then
        select case (jump_count)
        case (0)
            dNe = dNe*(1d0+(f_jump-1d0)* &
                  exp(-((RR-R_tr)*(RR-R_tr))/(2d0*(f_wide*R_tr)*(f_wide*R_tr))))
        case default
            dNe_base = dNe
            enhancement = 1d0
            do i = 1, jump_count
                width_cm = jump_width(i)*jump_radius(i)
                enhancement = enhancement+(jump_factor(i)-1d0)* &
                              exp(-((RR-jump_radius(i))*(RR-jump_radius(i)))/ &
                              (2d0*width_cm*width_cm))
            end do
            dNe = dNe_base*enhancement
        end select
    end if

    if (RR < R0) then
        select case (profile_count)
        case (0)
            if (A_star > 0d0) then
                dNe = A_star*3.0d35/R0**2
            else
                dNe = dNe_ISM
            end if
        case default
            call tab_density(R0, dNe)
        end select
    end if
end subroutine density_profile

subroutine set_density_profile(Boundary, n)
    implicit none
    integer, intent(in) :: n
    integer :: i, radius_index, factor_index, width_index
    real(8), intent(in), dimension(n) :: Boundary

    jump_count = 0
    profile_count = 0
    jump_radius = 0d0
    jump_factor = 1d0
    jump_width = 1d0
    profile_radius = 0d0
    profile_density = 0d0
    if (n < jump_slot) return
    jump_count = nint(Boundary(jump_slot))
    if (jump_count < 0 .or. jump_count > jump_max) &
        error stop 'density jump count outside supported range'
    radius_index = jump_slot+1
    factor_index = radius_index+jump_max
    width_index = factor_index+jump_max
    if (jump_count > 0 .and. n < width_index+jump_max-1) &
        error stop 'boundary is missing density jump arrays'
    do i = 1, jump_count
        jump_radius(i) = Boundary(radius_index+i-1)
        jump_factor(i) = Boundary(factor_index+i-1)
        jump_width(i) = Boundary(width_index+i-1)
        if (jump_radius(i) <= 0d0 .or. jump_factor(i) <= 0d0 .or. &
            jump_width(i) <= 0d0) &
            error stop 'density jump radii, factors, and widths must be positive'
    end do
    if (n < profile_slot) return
    profile_count = nint(Boundary(profile_slot))
    if (profile_count < 0 .or. profile_count > profile_max) &
        error stop 'density profile point count outside supported range'
    if (profile_count == 1) &
        error stop 'density profile requires at least 2d0 points'
    radius_index = profile_slot+1
    factor_index = radius_index+profile_max
    if (profile_count > 0 .and. n < factor_index+profile_max-1) &
        error stop 'boundary is missing density profile arrays'
    do i = 1, profile_count
        profile_radius(i) = Boundary(radius_index+i-1)
        profile_density(i) = Boundary(factor_index+i-1)
        if (profile_radius(i) <= 0d0 .or. profile_density(i) <= 0d0) &
            error stop 'density profile radii and densities must be positive'
    end do
    do i = 2, profile_count
        if (profile_radius(i) <= profile_radius(i-1)) &
            error stop 'density profile radii must be strictly increasing'
    end do
end subroutine set_density_profile

subroutine tab_density(RR, dNe)
    implicit none
    integer :: lo, hi, mid
    real(8), intent(in) :: RR
    real(8), intent(out) :: dNe
    real(8) :: x, x0, x1, y0, y1, weight

    if (RR <= profile_radius(1)) then
        x0 = log(profile_radius(1)); x1 = log(profile_radius(2))
        y0 = log(profile_density(1)); y1 = log(profile_density(2))
    else if (RR >= profile_radius(profile_count)) then
        x0 = log(profile_radius(profile_count-1))
        x1 = log(profile_radius(profile_count))
        y0 = log(profile_density(profile_count-1))
        y1 = log(profile_density(profile_count))
    else
        lo = 1; hi = profile_count
        do while (hi-lo > 1)
            mid = (lo+hi)/2
            if (profile_radius(mid) <= RR) then
                lo = mid
            else
                hi = mid
            end if
        end do
        x0 = log(profile_radius(lo)); x1 = log(profile_radius(hi))
        y0 = log(profile_density(lo)); y1 = log(profile_density(hi))
    end if
    x = log(RR)
    weight = (x-x0)/(x1-x0)
    dNe = exp(y0+weight*(y1-y0))
end subroutine tab_density

end module dynamics_density_profile
