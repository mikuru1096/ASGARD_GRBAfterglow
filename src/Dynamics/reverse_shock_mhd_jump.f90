module reverse_shock_mhd_jump
    use constants
    implicit none
    private
    public :: rs_vegas_ud, rs_vegas_comp, rs_mag_internal

contains

! Vegas 磁化跳跃条件的解析根：给定相对 Lorentz 因子和上游 sigma，返回下游四速度。
real(8) function rs_vegas_ud(gamma_rel, sigma)
    implicit none
    real(8), intent(in) :: gamma_rel, sigma
    real(8), parameter :: sigma_cut = 1d-6
    real(8) :: gamma_eff, ad, gamma_sq, gm1, gp1, hm1, hm2, term1, term2, u2_down
    real(8) :: a_cub, b_cub, c_cub, d_cub, b_red, c_red, d_red, p_cub, q_cub, amp, arg

    gamma_eff = max(1d0, gamma_rel)
    ad = 4d0/3d0+1d0/(3d0*gamma_eff)
    gamma_sq = gamma_eff*gamma_eff
    gm1 = gamma_eff-1d0; gp1 = gamma_eff+1d0; hm1 = ad-1d0; hm2 = ad-2d0
    term1 = -ad*hm2; term2 = gamma_sq-1d0
    if (sigma <= sigma_cut) then
        u2_down = gm1*hm1*hm1/(term1*gm1+2d0)
    else
        a_cub = term1*gm1+2d0
        b_cub = -gp1*(-hm2*(ad*gamma_sq+1d0)+ad*hm1*gamma_rel)*sigma &
                - gm1*(term1*(gamma_sq-2d0)+2d0*gamma_rel+3d0)
        c_cub = gp1*(ad*(1d0-ad/4d0)*term2+1d0)*sigma*sigma &
                + term2*(2d0*gamma_rel+hm2*(ad*gamma_rel-1d0))*sigma &
                + gp1*gm1*gm1*hm1*hm1
        d_cub = -gm1*gp1*gp1*hm2*hm2*sigma*sigma/4d0
        b_red = b_cub/a_cub; c_red = c_cub/a_cub; d_red = d_cub/a_cub
        p_cub = c_red-b_red*b_red/3d0
        q_cub = 2d0*b_red*b_red*b_red/27d0-b_red*c_red/3d0+d_red
        amp = dsqrt(-p_cub/3d0)
        arg = 3d0*q_cub/(2d0*p_cub*amp)
        arg = max(-1d0, min(1d0, arg))
        u2_down = 2d0*amp*dcos((dacos(arg)-2d0*pi)/3d0)-b_red/3d0
    end if
    rs_vegas_ud = dsqrt(u2_down)
end function rs_vegas_ud

! Vegas 压缩比 u_up/u_down；上游四速度由相对速度合成关系给出。
real(8) function rs_vegas_comp(gamma_rel, sigma)
    implicit none
    real(8), intent(in) :: gamma_rel, sigma
    real(8), parameter :: sigma_cut = 1d-6
    real(8) :: gamma_eff, u_down, u_up, gsq1

    gamma_eff = max(1d0, gamma_rel)
    if (gamma_eff <= 1d0) then
        rs_vegas_comp = 4d0*gamma_eff
        return
    end if
    if (sigma <= sigma_cut) then
        rs_vegas_comp = 4d0*gamma_eff
        return
    end if
    u_down = rs_vegas_ud(gamma_eff, sigma)
    gsq1 = (gamma_eff-1d0)*(gamma_eff+1d0)
    u_up = dsqrt((1d0+u_down*u_down)*gsq1)+u_down*gamma_eff
    rs_vegas_comp = u_up/u_down
end function rs_vegas_comp

! MHD jump 给出的下游热比内能；sigma=0 精确回到 hydrodynamic (gamma_rel-1)。
real(8) function rs_mag_internal(gamma_rel, sigma)
    implicit none
    real(8), intent(in) :: gamma_rel, sigma
    real(8) :: gamma_eff, ad, u_down, u_up, gamma_down, gamma_up, comp_ratio, h_down, gsq1

    gamma_eff = max(1d0, gamma_rel)
    if (sigma <= 0d0) then
        rs_mag_internal = gamma_eff-1d0
    else
        ad = 4d0/3d0+1d0/(3d0*gamma_eff)
        u_down = rs_vegas_ud(gamma_eff, sigma)
        gsq1 = (gamma_eff-1d0)*(gamma_eff+1d0)
        u_up = dsqrt((1d0+u_down*u_down)*gsq1)+u_down*gamma_eff
        gamma_down = dsqrt(1d0+u_down*u_down)
        gamma_up = dsqrt(1d0+u_up*u_up)
        comp_ratio = u_up/u_down
        h_down = (1d0+sigma)*gamma_up/gamma_down-comp_ratio*sigma
        rs_mag_internal = (h_down-1d0)/ad
    end if
end function rs_mag_internal

end module reverse_shock_mhd_jump
