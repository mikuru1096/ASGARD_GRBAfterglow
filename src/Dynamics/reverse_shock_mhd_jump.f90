module reverse_shock_mhd_jump
    use constants
    implicit none
    private
    public :: rs_vegas_ud, rs_mhd_state

contains

! Vegas 磁化跳跃条件的解析根：给定相对 Lorentz 因子和上游 sigma，返回下游四速度。
real(8) function rs_vegas_ud(gamma_rel, sigma)
    implicit none
    real(8), intent(in) :: gamma_rel, sigma
    real(8), parameter :: sigma_cut = 1d-6
    real(8) :: gamma_eff, ad, gm1, hm1, hm2, term1, u2_down
    real(16) :: gq, sq, adq, gsq, gm1q, gp1q, hm1q, hm2q, term1q, term2q
    real(16) :: aq, bq, cq, dq, discq, rmax, ymax, u2q

    gamma_eff = max(1d0, gamma_rel)
    if (gamma_eff == 1d0) then
        u2_down = sigma
    else if (sigma <= sigma_cut) then
        ad = 4d0/3d0+1d0/(3d0*gamma_eff)
        gm1 = gamma_eff-1d0
        hm1 = ad-1d0
        hm2 = ad-2d0
        term1 = -ad*hm2
        u2_down = gm1*hm1*hm1/(term1*gm1+2d0)
    else
        gq = real(gamma_eff, 16)
        sq = real(sigma, 16)
        adq = 4.0_16/3.0_16+1.0_16/(3.0_16*gq)
        gsq = gq*gq
        gm1q = gq-1.0_16
        gp1q = gq+1.0_16
        hm1q = adq-1.0_16
        hm2q = adq-2.0_16
        term1q = -adq*hm2q
        term2q = gsq-1.0_16
        aq = term1q*gm1q+2.0_16
        bq = -gp1q*(-hm2q*(adq*gsq+1.0_16)+adq*hm1q*gq)*sq &
             -gm1q*(term1q*(gsq-2.0_16)+2.0_16*gq+3.0_16)
        cq = gp1q*(adq*(1.0_16-adq/4.0_16)*term2q+1.0_16)*sq*sq &
             +term2q*(2.0_16*gq+hm2q*(adq*gq-1.0_16))*sq &
             +gp1q*gm1q*gm1q*hm1q*hm1q
        dq = -gm1q*gp1q*gp1q*hm2q*hm2q*sq*sq/4.0_16
        discq = cubic_disc(gq, sq)
        rmax = cubic_max(aq, bq, cq, dq, discq)
        ymax = cubic_max(dq, cq, bq, aq, discq)
        u2q = -dq*ymax/(aq*rmax)
        u2_down = real(u2q, 8)
    end if
    rs_vegas_ud = dsqrt(u2_down)
end function rs_vegas_ud

! 三实根多项式的最大根；正判别式因子通过 atan2 给出无定义域损失的三角角度。
! Largest root of a three-real-root cubic; the positive discriminant gives a domain-safe atan2 angle.
real(16) function cubic_max(a, b, c, d, disc)
    implicit none
    real(16), intent(in) :: a, b, c, d, disc
    real(16) :: delta0, delta1, amp, angle

    delta0 = b*b-3.0_16*a*c
    delta1 = 2.0_16*b*b*b-9.0_16*a*b*c+27.0_16*a*a*d
    amp = sqrt(delta0)/(3.0_16*abs(a))
    angle = atan2(sqrt(27.0_16)*abs(a)*sqrt(disc), -sign(1.0_16, a)*delta1)
    cubic_max = 2.0_16*amp*cos(angle/3.0_16)-b/(3.0_16*a)
end function cubic_max

! 判别式在 t=gamma_rel-1 和 sigma 上写成正系数 Horner 形式，避免近双根的大项消减。
! Positive-coefficient Horner form of the discriminant in t=gamma_rel-1 and sigma.
real(16) function cubic_disc(g, s)
    implicit none
    real(16), intent(in) :: g, s
    real(16) :: t, a0, a1, a2, a3, a4, a5, a6, a7, a8, f, lin

    t = g-1.0_16
    a8 = (((1024.0_16*s+3072.0_16)*s+3136.0_16)*s+1152.0_16)*s+64.0_16
    a7 = (((8448.0_16*s+24768.0_16)*s+24736.0_16)*s+8928.0_16)*s+512.0_16
    a6 = (((30672.0_16*s+85056.0_16)*s+80692.0_16)*s+27872.0_16)*s+1600.0_16
    a5 = (((63736.0_16*s+162420.0_16)*s+141344.0_16)*s+44840.0_16)*s+2432.0_16
    a4 = (((82641.0_16*s+188496.0_16)*s+143996.0_16)*s+39240.0_16)*s+1792.0_16
    a3 = (((68292.0_16*s+135696.0_16)*s+85556.0_16)*s+17728.0_16)*s+512.0_16
    a2 = s*((35068.0_16*s+58716.0_16)*s+27501.0_16)*s+3232.0_16*s
    a1 = s*s*((10224.0_16*s+13752.0_16)*s+3690.0_16)
    a0 = 1296.0_16*s*s*s*(s+1.0_16)
    f = ((((((((a8*t+a7)*t+a6)*t+a5)*t+a4)*t+a3)*t+a2)*t+a1)*t+a0)
    lin = (((8.0_16*s+4.0_16)*t+26.0_16*s+12.0_16)*t+26.0_16*s+8.0_16)*t+9.0_16*s
    cubic_disc = t*(g+1.0_16)**5*lin*lin*f/(104976.0_16*g**8)
end function cubic_disc

! One MHD state evaluation shared by dynamics and prompt shocks.
subroutine rs_mhd_state(gamma_rel, sigma, u_down, compression, specific_internal, shock_allowed)
    implicit none
    real(8), intent(in) :: gamma_rel, sigma
    real(8), intent(out) :: u_down, compression, specific_internal
    logical, intent(out) :: shock_allowed
    real(8), parameter :: sigma_cut = 1d-6
    real(8) :: gamma_eff, ad, u_up, gamma_down, h_down, gsq1, rel_u

    gamma_eff = max(1d0, gamma_rel)
    ad = 4d0/3d0+1d0/(3d0*gamma_eff)
    u_down = rs_vegas_ud(gamma_eff, sigma)
    gsq1 = (gamma_eff-1d0)*(gamma_eff+1d0)
    rel_u = dsqrt(gsq1)
    gamma_down = dsqrt(1d0+u_down*u_down)
    u_up = gamma_down*rel_u+u_down*gamma_eff
    if (gamma_eff <= 1d0 .and. sigma > 0d0) then
        compression = 1d0
    else if (sigma <= sigma_cut) then
        compression = 4d0*gamma_eff
    else
        compression = u_up/u_down
    end if
    if (gamma_eff <= 1d0) then
        specific_internal = 0d0
    else if (sigma <= 0d0) then
        specific_internal = gamma_eff-1d0
    else
        h_down = gamma_eff+rel_u*(u_down*u_down-sigma)/(u_down*gamma_down)
        specific_internal = (h_down-1d0)/ad
    end if
    if (gamma_eff <= 1d0) then
        shock_allowed = .false.
    else
        shock_allowed = u_up*u_up > sigma
    end if
end subroutine rs_mhd_state

end module reverse_shock_mhd_jump
