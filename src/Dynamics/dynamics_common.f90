module dynamics_common
    use constants
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none

    abstract interface
        subroutine dynamics_forward_rhs_iface(T, Y, M, D, E_e, E_iso, Eta_0, dNe_ISM, A_star, E_b, p, z, f_e, &
                                              E_inj_t1, E_inj_t2, E_inj, E_inj_q, R_tr, f_jump, f_wide, R0, index_dyn)
            implicit none
            integer, intent(in) :: M, index_dyn
            real(8), intent(in) :: T, E_e, E_iso, Eta_0, dNe_ISM, A_star, E_b, p, z, f_e
            real(8), intent(in) :: E_inj_t1, E_inj_t2, E_inj, E_inj_q, R_tr, f_jump, f_wide, R0
            real(8), intent(in) :: Y(M)
            real(8), intent(out) :: D(M)
        end subroutine dynamics_forward_rhs_iface
    end interface

    abstract interface
        subroutine dynamics_reverse_rhs_iface(dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross, &
                                              gam_m_cross, B3_ordered_cross, T, Y, D, para_m_ej, Delta_0, eta_0, &
                                              A_star, dNe_ISM, Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, &
                                              sigma_r)
            implicit none
            real(8), intent(inout) :: dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross
            real(8), intent(inout) :: gam_m_cross, B3_ordered_cross
            real(8), intent(in) :: T, para_m_ej, Delta_0, eta_0, A_star, dNe_ISM, Epsilon_b, Epsilon_e
            real(8), intent(in) :: p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r
            real(8), intent(in) :: Y(6)
            real(8), intent(out) :: D(6)
        end subroutine dynamics_reverse_rhs_iface
    end interface

contains

subroutine dynamics_deceleration_radius(A_star, dNe_ISM, Eta_0, DM_0, R_dec)
    implicit none
    real(8), intent(in) :: A_star, dNe_ISM, Eta_0, DM_0
    real(8), intent(out) :: R_dec
    real(8) :: R_dec_ism, R_dec_wind

    R_dec_ism = (dNe_ISM*Eta_0/DM_0)**(-one/3d0)
    if (A_star > zero) then
        R_dec_wind = DM_0/(2.0d35*A_star*Eta_0)
        R_dec = min(R_dec_wind, R_dec_ism)
    else
        R_dec = R_dec_ism
    end if
end subroutine dynamics_deceleration_radius

subroutine dynamics_external_density_base(A_star, dNe_ISM, RR, dNe)
    implicit none
    real(8), intent(in) :: A_star, dNe_ISM, RR
    real(8), intent(out) :: dNe
    real(8) :: dNe_wind

    if (A_star > zero) then
        dNe_wind = A_star*3.0d35/RR**2
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe = dNe_ISM
        else
            dNe = dNe_wind
        end if
    else
        dNe = dNe_ISM
    end if
end subroutine dynamics_external_density_base

subroutine dynamics_external_density_profile(A_star, dNe_ISM, RR, R0, apply_jump, R_tr, f_jump, f_wide, dNe)
    implicit none
    integer, intent(in) :: apply_jump
    real(8), intent(in) :: A_star, dNe_ISM, RR, R0, R_tr, f_jump, f_wide
    real(8), intent(out) :: dNe

    call dynamics_external_density_base(A_star, dNe_ISM, RR, dNe)

    if (A_star <= zero .and. apply_jump /= 0) then
        dNe = dNe_ISM*(one+(f_jump-one)*exp(-(log10(RR)-log10(R_tr))**2/(two*f_wide*f_wide)))
    end if

    if (RR < R0) then
        if (A_star > zero) then
            dNe = A_star*3.0d35/R0**2
        else
            dNe = dNe_ISM
        end if
    end if
end subroutine dynamics_external_density_profile

subroutine dynamics_log_time_step(T_base, Grid_Tobs_bin, T_log10, Num_R1, I_tobs, T, H)
    implicit none
    integer, intent(in) :: Num_R1, I_tobs
    real(8), intent(in) :: T_base, Grid_Tobs_bin, T_log10
    real(8), intent(out) :: T, H
    real(8) :: dT

    dT = ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
    T = T_base+dT
    if (I_tobs == 1) then
        H = dT
    else
        H = ten**(Grid_Tobs_bin+T_log10*I_tobs/Num_R1)-dT
    end if
end subroutine dynamics_log_time_step

subroutine dynamics_reverse_gamma_extrema(dB, gamma34, factor2, f_e_r, Gam_e_max, Gam_e_m)
    implicit none
    real(8), intent(in) :: dB, gamma34, factor2, f_e_r
    real(8), intent(out) :: Gam_e_max, Gam_e_m

    Gam_e_max = 3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
    Gam_e_m = factor2*(gamma34-one)/f_e_r+one
end subroutine dynamics_reverse_gamma_extrema

! Vegas 磁化跳跃条件的解析根：给定相对 Lorentz 因子和上游 sigma，返回下游四速度。
real(8) function rs_vegas_ud(gamma_rel, sigma)
    implicit none
    real(8), intent(in) :: gamma_rel, sigma
    real(8), parameter :: sigma_cut = 1d-6
    real(8) :: ad, gamma_sq, gm1, gp1, hm1, hm2, term1, term2, u2_down
    real(8) :: a_cub, b_cub, c_cub, d_cub, b_red, c_red, d_red, p_cub, q_cub, amp, arg

    ad = 4d0/3d0+one/(3d0*gamma_rel)
    gamma_sq = gamma_rel*gamma_rel
    gm1 = gamma_rel-one; gp1 = gamma_rel+one; hm1 = ad-one; hm2 = ad-two
    term1 = -ad*hm2; term2 = gamma_sq-one
    if (sigma <= sigma_cut) then
        u2_down = gm1*hm1*hm1/(term1*gm1+two)
    else
        a_cub = term1*gm1+two
        b_cub = -gp1*(-hm2*(ad*gamma_sq+one)+ad*hm1*gamma_rel)*sigma &
                - gm1*(term1*(gamma_sq-two)+two*gamma_rel+3d0)
        c_cub = gp1*(ad*(one-ad/4d0)*term2+one)*sigma*sigma &
                + term2*(two*gamma_rel+hm2*(ad*gamma_rel-one))*sigma &
                + gp1*gm1*gm1*hm1*hm1
        d_cub = -gm1*gp1*gp1*hm2*hm2*sigma*sigma/4d0
        b_red = b_cub/a_cub; c_red = c_cub/a_cub; d_red = d_cub/a_cub
        p_cub = c_red-b_red*b_red/3d0
        q_cub = two*b_red*b_red*b_red/27d0-b_red*c_red/3d0+d_red
        amp = dsqrt(-p_cub/3d0)
        arg = 3d0*q_cub/(two*p_cub*amp)
        arg = max(-one, min(one, arg))
        u2_down = two*amp*dcos((dacos(arg)-two*pi)/3d0)-b_red/3d0
    end if
    rs_vegas_ud = dsqrt(u2_down)
end function rs_vegas_ud

! Vegas 压缩比 u_up/u_down；上游四速度由相对速度合成关系给出。
real(8) function rs_vegas_comp(gamma_rel, sigma)
    implicit none
    real(8), intent(in) :: gamma_rel, sigma
    real(8) :: u_down, u_up

    u_down = rs_vegas_ud(gamma_rel, sigma)
    if (u_down == zero) then
        rs_vegas_comp = 4d0*gamma_rel
    else
        u_up = dsqrt((one+u_down*u_down)*(gamma_rel*gamma_rel-one))+u_down*gamma_rel
        rs_vegas_comp = u_up/u_down
    end if
end function rs_vegas_comp

! 磁化 RS 压缩比：Vegas 给出 sigma 依赖，ASGARD 非磁化 jump condition 固定零磁化极限。
real(8) function rs_mag_comp(gamma_rel, sigma)
    implicit none
    real(8), intent(in) :: gamma_rel, sigma
    real(8) :: comp_asgard0, comp_vegas0, comp_vegas_sigma

    comp_asgard0 = 4d0*gamma_rel+3d0
    if (sigma <= zero) then
        rs_mag_comp = comp_asgard0
    else
        comp_vegas0 = rs_vegas_comp(gamma_rel, zero)
        comp_vegas_sigma = rs_vegas_comp(gamma_rel, sigma)
        rs_mag_comp = comp_asgard0*comp_vegas_sigma/comp_vegas0
    end if
end function rs_mag_comp

! 上游喷流有序磁场强度，采用 Vegas 的 B4=sqrt(4*pi*c^2*sigma*rho4) 定义。
real(8) function rs_b4_up(rho_up, sigma)
    implicit none
    real(8), intent(in) :: rho_up, sigma

    if (sigma <= zero) then
        rs_b4_up = zero
    else
        rs_b4_up = dsqrt(4d0*pi*para_c*para_c*sigma*rho_up)
    end if
end function rs_b4_up

subroutine dynamics_rk4_error_n(Y, G, M, P)
    implicit none
    integer, intent(in) :: M
    integer :: I
    real(8), intent(in) :: Y(M), G(M)
    real(8), intent(out) :: P
    real(8) :: Q

    P = 0.0d0
    do I = 1, M
        if (.not. ieee_is_finite(Y(I)) .or. .not. ieee_is_finite(G(I))) then
            P = huge(1.0d0)
            return
        end if
        if (Y(I)+G(I) <= zero) then
            P = huge(1.0d0)
            return
        end if
        Q = 2.0d0*abs(Y(I)-G(I))/(Y(I)+G(I))
        if (Q > P) P = Q
    end do
end subroutine dynamics_rk4_error_n

subroutine dynamics_rk4_forward(rhs, T, H, Y, M, EPS, D, B, C, G, E, Epsilon_e, E_iso, Eta_0, dNe_ISM, A_star, Eb, pp, z0, f_e, &
                                E_inj_t1, E_inj_t2, E_inj, E_inj_q, R_tr, f_jump, f_wide, R0, index_dyn)
    implicit none
    procedure(dynamics_forward_rhs_iface) :: rhs
    integer, intent(in) :: M, index_dyn
    integer :: N, J, K
    real(8), intent(inout) :: T, Y(M), D(M), B(M), C(M), G(M), E(M)
    real(8), intent(in) :: H, EPS, Epsilon_e, E_iso, Eta_0, dNe_ISM, A_star, Eb, pp, z0, f_e
    real(8), intent(in) :: E_inj_t1, E_inj_t2, E_inj, E_inj_q, R_tr, f_jump, f_wide, R0
    real(8) :: HH, DT, P, X, TT
    real(8) :: A(4)

    HH = H
    N = 1
    if (index_dyn == 3 .and. T > zero) then
        do while (HH > 0.5d0*T)
            HH = HH/2.0d0
            N = N+N
        end do
    end if
    P = one+EPS
    X = T
    C = Y
    do while (P >= EPS)
        A = (/0.5d0*HH, 0.5d0*HH, HH, HH/)
        G = Y
        Y = C
        T = X
        DT = HH
        do J = 1, N
            call rhs(T, Y, M, D, Epsilon_e, E_iso, Eta_0, dNe_ISM, A_star, Eb, pp, z0, f_e, &
                     E_inj_t1, E_inj_t2, E_inj, E_inj_q, R_tr, f_jump, f_wide, R0, index_dyn)
            E = Y
            B = Y
            do K = 1, 3
                Y = E+A(K)*D
                B = B+A(K+1)*D/3.0d0
                TT = T+A(K)
                call rhs(TT, Y, M, D, Epsilon_e, E_iso, Eta_0, dNe_ISM, A_star, Eb, pp, z0, f_e, &
                         E_inj_t1, E_inj_t2, E_inj, E_inj_q, R_tr, f_jump, f_wide, R0, index_dyn)
            end do
            Y = B+HH*D/6.0d0
            T = T+DT
        end do
        call dynamics_rk4_error_n(Y, G, M, P)
        HH = HH/2.0d0
        N = N+N
    end do
    T = X
end subroutine dynamics_rk4_forward

subroutine dynamics_rk4_reverse(rhs, dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross, &
                                gam_m_cross, B3_ordered_cross, T, H, Y, para_m_ej, Delta_0, eta_0, A_star, &
                                dNe_ISM, Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
    implicit none
    procedure(dynamics_reverse_rhs_iface) :: rhs
    integer :: N, J, K
    real(8), intent(inout) :: dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross
    real(8), intent(inout) :: gam_m_cross, B3_ordered_cross, T, Y(6)
    real(8), intent(in) :: H, para_m_ej, Delta_0, eta_0, A_star, dNe_ISM, Epsilon_b, Epsilon_e
    real(8), intent(in) :: p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r
    real(8) :: D(6), A(4), B(6), C(6), G(6), E(6)
    real(8) :: EPS, HH, P, X, DT, TT, dB3_try, T_cross_try, R_cross_try, e3_cross_try, gam20_try
    real(8) :: U3_cross_try, V3_cross_try, M3_cross_try, gam_m_cross_try, B3_ordered_cross_try

    EPS = 1d-5
    HH = H
    N = 1
    P = one+EPS
    X = T
    C = Y
    do while (P >= EPS)
        A = (/0.5d0*HH, 0.5d0*HH, HH, HH/)
        G = Y
        Y = C
        DT = HH
        T = X
        dB3_try=dB3; T_cross_try=T_cross; R_cross_try=R_cross
        e3_cross_try=e3_cross; gam20_try=gam20
        U3_cross_try=U3_cross; V3_cross_try=V3_cross; M3_cross_try=M3_cross; gam_m_cross_try=gam_m_cross
        B3_ordered_cross_try=B3_ordered_cross
        do J = 1, N
            call rhs(dB3_try, T_cross_try, R_cross_try, e3_cross_try, gam20_try, U3_cross_try, &
                     V3_cross_try, M3_cross_try, gam_m_cross_try, B3_ordered_cross_try, T, Y, D, &
                     para_m_ej, Delta_0, eta_0, A_star, dNe_ISM, &
                     Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
            E = Y
            B = Y
            do K = 1, 3
                Y = E+A(K)*D
                B = B+A(K+1)*D/3.0d0
                TT = T+A(K)
                call rhs(dB3_try, T_cross_try, R_cross_try, e3_cross_try, gam20_try, U3_cross_try, &
                         V3_cross_try, M3_cross_try, gam_m_cross_try, B3_ordered_cross_try, TT, Y, D, &
                         para_m_ej, Delta_0, eta_0, A_star, dNe_ISM, &
                         Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
            end do
            Y = B+HH*D/6.0d0
            T = T+DT
        end do
        call dynamics_rk4_error_n(Y, G, 6, P)
        if (P < EPS) then
            dB3=dB3_try; T_cross=T_cross_try; R_cross=R_cross_try
            e3_cross=e3_cross_try; gam20=gam20_try
            U3_cross=U3_cross_try; V3_cross=V3_cross_try; M3_cross=M3_cross_try; gam_m_cross=gam_m_cross_try
            B3_ordered_cross=B3_ordered_cross_try
        end if
        HH = 0.5d0*HH
        N = N+N
    end do
    T = X
end subroutine dynamics_rk4_reverse

end module dynamics_common
