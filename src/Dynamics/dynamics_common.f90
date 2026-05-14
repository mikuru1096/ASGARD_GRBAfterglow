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
                                              gam_m_cross, T, Y, D, para_m_ej, Delta_0, eta_0, &
                                              A_star, dNe_ISM, Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r)
            implicit none
            real(8), intent(inout) :: dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross, gam_m_cross
            real(8), intent(in) :: T, para_m_ej, Delta_0, eta_0, A_star, dNe_ISM, Epsilon_b, Epsilon_e
            real(8), intent(in) :: p_f, f_e, e_r, b_r, p_r, f_e_r
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
                                gam_m_cross, T, H, Y, para_m_ej, Delta_0, eta_0, A_star, &
                                dNe_ISM, Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r)
    implicit none
    procedure(dynamics_reverse_rhs_iface) :: rhs
    integer :: N, J, K
    real(8), intent(inout) :: dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross, gam_m_cross, T, Y(6)
    real(8), intent(in) :: H, para_m_ej, Delta_0, eta_0, A_star, dNe_ISM, Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r
    real(8) :: D(6), A(4), B(6), C(6), G(6), E(6)
    real(8) :: EPS, HH, P, X, DT, TT, dB3_try, T_cross_try, R_cross_try, e3_cross_try, gam20_try
    real(8) :: U3_cross_try, V3_cross_try, M3_cross_try, gam_m_cross_try

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
        do J = 1, N
            call rhs(dB3_try, T_cross_try, R_cross_try, e3_cross_try, gam20_try, U3_cross_try, &
                     V3_cross_try, M3_cross_try, gam_m_cross_try, T, Y, D, &
                     para_m_ej, Delta_0, eta_0, A_star, dNe_ISM, &
                     Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r)
            E = Y
            B = Y
            do K = 1, 3
                Y = E+A(K)*D
                B = B+A(K+1)*D/3.0d0
                TT = T+A(K)
                call rhs(dB3_try, T_cross_try, R_cross_try, e3_cross_try, gam20_try, U3_cross_try, &
                         V3_cross_try, M3_cross_try, gam_m_cross_try, TT, Y, D, &
                         para_m_ej, Delta_0, eta_0, A_star, dNe_ISM, &
                         Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r)
            end do
            Y = B+HH*D/6.0d0
            T = T+DT
        end do
        call dynamics_rk4_error_n(Y, G, 6, P)
        if (P < EPS) then
            dB3=dB3_try; T_cross=T_cross_try; R_cross=R_cross_try
            e3_cross=e3_cross_try; gam20=gam20_try
            U3_cross=U3_cross_try; V3_cross=V3_cross_try; M3_cross=M3_cross_try; gam_m_cross=gam_m_cross_try
        end if
        HH = 0.5d0*HH
        N = N+N
    end do
    T = X
end subroutine dynamics_rk4_reverse

end module dynamics_common
