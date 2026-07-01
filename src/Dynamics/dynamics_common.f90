module dynamics_common
    use constants
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private
    public :: dynamics_forward_rhs_iface, dynamics_reverse_rhs_iface
    public :: density_jump_max, active_density_jump_count, active_density_profile_count, active_density_jump_r
    public :: active_density_jump_factor, active_density_jump_width, reverse_rhs_phase
    public :: dynamics_external_density_profile, dynamics_set_density_jump_profile
    public :: dynamics_rk4_forward
    public :: dynamics_rk4_reverse, dynamics_rk4_reverse_pre_m3, dynamics_rk4_reverse_waiting
    public :: rs_vegas_ud, rs_vegas_comp, rs_mag_specific_internal
    integer :: reverse_rhs_phase = 0
    integer, parameter :: density_jump_max = 8, density_profile_max = 96
    integer, parameter :: density_boundary_r0_index = 27, density_boundary_count_index = 28
    integer, parameter :: density_profile_count_index = density_boundary_count_index+1+3*density_jump_max
    integer :: active_density_jump_count = 0, active_density_profile_count = 0
    real(8) :: active_density_jump_r(density_jump_max) = zero
    real(8) :: active_density_jump_factor(density_jump_max) = one
    real(8) :: active_density_jump_width(density_jump_max) = one
    real(8) :: active_density_profile_r(density_profile_max) = zero
    real(8) :: active_density_profile_n(density_profile_max) = zero

    abstract interface
        subroutine dynamics_forward_rhs_iface(T, Y, M, D, &
                                              E_e, E_iso, Eta_0, dNe_ISM, A_star, E_b, p, z, f_e, &
                                              E_inj_t1, E_inj_t2, E_inj, E_inj_q, &
                                              R_tr, f_jump, f_wide, R0, &
                                              index_dyn)
            implicit none
            integer, intent(in) :: M, index_dyn
            real(8), intent(in) :: T, E_e, E_iso, Eta_0, dNe_ISM, A_star, E_b, p, z, f_e
            real(8), intent(in) :: E_inj_t1, E_inj_t2, E_inj, E_inj_q, R_tr, f_jump, f_wide, R0
            real(8), intent(in) :: Y(M)
            real(8), intent(out) :: D(M)
        end subroutine dynamics_forward_rhs_iface
    end interface

    abstract interface
        subroutine dynamics_reverse_rhs_iface(dB3, T_cross, R_cross, e3_cross, gam20, &
                                              U3_cross, V3_cross, M3_cross, &
                                              gam_m_cross, B3_ordered_cross, &
                                              T, Y, D, M, para_m_ej, V3_scale, Delta_0, &
                                              eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                              Epsilon_b, Epsilon_e, p_f, f_e, &
                                              e_r, b_r, p_r, f_e_r, sigma_r)
            implicit none
            integer, intent(in) :: M
            real(8), intent(inout) :: dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross
            real(8), intent(inout) :: gam_m_cross, B3_ordered_cross
            real(8), intent(in) :: T, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM
            real(8), intent(in) :: R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e
            real(8), intent(in) :: p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r
            real(8), intent(in) :: Y(M)
            real(8), intent(out) :: D(M)
        end subroutine dynamics_reverse_rhs_iface
    end interface

contains

subroutine dynamics_external_density_profile(A_star, dNe_ISM, RR, R0, apply_jump, R_tr, f_jump, f_wide, dNe)
    implicit none
    integer, intent(in) :: apply_jump
    integer :: i
    real(8), intent(in) :: A_star, dNe_ISM, RR, R0, R_tr, f_jump, f_wide
    real(8), intent(out) :: dNe
    real(8) :: dNe_base, dNe_wind, enhancement, width_cm

    select case (active_density_profile_count)
    case (0)
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
    case default
        call dynamics_density_tabulated_profile(RR, dNe)
    end select

    if (apply_jump /= 0 .and. active_density_profile_count == 0) then
        select case (active_density_jump_count)
        case (0)
            dNe = dNe*(one+(f_jump-one)* &
                  exp(-((RR-R_tr)*(RR-R_tr))/(two*(f_wide*R_tr)*(f_wide*R_tr))))
        case default
            dNe_base = dNe
            enhancement = one
            do i = 1, active_density_jump_count
                width_cm = active_density_jump_width(i)*active_density_jump_r(i)
                enhancement = enhancement+(active_density_jump_factor(i)-one)* &
                              exp(-((RR-active_density_jump_r(i))*(RR-active_density_jump_r(i)))/ &
                              (two*width_cm*width_cm))
            end do
            dNe = dNe_base*enhancement
        end select
    end if

    if (RR < R0) then
        select case (active_density_profile_count)
        case (0)
            if (A_star > zero) then
                dNe = A_star*3.0d35/R0**2
            else
                dNe = dNe_ISM
            end if
        case default
            call dynamics_density_tabulated_profile(R0, dNe)
        end select
    end if
end subroutine dynamics_external_density_profile

subroutine dynamics_set_density_jump_profile(Boundary, n)
    implicit none
    integer, intent(in) :: n
    integer :: i, radius_index, factor_index, width_index
    real(8), intent(in) :: Boundary(n)

    active_density_jump_count = 0
    active_density_profile_count = 0
    active_density_jump_r = zero
    active_density_jump_factor = one
    active_density_jump_width = one
    active_density_profile_r = zero
    active_density_profile_n = zero
    if (n < density_boundary_count_index) return
    active_density_jump_count = nint(Boundary(density_boundary_count_index))
    if (active_density_jump_count < 0 .or. active_density_jump_count > density_jump_max) &
        error stop 'density jump count outside supported range'
    radius_index = density_boundary_count_index+1
    factor_index = radius_index+density_jump_max
    width_index = factor_index+density_jump_max
    if (active_density_jump_count > 0 .and. n < width_index+density_jump_max-1) &
        error stop 'boundary is missing density jump arrays'
    do i = 1, active_density_jump_count
        active_density_jump_r(i) = Boundary(radius_index+i-1)
        active_density_jump_factor(i) = Boundary(factor_index+i-1)
        active_density_jump_width(i) = Boundary(width_index+i-1)
        if (active_density_jump_r(i) <= zero .or. active_density_jump_factor(i) <= zero .or. &
            active_density_jump_width(i) <= zero) &
            error stop 'density jump radii, factors, and widths must be positive'
    end do
    if (n < density_profile_count_index) return
    active_density_profile_count = nint(Boundary(density_profile_count_index))
    if (active_density_profile_count < 0 .or. active_density_profile_count > density_profile_max) &
        error stop 'density profile point count outside supported range'
    if (active_density_profile_count == 1) &
        error stop 'density profile requires at least two points'
    radius_index = density_profile_count_index+1
    factor_index = radius_index+density_profile_max
    if (active_density_profile_count > 0 .and. n < factor_index+density_profile_max-1) &
        error stop 'boundary is missing density profile arrays'
    do i = 1, active_density_profile_count
        active_density_profile_r(i) = Boundary(radius_index+i-1)
        active_density_profile_n(i) = Boundary(factor_index+i-1)
        if (active_density_profile_r(i) <= zero .or. active_density_profile_n(i) <= zero) &
            error stop 'density profile radii and densities must be positive'
    end do
    do i = 2, active_density_profile_count
        if (active_density_profile_r(i) <= active_density_profile_r(i-1)) &
            error stop 'density profile radii must be strictly increasing'
    end do
end subroutine dynamics_set_density_jump_profile

subroutine dynamics_density_tabulated_profile(RR, dNe)
    implicit none
    integer :: lo, hi, mid
    real(8), intent(in) :: RR
    real(8), intent(out) :: dNe
    real(8) :: x, x0, x1, y0, y1, weight

    if (RR <= active_density_profile_r(1)) then
        x0 = log(active_density_profile_r(1)); x1 = log(active_density_profile_r(2))
        y0 = log(active_density_profile_n(1)); y1 = log(active_density_profile_n(2))
    else if (RR >= active_density_profile_r(active_density_profile_count)) then
        x0 = log(active_density_profile_r(active_density_profile_count-1))
        x1 = log(active_density_profile_r(active_density_profile_count))
        y0 = log(active_density_profile_n(active_density_profile_count-1))
        y1 = log(active_density_profile_n(active_density_profile_count))
    else
        lo = 1; hi = active_density_profile_count
        do while (hi-lo > 1)
            mid = (lo+hi)/2
            if (active_density_profile_r(mid) <= RR) then
                lo = mid
            else
                hi = mid
            end if
        end do
        x0 = log(active_density_profile_r(lo)); x1 = log(active_density_profile_r(hi))
        y0 = log(active_density_profile_n(lo)); y1 = log(active_density_profile_n(hi))
    end if
    x = log(RR)
    weight = (x-x0)/(x1-x0)
    dNe = exp(y0+weight*(y1-y0))
end subroutine dynamics_density_tabulated_profile

! Vegas 磁化跳跃条件的解析根：给定相对 Lorentz 因子和上游 sigma，返回下游四速度。
real(8) function rs_vegas_ud(gamma_rel, sigma)
    implicit none
    real(8), intent(in) :: gamma_rel, sigma
    real(8), parameter :: sigma_cut = 1d-6
    real(8) :: gamma_eff, ad, gamma_sq, gm1, gp1, hm1, hm2, term1, term2, u2_down
    real(8) :: a_cub, b_cub, c_cub, d_cub, b_red, c_red, d_red, p_cub, q_cub, amp, arg

    gamma_eff = max(one, gamma_rel)
    ad = 4d0/3d0+one/(3d0*gamma_eff)
    gamma_sq = gamma_eff*gamma_eff
    gm1 = gamma_eff-one; gp1 = gamma_eff+one; hm1 = ad-one; hm2 = ad-two
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
    real(8), parameter :: sigma_cut = 1d-6
    real(8) :: gamma_eff, u_down, u_up, gamma_sq_minus_one

    gamma_eff = max(one, gamma_rel)
    if (gamma_eff <= one) then
        rs_vegas_comp = 4d0*gamma_eff
        return
    end if
    if (sigma <= sigma_cut) then
        rs_vegas_comp = 4d0*gamma_eff
        return
    end if
    u_down = rs_vegas_ud(gamma_eff, sigma)
    gamma_sq_minus_one = (gamma_eff-one)*(gamma_eff+one)
    u_up = dsqrt((one+u_down*u_down)*gamma_sq_minus_one)+u_down*gamma_eff
    rs_vegas_comp = u_up/u_down
end function rs_vegas_comp

! MHD jump 给出的下游热比内能；sigma=0 精确回到 hydrodynamic (gamma_rel-1)。
real(8) function rs_mag_specific_internal(gamma_rel, sigma)
    implicit none
    real(8), intent(in) :: gamma_rel, sigma
    real(8) :: gamma_eff, ad, u_down, u_up, gamma_down, gamma_up, comp_ratio, h_down, gamma_sq_minus_one

    gamma_eff = max(one, gamma_rel)
    if (sigma <= zero) then
        rs_mag_specific_internal = gamma_eff-one
    else
        ad = 4d0/3d0+one/(3d0*gamma_eff)
        u_down = rs_vegas_ud(gamma_eff, sigma)
        gamma_sq_minus_one = (gamma_eff-one)*(gamma_eff+one)
        u_up = dsqrt((one+u_down*u_down)*gamma_sq_minus_one)+u_down*gamma_eff
        gamma_down = dsqrt(one+u_down*u_down)
        gamma_up = dsqrt(one+u_up*u_up)
        comp_ratio = u_up/u_down
        h_down = (one+sigma)*gamma_up/gamma_down-comp_ratio*sigma
        rs_mag_specific_internal = (h_down-one)/ad
    end if
end function rs_mag_specific_internal

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

subroutine dynamics_rk4_reverse_error(Y, G, M, P)
    implicit none
    integer, intent(in) :: M
    integer :: I
    real(8), intent(in) :: Y(M), G(M)
    real(8), intent(out) :: P
    real(8) :: Q

    P = 0.0d0
    do I = 1, 3
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
end subroutine dynamics_rk4_reverse_error

subroutine dynamics_rk4_forward(rhs, T, H, Y, M, EPS, D, B, C, G, E, &
                                Epsilon_e, E_iso, Eta_0, dNe_ISM, A_star, Eb, pp, z0, f_e, &
                                E_inj_t1, E_inj_t2, E_inj, E_inj_q, &
                                R_tr, f_jump, f_wide, R0, &
                                index_dyn)
    implicit none
    procedure(dynamics_forward_rhs_iface) :: rhs
    integer, intent(in) :: M, index_dyn
    integer :: N, J, K
    real(8), intent(inout) :: T, Y(M), D(M), B(M), C(M), G(M), E(M)
    real(8), intent(in) :: H, EPS, Epsilon_e, E_iso, Eta_0, dNe_ISM, A_star, Eb, pp, z0, f_e
    real(8), intent(in) :: E_inj_t1, E_inj_t2, E_inj, E_inj_q, R_tr, f_jump, f_wide, R0
    real(8) :: HH, DT, P, X, TT, S, HS, SS, T_phys
    real(8) :: A(4)

    X = T
    C = Y
    if (X > zero) then
        HS = log((X+H)/X)
        HH = HS
        N = 1
        P = one+EPS
        do while (P >= EPS)
            A = (/0.5d0*HH, 0.5d0*HH, HH, HH/)
            G = Y
            Y = C
            S = zero
            do J = 1, N
                T_phys = X*exp(S)
                call rhs(T_phys, Y, M, D, Epsilon_e, E_iso, Eta_0, dNe_ISM, A_star, Eb, pp, z0, f_e, &
                         E_inj_t1, E_inj_t2, E_inj, E_inj_q, R_tr, f_jump, f_wide, R0, index_dyn)
                D = T_phys*D
                E = Y
                B = Y
                do K = 1, 3
                    Y = E+A(K)*D
                    B = B+A(K+1)*D/3.0d0
                    SS = S+A(K)
                    T_phys = X*exp(SS)
                    call rhs(T_phys, Y, M, D, Epsilon_e, E_iso, Eta_0, dNe_ISM, A_star, Eb, pp, z0, f_e, &
                             E_inj_t1, E_inj_t2, E_inj, E_inj_q, R_tr, f_jump, f_wide, R0, index_dyn)
                    D = T_phys*D
                end do
                Y = B+HH*D/6.0d0
                S = S+HH
            end do
            call dynamics_rk4_error_n(Y, G, M, P)
            HH = 0.5d0*HH
            N = N+N
        end do
        T = X
        return
    end if

    HH = H
    N = 1
    P = one+EPS
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

subroutine dynamics_rk4_reverse_plain_step(rhs, dB3, T_cross, R_cross, &
                                           e3_cross, gam20, U3_cross, V3_cross, M3_cross, &
                                           gam_m_cross, B3_ordered_cross, &
                                           T_base, S, H, Y, M, para_m_ej, V3_scale, Delta_0, &
                                           eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                           Epsilon_b, Epsilon_e, p_f, f_e, &
                                           e_r, b_r, p_r, f_e_r, sigma_r, rs_phase)
    implicit none
    procedure(dynamics_reverse_rhs_iface) :: rhs
    integer, intent(in) :: M, rs_phase
    integer :: K
    real(8), intent(inout) :: dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross
    real(8), intent(inout) :: gam_m_cross, B3_ordered_cross, S, Y(M)
    real(8), intent(in) :: T_base, H, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM
    real(8), intent(in) :: R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e
    real(8), intent(in) :: p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r
    real(8) :: D(M), A(4), B(M), E(M), SS, T_phys

    A = (/0.5d0*H, 0.5d0*H, H, H/)
    reverse_rhs_phase = rs_phase
    T_phys = T_base*exp(S)
    call rhs(dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, &
             V3_cross, M3_cross, gam_m_cross, B3_ordered_cross, &
             T_phys, Y, D, M, para_m_ej, V3_scale, Delta_0, &
             eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e, p_f, f_e, &
             e_r, b_r, p_r, f_e_r, sigma_r)
    D = T_phys*D
    E = Y
    B = Y
    do K = 1, 3
        Y = E+A(K)*D
        B = B+A(K+1)*D/3.0d0
        SS = S+A(K)
        T_phys = T_base*exp(SS)
        call rhs(dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, &
                 V3_cross, M3_cross, gam_m_cross, B3_ordered_cross, &
                 T_phys, Y, D, M, para_m_ej, V3_scale, Delta_0, &
                 eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e, p_f, f_e, &
                 e_r, b_r, p_r, f_e_r, sigma_r)
        D = T_phys*D
    end do
    Y = B+H*D/6.0d0
    S = S+H
    reverse_rhs_phase = 0
end subroutine dynamics_rk4_reverse_plain_step

subroutine dynamics_rk4_reverse_event_step(rhs, dB3, T_cross, R_cross, &
                                           e3_cross, gam20, U3_cross, V3_cross, M3_cross, &
                                           gam_m_cross, B3_ordered_cross, &
                                           T_base, S, H, Y, M, para_m_ej, V3_scale, Delta_0, &
                                           eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                           Epsilon_b, Epsilon_e, p_f, f_e, &
                                           e_r, b_r, p_r, f_e_r, sigma_r)
    implicit none
    procedure(dynamics_reverse_rhs_iface) :: rhs
    integer, intent(in) :: M
    integer :: I
    real(8), intent(inout) :: dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross
    real(8), intent(inout) :: gam_m_cross, B3_ordered_cross, S, Y(M)
    real(8), intent(in) :: T_base, H, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM
    real(8), intent(in) :: R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e
    real(8), intent(in) :: p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r
    real(8) :: D(M), Y_trial(M), S_trial, H_lo, H_hi, H_mid, H_event, H_post, T_phys
    real(8), parameter :: event_mass_tol = 1d-7, event_time_tol = 1d-12
    real(8) :: dB3_trial, T_cross_trial, R_cross_trial, e3_cross_trial, gam20_trial
    real(8) :: U3_cross_trial, V3_cross_trial, M3_cross_trial, gam_m_cross_trial, B3_ordered_cross_trial

    if (T_cross >= zero .or. Y(4) >= one) then
        call dynamics_rk4_reverse_plain_step(rhs, dB3, T_cross, R_cross, e3_cross, gam20, &
                                             U3_cross, V3_cross, M3_cross, gam_m_cross, &
                                             B3_ordered_cross, T_base, S, H, Y, M, para_m_ej, V3_scale, Delta_0, &
                                             eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                             Epsilon_b, Epsilon_e, p_f, f_e, &
                                             e_r, b_r, p_r, f_e_r, sigma_r, 2)
        return
    end if

    call advance_phase1_trial(H)
    if (Y_trial(4) < one) then
        dB3=dB3_trial
        S = S_trial
        Y = Y_trial
        return
    end if

    H_lo = zero
    H_hi = H
    do I = 1, 60
        H_mid = 0.5d0*(H_lo+H_hi)
        call advance_phase1_trial(H_mid)
        if (Y_trial(4) >= one) then
            H_hi = H_mid
        else
            H_lo = H_mid
        end if
        if (abs(Y_trial(4)-one) <= event_mass_tol) exit
        if (H_hi-H_lo <= event_time_tol*abs(H)) exit
    end do

    H_event = H_hi
    call dynamics_rk4_reverse_plain_step(rhs, dB3, T_cross, R_cross, e3_cross, gam20, &
                                         U3_cross, V3_cross, M3_cross, gam_m_cross, &
                                         B3_ordered_cross, T_base, S, H_event, Y, M, para_m_ej, V3_scale, Delta_0, &
                                         eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                         Epsilon_b, Epsilon_e, p_f, f_e, &
                                         e_r, b_r, p_r, f_e_r, sigma_r, 1)
    Y(4) = one
    reverse_rhs_phase = 2
    T_phys = T_base*exp(S)
    call rhs(dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, &
             V3_cross, M3_cross, gam_m_cross, B3_ordered_cross, &
             T_phys, Y, D, M, para_m_ej, V3_scale, Delta_0, &
             eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e, p_f, f_e, &
             e_r, b_r, p_r, f_e_r, sigma_r)
    reverse_rhs_phase = 0
    H_post = H-H_event
    if (H_post > zero) then
        call dynamics_rk4_reverse_plain_step(rhs, dB3, T_cross, R_cross, e3_cross, gam20, &
                                             U3_cross, V3_cross, M3_cross, gam_m_cross, &
                                             B3_ordered_cross, T_base, S, H_post, Y, M, para_m_ej, V3_scale, Delta_0, &
                                             eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                             Epsilon_b, Epsilon_e, p_f, f_e, &
                                             e_r, b_r, p_r, f_e_r, sigma_r, 2)
    end if

contains

subroutine reset_trial_state()
    implicit none

    Y_trial = Y
    S_trial = S
    dB3_trial = dB3
    T_cross_trial = T_cross
    R_cross_trial = R_cross
    e3_cross_trial = e3_cross
    gam20_trial = gam20
    U3_cross_trial = U3_cross
    V3_cross_trial = V3_cross
    M3_cross_trial = M3_cross
    gam_m_cross_trial = gam_m_cross
    B3_ordered_cross_trial = B3_ordered_cross
end subroutine reset_trial_state

subroutine advance_phase1_trial(H_try)
    implicit none
    real(8), intent(in) :: H_try

    call reset_trial_state()
    call dynamics_rk4_reverse_plain_step(rhs, dB3_trial, T_cross_trial, R_cross_trial, &
                                         e3_cross_trial, gam20_trial, U3_cross_trial, V3_cross_trial, &
                                         M3_cross_trial, gam_m_cross_trial, B3_ordered_cross_trial, &
                                         T_base, S_trial, H_try, Y_trial, M, para_m_ej, V3_scale, Delta_0, &
                                         eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                         Epsilon_b, Epsilon_e, p_f, f_e, &
                                         e_r, b_r, p_r, f_e_r, sigma_r, 1)
end subroutine advance_phase1_trial
end subroutine dynamics_rk4_reverse_event_step

subroutine dynamics_rk4_reverse_pre_step(rhs, dB3, T, Y, M, H, para_m_ej, V3_scale, Delta_0, &
                                         eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                         Epsilon_b, Epsilon_e, p_f, f_e, &
                                         e_r, b_r, p_r, f_e_r, sigma_r)
    implicit none
    procedure(dynamics_reverse_rhs_iface) :: rhs
    integer, intent(in) :: M
    real(8), intent(inout) :: dB3, T, Y(M)
    real(8), intent(in) :: H, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM
    real(8), intent(in) :: R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e
    real(8), intent(in) :: p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r
    real(8) :: Y0(M), D(M), K1(M), K2(M), K3(M), K4(M), T0, L1, L2, L3, L4
    real(8) :: dummy(8), scale

    Y0 = Y
    T0 = T
    dummy = zero
    reverse_rhs_phase = 1
    call rhs(dB3, dummy(1), dummy(2), dummy(3), dummy(4), dummy(5), dummy(6), dummy(7), dummy(8), &
             dummy(8), T0, Y0, D, M, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM, &
             R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
    scale = one/D(4)
    K1 = D*scale
    K1(4) = zero
    L1 = scale

    Y = Y0+0.5d0*H*K1
    Y(4) = Y0(4)+0.5d0*H
    T = T0+0.5d0*H*L1
    call rhs(dB3, dummy(1), dummy(2), dummy(3), dummy(4), dummy(5), dummy(6), dummy(7), dummy(8), &
             dummy(8), T, Y, D, M, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM, &
             R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
    scale = one/D(4)
    K2 = D*scale
    K2(4) = zero
    L2 = scale

    Y = Y0+0.5d0*H*K2
    Y(4) = Y0(4)+0.5d0*H
    T = T0+0.5d0*H*L2
    call rhs(dB3, dummy(1), dummy(2), dummy(3), dummy(4), dummy(5), dummy(6), dummy(7), dummy(8), &
             dummy(8), T, Y, D, M, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM, &
             R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
    scale = one/D(4)
    K3 = D*scale
    K3(4) = zero
    L3 = scale

    Y = Y0+H*K3
    Y(4) = Y0(4)+H
    T = T0+H*L3
    call rhs(dB3, dummy(1), dummy(2), dummy(3), dummy(4), dummy(5), dummy(6), dummy(7), dummy(8), &
             dummy(8), T, Y, D, M, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM, &
             R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
    scale = one/D(4)
    K4 = D*scale
    K4(4) = zero
    L4 = scale

    Y = Y0+H*(K1+two*K2+two*K3+K4)/6.0d0
    Y(4) = Y0(4)+H
    T = T0+H*(L1+two*L2+two*L3+L4)/6.0d0
    reverse_rhs_phase = 0
end subroutine dynamics_rk4_reverse_pre_step

subroutine dynamics_rk4_reverse_pre_m3(rhs, dB3, T_cross, R_cross, e3_cross, gam20, &
                                       U3_cross, V3_cross, M3_cross, gam_m_cross, B3_ordered_cross, &
                                       T_state, T_target, Y, M, para_m_ej, V3_scale, Delta_0, &
                                       eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                       Epsilon_b, Epsilon_e, p_f, f_e, &
                                       e_r, b_r, p_r, f_e_r, sigma_r)
    implicit none
    procedure(dynamics_reverse_rhs_iface) :: rhs
    integer, intent(in) :: M
    integer :: I, J, N, N_bracket
    real(8), intent(inout) :: dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross
    real(8), intent(inout) :: gam_m_cross, B3_ordered_cross, T_state, Y(M)
    real(8), intent(in) :: T_target, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM
    real(8), intent(in) :: R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e
    real(8), intent(in) :: p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r
    real(8), parameter :: pre_m3_fraction_step_max = 5d-2
    logical :: crossing_first, has_reference
    real(8) :: H_bound, H_event, H_lo, H_hi, H_mid, H_target_est, HH, P, Q
    real(8) :: dB3_try, T_try, T_ref, Y_try(M), G(M), D(M), dummy(8)

    dummy = zero
    reverse_rhs_phase = 1
    call rhs(dB3, dummy(1), dummy(2), dummy(3), dummy(4), dummy(5), dummy(6), dummy(7), dummy(8), &
             dummy(8), T_state, Y, D, M, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM, &
             R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
    reverse_rhs_phase = 0
    H_bound = one-Y(4)
    H_target_est = D(4)*(T_target-T_state)
    H_hi = min(pre_m3_fraction_step_max, H_bound, H_target_est)
    T_try = T_state
    Y_try = Y
    dB3_try = dB3
    N_bracket = max(2, ceiling(H_hi/pre_m3_fraction_step_max))
    HH = H_hi/dble(N_bracket)
    do J = 1, N_bracket
        call dynamics_rk4_reverse_pre_step(rhs, dB3_try, T_try, Y_try, M, HH, para_m_ej, V3_scale, &
                                           Delta_0, eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                           Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
    end do
    do while (T_try < T_target .and. H_hi < H_bound)
        H_hi = min(H_bound, two*H_hi)
        T_try = T_state
        Y_try = Y
        dB3_try = dB3
        N_bracket = max(2, ceiling(H_hi/pre_m3_fraction_step_max))
        HH = H_hi/dble(N_bracket)
        do J = 1, N_bracket
            call dynamics_rk4_reverse_pre_step(rhs, dB3_try, T_try, Y_try, M, HH, para_m_ej, V3_scale, &
                                               Delta_0, eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                               Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
        end do
    end do
    crossing_first = (T_try < T_target)
    if (crossing_first) then
        H_event = H_bound
    else
        H_lo = zero
        do I = 1, 30
            H_mid = 0.5d0*(H_lo+H_hi)
            T_try = T_state
            Y_try = Y
            dB3_try = dB3
            N_bracket = max(2, ceiling(H_mid/pre_m3_fraction_step_max))
            HH = H_mid/dble(N_bracket)
            do J = 1, N_bracket
                call dynamics_rk4_reverse_pre_step(rhs, dB3_try, T_try, Y_try, M, HH, para_m_ej, V3_scale, &
                                                   Delta_0, eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                                   Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
            end do
            if (T_try >= T_target) then
                H_hi = H_mid
            else
                H_lo = H_mid
            end if
        end do
        H_event = H_hi
    end if

    N = max(1, ceiling(H_event/pre_m3_fraction_step_max))
    P = one+1d-5
    has_reference = .false.
    do while (P >= 1d-5)
        T_try = T_state
        Y_try = Y
        dB3_try = dB3
        HH = H_event/dble(N)
        do J = 1, N
            call dynamics_rk4_reverse_pre_step(rhs, dB3_try, T_try, Y_try, M, HH, para_m_ej, V3_scale, &
                                               Delta_0, eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                               Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
        end do
        if (has_reference) then
            call dynamics_rk4_reverse_error(Y_try, G, M, P)
            Q = two*abs(T_try-T_ref)/(T_try+T_ref)
            if (Q > P) P = Q
        else
            P = one+1d-5
            has_reference = .true.
        end if
        if (P < 1d-5) then
            dB3 = dB3_try
            T_state = T_try
            Y = Y_try
        else
            G = Y_try
            T_ref = T_try
            N = N+N
        end if
    end do

    if (crossing_first) then
        Y(4) = one
        dummy = zero
        reverse_rhs_phase = 2
        call rhs(dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross, gam_m_cross, &
                 B3_ordered_cross, T_state, Y, D, M, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM, &
                 R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e, p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r)
        reverse_rhs_phase = 0
    end if
end subroutine dynamics_rk4_reverse_pre_m3

subroutine dynamics_rk4_reverse_waiting(rhs, dB3, T_cross, R_cross, &
                                        e3_cross, gam20, U3_cross, V3_cross, M3_cross, &
                                        gam_m_cross, B3_ordered_cross, &
                                        T, H, Y, M, para_m_ej, V3_scale, Delta_0, &
                                        eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                        Epsilon_b, Epsilon_e, p_f, f_e, &
                                        e_r, b_r, p_r, f_e_r, sigma_r)
    implicit none
    procedure(dynamics_reverse_rhs_iface) :: rhs
    integer, intent(in) :: M
    integer :: N, J
    real(8), intent(inout) :: dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross
    real(8), intent(inout) :: gam_m_cross, B3_ordered_cross, T, Y(M)
    real(8), intent(in) :: H, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM
    real(8), intent(in) :: R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e
    real(8), intent(in) :: p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r
    logical :: has_reference
    real(8) :: C(M), G(M), EPS, HH, P, X, S, HS
    real(8) :: dB3_try, T_cross_try, R_cross_try, e3_cross_try, gam20_try
    real(8) :: U3_cross_try, V3_cross_try, M3_cross_try, gam_m_cross_try, B3_ordered_cross_try

    EPS = 1d-5
    N = 1
    HS = log((T+H)/T)
    HH = HS
    P = one+EPS
    X = T
    C = Y
    has_reference = .false.
    do while (P >= EPS)
        Y = C
        S = zero
        dB3_try=dB3; T_cross_try=T_cross; R_cross_try=R_cross
        e3_cross_try=e3_cross; gam20_try=gam20
        U3_cross_try=U3_cross; V3_cross_try=V3_cross
        M3_cross_try=M3_cross; gam_m_cross_try=gam_m_cross
        B3_ordered_cross_try=B3_ordered_cross
        do J = 1, N
            call dynamics_rk4_reverse_plain_step(rhs, dB3_try, T_cross_try, R_cross_try, &
                                                 e3_cross_try, gam20_try, U3_cross_try, V3_cross_try, &
                                                 M3_cross_try, gam_m_cross_try, B3_ordered_cross_try, &
                                                 X, S, HH, Y, M, para_m_ej, V3_scale, Delta_0, &
                                                 eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                                 Epsilon_b, Epsilon_e, p_f, f_e, &
                                                 e_r, b_r, p_r, f_e_r, sigma_r, -1)
        end do
        if (has_reference) then
            call dynamics_rk4_reverse_error(Y, G, M, P)
        else
            P = one+EPS
            has_reference = .true.
        end if
        if (P < EPS) then
            dB3=dB3_try
        else
            G = Y
            HH = 0.5d0*HH
            N = N+N
        end if
    end do
    T = X
end subroutine dynamics_rk4_reverse_waiting

subroutine dynamics_rk4_reverse(rhs, dB3, T_cross, R_cross, &
                                e3_cross, gam20, U3_cross, V3_cross, M3_cross, &
                                gam_m_cross, B3_ordered_cross, &
                                T, H, Y, M, para_m_ej, V3_scale, Delta_0, &
                                eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                Epsilon_b, Epsilon_e, p_f, f_e, &
                                e_r, b_r, p_r, f_e_r, sigma_r)
    implicit none
    procedure(dynamics_reverse_rhs_iface) :: rhs
    integer, intent(in) :: M
    integer :: N, J
    real(8), intent(inout) :: dB3, T_cross, R_cross, e3_cross, gam20, U3_cross, V3_cross, M3_cross
    real(8), intent(inout) :: gam_m_cross, B3_ordered_cross, T, Y(M)
    real(8), intent(in) :: H, para_m_ej, V3_scale, Delta_0, eta_0, A_star, dNe_ISM
    real(8), intent(in) :: R_tr, f_jump, f_wide, R0, Epsilon_b, Epsilon_e
    real(8), intent(in) :: p_f, f_e, e_r, b_r, p_r, f_e_r, sigma_r
    logical :: has_reference
    real(8) :: C(M), G(M), EPS, HH, P, X, S, HS
    real(8) :: dB3_try, T_cross_try, R_cross_try, e3_cross_try, gam20_try
    real(8) :: U3_cross_try, V3_cross_try, M3_cross_try, gam_m_cross_try, B3_ordered_cross_try

    EPS = 1d-5
    N = 1
    HS = log((T+H)/T)
    HH = HS
    P = one+EPS
    X = T
    C = Y
    has_reference = .false.
    do while (P >= EPS)
        Y = C
        S = zero
        dB3_try=dB3; T_cross_try=T_cross; R_cross_try=R_cross
        e3_cross_try=e3_cross; gam20_try=gam20
        U3_cross_try=U3_cross; V3_cross_try=V3_cross
        M3_cross_try=M3_cross; gam_m_cross_try=gam_m_cross
        B3_ordered_cross_try=B3_ordered_cross
        do J = 1, N
            call dynamics_rk4_reverse_event_step(rhs, dB3_try, T_cross_try, R_cross_try, &
                                                 e3_cross_try, gam20_try, U3_cross_try, V3_cross_try, &
                                                 M3_cross_try, gam_m_cross_try, B3_ordered_cross_try, &
                                                 X, S, HH, Y, M, para_m_ej, V3_scale, Delta_0, &
                                                 eta_0, A_star, dNe_ISM, R_tr, f_jump, f_wide, R0, &
                                                 Epsilon_b, Epsilon_e, p_f, f_e, &
                                                 e_r, b_r, p_r, f_e_r, sigma_r)
        end do
        if (has_reference) then
            call dynamics_rk4_reverse_error(Y, G, M, P)
        else
            P = one+EPS
            has_reference = .true.
        end if
        if (P < EPS) then
            dB3=dB3_try; T_cross=T_cross_try; R_cross=R_cross_try
            e3_cross=e3_cross_try; gam20=gam20_try
            U3_cross=U3_cross_try; V3_cross=V3_cross_try
            M3_cross=M3_cross_try; gam_m_cross=gam_m_cross_try
            B3_ordered_cross=B3_ordered_cross_try
        else
            G = Y
            HH = 0.5d0*HH
            N = N+N
        end if
    end do
    T = X
end subroutine dynamics_rk4_reverse

end module dynamics_common
