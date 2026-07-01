!f2py: skip
module hadronic_species_transport_kernel
    use constants
    implicit none
    private

    real(8), parameter :: neutron_mass_gev = Para_m_n_GeV
    real(8), parameter :: pion_mass_gev = Para_m_pi_charged_GeV
    real(8), parameter :: muon_mass_gev = Para_m_mu_GeV
    real(8), parameter :: neutron_beta_decay_s = 879.4d0
    real(8), parameter :: charged_pion_decay_s = 2.6033d-8
    real(8), parameter :: muon_decay_s = 2.1969811d-6
    real(8), parameter :: gev_c2_to_g = Para_GeV2erg / (Para_c * Para_c)
    real(8), parameter :: log_spacing_rtol = 1d-3

    public :: hadronic_species_spherical_divergence_rate
    public :: hadronic_species_synchrotron_dgamma_dt
    public :: hadronic_species_adiabatic_dgamma_dt
    public :: hadronic_species_advance_operator

contains

! 球对称膨胀发散率：3*v_exp/r，用于绝热冷却计算。
real(8) function hadronic_species_spherical_divergence_rate(radius_cm,expansion_speed_cm_s)
    real(8), intent(in) :: radius_cm, expansion_speed_cm_s

    if (radius_cm <= zero) then
        error stop "hadronic species transport: radius_cm must be positive."
    end if
    if (expansion_speed_cm_s < zero) then
        error stop "hadronic species transport: expansion_speed_cm_s must be non-negative."
    end if
    hadronic_species_spherical_divergence_rate = 3d0*expansion_speed_cm_s/radius_cm
end function hadronic_species_spherical_divergence_rate

! 计算粒子同步辐射冷却率 dγ/dt = -(4/3)σ_T,c U_B γ²/(m c)，中性粒子返回零。
subroutine hadronic_species_synchrotron_dgamma_dt(num_gamma,gamma,b_field_g,mass_gev,charge_number,dgamma_dt)
    integer, intent(in) :: num_gamma
    integer, intent(in) :: charge_number
    real(8), intent(in) :: gamma(num_gamma), b_field_g, mass_gev
    real(8), intent(out) :: dgamma_dt(num_gamma)
    integer :: i_gamma
    real(8) :: mass_g, u_b, sigma_t_species

    call hadronic_species_validate_gamma_grid(num_gamma,gamma)
    if (b_field_g < zero) then
        error stop "hadronic species transport: b_field_g must be non-negative."
    end if
    if (mass_gev <= zero) then
        error stop "hadronic species transport: mass_gev must be positive."
    end if
    if (charge_number == 0) then
        dgamma_dt = zero
        return
    end if

    mass_g = mass_gev*gev_c2_to_g
    u_b = b_field_g*b_field_g/(8d0*pi)
    sigma_t_species = Para_sigmaT*(charge_number**4)*(Para_m_e/mass_g)**2
    do i_gamma=1,num_gamma
        dgamma_dt(i_gamma) = -(4d0/3d0)*sigma_t_species*u_b*gamma(i_gamma)**2/(mass_g*Para_c)
    end do
end subroutine hadronic_species_synchrotron_dgamma_dt

! 计算绝热冷却率 dγ/dt = -(∇·v/3)γ。
subroutine hadronic_species_adiabatic_dgamma_dt(num_gamma,gamma,divergence_rate_s_inv,dgamma_dt)
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: gamma(num_gamma), divergence_rate_s_inv
    real(8), intent(out) :: dgamma_dt(num_gamma)
    integer :: i_gamma

    call hadronic_species_validate_gamma_grid(num_gamma,gamma)
    if (divergence_rate_s_inv < zero) then
        error stop "hadronic species transport: divergence_rate_s_inv must be non-negative."
    end if

    do i_gamma=1,num_gamma
        dgamma_dt(i_gamma) = -(divergence_rate_s_inv/3d0)*gamma(i_gamma)
    end do
end subroutine hadronic_species_adiabatic_dgamma_dt

! 粒子种类输运主算子：对中子/π±/μ子各分量推进冷却+衰变，使用迎风+Strang分裂格式。
subroutine hadronic_species_advance_operator(num_gamma,gamma,dt_s,b_field_g,divergence_rate_s_inv, &
                                             neutron_prev,pion_plus_prev,pion_minus_prev, &
                                             muon_minus_left_prev,muon_minus_right_prev, &
                                             muon_plus_left_prev,muon_plus_right_prev, &
                                             neutron_source_per_gamma_s,pion_plus_source_per_gamma_s, &
                                             pion_minus_source_per_gamma_s,muon_minus_left_source_per_gamma_s, &
                                             muon_minus_right_source_per_gamma_s, &
                                             muon_plus_left_source_per_gamma_s, &
                                             muon_plus_right_source_per_gamma_s, &
                                             neutron_next,pion_plus_next,pion_minus_next, &
                                             muon_minus_left_next,muon_minus_right_next, &
                                             muon_plus_left_next,muon_plus_right_next)
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: gamma(num_gamma), dt_s, b_field_g, divergence_rate_s_inv
    ! 上一步粒子谱（7分量）
    real(8), intent(in) :: neutron_prev(num_gamma), pion_plus_prev(num_gamma), pion_minus_prev(num_gamma), &
        muon_minus_left_prev(num_gamma), muon_minus_right_prev(num_gamma), &
        muon_plus_left_prev(num_gamma), muon_plus_right_prev(num_gamma)
    ! 源项（7分量）
    real(8), intent(in) :: neutron_source_per_gamma_s(num_gamma), pion_plus_source_per_gamma_s(num_gamma), &
        pion_minus_source_per_gamma_s(num_gamma), muon_minus_left_source_per_gamma_s(num_gamma), &
        muon_minus_right_source_per_gamma_s(num_gamma), muon_plus_left_source_per_gamma_s(num_gamma), &
        muon_plus_right_source_per_gamma_s(num_gamma)
    ! 下一步粒子谱（7分量）
    real(8), intent(out) :: neutron_next(num_gamma), pion_plus_next(num_gamma), pion_minus_next(num_gamma), &
        muon_minus_left_next(num_gamma), muon_minus_right_next(num_gamma), &
        muon_plus_left_next(num_gamma), muon_plus_right_next(num_gamma)
    real(8) :: dgamma_syn_pion(num_gamma), dgamma_syn_muon(num_gamma), dgamma_ad(num_gamma), dgamma_total(num_gamma)

    if (dt_s <= zero) then
        error stop "hadronic species transport: dt_s must be positive."
    end if

    call validate_species_transport_inputs
    call hadronic_species_synchrotron_dgamma_dt(num_gamma,gamma,b_field_g,neutron_mass_gev,0,dgamma_total)
    call hadronic_species_synchrotron_dgamma_dt(num_gamma,gamma,b_field_g,pion_mass_gev,1,dgamma_syn_pion)
    call hadronic_species_synchrotron_dgamma_dt(num_gamma,gamma,b_field_g,muon_mass_gev,1,dgamma_syn_muon)
    call hadronic_species_adiabatic_dgamma_dt(num_gamma,gamma,divergence_rate_s_inv,dgamma_ad)

    dgamma_total = dgamma_ad
    call hadronic_species_advance_one(num_gamma,gamma,neutron_prev,neutron_source_per_gamma_s,dt_s, &
                                      neutron_beta_decay_s,dgamma_total,neutron_next)

    dgamma_total = dgamma_syn_pion + dgamma_ad
    call hadronic_species_advance_one(num_gamma,gamma,pion_plus_prev,pion_plus_source_per_gamma_s,dt_s, &
                                      charged_pion_decay_s,dgamma_total,pion_plus_next)
    call hadronic_species_advance_one(num_gamma,gamma,pion_minus_prev,pion_minus_source_per_gamma_s,dt_s, &
                                      charged_pion_decay_s,dgamma_total,pion_minus_next)

    dgamma_total = dgamma_syn_muon + dgamma_ad
    call advance_muon_helicity_species

contains

    subroutine validate_species_transport_inputs
        call hadronic_species_validate_gamma_grid(num_gamma,gamma)
        call hadronic_species_validate_non_negative(num_gamma,neutron_prev,"neutron_prev")
        call hadronic_species_validate_non_negative(num_gamma,pion_plus_prev,"pion_plus_prev")
        call hadronic_species_validate_non_negative(num_gamma,pion_minus_prev,"pion_minus_prev")
        call hadronic_species_validate_non_negative(num_gamma,muon_minus_left_prev,"muon_minus_left_prev")
        call hadronic_species_validate_non_negative(num_gamma,muon_minus_right_prev,"muon_minus_right_prev")
        call hadronic_species_validate_non_negative(num_gamma,muon_plus_left_prev,"muon_plus_left_prev")
        call hadronic_species_validate_non_negative(num_gamma,muon_plus_right_prev,"muon_plus_right_prev")
        call hadronic_species_validate_non_negative(num_gamma,neutron_source_per_gamma_s, &
                                                    "neutron_source_per_gamma_s")
        call hadronic_species_validate_non_negative(num_gamma,pion_plus_source_per_gamma_s, &
                                                    "pion_plus_source_per_gamma_s")
        call hadronic_species_validate_non_negative(num_gamma,pion_minus_source_per_gamma_s, &
                                                    "pion_minus_source_per_gamma_s")
        call hadronic_species_validate_non_negative(num_gamma,muon_minus_left_source_per_gamma_s, &
                                                    "muon_minus_left_source_per_gamma_s")
        call hadronic_species_validate_non_negative(num_gamma,muon_minus_right_source_per_gamma_s, &
                                                    "muon_minus_right_source_per_gamma_s")
        call hadronic_species_validate_non_negative(num_gamma,muon_plus_left_source_per_gamma_s, &
                                                    "muon_plus_left_source_per_gamma_s")
        call hadronic_species_validate_non_negative(num_gamma,muon_plus_right_source_per_gamma_s, &
                                                    "muon_plus_right_source_per_gamma_s")
    end subroutine validate_species_transport_inputs

    subroutine advance_muon_helicity_species
        call hadronic_species_advance_one(num_gamma,gamma,muon_minus_left_prev, &
                                          muon_minus_left_source_per_gamma_s,dt_s, &
                                          muon_decay_s,dgamma_total,muon_minus_left_next)
        call hadronic_species_advance_one(num_gamma,gamma,muon_minus_right_prev, &
                                          muon_minus_right_source_per_gamma_s,dt_s, &
                                          muon_decay_s,dgamma_total,muon_minus_right_next)
        call hadronic_species_advance_one(num_gamma,gamma,muon_plus_left_prev, &
                                          muon_plus_left_source_per_gamma_s,dt_s, &
                                          muon_decay_s,dgamma_total,muon_plus_left_next)
        call hadronic_species_advance_one(num_gamma,gamma,muon_plus_right_prev, &
                                          muon_plus_right_source_per_gamma_s,dt_s, &
                                          muon_decay_s,dgamma_total,muon_plus_right_next)
    end subroutine advance_muon_helicity_species
end subroutine hadronic_species_advance_operator

! 单一种类输运推进：迎风格式+Strang分裂（半衰变步-输运-半衰变步），Courant自适应子步。
subroutine hadronic_species_advance_one(num_gamma,gamma,density_prev,source_per_gamma_s,dt_s, &
                                        decay_lifetime_s,dgamma_dt,density_next)
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: gamma(num_gamma), density_prev(num_gamma), source_per_gamma_s(num_gamma), dt_s
    real(8), intent(in) :: decay_lifetime_s, dgamma_dt(num_gamma)
    real(8), intent(out) :: density_next(num_gamma)
    integer :: i_gamma, i_face, i_sub, n_sub
    real(8) :: dln_gamma, max_courant, dt_sub
    real(8) :: tau_decay(num_gamma), n_xi_next(num_gamma), q_xi(num_gamma), a_center(num_gamma)
    real(8) :: a_face(num_gamma+1), n_ext(num_gamma+2), flux(num_gamma+1), div_flux(num_gamma)
    real(8) :: decay_factor_half(num_gamma)

    call hadronic_species_validate_non_negative(num_gamma,density_prev,"density_prev")
    call hadronic_species_validate_non_negative(num_gamma,source_per_gamma_s,"source_per_gamma_s")
    dln_gamma = hadronic_species_log_spacing(num_gamma,gamma)

    do i_gamma=1,num_gamma
        tau_decay(i_gamma) = decay_lifetime_s*gamma(i_gamma)
        n_xi_next(i_gamma) = gamma(i_gamma)*density_prev(i_gamma)
        q_xi(i_gamma) = gamma(i_gamma)*source_per_gamma_s(i_gamma)
        a_center(i_gamma) = dgamma_dt(i_gamma)/gamma(i_gamma)
    end do

    a_face(1) = a_center(1)
    a_face(num_gamma+1) = a_center(num_gamma)
    do i_face=2,num_gamma
        a_face(i_face) = 0.5d0*(a_center(i_face-1)+a_center(i_face))
    end do

    max_courant = zero
    do i_face=1,num_gamma+1
        max_courant = max(max_courant,abs(a_face(i_face))*dt_s/dln_gamma)
    end do
    n_sub = max(1,ceiling(max_courant))
    dt_sub = dt_s/dble(n_sub)

    do i_gamma=1,num_gamma
        decay_factor_half(i_gamma) = dexp(-0.5d0*dt_sub/tau_decay(i_gamma))
    end do

    do i_sub=1,n_sub
        n_xi_next = n_xi_next*decay_factor_half
        n_ext(1) = zero
        n_ext(num_gamma+2) = zero
        n_ext(2:num_gamma+1) = n_xi_next

        do i_face=1,num_gamma+1
            if (a_face(i_face) >= zero) then
                flux(i_face) = a_face(i_face)*n_ext(i_face)
            else
                flux(i_face) = a_face(i_face)*n_ext(i_face+1)
            end if
        end do

        do i_gamma=1,num_gamma
            div_flux(i_gamma) = (flux(i_gamma+1)-flux(i_gamma))/dln_gamma
            n_xi_next(i_gamma) = n_xi_next(i_gamma) + dt_sub*(-div_flux(i_gamma)+q_xi(i_gamma))
            if (n_xi_next(i_gamma) < zero) n_xi_next(i_gamma) = zero
        end do
        n_xi_next = n_xi_next*decay_factor_half
    end do

    do i_gamma=1,num_gamma
        density_next(i_gamma) = n_xi_next(i_gamma)/gamma(i_gamma)
    end do
end subroutine hadronic_species_advance_one

! 验证gamma网格至少2点、全正、严格递增。
subroutine hadronic_species_validate_gamma_grid(num_gamma,gamma)
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: gamma(num_gamma)
    integer :: i_gamma

    if (num_gamma < 2) then
        error stop "hadronic species transport: gamma grid must have at least two points."
    end if
    do i_gamma=1,num_gamma
        if (gamma(i_gamma) <= zero) then
            error stop "hadronic species transport: gamma grid must be positive."
        end if
    end do
    do i_gamma=2,num_gamma
        if (gamma(i_gamma) <= gamma(i_gamma-1)) then
            error stop "hadronic species transport: gamma grid must be strictly increasing."
        end if
    end do
end subroutine hadronic_species_validate_gamma_grid

! 验证数组所有值为非负。
subroutine hadronic_species_validate_non_negative(num_gamma,values,name)
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: values(num_gamma)
    character(len=*), intent(in) :: name
    integer :: i_gamma

    do i_gamma=1,num_gamma
        if (values(i_gamma) < zero) then
            error stop "hadronic species transport: "//trim(name)//" must be non-negative."
        end if
    end do
end subroutine hadronic_species_validate_non_negative

! 检验并返回gamma网格的平均对数间距（允许松弛容差）。
real(8) function hadronic_species_log_spacing(num_gamma,gamma)
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: gamma(num_gamma)
    integer :: i_gamma
    real(8) :: dlog_i, dlog_ref

    dlog_ref = zero
    do i_gamma=1,num_gamma-1
        dlog_ref = dlog_ref + dlog(gamma(i_gamma+1)/gamma(i_gamma))
    end do
    dlog_ref = dlog_ref/dble(num_gamma-1)

    do i_gamma=1,num_gamma-1
        dlog_i = dlog(gamma(i_gamma+1)/gamma(i_gamma))
        if (abs(dlog_i-dlog_ref) > log_spacing_rtol*abs(dlog_ref)) then
            error stop "hadronic species transport: gamma grid must be approximately log-spaced."
        end if
    end do
    hadronic_species_log_spacing = dlog_ref
end function hadronic_species_log_spacing

end module hadronic_species_transport_kernel
