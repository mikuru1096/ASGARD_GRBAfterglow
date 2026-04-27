!f2py: skip
module hadronic_pp_kernel
    use constants
    use hadronic_common, only: hadronic_validate_log_grid
    implicit none
    private

    real(8), parameter :: proton_mass_gev = 0.9382720813d0
    real(8), parameter :: pion0_mass_gev = 0.1349768d0
    real(8), parameter :: mbarn_to_cm2 = 1.0d-27

    public :: hadronic_pp_delta_operator
    public :: hadronic_pp_sigma_inelastic_kelner2006
    public :: hadronic_pp_threshold_kinetic_energy_gev

contains

! pp δ-近似算子：基于Kelner+2006非弹性截面，计算pp碰撞产生的γ、ν和e±源项。
subroutine hadronic_pp_delta_operator(num_p,proton_energy_gev,proton_density_per_gev, &
                                      target_proton_density_cm3,num_gamma,gamma_energy_gev, &
                                      num_nu,neutrino_energy_gev,num_pair,pair_energy_gev, &
                                      gamma_rate_per_gev,neutrino_rate_per_gev, &
                                      pair_rate_per_gev,proton_loss_rate, &
                                      kappa_inelastic,pion_energy_fraction, &
                                      neutral_pion_fraction)
    integer, intent(in) :: num_p,num_gamma,num_nu,num_pair
    real(8), intent(in) :: proton_energy_gev(num_p),proton_density_per_gev(num_p)
    real(8), intent(in) :: target_proton_density_cm3
    real(8), intent(in) :: gamma_energy_gev(num_gamma),neutrino_energy_gev(num_nu)
    real(8), intent(in) :: pair_energy_gev(num_pair)
    real(8), intent(out) :: gamma_rate_per_gev(num_gamma),neutrino_rate_per_gev(num_nu)
    real(8), intent(out) :: pair_rate_per_gev(num_pair),proton_loss_rate(num_p)
    real(8), intent(in), optional :: kappa_inelastic,pion_energy_fraction
    real(8), intent(in), optional :: neutral_pion_fraction
    real(8) :: sigma_inel(num_p),collision_rate(num_p),parent_rate(num_p)
    real(8) :: kappa_loc,pion_frac_loc,neutral_frac_loc,charged_frac_loc
    real(8) :: x_gamma,x_nu,x_pair

    call hadronic_pp_validate_grid(num_p,proton_energy_gev,"proton_energy_gev")
    call hadronic_pp_validate_density_grid(num_p,proton_energy_gev,proton_density_per_gev)
    call hadronic_pp_validate_grid(num_gamma,gamma_energy_gev,"gamma_energy_gev")
    call hadronic_pp_validate_grid(num_nu,neutrino_energy_gev,"neutrino_energy_gev")
    call hadronic_pp_validate_grid(num_pair,pair_energy_gev,"pair_energy_gev")

    if (target_proton_density_cm3 < zero) then
        error stop "hadronic_pp_delta_operator: target_proton_density_cm3 must be non-negative."
    end if

    kappa_loc = 0.5d0
    if (present(kappa_inelastic)) kappa_loc = kappa_inelastic
    if (kappa_loc <= zero .or. kappa_loc > one) then
        error stop "hadronic_pp_delta_operator: kappa_inelastic must be in (0, 1]."
    end if

    pion_frac_loc = 0.17d0
    if (present(pion_energy_fraction)) pion_frac_loc = pion_energy_fraction
    if (pion_frac_loc <= zero .or. pion_frac_loc >= one) then
        error stop "hadronic_pp_delta_operator: pion_energy_fraction must be in (0, 1)."
    end if

    neutral_frac_loc = one/3d0
    if (present(neutral_pion_fraction)) neutral_frac_loc = neutral_pion_fraction
    if (neutral_frac_loc < zero .or. neutral_frac_loc > one) then
        error stop "hadronic_pp_delta_operator: neutral_pion_fraction must be in [0, 1]."
    end if

    charged_frac_loc = one - neutral_frac_loc
    x_gamma = 0.5d0*pion_frac_loc
    x_nu = 0.25d0*pion_frac_loc
    x_pair = x_nu

    call hadronic_pp_sigma_inelastic_kelner2006(num_p,proton_energy_gev,sigma_inel)
    collision_rate = target_proton_density_cm3*Para_c*sigma_inel
    parent_rate = collision_rate*proton_density_per_gev

    proton_loss_rate = -kappa_loc*parent_rate

    call hadronic_pp_delta_secondary_source(num_gamma,gamma_energy_gev,num_p,proton_energy_gev,parent_rate, &
                                            x_gamma,two*neutral_frac_loc,gamma_rate_per_gev)
    call hadronic_pp_delta_secondary_source(num_nu,neutrino_energy_gev,num_p,proton_energy_gev,parent_rate, &
                                            x_nu,3d0*charged_frac_loc,neutrino_rate_per_gev)
    call hadronic_pp_delta_secondary_source(num_pair,pair_energy_gev,num_p,proton_energy_gev,parent_rate, &
                                            x_pair,charged_frac_loc,pair_rate_per_gev)
end subroutine hadronic_pp_delta_operator

! Kelner+2006 pp非弹性截面参数化：σ_inel(E_p) = (34.3 + 1.88L + 0.25L²)·(1-(E_th/E_k)⁴)² mb。
subroutine hadronic_pp_sigma_inelastic_kelner2006(num_p,proton_energy_gev,sigma_inel_cm2)
    integer, intent(in) :: num_p
    real(8), intent(in) :: proton_energy_gev(num_p)
    real(8), intent(out) :: sigma_inel_cm2(num_p)
    integer :: i_p
    real(8) :: kinetic_energy_gev,threshold_gev,ratio,log_term,cutoff

    call hadronic_pp_validate_grid(num_p,proton_energy_gev,"proton_energy_gev")
    if (minval(proton_energy_gev) <= proton_mass_gev) then
        error stop "hadronic_pp_sigma_inelastic_kelner2006: proton_energy_gev must exceed proton rest energy."
    end if

    threshold_gev = hadronic_pp_threshold_kinetic_energy_gev()
    sigma_inel_cm2 = zero
    do i_p=1,num_p
        kinetic_energy_gev = proton_energy_gev(i_p) - proton_mass_gev
        if (kinetic_energy_gev <= threshold_gev) cycle
        ratio = kinetic_energy_gev/threshold_gev
        log_term = dlog(ratio)
        cutoff = one - (threshold_gev/kinetic_energy_gev)**4
        sigma_inel_cm2(i_p) = (34.3d0 + 1.88d0*log_term + 0.25d0*log_term*log_term) * &
                              cutoff*cutoff*mbarn_to_cm2
    end do
end subroutine hadronic_pp_sigma_inelastic_kelner2006

! pp反应阈值动能：E_th = 2m_π⁰ + m_π⁰²/(2m_p) ≈ 0.28 GeV。
real(8) function hadronic_pp_threshold_kinetic_energy_gev()
    hadronic_pp_threshold_kinetic_energy_gev = two*pion0_mass_gev + &
                                               pion0_mass_gev*pion0_mass_gev/(two*proton_mass_gev)
end function hadronic_pp_threshold_kinetic_energy_gev

! δ-近似次级粒子源项：次级粒子能量 = 能量份额 × 母粒子能量，多重数加权。
subroutine hadronic_pp_delta_secondary_source(num_secondary,secondary_energy_gev,num_parent, &
                                              parent_energy_gev,parent_rate_per_gev, &
                                              energy_fraction,multiplicity,secondary_rate_per_gev)
    integer, intent(in) :: num_secondary,num_parent
    real(8), intent(in) :: secondary_energy_gev(num_secondary),parent_energy_gev(num_parent)
    real(8), intent(in) :: parent_rate_per_gev(num_parent),energy_fraction,multiplicity
    real(8), intent(out) :: secondary_rate_per_gev(num_secondary)
    real(8) :: parent_eval_gev(num_secondary),sampled_parent(num_secondary)

    if (energy_fraction <= zero) then
        error stop "hadronic_pp_delta_secondary_source: energy_fraction must be positive."
    end if
    if (multiplicity < zero) then
        error stop "hadronic_pp_delta_secondary_source: multiplicity must be non-negative."
    end if

    parent_eval_gev = secondary_energy_gev/energy_fraction
    call hadronic_pp_loglog_interp_positive(num_parent,parent_energy_gev,parent_rate_per_gev, &
                                            num_secondary,parent_eval_gev,sampled_parent)
    secondary_rate_per_gev = (multiplicity/energy_fraction)*sampled_parent
end subroutine hadronic_pp_delta_secondary_source

! 双对数空间中正值的线性插值，仅使用y>0的区间。
subroutine hadronic_pp_loglog_interp_positive(num_x,x,y,num_x_new,x_new,y_new)
    integer, intent(in) :: num_x,num_x_new
    real(8), intent(in) :: x(num_x),y(num_x),x_new(num_x_new)
    real(8), intent(out) :: y_new(num_x_new)
    integer :: positive_idx(num_x),num_positive
    real(8) :: xp(num_x),yp(num_x)
    integer :: i,ipos

    call hadronic_pp_validate_grid(num_x,x,"interpolation x")
    do i=1,num_x_new
        if (x_new(i) <= zero) then
            error stop "hadronic_pp_loglog_interp_positive: x_new must be strictly positive."
        end if
    end do

    num_positive = 0
    do i=1,num_x
        if (y(i) > zero) then
            num_positive = num_positive + 1
            positive_idx(num_positive) = i
        end if
    end do

    y_new = zero
    if (num_positive < 2) return

    do i=1,num_positive
        xp(i) = x(positive_idx(i))
        yp(i) = y(positive_idx(i))
    end do

    do i=1,num_x_new
        if (x_new(i) < xp(1) .or. x_new(i) > xp(num_positive)) cycle
        ipos = hadronic_pp_upper_bracket(num_positive,xp,x_new(i))
        y_new(i) = dexp(hadronic_pp_loglog_linear_eval(xp(ipos),xp(ipos+1),yp(ipos),yp(ipos+1),x_new(i)))
    end do
end subroutine hadronic_pp_loglog_interp_positive

! 二分查找单调递增数组x中x_eval所在的区间上界索引。
integer function hadronic_pp_upper_bracket(num_x,x,x_eval)
    integer, intent(in) :: num_x
    real(8), intent(in) :: x(num_x),x_eval
    integer :: ilo,ihi,imid

    if (x_eval <= x(1)) then
        hadronic_pp_upper_bracket = 1
        return
    end if
    if (x_eval >= x(num_x)) then
        hadronic_pp_upper_bracket = num_x - 1
        return
    end if

    ilo = 1
    ihi = num_x
    do while (ihi - ilo > 1)
        imid = (ilo + ihi)/2
        if (x(imid) <= x_eval) then
            ilo = imid
        else
            ihi = imid
        end if
    end do
    hadronic_pp_upper_bracket = ilo
end function hadronic_pp_upper_bracket

! 双对数空间中单次线性插值求值。
real(8) function hadronic_pp_loglog_linear_eval(x0,x1,y0,y1,x_eval)
    real(8), intent(in) :: x0,x1,y0,y1,x_eval
    real(8) :: lx0,lx1,ly0,ly1,lx_eval,frac

    lx0 = dlog(x0)
    lx1 = dlog(x1)
    ly0 = dlog(y0)
    ly1 = dlog(y1)
    lx_eval = dlog(x_eval)
    frac = (lx_eval - lx0)/(lx1 - lx0)
    hadronic_pp_loglog_linear_eval = ly0 + frac*(ly1 - ly0)
end function hadronic_pp_loglog_linear_eval

! 验证密度网格与能量网格大小一致。
subroutine hadronic_pp_validate_density_grid(num_p,proton_energy_gev,proton_density_per_gev)
    integer, intent(in) :: num_p
    real(8), intent(in) :: proton_energy_gev(num_p),proton_density_per_gev(num_p)

    call hadronic_pp_validate_grid(num_p,proton_energy_gev,"proton_energy_gev")
    if (size(proton_density_per_gev) /= num_p) then
        error stop "hadronic_pp_validate_density_grid: proton_density_per_gev shape mismatch."
    end if
end subroutine hadronic_pp_validate_density_grid

! 验证网格为非空且对数均匀（包装hadronic_validate_log_grid）。
subroutine hadronic_pp_validate_grid(num_x,x,name)
    integer, intent(in) :: num_x
    real(8), intent(in) :: x(num_x)
    character(*), intent(in) :: name

    if (len_trim(name) <= 0) error stop "hadronic_pp_validate_grid: internal name must be non-empty."
    call hadronic_validate_log_grid(num_x,x,name)
end subroutine hadronic_pp_validate_grid

end module hadronic_pp_kernel
