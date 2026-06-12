!f2py: skip
module hadronic_acceleration_kernel
    use constants
    use hadronic_common, only: hadronic_proton_mass_gev, hadronic_electron_mass_gev, &
        hadronic_neutron_mass_gev, hadronic_pion_charged_mass_gev, hadronic_muon_mass_gev
    implicit none
    private

    real(8), parameter :: gev_to_gram = 1d9 * Para_eV2erg / (Para_c * Para_c)

    public :: hadronic_species_properties
    public :: hadronic_acceleration_timescale_s
    public :: hadronic_synchrotron_cooling_timescale_s
    public :: hadronic_external_photon_cooling_timescale_s
    public :: hadronic_species_injection_operator
    public :: hadronic_estimate_max_gamma
    public :: hadronic_acceleration_operator

contains

! 查询粒子种类（质子、中子、π介子、μ子）的质量、电荷数等基本属性。
subroutine hadronic_species_properties(species,mass_gev,charge_number,mass_g,abs_charge_esu)
    character(len=*), intent(in) :: species
    real(8), intent(out) :: mass_gev,mass_g,abs_charge_esu
    integer, intent(out) :: charge_number

    select case (trim(species))
    case ("proton")
        mass_gev = hadronic_proton_mass_gev
        charge_number = 1
    case ("neutron")
        mass_gev = hadronic_neutron_mass_gev
        charge_number = 0
    case ("pion_plus")
        mass_gev = hadronic_pion_charged_mass_gev
        charge_number = 1
    case ("pion_minus")
        mass_gev = hadronic_pion_charged_mass_gev
        charge_number = -1
    case ("muon_plus")
        mass_gev = hadronic_muon_mass_gev
        charge_number = 1
    case ("muon_minus")
        mass_gev = hadronic_muon_mass_gev
        charge_number = -1
    case default
        error stop "hadronic_species_properties: unsupported species."
    end select

    mass_g = mass_gev*gev_to_gram
    abs_charge_esu = dabs(dble(charge_number))*Para_e
end subroutine hadronic_species_properties

! 计算粒子在磁场中的费米加速时标 t_acc = eta_acc * gamma * m * c / (|q| * B)。
subroutine hadronic_acceleration_timescale_s(num_gamma,species,gamma,b_field_g,eta_acc,t_acc)
    integer, intent(in) :: num_gamma
    character(len=*), intent(in) :: species
    real(8), intent(in) :: gamma(num_gamma),b_field_g,eta_acc
    real(8), intent(out) :: t_acc(num_gamma)
    integer :: i_gam,charge_number
    real(8) :: mass_gev,mass_g,abs_charge_esu

    call hadronic_species_properties(species,mass_gev,charge_number,mass_g,abs_charge_esu)
    if (b_field_g <= zero) error stop "hadronic_acceleration_timescale_s: b_field_g must be > 0."
    if (eta_acc <= zero) error stop "hadronic_acceleration_timescale_s: eta_acc must be > 0."
    if (abs_charge_esu <= zero) then
        error stop "hadronic_acceleration_timescale_s: species has zero charge."
    end if

    do i_gam=1,num_gamma
        if (gamma(i_gam) <= zero) error stop "hadronic_acceleration_timescale_s: gamma must be > 0."
        t_acc(i_gam) = eta_acc*gamma(i_gam)*mass_g*Para_c/(abs_charge_esu*b_field_g)
    end do
end subroutine hadronic_acceleration_timescale_s

! 计算粒子在磁场中的同步辐射冷却时标。
subroutine hadronic_synchrotron_cooling_timescale_s(num_gamma,species,gamma,b_field_g,t_syn)
    integer, intent(in) :: num_gamma
    character(len=*), intent(in) :: species
    real(8), intent(in) :: gamma(num_gamma),b_field_g
    real(8), intent(out) :: t_syn(num_gamma)
    integer :: i_gam,charge_number
    real(8) :: mass_gev,mass_g,abs_charge_esu

    call hadronic_species_properties(species,mass_gev,charge_number,mass_g,abs_charge_esu)
    if (b_field_g <= zero) then
        error stop "hadronic_synchrotron_cooling_timescale_s: b_field_g must be > 0."
    end if

    do i_gam=1,num_gamma
        if (gamma(i_gam) <= zero) then
            error stop "hadronic_synchrotron_cooling_timescale_s: gamma must be > 0."
        end if
        t_syn(i_gam) = 6d0*pi*(mass_g**3)*Para_c / &
                       (Para_SigmaT*(Para_m_e**2)*(b_field_g**2)*gamma(i_gam))
    end do
end subroutine hadronic_synchrotron_cooling_timescale_s

! 根据给定的外部冷却率计算外部光子场冷却时标 t_ext = gamma / cooling_rate。
subroutine hadronic_external_photon_cooling_timescale_s(num_gamma,gamma,cooling_rate,t_ext)
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: gamma(num_gamma),cooling_rate(num_gamma)
    real(8), intent(out) :: t_ext(num_gamma)
    integer :: i_gam

    do i_gam=1,num_gamma
        if (gamma(i_gam) <= zero) then
            error stop "hadronic_external_photon_cooling_timescale_s: gamma must be > 0."
        end if
        if (cooling_rate(i_gam) <= zero) then
            error stop "hadronic_external_photon_cooling_timescale_s: cooling_rate must be > 0."
        end if
        t_ext(i_gam) = gamma(i_gam)/cooling_rate(i_gam)
    end do
end subroutine hadronic_external_photon_cooling_timescale_s

! 计算粒子注入源项：幂律谱分布 × 指数截断，按光度归一化。
subroutine hadronic_species_injection_operator(num_gamma,gamma,species,luminosity_erg_s, &
                                               spectral_index,gamma_min,gamma_max,gamma_cut, &
                                               has_gamma_cut,q_inj)
    integer, intent(in) :: num_gamma
    character(len=*), intent(in) :: species
    real(8), intent(in) :: gamma(num_gamma),luminosity_erg_s,spectral_index
    real(8), intent(in) :: gamma_min,gamma_max,gamma_cut
    logical, intent(in) :: has_gamma_cut
    real(8), intent(out) :: q_inj(num_gamma)
    integer :: i_gam,charge_number
    real(8) :: mass_gev,mass_g,abs_charge_esu,mass_energy_erg,norm,q0
    real(8) :: shape(num_gamma),integrand(num_gamma)

    call hadronic_species_properties(species,mass_gev,charge_number,mass_g,abs_charge_esu)
    if (gamma_min <= zero .or. gamma_max <= zero .or. gamma_min >= gamma_max) then
        error stop "hadronic_species_injection_operator: require 0 < gamma_min < gamma_max."
    end if
    if (luminosity_erg_s <= zero) then
        error stop "hadronic_species_injection_operator: luminosity_erg_s must be > 0."
    end if
    if (has_gamma_cut .and. gamma_cut <= zero) then
        error stop "hadronic_species_injection_operator: gamma_cut must be > 0."
    end if

    do i_gam=1,num_gamma
        if (gamma(i_gam) <= zero) error stop "hadronic_species_injection_operator: gamma must be > 0."
        shape(i_gam) = zero
        if (gamma(i_gam) >= gamma_min .and. gamma(i_gam) <= gamma_max) then
            shape(i_gam) = gamma(i_gam)**(-spectral_index)
            if (has_gamma_cut) then
                shape(i_gam) = shape(i_gam)*dexp(-gamma(i_gam)/gamma_cut)
            end if
        end if
        integrand(i_gam) = shape(i_gam)*gamma(i_gam)
    end do

    norm = hadronic_trapezoid(num_gamma,gamma,integrand)
    if (norm <= zero) then
        error stop "hadronic_species_injection_operator: normalization integral must be > 0."
    end if

    mass_energy_erg = mass_gev*Para_GeV2erg
    q0 = luminosity_erg_s/(mass_energy_erg*norm)
    q_inj = q0*shape
end subroutine hadronic_species_injection_operator

! 通过平衡加速、同步冷却、动力学和外部冷却时标来估计粒子的最大洛伦兹因子。
subroutine hadronic_estimate_max_gamma(species,b_field_g,radius_cm,gamma_bulk,eta_acc, &
                                       num_gamma_scan,gamma_scan,external_cooling_rate, &
                                       has_external_cooling,gamma_max,gamma_dyn,gamma_syn, &
                                       gamma_ext,has_gamma_ext)
    character(len=*), intent(in) :: species
    real(8), intent(in) :: b_field_g,radius_cm,gamma_bulk,eta_acc
    integer, intent(in) :: num_gamma_scan
    real(8), intent(in) :: gamma_scan(num_gamma_scan),external_cooling_rate(num_gamma_scan)
    logical, intent(in) :: has_external_cooling
    real(8), intent(out) :: gamma_max,gamma_dyn,gamma_syn,gamma_ext
    logical, intent(out) :: has_gamma_ext
    integer :: charge_number,i_cross
    real(8) :: mass_gev,mass_g,abs_charge_esu
    real(8) :: t_acc(num_gamma_scan),t_ext(num_gamma_scan),ratio_lo,ratio_hi
    real(8) :: x0,x1,y0,y1,x_star

    call hadronic_species_properties(species,mass_gev,charge_number,mass_g,abs_charge_esu)
    if (abs_charge_esu <= zero) then
        error stop "hadronic_estimate_max_gamma: species has zero charge."
    end if
    if (b_field_g <= zero .or. radius_cm <= zero .or. gamma_bulk <= zero .or. eta_acc <= zero) then
        error stop "hadronic_estimate_max_gamma: physical inputs must be > 0."
    end if

    call initialize_acceleration_limits

    if (.not. has_external_cooling) return

    call build_external_cooling_timescales
    i_cross = find_external_cooling_crossing()
    if (i_cross == 0) return

    call apply_external_cooling_limit(i_cross)

contains

    subroutine initialize_acceleration_limits
        real(8) :: t_dyn

        t_dyn = radius_cm/(gamma_bulk*Para_c)
        gamma_dyn = abs_charge_esu*b_field_g*t_dyn/(eta_acc*mass_g*Para_c)
        gamma_syn = dsqrt(6d0*pi*abs_charge_esu*(mass_g**2) / &
                    (eta_acc*Para_SigmaT*(Para_m_e**2)*b_field_g))

        has_gamma_ext = .false.
        gamma_ext = zero
        gamma_max = dmin1(gamma_dyn,gamma_syn)
    end subroutine initialize_acceleration_limits

    subroutine build_external_cooling_timescales
        call hadronic_acceleration_timescale_s(num_gamma_scan,species,gamma_scan,b_field_g,eta_acc,t_acc)
        call hadronic_external_photon_cooling_timescale_s(num_gamma_scan,gamma_scan,external_cooling_rate,t_ext)
    end subroutine build_external_cooling_timescales

    integer function find_external_cooling_crossing()
        integer :: i_gam

        find_external_cooling_crossing = 0
        do i_gam=1,num_gamma_scan-1
            ratio_lo = t_acc(i_gam)-t_ext(i_gam)
            ratio_hi = t_acc(i_gam+1)-t_ext(i_gam+1)
            if (ratio_lo == zero .or. ratio_lo*ratio_hi <= zero .or. ratio_hi == zero) then
                find_external_cooling_crossing = i_gam
                exit
            end if
        end do
    end function find_external_cooling_crossing

    subroutine apply_external_cooling_limit(i_cross)
        integer, intent(in) :: i_cross

        if (t_acc(i_cross) == t_ext(i_cross)) then
            gamma_ext = gamma_scan(i_cross)
        else if (t_acc(i_cross+1) == t_ext(i_cross+1)) then
            gamma_ext = gamma_scan(i_cross+1)
        else
            x0 = dlog(gamma_scan(i_cross))
            x1 = dlog(gamma_scan(i_cross+1))
            y0 = dlog(t_acc(i_cross)/t_ext(i_cross))
            y1 = dlog(t_acc(i_cross+1)/t_ext(i_cross+1))
            x_star = x0-y0*(x1-x0)/(y1-y0)
            gamma_ext = dexp(x_star)
        end if

        has_gamma_ext = .true.
        gamma_max = dmin1(gamma_max,gamma_ext)
    end subroutine apply_external_cooling_limit
end subroutine hadronic_estimate_max_gamma

! 统一调用加速时标、同步冷却、注入源项和最大能量估计，完成粒子加速算子的完整计算。
subroutine hadronic_acceleration_operator(num_gamma,species,gamma,b_field_g,eta_acc, &
                                          luminosity_erg_s,spectral_index,gamma_min, &
                                          gamma_max_inj,gamma_cut,has_gamma_cut, &
                                          radius_cm,gamma_bulk,num_gamma_scan,gamma_scan, &
                                          external_cooling_rate,has_external_cooling, &
                                          t_acc,t_syn,q_inj,gamma_max,gamma_dyn,gamma_syn, &
                                          gamma_ext,has_gamma_ext)
    integer, intent(in) :: num_gamma,num_gamma_scan
    character(len=*), intent(in) :: species
    real(8), intent(in) :: gamma(num_gamma),b_field_g,eta_acc,luminosity_erg_s,spectral_index
    real(8), intent(in) :: gamma_min,gamma_max_inj,gamma_cut,radius_cm,gamma_bulk
    real(8), intent(in) :: gamma_scan(num_gamma_scan),external_cooling_rate(num_gamma_scan)
    logical, intent(in) :: has_gamma_cut,has_external_cooling
    real(8), intent(out) :: t_acc(num_gamma),t_syn(num_gamma),q_inj(num_gamma)
    real(8), intent(out) :: gamma_max,gamma_dyn,gamma_syn,gamma_ext
    logical, intent(out) :: has_gamma_ext

    call hadronic_acceleration_timescale_s(num_gamma,species,gamma,b_field_g,eta_acc,t_acc)
    call hadronic_synchrotron_cooling_timescale_s(num_gamma,species,gamma,b_field_g,t_syn)
    call hadronic_species_injection_operator(num_gamma,gamma,species,luminosity_erg_s, &
                                             spectral_index,gamma_min,gamma_max_inj, &
                                             gamma_cut,has_gamma_cut,q_inj)
    call hadronic_estimate_max_gamma(species,b_field_g,radius_cm,gamma_bulk,eta_acc, &
                                     num_gamma_scan,gamma_scan,external_cooling_rate, &
                                     has_external_cooling,gamma_max,gamma_dyn,gamma_syn, &
                                     gamma_ext,has_gamma_ext)
end subroutine hadronic_acceleration_operator

! 梯形法则数值积分。
real(8) function hadronic_trapezoid(num_x,x,y)
    integer, intent(in) :: num_x
    real(8), intent(in) :: x(num_x),y(num_x)
    integer :: i_x

    hadronic_trapezoid = zero
    if (num_x <= 1) return
    do i_x=1,num_x-1
        hadronic_trapezoid = hadronic_trapezoid + 5d-1*(y(i_x)+y(i_x+1))*(x(i_x+1)-x(i_x))
    end do
end function hadronic_trapezoid

end module hadronic_acceleration_kernel
