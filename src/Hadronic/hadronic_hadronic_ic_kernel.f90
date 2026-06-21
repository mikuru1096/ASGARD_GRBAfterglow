!f2py: skip
module hadronic_hadronic_ic_kernel
    use constants
    use hadronic_common, only: hadronic_validate_log_grid
    implicit none
    private

    real(8), parameter :: am3_c_cgs = Para_c
    real(8), parameter :: am3_sigma_t_cgs = Para_sigmaT
    real(8), parameter :: am3_mass_electron_gev = Para_m_energy*Para_erg2eV*1d-9
    real(8), parameter :: am3_mass_proton_gev = Para_m_p_E*Para_erg2eV*1d-9
    real(8), parameter :: am3_mass_pion_charged_gev = 1.396d-1
    real(8), parameter :: am3_mass_muon_gev = 1.0566d-1

    public :: hadronic_hadronic_ic_initialize_kernel
    public :: hadronic_hadronic_ic_operator_from_kernel

contains

! 初始化强子IC计算核：验证网格一致性，为质子/π介子/μ子预计算kernel映射索引。
subroutine hadronic_hadronic_ic_initialize_kernel(num_had,hadron_energy_gev,num_ph,photon_energy_gev, &
                                                  ind_min_energy_pho_hadgrid,dln_energy, &
                                                  delta_e_p,jmax_p,delta_e_pi,jmax_pi,delta_e_mu,jmax_mu)
    integer, intent(in) :: num_had,num_ph,ind_min_energy_pho_hadgrid
    real(8), intent(in) :: hadron_energy_gev(num_had),photon_energy_gev(num_ph)
    real(8), intent(out) :: dln_energy
    integer, intent(out) :: delta_e_p(num_had),jmax_p(num_had),delta_e_pi(num_had),jmax_pi(num_had)
    integer, intent(out) :: delta_e_mu(num_had),jmax_mu(num_had)
    real(8) :: dln_had,dln_ph

    dln_had = hadronic_hadronic_ic_log_spacing(num_had,hadron_energy_gev)
    dln_ph = hadronic_hadronic_ic_log_spacing(num_ph,photon_energy_gev)
    if (dabs(dln_had-dln_ph) > dmax1(1d-12,1d-10*dabs(dln_had))) then
        error stop "hadronic IC requires hadron/photon grids with the same logarithmic spacing."
    end if

    dln_energy = dln_had
    call hadronic_hadronic_ic_build_species_kernel(num_had,hadron_energy_gev,dln_energy,num_ph, &
                                                   ind_min_energy_pho_hadgrid,am3_mass_proton_gev,delta_e_p,jmax_p)
    call hadronic_hadronic_ic_build_species_kernel(num_had,hadron_energy_gev,dln_energy,num_ph, &
                                                   ind_min_energy_pho_hadgrid,am3_mass_pion_charged_gev,delta_e_pi,jmax_pi)
    call hadronic_hadronic_ic_build_species_kernel(num_had,hadron_energy_gev,dln_energy,num_ph, &
                                                   ind_min_energy_pho_hadgrid,am3_mass_muon_gev,delta_e_mu,jmax_mu)
end subroutine hadronic_hadronic_ic_initialize_kernel

! 使用预计算的kernel索引直接计算强子IC（跳过初始化步骤）。
subroutine hadronic_hadronic_ic_operator_from_kernel(num_had,num_ph,photons_on_had_grid_per_gev, &
                                                     protons_per_gev,pion_plus_per_gev,pion_minus_per_gev, &
                                                     muon_minus_left_per_gev,muon_minus_right_per_gev, &
                                                     muon_plus_left_per_gev,muon_plus_right_per_gev, &
                                                     dln_energy,delta_e_p,jmax_p,delta_e_pi,jmax_pi, &
                                                     delta_e_mu,jmax_mu,epsilon_p_ic,epsilon_pi_ic,epsilon_mu_ic, &
                                                     coeff_p_cgs,coeff_pi_cgs,coeff_mu_cgs)
    integer, intent(in) :: num_had,num_ph
    real(8), intent(in) :: photons_on_had_grid_per_gev(num_ph)
    real(8), intent(in) :: protons_per_gev(num_had),pion_plus_per_gev(num_had),pion_minus_per_gev(num_had)
    real(8), intent(in) :: muon_minus_left_per_gev(num_had),muon_minus_right_per_gev(num_had)
    real(8), intent(in) :: muon_plus_left_per_gev(num_had),muon_plus_right_per_gev(num_had),dln_energy
    integer, intent(in) :: delta_e_p(num_had),jmax_p(num_had),delta_e_pi(num_had),jmax_pi(num_had)
    integer, intent(in) :: delta_e_mu(num_had),jmax_mu(num_had)
    real(8), intent(out) :: epsilon_p_ic(num_ph),epsilon_pi_ic(num_ph),epsilon_mu_ic(num_ph)
    real(8), intent(out) :: coeff_p_cgs,coeff_pi_cgs,coeff_mu_cgs

    call hadronic_hadronic_ic_apply_kernel(num_had,num_ph,photons_on_had_grid_per_gev,protons_per_gev, &
                                           pion_plus_per_gev,pion_minus_per_gev,muon_minus_left_per_gev, &
                                           muon_minus_right_per_gev,muon_plus_left_per_gev,muon_plus_right_per_gev, &
                                           dln_energy,delta_e_p,jmax_p,delta_e_pi,jmax_pi,delta_e_mu,jmax_mu, &
                                           epsilon_p_ic,epsilon_pi_ic,epsilon_mu_ic,coeff_p_cgs, &
                                           coeff_pi_cgs,coeff_mu_cgs)
end subroutine hadronic_hadronic_ic_operator_from_kernel

! 应用IC kernel：分别对质子、π介子和μ子三个通道计算IC冷却率。
subroutine hadronic_hadronic_ic_apply_kernel(num_had,num_ph,photons_on_had_grid_per_gev,protons_per_gev, &
                                             pion_plus_per_gev,pion_minus_per_gev,muon_minus_left_per_gev, &
                                             muon_minus_right_per_gev,muon_plus_left_per_gev,muon_plus_right_per_gev, &
                                             dln_energy,delta_e_p,jmax_p,delta_e_pi,jmax_pi,delta_e_mu,jmax_mu, &
                                             epsilon_p_ic,epsilon_pi_ic,epsilon_mu_ic,coeff_p_cgs, &
                                             coeff_pi_cgs,coeff_mu_cgs)
    integer, intent(in) :: num_had,num_ph
    real(8), intent(in) :: photons_on_had_grid_per_gev(num_ph)
    real(8), intent(in) :: protons_per_gev(num_had),pion_plus_per_gev(num_had),pion_minus_per_gev(num_had)
    real(8), intent(in) :: muon_minus_left_per_gev(num_had),muon_minus_right_per_gev(num_had)
    real(8), intent(in) :: muon_plus_left_per_gev(num_had),muon_plus_right_per_gev(num_had),dln_energy
    integer, intent(in) :: delta_e_p(num_had),jmax_p(num_had),delta_e_pi(num_had),jmax_pi(num_had)
    integer, intent(in) :: delta_e_mu(num_had),jmax_mu(num_had)
    real(8), intent(out) :: epsilon_p_ic(num_ph),epsilon_pi_ic(num_ph),epsilon_mu_ic(num_ph)
    real(8), intent(out) :: coeff_p_cgs,coeff_pi_cgs,coeff_mu_cgs
    real(8) :: hadron_sum(num_had)

    call apply_hadronic_ic_species(am3_mass_proton_gev,protons_per_gev,delta_e_p,jmax_p,coeff_p_cgs,epsilon_p_ic)

    call sum_charged_pion_density(hadron_sum)
    call apply_hadronic_ic_species(am3_mass_pion_charged_gev,hadron_sum,delta_e_pi,jmax_pi,coeff_pi_cgs,epsilon_pi_ic)

    call sum_muon_helicity_density(hadron_sum)
    call apply_hadronic_ic_species(am3_mass_muon_gev,hadron_sum,delta_e_mu,jmax_mu,coeff_mu_cgs,epsilon_mu_ic)

contains

    subroutine apply_hadronic_ic_species(mass_gev,hadron_density_per_gev,delta_e,jmax,coeff_cgs,epsilon_ic)
        real(8), intent(in) :: mass_gev,hadron_density_per_gev(num_had)
        integer, intent(in) :: delta_e(num_had),jmax(num_had)
        real(8), intent(out) :: coeff_cgs,epsilon_ic(num_ph)

        coeff_cgs = hadronic_hadronic_ic_coeff(mass_gev)
        call hadronic_hadronic_ic_compute_channel(num_ph,photons_on_had_grid_per_gev,num_had,hadron_density_per_gev, &
                                                  delta_e,jmax,dln_energy,coeff_cgs,epsilon_ic)
    end subroutine apply_hadronic_ic_species

    subroutine sum_charged_pion_density(hadron_sum)
        real(8), intent(out) :: hadron_sum(num_had)

        hadron_sum = pion_plus_per_gev + pion_minus_per_gev
    end subroutine sum_charged_pion_density

    subroutine sum_muon_helicity_density(hadron_sum)
        real(8), intent(out) :: hadron_sum(num_had)

        hadron_sum = muon_minus_left_per_gev + muon_minus_right_per_gev + muon_plus_left_per_gev + muon_plus_right_per_gev
    end subroutine sum_muon_helicity_density
end subroutine hadronic_hadronic_ic_apply_kernel

! 为给定粒子种类构建IC映射kernel：计算delta_e（能量偏移）和jmax（最大光子索引）。
subroutine hadronic_hadronic_ic_build_species_kernel(num_had,hadron_energy_gev,dln_energy,num_ph, &
                                                     ind_min_energy_pho_hadgrid,mass_gev,delta_e,jmax)
    integer, intent(in) :: num_had,num_ph,ind_min_energy_pho_hadgrid
    real(8), intent(in) :: hadron_energy_gev(num_had),dln_energy,mass_gev
    integer, intent(out) :: delta_e(num_had),jmax(num_had)
    integer :: i,candidate,jmax1
    real(8) :: gamma

    do i=1,num_had
        if (hadron_energy_gev(i) <= zero) then
            error stop "hadronic IC requires positive hadron energies."
        end if
        gamma = hadron_energy_gev(i)/mass_gev
        delta_e(i) = int(dlog(gamma*gamma)/dln_energy)
        jmax1 = int(dlog(0.5d0*mass_gev/gamma)/dln_energy) + ind_min_energy_pho_hadgrid
        candidate = jmax1 + delta_e(i)
        if (candidate > num_ph) then
            jmax(i) = num_ph
        else
            jmax(i) = candidate
        end if
    end do
end subroutine hadronic_hadronic_ic_build_species_kernel

! 计算单个强子种类的IC通道：对光子能格进行卷积求和。
subroutine hadronic_hadronic_ic_compute_channel(num_ph,photons_on_had_grid_per_gev,num_had,hadron_density_per_gev, &
                                                delta_e,jmax,dln_energy,coeff_cgs,epsilon_ic)
    integer, intent(in) :: num_ph,num_had,delta_e(num_had),jmax(num_had)
    real(8), intent(in) :: photons_on_had_grid_per_gev(num_ph),hadron_density_per_gev(num_had),dln_energy,coeff_cgs
    real(8), intent(out) :: epsilon_ic(num_ph)
    integer :: i,j,src_idx,j0
    real(8) :: z

    epsilon_ic = zero
    !$OMP PARALLEL DO if(num_ph*num_had >= 8192) schedule(static) private(i,j,src_idx,j0,z)
    do j=1,num_ph
        j0 = j - 1
        z = zero
        do i=1,num_had
            if (j0 < delta_e(i) .or. j0 > jmax(i)) cycle
            src_idx = j - delta_e(i)
            if (src_idx < 1 .or. src_idx > num_ph) then
                error stop "hadronic IC kernel maps to an out-of-grid photon source index."
            end if
            z = z + photons_on_had_grid_per_gev(src_idx)*hadron_density_per_gev(i)
        end do
        epsilon_ic(j) = z*dln_energy*coeff_cgs
    end do
    !$OMP END PARALLEL DO
end subroutine hadronic_hadronic_ic_compute_channel

! IC前因子系数：σ_T * c / (mass_ratio)^2。
real(8) function hadronic_hadronic_ic_coeff(mass_gev)
    real(8), intent(in) :: mass_gev
    real(8) :: mass_ratio

    mass_ratio = mass_gev/am3_mass_electron_gev
    hadronic_hadronic_ic_coeff = am3_c_cgs*am3_sigma_t_cgs/(mass_ratio*mass_ratio)
end function hadronic_hadronic_ic_coeff

! 获取能量网格的对数间距（调用验证逻辑）。
real(8) function hadronic_hadronic_ic_log_spacing(num_grid,energy_grid)
    integer, intent(in) :: num_grid
    real(8), intent(in) :: energy_grid(num_grid)
    real(8) :: dln_local
    call hadronic_validate_log_grid(num_grid,energy_grid,"hadronic_IC_grid",dln_local)
    hadronic_hadronic_ic_log_spacing = dln_local
end function hadronic_hadronic_ic_log_spacing

end module hadronic_hadronic_ic_kernel
