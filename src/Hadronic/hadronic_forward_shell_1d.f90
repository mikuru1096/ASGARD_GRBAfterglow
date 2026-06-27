module hadronic_forward_shell_1d
    implicit none
    private
    public :: hadronic_forward_pp_delta_shell
    public :: hadronic_forward_hadronic_ic_shell, hadronic_forward_hic_projected
    public :: hadronic_forward_species_transport_step, hadronic_forward_injection_content
    public :: hadronic_forward_global_gamma_p_max, hadronic_forward_secondary_radiation_shell, hadronic_forward_secondary_radiation_projected
    public :: hadronic_forward_continuous_loss_rates, hadronic_forward_secondary_electron_sequence
    public :: hadronic_forward_photon_loss_closure, hadronic_forward_interaction_effective_time, hadronic_forward_pgamma_proton_update
    public :: hadronic_forward_proton_transport_step, hadronic_forward_exponential_sink, hadronic_forward_energy_luminosity_from_rate
    public :: hadronic_forward_project_luminosity_from_rate, hadronic_forward_project_hic_luminosity, hadronic_forward_pair_source_content
    public :: hadronic_forward_shell_density_per_gev, hadronic_forward_gamma_edges, hadronic_forward_photon_density_hz_to_gev
    public :: hadronic_forward_process_power, hadronic_forward_positive_loglog_interp, hadronic_forward_source_per_gamma
    public :: hadronic_forward_distribution_per_gev, hadronic_forward_aligned_photon_grid, hadronic_forward_shell_volumes
    public :: hadronic_forward_shell_comoving_dt, hadronic_forward_dynamical_time
    public :: hadronic_forward_quantum_syn_cooling_factor
contains

! Single-shell pp delta-approximation wrapper.
subroutine hadronic_forward_pp_delta_shell(Num_p,proton_energy_gev,proton_density_per_gev,target_proton_density_cm3, &
                                      Num_gamma,gamma_energy_gev,Num_nu,neutrino_energy_gev,Num_pair,pair_energy_gev, &
                                      kappa_inelastic,pion_energy_fraction,neutral_pion_fraction,gamma_rate_per_gev, &
                                      neutrino_rate_per_gev,pair_rate_per_gev,proton_loss_rate)
    use hadronic_pp_kernel, only: hadronic_pp_delta_operator
    implicit none
    integer, intent(in) :: Num_p,Num_gamma,Num_nu,Num_pair
    real(8), intent(in) :: proton_energy_gev(Num_p),proton_density_per_gev(Num_p),target_proton_density_cm3
    real(8), intent(in) :: gamma_energy_gev(Num_gamma),neutrino_energy_gev(Num_nu),pair_energy_gev(Num_pair)
    real(8), intent(in) :: kappa_inelastic,pion_energy_fraction,neutral_pion_fraction
    real(8), intent(out) :: gamma_rate_per_gev(Num_gamma),neutrino_rate_per_gev(Num_nu),pair_rate_per_gev(Num_pair)
    real(8), intent(out) :: proton_loss_rate(Num_p)

    call hadronic_pp_delta_operator(num_p=Num_p,proton_energy_gev=proton_energy_gev, &
                                    proton_density_per_gev=proton_density_per_gev, &
                                    target_proton_density_cm3=target_proton_density_cm3, &
                                    num_gamma=Num_gamma,gamma_energy_gev=gamma_energy_gev, &
                                    num_nu=Num_nu,neutrino_energy_gev=neutrino_energy_gev, &
                                    num_pair=Num_pair,pair_energy_gev=pair_energy_gev, &
                                    gamma_rate_per_gev=gamma_rate_per_gev, &
                                    neutrino_rate_per_gev=neutrino_rate_per_gev, &
                                    pair_rate_per_gev=pair_rate_per_gev,proton_loss_rate=proton_loss_rate, &
                                    kappa_inelastic=kappa_inelastic,pion_energy_fraction=pion_energy_fraction, &
                                    neutral_pion_fraction=neutral_pion_fraction)
end subroutine hadronic_forward_pp_delta_shell

! Single-shell hadronic inverse-Compton wrapper.
subroutine hadronic_forward_hadronic_ic_shell(Num_had,hadron_energy_gev,Num_ph,photon_energy_gev,photons_on_had_grid_per_gev, &
                                         protons_per_gev,pion_plus_per_gev,pion_minus_per_gev,muon_minus_left_per_gev, &
                                         muon_minus_right_per_gev,muon_plus_left_per_gev,muon_plus_right_per_gev, &
                                         ind_min_energy_pho_hadgrid,epsilon_p_ic,epsilon_pi_ic,epsilon_mu_ic, &
                                         coeff_p_cgs,coeff_pi_cgs,coeff_mu_cgs,dln_energy,delta_e_p,jmax_p,delta_e_pi, &
                                         jmax_pi,delta_e_mu,jmax_mu)
    use hadronic_hadronic_ic_kernel, only: hadronic_hadronic_ic_initialize_kernel, &
                                           hadronic_hadronic_ic_operator_from_kernel
    implicit none
    integer, intent(in) :: Num_had,Num_ph,ind_min_energy_pho_hadgrid
    real(8), intent(in) :: hadron_energy_gev(Num_had),photon_energy_gev(Num_ph),photons_on_had_grid_per_gev(Num_ph)
    real(8), intent(in) :: protons_per_gev(Num_had),pion_plus_per_gev(Num_had),pion_minus_per_gev(Num_had)
    real(8), intent(in) :: muon_minus_left_per_gev(Num_had),muon_minus_right_per_gev(Num_had)
    real(8), intent(in) :: muon_plus_left_per_gev(Num_had),muon_plus_right_per_gev(Num_had)
    real(8), intent(out) :: epsilon_p_ic(Num_ph),epsilon_pi_ic(Num_ph),epsilon_mu_ic(Num_ph)
    real(8), intent(out) :: coeff_p_cgs,coeff_pi_cgs,coeff_mu_cgs,dln_energy
    integer, intent(out) :: delta_e_p(Num_had),jmax_p(Num_had),delta_e_pi(Num_had),jmax_pi(Num_had)
    integer, intent(out) :: delta_e_mu(Num_had),jmax_mu(Num_had)

    call hadronic_hadronic_ic_initialize_kernel(Num_had,hadron_energy_gev,Num_ph,photon_energy_gev, &
                                                ind_min_energy_pho_hadgrid,dln_energy, &
                                                delta_e_p,jmax_p,delta_e_pi,jmax_pi,delta_e_mu,jmax_mu)
    call hadronic_hadronic_ic_operator_from_kernel(Num_had,Num_ph,photons_on_had_grid_per_gev,protons_per_gev, &
                                                   pion_plus_per_gev,pion_minus_per_gev,muon_minus_left_per_gev, &
                                                   muon_minus_right_per_gev,muon_plus_left_per_gev, &
                                                   muon_plus_right_per_gev,dln_energy,delta_e_p,jmax_p, &
                                                   delta_e_pi,jmax_pi,delta_e_mu,jmax_mu,epsilon_p_ic, &
                                                   epsilon_pi_ic,epsilon_mu_ic,coeff_p_cgs,coeff_pi_cgs, &
                                                   coeff_mu_cgs)
end subroutine hadronic_forward_hadronic_ic_shell

! proton-only hadronic IC：目标光子插值、IC operator 和 luminosity 投影在 Fortran 内闭合。
subroutine hadronic_forward_hic_projected(num_had,num_ph,num_align,hadron_energy_gev,photon_energy_gev, &
                                     photon_density_per_gev,protons_per_gev,shell_volume_cm3,luminosity)
    implicit none
    integer, intent(in) :: num_had,num_ph,num_align
    real(8), intent(in) :: hadron_energy_gev(num_had),photon_energy_gev(num_ph),photon_density_per_gev(num_ph)
    real(8), intent(in) :: protons_per_gev(num_had),shell_volume_cm3
    real(8), intent(out) :: luminosity(num_ph)
    real(8) :: aligned_photon(num_align),photon_density_aligned(num_align),zero_had(num_had)
    real(8) :: epsilon_p(num_align),epsilon_pi(num_align),epsilon_mu(num_align)
    real(8) :: coeff_p,coeff_pi,coeff_mu,dln_energy
    integer :: delta_p(num_had),jmax_p(num_had),delta_pi(num_had),jmax_pi(num_had)
    integer :: delta_mu(num_had),jmax_mu(num_had)

    zero_had=0d0
    call hadronic_forward_aligned_photon_grid(num_had,num_ph,num_align,hadron_energy_gev,photon_energy_gev,aligned_photon)
    call hadronic_forward_positive_loglog_interp(num_ph,num_align,photon_energy_gev,photon_density_per_gev, &
                                            aligned_photon,photon_density_aligned)
    call hadronic_forward_hadronic_ic_shell(num_had,hadron_energy_gev,num_align,aligned_photon, &
                                       photon_density_aligned,protons_per_gev,zero_had,zero_had,zero_had,zero_had, &
                                       zero_had,zero_had,0,epsilon_p,epsilon_pi,epsilon_mu,coeff_p,coeff_pi, &
                                       coeff_mu,dln_energy,delta_p,jmax_p,delta_pi,jmax_pi,delta_mu,jmax_mu)
    call hadronic_forward_project_hic_luminosity(num_align,num_ph,aligned_photon,epsilon_p,epsilon_pi,epsilon_mu, &
                                            shell_volume_cm3,photon_energy_gev,luminosity)
end subroutine hadronic_forward_hic_projected

! 单壳层二级强子输运：源项映射、π/μ 推进和 neutron sink 在 Fortran 内闭合。
subroutine hadronic_forward_species_transport_step(num_gamma,num_src,gamma,source_energy_gev,neutron_source_per_gev_s, &
                                              pion_plus_source_per_gev_s,pion_minus_source_per_gev_s, &
                                              muon_minus_left_source_per_gev_s,muon_minus_right_source_per_gev_s, &
                                              muon_plus_left_source_per_gev_s,muon_plus_right_source_per_gev_s, &
                                              neutron_loss_src_s_inv,dt_s,b_field_g,divergence_rate_s_inv, &
                                              shell_volume_cm3,neutron_prev,pion_plus_prev,pion_minus_prev, &
                                              muon_minus_left_prev,muon_minus_right_prev,muon_plus_left_prev, &
                                              muon_plus_right_prev,neutron_next,pion_plus_next,pion_minus_next, &
                                              muon_minus_left_next,muon_minus_right_next,muon_plus_left_next, &
                                              muon_plus_right_next)
    use constants
    use hadronic_species_transport_kernel, only: hadronic_species_advance_operator
    implicit none
    integer, intent(in) :: num_gamma,num_src
    real(8), intent(in) :: gamma(num_gamma),source_energy_gev(num_src),dt_s,b_field_g,divergence_rate_s_inv
    real(8), intent(in) :: shell_volume_cm3,neutron_source_per_gev_s(num_src),pion_plus_source_per_gev_s(num_src)
    real(8), intent(in) :: pion_minus_source_per_gev_s(num_src),muon_minus_left_source_per_gev_s(num_src)
    real(8), intent(in) :: muon_minus_right_source_per_gev_s(num_src),muon_plus_left_source_per_gev_s(num_src)
    real(8), intent(in) :: muon_plus_right_source_per_gev_s(num_src),neutron_loss_src_s_inv(num_src)
    real(8), intent(in) :: neutron_prev(num_gamma),pion_plus_prev(num_gamma),pion_minus_prev(num_gamma)
    real(8), intent(in) :: muon_minus_left_prev(num_gamma),muon_minus_right_prev(num_gamma)
    real(8), intent(in) :: muon_plus_left_prev(num_gamma),muon_plus_right_prev(num_gamma)
    real(8), intent(out) :: neutron_next(num_gamma),pion_plus_next(num_gamma),pion_minus_next(num_gamma)
    real(8), intent(out) :: muon_minus_left_next(num_gamma),muon_minus_right_next(num_gamma)
    real(8), intent(out) :: muon_plus_left_next(num_gamma),muon_plus_right_next(num_gamma)
    integer :: i
    real(8) :: neutron_energy(num_gamma),pion_energy(num_gamma),muon_energy(num_gamma)
    real(8) :: neutron_source(num_gamma),pion_plus_source(num_gamma),pion_minus_source(num_gamma)
    real(8) :: muon_ml_source(num_gamma),muon_mr_source(num_gamma),muon_pl_source(num_gamma),muon_pr_source(num_gamma)
    real(8) :: neutron_loss(num_gamma),neutron_transport(num_gamma)

    do i=1,num_gamma
        neutron_energy(i)=gamma(i)*Para_m_n_GeV
        pion_energy(i)=gamma(i)*Para_m_pi_charged_GeV
        muon_energy(i)=gamma(i)*Para_m_mu_GeV
    end do
    call hadronic_forward_source_per_gamma(num_src,num_gamma,source_energy_gev,neutron_source_per_gev_s, &
                                      neutron_energy,Para_m_n_GeV,shell_volume_cm3,neutron_source)
    call hadronic_forward_source_per_gamma(num_src,num_gamma,source_energy_gev,pion_plus_source_per_gev_s, &
                                      pion_energy,Para_m_pi_charged_GeV,shell_volume_cm3,pion_plus_source)
    call hadronic_forward_source_per_gamma(num_src,num_gamma,source_energy_gev,pion_minus_source_per_gev_s, &
                                      pion_energy,Para_m_pi_charged_GeV,shell_volume_cm3,pion_minus_source)
    call hadronic_forward_source_per_gamma(num_src,num_gamma,source_energy_gev,muon_minus_left_source_per_gev_s, &
                                      muon_energy,Para_m_mu_GeV,shell_volume_cm3,muon_ml_source)
    call hadronic_forward_source_per_gamma(num_src,num_gamma,source_energy_gev,muon_minus_right_source_per_gev_s, &
                                      muon_energy,Para_m_mu_GeV,shell_volume_cm3,muon_mr_source)
    call hadronic_forward_source_per_gamma(num_src,num_gamma,source_energy_gev,muon_plus_left_source_per_gev_s, &
                                      muon_energy,Para_m_mu_GeV,shell_volume_cm3,muon_pl_source)
    call hadronic_forward_source_per_gamma(num_src,num_gamma,source_energy_gev,muon_plus_right_source_per_gev_s, &
                                      muon_energy,Para_m_mu_GeV,shell_volume_cm3,muon_pr_source)
    call hadronic_species_advance_operator(num_gamma,gamma,dt_s,b_field_g,divergence_rate_s_inv,neutron_prev, &
                                           pion_plus_prev,pion_minus_prev,muon_minus_left_prev, &
                                           muon_minus_right_prev,muon_plus_left_prev,muon_plus_right_prev, &
                                           neutron_source,pion_plus_source,pion_minus_source,muon_ml_source, &
                                           muon_mr_source,muon_pl_source,muon_pr_source,neutron_transport, &
                                           pion_plus_next,pion_minus_next,muon_minus_left_next, &
                                           muon_minus_right_next,muon_plus_left_next,muon_plus_right_next)
    call hadronic_forward_positive_loglog_interp(num_src,num_gamma,source_energy_gev,neutron_loss_src_s_inv, &
                                            neutron_energy,neutron_loss)
    call hadronic_forward_exponential_sink(num_gamma,neutron_transport,neutron_loss,dt_s,neutron_next)
end subroutine hadronic_forward_species_transport_step

! 单壳层注入 content：按 shell energy/dt 得到源项，再乘回本步 dt。
subroutine hadronic_forward_injection_content(num_gamma,species,gamma,shell_energy_erg,dt_s,spectral_index, &
                                         gamma_min,gamma_max_inj,gamma_cut,has_gamma_cut,q_content)
    use hadronic_acceleration_kernel, only: hadronic_species_injection_operator
    implicit none
    integer, intent(in) :: num_gamma
    character(len=*), intent(in) :: species
    real(8), intent(in) :: gamma(num_gamma),shell_energy_erg,dt_s,spectral_index
    real(8), intent(in) :: gamma_min,gamma_max_inj,gamma_cut
    logical, intent(in) :: has_gamma_cut
    real(8), intent(out) :: q_content(num_gamma)

    if (dt_s <= 0d0) error stop "hadronic injection content requires dt_s > 0."
    if (shell_energy_erg <= 0d0) then
        q_content=0d0
        return
    end if
    call hadronic_species_injection_operator(num_gamma,gamma,species,shell_energy_erg/dt_s,spectral_index, &
                                             gamma_min,gamma_max_inj,gamma_cut,has_gamma_cut,q_content)
    q_content(1:num_gamma)=dt_s*q_content(1:num_gamma)
end subroutine hadronic_forward_injection_content

! 全局质子最大 Lorentz 因子：逐半径壳层估计并取最大值作为输运网格上界。
subroutine hadronic_forward_global_gamma_p_max(num_r,radius_cm,gamma_bulk,b_field_g,eta_acc,gamma_max_global)
    use hadronic_acceleration_kernel, only: hadronic_estimate_max_gamma
    implicit none
    integer, intent(in) :: num_r
    real(8), intent(in) :: radius_cm(num_r),gamma_bulk(num_r),b_field_g(num_r),eta_acc
    real(8), intent(out) :: gamma_max_global
    integer :: i
    real(8) :: gamma_scan(2),external_rate(2),gamma_max,gamma_dyn,gamma_syn,gamma_ext
    logical :: has_gamma_ext

    if (eta_acc <= 0d0) error stop "hadronic eta_acc must be positive."
    gamma_scan=(/1d0,2d0/)
    external_rate=(/1d0,2d0/)
    gamma_max_global=0d0
    do i=1,num_r
        call hadronic_estimate_max_gamma("proton",b_field_g(i),radius_cm(i),gamma_bulk(i),eta_acc, &
                                         2,gamma_scan,external_rate,.false.,gamma_max,gamma_dyn, &
                                         gamma_syn,gamma_ext,has_gamma_ext)
        gamma_max_global=max(gamma_max_global,gamma_max)
    end do
    if (gamma_max_global <= 1d0) error stop "hadronic maximum proton Lorentz factor must exceed unity."
end subroutine hadronic_forward_global_gamma_p_max

! Single-shell secondary-radiation wrapper.
subroutine hadronic_forward_secondary_radiation_shell(Num_had,hadron_energy_gev,Num_ph,photon_energy_gev,pion_plus_per_gev, &
                                                 pion_minus_per_gev,muon_minus_left_per_gev,muon_minus_right_per_gev, &
                                                 muon_plus_left_per_gev,muon_plus_right_per_gev,photons_on_had_grid_per_gev, &
                                                 ind_min_energy_pho_hadgrid,magnetic_field_g,pion_synch_rate_per_gev, &
                                                 muon_synch_rate_per_gev,pion_ic_rate_per_gev,muon_ic_rate_per_gev)
    use hadronic_secondary_radiation_kernel, only: hadronic_secondary_radiation_operator
    implicit none
    integer, intent(in) :: Num_had,Num_ph,ind_min_energy_pho_hadgrid
    real(8), intent(in) :: hadron_energy_gev(Num_had),photon_energy_gev(Num_ph)
    real(8), intent(in) :: pion_plus_per_gev(Num_had),pion_minus_per_gev(Num_had)
    real(8), intent(in) :: muon_minus_left_per_gev(Num_had),muon_minus_right_per_gev(Num_had)
    real(8), intent(in) :: muon_plus_left_per_gev(Num_had),muon_plus_right_per_gev(Num_had)
    real(8), intent(in) :: photons_on_had_grid_per_gev(Num_ph),magnetic_field_g
    real(8), intent(out) :: pion_synch_rate_per_gev(Num_ph),muon_synch_rate_per_gev(Num_ph)
    real(8), intent(out) :: pion_ic_rate_per_gev(Num_ph),muon_ic_rate_per_gev(Num_ph)
    real(8) :: synch_dln_energy,ic_dln_energy
    real(8) :: kernel_pion(Num_ph,Num_had),kernel_muon(Num_ph,Num_had)
    integer :: delta_e_pi(Num_had),jmax_pi(Num_had),delta_e_mu(Num_had),jmax_mu(Num_had)

    call hadronic_secondary_radiation_operator(Num_had,hadron_energy_gev,Num_ph,photon_energy_gev, &
                                               pion_plus_per_gev,pion_minus_per_gev,muon_minus_left_per_gev, &
                                               muon_minus_right_per_gev,muon_plus_left_per_gev,muon_plus_right_per_gev, &
                                               photons_on_had_grid_per_gev,ind_min_energy_pho_hadgrid,magnetic_field_g, &
                                               pion_synch_rate_per_gev,muon_synch_rate_per_gev,pion_ic_rate_per_gev, &
                                               muon_ic_rate_per_gev,synch_dln_energy,kernel_pion,kernel_muon, &
                                               ic_dln_energy,delta_e_pi,jmax_pi,delta_e_mu,jmax_mu)
end subroutine hadronic_forward_secondary_radiation_shell

! 投影后的二级 π/μ 辐射：分布映射、目标光子插值、辐射和 luminosity 投影在 Fortran 内闭合。
subroutine hadronic_forward_secondary_radiation_projected(num_had,num_ph,num_align,hadron_energy_gev, &
                                                     photon_energy_gev,photon_density_per_gev, &
                                                     pion_plus_density_per_gamma,pion_minus_density_per_gamma, &
                                                     muon_minus_left_density_per_gamma, &
                                                     muon_minus_right_density_per_gamma, &
                                                     muon_plus_left_density_per_gamma, &
                                                     muon_plus_right_density_per_gamma, &
                                                     shell_volume_cm3,magnetic_field_g,pion_synch_luminosity, &
                                                     muon_synch_luminosity,pion_ic_luminosity,muon_ic_luminosity)
    use constants
    implicit none
    integer, intent(in) :: num_had,num_ph,num_align
    real(8), intent(in) :: hadron_energy_gev(num_had),photon_energy_gev(num_ph),photon_density_per_gev(num_ph)
    real(8), intent(in) :: pion_plus_density_per_gamma(num_had),pion_minus_density_per_gamma(num_had)
    real(8), intent(in) :: muon_minus_left_density_per_gamma(num_had),muon_minus_right_density_per_gamma(num_had)
    real(8), intent(in) :: muon_plus_left_density_per_gamma(num_had),muon_plus_right_density_per_gamma(num_had)
    real(8), intent(in) :: shell_volume_cm3,magnetic_field_g
    real(8), intent(out) :: pion_synch_luminosity(num_ph),muon_synch_luminosity(num_ph)
    real(8), intent(out) :: pion_ic_luminosity(num_ph),muon_ic_luminosity(num_ph)
    integer :: i
    real(8) :: gamma_species(num_had),pion_energy(num_had),muon_energy(num_had),aligned_photon(num_align)
    real(8) :: photon_density_aligned(num_align),pion_plus_per_gev(num_had),pion_minus_per_gev(num_had)
    real(8) :: muon_ml_per_gev(num_had),muon_mr_per_gev(num_had),muon_pl_per_gev(num_had),muon_pr_per_gev(num_had)
    real(8) :: pion_synch_rate(num_align),muon_synch_rate(num_align),pion_ic_rate(num_align),muon_ic_rate(num_align)

    do i=1,num_had
        gamma_species(i)=hadron_energy_gev(i)/Para_m_p_GeV
        pion_energy(i)=gamma_species(i)*Para_m_pi_charged_GeV
        muon_energy(i)=gamma_species(i)*Para_m_mu_GeV
    end do
    call hadronic_forward_distribution_per_gev(num_had,num_had,pion_energy,pion_plus_density_per_gamma, &
                                          hadron_energy_gev,Para_m_pi_charged_GeV,shell_volume_cm3,pion_plus_per_gev)
    call hadronic_forward_distribution_per_gev(num_had,num_had,pion_energy,pion_minus_density_per_gamma, &
                                          hadron_energy_gev,Para_m_pi_charged_GeV,shell_volume_cm3,pion_minus_per_gev)
    call hadronic_forward_distribution_per_gev(num_had,num_had,muon_energy,muon_minus_left_density_per_gamma, &
                                          hadron_energy_gev,Para_m_mu_GeV,shell_volume_cm3,muon_ml_per_gev)
    call hadronic_forward_distribution_per_gev(num_had,num_had,muon_energy,muon_minus_right_density_per_gamma, &
                                          hadron_energy_gev,Para_m_mu_GeV,shell_volume_cm3,muon_mr_per_gev)
    call hadronic_forward_distribution_per_gev(num_had,num_had,muon_energy,muon_plus_left_density_per_gamma, &
                                          hadron_energy_gev,Para_m_mu_GeV,shell_volume_cm3,muon_pl_per_gev)
    call hadronic_forward_distribution_per_gev(num_had,num_had,muon_energy,muon_plus_right_density_per_gamma, &
                                          hadron_energy_gev,Para_m_mu_GeV,shell_volume_cm3,muon_pr_per_gev)
    call hadronic_forward_aligned_photon_grid(num_had,num_ph,num_align,hadron_energy_gev,photon_energy_gev,aligned_photon)
    call hadronic_forward_positive_loglog_interp(num_ph,num_align,photon_energy_gev,photon_density_per_gev, &
                                            aligned_photon,photon_density_aligned)
    call hadronic_forward_secondary_radiation_shell(num_had,hadron_energy_gev,num_align,aligned_photon,pion_plus_per_gev, &
                                               pion_minus_per_gev,muon_ml_per_gev,muon_mr_per_gev,muon_pl_per_gev, &
                                               muon_pr_per_gev,photon_density_aligned,0,magnetic_field_g, &
                                               pion_synch_rate,muon_synch_rate,pion_ic_rate,muon_ic_rate)
    call hadronic_forward_project_luminosity_from_rate(num_align,num_ph,aligned_photon,pion_synch_rate, &
                                                  shell_volume_cm3,photon_energy_gev,pion_synch_luminosity)
    call hadronic_forward_project_luminosity_from_rate(num_align,num_ph,aligned_photon,muon_synch_rate, &
                                                  shell_volume_cm3,photon_energy_gev,muon_synch_luminosity)
    call hadronic_forward_project_luminosity_from_rate(num_align,num_ph,aligned_photon,pion_ic_rate, &
                                                  shell_volume_cm3,photon_energy_gev,pion_ic_luminosity)
    call hadronic_forward_project_luminosity_from_rate(num_align,num_ph,aligned_photon,muon_ic_rate, &
                                                  shell_volume_cm3,photon_energy_gev,muon_ic_luminosity)
end subroutine hadronic_forward_secondary_radiation_projected

! 连续冷却率 wrapper：绝热项加同步项，可选 quantum synch 修正。
subroutine hadronic_forward_continuous_loss_rates(num_gamma,gamma,b_field_g,t_dyn_s,mass_gev,quantum_syn,loss_total)
    use constants
    use hadronic_common, only: hadronic_quantum_syn_cooling_factor
    implicit none
    integer, intent(in) :: num_gamma,quantum_syn
    real(8), intent(in) :: gamma(num_gamma),b_field_g,t_dyn_s,mass_gev
    real(8), intent(out) :: loss_total(num_gamma)
    integer :: i
    real(8) :: coeff_syn,mass_ratio,syn_loss

    if (b_field_g < zero) error stop "hadronic continuous loss rates require b_field_g >= 0."
    if (t_dyn_s <= zero) error stop "hadronic continuous loss rates require t_dyn_s > 0."
    mass_ratio=mass_gev/Para_m_e_GeV
    coeff_syn=Para_sigmaT*b_field_g*b_field_g/(6d0*pi*Para_m_e*Para_c)/(mass_ratio**3)
    do i=1,num_gamma
        syn_loss=coeff_syn*gamma(i)*gamma(i)
        if (quantum_syn /= 0) syn_loss=syn_loss*hadronic_quantum_syn_cooling_factor(gamma(i),b_field_g,mass_gev)
        loss_total(i)=gamma(i)/t_dyn_s+syn_loss
    end do
end subroutine hadronic_forward_continuous_loss_rates

! BH/pp 二级电子序列：按半径壳层推进冷却谱，并输出同步辐射源项。
subroutine hadronic_forward_secondary_electron_sequence(num_e,num_nu,num_r,gamma_e,radius_cm,gamma_bulk,b_field_g, &
                                                   frequency_hz,source_content,index_syn_integr,n_threads,quantum_syn, &
                                                   electron_density,luminosity_syn,seed_syn,source_radius)
    use constants
    use electron_radiation_kernel, only: get_syn_selected
    use hadronic_transport_remap_kernel, only: hadronic_advance_energy_loggamma_remap
    implicit none
    integer, intent(in) :: num_e,num_nu,num_r,index_syn_integr,n_threads,quantum_syn
    real(8), intent(in) :: gamma_e(num_e),radius_cm(num_r),gamma_bulk(num_r),b_field_g(num_r)
    real(8), intent(in) :: frequency_hz(num_nu),source_content(num_e,num_r)
    real(8), intent(out) :: electron_density(num_e,num_r),luminosity_syn(num_nu,num_r),seed_syn(num_nu,num_r)
    real(8), intent(out) :: source_radius(num_e,num_r)
    integer :: i_r
    real(8) :: dr,dt_s,t_dyn_s,loss_total(num_e),prev_density(num_e),next_density(num_e)

    electron_density=zero; luminosity_syn=zero; seed_syn=zero; source_radius=zero
    prev_density=zero
    do i_r=1,num_r
        call hadronic_sequence_shell_geometry(num_r,radius_cm,gamma_bulk,i_r,dr,dt_s)
        t_dyn_s=radius_cm(i_r)/(gamma_bulk(i_r)*Para_c)
        call hadronic_forward_continuous_loss_rates(num_e,gamma_e,b_field_g(i_r),t_dyn_s, &
                                               Para_m_e_GeV,quantum_syn,loss_total)
        call hadronic_advance_energy_loggamma_remap(num_e,gamma_e,prev_density, &
                                                    source_content(:,i_r),loss_total,dt_s,next_density)
        call get_syn_selected(index_syn_integr,radius_cm(i_r),b_field_g(i_r),num_e,num_nu,n_threads, &
                              gamma_e,next_density,frequency_hz,luminosity_syn(:,i_r),seed_syn(:,i_r))
        source_radius(:,i_r)=source_content(:,i_r)*dt_s/dr*gamma_e(:)*dlog(ten)
        electron_density(:,i_r)=next_density
        prev_density=next_density
    end do
end subroutine hadronic_forward_secondary_electron_sequence

! 光子损失闭合：由壳层 comoving path time 得到 tau，并返回局部 survival 平均因子。
subroutine hadronic_forward_photon_loss_closure(num_ph,num_r,radius_cm,gamma_bulk,shell_index,loss_rate,tau, &
                                           survival)
    use constants
    implicit none
    integer, intent(in) :: num_ph,num_r,shell_index
    real(8), intent(in) :: radius_cm(num_r),gamma_bulk(num_r),loss_rate(num_ph)
    real(8), intent(out) :: tau(num_ph),survival(num_ph)
    integer :: i_ph
    real(8) :: dr,dt_s

    if (any(loss_rate < zero)) error stop "hadronic photon loss closure requires non-negative loss rate."
    call hadronic_sequence_shell_geometry(num_r,radius_cm,gamma_bulk,shell_index,dr,dt_s)
    do i_ph=1,num_ph
        tau(i_ph)=loss_rate(i_ph)*dt_s
        if (tau(i_ph) > 1d-6) then
            survival(i_ph)=(one-dexp(-tau(i_ph)))/tau(i_ph)
        else if (tau(i_ph) > zero) then
            survival(i_ph)=one-tau(i_ph)/two+tau(i_ph)*tau(i_ph)/6d0
        else
            survival(i_ph)=one
        end if
    end do
end subroutine hadronic_forward_photon_loss_closure

! 相互作用有效时间：对指数 sink 的同一步 reinjection 积分。
subroutine hadronic_forward_interaction_effective_time(num_rate,rate_s_inv,dt_s,effective_time_s)
    use constants
    implicit none
    integer, intent(in) :: num_rate
    real(8), intent(in) :: rate_s_inv(num_rate),dt_s
    real(8), intent(out) :: effective_time_s(num_rate)
    integer :: i
    real(8) :: tau

    if (dt_s <= zero) error stop "hadronic interaction effective time requires dt_s > 0."
    if (any(rate_s_inv < zero)) error stop "hadronic interaction effective time requires non-negative rates."
    do i=1,num_rate
        if (rate_s_inv(i) > zero) then
            tau=rate_s_inv(i)*dt_s
            if (tau < 1d-4) then
                effective_time_s(i)=dt_s*(one-tau/two+tau*tau/6d0)
            else
                effective_time_s(i)=(one-dexp(-tau))/rate_s_inv(i)
            end if
        else
            effective_time_s(i)=dt_s
        end if
    end do
end subroutine hadronic_forward_interaction_effective_time

! pγ 后质子更新：指数 sink 与同一步 reinjection 在同一壳层内闭合。
subroutine hadronic_forward_pgamma_proton_update(num_gamma,dn_transport,loss_rate_s_inv,reinj_rate_per_gev, &
                                            shell_volume_cm3,dt_s,dn_next)
    use constants
    implicit none
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: dn_transport(num_gamma),loss_rate_s_inv(num_gamma),reinj_rate_per_gev(num_gamma)
    real(8), intent(in) :: shell_volume_cm3,dt_s
    real(8), intent(out) :: dn_next(num_gamma)
    integer :: i
    real(8) :: effective_time(num_gamma)

    if (shell_volume_cm3 <= zero) error stop "hadronic p-gamma proton update requires positive shell volume."
    call hadronic_forward_interaction_effective_time(num_gamma,loss_rate_s_inv,dt_s,effective_time)
    do i=1,num_gamma
        dn_next(i)=dn_transport(i)*dexp(-loss_rate_s_inv(i)*dt_s)+ &
                   effective_time(i)*shell_volume_cm3*reinj_rate_per_gev(i)*Para_m_p_GeV
    end do
end subroutine hadronic_forward_pgamma_proton_update

! 单壳层质子输运：连续冷却 + BH/pp 损失 + pγ sink/reinjection 在 Fortran 内闭合。
subroutine hadronic_forward_proton_transport_step(num_gamma,gamma,dn_prev,q_inj,b_field_g,t_dyn_s,mass_gev, &
                                             quantum_syn,bh_loss,pp_loss,pg_loss_rate,pg_reinj_rate_per_gev, &
                                             shell_volume_cm3,dt_s,dn_next)
    use hadronic_transport_remap_kernel, only: hadronic_advance_energy_loggamma_remap
    implicit none
    integer, intent(in) :: num_gamma,quantum_syn
    real(8), intent(in) :: gamma(num_gamma),dn_prev(num_gamma),q_inj(num_gamma)
    real(8), intent(in) :: b_field_g,t_dyn_s,mass_gev,shell_volume_cm3,dt_s
    real(8), intent(in) :: bh_loss(num_gamma),pp_loss(num_gamma),pg_loss_rate(num_gamma)
    real(8), intent(in) :: pg_reinj_rate_per_gev(num_gamma)
    real(8), intent(out) :: dn_next(num_gamma)
    real(8) :: loss_total(num_gamma),dn_transport(num_gamma)

    call hadronic_forward_continuous_loss_rates(num_gamma,gamma,b_field_g,t_dyn_s,mass_gev,quantum_syn,loss_total)
    loss_total(1:num_gamma)=loss_total(1:num_gamma)+bh_loss(1:num_gamma)+pp_loss(1:num_gamma)
    call hadronic_advance_energy_loggamma_remap(num_gamma,gamma,dn_prev,q_inj,loss_total,dt_s,dn_transport)
    call hadronic_forward_pgamma_proton_update(num_gamma,dn_transport,pg_loss_rate,pg_reinj_rate_per_gev, &
                                          shell_volume_cm3,dt_s,dn_next)
end subroutine hadronic_forward_proton_transport_step

! 指数 sink：用于单步粒子损失项 N -> N exp(-rate dt)。
subroutine hadronic_forward_exponential_sink(num_value,values,loss_rate_s_inv,dt_s,values_next)
    use constants
    implicit none
    integer, intent(in) :: num_value
    real(8), intent(in) :: values(num_value),loss_rate_s_inv(num_value),dt_s
    real(8), intent(out) :: values_next(num_value)
    integer :: i

    if (dt_s <= zero) error stop "hadronic exponential sink requires dt_s > 0."
    if (any(loss_rate_s_inv < zero)) error stop "hadronic exponential sink requires non-negative rates."
    do i=1,num_value
        values_next(i)=values(i)*dexp(-loss_rate_s_inv(i)*dt_s)
    end do
end subroutine hadronic_forward_exponential_sink

! 将每能量反应率谱转换为 shell luminosity 谱。
subroutine hadronic_forward_energy_luminosity_from_rate(num_energy,energy_gev,rate_per_gev,shell_volume_cm3, &
                                                   luminosity)
    use constants
    implicit none
    integer, intent(in) :: num_energy
    real(8), intent(in) :: energy_gev(num_energy),rate_per_gev(num_energy),shell_volume_cm3
    real(8), intent(out) :: luminosity(num_energy)

    if (shell_volume_cm3 <= zero) error stop "hadronic luminosity conversion requires positive shell volume."
    luminosity(1:num_energy)=shell_volume_cm3*rate_per_gev(1:num_energy)*energy_gev(1:num_energy)*Para_h_GeV*Para_GeV2erg
end subroutine hadronic_forward_energy_luminosity_from_rate

! rate 谱投影：先转壳层 luminosity，再映射到目标 photon energy grid。
subroutine hadronic_forward_project_luminosity_from_rate(num_src,num_dst,energy_src_gev,rate_src_per_gev, &
                                                    shell_volume_cm3,energy_dst_gev,luminosity_dst)
    use constants
    implicit none
    integer, intent(in) :: num_src,num_dst
    real(8), intent(in) :: energy_src_gev(num_src),rate_src_per_gev(num_src),energy_dst_gev(num_dst)
    real(8), intent(in) :: shell_volume_cm3
    real(8), intent(out) :: luminosity_dst(num_dst)
    real(8) :: luminosity_src(num_src)

    if (shell_volume_cm3 <= 0d0) error stop "hadronic luminosity projection requires positive shell volume."
    luminosity_src(1:num_src)=shell_volume_cm3*rate_src_per_gev(1:num_src)*energy_src_gev(1:num_src)* &
                              Para_h_GeV*Para_GeV2erg
    call hadronic_forward_positive_loglog_interp(num_src,num_dst,energy_src_gev,luminosity_src, &
                                            energy_dst_gev,luminosity_dst)
end subroutine hadronic_forward_project_luminosity_from_rate

! hadronic IC luminosity 投影：合并 p/pi/mu IC 能量源项并映射到目标 photon grid。
subroutine hadronic_forward_project_hic_luminosity(num_src,num_dst,energy_src_gev,epsilon_p_ic,epsilon_pi_ic, &
                                              epsilon_mu_ic,shell_volume_cm3,energy_dst_gev,luminosity_dst)
    use constants
    implicit none
    integer, intent(in) :: num_src,num_dst
    real(8), intent(in) :: energy_src_gev(num_src),epsilon_p_ic(num_src),epsilon_pi_ic(num_src),epsilon_mu_ic(num_src)
    real(8), intent(in) :: shell_volume_cm3,energy_dst_gev(num_dst)
    real(8), intent(out) :: luminosity_dst(num_dst)
    real(8) :: luminosity_src(num_src)

    if (shell_volume_cm3 <= 0d0) error stop "hadronic IC luminosity projection requires positive shell volume."
    luminosity_src(1:num_src)=shell_volume_cm3*(epsilon_p_ic(1:num_src)+epsilon_pi_ic(1:num_src)+ &
                              epsilon_mu_ic(1:num_src))*Para_h_GeV*Para_GeV2erg
    call hadronic_forward_positive_loglog_interp(num_src,num_dst,energy_src_gev,luminosity_src, &
                                            energy_dst_gev,luminosity_dst)
end subroutine hadronic_forward_project_hic_luminosity

! BH/pp 二级电子源项：把每 GeV 产生率合并为壳层内每 gamma 注入 content。
subroutine hadronic_forward_pair_source_content(num_e,pp_pair_rate_per_gev,bh_pair_rate_per_gev,include_pp, &
                                           include_bh,shell_volume_cm3,source_content)
    use constants
    implicit none
    integer, intent(in) :: num_e,include_pp,include_bh
    real(8), intent(in) :: pp_pair_rate_per_gev(num_e),bh_pair_rate_per_gev(num_e),shell_volume_cm3
    real(8), intent(out) :: source_content(num_e)

    if (shell_volume_cm3 <= 0d0) error stop "hadronic pair source content requires positive shell volume."
    source_content=0d0
    if (include_pp /= 0) source_content(1:num_e)=source_content(1:num_e)+pp_pair_rate_per_gev(1:num_e)
    if (include_bh /= 0) source_content(1:num_e)=source_content(1:num_e)+bh_pair_rate_per_gev(1:num_e)
    source_content(1:num_e)=shell_volume_cm3*Para_m_e_GeV*source_content(1:num_e)
end subroutine hadronic_forward_pair_source_content

! 壳层粒子密度归一化：从每 gamma 数量转换为每 GeV 体密度。
subroutine hadronic_forward_shell_density_per_gev(num_gamma,density_per_gamma,mass_gev,shell_volume_cm3, &
                                             density_per_gev)
    implicit none
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: density_per_gamma(num_gamma),mass_gev,shell_volume_cm3
    real(8), intent(out) :: density_per_gev(num_gamma)

    if (mass_gev <= 0d0) error stop "hadronic shell density per GeV requires positive particle mass."
    if (shell_volume_cm3 <= 0d0) error stop "hadronic shell density per GeV requires positive shell volume."
    density_per_gev(1:num_gamma)=density_per_gamma(1:num_gamma)/(shell_volume_cm3*mass_gev)
end subroutine hadronic_forward_shell_density_per_gev

! gamma 网格边界：对数中心网格的几何中点边界。
subroutine hadronic_forward_gamma_edges(num_gamma,gamma,gamma_edge)
    implicit none
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: gamma(num_gamma)
    real(8), intent(out) :: gamma_edge(num_gamma+1)
    integer :: i

    if (num_gamma < 1) error stop "hadronic gamma edges require at least one grid point."
    if (num_gamma == 1) then
        if (gamma(1) <= 1d0) error stop "single-point hadronic gamma grid must exceed unity."
        gamma_edge(1)=0.5d0*gamma(1)
        gamma_edge(2)=2d0*gamma(1)
        return
    end if
    gamma_edge(1)=gamma(1)*dsqrt(gamma(1)/gamma(2))
    do i=2,num_gamma
        gamma_edge(i)=dsqrt(gamma(i-1)*gamma(i))
    end do
    gamma_edge(num_gamma+1)=gamma(num_gamma)*dsqrt(gamma(num_gamma)/gamma(num_gamma-1))
end subroutine hadronic_forward_gamma_edges

! photon density 单位变换：E=h nu，n_E=n_nu/h。
subroutine hadronic_forward_photon_density_hz_to_gev(num_ph,photon_nu_hz,photon_density_per_hz, &
                                                photon_energy_gev,photon_density_per_gev)
    use constants
    implicit none
    integer, intent(in) :: num_ph
    real(8), intent(in) :: photon_nu_hz(num_ph),photon_density_per_hz(num_ph)
    real(8), intent(out) :: photon_energy_gev(num_ph),photon_density_per_gev(num_ph)

    photon_energy_gev(1:num_ph)=Para_h_GeV*photon_nu_hz(1:num_ph)
    photon_density_per_gev(1:num_ph)=photon_density_per_hz(1:num_ph)/Para_h_GeV
end subroutine hadronic_forward_photon_density_hz_to_gev

! AM3 分过程功率归并：积分每个过程 luminosity，并按质子能量分布投到 hadron grid。
subroutine hadronic_forward_process_power(num_had,num_proc_energy,num_process,hadron_energy_gev,dn_had, &
                                     process_energy_gev,process_rate_per_gev,shell_volume_cm3,process_power)
    use constants
    implicit none
    integer, intent(in) :: num_had,num_proc_energy,num_process
    real(8), intent(in) :: hadron_energy_gev(num_had),dn_had(num_had),process_energy_gev(num_proc_energy)
    real(8), intent(in) :: process_rate_per_gev(num_process,num_proc_energy),shell_volume_cm3
    real(8), intent(out) :: process_power(num_process,num_had)
    integer :: i,j
    real(8) :: proton_weight(num_had),total_weight,luminosity(num_proc_energy),proc_total

    if (shell_volume_cm3 <= 0d0) error stop "hadronic process power requires positive shell volume."
    do i=1,num_had
        proton_weight(i)=dn_had(i)*hadron_energy_gev(i)
    end do
    total_weight=hadronic_trapezoid(num_had,hadron_energy_gev,proton_weight)
    process_power=0d0
    if (total_weight <= 0d0) return
    do j=1,num_process
        luminosity(1:num_proc_energy)=shell_volume_cm3*process_rate_per_gev(j,1:num_proc_energy)* &
                                      process_energy_gev(1:num_proc_energy)*Para_h_GeV*Para_GeV2erg
        proc_total=hadronic_trapezoid(num_proc_energy,process_energy_gev,luminosity)
        process_power(j,1:num_had)=proton_weight(1:num_had)/total_weight*proc_total
    end do
contains
    real(8) function hadronic_trapezoid(num_x,x,y)
        implicit none
        integer, intent(in) :: num_x
        real(8), intent(in) :: x(num_x),y(num_x)
        integer :: k

        hadronic_trapezoid=0d0
        do k=1,num_x-1
            hadronic_trapezoid=hadronic_trapezoid+0.5d0*(y(k)+y(k+1))*(x(k+1)-x(k))
        end do
    end function hadronic_trapezoid
end subroutine hadronic_forward_process_power

! 正值 log-log 插值：只使用 finite 且正的源点，范围外输出零。
subroutine hadronic_forward_positive_loglog_interp(num_src,num_dst,x_src,y_src,x_dst,y_dst)
    use ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer, intent(in) :: num_src,num_dst
    real(8), intent(in) :: x_src(num_src),y_src(num_src),x_dst(num_dst)
    real(8), intent(out) :: y_dst(num_dst)
    integer :: i,j,n_valid
    real(8) :: xv(num_src),yv(num_src),frac

    y_dst=0d0; n_valid=0
    do i=1,num_src
        if (ieee_is_finite(x_src(i)) .and. ieee_is_finite(y_src(i)) .and. x_src(i) > 0d0 .and. y_src(i) > 0d0) then
            n_valid=n_valid+1
            xv(n_valid)=x_src(i)
            yv(n_valid)=y_src(i)
        end if
    end do
    if (n_valid < 2) return
    do i=1,num_dst
        if (x_dst(i) < xv(1) .or. x_dst(i) > xv(n_valid) .or. .not. ieee_is_finite(x_dst(i))) cycle
        do j=1,n_valid-1
            if (x_dst(i) >= xv(j) .and. x_dst(i) <= xv(j+1)) then
                frac=dlog(x_dst(i)/xv(j))/dlog(xv(j+1)/xv(j))
                y_dst(i)=dexp(dlog(yv(j))+frac*dlog(yv(j+1)/yv(j)))
                exit
            end if
        end do
    end do
end subroutine hadronic_forward_positive_loglog_interp

! pγ 源项映射：从每 GeV 反应率谱变换到二级粒子的每 gamma 注入率。
subroutine hadronic_forward_source_per_gamma(num_src,num_dst,energy_src_gev,source_src_per_gev_s,energy_dst_gev, &
                                        mass_gev,shell_volume_cm3,source_dst_per_gamma_s)
    implicit none
    integer, intent(in) :: num_src,num_dst
    real(8), intent(in) :: energy_src_gev(num_src),source_src_per_gev_s(num_src),energy_dst_gev(num_dst)
    real(8), intent(in) :: mass_gev,shell_volume_cm3
    real(8), intent(out) :: source_dst_per_gamma_s(num_dst)

    if (mass_gev <= 0d0) error stop "hadronic source per gamma requires positive particle mass."
    if (shell_volume_cm3 <= 0d0) error stop "hadronic source per gamma requires positive shell volume."
    call hadronic_forward_positive_loglog_interp(num_src,num_dst,energy_src_gev,source_src_per_gev_s, &
                                            energy_dst_gev,source_dst_per_gamma_s)
    source_dst_per_gamma_s(1:num_dst)=shell_volume_cm3*mass_gev*source_dst_per_gamma_s(1:num_dst)
end subroutine hadronic_forward_source_per_gamma

! 二级粒子分布映射：从壳层内每 gamma 数量变换为局部每 GeV 数密度谱。
subroutine hadronic_forward_distribution_per_gev(num_src,num_dst,energy_src_gev,density_src_per_gamma,energy_dst_gev, &
                                            mass_gev,shell_volume_cm3,density_dst_per_gev)
    implicit none
    integer, intent(in) :: num_src,num_dst
    real(8), intent(in) :: energy_src_gev(num_src),density_src_per_gamma(num_src),energy_dst_gev(num_dst)
    real(8), intent(in) :: mass_gev,shell_volume_cm3
    real(8), intent(out) :: density_dst_per_gev(num_dst)

    if (mass_gev <= 0d0) error stop "hadronic distribution per GeV requires positive particle mass."
    if (shell_volume_cm3 <= 0d0) error stop "hadronic distribution per GeV requires positive shell volume."
    call hadronic_forward_positive_loglog_interp(num_src,num_dst,energy_src_gev,density_src_per_gamma, &
                                            energy_dst_gev,density_dst_per_gev)
    density_dst_per_gev(1:num_dst)=density_dst_per_gev(1:num_dst)/(shell_volume_cm3*mass_gev)
end subroutine hadronic_forward_distribution_per_gev

! 按 hadron log spacing 构造对齐 photon energy grid。
subroutine hadronic_forward_aligned_photon_grid(num_had,num_ph,num_out,hadron_energy_gev,photon_energy_gev, &
                                           aligned_photon_gev)
    implicit none
    integer, intent(in) :: num_had,num_ph,num_out
    real(8), intent(in) :: hadron_energy_gev(num_had),photon_energy_gev(num_ph)
    real(8), intent(out) :: aligned_photon_gev(num_out)
    integer :: i
    real(8) :: dln_had,log_min

    if (num_had < 2 .or. num_ph < 2 .or. num_out < 1) error stop "hadronic aligned photon grid needs valid grids."
    dln_had=dlog(hadron_energy_gev(2)/hadron_energy_gev(1))
    log_min=dlog(photon_energy_gev(1))
    do i=1,num_out
        aligned_photon_gev(i)=dexp(log_min+dln_had*dble(i-1))
    end do
end subroutine hadronic_forward_aligned_photon_grid

subroutine hadronic_sequence_shell_geometry(num_r,radius_cm,gamma_bulk,i_r,dr,dt_s)
    use constants
    implicit none
    integer, intent(in) :: num_r,i_r
    real(8), intent(in) :: radius_cm(num_r),gamma_bulk(num_r)
    real(8), intent(out) :: dr,dt_s
    real(8) :: beta

    if (gamma_bulk(i_r) <= one) error stop "hadronic sequence shell dt requires gamma_bulk > 1."
    if (i_r == 1) then
        if (num_r < 2) error stop "hadronic sequence shell dt requires at least two radii."
        dr=radius_cm(2)-radius_cm(1)
    else
        dr=radius_cm(i_r)-radius_cm(i_r-1)
    end if
    if (dr <= zero) error stop "hadronic sequence shell radii must be strictly increasing."
    beta=dsqrt(one-one/(gamma_bulk(i_r)*gamma_bulk(i_r)))
    dt_s=dr/(beta*gamma_bulk(i_r)*Para_c)
end subroutine hadronic_sequence_shell_geometry

! 半径壳层体积：第一个壳层以内边界 R=0，之后使用相邻半径。
subroutine hadronic_forward_shell_volumes(num_r,radius_cm,shell_volume_cm3)
    use constants
    implicit none
    integer, intent(in) :: num_r
    real(8), intent(in) :: radius_cm(num_r)
    real(8), intent(out) :: shell_volume_cm3(num_r)
    integer :: i
    real(8) :: r_prev

    if (num_r < 1) error stop "hadronic shell volumes require at least one radius."
    if (radius_cm(1) <= zero) error stop "hadronic shell radii must be positive."
    shell_volume_cm3(1)=(4d0/3d0)*pi*radius_cm(1)**3
    r_prev=radius_cm(1)
    do i=2,num_r
        if (radius_cm(i) <= zero) error stop "hadronic shell radii must be positive."
        if (radius_cm(i) <= radius_cm(i-1)) error stop "hadronic shell radii must be strictly increasing."
        shell_volume_cm3(i)=(4d0/3d0)*pi*(radius_cm(i)**3-r_prev**3)
        r_prev=radius_cm(i)
    end do
end subroutine hadronic_forward_shell_volumes

! 壳层 comoving 时间步长 wrapper：复用统一的半径壳层几何定义。
subroutine hadronic_forward_shell_comoving_dt(num_r,radius_cm,gamma_bulk,shell_index,dt_s)
    implicit none
    integer, intent(in) :: num_r,shell_index
    real(8), intent(in) :: radius_cm(num_r),gamma_bulk(num_r)
    real(8), intent(out) :: dt_s
    real(8) :: dr

    call hadronic_sequence_shell_geometry(num_r,radius_cm,gamma_bulk,shell_index,dr,dt_s)
end subroutine hadronic_forward_shell_comoving_dt

! 动力学时间 wrapper：t_dyn=R/(Gamma c)。
subroutine hadronic_forward_dynamical_time(radius_cm,gamma_bulk,t_dyn_s)
    use constants
    implicit none
    real(8), intent(in) :: radius_cm,gamma_bulk
    real(8), intent(out) :: t_dyn_s

    if (radius_cm <= zero) error stop "hadronic dynamical time requires positive radius."
    if (gamma_bulk < one) error stop "hadronic dynamical time requires gamma_bulk >= 1."
    t_dyn_s=radius_cm/(gamma_bulk*Para_c)
end subroutine hadronic_forward_dynamical_time

! Quantum synchrotron cooling-factor array wrapper.
subroutine hadronic_forward_quantum_syn_cooling_factor(num_gamma,gamma,b_field_g,mass_gev,factor)
    use hadronic_common, only: hadronic_quantum_syn_cooling_factor
    implicit none
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: gamma(num_gamma),b_field_g,mass_gev
    real(8), intent(out) :: factor(num_gamma)
    integer :: i

    do i=1,num_gamma
        factor(i) = hadronic_quantum_syn_cooling_factor(gamma(i),b_field_g,mass_gev)
    end do
end subroutine hadronic_forward_quantum_syn_cooling_factor

end module hadronic_forward_shell_1d
