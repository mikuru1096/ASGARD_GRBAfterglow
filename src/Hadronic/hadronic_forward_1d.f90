! Single-shell proton synchrotron wrapper.
subroutine fs_hadronic_proton_syn_shell(R_loc,B_field_g,Num_gam_p,Num_nu,gam_p,dN_gam_p,V_seed,P_had_syn,Seed_had_syn)
    use hadronic_radiation_kernel
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p,Num_nu
    real(8), intent(in) :: R_loc,B_field_g,gam_p(Num_gam_p),dN_gam_p(Num_gam_p),V_seed(Num_nu)
    real(8), intent(out) :: P_had_syn(Num_nu),Seed_had_syn(Num_nu)

    call hadronic_get_proton_syn_state(R_loc,B_field_g,Num_gam_p,Num_nu,gam_p,dN_gam_p,V_seed,P_had_syn,Seed_had_syn)
end subroutine fs_hadronic_proton_syn_shell

! Single-shell charged-hadron synchrotron polarization wrapper.
subroutine fs_hadronic_syn_polarization_shell(Num_had,hadron_energy_gev,density_per_gev,Num_ph, &
                                              photon_frequency_hz,particle_mass_gev,magnetic_field_g,p_index,Pi_nu)
    use hadronic_radiation_kernel, only: hadronic_syn_polarization_fraction
    implicit none
    integer, intent(in) :: Num_had,Num_ph
    real(8), intent(in) :: hadron_energy_gev(Num_had),density_per_gev(Num_had),photon_frequency_hz(Num_ph)
    real(8), intent(in) :: particle_mass_gev,magnetic_field_g,p_index
    real(8), intent(out) :: Pi_nu(Num_ph)

    call hadronic_syn_polarization_fraction(Num_had,hadron_energy_gev,density_per_gev,Num_ph, &
                                            photon_frequency_hz,particle_mass_gev,magnetic_field_g,p_index,Pi_nu)
end subroutine fs_hadronic_syn_polarization_shell

! Single-shell Hummer 2010 p-gamma operator wrapper.
subroutine fs_hadronic_pgamma_operator_shell(Num_gam_p,Num_nu,hadron_energy_gev,hadron_density_per_gev,photon_energy_gev, &
                                             photon_density_per_gev,neutron_density_per_gev,pion0_source_rate_per_gev, &
                                             pion_plus_source_rate_per_gev,pion_minus_source_rate_per_gev, &
                                             proton_reinjection_rate_per_gev,neutron_reinjection_rate_per_gev, &
                                             proton_loss_rate,neutron_loss_rate,photon_loss_rate)
    use hadronic_interaction_kernel
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p,Num_nu
    real(8), intent(in) :: hadron_energy_gev(Num_gam_p),hadron_density_per_gev(Num_gam_p),photon_energy_gev(Num_nu)
    real(8), intent(in) :: photon_density_per_gev(Num_nu),neutron_density_per_gev(Num_gam_p)
    real(8), intent(out) :: pion0_source_rate_per_gev(Num_gam_p),pion_plus_source_rate_per_gev(Num_gam_p)
    real(8), intent(out) :: pion_minus_source_rate_per_gev(Num_gam_p),proton_reinjection_rate_per_gev(Num_gam_p)
    real(8), intent(out) :: neutron_reinjection_rate_per_gev(Num_gam_p),proton_loss_rate(Num_gam_p)
    real(8), intent(out) :: neutron_loss_rate(Num_gam_p),photon_loss_rate(Num_nu)

    call hadronic_pg_hummer2010_operator(Num_gam_p,Num_nu,hadron_energy_gev,hadron_density_per_gev,photon_energy_gev, &
                                         photon_density_per_gev,neutron_density_per_gev,pion0_source_rate_per_gev, &
                                         pion_plus_source_rate_per_gev,pion_minus_source_rate_per_gev, &
                                         proton_reinjection_rate_per_gev,neutron_reinjection_rate_per_gev, &
                                         proton_loss_rate,neutron_loss_rate,photon_loss_rate)
end subroutine fs_hadronic_pgamma_operator_shell

! Single-shell gamma-gamma pair-production wrapper.
subroutine fs_hadronic_pair_production_shell(Num_gamma,photon_energy_gev,photon_density_per_gev,Num_e,electron_energy_gev, &
                                             max_com_energy_factor,photon_loss_rate,pair_injection_rate_per_gev_per_species, &
                                             pair_injection_rate_per_gev_total,absorbed_power_gev_per_cm3_s, &
                                             injected_power_gev_per_cm3_s)
    use hadronic_pair_production_kernel, only: hadronic_pair_production_operator
    implicit none
    integer, intent(in) :: Num_gamma,Num_e,max_com_energy_factor
    real(8), intent(in) :: photon_energy_gev(Num_gamma),photon_density_per_gev(Num_gamma),electron_energy_gev(Num_e)
    real(8), intent(out) :: photon_loss_rate(Num_gamma),pair_injection_rate_per_gev_per_species(Num_e)
    real(8), intent(out) :: pair_injection_rate_per_gev_total(Num_e),absorbed_power_gev_per_cm3_s
    real(8), intent(out) :: injected_power_gev_per_cm3_s

    call hadronic_pair_production_operator(Num_gamma,photon_energy_gev,photon_density_per_gev,Num_e,electron_energy_gev, &
                                           max_com_energy_factor,photon_loss_rate,pair_injection_rate_per_gev_per_species, &
                                           pair_injection_rate_per_gev_total,absorbed_power_gev_per_cm3_s, &
                                           injected_power_gev_per_cm3_s)
end subroutine fs_hadronic_pair_production_shell

! Single-shell pp delta-approximation wrapper.
subroutine fs_hadronic_pp_delta_shell(Num_p,proton_energy_gev,proton_density_per_gev,target_proton_density_cm3, &
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
end subroutine fs_hadronic_pp_delta_shell

! Single-shell Bethe-Heitler wrapper.
subroutine fs_hadronic_bethe_heitler_shell(Num_p,proton_energy_gev,proton_density_per_gev,Num_ph,photon_energy_gev, &
                                           photon_density_per_gev,Num_e,electron_energy_gev,pair_rate_per_gev, &
                                           proton_loss_rate,photon_loss_rate)
    use hadronic_bethe_heitler_kernel, only: hadronic_bethe_heitler_operator
    implicit none
    integer, intent(in) :: Num_p,Num_ph,Num_e
    real(8), intent(in) :: proton_energy_gev(Num_p),proton_density_per_gev(Num_p)
    real(8), intent(in) :: photon_energy_gev(Num_ph),photon_density_per_gev(Num_ph)
    real(8), intent(in) :: electron_energy_gev(Num_e)
    real(8), intent(out) :: pair_rate_per_gev(Num_e),proton_loss_rate(Num_p),photon_loss_rate(Num_ph)

    call hadronic_bethe_heitler_operator(Num_p,proton_energy_gev,proton_density_per_gev,Num_ph,photon_energy_gev, &
                                         photon_density_per_gev,Num_e,electron_energy_gev,pair_rate_per_gev, &
                                         proton_loss_rate,photon_loss_rate)
end subroutine fs_hadronic_bethe_heitler_shell

! Single-shell hadronic inverse-Compton wrapper.
subroutine fs_hadronic_hadronic_ic_shell(Num_had,hadron_energy_gev,Num_ph,photon_energy_gev,photons_on_had_grid_per_gev, &
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
end subroutine fs_hadronic_hadronic_ic_shell

! Single-shell secondary-species transport wrapper.
subroutine fs_hadronic_species_transport_shell(Num_gamma,gamma,dt_s,b_field_g,divergence_rate_s_inv, &
                                               neutron_prev,pion_plus_prev,pion_minus_prev,muon_minus_left_prev, &
                                               muon_minus_right_prev,muon_plus_left_prev,muon_plus_right_prev, &
                                               neutron_source_per_gamma_s,pion_plus_source_per_gamma_s, &
                                               pion_minus_source_per_gamma_s,muon_minus_left_source_per_gamma_s, &
                                               muon_minus_right_source_per_gamma_s,muon_plus_left_source_per_gamma_s, &
                                               muon_plus_right_source_per_gamma_s,neutron_next,pion_plus_next, &
                                               pion_minus_next,muon_minus_left_next,muon_minus_right_next, &
                                               muon_plus_left_next,muon_plus_right_next)
    use hadronic_species_transport_kernel, only: hadronic_species_advance_operator
    implicit none
    integer, intent(in) :: Num_gamma
    real(8), intent(in) :: gamma(Num_gamma),dt_s,b_field_g,divergence_rate_s_inv
    real(8), intent(in) :: neutron_prev(Num_gamma),pion_plus_prev(Num_gamma),pion_minus_prev(Num_gamma)
    real(8), intent(in) :: muon_minus_left_prev(Num_gamma),muon_minus_right_prev(Num_gamma)
    real(8), intent(in) :: muon_plus_left_prev(Num_gamma),muon_plus_right_prev(Num_gamma)
    real(8), intent(in) :: neutron_source_per_gamma_s(Num_gamma),pion_plus_source_per_gamma_s(Num_gamma)
    real(8), intent(in) :: pion_minus_source_per_gamma_s(Num_gamma),muon_minus_left_source_per_gamma_s(Num_gamma)
    real(8), intent(in) :: muon_minus_right_source_per_gamma_s(Num_gamma),muon_plus_left_source_per_gamma_s(Num_gamma)
    real(8), intent(in) :: muon_plus_right_source_per_gamma_s(Num_gamma)
    real(8), intent(out) :: neutron_next(Num_gamma),pion_plus_next(Num_gamma),pion_minus_next(Num_gamma)
    real(8), intent(out) :: muon_minus_left_next(Num_gamma),muon_minus_right_next(Num_gamma)
    real(8), intent(out) :: muon_plus_left_next(Num_gamma),muon_plus_right_next(Num_gamma)

    call hadronic_species_advance_operator(Num_gamma,gamma,dt_s,b_field_g,divergence_rate_s_inv, &
                                           neutron_prev,pion_plus_prev,pion_minus_prev,muon_minus_left_prev, &
                                           muon_minus_right_prev,muon_plus_left_prev,muon_plus_right_prev, &
                                           neutron_source_per_gamma_s,pion_plus_source_per_gamma_s, &
                                           pion_minus_source_per_gamma_s,muon_minus_left_source_per_gamma_s, &
                                           muon_minus_right_source_per_gamma_s,muon_plus_left_source_per_gamma_s, &
                                           muon_plus_right_source_per_gamma_s,neutron_next,pion_plus_next,pion_minus_next, &
                                           muon_minus_left_next,muon_minus_right_next,muon_plus_left_next,muon_plus_right_next)
end subroutine fs_hadronic_species_transport_shell

! Single-shell hadronic acceleration wrapper.
subroutine fs_hadronic_acceleration_shell(Num_gamma,species,gamma,b_field_g,eta_acc,luminosity_erg_s,spectral_index, &
                                          gamma_min,gamma_max_inj,gamma_cut,has_gamma_cut,radius_cm,gamma_bulk, &
                                          Num_gamma_scan,gamma_scan,external_cooling_rate,has_external_cooling, &
                                          t_acc,t_syn,q_inj,gamma_max,gamma_dyn,gamma_syn,gamma_ext,has_gamma_ext)
    use hadronic_acceleration_kernel, only: hadronic_acceleration_operator
    implicit none
    integer, intent(in) :: Num_gamma,Num_gamma_scan
    character(len=*), intent(in) :: species
    real(8), intent(in) :: gamma(Num_gamma),b_field_g,eta_acc,luminosity_erg_s,spectral_index
    real(8), intent(in) :: gamma_min,gamma_max_inj,gamma_cut,radius_cm,gamma_bulk
    real(8), intent(in) :: gamma_scan(Num_gamma_scan),external_cooling_rate(Num_gamma_scan)
    logical, intent(in) :: has_gamma_cut,has_external_cooling
    real(8), intent(out) :: t_acc(Num_gamma),t_syn(Num_gamma),q_inj(Num_gamma)
    real(8), intent(out) :: gamma_max,gamma_dyn,gamma_syn,gamma_ext
    logical, intent(out) :: has_gamma_ext

    call hadronic_acceleration_operator(Num_gamma,species,gamma,b_field_g,eta_acc,luminosity_erg_s,spectral_index, &
                                        gamma_min,gamma_max_inj,gamma_cut,has_gamma_cut,radius_cm,gamma_bulk, &
                                        Num_gamma_scan,gamma_scan,external_cooling_rate,has_external_cooling, &
                                        t_acc,t_syn,q_inj,gamma_max,gamma_dyn,gamma_syn,gamma_ext,has_gamma_ext)
end subroutine fs_hadronic_acceleration_shell

! Single-shell secondary-radiation wrapper.
subroutine fs_hadronic_secondary_radiation_shell(Num_had,hadron_energy_gev,Num_ph,photon_energy_gev,pion_plus_per_gev, &
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
end subroutine fs_hadronic_secondary_radiation_shell

! Single-shell Hummer 2010 decay wrapper.
subroutine fs_hadronic_decay_operator_shell(Num_gam_p,hadron_energy_gev,pion0_source_rate_per_gev,pion_plus_source_rate_per_gev, &
                                            pion_minus_source_rate_per_gev,Num_gamma,gamma_energy_gev,Num_nu,neutrino_energy_gev, &
                                            Num_proc,process_energy_gev,gamma_rate_per_gev,process_rate_per_gev, &
                                            muon_plus_right_source_rate_per_gev,muon_plus_left_source_rate_per_gev, &
                                            muon_minus_left_source_rate_per_gev,muon_minus_right_source_rate_per_gev, &
                                            prompt_pion_neutrino_rate_per_gev,muon_neutrino_rate_per_gev, &
                                            muon_electron_rate_per_gev,neutrino_rate_per_gev)
    use hadronic_decay_kernel, only: hadronic_hummer2010_decay_operator
    implicit none
    integer, intent(in) :: Num_gam_p,Num_gamma,Num_nu,Num_proc
    real(8), intent(in) :: hadron_energy_gev(Num_gam_p),pion0_source_rate_per_gev(Num_gam_p)
    real(8), intent(in) :: pion_plus_source_rate_per_gev(Num_gam_p),pion_minus_source_rate_per_gev(Num_gam_p)
    real(8), intent(in) :: gamma_energy_gev(Num_gamma),neutrino_energy_gev(Num_nu)
    real(8), intent(in) :: process_energy_gev(Num_proc)
    real(8), intent(out) :: gamma_rate_per_gev(Num_gamma),process_rate_per_gev(3,Num_proc)
    real(8), intent(out) :: muon_plus_right_source_rate_per_gev(Num_gam_p)
    real(8), intent(out) :: muon_plus_left_source_rate_per_gev(Num_gam_p)
    real(8), intent(out) :: muon_minus_left_source_rate_per_gev(Num_gam_p)
    real(8), intent(out) :: muon_minus_right_source_rate_per_gev(Num_gam_p)
    real(8), intent(out) :: prompt_pion_neutrino_rate_per_gev(Num_nu),muon_neutrino_rate_per_gev(Num_nu)
    real(8), intent(out) :: muon_electron_rate_per_gev(Num_proc),neutrino_rate_per_gev(Num_nu)

    call hadronic_hummer2010_decay_operator(Num_gam_p,hadron_energy_gev,pion0_source_rate_per_gev,pion_plus_source_rate_per_gev, &
                                            pion_minus_source_rate_per_gev,Num_gamma,gamma_energy_gev,Num_nu,neutrino_energy_gev, &
                                            Num_proc,process_energy_gev,gamma_rate_per_gev,process_rate_per_gev, &
                                            muon_plus_right_source_rate_per_gev,muon_plus_left_source_rate_per_gev, &
                                            muon_minus_left_source_rate_per_gev,muon_minus_right_source_rate_per_gev, &
                                            prompt_pion_neutrino_rate_per_gev,muon_neutrino_rate_per_gev, &
                                            muon_electron_rate_per_gev,neutrino_rate_per_gev)
end subroutine fs_hadronic_decay_operator_shell

! Forward-shock 1D hadronic driver.
subroutine fs_hadronic_1d(R_Tobs,R_Gamma,R,shell_energy_inj_erg,B_field_g,V_seed,Seed_target,p_p,epsilon_p,eta_acc, &
                          include_proton_synch,include_pg,include_neutrino,Num_nu,Num_R,num_gam_p,num_nu_nu,n_threads, &
                          gam_p,dN_gam_p,P_had_syn,Seed_had_syn,P_had_pg_gamma,V_nu,P_nu_all)
    use constants
    use hadronic_common
    use hadronic_transport_kernel
    use hadronic_radiation_kernel
    use hadronic_interaction_kernel
    use hadronic_decay_kernel
    use hadronic_acceleration_kernel, only: hadronic_species_injection_operator, hadronic_estimate_max_gamma
    use hadronic_pgamma_hummer_1d, only: hadronic_pgamma_hummer_shell
    use hadronic_transport_remap_kernel, only: hadronic_advance_energy_loggamma_remap
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: include_proton_synch,include_pg,include_neutrino,Num_nu,Num_R,num_gam_p,num_nu_nu,n_threads
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),shell_energy_inj_erg(Num_R),B_field_g(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: Seed_target(Num_nu,Num_R)
    real(8), intent(in) :: p_p,epsilon_p,eta_acc
    real(8), intent(out) :: gam_p(num_gam_p),dN_gam_p(num_gam_p,Num_R),P_had_syn(Num_nu,Num_R)
    real(8), intent(out) :: Seed_had_syn(Num_nu,Num_R),P_had_pg_gamma(Num_nu,Num_R)
    real(8), intent(out) :: V_nu(num_nu_nu),P_nu_all(num_nu_nu,Num_R)
    integer :: I_R
    real(8) :: gam_p_max_global,t_dyn_s,dt_s,gam_p_min,energy_budget_erg,shell_volume
    real(8) :: gamma_max_loc,gamma_dyn,gamma_syn,gamma_ext,gamma_scan(2),external_rate(2)
    logical :: has_gamma_ext
    real(8), allocatable :: dN_prev(:),dN_next(:),Q_inj(:),loss_ad(:),loss_syn(:),loss_total(:),dN_trial(:),surv_pg(:)
    real(8), allocatable :: n_prev(:),n_next(:),pip_prev(:),pip_next(:),pim_prev(:),pim_next(:)
    real(8), allocatable :: muml_prev(:),muml_next(:),mumr_prev(:),mumr_next(:),mupl_prev(:),mupl_next(:)
    real(8), allocatable :: mupr_prev(:),mupr_next(:)
    logical :: use_hummer

    allocate(dN_prev(num_gam_p),dN_next(num_gam_p),dN_trial(num_gam_p),Q_inj(num_gam_p),loss_ad(num_gam_p), & 
             loss_syn(num_gam_p),loss_total(num_gam_p),surv_pg(Num_nu))
    use_hummer=(include_pg /= 0 .or. include_neutrino /= 0)
    if (use_hummer) then
        allocate(n_prev(num_gam_p),n_next(num_gam_p),pip_prev(num_gam_p),pip_next(num_gam_p), &
                 pim_prev(num_gam_p),pim_next(num_gam_p),muml_prev(num_gam_p),muml_next(num_gam_p), &
                 mumr_prev(num_gam_p),mumr_next(num_gam_p),mupl_prev(num_gam_p),mupl_next(num_gam_p), &
                 mupr_prev(num_gam_p),mupr_next(num_gam_p))
        n_prev=zero; pip_prev=zero; pim_prev=zero; muml_prev=zero; mumr_prev=zero; mupl_prev=zero; mupr_prev=zero
    end if

    call initialize_proton_gamma_grid
    dN_prev=zero

    call initialize_output_grids

    do I_R=1,Num_R
        dt_s=hadronic_shell_dt(R_Tobs,I_R)
        t_dyn_s=hadronic_dynamical_time(R(I_R),R_Gamma(I_R))
        if (shell_energy_inj_erg(I_R) < zero) error stop "hadronic shell injection energy must be non-negative."
        energy_budget_erg=shell_energy_inj_erg(I_R)
        gam_p_min=max(gam_p(1),R_Gamma(I_R))
        call inject_protons_for_shell
        call advance_proton_transport_for_shell
        dN_next=dN_trial

        surv_pg=one
        if (use_hummer) then
            call advance_hummer_secondary_chain
        end if
        dN_gam_p(:,I_R)=dN_next

        if (include_proton_synch /= 0) then
            call emit_proton_synchrotron_for_shell
        end if

        dN_prev=dN_next
    end do

    if (use_hummer) deallocate(n_prev,n_next,pip_prev,pip_next,pim_prev,pim_next,muml_prev,muml_next, &
                               mumr_prev,mumr_next,mupl_prev,mupl_next,mupr_prev,mupr_next)
    deallocate(dN_prev,dN_next,dN_trial,Q_inj,loss_ad,loss_syn,loss_total,surv_pg)

contains

    subroutine initialize_proton_gamma_grid
        implicit real(8)(A-H,O-Z)
        integer :: I_scan

        gamma_scan=(/one,two/)
        external_rate=(/one,two/)
        t_dyn_s=hadronic_dynamical_time(R(1),R_Gamma(1))
        if (use_hummer) then
            call hadronic_estimate_max_gamma("proton",B_field_g(1),R(1),R_Gamma(1),eta_acc,2,gamma_scan,external_rate, &
                                             .false.,gam_p_max_global,gamma_dyn,gamma_syn,gamma_ext,has_gamma_ext)
        else
            gam_p_max_global=hadronic_gamma_p_max(B_field_g(1),t_dyn_s,eta_acc)
        end if
        do I_scan=2,Num_R
            t_dyn_s=hadronic_dynamical_time(R(I_scan),R_Gamma(I_scan))
            if (use_hummer) then
                call hadronic_estimate_max_gamma("proton",B_field_g(I_scan),R(I_scan),R_Gamma(I_scan),eta_acc,2,gamma_scan, &
                                                 external_rate,.false.,gamma_max_loc,gamma_dyn,gamma_syn,gamma_ext, &
                                                 has_gamma_ext)
            else
                gamma_max_loc=hadronic_gamma_p_max(B_field_g(I_scan),t_dyn_s,eta_acc)
            end if
            gam_p_max_global=max(gam_p_max_global,gamma_max_loc)
        end do
        if (gam_p_max_global <= one+1d-3) error stop "forward hadronic gamma_p_max must exceed the injection grid minimum."
        call hadronic_build_gamma_p_grid(num_gam_p,one+1d-3,gam_p_max_global,gam_p)
    end subroutine initialize_proton_gamma_grid

    subroutine initialize_output_grids
        implicit real(8)(A-H,O-Z)
        integer :: I_nu_out

        V_nu=zero
        if (num_nu_nu > 1) then
            do I_nu_out=1,num_nu_nu
                V_nu(I_nu_out)=ten**(dlog10(1d-3*Para_GeV2Hz)+dble(I_nu_out-1)* &
                                (dlog10(1d8*Para_GeV2Hz)-dlog10(1d-3*Para_GeV2Hz))/dble(num_nu_nu-1))
            end do
        else
            V_nu(1)=Para_GeV2Hz
        end if

        dN_gam_p=zero
        P_had_syn=zero
        Seed_had_syn=zero
        P_had_pg_gamma=zero
        P_nu_all=zero
    end subroutine initialize_output_grids

    subroutine inject_protons_for_shell
        implicit real(8)(A-H,O-Z)

        if (use_hummer .and. energy_budget_erg > zero) then
            call hadronic_species_injection_operator(num_gam_p,gam_p,"proton",energy_budget_erg/dt_s,p_p, &
                                                     gam_p_min,gam_p(num_gam_p),one,.false.,Q_inj)
            Q_inj=dt_s*Q_inj
        else
            call hadronic_proton_injection_powerlaw(num_gam_p,gam_p,p_p,energy_budget_erg,gam_p_min,gam_p(num_gam_p),Q_inj)
        end if
    end subroutine inject_protons_for_shell

    subroutine advance_proton_transport_for_shell
        implicit real(8)(A-H,O-Z)

        call hadronic_proton_loss_rates(num_gam_p,gam_p,B_field_g(I_R),t_dyn_s,loss_ad,loss_syn,loss_total)
        if (use_hummer) then
            call hadronic_advance_energy_loggamma_remap(num_gam_p,gam_p,dN_prev,Q_inj,loss_total,dt_s,dN_trial)
        else
            call hadronic_advance_energy_loggamma(num_gam_p,gam_p,dN_prev,Q_inj,loss_total,dt_s,dN_trial)
        end if
    end subroutine advance_proton_transport_for_shell

    subroutine advance_hummer_secondary_chain
        implicit real(8)(A-H,O-Z)

        if (I_R == 1) then
            shell_volume=(4d0/3d0)*pi*R(I_R)**3
        else
            shell_volume=(4d0/3d0)*pi*(R(I_R)**3-R(I_R-1)**3)
        end if
        call hadronic_pgamma_hummer_shell(num_gam_p,Num_nu,num_nu_nu,dt_s,R(I_R),R_Gamma(I_R),B_field_g(I_R), &
                                          shell_volume,gam_p,dN_trial,V_seed,Seed_target(:,I_R),V_nu, &
                                          include_pg,include_neutrino, &
                                          n_prev,pip_prev,pim_prev,muml_prev,mumr_prev,mupl_prev,mupr_prev, &
                                          dN_next,n_next,pip_next,pim_next,muml_next,mumr_next,mupl_next,mupr_next, &
                                          surv_pg,P_had_pg_gamma(:,I_R),P_nu_all(:,I_R))
        n_prev=n_next; pip_prev=pip_next; pim_prev=pim_next
        muml_prev=muml_next; mumr_prev=mumr_next; mupl_prev=mupl_next; mupr_prev=mupr_next
    end subroutine advance_hummer_secondary_chain

    subroutine emit_proton_synchrotron_for_shell
        implicit real(8)(A-H,O-Z)

        call hadronic_get_proton_syn_state(R(I_R),B_field_g(I_R),num_gam_p,Num_nu,gam_p,dN_next, &
                                           V_seed,P_had_syn(:,I_R),Seed_had_syn(:,I_R))
        if (use_hummer) then
            P_had_syn(:,I_R)=P_had_syn(:,I_R)*surv_pg
            Seed_had_syn(:,I_R)=Seed_had_syn(:,I_R)*surv_pg
        end if
    end subroutine emit_proton_synchrotron_for_shell
end subroutine fs_hadronic_1d

! Single-step electromagnetic pair-cascade wrapper.
subroutine fs_hadronic_pair_cascade_step(num_ph,photon_energy_gev,photon_density, &
                                          num_e,electron_energy_gev,b_field_g,path_time_s, &
                                          cascade_syn_spec,photon_loss_rate,absorbed_power)
    use hadronic_pair_cascade_kernel, only: hadronic_cascade_step
    implicit none
    integer, intent(in) :: num_ph,num_e
    real(8), intent(in) :: photon_energy_gev(num_ph),photon_density(num_ph)
    real(8), intent(in) :: electron_energy_gev(num_e),b_field_g,path_time_s
    real(8), intent(out) :: cascade_syn_spec(num_ph),photon_loss_rate(num_ph),absorbed_power

    call hadronic_cascade_step(num_ph,photon_energy_gev,photon_density, &
                                num_e,electron_energy_gev,b_field_g,path_time_s, &
                                cascade_syn_spec,photon_loss_rate,absorbed_power)
end subroutine fs_hadronic_pair_cascade_step

subroutine fs_hadronic_pair_cascade_sequence(num_ph,num_e,num_shell,photon_energy_gev,primary_photon_density, &
                                             electron_energy_gev,frequency_hz,radius_cm,gamma_bulk,observer_time_s, &
                                             b_field_g,num_threads,index_syn_integr,substeps_per_shell, &
                                             photon_loss_rate,tau_pair,pair_density,pair_luminosity, &
                                             pair_seed,cascade_photon_density,absorbed_power,injected_power)
    use hadronic_pair_cascade_kernel, only: hadronic_cascade_sequence
    implicit none
    integer, intent(in) :: num_ph,num_e,num_shell,num_threads,index_syn_integr,substeps_per_shell
    real(8), intent(in) :: photon_energy_gev(num_ph),primary_photon_density(num_ph,num_shell)
    real(8), intent(in) :: electron_energy_gev(num_e),frequency_hz(num_ph)
    real(8), intent(in) :: radius_cm(num_shell),gamma_bulk(num_shell),observer_time_s(num_shell),b_field_g(num_shell)
    real(8), intent(out) :: photon_loss_rate(num_ph,num_shell),tau_pair(num_ph,num_shell)
    real(8), intent(out) :: pair_density(num_e,num_shell),pair_luminosity(num_ph,num_shell)
    real(8), intent(out) :: pair_seed(num_ph,num_shell),cascade_photon_density(num_ph,num_shell)
    real(8), intent(out) :: absorbed_power(num_shell),injected_power(num_shell)

    call hadronic_cascade_sequence(num_ph,num_e,num_shell,photon_energy_gev,primary_photon_density, &
                                   electron_energy_gev,frequency_hz,radius_cm,gamma_bulk,observer_time_s, &
                                   b_field_g,num_threads,index_syn_integr,substeps_per_shell, &
                                   photon_loss_rate,tau_pair,pair_density,pair_luminosity, &
                                   pair_seed,cascade_photon_density,absorbed_power,injected_power)
end subroutine fs_hadronic_pair_cascade_sequence

! pp spectral source model: SIBYLL=0, QGSJET=1, Geant4=2, Pythia8=3.
subroutine fs_hadronic_pp_spectral_source(num_p,proton_kinetic_energy_gev, &
    proton_density_per_gev,num_g,gamma_energy_gev,target_density_cm3,model, &
    pi0_gamma_rate)
    use hadronic_pp_models_kernel, only: hadronic_pp_pi0_source_spectrum
    implicit none
    integer, intent(in) :: num_p,num_g,model
    real(8), intent(in) :: proton_kinetic_energy_gev(num_p),proton_density_per_gev(num_p)
    real(8), intent(in) :: gamma_energy_gev(num_g),target_density_cm3
    real(8), intent(out) :: pi0_gamma_rate(num_g)

    call hadronic_pp_pi0_source_spectrum(num_p,proton_kinetic_energy_gev, &
        proton_density_per_gev,num_g,gamma_energy_gev,target_density_cm3,model, &
        pi0_gamma_rate)
end subroutine fs_hadronic_pp_spectral_source

! Quantum synchrotron cooling-factor array wrapper.
subroutine fs_hadronic_quantum_syn_cooling_factor(num_gamma,gamma,b_field_g,mass_gev,factor)
    use hadronic_common, only: hadronic_quantum_syn_cooling_factor
    implicit none
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: gamma(num_gamma),b_field_g,mass_gev
    real(8), intent(out) :: factor(num_gamma)
    integer :: i

    do i=1,num_gamma
        factor(i) = hadronic_quantum_syn_cooling_factor(gamma(i),b_field_g,mass_gev)
    end do
end subroutine fs_hadronic_quantum_syn_cooling_factor
