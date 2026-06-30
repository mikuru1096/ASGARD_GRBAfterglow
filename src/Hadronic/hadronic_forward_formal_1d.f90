module hadronic_forward_formal_1d
    use hadronic_forward_shell_1d
    implicit none
    private
    public :: hadronic_forward_formal_transport_1d_impl
contains

! Formal 1D hadronic shell sequence: Python passes arrays; Fortran owns the shell transport order.
subroutine hadronic_forward_formal_transport_1d_impl(R_Tobs,R_Gamma,R,B_field_g,V_seed,Seed_target,gamma_e, &
        shell_energy_inj_erg,pp_target_density_cm3,p_p,eta_acc,index_syn_integr,include_proton_synch, &
        include_pg,include_neutrino,include_bethe_heitler,include_hadronic_inverse_compton,include_pp, &
        quantum_syn,n_threads,Num_nu,Num_R,Num_e,num_gam_p,num_nu_nu,gam_p,gam_secondary,dN_gam_p, &
        P_had_syn,Seed_had_syn,P_had_pg_gamma,V_nu,P_nu_all,P_had_bh,Seed_had_bh,dN_gam_e_bh, &
        secondary_electron_source_r,tau_bh,bh_photon_loss_rate,P_had_hic,dN_gam_n,dN_gam_pi_plus, &
        dN_gam_pi_minus,dN_gam_mu_minus_left,dN_gam_mu_minus_right,dN_gam_mu_plus_left, &
        dN_gam_mu_plus_right,P_had_pion_synch,P_had_muon_synch,P_had_pion_ic,P_had_muon_ic, &
        tau_pg,pg_photon_survival,am3_process_power)
    use constants
    use hadronic_common, only: hadronic_build_gamma_p_grid
    use hadronic_bethe_heitler_kernel, only: hadronic_bethe_heitler_operator
    use hadronic_decay_kernel, only: hadronic_hummer2010_decay_operator
    use hadronic_interaction_kernel, only: hadronic_pg_hummer2010_operator
    use hadronic_radiation_kernel, only: hadronic_get_proton_syn_state
    implicit none
    integer, intent(in) :: index_syn_integr,include_proton_synch,include_pg,include_neutrino
    integer, intent(in) :: include_bethe_heitler,include_hadronic_inverse_compton,include_pp,quantum_syn,n_threads
    integer, intent(in) :: Num_nu,Num_R,Num_e,num_gam_p,num_nu_nu
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),B_field_g(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: Seed_target(Num_nu,Num_R),gamma_e(Num_e),shell_energy_inj_erg(Num_R)
    real(8), intent(in) :: pp_target_density_cm3(Num_R),p_p,eta_acc
    real(8), intent(out) :: gam_p(num_gam_p),gam_secondary(num_gam_p),dN_gam_p(num_gam_p,Num_R)
    real(8), intent(out) :: P_had_syn(Num_nu,Num_R),Seed_had_syn(Num_nu,Num_R),P_had_pg_gamma(Num_nu,Num_R)
    real(8), intent(out) :: V_nu(num_nu_nu),P_nu_all(num_nu_nu,Num_R),P_had_bh(Num_nu,Num_R)
    real(8), intent(out) :: Seed_had_bh(Num_nu,Num_R),dN_gam_e_bh(Num_e,Num_R)
    real(8), intent(out) :: secondary_electron_source_r(Num_e,Num_R),tau_bh(Num_nu,Num_R)
    real(8), intent(out) :: bh_photon_loss_rate(Num_nu,Num_R),P_had_hic(Num_nu,Num_R)
    real(8), intent(out) :: dN_gam_n(num_gam_p,Num_R),dN_gam_pi_plus(num_gam_p,Num_R)
    real(8), intent(out) :: dN_gam_pi_minus(num_gam_p,Num_R),dN_gam_mu_minus_left(num_gam_p,Num_R)
    real(8), intent(out) :: dN_gam_mu_minus_right(num_gam_p,Num_R),dN_gam_mu_plus_left(num_gam_p,Num_R)
    real(8), intent(out) :: dN_gam_mu_plus_right(num_gam_p,Num_R),P_had_pion_synch(Num_nu,Num_R)
    real(8), intent(out) :: P_had_muon_synch(Num_nu,Num_R),P_had_pion_ic(Num_nu,Num_R)
    real(8), intent(out) :: P_had_muon_ic(Num_nu,Num_R),tau_pg(Num_nu,Num_R)
    real(8), intent(out) :: pg_photon_survival(Num_nu,Num_R),am3_process_power(3,num_gam_p,Num_R)
    integer :: i_r,i_nu_out,num_align
    real(8) :: dt_s,t_dyn_s,dr_shell,gam_p_min,gam_p_max_global,shell_volume,divergence_rate,dln_had
    real(8) :: shell_volumes(Num_R)
    real(8) :: hadron_energy(num_gam_p),photon_energy(Num_nu),electron_energy(Num_e),zero_rate(num_gam_p)
    real(8) :: neutrino_energy(num_nu_nu)
    real(8) :: dN_prev(num_gam_p),dN_next(num_gam_p),dN_trial(num_gam_p),Q_inj(num_gam_p)
    real(8) :: proton_density(num_gam_p),neutron_density(num_gam_p),photon_density_tau(Num_nu)
    real(8) :: photon_density(Num_nu),pion0_source(num_gam_p),pip_source(num_gam_p),pim_source(num_gam_p)
    real(8) :: pg_reinj(num_gam_p),neutron_reinj(num_gam_p),pg_loss(num_gam_p),neutron_loss(num_gam_p)
    real(8) :: pg_photon_loss(Num_nu),gamma_rate(Num_nu),process_rate(3,Num_nu),nu_rate(num_nu_nu)
    real(8) :: mupr_source(num_gam_p),mupl_source(num_gam_p),muml_source(num_gam_p),mumr_source(num_gam_p)
    real(8) :: prompt_nu_rate(num_nu_nu),muon_nu_rate(num_nu_nu),muon_e_rate(Num_nu)
    real(8) :: bh_pair_rate(Num_e),bh_proton_loss(num_gam_p),bh_loss(num_gam_p)
    real(8) :: pp_gamma_rate(Num_nu),pp_nu_rate(num_nu_nu),pp_pair_rate(Num_e),pp_proton_loss(num_gam_p)
    real(8) :: pp_loss(num_gam_p),pp_gamma_lum(Num_nu),pp_nu_lum(num_nu_nu),q_secondary_e(Num_e,Num_R)
    real(8) :: n_prev(num_gam_p),n_next(num_gam_p),pip_prev(num_gam_p),pip_next(num_gam_p)
    real(8) :: pim_prev(num_gam_p),pim_next(num_gam_p),muml_prev(num_gam_p),muml_next(num_gam_p)
    real(8) :: mumr_prev(num_gam_p),mumr_next(num_gam_p),mupl_prev(num_gam_p),mupl_next(num_gam_p)
    real(8) :: mupr_prev(num_gam_p),mupr_next(num_gam_p)

    call hadronic_forward_global_gamma_p_max(Num_R,R,R_Gamma,B_field_g,eta_acc,gam_p_max_global)
    call hadronic_build_gamma_p_grid(num_gam_p,1d0+1d-3,gam_p_max_global,gam_p)
    gam_secondary=gam_p
    hadron_energy=gam_p*Para_m_p_GeV
    electron_energy=gamma_e*Para_m_e_GeV
    photon_energy=V_seed*Para_h_GeV
    zero_rate=0d0
    dN_prev=0d0
    n_prev=0d0; pip_prev=0d0; pim_prev=0d0
    muml_prev=0d0; mumr_prev=0d0; mupl_prev=0d0; mupr_prev=0d0
    dN_gam_p=0d0; P_had_syn=0d0; Seed_had_syn=0d0; P_had_pg_gamma=0d0; P_nu_all=0d0
    P_had_bh=0d0; Seed_had_bh=0d0; dN_gam_e_bh=0d0; secondary_electron_source_r=0d0
    tau_bh=0d0; bh_photon_loss_rate=0d0; P_had_hic=0d0; q_secondary_e=0d0
    dN_gam_n=0d0; dN_gam_pi_plus=0d0; dN_gam_pi_minus=0d0
    dN_gam_mu_minus_left=0d0; dN_gam_mu_minus_right=0d0; dN_gam_mu_plus_left=0d0
    dN_gam_mu_plus_right=0d0; P_had_pion_synch=0d0; P_had_muon_synch=0d0
    P_had_pion_ic=0d0; P_had_muon_ic=0d0; tau_pg=0d0; pg_photon_survival=1d0
    am3_process_power=0d0

    if (num_nu_nu > 1) then
        do i_nu_out=1,num_nu_nu
            V_nu(i_nu_out)=ten**(dlog10(1d-3*Para_GeV2Hz)+dble(i_nu_out-1)* &
                           (dlog10(1d8*Para_GeV2Hz)-dlog10(1d-3*Para_GeV2Hz))/dble(num_nu_nu-1))
        end do
    else
        V_nu(1)=1d-3*Para_GeV2Hz
    end if
    neutrino_energy=Para_h_GeV*V_nu

    dln_had=dlog(hadron_energy(2)/hadron_energy(1))
    num_align=int(ceiling((dlog(photon_energy(Num_nu))-dlog(photon_energy(1)))/dln_had))+1
    call hadronic_forward_shell_volumes(Num_R,R,shell_volumes)

    do i_r=1,Num_R
        call hadronic_sequence_shell_geometry(Num_R,R,R_Gamma,i_r,dr_shell,dt_s)
        call hadronic_forward_dynamical_time(R(i_r),R_Gamma(i_r),t_dyn_s)
        gam_p_min=max(gam_p(1),R_Gamma(i_r))
        shell_volume=shell_volumes(i_r)
        call hadronic_forward_injection_content(num_gam_p,"proton",gam_p,shell_energy_inj_erg(i_r),dt_s, &
                                           p_p,gam_p_min,gam_p(num_gam_p),1d0,.false.,Q_inj)
        call hadronic_forward_proton_transport_step(num_gam_p,gam_p,dN_prev,Q_inj,B_field_g(i_r), &
                                               t_dyn_s,Para_m_p_GeV,quantum_syn,zero_rate,zero_rate,zero_rate, &
                                               zero_rate,shell_volume,dt_s,dN_trial)
        call hadronic_forward_shell_density_per_gev(num_gam_p,dN_trial,Para_m_p_GeV,shell_volume,proton_density)
        call hadronic_forward_shell_density_per_gev(num_gam_p,n_prev,Para_m_n_GeV,shell_volume,neutron_density)

        photon_energy(1:Num_nu)=Para_h_GeV*V_seed(1:Num_nu)
        photon_density_tau(1:Num_nu)=Seed_target(1:Num_nu,i_r)/Para_h_GeV
        call hadronic_pg_hummer2010_operator(num_gam_p,Num_nu,hadron_energy,proton_density, &
                                             photon_energy,photon_density_tau,neutron_density,pion0_source, &
                                             pip_source,pim_source,pg_reinj,neutron_reinj,pg_loss, &
                                             neutron_loss,pg_photon_loss)
        call hadronic_hummer2010_decay_operator(num_gam_p,hadron_energy,pion0_source,pip_source,pim_source, &
                                                Num_nu,photon_energy,num_nu_nu,neutrino_energy,Num_nu, &
                                                photon_energy,gamma_rate,process_rate,mupr_source,mupl_source, &
                                                muml_source,mumr_source,prompt_nu_rate,muon_nu_rate, &
                                                muon_e_rate,nu_rate)
        call hadronic_forward_photon_loss_closure(Num_nu,Num_R,R,R_Gamma,i_r,pg_photon_loss, &
                                             tau_pg(:,i_r),pg_photon_survival(:,i_r))
        photon_density=photon_density_tau*pg_photon_survival(:,i_r)
        call hadronic_pg_hummer2010_operator(num_gam_p,Num_nu,hadron_energy,proton_density, &
                                             photon_energy,photon_density,neutron_density,pion0_source, &
                                             pip_source,pim_source,pg_reinj,neutron_reinj,pg_loss, &
                                             neutron_loss,pg_photon_loss)
        call hadronic_hummer2010_decay_operator(num_gam_p,hadron_energy,pion0_source,pip_source,pim_source, &
                                                Num_nu,photon_energy,num_nu_nu,neutrino_energy,Num_nu, &
                                                photon_energy,gamma_rate,process_rate,mupr_source,mupl_source, &
                                                muml_source,mumr_source,prompt_nu_rate,muon_nu_rate, &
                                                muon_e_rate,nu_rate)

        bh_loss=0d0; bh_pair_rate=0d0; bh_photon_loss_rate(:,i_r)=0d0
        if (include_bethe_heitler /= 0) then
            call hadronic_bethe_heitler_operator(num_gam_p,hadron_energy,proton_density,Num_nu, &
                                                 photon_energy,photon_density,Num_e,electron_energy, &
                                                 bh_pair_rate,bh_proton_loss,bh_photon_loss_rate(:,i_r))
            if (any(bh_proton_loss > 0d0)) error stop "Bethe-Heitler proton loss rate must be non-positive."
            bh_loss=-bh_proton_loss
            call hadronic_forward_photon_loss_closure(Num_nu,Num_R,R,R_Gamma,i_r,bh_photon_loss_rate(:,i_r), &
                                                 tau_bh(:,i_r),photon_density_tau)
        end if

        pp_loss=0d0; pp_pair_rate=0d0; pp_gamma_lum=0d0; pp_nu_lum=0d0
        if (include_pp /= 0) then
            call hadronic_forward_pp_delta_shell(num_gam_p,hadron_energy,proton_density,pp_target_density_cm3(i_r), &
                                            Num_nu,photon_energy,num_nu_nu,neutrino_energy,Num_e,electron_energy, &
                                            0.5d0,0.17d0,1d0/3d0,pp_gamma_rate,pp_nu_rate,pp_pair_rate, &
                                            pp_proton_loss)
            if (any(pp_proton_loss > 0d0)) error stop "pp proton loss rate must be non-positive."
            pp_loss=-pp_proton_loss
            call hadronic_forward_energy_luminosity_from_rate(Num_nu,photon_energy,pp_gamma_rate,shell_volume,pp_gamma_lum)
            call hadronic_forward_energy_luminosity_from_rate(num_nu_nu,neutrino_energy,pp_nu_rate,shell_volume,pp_nu_lum)
        end if

        call hadronic_forward_proton_transport_step(num_gam_p,gam_p,dN_prev,Q_inj,B_field_g(i_r), &
                                               t_dyn_s,Para_m_p_GeV,quantum_syn,bh_loss,pp_loss,pg_loss, &
                                               pg_reinj,shell_volume,dt_s,dN_next)
        dN_gam_p(:,i_r)=dN_next

        if (include_proton_synch /= 0) then
            call hadronic_get_proton_syn_state(R(i_r),B_field_g(i_r),num_gam_p,Num_nu,gam_p,dN_next,V_seed, &
                                               P_had_syn(:,i_r),Seed_had_syn(:,i_r))
        end if

        divergence_rate=3d0/t_dyn_s
        call hadronic_forward_species_transport_step(num_gam_p,num_gam_p,gam_p,hadron_energy,neutron_reinj, &
                                                pip_source,pim_source,muml_source,mumr_source,mupl_source, &
                                                mupr_source,neutron_loss,dt_s,B_field_g(i_r),divergence_rate, &
                                                shell_volume,n_prev,pip_prev,pim_prev,muml_prev,mumr_prev, &
                                                mupl_prev,mupr_prev,n_next,pip_next,pim_next,muml_next, &
                                                mumr_next,mupl_next,mupr_next)
        dN_gam_n(:,i_r)=n_next; dN_gam_pi_plus(:,i_r)=pip_next; dN_gam_pi_minus(:,i_r)=pim_next
        dN_gam_mu_minus_left(:,i_r)=muml_next; dN_gam_mu_minus_right(:,i_r)=mumr_next
        dN_gam_mu_plus_left(:,i_r)=mupl_next; dN_gam_mu_plus_right(:,i_r)=mupr_next

        call hadronic_forward_secondary_radiation_projected(num_gam_p,Num_nu,num_align,hadron_energy, &
                                                       photon_energy,photon_density,pip_next,pim_next,muml_next, &
                                                       mumr_next,mupl_next,mupr_next,shell_volume,B_field_g(i_r), &
                                                       P_had_pion_synch(:,i_r),P_had_muon_synch(:,i_r), &
                                                       P_had_pion_ic(:,i_r),P_had_muon_ic(:,i_r))
        if (include_pg /= 0) then
            call hadronic_forward_energy_luminosity_from_rate(Num_nu,photon_energy,gamma_rate,shell_volume, &
                                                         P_had_pg_gamma(:,i_r))
            if (include_pp /= 0) P_had_pg_gamma(:,i_r)=P_had_pg_gamma(:,i_r)+pp_gamma_lum
            call hadronic_forward_process_power(num_gam_p,Num_nu,3,hadron_energy,dN_next,photon_energy, &
                                           process_rate,shell_volume,am3_process_power(:,:,i_r))
        end if
        if (include_neutrino /= 0) then
            call hadronic_forward_energy_luminosity_from_rate(num_nu_nu,neutrino_energy,nu_rate,shell_volume,P_nu_all(:,i_r))
            if (include_pp /= 0) P_nu_all(:,i_r)=P_nu_all(:,i_r)+pp_nu_lum
        end if
        if (include_hadronic_inverse_compton /= 0) then
            call hadronic_forward_hic_projected(num_gam_p,Num_nu,num_align,hadron_energy,photon_energy, &
                                           photon_density,proton_density,shell_volume,P_had_hic(:,i_r))
        end if
        if (include_bethe_heitler /= 0 .or. include_pp /= 0) then
            call hadronic_forward_pair_source_content(Num_e,pp_pair_rate,bh_pair_rate,include_pp, &
                                                 include_bethe_heitler,shell_volume,q_secondary_e(:,i_r))
        end if

        dN_prev=dN_next
        n_prev=n_next; pip_prev=pip_next; pim_prev=pim_next
        muml_prev=muml_next; mumr_prev=mumr_next; mupl_prev=mupl_next; mupr_prev=mupr_next
    end do

    if (include_bethe_heitler /= 0 .or. include_pp /= 0) then
        call hadronic_forward_secondary_electron_sequence(Num_e,Num_nu,Num_R,gamma_e,R,R_Gamma,B_field_g, &
                                                     V_seed,q_secondary_e,index_syn_integr,n_threads,quantum_syn, &
                                                     dN_gam_e_bh,P_had_bh,Seed_had_bh,secondary_electron_source_r)
    end if
    P_had_syn=P_had_syn*pg_photon_survival
    Seed_had_syn=Seed_had_syn*pg_photon_survival
    P_had_pg_gamma=P_had_pg_gamma*pg_photon_survival
    if (include_bethe_heitler /= 0) then
        P_had_bh=P_had_bh*pg_photon_survival
        Seed_had_bh=Seed_had_bh*pg_photon_survival
    end if
    if (include_hadronic_inverse_compton /= 0) P_had_hic=P_had_hic*pg_photon_survival
    P_had_pion_synch=P_had_pion_synch*pg_photon_survival
    P_had_muon_synch=P_had_muon_synch*pg_photon_survival
    P_had_pion_ic=P_had_pion_ic*pg_photon_survival
    P_had_muon_ic=P_had_muon_ic*pg_photon_survival
end subroutine hadronic_forward_formal_transport_1d_impl

end module hadronic_forward_formal_1d
