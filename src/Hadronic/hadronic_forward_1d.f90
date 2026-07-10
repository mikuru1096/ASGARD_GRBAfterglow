! Python process wrappers kept only for runtime-facing hadronic glue.
! Python/f2py ABI wrapper; the shell implementation lives in hadronic_shell.
subroutine had_syn_pol(Num_had,ehad,density_per_gev,Num_ph, &
                                              photon_frequency_hz,particle_mass_gev,magnetic_field_g,p_index,Pi_nu)
    use hadronic_rad, only: syn_polarization
    implicit none
    integer, intent(in) :: Num_had,Num_ph
    real(8), intent(in), dimension(Num_had) :: ehad,density_per_gev
    real(8), intent(in), dimension(Num_ph) :: photon_frequency_hz
    real(8), intent(in) :: particle_mass_gev,magnetic_field_g,p_index
    real(8), intent(out), dimension(Num_ph) :: Pi_nu

    call syn_polarization(Num_had, ehad, density_per_gev, Num_ph, &
        photon_frequency_hz, particle_mass_gev, magnetic_field_g, p_index, Pi_nu)
end subroutine had_syn_pol
! Python/f2py ABI wrapper; the shell implementation lives in hadronic_shell.
subroutine pg_operator(ngp,nnu,ehad,hden,eph, &
                                             phden,nden,qpi0, &
                                             qpip,qpim, &
                                             qpro,qntr, &
                                             ploss,nloss,phloss)
    use hadronic_pg, only: pg_hummer
    implicit none
    integer, intent(in) :: ngp,nnu
    real(8), intent(in), dimension(ngp) :: ehad,hden
    real(8), intent(in), dimension(nnu) :: eph
    real(8), intent(in), dimension(nnu) :: phden
    real(8), intent(in), dimension(ngp) :: nden
    real(8), intent(out), dimension(ngp) :: qpi0,qpip
    real(8), intent(out), dimension(ngp) :: qpim,qpro
    real(8), intent(out), dimension(ngp) :: qntr,ploss
    real(8), intent(out), dimension(ngp) :: nloss
    real(8), intent(out), dimension(nnu) :: phloss

    call pg_hummer(ngp, nnu, ehad, hden, &
        eph, phden, nden, qpi0, &
        qpip, qpim, qpro, &
        qntr, ploss, nloss, phloss)
end subroutine pg_operator
! Python/f2py ABI wrapper; the shell implementation lives in hadronic_shell.
subroutine pair_production(Num_gamma,eph,phden,Num_e,electron_energy_gev, &
                                             max_com_energy_factor,phloss,pair_injection_rate_per_gev_per_species, &
                                             pair_injection_rate_per_gev_total,absorbed_power_gev_per_cm3_s, &
                                             injected_power_gev_per_cm3_s)
    use hadronic_pair, only: pair_operator
    implicit none
    integer, intent(in) :: Num_gamma,Num_e,max_com_energy_factor
    real(8), intent(in), dimension(Num_gamma) :: eph,phden
    real(8), intent(in), dimension(Num_e) :: electron_energy_gev
    real(8), intent(out), dimension(Num_gamma) :: phloss
    real(8), intent(out), dimension(Num_e) :: pair_injection_rate_per_gev_per_species
    real(8), intent(out), dimension(Num_e) :: pair_injection_rate_per_gev_total
    real(8), intent(out) :: absorbed_power_gev_per_cm3_s
    real(8), intent(out) :: injected_power_gev_per_cm3_s

    call pair_operator(Num_gamma, eph, phden, Num_e, &
        electron_energy_gev, max_com_energy_factor, phloss, pair_injection_rate_per_gev_per_species, &
        pair_injection_rate_per_gev_total, absorbed_power_gev_per_cm3_s, injected_power_gev_per_cm3_s)
end subroutine pair_production
! Python/f2py ABI wrapper; the shell implementation lives in hadronic_shell.
subroutine pp_shell(Num_p,proton_energy_gev,proton_density_per_gev,target_proton_density_cm3, &
                                      Num_gamma,gamma_energy_gev,nnu,neutrino_energy_gev,Num_pair,pair_energy_gev, &
                                      kappa_inelastic,pion_energy_fraction,neutral_pion_fraction,gamma_rate_per_gev, &
                                      neutrino_rate_per_gev,pair_rate_per_gev,ploss)
    use hadronic_shell, only: pp_delta
    implicit none
    integer, intent(in) :: Num_p,Num_gamma,nnu,Num_pair
    real(8), intent(in), dimension(Num_p) :: proton_energy_gev,proton_density_per_gev
    real(8), intent(in) :: target_proton_density_cm3
    real(8), intent(in), dimension(Num_gamma) :: gamma_energy_gev
    real(8), intent(in), dimension(nnu) :: neutrino_energy_gev
    real(8), intent(in), dimension(Num_pair) :: pair_energy_gev
    real(8), intent(in) :: kappa_inelastic,pion_energy_fraction,neutral_pion_fraction
    real(8), intent(out), dimension(Num_gamma) :: gamma_rate_per_gev
    real(8), intent(out), dimension(nnu) :: neutrino_rate_per_gev
    real(8), intent(out), dimension(Num_pair) :: pair_rate_per_gev
    real(8), intent(out), dimension(Num_p) :: ploss

    call pp_delta(Num_p, proton_energy_gev, proton_density_per_gev, &
        target_proton_density_cm3, Num_gamma, gamma_energy_gev, nnu, neutrino_energy_gev, Num_pair, &
        pair_energy_gev, kappa_inelastic, pion_energy_fraction, neutral_pion_fraction, gamma_rate_per_gev, &
        neutrino_rate_per_gev, pair_rate_per_gev, ploss)
end subroutine pp_shell
! Python/f2py ABI wrapper; the shell implementation lives in hadronic_shell.
subroutine bethe_heitler(Num_p,proton_energy_gev,proton_density_per_gev,Num_ph,eph, &
                                           phden,Num_e,electron_energy_gev,pair_rate_per_gev, &
                                           ploss,phloss)
    use hadronic_bh, only: bh_calc
    implicit none
    integer, intent(in) :: Num_p,Num_ph,Num_e
    real(8), intent(in), dimension(Num_p) :: proton_energy_gev,proton_density_per_gev
    real(8), intent(in), dimension(Num_ph) :: eph,phden
    real(8), intent(in), dimension(Num_e) :: electron_energy_gev
    real(8), intent(out), dimension(Num_e) :: pair_rate_per_gev
    real(8), intent(out), dimension(Num_p) :: ploss
    real(8), intent(out), dimension(Num_ph) :: phloss

    call bh_calc(Num_p, proton_energy_gev, proton_density_per_gev, Num_ph, &
        eph, phden, Num_e, electron_energy_gev, pair_rate_per_gev, &
        ploss, phloss)
end subroutine bethe_heitler
! Python/f2py ABI wrapper; the shell implementation lives in hadronic_shell.
subroutine hadronic_ic(Num_had,ehad,Num_ph,eph,photons_on_had_grid_per_gev, &
                                         protons_per_gev,pion_plus_per_gev,pion_minus_per_gev,muon_minus_left_per_gev, &
                                         muon_minus_right_per_gev,muon_plus_left_per_gev,muon_plus_right_per_gev, &
                                         ind_min_energy_pho_hadgrid,epsilon_p_ic,epsilon_pi_ic,epsilon_mu_ic, &
                                         coeff_p_cgs,coeff_pi_cgs,coeff_mu_cgs,dln_energy,delta_e_p,jmax_p,delta_e_pi, &
                                         jmax_pi,delta_e_mu,jmax_mu)
    use hadronic_shell, only: hic_shell
    implicit none
    integer, intent(in) :: Num_had,Num_ph,ind_min_energy_pho_hadgrid
    real(8), intent(in), dimension(Num_had) :: ehad
    real(8), intent(in), dimension(Num_ph) :: eph,photons_on_had_grid_per_gev
    real(8), intent(in), dimension(Num_had) :: protons_per_gev,pion_plus_per_gev,pion_minus_per_gev
    real(8), intent(in), dimension(Num_had) :: muon_minus_left_per_gev,muon_minus_right_per_gev
    real(8), intent(in), dimension(Num_had) :: muon_plus_left_per_gev,muon_plus_right_per_gev
    real(8), intent(out), dimension(Num_ph) :: epsilon_p_ic,epsilon_pi_ic,epsilon_mu_ic
    real(8), intent(out) :: coeff_p_cgs,coeff_pi_cgs,coeff_mu_cgs,dln_energy
    integer, intent(out), dimension(Num_had) :: delta_e_p,jmax_p,delta_e_pi,jmax_pi
    integer, intent(out), dimension(Num_had) :: delta_e_mu,jmax_mu

    call hic_shell(Num_had, ehad, Num_ph, eph, &
        photons_on_had_grid_per_gev, protons_per_gev, pion_plus_per_gev, pion_minus_per_gev, &
        muon_minus_left_per_gev, muon_minus_right_per_gev, muon_plus_left_per_gev, muon_plus_right_per_gev, &
        ind_min_energy_pho_hadgrid, epsilon_p_ic, epsilon_pi_ic, epsilon_mu_ic, coeff_p_cgs, coeff_pi_cgs, &
        coeff_mu_cgs, dln_energy, delta_e_p, jmax_p, delta_e_pi, jmax_pi, delta_e_mu, jmax_mu)
end subroutine hadronic_ic
! Python/f2py ABI wrapper; the shell implementation lives in hadronic_shell.
subroutine decay_operator(ngp,ehad,qpi0,qpip, &
                                            qpim,Num_gamma,gamma_energy_gev,nnu,neutrino_energy_gev, &
                                            Num_proc,process_energy_gev,gamma_rate_per_gev,process_rate_per_gev, &
                                            muon_plus_right_source_rate_per_gev,muon_plus_left_source_rate_per_gev, &
                                            muon_minus_left_source_rate_per_gev,muon_minus_right_source_rate_per_gev, &
                                            prompt_pion_neutrino_rate_per_gev,muon_neutrino_rate_per_gev, &
                                            muon_electron_rate_per_gev,neutrino_rate_per_gev)
    use hadronic_decay, only: decay_hummer
    implicit none
    integer, intent(in) :: ngp,Num_gamma,nnu,Num_proc
    real(8), intent(in), dimension(ngp) :: ehad,qpi0
    real(8), intent(in), dimension(ngp) :: qpip,qpim
    real(8), intent(in), dimension(Num_gamma) :: gamma_energy_gev
    real(8), intent(in), dimension(nnu) :: neutrino_energy_gev
    real(8), intent(in), dimension(Num_proc) :: process_energy_gev
    real(8), intent(out), dimension(Num_gamma) :: gamma_rate_per_gev
    real(8), intent(out), dimension(3,Num_proc) :: process_rate_per_gev
    real(8), intent(out), dimension(ngp) :: muon_plus_right_source_rate_per_gev
    real(8), intent(out), dimension(ngp) :: muon_plus_left_source_rate_per_gev
    real(8), intent(out), dimension(ngp) :: muon_minus_left_source_rate_per_gev
    real(8), intent(out), dimension(ngp) :: muon_minus_right_source_rate_per_gev
    real(8), intent(out), dimension(nnu) :: prompt_pion_neutrino_rate_per_gev,muon_neutrino_rate_per_gev
    real(8), intent(out), dimension(Num_proc) :: muon_electron_rate_per_gev
    real(8), intent(out), dimension(nnu) :: neutrino_rate_per_gev

    call decay_hummer(ngp, ehad, qpi0, &
        qpip, qpim, Num_gamma, gamma_energy_gev, nnu, &
        neutrino_energy_gev, Num_proc, process_energy_gev, gamma_rate_per_gev, process_rate_per_gev, &
        muon_plus_right_source_rate_per_gev, muon_plus_left_source_rate_per_gev, &
        muon_minus_left_source_rate_per_gev, muon_minus_right_source_rate_per_gev, &
        prompt_pion_neutrino_rate_per_gev, muon_neutrino_rate_per_gev, muon_electron_rate_per_gev, &
        neutrino_rate_per_gev)
end subroutine decay_operator

! Python/f2py ABI wrapper for the shell-sequence pair cascade used by runtime pair-cascade mode.
subroutine cascade_sequence(nph,ne,num_shell,eph,primary_photon_density, &
                                             electron_energy_gev,frequency_hz,radius_cm,gamma_bulk,observer_time_s, &
                                             B_field_g,num_threads,index_syn_integr,substeps_per_shell, &
                                             phloss,tau_pair,pair_density,pair_luminosity, &
                                             pair_seed,cascade_photon_density,absorbed_power,injected_power)
    use hadronic_cascade, only: cascade_seq
    implicit none
    integer, intent(in) :: nph,ne,num_shell,num_threads,index_syn_integr,substeps_per_shell
    real(8), intent(in), dimension(nph) :: eph
    real(8), intent(in), dimension(nph,num_shell) :: primary_photon_density
    real(8), intent(in), dimension(ne) :: electron_energy_gev
    real(8), intent(in), dimension(nph) :: frequency_hz
    real(8), intent(in), dimension(num_shell) :: radius_cm,gamma_bulk,observer_time_s,B_field_g
    real(8), intent(out), dimension(nph,num_shell) :: phloss,tau_pair
    real(8), intent(out), dimension(ne,num_shell) :: pair_density
    real(8), intent(out), dimension(nph,num_shell) :: pair_luminosity
    real(8), intent(out), dimension(nph,num_shell) :: pair_seed,cascade_photon_density
    real(8), intent(out), dimension(num_shell) :: absorbed_power,injected_power

    call cascade_seq(nph, ne, num_shell, eph, &
        primary_photon_density, electron_energy_gev, frequency_hz, radius_cm, gamma_bulk, observer_time_s, &
        B_field_g, num_threads, index_syn_integr, substeps_per_shell, phloss, tau_pair, pair_density, &
        pair_luminosity, pair_seed, cascade_photon_density, absorbed_power, injected_power)
end subroutine cascade_sequence

! Public shell-sequence drivers.
! Forward-shock 1D hadronic light driver.
! 顺序: build proton grid -> loop shells: inject protons -> transport losses
!       -> optional Hummer p-gamma secondary chain -> proton synchrotron emission.
subroutine hadronic_1d(R_Tobs,R_Gamma,R,shell_energy_inj_erg,B_field_g,V_seed,Seed_target,p_p,epsilon_p,eta_acc, &
                          include_proton_synch,include_pg,include_neutrino,nnu,Num_R,num_gam_p,num_nu_nu,n_threads, &
                          gam_p,dN_gam_p,P_had_syn,Seed_had_syn,P_had_pg_gamma,V_nu,P_nu_all)
    use constants
    use hadronic_base
    use hadronic_transport
    use hadronic_rad, only: proton_syn
    use hadronic_pg, only: pg_hummer
    use hadronic_decay
    use hadronic_accel, only: gamma_limit
    use hadronic_hummer, only: hummer_shell
    use hadronic_remap, only: remap_loggamma
    use hadronic_shell, only: shell_geom
    implicit none
    integer, intent(in) :: include_proton_synch,include_pg,include_neutrino,nnu,Num_R,num_gam_p,num_nu_nu,n_threads
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R,shell_energy_inj_erg,B_field_g
    real(8), intent(in), dimension(nnu) :: V_seed
    real(8), intent(in), dimension(nnu,Num_R) :: Seed_target
    real(8), intent(in) :: p_p,epsilon_p,eta_acc
    real(8), intent(out), dimension(num_gam_p) :: gam_p
    real(8), intent(out), dimension(num_gam_p,Num_R) :: dN_gam_p
    real(8), intent(out), dimension(nnu,Num_R) :: P_had_syn
    real(8), intent(out), dimension(nnu,Num_R) :: Seed_had_syn,P_had_pg_gamma
    real(8), intent(out), dimension(num_nu_nu) :: V_nu
    real(8), intent(out), dimension(num_nu_nu,Num_R) :: P_nu_all
    integer :: ir
    real(8) :: gmax,tdyn,dt,dr,gmin,ebudget,volume
    real(8), dimension(2) :: gscan,xrate
    real(8) :: gloc,gdyn,gsyn,gext
    logical :: has_extlim
    real(8), allocatable, dimension(:) :: dnprev,dnnext,qinj,lad,lsyn,ltot,dntry,surv
    real(8), allocatable, dimension(:) :: nprev,nnext,pipp,pipn,pimp,pimn
    real(8), allocatable, dimension(:) :: mumlp,mumln,mumrp,mumrn,muplp,mupln
    real(8), allocatable, dimension(:) :: muprp,muprn
    logical :: usehum

    allocate(dnprev(num_gam_p),dnnext(num_gam_p),dntry(num_gam_p),qinj(num_gam_p),lad(num_gam_p), &
             lsyn(num_gam_p),ltot(num_gam_p),surv(nnu))
    usehum=(include_pg /= 0 .or. include_neutrino /= 0)
    if (usehum) then
        allocate(nprev(num_gam_p),nnext(num_gam_p),pipp(num_gam_p),pipn(num_gam_p), &
                 pimp(num_gam_p),pimn(num_gam_p),mumlp(num_gam_p),mumln(num_gam_p), &
                 mumrp(num_gam_p),mumrn(num_gam_p),muplp(num_gam_p),mupln(num_gam_p), &
                 muprp(num_gam_p),muprn(num_gam_p))
        nprev=0d0; pipp=0d0; pimp=0d0; mumlp=0d0; mumrp=0d0; muplp=0d0; muprp=0d0
    end if

    ! 初始化质子能量网格和输出网格。
    ! Initialize the proton energy grid and output arrays.
    call init_grid
    dnprev=0d0

    call init_out

    ! 壳层循环：注入 -> 质子输运 -> p-gamma secondary -> 同步辐射。
    ! Shell loop: injection -> proton transport -> p-gamma secondaries -> synchrotron.
    do ir=1,Num_R
        call shell_geom(Num_R,R,R_Gamma,ir,dr,dt)
        tdyn=dyn_time(R(ir),R_Gamma(ir))
        if (shell_energy_inj_erg(ir) < 0d0) error stop "hadronic shell injection energy must be non-negative."
        ebudget=shell_energy_inj_erg(ir)
        gmin=max(gam_p(1),R_Gamma(ir))
        call inject_p
        call advance_p
        dnnext=dntry

        surv=1d0
        if (usehum) then
            call advance_sec
        end if
        dN_gam_p(:,ir)=dnnext

        if (include_proton_synch /= 0) then
            call emit_syn
        end if

        dnprev=dnnext
    end do

    ! 释放只在 Hummer secondary 链中需要的状态缓存。
    ! Release state caches used only by the Hummer secondary chain.
    if (usehum) deallocate(nprev,nnext,pipp,pipn,pimp,pimn,mumlp,mumln, &
                               mumrp,mumrn,muplp,mupln,muprp,muprn)
    deallocate(dnprev,dnnext,dntry,qinj,lad,lsyn,ltot,surv)

contains

    subroutine init_grid
        integer :: iscan

        gscan=[1d0,2d0]
        xrate=[1d0,2d0]
        tdyn=dyn_time(R(1),R_Gamma(1))
        if (usehum) then
            call gamma_limit("proton",B_field_g(1),R(1),R_Gamma(1),eta_acc,2,gscan,xrate, &
                                             .false.,gmax,gdyn,gsyn,gext,has_extlim)
        else
            gmax=proton_limit(B_field_g(1),tdyn,eta_acc)
        end if
        do iscan=2,Num_R
            tdyn=dyn_time(R(iscan),R_Gamma(iscan))
            if (usehum) then
                call gamma_limit("proton",B_field_g(iscan),R(iscan),R_Gamma(iscan),eta_acc,2,gscan, &
                                                 xrate,.false.,gloc,gdyn,gsyn,gext, &
                                                 has_extlim)
            else
                gloc=proton_limit(B_field_g(iscan),tdyn,eta_acc)
            end if
            gmax=max(gmax,gloc)
        end do
        if (gmax <= 1d0+1d-3) error stop "forward hadronic gp_max must exceed the injection grid minimum."
        call build_grid(num_gam_p,1d0+1d-3,gmax,gam_p)
    end subroutine init_grid

    subroutine init_out
        integer :: inu

        V_nu=0d0
        if (num_nu_nu > 1) then
            do inu=1,num_nu_nu
                V_nu(inu)=1d1**(dlog10(1d-3*Para_GeV2Hz)+dble(inu-1)* &
                                (dlog10(1d8*Para_GeV2Hz)-dlog10(1d-3*Para_GeV2Hz))/dble(num_nu_nu-1))
            end do
        else
            V_nu(1)=Para_GeV2Hz
        end if

        dN_gam_p=0d0
        P_had_syn=0d0
        Seed_had_syn=0d0
        P_had_pg_gamma=0d0
        P_nu_all=0d0
    end subroutine init_out

    subroutine inject_p
        call proton_inject(num_gam_p,gam_p,p_p,ebudget,gmin,gam_p(num_gam_p),qinj)
    end subroutine inject_p

    subroutine advance_p
        call proton_loss(num_gam_p,gam_p,B_field_g(ir),tdyn,lad,lsyn,ltot)
        if (usehum) then
            call remap_loggamma(num_gam_p,gam_p,dnprev,qinj,ltot,dt,dntry)
        else
            call advance_loggamma(num_gam_p,gam_p,dnprev,qinj,ltot,dt,dntry)
        end if
    end subroutine advance_p

    subroutine advance_sec
        if (ir == 1) then
            volume=(4d0/3d0)*pi*R(ir)**3
        else
            volume=(4d0/3d0)*pi*(R(ir)**3-R(ir-1)**3)
        end if
        call hummer_shell(num_gam_p,nnu,num_nu_nu,dt,R(ir),R_Gamma(ir),B_field_g(ir), &
                                          volume,gam_p,dntry,V_seed,Seed_target(:,ir),V_nu, &
                                          include_pg,include_neutrino, &
                                          nprev,pipp,pimp,mumlp,mumrp,muplp,muprp, &
                                          dnnext,nnext,pipn,pimn,mumln,mumrn,mupln,muprn, &
                                          surv,P_had_pg_gamma(:,ir),P_nu_all(:,ir))
        nprev=nnext; pipp=pipn; pimp=pimn
        mumlp=mumln; mumrp=mumrn; muplp=mupln; muprp=muprn
    end subroutine advance_sec

    subroutine emit_syn
        call proton_syn(R(ir),B_field_g(ir),num_gam_p,nnu,gam_p,dnnext, &
                                           V_seed,P_had_syn(:,ir),Seed_had_syn(:,ir))
        if (usehum) then
            P_had_syn(:,ir)=P_had_syn(:,ir)*surv
            Seed_had_syn(:,ir)=Seed_had_syn(:,ir)*surv
        end if
    end subroutine emit_syn
end subroutine hadronic_1d

! Formal 1D hadronic shell-sequence ABI wrapper.
! The implementation advances proton transport, p-gamma/BH/pp interactions, secondary species,
! secondary radiation, photon survival, and secondary electron source on the shell grid.
subroutine formal_transport_1d(R_Tobs,R_Gamma,R,B_field_g,V_seed,Seed_target,gamma_e, &
        shell_energy_inj_erg,pp_target_density_cm3,p_p,eta_acc,index_syn_integr,include_proton_synch, &
        include_pg,include_neutrino,include_bethe_heitler,include_hadronic_inverse_compton,include_pp, &
        quantum_syn,n_threads,nnu,Num_R,Num_e,num_gam_p,num_nu_nu,gam_p,gam_secondary,dN_gam_p, &
        P_had_syn,Seed_had_syn,P_had_pg_gamma,V_nu,P_nu_all,P_had_bh,Seed_had_bh,dN_gam_e_bh, &
        secondary_electron_source_r,tau_bh,bh_phloss,P_had_hic,dN_gam_n,dN_gam_pi_plus, &
        dN_gam_pi_minus,dN_gam_mu_minus_left,dN_gam_mu_minus_right,dN_gam_mu_plus_left, &
        dN_gam_mu_plus_right,P_had_pion_synch,P_had_muon_synch,P_had_pion_ic,P_had_muon_ic, &
        tau_pg,pg_photon_survival,am3_process_power)
    use hadronic_formal, only: formal_transport
    implicit none
    integer, intent(in) :: index_syn_integr,include_proton_synch,include_pg,include_neutrino
    integer, intent(in) :: include_bethe_heitler,include_hadronic_inverse_compton,include_pp,quantum_syn,n_threads
    integer, intent(in) :: nnu,Num_R,Num_e,num_gam_p,num_nu_nu
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R,B_field_g
    real(8), intent(in), dimension(nnu) :: V_seed
    real(8), intent(in), dimension(nnu,Num_R) :: Seed_target
    real(8), intent(in), dimension(Num_e) :: gamma_e
    real(8), intent(in), dimension(Num_R) :: shell_energy_inj_erg
    real(8), intent(in), dimension(Num_R) :: pp_target_density_cm3
    real(8), intent(in) :: p_p,eta_acc
    real(8), intent(out), dimension(num_gam_p) :: gam_p,gam_secondary
    real(8), intent(out), dimension(num_gam_p,Num_R) :: dN_gam_p
    real(8), intent(out), dimension(nnu,Num_R) :: P_had_syn,Seed_had_syn,P_had_pg_gamma
    real(8), intent(out), dimension(num_nu_nu) :: V_nu
    real(8), intent(out), dimension(num_nu_nu,Num_R) :: P_nu_all
    real(8), intent(out), dimension(nnu,Num_R) :: P_had_bh
    real(8), intent(out), dimension(nnu,Num_R) :: Seed_had_bh
    real(8), intent(out), dimension(Num_e,Num_R) :: dN_gam_e_bh
    real(8), intent(out), dimension(Num_e,Num_R) :: secondary_electron_source_r
    real(8), intent(out), dimension(nnu,Num_R) :: tau_bh
    real(8), intent(out), dimension(nnu,Num_R) :: bh_phloss,P_had_hic
    real(8), intent(out), dimension(num_gam_p,Num_R) :: dN_gam_n,dN_gam_pi_plus
    real(8), intent(out), dimension(num_gam_p,Num_R) :: dN_gam_pi_minus,dN_gam_mu_minus_left
    real(8), intent(out), dimension(num_gam_p,Num_R) :: dN_gam_mu_minus_right,dN_gam_mu_plus_left
    real(8), intent(out), dimension(num_gam_p,Num_R) :: dN_gam_mu_plus_right
    real(8), intent(out), dimension(nnu,Num_R) :: P_had_pion_synch
    real(8), intent(out), dimension(nnu,Num_R) :: P_had_muon_synch,P_had_pion_ic
    real(8), intent(out), dimension(nnu,Num_R) :: P_had_muon_ic,tau_pg
    real(8), intent(out), dimension(nnu,Num_R) :: pg_photon_survival
    real(8), intent(out), dimension(3,num_gam_p,Num_R) :: am3_process_power

    call formal_transport(R_Tobs,R_Gamma,R,B_field_g,V_seed,Seed_target,gamma_e, &
        shell_energy_inj_erg,pp_target_density_cm3,p_p,eta_acc,index_syn_integr,include_proton_synch,include_pg, &
        include_neutrino,include_bethe_heitler,include_hadronic_inverse_compton,include_pp,quantum_syn, &
        n_threads,nnu,Num_R,Num_e,num_gam_p,num_nu_nu,gam_p,gam_secondary,dN_gam_p,P_had_syn,Seed_had_syn, &
        P_had_pg_gamma,V_nu,P_nu_all,P_had_bh,Seed_had_bh,dN_gam_e_bh,secondary_electron_source_r,tau_bh, &
        bh_phloss,P_had_hic,dN_gam_n,dN_gam_pi_plus,dN_gam_pi_minus,dN_gam_mu_minus_left, &
        dN_gam_mu_minus_right,dN_gam_mu_plus_left,dN_gam_mu_plus_right,P_had_pion_synch,P_had_muon_synch, &
        P_had_pion_ic,P_had_muon_ic,tau_pg,pg_photon_survival,am3_process_power)
end subroutine formal_transport_1d
