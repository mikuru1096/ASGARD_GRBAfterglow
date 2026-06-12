!f2py: skip
module hadronic_pgamma_hummer_1d
    use constants
    use hadronic_common
    use hadronic_interaction_kernel, only: hadronic_pg_hummer2010_operator
    use hadronic_decay_kernel, only: hadronic_hummer2010_decay_operator
    use hadronic_species_transport_kernel, only: hadronic_species_advance_operator
    implicit none

contains

subroutine hadronic_pgamma_hummer_shell(num_gam_p,num_ph,num_nu,dt_s,radius_cm,gamma_bulk,b_field_g,shell_volume_cm3, &
                                        gam_p,dn_trial,photon_frequency_hz,seed_target_hz,neutrino_frequency_hz, &
                                        include_pg,include_neutrino,neutron_prev,pip_prev,pim_prev,muml_prev,mumr_prev, &
                                        mupl_prev,mupr_prev,dn_next,neutron_next,pip_next,pim_next,muml_next,mumr_next, &
                                        mupl_next,mupr_next,pg_survival,pg_gamma_lum,nu_lum)
    implicit none
    integer, intent(in) :: num_gam_p,num_ph,num_nu,include_pg,include_neutrino
    real(8), intent(in) :: dt_s,radius_cm,gamma_bulk,b_field_g,shell_volume_cm3
    real(8), intent(in) :: gam_p(num_gam_p),dn_trial(num_gam_p),photon_frequency_hz(num_ph),seed_target_hz(num_ph)
    real(8), intent(in) :: neutrino_frequency_hz(num_nu)
    real(8), intent(in) :: neutron_prev(num_gam_p),pip_prev(num_gam_p),pim_prev(num_gam_p)
    real(8), intent(in) :: muml_prev(num_gam_p),mumr_prev(num_gam_p),mupl_prev(num_gam_p),mupr_prev(num_gam_p)
    real(8), intent(out) :: dn_next(num_gam_p),neutron_next(num_gam_p),pip_next(num_gam_p),pim_next(num_gam_p)
    real(8), intent(out) :: muml_next(num_gam_p),mumr_next(num_gam_p),mupl_next(num_gam_p),mupr_next(num_gam_p)
    real(8), intent(out) :: pg_survival(num_ph),pg_gamma_lum(num_ph),nu_lum(num_nu)

    dn_next=dn_trial
    neutron_next=neutron_prev; pip_next=pip_prev; pim_next=pim_prev
    muml_next=muml_prev; mumr_next=mumr_prev; mupl_next=mupl_prev; mupr_next=mupr_prev
    pg_survival=one
    pg_gamma_lum=zero
    nu_lum=zero

    if (include_pg == 0 .and. include_neutrino == 0) return

    block
        real(8) :: hadron_energy(num_gam_p),proton_density(num_gam_p),neutron_density(num_gam_p)
        real(8) :: photon_energy(num_ph),photon_density_tau(num_ph),photon_density_att(num_ph),neutrino_energy(num_nu)
        real(8) :: qpi0(num_gam_p),qpip(num_gam_p),qpim(num_gam_p),preinj(num_gam_p),nreinj(num_gam_p)
        real(8) :: ploss(num_gam_p),nloss(num_gam_p),phloss(num_ph),tau_pg(num_ph),gamma_rate(num_ph)
        real(8) :: proc_rate(3,num_ph),mupr(num_gam_p),mupl(num_gam_p),mumnl(num_gam_p),mumnr(num_gam_p)
        real(8) :: prompt_nu(num_nu),muon_nu(num_nu),muon_e(num_ph),path_time_s,rate_dt,tau
        real(8) :: neutron_source(num_gam_p),pip_source(num_gam_p),pim_source(num_gam_p)
        real(8) :: muml_source(num_gam_p),mumr_source(num_gam_p),mupl_source(num_gam_p),mupr_source(num_gam_p)
        real(8) :: divergence_rate
        integer :: i

        hadron_energy=gam_p*hadronic_proton_mass_gev
        proton_density=dn_trial/(shell_volume_cm3*hadronic_proton_mass_gev)
        neutron_density=neutron_prev/(shell_volume_cm3*hadronic_neutron_mass_gev)
        photon_energy=photon_frequency_hz*Para_h_GeV
        photon_density_tau=seed_target_hz/Para_h_GeV
        neutrino_energy=neutrino_frequency_hz*Para_h_GeV

        call hadronic_pg_hummer2010_operator(num_gam_p,num_ph,hadron_energy,proton_density,photon_energy, &
                                             photon_density_tau,neutron_density,qpi0,qpip,qpim,preinj,nreinj, &
                                             ploss,nloss,phloss)
        path_time_s=radius_cm/(12d0*gamma_bulk*Para_c)
        tau_pg=phloss*path_time_s
        do i=1,num_ph
            tau=tau_pg(i)
            if (tau > 1d-6) then
                pg_survival(i)=(one-exp(-tau))/tau
            else if (tau > zero) then
                pg_survival(i)=one-tau/two+tau*tau/6d0
            end if
        end do

        photon_density_att=photon_density_tau*pg_survival
        call hadronic_pg_hummer2010_operator(num_gam_p,num_ph,hadron_energy,proton_density,photon_energy, &
                                             photon_density_att,neutron_density,qpi0,qpip,qpim,preinj,nreinj, &
                                             ploss,nloss,phloss)

        do i=1,num_gam_p
            associate(loss_rate => ploss(i), injection_rate => preinj(i), &
                      proton_mass_volume_factor => shell_volume_cm3*hadronic_proton_mass_gev)
                rate_dt=loss_rate*dt_s
                if (loss_rate > zero) then
                    dn_next(i)=dn_trial(i)*exp(-rate_dt)+(one-exp(-rate_dt))/loss_rate &
                               *proton_mass_volume_factor*injection_rate
                else
                    dn_next(i)=dn_trial(i)+dt_s*proton_mass_volume_factor*injection_rate
                end if
            end associate
        end do

        call hadronic_hummer2010_decay_operator(num_gam_p,hadron_energy,qpi0,qpip,qpim,num_ph,photon_energy, &
                                                num_nu,neutrino_energy,num_ph,photon_energy,gamma_rate,proc_rate, &
                                                mupr,mupl,mumnl,mumnr,prompt_nu,muon_nu,muon_e,nu_lum)
        call map_pgamma_secondary_sources(num_gam_p,gam_p,hadron_energy,shell_volume_cm3,nreinj,qpip,qpim, &
                                          mumnl,mumnr,mupl,mupr,neutron_source,pip_source,pim_source, &
                                          muml_source,mumr_source,mupl_source,mupr_source)
        divergence_rate=3d0*gamma_bulk*Para_c/radius_cm
        call hadronic_species_advance_operator(num_gam_p,gam_p,dt_s,b_field_g,divergence_rate, &
                                               neutron_prev,pip_prev,pim_prev,muml_prev,mumr_prev,mupl_prev,mupr_prev, &
                                               neutron_source,pip_source,pim_source,muml_source,mumr_source,mupl_source, &
                                               mupr_source,neutron_next,pip_next,pim_next,muml_next,mumr_next,mupl_next,mupr_next)
        call apply_neutron_pgamma_loss(num_gam_p,gam_p,dt_s,hadron_energy,nloss,neutron_next)
        if (include_pg /= 0) pg_gamma_lum=shell_volume_cm3*gamma_rate*photon_energy*Para_h_GeV*Para_GeV2erg*pg_survival
        if (include_neutrino /= 0) nu_lum=shell_volume_cm3*nu_lum*neutrino_energy*Para_h_GeV*Para_GeV2erg
    end block
end subroutine hadronic_pgamma_hummer_shell

subroutine map_pgamma_secondary_sources(num_gam_p,gam_p,hadron_energy,shell_volume_cm3,nreinj,qpip,qpim, &
                                        mumnl,mumnr,mupl,mupr,neutron_source,pip_source,pim_source, &
                                        muml_source,mumr_source,mupl_source,mupr_source)
    implicit none
    integer, intent(in) :: num_gam_p
    real(8), intent(in) :: gam_p(num_gam_p),hadron_energy(num_gam_p),shell_volume_cm3
    real(8), intent(in) :: nreinj(num_gam_p),qpip(num_gam_p),qpim(num_gam_p)
    real(8), intent(in) :: mumnl(num_gam_p),mumnr(num_gam_p),mupl(num_gam_p),mupr(num_gam_p)
    real(8), intent(out) :: neutron_source(num_gam_p),pip_source(num_gam_p),pim_source(num_gam_p)
    real(8), intent(out) :: muml_source(num_gam_p),mumr_source(num_gam_p),mupl_source(num_gam_p),mupr_source(num_gam_p)
    real(8) :: neutron_energy(num_gam_p),pion_energy(num_gam_p),muon_energy(num_gam_p)

    neutron_energy=gam_p*hadronic_neutron_mass_gev
    pion_energy=gam_p*hadronic_pion_charged_mass_gev
    muon_energy=gam_p*hadronic_muon_mass_gev
    call interp_source_per_gamma(num_gam_p,hadron_energy,nreinj,neutron_energy,hadronic_neutron_mass_gev, &
                                 shell_volume_cm3,neutron_source)
    call interp_source_per_gamma(num_gam_p,hadron_energy,qpip,pion_energy,hadronic_pion_charged_mass_gev, &
                                 shell_volume_cm3,pip_source)
    call interp_source_per_gamma(num_gam_p,hadron_energy,qpim,pion_energy,hadronic_pion_charged_mass_gev, &
                                 shell_volume_cm3,pim_source)
    call interp_source_per_gamma(num_gam_p,hadron_energy,mumnl,muon_energy,hadronic_muon_mass_gev, &
                                 shell_volume_cm3,muml_source)
    call interp_source_per_gamma(num_gam_p,hadron_energy,mumnr,muon_energy,hadronic_muon_mass_gev, &
                                 shell_volume_cm3,mumr_source)
    call interp_source_per_gamma(num_gam_p,hadron_energy,mupl,muon_energy,hadronic_muon_mass_gev, &
                                 shell_volume_cm3,mupl_source)
    call interp_source_per_gamma(num_gam_p,hadron_energy,mupr,muon_energy,hadronic_muon_mass_gev, &
                                 shell_volume_cm3,mupr_source)
end subroutine map_pgamma_secondary_sources

subroutine apply_neutron_pgamma_loss(num_gam_p,gam_p,dt_s,hadron_energy,nloss,neutron_next)
    implicit none
    integer, intent(in) :: num_gam_p
    real(8), intent(in) :: gam_p(num_gam_p),dt_s,hadron_energy(num_gam_p),nloss(num_gam_p)
    real(8), intent(inout) :: neutron_next(num_gam_p)
    real(8) :: neutron_energy(num_gam_p),neutron_loss_interp(num_gam_p)
    integer :: i

    neutron_energy=gam_p*hadronic_neutron_mass_gev
    call interp_positive_loglog(num_gam_p,hadron_energy,nloss,neutron_energy,neutron_loss_interp)
    do i=1,num_gam_p
        if (one-dt_s*neutron_loss_interp(i) < zero) error stop "hadronic p-gamma neutron survival became negative."
        neutron_next(i)=neutron_next(i)*(one-dt_s*neutron_loss_interp(i))
    end do
end subroutine apply_neutron_pgamma_loss

subroutine interp_source_per_gamma(num_src,energy_src,source_src,energy_dst,mass_gev,shell_volume_cm3,source_dst)
    implicit none
    integer, intent(in) :: num_src
    real(8), intent(in) :: energy_src(num_src),source_src(num_src),energy_dst(num_src),mass_gev,shell_volume_cm3
    real(8), intent(out) :: source_dst(num_src)
    real(8) :: source_per_gev(num_src)

    call interp_positive_loglog(num_src,energy_src,source_src,energy_dst,source_per_gev)
    source_dst=shell_volume_cm3*source_per_gev*mass_gev
end subroutine interp_source_per_gamma

subroutine interp_positive_loglog(num_src,x_src,y_src,x_dst,y_dst)
    implicit none
    integer, intent(in) :: num_src
    real(8), intent(in) :: x_src(num_src),y_src(num_src),x_dst(num_src)
    real(8), intent(out) :: y_dst(num_src)
    integer :: i,j
    real(8) :: lx,frac

    y_dst=zero
    do i=1,num_src
        if (x_dst(i) < x_src(1) .or. x_dst(i) > x_src(num_src)) cycle
        do j=1,num_src-1
            if (x_dst(i) >= x_src(j) .and. x_dst(i) <= x_src(j+1)) then
                if (y_src(j) > zero .and. y_src(j+1) > zero) then
                    lx=log(x_dst(i)/x_src(j))/log(x_src(j+1)/x_src(j))
                    frac=log(y_src(j+1)/y_src(j))*lx
                    y_dst(i)=y_src(j)*exp(frac)
                end if
                exit
            end if
        end do
    end do
end subroutine interp_positive_loglog

end module hadronic_pgamma_hummer_1d
