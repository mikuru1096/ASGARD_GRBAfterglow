! 电磁对级联 Fortran 核心：调用已有对产生核，执行冷却+同步辐射。
module hadronic_pair_cascade_kernel
    use constants
    use hadronic_common, only: hadronic_electron_mass_gev, hadronic_validate_log_grid
    use hadronic_pair_production_kernel, only: hadronic_pair_production_operator
    use electron_radiation_kernel, only: get_syn_selected_state
    implicit none
    private

    public :: hadronic_cascade_step, hadronic_cascade_sequence

contains

! 单次级联步：γγ→e⁺e⁻ (已有核) → 同步冷却演化 → 同步辐射光子。
subroutine hadronic_cascade_step(num_ph,photon_energy_gev,photon_density, &
                                  num_e,electron_energy_gev,b_field_g,path_time_s, &
                                  cascade_syn_spec,photon_loss_rate,absorbed_power)
    integer, intent(in) :: num_ph,num_e
    real(8), intent(in) :: photon_energy_gev(num_ph),photon_density(num_ph)
    real(8), intent(in) :: electron_energy_gev(num_e),b_field_g,path_time_s
    real(8), intent(out) :: cascade_syn_spec(num_ph),photon_loss_rate(num_ph),absorbed_power

    integer :: max_com_factor,n_lep,np
    real(8) :: gm_e(num_e),dln_e,coeff_syn,gm_inj,gm_c
    real(8) :: pair_inj_per_species(num_e),pair_inj_total(num_e)
    real(8) :: photon_loss(num_ph),pair_steady(num_e)
    real(8) :: injected_power_local,nu_crit,x_arg,fx,temp

    n_lep = num_e
    np = num_ph
    max_com_factor = 138  ! 与 AM3 一致
    call validate_pair_cascade_inputs

    call run_pair_production_stage

    call build_pair_electron_gamma_grid
    if (b_field_g == zero) then
        cascade_syn_spec(1:np) = zero
        return
    end if

    call evolve_pair_cooling_stage

    call emit_pair_synchrotron_stage

contains

    subroutine validate_pair_cascade_inputs
        if (b_field_g < zero) error stop "pair cascade requires b_field_g >= 0."
        if (path_time_s < zero) error stop "pair cascade requires path_time_s >= 0."
        if (any(photon_density(1:np) < zero)) error stop "pair cascade requires non-negative photon density."
        call hadronic_validate_log_grid(np,photon_energy_gev,"pair cascade photon_energy_gev")
        call hadronic_validate_log_grid(n_lep,electron_energy_gev,"pair cascade electron_energy_gev",dln_e)
    end subroutine validate_pair_cascade_inputs

    subroutine run_pair_production_stage
        call hadronic_pair_production_operator(np,photon_energy_gev,photon_density, &
                                                n_lep,electron_energy_gev, &
                                                max_com_factor,photon_loss, &
                                                pair_inj_per_species, &
                                                pair_inj_total,absorbed_power, &
                                                injected_power_local)
        photon_loss_rate(1:np) = photon_loss(1:np)
    end subroutine run_pair_production_stage

    subroutine build_pair_electron_gamma_grid
        gm_e(1:n_lep) = electron_energy_gev(1:n_lep) / hadronic_electron_mass_gev
        if (gm_e(1) < one) error stop "pair cascade electron grid must start at gamma >= 1."
    end subroutine build_pair_electron_gamma_grid

    subroutine evolve_pair_cooling_stage
        integer :: i_e

        pair_steady(1:n_lep) = zero
        coeff_syn = Para_SigmaT * b_field_g**2 / (6d0*pi*Para_m_e*Para_c)

        do i_e=2,n_lep-1
            if (pair_inj_total(i_e) <= zero) cycle
            gm_inj = gm_e(i_e)
            ! 冷却后 Lorentz 因子: 1/γ_c = 1/γ_inj + coeff_syn * path_time_s
            gm_c = gm_inj / (one + coeff_syn * gm_inj * path_time_s)
            call distribute_cooled_power(n_lep,gm_e,dln_e,gm_inj,gm_c, &
                                          pair_inj_total(i_e),coeff_syn,pair_steady)
        end do
    end subroutine evolve_pair_cooling_stage

    subroutine emit_pair_synchrotron_stage
        integer :: i_e,i_ph

        cascade_syn_spec(1:np) = zero
        temp = dsqrt(3d0)*Para_e**3 / Para_m_p_e

        do i_ph=1,np
            do i_e=2,n_lep-1
                if (pair_steady(i_e) <= zero) cycle
                nu_crit = 4.2d6 * b_field_g * gm_e(i_e)*gm_e(i_e)
                x_arg = photon_energy_gev(i_ph)*Para_GeV2Hz / nu_crit
                if (x_arg > 1d2) cycle
                fx = 1.80842d0 * x_arg**(one/3d0) * dexp(-x_arg)
                cascade_syn_spec(i_ph) = cascade_syn_spec(i_ph) + &
                    temp*b_field_g*pair_steady(i_e)*fx*gm_e(i_e)*dln_e
            end do
        end do
    end subroutine emit_pair_synchrotron_stage
end subroutine hadronic_cascade_step

subroutine hadronic_cascade_sequence(num_ph,num_e,num_shell,photon_energy_gev,primary_photon_density, &
                                     electron_energy_gev,frequency_hz,radius_cm,gamma_bulk,observer_time_s, &
                                     b_field_g,num_threads,index_syn_integr,substeps_per_shell, &
                                     photon_loss_rate,tau_pair,pair_density,pair_luminosity, &
                                     pair_seed,cascade_photon_density,absorbed_power,injected_power)
    integer, intent(in) :: num_ph,num_e,num_shell,num_threads,index_syn_integr,substeps_per_shell
    real(8), intent(in) :: photon_energy_gev(num_ph),primary_photon_density(num_ph,num_shell)
    real(8), intent(in) :: electron_energy_gev(num_e),frequency_hz(num_ph)
    real(8), intent(in) :: radius_cm(num_shell),gamma_bulk(num_shell),observer_time_s(num_shell),b_field_g(num_shell)
    real(8), intent(out) :: photon_loss_rate(num_ph,num_shell),tau_pair(num_ph,num_shell)
    real(8), intent(out) :: pair_density(num_e,num_shell),pair_luminosity(num_ph,num_shell)
    real(8), intent(out) :: pair_seed(num_ph,num_shell),cascade_photon_density(num_ph,num_shell)
    real(8), intent(out) :: absorbed_power(num_shell),injected_power(num_shell)

    integer :: i_shell,i_sub,max_com_factor
    real(8) :: gm_e(num_e),gamma_edge(num_e+1),pair_prev(num_e),pair_current(num_e)
    real(8) :: photon_density(num_ph),source_rate(num_ph),loss_rate(num_ph),q_syn_density(num_ph)
    real(8) :: pair_inj_per_species(num_e),pair_inj_total(num_e),pair_source(num_e),pair_loss(num_e)
    real(8) :: pair_next(num_e),photon_next(num_ph),dt_shell,dt_sub,t_escape,shell_volume
    real(8) :: absorbed_local,injected_local

    max_com_factor=138
    call validate_sequence_inputs
    call build_sequence_grids

    photon_loss_rate=zero; tau_pair=zero; pair_density=zero
    pair_luminosity=zero; pair_seed=zero; cascade_photon_density=zero
    absorbed_power=zero; injected_power=zero; pair_prev=zero

    do i_shell=1,num_shell
        call shell_geometry(i_shell,dt_shell,t_escape,shell_volume)
        dt_sub=dt_shell/dble(substeps_per_shell)
        photon_density(1:num_ph)=primary_photon_density(1:num_ph,i_shell)
        pair_current(1:num_e)=pair_prev(1:num_e)

        do i_sub=1,substeps_per_shell
            call hadronic_pair_production_operator(num_ph,photon_energy_gev,photon_density, &
                                                   num_e,electron_energy_gev,max_com_factor,loss_rate, &
                                                   pair_inj_per_species,pair_inj_total, &
                                                   absorbed_local,injected_local)
            pair_source(1:num_e)=shell_volume*pair_inj_total(1:num_e)*hadronic_electron_mass_gev*dt_sub
            call electron_loss_rates(num_e,gm_e,b_field_g(i_shell),dynamical_time(i_shell),pair_loss)
            call advance_energy_loggamma(num_e,gm_e,gamma_edge,pair_current,pair_source,pair_loss,dt_sub,pair_next)
            pair_current(1:num_e)=pair_next(1:num_e)

            call pair_synchrotron_state(i_shell,pair_current,pair_luminosity(:,i_shell),pair_seed(:,i_shell))
            q_syn_density(1:num_ph)=pair_seed(1:num_ph,i_shell)/Para_h_GeV/t_escape
            source_rate(1:num_ph)=primary_photon_density(1:num_ph,i_shell)/t_escape+q_syn_density(1:num_ph)
            call advance_photon_density(num_ph,photon_density,source_rate,loss_rate+one/t_escape,dt_sub,photon_next)
            photon_density(1:num_ph)=photon_next(1:num_ph)
            absorbed_power(i_shell)=absorbed_local
            injected_power(i_shell)=injected_local
        end do

        call hadronic_pair_production_operator(num_ph,photon_energy_gev,photon_density, &
                                               num_e,electron_energy_gev,max_com_factor,loss_rate, &
                                               pair_inj_per_species,pair_inj_total,absorbed_local,injected_local)
        call pair_synchrotron_state(i_shell,pair_current,pair_luminosity(:,i_shell),pair_seed(:,i_shell))
        photon_loss_rate(1:num_ph,i_shell)=loss_rate(1:num_ph)
        tau_pair(1:num_ph,i_shell)=loss_rate(1:num_ph)*t_escape
        pair_density(1:num_e,i_shell)=pair_current(1:num_e)
        cascade_photon_density(1:num_ph,i_shell)=pair_seed(1:num_ph,i_shell)/Para_h_GeV
        pair_prev(1:num_e)=pair_current(1:num_e)
    end do

contains

    subroutine validate_sequence_inputs
        if (num_ph < 2 .or. num_e < 2 .or. num_shell < 1) error stop "pair cascade sequence grid is too small."
        if (substeps_per_shell < 1) error stop "pair cascade sequence requires positive substeps."
        if (any(primary_photon_density < zero)) error stop "pair cascade sequence requires non-negative photons."
        if (any(gamma_bulk < one)) error stop "pair cascade sequence requires gamma_bulk >= 1."
        if (any(radius_cm <= zero) .or. any(b_field_g < zero)) error stop "pair cascade sequence received bad shell state."
        call hadronic_validate_log_grid(num_ph,photon_energy_gev,"pair cascade photon_energy_gev")
        call hadronic_validate_log_grid(num_e,electron_energy_gev,"pair cascade electron_energy_gev")
        call hadronic_validate_log_grid(num_ph,frequency_hz,"pair cascade frequency_hz")
    end subroutine validate_sequence_inputs

    subroutine build_sequence_grids
        integer :: i_e

        gm_e(1:num_e)=electron_energy_gev(1:num_e)/hadronic_electron_mass_gev
        if (gm_e(1) < one) error stop "pair cascade sequence electron gamma grid must start at gamma >= 1."
        gamma_edge(1)=gm_e(1)*dsqrt(gm_e(1)/gm_e(2))
        do i_e=2,num_e
            gamma_edge(i_e)=dsqrt(gm_e(i_e-1)*gm_e(i_e))
        end do
        gamma_edge(num_e+1)=gm_e(num_e)*dsqrt(gm_e(num_e)/gm_e(num_e-1))
    end subroutine build_sequence_grids

    subroutine shell_geometry(i_shell_in,dt_s,t_escape_s,volume_cm3)
        integer, intent(in) :: i_shell_in
        real(8), intent(out) :: dt_s,t_escape_s,volume_cm3

        if (i_shell_in == 1) then
            dt_s=observer_time_s(1)
            volume_cm3=(4d0/3d0)*pi*radius_cm(1)**3
        else
            dt_s=observer_time_s(i_shell_in)-observer_time_s(i_shell_in-1)
            volume_cm3=(4d0/3d0)*pi*(radius_cm(i_shell_in)**3-radius_cm(i_shell_in-1)**3)
        end if
        if (dt_s <= zero) error stop "pair cascade sequence times must be strictly increasing."
        t_escape_s=radius_cm(i_shell_in)/(12d0*gamma_bulk(i_shell_in)*Para_c)
    end subroutine shell_geometry

    real(8) function dynamical_time(i_shell_in)
        integer, intent(in) :: i_shell_in
        dynamical_time=radius_cm(i_shell_in)/(gamma_bulk(i_shell_in)*Para_c)
    end function dynamical_time

    subroutine pair_synchrotron_state(i_shell_in,pair_state,p_syn,seed_syn)
        integer, intent(in) :: i_shell_in
        real(8), intent(in) :: pair_state(num_e)
        real(8), intent(out) :: p_syn(num_ph),seed_syn(num_ph)
        real(8) :: p_emit_tmp(num_ph),tau_syn_tmp(num_ph)

        if (b_field_g(i_shell_in) == zero) then
            p_syn=zero; seed_syn=zero
            return
        end if
        call get_syn_selected_state(index_syn_integr,radius_cm(i_shell_in),b_field_g(i_shell_in), &
                                    num_e,num_ph,num_threads,gm_e,pair_state,frequency_hz,p_emit_tmp,p_syn, &
                                    seed_syn,tau_syn_tmp)
    end subroutine pair_synchrotron_state
end subroutine hadronic_cascade_sequence

! ------------------------------------------------------------
! 将注入对按冷却演化分配: N(γ) dγ = (注入数) * (冷却时间) / dln_γ
! 同步主导: dγ/dt ∝ γ², 稳态 N(γ) ∝ γ⁻² 在 γ_c < γ < γ_inj
subroutine distribute_cooled_power(n_lep,gm_e,dln_e,gm_inj,gm_c, &
                                    inj_total,coeff_syn,pair_steady)
    integer, intent(in) :: n_lep
    real(8), intent(in) :: gm_e(n_lep),dln_e,gm_inj,gm_c,inj_total,coeff_syn
    real(8), intent(inout) :: pair_steady(n_lep)
    integer :: k
    real(8) :: weight

    do k=1,n_lep
        if (gm_e(k) > gm_c .and. gm_e(k) <= gm_inj) then
            weight = one / (coeff_syn * gm_e(k)*gm_e(k))  ! t_cool(γ)
            pair_steady(k) = pair_steady(k) + inj_total * weight * dln_e
        end if
    end do
    ! 不需要额外归一化，weight 已经直接是 dN/dγ 的贡献

end subroutine distribute_cooled_power

subroutine electron_loss_rates(num_e,gm_e,b_field_g,t_dyn_s,loss_rate)
    integer, intent(in) :: num_e
    real(8), intent(in) :: gm_e(num_e),b_field_g,t_dyn_s
    real(8), intent(out) :: loss_rate(num_e)
    real(8) :: coeff_syn

    coeff_syn=Para_SigmaT*b_field_g*b_field_g/(6d0*pi*Para_m_e*Para_c)
    loss_rate(1:num_e)=coeff_syn*gm_e(1:num_e)*gm_e(1:num_e)+gm_e(1:num_e)/t_dyn_s
end subroutine electron_loss_rates

subroutine advance_energy_loggamma(num_e,gm_e,gamma_edge,density_prev,source_content,loss_rate,dt_s,density_next)
    integer, intent(in) :: num_e
    real(8), intent(in) :: gm_e(num_e),gamma_edge(num_e+1),density_prev(num_e),source_content(num_e)
    real(8), intent(in) :: loss_rate(num_e),dt_s
    real(8), intent(out) :: density_next(num_e)
    integer :: i_e,target
    real(8) :: content(num_e),next_content(num_e),cooled_gamma

    content(1:num_e)=density_prev(1:num_e)*(gamma_edge(2:num_e+1)-gamma_edge(1:num_e))
    next_content=zero
    do i_e=1,num_e
        cooled_gamma=gm_e(i_e)-loss_rate(i_e)*dt_s
        target=gamma_bin_index(num_e,gamma_edge,cooled_gamma)
        if (target >= 1 .and. target <= num_e) next_content(target)=next_content(target)+content(i_e)
    end do
    density_next(1:num_e)=next_content(1:num_e)/(gamma_edge(2:num_e+1)-gamma_edge(1:num_e))+source_content(1:num_e)
end subroutine advance_energy_loggamma

integer function gamma_bin_index(num_e,gamma_edge,gamma_value)
    integer, intent(in) :: num_e
    real(8), intent(in) :: gamma_edge(num_e+1),gamma_value
    integer :: i_e

    gamma_bin_index=0
    do i_e=1,num_e
        if (gamma_value >= gamma_edge(i_e) .and. gamma_value < gamma_edge(i_e+1)) then
            gamma_bin_index=i_e
            return
        end if
    end do
end function gamma_bin_index

subroutine advance_photon_density(num_ph,density_prev,source_rate,loss_rate,dt_s,density_next)
    integer, intent(in) :: num_ph
    real(8), intent(in) :: density_prev(num_ph),source_rate(num_ph),loss_rate(num_ph),dt_s
    real(8), intent(out) :: density_next(num_ph)
    integer :: i_ph
    real(8) :: attenuation

    do i_ph=1,num_ph
        if (loss_rate(i_ph) > zero) then
            attenuation=dexp(-loss_rate(i_ph)*dt_s)
            density_next(i_ph)=density_prev(i_ph)*attenuation+source_rate(i_ph)*(one-attenuation)/loss_rate(i_ph)
        else
            density_next(i_ph)=density_prev(i_ph)+source_rate(i_ph)*dt_s
        end if
    end do
end subroutine advance_photon_density

end module hadronic_pair_cascade_kernel
