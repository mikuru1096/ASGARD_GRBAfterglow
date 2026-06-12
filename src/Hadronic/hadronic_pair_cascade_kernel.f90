! 电磁对级联 Fortran 核心：调用已有对产生核，执行冷却+同步辐射。
module hadronic_pair_cascade_kernel
    use constants
    use hadronic_common, only: hadronic_electron_mass_gev, hadronic_validate_log_grid
    use hadronic_pair_production_kernel, only: hadronic_pair_production_operator
    implicit none
    private

    public :: hadronic_cascade_step

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

end module hadronic_pair_cascade_kernel
