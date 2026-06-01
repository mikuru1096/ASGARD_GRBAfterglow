! 反激波强子核心：质子注入-冷却-输运-同步辐射，复用 hadronic 核。
subroutine fs_hadronic_reverse_1d(R_Tobs,R_Gamma,R,shell_energy_inj_erg,B_field_g,V_seed, &
                                  include_proton_synch,Num_nu,Num_R,num_gam_p, &
                                  gam_p,dN_gam_p,P_had_syn,Seed_had_syn)
    use constants
    use hadronic_common
    use hadronic_transport_kernel
    use hadronic_radiation_kernel
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: include_proton_synch,Num_nu,Num_R,num_gam_p
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),shell_energy_inj_erg(Num_R)
    real(8), intent(in) :: B_field_g(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: gam_p(num_gam_p),dN_gam_p(num_gam_p,Num_R)
    real(8), intent(out) :: P_had_syn(Num_nu,Num_R),Seed_had_syn(Num_nu,Num_R)
    integer :: I_R
    real(8) :: gam_p_max_global,t_dyn_s,dt_s,gam_p_min,energy_budget_erg
    real(8) :: dN_prev(num_gam_p),dN_next(num_gam_p),Q_inj(num_gam_p)
    real(8) :: loss_ad(num_gam_p),loss_syn(num_gam_p),loss_total(num_gam_p)

    ! 估算质子最高 Lorentz 因子（取所有壳层的最大值）
    t_dyn_s=hadronic_dynamical_time(R(1),R_Gamma(1))
    gam_p_max_global=hadronic_gamma_p_max(B_field_g(1),t_dyn_s,ten)
    do I_R=2,Num_R
        t_dyn_s=hadronic_dynamical_time(R(I_R),R_Gamma(I_R))
        gam_p_max_global=max(gam_p_max_global,hadronic_gamma_p_max( &
            B_field_g(I_R),t_dyn_s,ten))
    end do
    if (gam_p_max_global <= one+1d-3) error stop "reverse hadronic gamma_p_max must exceed the injection grid minimum."
    call hadronic_build_gamma_p_grid(num_gam_p,one+1d-3,gam_p_max_global,gam_p)
    dN_prev=zero

    dN_gam_p=zero
    P_had_syn=zero
    Seed_had_syn=zero

    do I_R=1,Num_R
        dt_s=hadronic_shell_dt(R_Tobs,I_R)
        t_dyn_s=hadronic_dynamical_time(R(I_R),R_Gamma(I_R))
        if (shell_energy_inj_erg(I_R) < zero) error stop "reverse hadronic shell injection energy must be non-negative."
        energy_budget_erg=shell_energy_inj_erg(I_R)
        gam_p_min=max(gam_p(1),R_Gamma(I_R))
        call hadronic_proton_injection_powerlaw(num_gam_p,gam_p,2.2d0,energy_budget_erg, &
                                                gam_p_min,gam_p(num_gam_p),Q_inj)
        call hadronic_proton_loss_rates(num_gam_p,gam_p,B_field_g(I_R), &
                                        t_dyn_s,loss_ad,loss_syn,loss_total)
        call hadronic_advance_energy_loggamma(num_gam_p,gam_p,dN_prev,Q_inj,loss_total,dt_s,dN_next)
        dN_gam_p(:,I_R)=dN_next

        if (include_proton_synch /= 0) then
            call hadronic_get_proton_syn_state(R(I_R),B_field_g(I_R), &
                                               num_gam_p,Num_nu,gam_p,dN_next, &
                                               V_seed,P_had_syn(:,I_R),Seed_had_syn(:,I_R))
        end if

        dN_prev=dN_next
    end do
end subroutine fs_hadronic_reverse_1d
