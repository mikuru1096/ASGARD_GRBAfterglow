! Calculate the electron distributions of forward shock.
! Internal conserved variable:
!   U(log10(gamma_e), log10(chi), R) = dN / (dlog10(gamma_e) dlog10(chi))
! Coordinate definition:
!   chi = 1 + 8 * Gamma_sh^2 * x / R
!****************************************************************************************
subroutine fs_electron_fullhide_2d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                   Num_chi,index_Y,index_syn_intger,n_threads, &
                                   gam_e,dN_gam_e,dN_gam_e_total,P_syn,Seed_syn,V_m,V_c,V_a)
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,Num_chi,index_Y,index_syn_intger,n_threads
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e, Num_chi, Num_R)
    real(8), intent(out) :: dN_gam_e_total(Num_gam_e, Num_R)
    real(8), intent(out) :: gam_e(Num_gam_e)
    real(8), intent(out) :: P_syn(Num_nu, Num_R)
    real(8), intent(out) :: Seed_syn(Num_nu, Num_R)
    real(8), intent(out) :: V_m(Num_R), V_c(Num_R), V_a(Num_R)
    call fs_electron_transport_2d_core(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                       Num_chi,index_Y,index_syn_intger,n_threads, &
                                       gam_e,dN_gam_e,dN_gam_e_total,P_syn,Seed_syn,V_m,V_c,V_a, &
                                       .false., 'fullhide_2d')
end subroutine fs_electron_fullhide_2d

subroutine fs_electron_transport_2d_core(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                         Num_chi,index_Y,index_syn_intger,n_threads, &
                                         gam_e,dN_gam_e,dN_gam_e_total,P_syn,Seed_syn,V_m,V_c,V_a, &
                                         use_charint_transport, profile_tag)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    use electron_common, only: electron_initial_density, electron_initialize_spectrum, &
                               electron_initial_profile_exp_cutoff, electron_initial_grid_gamma, &
                               electron_gamma_m_exact, electron_injection_prefactor, &
                               electron_gamma_c_from_loss_mean, electron_source_bounds
    use electron_injection_profiles, only: electron_build_source_term_exp_cutoff
    use electron_forward_kernel, only: prepare_forward_cooling_aux_batch, assemble_forward_cooling_split, &
                                       assemble_forward_cooling_split_batch
    use electron_radiation_kernel, only: get_syn_state, get_nu_a_2d_cell_path, reduce_syn_shell_from_chi, &
                                         build_reduced_log_grid, project_syn_state_logbands
    use electron_seed_history_kernel, only: integrate_downstream_proper_time, accumulate_comoving_history_fields
    use radiation_common, only: radiation_pair_tau_headon_segment
    use electron_transport_2d_kernel, only: compute_log_chi_geometry, get_shock_transport_state, &
                                             compute_downstream_comoving_grid, bm_beta2_lab, bm_beta2_shock
    use electron_transport_2d_kernel, only: advance_eta_logchi_implicit, advance_eta_logchi_advection_charint, &
                                             advance_eta_logchi_diffusion_implicit
    use electron_transport_2d_kernel, only: advance_energy_loggamma_chi, advance_energy_loggamma_chi_charint
    implicit real(8)(a-h,o-z)

    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,Num_chi,index_Y,index_syn_intger,n_threads
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    logical, intent(in) :: use_charint_transport
    character(len=*), intent(in) :: profile_tag
    real(8), intent(out) :: dN_gam_e(Num_gam_e, Num_chi, Num_R)
    real(8), intent(out) :: dN_gam_e_total(Num_gam_e, Num_R)
    real(8), intent(out) :: gam_e(Num_gam_e)
    real(8), intent(out) :: P_syn(Num_nu, Num_R)
    real(8), intent(out) :: Seed_syn(Num_nu, Num_R)
    real(8), intent(out) :: V_m(Num_R), V_c(Num_R), V_a(Num_R)

    real(8), allocatable :: dEl(:), dEL_mean(:), kappa2_arr(:)
    real(8), allocatable :: dN_init(:), dN_init_log(:), U_log(:,:), source_eta1(:)
    real(8), allocatable :: eta_grid(:), chi_grid(:), eta_face(:), chi_face(:)
    real(8), allocatable :: a_arr(:), dln_a_dR_arr(:), dF1(:), shell_population(:), chi_population(:), dEL_mean_shell(:)
    real(8), allocatable :: V_m_chi(:), V_c_chi(:), V_a_chi(:), chi_weight(:), Epsilon_b_chi(:), DB_chi(:), t_decay_chi(:)
    real(8), allocatable :: proper_time_arr(:), x_face_hist(:,:), x_comov_face_hist(:,:), &
                             x_comov_hist(:,:), dx_comov_hist(:,:), dN_cell(:)
    real(8), allocatable :: beta_hist(:,:)
    real(8), allocatable :: P_local(:,:), kappa2_chi(:,:)
    real(8), allocatable :: V_cool(:)
    real(8), allocatable :: P_hist(:,:,:), Seed_hist(:,:,:), Tau_hist(:,:,:)
    real(8), allocatable :: P_hist_cool(:,:,:), Seed_hist_cool(:,:,:), Tau_hist_cool(:,:,:), &
                             Tau_pair_hist_cool(:,:,:), Tau_prop_hist_cool(:,:,:)
    real(8), allocatable :: P_eff_cool_chi(:,:), Seed_eff_cool_chi(:,:)
    real(8), allocatable :: cooling_aux_chi(:,:), dEl_chi(:,:), dEL_mean_chi(:,:)

    real(8) :: temp, chi_max_global, deta, d_x_E, ln10
    real(8) :: R_loc, R_Gamma_loc, dNe, Para_N_e_ini, DB, DB_min
    real(8) :: Epsilon_b_floor, magnetic_decay_alpha_t, magnetic_decay_t0_s
    real(8) :: Gam_e_max, Gam_e_max_max, Gam_e_m, Gam_e_c, Gam_e_c_diag, temp_gam
    real(8) :: Gam_e_max_cell, Gam_e_m_cell, Gam_e_c_cell, beta_2_sh_loc
    real(8) :: beta_sh, beta_2, beta_2_sh
    real(8) :: dDR_try, dDR_xi, dDR_eta, dDD, dDR, max_xi_coeff, max_eta_coeff
    real(8) :: frac_sub, R_sub, Gamma_sh_sub, a_sub, dln_a_dR_sub
    real(8) :: Q, Gam_e_m_p
    real(8) :: support_floor, shell_peak
    real(8) :: x_l, x_c, x_r, y_l, y_c, y_r, x_peak, denom_peak
    real(8) :: t_start, t_stop
    real(8) :: t_hist_accum, t_syn_state, t_prepare_aux, t_cooling, t_eta, t_xi
    integer :: I_tobs, I_chi, I_gam_e, L1, L, src_lo, src_hi, active_hi, active_chi_hi, Num_nu_cool
    integer :: I_nu
    integer :: total_substeps, max_shell_substeps, shell_cooling_calls, substep_cooling_calls
    integer :: prepare_aux_calls, history_calls, syn_state_calls, eta_calls, xi_calls
    integer :: env_len, env_status
    logical :: profile_enabled, magnetic_decay_active
    character(len=32) :: profile_env

    allocate(dEl(Num_gam_e), dEL_mean(Num_gam_e-1), dEL_mean_shell(Num_gam_e-1), kappa2_arr(Num_gam_e), &
             dN_init(Num_gam_e), dN_init_log(Num_gam_e), dF1(Num_gam_e), shell_population(Num_gam_e), chi_population(Num_chi), &
             V_m_chi(Num_chi), V_c_chi(Num_chi), V_a_chi(Num_chi), chi_weight(Num_chi), &
             Epsilon_b_chi(Num_chi), DB_chi(Num_chi), t_decay_chi(Num_chi), &
             source_eta1(Num_gam_e), U_log(Num_gam_e, Num_chi), eta_grid(Num_chi), &
             chi_grid(Num_chi), eta_face(0:Num_chi), chi_face(0:Num_chi), &
             a_arr(Num_R), dln_a_dR_arr(Num_R), proper_time_arr(Num_R), x_face_hist(0:Num_chi, Num_R), &
             x_comov_face_hist(0:Num_chi, Num_R), x_comov_hist(Num_chi, Num_R), dx_comov_hist(Num_chi, Num_R), &
             beta_hist(Num_chi, Num_R), dN_cell(Num_gam_e), P_local(Num_nu, Num_chi), kappa2_chi(Num_gam_e, Num_chi), &
             cooling_aux_chi(Num_gam_e, Num_chi), dEl_chi(Num_gam_e, Num_chi), dEL_mean_chi(Num_gam_e-1, Num_chi))

    Eta_0            = Boundary(1)
    R_ini            = Boundary(4)
    Epsilon_e        = Boundary(5)
    Epsilon_b        = Boundary(6)
    p                = Boundary(7)
    z                = Boundary(8)
    dNe_ISM          = Boundary(11)
    A_star           = Boundary(12)
    E_iso            = Boundary(14)
    T_log10_duration = Boundary(15)
    f_e              = Boundary(16)
    R_tr             = Boundary(21)
    f_jump           = Boundary(22)
    f_wide           = Boundary(23)
    Epsilon_b_floor  = Epsilon_b
    magnetic_decay_alpha_t = zero
    magnetic_decay_t0_s = one
    if (n >= 27) then
        Epsilon_b_floor = Boundary(24)
        magnetic_decay_alpha_t = Boundary(25)
        magnetic_decay_t0_s = Boundary(26)
    end if
    R0               = Boundary(n)
    magnetic_decay_active = magnetic_decay_alpha_t < zero .and. magnetic_decay_t0_s > zero .and. &
                            Epsilon_b_floor > zero .and. Epsilon_b_floor < Epsilon_b

    ln10 = dlog(ten)
    profile_enabled = .false.
    profile_env = ''
    call get_environment_variable('ASGARD_PROFILE_2D', profile_env, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
        if (profile_env(1:1) /= '0') profile_enabled = .true.
    end if
    t_hist_accum = zero
    t_syn_state = zero
    t_prepare_aux = zero
    t_cooling = zero
    t_eta = zero
    t_xi = zero
    total_substeps = 0
    max_shell_substeps = 0
    shell_cooling_calls = 0
    substep_cooling_calls = 0
    prepare_aux_calls = 0
    history_calls = 0
    syn_state_calls = 0
    eta_calls = 0
    xi_calls = 0
    Num_nu_cool = min(6, Num_nu)
    allocate(V_cool(Num_nu_cool), P_hist(Num_nu, Num_chi, Num_R), &
             Seed_hist(Num_nu, Num_chi, Num_R), Tau_hist(Num_nu, Num_chi, Num_R), &
             P_hist_cool(Num_nu_cool, Num_chi, Num_R), Seed_hist_cool(Num_nu_cool, Num_chi, Num_R), &
             Tau_hist_cool(Num_nu_cool, Num_chi, Num_R), &
             Tau_pair_hist_cool(Num_nu_cool, Num_chi, Num_R), Tau_prop_hist_cool(Num_nu_cool, Num_chi, Num_R), &
             P_eff_cool_chi(Num_nu_cool, Num_chi), Seed_eff_cool_chi(Num_nu_cool, Num_chi))
    call build_reduced_log_grid(Num_nu,V_seed,Num_nu_cool,V_cool)

    P_syn          = zero
    Seed_syn       = zero
    V_m            = zero
    V_c            = zero
    V_a            = zero
    dN_gam_e       = zero
    dN_gam_e_total = zero
    U_log          = zero
    dF1            = zero

    call compute_log_chi_geometry(Num_R, Num_chi, R, R_Gamma, a_arr, dln_a_dR_arr, &
                                  chi_max_global, deta, eta_face, eta_grid, chi_face, chi_grid)
    call integrate_downstream_proper_time(Num_R,R,R_Gamma,proper_time_arr)

    call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)

    DB            = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    DB_min        = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max     = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    Gam_e_max_max = 3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam      = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-one)
    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c       = 7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)

    call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                      electron_initial_profile_exp_cutoff,electron_initial_grid_gamma,gam_e,dN_init)
    d_x_E = dlog10(gam_e(2)/gam_e(1))

    dN_init_log = dN_init * gam_e * ln10

    U_log(:,1) = dN_init_log / deta

    do I_chi = 1, Num_chi
        dN_gam_e(:, I_chi, 1) = U_log(:, I_chi) / (ln10*ln10*gam_e*chi_grid(I_chi))
    end do
    dN_gam_e_total(:,1) = zero
    do I_chi = 1, Num_chi
        dN_gam_e_total(:,1) = dN_gam_e_total(:,1) + U_log(:, I_chi) * deta / (gam_e*ln10)
    end do
    call compute_downstream_comoving_grid(Num_chi,R(1),R_Gamma(1),chi_face,chi_grid, &
                                          x_face_hist(:,1),x_comov_face_hist(:,1),x_comov_hist(:,1),dx_comov_hist(:,1))
    do I_chi = 1, Num_chi
        beta_hist(I_chi,1) = bm_beta2_lab(R_Gamma(1),chi_grid(I_chi))
        beta_2_sh_loc = bm_beta2_shock(R_Gamma(1),chi_grid(I_chi))
        t_decay_chi(I_chi) = x_comov_hist(I_chi,1)/max(beta_2_sh_loc*para_c,tiny(one))
        if (magnetic_decay_active) then
            Epsilon_b_chi(I_chi) = Epsilon_b_floor + (Epsilon_b-Epsilon_b_floor) * &
                                   (one+t_decay_chi(I_chi)/magnetic_decay_t0_s)**magnetic_decay_alpha_t
        else
            Epsilon_b_chi(I_chi) = Epsilon_b
        end if
        DB_chi(I_chi) = 0.39d0*dsqrt(Epsilon_b_chi(I_chi)*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    end do
    do I_chi = 1, Num_chi
        dN_cell = dN_gam_e(:,I_chi,1)
        if (profile_enabled) call cpu_time(t_start)
        call get_syn_state(R(1),DB_chi(I_chi),Num_gam_e,Num_nu,n_threads,gam_e,dN_cell,V_seed, &
                       P_local(:,I_chi),P_hist(:,I_chi,1),Seed_hist(:,I_chi,1),Tau_hist(:,I_chi,1))
        call project_syn_state_logbands(Num_nu,V_seed,P_hist(:,I_chi,1),Seed_hist(:,I_chi,1),Tau_hist(:,I_chi,1), &
                                        Num_nu_cool,V_cool,P_hist_cool(:,I_chi,1), &
                                        Seed_hist_cool(:,I_chi,1),Tau_hist_cool(:,I_chi,1))
        Tau_pair_hist_cool(:,I_chi,1) = zero
        call radiation_pair_tau_headon_segment(V_cool,Num_nu_cool,Seed_hist_cool(:,I_chi,1),dx_comov_hist(I_chi,1), &
                                               Tau_pair_hist_cool(:,I_chi,1))
        Tau_prop_hist_cool(:,I_chi,1) = Tau_hist_cool(:,I_chi,1) + Tau_pair_hist_cool(:,I_chi,1)
        if (profile_enabled) then
            call cpu_time(t_stop)
            t_syn_state = t_syn_state + (t_stop-t_start)
        end if
        syn_state_calls = syn_state_calls + 1
    end do

    do I_tobs = 2, Num_R
        R_loc       = R(I_tobs-1)
        R_Gamma_loc = (R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two

        call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)

        DB        = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam  = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_c   = 7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)

        call get_shock_transport_state(R_Gamma_loc, beta_sh, beta_2, beta_2_sh)

        call reduce_syn_shell_from_chi(Num_nu,Num_chi,deta,chi_grid,P_hist(:,:,I_tobs-1),Seed_hist(:,:,I_tobs-1), &
                                       P_syn(:,I_tobs),Seed_syn(:,I_tobs))

        if (profile_enabled) call cpu_time(t_start)
        call accumulate_comoving_history_fields(I_tobs-1,Num_R,Num_chi,Num_nu_cool,proper_time_arr,V_cool, &
                                                x_comov_face_hist,x_comov_hist,dx_comov_hist,beta_hist, &
                                                Tau_prop_hist_cool,P_hist_cool,Seed_hist_cool, &
                                                P_eff_cool_chi,Seed_eff_cool_chi)
        history_calls = history_calls + 1
        if (profile_enabled) then
            call cpu_time(t_stop)
            t_hist_accum = t_hist_accum + (t_stop-t_start)
        end if
        do I_chi = 1, Num_chi
            beta_2_sh_loc = bm_beta2_shock(R_Gamma_loc,chi_grid(I_chi))
            t_decay_chi(I_chi) = x_comov_hist(I_chi,I_tobs-1)/max(beta_2_sh_loc*para_c,tiny(one))
            if (magnetic_decay_active) then
                Epsilon_b_chi(I_chi) = Epsilon_b_floor + (Epsilon_b-Epsilon_b_floor) * &
                                       (one+t_decay_chi(I_chi)/magnetic_decay_t0_s)**magnetic_decay_alpha_t
            else
                Epsilon_b_chi(I_chi) = Epsilon_b
            end if
            DB_chi(I_chi) = 0.39d0*dsqrt(Epsilon_b_chi(I_chi)*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        end do

        if (profile_enabled) call cpu_time(t_start)
        call prepare_forward_cooling_aux_batch(index_Y,Num_gam_e,Num_nu_cool,Num_chi,n_threads,gam_e,V_cool, &
                                               P_eff_cool_chi,Seed_eff_cool_chi,cooling_aux_chi)
        if (profile_enabled) then
            call cpu_time(t_stop)
            t_prepare_aux = t_prepare_aux + (t_stop-t_start)
        end if
        prepare_aux_calls = prepare_aux_calls + 1
        if (profile_enabled) call cpu_time(t_start)
        if (magnetic_decay_active) then
            do I_chi = 1, Num_chi
                Gam_e_max_cell = 3d0*Para_m_energy/dsqrt(8d0*DB_chi(I_chi)*Para_e**3)
                temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
                call electron_gamma_m_exact(p,temp_gam,Gam_e_max_cell,Gam_e_m_cell)
                Gam_e_c_cell = 7.7d8*(one+z)/R_Gamma_loc/DB_chi(I_chi)**2/max(R_Tobs(I_tobs),tiny(one))
                call assemble_forward_cooling_split(index_Y,Epsilon_e,Epsilon_b_chi(I_chi),p,DB_chi(I_chi), &
                                                    Gam_e_m_cell,Gam_e_c_cell, &
                                                    Gam_e_max_cell,R_loc,R_Gamma_loc,beta_sh,dNe,Num_gam_e, &
                                                    Num_nu_cool,n_threads,gam_e,V_cool, &
                                                    P_eff_cool_chi(:,I_chi),Seed_eff_cool_chi(:,I_chi),Seed_eff_cool_chi(:,I_chi), &
                                                    cooling_aux_chi(:,I_chi),dEl_chi(:,I_chi))
                dEL_mean_chi(:,I_chi)=(dEl_chi(2:Num_gam_e,I_chi)+dEl_chi(1:Num_gam_e-1,I_chi))/two/dlog(ten)
                kappa2_chi(:,I_chi) = gam_e*Para_m_energy*para_c/(3d0*Para_e*DB_chi(I_chi))
            end do
        else
            call assemble_forward_cooling_split_batch(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c, &
                                                      Gam_e_max,R_loc,R_Gamma_loc, &
                                                      beta_sh,dNe,Num_gam_e,Num_nu_cool,Num_chi,n_threads,gam_e,V_cool, &
                                                      P_eff_cool_chi,Seed_eff_cool_chi,Seed_eff_cool_chi,cooling_aux_chi,dEl_chi)
            do I_chi = 1, Num_chi
                dEL_mean_chi(:,I_chi)=(dEl_chi(2:Num_gam_e,I_chi)+dEl_chi(1:Num_gam_e-1,I_chi))/two/dlog(ten)
            end do
            do I_chi = 1, Num_chi
                kappa2_chi(:,I_chi) = gam_e*Para_m_energy*para_c/(3d0*Para_e*DB)
            end do
        end if
        if (profile_enabled) then
            call cpu_time(t_stop)
            t_cooling = t_cooling + (t_stop-t_start)
        end if
        shell_cooling_calls = shell_cooling_calls + 1
        call get_nu_a_2d_cell_path(Num_nu,Num_chi,V_seed,Tau_hist(:,:,I_tobs-1),V_a_chi)
        temp = zero
        Q = zero
        do I_chi = 1, Num_chi
            chi_weight(I_chi) = max(sum(U_log(:,I_chi)), tiny(one))
            I_nu = maxloc(P_hist(:,I_chi,I_tobs-1), dim=1)
            V_m_chi(I_chi) = max(V_seed(I_nu), tiny(one))
            if (I_nu > 1 .and. I_nu < Num_nu) then
                if (P_hist(I_nu-1,I_chi,I_tobs-1) > zero .and. P_hist(I_nu,I_chi,I_tobs-1) > zero .and. &
                    P_hist(I_nu+1,I_chi,I_tobs-1) > zero) then
                    x_l = dlog(V_seed(I_nu-1))
                    x_c = dlog(V_seed(I_nu))
                    x_r = dlog(V_seed(I_nu+1))
                    y_l = dlog(P_hist(I_nu-1,I_chi,I_tobs-1))
                    y_c = dlog(P_hist(I_nu,I_chi,I_tobs-1))
                    y_r = dlog(P_hist(I_nu+1,I_chi,I_tobs-1))
                    denom_peak = y_l - two*y_c + y_r
                    if (dabs(denom_peak) > tiny(one)) then
                        x_peak = x_c + 0.5d0*(y_l-y_r)*(x_c-x_l)/denom_peak
                        x_peak = min(max(x_peak, x_l), x_r)
                        V_m_chi(I_chi) = dexp(x_peak)
                    end if
                end if
            end if
            dEL_mean_shell = dEL_mean_chi(:,I_chi)
            call electron_gamma_c_from_loss_mean(Num_gam_e,gam_e,dEL_mean_shell,R_loc,Gam_e_c_diag)
            V_c_chi(I_chi) = max(4.2d6*DB_chi(I_chi)*Gam_e_c_diag*Gam_e_c_diag, tiny(one))
            temp = temp + chi_weight(I_chi)*dlog(V_m_chi(I_chi))
            Q = Q + chi_weight(I_chi)
        end do
        V_m(I_tobs-1) = dexp(temp/max(Q,tiny(one)))/(R_Gamma_loc*(1d0-beta_sh)*(one+z))
        temp = zero
        do I_chi = 1, Num_chi
            temp = temp + chi_weight(I_chi)*dlog(V_c_chi(I_chi))
        end do
        V_c(I_tobs-1) = dexp(temp/max(Q,tiny(one)))/(R_Gamma_loc*(1d0-beta_sh)*(one+z))
        V_a(I_tobs-1) = V_a_chi(Num_chi)/(R_Gamma_loc*(1d0-beta_sh)*(one+z))
        call electron_source_bounds(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,src_lo,src_hi)
        call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,one,p,dF1)
        temp = maxval(dF1)
        if (temp > zero) then
            src_hi = max(2, src_lo)
            do I_gam_e = Num_gam_e, 1, -1
                if (dF1(I_gam_e) > 1d-12*temp) then
                    src_hi = max(src_hi, min(Num_gam_e, I_gam_e+1))
                    exit
                end if
            end do
        end if
        shell_population = sum(U_log, dim=2)
        chi_population = sum(U_log, dim=1)
        chi_peak = maxval(shell_population)
        shell_peak = max(chi_peak, maxval(dF1))
        support_floor = 1d-12*shell_peak
        active_hi = max(2, src_hi)
        if (shell_peak > zero) then
            do I_gam_e = Num_gam_e, 1, -1
                if (shell_population(I_gam_e) > support_floor) then
                    active_hi = max(active_hi, min(Num_gam_e, I_gam_e+1))
                    exit
                end if
            end do
        end if
        max_xi_coeff = zero
        if (active_hi > 1) then
            do I_chi = 1, Num_chi
                if (chi_peak > zero) then
                    if (chi_population(I_chi) <= 1d-10*chi_peak) cycle
                end if
                max_xi_coeff = max(max_xi_coeff, maxval(dabs(dEL_mean_chi(1:active_hi-1,I_chi) + one/R_loc/ln10)))
            end do
            if (max_xi_coeff <= zero) then
                max_xi_coeff = maxval(dabs(dEL_mean_chi(1:active_hi-1,:) + one/R_loc/ln10))
            end if
        end if

        dDR_xi = huge(one)
        if (max_xi_coeff > zero) dDR_xi = 0.4d0*d_x_E/max_xi_coeff
        active_chi_hi = Num_chi
        if (chi_peak > zero) then
            active_chi_hi = 1
            do I_chi = Num_chi, 1, -1
                if (chi_population(I_chi) > 1d-10*chi_peak) then
                    active_chi_hi = min(Num_chi, I_chi+1)
                    exit
                end if
            end do
        end if
        dDR_eta = huge(one)
        if (use_charint_transport) then
            max_eta_coeff = zero
            do I_chi = 1, max(1, active_chi_hi)
                beta_2_sh_loc = bm_beta2_shock(R_Gamma_loc,chi_face(I_chi))
                max_eta_coeff = max(max_eta_coeff, dabs((8d0*R_Gamma_loc*R_Gamma_loc/R_loc) * &
                                                        beta_2_sh_loc/(chi_face(I_chi)*beta_sh) + &
                                                        ((chi_face(I_chi)-one)/chi_face(I_chi))*dln_a_dR_arr(I_tobs-1)) / ln10)
            end do
            if (max_eta_coeff > zero) dDR_eta = 0.4d0*deta/max_eta_coeff
        end if

        dDD     = R(I_tobs)-R(I_tobs-1)
        dDR_try = min(dDD, min(dDR_xi, dDR_eta))
        L1      = max(1, ceiling(dDD/max(dDR_try, tiny(one))))
        if (use_charint_transport) L1 = min(L1, 2048)
        dDR     = dDD/dble(L1)
        total_substeps = total_substeps + L1
        max_shell_substeps = max(max_shell_substeps, L1)

        do L = 1, L1
            frac_sub = (dble(L)-0.5d0)/dble(L1)
            R_sub = R(I_tobs-1) + frac_sub*dDD
            Gamma_sh_sub = (one-frac_sub)*R_Gamma(I_tobs-1) + frac_sub*R_Gamma(I_tobs)
            a_sub = 8d0*Gamma_sh_sub*Gamma_sh_sub/R_sub
            dln_a_dR_sub = (one-frac_sub)*dln_a_dR_arr(I_tobs-1) + frac_sub*dln_a_dR_arr(I_tobs)

            call dynamics_external_density_profile(A_star,dNe_ISM,R_sub,R0,1,R_tr,f_jump,f_wide,dNe)

            DB = 0.39d0*dsqrt(Epsilon_b*dNe*(Gamma_sh_sub*(Gamma_sh_sub-one)))
            Gam_e_max = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
            temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(Gamma_sh_sub-one)
            call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
            Gam_e_c = 7.7d8*(one+z)/Gamma_sh_sub/DB**2/max(R_Tobs(I_tobs),tiny(one))
            call get_shock_transport_state(Gamma_sh_sub, beta_sh, beta_2, beta_2_sh)

            Gam_e_m_p = (one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
            call electron_injection_prefactor(R_sub,dDR,dNe,f_e,Gam_e_m_p,Q)
            call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,Q,p,dF1)
            source_eta1 = dF1/deta

            if (use_charint_transport) then
                if (profile_enabled) call cpu_time(t_start)
                call advance_eta_logchi_advection_charint(U_log, Num_gam_e, Num_chi, active_hi, deta, eta_face, chi_face, &
                                                          Gamma_sh_sub, a_sub, dln_a_dR_sub, beta_sh, source_eta1, dDR)
                call advance_eta_logchi_diffusion_implicit(U_log, Num_gam_e, Num_chi, active_hi, deta, chi_face, Gamma_sh_sub, &
                                                           a_sub, dln_a_dR_sub, beta_sh, kappa2_chi, dDR, n_threads)
                eta_calls = eta_calls + 1
                if (profile_enabled) then
                    call cpu_time(t_stop)
                    t_eta = t_eta + (t_stop-t_start)
                end if

                if (profile_enabled) call cpu_time(t_start)
                call advance_energy_loggamma_chi_charint(U_log, Num_gam_e, Num_chi, gam_e, DB_chi, dEl_chi, R_sub, &
                                                         Gamma_sh_sub, beta_sh, index_Y, dDR, active_chi_hi, n_threads)
                xi_calls = xi_calls + 1
                if (profile_enabled) then
                    call cpu_time(t_stop)
                    t_xi = t_xi + (t_stop-t_start)
                end if
            else
                if (profile_enabled) call cpu_time(t_start)
                call advance_eta_logchi_implicit(U_log, Num_gam_e, Num_chi, active_hi, deta, chi_face, Gamma_sh_sub, &
                                                 a_sub, dln_a_dR_sub, beta_sh, kappa2_chi, source_eta1, dDR, n_threads)
                eta_calls = eta_calls + 1
                if (profile_enabled) then
                    call cpu_time(t_stop)
                    t_eta = t_eta + (t_stop-t_start)
                end if

                if (profile_enabled) call cpu_time(t_start)
                call advance_energy_loggamma_chi(U_log, Num_gam_e, Num_chi, dEL_mean_chi, R_sub, d_x_E, dDR, n_threads)
                xi_calls = xi_calls + 1
                if (profile_enabled) then
                    call cpu_time(t_stop)
                    t_xi = t_xi + (t_stop-t_start)
                end if
            end if
        end do

        do I_chi = 1, Num_chi
            dN_gam_e(:, I_chi, I_tobs) = U_log(:, I_chi) / (ln10*ln10*gam_e*chi_grid(I_chi))
        end do

        dN_gam_e_total(:, I_tobs) = zero
        do I_chi = 1, Num_chi
            dN_gam_e_total(:, I_tobs) = dN_gam_e_total(:, I_tobs) + U_log(:, I_chi)*deta/(gam_e*ln10)
        end do
        call dynamics_external_density_profile(A_star,dNe_ISM,R(I_tobs),R0,1,R_tr,f_jump,f_wide,dNe)
        call compute_downstream_comoving_grid(Num_chi,R(I_tobs),R_Gamma(I_tobs),chi_face,chi_grid, &
                                              x_face_hist(:,I_tobs),x_comov_face_hist(:,I_tobs),x_comov_hist(:,I_tobs), &
                                              dx_comov_hist(:,I_tobs))
        do I_chi = 1, Num_chi
            beta_hist(I_chi,I_tobs) = bm_beta2_lab(R_Gamma(I_tobs),chi_grid(I_chi))
            beta_2_sh_loc = bm_beta2_shock(R_Gamma(I_tobs),chi_grid(I_chi))
            t_decay_chi(I_chi) = x_comov_hist(I_chi,I_tobs)/max(beta_2_sh_loc*para_c,tiny(one))
            if (magnetic_decay_active) then
                Epsilon_b_chi(I_chi) = Epsilon_b_floor + (Epsilon_b-Epsilon_b_floor) * &
                                       (one+t_decay_chi(I_chi)/magnetic_decay_t0_s)**magnetic_decay_alpha_t
            else
                Epsilon_b_chi(I_chi) = Epsilon_b
            end if
            DB_chi(I_chi) = 0.39d0*dsqrt(Epsilon_b_chi(I_chi)*dNe*(R_Gamma(I_tobs)*(R_Gamma(I_tobs)-one)))
            dN_cell = dN_gam_e(:,I_chi,I_tobs)
            if (profile_enabled) call cpu_time(t_start)
            call get_syn_state(R(I_tobs),DB_chi(I_chi),Num_gam_e,Num_nu,n_threads,gam_e,dN_cell,V_seed, &
                               P_local(:,I_chi),P_hist(:,I_chi,I_tobs),Seed_hist(:,I_chi,I_tobs),Tau_hist(:,I_chi,I_tobs))
            call project_syn_state_logbands(Num_nu,V_seed,P_hist(:,I_chi,I_tobs), &
                                            Seed_hist(:,I_chi,I_tobs),Tau_hist(:,I_chi,I_tobs), &
                                            Num_nu_cool,V_cool,P_hist_cool(:,I_chi,I_tobs), &
                                            Seed_hist_cool(:,I_chi,I_tobs),Tau_hist_cool(:,I_chi,I_tobs))
            Tau_pair_hist_cool(:,I_chi,I_tobs) = zero
            call radiation_pair_tau_headon_segment(V_cool,Num_nu_cool,Seed_hist_cool(:,I_chi,I_tobs),dx_comov_hist(I_chi,I_tobs), &
                                                   Tau_pair_hist_cool(:,I_chi,I_tobs))
            Tau_prop_hist_cool(:,I_chi,I_tobs) = Tau_hist_cool(:,I_chi,I_tobs) + Tau_pair_hist_cool(:,I_chi,I_tobs)
            if (profile_enabled) then
                call cpu_time(t_stop)
                t_syn_state = t_syn_state + (t_stop-t_start)
            end if
            syn_state_calls = syn_state_calls + 1
        end do
    end do

    R_loc = R(Num_R)
    R_Gamma_loc = R_Gamma(Num_R)
    call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
    DB = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
    Gam_e_max = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c = 7.7d8*(one+z)/R_Gamma_loc/DB**2/max(R_Tobs(Num_R),tiny(one))
    call get_shock_transport_state(R_Gamma_loc, beta_sh, beta_2, beta_2_sh)

    call reduce_syn_shell_from_chi(Num_nu,Num_chi,deta,chi_grid,P_hist(:,:,Num_R),Seed_hist(:,:,Num_R), &
                                   P_syn(:,Num_R),Seed_syn(:,Num_R))
    if (profile_enabled) call cpu_time(t_start)
    call accumulate_comoving_history_fields(Num_R,Num_R,Num_chi,Num_nu_cool,proper_time_arr,V_cool, &
                                            x_comov_face_hist,x_comov_hist,dx_comov_hist,beta_hist, &
                                            Tau_prop_hist_cool,P_hist_cool,Seed_hist_cool, &
                                            P_eff_cool_chi,Seed_eff_cool_chi)
    history_calls = history_calls + 1
    if (profile_enabled) then
        call cpu_time(t_stop)
        t_hist_accum = t_hist_accum + (t_stop-t_start)
    end if
    if (profile_enabled) call cpu_time(t_start)
    call prepare_forward_cooling_aux_batch(index_Y,Num_gam_e,Num_nu_cool,Num_chi,n_threads,gam_e,V_cool, &
                                           P_eff_cool_chi,Seed_eff_cool_chi,cooling_aux_chi)
    if (profile_enabled) then
        call cpu_time(t_stop)
        t_prepare_aux = t_prepare_aux + (t_stop-t_start)
    end if
    prepare_aux_calls = prepare_aux_calls + 1
    if (profile_enabled) call cpu_time(t_start)
    do I_chi = 1, Num_chi
        beta_2_sh_loc = bm_beta2_shock(R_Gamma_loc,chi_grid(I_chi))
        t_decay_chi(I_chi) = x_comov_hist(I_chi,Num_R)/max(beta_2_sh_loc*para_c,tiny(one))
        if (magnetic_decay_active) then
            Epsilon_b_chi(I_chi) = Epsilon_b_floor + (Epsilon_b-Epsilon_b_floor) * &
                                   (one+t_decay_chi(I_chi)/magnetic_decay_t0_s)**magnetic_decay_alpha_t
        else
            Epsilon_b_chi(I_chi) = Epsilon_b
        end if
        DB_chi(I_chi) = 0.39d0*dsqrt(Epsilon_b_chi(I_chi)*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
    end do
    if (magnetic_decay_active) then
        do I_chi = 1, Num_chi
            Gam_e_max_cell = 3d0*Para_m_energy/dsqrt(8d0*DB_chi(I_chi)*Para_e**3)
            temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
            call electron_gamma_m_exact(p,temp_gam,Gam_e_max_cell,Gam_e_m_cell)
            Gam_e_c_cell = 7.7d8*(one+z)/R_Gamma_loc/DB_chi(I_chi)**2/max(R_Tobs(Num_R),tiny(one))
            call assemble_forward_cooling_split(index_Y,Epsilon_e,Epsilon_b_chi(I_chi),p,DB_chi(I_chi),Gam_e_m_cell,Gam_e_c_cell, &
                                                Gam_e_max_cell,R_loc,R_Gamma_loc,beta_sh,dNe,Num_gam_e, &
                                                Num_nu_cool,n_threads,gam_e,V_cool, &
                                                P_eff_cool_chi(:,I_chi),Seed_eff_cool_chi(:,I_chi),Seed_eff_cool_chi(:,I_chi), &
                                                cooling_aux_chi(:,I_chi),dEl_chi(:,I_chi))
            dEL_mean_chi(:,I_chi)=(dEl_chi(2:Num_gam_e,I_chi)+dEl_chi(1:Num_gam_e-1,I_chi))/two/dlog(ten)
        end do
    else
        call assemble_forward_cooling_split_batch(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                                  beta_sh,dNe,Num_gam_e,Num_nu_cool,Num_chi,n_threads,gam_e,V_cool, &
                                                  P_eff_cool_chi,Seed_eff_cool_chi,Seed_eff_cool_chi,cooling_aux_chi,dEl_chi)
        do I_chi = 1, Num_chi
            dEL_mean_chi(:,I_chi)=(dEl_chi(2:Num_gam_e,I_chi)+dEl_chi(1:Num_gam_e-1,I_chi))/two/dlog(ten)
        end do
    end if
    if (profile_enabled) then
        call cpu_time(t_stop)
        t_cooling = t_cooling + (t_stop-t_start)
    end if
    call get_nu_a_2d_cell_path(Num_nu,Num_chi,V_seed,Tau_hist(:,:,Num_R),V_a_chi)
    temp = zero
    Q = zero
    do I_chi = 1, Num_chi
        chi_weight(I_chi) = max(sum(U_log(:,I_chi)), tiny(one))
        I_nu = maxloc(P_hist(:,I_chi,Num_R), dim=1)
        V_m_chi(I_chi) = max(V_seed(I_nu), tiny(one))
        if (I_nu > 1 .and. I_nu < Num_nu) then
            if (P_hist(I_nu-1,I_chi,Num_R) > zero .and. P_hist(I_nu,I_chi,Num_R) > zero .and. &
                P_hist(I_nu+1,I_chi,Num_R) > zero) then
                x_l = dlog(V_seed(I_nu-1))
                x_c = dlog(V_seed(I_nu))
                x_r = dlog(V_seed(I_nu+1))
                y_l = dlog(P_hist(I_nu-1,I_chi,Num_R))
                y_c = dlog(P_hist(I_nu,I_chi,Num_R))
                y_r = dlog(P_hist(I_nu+1,I_chi,Num_R))
                denom_peak = y_l - two*y_c + y_r
                if (dabs(denom_peak) > tiny(one)) then
                    x_peak = x_c + 0.5d0*(y_l-y_r)*(x_c-x_l)/denom_peak
                    x_peak = min(max(x_peak, x_l), x_r)
                    V_m_chi(I_chi) = dexp(x_peak)
                end if
            end if
        end if
        dEL_mean_shell = dEL_mean_chi(:,I_chi)
        call electron_gamma_c_from_loss_mean(Num_gam_e,gam_e,dEL_mean_shell,R_loc,Gam_e_c_diag)
        V_c_chi(I_chi) = max(4.2d6*DB_chi(I_chi)*Gam_e_c_diag*Gam_e_c_diag, tiny(one))
        temp = temp + chi_weight(I_chi)*dlog(V_m_chi(I_chi))
        Q = Q + chi_weight(I_chi)
    end do
    V_m(Num_R) = dexp(temp/max(Q,tiny(one)))/(R_Gamma_loc*(1d0-beta_sh)*(one+z))
    temp = zero
    do I_chi = 1, Num_chi
        temp = temp + chi_weight(I_chi)*dlog(V_c_chi(I_chi))
    end do
    V_c(Num_R) = dexp(temp/max(Q,tiny(one)))/(R_Gamma_loc*(1d0-beta_sh)*(one+z))
    V_a(Num_R) = V_a_chi(Num_chi)/(R_Gamma_loc*(1d0-beta_sh)*(one+z))
    if (profile_enabled) then
        print '(A,1X,F10.4)', 'PROFILE '//trim(profile_tag)//' history_s', t_hist_accum
        print '(A,1X,F10.4)', 'PROFILE '//trim(profile_tag)//' syn_state_s', t_syn_state
        print '(A,1X,F10.4)', 'PROFILE '//trim(profile_tag)//' prepare_aux_s', t_prepare_aux
        print '(A,1X,F10.4)', 'PROFILE '//trim(profile_tag)//' cooling_s', t_cooling
        print '(A,1X,F10.4)', 'PROFILE '//trim(profile_tag)//' eta_s', t_eta
        print '(A,1X,F10.4)', 'PROFILE '//trim(profile_tag)//' xi_s', t_xi
        print '(A,1X,I12)', 'PROFILE '//trim(profile_tag)//' total_substeps', total_substeps
        print '(A,1X,I12)', 'PROFILE '//trim(profile_tag)//' max_shell_substeps', max_shell_substeps
        print '(A,1X,I12)', 'PROFILE '//trim(profile_tag)//' shell_cooling_calls', shell_cooling_calls
        print '(A,1X,I12)', 'PROFILE '//trim(profile_tag)//' substep_cooling_calls', substep_cooling_calls
        print '(A,1X,I12)', 'PROFILE '//trim(profile_tag)//' prepare_aux_calls', prepare_aux_calls
        print '(A,1X,I12)', 'PROFILE '//trim(profile_tag)//' history_calls', history_calls
        print '(A,1X,I12)', 'PROFILE '//trim(profile_tag)//' syn_state_calls', syn_state_calls
        print '(A,1X,I12)', 'PROFILE '//trim(profile_tag)//' eta_calls', eta_calls
        print '(A,1X,I12)', 'PROFILE '//trim(profile_tag)//' xi_calls', xi_calls
    end if
    deallocate(dEl, dEL_mean, dEL_mean_shell, kappa2_arr, dN_init, dN_init_log, dF1, &
               shell_population, chi_population, U_log, source_eta1, &
               V_m_chi, V_c_chi, V_a_chi, chi_weight, Epsilon_b_chi, DB_chi, t_decay_chi, &
               eta_grid, chi_grid, eta_face, chi_face, a_arr, dln_a_dR_arr, proper_time_arr, &
               x_face_hist, x_comov_face_hist, &
               x_comov_hist, dx_comov_hist, beta_hist, dN_cell, P_local, V_cool, P_hist, Seed_hist, Tau_hist, &
               P_hist_cool, Seed_hist_cool, Tau_pair_hist_cool, P_eff_cool_chi, Seed_eff_cool_chi, &
               cooling_aux_chi, dEl_chi, dEL_mean_chi, kappa2_chi)
end subroutine fs_electron_transport_2d_core
