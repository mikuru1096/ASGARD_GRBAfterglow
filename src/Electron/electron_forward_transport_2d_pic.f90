! PIC 2D输运核心：复用fullhide_2d历史辐射框架，替换PIC γmax、χ输运、3/5绝热冷却和patchy B/P_hit。
subroutine fs_electron_transport_2d_pic_core(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                             Num_chi,index_Y,index_syn_intger,n_threads, &
                                             gam_e,dN_gam_e,dN_gam_e_total,P_syn,Seed_syn,V_m,V_c,V_a, &
                                             use_charint_transport,profile_tag,electron_pic_uniform_b, &
                                             electron_pic_eta_acc,electron_pic_kappa_diff_scale,electron_pic_bw_factor)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    use electron_common, only: electron_initial_density, electron_initialize_spectrum, electron_unpack_boundary, &
                               electron_initial_grid_gamma, &
                               electron_gamma_m_exact, electron_injection_prefactor, &
                               electron_gamma_c_from_loss_mean, electron_source_bounds
    use electron_injection_profiles, only: electron_build_source_term_exp_cutoff
    use electron_cooling_kernel, only: prepare_forward_cooling_aux_batch
    use electron_radiation_kernel, only: get_syn_state, get_nu_a_2d_cell_path, reduce_syn_shell_from_chi, &
                                         build_reduced_log_grid, project_syn_state_logbands
    use electron_seed_history_kernel, only: integrate_downstream_proper_time, accumulate_comoving_history_fields
    use radiation_common, only: radiation_pair_tau_headon_segment
    use electron_transport_2d_kernel, only: compute_log_chi_geometry, get_shock_transport_state, &
                                             compute_downstream_comoving_grid, bm_beta2_lab, bm_beta2_shock, &
                                             compute_logchi_eta_step_limit
    use electron_transport_common, only: electron_logparabola_peak_frequency, &
                                         electron_active_gamma_hi, electron_active_chi_hi, &
                                         electron_max_xi_coeff_uniform
    implicit real(8)(a-h,o-z)

    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,Num_chi,index_Y,index_syn_intger,n_threads
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    logical, intent(in) :: use_charint_transport, electron_pic_uniform_b
    character(len=*), intent(in) :: profile_tag
    real(8), intent(in) :: electron_pic_eta_acc,electron_pic_kappa_diff_scale,electron_pic_bw_factor
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
    real(8), allocatable :: B_bg_chi(:), B_strong_chi(:), lambda_chi(:), UBcool_chi(:), P_hit(:,:), Bcool_eff2(:,:)
    real(8), allocatable :: dN_bg(:), dN_strong(:), P_emit_bg(:), P_emit_strong(:), P_syn_bg(:), P_syn_strong(:)
    real(8), allocatable :: Seed_bg(:), Seed_strong(:), Tau_bg(:), Tau_strong(:)
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
    real(8) :: R_loc, R_Gamma_loc, dNe, Para_N_e_ini, DB
    real(8) :: Epsilon_b_floor, magnetic_decay_alpha_t, magnetic_decay_t0_s
    real(8) :: Gam_e_max, Gam_e_m, Gam_e_c, Gam_e_c_diag, temp_gam
    real(8) :: Gam_e_max_cell, Gam_e_m_cell, Gam_e_c_cell
    real(8) :: beta_sh, beta_2, beta_2_sh
    real(8) :: Bacc_loc, omega_p_up, Lp_loc, x_max_loc, gamma_ref_hit, B_for_cooling
    real(8) :: Q
    real(8) :: t_start, t_stop
    real(8) :: t_hist_accum, t_syn_state, t_prepare_aux, t_cooling, t_eta, t_xi
    integer :: I_tobs, I_chi, Num_nu_cool
    integer :: total_substeps, max_shell_substeps, shell_cooling_calls, substep_cooling_calls
    integer :: prepare_aux_calls, history_calls, syn_state_calls, eta_calls, xi_calls
    logical :: profile_enabled, magnetic_decay_active

    allocate(dEl(Num_gam_e), dEL_mean(Num_gam_e-1), dEL_mean_shell(Num_gam_e-1), kappa2_arr(Num_gam_e), &
             dN_init(Num_gam_e), dN_init_log(Num_gam_e), dF1(Num_gam_e), shell_population(Num_gam_e), chi_population(Num_chi), &
             V_m_chi(Num_chi), V_c_chi(Num_chi), V_a_chi(Num_chi), chi_weight(Num_chi), &
             Epsilon_b_chi(Num_chi), DB_chi(Num_chi), t_decay_chi(Num_chi), &
             B_bg_chi(Num_chi), B_strong_chi(Num_chi), lambda_chi(Num_chi), UBcool_chi(Num_chi), &
             P_hit(Num_gam_e,Num_chi), Bcool_eff2(Num_gam_e,Num_chi), &
             dN_bg(Num_gam_e), dN_strong(Num_gam_e), P_emit_bg(Num_nu), P_emit_strong(Num_nu), &
             P_syn_bg(Num_nu), P_syn_strong(Num_nu), Seed_bg(Num_nu), Seed_strong(Num_nu), &
             Tau_bg(Num_nu), Tau_strong(Num_nu), &
             source_eta1(Num_gam_e), U_log(Num_gam_e, Num_chi), eta_grid(Num_chi), &
             chi_grid(Num_chi), eta_face(0:Num_chi), chi_face(0:Num_chi), &
             a_arr(Num_R), dln_a_dR_arr(Num_R), proper_time_arr(Num_R), x_face_hist(0:Num_chi, Num_R), &
             x_comov_face_hist(0:Num_chi, Num_R), x_comov_hist(Num_chi, Num_R), dx_comov_hist(Num_chi, Num_R), &
             beta_hist(Num_chi, Num_R), dN_cell(Num_gam_e), P_local(Num_nu, Num_chi), kappa2_chi(Num_gam_e, Num_chi), &
             cooling_aux_chi(Num_gam_e, Num_chi), dEl_chi(Num_gam_e, Num_chi), dEL_mean_chi(Num_gam_e-1, Num_chi))

    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    Epsilon_b_floor  = Epsilon_b
    magnetic_decay_alpha_t = zero
    magnetic_decay_t0_s = one
    if (n >= 27) then
        Epsilon_b_floor = Boundary(24)
        magnetic_decay_alpha_t = Boundary(25)
        magnetic_decay_t0_s = Boundary(26)
    end if
    magnetic_decay_active = magnetic_decay_alpha_t < zero .and. magnetic_decay_t0_s > zero .and. &
                            Epsilon_b_floor > zero .and. Epsilon_b_floor < Epsilon_b

    ln10 = dlog(ten)
    profile_enabled = .false.
    block
        integer :: env_len, env_status
        character(len=32) :: profile_env

        profile_env = ''
        call get_environment_variable('ASGARD_PROFILE_2D', profile_env, length=env_len, status=env_status)
        if (env_status == 0 .and. env_len > 0) then
            if (profile_env(1:1) /= '0') profile_enabled = .true.
        end if
    end block
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

    call pic_build_field_state(Num_chi,chi_grid,R(1),R_Gamma(1),dNe,Epsilon_b,electron_pic_uniform_b, &
                               electron_pic_bw_factor,Epsilon_b_chi,DB_chi,B_bg_chi,B_strong_chi,lambda_chi, &
                               UBcool_chi,omega_p_up,Lp_loc,x_max_loc)
    Bacc_loc      = B_bg_chi(1)
    DB            = DB_chi(1)
    call pic_gamma_e_max(R_Gamma(1),dNe,Bacc_loc,electron_pic_eta_acc,Gam_e_max)
    temp_gam      = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-one)
    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c       = 7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)

    block
        real(8) :: dNe_end, DB_min, Gam_e_max_max

        call dynamics_external_density_profile(A_star,dNe_ISM,R(Num_R),R0,1,R_tr,f_jump,f_wide,dNe_end)
        DB_min = pic_background_B0(dNe_end,Epsilon_b,R_Gamma(Num_R))
        call pic_gamma_e_max(R_Gamma(Num_R),dNe_end,DB_min,electron_pic_eta_acc,Gam_e_max_max)
        call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                          electron_initial_grid_gamma,gam_e,dN_init)
    end block
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
    call pic_build_field_state(Num_chi,chi_grid,R(1),R_Gamma(1),dNe,Epsilon_b,electron_pic_uniform_b, &
                               electron_pic_bw_factor,Epsilon_b_chi,DB_chi,B_bg_chi,B_strong_chi,lambda_chi, &
                               UBcool_chi,omega_p_up,Lp_loc,x_max_loc)
    call pic_build_hit_probability(Num_gam_e,Num_chi,gam_e,B_bg_chi(1),electron_pic_bw_factor,omega_p_up, &
                                   electron_pic_uniform_b,B_bg_chi,B_strong_chi,P_hit,Bcool_eff2)
    do I_chi = 1, Num_chi
        beta_hist(I_chi,1) = bm_beta2_lab(R_Gamma(1),chi_grid(I_chi))
    end do
    do I_chi = 1, Num_chi
        dN_cell = dN_gam_e(:,I_chi,1)
        if (profile_enabled) call cpu_time(t_start)
        call pic_syn_state_with_hit(electron_pic_uniform_b,R(1),B_bg_chi(I_chi),B_strong_chi(I_chi),Num_gam_e,Num_nu, &
                                    n_threads,gam_e,dN_cell,V_seed,P_hit(:,I_chi),dN_bg,dN_strong,P_emit_bg,P_emit_strong, &
                                    P_syn_bg,P_syn_strong,Seed_bg,Seed_strong,Tau_bg,Tau_strong, &
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

        call pic_build_field_state(Num_chi,chi_grid,R_loc,R_Gamma_loc,dNe,Epsilon_b,electron_pic_uniform_b, &
                                   electron_pic_bw_factor,Epsilon_b_chi,DB_chi,B_bg_chi,B_strong_chi,lambda_chi, &
                                   UBcool_chi,omega_p_up,Lp_loc,x_max_loc)
        call pic_build_hit_probability(Num_gam_e,Num_chi,gam_e,B_bg_chi(1),electron_pic_bw_factor,omega_p_up, &
                                       electron_pic_uniform_b,B_bg_chi,B_strong_chi,P_hit,Bcool_eff2)
        DB        = DB_chi(1)
        Bacc_loc  = B_bg_chi(1)
        call pic_gamma_e_max(R_Gamma_loc,dNe,Bacc_loc,electron_pic_eta_acc,Gam_e_max)
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
            B_for_cooling = dsqrt(max(sum(Bcool_eff2(:,I_chi)*max(U_log(:,I_chi),zero)) / &
                                      max(sum(max(U_log(:,I_chi),zero)),tiny(one)),tiny(one)))
            Gam_e_max_cell = Gam_e_max
            temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
            call electron_gamma_m_exact(p,temp_gam,Gam_e_max_cell,Gam_e_m_cell)
            Gam_e_c_cell = 7.7d8*(one+z)/R_Gamma_loc/B_for_cooling**2/max(R_Tobs(I_tobs),tiny(one))
            call pic_assemble_cooling_column(index_Y,Epsilon_e,Epsilon_b_chi(I_chi),p,B_for_cooling, &
                                             Gam_e_m_cell,Gam_e_c_cell,Gam_e_max_cell,R_loc,R_Gamma_loc,beta_sh,dNe, &
                                             Num_gam_e,Num_nu_cool,n_threads,gam_e,V_cool,Seed_eff_cool_chi(:,I_chi), &
                                             cooling_aux_chi(:,I_chi),Bcool_eff2(:,I_chi),dEl_chi(:,I_chi))
            dEL_mean_chi(:,I_chi) = (dEl_chi(2:Num_gam_e,I_chi)+dEl_chi(1:Num_gam_e-1,I_chi))/two/dlog(ten)
        end do
        call pic_kappa2_patchy(Num_gam_e,Num_chi,gam_e,R_Gamma_loc,omega_p_up,electron_pic_kappa_diff_scale,kappa2_chi)
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
            V_m_chi(I_chi) = electron_logparabola_peak_frequency(Num_nu,V_seed,P_hist(:,I_chi,I_tobs-1))
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
        block
            real(8) :: dDD, dDR, dDR_try, dDR_xi, dDR_eta, max_xi_coeff
            real(8) :: chi_peak
            integer :: L, L1, src_lo, src_hi, active_hi, active_chi_hi

            call electron_source_bounds(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,src_lo,src_hi)
            call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,one,p,dF1)
            shell_population = sum(U_log, dim=2)
            chi_population = sum(U_log, dim=1)
            chi_peak = maxval(shell_population)
            active_hi = electron_active_gamma_hi(Num_gam_e,dF1,shell_population,src_lo,src_hi,chi_peak)
            max_xi_coeff = electron_max_xi_coeff_uniform(Num_gam_e,Num_chi,dEL_mean_chi, &
                                                         (3d0/5d0)/R_loc/ln10,chi_population,chi_peak,active_hi)

            dDR_xi = huge(one)
            if (max_xi_coeff > zero) dDR_xi = 0.4d0*d_x_E/max_xi_coeff
            active_chi_hi = electron_active_chi_hi(Num_chi,chi_population,chi_peak)
            dDR_eta = compute_logchi_eta_step_limit(Num_chi,active_chi_hi,R_loc,R_Gamma_loc,beta_sh, &
                                                    dln_a_dR_arr(I_tobs-1),deta,chi_face,0.4d0)

            dDD     = R(I_tobs)-R(I_tobs-1)
            dDR_try = min(dDD, min(dDR_xi, dDR_eta))
            L1      = max(100, min(1000, ceiling(dDD/max(dDR_try, tiny(one)))))
            dDR     = dDD/dble(L1)
            total_substeps = total_substeps + L1
            max_shell_substeps = max(max_shell_substeps, L1)

            do L = 1, L1
                block
                    real(8) :: frac_sub, R_sub, Gamma_sh_sub, a_sub, dln_a_dR_sub, Gam_e_m_p

                    frac_sub = (dble(L)-0.5d0)/dble(L1)
                    R_sub = R(I_tobs-1) + frac_sub*dDD
                    Gamma_sh_sub = (one-frac_sub)*R_Gamma(I_tobs-1) + frac_sub*R_Gamma(I_tobs)
                    a_sub = 8d0*Gamma_sh_sub*Gamma_sh_sub/R_sub
                    dln_a_dR_sub = (one-frac_sub)*dln_a_dR_arr(I_tobs-1) + frac_sub*dln_a_dR_arr(I_tobs)

                    call dynamics_external_density_profile(A_star,dNe_ISM,R_sub,R0,1,R_tr,f_jump,f_wide,dNe)

                    call pic_build_field_state(Num_chi,chi_grid,R_sub,Gamma_sh_sub,dNe,Epsilon_b,electron_pic_uniform_b, &
                                               electron_pic_bw_factor,Epsilon_b_chi,DB_chi,B_bg_chi,B_strong_chi,lambda_chi, &
                                               UBcool_chi,omega_p_up,Lp_loc,x_max_loc)
                    DB = DB_chi(1)
                    call pic_gamma_e_max(Gamma_sh_sub,dNe,B_bg_chi(1),electron_pic_eta_acc,Gam_e_max)
                    temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(Gamma_sh_sub-one)
                    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
                    Gam_e_c = 7.7d8*(one+z)/Gamma_sh_sub/DB**2/max(R_Tobs(I_tobs),tiny(one))
                    call get_shock_transport_state(Gamma_sh_sub, beta_sh, beta_2, beta_2_sh)

                    Gam_e_m_p = (one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
                    call electron_injection_prefactor(R_sub,dDR,dNe,f_e,Gam_e_m_p,Q)
                    call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,Q,p,dF1)
                    source_eta1 = dF1/deta

                    if (profile_enabled) call cpu_time(t_start)
                    call pic_advance_eta_advection_implicit(U_log,Num_gam_e,Num_chi,active_hi,deta,chi_face,Gamma_sh_sub, &
                                                            a_sub,dln_a_dR_sub,beta_sh,dDR,n_threads)
                    call pic_advance_eta_diffusion_implicit(U_log,Num_gam_e,Num_chi,active_hi,deta,chi_face,Gamma_sh_sub, &
                                                            a_sub,beta_sh,kappa2_chi,dDR,n_threads)
                    eta_calls = eta_calls + 1
                    if (profile_enabled) then
                        call cpu_time(t_stop)
                        t_eta = t_eta + (t_stop-t_start)
                    end if

                    if (profile_enabled) call cpu_time(t_start)
                    call pic_advance_energy_loggamma(U_log,Num_gam_e,Num_chi,dEL_mean_chi,R_sub,d_x_E,dDR,source_eta1,n_threads)
                    xi_calls = xi_calls + 1
                    if (profile_enabled) then
                        call cpu_time(t_stop)
                        t_xi = t_xi + (t_stop-t_start)
                    end if
                end block
            end do
        end block

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
        call pic_build_field_state(Num_chi,chi_grid,R(I_tobs),R_Gamma(I_tobs),dNe,Epsilon_b,electron_pic_uniform_b, &
                                   electron_pic_bw_factor,Epsilon_b_chi,DB_chi,B_bg_chi,B_strong_chi,lambda_chi, &
                                   UBcool_chi,omega_p_up,Lp_loc,x_max_loc)
        call pic_build_hit_probability(Num_gam_e,Num_chi,gam_e,B_bg_chi(1),electron_pic_bw_factor,omega_p_up, &
                                       electron_pic_uniform_b,B_bg_chi,B_strong_chi,P_hit,Bcool_eff2)
        do I_chi = 1, Num_chi
            beta_hist(I_chi,I_tobs) = bm_beta2_lab(R_Gamma(I_tobs),chi_grid(I_chi))
            dN_cell = dN_gam_e(:,I_chi,I_tobs)
            if (profile_enabled) call cpu_time(t_start)
            call pic_syn_state_with_hit(electron_pic_uniform_b,R(I_tobs),B_bg_chi(I_chi),B_strong_chi(I_chi),Num_gam_e, &
                                        Num_nu,n_threads,gam_e,dN_cell,V_seed,P_hit(:,I_chi),dN_bg,dN_strong, &
                                        P_emit_bg,P_emit_strong,P_syn_bg,P_syn_strong,Seed_bg,Seed_strong,Tau_bg,Tau_strong, &
                                        P_local(:,I_chi),P_hist(:,I_chi,I_tobs),Seed_hist(:,I_chi,I_tobs), &
                                        Tau_hist(:,I_chi,I_tobs))
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
    call pic_build_field_state(Num_chi,chi_grid,R_loc,R_Gamma_loc,dNe,Epsilon_b,electron_pic_uniform_b, &
                               electron_pic_bw_factor,Epsilon_b_chi,DB_chi,B_bg_chi,B_strong_chi,lambda_chi, &
                               UBcool_chi,omega_p_up,Lp_loc,x_max_loc)
    call pic_build_hit_probability(Num_gam_e,Num_chi,gam_e,B_bg_chi(1),electron_pic_bw_factor,omega_p_up, &
                                   electron_pic_uniform_b,B_bg_chi,B_strong_chi,P_hit,Bcool_eff2)
    DB = DB_chi(1)
    call pic_gamma_e_max(R_Gamma_loc,dNe,B_bg_chi(1),electron_pic_eta_acc,Gam_e_max)
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
        B_for_cooling = dsqrt(max(sum(Bcool_eff2(:,I_chi)*max(U_log(:,I_chi),zero)) / &
                                  max(sum(max(U_log(:,I_chi),zero)),tiny(one)),tiny(one)))
        Gam_e_max_cell = Gam_e_max
        temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max_cell,Gam_e_m_cell)
        Gam_e_c_cell = 7.7d8*(one+z)/R_Gamma_loc/B_for_cooling**2/max(R_Tobs(Num_R),tiny(one))
        call pic_assemble_cooling_column(index_Y,Epsilon_e,Epsilon_b_chi(I_chi),p,B_for_cooling, &
                                         Gam_e_m_cell,Gam_e_c_cell,Gam_e_max_cell,R_loc,R_Gamma_loc,beta_sh,dNe, &
                                         Num_gam_e,Num_nu_cool,n_threads,gam_e,V_cool,Seed_eff_cool_chi(:,I_chi), &
                                         cooling_aux_chi(:,I_chi),Bcool_eff2(:,I_chi),dEl_chi(:,I_chi))
        dEL_mean_chi(:,I_chi) = (dEl_chi(2:Num_gam_e,I_chi)+dEl_chi(1:Num_gam_e-1,I_chi))/two/dlog(ten)
    end do
    if (profile_enabled) then
        call cpu_time(t_stop)
        t_cooling = t_cooling + (t_stop-t_start)
    end if
    call get_nu_a_2d_cell_path(Num_nu,Num_chi,V_seed,Tau_hist(:,:,Num_R),V_a_chi)
    temp = zero
    Q = zero
    do I_chi = 1, Num_chi
        chi_weight(I_chi) = max(sum(U_log(:,I_chi)), tiny(one))
        V_m_chi(I_chi) = electron_logparabola_peak_frequency(Num_nu,V_seed,P_hist(:,I_chi,Num_R))
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
               B_bg_chi, B_strong_chi, lambda_chi, UBcool_chi, P_hit, Bcool_eff2, &
               dN_bg, dN_strong, P_emit_bg, P_emit_strong, P_syn_bg, P_syn_strong, &
               Seed_bg, Seed_strong, Tau_bg, Tau_strong, &
               eta_grid, chi_grid, eta_face, chi_face, a_arr, dln_a_dR_arr, proper_time_arr, &
               x_face_hist, x_comov_face_hist, &
               x_comov_hist, dx_comov_hist, beta_hist, dN_cell, P_local, V_cool, P_hist, Seed_hist, Tau_hist, &
               P_hist_cool, Seed_hist_cool, Tau_pair_hist_cool, P_eff_cool_chi, Seed_eff_cool_chi, &
               cooling_aux_chi, dEl_chi, dEL_mean_chi, kappa2_chi)
end subroutine fs_electron_transport_2d_pic_core


subroutine pic_build_field_state(Num_chi,chi_grid,R_loc,Gamma_sh,dNe_up,epsB_uniform,uniform_b,Bw_factor, &
                                 epsB_chi,Bemit_chi,B_bg_chi,B_strong_chi,lambda_chi,UBcool_chi, &
                                 omega_p_up,Lp_loc_out,x_max_loc_out)
    use constants
    implicit real(8)(a-h,o-z)

    integer, intent(in) :: Num_chi
    real(8), intent(in) :: chi_grid(Num_chi), R_loc, Gamma_sh, dNe_up, epsB_uniform, Bw_factor
    logical, intent(in) :: uniform_b
    real(8), intent(out) :: epsB_chi(Num_chi), Bemit_chi(Num_chi), B_bg_chi(Num_chi), B_strong_chi(Num_chi)
    real(8), intent(out) :: lambda_chi(Num_chi), UBcool_chi(Num_chi), omega_p_up, Lp_loc_out, x_max_loc_out

    real(8) :: e2_sh, lambda_scatt, x_profile(Num_chi), x_max_loc, B0_loc, Bw_loc, Lp_loc, B_loc
    integer :: I_chi

    omega_p_up = dsqrt(4d0*pi*Para_e*Para_e*max(dNe_up,tiny(one))/Para_m_p/max(Gamma_sh,one))
    lambda_scatt = 100d0*Para_c/max(omega_p_up,tiny(one))
    B0_loc = pic_background_B0(dNe_up,epsB_uniform,Gamma_sh)
    e2_sh = 4d0*max(Gamma_sh,one)*(max(Gamma_sh,one)-one)*max(dNe_up,zero)*Para_m_p_E

    do I_chi = 1, Num_chi
        x_profile(I_chi) = R_loc*(chi_grid(I_chi)-one)/(8d0*max(Gamma_sh*Gamma_sh,tiny(one)))
    end do
    x_max_loc = max(maxval(x_profile),tiny(one))
    Bw_loc = max(Bw_factor,tiny(one))*B0_loc
    Lp_loc = min(lambda_scatt,x_max_loc/(Bw_loc/B0_loc)**2)
    Lp_loc = max(Lp_loc,tiny(one))

    do I_chi = 1, Num_chi
        if (uniform_b) then
            B_loc = B0_loc
        else
            B_loc = Bw_loc*(one+x_profile(I_chi)/Lp_loc)**(-0.6d0) + B0_loc
        end if
        B_bg_chi(I_chi) = B0_loc
        B_strong_chi(I_chi) = max(B_loc,tiny(one))
        Bemit_chi(I_chi) = B_strong_chi(I_chi)
        UBcool_chi(I_chi) = Bemit_chi(I_chi)*Bemit_chi(I_chi)/(8d0*pi)
        epsB_chi(I_chi) = max(UBcool_chi(I_chi)/max(e2_sh,tiny(one)),tiny(one))
        lambda_chi(I_chi) = lambda_scatt
        if (uniform_b) epsB_chi(I_chi) = max(epsB_uniform,tiny(one))
    end do

    Lp_loc_out = Lp_loc
    x_max_loc_out = x_max_loc
end subroutine pic_build_field_state


real(8) function pic_background_B0(dNe_up,epsB_bg,Gamma_ref)
    use constants
    implicit real(8)(a-h,o-z)
    real(8), intent(in) :: dNe_up,epsB_bg,Gamma_ref

    pic_background_B0 = 0.39d0*dsqrt(max(epsB_bg,tiny(one))*max(dNe_up,zero)* &
                                     (max(Gamma_ref,one)*(max(Gamma_ref,one)-one)))
    pic_background_B0 = max(pic_background_B0,tiny(one))
end function pic_background_B0


subroutine pic_gamma_e_max(Gamma_bulk,dNe_up,Bprime_accel,eta_acc,gamma_e_max_pic)
    use constants
    implicit real(8)(a-h,o-z)
    real(8), intent(in) :: Gamma_bulk,dNe_up,Bprime_accel,eta_acc
    real(8), intent(out) :: gamma_e_max_pic
    real(8) :: omega_p_up,accel_prefactor

    omega_p_up = dsqrt(4d0*pi*Para_e*Para_e*max(dNe_up,tiny(one))/Para_m_p/max(Gamma_bulk,one))
    accel_prefactor = 0.6d0*pi*Para_m_p*Para_m_p*omega_p_up*max(Gamma_bulk,one)*max(Gamma_bulk,one)*Para_c
    gamma_e_max_pic = (accel_prefactor/(Para_SigmaT*Para_m_e*max(Bprime_accel,tiny(one))**2* &
                       max(eta_acc,tiny(one))))**(one/3d0)
end subroutine pic_gamma_e_max


subroutine pic_build_hit_probability(Num_gam_e,Num_chi,gam_e,Bacc_loc,Bw_factor,omega_p_up,uniform_b, &
                                     B_bg_chi,B_strong_chi,P_hit,Bcool_eff2)
    use constants
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: Num_gam_e,Num_chi
    real(8), intent(in) :: gam_e(Num_gam_e),Bacc_loc,Bw_factor,omega_p_up
    real(8), intent(in) :: B_bg_chi(Num_chi),B_strong_chi(Num_chi)
    logical, intent(in) :: uniform_b
    real(8), intent(out) :: P_hit(Num_gam_e,Num_chi),Bcool_eff2(Num_gam_e,Num_chi)
    real(8), parameter :: P_max_hit = 0.01d0, q_cov_hit = 2d0, alpha_hit = 3d0
    real(8) :: gamma_ref_hit, spatial_norm, spatial_factor, energy_factor, b_bg2, b_strong2
    integer :: I_gam_e,I_chi

    P_hit = zero
    do I_chi = 1, Num_chi
        Bcool_eff2(:,I_chi) = B_bg_chi(I_chi)*B_bg_chi(I_chi)
    end do
    if (uniform_b) return

    gamma_ref_hit = 100d0*Para_e*Bacc_loc*Bacc_loc / &
                    max(two*pi*Para_m_e*Para_c*omega_p_up*max(Bw_factor,tiny(one))*Bacc_loc,tiny(one))
    gamma_ref_hit = max(gamma_ref_hit,one)
    spatial_norm = max(B_strong_chi(1)-B_bg_chi(1),tiny(one))
    do I_chi = 1, Num_chi
        b_bg2 = B_bg_chi(I_chi)*B_bg_chi(I_chi)
        b_strong2 = B_strong_chi(I_chi)*B_strong_chi(I_chi)
        spatial_factor = ((B_strong_chi(I_chi)-B_bg_chi(I_chi))/spatial_norm)**q_cov_hit
        spatial_factor = min(one,max(zero,spatial_factor))
        do I_gam_e = 1, Num_gam_e
            energy_factor = (gam_e(I_gam_e)/gamma_ref_hit)**alpha_hit
            energy_factor = energy_factor/(one+energy_factor)
            P_hit(I_gam_e,I_chi) = min(P_max_hit,max(zero,P_max_hit*spatial_factor*energy_factor))
            Bcool_eff2(I_gam_e,I_chi) = (one-P_hit(I_gam_e,I_chi))*b_bg2 + P_hit(I_gam_e,I_chi)*b_strong2
        end do
    end do
end subroutine pic_build_hit_probability


subroutine pic_syn_state_with_hit(uniform_b,R_loc,B_bg,B_strong,Num_gam_e,Num_nu,n_threads,gam_e,dN_cell,V_seed,P_hit_col, &
                                  dN_bg,dN_strong,P_emit_bg,P_emit_strong,P_syn_bg,P_syn_strong,Seed_bg,Seed_strong, &
                                  Tau_bg,Tau_strong,P_emit,P_syn,Seed_syn,Tau_syn)
    use constants
    use electron_radiation_kernel, only: get_syn_state
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: Num_gam_e,Num_nu,n_threads
    real(8), intent(in) :: R_loc,B_bg,B_strong,gam_e(Num_gam_e),dN_cell(Num_gam_e),V_seed(Num_nu),P_hit_col(Num_gam_e)
    logical, intent(in) :: uniform_b
    real(8), intent(inout) :: dN_bg(Num_gam_e),dN_strong(Num_gam_e),P_emit_bg(Num_nu),P_emit_strong(Num_nu)
    real(8), intent(inout) :: P_syn_bg(Num_nu),P_syn_strong(Num_nu),Seed_bg(Num_nu),Seed_strong(Num_nu)
    real(8), intent(inout) :: Tau_bg(Num_nu),Tau_strong(Num_nu)
    real(8), intent(out) :: P_emit(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu),Tau_syn(Num_nu)

    if (uniform_b) then
        call get_syn_state(R_loc,B_bg,Num_gam_e,Num_nu,n_threads,gam_e,dN_cell,V_seed,P_emit,P_syn,Seed_syn,Tau_syn)
        return
    end if

    dN_bg = (one-P_hit_col)*dN_cell
    dN_strong = P_hit_col*dN_cell
    call get_syn_state(R_loc,B_bg,Num_gam_e,Num_nu,n_threads,gam_e,dN_bg,V_seed,P_emit_bg,P_syn_bg,Seed_bg,Tau_bg)
    call get_syn_state(R_loc,B_strong,Num_gam_e,Num_nu,n_threads,gam_e,dN_strong,V_seed, &
                       P_emit_strong,P_syn_strong,Seed_strong,Tau_strong)
    P_emit = P_emit_bg + P_emit_strong
    P_syn = P_syn_bg + P_syn_strong
    Seed_syn = Seed_bg + Seed_strong
    Tau_syn = Tau_bg + Tau_strong
end subroutine pic_syn_state_with_hit


subroutine pic_assemble_cooling_column(index_Y,Epsilon_e,Epsilon_b,p,B_for_cooling,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                       R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn_ssa, &
                                       cooling_aux,Bcool_eff2_col,dEl)
    use constants
    use electron_cooling_kernel, only: electron_cooling_ssa_loss, electron_cooling_y_fan
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: index_Y,Num_gam_e,Num_nu,n_threads
    real(8), intent(in) :: Epsilon_e,Epsilon_b,p,B_for_cooling,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
    real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn_ssa(Num_nu),cooling_aux(Num_gam_e)
    real(8), intent(in) :: Bcool_eff2_col(Num_gam_e)
    real(8), intent(inout) :: Gam_e_max
    real(8), intent(out) :: dEl(Num_gam_e)
    real(8) :: Compton(Num_gam_e),dot_gam_e_SSA(Num_gam_e)
    real(8) :: cooling_scale,ssa_scale,f_r(Num_gam_e),Q

    cooling_scale = one/(beta_Gam*R_Gamma_loc)
    ssa_scale = cooling_scale/Para_c
    f_r = 1.35d-19*Bcool_eff2_col*cooling_scale/pi
    call electron_cooling_ssa_loss(B_for_cooling,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn_ssa,dot_gam_e_SSA)

    select case(index_Y)
    case(0)
        dEl = (f_r-dot_gam_e_SSA*ssa_scale)*gam_e
    case(1)
        dEl = (f_r+(cooling_aux-dot_gam_e_SSA)*ssa_scale)*gam_e
    case(2)
        Q = 4d0*pi*R_loc*R_loc*Para_c
        Compton = one+cooling_aux/Q/(4d0*R_Gamma_loc*R_Gamma_loc*dNe*Para_m_p_E)
        Gam_e_max = Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl = (f_r*Compton-dot_gam_e_SSA*ssa_scale)*gam_e
    case(3)
        call electron_cooling_y_fan(Epsilon_e,Epsilon_b,p,B_for_cooling,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,Compton)
        Compton = one+Compton
        Gam_e_max = Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl = (f_r*Compton-dot_gam_e_SSA*ssa_scale)*gam_e
    case default
        print*, 'invalid Compton case, check your chosen model!'
        stop
    end select
end subroutine pic_assemble_cooling_column


subroutine pic_kappa2_patchy(Num_gam_e,Num_chi,gam_e,Gamma_bulk,omega_p_up,kappa_diff_scale,kappa2_arr)
    use constants
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: Num_gam_e,Num_chi
    real(8), intent(in) :: gam_e(Num_gam_e),Gamma_bulk,omega_p_up,kappa_diff_scale
    real(8), intent(out) :: kappa2_arr(Num_gam_e,Num_chi)
    real(8) :: gamma0_bulk,nu_scatt,kappa_gamma
    integer :: I_gam_e,I_chi

    gamma0_bulk = max(Gamma_bulk,one)
    do I_gam_e = 1, Num_gam_e
        nu_scatt = max(0.4d0*(gam_e(I_gam_e)/gamma0_bulk)**(-2d0)*omega_p_up,tiny(one))
        kappa_gamma = max(kappa_diff_scale,tiny(one))*Para_c*Para_c/(3d0*nu_scatt)
        do I_chi = 1, Num_chi
            kappa2_arr(I_gam_e,I_chi) = kappa_gamma
        end do
    end do
end subroutine pic_kappa2_patchy


subroutine pic_solve_tridiagonal(n_eq,lower,diag,upper,rhs,sol)
    use constants
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: n_eq
    real(8), intent(in) :: lower(n_eq),diag(n_eq),upper(n_eq),rhs(n_eq)
    real(8), intent(out) :: sol(n_eq)
    real(8) :: diag_work(n_eq),rhs_work(n_eq),pivot
    integer :: I_eq

    diag_work = diag
    rhs_work = rhs
    do I_eq = 2, n_eq
        pivot = lower(I_eq)/max(diag_work(I_eq-1),tiny(one))
        diag_work(I_eq) = diag_work(I_eq)-pivot*upper(I_eq-1)
        rhs_work(I_eq) = rhs_work(I_eq)-pivot*rhs_work(I_eq-1)
    end do
    sol(n_eq) = rhs_work(n_eq)/max(diag_work(n_eq),tiny(one))
    do I_eq = n_eq-1, 1, -1
        sol(I_eq) = (rhs_work(I_eq)-upper(I_eq)*sol(I_eq+1))/max(diag_work(I_eq),tiny(one))
    end do
end subroutine pic_solve_tridiagonal


subroutine pic_advance_eta_advection_implicit(U_log,Num_gam_e,Num_chi,active_hi,deta,chi_face,Gamma_sh,a_loc,dln_a_dR_loc, &
                                             beta_sh,dR_step,n_threads)
    use constants
    use electron_transport_2d_kernel, only: bm_beta2_shock
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: Num_gam_e,Num_chi,active_hi,n_threads
    real(8), intent(inout) :: U_log(Num_gam_e,Num_chi)
    real(8), intent(in) :: deta,chi_face(0:Num_chi),Gamma_sh,a_loc,dln_a_dR_loc,beta_sh,dR_step
    real(8) :: A_face(0:Num_chi),flux_coeff(0:Num_chi),lower(Num_chi),diag(Num_chi),upper(Num_chi)
    real(8) :: rhs(Num_chi),sol(Num_chi),val
    integer :: I_face,I_chi,I_gam_e,hi

    A_face = zero
    flux_coeff = zero
    do I_face = 1, Num_chi-1
        A_face(I_face) = (a_loc*bm_beta2_shock(Gamma_sh,chi_face(I_face))/(chi_face(I_face)*beta_sh) + &
                         ((chi_face(I_face)-one)/chi_face(I_face))*dln_a_dR_loc)/dlog(ten)
        flux_coeff(I_face) = dR_step*A_face(I_face)/deta
    end do
    hi = min(max(active_hi,2),Num_gam_e)
    !$omp parallel do private(I_gam_e,I_chi,lower,diag,upper,rhs,sol,val) if(n_threads>1) num_threads(n_threads)
    do I_gam_e = 1, hi
        lower = zero
        diag = one
        upper = zero
        rhs = U_log(I_gam_e,:)
        do I_chi = 1, Num_chi-1
            val = flux_coeff(I_chi)
            if (val >= zero) then
                diag(I_chi) = diag(I_chi) + val
                lower(I_chi+1) = lower(I_chi+1) - val
            else
                upper(I_chi) = upper(I_chi) + val
                diag(I_chi+1) = diag(I_chi+1) - val
            end if
        end do
        call pic_solve_tridiagonal(Num_chi,lower,diag,upper,rhs,sol)
        U_log(I_gam_e,:) = max(sol,zero)
    end do
    !$omp end parallel do
end subroutine pic_advance_eta_advection_implicit


subroutine pic_advance_eta_diffusion_implicit(U_log,Num_gam_e,Num_chi,active_hi,deta,chi_face,Gamma_sh,a_loc,beta_sh, &
                                             kappa2_chi,dR_step,n_threads)
    use constants
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: Num_gam_e,Num_chi,active_hi,n_threads
    real(8), intent(inout) :: U_log(Num_gam_e,Num_chi)
    real(8), intent(in) :: deta,chi_face(0:Num_chi),Gamma_sh,a_loc,beta_sh,kappa2_chi(Num_gam_e,Num_chi),dR_step
    real(8) :: lower(Num_chi),diag(Num_chi),upper(Num_chi),rhs(Num_chi),sol(Num_chi)
    real(8) :: kappa_face,diff_coeff
    integer :: I_face,I_chi,I_gam_e,hi

    hi = min(max(active_hi,2),Num_gam_e)
    !$omp parallel do private(I_gam_e,I_face,I_chi,lower,diag,upper,rhs,sol,kappa_face,diff_coeff) &
    !$omp if(n_threads>1) num_threads(n_threads)
    do I_gam_e = 1, hi
        lower = zero
        diag = one
        upper = zero
        rhs = U_log(I_gam_e,:)
        do I_face = 1, Num_chi-1
            kappa_face = 0.5d0*(kappa2_chi(I_gam_e,I_face)+kappa2_chi(I_gam_e,I_face+1))
            diff_coeff = dR_step*a_loc*a_loc*kappa_face/(chi_face(I_face)*chi_face(I_face)* &
                         beta_sh*Para_c*dlog(ten)*dlog(ten)*deta*deta)
            diag(I_face) = diag(I_face) + diff_coeff
            upper(I_face) = upper(I_face) - diff_coeff
            lower(I_face+1) = lower(I_face+1) - diff_coeff
            diag(I_face+1) = diag(I_face+1) + diff_coeff
        end do
        call pic_solve_tridiagonal(Num_chi,lower,diag,upper,rhs,sol)
        U_log(I_gam_e,:) = max(sol,zero)
    end do
    !$omp end parallel do
end subroutine pic_advance_eta_diffusion_implicit


subroutine pic_advance_energy_loggamma(U_log,Num_gam_e,Num_chi,dEL_mean_2d,R_loc,d_x_E,dR_step,source_eta1,n_threads)
    use constants
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: Num_gam_e,Num_chi,n_threads
    real(8), intent(inout) :: U_log(Num_gam_e,Num_chi)
    real(8), intent(in) :: dEL_mean_2d(Num_gam_e-1,Num_chi),R_loc,d_x_E,dR_step,source_eta1(Num_gam_e)
    real(8) :: coeff_xi(Num_gam_e-1),up(Num_gam_e-1),principal(Num_gam_e),temp1(Num_gam_e-1)
    real(8) :: lower(Num_gam_e),upper(Num_gam_e)
    real(8) :: rhs(Num_gam_e),sol(Num_gam_e),CFL
    integer :: I_chi

    CFL = dR_step/d_x_E
    !$omp parallel do private(I_chi,coeff_xi,up,principal,temp1,lower,upper,rhs,sol) if(n_threads>1) num_threads(n_threads)
    do I_chi = 1, Num_chi
        coeff_xi = dEL_mean_2d(:,I_chi) + (3d0/5d0)/R_loc/dlog(ten)
        up = -CFL*coeff_xi
        principal(2:Num_gam_e) = one - up
        principal(1) = principal(2)
        temp1 = up/(principal(2:Num_gam_e)+principal(1:Num_gam_e-1))*two
        lower = zero
        upper = zero
        lower(2:Num_gam_e) = -temp1
        upper(1:Num_gam_e-1) = temp1
        rhs = U_log(:,I_chi)
        if (I_chi == 1) rhs = rhs + dR_step*source_eta1
        call pic_solve_tridiagonal(Num_gam_e,lower,principal,upper,rhs,sol)
        U_log(:,I_chi) = max(sol,zero)
    end do
    !$omp end parallel do
end subroutine pic_advance_energy_loggamma
