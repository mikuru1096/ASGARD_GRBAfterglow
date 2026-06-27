! 电子2D输运核心：χ网格构建→壳层循环（历史场叠加+冷却+η对流/扩散+能量维冷却），支持charint/fullhide双模式。
subroutine fs_electron_transport_2d_core(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                         Num_chi,index_Y,index_syn_intger,n_threads,emit_full_chi_spectrum, &
                                         gam_e,dN_gam_e,dN_gam_e_total,P_syn,Seed_syn,V_m,V_c,V_a, &
                                         P_syn_chi,Seed_syn_chi,Tau_syn_chi,chi_radius,chi_gamma_bulk,chi_weight_out, &
                                         B_chi_out,substep_max,use_charint_transport, profile_tag)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    use electron_common, only: electron_initial_density, electron_initialize_spectrum, electron_unpack_boundary, &
                               electron_initial_grid_log_edges, electron_exp_tail_grid_factor, &
                               electron_gamma_m_exact, electron_injection_prefactor, &
                               electron_gamma_c_from_loss_mean, electron_source_bounds
    use electron_injection_profiles, only: electron_build_source_term_exp_cutoff_edges, electron_profile_log_cell_edges, &
                                           electron_initial_powerlaw_exp_cutoff_coord_edges, &
                                           electron_build_source_term_exp_cutoff_coord_edges
    use electron_cooling_kernel, only: prepare_forward_cooling_aux_batch, assemble_forward_cooling_split, &
                                       assemble_forward_cooling_split_batch
    use electron_radiation_kernel, only: get_syn_state, get_syn_selected_state, get_syn_chi_batch_state, &
                                         get_nu_a_2d_cell_path, &
                                         build_reduced_log_grid, project_syn_state_logbands
    use electron_seed_history_kernel, only: integrate_downstream_proper_time, advance_comoving_history_stream
    use radiation_common, only: radiation_pair_tau_headon_segment
    use electron_transport_2d_kernel, only: compute_q_geometry, compute_q_cell_geometry, get_shock_transport_state, &
                                             compute_downstream_comoving_grid, compute_q_divergence, compute_q_step_limit
    use electron_transport_2d_kernel, only: advance_q_implicit, advance_q_advection_charint, &
                                             advance_q_diffusion_implicit, advance_q_pwncr_implicit
    use electron_transport_2d_kernel, only: advance_energy_loggamma_chi, advance_energy_loggamma_chi_charint, &
                                             advance_energy_loggamma_chi_pwncr, advance_energy_stochastic_loggamma_chi
    use electron_transport_common, only: electron_logparabola_peak_frequency, &
                                         electron_active_gamma_hi, electron_active_chi_hi, &
                                         electron_max_xi_coeff_chi, &
                                         electron_dnx_to_dndgamma_exp_centers
    use electron_energy_coordinate_common, only: electron_build_four_velocity_grid, electron_coord_log_four_velocity_sq, &
                                                 electron_dxgamma_dcoord, electron_four_velocity_grid_gamma_scale
    use electron_shell_transport_common, only: electron_shell_flux_split_coord_sequence, &
                                               electron_shell_dcoord_to_dndgamma_exp_centers
    implicit real(8)(a-h,o-z)

    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,Num_chi,index_Y,index_syn_intger,n_threads,substep_max
    integer, intent(in) :: emit_full_chi_spectrum
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    logical, intent(in) :: use_charint_transport
    character(len=*), intent(in) :: profile_tag
    real(8), intent(out) :: dN_gam_e(Num_gam_e, Num_chi, Num_R)
    real(8), intent(out) :: dN_gam_e_total(Num_gam_e, Num_R)
    real(8), intent(out) :: gam_e(Num_gam_e)
    real(8), intent(out) :: P_syn(Num_nu, Num_R)
    real(8), intent(out) :: Seed_syn(Num_nu, Num_R)
    real(8), intent(out) :: V_m(Num_R), V_c(Num_R), V_a(Num_R)
    real(8), intent(out) :: P_syn_chi(Num_nu, Num_chi, Num_R), Seed_syn_chi(Num_nu, Num_chi, Num_R)
    real(8), intent(out) :: Tau_syn_chi(Num_nu, Num_chi, Num_R)
    real(8), intent(out) :: chi_radius(Num_chi, Num_R), chi_gamma_bulk(Num_chi, Num_R), chi_weight_out(Num_chi, Num_R)
    real(8), intent(out) :: B_chi_out(Num_chi, Num_R)

    real(8), allocatable :: dEl(:), dEL_mean(:)
    real(8), allocatable :: dN_init(:), dN_init_log(:), U_log(:,:), U_shell(:,:), source_q1(:)
    real(8), allocatable :: q_grid(:), q_face(:)
    real(8), allocatable :: dF1(:), shell_population(:), chi_population(:), dEL_mean_shell(:)
    real(8), allocatable :: V_m_chi(:), V_c_chi(:), V_a_chi(:), chi_weight(:), Epsilon_b_chi(:), DB_chi(:), t_decay_chi(:)
    real(8), allocatable :: proper_time_arr(:), x_face_hist(:,:), x_comov_face_hist(:,:), &
                             x_comov_hist(:,:), dx_comov_hist(:,:), dN_cell(:)
    real(8), allocatable :: radius_cell_hist(:,:), gamma_cell_hist(:,:), beta_hist(:,:)
    real(8), allocatable :: radius_cell_chi(:), gamma_cell_chi(:), beta_cell_chi(:), beta_rel_sh_chi(:)
    real(8), allocatable :: P_local(:,:), kappa2_chi(:,:)
    real(8), allocatable :: V_cool(:), P_emit_cool(:), P_emit_shell(:), Tau_shell(:)
    real(8), allocatable :: P_hist(:,:,:), Seed_hist(:,:,:), Tau_hist(:,:,:)
    real(8), allocatable :: P_hist_cool(:,:,:), Seed_hist_cool(:,:,:), Tau_hist_cool(:,:,:), &
                             Tau_pair_hist_cool(:,:,:), Tau_prop_hist_cool(:,:,:)
    real(8), allocatable :: P_eff_cool_chi(:,:), Seed_eff_cool_chi(:,:), P_stream_cool(:,:), Seed_stream_cool(:,:)
    real(8), allocatable :: cooling_aux_chi(:,:), dEl_chi(:,:), dEL_mean_chi(:,:), adiabatic_log_coeff_chi(:)

    real(8) :: temp, dq, d_x_E, ln10, coord_scale, x_edge_E(Num_gam_e+1), coord_edge_E(Num_gam_e+1)
    real(8) :: R_loc, R_Gamma_loc, dNe, Para_N_e_ini, DB, DB_min
    real(8) :: Epsilon_b_floor, magnetic_decay_alpha_t, magnetic_decay_t0_s
    real(8) :: stochastic_accel_norm
    real(8) :: Gam_e_max, Gam_e_max_max, Gam_e_m, Gam_e_c, Gam_e_c_diag, temp_gam
    real(8) :: Gam_e_max_cell, Gam_e_m_cell, Gam_e_c_cell
    real(8) :: beta_sh, beta_2, beta_2_sh
    real(8) :: Q, Q_rate
    real(8) :: t_start, t_stop
    real(8) :: t_hist_accum, t_syn_state, t_prepare_aux, t_cooling, t_eta, t_xi
    integer :: I_tobs, I_chi, Num_nu_cool, k_medium, substep_limit
    integer :: total_substeps, max_shell_substeps, shell_cooling_calls, substep_cooling_calls
    integer :: prepare_aux_calls, history_calls, syn_state_calls, eta_calls, xi_calls
    logical :: profile_enabled, magnetic_decay_active, pwn_cr_transport, free_outer_escape, emit_full_spectrum
    logical :: four_velocity_coord

    allocate(dEl(Num_gam_e), dEL_mean(Num_gam_e-1), dEL_mean_shell(Num_gam_e-1), &
             dN_init(Num_gam_e), dN_init_log(Num_gam_e), dF1(Num_gam_e), shell_population(Num_gam_e), chi_population(Num_chi), &
             V_m_chi(Num_chi), V_c_chi(Num_chi), V_a_chi(Num_chi), chi_weight(Num_chi), &
             Epsilon_b_chi(Num_chi), DB_chi(Num_chi), t_decay_chi(Num_chi), &
             source_q1(Num_gam_e), U_log(Num_gam_e, Num_chi), U_shell(Num_gam_e, Num_chi), &
             q_grid(Num_chi), q_face(0:Num_chi), &
             proper_time_arr(Num_R), x_face_hist(0:Num_chi, Num_R), &
             x_comov_face_hist(0:Num_chi, Num_R), x_comov_hist(Num_chi, Num_R), dx_comov_hist(Num_chi, Num_R), &
             radius_cell_hist(Num_chi, Num_R), gamma_cell_hist(Num_chi, Num_R), beta_hist(Num_chi, Num_R), &
             radius_cell_chi(Num_chi), gamma_cell_chi(Num_chi), beta_cell_chi(Num_chi), beta_rel_sh_chi(Num_chi), &
             dN_cell(Num_gam_e), P_local(Num_nu, Num_chi), kappa2_chi(Num_gam_e, Num_chi), &
             cooling_aux_chi(Num_gam_e, Num_chi), dEl_chi(Num_gam_e, Num_chi), dEL_mean_chi(Num_gam_e-1, Num_chi), &
             adiabatic_log_coeff_chi(Num_chi))

    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    if (A_star > zero) then
        k_medium = 2
    else
        k_medium = 0
    end if
    Epsilon_b_floor  = Epsilon_b
    magnetic_decay_alpha_t = zero
    magnetic_decay_t0_s = one
    if (n >= 27) then
        Epsilon_b_floor = Boundary(24)
        magnetic_decay_alpha_t = Boundary(25)
        magnetic_decay_t0_s = Boundary(26)
    end if
    stochastic_accel_norm = zero
    block
        real(8) :: transport_model_selector, escape_mode_selector

        transport_model_selector = zero
        escape_mode_selector = zero
        if (n >= 30) then
            transport_model_selector = Boundary(28)
            stochastic_accel_norm = Boundary(29)
            escape_mode_selector = Boundary(30)
        end if
        pwn_cr_transport = nint(transport_model_selector) == 1
        free_outer_escape = nint(escape_mode_selector) == 1
    end block
    magnetic_decay_active = magnetic_decay_alpha_t < zero .and. magnetic_decay_t0_s > zero .and. &
                            Epsilon_b_floor > zero .and. Epsilon_b_floor < Epsilon_b
    four_velocity_coord = .not. use_charint_transport .and. .not. pwn_cr_transport

    ln10 = dlog(ten)
    coord_scale = electron_four_velocity_grid_gamma_scale*electron_four_velocity_grid_gamma_scale-one
    profile_enabled = .false.
    emit_full_spectrum = emit_full_chi_spectrum /= 0
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
    substep_limit = max(1, substep_max)
    allocate(V_cool(Num_nu_cool), P_emit_cool(Num_nu_cool), P_emit_shell(Num_nu), Tau_shell(Num_nu), &
             P_hist(Num_nu, Num_chi, Num_R), &
             Seed_hist(Num_nu, Num_chi, Num_R), Tau_hist(Num_nu, Num_chi, Num_R), &
             P_hist_cool(Num_nu_cool, Num_chi, Num_R), Seed_hist_cool(Num_nu_cool, Num_chi, Num_R), &
             Tau_hist_cool(Num_nu_cool, Num_chi, Num_R), &
             Tau_pair_hist_cool(Num_nu_cool, Num_chi, Num_R), Tau_prop_hist_cool(Num_nu_cool, Num_chi, Num_R), &
             P_eff_cool_chi(Num_nu_cool, Num_chi), Seed_eff_cool_chi(Num_nu_cool, Num_chi), &
             P_stream_cool(Num_nu_cool, Num_chi), Seed_stream_cool(Num_nu_cool, Num_chi))
    call build_reduced_log_grid(Num_nu,V_seed,Num_nu_cool,V_cool)

    P_syn          = zero
    Seed_syn       = zero
    P_syn_chi      = zero
    Seed_syn_chi   = zero
    Tau_syn_chi    = zero
    chi_radius     = zero
    chi_gamma_bulk = one
    chi_weight_out = zero
    B_chi_out      = zero
    V_m            = zero
    V_c            = zero
    V_a            = zero
    dN_gam_e       = zero
    dN_gam_e_total = zero
    U_log          = zero
    dF1            = zero
    P_stream_cool  = zero
    Seed_stream_cool = zero

    call compute_q_geometry(Num_chi, dq, q_face, q_grid)
    call integrate_downstream_proper_time(Num_R,R,R_Gamma,proper_time_arr)

    call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)

    DB            = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    DB_min        = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max     = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    Gam_e_max_max = 3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam      = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-one)
    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c       = 7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)

    if (four_velocity_coord) then
        call electron_build_four_velocity_grid(Num_gam_e,one,electron_exp_tail_grid_factor*Gam_e_max_max, &
                                               electron_four_velocity_grid_gamma_scale,gam_e,coord_edge_E,x_edge_E)
        call electron_initial_powerlaw_exp_cutoff_coord_edges(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                                              Num_gam_e,coord_edge_E,coord_scale,dN_init)
    else
        call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                          electron_initial_grid_log_edges,gam_e,dN_init,x_edge_E)
        coord_edge_E = x_edge_E
    end if
    d_x_E = dlog10(gam_e(2)/gam_e(1))
    dN_init_log = dN_init

    U_log(:,1) = dN_init_log / dq

    do I_chi = 1, Num_chi
        if (four_velocity_coord) then
            call electron_shell_dcoord_to_dndgamma_exp_centers(Num_gam_e,coord_edge_E,coord_scale,gam_e, &
                                                               U_log(:,I_chi),dN_gam_e(:,I_chi,1))
        else
            call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge_E,gam_e,U_log(:,I_chi),dN_gam_e(:,I_chi,1))
        end if
    end do
    dN_gam_e_total(:,1) = zero
    do I_chi = 1, Num_chi
        dN_gam_e_total(:,1) = dN_gam_e_total(:,1) + dN_gam_e(:,I_chi,1)*dq
    end do
    call compute_downstream_comoving_grid(Num_chi,k_medium,R(1),R_Gamma(1),q_face,q_grid, &
                                          x_face_hist(:,1),x_comov_face_hist(:,1),x_comov_hist(:,1),dx_comov_hist(:,1), &
                                          radius_cell_hist(:,1),gamma_cell_hist(:,1),beta_hist(:,1),beta_rel_sh_chi)
    call update_epsilon_b_db_chi(x_comov_hist(:,1), R_Gamma(1), dNe, beta_rel_sh_chi)
    do I_chi = 1, Num_chi
        B_chi_out(I_chi,1) = DB_chi(I_chi)
    end do
    if (profile_enabled) call cpu_time(t_start)
    if (emit_full_spectrum) then
        call get_syn_chi_batch_state(R(1),Num_gam_e,Num_nu,Num_chi,gam_e,dN_gam_e(:,:,1),V_seed,DB_chi, &
                                     P_local,P_hist(:,:,1),Seed_hist(:,:,1),Tau_hist(:,:,1))
        do I_chi = 1, Num_chi
            call project_syn_state_logbands(Num_nu,V_seed,P_hist(:,I_chi,1),Seed_hist(:,I_chi,1),Tau_hist(:,I_chi,1), &
                                            Num_nu_cool,V_cool,P_hist_cool(:,I_chi,1), &
                                            Seed_hist_cool(:,I_chi,1),Tau_hist_cool(:,I_chi,1))
            Tau_pair_hist_cool(:,I_chi,1) = zero
            call radiation_pair_tau_headon_segment(V_cool,Num_nu_cool,Seed_hist_cool(:,I_chi,1),dx_comov_hist(I_chi,1), &
                                                   Tau_pair_hist_cool(:,I_chi,1))
            Tau_prop_hist_cool(:,I_chi,1) = Tau_hist_cool(:,I_chi,1) + Tau_pair_hist_cool(:,I_chi,1)
        end do
        syn_state_calls = syn_state_calls + Num_chi
    else
        do I_chi = 1, Num_chi
            dN_cell = dN_gam_e(:,I_chi,1)
            call get_syn_state(R(1),DB_chi(I_chi),Num_gam_e,Num_nu_cool,n_threads,gam_e,dN_cell,V_cool, &
                               P_emit_cool,P_hist_cool(:,I_chi,1),Seed_hist_cool(:,I_chi,1),Tau_hist_cool(:,I_chi,1))
            Tau_pair_hist_cool(:,I_chi,1) = zero
            call radiation_pair_tau_headon_segment(V_cool,Num_nu_cool,Seed_hist_cool(:,I_chi,1),dx_comov_hist(I_chi,1), &
                                                   Tau_pair_hist_cool(:,I_chi,1))
            Tau_prop_hist_cool(:,I_chi,1) = Tau_hist_cool(:,I_chi,1) + Tau_pair_hist_cool(:,I_chi,1)
        end do
        syn_state_calls = syn_state_calls + Num_chi
    end if
    if (profile_enabled) then
        call cpu_time(t_stop)
        t_syn_state = t_syn_state + (t_stop-t_start)
    end if

    do I_tobs = 2, Num_R
        R_loc       = R(I_tobs-1)
        R_Gamma_loc = (R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two

        call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)

        call update_shock_cooling_scales(R_Gamma_loc, R_Tobs(I_tobs), beta_sh, beta_2, beta_2_sh)

        if (emit_full_spectrum) then
            call get_syn_selected_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e, &
                                        dN_gam_e_total(:,I_tobs-1),V_seed,P_emit_shell, &
                                        P_syn(:,I_tobs),Seed_syn(:,I_tobs),Tau_shell)
        end if

        if (profile_enabled) call cpu_time(t_start)
        P_eff_cool_chi = P_hist_cool(:,:,I_tobs-1) + P_stream_cool
        Seed_eff_cool_chi = Seed_hist_cool(:,:,I_tobs-1) + Seed_stream_cool
        history_calls = history_calls + 1
        if (profile_enabled) then
            call cpu_time(t_stop)
            t_hist_accum = t_hist_accum + (t_stop-t_start)
        end if
        call compute_q_cell_geometry(Num_chi,k_medium,R_loc,R_Gamma_loc,q_grid, &
                                     radius_cell_chi,gamma_cell_chi,beta_cell_chi,beta_rel_sh_chi)
        call update_epsilon_b_db_chi(x_comov_hist(:,I_tobs-1), R_Gamma_loc, dNe, beta_rel_sh_chi)

        if (profile_enabled) call cpu_time(t_start)
        call prepare_forward_cooling_aux_batch(index_Y,Num_gam_e,Num_nu_cool,Num_chi,n_threads,gam_e,V_cool, &
                                               P_eff_cool_chi,Seed_eff_cool_chi,cooling_aux_chi)
        if (profile_enabled) then
            call cpu_time(t_stop)
            t_prepare_aux = t_prepare_aux + (t_stop-t_start)
        end if
        prepare_aux_calls = prepare_aux_calls + 1
        if (profile_enabled) call cpu_time(t_start)
        call assemble_cooling_chi(R_loc, R_Gamma_loc, beta_sh, Seed_eff_cool_chi, R_Tobs(I_tobs))
        if (profile_enabled) then
            call cpu_time(t_stop)
            t_cooling = t_cooling + (t_stop-t_start)
        end if
        shell_cooling_calls = shell_cooling_calls + 1
        call compute_vm_vc_va(R_Gamma_loc, beta_sh, I_tobs-1, R_loc)
        block
            real(8) :: dDD, dDR, dDR_try, dDR_xi, dDR_q, max_xi_coeff
            real(8) :: chi_peak
            integer :: L, L1, src_lo, src_hi, active_hi, active_chi_hi

            call electron_source_bounds(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,src_lo,src_hi)
            if (four_velocity_coord) then
                call electron_build_source_term_exp_cutoff_coord_edges(Num_gam_e,coord_edge_E,coord_scale, &
                                                                       Gam_e_m,Gam_e_max,one,p,dF1)
            else
                call electron_build_source_term_exp_cutoff_edges(Num_gam_e,x_edge_E,Gam_e_m,Gam_e_max,one,p,dF1)
            end if
            shell_population = sum(U_log, dim=2)
            chi_population = sum(U_log, dim=1)
            chi_peak = maxval(chi_population)
            active_hi = electron_active_gamma_hi(Num_gam_e,dF1,shell_population,src_lo,src_hi,chi_peak)
            if (pwn_cr_transport) then
                call compute_q_divergence(Num_chi,k_medium,R_loc,R_Gamma_loc,beta_sh,q_grid,adiabatic_log_coeff_chi)
            else
                adiabatic_log_coeff_chi = one/(R_loc*ln10)
            end if
            max_xi_coeff = electron_max_xi_coeff_chi(Num_gam_e,Num_chi,dEL_mean_chi, &
                                                     adiabatic_log_coeff_chi,chi_population,chi_peak,active_hi)

            dDR_xi = huge(one)
            if (max_xi_coeff > zero) dDR_xi = 4d0*d_x_E/max_xi_coeff
            active_chi_hi = electron_active_chi_hi(Num_chi,chi_population,chi_peak)
            dDR_q = huge(one)
            if (use_charint_transport) then
                dDR_q = compute_q_step_limit(Num_chi,k_medium,R_loc,dq,q_face,4d0)
            end if

            dDD     = R(I_tobs)-R(I_tobs-1)
            dDR_try = min(dDD, min(dDR_xi, dDR_q))
            L1      = max(1, min(substep_limit, ceiling(dDD/dDR_try)))
            if (.not. use_charint_transport .and. .not. pwn_cr_transport) then
                call transport_step_fullhide(R(I_tobs-1), R(I_tobs), R_Gamma(I_tobs-1), R_Gamma(I_tobs), &
                                              dDD, active_hi, active_chi_hi)
            else
                dDR     = dDD/dble(L1)
                total_substeps = total_substeps + L1
                max_shell_substeps = max(max_shell_substeps, L1)

                do L = 1, L1
                    block
                        real(8) :: frac_sub, R_sub, Gamma_sh_sub, Gam_e_m_p

                        frac_sub = (dble(L)-0.5d0)/dble(L1)
                        R_sub = R(I_tobs-1) + frac_sub*dDD
                        Gamma_sh_sub = (one-frac_sub)*R_Gamma(I_tobs-1) + frac_sub*R_Gamma(I_tobs)

                        call dynamics_external_density_profile(A_star,dNe_ISM,R_sub,R0,1,R_tr,f_jump,f_wide,dNe)

                        call update_shock_cooling_scales(Gamma_sh_sub, R_Tobs(I_tobs), beta_sh, beta_2, beta_2_sh)
                        if (pwn_cr_transport) then
                            call compute_q_divergence(Num_chi,k_medium,R_sub,Gamma_sh_sub,beta_sh,q_grid,adiabatic_log_coeff_chi)
                        else
                            adiabatic_log_coeff_chi = one/(R_sub*ln10)
                        end if

                        Gam_e_m_p = (one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
                        call electron_injection_prefactor(R_sub,dDR,dNe,f_e,Gam_e_m_p,Q)
                        Q_rate = Q/dDR
                        call electron_build_source_term_exp_cutoff_edges(Num_gam_e,x_edge_E,Gam_e_m,Gam_e_max,Q_rate,p,dF1)
                        source_q1 = dF1/dq

                        if (use_charint_transport) then
                            if (profile_enabled) call cpu_time(t_start)
                            call advance_q_advection_charint(U_log, Num_gam_e, Num_chi, active_hi, dq, q_face, &
                                                             k_medium, R_sub, zero*source_q1, dDR)
                            call advance_q_diffusion_implicit(U_log, Num_gam_e, Num_chi, active_hi, dq, q_face, &
                                                              k_medium, R_sub, Gamma_sh_sub, beta_sh, &
                                                              kappa2_chi, dDR, n_threads)
                            eta_calls = eta_calls + 1
                            if (profile_enabled) then
                                call cpu_time(t_stop)
                                t_eta = t_eta + (t_stop-t_start)
                            end if

                            if (profile_enabled) call cpu_time(t_start)
                            call advance_energy_loggamma_chi_charint(U_log, Num_gam_e, Num_chi, gam_e, DB_chi, &
                                                                     dEl_chi, R_sub, Gamma_sh_sub, beta_sh, index_Y, &
                                                                     dDR, active_chi_hi, n_threads, source_q1)
                            xi_calls = xi_calls + 1
                            if (profile_enabled) then
                                call cpu_time(t_stop)
                                t_xi = t_xi + (t_stop-t_start)
                            end if
                        else
                            if (profile_enabled) call cpu_time(t_start)
                            if (pwn_cr_transport) then
                                call advance_q_pwncr_implicit(U_log, Num_gam_e, Num_chi, active_hi, dq, q_face, &
                                                              k_medium, R_sub, Gamma_sh_sub, beta_sh, kappa2_chi, &
                                                              source_q1, dDR, free_outer_escape, n_threads)
                            else
                                call advance_q_implicit(U_log, Num_gam_e, Num_chi, active_hi, dq, q_face, &
                                                        k_medium, R_sub, Gamma_sh_sub, beta_sh, kappa2_chi, &
                                                        source_q1, dDR, n_threads)
                            end if
                            eta_calls = eta_calls + 1
                            if (profile_enabled) then
                                call cpu_time(t_stop)
                                t_eta = t_eta + (t_stop-t_start)
                            end if

                            if (profile_enabled) call cpu_time(t_start)
                            if (pwn_cr_transport) then
                                if (stochastic_accel_norm > zero) then
                                    call advance_energy_stochastic_loggamma_chi(U_log, Num_gam_e, Num_chi, &
                                                                                stochastic_accel_norm, R_sub, &
                                                                                d_x_E, 0.5d0*dDR, n_threads)
                                end if
                                call advance_energy_loggamma_chi_pwncr(U_log, Num_gam_e, Num_chi, dEL_mean_chi, &
                                                                       adiabatic_log_coeff_chi, d_x_E, dDR, n_threads)
                                if (stochastic_accel_norm > zero) then
                                    call advance_energy_stochastic_loggamma_chi(U_log, Num_gam_e, Num_chi, &
                                                                                stochastic_accel_norm, R_sub, &
                                                                                d_x_E, 0.5d0*dDR, n_threads)
                                end if
                            else
                                call advance_energy_loggamma_chi(U_log, Num_gam_e, Num_chi, dEL_mean_chi, &
                                                                 R_sub, d_x_E, dDR, n_threads)
                            end if
                            xi_calls = xi_calls + 1
                            if (profile_enabled) then
                                call cpu_time(t_stop)
                                t_xi = t_xi + (t_stop-t_start)
                            end if
                        end if
                    end block
                end do
            end if
        end block

        call dynamics_external_density_profile(A_star,dNe_ISM,R(I_tobs),R0,1,R_tr,f_jump,f_wide,dNe)
        call compute_downstream_comoving_grid(Num_chi,k_medium,R(I_tobs),R_Gamma(I_tobs),q_face,q_grid, &
                                              x_face_hist(:,I_tobs),x_comov_face_hist(:,I_tobs),x_comov_hist(:,I_tobs), &
                                              dx_comov_hist(:,I_tobs),radius_cell_hist(:,I_tobs),gamma_cell_hist(:,I_tobs), &
                                              beta_hist(:,I_tobs),beta_rel_sh_chi)
        call update_epsilon_b_db_chi(x_comov_hist(:,I_tobs), R_Gamma(I_tobs), dNe, beta_rel_sh_chi)
        if (profile_enabled) call cpu_time(t_start)
        !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. Num_chi*Num_nu >= 128) schedule(static) &
        !$OMP& private(I_chi)
        do I_chi = 1, Num_chi
            if (four_velocity_coord) then
                call electron_shell_dcoord_to_dndgamma_exp_centers(Num_gam_e,coord_edge_E,coord_scale,gam_e, &
                                                                   U_log(:,I_chi),dN_gam_e(:,I_chi,I_tobs))
            else
                call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge_E,gam_e,U_log(:,I_chi),dN_gam_e(:,I_chi,I_tobs))
            end if
            B_chi_out(I_chi,I_tobs) = DB_chi(I_chi)
            if (.not. emit_full_spectrum) then
                block
                    real(8) :: P_emit_tmp(Num_nu_cool)
                    call get_syn_state(R(I_tobs),DB_chi(I_chi),Num_gam_e,Num_nu_cool,1, &
                                       gam_e,dN_gam_e(:,I_chi,I_tobs),V_cool,P_emit_tmp, &
                                       P_hist_cool(:,I_chi,I_tobs),Seed_hist_cool(:,I_chi,I_tobs), &
                                       Tau_hist_cool(:,I_chi,I_tobs))
                end block
                Tau_pair_hist_cool(:,I_chi,I_tobs) = zero
                call radiation_pair_tau_headon_segment(V_cool,Num_nu_cool,Seed_hist_cool(:,I_chi,I_tobs), &
                                                       dx_comov_hist(I_chi,I_tobs),Tau_pair_hist_cool(:,I_chi,I_tobs))
                Tau_prop_hist_cool(:,I_chi,I_tobs) = Tau_hist_cool(:,I_chi,I_tobs) + Tau_pair_hist_cool(:,I_chi,I_tobs)
            end if
        end do
        !$OMP END PARALLEL DO
        if (emit_full_spectrum) then
            call get_syn_chi_batch_state(R(I_tobs),Num_gam_e,Num_nu,Num_chi,gam_e,dN_gam_e(:,:,I_tobs),V_seed,DB_chi, &
                                         P_local,P_hist(:,:,I_tobs),Seed_hist(:,:,I_tobs),Tau_hist(:,:,I_tobs))
            do I_chi = 1, Num_chi
                call project_syn_state_logbands(Num_nu,V_seed,P_hist(:,I_chi,I_tobs), &
                                                Seed_hist(:,I_chi,I_tobs),Tau_hist(:,I_chi,I_tobs), &
                                                Num_nu_cool,V_cool,P_hist_cool(:,I_chi,I_tobs), &
                                                Seed_hist_cool(:,I_chi,I_tobs),Tau_hist_cool(:,I_chi,I_tobs))
                Tau_pair_hist_cool(:,I_chi,I_tobs) = zero
                call radiation_pair_tau_headon_segment(V_cool,Num_nu_cool,Seed_hist_cool(:,I_chi,I_tobs), &
                                                       dx_comov_hist(I_chi,I_tobs),Tau_pair_hist_cool(:,I_chi,I_tobs))
                Tau_prop_hist_cool(:,I_chi,I_tobs) = Tau_hist_cool(:,I_chi,I_tobs) + Tau_pair_hist_cool(:,I_chi,I_tobs)
            end do
        end if
        if (profile_enabled) then
            call cpu_time(t_stop)
            t_syn_state = t_syn_state + (t_stop-t_start)
        end if
        syn_state_calls = syn_state_calls + Num_chi

        if (profile_enabled) call cpu_time(t_start)
        call advance_comoving_history_stream(I_tobs-1,I_tobs,Num_R,Num_chi,Num_nu_cool,proper_time_arr,V_cool, &
                                             x_comov_face_hist,x_comov_hist,dx_comov_hist,beta_hist, &
                                             Tau_prop_hist_cool,P_hist_cool,Seed_hist_cool,P_stream_cool,Seed_stream_cool)
        if (profile_enabled) then
            call cpu_time(t_stop)
            t_hist_accum = t_hist_accum + (t_stop-t_start)
        end if

        dN_gam_e_total(:, I_tobs) = zero
        do I_chi = 1, Num_chi
            dN_gam_e_total(:, I_tobs) = dN_gam_e_total(:, I_tobs) + dN_gam_e(:,I_chi,I_tobs)*dq
        end do
    end do

    R_loc = R(Num_R)
    R_Gamma_loc = R_Gamma(Num_R)
    call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
    call update_shock_cooling_scales(R_Gamma_loc, R_Tobs(Num_R), beta_sh, beta_2, beta_2_sh)

    if (emit_full_spectrum) then
        call get_syn_selected_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e, &
                                    dN_gam_e_total(:,Num_R),V_seed,P_emit_shell,P_syn(:,Num_R), &
                                    Seed_syn(:,Num_R),Tau_shell)
    end if
    if (profile_enabled) call cpu_time(t_start)
    P_eff_cool_chi = P_hist_cool(:,:,Num_R) + P_stream_cool
    Seed_eff_cool_chi = Seed_hist_cool(:,:,Num_R) + Seed_stream_cool
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
    call compute_q_cell_geometry(Num_chi,k_medium,R_loc,R_Gamma_loc,q_grid, &
                                 radius_cell_chi,gamma_cell_chi,beta_cell_chi,beta_rel_sh_chi)
    call update_epsilon_b_db_chi(x_comov_hist(:,Num_R), R_Gamma_loc, dNe, beta_rel_sh_chi)
    if (profile_enabled) call cpu_time(t_start)
    call assemble_cooling_chi(R_loc, R_Gamma_loc, beta_sh, Seed_eff_cool_chi, R_Tobs(Num_R))
    if (profile_enabled) then
        call cpu_time(t_stop)
        t_cooling = t_cooling + (t_stop-t_start)
    end if
    call compute_vm_vc_va(R_Gamma_loc, beta_sh, Num_R, R_loc)
    if (emit_full_spectrum) then
        do I_tobs = 1, Num_R
            call project_q_projection_shell(Num_nu,Num_chi,dq,P_hist(:,:,I_tobs),Seed_hist(:,:,I_tobs), &
                                            Tau_hist(:,:,I_tobs),radius_cell_hist(:,I_tobs),gamma_cell_hist(:,I_tobs), &
                                            P_syn_chi(:,:,I_tobs),Seed_syn_chi(:,:,I_tobs),Tau_syn_chi(:,:,I_tobs), &
                                            chi_radius(:,I_tobs),chi_gamma_bulk(:,I_tobs),chi_weight_out(:,I_tobs))
        end do
    else
        do I_tobs = 1, Num_R
            call project_q_projection_geometry(Num_chi,dq,radius_cell_hist(:,I_tobs),gamma_cell_hist(:,I_tobs), &
                                               chi_radius(:,I_tobs),chi_gamma_bulk(:,I_tobs),chi_weight_out(:,I_tobs))
        end do
    end if
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
    deallocate(dEl, dEL_mean, dEL_mean_shell, dN_init, dN_init_log, dF1, &
               shell_population, chi_population, U_log, U_shell, source_q1, &
               V_m_chi, V_c_chi, V_a_chi, chi_weight, Epsilon_b_chi, DB_chi, t_decay_chi, &
               q_grid, q_face, proper_time_arr, &
               x_face_hist, x_comov_face_hist, &
               x_comov_hist, dx_comov_hist, radius_cell_hist, gamma_cell_hist, beta_hist, &
               radius_cell_chi, gamma_cell_chi, beta_cell_chi, beta_rel_sh_chi, &
               dN_cell, P_local, V_cool, P_emit_cool, P_emit_shell, Tau_shell, P_hist, Seed_hist, Tau_hist, &
               P_hist_cool, Seed_hist_cool, Tau_hist_cool, Tau_pair_hist_cool, Tau_prop_hist_cool, &
               P_eff_cool_chi, Seed_eff_cool_chi, P_stream_cool, Seed_stream_cool, &
               cooling_aux_chi, dEl_chi, dEL_mean_chi, adiabatic_log_coeff_chi, kappa2_chi)

contains

subroutine transport_step_fullhide(R_prev, R_curr, Gamma_prev, Gamma_curr, dDR_step, &
                                    active_hi, active_chi_hi)
    real(8), intent(in) :: R_prev, R_curr, Gamma_prev, Gamma_curr, dDR_step
    integer, intent(in) :: active_hi, active_chi_hi
    real(8) :: R_sub, R_eff, half_dR, Gamma_sh_sub, Gam_e_m_p
    real(8) :: adiabatic_integral, face_coord, face_jac
    real(8) :: DB_loc, Gam_e_max_loc, Gam_e_m_loc, temp_gam_loc
    real(8) :: Q_rate_loc, source_q1_loc(Num_gam_e), dF1_zero(Num_gam_e)
    real(8) :: coord_face_step_loc(Num_gam_e-1)
    integer :: I

    R_sub = 0.5d0*(R_prev+R_curr)
    Gamma_sh_sub = 0.5d0*(Gamma_prev+Gamma_curr)
    half_dR = 0.5d0*dDR_step
    if (R_curr > R_prev) then
        R_eff = dDR_step/dlog(R_curr/R_prev)
    else
        R_eff = R_sub
    end if

    call get_shock_transport_state(Gamma_sh_sub, beta_sh, beta_2, beta_2_sh)
    if (pwn_cr_transport) then
        call compute_q_divergence(Num_chi,k_medium,R_sub,Gamma_sh_sub,beta_sh,q_grid,adiabatic_log_coeff_chi)
    else
        adiabatic_log_coeff_chi = one/(R_sub*ln10)
    end if

    call dynamics_external_density_profile(A_star,dNe_ISM,R_sub,R0,1,R_tr,f_jump,f_wide,dNe)
    DB_loc = 0.39d0*dsqrt(Epsilon_b*dNe*(Gamma_sh_sub*(Gamma_sh_sub-one)))
    Gam_e_max_loc = 3d0*Para_m_energy/dsqrt(8d0*DB_loc*Para_e**3)
    temp_gam_loc = Epsilon_e/f_e*para_m_p/para_m_e*(Gamma_sh_sub-one)
    call electron_gamma_m_exact(p,temp_gam_loc,Gam_e_max_loc,Gam_e_m_loc)
    Gam_e_m_p = (one-p)/(Gam_e_max_loc**(one-p)-Gam_e_m_loc**(one-p))
    call electron_injection_prefactor(R_prev,dDR_step,dNe,f_e,Gam_e_m_p,Q_rate_loc)
    call electron_build_source_term_exp_cutoff_coord_edges(Num_gam_e,coord_edge_E,coord_scale, &
                                                           Gam_e_m_loc,Gam_e_max_loc,Q_rate_loc*dDR_step/dq,p,source_q1_loc)

    U_shell = U_log
    dF1_zero = zero
    if (profile_enabled) call cpu_time(t_start)
    call advance_q_advection_charint(U_shell, Num_gam_e, Num_chi, active_hi, dq, q_face, &
                                     k_medium, R_eff, dF1_zero, half_dR)
    call advance_q_diffusion_implicit(U_shell, Num_gam_e, Num_chi, active_hi, dq, q_face, &
                                      k_medium, R_eff, Gamma_sh_sub, beta_sh, &
                                      kappa2_chi, half_dR, n_threads)
    eta_calls = eta_calls + 1
    if (profile_enabled) then
        call cpu_time(t_stop)
        t_eta = t_eta + (t_stop-t_start)
    end if

    if (profile_enabled) call cpu_time(t_start)
    adiabatic_integral = dlog(R_curr/R_prev)
    do I_chi = 1, active_chi_hi
        do I = 1, Num_gam_e-1
            face_coord = coord_edge_E(I+1)
            face_jac = ln10*electron_dxgamma_dcoord(electron_coord_log_four_velocity_sq,coord_scale,face_coord)
            coord_face_step_loc(I) = &
                (dDR_step*(dEl_chi(I,I_chi)+dEl_chi(I+1,I_chi))/two+adiabatic_integral)/face_jac
        end do
        if (I_chi == 1) then
            call electron_shell_flux_split_coord_sequence(Num_gam_e,coord_edge_E,coord_face_step_loc, &
                                                          source_q1_loc,U_shell(:,I_chi),U_log(:,I_chi))
        else
            call electron_shell_flux_split_coord_sequence(Num_gam_e,coord_edge_E,coord_face_step_loc, &
                                                          dF1_zero,U_shell(:,I_chi),U_log(:,I_chi))
        end if
    end do
    if (active_chi_hi < Num_chi) U_log(:,active_chi_hi+1:Num_chi) = U_shell(:,active_chi_hi+1:Num_chi)
    xi_calls = xi_calls + 1
    if (profile_enabled) then
        call cpu_time(t_stop)
        t_xi = t_xi + (t_stop-t_start)
    end if

    if (profile_enabled) call cpu_time(t_start)
    U_shell = U_log
    call advance_q_advection_charint(U_shell, Num_gam_e, Num_chi, active_hi, dq, q_face, &
                                     k_medium, R_eff, dF1_zero, half_dR)
    call advance_q_diffusion_implicit(U_shell, Num_gam_e, Num_chi, active_hi, dq, q_face, &
                                      k_medium, R_eff, Gamma_sh_sub, beta_sh, &
                                      kappa2_chi, half_dR, n_threads)
    eta_calls = eta_calls + 1
    if (profile_enabled) then
        call cpu_time(t_stop)
        t_eta = t_eta + (t_stop-t_start)
    end if
    U_log = U_shell
    total_substeps = total_substeps + 1
    max_shell_substeps = max(max_shell_substeps, 1)
end subroutine transport_step_fullhide

subroutine update_epsilon_b_db_chi(x_comov_col, Gamma_f, dNe_val, beta_rel_sh)
    real(8), intent(in) :: x_comov_col(Num_chi), Gamma_f, dNe_val, beta_rel_sh(Num_chi)
    real(8) :: Gam_m1
    integer :: I
    Gam_m1 = Gamma_f*(Gamma_f-one)
    do I = 1, Num_chi
        t_decay_chi(I) = x_comov_col(I)/max(beta_rel_sh(I)*para_c, tiny(one))
        if (magnetic_decay_active) then
            Epsilon_b_chi(I) = Epsilon_b_floor + (Epsilon_b-Epsilon_b_floor) * &
                               (one+t_decay_chi(I)/magnetic_decay_t0_s)**magnetic_decay_alpha_t
        else
            Epsilon_b_chi(I) = Epsilon_b
        end if
        DB_chi(I) = 0.39d0*dsqrt(Epsilon_b_chi(I)*dNe_val*Gam_m1)
    end do
end subroutine update_epsilon_b_db_chi

subroutine update_shock_cooling_scales(Gamma_f, R_Tobs, beta_sh, beta_2, beta_2_sh)
    real(8), intent(in) :: Gamma_f, R_Tobs
    real(8), intent(out) :: beta_sh, beta_2, beta_2_sh
    DB = 0.39d0*dsqrt(Epsilon_b*dNe*(Gamma_f*(Gamma_f-one)))
    Gam_e_max = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(Gamma_f-one)
    call electron_gamma_m_exact(p, temp_gam, Gam_e_max, Gam_e_m)
    Gam_e_c = 7.7d8*(one+z)/Gamma_f/DB**2/R_Tobs
    call get_shock_transport_state(Gamma_f, beta_sh, beta_2, beta_2_sh)
end subroutine update_shock_cooling_scales

subroutine assemble_cooling_chi(Radius, Gamma_f, beta_sh, Seed_eff_chi, R_Tobs_val)
    real(8), intent(in) :: Radius, Gamma_f, beta_sh, Seed_eff_chi(Num_nu_cool, Num_chi), R_Tobs_val
    if (magnetic_decay_active) then
        do I_chi = 1, Num_chi
            Gam_e_max_cell = 3d0*Para_m_energy/dsqrt(8d0*DB_chi(I_chi)*Para_e**3)
            temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(Gamma_f-one)
            call electron_gamma_m_exact(p, temp_gam, Gam_e_max_cell, Gam_e_m_cell)
            Gam_e_c_cell = 7.7d8*(one+z)/Gamma_f/DB_chi(I_chi)**2/R_Tobs_val
            call assemble_forward_cooling_split(index_Y, Epsilon_e, Epsilon_b_chi(I_chi), p, DB_chi(I_chi), &
                                                Gam_e_m_cell, Gam_e_c_cell, &
                                                Gam_e_max_cell, Radius, Gamma_f, beta_sh, dNe, Num_gam_e, &
                                                Num_nu_cool, n_threads, gam_e, V_cool, Seed_eff_chi(:,I_chi), &
                                                cooling_aux_chi(:,I_chi), dEl_chi(:,I_chi))
            dEL_mean_chi(:,I_chi)=(dEl_chi(2:Num_gam_e,I_chi)+dEl_chi(1:Num_gam_e-1,I_chi))/two/dlog(ten)
            kappa2_chi(:,I_chi) = gam_e*Para_m_energy*para_c/(3d0*Para_e*DB_chi(I_chi))
        end do
    else
        call assemble_forward_cooling_split_batch(index_Y, Epsilon_e, Epsilon_b, p, DB, Gam_e_m, Gam_e_c, &
                                                  Gam_e_max, Radius, Gamma_f, &
                                                  beta_sh, dNe, Num_gam_e, Num_nu_cool, Num_chi, n_threads, gam_e, V_cool, &
                                                  Seed_eff_chi, cooling_aux_chi, dEl_chi)
        do I_chi = 1, Num_chi
            dEL_mean_chi(:,I_chi)=(dEl_chi(2:Num_gam_e,I_chi)+dEl_chi(1:Num_gam_e-1,I_chi))/two/dlog(ten)
        end do
        do I_chi = 1, Num_chi
            kappa2_chi(:,I_chi) = gam_e*Para_m_energy*para_c/(3d0*Para_e*DB)
        end do
    end if
end subroutine assemble_cooling_chi

subroutine compute_vm_vc_va(Gamma_f, beta_sh, I_tobs_out, Radius_val)
    real(8), intent(in) :: Gamma_f, beta_sh
    integer, intent(in) :: I_tobs_out
    real(8), intent(in) :: Radius_val
    if (emit_full_spectrum) then
        call get_nu_a_2d_cell_path(Num_nu, Num_chi, V_seed, Tau_hist(:,:,I_tobs_out), V_a_chi)
    else
        call get_nu_a_2d_cell_path(Num_nu_cool, Num_chi, V_cool, Tau_hist_cool(:,:,I_tobs_out), V_a_chi)
    end if
    temp = zero
    Q = zero
    do I_chi = 1, Num_chi
        chi_weight(I_chi) = sum(U_log(:,I_chi))
        if (chi_weight(I_chi) > zero) then
            if (emit_full_spectrum) then
                V_m_chi(I_chi) = electron_logparabola_peak_frequency(Num_nu, V_seed, P_hist(:,I_chi,I_tobs_out))
            else
                V_m_chi(I_chi) = electron_logparabola_peak_frequency(Num_nu_cool, V_cool, P_hist_cool(:,I_chi,I_tobs_out))
            end if
            dEL_mean_shell = dEL_mean_chi(:,I_chi)
            call electron_gamma_c_from_loss_mean(Num_gam_e, gam_e, dEL_mean_shell, Radius_val, Gam_e_c_diag)
            V_c_chi(I_chi) = 4.2d6*DB_chi(I_chi)*Gam_e_c_diag*Gam_e_c_diag
            temp = temp + chi_weight(I_chi)*dlog(V_m_chi(I_chi))
            Q = Q + chi_weight(I_chi)
        end if
    end do
    V_m(I_tobs_out) = dexp(temp/Q)/(Gamma_f*(1d0-beta_sh)*(one+z))
    temp = zero
    do I_chi = 1, Num_chi
        if (chi_weight(I_chi) > zero) temp = temp + chi_weight(I_chi)*dlog(V_c_chi(I_chi))
    end do
    V_c(I_tobs_out) = dexp(temp/Q)/(Gamma_f*(1d0-beta_sh)*(one+z))
    V_a(I_tobs_out) = V_a_chi(Num_chi)/(Gamma_f*(1d0-beta_sh)*(one+z))
end subroutine compute_vm_vc_va

end subroutine fs_electron_transport_2d_core

subroutine reduce_syn_shell_from_q(Num_nu,Num_chi,dq,P_q,Seed_q,P_shell,Seed_shell)
    use constants
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: Num_nu,Num_chi
    real(8), intent(in) :: dq,P_q(Num_nu,Num_chi),Seed_q(Num_nu,Num_chi)
    real(8), intent(out) :: P_shell(Num_nu),Seed_shell(Num_nu)
    integer :: I_chi

    P_shell = zero
    Seed_shell = zero
    do I_chi = 1, Num_chi
        P_shell = P_shell + dq*P_q(:,I_chi)
        Seed_shell = Seed_shell + dq*Seed_q(:,I_chi)
    end do
end subroutine reduce_syn_shell_from_q

subroutine project_q_projection_shell(Num_nu,Num_chi,dq,P_src,Seed_src,Tau_src,radius_cell,gamma_cell, &
                                      P_dst,Seed_dst,Tau_dst,chi_radius_dst,chi_gamma_dst,chi_weight_dst)
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: Num_nu,Num_chi
    real(8), intent(in) :: dq
    real(8), intent(in) :: P_src(Num_nu,Num_chi),Seed_src(Num_nu,Num_chi),Tau_src(Num_nu,Num_chi)
    real(8), intent(in) :: radius_cell(Num_chi),gamma_cell(Num_chi)
    real(8), intent(out) :: P_dst(Num_nu,Num_chi),Seed_dst(Num_nu,Num_chi),Tau_dst(Num_nu,Num_chi)
    real(8), intent(out) :: chi_radius_dst(Num_chi),chi_gamma_dst(Num_chi),chi_weight_dst(Num_chi)
    integer :: I_chi

    P_dst = P_src
    Seed_dst = Seed_src
    Tau_dst = Tau_src
    do I_chi = 1, Num_chi
        chi_radius_dst(I_chi) = radius_cell(I_chi)
        chi_gamma_dst(I_chi) = gamma_cell(I_chi)
        chi_weight_dst(I_chi) = dq
    end do
end subroutine project_q_projection_shell

subroutine project_q_projection_geometry(Num_chi,dq,radius_cell,gamma_cell,chi_radius_dst,chi_gamma_dst,chi_weight_dst)
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: Num_chi
    real(8), intent(in) :: dq
    real(8), intent(in) :: radius_cell(Num_chi),gamma_cell(Num_chi)
    real(8), intent(out) :: chi_radius_dst(Num_chi),chi_gamma_dst(Num_chi),chi_weight_dst(Num_chi)
    integer :: I_chi

    do I_chi = 1, Num_chi
        chi_radius_dst(I_chi) = radius_cell(I_chi)
        chi_gamma_dst(I_chi) = gamma_cell(I_chi)
        chi_weight_dst(I_chi) = dq
    end do
end subroutine project_q_projection_geometry
