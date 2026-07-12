! 电子2D输运核心。
! 顺序: unpack config -> construct finite-q shell geometry -> initialize shock-front electrons
!       -> loop shells/substeps: shock state -> chi-local cooling/radiation history
!       -> energy advance -> q advection/diffusion -> shell/chi spectrum accumulation
!       -> return shell-level outputs plus projection geometry.  The chi_* arrays are observer
!       projection geometry, not a chi-local hadronic contract.
subroutine fs_transport_2d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,ng, &
                                         nchi,index_Y,index_syn_intger,n_threads,emit_full_chi_spectrum, &
                                         gam_e,dN_gam_e,dN_gam_e_total,P_syn,Seed_syn,V_m,V_c,V_a, &
                                         P_syn_chi,Seed_syn_chi,Tau_syn_chi,chi_radius,chi_gamma_bulk,chi_weight_out, &
                                         B_chi_out,substep_max,use_charint_transport, profile_tag)
    !$ use omp_lib
    use constants
    use dynamics_density_profile, only: density_profile, profile_count
    use electron_common, only: electron_initial_density, electron_initialize_spectrum, electron_unpack_boundary, &
                               imodelog, tail_factor, &
                               electron_gm_exact, electron_injection_prefactor, &
                               electron_gc_loss, electron_source_bounds
    use electron_injection_profiles, only: source_edges, &
                                           log_edges, &
                                           init_coord, &
                                           source_coord
    use electron_cooling_kernel, only: forward_cooling
    use electron_radiation_kernel, only: syn_state, nua_path, &
                                         reduce_grid, project_syn
    use electron_seed_history_kernel, only: integrate_proper_time, advance_history_stream
    use rad_common, only: pair_tau, syn_seed_chi, transfer_factor
    use electron_transport_2d, only: q_geometry, q_cell_geometry, shock_state, &
                                             downstream_grid, q_divergence, q_step_limit
    use electron_transport_2d, only: advance_q_implicit, advance_q_charint, &
                                             advance_q_diffusion, advance_pwncr_q
    use electron_transport_2d, only: advance_energy_chi, advance_energy_charint, &
                                             advance_pwncr_energy, advance_stoch_chi
    use electron_transport_common, only: logparabola_peak, &
                                         gamma_active_hi, chi_active_hi, &
                                         max_xi_chi, &
                                         dnx_dgamma, &
                                         flux_seq_nonuniform
    use electron_coord_common, only: build_fourvel_grid, coord_fourvel, &
                                                 dxg_dcoord, fourvel_scale
    use electron_shell_transport, only: coord_to_dgamma
    implicit real(8)(a-h,o-z)

    integer, intent(in) :: n,Num_nu,Num_R,ng,nchi,index_Y,index_syn_intger,n_threads,substep_max
    integer, intent(in) :: emit_full_chi_spectrum
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    logical, intent(in) :: use_charint_transport
    character(len=*), intent(in) :: profile_tag
    real(8), intent(out), dimension(ng,nchi,Num_R) :: dN_gam_e
    real(8), intent(out), dimension(ng,Num_R) :: dN_gam_e_total
    real(8), intent(out), dimension(ng) :: gam_e
    real(8), intent(out), dimension(Num_nu,Num_R) :: P_syn, Seed_syn
    real(8), intent(out), dimension(Num_R) :: V_m,V_c,V_a
    real(8), intent(out), dimension(Num_nu,nchi,Num_R) :: P_syn_chi, Seed_syn_chi, Tau_syn_chi
    real(8), intent(out), dimension(nchi,Num_R) :: chi_radius, chi_gamma_bulk, chi_weight_out, B_chi_out

    real(8), allocatable, dimension(:) :: dEl,dEL_mean
    real(8), allocatable, dimension(:) :: dN_init,dN_init_log,source_q1
    real(8), allocatable, dimension(:,:) :: ulog,U_shell
    real(8), allocatable, dimension(:) :: q_grid,q_face,q_weight
    real(8), allocatable, dimension(:) :: dF1,shell_population,chi_population,dEL_mean_shell
    real(8), allocatable, dimension(:) :: V_m_chi,V_c_chi,V_a_chi,chi_weight,Epsilon_b_chi,DB_chi,t_decay_chi
    real(8), allocatable, dimension(:) :: tprop
    real(8), allocatable, dimension(:,:) :: xface,x_comov_face_hist,x_comov_hist,dx_comov_hist
    real(8), allocatable, dimension(:,:) :: rcell_hist,gcell_hist,beta
    real(8), allocatable, dimension(:) :: rcell_chi,gcell_chi,bcell_chi,brel_sh_chi
    real(8), allocatable, dimension(:,:) :: P_local,kappa2_chi
    real(8), allocatable, dimension(:) :: V_cool,P_emit_shell,Taushell,hist_inv
    real(8), allocatable, dimension(:,:) :: shell_emit,shell_tau
    real(8), allocatable, dimension(:,:,:) :: pemit,seed,Tau_hist
    real(8), allocatable, dimension(:,:,:) :: pemit_cool,seed_cool,Tau_hist_cool,Tau_pair_hist_cool &
        & ,Tau_prop_hist_cool
    real(8), allocatable, dimension(:,:) :: peff_cool_chi,seeff_cool_chi,pstream_cool,sstream_cool,hist_prefix
    real(8), allocatable, dimension(:,:) :: cooling_aux_chi,dEl_chi,dEL_mean_chi
    real(8), allocatable, dimension(:) :: adiabatic_log_coeff_chi

    real(8), dimension(ng+1) :: x_edge_E,coord_edge_E
    real(8) :: temp, dq, d_x_E, coord_scale, rloc, R_Gamma_loc, dNe, Para_N_e_ini, DB, DB_min
    real(8) :: Epsilon_b_floor, magnetic_decay_alpha_t, magnetic_decay_t0_s, stochastic_accel_norm
    real(8) :: Gam_e_max, Gam_e_max_max, Gam_e_m, Gam_e_c, Gam_e_c_diag, temp_gam, Gam_e_max_cell
    real(8) :: Gam_e_m_cell, Gam_e_c_cell, bsh, b2, b2sh, Q
    integer :: I_tobs, ichi, Num_nu_cool, k_medium, substep_limit
    logical :: magnetic_decay_active, pwn_cr_transport, free_outer_escape, emit_full_spectrum
    logical :: four_velocity_coord

    call electron_unpack_boundary(Boundary,n,Eta_0,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    if (profile_count > 0) error stop "fs_transport_2d does not support tabulated density profiles"

    allocate(dEl(ng), dEL_mean(ng-1), dEL_mean_shell(ng-1), &
             dN_init(ng), dN_init_log(ng), dF1(ng), shell_population(ng), chi_population(nchi), &
             V_m_chi(nchi), V_c_chi(nchi), V_a_chi(nchi), chi_weight(nchi), &
             Epsilon_b_chi(nchi), DB_chi(nchi), t_decay_chi(nchi), &
             source_q1(ng), ulog(ng, nchi), U_shell(ng, nchi), &
             q_grid(nchi), q_face(0:nchi), q_weight(nchi), &
             tprop(Num_R), xface(0:nchi, Num_R), &
             x_comov_face_hist(0:nchi, Num_R), x_comov_hist(nchi, Num_R), dx_comov_hist(nchi, Num_R), &
             rcell_hist(nchi, Num_R), gcell_hist(nchi, Num_R), beta(nchi, Num_R), &
             rcell_chi(nchi), gcell_chi(nchi), bcell_chi(nchi), brel_sh_chi(nchi), &
             P_local(Num_nu, nchi), kappa2_chi(ng, nchi), shell_emit(Num_nu, Num_R), &
             shell_tau(Num_nu, Num_R), &
             cooling_aux_chi(ng, nchi), dEl_chi(ng, nchi), dEL_mean_chi(ng-1, nchi), &
             adiabatic_log_coeff_chi(nchi))

    if (A_star > 0d0) then
        k_medium = 2
    else
        k_medium = 0
    end if
    Epsilon_b_floor  = Epsilon_b
    magnetic_decay_alpha_t = 0d0
    magnetic_decay_t0_s = 1d0
    if (n >= 27) then
        Epsilon_b_floor = Boundary(24)
        magnetic_decay_alpha_t = Boundary(25)
        magnetic_decay_t0_s = Boundary(26)
    end if
    block
        real(8) :: transport_model_selector, escape_mode_selector

        transport_model_selector = Boundary(n-2)
        stochastic_accel_norm = Boundary(n-1)
        escape_mode_selector = Boundary(n)
        pwn_cr_transport = nint(transport_model_selector) == 1
        free_outer_escape = nint(escape_mode_selector) == 1
    end block
    magnetic_decay_active = magnetic_decay_alpha_t < 0d0 .and. magnetic_decay_t0_s > 0d0 .and. &
                            Epsilon_b_floor > 0d0 .and. Epsilon_b_floor < Epsilon_b
    four_velocity_coord = .not. use_charint_transport .and. .not. pwn_cr_transport

    coord_scale = fourvel_scale*fourvel_scale-1d0
    emit_full_spectrum = emit_full_chi_spectrum /= 0
    Num_nu_cool = min(6, Num_nu)
    substep_limit = max(1, substep_max)
    allocate(V_cool(Num_nu_cool), P_emit_shell(Num_nu), Taushell(Num_nu), hist_inv(nchi), &
             pemit(Num_nu, nchi, Num_R), &
             seed(Num_nu, nchi, Num_R), Tau_hist(Num_nu, nchi, Num_R), &
             pemit_cool(Num_nu_cool, nchi, Num_R), seed_cool(Num_nu_cool, nchi, Num_R), &
             Tau_hist_cool(Num_nu_cool, nchi, Num_R), &
             Tau_pair_hist_cool(Num_nu_cool, nchi, Num_R), Tau_prop_hist_cool(Num_nu_cool, nchi, Num_R), &
             peff_cool_chi(Num_nu_cool, nchi), seeff_cool_chi(Num_nu_cool, nchi), &
             pstream_cool(Num_nu_cool, nchi), sstream_cool(Num_nu_cool, nchi), &
             hist_prefix(Num_nu_cool,0:nchi))
    call reduce_grid(Num_nu,V_seed,Num_nu_cool,V_cool)

    P_syn          = 0d0
    Seed_syn       = 0d0
    P_syn_chi      = 0d0
    Seed_syn_chi   = 0d0
    Tau_syn_chi    = 0d0
    shell_emit     = 0d0
    shell_tau      = 0d0
    chi_radius     = 0d0
    chi_gamma_bulk = 1d0
    chi_weight_out = 0d0
    B_chi_out      = 0d0
    V_m            = 0d0
    V_c            = 0d0
    V_a            = 0d0
    dN_gam_e       = 0d0
    dN_gam_e_total = 0d0
    ulog          = 0d0
    dF1            = 0d0
    pstream_cool  = 0d0
    sstream_cool = 0d0

    call q_geometry(nchi, dq, q_face, q_grid)
    q_weight = dq
    call integrate_proper_time(Num_R,R,R_Gamma,tprop)
    call init_front()

    do I_tobs = 2, Num_R
        rloc       = R(I_tobs-1)
        R_Gamma_loc = (R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/2d0

        call density_profile(A_star,dNe_ISM,rloc,R0,1,R_tr,f_jump,f_wide,dNe)

        call shock_scales(R_Gamma_loc, R_Tobs(I_tobs), bsh, b2, b2sh)

        if (emit_full_spectrum) then
            call syn_state(index_syn_intger,rloc,DB,ng,Num_nu,n_threads,gam_e, &
                                        dN_gam_e_total(:,I_tobs-1),V_seed,P_emit_shell, &
                                        P_syn(:,I_tobs),Seed_syn(:,I_tobs),Taushell)
            shell_emit(:,I_tobs) = P_emit_shell
            shell_tau(:,I_tobs) = Taushell
        end if

        peff_cool_chi = pemit_cool(:,:,I_tobs-1) + pstream_cool
        seeff_cool_chi = seed_cool(:,:,I_tobs-1) + sstream_cool
        call q_cell_geometry(nchi,k_medium,rloc,R_Gamma_loc,q_grid, &
                                     rcell_chi,gcell_chi,bcell_chi,brel_sh_chi)
        call update_bchi(x_comov_hist(:,I_tobs-1), R_Gamma_loc, dNe, brel_sh_chi)

        call forward_cooling(0,index_Y,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0, &
                             ng,Num_nu_cool,nchi,n_threads,gam_e,V_cool, &
                             peff_cool_chi,seeff_cool_chi,seeff_cool_chi,cooling_aux_chi,dEl_chi)
        call assemble_cooling_chi(rloc, R_Gamma_loc, bsh, seeff_cool_chi, R_Tobs(I_tobs))
        call break_freqs(R_Gamma_loc, bsh, I_tobs-1, rloc)
        call advance_shell(I_tobs)

        call store_shell(I_tobs)

        call advance_history_stream(I_tobs-1,I_tobs,Num_R,nchi,Num_nu_cool,tprop,V_cool, &
                                             x_comov_face_hist,x_comov_hist,dx_comov_hist,beta, &
                                             Tau_prop_hist_cool,pemit_cool,seed_cool,pstream_cool,sstream_cool, &
                                             hist_inv,hist_prefix)

        dN_gam_e_total(:, I_tobs) = 0d0
        do ichi = 1, nchi
            dN_gam_e_total(:, I_tobs) = dN_gam_e_total(:, I_tobs) + dN_gam_e(:,ichi,I_tobs)*dq
        end do
    end do

    call finish_output()
    deallocate(dEl, dEL_mean, dEL_mean_shell, dN_init, dN_init_log, dF1, &
               shell_population, chi_population, ulog, U_shell, source_q1, &
               V_m_chi, V_c_chi, V_a_chi, chi_weight, Epsilon_b_chi, DB_chi, t_decay_chi, &
               q_grid, q_face, q_weight, tprop, &
               xface, x_comov_face_hist, &
               x_comov_hist, dx_comov_hist, rcell_hist, gcell_hist, beta, &
               rcell_chi, gcell_chi, bcell_chi, brel_sh_chi, &
               P_local, V_cool, P_emit_shell, Taushell, hist_inv, hist_prefix, pemit, seed, Tau_hist, &
               shell_emit, shell_tau, &
               pemit_cool, seed_cool, Tau_hist_cool, Tau_pair_hist_cool, Tau_prop_hist_cool, &
               peff_cool_chi, seeff_cool_chi, pstream_cool, sstream_cool, &
               cooling_aux_chi, dEl_chi, dEL_mean_chi, adiabatic_log_coeff_chi, kappa2_chi)

contains

subroutine init_front()
    call electron_initial_density(A_star,dNe_ISM,R(1),R0,R_tr,f_jump,f_wide,dNe,Para_N_e_ini)

    DB            = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-1d0)))
    DB_min        = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-1d0)))
    Gam_e_max     = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    Gam_e_max_max = 3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam      = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-1d0)
    call electron_gm_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c       = 7.7d8/(1d0+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/2d0)

    if (four_velocity_coord) then
        call build_fourvel_grid(ng,1d0,tail_factor*Gam_e_max_max, &
                                               fourvel_scale,gam_e,coord_edge_E,x_edge_E)
        call init_coord(f_e*Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                                              ng,coord_edge_E,coord_scale,dN_init)
    else
        call electron_initialize_spectrum(ng,Gam_e_max_max,f_e*Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                          imodelog,gam_e,dN_init,x_edge_E)
        coord_edge_E = x_edge_E
    end if
    d_x_E = dlog(gam_e(2)/gam_e(1))
    dN_init_log = dN_init

    ulog(:,1) = dN_init_log / dq

    do ichi = 1, nchi
        if (four_velocity_coord) then
            call coord_to_dgamma(ng,coord_edge_E,coord_scale,gam_e, &
                                                               ulog(:,ichi),dN_gam_e(:,ichi,1))
        else
            call dnx_dgamma(ng,x_edge_E,gam_e,ulog(:,ichi),dN_gam_e(:,ichi,1))
        end if
    end do
    dN_gam_e_total(:,1) = 0d0
    do ichi = 1, nchi
        dN_gam_e_total(:,1) = dN_gam_e_total(:,1) + dN_gam_e(:,ichi,1)*dq
    end do
    call downstream_grid(nchi,k_medium,R(1),R_Gamma(1),q_face,q_grid, &
                                          xface(:,1),x_comov_face_hist(:,1),x_comov_hist(:,1),dx_comov_hist(:,1), &
                                          rcell_hist(:,1),gcell_hist(:,1),beta(:,1),brel_sh_chi)
    call update_bchi(x_comov_hist(:,1), R_Gamma(1), dNe, brel_sh_chi)
    if (nchi == 1) then
        B_chi_out(1,1) = DB
    else
        do ichi = 1, nchi
            B_chi_out(ichi,1) = DB_chi(ichi)
        end do
    end if
    if (emit_full_spectrum) then
        call syn_seed_chi(R(1),ng,Num_nu,nchi,gam_e,dN_gam_e(:,:,1),V_seed, &
                                               DB_chi,q_weight,1.046d4, &
                                               P_local,pemit(:,:,1),seed(:,:,1),Tau_hist(:,:,1))
        do ichi = 1, nchi
            call project_syn(Num_nu,V_seed,pemit(:,ichi,1),seed(:,ichi,1),Tau_hist(:,ichi,1), &
                                            Num_nu_cool,V_cool,pemit_cool(:,ichi,1), &
                                            seed_cool(:,ichi,1),Tau_hist_cool(:,ichi,1))
            Tau_pair_hist_cool(:,ichi,1) = 0d0
            call pair_tau(V_cool,Num_nu_cool,seed_cool(:,ichi,1),dx_comov_hist(ichi,1), &
                                                   Tau_pair_hist_cool(:,ichi,1))
            Tau_prop_hist_cool(:,ichi,1) = Tau_hist_cool(:,ichi,1) + Tau_pair_hist_cool(:,ichi,1)
        end do
    else
        if (magnetic_decay_active) then
            do ichi = 1, nchi
                call coolscalar(R(1),1,ichi)
            end do
        else
            call coolbatch(R(1),1)
        end if
    end if
end subroutine init_front

subroutine advance_shell(I_tobs)
    integer, intent(in) :: I_tobs
    real(8) :: dDD, dDR, dDR_try, dDR_xi, dDR_q, max_xi_coeff, chi_peak
    integer :: L, L1, src_lo, src_hi, active_hi, active_chi

    call electron_source_bounds(ng,gam_e,Gam_e_m,Gam_e_max,src_lo,src_hi)
    if (four_velocity_coord) then
        call source_coord(ng,coord_edge_E,coord_scale, &
                                                               Gam_e_m,Gam_e_max,1d0,p,dF1)
    else
        call source_edges(ng,x_edge_E,Gam_e_m,Gam_e_max,1d0,p,dF1)
    end if
    shell_population = sum(ulog, dim=2)
    chi_population = sum(ulog, dim=1)
    chi_peak = maxval(chi_population)
    active_hi = gamma_active_hi(ng,dF1,shell_population,src_lo,src_hi,chi_peak)
    if (pwn_cr_transport) then
        call q_divergence(nchi,k_medium,rloc,R_Gamma_loc,bsh,q_grid,adiabatic_log_coeff_chi)
    else
        adiabatic_log_coeff_chi = 1d0/rloc
    end if
    max_xi_coeff = max_xi_chi(ng,nchi,dEL_mean_chi, &
                                             adiabatic_log_coeff_chi,chi_population,chi_peak,active_hi)

    dDR_xi = huge(1d0)
    if (max_xi_coeff > 0d0) dDR_xi = 4d0*d_x_E/max_xi_coeff
    active_chi = chi_active_hi(nchi,chi_population,chi_peak)
    dDR_q = huge(1d0)
    if (use_charint_transport) then
        dDR_q = q_step_limit(nchi,k_medium,rloc,dq,q_face,4d0)
    end if

    dDD     = R(I_tobs)-R(I_tobs-1)
    dDR_try = min(dDD, min(dDR_xi, dDR_q))
    L1      = max(1, min(substep_limit, ceiling(dDD/dDR_try)))
    if (.not. use_charint_transport .and. .not. pwn_cr_transport) then
        call transport_step_fullhide(R(I_tobs-1), R(I_tobs), R_Gamma(I_tobs-1), R_Gamma(I_tobs), &
                                      dDD, active_hi, active_chi)
    else
        dDR     = dDD/dble(L1)

        do L = 1, L1
            block
                real(8) :: frac_sub, R_sub, gsh_sub, Gam_e_m_p

                frac_sub = (dble(L)-0.5d0)/dble(L1)
                R_sub = R(I_tobs-1) + frac_sub*dDD
                gsh_sub = (1d0-frac_sub)*R_Gamma(I_tobs-1) + frac_sub*R_Gamma(I_tobs)

                call density_profile(A_star,dNe_ISM,R_sub,R0,1,R_tr,f_jump,f_wide,dNe)

                call shock_scales(gsh_sub, R_Tobs(I_tobs), bsh, b2, b2sh)
                if (pwn_cr_transport) then
                    call q_divergence(nchi,k_medium,R_sub,gsh_sub,bsh,q_grid,adiabatic_log_coeff_chi)
                else
                    adiabatic_log_coeff_chi = 1d0/R_sub
                end if

                Gam_e_m_p = (1d0-p)/(Gam_e_max**(1d0-p)-Gam_e_m**(1d0-p))
                call electron_injection_prefactor(R_sub-0.5d0*dDR,dDR,dNe,f_e,Gam_e_m_p,Q)
                call source_edges(ng,x_edge_E,Gam_e_m,Gam_e_max,Q,p,dF1)
                source_q1 = dF1/dq

                if (use_charint_transport) then
                    call advance_q_charint(ulog, ng, nchi, active_hi, dq, q_face, &
                                                     k_medium, R_sub, 0d0*source_q1, dDR, n_threads)
                    call advance_q_diffusion(ulog, ng, nchi, active_hi, dq, q_face, &
                                                      k_medium, R_sub, gsh_sub, bsh, &
                                                      kappa2_chi, dDR, n_threads)

                    call advance_energy_charint(ulog, ng, nchi, gam_e, DB_chi, &
                                                             dEl_chi, R_sub, gsh_sub, bsh, index_Y, &
                                                             dDR, active_chi, n_threads, source_q1)
                else
                    if (pwn_cr_transport) then
                        call advance_pwncr_q(ulog, ng, nchi, active_hi, dq, q_face, &
                                                      k_medium, R_sub, gsh_sub, bsh, kappa2_chi, &
                                                      source_q1, dDR, free_outer_escape, n_threads)
                    else
                        call advance_q_implicit(ulog, ng, nchi, active_hi, dq, q_face, &
                                                k_medium, R_sub, gsh_sub, bsh, kappa2_chi, &
                                                source_q1, dDR, n_threads)
                    end if

                    if (pwn_cr_transport) then
                        if (stochastic_accel_norm > 0d0) then
                            call advance_stoch_chi(ulog, ng, nchi, &
                                                                        stochastic_accel_norm, R_sub, &
                                                                        d_x_E, 0.5d0*dDR, n_threads)
                        end if
                        call advance_pwncr_energy(ulog, ng, nchi, dEL_mean_chi, &
                                                               adiabatic_log_coeff_chi, d_x_E, dDR, n_threads)
                        if (stochastic_accel_norm > 0d0) then
                            call advance_stoch_chi(ulog, ng, nchi, &
                                                                        stochastic_accel_norm, R_sub, &
                                                                        d_x_E, 0.5d0*dDR, n_threads)
                        end if
                    else
                        call advance_energy_chi(ulog, ng, nchi, dEL_mean_chi, &
                                                         R_sub, d_x_E, dDR, n_threads)
                    end if
                end if
            end block
        end do
    end if
end subroutine advance_shell

subroutine store_shell(I_tobs)
    integer, intent(in) :: I_tobs

    call density_profile(A_star,dNe_ISM,R(I_tobs),R0,1,R_tr,f_jump,f_wide,dNe)
    call downstream_grid(nchi,k_medium,R(I_tobs),R_Gamma(I_tobs),q_face,q_grid, &
                                          xface(:,I_tobs),x_comov_face_hist(:,I_tobs),x_comov_hist(:,I_tobs), &
                                          dx_comov_hist(:,I_tobs),rcell_hist(:,I_tobs),gcell_hist(:,I_tobs), &
                                          beta(:,I_tobs),brel_sh_chi)
    call update_bchi(x_comov_hist(:,I_tobs), R_Gamma(I_tobs), dNe, brel_sh_chi)
    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. nchi*Num_nu >= 128) schedule(static) &
    !$OMP& private(ichi)
    do ichi = 1, nchi
        if (four_velocity_coord) then
            call coord_to_dgamma(ng,coord_edge_E,coord_scale,gam_e, &
                                                               ulog(:,ichi),dN_gam_e(:,ichi,I_tobs))
        else
            call dnx_dgamma(ng,x_edge_E,gam_e,ulog(:,ichi),dN_gam_e(:,ichi,I_tobs))
        end if
        B_chi_out(ichi,I_tobs) = DB_chi(ichi)
    end do
    !$OMP END PARALLEL DO
    if (.not. emit_full_spectrum) then
        call coolbatch(R(I_tobs),I_tobs)
    end if
    if (emit_full_spectrum) then
        call syn_seed_chi(R(I_tobs),ng,Num_nu,nchi,gam_e,dN_gam_e(:,:,I_tobs),V_seed, &
                                               DB_chi,q_weight,1.046d4, &
                                               P_local,pemit(:,:,I_tobs),seed(:,:,I_tobs), &
                                               Tau_hist(:,:,I_tobs))
        do ichi = 1, nchi
            call project_syn(Num_nu,V_seed,pemit(:,ichi,I_tobs), &
                                            seed(:,ichi,I_tobs),Tau_hist(:,ichi,I_tobs), &
                                            Num_nu_cool,V_cool,pemit_cool(:,ichi,I_tobs), &
                                            seed_cool(:,ichi,I_tobs),Tau_hist_cool(:,ichi,I_tobs))
        end do
        call pairall(I_tobs)
    end if
end subroutine store_shell

subroutine coolscalar(Rad,Iout,Ichi)
    real(8), intent(in) :: Rad
    integer, intent(in) :: Iout,Ichi
    real(8), dimension(1) :: DBcell, Wcell
    real(8), dimension(ng,1) :: DNcell
    real(8), dimension(Num_nu_cool,1) :: Pemcell,Psyncell,Seedcell,Taucell

    DBcell(1) = DB_chi(Ichi)
    Wcell(1) = q_weight(Ichi)
    DNcell(:,1) = dN_gam_e(:,Ichi,Iout)
    call syn_seed_chi(Rad,ng,Num_nu_cool,1,gam_e,DNcell,V_cool,DBcell, &
                                           Wcell,1.046d4, &
                                           Pemcell,Psyncell,Seedcell,Taucell)
    pemit_cool(:,Ichi,Iout) = Psyncell(:,1)
    seed_cool(:,Ichi,Iout) = Seedcell(:,1)
    Tau_hist_cool(:,Ichi,Iout) = Taucell(:,1)
    call pairone(Iout,Ichi)
end subroutine coolscalar

subroutine coolbatch(Rad,Iout)
    real(8), intent(in) :: Rad
    integer, intent(in) :: Iout

    call syn_seed_chi(Rad,ng,Num_nu_cool,nchi,gam_e,dN_gam_e(:,:,Iout),V_cool, &
                                           DB_chi,q_weight,1.046d4, &
                                           peff_cool_chi,pemit_cool(:,:,Iout), &
                                           seed_cool(:,:,Iout),Tau_hist_cool(:,:,Iout))
    call pairall(Iout)
end subroutine coolbatch

subroutine pairone(Iout,Ichi)
    integer, intent(in) :: Iout,Ichi

    Tau_pair_hist_cool(:,Ichi,Iout) = 0d0
    call pair_tau(V_cool,Num_nu_cool,seed_cool(:,Ichi,Iout), &
                                           dx_comov_hist(Ichi,Iout),Tau_pair_hist_cool(:,Ichi,Iout))
    Tau_prop_hist_cool(:,Ichi,Iout) = Tau_hist_cool(:,Ichi,Iout) + Tau_pair_hist_cool(:,Ichi,Iout)
end subroutine pairone

subroutine pairall(Iout)
    integer, intent(in) :: Iout
    integer :: J

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. nchi*Num_nu_cool >= 128) schedule(static) &
    !$OMP& private(J)
    do J = 1, nchi
        call pairone(Iout,J)
    end do
    !$OMP END PARALLEL DO
end subroutine pairall

subroutine finish_output()
    integer :: I_proj, I_proj_chi, I_nu
    real(8) :: trans

    rloc = R(Num_R)
    R_Gamma_loc = R_Gamma(Num_R)
    call density_profile(A_star,dNe_ISM,rloc,R0,1,R_tr,f_jump,f_wide,dNe)
    call shock_scales(R_Gamma_loc, R_Tobs(Num_R), bsh, b2, b2sh)

    if (emit_full_spectrum) then
        call syn_state(index_syn_intger,rloc,DB,ng,Num_nu,n_threads,gam_e, &
                                    dN_gam_e_total(:,Num_R),V_seed,P_emit_shell,P_syn(:,Num_R), &
                                    Seed_syn(:,Num_R),Taushell)
        shell_emit(:,Num_R) = P_emit_shell
        shell_tau(:,Num_R) = Taushell
    end if
    peff_cool_chi = pemit_cool(:,:,Num_R) + pstream_cool
    seeff_cool_chi = seed_cool(:,:,Num_R) + sstream_cool
    call forward_cooling(0,index_Y,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0, &
                         ng,Num_nu_cool,nchi,n_threads,gam_e,V_cool, &
                         peff_cool_chi,seeff_cool_chi,seeff_cool_chi,cooling_aux_chi,dEl_chi)
    call q_cell_geometry(nchi,k_medium,rloc,R_Gamma_loc,q_grid, &
                                 rcell_chi,gcell_chi,bcell_chi,brel_sh_chi)
    call update_bchi(x_comov_hist(:,Num_R), R_Gamma_loc, dNe, brel_sh_chi)
    call assemble_cooling_chi(rloc, R_Gamma_loc, bsh, seeff_cool_chi, R_Tobs(Num_R))
    call break_freqs(R_Gamma_loc, bsh, Num_R, rloc)
    if (emit_full_spectrum) then
        do I_proj = 1, Num_R
            if (nchi == 1) then
                chi_radius(1,I_proj) = R(I_proj)
                chi_gamma_bulk(1,I_proj) = R_Gamma(I_proj)
                chi_weight_out(1,I_proj) = dq
                P_syn_chi(:,1,I_proj) = shell_emit(:,I_proj)/dq
                Tau_syn_chi(:,1,I_proj) = shell_tau(:,I_proj)
                Seed_syn_chi(:,1,I_proj) = Seed_syn(:,I_proj)
            else
                Seed_syn_chi(:,:,I_proj) = seed(:,:,I_proj)
                Tau_syn_chi(:,:,I_proj) = Tau_hist(:,:,I_proj)
                do I_proj_chi = 1, nchi
                    chi_radius(I_proj_chi,I_proj) = rcell_hist(I_proj_chi,I_proj)
                    chi_gamma_bulk(I_proj_chi,I_proj) = gcell_hist(I_proj_chi,I_proj)
                    chi_weight_out(I_proj_chi,I_proj) = dq
                    do I_nu = 1, Num_nu
                        call transfer_factor(Tau_hist(I_nu,I_proj_chi,I_proj), trans)
                        P_syn_chi(I_nu,I_proj_chi,I_proj) = pemit(I_nu,I_proj_chi,I_proj)/trans
                    end do
                end do
            end if
        end do
    else
        do I_proj = 1, Num_R
            do I_proj_chi = 1, nchi
                if (nchi == 1) then
                    chi_radius(I_proj_chi,I_proj) = R(I_proj)
                    chi_gamma_bulk(I_proj_chi,I_proj) = R_Gamma(I_proj)
                else
                    chi_radius(I_proj_chi,I_proj) = rcell_hist(I_proj_chi,I_proj)
                    chi_gamma_bulk(I_proj_chi,I_proj) = gcell_hist(I_proj_chi,I_proj)
                end if
                chi_weight_out(I_proj_chi,I_proj) = dq
            end do
        end do
    end if
end subroutine finish_output

subroutine transport_step_fullhide(R_prev, R_curr, Gamma_prev, Gamma_curr, dDR_step, &
                                    active_hi, active_chi)
    real(8), intent(in) :: R_prev, R_curr, Gamma_prev, Gamma_curr, dDR_step
    integer, intent(in) :: active_hi, active_chi
    real(8) :: R_sub, R_eff, half_dR, gsh_sub, Gam_e_m_p, adiabatic_integral, face_coord
    real(8) :: DB_loc, Gam_e_max_loc, Gam_e_m_loc, temp_gam_loc
    real(8), dimension(ng) :: source_q1_loc,dF1_zero
    real(8) :: Q_rate_loc
    real(8), dimension(ng-1) :: coord_face_step_loc, face_invjac
    integer :: I

    R_sub = 0.5d0*(R_prev+R_curr)
    gsh_sub = 0.5d0*(Gamma_prev+Gamma_curr)
    half_dR = 0.5d0*dDR_step
    if (R_curr > R_prev) then
        R_eff = dDR_step/dlog(R_curr/R_prev)
    else
        R_eff = R_sub
    end if

    call shock_state(gsh_sub, bsh, b2, b2sh)
    if (pwn_cr_transport) then
        call q_divergence(nchi,k_medium,R_sub,gsh_sub,bsh,q_grid,adiabatic_log_coeff_chi)
    else
        adiabatic_log_coeff_chi = 1d0/R_sub
    end if

    call density_profile(A_star,dNe_ISM,R_sub,R0,1,R_tr,f_jump,f_wide,dNe)
    DB_loc = 0.39d0*dsqrt(Epsilon_b*dNe*(gsh_sub*(gsh_sub-1d0)))
    Gam_e_max_loc = 3d0*Para_m_energy/dsqrt(8d0*DB_loc*Para_e**3)
    temp_gam_loc = Epsilon_e/f_e*para_m_p/para_m_e*(gsh_sub-1d0)
    call electron_gm_exact(p,temp_gam_loc,Gam_e_max_loc,Gam_e_m_loc)
    Gam_e_m_p = (1d0-p)/(Gam_e_max_loc**(1d0-p)-Gam_e_m_loc**(1d0-p))
    call electron_injection_prefactor(R_prev,dDR_step,dNe,f_e,Gam_e_m_p,Q_rate_loc)
    call source_coord(ng,coord_edge_E,coord_scale, &
                                                           Gam_e_m_loc,Gam_e_max_loc, &
                                                           Q_rate_loc*dDR_step/dq,p,source_q1_loc)

    U_shell = ulog
    dF1_zero = 0d0
    call advance_q_charint(U_shell, ng, nchi, active_hi, dq, q_face, &
                                     k_medium, R_eff, dF1_zero, half_dR, n_threads)
    call advance_q_diffusion(U_shell, ng, nchi, active_hi, dq, q_face, &
                                      k_medium, R_eff, gsh_sub, bsh, &
                                      kappa2_chi, half_dR, n_threads)

    adiabatic_integral = dlog(R_curr/R_prev)
    do I = 1, ng-1
        face_coord = coord_edge_E(I+1)
        face_invjac(I) = 1d0/dxg_dcoord(coord_fourvel,coord_scale,face_coord)
    end do
    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. active_chi*ng >= 512) schedule(static) &
    !$OMP& private(I,ichi,coord_face_step_loc)
    do ichi = 1, active_chi
        do I = 1, ng-1
            coord_face_step_loc(I) = &
                (dDR_step*(dEl_chi(I,ichi)+dEl_chi(I+1,ichi))/2d0+adiabatic_integral)*face_invjac(I)
        end do
        if (ichi == 1) then
            call flux_seq_nonuniform(ng,coord_edge_E,coord_face_step_loc, &
                                                                  source_q1_loc,U_shell(:,ichi),ulog(:,ichi),.true.)
        else
            call flux_seq_nonuniform(ng,coord_edge_E,coord_face_step_loc, &
                                                                  dF1_zero,U_shell(:,ichi),ulog(:,ichi),.true.)
        end if
    end do
    !$OMP END PARALLEL DO
    if (active_chi < nchi) ulog(:,active_chi+1:nchi) = U_shell(:,active_chi+1:nchi)

    U_shell = ulog
    call advance_q_charint(U_shell, ng, nchi, active_hi, dq, q_face, &
                                     k_medium, R_eff, dF1_zero, half_dR, n_threads)
    call advance_q_diffusion(U_shell, ng, nchi, active_hi, dq, q_face, &
                                      k_medium, R_eff, gsh_sub, bsh, &
                                      kappa2_chi, half_dR, n_threads)
    ulog = U_shell
end subroutine transport_step_fullhide

subroutine update_bchi(x_comov_col, gf, dNe_val, brel_sh)
    real(8), intent(in), dimension(nchi) :: x_comov_col,brel_sh
    real(8), intent(in) :: gf,dNe_val
    real(8) :: Gam_m1
    integer :: I
    Gam_m1 = gf*(gf-1d0)
    do I = 1, nchi
        if (brel_sh(I) <= 0d0) error stop 'electron 2d magnetic decay requires positive shock-relative speed'
        t_decay_chi(I) = x_comov_col(I)/(brel_sh(I)*para_c)
        if (magnetic_decay_active) then
            Epsilon_b_chi(I) = Epsilon_b_floor + (Epsilon_b-Epsilon_b_floor) * &
                               (1d0+t_decay_chi(I)/magnetic_decay_t0_s)**magnetic_decay_alpha_t
        else
            Epsilon_b_chi(I) = Epsilon_b
        end if
        DB_chi(I) = 0.39d0*dsqrt(Epsilon_b_chi(I)*dNe_val*Gam_m1)
    end do
end subroutine update_bchi

subroutine shock_scales(gf, R_Tobs, bsh, b2, b2sh)
    real(8), intent(in) :: gf, R_Tobs
    real(8), intent(out) :: bsh, b2, b2sh
    DB = 0.39d0*dsqrt(Epsilon_b*dNe*(gf*(gf-1d0)))
    Gam_e_max = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(gf-1d0)
    call electron_gm_exact(p, temp_gam, Gam_e_max, Gam_e_m)
    Gam_e_c = 7.7d8*(1d0+z)/gf/DB**2/R_Tobs
    call shock_state(gf, bsh, b2, b2sh)
end subroutine shock_scales

subroutine assemble_cooling_chi(Rad, gf, bsh, seeff_chi, R_Tobs_val)
    real(8), intent(in), dimension(Num_nu_cool,nchi) :: seeff_chi
    real(8), intent(in) :: Rad,gf,bsh,R_Tobs_val
    real(8), dimension(Num_nu_cool,1) :: seeff_column
    real(8), dimension(ng,1) :: cooling_aux_column,dEl_column
    if (magnetic_decay_active) then
        !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. nchi*ng >= 512) schedule(static) &
        !$OMP& private(ichi,Gam_e_max_cell,temp_gam,Gam_e_m_cell,Gam_e_c_cell, &
        !$OMP&         seeff_column,cooling_aux_column,dEl_column)
        do ichi = 1, nchi
            Gam_e_max_cell = 3d0*Para_m_energy/dsqrt(8d0*DB_chi(ichi)*Para_e**3)
            temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(gf-1d0)
            call electron_gm_exact(p, temp_gam, Gam_e_max_cell, Gam_e_m_cell)
            Gam_e_c_cell = 7.7d8*(1d0+z)/gf/DB_chi(ichi)**2/R_Tobs_val
            seeff_column(:,1)=seeff_chi(:,ichi)
            cooling_aux_column(:,1)=cooling_aux_chi(:,ichi)
            call forward_cooling(2,index_Y, Epsilon_e, Epsilon_b_chi(ichi), p, DB_chi(ichi), &
                                 Gam_e_m_cell, Gam_e_c_cell, Gam_e_max_cell, Rad, gf, &
                                 bsh, dNe, ng, Num_nu_cool, 1, 1, gam_e, &
                                 V_cool, seeff_column, seeff_column, seeff_column, &
                                 cooling_aux_column, dEl_column)
            dEl_chi(:,ichi)=dEl_column(:,1)
            dEL_mean_chi(:,ichi)=(dEl_chi(2:ng,ichi)+dEl_chi(1:ng-1,ichi))/2d0
            kappa2_chi(:,ichi) = gam_e*Para_m_energy*para_c/(3d0*Para_e*DB_chi(ichi))
        end do
        !$OMP END PARALLEL DO
    else
        call forward_cooling(2,index_Y, Epsilon_e, Epsilon_b, p, DB, Gam_e_m, Gam_e_c, &
                             Gam_e_max, Rad, gf, &
                             bsh, dNe, ng, Num_nu_cool, nchi, n_threads, gam_e, V_cool, &
                             seeff_chi, seeff_chi, seeff_chi, cooling_aux_chi, dEl_chi)
        do ichi = 1, nchi
            dEL_mean_chi(:,ichi)=(dEl_chi(2:ng,ichi)+dEl_chi(1:ng-1,ichi))/2d0
        end do
        do ichi = 1, nchi
            kappa2_chi(:,ichi) = gam_e*Para_m_energy*para_c/(3d0*Para_e*DB)
        end do
    end if
end subroutine assemble_cooling_chi

subroutine break_freqs(gf, bsh, I_tobs_out, Rad_val)
    real(8), intent(in) :: gf, bsh
    integer, intent(in) :: I_tobs_out
    real(8), intent(in) :: Rad_val
    if (emit_full_spectrum) then
        call nua_path(Num_nu, nchi, V_seed, Tau_hist(:,:,I_tobs_out), V_a_chi)
    else
        call nua_path(Num_nu_cool, nchi, V_cool, Tau_hist_cool(:,:,I_tobs_out), V_a_chi)
    end if
    temp = 0d0
    Q = 0d0
    do ichi = 1, nchi
        chi_weight(ichi) = sum(ulog(:,ichi))
        if (chi_weight(ichi) > 0d0) then
            if (emit_full_spectrum) then
                V_m_chi(ichi) = logparabola_peak(Num_nu, V_seed, pemit(:,ichi,I_tobs_out))
            else
                V_m_chi(ichi) = logparabola_peak(Num_nu_cool, V_cool, pemit_cool(:,ichi,I_tobs_out))
            end if
            dEL_mean_shell = dEL_mean_chi(:,ichi)
            call electron_gc_loss(ng, gam_e, dEL_mean_shell, Rad_val, Gam_e_c_diag)
            V_c_chi(ichi) = 4.2d6*DB_chi(ichi)*Gam_e_c_diag*Gam_e_c_diag
            temp = temp + chi_weight(ichi)*dlog(V_m_chi(ichi))
            Q = Q + chi_weight(ichi)
        end if
    end do
    V_m(I_tobs_out) = dexp(temp/Q)/(gf*(1d0-bsh)*(1d0+z))
    temp = 0d0
    do ichi = 1, nchi
        if (chi_weight(ichi) > 0d0) temp = temp + chi_weight(ichi)*dlog(V_c_chi(ichi))
    end do
    V_c(I_tobs_out) = dexp(temp/Q)/(gf*(1d0-bsh)*(1d0+z))
    V_a(I_tobs_out) = V_a_chi(nchi)/(gf*(1d0-bsh)*(1d0+z))
end subroutine break_freqs

end subroutine fs_transport_2d
