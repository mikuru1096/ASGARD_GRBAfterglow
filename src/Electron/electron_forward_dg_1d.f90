! Forward-shock electron solver using a moving multi-domain LGL DG discretization.
subroutine fs_electron_dg_1d(Boundary, R_Tobs, R_Gamma, R, V_seed, n, Num_nu, Num_R, Num_gam_e, &
                             index_Y, index_syn_intger, n_threads, gam_e, dN_gam_e, P_syn, Seed_syn, V_m, V_c, V_a)
    use constants
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use dynamics_common, only: dynamics_external_density_profile, active_density_jump_count, &
                               active_density_jump_factor, active_density_jump_r, active_density_jump_width, &
                               active_density_profile_count
    use electron_common, only: electron_initial_density, electron_unpack_boundary, electron_gamma_m_exact, &
                               electron_exp_tail_grid_factor, electron_injection_prefactor
    use electron_cooling_kernel, only: get_forward_cooling
    use electron_energy_coordinate_common, only: electron_build_four_velocity_grid, electron_four_velocity_grid_gamma_scale
    use electron_transport_dg_1d_kernel, only: electron_dg1d_mesh, electron_dg1d_build_four_velocity_mesh, &
                                               electron_dg1d_initial_state, electron_dg1d_project_state, &
                                               electron_dg1d_gamma_nodes, electron_dg1d_project_source, &
                                               electron_dg1d_advance_step, electron_dg1d_project_to_coord_cells, &
                                               electron_dg1d_limit_positive_cell_preserving, electron_dg1d_integral, &
                                               electron_dg1d_tail_moment_fraction, &
                                               electron_dg1d_apply_positive_kernel_filter
    use electron_shell_transport_common, only: electron_shell_dcoord_to_dndgamma_exp_centers
    use electron_injection_profiles, only: electron_initial_powerlaw_exp_cutoff_coord_edges, &
                                           electron_build_source_term_exp_cutoff_coord_edges
    use electron_radiation_kernel, only: get_syn_selected_state, get_nu_a_from_tau_grid
    implicit none

    integer, intent(in) :: n, Num_nu, Num_R, Num_gam_e, index_Y, index_syn_intger, n_threads
    real(8), intent(in) :: Boundary(n), R_Tobs(Num_R), R_Gamma(Num_R), R(Num_R), V_seed(Num_nu)
    real(8), intent(out) :: gam_e(Num_gam_e), dN_gam_e(Num_gam_e,Num_R), P_syn(Num_nu,Num_R)
    real(8), intent(out) :: Seed_syn(Num_nu,Num_R), V_m(Num_R), V_c(Num_R), V_a(Num_R)

    type(electron_dg1d_mesh) :: mesh, new_mesh
    real(8), allocatable :: state(:), projected(:), gamma_dg(:), dEl_dg(:), dEl_dg_base(:), source_dg(:), source_template(:)
    real(8), allocatable :: dN_x_init(:), x_edge(:), coord_edge(:), source_grid(:), P_emit_shell(:), Tau_syn_shell(:)
    real(8) :: Eta_0, R_ini, Epsilon_e, Epsilon_b, p, z, dNe_ISM, A_star, E_iso, T_log10_duration
    real(8) :: f_e, R_tr, f_jump, f_wide, R0, dNe, Para_N_e_ini, DB, DB_min, Gam_e_max_max
    real(8) :: Gam_e_max, Gam_e_m, Gam_e_c, temp_gam, beta_Gam, dDD, R_loc, R_Gamma_loc
    real(8) :: dNe_shell, dNe_step, DB_step, Gam_e_max_step, Gam_e_m_step, Gam_e_m_p_step, source_norm, temp
    real(8) :: Gam_e_max_inj_shell, Gam_e_m_inj_shell, Gam_e_m_p_shell
    real(8) :: dR_base, dR_step, R_end, R_mid, dg_gamma_scale, coord_scale
    real(8) :: grid_content, dg_content, source_cache_gamma_m, source_cache_gamma_max
    real(8), parameter :: dg_base_substeps = 10d0, dg_jump_window_sigma = 4d0
    real(8), parameter :: dg_jump_substeps_per_sigma = 8d0, dg_jump_log_density_step = 5d-2
    real(8), parameter :: dg_tail_moment_threshold = 1d-10, dg_tail_moment_power = 2d0
    integer :: I_tobs, source_cache_ntot
    logical :: source_cache_ready, uniform_density_shell

    allocate(dN_x_init(Num_gam_e), x_edge(Num_gam_e+1), coord_edge(Num_gam_e+1), source_grid(Num_gam_e), &
             P_emit_shell(Num_nu), Tau_syn_shell(Num_nu))

    call electron_unpack_boundary(Boundary, n, Eta_0, R_ini, Epsilon_e, Epsilon_b, p, z, dNe_ISM, A_star, &
                                  E_iso, T_log10_duration, f_e, R_tr, f_jump, f_wide, R0)
    P_syn = zero
    Seed_syn = zero
    V_m = zero
    V_c = zero
    V_a = zero
    source_cache_ready = .false.
    source_cache_ntot = 0
    source_cache_gamma_m = -one
    source_cache_gamma_max = -one

    call electron_initial_density(A_star, dNe_ISM, R_ini, R(1), R0, dNe, Para_N_e_ini)
    DB = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1) - one)))
    Gam_e_max = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R) - one)))
    Gam_e_max_max = 3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1) - one)
    call electron_gamma_m_exact(p, temp_gam, Gam_e_max, Gam_e_m)
    Gam_e_c = 7.7d8/(one + dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    call initialize_forward_four_velocity_grid()

    do I_tobs = 2, Num_R
        call prepare_shell(I_tobs)
        call write_radiation_and_breaks(I_tobs)
        call remesh_shell(Gam_e_max, Gam_e_m, Gam_e_c, Gam_e_max)

        R_end = R(I_tobs)
        dR_base = dDD/dg_base_substeps
        do while (R_loc < R_end)
            dR_step = min(dR_base, R_end - R_loc)
            call limit_density_jump_step(R_loc, R_end, dR_step)
            R_mid = R_loc + 0.5d0*dR_step
            call advance_substep(dR_step, R_mid)
            R_loc = R_loc + dR_step
        enddo
        call write_positive_output(I_tobs)
    enddo

    deallocate(state, dN_x_init, x_edge, coord_edge, source_grid, P_emit_shell, Tau_syn_shell)
    if (allocated(projected)) deallocate(projected)
    if (allocated(gamma_dg)) deallocate(gamma_dg, dEl_dg, dEl_dg_base, source_dg, source_template)

    contains

    subroutine initialize_forward_four_velocity_grid()
        dg_gamma_scale = electron_four_velocity_grid_gamma_scale
        coord_scale = dg_gamma_scale*dg_gamma_scale - one
        call electron_build_four_velocity_grid(Num_gam_e, one, electron_exp_tail_grid_factor*Gam_e_max_max, &
                                               dg_gamma_scale, gam_e, coord_edge, x_edge)
        call electron_dg1d_build_four_velocity_mesh(x_edge(1), forward_dg_active_xmax(Gam_e_max), dlog10(Gam_e_m), &
                                                    dlog10(Gam_e_c), dlog10(Gam_e_max), dg_gamma_scale, mesh)
        allocate(state(mesh%ntot))
        call electron_dg1d_initial_state(mesh, Para_N_e_ini, p, Gam_e_m, Gam_e_c, Gam_e_max, state)
        call electron_initial_powerlaw_exp_cutoff_coord_edges(Para_N_e_ini, p, Gam_e_m, Gam_e_c, Gam_e_max, &
                                                              Num_gam_e, coord_edge, coord_scale, dN_x_init)
        call scale_dg_state_to_grid_content(state, dN_x_init)
        call electron_dg1d_project_to_coord_cells(mesh, state, Num_gam_e, coord_edge, source_grid)
        call electron_shell_dcoord_to_dndgamma_exp_centers(Num_gam_e, coord_edge, coord_scale, gam_e, &
                                                           source_grid, dN_gam_e(:,1))
        call enforce_output_positivity(1)
    end subroutine initialize_forward_four_velocity_grid

    subroutine prepare_shell(I_tobs_local)
        integer, intent(in) :: I_tobs_local

        R_loc = R(I_tobs_local - 1)
        R_Gamma_loc = (R_Gamma(I_tobs_local) + R_Gamma(I_tobs_local - 1))/two
        if (R_Gamma_loc < one) error stop 'fs_electron_dg_1d requires Gamma >= 1'
        beta_Gam = dsqrt(one - one/R_Gamma_loc**2)
        call dynamics_external_density_profile(A_star, dNe_ISM, R_loc, R0, 1, R_tr, f_jump, f_wide, dNe)
        dNe_shell = dNe
        DB = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc - one)))
        Gam_e_max = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc - one)
        call electron_gamma_m_exact(p, temp_gam, Gam_e_max, Gam_e_m)
        Gam_e_max_inj_shell = Gam_e_max
        Gam_e_m_inj_shell = Gam_e_m
        Gam_e_m_p_shell = (one - p)/(Gam_e_max_inj_shell**(one - p) - Gam_e_m_inj_shell**(one - p))
        uniform_density_shell = (A_star <= zero .and. active_density_profile_count == 0 .and. &
                                 active_density_jump_count == 0 .and. abs(f_jump - one) <= zero)
        Gam_e_c = 7.7d8*(one + z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs_local)
        dDD = R(I_tobs_local) - R(I_tobs_local - 1)
    end subroutine prepare_shell

    subroutine write_radiation_and_breaks(I_tobs_local)
        integer, intent(in) :: I_tobs_local

        V_m(I_tobs_local - 1) = 4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(one - beta_Gam)*(one + z))
        V_c(I_tobs_local - 1) = 4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(one - beta_Gam)*(one + z))
        call get_syn_selected_state(index_syn_intger, R_loc, DB, Num_gam_e, Num_nu, n_threads, &
                                    gam_e, dN_gam_e(:,I_tobs_local - 1), V_seed, P_emit_shell, &
                                    P_syn(:,I_tobs_local), Seed_syn(:,I_tobs_local), Tau_syn_shell)
        call get_nu_a_from_tau_grid(Num_nu, V_seed, Tau_syn_shell, temp)
        V_a(I_tobs_local - 1) = temp/(R_Gamma_loc*(one - beta_Gam)*(one + z))
    end subroutine write_radiation_and_breaks

    subroutine remesh_shell(gamma_upper, gamma_m_break, gamma_c_break, gamma_max_break)
        real(8), intent(in) :: gamma_upper, gamma_m_break, gamma_c_break, gamma_max_break

        call electron_dg1d_build_four_velocity_mesh(x_edge(1), forward_dg_active_xmax(gamma_upper), &
                                                    dlog10(gamma_m_break), dlog10(gamma_c_break), &
                                                    dlog10(gamma_max_break), dg_gamma_scale, new_mesh)
        call ensure_dg_work(new_mesh%ntot)
        call electron_dg1d_project_state(mesh, state, new_mesh, projected)
        call electron_dg1d_apply_positive_kernel_filter(new_mesh, projected)
        call electron_dg1d_limit_positive_cell_preserving(new_mesh, projected)
        if (size(state) /= new_mesh%ntot) then
            deallocate(state)
            allocate(state(new_mesh%ntot))
        endif
        state = projected
        mesh = new_mesh
        source_cache_ready = .false.
        call prepare_dg_cooling()
    end subroutine remesh_shell

    subroutine ensure_dg_work(ntot)
        integer, intent(in) :: ntot

        if (allocated(projected)) then
            if (size(projected) /= ntot) deallocate(projected, gamma_dg, dEl_dg, dEl_dg_base, source_dg, source_template)
        endif
        if (.not. allocated(projected)) allocate(projected(ntot), gamma_dg(ntot), dEl_dg(ntot), dEl_dg_base(ntot), &
                                                 source_dg(ntot), source_template(ntot))
    end subroutine ensure_dg_work

    subroutine prepare_dg_cooling()
        call electron_dg1d_gamma_nodes(mesh, gamma_dg)
        call get_forward_cooling(index_Y, Epsilon_e, Epsilon_b, p, DB, Gam_e_m, Gam_e_c, Gam_e_max, &
                                 R_loc, R_Gamma_loc, beta_Gam, dNe_shell, mesh%ntot, Num_nu, n_threads, &
                                 gamma_dg, V_seed, P_syn(:,I_tobs), Seed_syn(:,I_tobs), dEl_dg_base)
    end subroutine prepare_dg_cooling

    subroutine advance_substep(dR_local, R_step)
        real(8), intent(in) :: dR_local, R_step

        if (uniform_density_shell) then
            dNe_step = dNe_shell
            DB_step = DB
            Gam_e_max_step = Gam_e_max_inj_shell
            Gam_e_m_step = Gam_e_m_inj_shell
            Gam_e_m_p_step = Gam_e_m_p_shell
        else
            call dynamics_external_density_profile(A_star, dNe_ISM, R_step, R0, 1, R_tr, f_jump, f_wide, dNe_step)
            DB_step = 0.39d0*dsqrt(Epsilon_b*dNe_step*(R_Gamma_loc*(R_Gamma_loc - one)))
            Gam_e_max_step = 3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
            temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc - one)
            call electron_gamma_m_exact(p, temp_gam, Gam_e_max_step, Gam_e_m_step)
            Gam_e_m_p_step = (one - p)/(Gam_e_max_step**(one - p) - Gam_e_m_step**(one - p))
        endif
        call electron_injection_prefactor(R_step, dR_local, dNe_step, f_e, Gam_e_m_p_step, source_norm)
        if (forward_dg_source_upper_xmax(Gam_e_max_step) > mesh%x_gamma(mesh%ntot)) &
            call remesh_shell(Gam_e_max_step, Gam_e_m_step, Gam_e_c, Gam_e_max_step)

        call prepare_forward_dg_source()
        if (dNe_shell > zero) then
            dEl_dg = dEl_dg_base*(dNe_step/dNe_shell)
        else
            dEl_dg = dEl_dg_base
        endif
        call electron_dg1d_advance_step(mesh, one/R_step, dR_local, dEl_dg, source_dg, state, projected)
        call electron_dg1d_apply_positive_kernel_filter(mesh, projected)
        call electron_dg1d_limit_positive_cell_preserving(mesh, projected)
        state = projected
    end subroutine advance_substep

    subroutine prepare_forward_dg_source()
        if (source_norm <= zero) then
            source_dg = zero
            return
        endif
        if ((.not. source_cache_ready) .or. source_cache_ntot /= mesh%ntot .or. &
            source_cache_gamma_m /= Gam_e_m_step .or. source_cache_gamma_max /= Gam_e_max_step) then
            call electron_dg1d_project_source(mesh, one, p, Gam_e_m_step, Gam_e_max_step, source_template)
            call electron_build_source_term_exp_cutoff_coord_edges(Num_gam_e, coord_edge, coord_scale, &
                                                                   Gam_e_m_step, Gam_e_max_step, one, p, source_grid)
            call scale_dg_state_to_grid_content(source_template, source_grid)
            source_cache_ready = .true.
            source_cache_ntot = mesh%ntot
            source_cache_gamma_m = Gam_e_m_step
            source_cache_gamma_max = Gam_e_max_step
        endif
        source_dg = source_norm*source_template
    end subroutine prepare_forward_dg_source

    subroutine limit_density_jump_step(R_left, R_stop, dR_limited)
        real(8), intent(in) :: R_left, R_stop
        real(8), intent(inout) :: dR_limited
        integer :: j

        if (active_density_jump_count > 0) then
            do j = 1, active_density_jump_count
                call limit_one_density_jump(active_density_jump_r(j), active_density_jump_factor(j), &
                                            active_density_jump_width(j), R_left, R_stop, dR_limited)
            enddo
        else if (abs(f_jump - one) > zero) then
            call limit_one_density_jump(R_tr, f_jump, f_wide, R_left, R_stop, dR_limited)
        endif
        call limit_density_log_slope(R_left, dR_limited)
    end subroutine limit_density_jump_step

    subroutine limit_one_density_jump(R_jump, jump_factor, width_frac, R_left, R_stop, dR_limited)
        real(8), intent(in) :: R_jump, jump_factor, width_frac, R_left, R_stop
        real(8), intent(inout) :: dR_limited
        real(8) :: sigma_r, marks(7), R_trial, eps_r
        integer :: i

        if (abs(jump_factor - one) <= zero) return
        sigma_r = width_frac*R_jump
        marks = R_jump + sigma_r*(/-4d0, -2d0, -one, zero, one, 2d0, 4d0/)
        R_trial = min(R_left + dR_limited, R_stop)
        eps_r = 1d-12*max(R_stop, one)
        do i = 1, 7
            if (marks(i) > R_left + eps_r .and. marks(i) < R_trial + eps_r) then
                dR_limited = min(dR_limited, marks(i) - R_left)
                R_trial = R_left + dR_limited
            endif
        enddo
        if (R_left < R_jump + dg_jump_window_sigma*sigma_r .and. &
            R_trial > R_jump - dg_jump_window_sigma*sigma_r) then
            dR_limited = min(dR_limited, sigma_r/dg_jump_substeps_per_sigma)
        endif
    end subroutine limit_one_density_jump

    subroutine limit_density_log_slope(R_left, dR_limited)
        real(8), intent(in) :: R_left
        real(8), intent(inout) :: dR_limited
        real(8) :: R_probe, slope

        R_probe = R_left + 0.5d0*dR_limited
        slope = density_jump_log_slope(R_probe)
        if (abs(slope) > zero) dR_limited = min(dR_limited, dg_jump_log_density_step/abs(slope))
    end subroutine limit_density_log_slope

    real(8) function density_jump_log_slope(radius) result(slope)
        real(8), intent(in) :: radius
        real(8) :: enhancement, derivative
        integer :: j

        enhancement = one
        derivative = zero
        if (active_density_jump_count > 0) then
            do j = 1, active_density_jump_count
                call add_density_jump_derivative(radius, active_density_jump_r(j), active_density_jump_factor(j), &
                                                 active_density_jump_width(j), enhancement, derivative)
            enddo
        else if (abs(f_jump - one) > zero) then
            call add_density_jump_derivative(radius, R_tr, f_jump, f_wide, enhancement, derivative)
        endif
        slope = derivative/enhancement
    end function density_jump_log_slope

    subroutine add_density_jump_derivative(radius, R_jump, jump_factor, width_frac, enhancement, derivative)
        real(8), intent(in) :: radius, R_jump, jump_factor, width_frac
        real(8), intent(inout) :: enhancement, derivative
        real(8) :: sigma_r, profile

        if (abs(jump_factor - one) <= zero) return
        sigma_r = width_frac*R_jump
        profile = (jump_factor - one)*dexp(-((radius - R_jump)*(radius - R_jump))/(two*sigma_r*sigma_r))
        enhancement = enhancement + profile
        derivative = derivative - profile*(radius - R_jump)/(sigma_r*sigma_r)
    end subroutine add_density_jump_derivative

    real(8) function forward_dg_source_upper_xmax(gamma_upper) result(x_max)
        real(8), intent(in) :: gamma_upper
        real(8) :: gamma_grid_max

        gamma_grid_max = ten**x_edge(Num_gam_e+1)
        x_max = dlog10(min(gamma_grid_max, electron_exp_tail_grid_factor*gamma_upper))
    end function forward_dg_source_upper_xmax

    real(8) function forward_dg_active_xmax(gamma_upper) result(x_max)
        real(8), intent(in) :: gamma_upper
        real(8) :: tail_fraction

        x_max = forward_dg_source_upper_xmax(gamma_upper)
        if (allocated(state)) then
            if (x_max < mesh%x_gamma(mesh%ntot)) then
                call electron_dg1d_tail_moment_fraction(mesh, state, ten**x_max, &
                                                        dg_tail_moment_power, tail_fraction)
                if (tail_fraction > dg_tail_moment_threshold) x_max = mesh%x_gamma(mesh%ntot)
            endif
        endif
    end function forward_dg_active_xmax

    subroutine scale_dg_state_to_grid_content(dg_state, grid_state)
        real(8), intent(inout) :: dg_state(:)
        real(8), intent(in) :: grid_state(Num_gam_e)

        grid_content = sum(grid_state*(coord_edge(2:Num_gam_e+1) - coord_edge(1:Num_gam_e)))
        call electron_dg1d_integral(mesh, dg_state, dg_content)
        if (grid_content > zero) then
            if (dg_content <= zero) error stop 'fs_electron_dg_1d source projection has non-positive DG content'
            dg_state = dg_state*(grid_content/dg_content)
        endif
    end subroutine scale_dg_state_to_grid_content

    subroutine write_positive_output(I_tobs_local)
        integer, intent(in) :: I_tobs_local
        real(8) :: projected_content, dg_content_local

        call electron_dg1d_project_to_coord_cells(mesh, state, Num_gam_e, coord_edge, source_grid)
        call electron_shell_dcoord_to_dndgamma_exp_centers(Num_gam_e, coord_edge, coord_scale, gam_e, &
                                                           source_grid, dN_gam_e(:,I_tobs_local))
        call enforce_output_positivity(I_tobs_local)
        call electron_dg1d_integral(mesh, state, dg_content_local)
        projected_content = sum(dN_gam_e(:,I_tobs_local)*gam_e*dlog(ten)*(x_edge(2:Num_gam_e+1) - x_edge(1:Num_gam_e)))
        if (.not. (dg_content_local > zero .and. ieee_is_finite(dg_content_local))) &
            error stop 'fs_electron_dg_1d output projection has non-positive DG content'
        if (.not. (projected_content > zero .and. ieee_is_finite(projected_content))) &
            error stop 'fs_electron_dg_1d output projection lost positive content'
    end subroutine write_positive_output

    subroutine enforce_output_positivity(I_tobs_local)
        integer, intent(in) :: I_tobs_local

        where (dN_gam_e(:,I_tobs_local) > zero .and. ieee_is_finite(dN_gam_e(:,I_tobs_local)))
            dN_gam_e(:,I_tobs_local) = dN_gam_e(:,I_tobs_local)
        elsewhere
            dN_gam_e(:,I_tobs_local) = zero
        end where
    end subroutine enforce_output_positivity

end subroutine fs_electron_dg_1d
