! Forward-shock electron solver using a moving multi-domain LGL DG discretization.
subroutine fs_electron_dg_1d(Boundary, R_Tobs, R_Gamma, R, V_seed, n, Num_nu, Num_R, Num_gam_e, &
                             index_Y, index_syn_intger, n_threads, gam_e, dN_gam_e, P_syn, Seed_syn, V_m, V_c, V_a)
    use constants
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use dynamics_common, only: dynamics_external_density_profile, active_density_jump_count, &
                               active_density_jump_factor, active_density_jump_r, active_density_jump_width
    use electron_common, only: electron_initial_density, electron_unpack_boundary, electron_gamma_m_exact, &
                               electron_initialize_spectrum, electron_initial_grid_log_edges, &
                               electron_injection_prefactor, electron_prepare_radiation_spectrum
    use electron_cooling_kernel, only: get_forward_cooling
    use electron_transport_dg_1d_kernel, only: electron_dg1d_mesh, electron_dg1d_build_mesh, &
                                              electron_dg1d_initial_state, electron_dg1d_project_state, &
                                              electron_dg1d_gamma_nodes, electron_dg1d_project_source, &
                                              electron_dg1d_advance_step, electron_dg1d_project_to_log_cells, &
                                              electron_dg1d_limit_positive_cell_preserving, electron_dg1d_integral
    use electron_injection_profiles, only: electron_profile_log_cell_edges, electron_build_source_term_exp_cutoff_edges
    use electron_radiation_kernel, only: get_nu_a, get_syn_selected
    use electron_transport_common, only: electron_dnx_to_dndgamma_exp_centers
    implicit none

    integer, intent(in) :: n, Num_nu, Num_R, Num_gam_e, index_Y, index_syn_intger, n_threads
    real(8), intent(in) :: Boundary(n), R_Tobs(Num_R), R_Gamma(Num_R), R(Num_R), V_seed(Num_nu)
    real(8), intent(out) :: gam_e(Num_gam_e), dN_gam_e(Num_gam_e,Num_R), P_syn(Num_nu,Num_R)
    real(8), intent(out) :: Seed_syn(Num_nu,Num_R), V_m(Num_R), V_c(Num_R), V_a(Num_R)

    type(electron_dg1d_mesh) :: mesh, new_mesh
    real(8), allocatable :: state(:), projected(:), gamma_dg(:), dEl_dg(:), dEl_dg_base(:), source_dg(:)
    real(8), allocatable :: dN_x_init(:), x_edge(:), gam_e_rad(:), dN_gam_e_rad(:), source_grid(:)
    real(8) :: Eta_0, R_ini, Epsilon_e, Epsilon_b, p, z, dNe_ISM, A_star, E_iso, T_log10_duration
    real(8) :: f_e, R_tr, f_jump, f_wide, R0, dNe, Para_N_e_ini, DB, DB_min, Gam_e_max_max
    real(8) :: Gam_e_max, Gam_e_m, Gam_e_c, temp_gam, beta_Gam, dDD, R_loc, R_Gamma_loc
    real(8) :: dNe_shell, dNe_step, DB_step, Gam_e_max_step, Gam_e_m_step, Gam_e_m_p_step, source_norm, temp
    real(8) :: dR_base, dR_step, R_end, R_mid
    real(8) :: grid_content, dg_content
    real(8), parameter :: dg_base_substeps = 10d0, dg_jump_window_sigma = 4d0
    real(8), parameter :: dg_jump_substeps_per_sigma = 8d0, dg_jump_log_density_step = 5d-2
    integer :: I_tobs, Num_gam_rad

    allocate(dN_x_init(Num_gam_e), x_edge(Num_gam_e+1), gam_e_rad(Num_gam_e), dN_gam_e_rad(Num_gam_e), &
             source_grid(Num_gam_e))

    call electron_unpack_boundary(Boundary, n, Eta_0, R_ini, Epsilon_e, Epsilon_b, p, z, dNe_ISM, A_star, &
                                  E_iso, T_log10_duration, f_e, R_tr, f_jump, f_wide, R0)
    P_syn = zero
    Seed_syn = zero
    V_m = zero
    V_c = zero
    V_a = zero

    call electron_initial_density(A_star, dNe_ISM, R_ini, R(1), R0, dNe, Para_N_e_ini)
    DB = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1) - one)))
    Gam_e_max = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R) - one)))
    Gam_e_max_max = 3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1) - one)
    call electron_gamma_m_exact(p, temp_gam, Gam_e_max, Gam_e_m)
    Gam_e_c = 7.7d8/(one + dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    call electron_initialize_spectrum(Num_gam_e, Gam_e_max_max, Para_N_e_ini, p, Gam_e_m, Gam_e_c, Gam_e_max, &
                                      electron_initial_grid_log_edges, gam_e, dN_x_init, x_edge)
    call electron_profile_log_cell_edges(Num_gam_e, gam_e, x_edge)
    call electron_dnx_to_dndgamma_exp_centers(Num_gam_e, x_edge, gam_e, dN_x_init, dN_gam_e(:,1))

    call electron_dg1d_build_mesh(x_edge(1), x_edge(Num_gam_e+1), dlog10(Gam_e_m), dlog10(Gam_e_c), dlog10(Gam_e_max), mesh)
    allocate(state(mesh%ntot))
    call electron_dg1d_initial_state(mesh, Para_N_e_ini, p, Gam_e_m, Gam_e_c, Gam_e_max, state)
    call scale_dg_state_to_grid_content(state, dN_x_init)

    do I_tobs = 2, Num_R
        call prepare_shell(I_tobs)
        call write_radiation_and_breaks(I_tobs)
        call remesh_shell()

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

    deallocate(state, dN_x_init, x_edge, gam_e_rad, dN_gam_e_rad, source_grid)
    if (allocated(projected)) deallocate(projected)
    if (allocated(gamma_dg)) deallocate(gamma_dg, dEl_dg, dEl_dg_base, source_dg)

contains

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
        Gam_e_c = 7.7d8*(one + z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs_local)
        dDD = R(I_tobs_local) - R(I_tobs_local - 1)
    end subroutine prepare_shell

    subroutine write_radiation_and_breaks(I_tobs_local)
        integer, intent(in) :: I_tobs_local

        call electron_prepare_radiation_spectrum(Num_gam_e, gam_e, dN_gam_e(:,I_tobs_local - 1), &
                                                 Num_gam_rad, gam_e_rad, dN_gam_e_rad)
        V_m(I_tobs_local - 1) = 4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(one - beta_Gam)*(one + z))
        V_c(I_tobs_local - 1) = 4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(one - beta_Gam)*(one + z))
        call get_nu_a(R_loc, DB, Num_gam_rad, gam_e_rad(1:Num_gam_rad), dN_gam_e_rad(1:Num_gam_rad), temp)
        V_a(I_tobs_local - 1) = temp/(R_Gamma_loc*(one - beta_Gam)*(one + z))
        call get_syn_selected(index_syn_intger, R_loc, DB, Num_gam_e, Num_nu, n_threads, &
                              gam_e, dN_gam_e(:,I_tobs_local - 1), V_seed, P_syn(:,I_tobs_local), &
                              Seed_syn(:,I_tobs_local))
    end subroutine write_radiation_and_breaks

    subroutine remesh_shell()
        call electron_dg1d_build_mesh(x_edge(1), x_edge(Num_gam_e+1), dlog10(Gam_e_m), dlog10(Gam_e_c), &
                                      dlog10(Gam_e_max), new_mesh)
        call ensure_dg_work(new_mesh%ntot)
        call electron_dg1d_project_state(mesh, state, new_mesh, projected)
        deallocate(state)
        allocate(state(new_mesh%ntot))
        state = projected
        mesh = new_mesh
        call prepare_dg_cooling()
    end subroutine remesh_shell

    subroutine ensure_dg_work(ntot)
        integer, intent(in) :: ntot

        if (allocated(projected)) then
            if (size(projected) /= ntot) deallocate(projected, gamma_dg, dEl_dg, dEl_dg_base, source_dg)
        endif
        if (.not. allocated(projected)) allocate(projected(ntot), gamma_dg(ntot), dEl_dg(ntot), dEl_dg_base(ntot), source_dg(ntot))
    end subroutine ensure_dg_work

    subroutine prepare_dg_cooling()
        call electron_dg1d_gamma_nodes(mesh, gamma_dg)
        call get_forward_cooling(index_Y, Epsilon_e, Epsilon_b, p, DB, Gam_e_m, Gam_e_c, Gam_e_max, &
                                 R_loc, R_Gamma_loc, beta_Gam, dNe_shell, mesh%ntot, Num_nu, n_threads, &
                                 gamma_dg, V_seed, P_syn(:,I_tobs), Seed_syn(:,I_tobs), dEl_dg_base)
    end subroutine prepare_dg_cooling

    subroutine advance_substep(dR_local, R_step)
        real(8), intent(in) :: dR_local, R_step

        call dynamics_external_density_profile(A_star, dNe_ISM, R_step, R0, 1, R_tr, f_jump, f_wide, dNe_step)
        DB_step = 0.39d0*dsqrt(Epsilon_b*dNe_step*(R_Gamma_loc*(R_Gamma_loc - one)))
        Gam_e_max_step = 3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
        temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc - one)
        call electron_gamma_m_exact(p, temp_gam, Gam_e_max_step, Gam_e_m_step)
        Gam_e_m_p_step = (one - p)/(Gam_e_max_step**(one - p) - Gam_e_m_step**(one - p))
        call electron_injection_prefactor(R_step, dR_local, dNe_step, f_e, Gam_e_m_p_step, source_norm)

        call electron_dg1d_project_source(mesh, source_norm, p, Gam_e_m_step, Gam_e_max_step, source_dg)
        call electron_build_source_term_exp_cutoff_edges(Num_gam_e, x_edge, Gam_e_m_step, Gam_e_max_step, &
                                                         source_norm, p, source_grid)
        call scale_dg_state_to_grid_content(source_dg, source_grid)
        if (dNe_shell > zero) then
            dEl_dg = dEl_dg_base*(dNe_step/dNe_shell)
        else
            dEl_dg = dEl_dg_base
        endif
        call electron_dg1d_advance_step(mesh, one/R_step, dR_local, dEl_dg, source_dg, state, projected)
        call electron_dg1d_limit_positive_cell_preserving(mesh, projected)
        state = projected
    end subroutine advance_substep

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

    subroutine scale_dg_state_to_grid_content(dg_state, grid_state)
        real(8), intent(inout) :: dg_state(:)
        real(8), intent(in) :: grid_state(Num_gam_e)

        grid_content = sum(grid_state*(x_edge(2:Num_gam_e+1) - x_edge(1:Num_gam_e)))
        call electron_dg1d_integral(mesh, dg_state, dg_content)
        if (grid_content > zero) then
            if (dg_content <= zero) error stop 'fs_electron_dg_1d source projection has non-positive DG content'
            dg_state = dg_state*(grid_content/dg_content)
        endif
    end subroutine scale_dg_state_to_grid_content

    subroutine write_positive_output(I_tobs_local)
        integer, intent(in) :: I_tobs_local
        real(8) :: projected_content, dg_content_local

        call electron_dg1d_project_to_log_cells(mesh, state, Num_gam_e, x_edge, gam_e, dN_gam_e(:,I_tobs_local))
        ! Radiation-boundary positivity only: do not renormalize or feed this projection back into DG transport.
        where (dN_gam_e(:,I_tobs_local) > zero .and. ieee_is_finite(dN_gam_e(:,I_tobs_local)))
            dN_gam_e(:,I_tobs_local) = dN_gam_e(:,I_tobs_local)
        elsewhere
            dN_gam_e(:,I_tobs_local) = zero
        end where
        call electron_dg1d_integral(mesh, state, dg_content_local)
        projected_content = sum(dN_gam_e(:,I_tobs_local)*gam_e*dlog(ten)*(x_edge(2:Num_gam_e+1) - x_edge(1:Num_gam_e)))
        if (.not. (dg_content_local > zero .and. ieee_is_finite(dg_content_local))) &
            error stop 'fs_electron_dg_1d output projection has non-positive DG content'
        if (.not. (projected_content > zero .and. ieee_is_finite(projected_content))) &
            error stop 'fs_electron_dg_1d output projection lost positive content'
    end subroutine write_positive_output

end subroutine fs_electron_dg_1d
