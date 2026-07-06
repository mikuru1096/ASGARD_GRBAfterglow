! 正向激波电子 DG 求解器：移动多域 LGL 网格跟踪注入、冷却和高能尾部。
! Forward-shock electron DG solver: a moving multi-domain LGL mesh tracks injection, cooling, and the high-energy tail.
subroutine fs_dg_1d(Boundary, R_Tobs, R_Gamma, R, V_seed, n, Num_nu, Num_R, Num_gam_e, &
                             index_Y, index_syn_intger, n_threads, thermal_electrons, &
                             gam_e, dN_gam_e, P_syn, Seed_syn, V_m, V_c, V_a)
    use constants
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use dynamics_density_profile, only: density_profile, jump_count, &
                                        jump_factor, jump_radius, jump_width, &
                                        profile_count
    use electron_common, only: electron_initial_density, electron_unpack_boundary, electron_gm_exact, &
                               tail_factor, electron_injection_prefactor
    use electron_cooling_kernel, only: get_forward_cooling
    use electron_coord_common, only: build_fourvel_grid, fourvel_scale
    use electron_dg_transport, only: dg_mesh, dg_build_mesh, &
                                               dg_initial_state, dg_project_state, &
                                               dg_project_source, dg_project_cell_density, &
                                               dg_advance_step, dg_project_cells, &
                                               dg_limit_positive, dg_integral, &
                                               dg_tail_fraction, &
                                               dg_filter_positive
    use electron_shell_transport, only: coord_to_dgamma
    use electron_injection_profiles, only: init_coord, &
                                           source_coord
    use hybrid_spectrum, only: hybrid_coord
    use electron_radiation_kernel, only: syn_state, nua_fromtau
    implicit none

    integer, intent(in) :: n, Num_nu, Num_R, Num_gam_e, index_Y, index_syn_intger, n_threads
    integer, intent(in) :: thermal_electrons
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(out), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_gam_e,Num_R) :: dN_gam_e
    real(8), intent(out), dimension(Num_nu,Num_R) :: P_syn
    real(8), intent(out), dimension(Num_nu,Num_R) :: Seed_syn
    real(8), intent(out), dimension(Num_R) :: V_m,V_c,V_a

    type(dg_mesh) :: mesh, new_mesh
    real(8), allocatable, dimension(:) :: state,projected,gdg,deldg,delbase,srcdg,srctpl
    real(8), allocatable, dimension(:) :: nxinit,xedge,yedge,srcgrid,pemit,tau
    real(8) :: Eta_0, R_ini, Epsilon_e, Epsilon_b, p, z, nism, A_star, E_iso, tdur_log
    real(8) :: f_e, R_tr, f_jump, f_wide, R0, dNe, ninit, DB, DB_min, gemax0
    real(8) :: gemax, gm, gc, temp_gam, beta_Gam, dDD, R_loc, gloc
    real(8) :: dNe_shell, dNe_step, DB_step, gmax_step, gm_step, gmp_step, source_norm, temp
    real(8) :: gmax_inj, gm_inj, gmp_shell
    real(8) :: dR_base, dR_step, R_end, R_mid, dgscale, coord_scale
    real(8) :: cache_gm, cache_gmax
    real(8), parameter :: dg_substeps = 10d0, jump_sigma = 4d0
    real(8), parameter :: jump_nsigma = 8d0, jump_logstep = 5d-2
    real(8), parameter :: tail_thresh = 1d-10, tail_power = 2d0
    integer :: I_tobs, cache_n
    logical :: cache_ready, uniform_shell

    allocate(nxinit(Num_gam_e), xedge(Num_gam_e+1), yedge(Num_gam_e+1), srcgrid(Num_gam_e), &
             pemit(Num_nu), tau(Num_nu))

    call electron_unpack_boundary(Boundary, n, Eta_0, R_ini, Epsilon_e, Epsilon_b, p, z, nism, A_star, &
                                  E_iso, tdur_log, f_e, R_tr, f_jump, f_wide, R0)
    if (thermal_electrons /= 0) then
        if (f_e <= 0d0 .or. f_e > 1d0) error stop 'thermal electrons require 0 < f_e <= 1'
    endif
    P_syn = 0d0
    Seed_syn = 0d0
    V_m = 0d0
    V_c = 0d0
    V_a = 0d0
    cache_ready = .false.
    cache_n = 0
    cache_gm = -1d0
    cache_gmax = -1d0

    call electron_initial_density(A_star, nism, R_ini, R(1), R0, dNe, ninit)
    DB = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1) - 1d0)))
    gemax = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min = 0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R) - 1d0)))
    gemax0 = 3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1) - 1d0)
    call electron_gm_exact(p, temp_gam, gemax, gm)
    gc = 7.7d8/(1d0 + dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/2d0)
    call init_fourvel_grid()

    do I_tobs = 2, Num_R
        call prepare_shell(I_tobs)
        call write_rad_breaks(I_tobs)
        call remesh_shell(gemax, gm, gc, gemax)

        R_end = R(I_tobs)
        dR_base = dDD/dg_substeps
        do while (R_loc < R_end)
            dR_step = min(dR_base, R_end - R_loc)
            call limit_jump_step(R_loc, R_end, dR_step)
            R_mid = R_loc + 0.5d0*dR_step
            call advance_substep(dR_step, R_mid)
            R_loc = R_loc + dR_step
        enddo
        call write_positive_output(I_tobs)
    enddo
    call write_final()

    deallocate(state, nxinit, xedge, yedge, srcgrid, pemit, tau)
    if (allocated(projected)) deallocate(projected)
    if (allocated(gdg)) deallocate(gdg, deldg, delbase, srcdg, srctpl)

    contains

    ! 初始网格使用 four-velocity 坐标，低能端贴近 gamma=1，高能端覆盖初始 cutoff。
    ! The initial grid uses the four-velocity coordinate, starting near gamma=1 and covering the initial cutoff.
    subroutine init_fourvel_grid()
        dgscale = fourvel_scale
        coord_scale = dgscale*dgscale - 1d0
        call build_fourvel_grid(Num_gam_e, 1d0, tail_factor*gemax0, &
                                               dgscale, gam_e, yedge, xedge)
        call dg_build_mesh(xedge(1), dg_active_xmax(gemax), dlog(gm), &
                                                    dlog(gc), dlog(gemax), dgscale, mesh)
        allocate(state(mesh%ntot))
        if (thermal_electrons == 0) then
            call dg_initial_state(mesh, ninit, p, gm, gc, gemax, state)
            call init_coord(ninit, p, gm, gc, gemax, &
                                                                  Num_gam_e, yedge, coord_scale, nxinit)
        else
            call hybrid_coord(Num_gam_e, yedge, coord_scale, p, gm, gemax, f_e, nxinit)
            nxinit = nxinit*ninit
            call dg_project_cell_density(mesh, Num_gam_e, yedge, nxinit, state)
        endif
        call scale_dg_content(state, nxinit)
        call write_positive_output(1)
    end subroutine init_fourvel_grid

    subroutine prepare_shell(it)
        integer, intent(in) :: it

        R_loc = R(it - 1)
        gloc = (R_Gamma(it) + R_Gamma(it - 1))/2d0
        if (gloc < 1d0) error stop 'fs_dg_1d requires Gamma >= 1'
        beta_Gam = dsqrt(1d0 - 1d0/gloc**2)
        call density_profile(A_star, nism, R_loc, R0, 1, R_tr, f_jump, f_wide, dNe_shell)
        DB = 0.39d0*dsqrt(Epsilon_b*dNe_shell*(gloc*(gloc - 1d0)))
        gmax_inj = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(gloc - 1d0)
        call electron_gm_exact(p, temp_gam, gmax_inj, gm_inj)
        gmp_shell = (1d0 - p)/(gmax_inj**(1d0 - p) - gm_inj**(1d0 - p))
        uniform_shell = (A_star <= 0d0 .and. profile_count == 0 .and. &
                                 jump_count == 0 .and. abs(f_jump - 1d0) <= 0d0)
        gc = 7.7d8*(1d0 + z)/gloc/DB**2/R_Tobs(it)
        dDD = R(it) - R(it - 1)
    end subroutine prepare_shell

    ! 写当前壳层的同步辐射谱和三个 break frequency 诊断。
    ! Write the current shell synchrotron spectrum and the 3 break-frequency diagnostics.
    subroutine write_rad_breaks(it)
        integer, intent(in) :: it

        V_m(it - 1) = 4.2d6*DB*gm_inj*gm_inj/(gloc*(1d0 - beta_Gam)*(1d0 + z))
        V_c(it - 1) = 4.2d6*DB*gc*gc/(gloc*(1d0 - beta_Gam)*(1d0 + z))
        call syn_state(index_syn_intger, R_loc, DB, Num_gam_e, Num_nu, n_threads, &
                                    gam_e, dN_gam_e(:,it - 1), V_seed, pemit, &
                                    P_syn(:,it), Seed_syn(:,it), tau)
        call nua_fromtau(Num_nu, V_seed, tau, temp)
        V_a(it - 1) = temp/(gloc*(1d0 - beta_Gam)*(1d0 + z))
    end subroutine write_rad_breaks

    ! 最后一个输出点没有后续推进，只刷新与最终电子谱一致的辐射诊断。
    ! The final output point has no following advance, so refresh diagnostics from the final electron spectrum.
    subroutine write_final()
        R_loc = R(Num_R)
        gloc = R_Gamma(Num_R)
        if (gloc < 1d0) error stop 'fs_dg_1d requires Gamma >= 1'
        beta_Gam = dsqrt(1d0 - 1d0/gloc**2)
        call density_profile(A_star, nism, R_loc, R0, 1, R_tr, f_jump, f_wide, dNe_shell)
        DB = 0.39d0*dsqrt(Epsilon_b*dNe_shell*(gloc*(gloc - 1d0)))
        gemax = 3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(gloc - 1d0)
        call electron_gm_exact(p, temp_gam, gemax, gm)
        gc = 7.7d8*(1d0 + z)/gloc/DB**2/R_Tobs(Num_R)
        V_m(Num_R) = 4.2d6*DB*gm*gm/(gloc*(1d0 - beta_Gam)*(1d0 + z))
        V_c(Num_R) = 4.2d6*DB*gc*gc/(gloc*(1d0 - beta_Gam)*(1d0 + z))
        call syn_state(index_syn_intger, R_loc, DB, Num_gam_e, Num_nu, n_threads, &
                                    gam_e, dN_gam_e(:,Num_R), V_seed, pemit, &
                                    P_syn(:,Num_R), Seed_syn(:,Num_R), tau)
        call nua_fromtau(Num_nu, V_seed, tau, temp)
        V_a(Num_R) = temp/(gloc*(1d0 - beta_Gam)*(1d0 + z))
    end subroutine write_final

    ! 根据当前 cutoff 和尾部矩判断是否扩展 DG 网格，并守恒投影旧状态。
    ! Resize the DG mesh from the current cutoff and tail moment, then conservatively project the old state.
    subroutine remesh_shell(gup, gmbr, gcbr, gxbr)
        real(8), intent(in) :: gup, gmbr, gcbr, gxbr

        call dg_build_mesh(xedge(1), dg_active_xmax(gup), &
                                                    dlog(gmbr), dlog(gcbr), &
                                                    dlog(gxbr), dgscale, new_mesh)
        if (allocated(projected)) then
            if (size(projected) /= new_mesh%ntot) &
                deallocate(projected, gdg, deldg, delbase, srcdg, srctpl)
        endif
        if (.not. allocated(projected)) &
            allocate(projected(new_mesh%ntot), gdg(new_mesh%ntot), deldg(new_mesh%ntot), &
                     delbase(new_mesh%ntot), srcdg(new_mesh%ntot), srctpl(new_mesh%ntot))
        call dg_project_state(mesh, state, new_mesh, projected)
        call dg_filter_positive(new_mesh, projected)
        call dg_limit_positive(new_mesh, projected)
        if (size(state) /= new_mesh%ntot) then
            deallocate(state)
            allocate(state(new_mesh%ntot))
        endif
        state = projected
        mesh = new_mesh
        cache_ready = .false.
        gdg = mesh%gamma
        call get_forward_cooling(index_Y, Epsilon_e, Epsilon_b, p, DB, gm, gc, gemax, &
                                 R_loc, gloc, beta_Gam, dNe_shell, mesh%ntot, Num_nu, n_threads, &
                                 gdg, V_seed, P_syn(:,I_tobs), Seed_syn(:,I_tobs), delbase)

    end subroutine remesh_shell

    ! 单个半径子步：更新局部介质、注入源、冷却率，再推进 DG 状态。
    ! One radial substep: update local medium, source, cooling rate, then advance the DG state.
    subroutine advance_substep(dR_local, R_step)
        real(8), intent(in) :: dR_local, R_step

        if (uniform_shell) then
            dNe_step = dNe_shell
            DB_step = DB
            gmax_step = gmax_inj
            gm_step = gm_inj
            gmp_step = gmp_shell
        else
            call density_profile(A_star, nism, R_step, R0, 1, R_tr, f_jump, f_wide, dNe_step)
            DB_step = 0.39d0*dsqrt(Epsilon_b*dNe_step*(gloc*(gloc - 1d0)))
            gmax_step = 3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
            temp_gam = Epsilon_e/f_e*para_m_p/para_m_e*(gloc - 1d0)
            call electron_gm_exact(p, temp_gam, gmax_step, gm_step)
            gmp_step = (1d0 - p)/(gmax_step**(1d0 - p) - gm_step**(1d0 - p))
        endif
        if (thermal_electrons == 0) then
            call electron_injection_prefactor(R_step, dR_local, dNe_step, f_e, gmp_step, source_norm)
        else
            call electron_injection_prefactor(R_step, dR_local, dNe_step, 1d0, 1d0, source_norm)
        endif
        if (dg_source_xmax(gmax_step) > mesh%x_gamma(mesh%ntot)) &
            call remesh_shell(gmax_step, gm_step, gc, gmax_step)

        call prepare_dg_source()
        if (dNe_shell > 0d0) then
            deldg = delbase*(dNe_step/dNe_shell)
        else
            deldg = delbase
        endif
        call dg_advance_step(mesh, 1d0/R_step, dR_local, deldg, srcdg, state, projected)
        call dg_filter_positive(mesh, projected)
        call dg_limit_positive(mesh, projected)
        state = projected
    end subroutine advance_substep

    subroutine prepare_dg_source()
        if (source_norm <= 0d0) then
            srcdg = 0d0
            return
        endif
        if ((.not. cache_ready) .or. cache_n /= mesh%ntot .or. &
            cache_gm /= gm_step .or. cache_gmax /= gmax_step) then
            if (thermal_electrons == 0) then
                call dg_project_source(mesh, 1d0, p, gm_step, gmax_step, srctpl)
                call source_coord(Num_gam_e, yedge, coord_scale, &
                                                                       gm_step, gmax_step, 1d0, p, srcgrid)
            else
                call hybrid_coord(Num_gam_e, yedge, coord_scale, p, gm_step, gmax_step, f_e, srcgrid)
                call dg_project_cell_density(mesh, Num_gam_e, yedge, srcgrid, srctpl)
            endif
            call scale_dg_content(srctpl, srcgrid)
            cache_ready = .true.
            cache_n = mesh%ntot
            cache_gm = gm_step
            cache_gmax = gmax_step
        endif
        srcdg = source_norm*srctpl
    end subroutine prepare_dg_source

    subroutine limit_jump_step(R_left, R_stop, dR_limited)
        real(8), intent(in) :: R_left, R_stop
        real(8), intent(inout) :: dR_limited
        real(8) :: R_probe, slope
        integer :: j

        if (jump_count > 0) then
            do j = 1, jump_count
                call limit_one_jump(jump_radius(j), jump_factor(j), &
                                            jump_width(j), R_left, R_stop, dR_limited)
            enddo
        else if (abs(f_jump - 1d0) > 0d0) then
            call limit_one_jump(R_tr, f_jump, f_wide, R_left, R_stop, dR_limited)
        endif
        ! 用局部 Gaussian density-jump 的对数斜率限制子步，避免跨过跳变结构。
        ! Bound the substep by the local Gaussian density-jump logarithmic slope.
        R_probe = R_left + 0.5d0*dR_limited
        slope = jump_log_slope(R_probe)
        if (abs(slope) > 0d0) dR_limited = min(dR_limited, jump_logstep/abs(slope))
    end subroutine limit_jump_step

    subroutine limit_one_jump(R_jump, jump_factor, width_frac, R_left, R_stop, dR_limited)
        real(8), intent(in) :: R_jump, jump_factor, width_frac, R_left, R_stop
        real(8), intent(inout) :: dR_limited
        real(8), dimension(7) :: marks
        real(8) :: sigma_r,R_trial,eps_r
        integer :: i

        if (abs(jump_factor - 1d0) <= 0d0) return
        sigma_r = width_frac*R_jump
        marks = R_jump + sigma_r*[-4d0, -2d0, -1d0, 0d0, 1d0, 2d0, 4d0]
        R_trial = min(R_left + dR_limited, R_stop)
        eps_r = 1d-12*max(R_stop, 1d0)
        do i = 1, 7
            if (marks(i) > R_left + eps_r .and. marks(i) < R_trial + eps_r) then
                dR_limited = min(dR_limited, marks(i) - R_left)
                R_trial = R_left + dR_limited
            endif
        enddo
        if (R_left < R_jump + jump_sigma*sigma_r .and. &
            R_trial > R_jump - jump_sigma*sigma_r) then
            dR_limited = min(dR_limited, sigma_r/jump_nsigma)
        endif
    end subroutine limit_one_jump

    real(8) function jump_log_slope(radius) result(slope)
        real(8), intent(in) :: radius
        real(8) :: enhancement, derivative
        integer :: j

        enhancement = 1d0
        derivative = 0d0
        if (jump_count > 0) then
            do j = 1, jump_count
                call add_jump_deriv(radius, jump_radius(j), jump_factor(j), &
                                                 jump_width(j), enhancement, derivative)
            enddo
        else if (abs(f_jump - 1d0) > 0d0) then
            call add_jump_deriv(radius, R_tr, f_jump, f_wide, enhancement, derivative)
        endif
        slope = derivative/enhancement
    end function jump_log_slope

    subroutine add_jump_deriv(radius, R_jump, jump_factor, width_frac, enhancement, derivative)
        real(8), intent(in) :: radius, R_jump, jump_factor, width_frac
        real(8), intent(inout) :: enhancement, derivative
        real(8) :: sigma_r, profile

        if (abs(jump_factor - 1d0) <= 0d0) return
        sigma_r = width_frac*R_jump
        profile = (jump_factor - 1d0)*dexp(-((radius - R_jump)*(radius - R_jump))/(2d0*sigma_r*sigma_r))
        enhancement = enhancement + profile
        derivative = derivative - profile*(radius - R_jump)/(sigma_r*sigma_r)
    end subroutine add_jump_deriv

    real(8) function dg_source_xmax(gup) result(x_max)
        real(8), intent(in) :: gup
        real(8) :: grid_gmax

        grid_gmax = dexp(xedge(Num_gam_e+1))
        x_max = dlog(min(grid_gmax, tail_factor*gup))
    end function dg_source_xmax

    real(8) function dg_active_xmax(gup) result(x_max)
        real(8), intent(in) :: gup
        real(8) :: tail_fraction

        x_max = dg_source_xmax(gup)
        if (allocated(state)) then
            if (x_max < mesh%x_gamma(mesh%ntot)) then
                call dg_tail_fraction(mesh, state, dexp(x_max), &
                                                        tail_power, tail_fraction)
                if (tail_fraction > tail_thresh) x_max = mesh%x_gamma(mesh%ntot)
            endif
        endif
    end function dg_active_xmax

    subroutine scale_dg_content(dg_state, grid_state)
        real(8), intent(inout), dimension(:) :: dg_state
        real(8), intent(in), dimension(Num_gam_e) :: grid_state
        real(8) :: grid_content, dg_content

        grid_content = sum(grid_state*(yedge(2:Num_gam_e+1) - yedge(1:Num_gam_e)))
        call dg_integral(mesh, dg_state, dg_content)
        if (grid_content > 0d0) then
            if (dg_content <= 0d0) error stop 'fs_dg_1d source projection has non-positive DG content'
            dg_state = dg_state*(grid_content/dg_content)
        endif
    end subroutine scale_dg_content

    subroutine write_positive_output(it)
        integer, intent(in) :: it
        real(8) :: proj, dgcontent

        call dg_project_cells(mesh, state, Num_gam_e, yedge, srcgrid)
        call coord_to_dgamma(Num_gam_e, yedge, coord_scale, gam_e, &
                                                           srcgrid, dN_gam_e(:,it))
        where (dN_gam_e(:,it) <= 0d0 .or. .not. ieee_is_finite(dN_gam_e(:,it)))
            dN_gam_e(:,it) = 0d0
        end where
        call dg_integral(mesh, state, dgcontent)
        proj = sum(dN_gam_e(:,it)*gam_e*(xedge(2:Num_gam_e+1) - xedge(1:Num_gam_e)))
        if (.not. (dgcontent > 0d0 .and. ieee_is_finite(dgcontent))) &
            error stop 'fs_dg_1d output projection has non-positive DG content'
        if (.not. (proj > 0d0 .and. ieee_is_finite(proj))) &
            error stop 'fs_dg_1d output projection lost positive content'
    end subroutine write_positive_output

end subroutine fs_dg_1d
