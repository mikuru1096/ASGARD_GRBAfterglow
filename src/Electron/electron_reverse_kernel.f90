module electron_reverse_kernel
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    use electron_energy_coordinate_common, only: electron_build_four_velocity_grid, electron_coord_from_xgamma, &
                                                 electron_dxgamma_dcoord, &
                                                 electron_coord_log_four_velocity_sq, &
                                                 electron_four_velocity_grid_gamma_scale
    use electron_injection_profiles, only: electron_build_kinetic_source_term_exp_cutoff_edges, &
                                           electron_build_kinetic_source_term_exp_cutoff_coord_edges, &
                                           electron_profile_log_cell_edges
    use electron_common, only: electron_exp_tail_grid_factor
    use electron_shell_transport_common, only: electron_resolve_1d_solver_id, &
                                               electron_shell_flux_split_coord_step, &
                                               electron_shell_dcoord_to_dndgamma_exp_centers, &
                                               electron_solver_dg_1d, electron_solver_fullhide_1d
    use electron_transport_common, only: electron_build_piecewise_affine_u, electron_characteristic_remap_edges, &
                                         electron_dnx_to_dndgamma_exp_centers, electron_trace_affine_u_edges_batch, &
                                         electron_trace_piecewise_affine_u_edges_batch, electron_u_edges_from_x, &
                                         electron_fullhide_flux_split_step
    use electron_transport_dg_1d_kernel, only: electron_dg1d_mesh, electron_dg1d_advance_characteristic_step, &
                                               electron_dg1d_advance_step, &
                                               electron_dg1d_build_four_velocity_mesh, electron_dg1d_integral, &
                                               electron_dg1d_limit_positive_cell_preserving, &
                                               electron_dg1d_project_kinetic_source, &
                                               electron_dg1d_project_state, electron_dg1d_project_to_coord_cells, &
                                               electron_dg1d_scale_to_content, electron_dg1d_tail_moment_fraction, &
                                               electron_dg1d_apply_positive_kernel_filter
    use electron_radiation_kernel, only: get_syn_selected_state, get_nu_a
    use electron_cooling_kernel, only: prepare_forward_cooling_aux_batch
    implicit none
    integer, parameter :: reverse_dg_base_substeps = 10
contains

! 反向激波电子演化主驱动：注入→同步+IC冷却→隐式输运推进，IC模式遵循公开标准选择。
    subroutine electron_reverse_evolve(Delta_0,e_r,b_r,p_r,f_e_r,eta_0,Epsilon_e,Epsilon_b,z,A_star,dNe_ISM,para_m_ej, &
                                       R_tr,f_jump,f_wide,R0, &
                                       T_cross,R_cross,U3_cross,V3_cross,M3_cross,R_Tobs,R_Gamma,R,B3,M3_shell, &
                                       U3_shell,V3_shell,V_seed, &
                                       Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads,gam_e,dN_gam_e, &
                                       solver_id)
    implicit none
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in), optional :: solver_id
    integer :: I_tobs,I_gam_e,L1,L,active_solver,i_empty,i_edge,i_coord
    real(8), intent(in) :: Delta_0,e_r,b_r,p_r,f_e_r,eta_0,Epsilon_e,Epsilon_b,z,A_star,dNe_ISM,para_m_ej
    real(8), intent(in) :: R_tr,f_jump,f_wide,R0
    real(8), intent(in) :: T_cross,R_cross,U3_cross,V3_cross,M3_cross
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),B3(Num_R),M3_shell(Num_R)
    real(8), intent(in) :: U3_shell(Num_R),V3_shell(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R)
    real(8), parameter :: reverse_gamma_c_coeff=7.7d8, reverse_synch_b_coeff=0.39d0, reverse_adv_coeff=1.35d-19
    real(8), parameter :: reverse_transrel_gamma_break=2d0
    real(8), parameter :: reverse_dg_tail_moment_threshold=1d-10, reverse_dg_tail_moment_power=2d0
    real(8) :: factor2,dB,gamma34,Gam_e_max,Gam_e_m,Gam_e_c,dNe,DB_min,Gam_e_max_max,Gam_e_min_global,d_x,R_loc,R_Gamma_loc,Delta
    real(8) :: R_n4,beta4,beta2,u2,u4,f_r,dDR,dDD,Qshell,cooling_scale,R_step_mid
    real(8) :: thermal_scale_lo,thermal_scale_hi,thermal_loss_rate,adiabatic_rate,dg_gamma_scale,coord_scale
    real(8) :: coord_mid,dxdy
    real(8) :: dg_gamma_low,dg_gamma_m_front,dg_gamma_high
    real(8) :: injection_rate,inj_hi,inj_width,mass_lo,mass_hi
    real(8), allocatable :: dEl(:),x(:),dF1(:),temp3(:),dN_x(:),x_edge(:),coord_edge(:)
    real(8), allocatable :: post_cross_log_state(:),post_cross_edge_map(:),post_cross_edge_work(:)
    real(8), allocatable :: post_cross_step_back(:),post_cross_u_edge(:),post_cross_a_cell(:),post_cross_b_cell(:)
    real(8), allocatable :: post_cross_affine_back(:,:)
    real(8), allocatable :: dB3_serial(:),P_syn(:),Seed_syn(:),cooling_aux(:),Compton(:)
    type(electron_dg1d_mesh) :: dg_mesh,dg_new_mesh
    real(8), allocatable :: dg_state(:),dg_work(:),dg_dEl(:),dg_source(:)
    logical :: crossed_reverse,post_cross_cache_ready

    allocate(dEl(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e),temp3(Num_gam_e-1),dN_x(Num_gam_e), &
             x_edge(Num_gam_e+1),coord_edge(Num_gam_e+1), &
             post_cross_log_state(Num_gam_e),post_cross_edge_map(Num_gam_e+1), &
             post_cross_edge_work(Num_gam_e+1),post_cross_step_back(Num_gam_e+1), &
             post_cross_u_edge(Num_gam_e+1),post_cross_a_cell(Num_gam_e),post_cross_b_cell(Num_gam_e), &
             post_cross_affine_back(Num_gam_e+1,1), &
             dB3_serial(Num_R),P_syn(Num_nu),Seed_syn(Num_nu), &
             cooling_aux(Num_gam_e),Compton(Num_gam_e))
    dB3_serial=B3

    active_solver=electron_resolve_1d_solver_id(solver_id)
    crossed_reverse=(T_cross > zero .and. R_cross > zero .and. M3_cross > zero .and. V3_cross > zero)
    post_cross_cache_ready=.false.
    if (maxval(B3) <= zero .or. maxval(M3_shell) <= M3_shell(1)) then
        do i_empty=1,Num_gam_e
            if (Num_gam_e == 1) then
                gam_e(i_empty)=one
            else
                gam_e(i_empty)=ten**(dble(i_empty-1)/dble(Num_gam_e-1))
            end if
        end do
        dN_gam_e=zero
        deallocate(dEl,x,dF1,temp3,dN_x,x_edge,coord_edge,post_cross_log_state,post_cross_edge_map, &
                   post_cross_edge_work,post_cross_step_back,post_cross_u_edge,post_cross_a_cell, &
                   post_cross_b_cell,post_cross_affine_back,dB3_serial,P_syn,Seed_syn,cooling_aux,Compton)
        return
    end if

    factor2=(p_r-two)/(p_r-one)*e_r*Para_m_p_div_m_e
    if (p_r < 2.05d0) factor2=0.05d0/1.05d0*e_r*Para_m_p_div_m_e
    beta4=dsqrt(one-one/eta_0**2); u4=dsqrt(eta_0*eta_0-one)
    dB3_serial(1)=dB3_serial(min(2,Num_R))
    dB=dB3_serial(1); gamma34=1.001d0
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
    Gam_e_m=factor2*(gamma34-one)/f_e_r+one
    dg_gamma_scale=electron_four_velocity_grid_gamma_scale
    dg_gamma_low=one
    dg_gamma_m_front=Gam_e_m
    dg_gamma_high=Gam_e_max
    Gam_e_c=reverse_gamma_c_coeff/(one+dsqrt(e_r/b_r))/R_Gamma(1)/dB**2/(R_Tobs(1)/two)

    call dynamics_external_density_profile(A_star,dNe_ISM,R(1),R0,1,R_tr,f_jump,f_wide,dNe)
    DB_min=reverse_synch_b_coeff*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    Gam_e_min_global=one

    do I_tobs=2,Num_R
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        beta2=dsqrt(one-one/R_Gamma_loc**2)
        u2=dsqrt(R_Gamma_loc*R_Gamma_loc-one)
        gamma34=(R_Gamma_loc*R_Gamma_loc+eta_0*eta_0-one)/(eta_0*R_Gamma_loc+u2*u4)
        dB=(dB3_serial(I_tobs)+dB3_serial(I_tobs-1))/two
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        Gam_e_m=factor2*(gamma34-one)/f_e_r+one
    end do
    if (Gam_e_max_max <= Gam_e_min_global) error stop "electron_reverse_evolve: reverse electron grid maximum must exceed minimum."

    do I_gam_e=1,Num_gam_e
        if (Num_gam_e == 1) then
            gam_e(I_gam_e)=Gam_e_min_global
        else
            gam_e(I_gam_e)=Gam_e_min_global*ten**(dlog10(Gam_e_max_max/Gam_e_min_global)*(I_gam_e-1)/(Num_gam_e-1))
        end if
        dN_gam_e(I_gam_e,1)=zero
    end do

    coord_scale=dg_gamma_scale*dg_gamma_scale-one
    call electron_build_four_velocity_grid(Num_gam_e,Gam_e_min_global,Gam_e_max_max,dg_gamma_scale,gam_e,coord_edge,x_edge)
    dN_x=zero
    d_x=coord_edge(2)-coord_edge(1)
    R_loc=R(1)
    if (active_solver == electron_solver_dg_1d) call initialize_reverse_dg_state()

    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        Delta=max(Delta_0,R_loc/Eta_0**2)
        R_n4=para_m_ej/(4d0*pi*Para_m_p*R_loc*R_loc*Eta_0*Delta)
        beta2=dsqrt(one-one/R_Gamma_loc**2)
        u2=dsqrt(R_Gamma_loc*R_Gamma_loc-one)
        gamma34=(R_Gamma_loc*R_Gamma_loc+eta_0*eta_0-one)/(eta_0*R_Gamma_loc+u2*u4)
        dB=(dB3_serial(I_tobs)+dB3_serial(I_tobs-1))/two
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        Gam_e_m=factor2*(gamma34-one)/f_e_r+one
        Gam_e_c=reverse_gamma_c_coeff*(one+z)/R_Gamma_loc/dB**2/R_Tobs(I_tobs)
        f_r=reverse_adv_coeff/beta2/R_Gamma_loc*dB**2/pi
        dDR=0.7d0/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        if (dB <= zero) then
            dN_gam_e(:,I_tobs)=dN_gam_e(:,I_tobs-1)
            cycle
        end if
        thermal_loss_rate=zero
        if (crossed_reverse .and. R(I_tobs) > R_cross) then
            if (R(I_tobs-1) < R_cross) then
                thermal_scale_lo=V3_cross
                thermal_scale_hi=V3_shell(I_tobs)
            else
                thermal_scale_lo=V3_shell(I_tobs-1)
                thermal_scale_hi=V3_shell(I_tobs)
            end if
            if (thermal_scale_lo <= zero .or. thermal_scale_hi <= zero) &
                error stop "electron_reverse_evolve: post-crossing comoving volume must be positive."
            if (R(I_tobs-1) < R_cross) then
                thermal_loss_rate=dlog(thermal_scale_hi/thermal_scale_lo)/(3d0*(R(I_tobs)-R_cross))
            else
                thermal_loss_rate=dlog(thermal_scale_hi/thermal_scale_lo)/(3d0*dDD)
            end if
            if (thermal_loss_rate < zero) &
                error stop "electron_reverse_evolve: post-crossing comoving volume must not decrease."
        end if
        injection_rate=zero
        if ((.not. crossed_reverse) .or. R(I_tobs-1) < R_cross) then
            if (crossed_reverse) then
                inj_hi=min(R(I_tobs),R_cross)
            else
                inj_hi=R(I_tobs)
            end if
            if (inj_hi > R(I_tobs-1)) then
                mass_lo=M3_shell(I_tobs-1)
                mass_hi=M3_shell(I_tobs)
                if (crossed_reverse .and. R(I_tobs) > R_cross) mass_hi=M3_cross
                if (mass_hi < mass_lo) error stop "electron_reverse_evolve: reverse swept mass must not decrease."
                inj_width=inj_hi-R(I_tobs-1)
                injection_rate=f_e_r*(mass_hi-mass_lo)/(Para_m_p*inj_width)
            end if
        end if
        L1=reverse_transport_substeps(dDR,dDD,active_solver)
        dDR=dDD/L1
        do i_coord=1,Num_gam_e
            coord_mid=0.5d0*(coord_edge(i_coord)+coord_edge(i_coord+1))
            dxdy=electron_dxgamma_dcoord(electron_coord_log_four_velocity_sq,coord_scale,coord_mid)
            dN_x(i_coord)=dN_gam_e(i_coord,I_tobs-1)*gam_e(i_coord)*dlog(ten)*dxdy
        end do

        call compute_reverse_cooling_loss(I_tobs)
        call advance_reverse_transport_shell(I_tobs)
    end do

    deallocate(dEl,x,dF1,temp3,dN_x,x_edge,coord_edge,post_cross_log_state,post_cross_edge_map, &
               post_cross_edge_work,post_cross_step_back,post_cross_u_edge,post_cross_a_cell, &
               post_cross_b_cell,post_cross_affine_back,dB3_serial,P_syn,Seed_syn,cooling_aux,Compton)
    if (allocated(dg_state)) deallocate(dg_state)
    if (allocated(dg_work)) deallocate(dg_work,dg_dEl,dg_source)

contains

    subroutine compute_reverse_cooling_loss(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs
    real(8) :: P_syn_column(Num_nu,1),Seed_syn_column(Num_nu,1),cooling_aux_column(Num_gam_e,1)
    real(8) :: P_emit_tmp(Num_nu),Tau_syn_tmp(Num_nu)

        call get_syn_selected_state(index_syn_intger,R(I_tobs-1),dB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_emit_tmp,P_syn,Seed_syn,Tau_syn_tmp)
        P_syn_column(:,1)=P_syn
        Seed_syn_column(:,1)=Seed_syn
        call prepare_forward_cooling_aux_batch(index_Y,Num_gam_e,Num_nu,1,n_threads,gam_e,V_seed, &
                                               P_syn_column,Seed_syn_column,cooling_aux_column)
        cooling_aux=cooling_aux_column(:,1)
        select case(index_Y)
        case(0)
            dEl=f_r*gam_e
        case(1)
            cooling_scale=one/(beta2*R_Gamma_loc*Para_c)
            dEl=(f_r+cooling_aux*cooling_scale)*gam_e
        case(2)
            Qshell=4d0*pi*R(I_tobs-1)*R(I_tobs-1)*Para_c
            Compton=one+cooling_aux/Qshell/(4d0*R_Gamma_loc*R_Gamma_loc*R_n4*Para_m_p_E)
            Gam_e_max=Gam_e_max/dsqrt(Compton(Num_gam_e))
            dEl=f_r*Compton*gam_e
        case default
            error stop 'electron_reverse_evolve: index_Y must be 0, 1, or 2.'
        end select
    end subroutine compute_reverse_cooling_loss

    subroutine advance_reverse_transport_shell(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs
    real(8) :: source_norm_step
    logical :: pure_post_cross_shell

        pure_post_cross_shell = (crossed_reverse .and. R(I_tobs - 1) >= R_cross)
        if (pure_post_cross_shell) then
            call advance_reverse_post_cross_analytic(I_tobs)
            return
        end if
        if (active_solver == electron_solver_dg_1d) then
            call remesh_reverse_dg_state()
            if (index_Y /= 0) call prepare_reverse_dg_cooling()
            do L=1,L1
                if (index_Y == 0) then
                    R_step_mid=R_loc+0.5d0*dDR
                    call prepare_reverse_transport_substep_state(I_tobs,R_step_mid)
                else
                    R_step_mid=R_loc+dDR
                end if
                if ((.not. crossed_reverse) .or. R_cross >= R_step_mid) then
                    source_norm_step=injection_rate
                else
                    source_norm_step=zero
                end if
                if ((.not. crossed_reverse) .or. R_step_mid <= R_cross) then
                    adiabatic_rate=one/R_step_mid
                else
                    adiabatic_rate=thermal_loss_rate
                end if
                if (source_norm_step > zero) then
                    call remesh_reverse_dg_state()
                    if (index_Y /= 0) call prepare_reverse_dg_cooling()
                end if
                if (source_norm_step > zero) then
                    call electron_dg1d_project_kinetic_source(dg_mesh,source_norm_step,p_r,Gam_e_m,Gam_e_max,dg_source)
                else
                    dg_source=zero
                end if
                call electron_dg1d_scale_to_content(dg_mesh,source_norm_step,dg_source)
                if (index_Y == 0) then
                    call electron_dg1d_advance_characteristic_step(dg_mesh,dDR,f_r,adiabatic_rate,dg_source, &
                                                                   dg_state,dg_work)
                else
                    call electron_dg1d_advance_step(dg_mesh,adiabatic_rate,dDR,dg_dEl,dg_source,dg_state,dg_work)
                    call electron_dg1d_limit_positive_cell_preserving(dg_mesh,dg_work)
                end if
                call electron_dg1d_apply_positive_kernel_filter(dg_mesh,dg_work)
                call electron_dg1d_limit_positive_cell_preserving(dg_mesh,dg_work)
                dg_state=dg_work
                call advance_reverse_dg_low_front(source_norm_step)
                call advance_reverse_dg_injection_front(source_norm_step)
                call advance_reverse_dg_high_front(source_norm_step)
                R_loc=R_loc+dDR
            end do
            call electron_dg1d_project_to_coord_cells(dg_mesh,dg_state,Num_gam_e,coord_edge,dF1)
            call electron_shell_dcoord_to_dndgamma_exp_centers(Num_gam_e,coord_edge,coord_scale, &
                                                               gam_e,dF1,dN_gam_e(:,I_tobs))
            where (dN_gam_e(:,I_tobs) > zero)
                dN_gam_e(:,I_tobs)=dN_gam_e(:,I_tobs)
            elsewhere
                dN_gam_e(:,I_tobs)=zero
            end where
            dN_x=dN_gam_e(:,I_tobs)*gam_e*dlog(ten)
            return
        end if

        do L=1,L1
            if (index_Y == 0) then
                R_step_mid=R_loc+0.5d0*dDR
                call prepare_reverse_transport_substep_state(I_tobs,R_step_mid)
            else
                R_loc=R_loc+dDR
                R_step_mid=R_loc
            end if
            if ((.not. crossed_reverse) .or. R_cross >= R_step_mid) then
                call electron_build_kinetic_source_term_exp_cutoff_coord_edges(Num_gam_e,coord_edge,coord_scale, &
                                                                               Gam_e_m,Gam_e_max,injection_rate,p_r,dF1)
            else
                dF1=zero
            end if
            if ((.not. crossed_reverse) .or. R_step_mid <= R_cross) then
                adiabatic_rate=one/R_step_mid
            else
                adiabatic_rate=thermal_loss_rate
            end if
            call electron_shell_flux_split_coord_step(Num_gam_e,dDR,coord_edge,coord_scale,dEl,adiabatic_rate,dF1,dN_x,x)
            dN_x=x
            if (index_Y == 0) R_loc=R_loc+dDR
            if (L == L1) call electron_shell_dcoord_to_dndgamma_exp_centers(Num_gam_e,coord_edge,coord_scale, &
                                                                            gam_e,dN_x,dN_gam_e(:,I_tobs))
        end do
    end subroutine advance_reverse_transport_shell

    subroutine advance_reverse_post_cross_analytic(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs

        dF1=zero
        if (.not. post_cross_cache_ready) then
            post_cross_log_state=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)
            post_cross_edge_map=x_edge
            post_cross_cache_ready=.true.
        end if
        do L=1,L1
            R_step_mid=R_loc+0.5d0*dDR
            if (index_Y == 0) then
                call prepare_reverse_transport_substep_state(I_tobs,R_step_mid)
                call electron_u_edges_from_x(Num_gam_e,x_edge,post_cross_u_edge)
                call electron_trace_affine_u_edges_batch(Num_gam_e,1,post_cross_u_edge,(/dDR/), &
                                                         f_r,thermal_loss_rate,post_cross_affine_back)
                post_cross_step_back=post_cross_affine_back(:,1)
            else
                call electron_build_piecewise_affine_u(Num_gam_e,x_edge,gam_e,dEl,thermal_loss_rate, &
                                                       post_cross_u_edge,post_cross_a_cell,post_cross_b_cell)
                call electron_trace_piecewise_affine_u_edges_batch(Num_gam_e,1,post_cross_u_edge,post_cross_u_edge, &
                                                                   post_cross_a_cell,post_cross_b_cell,(/dDR/), &
                                                                   post_cross_affine_back)
                post_cross_step_back=post_cross_affine_back(:,1)
            end if
            do i_edge=1,Num_gam_e+1
                post_cross_edge_work(i_edge)=reverse_post_cross_map_value(post_cross_step_back(i_edge))
            end do
            post_cross_edge_map=post_cross_edge_work
            R_loc=R_loc+dDR
        end do
        call electron_characteristic_remap_edges(Num_gam_e,x_edge,post_cross_edge_map,post_cross_log_state,dN_x,.true.)
        call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
        where (dN_gam_e(:,I_tobs) > zero)
            dN_gam_e(:,I_tobs)=dN_gam_e(:,I_tobs)
        elsewhere
            dN_gam_e(:,I_tobs)=zero
        end where
    end subroutine advance_reverse_post_cross_analytic

    real(8) function reverse_post_cross_map_value(x_eval) result(x_cross)
    implicit none
    integer :: left,right,mid
    real(8), intent(in) :: x_eval
    real(8) :: frac

        if (x_eval <= x_edge(1)) then
            frac=(x_eval-x_edge(1))/(x_edge(2)-x_edge(1))
            x_cross=post_cross_edge_map(1)+frac*(post_cross_edge_map(2)-post_cross_edge_map(1))
            return
        end if
        if (x_eval >= x_edge(Num_gam_e+1)) then
            frac=(x_eval-x_edge(Num_gam_e+1))/(x_edge(Num_gam_e+1)-x_edge(Num_gam_e))
            x_cross=post_cross_edge_map(Num_gam_e+1)+ &
                    frac*(post_cross_edge_map(Num_gam_e+1)-post_cross_edge_map(Num_gam_e))
            return
        end if
        left=1
        right=Num_gam_e
        do while (left < right)
            mid=(left+right)/2
            if (x_edge(mid+1) >= x_eval) then
                right=mid
            else
                left=mid+1
            end if
        end do
        frac=(x_eval-x_edge(left))/(x_edge(left+1)-x_edge(left))
        x_cross=post_cross_edge_map(left)+frac*(post_cross_edge_map(left+1)-post_cross_edge_map(left))
    end function reverse_post_cross_map_value

    subroutine initialize_reverse_dg_state()
    implicit none

        call electron_dg1d_build_four_velocity_mesh(x_edge(1),reverse_dg_active_xmax(),dlog10(reverse_dg_low_break()), &
                                                    dlog10(reverse_dg_injection_break()), &
                                                    dlog10(reverse_dg_upper_break()),dg_gamma_scale,dg_mesh)
        allocate(dg_state(dg_mesh%ntot))
        dg_state=zero
    end subroutine initialize_reverse_dg_state

    subroutine remesh_reverse_dg_state()
    implicit none

        call electron_dg1d_build_four_velocity_mesh(x_edge(1),reverse_dg_active_xmax(),dlog10(reverse_dg_low_break()), &
                                                    dlog10(reverse_dg_injection_break()), &
                                                    dlog10(reverse_dg_upper_break()),dg_gamma_scale,dg_new_mesh)
        call ensure_reverse_dg_work(dg_new_mesh%ntot)
        call electron_dg1d_project_state(dg_mesh,dg_state,dg_new_mesh,dg_work)
        call electron_dg1d_apply_positive_kernel_filter(dg_new_mesh,dg_work)
        call electron_dg1d_limit_positive_cell_preserving(dg_new_mesh,dg_work)
        if (size(dg_state) /= dg_new_mesh%ntot) then
            deallocate(dg_state)
            allocate(dg_state(dg_new_mesh%ntot))
        end if
        dg_state=dg_work
        dg_mesh=dg_new_mesh
    end subroutine remesh_reverse_dg_state

    subroutine ensure_reverse_dg_work(ntot)
    implicit none
    integer, intent(in) :: ntot

        if (allocated(dg_work)) then
            if (size(dg_work) /= ntot) deallocate(dg_work,dg_dEl,dg_source)
        end if
        if (.not. allocated(dg_work)) allocate(dg_work(ntot),dg_dEl(ntot),dg_source(ntot))
    end subroutine ensure_reverse_dg_work

    subroutine prepare_reverse_dg_cooling()
    implicit none
    integer :: i_node

        do i_node=1,dg_mesh%ntot
            dg_dEl(i_node)=reverse_interp_log_grid(Num_gam_e,x_edge,dEl,dg_mesh%x_gamma(i_node))
        end do
    end subroutine prepare_reverse_dg_cooling

    subroutine prepare_reverse_transport_substep_state(I_tobs,radius_eval)
    implicit none
    integer, intent(in) :: I_tobs
    real(8), intent(in) :: radius_eval

        R_Gamma_loc=reverse_shell_linear_value(I_tobs,R_Gamma,radius_eval)
        beta2=dsqrt(one-one/R_Gamma_loc**2)
        u2=dsqrt(R_Gamma_loc*R_Gamma_loc-one)
        gamma34=(R_Gamma_loc*R_Gamma_loc+eta_0*eta_0-one)/(eta_0*R_Gamma_loc+u2*u4)
        dB=reverse_shell_linear_value(I_tobs,dB3_serial,radius_eval)
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        Gam_e_m=factor2*(gamma34-one)/f_e_r+one
        f_r=reverse_adv_coeff/beta2/R_Gamma_loc*dB**2/pi
        dEl=f_r*gam_e
    end subroutine prepare_reverse_transport_substep_state

    real(8) function reverse_shell_linear_value(i_shell,values,radius_eval) result(value)
    implicit none
    integer, intent(in) :: i_shell
    real(8), intent(in) :: values(Num_R),radius_eval
    real(8) :: frac

        frac=(radius_eval-R(i_shell-1))/(R(i_shell)-R(i_shell-1))
        value=values(i_shell-1)+frac*(values(i_shell)-values(i_shell-1))
    end function reverse_shell_linear_value

    real(8) function reverse_dg_upper_break() result(gamma_break)
    implicit none

        if ((.not. crossed_reverse) .or. R_loc < R_cross) then
            gamma_break=Gam_e_max
        else
            gamma_break=dg_gamma_high
        end if
    end function reverse_dg_upper_break

    real(8) function reverse_dg_source_upper_xmax() result(x_max)
    implicit none
    real(8) :: gamma_grid_max

        gamma_grid_max=ten**x_edge(Num_gam_e+1)
        x_max=dlog10(min(gamma_grid_max,electron_exp_tail_grid_factor*reverse_dg_upper_break()))
    end function reverse_dg_source_upper_xmax

    real(8) function reverse_dg_active_xmax() result(x_max)
    implicit none
    real(8) :: tail_fraction

        x_max=reverse_dg_source_upper_xmax()
        if (allocated(dg_state)) then
            if (x_max < dg_mesh%x_gamma(dg_mesh%ntot)) then
                call electron_dg1d_tail_moment_fraction(dg_mesh,dg_state,ten**x_max, &
                                                        reverse_dg_tail_moment_power,tail_fraction)
                if (tail_fraction > reverse_dg_tail_moment_threshold) x_max=dg_mesh%x_gamma(dg_mesh%ntot)
            end if
        end if
    end function reverse_dg_active_xmax

    subroutine advance_reverse_dg_high_front(source_norm)
    implicit none
    real(8), intent(in) :: source_norm

        if (dg_gamma_high > one) call advance_reverse_dg_front_value(dg_gamma_high)
        if (source_norm > zero) dg_gamma_high=max(dg_gamma_high,Gam_e_max)
    end subroutine advance_reverse_dg_high_front

    subroutine advance_reverse_dg_low_front(source_norm)
    implicit none
    real(8), intent(in) :: source_norm

        if (source_norm > zero .and. dg_gamma_low <= one) dg_gamma_low=Gam_e_m
        if (dg_gamma_low > one) call advance_reverse_dg_front_value(dg_gamma_low)
        if (source_norm > zero) dg_gamma_low=min(dg_gamma_low,Gam_e_m)
        dg_gamma_low=max(one,dg_gamma_low)
    end subroutine advance_reverse_dg_low_front

    real(8) function reverse_dg_low_break() result(gamma_break)
    implicit none

        gamma_break=one
        if (dg_gamma_low > one) gamma_break=min(dg_gamma_low,Gam_e_m)
        ! Resolve the trans-relativistic kinetic-source singularity when gamma_m is nearly 1.
        if (Gam_e_m < reverse_transrel_gamma_break .and. Gam_e_max > reverse_transrel_gamma_break) &
            gamma_break=reverse_transrel_gamma_break
    end function reverse_dg_low_break

    subroutine advance_reverse_dg_injection_front(source_norm)
    implicit none
    real(8), intent(in) :: source_norm

        if (source_norm > zero) then
            dg_gamma_m_front=Gam_e_m
        else if (dg_gamma_m_front > one) then
            call advance_reverse_dg_front_value(dg_gamma_m_front)
        end if
        dg_gamma_m_front=max(one,dg_gamma_m_front)
    end subroutine advance_reverse_dg_injection_front

    real(8) function reverse_dg_injection_break() result(gamma_break)
    implicit none

        gamma_break=Gam_e_m
        if (crossed_reverse .and. R_loc >= R_cross .and. dg_gamma_m_front > one) gamma_break=dg_gamma_m_front
    end function reverse_dg_injection_break

    subroutine advance_reverse_dg_front_value(gamma_front)
    implicit none
    real(8), intent(inout) :: gamma_front
    real(8) :: exp_b,loss_front

        if (index_Y == 0) then
            if (abs(adiabatic_rate*dDR) <= 1d-10) then
                gamma_front=gamma_front/(one+f_r*gamma_front*dDR)
            else
                exp_b=dexp(-adiabatic_rate*dDR)
                gamma_front=gamma_front*exp_b/(one+(f_r/adiabatic_rate)*gamma_front*(one-exp_b))
            end if
        else
            loss_front=reverse_interp_log_grid(Num_gam_e,x_edge,dEl,dlog10(gamma_front))
            gamma_front=gamma_front*dexp(-(loss_front+adiabatic_rate)*dDR)
        end if
    end subroutine advance_reverse_dg_front_value
end subroutine electron_reverse_evolve

subroutine electron_secondary_reverse_evolve(e_r,b_r,p_r,f_e_r,z,R_Tobs,R_Gamma,R,B3,M3_shell,U3_shell,V3_shell,Gam_m_shell, &
                                             V_seed,Num_nu,Num_R,Num_gam_e,index_syn_intger,n_threads,gam_e,dN_gam_e, &
                                             solver_id)
    implicit none
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,index_syn_intger,n_threads
    integer, intent(in), optional :: solver_id
    integer :: I_tobs,I_gam_e,L1,L,active_solver
    real(8), intent(in) :: e_r,b_r,p_r,f_e_r,z
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),B3(Num_R),M3_shell(Num_R),U3_shell(Num_R)
    real(8), intent(in) :: V3_shell(Num_R),Gam_m_shell(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R)
    real(8), parameter :: secondary_adv_coeff=1.35d-19
    real(8) :: dB,Gam_e_max,Gam_e_m,Gam_e_max_max,Gam_e_min_global,d_x,R_loc,R_Gamma_loc,beta2
    real(8) :: f_r,dDR,dDD,injection_rate,mass_lo,mass_hi,inj_width,adiabatic_rate
    real(8), allocatable :: dEl(:),x(:),dF1(:),temp3(:),dN_x(:),x_edge(:)

    active_solver=electron_resolve_1d_solver_id(solver_id)
    allocate(dEl(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e),temp3(Num_gam_e-1),dN_x(Num_gam_e),x_edge(Num_gam_e+1))

    Gam_e_min_global=one
    Gam_e_max_max=zero
    do I_tobs=2,Num_R
        dB=(B3(I_tobs)+B3(I_tobs-1))/two
        if (dB > zero .and. Gam_m_shell(I_tobs) > one) then
            Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
            Gam_e_max_max=max(Gam_e_max_max,Gam_e_max)
        end if
    end do
    if (Gam_e_max_max <= Gam_e_min_global) error stop "electron_secondary_reverse_evolve: empty secondary electron grid."

    do I_gam_e=1,Num_gam_e
        if (Num_gam_e == 1) then
            gam_e(I_gam_e)=Gam_e_min_global
        else
            gam_e(I_gam_e)=Gam_e_min_global*ten**(dlog10(Gam_e_max_max/Gam_e_min_global)*(I_gam_e-1)/(Num_gam_e-1))
        end if
        dN_gam_e(I_gam_e,1)=zero
    end do

    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
    d_x=dlog10(gam_e(2)/gam_e(1))
    call electron_profile_log_cell_edges(Num_gam_e,gam_e,x_edge)

    do I_tobs=2,Num_R
        if (M3_shell(I_tobs) <= zero .and. M3_shell(I_tobs-1) <= zero) then
            dN_gam_e(:,I_tobs)=dN_gam_e(:,I_tobs-1)
            cycle
        end if
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        beta2=dsqrt(one-one/R_Gamma_loc**2)
        dB=(B3(I_tobs)+B3(I_tobs-1))/two
        if (dB <= zero) error stop "electron_secondary_reverse_evolve: secondary reservoir requires B3 > 0."
        if (Gam_m_shell(I_tobs-1) > one) then
            Gam_e_m=(Gam_m_shell(I_tobs)+Gam_m_shell(I_tobs-1))/two
        else
            Gam_e_m=Gam_m_shell(I_tobs)
        end if
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        f_r=secondary_adv_coeff/beta2/R_Gamma_loc*dB**2/pi
        dDD=R(I_tobs)-R(I_tobs-1)
        dDR=0.7d0/(f_r*Gam_e_max+one/R(I_tobs-1))
        L1=reverse_transport_substeps(dDR,dDD,active_solver)
        dDR=dDD/L1
        dEl=f_r*gam_e
        mass_lo=M3_shell(I_tobs-1)
        mass_hi=M3_shell(I_tobs)
        if (mass_hi < mass_lo) error stop "electron_secondary_reverse_evolve: secondary swept mass must not decrease."
        inj_width=R(I_tobs)-R(I_tobs-1)
        injection_rate=f_e_r*(mass_hi-mass_lo)/(Para_m_p*inj_width)
        if (V3_shell(I_tobs) <= zero .or. V3_shell(I_tobs-1) <= zero) then
            adiabatic_rate=zero
        else
            adiabatic_rate=dlog(V3_shell(I_tobs)/V3_shell(I_tobs-1))/(3d0*(R(I_tobs)-R(I_tobs-1)))
        end if
        call advance_secondary_transport_shell(I_tobs)
    end do

    deallocate(dEl,x,dF1,temp3,dN_x,x_edge)

contains

    subroutine advance_secondary_transport_shell(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs
    real(8) :: dg_adiabatic(L1), dg_source_norm(L1)
    real(8) :: face_speed(Num_gam_e-1)

        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)
        if (active_solver == electron_solver_dg_1d) then
            do L=1,L1
                R_loc=R_loc+dDR
                dg_adiabatic(L)=adiabatic_rate
                if (injection_rate > zero) then
                    dg_source_norm(L)=injection_rate
                else
                    dg_source_norm(L)=zero
                end if
            end do
            call reverse_dg_grid_sequence(Num_gam_e,x_edge,gam_e,L1,dDR,f_r,dg_adiabatic,dg_source_norm, &
                                          p_r,Gam_e_m,Gam_e_max,dN_x,x)
            dN_x=x
            call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
            return
        end if

        do L=1,L1
            R_loc=R_loc+dDR
            if (injection_rate > zero) then
                call electron_build_kinetic_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m,Gam_e_max, &
                                                                          injection_rate,p_r,dF1)
            else
                dF1=zero
            end if
            face_speed=((dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two+adiabatic_rate)/dlog(ten)
            call electron_fullhide_flux_split_step(Num_gam_e,dDR,d_x,face_speed,dF1,dN_x,x,.true.)
            dN_x=x
            if (L == L1) call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
        end do
    end subroutine advance_secondary_transport_shell
end subroutine electron_secondary_reverse_evolve

! Secondary RS 源项历史：逐半径壳层调用同步辐射和SSA核，返回给统一EATS投影使用。
subroutine electron_secondary_reverse_synchrotron(index_syn_intger,Num_nu,Num_R,Num_gam_e,n_threads,R,R_Gamma,B3,gam_e, &
                                                  dN_gam_e,V_seed,z,L_syn_spec,Seed_syn,Nu_a)
    implicit none
    integer, intent(in) :: index_syn_intger,Num_nu,Num_R,Num_gam_e,n_threads
    integer :: I_tobs
    real(8), intent(in) :: R(Num_R),R_Gamma(Num_R),B3(Num_R),gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R)
    real(8), intent(in) :: V_seed(Num_nu),z
    real(8), intent(out) :: L_syn_spec(Num_nu,Num_R),Seed_syn(Num_nu,Num_R),Nu_a(Num_R)
    real(8) :: doppler_den,P_emit_tmp(Num_nu),Tau_syn_tmp(Num_nu)

    L_syn_spec=zero
    Seed_syn=zero
    Nu_a=zero
    do I_tobs=1,Num_R
        if (B3(I_tobs) <= zero) cycle
        call get_syn_selected_state(index_syn_intger,R(I_tobs),B3(I_tobs),Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,I_tobs),V_seed,P_emit_tmp,L_syn_spec(:,I_tobs), &
                                    Seed_syn(:,I_tobs),Tau_syn_tmp)
        doppler_den=R_Gamma(I_tobs)*(one-dsqrt(one-R_Gamma(I_tobs)**(-2)))*(one+z)
        call get_nu_a(R(I_tobs),B3(I_tobs),Num_gam_e,gam_e,dN_gam_e(:,I_tobs),Nu_a(I_tobs))
        Nu_a(I_tobs)=Nu_a(I_tobs)/doppler_den
    end do
end subroutine electron_secondary_reverse_synchrotron

! Secondary RS 分支辐射归并：每个 reservoir 独立输运电子谱，再在源项半径网格上求和。
subroutine electron_secondary_reverse_branch_synchrotron(e_r,b_r,p_r,f_e_r,z,R_Tobs,R_Gamma,R,B3_branch, &
                                                         M3_branch,U3_branch,V3_branch,Gam_m_branch,V_seed, &
                                                         Num_jump,Num_nu,Num_R,Num_gam_e,index_syn_intger,n_threads, &
                                                         Branch_L_syn_spec,L_syn_spec,solver_id)
    implicit none
    integer, intent(in) :: Num_jump,Num_nu,Num_R,Num_gam_e,index_syn_intger,n_threads
    integer, intent(in), optional :: solver_id
    integer :: I_jump
    real(8), intent(in) :: e_r,b_r,p_r,f_e_r,z
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: B3_branch(Num_jump,Num_R),M3_branch(Num_jump,Num_R),U3_branch(Num_jump,Num_R)
    real(8), intent(in) :: V3_branch(Num_jump,Num_R),Gam_m_branch(Num_jump,Num_R)
    real(8), intent(out) :: Branch_L_syn_spec(Num_jump,Num_nu,Num_R),L_syn_spec(Num_nu,Num_R)
    real(8), allocatable :: gam_e_branch(:),dN_branch(:,:),seed_dummy(:,:),nu_a_dummy(:)

    allocate(gam_e_branch(Num_gam_e),dN_branch(Num_gam_e,Num_R),seed_dummy(Num_nu,Num_R),nu_a_dummy(Num_R))
    Branch_L_syn_spec=zero
    L_syn_spec=zero
    do I_jump=1,Num_jump
        if (.not. any(M3_branch(I_jump,:) > zero)) cycle
        call electron_secondary_reverse_evolve(e_r,b_r,p_r,f_e_r,z,R_Tobs,R_Gamma,R,B3_branch(I_jump,:), &
                                               M3_branch(I_jump,:),U3_branch(I_jump,:),V3_branch(I_jump,:), &
                                               Gam_m_branch(I_jump,:),V_seed,Num_nu,Num_R,Num_gam_e, &
                                               index_syn_intger,n_threads,gam_e_branch,dN_branch,solver_id)
        call electron_secondary_reverse_synchrotron(index_syn_intger,Num_nu,Num_R,Num_gam_e,n_threads,R,R_Gamma, &
                                                    B3_branch(I_jump,:),gam_e_branch,dN_branch,V_seed,z, &
                                                    Branch_L_syn_spec(I_jump,:,:),seed_dummy,nu_a_dummy)
        L_syn_spec=L_syn_spec+Branch_L_syn_spec(I_jump,:,:)
    end do
    deallocate(gam_e_branch,dN_branch,seed_dummy,nu_a_dummy)
end subroutine electron_secondary_reverse_branch_synchrotron

subroutine electron_secondary_reverse_branch_reaccelerated(e_r,b_r,p_r,f_e_r,z,R_Tobs,R_Gamma,R,B3_branch, &
                                                           M3_branch,U3_branch,V3_branch,Gam_m_branch,Gamma43_branch, &
                                                           Comp_branch,Parent_branch,V_seed,Num_jump,Num_nu,Num_R,Num_gam_e, &
                                                           index_syn_intger,n_threads,gam_e,dN_total,Branch_L_syn_spec, &
                                                           L_syn_spec,Branch_seed_energy,Branch_reaccel_energy,solver_id)
    implicit none
    integer, intent(in) :: Num_jump,Num_nu,Num_R,Num_gam_e,index_syn_intger,n_threads
    integer, intent(in) :: Parent_branch(Num_jump)
    integer, intent(in), optional :: solver_id
    integer :: I_tobs,I_jump,I_gam_e,L1,L,parent,active_solver
    real(8), intent(in) :: e_r,b_r,p_r,f_e_r,z
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: B3_branch(Num_jump,Num_R),M3_branch(Num_jump,Num_R),U3_branch(Num_jump,Num_R)
    real(8), intent(in) :: V3_branch(Num_jump,Num_R),Gam_m_branch(Num_jump,Num_R),Gamma43_branch(Num_jump,Num_R)
    real(8), intent(in) :: Comp_branch(Num_jump,Num_R)
    real(8), intent(out) :: gam_e(Num_gam_e),dN_total(Num_gam_e,Num_R)
    real(8), intent(out) :: Branch_L_syn_spec(Num_jump,Num_nu,Num_R),L_syn_spec(Num_nu,Num_R)
    real(8), intent(out) :: Branch_seed_energy(Num_jump,Num_R),Branch_reaccel_energy(Num_jump,Num_R)
    real(8), parameter :: secondary_adv_coeff=1.35d-19
    real(8) :: dB,Gam_e_max,Gam_e_m,Gam_e_max_max,Gam_e_min_global,d_x,R_loc,R_Gamma_loc,beta2
    real(8) :: f_r,dDR,dDD,injection_rate,mass_lo,mass_hi,mass_delta,fresh_mass,adiabatic_rate
    real(8) :: parent_mass,transfer_mass,transfer_fraction,boost_factor,seed_energy,out_energy
    logical :: dissipative_shell
    real(8), allocatable :: dEl(:),x(:),dF1(:),temp3(:),x_edge(:),dN_x(:,:),dN_work(:,:)
    real(8), allocatable :: dN_seed(:),dN_boost(:),dN_reaccel(:),dN_gamma_branch(:,:,:),seed_dummy(:,:),nu_a_dummy(:)
    real(8), allocatable :: branch_mass_available(:),fresh_mass_branch(:)

    active_solver=electron_resolve_1d_solver_id(solver_id)
    allocate(dEl(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e),temp3(Num_gam_e-1),x_edge(Num_gam_e+1), &
             dN_x(Num_jump,Num_gam_e),dN_work(Num_jump,Num_gam_e),dN_seed(Num_gam_e),dN_boost(Num_gam_e), &
             dN_reaccel(Num_gam_e),dN_gamma_branch(Num_jump,Num_gam_e,Num_R),seed_dummy(Num_nu,Num_R), &
             nu_a_dummy(Num_R),branch_mass_available(Num_jump),fresh_mass_branch(Num_jump))

    call build_secondary_reaccel_gamma_grid()
    call electron_profile_log_cell_edges(Num_gam_e,gam_e,x_edge)
    d_x=dlog10(gam_e(2)/gam_e(1))
    dN_gamma_branch=zero; dN_total=zero; dN_x=zero; branch_mass_available=zero
    fresh_mass_branch=zero
    Branch_L_syn_spec=zero; L_syn_spec=zero; Branch_seed_energy=zero; Branch_reaccel_energy=zero

    do I_tobs=2,Num_R
        dN_work=dN_gamma_branch(:,:,I_tobs-1)*spread(gam_e*dlog(ten),1,Num_jump)
        call transfer_reaccelerated_parent_electrons(I_tobs)
        do I_jump=1,Num_jump
            call advance_reaccelerated_branch_shell(I_tobs,I_jump)
        end do
    end do
    do I_jump=1,Num_jump
        dN_total=dN_total+dN_gamma_branch(I_jump,:,:)
        if (.not. any(dN_gamma_branch(I_jump,:,:) > zero)) cycle
        call electron_secondary_reverse_synchrotron(index_syn_intger,Num_nu,Num_R,Num_gam_e,n_threads,R,R_Gamma, &
                                                    B3_branch(I_jump,:),gam_e,dN_gamma_branch(I_jump,:,:),V_seed,z, &
                                                    Branch_L_syn_spec(I_jump,:,:),seed_dummy,nu_a_dummy)
        L_syn_spec=L_syn_spec+Branch_L_syn_spec(I_jump,:,:)
    end do
    deallocate(dEl,x,dF1,temp3,x_edge,dN_x,dN_work,dN_seed,dN_boost,dN_reaccel,dN_gamma_branch, &
               seed_dummy,nu_a_dummy,branch_mass_available,fresh_mass_branch)

contains

    subroutine build_secondary_reaccel_gamma_grid()
    implicit none

        Gam_e_min_global=one
        Gam_e_max_max=zero
        do I_tobs=2,Num_R
            do I_jump=1,Num_jump
                dB=(B3_branch(I_jump,I_tobs)+B3_branch(I_jump,I_tobs-1))/two
                if (dB > zero .and. Gam_m_branch(I_jump,I_tobs) > one) then
                    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
                    Gam_e_max_max=max(Gam_e_max_max,Gam_e_max)
                end if
            end do
        end do
        if (Gam_e_max_max <= Gam_e_min_global) &
            error stop "electron_secondary_reverse_branch_reaccelerated: empty electron grid."
        do I_gam_e=1,Num_gam_e
            if (Num_gam_e == 1) then
                gam_e(I_gam_e)=Gam_e_min_global
            else
                gam_e(I_gam_e)=Gam_e_min_global*ten**(dlog10(Gam_e_max_max/Gam_e_min_global)*(I_gam_e-1)/(Num_gam_e-1))
            end if
        end do
    end subroutine build_secondary_reaccel_gamma_grid

    subroutine transfer_reaccelerated_parent_electrons(i_shell)
    implicit none
    integer, intent(in) :: i_shell

        fresh_mass_branch=zero
        do I_jump=1,Num_jump
            mass_delta=M3_branch(I_jump,i_shell)-M3_branch(I_jump,i_shell-1)
            if (mass_delta < zero) error stop "electron_secondary_reverse_branch_reaccelerated: branch mass decreased."
            dissipative_shell=Gam_m_branch(I_jump,i_shell) > one .and. Gamma43_branch(I_jump,i_shell) > one
            fresh_mass=zero
            if (dissipative_shell) fresh_mass=mass_delta
            if (Parent_branch(I_jump) > 0 .and. dissipative_shell .and. mass_delta > zero) then
                parent=Parent_branch(I_jump)
                parent_mass=branch_mass_available(parent)
                if (parent_mass > zero) then
                    transfer_mass=min(mass_delta,parent_mass)
                    transfer_fraction=transfer_mass/parent_mass
                    dN_seed=dN_work(parent,:)*transfer_fraction
                    dN_work(parent,:)=dN_work(parent,:)-dN_seed
                    boost_factor=Gamma43_branch(I_jump,i_shell)
                    if (Comp_branch(I_jump,i_shell) > one) boost_factor=boost_factor*Comp_branch(I_jump,i_shell)**(one/3d0)
                    call boost_log_distribution(boost_factor,dN_seed,dN_boost)
                    call dsa_reaccelerate_distribution(p_r,dN_boost,dN_reaccel)
                    seed_energy=distribution_energy_from_log(dN_seed)
                    out_energy=distribution_energy_from_log(dN_reaccel)
                    Branch_seed_energy(I_jump,i_shell)=Branch_seed_energy(I_jump,i_shell)+seed_energy
                    Branch_reaccel_energy(I_jump,i_shell)=Branch_reaccel_energy(I_jump,i_shell)+out_energy
                    dN_work(I_jump,:)=dN_work(I_jump,:)+dN_reaccel
                    branch_mass_available(parent)=branch_mass_available(parent)-transfer_mass
                    branch_mass_available(I_jump)=branch_mass_available(I_jump)+transfer_mass
                    fresh_mass=mass_delta-transfer_mass
                end if
            end if
            fresh_mass_branch(I_jump)=fresh_mass
            branch_mass_available(I_jump)=branch_mass_available(I_jump)+fresh_mass
            dN_x(I_jump,:)=dN_work(I_jump,:)
        end do
    end subroutine transfer_reaccelerated_parent_electrons

    subroutine advance_reaccelerated_branch_shell(i_shell,jump_index)
    implicit none
    integer, intent(in) :: i_shell,jump_index
    real(8) :: dg_adiabatic(L1), dg_source_norm(L1)
    real(8) :: face_speed(Num_gam_e-1)

        if (M3_branch(jump_index,i_shell) <= zero .and. M3_branch(jump_index,i_shell-1) <= zero) then
            dN_gamma_branch(jump_index,:,i_shell)=dN_x(jump_index,:)/gam_e/dlog(ten)
            return
        end if
        call prepare_branch_shell(i_shell,jump_index)
        call compute_branch_injection(i_shell,jump_index)
        if (active_solver == electron_solver_dg_1d) then
            do L=1,L1
                R_loc=R_loc+dDR
                dg_adiabatic(L)=adiabatic_rate
                if (injection_rate > zero) then
                    dg_source_norm(L)=injection_rate
                else
                    dg_source_norm(L)=zero
                end if
            end do
            call reverse_dg_grid_sequence(Num_gam_e,x_edge,gam_e,L1,dDR,f_r,dg_adiabatic,dg_source_norm, &
                                          p_r,Gam_e_m,Gam_e_max,dN_x(jump_index,:),x)
            dN_x(jump_index,:)=x
            dN_gamma_branch(jump_index,:,i_shell)=dN_x(jump_index,:)/gam_e/dlog(ten)
            return
        end if

        do L=1,L1
            R_loc=R_loc+dDR
            if (injection_rate > zero) then
                call electron_build_kinetic_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m,Gam_e_max, &
                                                                          injection_rate,p_r,dF1)
            else
                dF1=zero
            end if
            face_speed=((dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two+adiabatic_rate)/dlog(ten)
            call electron_fullhide_flux_split_step(Num_gam_e,dDR,d_x,face_speed,dF1,dN_x(jump_index,:),x,.true.)
            dN_x(jump_index,:)=x
        end do
        dN_gamma_branch(jump_index,:,i_shell)=dN_x(jump_index,:)/gam_e/dlog(ten)
    end subroutine advance_reaccelerated_branch_shell

    subroutine prepare_branch_shell(i_shell,jump_index)
    implicit none
    integer, intent(in) :: i_shell,jump_index

        R_loc=R(i_shell-1)
        R_Gamma_loc=(R_Gamma(i_shell)+R_Gamma(i_shell-1))/two
        beta2=dsqrt(one-one/R_Gamma_loc**2)
        dB=(B3_branch(jump_index,i_shell)+B3_branch(jump_index,i_shell-1))/two
        if (dB <= zero) error stop "electron_secondary_reverse_branch_reaccelerated: active branch requires B3 > 0."
        if (Gam_m_branch(jump_index,i_shell-1) > one) then
            Gam_e_m=(Gam_m_branch(jump_index,i_shell)+Gam_m_branch(jump_index,i_shell-1))/two
        else
            Gam_e_m=Gam_m_branch(jump_index,i_shell)
        end if
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        f_r=secondary_adv_coeff/beta2/R_Gamma_loc*dB**2/pi
        dDD=R(i_shell)-R(i_shell-1)
        dDR=0.7d0/(f_r*Gam_e_max+one/R(i_shell-1))
        L1=reverse_transport_substeps(dDR,dDD,active_solver)
        dDR=dDD/L1
        dEl=f_r*gam_e
        if (V3_branch(jump_index,i_shell) <= zero .or. V3_branch(jump_index,i_shell-1) <= zero) then
            adiabatic_rate=zero
        else
            adiabatic_rate=dlog(V3_branch(jump_index,i_shell)/V3_branch(jump_index,i_shell-1))/(3d0*dDD)
        end if
    end subroutine prepare_branch_shell

    subroutine compute_branch_injection(i_shell,jump_index)
    implicit none
    integer, intent(in) :: i_shell,jump_index

        mass_lo=M3_branch(jump_index,i_shell-1)
        mass_hi=M3_branch(jump_index,i_shell)
        if (mass_hi < mass_lo) error stop "electron_secondary_reverse_branch_reaccelerated: swept mass must not decrease."
        injection_rate=f_e_r*fresh_mass_branch(jump_index)/(Para_m_p*dDD)
    end subroutine compute_branch_injection

    subroutine boost_log_distribution(boost,dN_in,dN_out)
    implicit none
    real(8), intent(in) :: boost,dN_in(Num_gam_e)
    real(8), intent(out) :: dN_out(Num_gam_e)
    integer :: i_src,i_hi
    real(8) :: x_src,pos,frac

        if (boost <= one) error stop "electron_secondary_reverse_branch_reaccelerated: boost must exceed unity."
        dN_out=zero
        do I_gam_e=1,Num_gam_e
            x_src=dlog10(gam_e(I_gam_e)/boost)
            if (x_src < dlog10(gam_e(1)) .or. x_src > dlog10(gam_e(Num_gam_e))) cycle
            pos=(x_src-dlog10(gam_e(1)))/d_x+one
            i_src=int(pos)
            if (i_src < 1) cycle
            if (i_src >= Num_gam_e) then
                dN_out(I_gam_e)=dN_in(Num_gam_e)
            else
                i_hi=i_src+1
                frac=pos-dble(i_src)
                dN_out(I_gam_e)=(one-frac)*dN_in(i_src)+frac*dN_in(i_hi)
            end if
        end do
    end subroutine boost_log_distribution

    subroutine dsa_reaccelerate_distribution(p,dN_seed_log,dN_out_log)
    implicit none
    real(8), intent(in) :: p,dN_seed_log(Num_gam_e)
    real(8), intent(out) :: dN_out_log(Num_gam_e)
    integer :: i
    real(8) :: integral,dN_dgamma

        if (p <= two) error stop "electron_secondary_reverse_branch_reaccelerated: DSA requires p > 2."
        integral=zero
        dN_out_log=zero
        do i=1,Num_gam_e
            dN_dgamma=dN_seed_log(i)/(gam_e(i)*dlog(ten))
            integral=integral+dN_dgamma*gam_e(i)**(p-one)*(gam_e(i)*dlog(ten)*d_x)
            dN_out_log(i)=(p-one)*gam_e(i)**(-p)*integral*gam_e(i)*dlog(ten)
        end do
    end subroutine dsa_reaccelerate_distribution

    real(8) function distribution_energy_from_log(dN_log)
    implicit none
    real(8), intent(in) :: dN_log(Num_gam_e)

        distribution_energy_from_log=sum((gam_e-one)*Para_m_e*Para_c**2*dN_log)*d_x
    end function distribution_energy_from_log
end subroutine electron_secondary_reverse_branch_reaccelerated

integer function reverse_transport_substeps(candidate_dr,shell_dr,solver) result(nstep)
    implicit none
    integer, intent(in) :: solver
    real(8), intent(in) :: candidate_dr,shell_dr

    if (solver == electron_solver_dg_1d) then
        nstep=reverse_dg_base_substeps
    else
        nstep=max(100,min(1000,int(shell_dr/candidate_dr)))
    end if
end function reverse_transport_substeps

real(8) function reverse_dg_kinetic_break(gamma_m,gamma_max) result(gamma_break)
    implicit none
    real(8), intent(in) :: gamma_m,gamma_max

    gamma_break=gamma_m
    if (gamma_max > gamma_m) gamma_break=min(gamma_max,20d0*max(gamma_m,one))
end function reverse_dg_kinetic_break

subroutine reverse_dg_grid_sequence(Num_gam_e,x_edge,gam_e,num_step,dR,rad_coeff,adiabatic_rate_step,source_norm_step, &
                                    p,gamma_m,gamma_max,dN_x_in,dN_x_out)
    implicit none
    integer, intent(in) :: Num_gam_e,num_step
    integer :: i,step
    real(8), intent(in) :: x_edge(Num_gam_e+1),gam_e(Num_gam_e),dR,rad_coeff
    real(8), intent(in) :: adiabatic_rate_step(num_step),source_norm_step(num_step),p,gamma_m,gamma_max,dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    type(electron_dg1d_mesh) :: mesh
    real(8), allocatable :: state(:),advanced(:),source_nodes(:),source_template(:),dN_coord(:),dN_dgamma(:),coord_edge(:)
    real(8) :: input_content,dg_content,projected_content
    logical :: has_source

    call electron_dg1d_build_four_velocity_mesh(x_edge(1),x_edge(Num_gam_e+1),dlog10(gamma_m), &
                                                dlog10(reverse_dg_kinetic_break(gamma_m,gamma_max)), &
                                                dlog10(gamma_max),electron_four_velocity_grid_gamma_scale,mesh)
    allocate(state(mesh%ntot),advanced(mesh%ntot),source_nodes(mesh%ntot),source_template(mesh%ntot),dN_coord(Num_gam_e), &
             dN_dgamma(Num_gam_e),coord_edge(Num_gam_e+1))
    do i=1,Num_gam_e+1
        coord_edge(i)=electron_coord_from_xgamma(mesh%coord_kind,mesh%coord_scale,x_edge(i))
    end do
    do i=1,mesh%ntot
        state(i)=reverse_interp_log_grid(Num_gam_e,x_edge,dN_x_in,mesh%x_gamma(i))*mesh%dxgamma_dcoord(i)
    end do
    input_content=sum(dN_x_in*(x_edge(2:Num_gam_e+1)-x_edge(1:Num_gam_e)))
    call electron_dg1d_scale_to_content(mesh,input_content,state)
    has_source = maxval(source_norm_step) > zero
    if (has_source) then
        call electron_dg1d_project_kinetic_source(mesh,one,p,gamma_m,gamma_max,source_template)
        call electron_dg1d_scale_to_content(mesh,one,source_template)
    endif
    do step=1,num_step
        if (source_norm_step(step) > zero) then
            source_nodes=source_norm_step(step)*source_template
        else
            source_nodes=zero
        end if
        call electron_dg1d_advance_characteristic_step(mesh,dR,rad_coeff,adiabatic_rate_step(step), &
                                                       source_nodes,state,advanced)
        call electron_dg1d_apply_positive_kernel_filter(mesh,advanced)
        call electron_dg1d_limit_positive_cell_preserving(mesh,advanced)
        state=advanced
    end do
    call electron_dg1d_integral(mesh,state,dg_content)
    call electron_dg1d_project_to_coord_cells(mesh,state,Num_gam_e,coord_edge,dN_coord)
    call electron_shell_dcoord_to_dndgamma_exp_centers(Num_gam_e,coord_edge,mesh%coord_scale,gam_e,dN_coord,dN_dgamma)
    dN_x_out=dN_dgamma*gam_e*dlog(ten)
    projected_content=sum(dN_x_out*(x_edge(2:Num_gam_e+1)-x_edge(1:Num_gam_e)))
    if (dg_content > zero) then
        if (projected_content <= zero) error stop "reverse_dg_grid_sequence: projection lost positive content."
        dN_x_out=dN_x_out*(dg_content/projected_content)
    end if
    deallocate(state,advanced,source_nodes,source_template,dN_coord,dN_dgamma,coord_edge)
end subroutine reverse_dg_grid_sequence

real(8) function reverse_interp_log_grid(Num_gam_e,x_edge,values,x_eval) result(value)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: i_cell,lo,hi,mid
    real(8), intent(in) :: x_edge(Num_gam_e+1),values(Num_gam_e),x_eval
    real(8) :: x_left,x_right

    x_left=0.5d0*(x_edge(1)+x_edge(2))
    if (x_eval <= x_left) then
        value=values(1)
        return
    end if
    x_right=0.5d0*(x_edge(Num_gam_e)+x_edge(Num_gam_e+1))
    if (x_eval >= x_right) then
        value=values(Num_gam_e)
        return
    end if
    lo=1
    hi=Num_gam_e-1
    do while (lo < hi)
        mid=(lo+hi)/2
        x_right=0.5d0*(x_edge(mid+1)+x_edge(mid+2))
        if (x_eval <= x_right) then
            hi=mid
        else
            lo=mid+1
        end if
    end do
    i_cell=lo
    x_left=0.5d0*(x_edge(i_cell)+x_edge(i_cell+1))
    x_right=0.5d0*(x_edge(i_cell+1)+x_edge(i_cell+2))
    value=values(i_cell)+(values(i_cell+1)-values(i_cell))*(x_eval-x_left)/(x_right-x_left)
end function reverse_interp_log_grid

end module electron_reverse_kernel
