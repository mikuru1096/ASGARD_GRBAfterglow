!Calculate the electron distibutions of forward shock.
!New improvement at 11.29.2022
!****************************************************************************************
!******************************* main program *******************************************
!****************************************************************************************
! 电子1D全隐格式主驱动。
! 顺序: unpack boundary/config -> build log-four-velocity grid -> initialize electron distribution
!       -> loop shells: density/B/gamma_m/c -> injection/cooling -> transport -> synch/SSA diagnostics
!       -> return shell-level electron, synchrotron seed, and break frequencies.
subroutine fs_electron_fullhide_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                                   adaptive_substeps,substep_rtol,substep_min,substep_max,thermal_electrons, &
                                   gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    use electron_common
    use electron_injection_profiles, only: electron_build_source_term_exp_cutoff_edges, electron_add_thermal_source_term, &
                                           electron_initial_powerlaw_exp_cutoff, electron_initial_powerlaw_exp_cutoff_coord_edges, &
                                           electron_build_source_term_exp_cutoff_coord_edges, electron_add_thermal_population, &
                                           electron_profile_log_cell_edges
    use electron_energy_coordinate_common, only: electron_build_four_velocity_grid, electron_dxgamma_dcoord, &
                                                 electron_coord_log_four_velocity_sq, &
                                                 electron_four_velocity_grid_gamma_scale
    use electron_radiation_kernel, only: get_syn_selected_state, get_nu_a_from_tau_grid
    use electron_cooling_kernel, only: get_forward_cooling
    use electron_shell_transport_common, only: electron_shell_fullhide_step, electron_shell_fullhide_spacetime_sequence, &
                                               electron_shell_flux_split_coord_step, &
                                               electron_shell_flux_split_coord_sequence, &
                                               electron_shell_dcoord_to_dndgamma_exp_centers
    use electron_transport_common, only: electron_dnx_to_dndgamma_exp_centers
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in) :: adaptive_substeps,substep_min,substep_max,thermal_electrons
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: substep_rtol
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R), &
                            Seed_syn(Num_nu,Num_R), V_m(Num_R), V_c(Num_R), V_a(Num_R)
    
    real(8),allocatable,dimension (:) :: dEl,dEl_step,dEL_mean,x,dN_x,x_edge,coord_edge,dxdy_grid, &
                                         dN_full,dN_half,dN_half2,dF1,dEL_mean_step,P_emit_shell,Tau_syn_shell
    real(8),allocatable,dimension (:,:) :: dF_steps,face_coupling
    logical :: is_uniform_density,budget_diag_enabled
    integer :: env_len,env_status,I_face
    character(len=32) :: diag_env
    real(8) :: dDR_xi
    real(8) :: n_before_step,n_after_step,inj_step,n_budget,rel_loss_xi_max
    real(8) :: source_integral,adiabatic_integral,l_count_real,step_sum,step_sq_sum
    real(8) :: radius_sum,radius_sq_sum,source_prefactor,coord_scale,dg_gamma_scale,face_coord,face_jac
    allocate (dEl(Num_gam_e),dEl_step(Num_gam_e),dEL_mean(Num_gam_e-1),x(Num_gam_e),dN_x(Num_gam_e), &
              dN_full(Num_gam_e), &
              x_edge(Num_gam_e+1),coord_edge(Num_gam_e+1),dxdy_grid(Num_gam_e), &
              dN_half(Num_gam_e),dN_half2(Num_gam_e),dF1(Num_gam_e),dEL_mean_step(Num_gam_e-1), &
              P_emit_shell(Num_nu),Tau_syn_shell(Num_nu))
    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    if (thermal_electrons /= 0) then
        if (f_e <= zero .or. f_e > one) error stop 'thermal electrons require 0 < f_e <= 1'
    end if

    P_syn=zero
    Seed_syn=zero
    V_m=zero
    V_c=zero
    V_a=zero

    call initialize_forward_electron_state()
    d_x=dlog10(gam_e(2)/gam_e(1))
    is_uniform_density=(A_star <= zero .and. f_jump == one)
    budget_diag_enabled=.false.
    diag_env=''
    call get_environment_variable('ASGARD_DIAG_1D_BUDGET',diag_env,length=env_len,status=env_status)
    if (env_status == 0 .and. env_len > 0) then
        if (diag_env(1:1) /= '0') budget_diag_enabled=.true.
    end if
    rel_loss_xi_max=zero
!    factor_adv=Para_sigmaT/(6.0d0*pi*Para_m_energy)

    do I_tobs=2,Num_R
        call prepare_fullhide_shell(I_tobs)
        call prepare_shell_cooling_faces()
        if (adaptive_substeps == 0) then
            if (dDD <= zero) error stop 'fs_electron_fullhide_1d requires increasing radius grid'
            if (dDR <= zero) error stop 'fs_electron_fullhide_1d requires positive cooling substep width'
            L1=max(100,min(1000,int(dDD/dDR)))
            dDR=dDD/dble(L1)
            if (is_uniform_density .and. thermal_electrons == 0) then
                allocate(dF_steps(Num_gam_e,1),face_coupling(Num_gam_e-1,1))
                l_count_real=dble(L1)
                step_sum=l_count_real*(l_count_real+one)/two
                step_sq_sum=l_count_real*(l_count_real+one)*(two*l_count_real+one)/6d0
                radius_sum=l_count_real*R_loc+dDR*step_sum
                radius_sq_sum=l_count_real*R_loc*R_loc+two*R_loc*dDR*step_sum+dDR*dDR*step_sq_sum
                source_prefactor=4d0/3d0*pi*dNe*f_e*Gam_e_m_p
                source_integral=dDR*source_prefactor*(3d0*radius_sq_sum+dDR*(3d0*radius_sum+l_count_real*dDR))
                adiabatic_integral=zero
                do L=1,L1
                    adiabatic_integral=adiabatic_integral+dDR/(R_loc+dDR*dble(L))
                end do
                call electron_build_source_term_exp_cutoff_coord_edges(Num_gam_e,coord_edge,coord_scale, &
                                                                       Gam_e_m,Gam_e_max,source_integral,p, &
                                                                       dF_steps(:,1))
                do I_face=1,Num_gam_e-1
                    face_coord=coord_edge(I_face+1)
                    face_jac=dlog(ten)*electron_dxgamma_dcoord(electron_coord_log_four_velocity_sq,coord_scale,face_coord)
                    face_coupling(I_face,1)=(dDD*(dEl(I_face)+dEl(I_face+1))/two+adiabatic_integral)/face_jac
                end do
                if (budget_diag_enabled) then
                    n_before_step=sum(dN_x*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
                    inj_step=sum(dF_steps(:,1)*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
                end if
                call electron_shell_flux_split_coord_sequence(Num_gam_e,coord_edge,face_coupling(:,1),dF_steps(:,1),dN_x,x)
                if (budget_diag_enabled) then
                    n_after_step=sum(x*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
                    n_budget=n_before_step+inj_step
                    if (n_budget > zero) rel_loss_xi_max=max(rel_loss_xi_max,max(zero,(n_budget-n_after_step)/n_budget))
                end if
                dN_x=x
                deallocate(dF_steps,face_coupling)
            else
                do L=1,L1
                    R_loc=R_loc+dDR
                    call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
                    DB_step=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
                    Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
                    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
                    call electron_gamma_m_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
                    Gam_e_m_p_step=(one-p)/(Gam_e_max_step**(one-p)-Gam_e_m_step**(one-p))
                    Q=4d0/3d0*pi*(3d0*R_loc**2+dDR*(3d0*R_loc+dDR))*dNe*f_e*Gam_e_m_p_step
                    call electron_build_source_term_exp_cutoff_coord_edges(Num_gam_e,coord_edge,coord_scale, &
                                                                           Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
                    if (thermal_electrons /= 0) then
                        call electron_add_thermal_source_term(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam, &
                                                              Q*(one-f_e)/(f_e*Gam_e_m_p_step),dF1)
                    end if
                    if (dNe_shell > zero) then
                        dEl_step=dEl*(dNe/dNe_shell)
                    else
                        dEl_step=dEl
                    end if
                    if (budget_diag_enabled) then
                        n_before_step=sum(dN_x*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
                        inj_step=dDR*sum(dF1*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
                    end if
                    call electron_shell_flux_split_coord_step(Num_gam_e,dDR,coord_edge,coord_scale, &
                                                              dEl_step,one/R_loc,dF1,dN_x,x)
                    if (budget_diag_enabled) then
                        n_after_step=sum(x*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
                        n_budget=n_before_step+inj_step
                        if (n_budget > zero) rel_loss_xi_max=max(rel_loss_xi_max,max(zero,(n_budget-n_after_step)/n_budget))
                        if (I_tobs <= 6 .and. L == L1) then
                            print '(A,1X,I4,1X,ES12.4,1X,ES12.4,1X,ES12.4)', &
                                  'BUDGET1D shell', I_tobs, n_before_step, inj_step, n_after_step
                        end if
                    end if
                    dN_x=x
                end do
            end if
            call store_forward_shell_distribution(I_tobs)
        else
            dR_min=dDD/max(substep_max,1)
            dR_max=dDD/max(substep_min,1)
            dR_try=min(min(dDR_xi,dR_max),dDD)
            dR_left=dDD

            do while (dR_left > zero)
                dR_try=min(dR_try,dR_left)
                if (dR_left < 1.5d0*dR_try) dR_try=dR_left

                R_full=R_loc+dR_try
                if (is_uniform_density) then
                    dNe_full=dNe
                else
                    call dynamics_external_density_profile(A_star,dNe_ISM,R_full,R0,1,R_tr,f_jump,f_wide,dNe_full)
                end if
                DB_full=0.39d0*dsqrt(Epsilon_b*dNe_full*(R_Gamma_loc*(R_Gamma_loc-one)))
                Gam_e_max_full=3d0*Para_m_energy/dsqrt(8d0*DB_full*Para_e**3)
                temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
                call electron_gamma_m_exact(p,temp_gam,Gam_e_max_full,Gam_e_m_full)
                Gam_e_m_p_full=(one-p)/(Gam_e_max_full**(one-p)-Gam_e_m_full**(one-p))
                call electron_injection_prefactor(R_full,dR_try,dNe_full,f_e,Gam_e_m_p_full,Q)
                call electron_build_source_term_exp_cutoff_coord_edges(Num_gam_e,coord_edge,coord_scale, &
                                                                       Gam_e_m_full,Gam_e_max_full,Q,p,dF1)
                if (thermal_electrons /= 0) then
                    thermal_count=Q*(one-f_e)/(f_e*Gam_e_m_p_full)
                    call electron_add_thermal_source_term(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam,thermal_count,dF1)
                end if
                if (dNe_shell > zero) then
                    dEl_step=dEl*(dNe_full/dNe_shell)
                else
                    dEl_step=dEl
                end if
                call electron_shell_flux_split_coord_step(Num_gam_e,dR_try,coord_edge,coord_scale, &
                                                          dEl_step,one/R_full,dF1,dN_x,dN_full)

                dR_half=0.5d0*dR_try
                R_half=R_loc+dR_half
                if (is_uniform_density) then
                    dNe_half=dNe
                else
                    call dynamics_external_density_profile(A_star,dNe_ISM,R_half,R0,1,R_tr,f_jump,f_wide,dNe_half)
                end if
                DB_half=0.39d0*dsqrt(Epsilon_b*dNe_half*(R_Gamma_loc*(R_Gamma_loc-one)))
                Gam_e_max_half=3d0*Para_m_energy/dsqrt(8d0*DB_half*Para_e**3)
                temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
                call electron_gamma_m_exact(p,temp_gam,Gam_e_max_half,Gam_e_m_half)
                Gam_e_m_p_half=(one-p)/(Gam_e_max_half**(one-p)-Gam_e_m_half**(one-p))
                call electron_injection_prefactor(R_half,dR_half,dNe_half,f_e,Gam_e_m_p_half,Q)
                call electron_build_source_term_exp_cutoff_coord_edges(Num_gam_e,coord_edge,coord_scale, &
                                                                       Gam_e_m_half,Gam_e_max_half,Q,p,dF1)
                if (thermal_electrons /= 0) then
                    thermal_count=Q*(one-f_e)/(f_e*Gam_e_m_p_half)
                    call electron_add_thermal_source_term(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam,thermal_count,dF1)
                end if
                if (dNe_shell > zero) then
                    dEl_step=dEl*(dNe_half/dNe_shell)
                else
                    dEl_step=dEl
                end if
                call electron_shell_flux_split_coord_step(Num_gam_e,dR_half,coord_edge,coord_scale, &
                                                          dEl_step,one/R_half,dF1,dN_x,dN_half)

                call electron_injection_prefactor(R_full,dR_half,dNe_full,f_e,Gam_e_m_p_full,Q)
                call electron_build_source_term_exp_cutoff_coord_edges(Num_gam_e,coord_edge,coord_scale, &
                                                                       Gam_e_m_full,Gam_e_max_full,Q,p,dF1)
                if (thermal_electrons /= 0) then
                    thermal_count=Q*(one-f_e)/(f_e*Gam_e_m_p_full)
                    call electron_add_thermal_source_term(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam,thermal_count,dF1)
                end if
                if (dNe_shell > zero) then
                    dEl_step=dEl*(dNe_full/dNe_shell)
                else
                    dEl_step=dEl
                end if
                call electron_shell_flux_split_coord_step(Num_gam_e,dR_half,coord_edge,coord_scale, &
                                                          dEl_step,one/R_full,dF1,dN_half,dN_half2)

                call electron_max_relative_error(Num_gam_e,dN_half2,dN_full,step_error)
                if (step_error <= substep_rtol .or. dR_try <= dR_min) then
                    dN_x=dN_half2
                    R_loc=R_full
                    dNe=dNe_full
                    dR_left=dR_left-dR_try
                    if (dR_left > zero) then
                        if (step_error <= 0.25d0*substep_rtol) then
                            dR_try=min(two*dR_try,dR_max)
                        end if
                    end if
                else
                    dR_try=max(0.5d0*dR_try,dR_min)
                end if
            end do
            call store_forward_shell_distribution(I_tobs)
        end if
    end do

    deallocate (dEl,dEl_step,dEL_mean,x,dN_x,x_edge,coord_edge,dxdy_grid,dN_full,dN_half,dN_half2,dF1, &
                dEL_mean_step,P_emit_shell,Tau_syn_shell)
    if (budget_diag_enabled) then
        print '(A,1X,ES12.4)', 'BUDGET1D max_rel_loss', rel_loss_xi_max
    end if

    return

contains

    subroutine initialize_forward_electron_state()
    implicit real(8)(A-H,O-Z)

        call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
        Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
        if (R_Gamma(1) < one) error stop 'fs_electron_fullhide_1d requires initial Gamma >= 1'
        beta_Gam=dsqrt(one-one/R_Gamma(1)**2)

        call initialize_forward_four_velocity_grid(Gam_e_max_max)
        if (thermal_electrons == 0) then
            call electron_initial_powerlaw_exp_cutoff_coord_edges(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                                                  Num_gam_e,coord_edge,coord_scale,dN_x)
            call electron_shell_dcoord_to_dndgamma_exp_centers(Num_gam_e,coord_edge,coord_scale,gam_e,dN_x, &
                                                               dN_gam_e(:,1))
        else
            call electron_initial_powerlaw_exp_cutoff(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                                      Num_gam_e,gam_e,dN_gam_e(:,1))
            call electron_add_thermal_population(Num_gam_e,gam_e,R_Gamma(1)*beta_Gam,Para_N_e_ini*(one-f_e), &
                                                 dN_gam_e(:,1))
            dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)*dxdy_grid
        end if
    end subroutine initialize_forward_electron_state

    subroutine initialize_forward_four_velocity_grid(grid_gamma_max)
    implicit real(8)(A-H,O-Z)
    integer :: I_grid
    real(8), intent(in) :: grid_gamma_max

        dg_gamma_scale=electron_four_velocity_grid_gamma_scale
        coord_scale=dg_gamma_scale*dg_gamma_scale-one
        call electron_build_four_velocity_grid(Num_gam_e,one,electron_exp_tail_grid_factor*grid_gamma_max, &
                                               dg_gamma_scale,gam_e,coord_edge,x_edge)
        do I_grid=1,Num_gam_e
            dxdy_grid(I_grid)=electron_dxgamma_dcoord(electron_coord_log_four_velocity_sq,coord_scale, &
                                                      0.5d0*(coord_edge(I_grid)+coord_edge(I_grid+1)))
        end do
    end subroutine initialize_forward_four_velocity_grid

    subroutine prepare_shell_cooling_faces()
    implicit real(8)(A-H,O-Z)

        dEL_mean=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two/dlog(ten)
        dDR_xi=dDR
    end subroutine prepare_shell_cooling_faces

    subroutine store_forward_shell_distribution(I_tobs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs

        call electron_shell_dcoord_to_dndgamma_exp_centers(Num_gam_e,coord_edge,coord_scale,gam_e,dN_x, &
                                                           dN_gam_e(:,I_tobs))
    end subroutine store_forward_shell_distribution

    subroutine prepare_fullhide_shell(I_tobs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs

        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        if (R_Gamma_loc < one) error stop 'fs_electron_fullhide_1d requires Gamma >= 1'
        beta_Gam=dsqrt(one-one/R_Gamma_loc**2)
        call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        dNe_shell=dNe

        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1d0/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)*dxdy_grid
        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))

        call get_syn_selected_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_emit_shell, &
                                    P_syn(:,I_tobs),Seed_syn(:,I_tobs),Tau_syn_shell)
        call get_nu_a_from_tau_grid(Num_nu,V_seed,Tau_syn_shell,temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))

        call get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                 R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                 P_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl)
    end subroutine prepare_fullhide_shell
end subroutine fs_electron_fullhide_1d

! 1D joint-coupling electron pass: use externally closed photon and secondary pair source fields.
! Joint electron-photon pass: same fixed-substep fullhide shell transport, but cooling seed and
! secondary source are supplied by the shell-level photon/hadronic feedback stage.
subroutine fs_electron_fullhide_1d_coupled(Boundary,R_Tobs,R_Gamma,R,V_seed,Seed_cooling,Secondary_source,n,Num_nu,Num_R, &
                                           Num_gam_e,index_Y,index_syn_intger,n_threads,adaptive_substeps, &
                                           substep_rtol,substep_min,substep_max,thermal_electrons, &
                                           gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    use electron_common
    use electron_injection_profiles, only: electron_build_source_term_exp_cutoff_edges, electron_add_thermal_source_term, &
                                           electron_profile_log_cell_edges
    use electron_radiation_kernel, only: get_nu_a, get_syn_selected
    use electron_cooling_kernel, only: electron_cooling_ic_loss_emissivity_budget, assemble_forward_cooling_split
    use electron_shell_transport_common, only: electron_shell_fullhide_step
    use electron_transport_common, only: electron_dnx_to_dndgamma_exp_centers
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in) :: adaptive_substeps,substep_min,substep_max,thermal_electrons
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: Seed_cooling(Num_nu,Num_R),Secondary_source(Num_gam_e,Num_R),substep_rtol
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R), &
                            Seed_syn(Num_nu,Num_R),V_m(Num_R),V_c(Num_R),V_a(Num_R)
    real(8), allocatable :: dEl(:),dEL_mean(:),dEL_mean_step(:),cooling_aux(:),x(:),dN_x(:),x_edge(:),dF1(:), &
                            gam_e_rad(:),dN_gam_e_rad(:)
    logical :: budget_diag_enabled
    integer :: Num_gam_rad,env_len,env_status
    character(len=32) :: diag_env
    real(8) :: n_before_step,n_after_step,inj_step,n_budget,rel_loss_xi_max,thermal_count

    if (adaptive_substeps /= 0) error stop 'fs_electron_fullhide_1d_coupled requires fixed shell substeps'
    if (index_Y /= 1) error stop 'fs_electron_fullhide_1d_coupled requires index_Y=1'
    if (substep_min < 1 .or. substep_max < 1) error stop 'fs_electron_fullhide_1d_coupled requires positive substep bounds'
    if (substep_rtol <= zero) error stop 'fs_electron_fullhide_1d_coupled requires positive substep_rtol'
    allocate(dEl(Num_gam_e),dEL_mean(Num_gam_e-1),dEL_mean_step(Num_gam_e-1),cooling_aux(Num_gam_e), &
             x(Num_gam_e),dN_x(Num_gam_e),x_edge(Num_gam_e+1),dF1(Num_gam_e), &
             gam_e_rad(Num_gam_e),dN_gam_e_rad(Num_gam_e))

    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    if (thermal_electrons /= 0) then
        if (f_e <= zero .or. f_e > one) error stop 'thermal electrons require 0 < f_e <= 1'
    end if

    P_syn=zero
    Seed_syn=zero
    V_m=zero
    V_c=zero
    V_a=zero

    call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)
    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-one)
    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    if (R_Gamma(1) < one) error stop 'fs_electron_fullhide_1d_coupled requires initial Gamma >= 1'
    beta_Gam=dsqrt(one-one/R_Gamma(1)**2)
    if (thermal_electrons == 0) then
        call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                          electron_initial_grid_log_edges,gam_e,dN_x,x_edge)
        call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,1))
    else
        call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                          electron_initial_grid_gamma,gam_e,dN_gam_e(:,1),thermal_electrons=thermal_electrons, &
                                          f_e=f_e,four_v=R_Gamma(1)*beta_Gam)
        dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
        call electron_profile_log_cell_edges(Num_gam_e,gam_e,x_edge)
    end if

    d_x=dlog10(gam_e(2)/gam_e(1))
    budget_diag_enabled=.false.
    diag_env=''
    call get_environment_variable('ASGARD_DIAG_1D_BUDGET',diag_env,length=env_len,status=env_status)
    if (env_status == 0 .and. env_len > 0) then
        if (diag_env(1:1) /= '0') budget_diag_enabled=.true.
    end if
    rel_loss_xi_max=zero

    do I_tobs=2,Num_R
        call prepare_coupled_shell(I_tobs)
        dEL_mean=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two/dlog(ten)
        if (dDD <= zero) error stop 'fs_electron_fullhide_1d_coupled requires increasing radius grid'
        if (dDR <= zero) error stop 'fs_electron_fullhide_1d_coupled requires positive cooling substep width'
        L1=max(100,min(1000,int(dDD/dDR)))
        dDR=dDD/dble(L1)

        do L=1,L1
            R_loc=R_loc+dDR
            call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
            DB_step=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
            Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
            temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
            call electron_gamma_m_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
            Gam_e_m_p_step=(one-p)/(Gam_e_max_step**(one-p)-Gam_e_m_step**(one-p))
            Q=4d0/3d0*pi*(3d0*R_loc**2+dDR*(3d0*R_loc+dDR))*dNe*f_e*Gam_e_m_p_step
            call electron_build_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
            dF1=dF1+Secondary_source(:,I_tobs)
            if (thermal_electrons /= 0) then
                thermal_count=Q*(one-f_e)/(f_e*Gam_e_m_p_step)
                call electron_add_thermal_source_term(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam,thermal_count,dF1)
            end if
            if (dNe_shell > zero) then
                dEL_mean_step=dEL_mean*(dNe/dNe_shell)
            else
                dEL_mean_step=dEL_mean
            end if
            if (budget_diag_enabled) then
                n_before_step=sum(dN_x)*d_x
                inj_step=dDR*sum(dF1)*d_x
            end if
            call electron_shell_fullhide_step(Num_gam_e,R_loc,dDR,d_x,dEL_mean_step,dF1,dN_x,x)
            if (budget_diag_enabled) then
                n_after_step=sum(x)*d_x
                n_budget=n_before_step+inj_step
                if (n_budget > zero) rel_loss_xi_max=max(rel_loss_xi_max,max(zero,(n_budget-n_after_step)/n_budget))
            end if
            dN_x=x
        end do
        call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
    end do

    deallocate(dEl,dEL_mean,dEL_mean_step,cooling_aux,x,dN_x,x_edge,dF1,gam_e_rad,dN_gam_e_rad)
    if (budget_diag_enabled) print '(A,1X,ES12.4)', 'BUDGET1D coupled max_rel_loss', rel_loss_xi_max
    return

contains

    subroutine prepare_coupled_shell(I_tobs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs

        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        if (R_Gamma_loc < one) error stop 'fs_electron_fullhide_1d_coupled requires Gamma >= 1'
        beta_Gam=dsqrt(one-one/R_Gamma_loc**2)
        call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        dNe_shell=dNe
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1d0/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)
        call electron_prepare_radiation_spectrum(Num_gam_e,gam_e,dN_gam_e(:,I_tobs-1), &
                                                 Num_gam_rad,gam_e_rad,dN_gam_e_rad)
        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        call get_nu_a(R_loc,DB,Num_gam_rad,gam_e_rad(1:Num_gam_rad),dN_gam_e_rad(1:Num_gam_rad),temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        call get_syn_selected(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                              gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_syn(:,I_tobs),Seed_syn(:,I_tobs))
        Gam_e_max_cool=Gam_e_max
        call electron_cooling_ic_loss_emissivity_budget(Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                                        Seed_cooling(:,I_tobs),cooling_aux)
        call assemble_forward_cooling_split(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max_cool,R_loc, &
                                            R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                            Seed_syn(:,I_tobs),cooling_aux,dEl)
    end subroutine prepare_coupled_shell
end subroutine fs_electron_fullhide_1d_coupled
