! 电子 1D 全隐格式主驱动。
! 顺序：解包边界 -> 构建 log-four-velocity 网格 -> 初始化电子谱
!       -> 壳层循环：介质/磁场/注入能标 -> 注入/冷却 -> 输运 -> 同步辐射与 SSA 诊断。
! 1D fully implicit electron driver.
! Order: unpack boundary -> build the log-four-velocity grid -> initialize the electron spectrum
!        -> shell loop: medium/field/injection scales -> injection/cooling -> transport -> synchrotron and SSA diagnostics.
subroutine fs_fullhide_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                                   adaptive_substeps,substep_rtol,substep_min,substep_max,thermal_electrons, &
                                   gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    use constants
    use dynamics_density_profile, only: density_profile
    use electron_common
    use electron_injection_profiles, only: source_edges, add_thermal, &
                                           init_powerlaw, init_coord, &
                                           source_coord, thermal_pop, &
                                           log_edges
    use electron_coord_common, only: build_fourvel_grid, dxg_dcoord, &
                                                 coord_fourvel, &
                                                 fourvel_scale
    use electron_radiation_kernel, only: syn_state, nua_fromtau
    use electron_cooling_kernel, only: forward_cooling
    use electron_shell_transport, only: shell_coord_step, &
                                               coord_to_dgamma
    use electron_transport_common, only: dnx_dgamma, &
                                         flux_seq_nonuniform
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in) :: adaptive_substeps,substep_min,substep_max,thermal_electrons
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in) :: substep_rtol
    real(8), intent(out), dimension(Num_gam_e,Num_R) :: dN_gam_e
    real(8), intent(out), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_nu,Num_R) :: P_syn,Seed_syn
    real(8), intent(out), dimension(Num_R) :: V_m,V_c,V_a

    real(8),allocatable,dimension (:) :: dEl,dEl_step,dEL_mean,x,dN_x,x_edge,coord_edge,dxdy_grid, &
                                         dN_full,dN_half,dN_half2,dF1,dEL_mean_step,P_emit_shell,Tau_syn_shell
    real(8),allocatable,dimension (:,:) :: dF_steps,face_coupling
    logical :: is_uniform_density,budget_diag_enabled
    integer :: env_len,env_status,I_face
    character(len=32) :: diag_env
    real(8) :: dDR_xi, n_before_step, n_after_step, inj_step, n_budget, rel_loss_xi_max, source_integral
    real(8) :: adiabatic_integral, l_count_real, step_sum, step_sq_sum, radius_sum, radius_sq_sum
    real(8) :: source_prefactor, coord_scale, dg_gamma_scale, face_coord, face_jac, R_loc, R_Gamma_loc
    real(8) :: beta_Gam, dNe, DB, Gam_e_max, Gam_e_m, Gam_e_m_p, Gam_e_c, dNe_shell, dDR, dDD, f_r
    allocate (dEl(Num_gam_e),dEl_step(Num_gam_e),dEL_mean(Num_gam_e-1),x(Num_gam_e),dN_x(Num_gam_e), &
              dN_full(Num_gam_e), &
              x_edge(Num_gam_e+1),coord_edge(Num_gam_e+1),dxdy_grid(Num_gam_e), &
              dN_half(Num_gam_e),dN_half2(Num_gam_e),dF1(Num_gam_e),dEL_mean_step(Num_gam_e-1), &
              P_emit_shell(Num_nu),Tau_syn_shell(Num_nu))
    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    if (thermal_electrons /= 0) then
        if (f_e <= 0d0 .or. f_e > 1d0) error stop 'thermal electrons require 0 < f_e <= 1'
    end if

    P_syn=0d0
    Seed_syn=0d0
    V_m=0d0
    V_c=0d0
    V_a=0d0

    call init_electron_state()
    d_x=dlog(gam_e(2)/gam_e(1))
    is_uniform_density=(A_star <= 0d0 .and. f_jump == 1d0)
    budget_diag_enabled=.false.
    diag_env=''
    call get_environment_variable('ASGARD_DIAG_1D_BUDGET',diag_env,length=env_len,status=env_status)
    if (env_status == 0 .and. env_len > 0) then
        if (diag_env(1:1) /= '0') budget_diag_enabled=.true.
    end if
    rel_loss_xi_max=0d0
!    factor_adv=Para_sigmaT/(6.0d0*pi*Para_m_energy)

    do I_tobs=2,Num_R
        call prepare_fullhide_shell(I_tobs)
        dEL_mean=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/2d0
        dDR_xi=dDR
        if (adaptive_substeps == 0) then
            call advance_fixed_shell(I_tobs)
        else
            call advance_adaptive_shell()
        end if
        call coord_to_dgamma(Num_gam_e,coord_edge,coord_scale,gam_e,dN_x, &
                                                           dN_gam_e(:,I_tobs))
    end do
    call write_finaldiag()

    deallocate (dEl,dEl_step,dEL_mean,x,dN_x,x_edge,coord_edge,dxdy_grid,dN_full,dN_half,dN_half2,dF1, &
                dEL_mean_step,P_emit_shell,Tau_syn_shell)
    if (budget_diag_enabled) then
        print '(A,1X,ES12.4)', 'BUDGET1D max_rel_loss', rel_loss_xi_max
    end if

contains

    ! 初始化电子能谱和输运坐标；thermal 分支保留公开参数面。
    ! Initialize the electron spectrum and transport coordinate; the thermal branch keeps the public argument contract.
    subroutine init_electron_state()
    implicit real(8)(A-H,O-Z)

        call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-1d0)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-1d0)))
        Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-1d0)
        call electron_gm_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_c=7.7d8/(1d0+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/2d0)
        if (R_Gamma(1) < 1d0) error stop 'fs_fullhide_1d requires initial Gamma >= 1'
        beta_Gam=dsqrt(1d0-1d0/R_Gamma(1)**2)

        call init_fourvel_grid(Gam_e_max_max)
        if (thermal_electrons == 0) then
            call init_coord(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                                                  Num_gam_e,coord_edge,coord_scale,dN_x)
            call coord_to_dgamma(Num_gam_e,coord_edge,coord_scale,gam_e,dN_x, &
                                                               dN_gam_e(:,1))
        else
            call init_powerlaw(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                                      Num_gam_e,gam_e,dN_gam_e(:,1))
            call thermal_pop(Num_gam_e,gam_e,R_Gamma(1)*beta_Gam,Para_N_e_ini*(1d0-f_e), &
                                                 dN_gam_e(:,1))
            dN_x=dN_gam_e(:,1)*gam_e*dxdy_grid
        end if
    end subroutine init_electron_state

    ! four-velocity 坐标把 gamma~1 附近拉开，避免低能端网格过粗。
    ! The four-velocity coordinate resolves the gamma~1 region without over-coarsening the low-energy end.
    subroutine init_fourvel_grid(imodeg_max)
    implicit real(8)(A-H,O-Z)
    integer :: I_grid
    real(8), intent(in) :: imodeg_max

        dg_gamma_scale=fourvel_scale
        coord_scale=dg_gamma_scale*dg_gamma_scale-1d0
        call build_fourvel_grid(Num_gam_e,1d0,tail_factor*imodeg_max, &
                                               dg_gamma_scale,gam_e,coord_edge,x_edge)
        do I_grid=1,Num_gam_e
            dxdy_grid(I_grid)=dxg_dcoord(coord_fourvel,coord_scale, &
                                                      0.5d0*(coord_edge(I_grid)+coord_edge(I_grid+1)))
        end do
    end subroutine init_fourvel_grid

    ! 准备当前壳层的介质、磁场、break frequency 和冷却数组。
    ! Prepare the current shell medium, magnetic field, break frequencies, and cooling array.
    subroutine prepare_fullhide_shell(I_tobs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs

        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/2d0
        if (R_Gamma_loc < 1d0) error stop 'fs_fullhide_1d requires Gamma >= 1'
        beta_Gam=dsqrt(1d0-1d0/R_Gamma_loc**2)
        call density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-1d0)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
        call electron_gm_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(1d0-p)/(Gam_e_max**(1d0-p)-Gam_e_m**(1d0-p))
        Gam_e_c=7.7d8*(1d0+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        dNe_shell=dNe

        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1d0/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dxdy_grid
        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))

        call syn_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_emit_shell, &
                                    P_syn(:,I_tobs),Seed_syn(:,I_tobs),Tau_syn_shell)
        call nua_fromtau(Num_nu,V_seed,Tau_syn_shell,temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))

        call forward_cooling(1,index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                             R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,1,n_threads,gam_e,V_seed, &
                             P_syn(:,I_tobs),Seed_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl_step,dEl)
    end subroutine prepare_fullhide_shell

    ! 最后一个输出点没有后续推进，只刷新与最终电子谱一致的辐射诊断。
    ! The final output point has no following advance, so refresh diagnostics from the final electron spectrum.
    subroutine write_finaldiag()
    implicit real(8)(A-H,O-Z)

        R_loc=R(Num_R)
        R_Gamma_loc=R_Gamma(Num_R)
        if (R_Gamma_loc < 1d0) error stop 'fs_fullhide_1d requires Gamma >= 1'
        beta_Gam=dsqrt(1d0-1d0/R_Gamma_loc**2)
        call density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-1d0)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
        call electron_gm_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_c=7.7d8*(1d0+z)/R_Gamma_loc/DB**2/R_Tobs(Num_R)
        V_m(Num_R)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        V_c(Num_R)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        call syn_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,Num_R),V_seed,P_emit_shell, &
                                    P_syn(:,Num_R),Seed_syn(:,Num_R),Tau_syn_shell)
        call nua_fromtau(Num_nu,V_seed,Tau_syn_shell,temp)
        V_a(Num_R)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
    end subroutine write_finaldiag

    ! 固定子步路径：均匀介质可合并源项，非均匀介质逐子步更新介质。
    ! Fixed-substep path: uniform media can merge sources, while nonuniform media update the medium per substep.
    subroutine advance_fixed_shell(I_tobs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs

        if (dDD <= 0d0) error stop 'fs_fullhide_1d requires increasing radius grid'
        if (dDR <= 0d0) error stop 'fs_fullhide_1d requires positive cooling substep width'
        L1=max(100,min(1000,int(dDD/dDR)))
        dDR=dDD/dble(L1)
        if (is_uniform_density .and. thermal_electrons == 0) then
            call advance_uniform_shell(L1)
        else
            call advance_general_shell(I_tobs,L1)
        end if
    end subroutine advance_fixed_shell

    subroutine advance_uniform_shell(L1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: L1

        allocate(dF_steps(Num_gam_e,1),face_coupling(Num_gam_e-1,1))
        l_count_real=dble(L1)
        step_sum=l_count_real*(l_count_real+1d0)/2d0
        step_sq_sum=l_count_real*(l_count_real+1d0)*(2d0*l_count_real+1d0)/6d0
        radius_sum=l_count_real*R_loc+dDR*step_sum
        radius_sq_sum=l_count_real*R_loc*R_loc+2d0*R_loc*dDR*step_sum+dDR*dDR*step_sq_sum
        source_prefactor=4d0/3d0*pi*dNe*f_e*Gam_e_m_p
        source_integral=dDR*source_prefactor*(3d0*radius_sq_sum+dDR*(3d0*radius_sum+l_count_real*dDR))
        adiabatic_integral=0d0
        do L=1,L1
            adiabatic_integral=adiabatic_integral+dDR/(R_loc+dDR*dble(L))
        end do
        call source_coord(Num_gam_e,coord_edge,coord_scale, &
                                                               Gam_e_m,Gam_e_max,source_integral,p, &
                                                               dF_steps(:,1))
        do I_face=1,Num_gam_e-1
            face_coord=coord_edge(I_face+1)
            face_jac=dxg_dcoord(coord_fourvel,coord_scale,face_coord)
            face_coupling(I_face,1)=(dDD*(dEl(I_face)+dEl(I_face+1))/2d0+adiabatic_integral)/face_jac
        end do
        if (budget_diag_enabled) then
            n_before_step=sum(dN_x*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
            inj_step=sum(dF_steps(:,1)*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
        end if
        call flux_seq_nonuniform(Num_gam_e,coord_edge,face_coupling(:,1),dF_steps(:,1), &
                                                              dN_x,x,.true.)
        if (budget_diag_enabled) then
            n_after_step=sum(x*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
            n_budget=n_before_step+inj_step
            if (n_budget > 0d0) rel_loss_xi_max=max(rel_loss_xi_max,max(0d0,(n_budget-n_after_step)/n_budget))
        end if
        dN_x=x
        deallocate(dF_steps,face_coupling)
    end subroutine advance_uniform_shell

    subroutine advance_general_shell(I_tobs,L1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs,L1

        do L=1,L1
            R_loc=R_loc+dDR
            call density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
            DB_step=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-1d0)))
            Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
            temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
            call electron_gm_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
            Gam_e_m_p_step=(1d0-p)/(Gam_e_max_step**(1d0-p)-Gam_e_m_step**(1d0-p))
            Q=4d0/3d0*pi*(3d0*R_loc**2+dDR*(3d0*R_loc+dDR))*dNe*f_e*Gam_e_m_p_step
            call source_coord(Num_gam_e,coord_edge,coord_scale, &
                                                                   Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
            if (thermal_electrons /= 0) then
                call add_thermal(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam, &
                                                      Q*(1d0-f_e)/(f_e*Gam_e_m_p_step),dF1)
            end if
            if (dNe_shell > 0d0) then
                dEl_step=dEl*(dNe/dNe_shell)
            else
                dEl_step=dEl
            end if
            if (budget_diag_enabled) then
                n_before_step=sum(dN_x*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
                inj_step=dDR*sum(dF1*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
            end if
            call shell_coord_step(Num_gam_e,dDR,coord_edge,coord_scale, &
                                                      dEl_step,1d0/R_loc,dF1,dN_x,x)
            if (budget_diag_enabled) then
                n_after_step=sum(x*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
                n_budget=n_before_step+inj_step
                if (n_budget > 0d0) rel_loss_xi_max=max(rel_loss_xi_max,max(0d0,(n_budget-n_after_step)/n_budget))
                if (I_tobs <= 6 .and. L == L1) then
                    print '(A,1X,I4,1X,ES12.4,1X,ES12.4,1X,ES12.4)', &
                          'BUDGET1D shell', I_tobs, n_before_step, inj_step, n_after_step
                end if
            end if
            dN_x=x
        end do
    end subroutine advance_general_shell

    ! 自适应子步用整步/半步差估计局部误差，不做经验平滑。
    ! Adaptive substepping estimates local error from full-step versus half-step updates, with no heuristic smoothing.
    subroutine advance_adaptive_shell()
    implicit real(8)(A-H,O-Z)

        dR_min=dDD/max(substep_max,1)
        dR_max=dDD/max(substep_min,1)
        dR_try=min(min(dDR_xi,dR_max),dDD)
        dR_left=dDD

        do while (dR_left > 0d0)
            dR_try=min(dR_try,dR_left)
            if (dR_left < 1.5d0*dR_try) dR_try=dR_left

            R_full=R_loc+dR_try
            if (is_uniform_density) then
                dNe_full=dNe
            else
                call density_profile(A_star,dNe_ISM,R_full,R0,1,R_tr,f_jump,f_wide,dNe_full)
            end if
            DB_full=0.39d0*dsqrt(Epsilon_b*dNe_full*(R_Gamma_loc*(R_Gamma_loc-1d0)))
            Gam_e_max_full=3d0*Para_m_energy/dsqrt(8d0*DB_full*Para_e**3)
            temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
            call electron_gm_exact(p,temp_gam,Gam_e_max_full,Gam_e_m_full)
            Gam_e_m_p_full=(1d0-p)/(Gam_e_max_full**(1d0-p)-Gam_e_m_full**(1d0-p))
            call electron_injection_prefactor(R_full,dR_try,dNe_full,f_e,Gam_e_m_p_full,Q)
            call source_coord(Num_gam_e,coord_edge,coord_scale, &
                                                                   Gam_e_m_full,Gam_e_max_full,Q,p,dF1)
            if (thermal_electrons /= 0) then
                thermal_count=Q*(1d0-f_e)/(f_e*Gam_e_m_p_full)
                call add_thermal(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam,thermal_count,dF1)
            end if
            if (dNe_shell > 0d0) then
                dEl_step=dEl*(dNe_full/dNe_shell)
            else
                dEl_step=dEl
            end if
            call shell_coord_step(Num_gam_e,dR_try,coord_edge,coord_scale, &
                                                      dEl_step,1d0/R_full,dF1,dN_x,dN_full)

            dR_half=0.5d0*dR_try
            R_half=R_loc+dR_half
            if (is_uniform_density) then
                dNe_half=dNe
            else
                call density_profile(A_star,dNe_ISM,R_half,R0,1,R_tr,f_jump,f_wide,dNe_half)
            end if
            DB_half=0.39d0*dsqrt(Epsilon_b*dNe_half*(R_Gamma_loc*(R_Gamma_loc-1d0)))
            Gam_e_max_half=3d0*Para_m_energy/dsqrt(8d0*DB_half*Para_e**3)
            temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
            call electron_gm_exact(p,temp_gam,Gam_e_max_half,Gam_e_m_half)
            Gam_e_m_p_half=(1d0-p)/(Gam_e_max_half**(1d0-p)-Gam_e_m_half**(1d0-p))
            call electron_injection_prefactor(R_half,dR_half,dNe_half,f_e,Gam_e_m_p_half,Q)
            call source_coord(Num_gam_e,coord_edge,coord_scale, &
                                                                   Gam_e_m_half,Gam_e_max_half,Q,p,dF1)
            if (thermal_electrons /= 0) then
                thermal_count=Q*(1d0-f_e)/(f_e*Gam_e_m_p_half)
                call add_thermal(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam,thermal_count,dF1)
            end if
            if (dNe_shell > 0d0) then
                dEl_step=dEl*(dNe_half/dNe_shell)
            else
                dEl_step=dEl
            end if
            call shell_coord_step(Num_gam_e,dR_half,coord_edge,coord_scale, &
                                                      dEl_step,1d0/R_half,dF1,dN_x,dN_half)

            call electron_injection_prefactor(R_full,dR_half,dNe_full,f_e,Gam_e_m_p_full,Q)
            call source_coord(Num_gam_e,coord_edge,coord_scale, &
                                                                   Gam_e_m_full,Gam_e_max_full,Q,p,dF1)
            if (thermal_electrons /= 0) then
                thermal_count=Q*(1d0-f_e)/(f_e*Gam_e_m_p_full)
                call add_thermal(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam,thermal_count,dF1)
            end if
            if (dNe_shell > 0d0) then
                dEl_step=dEl*(dNe_full/dNe_shell)
            else
                dEl_step=dEl
            end if
            call shell_coord_step(Num_gam_e,dR_half,coord_edge,coord_scale, &
                                                      dEl_step,1d0/R_full,dF1,dN_half,dN_half2)

            call electron_relerr_max(Num_gam_e,dN_half2,dN_full,step_error)
            if (step_error <= substep_rtol .or. dR_try <= dR_min) then
                dN_x=dN_half2
                R_loc=R_full
                dNe=dNe_full
                dR_left=dR_left-dR_try
                if (dR_left > 0d0) then
                    if (step_error <= 0.25d0*substep_rtol) then
                        dR_try=min(2d0*dR_try,dR_max)
                    end if
                end if
            else
                dR_try=max(0.5d0*dR_try,dR_min)
            end if
        end do
    end subroutine advance_adaptive_shell
end subroutine fs_fullhide_1d

! 1D joint-coupling electron pass: use externally closed photon and secondary pair source fields.
! Joint electron-photon pass: same fixed-substep fullhide shell transport, but cooling seed and
! secondary source are supplied by the shell-level photon/hadronic feedback stage.
subroutine fs_fullhide_coupled(Boundary,R_Tobs,R_Gamma,R,V_seed,Seed_cooling,sec_source,n,Num_nu,Num_R, &
                                           Num_gam_e,index_Y,index_syn_intger,n_threads,adaptive_substeps, &
                                           substep_rtol,substep_min,substep_max,thermal_electrons, &
                                           gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    use constants
    use dynamics_density_profile, only: density_profile
    use electron_common
    use electron_injection_profiles, only: source_edges, add_thermal, &
                                           log_edges
    use electron_radiation_kernel, only: syn_state, nua_fromtau
    use electron_ic_kernel, only: electron_ic_budget
    use electron_cooling_kernel, only: forward_cooling
    use electron_transport_common, only: dnx_dgamma, fullhide_step
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in) :: adaptive_substeps,substep_min,substep_max,thermal_electrons
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu,Num_R) :: Seed_cooling
    real(8), intent(in), dimension(Num_gam_e,Num_R) :: sec_source
    real(8), intent(in) :: substep_rtol
    real(8), intent(out), dimension(Num_gam_e,Num_R) :: dN_gam_e
    real(8), intent(out), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_nu,Num_R) :: P_syn,Seed_syn
    real(8), intent(out), dimension(Num_R) :: V_m,V_c,V_a
    real(8), allocatable, dimension(:) :: dEl,dEL_mean,dEL_mean_step,cooling_aux,x,dN_x,x_edge,dF1
    logical :: budget_diag_enabled
    integer :: env_len,env_status
    character(len=32) :: diag_env
    real(8) :: n_before_step,n_after_step,inj_step,n_budget,rel_loss_xi_max,thermal_count
    real(8) :: R_loc,R_Gamma_loc,beta_Gam,dNe,DB,Gam_e_max,Gam_e_m,Gam_e_m_p,Gam_e_c,dNe_shell,dDR,dDD,f_r

    if (adaptive_substeps /= 0) error stop 'fs_fullhide_coupled requires fixed shell substeps'
    if (index_Y /= 1) error stop 'fs_fullhide_coupled requires index_Y=1'
    if (substep_min < 1 .or. substep_max < 1) error stop 'fs_fullhide_coupled requires positive substep bounds'
    if (substep_rtol <= 0d0) error stop 'fs_fullhide_coupled requires positive substep_rtol'
    allocate(dEl(Num_gam_e),dEL_mean(Num_gam_e-1),dEL_mean_step(Num_gam_e-1),cooling_aux(Num_gam_e), &
             x(Num_gam_e),dN_x(Num_gam_e),x_edge(Num_gam_e+1),dF1(Num_gam_e))

    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    if (thermal_electrons /= 0) then
        if (f_e <= 0d0 .or. f_e > 1d0) error stop 'thermal electrons require 0 < f_e <= 1'
    end if

    P_syn=0d0
    Seed_syn=0d0
    V_m=0d0
    V_c=0d0
    V_a=0d0

    call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)
    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-1d0)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-1d0)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-1d0)
    call electron_gm_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c=7.7d8/(1d0+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/2d0)
    if (R_Gamma(1) < 1d0) error stop 'fs_fullhide_coupled requires initial Gamma >= 1'
    beta_Gam=dsqrt(1d0-1d0/R_Gamma(1)**2)
    if (thermal_electrons == 0) then
        call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                          imodelog,gam_e,dN_x,x_edge)
        call dnx_dgamma(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,1))
    else
        call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                          imodeg,gam_e,dN_gam_e(:,1),thermal_electrons=thermal_electrons, &
                                          f_e=f_e,four_v=R_Gamma(1)*beta_Gam)
        dN_x=dN_gam_e(:,1)*gam_e
        call log_edges(Num_gam_e,gam_e,x_edge)
    end if

    d_x=dlog(gam_e(2)/gam_e(1))
    budget_diag_enabled=.false.
    diag_env=''
    call get_environment_variable('ASGARD_DIAG_1D_BUDGET',diag_env,length=env_len,status=env_status)
    if (env_status == 0 .and. env_len > 0) then
        if (diag_env(1:1) /= '0') budget_diag_enabled=.true.
    end if
    rel_loss_xi_max=0d0

    do I_tobs=2,Num_R
        call prepare_coupled_shell(I_tobs)
        dEL_mean=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/2d0
        call advance_coupled_shell(I_tobs)
        call dnx_dgamma(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
    end do
    call write_finaldiag()

    deallocate(dEl,dEL_mean,dEL_mean_step,cooling_aux,x,dN_x,x_edge,dF1)
    if (budget_diag_enabled) print '(A,1X,ES12.4)', 'BUDGET1D coupled max_rel_loss', rel_loss_xi_max

contains

    subroutine prepare_coupled_shell(I_tobs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs
    real(8), dimension(Num_nu,1) :: Seed_ssa_column
    real(8), dimension(Num_gam_e,1) :: cooling_aux_column,dEl_column
    real(8), dimension(Num_nu) :: P_emit_tmp,Tau_syn_tmp

        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/2d0
        if (R_Gamma_loc < 1d0) error stop 'fs_fullhide_coupled requires Gamma >= 1'
        beta_Gam=dsqrt(1d0-1d0/R_Gamma_loc**2)
        call density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-1d0)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
        call electron_gm_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(1d0-p)/(Gam_e_max**(1d0-p)-Gam_e_m**(1d0-p))
        Gam_e_c=7.7d8*(1d0+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        dNe_shell=dNe
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1d0/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e
        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        call syn_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_emit_tmp,P_syn(:,I_tobs), &
                                    Seed_syn(:,I_tobs),Tau_syn_tmp)
        call nua_fromtau(Num_nu,V_seed,Tau_syn_tmp,temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        Gam_e_max_cool=Gam_e_max
        call electron_ic_budget(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_cooling(:,I_tobs),cooling_aux)
        Seed_ssa_column(:,1)=Seed_syn(:,I_tobs)
        cooling_aux_column(:,1)=cooling_aux
        call forward_cooling(2,index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max_cool, &
                             R_loc,R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,1,n_threads, &
                             gam_e,V_seed,Seed_ssa_column,Seed_ssa_column,Seed_ssa_column, &
                             cooling_aux_column,dEl_column)
        dEl=dEl_column(:,1)
    end subroutine prepare_coupled_shell

    ! 最后一个输出点没有后续推进，只刷新与最终电子谱一致的辐射诊断。
    ! The final output point has no following advance, so refresh diagnostics from the final electron spectrum.
    subroutine write_finaldiag()
    implicit real(8)(A-H,O-Z)
    real(8), dimension(Num_nu) :: Ptmp,Tautmp

        R_loc=R(Num_R)
        R_Gamma_loc=R_Gamma(Num_R)
        if (R_Gamma_loc < 1d0) error stop 'fs_fullhide_coupled requires Gamma >= 1'
        beta_Gam=dsqrt(1d0-1d0/R_Gamma_loc**2)
        call density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-1d0)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
        call electron_gm_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_c=7.7d8*(1d0+z)/R_Gamma_loc/DB**2/R_Tobs(Num_R)
        V_m(Num_R)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        V_c(Num_R)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        call syn_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,Num_R),V_seed,Ptmp, &
                                    P_syn(:,Num_R),Seed_syn(:,Num_R),Tautmp)
        call nua_fromtau(Num_nu,V_seed,Tautmp,temp)
        V_a(Num_R)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
    end subroutine write_finaldiag

    subroutine advance_coupled_shell(I_tobs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs

        if (dDD <= 0d0) error stop 'fs_fullhide_coupled requires increasing radius grid'
        if (dDR <= 0d0) error stop 'fs_fullhide_coupled requires positive cooling substep width'
        L1=max(100,min(1000,int(dDD/dDR)))
        dDR=dDD/dble(L1)

        do L=1,L1
            R_loc=R_loc+dDR
            call density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
            DB_step=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-1d0)))
            Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
            temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
            call electron_gm_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
            Gam_e_m_p_step=(1d0-p)/(Gam_e_max_step**(1d0-p)-Gam_e_m_step**(1d0-p))
            Q=4d0/3d0*pi*(3d0*R_loc**2+dDR*(3d0*R_loc+dDR))*dNe*f_e*Gam_e_m_p_step
            call source_edges(Num_gam_e,x_edge,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
            dF1=dF1+sec_source(:,I_tobs)
            if (thermal_electrons /= 0) then
                thermal_count=Q*(1d0-f_e)/(f_e*Gam_e_m_p_step)
                call add_thermal(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam,thermal_count,dF1)
            end if
            if (dNe_shell > 0d0) then
                dEL_mean_step=dEL_mean*(dNe/dNe_shell)
            else
                dEL_mean_step=dEL_mean
            end if
            if (budget_diag_enabled) then
                n_before_step=sum(dN_x)*d_x
                inj_step=dDR*sum(dF1)*d_x
            end if
            call fullhide_step(Num_gam_e,R_loc,dDR,d_x,dEL_mean_step,dF1,dN_x,x)
            if (budget_diag_enabled) then
                n_after_step=sum(x)*d_x
                n_budget=n_before_step+inj_step
                if (n_budget > 0d0) rel_loss_xi_max=max(rel_loss_xi_max,max(0d0,(n_budget-n_after_step)/n_budget))
            end if
            dN_x=x
        end do
    end subroutine advance_coupled_shell
end subroutine fs_fullhide_coupled
