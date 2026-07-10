!f2py: skip
module fullhide_drivers
contains

! 电子 1D 全隐格式主驱动。
! 顺序：解包边界 -> 构建 log-four-velocity 网格 -> 初始化电子谱
!       -> 壳层循环：介质/磁场/注入能标 -> 注入/冷却 -> 输运 -> 同步辐射与 SSA 诊断。
! 1D fully implicit electron driver.
! Order: unpack boundary -> build the log-four-velocity grid -> initialize the electron spectrum
!        -> shell loop: medium/field/injection scales -> injection/cooling -> transport -> synchrotron and SSA diagnostics.
subroutine fullhide_core(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                                  adaptive_substeps,substep_rtol,substep_min,substep_max,thermal_electrons, &
                                  grid_mode,grid_top,gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    use constants
    use dynamics_density_profile, only: density_profile, uniform_density
    use electron_common
    use electron_injection_profiles, only: source_edges, add_thermal, add_thermal_edges, &
                                           init_powerlaw, init_edges, init_coord, &
                                           source_coord, thermal_pop, &
                                           log_edges
    use electron_coord_common, only: build_fourvel_grid, dxg_dcoord, &
                                                 coord_fourvel, &
                                                 fourvel_scale
    use electron_radiation_kernel, only: syn_state, nua_fromtau
    use electron_cooling_kernel, only: forward_cooling
    use electron_shell_transport, only: shellstep_cached, &
                                               coord_to_dgamma
    use electron_transport_common, only: dnx_dgamma, &
                                         flux_seq_nonuniform
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in) :: adaptive_substeps,substep_min,substep_max,thermal_electrons,grid_mode
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in) :: substep_rtol,grid_top
    real(8), intent(out), dimension(Num_gam_e,Num_R) :: dN_gam_e
    real(8), intent(out), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_nu,Num_R) :: P_syn,Seed_syn
    real(8), intent(out), dimension(Num_R) :: V_m,V_c,V_a

    real(8),allocatable,dimension (:) :: dEl,dEl_step,x,dN_x,x_edge,coord_edge,dxdy_grid,face_invjac, &
                                         dN_full,dN_half,dN_half2,dF1,face_disp,P_emit_shell,Tau_syn_shell
    logical :: is_uniform_density,is_log
    integer :: I_face
    real(8) :: dDR_xi,source_integral
    real(8) :: adiabatic_integral
    real(8) :: coord_scale, dg_gamma_scale, grid_max, R_loc, R_Gamma_loc
    real(8) :: beta_Gam, dNe, DB, Gam_e_max, Gam_e_m, Gam_e_m_p, Gam_e_c, dNe_shell, dDR, dDD, f_r,neinit
    allocate (dEl(Num_gam_e),dEl_step(Num_gam_e),x(Num_gam_e),dN_x(Num_gam_e), &
              dN_full(Num_gam_e), &
              x_edge(Num_gam_e+1),coord_edge(Num_gam_e+1),dxdy_grid(Num_gam_e), &
              face_invjac(Num_gam_e-1), &
              dN_half(Num_gam_e),dN_half2(Num_gam_e),dF1(Num_gam_e),face_disp(Num_gam_e-1), &
              P_emit_shell(Num_nu),Tau_syn_shell(Num_nu))
    call electron_unpack_boundary(Boundary,n,Eta_0,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    if (thermal_electrons /= 0) then
        if (f_e <= 0d0 .or. f_e > 1d0) error stop 'thermal electrons require 0 < f_e <= 1'
    end if

    P_syn=0d0
    Seed_syn=0d0
    V_m=0d0
    V_c=0d0
    V_a=0d0

    is_log=grid_mode /= 0
    call init_electron_state()
    is_uniform_density=uniform_density(A_star,f_jump)
    do I_tobs=2,Num_R
        call prepare_fullhide_shell(I_tobs)
        dDR_xi=dDR
        if (adaptive_substeps == 0) then
            call advance_fixed_shell()
        else
            call advance_adaptive_shell()
        end if
        call project_density(I_tobs)
    end do
    call write_finaldiag()

    deallocate (dEl,dEl_step,x,dN_x,x_edge,coord_edge,dxdy_grid,face_invjac,dN_full,dN_half,dN_half2,dF1, &
                face_disp,P_emit_shell,Tau_syn_shell)

contains

    ! 初始化电子能谱和输运坐标；thermal 分支保留公开参数面。
    ! Initialize the electron spectrum and transport coordinate; the thermal branch keeps the public argument contract.
    subroutine init_electron_state()
    implicit real(8)(A-H,O-Z)

        call electron_initial_density(A_star,dNe_ISM,R(1),R0,R_tr,f_jump,f_wide,dNe,neinit)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-1d0)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        if (is_log) then
            call electron_scanmax(Num_R,R_Gamma,R,Epsilon_b,A_star,dNe_ISM,R0,R_tr,f_jump,f_wide,Gam_e_max_max)
        else
            DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-1d0)))
            Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
        end if
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-1d0)
        call electron_gm_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_c=7.7d8/(1d0+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/2d0)
        if (R_Gamma(1) < 1d0) error stop 'fs_fullhide_1d requires initial Gamma >= 1'
        beta_Gam=dsqrt(1d0-1d0/R_Gamma(1)**2)

        if (is_log) then
            grid_max=max(grid_top,tail_factor*Gam_e_max_max)
            call electron_loggrid(Num_gam_e,grid_max,gam_e,x_edge)
            coord_edge=x_edge
            dxdy_grid=1d0
            face_invjac=1d0
        else
            call init_fourvel_grid(Gam_e_max_max)
        end if
        call init_transport()
    end subroutine init_electron_state

    ! 按输运坐标初始化保守 cell density；BH log 路径不经过点值投影。
    ! Initialize conservative cell density in the selected coordinate; the BH log path avoids point projection.
    subroutine init_transport()
    implicit real(8)(A-H,O-Z)

        if (is_log) then
            call init_edges(f_e*neinit,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                         Num_gam_e,x_edge,dN_x)
            if (thermal_electrons /= 0) then
                call transport_thermal(R_Gamma(1)*beta_Gam,neinit*(1d0-f_e),dN_x)
            end if
            call dnx_dgamma(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,1))
        else if (thermal_electrons == 0) then
            call init_coord(f_e*neinit,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                          Num_gam_e,coord_edge,coord_scale,dN_x)
            call coord_to_dgamma(Num_gam_e,coord_edge,coord_scale,gam_e,dN_x, &
                                                dN_gam_e(:,1))
        else
            call init_powerlaw(f_e*neinit,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                               Num_gam_e,gam_e,dN_gam_e(:,1))
            call thermal_pop(Num_gam_e,gam_e,R_Gamma(1)*beta_Gam, &
                             neinit*(1d0-f_e),dN_gam_e(:,1))
            dN_x=dN_gam_e(:,1)*gam_e*dxdy_grid
        end if
    end subroutine init_transport

    ! 在 log-gamma 与 four-velocity 坐标上分别构造同一物理注入源。
    ! Build the same physical injection source in log-gamma or four-velocity coordinates.
    subroutine transport_source(gmin,gmax,qsrc,dst)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: gmin,gmax,qsrc
    real(8), intent(out), dimension(Num_gam_e) :: dst

        if (is_log) then
            call source_edges(Num_gam_e,x_edge,gmin,gmax,qsrc,p,dst)
        else
            call source_coord(Num_gam_e,coord_edge,coord_scale,gmin,gmax,qsrc,p,dst)
        end if
    end subroutine transport_source

    ! 热电子源在 BH log cell 上按边界积分，旧 four-velocity 路径保留原语义。
    ! Integrate the thermal source over BH log cells while retaining legacy four-velocity semantics.
    subroutine transport_thermal(fourv,count,dst)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: fourv,count
    real(8), intent(inout), dimension(Num_gam_e) :: dst

        if (is_log) then
            call add_thermal_edges(Num_gam_e,x_edge,fourv,count,dst)
        else
            call add_thermal(Num_gam_e,gam_e,fourv,count,dst)
        end if
    end subroutine transport_thermal

    ! 将保守 cell density 只投影一次到当前输出壳层。
    ! Project conservative cell density once to the current output shell.
    subroutine project_density(ir)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: ir

        if (is_log) then
            call dnx_dgamma(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,ir))
        else
            call coord_to_dgamma(Num_gam_e,coord_edge,coord_scale,gam_e,dN_x, &
                                 dN_gam_e(:,ir))
        end if
    end subroutine project_density

    ! four-velocity 坐标把 gamma~1 附近拉开，避免低能端网格过粗。
    ! The four-velocity coordinate resolves the gamma~1 region without over-coarsening the low-energy end.
    subroutine init_fourvel_grid(imodeg_max)
    implicit real(8)(A-H,O-Z)
    integer :: I_grid,I_face
    real(8), intent(in) :: imodeg_max

        dg_gamma_scale=fourvel_scale
        coord_scale=dg_gamma_scale*dg_gamma_scale-1d0
        call build_fourvel_grid(Num_gam_e,1d0,tail_factor*imodeg_max, &
                                               dg_gamma_scale,gam_e,coord_edge,x_edge)
        do I_grid=1,Num_gam_e
            dxdy_grid(I_grid)=dxg_dcoord(coord_fourvel,coord_scale, &
                                                      0.5d0*(coord_edge(I_grid)+coord_edge(I_grid+1)))
        end do
        do I_face=1,Num_gam_e-1
            face_invjac(I_face)=1d0/dxg_dcoord(coord_fourvel,coord_scale,coord_edge(I_face+1))
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
        if (.not. is_log) dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dxdy_grid
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
    subroutine advance_fixed_shell()
    implicit real(8)(A-H,O-Z)

        if (dDD <= 0d0) error stop 'fs_fullhide_1d requires increasing radius grid'
        if (dDR <= 0d0) error stop 'fs_fullhide_1d requires positive cooling substep width'
        L1=max(100,min(1000,int(dDD/dDR)))
        dDR=dDD/dble(L1)
        if (is_uniform_density .and. thermal_electrons == 0) then
            call advance_uniform_shell(L1)
        else
            call advance_general_shell(L1)
        end if
    end subroutine advance_fixed_shell

    subroutine advance_uniform_shell(L1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: L1

        call electron_injection_prefactor(R_loc,dDD,dNe,f_e,Gam_e_m_p,source_integral)
        source_integral=dDD*source_integral
        adiabatic_integral=0d0
        do L=1,L1
            adiabatic_integral=adiabatic_integral+dDR/(R_loc+dDR*dble(L))
        end do
        call transport_source(Gam_e_m,Gam_e_max,source_integral,dF1)
        do I_face=1,Num_gam_e-1
            face_disp(I_face)=(dDD*(dEl(I_face)+dEl(I_face+1))/2d0+adiabatic_integral)*face_invjac(I_face)
        end do
        call flux_seq_nonuniform(Num_gam_e,coord_edge,face_disp,dF1, &
                                                              dN_x,x,.true.)
        dN_x=x
    end subroutine advance_uniform_shell

    subroutine advance_general_shell(L1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: L1

        do L=1,L1
            R_loc=R_loc+dDR
            call density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
            DB_step=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-1d0)))
            Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
            temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
            call electron_gm_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
            Gam_e_m_p_step=(1d0-p)/(Gam_e_max_step**(1d0-p)-Gam_e_m_step**(1d0-p))
            call electron_injection_prefactor(R_loc-dDR,dDR,dNe,f_e,Gam_e_m_p_step,Q)
            call transport_source(Gam_e_m_step,Gam_e_max_step,Q,dF1)
            if (thermal_electrons /= 0) then
                call transport_thermal(R_Gamma_loc*beta_Gam,Q*(1d0-f_e)/(f_e*Gam_e_m_p_step),dF1)
            end if
            if (dNe_shell > 0d0) then
                dEl_step=dEl*(dNe/dNe_shell)
            else
                dEl_step=dEl
            end if
            call shellstep_cached(Num_gam_e,dDR,coord_edge,face_invjac, &
                                                      dEl_step,1d0/R_loc,dF1,dN_x,x)
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
            call electron_injection_prefactor(R_loc,dR_try,dNe_full,f_e,Gam_e_m_p_full,Q)
            call transport_source(Gam_e_m_full,Gam_e_max_full,Q,dF1)
            if (thermal_electrons /= 0) then
                thermal_count=Q*(1d0-f_e)/(f_e*Gam_e_m_p_full)
                call transport_thermal(R_Gamma_loc*beta_Gam,thermal_count,dF1)
            end if
            if (dNe_shell > 0d0) then
                dEl_step=dEl*(dNe_full/dNe_shell)
            else
                dEl_step=dEl
            end if
            call shellstep_cached(Num_gam_e,dR_try,coord_edge,face_invjac, &
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
            call electron_injection_prefactor(R_loc,dR_half,dNe_half,f_e,Gam_e_m_p_half,Q)
            call transport_source(Gam_e_m_half,Gam_e_max_half,Q,dF1)
            if (thermal_electrons /= 0) then
                thermal_count=Q*(1d0-f_e)/(f_e*Gam_e_m_p_half)
                call transport_thermal(R_Gamma_loc*beta_Gam,thermal_count,dF1)
            end if
            if (dNe_shell > 0d0) then
                dEl_step=dEl*(dNe_half/dNe_shell)
            else
                dEl_step=dEl
            end if
            call shellstep_cached(Num_gam_e,dR_half,coord_edge,face_invjac, &
                                                      dEl_step,1d0/R_half,dF1,dN_x,dN_half)

            call electron_injection_prefactor(R_half,dR_half,dNe_full,f_e,Gam_e_m_p_full,Q)
            call transport_source(Gam_e_m_full,Gam_e_max_full,Q,dF1)
            if (thermal_electrons /= 0) then
                thermal_count=Q*(1d0-f_e)/(f_e*Gam_e_m_p_full)
                call transport_thermal(R_Gamma_loc*beta_Gam,thermal_count,dF1)
            end if
            if (dNe_shell > 0d0) then
                dEl_step=dEl*(dNe_full/dNe_shell)
            else
                dEl_step=dEl
            end if
            call shellstep_cached(Num_gam_e,dR_half,coord_edge,face_invjac, &
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
end subroutine fullhide_core

! 1D joint-coupling electron pass: use externally closed photon and secondary pair source fields.
! Joint electron-photon pass: same fixed-substep fullhide shell transport, but cooling seed and
! secondary source are supplied by the shell-level photon/hadronic feedback stage.
subroutine coupled_core(Boundary,R_Tobs,R_Gamma,R,V_seed,Seed_cooling,sec_source,n,Num_nu,Num_R, &
                                 Num_gam_e,index_Y,index_syn_intger,n_threads,adaptive_substeps, &
                                 substep_rtol,substep_min,substep_max,thermal_electrons,grid_mode,grid_top, &
                                 gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    use constants
    use dynamics_density_profile, only: density_profile
    use electron_common
    use electron_injection_profiles, only: source_edges, add_thermal, add_thermal_edges, init_edges, &
                                           log_edges
    use electron_radiation_kernel, only: syn_state, nua_fromtau
    use electron_ic_kernel, only: electron_ic_budget
    use electron_cooling_kernel, only: forward_cooling
    use electron_transport_common, only: dnx_dgamma, fullhide_step
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in) :: adaptive_substeps,substep_min,substep_max,thermal_electrons,grid_mode
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu,Num_R) :: Seed_cooling
    real(8), intent(in), dimension(Num_gam_e,Num_R) :: sec_source
    real(8), intent(in) :: substep_rtol,grid_top
    real(8), intent(out), dimension(Num_gam_e,Num_R) :: dN_gam_e
    real(8), intent(out), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_nu,Num_R) :: P_syn,Seed_syn
    real(8), intent(out), dimension(Num_R) :: V_m,V_c,V_a
    real(8), allocatable, dimension(:) :: dEl,dEL_mean,dEL_mean_step,cooling_aux,x,dN_x,x_edge,dF1
    real(8) :: thermal_count,grid_max
    real(8) :: R_loc,R_Gamma_loc,beta_Gam,dNe,DB,Gam_e_max,Gam_e_m,Gam_e_m_p,Gam_e_c,dNe_shell,dDR,dDD,f_r

    if (adaptive_substeps /= 0) error stop 'fs_fullhide_coupled requires fixed shell substeps'
    if (index_Y /= 1) error stop 'fs_fullhide_coupled requires index_Y=1'
    if (substep_min < 1 .or. substep_max < 1) error stop 'fs_fullhide_coupled requires positive substep bounds'
    if (substep_rtol <= 0d0) error stop 'fs_fullhide_coupled requires positive substep_rtol'
    allocate(dEl(Num_gam_e),dEL_mean(Num_gam_e-1),dEL_mean_step(Num_gam_e-1),cooling_aux(Num_gam_e), &
             x(Num_gam_e),dN_x(Num_gam_e),x_edge(Num_gam_e+1),dF1(Num_gam_e))

    call electron_unpack_boundary(Boundary,n,Eta_0,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    if (thermal_electrons /= 0) then
        if (f_e <= 0d0 .or. f_e > 1d0) error stop 'thermal electrons require 0 < f_e <= 1'
    end if

    P_syn=0d0
    Seed_syn=0d0
    V_m=0d0
    V_c=0d0
    V_a=0d0

    call electron_initial_density(A_star,dNe_ISM,R(1),R0,R_tr,f_jump,f_wide,dNe,Para_N_e_ini)
    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-1d0)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    if (grid_mode /= 0) then
        call electron_scanmax(Num_R,R_Gamma,R,Epsilon_b,A_star,dNe_ISM,R0,R_tr,f_jump,f_wide,Gam_e_max_max)
    else
        DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-1d0)))
        Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    end if
    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-1d0)
    call electron_gm_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c=7.7d8/(1d0+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/2d0)
    if (R_Gamma(1) < 1d0) error stop 'fs_fullhide_coupled requires initial Gamma >= 1'
    beta_Gam=dsqrt(1d0-1d0/R_Gamma(1)**2)
    if (grid_mode /= 0) then
        grid_max=max(grid_top,tail_factor*Gam_e_max_max)
        call electron_loggrid(Num_gam_e,grid_max,gam_e,x_edge)
        call init_edges(f_e*Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                     Num_gam_e,x_edge,dN_x)
        if (thermal_electrons /= 0) then
            call coupled_thermal(R_Gamma(1)*beta_Gam,Para_N_e_ini*(1d0-f_e),dN_x)
        end if
        call dnx_dgamma(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,1))
    else if (thermal_electrons == 0) then
        call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,f_e*Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                          imodelog,gam_e,dN_x,x_edge)
        call dnx_dgamma(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,1))
    else
        call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                          imodeg,gam_e,dN_gam_e(:,1),thermal_electrons=thermal_electrons, &
                                          f_e=f_e,four_v=R_Gamma(1)*beta_Gam)
        dN_x=dN_gam_e(:,1)*gam_e
        call log_edges(Num_gam_e,gam_e,x_edge)
    end if

    if (grid_mode /= 0) then
        d_x=x_edge(2)-x_edge(1)
    else
        d_x=dlog(gam_e(2)/gam_e(1))
    end if

    do I_tobs=2,Num_R
        call prepare_coupled_shell(I_tobs)
        dEL_mean=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/2d0
        call advance_coupled_shell(I_tobs)
        call dnx_dgamma(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
    end do
    call write_finaldiag()

    deallocate(dEl,dEL_mean,dEL_mean_step,cooling_aux,x,dN_x,x_edge,dF1)

contains

    ! joint 输运按网格模式分派热电子 cell 源。
    ! Route the joint thermal cell source by transport-grid mode.
    subroutine coupled_thermal(fourv,count,dst)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: fourv,count
    real(8), intent(inout), dimension(Num_gam_e) :: dst

        if (grid_mode /= 0) then
            call add_thermal_edges(Num_gam_e,x_edge,fourv,count,dst)
        else
            call add_thermal(Num_gam_e,gam_e,fourv,count,dst)
        end if
    end subroutine coupled_thermal

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
        if (grid_mode == 0) dN_x=dN_gam_e(:,I_tobs-1)*gam_e
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
            call electron_injection_prefactor(R_loc-dDR,dDR,dNe,f_e,Gam_e_m_p_step,Q)
            call source_edges(Num_gam_e,x_edge,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
            dF1=dF1+sec_source(:,I_tobs)
            if (thermal_electrons /= 0) then
                thermal_count=Q*(1d0-f_e)/(f_e*Gam_e_m_p_step)
                call coupled_thermal(R_Gamma_loc*beta_Gam,thermal_count,dF1)
            end if
            if (dNe_shell > 0d0) then
                dEL_mean_step=dEL_mean*(dNe/dNe_shell)
            else
                dEL_mean_step=dEL_mean
            end if
            call fullhide_step(Num_gam_e,R_loc,dDR,d_x,dEL_mean_step,dF1,dN_x,x)
            dN_x=x
        end do
    end subroutine advance_coupled_shell
end subroutine coupled_core

end module fullhide_drivers

! 保留原 fullhide public ABI；非 BH 路径仍使用 four-velocity 坐标。
subroutine fs_fullhide_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                                   adaptive_substeps,substep_rtol,substep_min,substep_max,thermal_electrons, &
                                   gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    use fullhide_drivers, only: fullhide_core
    implicit none
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in) :: adaptive_substeps,substep_min,substep_max,thermal_electrons
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu),substep_rtol
    real(8), intent(out) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R)
    real(8), intent(out) :: P_syn(Num_nu,Num_R),Seed_syn(Num_nu,Num_R),V_m(Num_R),V_c(Num_R),V_a(Num_R)

    call fullhide_core(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                       adaptive_substeps,substep_rtol,substep_min,substep_max,thermal_electrons,0,0d0, &
                       gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
end subroutine fs_fullhide_1d

! Forward BH 入口显式接收运动学支撑上边界，并复用同一输运核心。
subroutine fs_fullhide_bh_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                                      adaptive_substeps,substep_rtol,substep_min,substep_max,thermal_electrons,grid_top, &
                                      gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    use fullhide_drivers, only: fullhide_core
    implicit none
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in) :: adaptive_substeps,substep_min,substep_max,thermal_electrons
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu),substep_rtol,grid_top
    real(8), intent(out) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R)
    real(8), intent(out) :: P_syn(Num_nu,Num_R),Seed_syn(Num_nu,Num_R),V_m(Num_R),V_c(Num_R),V_a(Num_R)

    call fullhide_core(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                       adaptive_substeps,substep_rtol,substep_min,substep_max,thermal_electrons,1,grid_top, &
                       gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
end subroutine fs_fullhide_bh_1d

! 保留原 joint-coupling public ABI；旧调用仍使用原 log-gamma 中心网格。
subroutine fs_fullhide_coupled(Boundary,R_Tobs,R_Gamma,R,V_seed,Seed_cooling,sec_source,n,Num_nu,Num_R, &
                                           Num_gam_e,index_Y,index_syn_intger,n_threads,adaptive_substeps, &
                                           substep_rtol,substep_min,substep_max,thermal_electrons, &
                                           gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    use fullhide_drivers, only: coupled_core
    implicit none
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in) :: adaptive_substeps,substep_min,substep_max,thermal_electrons
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: Seed_cooling(Num_nu,Num_R),sec_source(Num_gam_e,Num_R),substep_rtol
    real(8), intent(out) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R)
    real(8), intent(out) :: P_syn(Num_nu,Num_R),Seed_syn(Num_nu,Num_R),V_m(Num_R),V_c(Num_R),V_a(Num_R)

    call coupled_core(Boundary,R_Tobs,R_Gamma,R,V_seed,Seed_cooling,sec_source,n,Num_nu,Num_R,Num_gam_e, &
                      index_Y,index_syn_intger,n_threads,adaptive_substeps,substep_rtol,substep_min,substep_max, &
                      thermal_electrons,0,0d0,gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
end subroutine fs_fullhide_coupled

! Forward BH joint 入口接收与 primary/formal 共用的同一支撑上边界。
subroutine fs_fullhide_coupled_bh(Boundary,R_Tobs,R_Gamma,R,V_seed,Seed_cooling,sec_source,n,Num_nu,Num_R, &
                                              Num_gam_e,index_Y,index_syn_intger,n_threads,adaptive_substeps, &
                                              substep_rtol,substep_min,substep_max,thermal_electrons,grid_top, &
                                              gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    use fullhide_drivers, only: coupled_core
    implicit none
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in) :: adaptive_substeps,substep_min,substep_max,thermal_electrons
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: Seed_cooling(Num_nu,Num_R),sec_source(Num_gam_e,Num_R),substep_rtol,grid_top
    real(8), intent(out) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R)
    real(8), intent(out) :: P_syn(Num_nu,Num_R),Seed_syn(Num_nu,Num_R),V_m(Num_R),V_c(Num_R),V_a(Num_R)

    call coupled_core(Boundary,R_Tobs,R_Gamma,R,V_seed,Seed_cooling,sec_source,n,Num_nu,Num_R,Num_gam_e, &
                      index_Y,index_syn_intger,n_threads,adaptive_substeps,substep_rtol,substep_min,substep_max, &
                      thermal_electrons,1,grid_top,gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
end subroutine fs_fullhide_coupled_bh
