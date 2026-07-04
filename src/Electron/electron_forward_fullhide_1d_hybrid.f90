! HZ hybrid 电子分布路径：热/非热源项共用 fullhide 输运，辐射诊断仍走 selected synchrotron。
! HZ hybrid electron path: thermal and non-thermal sources share fullhide transport, with selected synchrotron diagnostics.
subroutine fs_fullhide_hz(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                      index_Y,index_syn_intger,n_threads,adaptive_substeps,substep_rtol, &
                                      substep_min,substep_max,thermal_electrons,gam_e,dN_gam_e,P_syn, &
                                      Seed_syn,V_m,V_c,V_a)
    use constants
    use dynamics_density_profile, only: density_profile
    use electron_common
    use electron_injection_profiles, only: init_coord, source_coord
    use electron_coord_common, only: build_fourvel_grid, dxg_dcoord, coord_fourvel, fourvel_scale
    use electron_radiation_kernel, only: syn_state, nua_fromtau
    use electron_cooling_kernel, only: get_forward_cooling
    use electron_shell_transport, only: shell_coord_step, coord_to_dgamma
    use hybrid_spectrum, only: hybrid_coord
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

    real(8),allocatable,dimension (:) :: dEl,dEl_step,x,dN_x,x_edge,coord_edge,dxdy_grid, &
                                         dN_full,dN_half,dN_half2,dF1
    logical :: is_uniform_density,budget_diag_enabled
    integer :: env_len,env_status
    character(len=32) :: diag_env
    real(8) :: dDR_xi,ln10,coord_scale,dg_gamma_scale
    real(8) :: n_before_step,n_after_step,inj_step,rel_loss_xi_max
    real(8), dimension(Num_nu) :: P_emit_tmp,Tau_syn_tmp
    allocate (dEl(Num_gam_e),dEl_step(Num_gam_e),x(Num_gam_e),dN_x(Num_gam_e), &
              dN_full(Num_gam_e), &
              x_edge(Num_gam_e+1),coord_edge(Num_gam_e+1),dxdy_grid(Num_gam_e), &
              dN_half(Num_gam_e),dN_half2(Num_gam_e),dF1(Num_gam_e))

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
    if (R_Gamma(1) < 1d0) error stop 'fs_fullhide_1d requires initial Gamma >= 1'
    beta_Gam=dsqrt(1d0-1d0/R_Gamma(1)**2)
    ln10=dlog(1d1)
    call init_fourvel_grid(Gam_e_max_max)
    if (thermal_electrons == 0) then
        call init_coord(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,coord_edge,coord_scale,dN_x)
    else
        call hybrid_coord(Num_gam_e,coord_edge,coord_scale,p,Gam_e_m,Gam_e_max,f_e,dN_x)
        dN_x=dN_x*Para_N_e_ini
    end if
    call coord_to_dgamma(Num_gam_e,coord_edge,coord_scale,gam_e,dN_x,dN_gam_e(:,1))
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
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*ln10*dxdy_grid
        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))

        call syn_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_emit_tmp,P_syn(:,I_tobs), &
                                    Seed_syn(:,I_tobs),Tau_syn_tmp)
        call nua_fromtau(Num_nu,V_seed,Tau_syn_tmp,temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))

        call get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                 R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                 P_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl)

        dDR_xi=dDR
        if (adaptive_substeps == 0) then
            L1=max(100,min(1000,int(dDD/max(dDR,tiny(1d0)))))
            dDR=dDD/dble(L1)

            do L=1,L1
                    R_loc=R_loc+dDR

                    if (is_uniform_density .and. thermal_electrons == 0) then
                        call build_hybrid_source(R_loc,dDR,dNe,Gam_e_m,Gam_e_max,Gam_e_m_p)
                    else
                        call density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
                        DB_step=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-1d0)))
                        Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
                        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
                        call electron_gm_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
                        Gam_e_m_p_step=(1d0-p)/(Gam_e_max_step**(1d0-p)-Gam_e_m_step**(1d0-p))
                        call build_hybrid_source(R_loc,dDR,dNe,Gam_e_m_step,Gam_e_max_step,Gam_e_m_p_step)

                        if (dNe_shell > 0d0) then
                            dEl_step=dEl*(dNe/dNe_shell)
                        else
                            dEl_step=dEl
                        end if
                    end if
                    if (is_uniform_density .and. thermal_electrons == 0) dEl_step=dEl
                    if (budget_diag_enabled) then
                        n_before_step=sum(dN_x*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
                        inj_step=dDR*sum(dF1*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
                    end if

                    call shell_coord_step(Num_gam_e,dDR,coord_edge,coord_scale,dEl_step,1d0/R_loc,dF1,dN_x,x)
                    if (budget_diag_enabled) then
                        n_after_step=sum(x*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e)))
                        rel_loss_xi_max=max(rel_loss_xi_max, &
                            max(0d0,(n_before_step+inj_step-n_after_step)/max(n_before_step+inj_step,tiny(1d0))))
                        if (I_tobs <= 6 .and. L == L1) then
                            print '(A,1X,I4,1X,ES12.4,1X,ES12.4,1X,ES12.4)', &
                                  'BUDGET1D shell', I_tobs, n_before_step, inj_step, n_after_step
                        end if
                    end if
                    dN_x=x

                    if (L1 == L) then
                        call coord_to_dgamma(Num_gam_e,coord_edge,coord_scale,gam_e,dN_x,dN_gam_e(:,I_tobs))
                    end if
            end do
        else
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
                call build_hybrid_count(Q,Gam_e_m_full,Gam_e_max_full,Gam_e_m_p_full)

                if (dNe_shell > 0d0) then
                    dEl_step=dEl*(dNe_full/dNe_shell)
                else
                    dEl_step=dEl
                end if
                call shell_coord_step(Num_gam_e,dR_try,coord_edge,coord_scale,dEl_step,1d0/R_full,dF1,dN_x,dN_full)

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
                call build_hybrid_count(Q,Gam_e_m_half,Gam_e_max_half,Gam_e_m_p_half)

                if (dNe_shell > 0d0) then
                    dEl_step=dEl*(dNe_half/dNe_shell)
                else
                    dEl_step=dEl
                end if
                call shell_coord_step(Num_gam_e,dR_half,coord_edge,coord_scale,dEl_step,1d0/R_half,dF1,dN_x,dN_half)

                call electron_injection_prefactor(R_full,dR_half,dNe_full,f_e,Gam_e_m_p_full,Q)
                call build_hybrid_count(Q,Gam_e_m_full,Gam_e_max_full,Gam_e_m_p_full)

                if (dNe_shell > 0d0) then
                    dEl_step=dEl*(dNe_full/dNe_shell)
                else
                    dEl_step=dEl
                end if
                call shell_coord_step(Num_gam_e,dR_half,coord_edge,coord_scale,dEl_step,1d0/R_full,dF1,dN_half,dN_half2)

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
            call coord_to_dgamma(Num_gam_e,coord_edge,coord_scale,gam_e,dN_x,dN_gam_e(:,I_tobs))
        end if
    end do
    call write_finaldiag()

    deallocate (dEl,dEl_step,x,dN_x,x_edge,coord_edge,dxdy_grid,dN_full,dN_half,dN_half2,dF1)
    if (budget_diag_enabled) then
        print '(A,1X,ES12.4)', 'BUDGET1D max_rel_loss', rel_loss_xi_max
    end if

    return

contains

    subroutine init_fourvel_grid(imodeg_max)
    implicit real(8)(A-H,O-Z)
    integer :: I_grid
    real(8), intent(in) :: imodeg_max

        dg_gamma_scale=fourvel_scale
        coord_scale=dg_gamma_scale*dg_gamma_scale-1d0
        call build_fourvel_grid(Num_gam_e,1d0,tail_factor*imodeg_max,dg_gamma_scale,gam_e,coord_edge,x_edge)
        do I_grid=1,Num_gam_e
            dxdy_grid(I_grid)=dxg_dcoord(coord_fourvel,coord_scale, &
                                         0.5d0*(coord_edge(I_grid)+coord_edge(I_grid+1)))
        end do
    end subroutine init_fourvel_grid

    ! 最后一个输出点没有后续推进，只刷新与最终电子谱一致的辐射诊断。
    ! The final output point has no following advance, so refresh diagnostics from the final electron spectrum.
    subroutine write_finaldiag()
    implicit none

        R_loc=R(Num_R)
        R_Gamma_loc=R_Gamma(Num_R)
        if (R_Gamma_loc < 1d0) error stop 'fs_fullhide_hz requires Gamma >= 1'
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
                                    gam_e,dN_gam_e(:,Num_R),V_seed,P_emit_tmp, &
                                    P_syn(:,Num_R),Seed_syn(:,Num_R),Tau_syn_tmp)
        call nua_fromtau(Num_nu,V_seed,Tau_syn_tmp,temp)
        V_a(Num_R)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
    end subroutine write_finaldiag

    ! 从半径壳层几何直接构建本子步注入源。
    ! Build the substep injection source directly from shell geometry.
    subroutine build_hybrid_source(rsrc,drsrc,nsrc,gm_src,gmax_src, &
                                               gmp_src)
    implicit none
    real(8), intent(in) :: rsrc,drsrc,nsrc,gm_src,gmax_src,gmp_src

        if (thermal_electrons /= 0) then
            Q = 4d0/3d0*pi*(3d0*rsrc**2+drsrc*(3d0*rsrc+drsrc))*nsrc
            call hybrid_coord(Num_gam_e,coord_edge,coord_scale,p,gm_src,gmax_src,f_e,dF1)
            dF1 = dF1*Q
        else
            Q = 4d0/3d0*pi*(3d0*rsrc**2+drsrc*(3d0*rsrc+drsrc))*nsrc*f_e*gmp_src
            call source_coord(Num_gam_e,coord_edge,coord_scale,gm_src,gmax_src,Q,p,dF1)
        end if
    end subroutine build_hybrid_source

    ! 自适应半步已给出粒子数时，只按 hybrid/power-law 谱形分配源项。
    ! When adaptive half-steps already provide the particle count, only distribute it by hybrid/power-law shape.
    subroutine build_hybrid_count(nsource,gm_src,gmax_src,gmp_src)
    implicit none
    real(8), intent(in) :: nsource,gm_src,gmax_src,gmp_src

        if (thermal_electrons /= 0) then
            Q = nsource/(f_e*gmp_src)
            call hybrid_coord(Num_gam_e,coord_edge,coord_scale,p,gm_src,gmax_src,f_e,dF1)
            dF1 = dF1*Q
        else
            call source_coord(Num_gam_e,coord_edge,coord_scale,gm_src,gmax_src,nsource,p,dF1)
        end if
    end subroutine build_hybrid_count
end subroutine fs_fullhide_hz
