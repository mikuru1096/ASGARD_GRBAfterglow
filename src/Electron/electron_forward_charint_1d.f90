! 电子 1D 特征线积分主驱动：初始化电子谱，再逐壳层做辐射、冷却和特征线更新。
! 1D electron characteristic integrator: initialize the spectrum, then update radiation, cooling, and characteristics shell by shell.
subroutine fs_charint_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                                adaptive_substeps,substep_rtol,substep_min,substep_max,gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    use constants
    use dynamics_density_profile, only: density_profile
    use electron_common
    use electron_transport_common, only: cfl_relax, rtol_relax, &
        cooling_affine, cooling_piecewise, char_update, &
        dnx_dgamma
    use electron_injection_profiles, only: source_edges
    use electron_radiation_kernel, only: syn_state, nua_fromtau
    use electron_cooling_kernel, only: forward_cooling
    implicit none
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads,adaptive_substeps,substep_min,substep_max
    integer :: I_tobs,L,L1
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in) :: substep_rtol
    real(8), intent(out), dimension(Num_gam_e,Num_R) :: dN_gam_e
    real(8), intent(out), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_nu,Num_R) :: P_syn,Seed_syn
    real(8), intent(out), dimension(Num_R) :: V_m,V_c,V_a
    real(8) :: Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star,E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0
    real(8) :: dNe,Para_N_e_ini,DB,DB_min,DB_step,beta_Gam,f_r,cooling_scale,a_rad,b_ad,Q,temp
    real(8) :: Gam_e_max,Gam_e_max_max,Gam_e_max_step,Gam_e_m,Gam_e_m_p,Gam_e_m_step,Gam_e_m_p_step,Gam_e_c,temp_gam
    real(8) :: R_loc,R_left,R_right,R_mid,R_Gamma_loc,dDR,dDD,dR_min,dR_max,dR_try,dR_left,dR_cap,step_limit,step_error
    real(8) :: dNe_shell,dNe_mid,dNe_right
    real(8), allocatable, dimension(:) :: dEl_base,dEl_step,dN_x,dN_step,dF1,dF1_shape,x_edge
    real(8), dimension(Num_nu) :: P_emit_tmp,Tau_syn_tmp
    logical :: is_uniform_density

    allocate(dEl_base(Num_gam_e),dEl_step(Num_gam_e),dN_x(Num_gam_e),dN_step(Num_gam_e), &
             dF1(Num_gam_e),dF1_shape(Num_gam_e),x_edge(Num_gam_e+1))

    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)

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
    call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                      imodelog,gam_e,dN_x,x_edge)
    call dnx_dgamma(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,1))
    is_uniform_density=(A_star <= 0d0 .and. f_jump == 1d0)

    do I_tobs=2,Num_R
        call prepare_characteristic_shell(I_tobs)

        if (adaptive_substeps == 0) then
            if (is_uniform_density) then
                dNe_mid=dNe_shell
                DB_step=0.39d0*dsqrt(Epsilon_b*dNe_mid*(R_Gamma_loc*(R_Gamma_loc-1d0)))
                Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
                temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
                call electron_gm_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
                Gam_e_m_p_step=(1d0-p)/(Gam_e_max_step**(1d0-p)-Gam_e_m_step**(1d0-p))
                call source_edges(Num_gam_e,x_edge,Gam_e_m_step,Gam_e_max_step,1d0,p,dF1_shape)
                if (index_Y == 0) then
                    cooling_scale=1d0/(beta_Gam*R_Gamma_loc)
                    a_rad=1.35d-19*DB_step**2*cooling_scale/pi
                end if

                R_mid=R_loc+0.5d0*dDR
                do L=1,L1
                    call electron_injection_prefactor(R_mid,dDR,dNe_mid,f_e,Gam_e_m_p_step,Q)

                    if (index_Y == 0) then
                        b_ad=1d0/R_mid
                        call char_update(Num_gam_e,dDR,x_edge,cooling_affine, &
                            a_rad,b_ad,gam_e,dEl_step,R_mid,Q,dF1_shape,dN_x,dN_step)
                    else
                        dEl_step=dEl_base
                        call char_update(Num_gam_e,dDR,x_edge,cooling_piecewise, &
                            0d0,0d0,gam_e,dEl_step,R_mid,Q,dF1_shape,dN_x,dN_step)
                    end if

                    dN_x=dN_step
                    R_mid=R_mid+dDR
                end do
                R_loc=R(I_tobs)
                dNe=dNe_shell
            else
                do L=1,L1
                    R_left=R_loc
                    R_right=R_left+dDR
                    call density_profile(A_star,dNe_ISM,R_right,R0,1,R_tr,f_jump,f_wide,dNe_right)
                    R_mid=0.5d0*(R_left+R_right)
                    call density_profile(A_star,dNe_ISM,R_mid,R0,1,R_tr,f_jump,f_wide,dNe_mid)
                    DB_step=0.39d0*dsqrt(Epsilon_b*dNe_mid*(R_Gamma_loc*(R_Gamma_loc-1d0)))
                    Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
                    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
                    call electron_gm_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
                    Gam_e_m_p_step=(1d0-p)/(Gam_e_max_step**(1d0-p)-Gam_e_m_step**(1d0-p))
                    call electron_injection_prefactor(R_mid,dDR,dNe_mid,f_e,Gam_e_m_p_step,Q)
                    call source_edges(Num_gam_e,x_edge,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)

                    if (index_Y == 0) then
                        cooling_scale=1d0/(beta_Gam*R_Gamma_loc)
                        a_rad=1.35d-19*DB_step**2*cooling_scale/pi
                        b_ad=1d0/R_mid
                        call char_update(Num_gam_e,dDR,x_edge,cooling_affine, &
                            a_rad,b_ad,gam_e,dEl_step,R_mid,1d0,dF1,dN_x,dN_step)
                    else
                        if (dNe_shell > 0d0) then
                            dEl_step=dEl_base*(dNe_mid/dNe_shell)
                        else
                            dEl_step=dEl_base
                        end if
                        call char_update(Num_gam_e,dDR,x_edge,cooling_piecewise, &
                            0d0,0d0,gam_e,dEl_step,R_mid,1d0,dF1,dN_x,dN_step)
                    end if

                    dN_x=dN_step
                    R_loc=R_right
                    dNe=dNe_right
                end do
            end if
        else
            dR_min=dDD/max(substep_max,1)
            dR_max=dDD/max(substep_min,1)
            dR_try=min(dDR,dR_max)
            dR_left=dDD

            if (is_uniform_density) then
                dNe_mid=dNe_shell
                DB_step=0.39d0*dsqrt(Epsilon_b*dNe_mid*(R_Gamma_loc*(R_Gamma_loc-1d0)))
                Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
                temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
                call electron_gm_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
                Gam_e_m_p_step=(1d0-p)/(Gam_e_max_step**(1d0-p)-Gam_e_m_step**(1d0-p))
                call source_edges(Num_gam_e,x_edge,Gam_e_m_step,Gam_e_max_step,1d0,p,dF1_shape)
                if (index_Y == 0) then
                    cooling_scale=1d0/(beta_Gam*R_Gamma_loc)
                    a_rad=1.35d-19*DB_step**2*cooling_scale/pi
                end if
            end if

            do while (dR_left > 0d0)
                dR_try=min(dR_try,dR_left)
                if (dR_left < 1.5d0*dR_try) dR_try=dR_left
                step_limit=rtol_relax*substep_rtol
                if (step_limit > 0d0) then
                    dR_cap=step_limit*max(R_loc,1d-30)/max(1d0-0.5d0*step_limit,0.5d0)
                    if (dR_cap > dR_min) dR_try=min(dR_try,dR_cap)
                end if
                if (is_uniform_density) then
                    R_mid=R_loc+0.5d0*dR_try
                    step_error=dR_try/max(R_mid,1d-30)
                    if (step_error > rtol_relax*substep_rtol .and. dR_try > dR_min) then
                        dR_try=max(0.5d0*dR_try,dR_min)
                        cycle
                    end if

                    call electron_injection_prefactor(R_mid,dR_try,dNe_shell,f_e,Gam_e_m_p_step,Q)
                    if (index_Y == 0) then
                        b_ad=1d0/R_mid
                        call char_update(Num_gam_e,dR_try,x_edge,cooling_affine, &
                            a_rad,b_ad,gam_e,dEl_step,R_mid,Q,dF1_shape,dN_x,dN_step)
                    else
                        dEl_step=dEl_base
                        call char_update(Num_gam_e,dR_try,x_edge,cooling_piecewise, &
                            0d0,0d0,gam_e,dEl_step,R_mid,Q,dF1_shape,dN_x,dN_step)
                    end if
                else
                    R_right=R_loc+dR_try
                    call density_profile(A_star,dNe_ISM,R_right,R0,1,R_tr,f_jump,f_wide,dNe_right)
                    R_mid=R_loc+0.5d0*dR_try
                    call density_profile(A_star,dNe_ISM,R_mid,R0,1,R_tr,f_jump,f_wide,dNe_mid)
                    step_error=max(dR_try/max(R_mid,1d-30),abs(dNe_right-dNe)/max(abs(dNe),1d-99))
                    if (step_error > rtol_relax*substep_rtol .and. dR_try > dR_min) then
                        dR_try=max(0.9d0*step_limit*dR_try/max(step_error,1d-30),dR_min)
                        cycle
                    end if

                    DB_step=0.39d0*dsqrt(Epsilon_b*dNe_mid*(R_Gamma_loc*(R_Gamma_loc-1d0)))
                    Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
                    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
                    call electron_gm_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
                    Gam_e_m_p_step=(1d0-p)/(Gam_e_max_step**(1d0-p)-Gam_e_m_step**(1d0-p))
                    call electron_injection_prefactor(R_mid,dR_try,dNe_mid,f_e,Gam_e_m_p_step,Q)
                    call source_edges(Num_gam_e,x_edge,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
                    if (index_Y == 0) then
                        cooling_scale=1d0/(beta_Gam*R_Gamma_loc)
                        a_rad=1.35d-19*DB_step**2*cooling_scale/pi
                        b_ad=1d0/R_mid
                        call char_update(Num_gam_e,dR_try,x_edge,cooling_affine, &
                            a_rad,b_ad,gam_e,dEl_step,R_mid,1d0,dF1,dN_x,dN_step)
                    else
                        if (dNe_shell > 0d0) then
                            dEl_step=dEl_base*(dNe_mid/dNe_shell)
                        else
                            dEl_step=dEl_base
                        end if
                        call char_update(Num_gam_e,dR_try,x_edge,cooling_piecewise, &
                            0d0,0d0,gam_e,dEl_step,R_mid,1d0,dF1,dN_x,dN_step)
                    end if
                end if

                if (step_error <= rtol_relax*substep_rtol .or. dR_try <= dR_min) then
                    dN_x=dN_step
                    R_loc=R_loc+dR_try
                    if (is_uniform_density) then
                        dNe=dNe_shell
                    else
                        dNe=dNe_right
                    end if
                    dR_left=dR_left-dR_try
                    if (dR_left > 0d0) then
                        if (step_error <= 0.25d0*rtol_relax*substep_rtol) then
                            dR_try=min(2d0*dR_try,dR_max)
                        end if
                    end if
                else
                    dR_try=max(0.5d0*dR_try,dR_min)
                end if
            end do
        end if
        call dnx_dgamma(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
    end do

    call write_final_diag()

    deallocate(dEl_base,dEl_step,dN_x,dN_step,dF1,dF1_shape,x_edge)

contains

    ! 准备当前输出壳层的磁场、注入能标、辐射和基准冷却率。
    ! Prepare magnetic field, injection scales, radiation, and baseline cooling for the current output shell.
    subroutine prepare_characteristic_shell(I_tobs)
        implicit none
        integer, intent(in) :: I_tobs

        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/2d0
        call density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-1d0)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
        call electron_gm_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(1d0-p)/(Gam_e_max**(1d0-p)-Gam_e_m**(1d0-p))
        Gam_e_c=7.7d8*(1d0+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        dNe_shell=dNe

        beta_Gam=sqrt(1d0-1d0/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=cfl_relax*0.1d0/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=max(4,min(2048,ceiling(dDD/max(dDR,1d-30))))
        dDR=dDD/L1

        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))

        call syn_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_emit_tmp,P_syn(:,I_tobs), &
                                    Seed_syn(:,I_tobs),Tau_syn_tmp)
        call nua_fromtau(Num_nu,V_seed,Tau_syn_tmp,temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        if (index_Y == 0) then
            dEl_base=0d0
        else
            call forward_cooling(1,index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                 R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,1,n_threads,gam_e,V_seed, &
                                 P_syn(:,I_tobs),Seed_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl_step,dEl_base)
        end if
    end subroutine prepare_characteristic_shell

    ! 最后一个输出点没有后续推进，只刷新诊断辐射量。
    ! The final output point has no following advance, so only refresh diagnostic radiation quantities.
    subroutine write_final_diag()
        implicit none

        R_loc=R(Num_R)
        R_Gamma_loc=R_Gamma(Num_R)
        call density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-1d0)))
        beta_Gam=sqrt(1d0-1d0/R_Gamma_loc**2)
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-1d0)
        call electron_gm_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_c=7.7d8*(1d0+z)/R_Gamma_loc/DB**2/R_Tobs(Num_R)
        V_m(Num_R)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        V_c(Num_R)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
        call syn_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,Num_R),V_seed,P_emit_tmp,P_syn(:,Num_R), &
                                    Seed_syn(:,Num_R),Tau_syn_tmp)
        call nua_fromtau(Num_nu,V_seed,Tau_syn_tmp,temp)
        V_a(Num_R)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(1d0+z))
    end subroutine write_final_diag

end subroutine fs_charint_1d
