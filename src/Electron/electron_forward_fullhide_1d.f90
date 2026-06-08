!Calculate the electron distibutions of forward shock.
!New improvement at 11.29.2022
!****************************************************************************************
!******************************* main program *******************************************
!****************************************************************************************
! 电子1D全隐格式主驱动：自适应子步+隐式迎风冷却，支持均匀/非均匀介质和粒子数守恒诊断。
subroutine fs_electron_fullhide_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                                   adaptive_substeps,substep_rtol,substep_min,substep_max,thermal_electrons, &
                                   gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    use electron_common
    use electron_injection_profiles, only: electron_build_source_term_exp_cutoff, electron_add_thermal_source_term
    use electron_radiation_kernel, only: get_nu_a, get_syn_selected
    use electron_cooling_kernel, only: get_forward_cooling
    use electron_transport_common, only: electron_fullhide_step, electron_fullhide_spacetime_sequence
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in) :: adaptive_substeps,substep_min,substep_max,thermal_electrons
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: substep_rtol
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R), &
                            Seed_syn(Num_nu,Num_R), V_m(Num_R), V_c(Num_R), V_a(Num_R)
    
    real(8),allocatable,dimension (:) :: dEl,dEL_mean,x,dN_x, &
                                         dN_full,dN_half,dN_half2,dF1,dEL_mean_base,dEL_mean_step, &
                                         gam_e_rad,dN_gam_e_rad
    real(8),allocatable,dimension (:,:) :: dF_steps,face_coupling
    logical :: is_uniform_density,budget_diag_enabled
    integer :: Num_gam_rad,env_len,env_status
    character(len=32) :: diag_env
    real(8) :: dDR_xi
    real(8) :: n_before_step,n_after_step,inj_step,rel_loss_xi_max
    real(8) :: source_integral,adiabatic_integral,l_count_real,step_sum,step_sq_sum
    real(8) :: radius_sum,radius_sq_sum,source_prefactor
    allocate (dEl(Num_gam_e),dEL_mean(Num_gam_e-1),x(Num_gam_e),dN_x(Num_gam_e), &
              dN_full(Num_gam_e), &
              dN_half(Num_gam_e),dN_half2(Num_gam_e),dF1(Num_gam_e),dEL_mean_base(Num_gam_e-1), &
              dEL_mean_step(Num_gam_e-1),gam_e_rad(Num_gam_e),dN_gam_e_rad(Num_gam_e))
    
    !***********************[Parameter Initial]**********************
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

    !*****************Part 1: given the boundary consition [Using the analytical approximation]*********************
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
    call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                      electron_initial_grid_gamma,gam_e,dN_gam_e(:,1),thermal_electrons=thermal_electrons, &
                                      f_e=f_e,four_v=R_Gamma(1)*beta_Gam)
    !*******************Part 1 is completed [has been checked and there is no bug]**********************************
    !*******************Part 2: To calculate the electron distribution**********************************************
    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
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
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)
        call electron_prepare_radiation_spectrum(Num_gam_e,gam_e,dN_gam_e(:,I_tobs-1), &
                                                 Num_gam_rad,gam_e_rad,dN_gam_e_rad)
        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        call get_nu_a(R_loc,DB,Num_gam_rad,gam_e_rad(1:Num_gam_rad),dN_gam_e_rad(1:Num_gam_rad),temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))

        call get_syn_selected(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                              gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_syn(:,I_tobs),Seed_syn(:,I_tobs))
        
        call get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                 R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                 P_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl)
         
        dEL_mean=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two/dlog(ten)
        dEL_mean_base=dEL_mean
        dDR_xi=dDR
        if (adaptive_substeps == 0) then
            L1=max(100,min(1000,int(dDD/max(dDR,tiny(one)))))
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
                    adiabatic_integral=adiabatic_integral+dDR/((R_loc+dDR*dble(L))*dlog(ten))
                end do
                call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,source_integral,p, &
                                                           dF_steps(:,1))
                face_coupling(:,1)=(dDD*dEL_mean_base+adiabatic_integral)/d_x
                if (budget_diag_enabled) then
                    n_before_step=sum(dN_x)*d_x
                    inj_step=sum(dF_steps)*d_x
                end if
                call electron_fullhide_spacetime_sequence(Num_gam_e,1,face_coupling,dF_steps,dN_x,x)
                if (budget_diag_enabled) then
                    n_after_step=sum(x)*d_x
                    rel_loss_xi_max=max(rel_loss_xi_max, &
                        max(zero,(n_before_step+inj_step-n_after_step)/max(n_before_step+inj_step,tiny(one))))
                    if (I_tobs <= 6) then
                        print '(A,1X,I4,1X,ES12.4,1X,ES12.4,1X,ES12.4)', &
                              'BUDGET1D shell', I_tobs, n_before_step, inj_step, n_after_step
                    end if
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
                    call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
                    if (thermal_electrons /= 0) then
                        call electron_add_thermal_source_term(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam, &
                                                              Q*(one-f_e)/(f_e*Gam_e_m_p_step),dF1)
                    end if
                    if (dNe_shell > zero) then
                        dEL_mean_step=dEL_mean_base*(dNe/dNe_shell)
                    else
                        dEL_mean_step=dEL_mean_base
                    end if
                    if (budget_diag_enabled) then
                        n_before_step=sum(dN_x)*d_x
                        inj_step=dDR*sum(dF1)*d_x
                    end if
                    call electron_fullhide_step(Num_gam_e,R_loc,dDR,d_x,dEL_mean_step,dF1,dN_x,x)
                    if (budget_diag_enabled) then
                        n_after_step=sum(x)*d_x
                        rel_loss_xi_max=max(rel_loss_xi_max, &
                            max(zero,(n_before_step+inj_step-n_after_step)/max(n_before_step+inj_step,tiny(one))))
                        if (I_tobs <= 6 .and. L == L1) then
                            print '(A,1X,I4,1X,ES12.4,1X,ES12.4,1X,ES12.4)', &
                                  'BUDGET1D shell', I_tobs, n_before_step, inj_step, n_after_step
                        end if
                    end if
                    dN_x=x
                end do
            end if
            dN_gam_e(:,I_tobs)=dN_x/gam_e/dlog(ten)
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
                call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m_full,Gam_e_max_full,Q,p,dF1)
                if (thermal_electrons /= 0) then
                    thermal_count=Q*(one-f_e)/(f_e*Gam_e_m_p_full)
                    call electron_add_thermal_source_term(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam,thermal_count,dF1)
                end if
                if (dNe_shell > zero) then
                    dEL_mean_step=dEL_mean_base*(dNe_full/dNe_shell)
                else
                    dEL_mean_step=dEL_mean_base
                end if
                call electron_fullhide_step(Num_gam_e,R_full,dR_try,d_x,dEL_mean_step,dF1,dN_x,dN_full)

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
                call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m_half,Gam_e_max_half,Q,p,dF1)
                if (thermal_electrons /= 0) then
                    thermal_count=Q*(one-f_e)/(f_e*Gam_e_m_p_half)
                    call electron_add_thermal_source_term(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam,thermal_count,dF1)
                end if
                if (dNe_shell > zero) then
                    dEL_mean_step=dEL_mean_base*(dNe_half/dNe_shell)
                else
                    dEL_mean_step=dEL_mean_base
                end if
                call electron_fullhide_step(Num_gam_e,R_half,dR_half,d_x,dEL_mean_step,dF1,dN_x,dN_half)

                call electron_injection_prefactor(R_full,dR_half,dNe_full,f_e,Gam_e_m_p_full,Q)
                call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m_full,Gam_e_max_full,Q,p,dF1)
                if (thermal_electrons /= 0) then
                    thermal_count=Q*(one-f_e)/(f_e*Gam_e_m_p_full)
                    call electron_add_thermal_source_term(Num_gam_e,gam_e,R_Gamma_loc*beta_Gam,thermal_count,dF1)
                end if
                if (dNe_shell > zero) then
                    dEL_mean_step=dEL_mean_base*(dNe_full/dNe_shell)
                else
                    dEL_mean_step=dEL_mean_base
                end if
                call electron_fullhide_step(Num_gam_e,R_full,dR_half,d_x,dEL_mean_step,dF1,dN_half,dN_half2)

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
            dN_gam_e(:,I_tobs)=dN_x/gam_e/dlog(ten)
        end if
    end do

    deallocate (dEl,dEL_mean,x,dN_x,dN_full,dN_half,dN_half2,dF1, &
                dEL_mean_base,dEL_mean_step,gam_e_rad,dN_gam_e_rad)
    if (budget_diag_enabled) then
        print '(A,1X,ES12.4)', 'BUDGET1D max_rel_loss', rel_loss_xi_max
    end if

    return
end subroutine fs_electron_fullhide_1d
