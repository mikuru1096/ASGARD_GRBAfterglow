!Calculate the electron distributions of forward shock.
!Modified to use Implicit Three-Level Method (second-order time accuracy)
!Modified at 2024 by AI assistant based on original code from 11.29.2022
!****************************************************************************************
!******************************* main program *******************************************
!****************************************************************************************
! 电子1D三层隐式格式主驱动：二阶时间精度BDF2（启动用单步），迎风+隐式冷却输运。
subroutine fs_electron_t2g1_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                            gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    use electron_common
    use electron_injection_profiles, only: electron_build_source_term_exp_cutoff_edges, electron_profile_log_cell_edges
    use electron_radiation_kernel, only: get_nu_a, get_syn_selected_state
    use electron_cooling_kernel, only: get_forward_cooling
    use electron_transport_common, only: electron_prepare_implicit_coeffs_common, electron_backward_sweep_common, &
                                         electron_dnx_to_dndgamma_exp_centers
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer :: I_tobs,L,L1,Num_gam_rad
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R), &
                            Seed_syn(Num_nu,Num_R), V_m(Num_R), V_c(Num_R), V_a(Num_R)
    real(8) :: Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star,E_iso,T_log10_duration,f_e
    real(8) :: R_tr,f_jump,f_wide,R0,dNe,Para_N_e_ini,DB,Gam_e_max,DB_min,Gam_e_max_max
    real(8) :: temp_gam,Gam_e_m,Gam_e_c,d_x,R_loc,R_Gamma_loc,Gam_e_m_p,dNe_shell
    real(8) :: beta_Gam,f_r,dDR,dDD,CFL,temp,DB_step,Gam_e_max_step,Gam_e_m_step,Gam_e_m_p_step,Q
    real(8) :: P_emit_tmp(Num_nu),Tau_syn_tmp(Num_nu)
    
    real(8),allocatable,dimension (:) :: dEl,dEL_mean,dEL_mean_base,principal,x,dF1,up,dot_gam_e_SSA, &
                                         dN_x,dN_x_prev,x_edge,temp1,temp2,temp3,temp4,para_maxwell,Compton,Compton1,dot_gam_e, &
                                         gam_e_rad,dN_gam_e_rad
    allocate (dEl(Num_gam_e),dEL_mean(Num_gam_e-1),dEL_mean_base(Num_gam_e-1),principal(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e), &
              up(Num_gam_e-1),dN_x(Num_gam_e),dN_x_prev(Num_gam_e),x_edge(Num_gam_e+1),temp1(Num_gam_e-1), &
              temp2(Num_gam_e),para_maxwell(Num_gam_e),temp3(Num_gam_e-1),temp4(Num_gam_e-1), &
              Compton(Num_gam_e),dot_gam_e(Num_gam_e),dot_gam_e_SSA(Num_gam_e), &
              Compton1(Num_gam_e),gam_e_rad(Num_gam_e),dN_gam_e_rad(Num_gam_e))
    
    !***********************[Parameter Initial]**********************
    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    
    P_syn=zero
    Seed_syn=zero
    V_m=zero
    V_c=zero
    V_a=zero

    !*****************Part 1: given the boundary condition [Using the analytical approximation]*********************
    call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)

    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-one)
    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                      electron_initial_grid_log_edges,gam_e,dN_x,x_edge)
    call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,1))
    !*******************Part 1 is completed [has been checked and there is no bug]**********************************
    !*******************Part 2: To calculate the electron distribution**********************************************
    dN_x_prev = dN_x
    d_x=dlog10(gam_e(2)/gam_e(1))

    ! For the first few steps, we need to use single-step methods
    ! We'll use a startup procedure
    do I_tobs=2,Num_R
        call prepare_t2g1_shell(I_tobs)
        ! For the first step of each I_tobs, we use the solution from previous step
        if (I_tobs == 2) then
            dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
            dN_x_prev = dN_x
        else
            dN_x=dN_x  ! Keep the solution from previous I_tobs
        end if
        call write_t2g1_radiation_and_cooling(I_tobs)

        ! Main loop for sub-steps
        do L=1,L1
            call advance_t2g1_substep(I_tobs,L)
        end do
    end do

    deallocate (dEl,dEL_mean,dEL_mean_base,principal,x,dF1,up,dN_x,dN_x_prev,x_edge,temp1,temp2, &
                para_maxwell,temp3,temp4,Compton,Compton1,gam_e_rad,dN_gam_e_rad)

    return

contains

    subroutine prepare_t2g1_shell(I_tobs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs

        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,0,R_tr,f_jump,f_wide,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        dNe_shell=dNe

        beta_Gam=sqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1/(f_r*Gam_e_max+1.333/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=max(100,min(1000,Int(dDD/dDR)))
        dDR=dDD/L1
        CFL=dDR/d_x
    end subroutine prepare_t2g1_shell

    subroutine write_t2g1_radiation_and_cooling(I_tobs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs

        call electron_prepare_radiation_spectrum(Num_gam_e,gam_e,dN_gam_e(:,I_tobs-1), &
                                                 Num_gam_rad,gam_e_rad,dN_gam_e_rad)

        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        call get_nu_a(R_loc,DB,Num_gam_rad,gam_e_rad(1:Num_gam_rad),dN_gam_e_rad(1:Num_gam_rad),temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))

        call get_syn_selected_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_emit_tmp, &
                                    P_syn(:,I_tobs),Seed_syn(:,I_tobs),Tau_syn_tmp)

        call get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                 R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                 P_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl)

        dEL_mean=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two/dlog(ten)
        dEL_mean_base=dEL_mean
    end subroutine write_t2g1_radiation_and_cooling

    subroutine advance_t2g1_substep(I_tobs,L)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs,L

        R_loc=R_loc+dDR

        call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
        DB_step=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
        Gam_e_m_p_step=(one-p)/(Gam_e_max_step**(one-p)-Gam_e_m_step**(one-p))

        call electron_injection_prefactor(R_loc,dDR,dNe,f_e,Gam_e_m_p_step,Q)
        call electron_build_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
        if (dNe_shell > zero) then
            dEL_mean=dEL_mean_base*(dNe/dNe_shell)
        else
            dEL_mean=dEL_mean_base
        end if

        temp3=dEL_mean+one/R_loc/dlog(ten)
        up=-CFL*temp3

        if (I_tobs == 2 .and. L <= 2) then
            call electron_prepare_implicit_coeffs_common(Num_gam_e,one,up,principal,temp1)
            temp2 = (dN_x + dDR * dF1) / principal
        else
            call electron_prepare_implicit_coeffs_common(Num_gam_e,1.5d0,up,principal,temp1)
            temp2 = ( (2d0)*dN_x - 0.5d0*dN_x_prev + dF1 * dDR ) / principal
        end if

        call electron_backward_sweep_common(Num_gam_e,temp1,temp2,x)

        dN_x_prev = dN_x
        dN_x = x

        if (L1 == L) then
            call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
        end if
    end subroutine advance_t2g1_substep
end subroutine fs_electron_t2g1_1d
