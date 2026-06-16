! 电子1D半拉格朗日格式主驱动：初始化→壳层循环（辐射+冷却+半拉格朗日步进）。
subroutine fs_electron_slc1_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                            gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    use electron_common
    use electron_transport_common, only: electron_semi_lagrangian_step, electron_dnx_to_dndgamma_exp_centers
    use electron_injection_profiles, only: electron_build_source_term_exp_cutoff_edges
    use electron_radiation_kernel, only: get_nu_a, get_syn_selected
    use electron_cooling_kernel, only: get_forward_cooling
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R), &
                            Seed_syn(Num_nu,Num_R),V_m(Num_R),V_c(Num_R),V_a(Num_R)

    real(8), allocatable :: dEl(:),dEL_mean(:),dEL_mean_base(:),dN_x(:),dN_step(:),dF1(:),x_edge(:)
    real(8), allocatable :: gam_e_rad(:),dN_gam_e_rad(:)
    logical :: is_uniform_density,anchor_gamma_m,anchor_gamma_c
    integer :: Num_gam_rad

    allocate(dEl(Num_gam_e),dEL_mean(Num_gam_e-1),dEL_mean_base(Num_gam_e-1),dN_x(Num_gam_e), &
             dN_step(Num_gam_e),dF1(Num_gam_e),x_edge(Num_gam_e+1),gam_e_rad(Num_gam_e),dN_gam_e_rad(Num_gam_e))

    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)

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
    call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                      electron_initial_grid_log_edges,gam_e,dN_x,x_edge)
    call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,1))
    d_x=dlog10(gam_e(2)/gam_e(1))
    is_uniform_density=(A_star <= zero .and. f_jump == one)

    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        dNe_shell=dNe

        beta_Gam=sqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=max(100,min(1000,Int(dDD/dDR)))
        dDR=dDD/L1

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
        dEL_mean_base=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two/dlog(ten)

        do L=1,L1
            R_left=R_loc
            dNe_left=dNe
            R_right=R_left+dDR
            if (.not. is_uniform_density) then
                call dynamics_external_density_profile(A_star,dNe_ISM,R_right,R0,1,R_tr,f_jump,f_wide,dNe_right)
                R_mid=0.5d0*(R_left+R_right)
                call dynamics_external_density_profile(A_star,dNe_ISM,R_mid,R0,1,R_tr,f_jump,f_wide,dNe_mid)
            else
                dNe_right=dNe_left
                R_mid=0.5d0*(R_left+R_right)
                dNe_mid=dNe_left
            end if
            DB_step=0.39d0*dsqrt(Epsilon_b*dNe_mid*(R_Gamma_loc*(R_Gamma_loc-one)))
            Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
            temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
            call electron_gamma_m_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
            Gam_e_m_p_step=(one-p)/(Gam_e_max_step**(one-p)-Gam_e_m_step**(one-p))
            call electron_injection_prefactor(R_mid,dDR,dNe_mid,f_e,Gam_e_m_p_step,Q)
            call electron_build_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
            if (dNe_shell > zero) then
                dEL_mean=dEL_mean_base*(dNe_mid/dNe_shell)
            else
                dEL_mean=dEL_mean_base
            end if
            call electron_semi_lagrangian_step(Num_gam_e,dDR,d_x,dEL_mean+one/R_mid/dlog(ten),dF1,dN_x,dN_step)
            dN_x=dN_step
            R_loc=R_right
            dNe=dNe_right
        end do
        call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
    end do

    deallocate(dEl,dEL_mean,dEL_mean_base,dN_x,dN_step,dF1,x_edge,gam_e_rad,dN_gam_e_rad)
end subroutine fs_electron_slc1_1d
