subroutine fs_electron_charint_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                                adaptive_substeps,substep_rtol,substep_min,substep_max,gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    !$ use omp_lib
    use constants
    use electron_common
    use electron_injection_profiles, only: electron_profile_log_cell_edges, electron_initial_powerlaw_exp_cutoff_edges, &
                                           electron_build_source_term_exp_cutoff
    use get_Y
    use electron_forward_kernel, only: electron_forward_cooling_step
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads,adaptive_substeps,substep_min,substep_max
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: substep_rtol
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R), &
                            Seed_syn(Num_nu,Num_R),V_m(Num_R),V_c(Num_R),V_a(Num_R)

    real(8), allocatable :: dEl_base(:),dEl_step(:),dN_x(:),dN_step(:),dF1(:),dF1_shape(:),x_edge(:)
    real(8), allocatable :: q_left_source_shape(:),q_right_source_shape(:),prefix_source_lin_shape(:)
    real(8), allocatable :: gam_e_rad(:),dN_gam_e_rad(:)
    logical :: is_uniform_density
    integer :: Num_gam_rad

    allocate(dEl_base(Num_gam_e),dEl_step(Num_gam_e),dN_x(Num_gam_e),dN_step(Num_gam_e), &
             dF1(Num_gam_e),dF1_shape(Num_gam_e),x_edge(Num_gam_e+1),gam_e_rad(Num_gam_e),dN_gam_e_rad(Num_gam_e), &
             q_left_source_shape(Num_gam_e),q_right_source_shape(Num_gam_e),prefix_source_lin_shape(0:Num_gam_e))

    Eta_0=Boundary(1)
    R_ini=Boundary(4)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    f_e=Boundary(16)
    R_tr=Boundary(21)
    f_jump=Boundary(22)
    f_wide=Boundary(23)
    R0=Boundary(n)

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
    call electron_build_gamma_grid(Num_gam_e,Gam_e_max_max,gam_e)
    call electron_profile_log_cell_edges(Num_gam_e,gam_e,x_edge)
    call electron_initial_powerlaw_exp_cutoff_edges(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,x_edge,dN_x)
    dN_gam_e(:,1)=dN_x/gam_e/dlog(ten)
    is_uniform_density=(A_star <= zero .and. f_jump == one)

    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        call electron_external_density(A_star,dNe_ISM,R_loc,R0,R_tr,f_jump,f_wide,1,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        dNe_shell=dNe

        beta_Gam=sqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=charint_cfl_relax*0.1d0/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=max(4,min(2048,ceiling(dDD/max(dDR,1d-30))))
        dDR=dDD/L1

        call electron_prepare_radiation_spectrum(Num_gam_e,gam_e,dN_gam_e(:,I_tobs-1), &
                                                 Num_gam_rad,gam_e_rad,dN_gam_e_rad)
        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        call get_nu_a(R_loc,DB,Num_gam_rad,gam_e_rad(1:Num_gam_rad),dN_gam_e_rad(1:Num_gam_rad),temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))

        call get_syn_selected(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                              gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_syn(:,I_tobs),Seed_syn(:,I_tobs))
        if (index_Y == 0) then
            dEl_base=zero
        else
            call electron_forward_cooling_step(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                     R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                     P_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl_base)
        end if

        if (adaptive_substeps == 0) then
            if (is_uniform_density) then
                dNe_mid=dNe_shell
                DB_step=0.39d0*dsqrt(Epsilon_b*dNe_mid*(R_Gamma_loc*(R_Gamma_loc-one)))
                Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
                temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
                call electron_gamma_m_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
                Gam_e_m_p_step=(one-p)/(Gam_e_max_step**(one-p)-Gam_e_m_step**(one-p))
                call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m_step,Gam_e_max_step,one,p,dF1_shape)
                call electron_prepare_conservative_remap_nonuniform(Num_gam_e,x_edge,dF1_shape, &
                                                                   q_left_source_shape,q_right_source_shape,prefix_source_lin_shape)

                if (index_Y == 0) then
                    cooling_scale=one/(beta_Gam*R_Gamma_loc)
                    a_rad=1.35d-19*DB_step**2*cooling_scale/pi
                end if

                R_mid=R_loc+0.5d0*dDR
                do L=1,L1
                    call electron_injection_prefactor(R_mid,dDR,dNe_mid,f_e,Gam_e_m_p_step,Q)

                    if (index_Y == 0) then
                        b_ad=one/R_mid
                        call electron_characteristic_step_affine_u_prepared_source( &
                            Num_gam_e,dDR,x_edge,a_rad,b_ad,Q,q_left_source_shape,q_right_source_shape, &
                            prefix_source_lin_shape,dF1_shape,dN_x,dN_step)
                    else
                        dEl_step=dEl_base
                        call electron_characteristic_step_piecewise_u_prepared_source( &
                            Num_gam_e,dDR,x_edge,gam_e,dEl_step,R_mid,Q,q_left_source_shape,q_right_source_shape, &
                            prefix_source_lin_shape,dF1_shape,dN_x,dN_step)
                    end if

                    dN_x=dN_step
                    R_mid=R_mid+dDR
                end do
                R_loc=R(I_tobs)
                dNe=dNe_shell
            else
                do L=1,L1
                    R_left=R_loc
                    dNe_left=dNe
                    R_right=R_left+dDR
                    call electron_external_density(A_star,dNe_ISM,R_right,R0,R_tr,f_jump,f_wide,1,dNe_right)
                    R_mid=0.5d0*(R_left+R_right)
                    call electron_external_density(A_star,dNe_ISM,R_mid,R0,R_tr,f_jump,f_wide,1,dNe_mid)
                    DB_step=0.39d0*dsqrt(Epsilon_b*dNe_mid*(R_Gamma_loc*(R_Gamma_loc-one)))
                    Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
                    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
                    call electron_gamma_m_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
                    Gam_e_m_p_step=(one-p)/(Gam_e_max_step**(one-p)-Gam_e_m_step**(one-p))
                    call electron_injection_prefactor(R_mid,dDR,dNe_mid,f_e,Gam_e_m_p_step,Q)
                    call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)

                    if (index_Y == 0) then
                        cooling_scale=one/(beta_Gam*R_Gamma_loc)
                        a_rad=1.35d-19*DB_step**2*cooling_scale/pi
                        b_ad=one/R_mid
                        call electron_characteristic_step_affine_u(Num_gam_e,dDR,x_edge,a_rad,b_ad,dF1,dN_x,dN_step)
                    else
                        if (dNe_shell > zero) then
                            dEl_step=dEl_base*(dNe_mid/dNe_shell)
                        else
                            dEl_step=dEl_base
                        end if
                        call electron_characteristic_step_piecewise_u(Num_gam_e,dDR,x_edge,gam_e,dEl_step,R_mid,dF1,dN_x,dN_step)
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
                DB_step=0.39d0*dsqrt(Epsilon_b*dNe_mid*(R_Gamma_loc*(R_Gamma_loc-one)))
                Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
                temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
                call electron_gamma_m_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
                Gam_e_m_p_step=(one-p)/(Gam_e_max_step**(one-p)-Gam_e_m_step**(one-p))
                call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m_step,Gam_e_max_step,one,p,dF1_shape)
                call electron_prepare_conservative_remap_nonuniform(Num_gam_e,x_edge,dF1_shape, &
                                                                   q_left_source_shape,q_right_source_shape,prefix_source_lin_shape)
                if (index_Y == 0) then
                    cooling_scale=one/(beta_Gam*R_Gamma_loc)
                    a_rad=1.35d-19*DB_step**2*cooling_scale/pi
                end if
            end if

            do while (dR_left > zero)
                dR_try=min(dR_try,dR_left)
                if (dR_left < 1.5d0*dR_try) dR_try=dR_left
                step_limit=charint_substep_rtol_relax*substep_rtol
                if (step_limit > zero) then
                    dR_cap=step_limit*max(R_loc,1d-30)/max(one-0.5d0*step_limit,0.5d0)
                    if (dR_cap > dR_min) dR_try=min(dR_try,dR_cap)
                end if
                if (is_uniform_density) then
                    R_mid=R_loc+0.5d0*dR_try
                    step_error=dR_try/max(R_mid,1d-30)
                    if (step_error > charint_substep_rtol_relax*substep_rtol .and. dR_try > dR_min) then
                        dR_try=max(0.5d0*dR_try,dR_min)
                        cycle
                    end if

                    call electron_injection_prefactor(R_mid,dR_try,dNe_shell,f_e,Gam_e_m_p_step,Q)
                    if (index_Y == 0) then
                        b_ad=one/R_mid
                        call electron_characteristic_step_affine_u_prepared_source( &
                            Num_gam_e,dR_try,x_edge,a_rad,b_ad,Q,q_left_source_shape,q_right_source_shape, &
                            prefix_source_lin_shape,dF1_shape,dN_x,dN_step)
                    else
                        dEl_step=dEl_base
                        call electron_characteristic_step_piecewise_u_prepared_source( &
                            Num_gam_e,dR_try,x_edge,gam_e,dEl_step,R_mid,Q,q_left_source_shape,q_right_source_shape, &
                            prefix_source_lin_shape,dF1_shape,dN_x,dN_step)
                    end if
                else
                    R_right=R_loc+dR_try
                    call electron_external_density(A_star,dNe_ISM,R_right,R0,R_tr,f_jump,f_wide,1,dNe_right)
                    R_mid=R_loc+0.5d0*dR_try
                    call electron_external_density(A_star,dNe_ISM,R_mid,R0,R_tr,f_jump,f_wide,1,dNe_mid)
                    step_error=max(dR_try/max(R_mid,1d-30),abs(dNe_right-dNe)/max(abs(dNe),1d-99))
                    if (step_error > charint_substep_rtol_relax*substep_rtol .and. dR_try > dR_min) then
                        dR_try=max(0.9d0*step_limit*dR_try/max(step_error,1d-30),dR_min)
                        cycle
                    end if

                    DB_step=0.39d0*dsqrt(Epsilon_b*dNe_mid*(R_Gamma_loc*(R_Gamma_loc-one)))
                    Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
                    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
                    call electron_gamma_m_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
                    Gam_e_m_p_step=(one-p)/(Gam_e_max_step**(one-p)-Gam_e_m_step**(one-p))
                    call electron_injection_prefactor(R_mid,dR_try,dNe_mid,f_e,Gam_e_m_p_step,Q)
                    call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
                    if (index_Y == 0) then
                        cooling_scale=one/(beta_Gam*R_Gamma_loc)
                        a_rad=1.35d-19*DB_step**2*cooling_scale/pi
                        b_ad=one/R_mid
                        call electron_characteristic_step_affine_u(Num_gam_e,dR_try,x_edge,a_rad,b_ad,dF1,dN_x,dN_step)
                    else
                        if (dNe_shell > zero) then
                            dEl_step=dEl_base*(dNe_mid/dNe_shell)
                        else
                            dEl_step=dEl_base
                        end if
                        call electron_characteristic_step_piecewise_u(Num_gam_e,dR_try,x_edge,gam_e,dEl_step,R_mid,dF1,dN_x,dN_step)
                    end if
                end if

                if (step_error <= charint_substep_rtol_relax*substep_rtol .or. dR_try <= dR_min) then
                    dN_x=dN_step
                    R_loc=R_loc+dR_try
                    if (is_uniform_density) then
                        dNe=dNe_shell
                    else
                        dNe=dNe_right
                    end if
                    dR_left=dR_left-dR_try
                    if (dR_left > zero) then
                        if (step_error <= 0.25d0*charint_substep_rtol_relax*substep_rtol) then
                            dR_try=min(two*dR_try,dR_max)
                        end if
                    end if
                else
                    dR_try=max(0.5d0*dR_try,dR_min)
                end if
            end do
        end if
        dN_gam_e(:,I_tobs)=dN_x/gam_e/dlog(ten)
    end do

    R_loc=R(Num_R)
    R_Gamma_loc=R_Gamma(Num_R)
    call electron_external_density(A_star,dNe_ISM,R_loc,R0,R_tr,f_jump,f_wide,1,dNe)
    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
    beta_Gam=sqrt(one-one/R_Gamma_loc**2)
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(Num_R)
    V_m(Num_R)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
    V_c(Num_R)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
    call electron_prepare_radiation_spectrum(Num_gam_e,gam_e,dN_gam_e(:,Num_R), &
                                             Num_gam_rad,gam_e_rad,dN_gam_e_rad)
    call get_nu_a(R_loc,DB,Num_gam_rad,gam_e_rad(1:Num_gam_rad),dN_gam_e_rad(1:Num_gam_rad),temp)
    V_a(Num_R)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))

    deallocate(dEl_base,dEl_step,dN_x,dN_step,dF1,dF1_shape,x_edge,gam_e_rad,dN_gam_e_rad, &
               q_left_source_shape,q_right_source_shape,prefix_source_lin_shape)
end subroutine fs_electron_charint_1d

subroutine fs_electron_charint_1d_affine_step_test(Num_gam_e,dDR,x_edge,a_u,b_u,dF1,dN_x_in,dN_x_out)
    use electron_common
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: dDR,x_edge(Num_gam_e+1),a_u,b_u,dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)

    call electron_characteristic_step_affine_u(Num_gam_e,dDR,x_edge,a_u,b_u,dF1,dN_x_in,dN_x_out)
end subroutine fs_electron_charint_1d_affine_step_test
