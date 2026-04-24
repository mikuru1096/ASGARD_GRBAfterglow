!Calculate the electron distributions of forward shock.
!Modified to use Implicit Three-Level Method (second-order time accuracy)
!Modified at 2024 by AI assistant based on original code from 11.29.2022
!****************************************************************************************
!******************************* main program *******************************************
!****************************************************************************************
subroutine fs_electron_t2g1_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                            gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    !$ use omp_lib
    use constants
    use electron_common
    use electron_injection_profiles, only: electron_initial_powerlaw, electron_build_source_term_profile
    use get_Y
    use electron_forward_kernel, only: electron_forward_cooling_step
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R), &
                            Seed_syn(Num_nu,Num_R), V_m(Num_R), V_c(Num_R), V_a(Num_R)
    
    real(8),allocatable,dimension (:) :: dEl,dEL_mean,dEL_mean_base,principal,x,dF1,up,para_minus_gam_e_p,dot_gam_e_SSA, &
                                         dN_x,dN_x_prev,temp1,temp2,temp3,temp4,para_maxwell,Compton,Compton1,dot_gam_e, &
                                         gam_e_rad,dN_gam_e_rad
    integer :: Num_gam_rad
    allocate (dEl(Num_gam_e),dEL_mean(Num_gam_e-1),dEL_mean_base(Num_gam_e-1),principal(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e), &
              up(Num_gam_e-1),dN_x(Num_gam_e),dN_x_prev(Num_gam_e),temp1(Num_gam_e-1), &
              temp2(Num_gam_e),para_maxwell(Num_gam_e),temp3(Num_gam_e-1),temp4(Num_gam_e-1), &
              para_minus_gam_e_p(Num_gam_e),Compton(Num_gam_e),dot_gam_e(Num_gam_e),dot_gam_e_SSA(Num_gam_e), &
              Compton1(Num_gam_e),gam_e_rad(Num_gam_e),dN_gam_e_rad(Num_gam_e))
    
    !***********************[Parameter Initial]**********************
    Eta_0=Boundary(1)
    R_ini=Boundary(4)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    E_iso=Boundary(14)
    T_log10_duration=Boundary(15)
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

    !*****************Part 1: given the boundary condition [Using the analytical approximation]*********************
    call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)

    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam=Epsilon_e/f_e*1836d0*(R_Gamma(1)-one)
    call electron_gamma_m_near_two(p,2.01d0,0.01d0,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    call electron_build_gamma_grid(Num_gam_e,Gam_e_max_max,Gam_e)
    call electron_initial_powerlaw(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,Gam_e,dN_gam_e(:,1))
    !*******************Part 1 is completed [has been checked and there is no bug]**********************************
    !*******************Part 2: To calculate the electron distribution**********************************************
    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
    dN_x_prev = dN_x
    d_x=dlog10(gam_e(2)/gam_e(1))
    para_minus_gam_e_p=one/(gam_e-one)**p*gam_e*dlog(ten)

    ! For the first few steps, we need to use single-step methods
    ! We'll use a startup procedure
    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        call electron_external_density(A_star,dNe_ISM,R_loc,R0,R_tr,f_jump,f_wide,0,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*1836d0*(R_Gamma_loc-one)
        call electron_gamma_m_near_two(p,2.05d0,0.05d0,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(p-one)*(Gam_e_m-one)**(p-one)
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        dNe_shell=dNe

        beta_Gam=sqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1/(f_r*Gam_e_max+1.333/(R(I_tobs)+R(I_tobs-1)))
        !***********************[Here we have presented the choice on Delta_r]******************************************
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=max(100,min(1000,Int(dDD/dDR)))
        dDR=dDD/L1
        CFL=dDR/d_x
        ! For the first step of each I_tobs, we use the solution from previous step
        if (I_tobs == 2) then
            dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
            dN_x_prev = dN_x
        else
            dN_x=dN_x  ! Keep the solution from previous I_tobs
        end if
        call electron_prepare_radiation_spectrum(Num_gam_e,gam_e,dN_gam_e(:,I_tobs-1), &
                                                 Num_gam_rad,gam_e_rad,dN_gam_e_rad)
        
        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        call get_nu_a(R_loc,DB,Num_gam_rad,gam_e_rad(1:Num_gam_rad),dN_gam_e_rad(1:Num_gam_rad),temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))

        call get_syn_selected(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                              gam_e,dN_gam_e(:,I_tobs-1),V_seed, &
                              P_syn(:,I_tobs),Seed_syn(:,I_tobs))
        
        call electron_forward_cooling_step(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                 R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                 P_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl)
         
        call electron_loss_mean(Num_gam_e,dEl,dEL_mean)
        dEL_mean_base=dEL_mean

        ! Main loop for sub-steps
        do L=1,L1
            R_loc=R_loc+dDR
            
            call electron_external_density(A_star,dNe_ISM,R_loc,R0,R_tr,f_jump,f_wide,1,dNe)
            DB_step=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
            Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
            temp_gam=Epsilon_e/f_e*1836d0*(R_Gamma_loc-one)
            call electron_gamma_m_near_two(p,2.05d0,0.05d0,temp_gam,Gam_e_max_step,Gam_e_m_step)
            Gam_e_m_p_step=(p-one)*(Gam_e_m_step-one)**(p-one)
            
            call electron_injection_prefactor(R_loc,dDR,dNe,f_e,Gam_e_m_p_step,Q)
            call electron_build_source_term_profile(Num_gam_e,gam_e,Gam_e_m_step,Gam_e_max_step,Q,para_minus_gam_e_p,dF1)
            if (dNe_shell > zero) then
                dEL_mean=dEL_mean_base*(dNe/dNe_shell)
            else
                dEL_mean=dEL_mean_base
            end if

            temp3=dEL_mean+one/R_loc/dlog(ten)
            up=-CFL*temp3
            
            if (I_tobs == 2 .and. L <= 2) then
                call electron_prepare_implicit_coeffs(Num_gam_e,one,up,principal,temp1)
                temp2 = (dN_x + dDR * dF1) / principal
            else
                call electron_prepare_implicit_coeffs(Num_gam_e,1.5d0,up,principal,temp1)
                temp2 = ( (2d0)*dN_x - 0.5d0*dN_x_prev + dF1 * dDR ) / principal
            end if

            call electron_backward_sweep(Num_gam_e,temp1,temp2,x)
            
            dN_x_prev = dN_x
            dN_x = x

            if (L1 == L) then
                dN_gam_e(:,I_tobs)=dN_x/gam_e/dlog(ten)
            end if
        end do
    end do

    deallocate (dEl,dEL_mean,dEL_mean_base,principal,x,dF1,up,para_minus_gam_e_p,dN_x,dN_x_prev,temp1,temp2, &
                para_maxwell,temp3,temp4,Compton,Compton1,gam_e_rad,dN_gam_e_rad)

    return
end subroutine fs_electron_t2g1_1d
