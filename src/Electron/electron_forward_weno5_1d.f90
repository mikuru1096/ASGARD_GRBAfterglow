! 电子1D WENO5格式主驱动：三阶TVD Runge-Kutta + WENO5通量重构。
subroutine fs_electron_weno5_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,n_threads, &
                                gam_e,dN_gam_e,P_syn,Seed_syn)
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    use electron_common
    use electron_injection_profiles, only: electron_build_source_term_exp_cutoff_edges
    use radiation_common, only: radiation_syn_seed_chi_batch_core
    use electron_cooling_kernel, only: get_forward_cooling
    use electron_transport_common, only: electron_dnx_to_dndgamma_exp_centers
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,n_threads
    integer :: I_tobs,L,L1
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R),Seed_syn(Num_nu,Num_R)
    real(8) :: Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star,E_iso,T_log10_duration,f_e
    real(8) :: R_tr,f_jump,f_wide,R0,dNe,Para_N_e_ini,DB,Gam_e_max,DB_min,Gam_e_max_max
    real(8) :: temp_gam,Gam_e_m,Gam_e_c,d_x,CFL_target,R_loc,R_Gamma_loc,Gam_e_m_p
    real(8) :: beta_Gam,f_r,dDR,dDD,WENO_speed_max,CFL,Q

    real(8),allocatable,dimension (:) :: dEl,dEl1,dN_x,x_edge,dF1
    real(8),allocatable,dimension (:,:) :: temp_store_extended
    real(8),allocatable,dimension (:) :: dN_x_extended,fp_extended,flux_extended,dEl1_extended
    allocate (dEl(Num_gam_e),dEl1(Num_gam_e),dN_x(Num_gam_e),x_edge(Num_gam_e+1),dF1(Num_gam_e))
    allocate(dN_x_extended(1-3:Num_gam_e+3),temp_store_extended(3, 1-3:Num_gam_e+3),&
             fp_extended(1-3:Num_gam_e+3),flux_extended(0:Num_gam_e),dEl1_extended(1-3:Num_gam_e+3)) !ghost cells

    !***********************[Parameter Initial]**********************
    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    
    P_syn=zero
    Seed_syn=zero
    
    !*****************Part 1: given the boundary consition [Using the analytical approximation]*********************
    call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)

    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam=Epsilon_e/f_e*1836d0*(R_Gamma(1)-one)
    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    call electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                      electron_initial_grid_log_edges,gam_e,dN_x,x_edge)
    call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,1))
    !*******************Part 1 is completed [has been checked and there is no bug]**********************************
    !*******************Part 2: To calculate the electron distribution**********************************************
    d_x=dlog10(gam_e(2)/gam_e(1))
    
    CFL_target=0.8d0

    do I_tobs=2,Num_R
        call prepare_weno_shell(I_tobs)
        call write_weno_radiation_and_cooling(I_tobs)

        do L=1,L1
            call advance_weno_substep(I_tobs,L)
        end do
    end do

    deallocate (dEl,dEl1,dN_x,x_edge,dF1,dN_x_extended,temp_store_extended,fp_extended,flux_extended,dEl1_extended)

contains

    subroutine prepare_weno_shell(I_tobs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs

        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,0,R_tr,f_jump,f_wide,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*1836d0*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)

        beta_Gam=dsqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1/(f_r*Gam_e_max+1.333/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)
    end subroutine prepare_weno_shell

    subroutine write_weno_radiation_and_cooling(I_tobs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs
    real(8) :: DB_chi(1),DNe_chi(Num_gam_e,1),P_emit_shell(Num_nu,1),P_syn_shell(Num_nu,1)
    real(8) :: Seed_syn_shell(Num_nu,1),Tau_syn_shell(Num_nu,1)

        DB_chi(1)=DB
        DNe_chi(:,1)=dN_gam_e(:,I_tobs-1)
        call radiation_syn_seed_chi_batch_core(R_loc,Num_gam_e,Num_nu,1,gam_e,DNe_chi,V_seed,DB_chi,1.046d4, &
                                               P_emit_shell,P_syn_shell,Seed_syn_shell,Tau_syn_shell)
        P_syn(:,I_tobs)=P_syn_shell(:,1)
        Seed_syn(:,I_tobs)=Seed_syn_shell(:,1)

        call get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                 R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                 P_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl)

        dEl(1:Num_gam_e-1)=(dEl(1:Num_gam_e-1)+dEl(2:Num_gam_e))*0.5d0
        dEl(Num_gam_e)=dEl(Num_gam_e-1)*0.5d0
        dEl1=(dEl+one/R_loc)/dlog(ten)

        L1=Int(dDD/dDR)+10
        WENO_speed_max=maxval(abs(dEl1))
        L1=max(L1,ceiling(dDD*WENO_speed_max/(CFL_target*d_x)))
        dDR=dDD/L1
        CFL=dDR/d_x
    end subroutine write_weno_radiation_and_cooling

    subroutine advance_weno_substep(I_tobs,L)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_tobs,L

        R_loc=R_loc+dDR
        
        call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,1,R_tr,f_jump,f_wide,dNe)
        
        dEl1=(dEl+one/R_loc)/dlog(ten)
        call electron_injection_prefactor(R_loc,dDR,dNe,f_e,Gam_e_m_p,Q)
        call electron_build_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m,Gam_e_max,Q,p,dF1)

        call load_weno_extended_state()

        do j=1,3
            call compute_weno_fluxes()
            call advance_weno_rk_stage(j)
        end do

        dN_x = dN_x_extended(1:Num_gam_e)

        if (L1 == L) then
            call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
        end if
    end subroutine advance_weno_substep

    subroutine load_weno_extended_state()
    implicit real(8)(A-H,O-Z)

        dN_x_extended(1-3:0) = dN_x(1)
        dN_x_extended(1:Num_gam_e) = dN_x
        dN_x_extended(Num_gam_e+1:Num_gam_e+3) = dN_x(Num_gam_e)
        temp_store_extended(1,:) = dN_x_extended

        dEl1_extended(1-3:0) = dEl1(1)
        dEl1_extended(1:Num_gam_e) = dEl1
        dEl1_extended(Num_gam_e+1:Num_gam_e+3) = dEl1(Num_gam_e)
    end subroutine load_weno_extended_state

    subroutine compute_weno_fluxes()
    implicit real(8)(A-H,O-Z)

        call weno5_update_ghost_cells(dN_x_extended, Num_gam_e)

        fp_extended = dEl1_extended * dN_x_extended

        do i_gam_e = 0, Num_gam_e
            if (dEl1_extended(i_gam_e) <= 0.0d0) then
                flux_extended(i_gam_e) = weno5_positive_flux(fp_extended(i_gam_e-2:i_gam_e+2))
            else
                flux_extended(i_gam_e) = weno5_negative_flux(fp_extended(i_gam_e-1:i_gam_e+3))
            end if
        end do
    end subroutine compute_weno_fluxes

    subroutine advance_weno_rk_stage(j)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: j

        if(j==1) then
            do i = 1, Num_gam_e
                dN_x_extended(i) = temp_store_extended(1,i) + CFL*(flux_extended(i)-flux_extended(i-1)) + dF1(i)*dDR
            end do
            temp_store_extended(2,:) = dN_x_extended
        else if(j==2) then
            do i = 1, Num_gam_e
                dN_x_extended(i) = 0.75d0*temp_store_extended(1,i) + 0.25d0*(temp_store_extended(2,i) + &
                               CFL*(flux_extended(i)-flux_extended(i-1))) + 0.25d0*dF1(i)*dDR
            end do
            temp_store_extended(3,:) = dN_x_extended
        else if(j==3) then
            do i = 1, Num_gam_e
                dN_x_extended(i) = (temp_store_extended(1,i) + 2.0d0*(temp_store_extended(3,i) + &
                               CFL*(flux_extended(i)-flux_extended(i-1))))/3.0d0 + 2.0d0/3.0d0*dF1(i)*dDR
            end do
        end if
    end subroutine advance_weno_rk_stage

! 更新WENO5鬼点：零阶外推（复制边界值）。
subroutine weno5_update_ghost_cells(arr, n)
    implicit none
    integer, intent(in) :: n
    real(8), intent(inout) :: arr(1-3:n+3)
        
    arr(1-3:0) = arr(1)
    arr(n+1:n+3) = arr(n)
end subroutine weno5_update_ghost_cells

! WENO5正通量重构 f⁺(f_{i-2},...,f_{i+2})。
function weno5_positive_flux(fps)
    implicit none
    real(8), intent(in) :: fps(-2:2)
    real(8) :: weno5_positive_flux
    real(8) :: omega(3), fu(3), beta(3)
    real(8) :: tao5, totalpha, alpha(3), fomega(3), eps

    omega(1) = 0.1d0;   omega(2) = 0.6d0;   omega(3) = 0.3d0
    eps=1d-30
    
    fu(1) =  1.0d0/3.0d0*fps(-2) - 7.0d0/6.0d0*fps(-1) + 11.0d0/6.0d0*fps(0)
    fu(2) = -1.0d0/6.0d0*fps(-1) + 5.0d0/6.0d0*fps(0)  + 1.0d0/3.0d0*fps(1)
    fu(3) =  1.0d0/3.0d0*fps(0)  + 5.0d0/6.0d0*fps(1)  - 1.0d0/6.0d0*fps(2)
    
    beta(1) = 13.0d0/12.0d0*( fps(-2) - 2.0d0*fps(-1) + fps(0) )**2 &
            + 0.25d0*( fps(-2) - 4.0d0*fps(-1) + 3.0d0*fps(0) )**2
    beta(2) = 13.0d0/12.0d0*( fps(-1) - 2.0d0*fps(0) + fps(1) )**2 &
            + 0.25d0*( fps(1) - fps(-1) )**2
    beta(3) = 13.0d0/12.0d0*( fps(0) - 2.0d0*fps(1) + fps(2) )**2 &
            + 0.25d0*( 3.0d0*fps(0) - 4.0d0*fps(1) + fps(2) )**2
    
    tao5 = abs(beta(1) - beta(3))
    
    alpha(:) = omega(:)*( 1.0d0 + (tao5/(beta(:)+eps))**2 )
    totalpha = alpha(1) + alpha(2) + alpha(3)
    
    if(totalpha < eps) then
        weno5_positive_flux = fu(2)
    else
        fomega(:) = alpha(:)/totalpha
        weno5_positive_flux = fu(1)*fomega(1) + fu(2)*fomega(2) + fu(3)*fomega(3)
    end if
    
end function weno5_positive_flux

! WENO5负通量重构 f⁻(f_{i-1},...,f_{i+3})。
function weno5_negative_flux(fms)
    implicit none
    real(8), intent(in) :: fms(-1:3)
    real(8) :: weno5_negative_flux
    real(8) :: omega(3), fu(3), beta(3)
    real(8) :: tao5, totalpha, alpha(3), fomega(3), eps
    
    omega(1) = 0.1d0;   omega(2) = 0.6d0;   omega(3) = 0.3d0
    eps=1d-30
    
    fu(1) =  1.0d0/3.0d0*fms(3) - 7.0d0/6.0d0*fms(2) + 11.0d0/6.0d0*fms(1)
    fu(2) = -1.0d0/6.0d0*fms(2) + 5.0d0/6.0d0*fms(1) + 1.0d0/3.0d0*fms(0)
    fu(3) =  1.0d0/3.0d0*fms(1) + 5.0d0/6.0d0*fms(0) - 1.0d0/6.0d0*fms(-1)
    
    beta(1) = 13.0d0/12.0d0*( fms(3) - 2.0d0*fms(2) + fms(1) )**2 &
            + 0.25d0*( fms(3) - 4.0d0*fms(2) + 3.0d0*fms(1) )**2
    beta(2) = 13.0d0/12.0d0*( fms(2) - 2.0d0*fms(1) + fms(0) )**2 &
            + 0.25d0*( fms(2) - fms(0) )**2
    beta(3) = 13.0d0/12.0d0*( fms(1) - 2.0d0*fms(0) + fms(-1) )**2 &
            + 0.25d0*( 3.0d0*fms(1) - 4.0d0*fms(0) + fms(-1) )**2
    
    tao5 = abs(beta(1) - beta(3))
    
    alpha(:) = omega(:)*( 1.0d0 + (tao5/(beta(:)+eps))**2 )
    totalpha = alpha(1) + alpha(2) + alpha(3)
    
    if(totalpha < eps) then
        weno5_negative_flux = fu(2)
    else
        fomega(:) = alpha(:)/totalpha
        weno5_negative_flux = fu(1)*fomega(1) + fu(2)*fomega(2) + fu(3)*fomega(3)
    end if
    
end function weno5_negative_flux
end subroutine fs_electron_weno5_1d
