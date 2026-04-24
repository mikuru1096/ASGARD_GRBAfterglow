subroutine fs_electron_weno5_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,n_threads, &
                                gam_e,dN_gam_e,P_syn,Seed_syn)
    !$ use omp_lib
    use constants
    use electron_common
    use electron_injection_profiles, only: electron_initial_powerlaw_exp_cutoff, electron_build_source_term_exp_cutoff
    use get_Y
    use electron_forward_kernel, only: electron_forward_cooling_step
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_R,Num_gam_e
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R),Seed_syn(Num_nu,Num_R)

    real(8),allocatable,dimension (:) :: f_r_times_gam_e,dEl,para_minus_gam_e_p,dEl1,x, &
            dN_x,dF1,fp,flux,Compton,Compton1,dot_gam_e,dot_gam_e_SSA
    real(8),allocatable,dimension (:,:) :: temp_store,temp_store_extended
    real(8),allocatable,dimension (:) :: dN_x_extended,fp_extended,flux_extended,dEl1_extended
    allocate (f_r_times_gam_e(Num_gam_e),dEl(Num_gam_e),para_minus_gam_e_p(Num_gam_e), &
            dEl1(Num_gam_e),x(Num_gam_e),dN_x(Num_gam_e),dF1(Num_gam_e),fp(Num_gam_e), &
            flux(0:Num_gam_e), temp_store(3,Num_gam_e),Compton(Num_gam_e), &
              dot_gam_e(Num_gam_e),dot_gam_e_SSA(Num_gam_e),Compton1(Num_gam_e))
    allocate(dN_x_extended(1-3:Num_gam_e+3),temp_store_extended(3, 1-3:Num_gam_e+3),&
             fp_extended(1-3:Num_gam_e+3),flux_extended(0:Num_gam_e),dEl1_extended(1-3:Num_gam_e+3)) !ghost cells

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
    
    !*****************Part 1: given the boundary consition [Using the analytical approximation]*********************
    call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)

    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam=Epsilon_e/f_e*1836d0*(R_Gamma(1)-one)
    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    call electron_build_gamma_grid(Num_gam_e,Gam_e_max_max,Gam_e)
    call electron_initial_powerlaw_exp_cutoff(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,Gam_e,dN_gam_e(:,1))
    !*******************Part 1 is completed [has been checked and there is no bug]**********************************
    !*******************Part 2: To calculate the electron distribution**********************************************
    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
    d_x=dlog10(gam_e(2)/gam_e(1))
    factor_adv=Para_sigmaT/(6.0d0*pi*Para_m_energy)
    para_minus_gam_e_p=gam_e**(-p)*gam_e*dlog(ten)
    
    CFL_target=0.8d0

    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        call electron_external_density(A_star,dNe_ISM,R_loc,R0,R_tr,f_jump,f_wide,0,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*1836d0*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)

        beta_Gam=dsqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1/(f_r*Gam_e_max+1.333/(R(I_tobs)+R(I_tobs-1)))
        !***********************[Here we have presented the choice on Delta_r]******************************************
        dDD=R(I_tobs)-R(I_tobs-1)
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)
        
        call get_syn(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e(:,I_tobs-1),V_seed, &
                     P_syn(:,I_tobs),Seed_syn(:,I_tobs))
        
        call electron_forward_cooling_step(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
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
        
        do L=1,L1
            R_loc=R_loc+dDR
            
            call electron_external_density(A_star,dNe_ISM,R_loc,R0,R_tr,f_jump,f_wide,1,dNe)
            
            dEl1=(dEl+one/R_loc)/dlog(ten)
            call electron_injection_prefactor(R_loc,dDR,dNe,f_e,Gam_e_m_p,Q)
            call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,Q,p,dF1)
!            dF1=dF1+Q*para_maxwell/Gam_e_m_p*(one-f_e)
            
            
            dN_x_extended(1-3:0) = dN_x(1)
            dN_x_extended(1:Num_gam_e) = dN_x
            dN_x_extended(Num_gam_e+1:Num_gam_e+3) = dN_x(Num_gam_e)
            temp_store_extended(1,:) = dN_x_extended
            
            dEl1_extended(1-3:0) = dEl1(1)
            dEl1_extended(1:Num_gam_e) = dEl1
            dEl1_extended(Num_gam_e+1:Num_gam_e+3) = dEl1(Num_gam_e)

            do j=1,3

              call update_ghost_cells(dN_x_extended, Num_gam_e)
        
              fp_extended = dEl1_extended * dN_x_extended
        
              do i_gam_e = 0, Num_gam_e
                 if (dEl1_extended(i_gam_e) <= 0.0d0) then
                    flux_extended(i_gam_e) = fpx(fp_extended(i_gam_e-2:i_gam_e+2))
                 else
                    flux_extended(i_gam_e) = fmx(fp_extended(i_gam_e-1:i_gam_e+3))
                 end if
              end do
        
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
           end do
           
           dN_x = dN_x_extended(1:Num_gam_e)
            
           where(dN_x < 0.0d0) dN_x = 0.0d0
        
           if (L1 == L) then
              dN_gam_e(:,I_tobs)=dN_x/gam_e/dlog(ten)
           end if
        end do
    end do
 
    deallocate (f_r_times_gam_e,dEl,para_minus_gam_e_p,dEl1,x,dN_x,dF1,fp,flux,temp_store,Compton,Compton1,dot_gam_e, &
                dN_x_extended, temp_store_extended, fp_extended, flux_extended)
    
    return
end subroutine

subroutine update_ghost_cells(arr, n)
        integer, intent(in) :: n
        real(8), intent(inout) :: arr(1-3:n+3)
        
        arr(1-3:0) = arr(1)
        arr(n+1:n+3) = arr(n)
        
        ! Periodic boundary conditions, not used.
        ! arr(1-2) = arr(n-1)
        ! arr(1-1) = arr(n)
        ! arr(n+1) = arr(1)
        ! arr(n+2) = arr(2)
    end subroutine update_ghost_cells

function fpx(fps)
    real(8) :: fps(-2:2), fpx
    real(8) :: omega(3), fu(3), beta(3)
    real(8) :: tao5, totalpha, alpha(3), fomega(3), eps

    omega(1) = 0.1d0;   omega(2) = 0.6d0;   omega(3) = 0.3d0
    eps=1d-30
    
    if(any(isnan(fps))) then
        fpx = 0.0d0
        return
    end if
    
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
        fpx = fu(2)
    else
        fomega(:) = alpha(:)/totalpha
        fpx = fu(1)*fomega(1) + fu(2)*fomega(2) + fu(3)*fomega(3)
    end if
    
    return
end function fpx

function fmx(fms)
    real(8) :: fms(-1:3), fmx
    real(8) :: omega(3), fu(3), beta(3)
    real(8) :: tao5, totalpha, alpha(3), fomega(3), eps
    
    omega(1) = 0.1d0;   omega(2) = 0.6d0;   omega(3) = 0.3d0
    eps=1d-30
    
    if(any(isnan(fms))) then
        fmx = 0.0d0
        return
    end if
    
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
        fmx = fu(2)
    else
        fomega(:) = alpha(:)/totalpha
        fmx = fu(1)*fomega(1) + fu(2)*fomega(2) + fu(3)*fomega(3)
    end if
    
    return
end function fmx
