!Calculate the electron distibutions of forward shock.
!New improvement at 11.29.2022
include 'calling_modules.f90'
!****************************************************************************************
!******************************* main program *******************************************
!****************************************************************************************
subroutine fs_electron_fullhide(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                                gam_e,dN_gam_e,P_syn,Seed_syn)
    !$ use omp_lib
    use constants
    use get_Y
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R),Seed_syn(Num_nu,Num_R)
    
    real(8),allocatable,dimension (:) :: dEl,dEL_mean,principal,x,dF1,up,para_minus_gam_e_p,dot_gam_e_SSA, &
                                         dN_x,temp1,temp2,temp3,temp4,para_maxwell,Compton,Compton1,dot_gam_e
    allocate (dEl(Num_gam_e),dEL_mean(Num_gam_e-1),principal(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e), &
              up(Num_gam_e-1),dN_x(Num_gam_e),temp1(Num_gam_e-1),temp2(Num_gam_e),para_maxwell(Num_gam_e), &
              temp3(Num_gam_e-1),temp4(Num_gam_e-1),para_minus_gam_e_p(Num_gam_e),Compton(Num_gam_e), &
              dot_gam_e(Num_gam_e),dot_gam_e_SSA(Num_gam_e),Compton1(Num_gam_e))
    
    !***********************[Parameter Initial]**********************
    Eta_0=Boundary(1)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    E_iso=Boundary(14)
    f_e=Boundary(16)
    R0=Boundary(n)
    
    P_syn=zero
    Seed_syn=zero
    
    !*****************Part 1: given the boundary consition [Using the analytical approximation]*********************
    if (A_star > zero) then
        dNe_wind=A_star*3.0d35/R(1)**2
        Para_N_e_ini=4d0*pi*R_ini*A_star*3.0d35
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
        Para_N_e_ini=4d0/3d0*pi*R_ini**3*dNe_ISM
    end if
    
    if (R(1)<R0) then
        dNe=A_star*3.0d35/R0**2*4
        Para_N_e_ini=4d0/3d0*pi*R_ini**3*dNe_ISM
    end if

    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    Gam_e_m=(p-two)/(p-one)*Epsilon_e/f_e*1836d0*(R_Gamma(1)-one)+one
    if (p<2.01 .and. p>=2.0) then
        Gam_e_m=0.01d0/1.01d0*Epsilon_e/f_e*1836d0*(R_Gamma(1)-one)+one
    else if (p<2 .and. p>1) then
        Gam_e_m=((two-p)/(p-one)*Epsilon_e/f_e*1836d0*(R_Gamma(1)-one)*Gam_e_max**(p-two))**(one/(p-one))+one
    end if
    Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    do I_gam_e=1,Num_gam_e
        Gam_e(I_gam_e)=3d0*ten**(dlog10(Gam_e_max_max)*(I_gam_e-1)/(Num_gam_e-1))
        if (Gam_e_m > Gam_e_c) then
            if (Gam_e_c > Gam_e(I_gam_e) .or. Gam_e_max < Gam_e(I_gam_e)) then
                dN_gam_e(I_gam_e,1)=zero
            else
                Q1=Para_N_e_ini*Gam_e_c
                if (Gam_e_m > Gam_e(I_gam_e)) then
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e(I_gam_e)**(-2)
                else
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e_m**(p-one)*Gam_e(I_gam_e)**(-(p+one))
                end if
            end if
        else
            if (Gam_e_m > Gam_e(I_gam_e) .or. Gam_e_max < Gam_e(I_gam_e)) then
                dN_gam_e(I_gam_e,1)=zero
            else
                Q1=Para_N_e_ini*Gam_e_m**(p-one)
                if (Gam_e_c > Gam_e(I_gam_e)) then
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e(I_gam_e)**(-p)
                else
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e_c*Gam_e(I_gam_e)**(-(p+one))
                end if
            end if
        end if
    end do
    !*******************Part 1 is completed [has been checked and there is no bug]**********************************
    !*******************Part 2: To calculate the electron distribution**********************************************
    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
    d_x=dlog10(gam_e(2)/gam_e(1))
!    factor_adv=Para_sigmaT/(6.0d0*pi*Para_m_energy)
    para_minus_gam_e_p=one/(gam_e-one)**p*gam_e*dlog(ten)
    
    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        if (A_star > zero) then
            dNe_wind=A_star*3.0d35/R_loc**2
            if (dNe_wind <= dNe_ISM/4d0) then
                dNe=dNe_ISM
            else
                dNe=dNe_wind
            end if
        else
            dNe=dNe_ISM
        end if
        
        if (R_loc<R0) then
            dNe=A_star*3.0d35/R0**2
        end if

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        Gam_e_m=(p-two)/(p-one)*Epsilon_e*1836d0*(R_Gamma_loc-one)/f_e+one
        if (p<2.05 .and. p>=2.0) then
            Gam_e_m=0.05d0/1.05d0*Epsilon_e*1836d0*(R_Gamma_loc-one)/f_e+one
        else if (p<2 .and. p>1) then
            Gam_e_m=((two-p)/(p-one)*Epsilon_e/f_e*1836d0*(R_Gamma_loc-one)*Gam_e_max**(p-two))**(one/(p-one))+one
        end if
        Gam_e_m_p=(p-one)*(Gam_e_m-one)**(p-one)
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)

        beta_Gam=dsqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1/(f_r*Gam_e_max+1.333/(R(I_tobs)+R(I_tobs-1)))
        !***********************[Here we have presented the choice on Delta_r]******************************************
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=max(100,min(1000,Int(dDD/dDR)))
        dDR=dDD/L1
        CFL=dDR/d_x
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)
        
        !Compton_max=one+(-one+dsqrt(one+4d0*eta*Epsilon_e/Epsilon_b))/two
        !Gam_e_max=Gam_e_max!/sqrt(Compton_max)
!        Compton = zero

        select case(index_syn_intger)
        
        case(1)
        call get_syn(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e(:,I_tobs-1),V_seed, &
                     P_syn(:,I_tobs),Seed_syn(:,I_tobs))
                     
        case(2)
        call get_syn_simpson(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e(:,I_tobs-1),V_seed, &
                     P_syn(:,I_tobs),Seed_syn(:,I_tobs))
                     
        end select
        
        call get_SSA_numerical(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn(:,I_tobs), dot_gam_e_SSA)
        
        select case(index_Y)
        
        case(1)
        call get_IC_numerical(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn(:,I_tobs), &
                       dot_gam_e)
        
        dEl=(f_r+(dot_gam_e-dot_gam_e_SSA)/beta_Gam/R_Gamma_loc/para_c)*gam_e
        
        case(2)
        call get_Y_Nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn(:,I_tobs), &
                         Compton)

        Q=4d0*pi*R_loc*R_loc*para_c
        Compton=one+Compton/Q/(4d0*R_Gamma_loc*R_Gamma_loc*dNe*Para_m_p_E)
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA/beta_Gam/R_Gamma_loc/para_c)*gam_e
        
        case(3)
        call get_Y_Fan(Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e, &
                       Compton)
        Compton=one+Compton
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA/beta_Gam/R_Gamma_loc/para_c)*gam_e
        
        case default
         
        print*, 'invalid Compton case, check your chosen model!'
        stop
         
        end select
         
!        four_v=R_Gamma_loc*beta_Gam
!        theta=four_v/3d0*(four_v+1.07*four_v*four_v)/(one+four_v+1.07*four_v*four_v)
!        theta=max(theta,2d-1)
!        para_maxwell=gam_e*gam_e*dsqrt(one-one/gam_e**2)/theta/besselk(1d0/theta)* &
!                     dexp(-gam_e/theta)
!        para_normalize=sum((para_maxwell(2:Num_gam_e)+para_maxwell(1:Num_gam_e-1))* &
!                       (gam_e(2:Num_gam_e)-gam_e(1:Num_gam_e-1)))/two
!        para_maxwell=para_maxwell/para_normalize*gam_e*dlog(ten)


        dEL_mean=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two/dlog(ten)

        do L=1,L1
            R_loc=R_loc+dDR
            Q=4d0/3d0*pi*(3d0*R_loc**2+dDR*(3d0*R_loc+dDR))*dNe*f_e*Gam_e_m_p  !here Q is Q_0*\gamma_m**p
            dF1=zero
            where(gam_e<Gam_e_max .and. gam_e>Gam_e_m) dF1=Q*para_minus_gam_e_p!*f_e
!            dF1=dF1+Q*para_maxwell/Gam_e_m_p*(one-f_e)

            temp3=dEL_mean+one/R_loc/dlog(ten)
            up=-CFL*temp3 !up
            principal(2:Num_gam_e)=one-up !main
            principal(1)=principal(2)

            temp1=up/(principal(2:Num_gam_e)+principal(1:Num_gam_e-1))*two
            temp2=(dN_x+dDR*dF1)/principal
            x(Num_gam_e)=temp2(Num_gam_e)
            do i=Num_gam_e-1,1,-1
                x(i)=max(zero, temp2(i)-temp1(i)*x(i+1))
            end do
            dN_x=x

            if (L1 == L) then
                dN_gam_e(:,I_tobs)=dN_x/gam_e/dlog(ten)
            end if
        end do
    end do

    deallocate (dEl,dEL_mean,principal,x,dF1,up,dN_x,temp1,temp2,para_maxwell,temp3,temp4,para_minus_gam_e_p,Compton,Compton1)

    return
end subroutine fs_electron_fullhide
