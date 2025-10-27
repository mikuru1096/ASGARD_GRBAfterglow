include 'Constants.f90'

subroutine fs_electron_fullhide (Boundary,R_Tobs,R_Gamma,R,n,Num_R,Num_gam_e, gam_e,dN_gam_e)
    !$ use omp_lib
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_R,Num_gam_e
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e)
    
    real(8),allocatable,dimension (:) :: gam_e,gam_e2,dgam_e,temp_N_e,f_r_times_gam_e2, &
            dEl,dEl1,para_minus_gam_e_p,cal_para_minus_gam_e_p,dF1,dgam_e_times_dDR
    allocate (gam_e(Num_gam_e),gam_e2(Num_gam_e),dgam_e(Num_gam_e-1),temp_N_e(Num_gam_e), &
            f_r_times_gam_e2(Num_gam_e),dEl(Num_gam_e),dEl1(Num_gam_e),para_minus_gam_e_p(Num_gam_e), &
            cal_para_minus_gam_e_p(Num_gam_e),dF1(Num_gam_e),dgam_e_times_dDR(Num_gam_e-1))


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
    
    !*****************Part 1: given the boundary consition [Using the analytical approximation]*********************
    if (A_star >= 0.0) then
        dNe_wind=A_star*3.0D35/R(1)**2
        if (dNe_wind <= dNe_ISM/4.0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
    end if
    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    Gam_e_m=(p-two)/(p-one)*Epsilon_e*1836d0*(R_Gamma(1)-one)+one
    if (p-2.05 < 0.001) then
        Gam_e_m=0.05d0/1.05d0*Epsilon_e*1836d0*(R_Gamma(1)-one)+one
    end if
    Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    do I_gam_e=1,Num_gam_e
        Gam_e(I_gam_e)=3d0*ten**(dlog10(Gam_e_max_max)*(I_gam_e-1)/(Num_gam_e-1))
       if (Gam_e_m-Gam_e_c > 0.1) then
           Q1=R_m(1)*Gam_e_c
           if ((Gam_e_c-Gam_e(I_gam_e)) > 1.0) then
               dN_gam_e(1,I_gam_e)=0.0
           else
               if ((Gam_e_m-Gam_e(I_gam_e)) > 1.0) then
                   dN_gam_e(1,I_gam_e)=Q1*Gam_e(I_gam_e)**(-2)
               else
                   if ((Gam_e_max-Gam_e(I_gam_e)) > 1.0) then
                       dN_gam_e(1,I_gam_e)=Q1*Gam_e_m**(p-1.0)*Gam_e(I_gam_e)**(-(p+1.0))
                   else
                       dN_gam_e(1,I_gam_e)=0.0
                   end if
               end if
           end if
       else
           Q1=R_m(1)*Gam_e_m**(p-1.0)
           if ((Gam_e_m-Gam_e(I_gam_e)) > 1.0) then
               dN_gam_e(1,I_gam_e)=0.0
           else
               if ((Gam_e_c-Gam_e(I_gam_e)) > 1.0) then
                   dN_gam_e(1,I_gam_e)=Q1*Gam_e(I_gam_e)**(-p)
               else
                   if ((Gam_e_max-Gam_e(I_gam_e)) > 1.0) then
                       dN_gam_e(1,I_gam_e)=Q1*Gam_e_c*Gam_e(I_gam_e)**(-(p+1.0))
                   else
                       dN_gam_e(1,I_gam_e)=0.0
                   end if
               end if
           end if
       end if
    end do
!*******************Part 1 is completed [has been checked and there is no bug]**********************************	
!*******************Part 2: To calculate the electron distribution**********************************************
    gam_e2=gam_e**2
    dgam_e=1.0/(gam_e(2:Num_gam_e)-gam_e(1:Num_gam_e-1)+0.01)
    para_minus_gam_e_p=1.0/(gam_e-1.0)**p

    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/2.0
        if (A_star >= 0.0) then
            dNe_wind=A_star*3.0D35/R_loc**2
            if (dNe_wind <= dNe_ISM/4.0) then
                dNe=dNe_ISM
            else
                dNe=dNe_wind
            end if
        else
            dNe=dNe_ISM
        end if
        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        Gam_e_m=(p-two)/(p-one)*Epsilon_e*1836d0*(R_Gamma_loc-one)+one
        if (p-2.05 < 0.001) then
            Gam_e_m=0.05/1.05*Epsilon_e*1836d0*(R_Gamma_loc-one)+one
        end if
        Gam_e_m_p=(p-one)*(Gam_e_m-one)**(p-one)
        Gam_e_c=7.7D8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        eta=(Gam_e_m/Gam_e_c)**(p-2.0)
        if (eta-1.0 > 0.001) eta=1.0
        beta_Gam=sqrt(1.0-1.0/R_Gamma_loc**2)
        f_r=(1.35D-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        f_r_times_gam_e2=f_r*gam_e2
        dDR=0.1/(f_r*Gam_e_max+1.333/(R(I_tobs)+R(I_tobs-1)))
        !***********************[Here we have presented the choice on Delta_r]******************************************
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=Int(dDD/dDR)
        L1=L1+100
        dDR=dDD/L1
        temp_N_e=dN_gam_e(I_tobs-1,:)
        dgam_e_times_dDR=dgam_e*dDR
        !$OMP PARALLEL  NUM_THREADS(8)
        !$OMP DO
        do i_gam_e=1,Num_gam_e
            if (gam_e(i_gam_e) < Gam_e_max .and. gam_e(i_gam_e) > Gam_e_m) then
                cal_para_minus_gam_e_p(i_gam_e)=para_minus_gam_e_p(i_gam_e)
            else
                cal_para_minus_gam_e_p(i_gam_e)=0.0
            end if
            if (i_gam_e < Num_gam_e) then
                !*****************************[The general inverse Compton effect]******************************
                hat_gam=6.5D6/sqrt(DB*gam_e(i_gam_e+1))
                if (Gam_e_m-Gam_e_c > 0.0) then
                    if (hat_gam-Gam_e_c < 0.0) then
                        eta_NK=0.0
                    else
                        if (hat_gam-Gam_e_m < 0.0) then
                            Step1=(p-1.0)/(p-2.0)*Gam_e_m-Gam_e_c
                            eta_NK=(hat_gam-Gam_e_c)/Step1
                        else
                            Step2=Gam_e_m**(p-1.0)*hat_gam**(2.0-p)
                            Step3=(p-1.0)*Gam_e_m-(p-2.0)*Gam_e_c
                            eta_NK=1.0-Step2/Step3
                        end if
                    end if
                else
                    if (hat_gam-Gam_e_m < 0.0) then
                        eta_NK=0.0
                    else
                        if (hat_gam-Gam_e_c < 0.0) then
                            Step4=Gam_e_c**(3.0-p)/(p-2.0)-Gam_e_m**(3.0-p)
                            eta_NK=(hat_gam**(3.0-p)-Gam_e_m**(3.0-p))/Step4
                        else
                            Step5=(3.0-p)*Gam_e_c*hat_gam**(2.0-p)
                            Step6=Gam_e_c**(3.0-p)-(p-2.0)*Gam_e_m**(3.0-p)
                            eta_NK=1.0-Step5/Step6
                        end if
                    end if
                end if
                Compton=1.0+(-1.0+sqrt(1.0+4.0*eta*eta_NK*Epsilon_e/Epsilon_b))/2.0
                dEl(i_gam_e)=f_r_times_gam_e2(i_gam_e)*Compton
            else
                dEl(i_gam_e)=f_r_times_gam_e2(i_gam_e)*Compton
            end if
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        do L=1,L1
            R_loc=R_loc+dDR
            dEl1=dEl+gam_e/R_Loc
            if (A_star >= 0.0) then
                dNe_wind=A_star*3.0D35/R_loc**2
                if (dNe_wind <= dNe_ISM/4.0) then
                    dNe=dNe_ISM
                else
                    dNe=dNe_wind
                end if
            else
                dNe=dNe_ISM
            end if
            Q=4.0*pi*R_loc*R_loc*dNe*Gam_e_m_p  !here Q is Q_0*\gamma_m**p
            dF1=Q*dDR*cal_para_minus_gam_e_p
            temp_N_e=temp_N_e+dF1
            dF3=(temp_N_e(2:Num_gam_e)/gam_e(2:Num_gam_e)**2-temp_N_e(1:Num_gam_e-1)/gam_e(1:Num_gam_e-1)**2) &
                        *dgam_e_times_dDR(1:Num_gam_e-1)*gam_e(1:Num_gam_e-1)**2*para_ssa
            do i_gam_e=1,Num_gam_e
                if (i_gam_e < Num_gam_e) then
                    dF2=(temp_N_e(i_gam_e+1)*dEl1(i_gam_e+1)-temp_N_e(i_gam_e)*dEl1(i_gam_e) &
                        +dF3(i_gam_e)-dF3(i_gam_e))*dgam_e_times_dDR(i_gam_e)
               !need SSA
                    temp_N_e(i_gam_e)=temp_N_e(i_gam_e)+dF2
                else
                    temp_N_e(i_gam_e)=temp_N_e(i_gam_e)+dF2
                end if
                !           Compton=0.0
                !**********************************************************************
                if (temp_N_e(i_gam_e) < 0.0) then
                    temp_N_e(i_gam_e)=0.0
                end if
            end do
        end do

        do i_gam_e=1,Num_gam_e
            if (temp_N_e(i_gam_e) < 0.0) then
                dN_gam_e(I_tobs,i_gam_e)=0.0
            else
                dN_gam_e(I_tobs,i_gam_e)=temp_N_e(i_gam_e)
            end if
        end do
    end do

    return
end subroutine
