subroutine dynamics_reverse (Delta_t,e_r,b_r,p_r,f_e_r,Boundary,n,Num_R,Num_gam_e, &
                             T_cross,R_cross,e3_cross,gam20,R_Tobs,R_Gamma,R,M2,M3,B3,gam_e,dN_gam_e)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_deceleration_radius, dynamics_external_density_base
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_R,Num_gam_e
    real(8), intent(in) :: Boundary(n)
    real(8), intent(in) :: Delta_t,e_r,b_r,p_r,f_e_r
    real(8), intent(out) :: T_cross,R_cross,e3_cross,gam20,R_Tobs(Num_R),R(Num_R),M2(Num_R),M3(Num_R),B3(Num_R)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),R_Gamma(Num_R)
    real(8), parameter :: reverse_gamma_c_coeff=7.7d8
    real(8), parameter :: reverse_synch_b_coeff=0.39d0
    real(8), parameter :: reverse_adv_coeff=1.35d-19

    real(8),allocatable,dimension (:) :: dEl,principal,x,dF1,up,temp1,temp2,temp3,Y,dB3_serial, &
                                         dN_x,para_minus_gam_e_p
    allocate (dEl(Num_gam_e),principal(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e),up(Num_gam_e-1), &
              temp1(Num_gam_e-1),temp2(Num_gam_e),temp3(Num_gam_e-1),dN_x(Num_gam_e),Y(4), &
              para_minus_gam_e_p(Num_gam_e),dB3_serial(Num_R))


    Eta_0=Boundary(1)
    R(1)=Boundary(4)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p_f=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    E_iso=Boundary(14)
    T_log10_duration=Boundary(15)
    f_e=Boundary(16)
    
    Delta_0=Delta_t*para_c
    para_m_ej=E_iso/eta_0/para_c**2
    
    if (A_star > zero) then
        para_m2=4d0*pi*R(1)*A_star*3d35*para_m_p
    else
        para_m2=4d0/3d0*pi*R(1)**3*dNe_ISM*para_m_p
    end if
    para_m3=1d1
    R_Gamma(1)=Eta_0-0.001
    dB3_serial(1)=zero
    
    Y=[R_Gamma(1),R(1),para_m2,para_m3]
    
    !**********************[Time bin]**********************
    DM_0=E_iso/((Eta_0-one)*4d0*pi*Para_m_p_E)
    call dynamics_deceleration_radius(A_star,dNe_ISM,Eta_0,DM_0,R_dec)
    
    T_cross=-1d0
    T00=Y(4)*(one/dsqrt(one-one/Eta_0**2)-one)/Para_c

    t_dec=R_dec/(two*Eta_0*Eta_0*Para_c)
    Grid_Tobs_bin=min(log10(T00)-2.0,dlog10(t_dec*0.1))
    T_log10=T_log10_duration-Grid_Tobs_bin
    Num_R1=Num_R-1
    !   log time bin

    do I_tobs=1,Num_R
        T=T00+ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        if (I_tobs == 1) then
            H=ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        else
            H=ten**(Grid_Tobs_bin+T_log10*I_tobs/Num_R1)-ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        end if
        
        call GRKT4(dB3,T_cross,R_cross,e3_cross,gam20, T,H,Y,para_m_ej,Delta_0,eta_0,A_star, &
                 dNe_ISM,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r)
        R_Tobs(I_tobs)=T*(one+z)
        R_Gamma(I_tobs)=Y(1)
        R(I_tobs)=Y(2)
        M2(I_tobs)=Y(3)
        M3(I_tobs)=Y(4)
        dB3_serial(I_tobs)=dB3
    end do
    B3 = dB3_serial

    factor2=(p_r-two)/(p_r-one)*e_r*Para_m_p_div_m_e
    if (p_r<2.05) factor2=0.05/1.05*e_r*Para_m_p_div_m_e

    dB3_serial(1)=dB3_serial(min(2,Num_R))
    dB=dB3_serial(1)
    gamma34=1.001
    call reverse_gamma_extrema(dB,gamma34,factor2,f_e_r,Gam_e_max,Gam_e_m)
    Gam_e_c=reverse_gamma_c_coeff/(one+dsqrt(e_r/b_r))/R_Gamma(1)/dB**2/(R_Tobs(1)/two)
    
    
    call dynamics_external_density_base(A_star,dNe_ISM,R(1),dNe)

    DB_min=reverse_synch_b_coeff*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    
    do I_gam_e=1,Num_gam_e
        Gam_e(I_gam_e)=3d0*ten**(dlog10(Gam_e_max_max)*(I_gam_e-1)/(Num_gam_e-1))
        call reverse_initial_electron_bin(Gam_e(I_gam_e),Gam_e_m,Gam_e_c,Gam_e_max,p_r,dN_gam_e(I_gam_e,1))
    end do
    !*******************Part 1 is completed [has been checked and there is no bug]**********************************
    !*******************Part 2: To calculate the electron distribution**********************************************
    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
    d_x=dlog10(gam_e(2)/gam_e(1))
    factor_adv=Para_sigmaT/(6.0*pi*Para_m_energy)
    para_minus_gam_e_p=one/(gam_e-one)**p_r
    
    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        Delta=max(Delta_0,R_loc/Eta_0**2)
        R_n4=para_m_ej/(4d0*pi*Para_m_p*R_loc*R_loc*Eta_0*Delta)
        beta4=dsqrt(one-one/eta_0**2)
        beta2=dsqrt(one-one/R_Gamma_loc**2)
        gamma34=(one-beta2*beta4)*eta_0*R_Gamma_loc
!        dB=0.39*sqrt(b_r*R_n4)*gamma34
        dB=(dB3_serial(I_tobs)+dB3_serial(I_tobs-1))/two
        call reverse_gamma_extrema(dB,gamma34,factor2,f_e_r,Gam_e_max,Gam_e_m)
        Gam_e_c=reverse_gamma_c_coeff*(one+z)/R_Gamma_loc/dB**2/R_Tobs(I_tobs)
        eta=(Gam_e_m/Gam_e_c)**(p_r-two)
        if (eta-one > 0.001) eta=one
        f_r=reverse_adv_coeff/beta2/R_Gamma_loc*dB**2/pi
        dDR=0.7/(f_r*Gam_e_max+1.333/(R(I_tobs)+R(I_tobs-1)))
        !***********************[Here we have presented the choice on Delta_r]******************************************
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=max(100,min(1000,Int(dDD/dDR)))
        dDR=dDD/L1
        CFL=dDR/d_x
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)
        
        Compton_max=one+(-one+dsqrt(one+4d0*eta*Epsilon_e/Epsilon_b))/two
        Gam_e_max=Gam_e_max/sqrt(Compton_max)
        
        DO i_gam_e=1,Num_gam_e
            IF (Num_gam_e-i_gam_e > 0.1) THEN
!*****************************[The general inverse Compton effect]******************************
                hat_gam=5.4246D6*sqrt(R_Gamma_loc/(DB*gam_e(i_gam_e+1)))
                call reverse_compton_factor(eta,e_r,b_r,p_r,Gam_e_m,Gam_e_c,hat_gam,Compton)
!**********************************************************************
                f_r1=f_r*(1.0+Compton)
                dEl(i_gam_e)=f_r1*gam_e(i_gam_e)
            else
                dEl(i_gam_e)=f_r1*gam_e(i_gam_e)+0.1
            end if
        END DO
        
        Q0=4d0*pi*R_n4*(p_r-one)*(Gam_e_m-one)**(p_r-one)*f_e_r
        DO L=1,L1
            R_loc=R_loc+dDR
            
            if (R_cross >= R_loc) then
                Q=Q0*R_loc*R_loc  !here Q is Q_0*\gamma_m**p
            else
                Q=zero
            end if
            
            dF1=zero
            where(gam_e<Gam_e_max .and. gam_e>Gam_e_m) dF1=Q*para_minus_gam_e_p*gam_e*dlog(ten)
            
            temp3=((dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two+one/R_loc)/dlog(ten)
            up=-CFL*temp3 !up
            principal(2:Num_gam_e)=one-up !main
            principal(1)=principal(2)

            temp1=up/principal(2:Num_gam_e)!+principal(1:Num_gam_e-1))*two
            temp2=(dN_x+dDR*dF1)/principal
            x(Num_gam_e)=temp2(Num_gam_e)
            do i=Num_gam_e-1,1,-1
                x(i)=max(zero, temp2(i)-temp1(i)*x(i+1))
            end do
            dN_x=x
            
            IF (L1 == L) THEN
                dN_gam_e(:,I_tobs)=dN_x/gam_e/dlog(ten)
            END IF
        end do
    end do

    deallocate (dEl,principal,x,dF1,up,temp1,temp2,temp3,dN_x,Y,para_minus_gam_e_p,dB3_serial)


    return
end subroutine dynamics_reverse

subroutine reverse_gamma_extrema(dB,gamma34,factor2,f_e_r,Gam_e_max,Gam_e_m)
    use constants
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: dB,gamma34,factor2,f_e_r
    real(8), intent(out) :: Gam_e_max,Gam_e_m

    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
    Gam_e_m=factor2*(gamma34-one)/f_e_r+one
end subroutine reverse_gamma_extrema

subroutine reverse_initial_electron_bin(Gam_e_val,Gam_e_m,Gam_e_c,Gam_e_max,p_r,dN_val)
    use constants
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: Gam_e_val,Gam_e_m,Gam_e_c,Gam_e_max,p_r
    real(8), intent(out) :: dN_val

    if (Gam_e_m > Gam_e_c) then
        Q1=1d10*Gam_e_c
        if ((Gam_e_c-Gam_e_val) > one) then
            dN_val=zero
        else
            if (Gam_e_m > Gam_e_val) then
                dN_val=Q1*Gam_e_val**(-2)
            else
                if (Gam_e_max > Gam_e_val) then
                    dN_val=Q1*Gam_e_m**(p_r-one)*Gam_e_val**(-(p_r+one))
                else
                    dN_val=zero
                end if
            end if
        end if
    else
        Q1=1d10*Gam_e_m**(p_r-one)
        if (Gam_e_m > Gam_e_val) then
            dN_val=zero
        else
            if (Gam_e_c > Gam_e_val) then
                dN_val=Q1*Gam_e_val**(-p_r)
            else
                if (Gam_e_max > Gam_e_val) then
                    dN_val=Q1*Gam_e_c*Gam_e_val**(-(p_r+one))
                else
                    dN_val=zero
                end if
            end if
        end if
    end if
end subroutine reverse_initial_electron_bin

subroutine reverse_compton_factor(eta,e_r,b_r,p_r,Gam_e_m,Gam_e_c,hat_gam,Compton)
    use constants
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: eta,e_r,b_r,p_r,Gam_e_m,Gam_e_c,hat_gam
    real(8), intent(out) :: Compton

    IF (Gam_e_m-Gam_e_c > 0.0d0) THEN
        IF (hat_gam-Gam_e_c < 0.0d0) THEN
            eta_NK=0.0d0
        ELSE
            IF (hat_gam-Gam_e_m < 0.0d0) THEN
                Step1=(p_r-1d0)/(p_r-2d0)*Gam_e_m-Gam_e_c
                eta_NK=(hat_gam-Gam_e_c)/Step1
            ELSE
                Step2=Gam_e_m**(p_r-1d0)*hat_gam**(2d0-p_r)
                Step3=(p_r-1d0)*Gam_e_m-(p_r-2d0)*Gam_e_c
                eta_NK=1d0-Step2/Step3
            END IF
        END IF
    ELSE
        IF (hat_gam-Gam_e_m < 0.0d0) THEN
            eta_NK=0.0d0
        ELSE
            IF (hat_gam-Gam_e_c < 0.0d0) THEN
                Step4=Gam_e_c**(3d0-p_r)/(p_r-2.0d0)-Gam_e_m**(3d0-p_r)
                eta_NK=(hat_gam**(3d0-p_r)-Gam_e_m**(3d0-p_r))/Step4
            ELSE
                Step5=(3d0-p_r)*Gam_e_c*hat_gam**(2d0-p_r)
                Step6=Gam_e_c**(3.0d0-p_r)-(p_r-2d0)*Gam_e_m**(3.0d0-p_r)
                eta_NK=1d0-Step5/Step6
            END IF
        END IF
    END IF
    Compton=(-1.0d0+dsqrt(1.0d0+4.0d0*eta*eta_NK*e_r/b_r))/2.0d0
end subroutine reverse_compton_factor

subroutine F(dB3,T_cross,R_cross,e3_cross,gam20, &
             T,Y,D,para_m_ej,Delta_0,eta_0,A_star,dNe_ISM,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r)
    use constants
    use dynamics_common, only: dynamics_external_density_base
    IMPLICIT REAL(8)(A-H,O-Z)
    real(8), intent(inout) :: dB3,T_cross,R_cross,e3_cross,gam20
    real(8), intent(in) :: T,para_m_ej,Delta_0,eta_0,A_star,dNe_ISM,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r
    real(8), intent(in) :: Y(4)
    real(8), intent(out) :: D(4)
    real(8), parameter :: reverse_synch_b_coeff=0.39d0
    real(8), parameter :: reverse_gamma_c_precise_coeff=7.739d8
    
    gam2=Y(1)
    RR=Y(2)
    para_m2=Y(3)
    para_m3=Y(4)
    
    call dynamics_external_density_base(A_star,dNe_ISM,RR,dNe)
    
    u2=dsqrt(gam2*gam2-one)
    u4=dsqrt(eta_0*eta_0-one)
    Delta=max(Delta_0,RR/eta_0**2)
    para1=4d0*pi*Para_m_p*RR*RR
    para_n4=para_m_ej/(para1*eta_0*Delta)
    
    beta4=u4/eta_0
    beta2=u2/gam2
    gam34=eta_0*gam2-u2*u4
    
    para_n3=(4.0*gam34+3.0)*para_n4
    betars=(u2*para_n3-u4*para_n4)/(gam2*para_n3-eta_0*para_n4)

    dB2=reverse_synch_b_coeff*dsqrt((Epsilon_b*dNe)*(gam2*gam2-one))
    gam_c2=reverse_gamma_c_precise_coeff/(dB2*dB2*gam2*T)
    gam_m2=Epsilon_e/f_e*Para_m_p_div_m_e*(p_f-two)*(gam2-one)/(p_f-one)+one
    eps2=Epsilon_e*min(one,(gam_m2/gam_c2)**(p_f-two))

    e2=4d0*gam2*gam2*dNe*Para_m_p_E
    if (T_cross < zero) then
        e3=e2
        e3_cross=e2
        dB3=dB2
        gam_c3=gam_c2
    else
        e3=e3_cross*(R_cross/RR)**3*gam2/gam20
        dB3=dsqrt(8d0*pi*b_r*e3)
        gam_c3=reverse_gamma_c_precise_coeff/(dB3*dB3*gam2*T)
    end if
    
    gam_m3=e_r/f_e_r*Para_m_p_div_m_e*(p_r-two)*(gam34-one)/(p_r-one)+one
    eps3=e_r*min(one,(gam_m3/gam_c3)**(p_r-two))

    dgam2_1=-para1*((gam2*gam2-one)*dNe+(gam2*gam34-eta_0)*(beta4-betars)*eta_0*para_n4)
    dgam2_2=(para_m2+para_m3+(one-eps2)*(two*gam2-one)*para_m2+(one-eps3)*(gam34-one)*para_m3+ &
 &          (one-eps3)*gam2*para_m3*(eta_0*(one-beta4*beta2)-eta_0*beta4/(gam2**2*beta2)))
    dgam2=dgam2_1/dgam2_2

    dR=beta2/(one-beta2)*para_c
    dm2=para1*dNe*dR
    if (para_m_ej > para_m3) then
        dm3=para1*(beta4-betars)*eta_0*para_n4*dR
    else
        if (T_cross < zero) then
            T_cross = T
            R_cross = RR
            gam20 = gam2
        end if
        dm3=zero
        dgam2_1=-u2**2*dm2/dR
        dgam2_2=para_m_ej+(eps2+two*(one-eps2)*gam2)*para_m2
        dgam2=dgam2_1/dgam2_2
    end if
    dgam2=dgam2*dR
    D=[dgam2,dR,dm2,dm3]
    
end subroutine F

SUBROUTINE GRKT4(dB3,T_cross,R_cross,e3_cross,gam20, T,H,Y,para_m_ej,Delta_0,eta_0,A_star, &
                 dNe_ISM,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r)
    use dynamics_common, only: dynamics_rk4_coefficients, dynamics_rk4_error_n
    IMPLICIT REAL(8)(A-H,O-Z)
    real(8), intent(inout) :: dB3,T_cross,R_cross,e3_cross,gam20,T,Y(4)
    real(8), intent(in) :: H,para_m_ej,Delta_0,eta_0,A_star,dNe_ISM,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r
    real(8) :: D(4),A(4),B(4),C(4),G(4),E(4)
    
    EPS=1d-5
    HH=H
    N=1
    P=1d0+EPS
    X=T
    C=Y
    do while (P >= EPS)
        call dynamics_rk4_coefficients(HH,A)
        G=Y
        Y=C
        DT=H/N
        T=X
        do J=1,N
            call F(dB3,T_cross,R_cross,e3_cross,gam20, &
                   T,Y,D,para_m_ej,Delta_0,eta_0,A_star,dNe_ISM,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r)
            E=Y
            B=Y
            do K=1,3
                Y=E+A(K)*D
                B=B+A(K+1)*D/3.0
                TT=T+A(K)
                call F(dB3,T_cross,R_cross,e3_cross,gam20, &
                       TT,Y,D,para_m_ej,Delta_0,eta_0,A_star,dNe_ISM,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r)
            end do
            Y=B+HH*D/6.0
            T=T+DT
        end do
        call dynamics_rk4_error_n(Y,G,4,P)
        HH=0.5*HH
        N=N+N
    end do
    T=X

    return
end subroutine GRKT4
