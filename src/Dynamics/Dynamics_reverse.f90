subroutine dynamics_reverse (Delta_t,e_r,b_r,p_r,f_e_r,Boundary,n,Num_R,Num_gam_e, &
                             T_cross,R_cross,e3_cross,gam20,R_Tobs,R_Gamma,R,gam_e,dN_gam_e)
    !$ use omp_lib
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_R,Num_gam_e
    real(8), intent(in) :: Boundary(n)
    real(8), intent(in) :: Delta_t,e_r,b_r,p_r,f_e_r
    real(8), intent(out) :: T_cross,R_cross,e3_cross,gam20,R_Tobs(Num_R),R(Num_R)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),R_Gamma(Num_R)

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
    
    R_dec_ISM=(dNe_ISM*Eta_0/DM_0)**(-one/3d0)
    if (A_star > zero) then
        R_dec_wind=DM_0/(2.0d35*A_star*Eta_0)
        R_dec=min(R_dec_wind,R_dec_ISM)
    else
        R_dec=R_dec_ISM
    end if
    
    T_cross=-1d0
    T00=Y(4)*(one/dsqrt(one-one/Eta_0**2)-one)/Para_c

    t_dec=R_dec/(two*Eta_0*Eta_0*Para_c)
    Grid_Tobs_bin=min(log10(T00)-2.0,dlog10(t_dec*0.1))
    T_log10=T_log10_duration-Grid_Tobs_bin
    !   log time bin

    do I_tobs=1,Num_R
        T=T00+ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/(Num_R-1))
        if (I_tobs < one) then
            H=ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/(Num_R-1))
        else
            H=ten**(Grid_Tobs_bin+T_log10*I_tobs/(Num_R-1))-ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/(Num_R-1))
        end if
        
        call GRKT4(dB3,T_cross,R_cross,e3_cross,gam20, T,H,Y,para_m_ej,Delta_0,eta_0,A_star, &
                 dNe_ISM,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r)
        R_Tobs(I_tobs)=T*(one+z)
        R_Gamma(I_tobs)=Y(1)
        R(I_tobs)=Y(2)
        dB3_serial(I_tobs)=dB3
    end do

    factor2=(p_r-two)/(p_r-one)*e_r*1836.5
    if (p_r<2.05) factor2=0.05/1.05*e_r*1836.5

    dB3_serial(1)=dB3_serial(2)
    dB=dB3_serial(1)
    gamma34=1.001
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
    Gam_e_m=(p_r-two)/(p_r-one)*e_r/f_e_r*1836d0*(gamma34-one)+one
    if (p_r-2.05 < 0.001) then
        Gam_e_m=0.05d0/1.05d0*e_r*1836d0*(gamma34-one)+one
    end if
    Gam_e_c=7.7d8/(one+dsqrt(e_r/b_r))/R_Gamma(1)/dB**2/(R_Tobs(1)/two)
    
    
    if (A_star > zero) then
        dNe_wind=A_star*3.0d35/R(1)**2
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
    end if

    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    
    do I_gam_e=1,Num_gam_e
        Gam_e(I_gam_e)=3d0*ten**(dlog10(Gam_e_max_max)*(I_gam_e-1)/(Num_gam_e-1))
        if (Gam_e_m > Gam_e_c) then
            Q1=1d10*Gam_e_c
            if ((Gam_e_c-Gam_e(I_gam_e)) > one) then
                dN_gam_e(I_gam_e,1)=zero
            else
                if (Gam_e_m > Gam_e(I_gam_e)) then
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e(I_gam_e)**(-2)
                else
                    if (Gam_e_max > Gam_e(I_gam_e)) then
                        dN_gam_e(I_gam_e,1)=Q1*Gam_e_m**(p_r-one)*Gam_e(I_gam_e)**(-(p_r+one))
                    else
                        dN_gam_e(I_gam_e,1)=zero
                    end if
                end if
            end if
        else
            Q1=1d10*Gam_e_m**(p_r-one)
            if (Gam_e_m > Gam_e(I_gam_e)) then
                dN_gam_e(I_gam_e,1)=zero
            else
                if (Gam_e_c > Gam_e(I_gam_e)) then
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e(I_gam_e)**(-p_r)
                else
                    if (Gam_e_max > Gam_e(I_gam_e)) then
                        dN_gam_e(I_gam_e,1)=Q1*Gam_e_c*Gam_e(I_gam_e)**(-(p_r+one))
                    else
                        dN_gam_e(I_gam_e,1)=zero
                    end if
                end if
            end if
        end if
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
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        Gam_e_m=(p_r-two)/(p_r-one)*e_r*1836d0*(gamma34-one)/f_e_r+one
        if (p_r-2.05 < 0.001) then
            Gam_e_m=0.05d0/1.05d0*e_r*1836d0*(gamma34-one)/f_e_r+one
        end if
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/dB**2/R_Tobs(I_tobs)
        eta=(Gam_e_m/Gam_e_c)**(p_r-two)
        if (eta-one > 0.001) eta=one
        f_r=1.35d-19/beta2/R_Gamma_loc*dB**2/pi
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
                IF (Gam_e_m-Gam_e_c > 0.0) THEN
                    IF (hat_gam-Gam_e_c < 0.0) THEN
                        eta_NK=0.0
                    ELSE
                        IF (hat_gam-Gam_e_m < 0.0) THEN
                            Step1=(p_r-1)/(p_r-2)*Gam_e_m-Gam_e_c
                            eta_NK=(hat_gam-Gam_e_c)/Step1
                        ELSE
                            Step2=Gam_e_m**(p_r-1)*hat_gam**(2-p_r)
                            Step3=(p_r-1)*Gam_e_m-(p_r-2)*Gam_e_c
                            eta_NK=1-Step2/Step3
                        END IF
                    END IF
                ELSE
                    IF (hat_gam-Gam_e_m < 0.0) THEN
                        eta_NK=0.0
                    ELSE
                        IF (hat_gam-Gam_e_c < 0.0) THEN
                            Step4=Gam_e_c**(3-p_r)/(p_r-2.0)-Gam_e_m**(3-p_r)
                            eta_NK=(hat_gam**(3-p_r)-Gam_e_m**(3-p_r))/Step4
                        ELSE
                            Step5=(3-p_r)*Gam_e_c*hat_gam**(2-p_r)
                            Step6=Gam_e_c**(3.0-p_r)-(p_r-2)*Gam_e_m**(3.0-p_r)
                            eta_NK=1-Step5/Step6
                        END IF
                    END IF
                END IF
                Compton=(-1.0+dsqrt(1.0+4.0*eta*eta_NK*e_r/b_r))/2.0
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

subroutine F(dB3,T_cross,R_cross,e3_cross,gam20, &
             T,Y,D,para_m_ej,Delta_0,eta_0,A_star,dNe_ISM,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r)
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    DIMENSION Y(4),D(4)
    
    gam2=Y(1)
    RR=Y(2)
    para_m2=Y(3)
    para_m3=Y(4)
    
    if (A_star >= 0.0) then
        dNe_wind=A_star*3.0D35/RR**2
        if (dNe_wind <= dNe_ISM/4.0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
    end if
    
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

    dB2=0.39d0*dsqrt((Epsilon_b*dNe)*(gam2*gam2-one))
    gam_c2=7.739d8/(dB2*dB2*gam2*T)
    gam_m2=Epsilon_e/f_e*1836d0*(p_f-two)*(gam2-one)/(p_f-one)+one
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
        gam_c3=7.739d8/(dB3*dB3*gam2*T)
    end if
    
    gam_m3=e_r/f_e_r*1836d0*(p_r-two)*(gam2-one)/(p_r-one)+one
    eps3=e_r*min(one,(gam_m3/gam_c3)**(p_r-two))

    dgam2_1=-para1*((gam2*gam2-one)*dNe+(gam2*gam34-eta_0)*(beta4-betars)*eta_0*para_n4)
    dgam2_2=(para_m2+para_m3+(one-eps2)*(two*gam2-one)*para_m2+(one-eps3)*(gam34-one)*para_m3+ &
            (one-eps3)*gam2*para_m3*(eta_0*(one-beta4*beta2)-eta_0*beta4/(gam2**2*beta2)))
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
    IMPLICIT REAL(8)(A-H,O-Z)
    DIMENSION Y(4),D(4),A(4),B(4),C(4),G(4),E(4)
    
    EPS=1d-5
    HH=H
    N=1
    P=1d0+EPS
    X=T
    C=Y
    do while (P >= EPS)
        A(1)=0.5*HH
        A(2)=A(1)
        A(3)=HH
        A(4)=HH
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
        P=0.0
        do I=1,M
            Q=2.0*abs(Y(I)-G(I))/(Y(I)+G(I))
            if (Q > P) P=Q
        end do
        HH=0.5*HH
        N=N+N
    end do
    T=X

    return
end subroutine GRKT4

