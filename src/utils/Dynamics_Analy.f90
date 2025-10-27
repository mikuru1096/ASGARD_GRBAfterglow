include 'Constants.f90'
subroutine dynamics_analy(Boundary,Num_gam_e, R_Tobs,R_Gamma,R,gam_e,dN_gam_e)
    !$ use omp_lib
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    PARAMETER(Num_R=900,Num_R1=Num_R-1)
    PARAMETER(Num_RK4=4)
    DIMENSION Y(Num_RK4),R_Tobs(Num_R),R_Gamma(Num_R),R_m(Num_R), &
            R(Num_R),dN_gam_e(Num_R,Num_gam_e),gam_e(Num_gam_e),Boundary(15)

    !f2py intent(in) :: Boundary,Num_gam_e
    !f2py intent(out) :: R_Tobs,R_Gamma,R,gam_e,dN_gam_e
    !f2py depend(Num_gam_e) :: gam_e,dN_gam_e

    Eta_0=Boundary(1)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    E_iso=Boundary(14)
    T_log10_duration=Boundary(15)

    !***********************[Parameter Initial]**********************
    !   Y(1)=Gamma,Y(2)=m,Y(3)=thermal energy,Y(4)=R

    Y(1)=Eta_0-0.001
    Y(2)=Boundary(2)
    Y(3)=Boundary(3)
    Y(4)=Boundary(4)
    T00=Y(4)*((1.0-1.0/Eta_0**2)**(-0.5)-1.0)/Para_c
    M=4
    EPS=1.0D-5
    !**********************[Time bin]**********************
    DM_0=E_iso/((Eta_0-1.0)*4.0*pi*Para_m_p_E)
    R_dec_ISM=(dNe_ISM*Eta_0/DM_0)**(-1.0/3.0)
    if (A_star > 0.0) then
        R_dec_wind=DM_0/(2.0D35*A_star*Eta_0)
        R_dec=min(R_dec_wind,R_dec_ISM)
    else
        R_dec=R_dec_ISM
    end if
    t_dec=R_dec/(2.0*Eta_0*Eta_0*Para_c)
    Grid_Tobs_bin=min(-5.0,LOG10(t_dec*0.1))
    T_log10=T_log10_duration-Grid_Tobs_bin
    !   log time bin
    do I_tobs=1,Num_R
        T=T00+10.0**(Grid_Tobs_bin+T_log10*(I_tobs-1)/Num_R1)
        if (I_tobs-1 < 0.1) then
            H=10.0**(Grid_Tobs_bin+T_log10*(I_tobs-1)/Num_R1)
        else
            H=10.0**(Grid_Tobs_bin+T_log10*I_tobs/Num_R1)-10.0**(Grid_Tobs_bin+T_log10*(I_tobs-1.0)/Num_R1)
        end if

        call GRKT2(T,H,Y,M,EPS,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Epsilon_b,p,z)

        R_Tobs(I_tobs)=T*(1.0+z)
        R_Gamma(I_tobs)=Y(1)
        R_m(I_tobs)=Y(2)!/Para_m_p
        R(I_tobs)=Y(4)
    end do

    !*****************Part 1: given the boundary consition [Using the analytical approximation]*********************

    R_Tobs=R_Tobs-T00*(1+z)
    DO K=1,Num_gam_e
        gam_e(K)=3.0*10**(8.0*(K-1.0)/(Num_gam_e-1.0))
    ENDDO
    factor2=(p-2.0)/(p-1.0)*Epsilon_e*1836.5
    IF (p-2.05<0.001) factor2=0.05/1.05*Epsilon_e*1836.5

    DO J=1,Num_R
        DR0=R(J)
        Gamma0=R_Gamma(J)
        if (A_star >= 0.0) then
            dNe_wind=A_star*3.0d35/DR0**2
            if (dNe_wind <= dNe_ISM/4d0) then
                dNe=dNe_ISM
            else
                dNe=dNe_wind
            end if
        else
            dNe=dNe_ISM
        end if
        dB=0.39*sqrt(Epsilon_b*dNe)*sqrt(Gamma0*(Gamma0-1.0))
        Gam_e_max=3.0*Para_m_energy/sqrt(8.0*dB*Para_e**3)
        Gam_e_m=factor2*(Gamma0-1.0)+1.0
        Gam_e_c=7.739D8*(1.0+z)/Gamma0/dB**2/R_Tobs(J)
        if (J==1) then
            Gam_e_c=7.7D8/(1.0+sqrt(Epsilon_e/Epsilon_b))/Gamma0/dB**2/(R_Tobs(J)/(1.0+z))
        end if

        eta=(Gam_e_m/Gam_e_c)**(p-two)
        if (eta-one > 0.001) eta=one
        Compton=one+(-one+dsqrt(one+4d0*eta*Epsilon_e/Epsilon_b))/two
        Gam_e_c=Gam_e_c/Compton
        Gam_e_max=Gam_e_max
        
        q=R_m(J)*4.0*pi
        do K=1,Num_gam_e
            if (Gam_e_m<Gam_e_c) then
                if (gam_e(K)<=Gam_e_m.or.(gam_e(K)> Gam_e_max)) then
                    dN_gam_e(J,K)=0.0
                else if (Gam_e_m<gam_e(K) .and. gam_e(K)<Gam_e_c) then
                    dN_gam_e(J,K)=q*(1.0-p)*p/(Gam_e_c**(1.0-p)-p*Gam_e_m**(1.0-p))/gam_e(K)**p
                else
                    dN_gam_e(J,K)=q*(1.0-p)*p/(Gam_e_c**(1.0-p)-p*Gam_e_m**(1.0-p))*Gam_e_c*gam_e(K)**(-p-1.0)
                end if
            else !gam_m>gam_e_c
                if (gam_e(K)<=Gam_e_c .or. (gam_e(K)>Gam_e_max)) then
                    dN_gam_e(J,K)=0.0
                else if (Gam_e_c<gam_e(K) .and. gam_e(K)<Gam_e_m) then
                    dN_gam_e(J,K)=q*p/((1.0-p)/Gam_e_m+p/Gam_e_c)/gam_e(K)**2
                else
                    dN_gam_e(J,K)=q*p/((1.0-p)/Gam_e_m+p/Gam_e_c)*Gam_e_m**(p-1.0)*gam_e(K)**(-p-1.0)
                end if
            end if
        end do
    end do

    return
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE F(T, Y, M, D, E_e, E_iso, Eta_0, dNe_ISM, A_star, E_b, p, z)
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    DIMENSION Y(M),D(M)

    if (A_star >= 0.0) then
        dNe_wind=A_star*3.0d35/Y(4)**2
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
    end if

    Epe=E_e
    qq=0.0
    t_1=700.0/(1d0+z)
    t_2=1500.0/(1d0+z)
    if (T >= t_1 .and. T<= t_2) then
        A=1.07d50
    else
        A=0.0
    end if
    A=0.0
    dB=0.39d0*dsqrt((E_b*dNe)*(Y(1)**2-one))
    gam_c=7.739d8/(dB**2*Y(1)*T)
    gam_m=Epe*1836d0*(p-two)*(Y(1)-one)/(p-one)+one
    if ((gam_c-gam_m) > 0.001) Epe=Epe*(gam_m/gam_c)**(p-two)
    Bgam=dsqrt(one-one/Y(1)**2)

    ! **************************Huang's Dynamic model******************************
         DNe0=E_iso/(Eta_0-1.0)/1.5*1.0D3/(4d0*pi)
         D00=Bgam*2.997D10/(1.0-Bgam)
         D01=Y(1)**2-1.0
         D02=DNe0+Epe*Y(2)+2.0*(1.0-Epe)*Y(1)*Y(2)
         dM=dNe*Y(4)**2*D00
         D(1)=(-D01*dM+A/(4d0*pi)*(1.0+T/T_1)**qq/1.5*1.0D3)/D02
         D(2)=dM
         D(3)=(Y(1)-1.0)*dM
         D(4)=D00
    !! ************************Pe'er's Dynamic model*******************************
    ! *******************Caculate hat_gamma derived by Pe'er **********************
    !             four_v=Bgam*Y(1)
    !             theta=four_v/3.0*(four_v+1.07*four_v**2)
    !    $        /(1.0+four_v+1.07*four_v**2)
    !             cz=theta/(0.24+theta)
    !             hat_gam=(5.0-1.21937*cz+0.18203*cz**2-0.96583*cz**3
    !    $        +2.32513*cz**4-2.39332*cz**5+1.07136*cz**6)/3.0
    !!! *************************new equations by Pe'er*****************************
    !             DNe0=E_iso/(Eta_0-1.0)/1.5*1.0D3/12.5664
    !             D00=Bgam*2.997D10/(1.0-Bgam)
    !             D01=Y(1)**2-1.0
    !             D011=hat_gam*D01-(hat_gam-1.0)*Y(1)*Bgam**2
    !             D02=DNe0+Epe*Y(2)+2.0*(1.0-Epe)*Y(1)*Y(2)
    !             D022=DNe0+Epe*Y(2)+(1.0-Epe)*Y(2)*(2.0*hat_gam*Y(1)
    !    $         -(hat_gam-1.0)*(1.0+Y(1)**(-2)))
    !             dM=DNe*Y(4)**2*D00
    !             D(1)=-(D011/D022)*dM
    !             D(2)=dM
    !             D(3)=(Y(1)-1.0)*dM
    !             D(4)=D00
    !! ************************Zhang's Dynamic model*******************************
    ! *******************Caculate hat_gamma derived by Pe'er **********************
!    four_v=Bgam*Y(1)
!    theta=four_v/3d0*(four_v+1.07d0*four_v**2)/(one+four_v+1.07d0*four_v**2)
!    cz=theta/(0.24d0+theta)
!    hat_gam=(5d0-1.21937d0*cz+0.18203d0*cz**2-0.96583d0*cz**3+ &
!            2.32513d0*cz**4-2.39332d0*cz**5+1.07136d0*cz**6)/3d0
    ! *************************new equations by Zhang******************************
!    DNe0=E_iso/(Eta_0-one)/Para_c**2
!    D00=Bgam*Para_c/(one-Bgam)
!    D01=Y(1)**2-one
!    dM=4d0*pi*dNe*Y(4)**2*D00*Para_m_p
!    D011=-dM/D00*Para_c**2*Y(1)*D01*(hat_gam*Y(1)-hat_gam+one)- &
!            (hat_gam*D01+one)*Y(1)*(one-hat_gam)*3d0*Y(3)/Y(4)+Y(1)**2*A*(one+T/t_1)**qq/D00
!    D022=Y(1)**2*(DNe0+Y(2))*Para_c**2+(hat_gam**2*D01+3d0*hat_gam-two)*Y(3)

!    D(1)=D011/D022*D00
!    D(2)=dM
!    D(3)=((1d0-Epe)*(Y(1)-one)*dM/D00*Para_c**2- &
!            (hat_gam-one)*(3d0/Y(4)-one/Y(1)*D(1)/D00)*Y(3))*D00+A*(one+T/t_1)**qq/D00
!    D(4)=D00
    return
end subroutine

SUBROUTINE GRKT2(T,H,Y,M,EPS,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0)
    IMPLICIT REAL(8)(A-H,O-Z)
    DIMENSION Y(M),D(M),A(4),B(M),C(M),G(M),E(M)

    HH=H
    N=1
    P=1+EPS
    X=T
    C=Y
    do while (P >= EPS)
        A(1)=HH/2.0
        A(2)=A(1)
        A(3)=HH
        A(4)=HH
        G=Y
        Y=C
        DT=H/N
        T=X
        do J=1,N
            call F(T,Y,M,D,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0)
            E=Y
            B=Y
            do K=1,3
                Y=E+A(K)*D
                B=B+A(K+1)*D/3.0
                TT=T+A(K)
                call F(TT,Y,M,D,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0)
            end do
            Y=B+HH*D/6.0
            T=T+DT
        end do
        P=0.0
        do I=1,M
            Q=2.0*abs(Y(I)-G(I))/(Y(I)+G(I))
            if (Q > P) P=Q
        end do
        HH=HH/2.0
        N=N+N
    end do
    T=X

    return
end subroutine GRKT2

