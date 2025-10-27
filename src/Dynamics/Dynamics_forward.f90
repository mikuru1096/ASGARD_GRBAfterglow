subroutine dynamics_forward(Boundary,n,Num_R,index_dyn, R_Tobs,R_Gamma,R)
    !$ use omp_lib
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_R,index_dyn
    real(8), intent(in) :: Boundary(n)
    real(8), intent(out) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R)

    real(8),dimension (4) :: Y,D,B,C,G,E
    real(8),dimension (Num_R) :: R_m
    !***********************[Parameter Initial]**********************
    Eta_0=Boundary(1)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    E_iso=Boundary(14)
    T_log10_duration=Boundary(15)
    f_e=Boundary(16)
    R0=Boundary(n)
    
    Num_R1=Num_R-1
    
    !***********************[Parameter Initial]**********************
    !   Y(1)=Gamma,Y(2)=m,Y(3)=thermal energy,Y(4)=R
    Y(1)=Eta_0-0.001
    Y(2)=Boundary(2)
    Y(3)=Boundary(3)
    Y(4)=Boundary(4)
    T00=Y(4)*(one/dsqrt(one-one/Eta_0**2)-one)/Para_c
    M=4
    EPS=1.0D-5
    !**********************[Time bin]**********************
    DM_0=E_iso/((Eta_0-one)*4d0*pi*Para_m_p_E)
    R_dec_ISM=(dNe_ISM*Eta_0/DM_0)**(-one/3d0)
    if (A_star > zero) then
        R_dec_wind=DM_0/(2.0d35*A_star*Eta_0)
        R_dec=min(R_dec_wind,R_dec_ISM)
    else
        R_dec=R_dec_ISM
    end if
    t_dec=R_dec/(two*Eta_0*Eta_0*Para_c)
    Grid_Tobs_bin=min(-5d0,dlog10(t_dec*0.1))
    T_log10=T_log10_duration-Grid_Tobs_bin
    !   log time bin
    do I_tobs=1,Num_R
        T=ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        if (I_tobs < one) then
            H=ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        else
            H=ten**(Grid_Tobs_bin+T_log10*I_tobs/Num_R1)-ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        end if

        call GRKT4(T,H,Y,M,EPS,D,B,C,G,E,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Epsilon_b,p,z,f_e,R0,index_dyn)
        R_Tobs(I_tobs)=T*(one+z)
        R_Gamma(I_tobs)=Y(1)
        if (index_dyn == 3) then
            R_m(I_tobs)=Y(2)/Para_m_p
        else
            R_m(I_tobs)=Y(2)
        end if
        R(I_tobs)=Y(4)
    end do

    return
end subroutine dynamics_forward

!**********************[Dynamic]**********************
SUBROUTINE F(T, Y, M, D, E_e, E_iso, Eta_0, dNe_ISM, A_star, E_b, p, z, f_e, R0,index_dyn)
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    DIMENSION Y(M),D(M)

    if (A_star > 0.0) then
        dNe_wind=A_star*3.0d35/Y(4)**2
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
    end if
    
    if (Y(4)<R0) then
        dNe=A_star*3.0d35/R0**2
    end if

    Epe=E_e
    qq=-0.2
    t_1=5000.0/(one+z)
    t_2=9000.0/(one+z)
    if (T >= t_1 .and. T<= t_2) then
        A=zero
    else
        A=zero
    end if
!    A=zero
    dB=0.39d0*dsqrt((E_b*dNe)*(Y(1)**2-one))
    gam_c=7.739d8/(dB**2*Y(1)*T)
    gam_m=Epe/f_e*1836d0*(p-two)*(Y(1)-one)/(p-one)+one
    if ((gam_c-gam_m) > 0.001) Epe=Epe*(gam_m/gam_c)**(p-two)
    Bgam=dsqrt(one-one/Y(1)**2)
    
    select case(index_dyn)
    
    case(1)
    ! **************************Huang's Dynamic model******************************
    DNe0=E_iso/(Eta_0-one)/1.5*1.0D3/(4d0*pi)
    D00=Bgam*Para_c/(one-Bgam)
    D01=Y(1)**2-one
    D02=DNe0+Epe*Y(2)+2.0*(one-Epe)*Y(1)*Y(2)
    dM=dNe*Y(4)**2*D00
    D(1)=(-D01*dM+A/(4d0*pi)*(one+T/T_1)**qq/1.5*1.0D3)/D02
    D(2)=dM
    D(3)=(Y(1)-one)*dM
    D(4)=D00
    
    case(2)
    !! ************************Pe'er's Dynamic model*******************************
    ! *******************Caculate hat_gamma derived by Pe'er **********************
    four_v=Bgam*Y(1)
    theta=four_v/3d0*(four_v+1.07d0*four_v**2)/(one+four_v+1.07d0*four_v**2)
    cz=theta/(0.24d0+theta)
    hat_gam=(5d0-1.21937d0*cz+0.18203d0*cz**2-0.96583d0*cz**3+  &
            2.32513d0*cz**4-2.39332d0*cz**5+1.07136d0*cz**6)/3d0
    !!! *************************new equations by Pe'er****************************
    DNe0=E_iso/(Eta_0-one)/1.5*1.0D3/(4d0*pi)
    D00=Bgam*Para_c/(one-Bgam)
    D01=Y(1)**2-one
    D011=hat_gam*D01-(hat_gam-one)*Y(1)*Bgam**2
    D02=DNe0+Epe*Y(2)+2.0*(one-Epe)*Y(1)*Y(2)
    D022=DNe0+Epe*Y(2)+(one-Epe)*Y(2)*(2.0*hat_gam*Y(1)-(hat_gam-1.0)*(one+Y(1)**(-2)))
    dM=DNe*Y(4)**2*D00
    D(1)=-(D011/D022)*dM
    D(2)=dM
    D(3)=(Y(1)-one)*dM
    D(4)=D00
    
    case(3)
    !! ************************Zhang's Dynamic model*******************************
    ! *******************Caculate hat_gamma derived by Pe'er **********************
    four_v=Bgam*Y(1)
    theta=four_v/3d0*(four_v+1.07d0*four_v**2)/(one+four_v+1.07d0*four_v**2)
    cz=theta/(0.24d0+theta)
    hat_gam=(5d0-1.21937d0*cz+0.18203d0*cz**2-0.96583d0*cz**3+  &
            2.32513d0*cz**4-2.39332d0*cz**5+1.07136d0*cz**6)/3d0
    ! *************************new equations by Zhang******************************
    DNe0=E_iso/(Eta_0-one)/Para_c**2
    D00=Bgam*Para_c/(one-Bgam)
    D01=Y(1)**2-one
    dM=4d0*pi*dNe*Y(4)**2*D00*Para_m_p
    D011=-dM/D00*Para_c**2*Y(1)*D01*(hat_gam*Y(1)-hat_gam+one)- &
            (hat_gam*D01+one)*Y(1)*(one-hat_gam)*3d0*Y(3)/Y(4)+Y(1)**2*A*(T/t_1)**qq/D00
    D022=Y(1)**2*(DNe0+Y(2))*Para_c**2+(hat_gam**2*D01+3d0*hat_gam-two)*Y(3)

    D(1)=D011/D022*D00
    D(2)=dM
    D(3)=((1d0-Epe)*(Y(1)-one)*dM/D00*Para_c**2-                &
            (hat_gam-one)*(3d0/Y(4)-one/Y(1)*D(1)/D00)*Y(3))*D00+A*(T/t_1)**qq/D00
    D(4)=D00
    
    end select
    
    return
end subroutine F

SUBROUTINE GRKT4(T,H,Y,M,EPS,D,B,C,G,E,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0,f_e,R0,index_dyn)
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
            call F(T,Y,M,D,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0,f_e,R0,index_dyn)
            E=Y
            B=Y
            do K=1,3
                Y=E+A(K)*D
                B=B+A(K+1)*D/3.0
                TT=T+A(K)
                call F(TT,Y,M,D,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0,f_e,R0,index_dyn)
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
end subroutine GRKT4


