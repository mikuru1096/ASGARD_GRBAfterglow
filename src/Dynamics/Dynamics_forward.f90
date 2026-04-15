subroutine dynamics_forward(Boundary,n,Num_R,index_dyn, R_Tobs,R_Gamma,R,R_m)
    !$ use omp_lib
    use constants
    use dynamics_common
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_R,index_dyn
    real(8), intent(in) :: Boundary(n)
    real(8), intent(out) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),R_m(Num_R)

    real(8),dimension (4) :: Y,D,B,C,G,E

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
    E_inj_t1=Boundary(17)
    E_inj_t2=Boundary(18)
    E_inj=Boundary(19)
    E_inj_q=Boundary(20)
    R_tr=Boundary(21)
    f_jump=Boundary(22)
    f_wide=Boundary(23)
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
    call dynamics_deceleration_radius(A_star,dNe_ISM,Eta_0,DM_0,R_dec)
    t_dec=R_dec/(two*Eta_0*Eta_0*Para_c)
    Grid_Tobs_bin=min(-5d0,dlog10(t_dec*0.1))
    T_log10=T_log10_duration-Grid_Tobs_bin
    !   log time bin
    do I_tobs=1,Num_R
        T=0.0*T00+ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        if (I_tobs < one) then
            H=ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        else
            H=ten**(Grid_Tobs_bin+T_log10*I_tobs/Num_R1)-ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        end if

        call GRKT4(T,H,Y,M,EPS,D,B,C,G,E,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Epsilon_b,p,z,f_e,&
        E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0,index_dyn)
        R_Tobs(I_tobs)=T*(one+z)
        R_Gamma(I_tobs)=Y(1)
        R_m(I_tobs)=Y(2)/Para_m_p
        R(I_tobs)=Y(4)
    end do

    return
end subroutine dynamics_forward

subroutine forward_external_density(A_star,dNe_ISM,RR,f_jump,R_tr,f_wide,R0,dNe)
    use constants
    use dynamics_common
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: A_star,dNe_ISM,RR,f_jump,R_tr,f_wide,R0
    real(8), intent(out) :: dNe

    call dynamics_external_density_base(A_star,dNe_ISM,RR,dNe)
    if (A_star <= zero) then
        dNe=dNe_ISM*(1.0+(f_jump-1d0)*exp(-(log10(RR)-log10(R_tr))**2/(2*f_wide*f_wide)))
    end if

    if (RR < R0) dNe=A_star*3.0d35/R0**2
end subroutine forward_external_density

subroutine forward_rk4_cycle(T,Y,H,N,M,D,B,E,A,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0,f_e, &
                             E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0,index_dyn)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: N,M,index_dyn
    real(8), intent(inout) :: T,Y(M)
    real(8), intent(in) :: H,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0,f_e
    real(8), intent(in) :: E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0
    real(8), intent(inout) :: D(M),B(M),E(M)
    real(8), intent(in) :: A(4)

    DT=H
    do J=1,N
        call F(T,Y,M,D,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0,f_e, &
               E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0,index_dyn)
        E=Y
        B=Y
        do K=1,3
            Y=E+A(K)*D
            B=B+A(K+1)*D/3.0
            TT=T+A(K)
            call F(TT,Y,M,D,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0,f_e, &
                   E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0,index_dyn)
        end do
        Y=B+H*D/6.0
        T=T+DT
    end do
end subroutine forward_rk4_cycle

subroutine forward_energy_injection(T,E_inj_t1,E_inj_t2,E_inj,E_inj_q,z,t_1,t_2,A)
    use constants
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: T,E_inj_t1,E_inj_t2,E_inj,E_inj_q,z
    real(8), intent(out) :: t_1,t_2,A

    t_1=E_inj_t1/(one+z)
    t_2=E_inj_t2/(one+z)
    if (T >= t_1 .and. T <= t_2) then
        A=E_inj*(T/t_1)**E_inj_q
    else
        A=zero
    end if
end subroutine forward_energy_injection

subroutine forward_hat_gamma(Bgam,gamma,hat_gam)
    use constants
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: Bgam,gamma
    real(8), intent(out) :: hat_gam

    four_v=Bgam*gamma
    theta=four_v/3d0*(four_v+1.07d0*four_v**2)/(one+four_v+1.07d0*four_v**2)
    cz=theta/(0.24d0+theta)
    hat_gam=(5d0-1.21937d0*cz+0.18203d0*cz**2-0.96583d0*cz**3+  &
            2.32513d0*cz**4-2.39332d0*cz**5+1.07136d0*cz**6)/3d0
end subroutine forward_hat_gamma

subroutine forward_shock_kinematics(gamma,bgam,D00,D01)
    use constants
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: gamma,bgam
    real(8), intent(out) :: D00,D01

    D00=Bgam*Para_c/(one-Bgam)
    D01=gamma**2-one
end subroutine forward_shock_kinematics

subroutine forward_mass_step(index_dyn,dNe,RR,D00,dM)
    use constants
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: index_dyn
    real(8), intent(in) :: dNe,RR,D00
    real(8), intent(out) :: dM

    if (index_dyn == 3) then
        dM=4d0*pi*dNe*RR**2*D00*Para_m_p
    else
        dM=dNe*RR**2*D00
    end if
end subroutine forward_mass_step

subroutine forward_rhs_huang(T,gamma,mass,Epe,E_iso,Eta_0,dM,D01,A,t_1,E_inj_q,D)
    use constants
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: T,gamma,mass,Epe,E_iso,Eta_0,dM,D01,A,t_1,E_inj_q
    real(8), intent(inout) :: D(4)

    DNe0=E_iso/(Eta_0-one)/1.5*1.0D3/(4d0*pi)
    D02=DNe0+Epe*mass+2.0*(one-Epe)*gamma*mass
    D(1)=(-D01*dM+A/(4d0*pi)*(one+T/t_1)**E_inj_q/1.5*1.0D3)/D02
    D(3)=(gamma-one)*dM
end subroutine forward_rhs_huang

subroutine forward_rhs_peer(gamma,mass,Bgam,Epe,E_iso,Eta_0,dM,D01,hat_gam,D)
    use constants
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: gamma,mass,Bgam,Epe,E_iso,Eta_0,dM,D01,hat_gam
    real(8), intent(inout) :: D(4)

    DNe0=E_iso/(Eta_0-one)/1.5*1.0D3/(4d0*pi)
    D011=hat_gam*D01-(hat_gam-one)*gamma*Bgam**2
    D02=DNe0+Epe*mass+2.0*(one-Epe)*gamma*mass
    D022=DNe0+Epe*mass+(one-Epe)*mass*(2.0*hat_gam*gamma-(hat_gam-1.0)*(one+gamma**(-2)))
    D(1)=-(D011/D022)*dM
    D(3)=(gamma-one)*dM
end subroutine forward_rhs_peer

subroutine forward_rhs_zhang(gamma,mass,energy,radius,D00,Epe,E_iso,Eta_0,dM,D01,hat_gam,A,D)
    use constants
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: gamma,mass,energy,radius,D00,Epe,E_iso,Eta_0,dM,D01,hat_gam,A
    real(8), intent(inout) :: D(4)

    DNe0=E_iso/(Eta_0-one)/Para_c**2
    D011=-dM/D00*Para_c**2*gamma*D01*(hat_gam*gamma-hat_gam+one)- &
            (hat_gam*D01+one)*gamma*(one-hat_gam)*3d0*energy/radius+gamma**2*A/D00
    D022=gamma**2*(DNe0+mass)*Para_c**2+(hat_gam**2*D01+3d0*hat_gam-two)*energy

    D(1)=D011/D022*D00
    D(3)=((1d0-Epe)*(gamma-one)*dM/D00*Para_c**2-                &
            (hat_gam-one)*(3d0/radius-one/gamma*D(1)/D00)*energy)*D00+A/D00
end subroutine forward_rhs_zhang

!**********************[Dynamic]**********************
SUBROUTINE F(T, Y, M, D, E_e, E_iso, Eta_0, dNe_ISM, A_star, E_b, p, z, f_e, &
             E_inj_t1, E_inj_t2, E_inj, E_inj_q, R_tr, f_jump, f_wide, R0, index_dyn)
    use constants
    use dynamics_common
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: M,index_dyn
    real(8), intent(in) :: T,E_e,E_iso,Eta_0,dNe_ISM,A_star,E_b,p,z,f_e
    real(8), intent(in) :: E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0
    real(8), intent(in) :: Y(M)
    real(8), intent(out) :: D(M)

    call forward_external_density(A_star,dNe_ISM,Y(4),f_jump,R_tr,f_wide,R0,dNe)
    call forward_energy_injection(T,E_inj_t1,E_inj_t2,E_inj,E_inj_q,z,t_1,t_2,A)
    
    Epe=E_e
    dB=0.39d0*dsqrt((E_b*dNe)*(Y(1)**2-one))
    gam_c=7.739d8/(dB**2*Y(1)*T)
    gam_m=Epe/f_e*1836d0*(p-two)*(Y(1)-one)/(p-one)+one
    if ((gam_c-gam_m) > 0.001 .and. p>=2) Epe=Epe*(gam_m/gam_c)**(p-two)
    Bgam=dsqrt(one-one/Y(1)**2)
    call forward_shock_kinematics(Y(1),Bgam,D00,D01)
    call forward_mass_step(index_dyn,dNe,Y(4),D00,dM)
    if (index_dyn == 2 .or. index_dyn == 3) then
        call forward_hat_gamma(Bgam,Y(1),hat_gam)
    end if

    D(2)=dM
    D(4)=D00

    select case(index_dyn)
    
    case(1)
        call forward_rhs_huang(T,Y(1),Y(2),Epe,E_iso,Eta_0,dM,D01,A,t_1,E_inj_q,D)
    
    case(2)
        call forward_rhs_peer(Y(1),Y(2),Bgam,Epe,E_iso,Eta_0,dM,D01,hat_gam,D)
    
    case(3)
        call forward_rhs_zhang(Y(1),Y(2),Y(3),Y(4),D00,Epe,E_iso,Eta_0,dM,D01,hat_gam,A,D)
    end select
    
    return
end subroutine F

SUBROUTINE GRKT4(T,H,Y,M,EPS,D,B,C,G,E,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0,f_e,&
E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0,index_dyn)
    use dynamics_common
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: M,index_dyn
    real(8), intent(inout) :: T,Y(M),D(M),B(M),C(M),G(M),E(M)
    real(8), intent(in) :: H,EPS,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0,f_e
    real(8), intent(in) :: E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0
    real(8) :: A(4)

    HH=H
    N=1
    P=1+EPS
    X=T
    C=Y
    do while (P >= EPS)
        call dynamics_rk4_coefficients(HH,A)
        G=Y
        Y=C
        T=X
        call forward_rk4_cycle(T,Y,HH,N,M,D,B,E,A,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0,f_e, &
        E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0,index_dyn)
        call dynamics_rk4_error_n(Y,G,M,P)
        HH=HH/2.0
        N=N+N
    end do
    T=X

    return
end subroutine GRKT4
