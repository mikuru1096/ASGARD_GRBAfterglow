subroutine dynamics_forward(Boundary,n,Num_R,index_dyn,R_Tobs,R_Gamma,R,R_m)
    !$ use omp_lib
    use constants
    use dynamics_common
    implicit none
    integer, intent(in) :: n,Num_R,index_dyn
    integer :: I_tobs, Num_R1, M
    procedure(dynamics_forward_rhs_iface) :: forward_dynamics_rhs
    real(8), intent(in) :: Boundary(n)
    real(8), intent(out) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),R_m(Num_R)
    real(8) :: Eta_0,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star,E_iso,T_log10_duration,f_e
    real(8) :: E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0,T00,EPS,DM_0,R_dec,t_dec,Grid_Tobs_bin,T_log10,T,H
    real(8),dimension(4) :: Y,D,B,C,G,E

    Eta_0=Boundary(1); Epsilon_e=Boundary(5); Epsilon_b=Boundary(6); p=Boundary(7); z=Boundary(8)
    dNe_ISM=Boundary(11); A_star=Boundary(12); E_iso=Boundary(14); T_log10_duration=Boundary(15); f_e=Boundary(16)
    E_inj_t1=Boundary(17); E_inj_t2=Boundary(18); E_inj=Boundary(19); E_inj_q=Boundary(20)
    R_tr=Boundary(21); f_jump=Boundary(22); f_wide=Boundary(23)
    call dynamics_boundary_r0(Boundary,n,R0)
    call dynamics_set_density_jump_profile(Boundary,n)

    Num_R1=Num_R-1
    Y=[Eta_0-0.001d0,Boundary(2),Boundary(3),Boundary(4)]
    T00=Y(4)*(one/dsqrt(one-one/Eta_0**2)-one)/Para_c
    M=4; EPS=1.0D-5

    DM_0=E_iso/((Eta_0-one)*4d0*pi*Para_m_p_E)
    call dynamics_deceleration_radius(A_star,dNe_ISM,Eta_0,DM_0,R_dec)
    t_dec=R_dec/(two*Eta_0*Eta_0*Para_c)
    Grid_Tobs_bin=min(-5d0,dlog10(t_dec*0.1d0))
    T_log10=T_log10_duration-Grid_Tobs_bin

    do I_tobs=1,Num_R
        call dynamics_log_time_step(zero,Grid_Tobs_bin,T_log10,Num_R1,I_tobs,T,H)
        call dynamics_rk4_forward(forward_dynamics_rhs,T,H,Y,M,EPS,D,B,C,G,E,Epsilon_e,E_iso,Eta_0, &
                                  dNe_ISM,A_star,Epsilon_b,p,z,f_e, &
                                  E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0,index_dyn)
        R_Tobs(I_tobs)=T*(one+z); R_Gamma(I_tobs)=Y(1); R_m(I_tobs)=Y(2)/Para_m_p; R(I_tobs)=Y(4)
    end do
end subroutine dynamics_forward

subroutine forward_dynamics_rhs(T,Y,M,D,E_e,E_iso,Eta_0,dNe_ISM,A_star,E_b,p,z,f_e, &
                                E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0,index_dyn)
    use constants
    use dynamics_common
    implicit none
    integer, intent(in) :: M,index_dyn
    real(8), intent(in) :: T,E_e,E_iso,Eta_0,dNe_ISM,A_star,E_b,p,z,f_e,E_inj_t1,E_inj_t2,E_inj,E_inj_q,R_tr,f_jump,f_wide,R0
    real(8), intent(in) :: Y(M)
    real(8), intent(out) :: D(M)
    real(8), parameter :: forward_synch_b_coeff=0.39d0, forward_gamma_c_coeff=7.739d8
    real(8) :: dNe,t_1,t_2,A,Epe,dB,gam_c,gam_m,Bgam,D00,D01,dM,hat_gam,four_v,theta,cz,DNe0,D011,D02,D022

    call dynamics_external_density_profile(A_star,dNe_ISM,Y(4),R0,1,R_tr,f_jump,f_wide,dNe)
    t_1=E_inj_t1/(one+z); t_2=E_inj_t2/(one+z)
    if (T >= t_1 .and. T <= t_2) then
        A=E_inj*(T/t_1)**E_inj_q
    else
        A=zero
    end if

    Epe=E_e
    dB=forward_synch_b_coeff*dsqrt((E_b*dNe)*(Y(1)**2-one))
    gam_c=forward_gamma_c_coeff/(dB**2*Y(1)*T)
    gam_m=Epe/f_e*Para_m_p_div_m_e*(p-two)*(Y(1)-one)/(p-one)+one
    if ((gam_c-gam_m) > 0.001d0 .and. p >= 2d0) Epe=Epe*(gam_m/gam_c)**(p-two)

    Bgam=dsqrt(one-one/Y(1)**2); D00=Bgam*Para_c/(one-Bgam); D01=Y(1)**2-one
    if (index_dyn == 3) then
        dM=4d0*pi*dNe*Y(4)**2*D00*Para_m_p
    else
        dM=dNe*Y(4)**2*D00
    end if
    if (index_dyn == 2 .or. index_dyn == 3) then
        four_v=Bgam*Y(1)
        theta=four_v/3d0*(four_v+1.07d0*four_v**2)/(one+four_v+1.07d0*four_v**2)
        cz=theta/(0.24d0+theta)
        hat_gam=(5d0-1.21937d0*cz+0.18203d0*cz**2-0.96583d0*cz**3+2.32513d0*cz**4-2.39332d0*cz**5+1.07136d0*cz**6)/3d0
    end if
    D(2)=dM; D(4)=D00

    select case(index_dyn)
    case(1)
        DNe0=E_iso/(Eta_0-one)/1.5d0*1.0D3/(4d0*pi)
        D02=DNe0+Epe*Y(2)+2d0*(one-Epe)*Y(1)*Y(2)
        D(1)=(-D01*dM+A/(4d0*pi)*(one+T/t_1)**E_inj_q/1.5d0*1.0D3)/D02
        D(3)=(Y(1)-one)*dM
    case(2)
        DNe0=E_iso/(Eta_0-one)/1.5d0*1.0D3/(4d0*pi)
        D011=hat_gam*D01-(hat_gam-one)*Y(1)*Bgam**2
        D022=DNe0+Epe*Y(2)+(one-Epe)*Y(2)*(2d0*hat_gam*Y(1)-(hat_gam-1d0)*(one+Y(1)**(-2)))
        D(1)=-(D011/D022)*dM
        D(3)=(Y(1)-one)*dM
    case(3)
        DNe0=E_iso/(Eta_0-one)/Para_c**2
        D011=-dM/D00*Para_c**2*Y(1)*D01*(hat_gam*Y(1)-hat_gam+one)-(hat_gam*D01+one)*Y(1)*(one-hat_gam)*3d0*Y(3)/Y(4)+Y(1)**2*A/D00
        D022=Y(1)**2*(DNe0+Y(2))*Para_c**2+(hat_gam**2*D01+3d0*hat_gam-two)*Y(3)
        D(1)=D011/D022*D00
        D(3)=((1d0-Epe)*(Y(1)-one)*dM/D00*Para_c**2-(hat_gam-one)*(3d0/Y(4)-one/Y(1)*D(1)/D00)*Y(3))*D00+A/D00
    end select
end subroutine forward_dynamics_rhs
