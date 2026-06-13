subroutine dynamics_reverse(Delta_t,e_r,b_r,p_r,f_e_r,sigma_r,Boundary,n,Num_R, &
                            T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross, &
                            R_Tobs,R_Gamma,R,M2,M3,B3,U3_th,V3_comoving,Gamma34_inst)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_deceleration_radius, dynamics_rk4_reverse, &
                               dynamics_rk4_reverse_pre_m3, &
                               dynamics_reverse_rhs_iface, dynamics_log_time_step, &
                               dynamics_external_density_profile, rs_mag_comp, &
                               dynamics_boundary_r0, dynamics_set_density_jump_profile
    implicit none
    integer, intent(in) :: n,Num_R
    integer :: I_tobs, Num_R1
    procedure(dynamics_reverse_rhs_iface) :: reverse_dynamics_rhs
    real(8), intent(in) :: Boundary(n),Delta_t,e_r,b_r,p_r,f_e_r,sigma_r
    real(8), intent(out) :: T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross
    real(8), intent(out) :: R_Tobs(Num_R),R(Num_R),M2(Num_R),M3(Num_R),B3(Num_R),R_Gamma(Num_R)
    real(8), intent(out) :: U3_th(Num_R),V3_comoving(Num_R),Gamma34_inst(Num_R)
    real(8) :: Eta_0,Epsilon_e,Epsilon_b,p_f,z,dNe_ISM,A_star,E_iso,T_log10_duration,f_e
    real(8) :: R_tr,f_jump,f_wide,R0
    real(8) :: Delta_0,para_m_ej,V3_scale,para_m2,para_m3,DM_0,R_dec,T00,t_dec,Grid_Tobs_bin,T_log10,T,H,dB3
    real(8) :: T_state,T_target
    real(8) :: u2_init,u4_init,Delta_init,para_n4_init,gam34_init,para_n3_init,comp_init
    real(8),allocatable :: Y(:)

    allocate(Y(6))
    Eta_0=Boundary(1); R(1)=Boundary(4); Epsilon_e=Boundary(5); Epsilon_b=Boundary(6); p_f=Boundary(7); z=Boundary(8)
    dNe_ISM=Boundary(11); A_star=Boundary(12); E_iso=Boundary(14); T_log10_duration=Boundary(15); f_e=Boundary(16)
    R_tr=Boundary(21); f_jump=Boundary(22); f_wide=Boundary(23)
    call dynamics_boundary_r0(Boundary,n,R0)
    call dynamics_set_density_jump_profile(Boundary,n)
    Delta_0=Delta_t*para_c; para_m_ej=E_iso/eta_0/para_c**2

    if (A_star > zero) then
        para_m2=4d0*pi*R(1)*A_star*3d35*para_m_p
    else
        para_m2=4d0*pi*R(1)**3*dNe_ISM*para_m_p/3d0
    end if
    para_m3=1d1; R_Gamma(1)=Eta_0-0.001d0
    u2_init=dsqrt(R_Gamma(1)*R_Gamma(1)-one); u4_init=dsqrt(Eta_0*Eta_0-one)
    Delta_init=max(Delta_0,R(1)/Eta_0**2)
    para_n4_init=para_m_ej/(4d0*pi*Para_m_p*R(1)*R(1)*Eta_0*Delta_init)
    gam34_init=(R_Gamma(1)*R_Gamma(1)+Eta_0*Eta_0-one)/(Eta_0*R_Gamma(1)+u2_init*u4_init)
    comp_init=rs_mag_comp(gam34_init,sigma_r)
    para_n3_init=comp_init*para_n4_init
    V3_scale=para_m_ej/(para_n3_init*Para_m_p)
    Y=[R_Gamma(1),R(1),para_m2,para_m3/para_m_ej,(gam34_init-one)*para_m3/para_m_ej, &
       para_m3/(para_n3_init*Para_m_p)/V3_scale]

    DM_0=E_iso/((Eta_0-one)*4d0*pi*Para_m_p_E)
    call dynamics_deceleration_radius(A_star,dNe_ISM,Eta_0,DM_0,R_dec)
    T_cross=-1d0; U3_cross=zero; V3_cross=zero; M3_cross=zero; gam_m_cross=zero
    B3_ordered_cross=zero
    ! RS 状态向量中 Y(2) 是 shock 半径；初始到达时刻必须由半径而不是 swept mass 定义。
    T00=Y(2)*(one/dsqrt(one-one/Eta_0**2)-one)/Para_c
    t_dec=R_dec/(two*Eta_0*Eta_0*Para_c)
    Grid_Tobs_bin=min(log10(T00)-2d0,dlog10(t_dec*0.1d0))
    T_log10=T_log10_duration-Grid_Tobs_bin; Num_R1=Num_R-1
    T_state=T00

    do I_tobs=1,Num_R
        call dynamics_log_time_step(T00,Grid_Tobs_bin,T_log10,Num_R1,I_tobs,T_target,H)
        do while (T_state < T_target)
            if (T_cross < zero .and. Y(4) < one) then
                call dynamics_rk4_reverse_pre_m3(reverse_dynamics_rhs,dB3,T_cross,R_cross,e3_cross,gam20, &
                                                 U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross, &
                                                 T_state,T_target,Y,para_m_ej,V3_scale,Delta_0,eta_0,A_star,dNe_ISM, &
                                                 R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r,sigma_r)
            else
                H=T_target-T_state
                T=T_state
                call dynamics_rk4_reverse(reverse_dynamics_rhs,dB3,T_cross,R_cross,e3_cross,gam20, &
                                          U3_cross,V3_cross,M3_cross, &
                                          gam_m_cross,B3_ordered_cross,T,H,Y,para_m_ej,V3_scale,Delta_0, &
                                          eta_0,A_star,dNe_ISM,R_tr,f_jump,f_wide,R0, &
                                          Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r,sigma_r)
                T_state=T_target
            end if
        end do
        R_Tobs(I_tobs)=T_target*(one+z); R_Gamma(I_tobs)=Y(1); R(I_tobs)=Y(2); M2(I_tobs)=Y(3)
        M3(I_tobs)=Y(4)*para_m_ej; U3_th(I_tobs)=Y(5)*para_m_ej*para_c**2
        V3_comoving(I_tobs)=Y(6)*V3_scale; B3(I_tobs)=dB3
        Gamma34_inst(I_tobs)=(Y(1)*Y(1)+Eta_0*Eta_0-one)/(Eta_0*Y(1)+dsqrt(Y(1)*Y(1)-one)*u4_init)
    end do

    deallocate(Y)
end subroutine dynamics_reverse

subroutine secondary_reverse_contact_rh(gamma4,n1,n4,e4,p4,gamma_c,p3,gamma43,n3_over_n4)
    use constants
    implicit none
    integer :: I
    real(8), intent(in) :: gamma4,n1,n4,e4,p4
    real(8), intent(out) :: gamma_c,p3,gamma43,n3_over_n4
    real(8) :: lo,hi,mid,f_lo,f_hi,f_mid,beta_c,beta4

    if (gamma4 <= one .or. n1 <= zero .or. n4 <= zero .or. e4 <= zero .or. p4 <= zero) &
        error stop 'secondary_reverse_contact_rh requires positive hot upstream state'
    lo=one+1d-10; hi=gamma4*(one-1d-10)
    call pressure_difference(lo,f_lo)
    call pressure_difference(hi,f_hi)
    if (f_lo*f_hi > zero) error stop 'secondary_reverse_contact_rh has no physical bracket'
    do I=1,80
        mid=0.5d0*(lo+hi)
        call pressure_difference(mid,f_mid)
        if (f_lo*f_mid <= zero) then
            hi=mid; f_hi=f_mid
        else
            lo=mid; f_lo=f_mid
        end if
    end do
    gamma_c=0.5d0*(lo+hi)
    beta_c=dsqrt(one-gamma_c**(-2)); beta4=dsqrt(one-gamma4**(-2))
    gamma43=gamma_c*gamma4*(one-beta_c*beta4)
    p3=p4+(4d0/3d0)*(gamma43*gamma43-one)*(e4+p4)
    call solve_reverse_shock_compression(n3_over_n4)

contains

    subroutine pressure_difference(gamma_trial, diff)
    implicit none
    real(8), intent(in) :: gamma_trial
    real(8), intent(out) :: diff
    real(8) :: beta_trial,beta_up,gamma_rel,p2_trial,p3_trial

        beta_trial=dsqrt(one-gamma_trial**(-2)); beta_up=dsqrt(one-gamma4**(-2))
        gamma_rel=gamma_trial*gamma4*(one-beta_trial*beta_up)
        p2_trial=(4d0/3d0)*(gamma_trial*gamma_trial-one)*n1*Para_m_p*Para_c**2
        p3_trial=p4+(4d0/3d0)*(gamma_rel*gamma_rel-one)*(e4+p4)
        diff=p2_trial-p3_trial
    end subroutine pressure_difference

    subroutine solve_reverse_shock_compression(comp)
    implicit none
    real(8), intent(out) :: comp
    integer :: K
    real(8) :: beta_s_lo,beta_s_hi,beta_s_mid,g_lo,g_hi,g_mid,eps_beta,beta_comp_one

        if (beta4 <= beta_c) error stop 'secondary_reverse_contact_rh reverse shock speed bracket is empty'
        eps_beta=1d-14*(beta4-beta_c)
        beta_s_lo=beta_c+eps_beta
        beta_s_hi=beta4-eps_beta
        call solve_compression_unity_speed(beta_s_lo,beta_s_hi,beta_comp_one)
        beta_s_hi=beta_comp_one-eps_beta
        call shock_momentum_difference(beta_s_lo,g_lo,comp)
        call shock_momentum_difference(beta_s_hi,g_hi,comp)
        if (g_lo*g_hi > zero) then
            comp=-one
            return
        end if
        do K=1,80
            beta_s_mid=0.5d0*(beta_s_lo+beta_s_hi)
            call shock_momentum_difference(beta_s_mid,g_mid,comp)
            if (g_lo*g_mid <= zero) then
                beta_s_hi=beta_s_mid; g_hi=g_mid
            else
                beta_s_lo=beta_s_mid; g_lo=g_mid
            end if
        end do
        call shock_momentum_difference(0.5d0*(beta_s_lo+beta_s_hi),g_mid,comp)
    end subroutine solve_reverse_shock_compression

    subroutine solve_compression_unity_speed(beta_lo,beta_hi,beta_one)
    implicit none
    real(8), intent(in) :: beta_lo,beta_hi
    real(8), intent(out) :: beta_one
    integer :: K
    real(8) :: lo_c,hi_c,mid_c,c_lo,c_mid,comp_tmp,diff_tmp

        lo_c=beta_lo; hi_c=beta_hi
        call shock_momentum_difference(lo_c,diff_tmp,comp_tmp)
        c_lo=comp_tmp-one
        do K=1,80
            mid_c=0.5d0*(lo_c+hi_c)
            call shock_momentum_difference(mid_c,diff_tmp,comp_tmp)
            c_mid=comp_tmp-one
            if (c_lo*c_mid <= zero) then
                hi_c=mid_c
            else
                lo_c=mid_c; c_lo=c_mid
            end if
        end do
        beta_one=0.5d0*(lo_c+hi_c)
    end subroutine solve_compression_unity_speed

    subroutine shock_momentum_difference(beta_s,diff,comp)
    implicit none
    real(8), intent(in) :: beta_s
    real(8), intent(out) :: diff,comp
    real(8) :: beta4_s,beta3_s,gamma4_s,gamma3_s,u4_s,u3_s,w4,w3,e3_trial

        beta4_s=(beta4-beta_s)/(one-beta4*beta_s)
        beta3_s=(beta_c-beta_s)/(one-beta_c*beta_s)
        gamma4_s=one/dsqrt(one-beta4_s*beta4_s)
        gamma3_s=one/dsqrt(one-beta3_s*beta3_s)
        u4_s=gamma4_s*dabs(beta4_s)
        u3_s=gamma3_s*dabs(beta3_s)
        comp=u4_s/u3_s
        e3_trial=3d0*p3
        w4=n4*Para_m_p*Para_c**2+e4+p4
        w3=comp*n4*Para_m_p*Para_c**2+e3_trial+p3
        diff=(w4*u4_s*u4_s+p4)-(w3*u3_s*u3_s+p3)
    end subroutine shock_momentum_difference
end subroutine secondary_reverse_contact_rh

subroutine reverse_dynamics_rhs(dB3,T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross, &
             T,Y,D,para_m_ej,V3_scale,Delta_0,eta_0,A_star,dNe_ISM,R_tr,f_jump,f_wide,R0, &
             Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r,sigma_r)
    use constants
    use dynamics_common, only: dynamics_external_density_profile, rs_mag_comp, rs_b4_up, reverse_rhs_phase
    implicit none
    real(8), intent(inout) :: dB3,T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross
    real(8), intent(in) :: T,para_m_ej,V3_scale,Delta_0,eta_0,A_star,dNe_ISM,R_tr,f_jump,f_wide,R0
    real(8), intent(in) :: Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r,sigma_r
    real(8), intent(in) :: Y(6)
    real(8), intent(out) :: D(6)
    real(8), parameter :: reverse_synch_b_coeff=0.39d0, reverse_gamma_c_precise_coeff=7.739d8
    real(8) :: gam2,RR,para_m2,para_m3,U3,V3,dNe,u2,u4,Delta,para1,para_n4,beta4,beta2,gam34,para_n3,betars
    real(8) :: dB2,gam_c2,gam_m2,eps2,e3,gam_c3,gam_m3,eps3,dgam2_1,dgam2_2,dgam2,dR,dm2,dm3
    real(8) :: thermal_gamma3,ad3,dV3_exp,dV3_shock,dU3_shock,dU3_ad,dU3,dV3
    real(8) :: comp_ratio,rho4,B4_ordered,B3_ordered
    logical :: pre_crossing

    call decode_reverse_state()
    pre_crossing=(reverse_rhs_phase == 1 .or. (reverse_rhs_phase == 0 .and. para_m_ej > para_m3))

    call compute_region2_radiative_efficiency()
    e3=U3/V3
    call compute_region3_field()
    call compute_region3_radiative_efficiency()

    dgam2_1=-para1*((gam2*gam2-one)*dNe+(gam2*gam34-eta_0)*(beta4-betars)*eta_0*para_n4)
    dgam2_2=(para_m2+para_m3+(one-eps2)*(two*gam2-one)*para_m2+(one-eps3)*(gam34-one)*para_m3+ &
             (one-eps3)*gam2*para_m3*(eta_0*(one-beta4*beta2)-eta_0*beta4/(gam2**2*beta2)))
    dgam2=dgam2_1/dgam2_2
    dR=beta2/(one-beta2)*para_c; dm2=para1*dNe*dR

    if (pre_crossing) then
        dm3=para1*(beta4-betars)*eta_0*para_n4*dR
        thermal_gamma3=gam34
    else
        dm3=zero
        thermal_gamma3=one+U3/(para_m3*para_c**2)
        if (T_cross < zero .and. gam2 > one) then
            T_cross=T; R_cross=RR; gam20=gam2; e3_cross=e3
            U3_cross=U3; V3_cross=V3; M3_cross=para_m3; gam_m_cross=gam_m3
            B3_ordered_cross=B3_ordered
        end if
        dgam2_1=-u2**2*dm2/dR; dgam2_2=para_m_ej+(eps2+two*(one-eps2)*gam2)*para_m2; dgam2=dgam2_1/dgam2_2
    end if
    dgam2=dgam2*dR
    ad3=(4d0*thermal_gamma3+one)/(3d0*thermal_gamma3)
    dV3_exp=V3*(3d0*dR/RR-dgam2/gam2)
    if (pre_crossing) then
        dV3_shock=dm3/(para_n3*Para_m_p)
        dU3_shock=(gam34-one)*dm3*para_c**2
    else
        dV3_shock=zero; dU3_shock=zero
    end if
    dU3_ad=-(ad3-one)*U3*dV3_exp/V3
    dU3=dU3_shock+dU3_ad; dV3=dV3_shock+dV3_exp
    D=[dgam2,dR,dm2,dm3/para_m_ej,dU3/(para_m_ej*para_c**2),dV3/V3_scale]

contains

    subroutine decode_reverse_state()
    implicit none

        gam2=Y(1); RR=Y(2); para_m2=Y(3); para_m3=Y(4)*para_m_ej
        U3=Y(5)*para_m_ej*para_c**2; V3=Y(6)*V3_scale
        call dynamics_external_density_profile(A_star,dNe_ISM,RR,R0,1,R_tr,f_jump,f_wide,dNe)
        u2=dsqrt(gam2*gam2-one); u4=dsqrt(eta_0*eta_0-one); Delta=max(Delta_0,RR/eta_0**2)
        para1=4d0*pi*Para_m_p*RR*RR; para_n4=para_m_ej/(para1*eta_0*Delta)
        beta4=u4/eta_0; beta2=u2/gam2
        gam34=(gam2*gam2+eta_0*eta_0-one)/(eta_0*gam2+u2*u4)
        comp_ratio=rs_mag_comp(gam34,sigma_r)
        para_n3=comp_ratio*para_n4; betars=(u2*para_n3-u4*para_n4)/(gam2*para_n3-eta_0*para_n4)
    end subroutine decode_reverse_state

    subroutine compute_region2_radiative_efficiency()
    implicit none

        dB2=reverse_synch_b_coeff*dsqrt((Epsilon_b*dNe)*(gam2*gam2-one))
        gam_c2=reverse_gamma_c_precise_coeff/(dB2*dB2*gam2*T)
        gam_m2=Epsilon_e/f_e*Para_m_p_div_m_e*(p_f-two)*(gam2-one)/(p_f-one)+one
        eps2=Epsilon_e*min(one,(gam_m2/gam_c2)**(p_f-two))
    end subroutine compute_region2_radiative_efficiency

    subroutine compute_region3_field()
    implicit none

        B3_ordered=zero
        if (pre_crossing) then
            rho4=para_n4*Para_m_p
            B4_ordered=rs_b4_up(rho4,sigma_r)
            B3_ordered=B4_ordered*comp_ratio
        else
            if (T_cross < zero .and. gam2 > one) then
                rho4=para_n4*Para_m_p
                B4_ordered=rs_b4_up(rho4,sigma_r)
                B3_ordered=B4_ordered*comp_ratio
            else
                B3_ordered=B3_ordered_cross*V3_cross/V3*RR/R_cross
            end if
        end if
        dB3=dsqrt(8d0*pi*b_r*e3+B3_ordered*B3_ordered)
    end subroutine compute_region3_field

    subroutine compute_region3_radiative_efficiency()
    implicit none

        gam_c3=reverse_gamma_c_precise_coeff/(dB3*dB3*gam2*T)
        gam_m3=e_r/f_e_r*Para_m_p_div_m_e*(p_r-two)*(gam34-one)/(p_r-one)+one
        eps3=e_r*min(one,(gam_m3/gam_c3)**(p_r-two))
    end subroutine compute_region3_radiative_efficiency
end subroutine reverse_dynamics_rhs
