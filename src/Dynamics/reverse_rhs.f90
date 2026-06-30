! Reverse-shock dynamical RHS, including MHD jump state and secondary branch derivatives.
subroutine reverse_dynamics_rhs(dB3,T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross, &
             T,Y,D,M,para_m_ej,V3_scale,Delta_0,eta_0,A_star,dNe_ISM,R_tr,f_jump,f_wide,R0, &
             Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r,sigma_r)
    use constants
    use dynamics_common, only: dynamics_external_density_profile, &
                               rs_mag_specific_internal, rs_vegas_ud, reverse_rhs_phase, &
                               active_density_jump_count, active_density_jump_r, active_density_jump_factor, &
                               active_density_jump_width
    use reverse_jump_conditions, only: secondary_reverse_contact_rh
    implicit none
    integer, intent(in) :: M
    real(8), intent(inout) :: dB3,T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross
    real(8), intent(in) :: T,para_m_ej,V3_scale,Delta_0,eta_0,A_star,dNe_ISM,R_tr,f_jump,f_wide,R0
    real(8), intent(in) :: Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r,sigma_r
    real(8), intent(in) :: Y(M)
    real(8), intent(out) :: D(M)
    real(8), parameter :: reverse_synch_b_coeff=0.39d0, reverse_gamma_c_precise_coeff=7.739d8
    real(8) :: gam2,RR,para_m2,para_m3,U3,V3,dNe,u2,u4,Delta,para1,para_n4,beta4,beta2,gam34
    real(8) :: para_n3,betars,delta_beta_rs,u3s,u4s_shock,gamma4s
    real(8) :: dB2,gam_c2,gam_m2,eps2,e3,gam_c3,gam_m3,eps3,dgam2_1,dgam2_2,dgam2,dR,dm2,dm3
    real(8) :: thermal_specific3,thermal_response3,thermal_gamma3,ad3,dV3_exp,dV3_shock,dU3_shock,dU3_ad,dU3,dV3
    real(8) :: secondary_m_total,secondary_u_total,secondary_v_total,secondary_p_total,secondary_inertia_mass
    real(8) :: comp_ratio,rho4,B4_ordered,B3_ordered,sigma_inertia,para_n4_inertia
    real(8) :: magnetic_inertia_mass,magnetic_pressure,magnetic_energy,waiting_shell_inertia_mass,smooth_shell_fraction
    real(8) :: beta_fast,fast_wave_depth,gamma_sq_minus_one
    integer :: j_inertia,m_idx_inertia,u_idx_inertia,v_idx_inertia
    logical :: pre_crossing, waiting_reverse

    waiting_reverse=(reverse_rhs_phase == -1)
    gam2=Y(1); RR=Y(2); para_m2=Y(3); para_m3=Y(4)*para_m_ej
    U3=Y(5)*para_m_ej*para_c**2; V3=Y(6)*V3_scale
    call dynamics_external_density_profile(A_star,dNe_ISM,RR,R0,1,R_tr,f_jump,f_wide,dNe)
    u2=dsqrt(gam2*gam2-one); u4=dsqrt(eta_0*eta_0-one)
    beta4=u4/eta_0; beta2=u2/gam2
    Delta=max(Delta_0,RR/(eta_0*eta_0))
    para1=4d0*pi*Para_m_p*RR*RR; para_n4=para_m_ej/(para1*eta_0*Delta)
    if (sigma_r <= zero) then
        sigma_inertia=one
        smooth_shell_fraction=one
    else
        sigma_inertia=one+sigma_r
        beta_fast=dsqrt(sigma_r/(one+sigma_r))
        fast_wave_depth=RR*beta_fast/(eta_0*eta_0*(one-beta4*beta_fast))/beta4
        if (fast_wave_depth >= Delta) then
            smooth_shell_fraction=one
        else
            smooth_shell_fraction=fast_wave_depth/Delta
        end if
    end if
    para_n4_inertia=para_n4*sigma_inertia
    waiting_shell_inertia_mass=para_m_ej*sigma_inertia*eta_0/(eta_0-one)*smooth_shell_fraction
    gam34=(gam2*gam2+eta_0*eta_0-one)/(eta_0*gam2+u2*u4)
    if (waiting_reverse) then
        comp_ratio=one
        para_n3=para_n4
        delta_beta_rs=zero
        betars=beta4
    else
        u3s=rs_vegas_ud(gam34,sigma_r)
        gamma_sq_minus_one=(gam34-one)*(gam34+one)
        u4s_shock=dsqrt((one+u3s*u3s)*gamma_sq_minus_one)+u3s*gam34
        gamma4s=dsqrt(one+u4s_shock*u4s_shock)
        comp_ratio=u4s_shock/u3s
        para_n3=comp_ratio*para_n4
        delta_beta_rs=u4s_shock/(eta_0*(eta_0*gamma4s-u4*u4s_shock))
        betars=beta4-delta_beta_rs
    end if
    D=zero
    pre_crossing=(reverse_rhs_phase == 1 .or. (reverse_rhs_phase == 0 .and. para_m_ej > para_m3))

    dB2=reverse_synch_b_coeff*dsqrt((Epsilon_b*dNe)*(gam2*gam2-one))
    gam_c2=reverse_gamma_c_precise_coeff/(dB2*dB2*gam2*T)
    gam_m2=Epsilon_e/f_e*Para_m_p_div_m_e*(p_f-two)*(gam2-one)/(p_f-one)+one
    eps2=Epsilon_e*min(one,(gam_m2/gam_c2)**(p_f-two))
    e3=U3/V3
    if (waiting_reverse) then
        B3_ordered=zero
        dB3=zero
    else
        B3_ordered=zero
        if (pre_crossing) then
            rho4=para_n4*Para_m_p
            if (sigma_r <= zero) then
                B4_ordered=zero
            else
                B4_ordered=dsqrt(4d0*pi*para_c*para_c*sigma_r*rho4)
            end if
            B3_ordered=B4_ordered*comp_ratio
        else
            if (T_cross < zero .and. gam2 > one) then
                rho4=para_n4*Para_m_p
                if (sigma_r <= zero) then
                    B4_ordered=zero
                else
                    B4_ordered=dsqrt(4d0*pi*para_c*para_c*sigma_r*rho4)
                end if
                B3_ordered=B4_ordered*comp_ratio
            else
                B3_ordered=B3_ordered_cross*V3_cross/V3*RR/R_cross
            end if
        end if
        dB3=dsqrt(8d0*pi*b_r*e3+B3_ordered*B3_ordered)
    end if
    if (waiting_reverse) then
        thermal_specific3=zero
        thermal_response3=zero
    else if (pre_crossing) then
        thermal_specific3=rs_mag_specific_internal(gam34,sigma_r)
        if (sigma_r <= zero) then
            thermal_response3=one
        else
            thermal_response3=thermal_specific3/(gam34-one)
        end if
    else
        thermal_specific3=U3/(para_m3*para_c**2)
        thermal_response3=zero
    end if
    thermal_gamma3=one+thermal_specific3
    if (waiting_reverse) then
        gam_c3=zero
        gam_m3=one
        eps3=zero
    else
        gam_c3=reverse_gamma_c_precise_coeff/(dB3*dB3*gam2*T)
        gam_m3=e_r/f_e_r*Para_m_p_div_m_e*(p_r-two)*thermal_specific3/(p_r-one)+one
        eps3=e_r*min(one,(gam_m3/gam_c3)**(p_r-two))
    end if
    secondary_m_total=zero; secondary_u_total=zero; secondary_v_total=zero
    secondary_p_total=zero; secondary_inertia_mass=zero
    do j_inertia=1,active_density_jump_count
        m_idx_inertia=6+j_inertia
        u_idx_inertia=6+active_density_jump_count+j_inertia
        v_idx_inertia=6+2*active_density_jump_count+j_inertia
        if (v_idx_inertia > M) error stop 'secondary RS branch state exceeds reverse dynamics vector'
        secondary_m_total=secondary_m_total+Y(m_idx_inertia)*para_m_ej
        secondary_u_total=secondary_u_total+Y(u_idx_inertia)*para_m_ej*para_c**2
        secondary_v_total=secondary_v_total+Y(v_idx_inertia)*V3_scale
    end do
    if (secondary_v_total > zero) then
        secondary_p_total=secondary_u_total/(3d0*secondary_v_total)
        secondary_inertia_mass=(secondary_m_total*para_c**2+secondary_u_total+ &
                                secondary_p_total*secondary_v_total)/para_c**2
        secondary_inertia_mass=secondary_inertia_mass+secondary_p_total*secondary_v_total/(gam2*gam2*para_c**2)
    end if
    magnetic_pressure=B3_ordered*B3_ordered/(8d0*pi)
    magnetic_energy=magnetic_pressure*V3
    magnetic_inertia_mass=(magnetic_energy+magnetic_pressure*V3)/(para_c**2)
    magnetic_inertia_mass=magnetic_inertia_mass+magnetic_pressure*V3/(gam2*gam2*para_c**2)

    dgam2_1=-para1*((gam2*gam2-one)*dNe+(gam2*gam34-eta_0)*delta_beta_rs*eta_0*para_n4_inertia)
    dgam2_2=(para_m2+para_m3+magnetic_inertia_mass+(one-eps2)*(two*gam2-one)*para_m2+ &
             (one-eps3)*thermal_specific3*para_m3 + &
             (one-eps3)*gam2*thermal_response3*para_m3* &
             (eta_0*(one-beta4*beta2)-eta_0*beta4/(gam2**2*beta2)) + secondary_inertia_mass)
    dgam2=dgam2_1/dgam2_2
    dR=beta2/(one-beta2)*para_c; dm2=para1*dNe*dR

    if (waiting_reverse) then
        dm3=zero
        dgam2_1=-u2**2*dm2/dR
        dgam2_2=waiting_shell_inertia_mass+(eps2+two*(one-eps2)*gam2)*para_m2
        dgam2=dgam2_1/dgam2_2
    else if (pre_crossing) then
        dm3=para1*delta_beta_rs*eta_0*para_n4*dR
    else
        dm3=zero
        if (T_cross < zero .and. gam2 > one) then
            T_cross=T; R_cross=RR; gam20=gam2; e3_cross=e3
            U3_cross=U3; V3_cross=V3; M3_cross=para_m3; gam_m_cross=gam_m3
            B3_ordered_cross=B3_ordered
        end if
        dgam2_1=-u2**2*dm2/dR
        dgam2_2=para_m3+magnetic_inertia_mass+secondary_inertia_mass+(eps2+two*(one-eps2)*gam2)*para_m2
        dgam2=dgam2_1/dgam2_2
    end if
    dgam2=dgam2*dR
    ad3=(4d0*thermal_gamma3+one)/(3d0*thermal_gamma3)
    if (waiting_reverse) then
        dV3_exp=zero
    else
        dV3_exp=V3*(3d0*dR/RR-dgam2/gam2)
    end if
    if (pre_crossing) then
        dV3_shock=dm3/(para_n3*Para_m_p)
        dU3_shock=thermal_specific3*dm3*para_c**2
    else
        dV3_shock=zero; dU3_shock=zero
    end if
    if (waiting_reverse) then
        dU3_ad=zero
    else
        dU3_ad=-(ad3-one)*U3*dV3_exp/V3
    end if
    dU3=dU3_shock+dU3_ad; dV3=dV3_shock+dV3_exp
    D(1:6)=[dgam2,dR,dm2,dm3/para_m_ej,dU3/(para_m_ej*para_c**2),dV3/V3_scale]
    call compute_secondary_branch_derivatives()

contains

    logical function secondary_parent_upstream_available(jump_index,parent_mass,parent_energy,parent_volume)
    implicit none
    integer, intent(in) :: jump_index
    real(8), intent(in) :: parent_mass,parent_energy,parent_volume
        secondary_parent_upstream_available=.false.
        if (jump_index <= 1) return
        if (parent_mass <= zero .or. parent_energy <= zero .or. parent_volume <= zero) return
        secondary_parent_upstream_available=.true.
    end function secondary_parent_upstream_available

    subroutine compute_secondary_branch_derivatives()
    implicit none
    integer :: j, m_idx, u_idx, v_idx, e_idx, g_idx, parent_m_idx, parent_u_idx, parent_v_idx
    real(8) :: density_factor,branch_weight,n1,n_excess,n_pre,n4,e4,p4,p3,e3_sec,e_ad,comp,beta_s
    real(8) :: gamma_c_j,gamma43_j,beta_c_j,n3,dm_dR,dM_shock,dV_shock,dU_shock
    real(8) :: u_sec,v_sec,dV_exp_sec,dU_ad_sec,gamma_m_sec,b_i,gam_e_max

        do j=1,active_density_jump_count
            m_idx=6+j
            u_idx=6+active_density_jump_count+j
            v_idx=6+2*active_density_jump_count+j
            e_idx=6+3*active_density_jump_count+j
            g_idx=6+4*active_density_jump_count+j
            if (g_idx > M) error stop 'secondary RS branch state exceeds reverse dynamics vector'
            u_sec=Y(u_idx)*para_m_ej*para_c**2
            v_sec=Y(v_idx)*V3_scale
            D(e_idx)=zero; D(g_idx)=zero
            dV_exp_sec=zero; dU_ad_sec=zero
            if (v_sec > zero) then
                dV_exp_sec=v_sec*(3d0*dR/RR-dgam2/gam2)
                dU_ad_sec=-(one/3d0)*u_sec*dV_exp_sec/v_sec
            end if
            call secondary_reverse_density_branch_rhs(RR,j,density_factor,branch_weight)
            if (branch_weight <= zero .or. gam2 <= one) then
                D(m_idx)=zero
                D(u_idx)=dU_ad_sec/(para_m_ej*para_c**2)
                D(v_idx)=dV_exp_sec/V3_scale
                cycle
            end if
            n1=dNe_ISM*density_factor
            n_excess=dNe_ISM*branch_weight
            n_pre=n1-n_excess
            if (n_pre <= zero) error stop 'secondary branch RHS found non-positive pre-bump density'
            n4=4d0*gam2*n_pre
            e4=4d0*gam2*(gam2-one)*n_pre*Para_m_p*Para_c**2
            p4=e4/3d0
            if (j > 1) then
                parent_m_idx=6+j-1
                parent_u_idx=6+active_density_jump_count+j-1
                parent_v_idx=6+2*active_density_jump_count+j-1
                if (secondary_parent_upstream_available(j,Y(parent_m_idx),Y(parent_u_idx),Y(parent_v_idx))) then
                    n4=Y(parent_m_idx)*para_m_ej/(Para_m_p*Y(parent_v_idx)*V3_scale)
                    e4=Y(parent_u_idx)*para_m_ej*para_c**2/(Y(parent_v_idx)*V3_scale)
                    p4=e4/3d0
                end if
            end if
            call secondary_reverse_contact_rh(gam2,n1,n4,e4,p4,gamma_c_j,p3,gamma43_j,comp,beta_s)
            if (comp <= zero) then
                D(m_idx)=zero
                D(u_idx)=dU_ad_sec/(para_m_ej*para_c**2)
                D(v_idx)=dV_exp_sec/V3_scale
                cycle
            end if
            e3_sec=3d0*p3
            e_ad=e4*comp**(4d0/3d0)
            if (e3_sec <= e_ad) then
                D(m_idx)=zero
                D(u_idx)=dU_ad_sec/(para_m_ej*para_c**2)
                D(v_idx)=dV_exp_sec/V3_scale
                cycle
            end if
            beta_c_j=dsqrt(one-gamma_c_j**(-2))
            n3=comp*n4
            if (p_r <= two) error stop 'secondary branch RHS requires p_r > 2'
            b_i=dsqrt(8d0*pi*b_r*e3_sec)
            gamma_m_sec=one+e_r/f_e_r*(p_r-two)/(p_r-one)*(e3_sec-e_ad)/(n3*Para_m_e*Para_c**2)
            gam_e_max=3d0*Para_m_energy/dsqrt(8d0*b_i*Para_e**3)
            if (gamma_m_sec >= gam_e_max) error stop 'secondary branch RHS electron injection exceeds maximum'
            dm_dR=4d0*pi*RR*RR*n4*Para_m_p*gam2*(beta2-beta_s)/beta_c_j
            dM_shock=dm_dR*dR
            dV_shock=dM_shock/(n3*Para_m_p)
            dU_shock=(e3_sec-e_ad)*dV_shock
            D(m_idx)=dM_shock/para_m_ej
            D(u_idx)=(dU_shock+dU_ad_sec)/(para_m_ej*para_c**2)
            D(v_idx)=(dV_shock+dV_exp_sec)/V3_scale
            D(e_idx)=dU_shock/(para_m_ej*para_c**2)
            D(g_idx)=gamma_m_sec*dU_shock/(para_m_ej*para_c**2)
        end do
    end subroutine compute_secondary_branch_derivatives

    subroutine secondary_reverse_density_branch_rhs(radius,jump_index,density_factor,branch_weight)
    implicit none
    integer, intent(in) :: jump_index
    integer :: k
    real(8), intent(in) :: radius
    real(8), intent(out) :: density_factor,branch_weight
    real(8) :: x,width,profile,total_density,base_density,enhancement

        call dynamics_external_density_profile(A_star,dNe_ISM,radius,R0,1,R_tr,f_jump,f_wide,total_density)
        enhancement=one; branch_weight=zero
        do k=1,active_density_jump_count
            x=radius-active_density_jump_r(k)
            width=active_density_jump_width(k)*active_density_jump_r(k)
            profile=(active_density_jump_factor(k)-one)*dexp(-(x*x)/(2d0*width*width))
            enhancement=enhancement+profile
        end do
        base_density=total_density/enhancement
        do k=1,active_density_jump_count
            x=radius-active_density_jump_r(k)
            width=active_density_jump_width(k)*active_density_jump_r(k)
            profile=(active_density_jump_factor(k)-one)*dexp(-(x*x)/(2d0*width*width))
            if (k == jump_index .and. x >= -4d0*width .and. x < zero) branch_weight=base_density*profile
        end do
        density_factor=total_density
    end subroutine secondary_reverse_density_branch_rhs
end subroutine reverse_dynamics_rhs
