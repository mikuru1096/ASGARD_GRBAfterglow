module reverse_rhs_module
    implicit none
    private
    public :: reverse_dynamics_rhs

contains

! Reverse-shock dynamical RHS, including MHD jump state and secondary branch derivatives.
subroutine reverse_dynamics_rhs(phase,rs_state,T,Y,D,M,mej,V3_scale,Delta_0,eta_0,A_star,dNe_ISM,R_tr,f_jump,f_wide,R0, &
             Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,fer,sigma_r)
    use constants
    use dynamics_density_profile, only: density_profile, jump_max, &
                                        secondary_branch_density
    use reverse_shock_mhd_jump, only: rs_mhd_state
    use reverse_shock_state, only: wait_phase, precross_phase, &
                                   rs_db3, rs_tcross, rs_rcross, rs_e3cross, rs_gam20, &
                                   rs_u3cross, rs_v3cross, rs_m3cross, rs_gammcross, rs_b3ordered, rs_nstate
    use reverse_jump_conditions, only: reverse_contact
    implicit none
    integer, intent(in) :: phase,M
    real(8), intent(inout), dimension(rs_nstate) :: rs_state
    real(8), intent(in) :: T,mej,V3_scale,Delta_0,eta_0,A_star,dNe_ISM,R_tr,f_jump,f_wide,R0
    real(8), intent(in) :: Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,fer,sigma_r
    real(8), intent(in), dimension(M) :: Y
    real(8), intent(out), dimension(M) :: D
    real(8), parameter :: rs_bcoeff=0.39d0, cool_coeff=7.739d8
    real(8) :: gam2,RR,para_m2,para_m3,U3,V3,dNe,u2,u4,Delta,para1,para_n4,beta4,beta2,gam34
    real(8) :: para_n3,betars,dbrs,u3s,u4s_shock,gamma4s,jump_internal
    real(8) :: dB2,gam_c2,gam_m2,eps2,e3,gam_c3,gam_m3,eps3,dgam2_1,dgam2_2,dgam2,dR,dm2,dm3
    real(8) :: thermal_specific3,thermal_response3,thermal_gamma3,ad3,dV3_exp,dV3_shock,dU3_shock,dU3_ad,dU3,dV3
    real(8) :: sec_m,sec_u,sec_v,sec_p,sec_inertia
    real(8) :: comp_ratio,rho4,B4_ordered,B3_ordered,sigma_inertia,n4_inertia
    real(8) :: mag_inertia,magnetic_pressure,magnetic_energy,wait_inertia,shell_frac
    real(8) :: beta_fast,fast_depth,gsq1
    real(8) :: base_density
    real(8), dimension(jump_max) :: branch_weight
    integer :: j_inertia,mi_idx,ui_idx,vi_idx,nbranch
    integer, parameter :: idb=rs_db3,itc=rs_tcross,irc=rs_rcross,iec=rs_e3cross,ig20=rs_gam20
    integer, parameter :: iuc=rs_u3cross,ivc=rs_v3cross,imc=rs_m3cross,igm=rs_gammcross,ib3=rs_b3ordered
    logical :: pre_crossing, waiting_reverse, shock_allowed


    waiting_reverse=(phase == wait_phase)
    nbranch=(M-6)/7
    gam2=Y(1); RR=Y(2); para_m2=Y(3); para_m3=Y(4)*mej
    U3=Y(5)*mej*para_c**2; V3=Y(6)*V3_scale
    call density_profile(A_star,dNe_ISM,RR,R0,1,R_tr,f_jump,f_wide,dNe)
    u2=dsqrt(gam2*gam2-1d0); u4=dsqrt(eta_0*eta_0-1d0)
    beta4=u4/eta_0; beta2=u2/gam2
    Delta=max(Delta_0,RR/(eta_0*eta_0))
    para1=4d0*pi*Para_m_p*RR*RR; para_n4=mej/(para1*eta_0*Delta)
    if (sigma_r <= 0d0) then
        sigma_inertia=1d0
        shell_frac=1d0
    else
        sigma_inertia=1d0+sigma_r
        beta_fast=dsqrt(sigma_r/(1d0+sigma_r))
        fast_depth=RR*beta_fast/(eta_0*eta_0*(1d0-beta4*beta_fast))/beta4
        if (fast_depth >= Delta) then
            shell_frac=1d0
        else
            shell_frac=fast_depth/Delta
        end if
    end if
    n4_inertia=para_n4*sigma_inertia
    wait_inertia=mej*sigma_inertia*eta_0/(eta_0-1d0)*shell_frac
    gam34=(gam2*gam2+eta_0*eta_0-1d0)/(eta_0*gam2+u2*u4)
    if (waiting_reverse) then
        comp_ratio=1d0
        para_n3=para_n4
        dbrs=0d0
        betars=beta4
    else
        call rs_mhd_state(gam34,sigma_r,u3s,comp_ratio,jump_internal,shock_allowed)
        gsq1=(gam34-1d0)*(gam34+1d0)
        u4s_shock=dsqrt((1d0+u3s*u3s)*gsq1)+u3s*gam34
        gamma4s=dsqrt(1d0+u4s_shock*u4s_shock)
        para_n3=comp_ratio*para_n4
        dbrs=u4s_shock/(eta_0*(eta_0*gamma4s-u4*u4s_shock))
        betars=beta4-dbrs
    end if
    pre_crossing=(phase == precross_phase)

    dB2=rs_bcoeff*dsqrt((Epsilon_b*dNe)*(gam2*gam2-1d0))
    gam_c2=cool_coeff/(dB2*dB2*gam2*T)
    gam_m2=Epsilon_e/f_e*Para_m_p_DIV_m_e*(p_f-2d0)*(gam2-1d0)/(p_f-1d0)+1d0
    eps2=Epsilon_e*min(1d0,(gam_m2/gam_c2)**(p_f-2d0))
    e3=U3/V3
    if (waiting_reverse) then
        B3_ordered=0d0
        rs_state(idb)=0d0
    else
        B3_ordered=0d0
        if (pre_crossing) then
            rho4=para_n4*Para_m_p
            if (sigma_r <= 0d0) then
                B4_ordered=0d0
            else
                B4_ordered=dsqrt(4d0*pi*para_c*para_c*sigma_r*rho4)
            end if
            B3_ordered=B4_ordered*comp_ratio
        else
            if (rs_state(itc) < 0d0 .and. gam2 > 1d0) then
                rho4=para_n4*Para_m_p
                if (sigma_r <= 0d0) then
                    B4_ordered=0d0
                else
                    B4_ordered=dsqrt(4d0*pi*para_c*para_c*sigma_r*rho4)
                end if
                B3_ordered=B4_ordered*comp_ratio
            else
                B3_ordered=rs_state(ib3)*rs_state(ivc)/V3*RR/rs_state(irc)
            end if
        end if
        rs_state(idb)=dsqrt(8d0*pi*b_r*e3+B3_ordered*B3_ordered)
    end if
    if (waiting_reverse) then
        thermal_specific3=0d0
        thermal_response3=0d0
    else if (pre_crossing) then
        thermal_specific3=jump_internal
        if (sigma_r <= 0d0) then
            thermal_response3=1d0
        else
            thermal_response3=thermal_specific3/(gam34-1d0)
        end if
    else
        thermal_specific3=U3/(para_m3*para_c**2)
        thermal_response3=0d0
    end if
    thermal_gamma3=1d0+thermal_specific3
    if (waiting_reverse) then
        gam_c3=0d0
        gam_m3=1d0
        eps3=0d0
    else
        gam_c3=cool_coeff/(rs_state(idb)*rs_state(idb)*gam2*T)
        gam_m3=e_r/fer*Para_m_p_DIV_m_e*(p_r-2d0)*thermal_specific3/(p_r-1d0)+1d0
        eps3=e_r*min(1d0,(gam_m3/gam_c3)**(p_r-2d0))
    end if
    sec_m=0d0; sec_u=0d0; sec_v=0d0
    sec_p=0d0; sec_inertia=0d0
    do j_inertia=1,nbranch
        mi_idx=6+j_inertia
        ui_idx=6+nbranch+j_inertia
        vi_idx=6+2*nbranch+j_inertia
        sec_m=sec_m+Y(mi_idx)*mej
        sec_u=sec_u+Y(ui_idx)*mej*para_c**2
        sec_v=sec_v+Y(vi_idx)*V3_scale
    end do
    if (sec_v > 0d0) then
        sec_p=sec_u/(3d0*sec_v)
        sec_inertia=(sec_m*para_c**2+sec_u+ &
                                sec_p*sec_v)/para_c**2
        sec_inertia=sec_inertia+sec_p*sec_v/(gam2*gam2*para_c**2)
    end if
    magnetic_pressure=B3_ordered*B3_ordered/(8d0*pi)
    magnetic_energy=magnetic_pressure*V3
    mag_inertia=(magnetic_energy+magnetic_pressure*V3)/(para_c**2)
    mag_inertia=mag_inertia+magnetic_pressure*V3/(gam2*gam2*para_c**2)

    dgam2_1=-para1*((gam2*gam2-1d0)*dNe+(gam2*gam34-eta_0)*dbrs*eta_0*n4_inertia)
    dgam2_2=(para_m2+para_m3+mag_inertia+(1d0-eps2)*(2d0*gam2-1d0)*para_m2+ &
             (1d0-eps3)*thermal_specific3*para_m3 + &
             (1d0-eps3)*gam2*thermal_response3*para_m3* &
             (eta_0*(1d0-beta4*beta2)-eta_0*beta4/(gam2**2*beta2)) + sec_inertia)
    dgam2=dgam2_1/dgam2_2
    dR=beta2/(1d0-beta2)*para_c; dm2=para1*dNe*dR

    if (waiting_reverse) then
        dm3=0d0
        dgam2_1=-u2**2*dm2/dR
        dgam2_2=wait_inertia+(eps2+2d0*(1d0-eps2)*gam2)*para_m2
        dgam2=dgam2_1/dgam2_2
    else if (pre_crossing) then
        dm3=para1*dbrs*eta_0*para_n4*dR
    else
        dm3=0d0
        if (rs_state(itc) < 0d0 .and. gam2 > 1d0) then
            rs_state(itc)=T; rs_state(irc)=RR; rs_state(ig20)=gam2; rs_state(iec)=e3
            rs_state(iuc)=U3; rs_state(ivc)=V3; rs_state(imc)=para_m3; rs_state(igm)=gam_m3
            rs_state(ib3)=B3_ordered
        end if
        dgam2_1=-u2**2*dm2/dR
        dgam2_2=para_m3+mag_inertia+sec_inertia+(eps2+2d0*(1d0-eps2)*gam2)*para_m2
        dgam2=dgam2_1/dgam2_2
    end if
    dgam2=dgam2*dR
    ad3=(4d0*thermal_gamma3+1d0)/(3d0*thermal_gamma3)
    if (waiting_reverse) then
        dV3_exp=0d0
    else
        dV3_exp=V3*(3d0*dR/RR-dgam2/gam2)
    end if
    if (pre_crossing) then
        dV3_shock=dm3/(para_n3*Para_m_p)
        dU3_shock=thermal_specific3*dm3*para_c**2
    else
        dV3_shock=0d0; dU3_shock=0d0
    end if
    if (waiting_reverse) then
        dU3_ad=0d0
    else
        dU3_ad=-(ad3-1d0)*U3*dV3_exp/V3
    end if
    dU3=dU3_shock+dU3_ad; dV3=dV3_shock+dV3_exp
    D(1)=dgam2
    D(2)=dR
    D(3)=dm2
    D(4)=dm3/mej
    D(5)=dU3/(mej*para_c**2)
    D(6)=dV3/V3_scale

    ! Evaluate exact explicit-jump or tabulated-profile branch excess once per RHS.
    do j_inertia=1,nbranch
        call secondary_branch_density(A_star,dNe_ISM,RR,R0,1,R_tr,f_jump,f_wide, &
                                      j_inertia,base_density,branch_weight(j_inertia),dNe)
    end do
    call branch_derivs()

contains

    logical function parent_ready(jump_index,parent_mass,parent_energy,parent_volume)
    implicit none
    integer, intent(in) :: jump_index
    real(8), intent(in) :: parent_mass,parent_energy,parent_volume
        parent_ready=.false.
        if (jump_index <= 1) return
        if (parent_mass <= 0d0 .or. parent_energy <= 0d0 .or. parent_volume <= 0d0) return
        parent_ready=.true.
    end function parent_ready

    subroutine branch_derivs()
    implicit none
    integer :: j, m_idx, u_idx, v_idx, e_idx, g_idx, h_idx, c_idx, pm_idx, pu_idx, pv_idx
    real(8) :: n1,n_excess,n_pre,n4,e4,p4,p3,e3_sec,e_ad,comp,beta_s
    real(8) :: gc_j,gamma43_j,bc_j,n3,dm_dR,dM_shock,dV_shock,dU_shock
    real(8) :: u_sec,v_sec,dv_exp,du_ad,gm_sec,b_i,gemax

        do j=1,nbranch
            m_idx=6+j
            u_idx=6+nbranch+j
            v_idx=6+2*nbranch+j
            e_idx=6+3*nbranch+j
            g_idx=6+4*nbranch+j
            h_idx=6+5*nbranch+j
            c_idx=6+6*nbranch+j
            u_sec=Y(u_idx)*mej*para_c**2
            v_sec=Y(v_idx)*V3_scale
            D(e_idx)=0d0; D(g_idx)=0d0
            D(h_idx)=0d0; D(c_idx)=0d0
            dv_exp=0d0; du_ad=0d0
            if (v_sec > 0d0) then
                dv_exp=v_sec*(3d0*dR/RR-dgam2/gam2)
                du_ad=-(1d0/3d0)*u_sec*dv_exp/v_sec
            end if
            if (branch_weight(j) <= 0d0 .or. gam2 <= 1d0) then
                D(m_idx)=0d0
                D(u_idx)=du_ad/(mej*para_c**2)
                D(v_idx)=dv_exp/V3_scale
                cycle
            end if
            n1=dNe
            n_excess=branch_weight(j)
            n_pre=n1-n_excess
            if (n_pre <= 0d0) error stop 'secondary branch RHS found non-positive pre-bump density'
            n4=4d0*gam2*n_pre
            e4=4d0*gam2*(gam2-1d0)*n_pre*Para_m_p*Para_c**2
            p4=e4/3d0
            if (j > 1) then
                pm_idx=6+j-1
                pu_idx=6+nbranch+j-1
                pv_idx=6+2*nbranch+j-1
                if (parent_ready(j,Y(pm_idx),Y(pu_idx),Y(pv_idx))) then
                    n4=Y(pm_idx)*mej/(Para_m_p*Y(pv_idx)*V3_scale)
                    e4=Y(pu_idx)*mej*para_c**2/(Y(pv_idx)*V3_scale)
                    p4=e4/3d0
                end if
            end if
            call reverse_contact(gam2,n1,n4,e4,p4,gc_j,p3,gamma43_j,comp,beta_s)
            if (comp <= 0d0) then
                D(m_idx)=0d0
                D(u_idx)=du_ad/(mej*para_c**2)
                D(v_idx)=dv_exp/V3_scale
                cycle
            end if
            e3_sec=3d0*p3
            e_ad=e4*comp**(4d0/3d0)
            if (e3_sec <= e_ad) then
                D(m_idx)=0d0
                D(u_idx)=du_ad/(mej*para_c**2)
                D(v_idx)=dv_exp/V3_scale
                cycle
            end if
            bc_j=dsqrt(1d0-gc_j**(-2))
            n3=comp*n4
            if (p_r <= 2d0) error stop 'secondary branch RHS requires p_r > 2'
            b_i=dsqrt(8d0*pi*b_r*e3_sec)
            gm_sec=1d0+e_r/fer*(p_r-2d0)/(p_r-1d0)*(e3_sec-e_ad)/(n3*Para_m_e*Para_c**2)
            gemax=3d0*Para_m_energy/dsqrt(8d0*b_i*Para_e**3)
            if (gm_sec >= gemax) error stop 'secondary branch RHS electron injection exceeds maximum'
            dm_dR=4d0*pi*RR*RR*n4*Para_m_p*gam2*(beta2-beta_s)/bc_j
            dM_shock=dm_dR*dR
            dV_shock=dM_shock/(n3*Para_m_p)
            dU_shock=(e3_sec-e_ad)*dV_shock
            D(m_idx)=dM_shock/mej
            D(u_idx)=(dU_shock+du_ad)/(mej*para_c**2)
            D(v_idx)=(dV_shock+dv_exp)/V3_scale
            D(e_idx)=dU_shock/(mej*para_c**2)
            D(g_idx)=gm_sec*dU_shock/(mej*para_c**2)
            D(h_idx)=gamma43_j*dM_shock/mej
            D(c_idx)=comp*dM_shock/mej
        end do
    end subroutine branch_derivs

end subroutine reverse_dynamics_rhs

end module reverse_rhs_module
