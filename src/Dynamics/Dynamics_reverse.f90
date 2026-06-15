subroutine dynamics_reverse(Delta_t,e_r,b_r,p_r,f_e_r,sigma_r,Boundary,n,Num_R, &
                            T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross, &
                            R_Tobs,R_Gamma,R,M2,M3,B3,U3_th,V3_comoving,Gamma34_inst, &
                            Secondary_M3,Secondary_U3,Secondary_V3,Secondary_B3, &
                            Secondary_M3_total,Secondary_U3_total,Secondary_V3_total,Secondary_B3_total, &
                            Secondary_pressure_total,Secondary_enthalpy_density_total, &
                            Secondary_gamma_contact,Secondary_pressure_3,Secondary_gamma_43,Secondary_beta_rs, &
                            Secondary_u_diss,Secondary_dissipated_energy,Secondary_electron_injected_energy, &
                            Secondary_branch_gamma_m,Secondary_branch_gamma_contact,Secondary_branch_gamma_43, &
                            Secondary_branch_compression,Secondary_branch_beta_rs,Secondary_branch_u_diss, &
                            Secondary_nu_m,Secondary_nu_c, &
                            Secondary_event_active,Secondary_start_radius,Secondary_end_radius, &
                            Secondary_start_tobs_axis,Secondary_end_tobs_axis)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_deceleration_radius, dynamics_rk4_reverse, &
                               dynamics_rk4_reverse_pre_m3, &
                               dynamics_reverse_rhs_iface, dynamics_log_time_step, &
                               dynamics_external_density_base, dynamics_external_density_profile, rs_mag_comp, &
                               dynamics_boundary_r0, dynamics_set_density_jump_profile, &
                               active_density_jump_count, density_jump_max, active_density_jump_r, &
                               active_density_jump_factor, active_density_jump_width
    implicit none
    integer, intent(in) :: n,Num_R
    integer :: I_tobs, Num_R1, Num_state
    procedure(dynamics_reverse_rhs_iface) :: reverse_dynamics_rhs
    real(8), intent(in) :: Boundary(n),Delta_t,e_r,b_r,p_r,f_e_r,sigma_r
    real(8), intent(out) :: T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross
    real(8), intent(out) :: R_Tobs(Num_R),R(Num_R),M2(Num_R),M3(Num_R),B3(Num_R),R_Gamma(Num_R)
    real(8), intent(out) :: U3_th(Num_R),V3_comoving(Num_R),Gamma34_inst(Num_R)
    real(8), intent(out) :: Secondary_M3(density_jump_max,Num_R),Secondary_U3(density_jump_max,Num_R)
    real(8), intent(out) :: Secondary_V3(density_jump_max,Num_R),Secondary_B3(density_jump_max,Num_R)
    real(8), intent(out) :: Secondary_M3_total(Num_R),Secondary_U3_total(Num_R)
    real(8), intent(out) :: Secondary_V3_total(Num_R),Secondary_B3_total(Num_R)
    real(8), intent(out) :: Secondary_pressure_total(Num_R),Secondary_enthalpy_density_total(Num_R)
    real(8), intent(out) :: Secondary_gamma_contact(Num_R),Secondary_pressure_3(Num_R),Secondary_gamma_43(Num_R)
    real(8), intent(out) :: Secondary_beta_rs(Num_R),Secondary_u_diss(Num_R),Secondary_dissipated_energy(Num_R)
    real(8), intent(out) :: Secondary_electron_injected_energy(Num_R),Secondary_branch_gamma_m(density_jump_max,Num_R)
    real(8), intent(out) :: Secondary_branch_gamma_contact(density_jump_max,Num_R)
    real(8), intent(out) :: Secondary_branch_gamma_43(density_jump_max,Num_R)
    real(8), intent(out) :: Secondary_branch_compression(density_jump_max,Num_R),Secondary_branch_beta_rs(density_jump_max,Num_R)
    real(8), intent(out) :: Secondary_branch_u_diss(density_jump_max,Num_R)
    real(8), intent(out) :: Secondary_nu_m(Num_R),Secondary_nu_c(Num_R)
    logical, intent(out) :: Secondary_event_active(density_jump_max)
    real(8), intent(out) :: Secondary_start_radius(density_jump_max),Secondary_end_radius(density_jump_max)
    real(8), intent(out) :: Secondary_start_tobs_axis(density_jump_max),Secondary_end_tobs_axis(density_jump_max)
    real(8) :: Eta_0,Epsilon_e,Epsilon_b,p_f,z,dNe_ISM,A_star,E_iso,T_log10_duration,f_e
    real(8) :: R_tr,f_jump,f_wide,R0
    real(8) :: Delta_0,para_m_ej,V3_scale,para_m2,para_m3,DM_0,R_dec,T00,t_dec,Grid_Tobs_bin,T_log10,T,H,dB3
    real(8) :: T_state,T_target
    real(8) :: u2_init,u4_init,Delta_init,para_n4_init,gam34_init,para_n3_init,comp_init
    real(8) :: event_prev_radius,event_prev_gamma,event_prev_tobs,event_curr_radius,event_curr_gamma,event_curr_tobs
    real(8) :: event_prev_source(density_jump_max),event_curr_source(density_jump_max)
    real(8) :: secondary_dissipated_prev(density_jump_max),secondary_gammam_prev(density_jump_max)
    real(8),allocatable :: Y(:),event_prev_state(:),event_curr_state(:)
    logical :: Secondary_event_closed(density_jump_max)

    Eta_0=Boundary(1); R(1)=Boundary(4); Epsilon_e=Boundary(5); Epsilon_b=Boundary(6); p_f=Boundary(7); z=Boundary(8)
    dNe_ISM=Boundary(11); A_star=Boundary(12); E_iso=Boundary(14); T_log10_duration=Boundary(15); f_e=Boundary(16)
    R_tr=Boundary(21); f_jump=Boundary(22); f_wide=Boundary(23)
    call dynamics_boundary_r0(Boundary,n,R0)
    call dynamics_set_density_jump_profile(Boundary,n)
    Num_state=6+5*active_density_jump_count
    allocate(Y(Num_state),event_prev_state(Num_state),event_curr_state(Num_state))
    Y=zero
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
    Y(1:6)=[R_Gamma(1),R(1),para_m2,para_m3/para_m_ej,(gam34_init-one)*para_m3/para_m_ej, &
            para_m3/(para_n3_init*Para_m_p)/V3_scale]

    DM_0=E_iso/((Eta_0-one)*4d0*pi*Para_m_p_E)
    call dynamics_deceleration_radius(A_star,dNe_ISM,Eta_0,DM_0,R_dec)
    T_cross=-1d0; U3_cross=zero; V3_cross=zero; M3_cross=zero; gam_m_cross=zero
    Secondary_M3=zero; Secondary_U3=zero; Secondary_V3=zero; Secondary_B3=zero
    Secondary_M3_total=zero; Secondary_U3_total=zero; Secondary_V3_total=zero; Secondary_B3_total=zero
    Secondary_pressure_total=zero; Secondary_enthalpy_density_total=zero
    Secondary_gamma_contact=zero; Secondary_pressure_3=zero; Secondary_gamma_43=one; Secondary_beta_rs=zero
    Secondary_u_diss=zero; Secondary_dissipated_energy=zero; Secondary_electron_injected_energy=zero
    Secondary_branch_gamma_m=zero; Secondary_branch_gamma_contact=zero; Secondary_branch_gamma_43=one
    Secondary_branch_compression=one
    Secondary_branch_beta_rs=zero; Secondary_branch_u_diss=zero; Secondary_nu_m=zero; Secondary_nu_c=zero
    secondary_dissipated_prev=zero; secondary_gammam_prev=zero
    Secondary_event_active=.false.; Secondary_event_closed=.false.
    Secondary_start_radius=zero; Secondary_end_radius=zero
    Secondary_start_tobs_axis=zero; Secondary_end_tobs_axis=zero
    B3_ordered_cross=zero
    ! RS 状态向量中 Y(2) 是 shock 半径；初始到达时刻必须由半径而不是 swept mass 定义。
    T00=Y(2)*(one/dsqrt(one-one/Eta_0**2)-one)/Para_c
    t_dec=R_dec/(two*Eta_0*Eta_0*Para_c)
    Grid_Tobs_bin=min(log10(T00)-2d0,dlog10(t_dec*0.1d0))
    T_log10=T_log10_duration-Grid_Tobs_bin; Num_R1=Num_R-1
    T_state=T00
    event_prev_radius=Y(2); event_prev_gamma=Y(1); event_prev_tobs=T00*(one+z)
    event_prev_state=Y
    call secondary_reverse_event_sources(event_prev_radius,event_prev_gamma,event_prev_state,event_prev_source)

    do I_tobs=1,Num_R
        call dynamics_log_time_step(T00,Grid_Tobs_bin,T_log10,Num_R1,I_tobs,T_target,H)
        do while (T_state < T_target)
            if (T_cross < zero .and. Y(4) < one) then
                call dynamics_rk4_reverse_pre_m3(reverse_dynamics_rhs,dB3,T_cross,R_cross,e3_cross,gam20, &
                                                 U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross, &
                                                 T_state,T_target,Y,Num_state,para_m_ej,V3_scale,Delta_0,eta_0,A_star,dNe_ISM, &
                                                 R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r,sigma_r)
            else
                H=T_target-T_state
                T=T_state
                call dynamics_rk4_reverse(reverse_dynamics_rhs,dB3,T_cross,R_cross,e3_cross,gam20, &
                                          U3_cross,V3_cross,M3_cross, &
                                          gam_m_cross,B3_ordered_cross,T,H,Y,Num_state,para_m_ej,V3_scale,Delta_0, &
                                          eta_0,A_star,dNe_ISM,R_tr,f_jump,f_wide,R0, &
                                          Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r,sigma_r)
                T_state=T_target
            end if
        end do
        R_Tobs(I_tobs)=T_target*(one+z); R_Gamma(I_tobs)=Y(1); R(I_tobs)=Y(2); M2(I_tobs)=Y(3)
        event_curr_radius=Y(2); event_curr_gamma=Y(1); event_curr_tobs=R_Tobs(I_tobs)
        event_curr_state=Y
        call secondary_reverse_event_sources(event_curr_radius,event_curr_gamma,event_curr_state,event_curr_source)
        call update_secondary_reverse_events()
        event_prev_radius=event_curr_radius; event_prev_gamma=event_curr_gamma
        event_prev_tobs=event_curr_tobs; event_prev_source=event_curr_source
        event_prev_state=event_curr_state
        M3(I_tobs)=Y(4)*para_m_ej; U3_th(I_tobs)=Y(5)*para_m_ej*para_c**2
        V3_comoving(I_tobs)=Y(6)*V3_scale; B3(I_tobs)=dB3
        call store_secondary_branch_state(I_tobs)
        Gamma34_inst(I_tobs)=(Y(1)*Y(1)+Eta_0*Eta_0-one)/(Eta_0*Y(1)+dsqrt(Y(1)*Y(1)-one)*u4_init)
    end do

    call close_open_secondary_reverse_events()
    deallocate(Y,event_prev_state,event_curr_state)

contains

    ! Secondary RS 事件诊断：用动力学积分接受态的 mechanical source 符号定位 start/end。
    subroutine update_secondary_reverse_events()
    implicit none
    integer :: j,k,n_scan
    real(8) :: radius_lo,radius_hi,gamma_lo,gamma_hi,tobs_lo,tobs_hi,source_lo,source_hi
    real(8) :: width,overlap_lo,overlap_hi,scan_length
    real(8) :: state_lo(Num_state),state_hi(Num_state),frac_scan

        do j=1,active_density_jump_count
            width=active_density_jump_width(j)*active_density_jump_r(j)
            overlap_lo=max(event_prev_radius,active_density_jump_r(j)-4d0*width)
            overlap_hi=min(event_curr_radius,active_density_jump_r(j))
            scan_length=max(zero,overlap_hi-overlap_lo)
            n_scan=max(1,ceiling(scan_length/(width/16d0)))
            radius_lo=event_prev_radius; gamma_lo=event_prev_gamma
            tobs_lo=event_prev_tobs; source_lo=event_prev_source(j)
            state_lo=event_prev_state
            do k=1,n_scan
                frac_scan=dble(k)/dble(n_scan)
                radius_hi=event_prev_radius+(event_curr_radius-event_prev_radius)*frac_scan
                gamma_hi=event_prev_gamma+(event_curr_gamma-event_prev_gamma)*frac_scan
                tobs_hi=event_prev_tobs+(event_curr_tobs-event_prev_tobs)*frac_scan
                state_hi=event_prev_state+frac_scan*(event_curr_state-event_prev_state)
                call secondary_reverse_event_source(j,radius_hi,gamma_hi,state_hi,source_hi)
                call record_secondary_reverse_event_segment(j,radius_lo,radius_hi,gamma_lo,gamma_hi,tobs_lo,tobs_hi, &
                                                            state_lo,state_hi,source_lo,source_hi)
                radius_lo=radius_hi; gamma_lo=gamma_hi; tobs_lo=tobs_hi; source_lo=source_hi; state_lo=state_hi
            end do
        end do
    end subroutine update_secondary_reverse_events

    subroutine record_secondary_reverse_event_segment(j,radius_lo,radius_hi,gamma_lo,gamma_hi,tobs_lo,tobs_hi, &
                                                      state_lo,state_hi,source_lo,source_hi)
    implicit none
    integer, intent(in) :: j
    real(8), intent(in) :: radius_lo,radius_hi,gamma_lo,gamma_hi,tobs_lo,tobs_hi,source_lo,source_hi
    real(8), intent(in) :: state_lo(Num_state),state_hi(Num_state)
    real(8) :: root_radius,root_tobs

        if (.not. Secondary_event_active(j)) then
            if (source_lo > zero) then
                Secondary_event_active(j)=.true.
                Secondary_start_radius(j)=radius_lo
                Secondary_start_tobs_axis(j)=tobs_lo
            else if (source_lo <= zero .and. source_hi > zero) then
                call secondary_reverse_event_root_between(j,radius_lo,radius_hi,gamma_lo,gamma_hi, &
                                                          tobs_lo,tobs_hi,state_lo,state_hi,root_radius,root_tobs)
                Secondary_event_active(j)=.true.
                Secondary_start_radius(j)=root_radius
                Secondary_start_tobs_axis(j)=root_tobs
            end if
        end if
        if (Secondary_event_active(j) .and. .not. Secondary_event_closed(j)) then
            if (source_lo > zero .and. source_hi <= zero) then
                call secondary_reverse_event_root_between(j,radius_lo,radius_hi,gamma_lo,gamma_hi, &
                                                          tobs_lo,tobs_hi,state_lo,state_hi,root_radius,root_tobs)
                Secondary_event_closed(j)=.true.
                if (root_radius >= active_density_jump_r(j)) root_radius=nearest(active_density_jump_r(j),-one)
                Secondary_end_radius(j)=root_radius
                Secondary_end_tobs_axis(j)=root_tobs
            end if
        end if
    end subroutine record_secondary_reverse_event_segment

    subroutine secondary_reverse_event_root_between(jump_index,r_lo_in,r_hi_in,g_lo_in,g_hi_in,t_lo_in,t_hi_in, &
                                                    state_lo_in,state_hi_in,root_r,root_t)
    implicit none
    integer, intent(in) :: jump_index
    integer :: iter
    real(8), intent(in) :: r_lo_in,r_hi_in,g_lo_in,g_hi_in,t_lo_in,t_hi_in
    real(8), intent(in) :: state_lo_in(Num_state),state_hi_in(Num_state)
    real(8), intent(out) :: root_r,root_t
    real(8) :: r_lo,r_hi,r_mid,gamma_lo,gamma_hi,gamma_mid,t_lo,t_hi,t_mid,source_lo,source_mid,frac
    real(8) :: state_mid(Num_state)

        r_lo=r_lo_in; r_hi=r_hi_in
        gamma_lo=g_lo_in; gamma_hi=g_hi_in
        t_lo=t_lo_in; t_hi=t_hi_in
        call secondary_reverse_event_source(jump_index,r_lo,gamma_lo,state_lo_in,source_lo)
        do iter=1,80
            r_mid=0.5d0*(r_lo+r_hi)
            frac=(r_mid-r_lo_in)/(r_hi_in-r_lo_in)
            gamma_mid=g_lo_in+frac*(g_hi_in-g_lo_in)
            t_mid=t_lo_in+frac*(t_hi_in-t_lo_in)
            state_mid=state_lo_in+frac*(state_hi_in-state_lo_in)
            call secondary_reverse_event_source(jump_index,r_mid,gamma_mid,state_mid,source_mid)
            if (source_lo*source_mid <= zero) then
                r_hi=r_mid; gamma_hi=gamma_mid; t_hi=t_mid
            else
                r_lo=r_mid; gamma_lo=gamma_mid; t_lo=t_mid; source_lo=source_mid
            end if
        end do
        root_r=0.5d0*(r_lo+r_hi)
        root_t=0.5d0*(t_lo+t_hi)
    end subroutine secondary_reverse_event_root_between

    subroutine close_open_secondary_reverse_events()
    implicit none
    integer :: j

        do j=1,active_density_jump_count
            if (Secondary_event_active(j) .and. .not. Secondary_event_closed(j)) then
                Secondary_event_closed(j)=.true.
                Secondary_end_radius(j)=event_prev_radius
                Secondary_end_tobs_axis(j)=event_prev_tobs
            end if
        end do
    end subroutine close_open_secondary_reverse_events

    subroutine secondary_reverse_event_sources(radius,gamma_bulk,state,sources)
    implicit none
    integer :: j
    real(8), intent(in) :: radius,gamma_bulk
    real(8), intent(in) :: state(Num_state)
    real(8), intent(out) :: sources(density_jump_max)

        sources=-one
        do j=1,active_density_jump_count
            call secondary_reverse_event_source(j,radius,gamma_bulk,state,sources(j))
        end do
    end subroutine secondary_reverse_event_sources

    subroutine secondary_reverse_event_source(jump_index,radius,gamma_bulk,state,source)
    implicit none
    integer, intent(in) :: jump_index
    integer :: parent_m_idx,parent_u_idx,parent_v_idx
    real(8), intent(in) :: radius,gamma_bulk
    real(8), intent(in) :: state(Num_state)
    real(8), intent(out) :: source
    real(8) :: density_factor,branch_weight,n1,n_excess,n_pre,n4,e4,p4,p3,e3_sec,e_ad
    real(8) :: gamma_c_j,gamma43_j,comp,beta_s

        call secondary_reverse_density_branch_state(radius,jump_index,density_factor,branch_weight)
        if (branch_weight <= zero .or. gamma_bulk <= one) then
            source=-one
            return
        end if
        n1=density_factor
        n_excess=branch_weight
        n_pre=n1-n_excess
        if (n_pre <= zero) error stop 'secondary reverse event source found non-positive pre-bump density'
        n4=4d0*gamma_bulk*n_pre
        e4=4d0*gamma_bulk*(gamma_bulk-one)*n_pre*Para_m_p*Para_c**2
        p4=e4/3d0
        if (jump_index > 1) then
            parent_m_idx=6+jump_index-1
            parent_u_idx=6+active_density_jump_count+jump_index-1
            parent_v_idx=6+2*active_density_jump_count+jump_index-1
            if (secondary_parent_upstream_available(jump_index,state(parent_m_idx),state(parent_u_idx),state(parent_v_idx))) then
                n4=state(parent_m_idx)*para_m_ej/(Para_m_p*state(parent_v_idx)*V3_scale)
                e4=state(parent_u_idx)*para_m_ej*para_c**2/(state(parent_v_idx)*V3_scale)
                p4=e4/3d0
            end if
        end if
        call secondary_reverse_contact_rh(gamma_bulk,n1,n4,e4,p4,gamma_c_j,p3,gamma43_j,comp,beta_s)
        if (comp <= zero) then
            source=-one
            return
        end if
        e3_sec=3d0*p3
        e_ad=e4*comp**(4d0/3d0)
        source=e3_sec-e_ad
    end subroutine secondary_reverse_event_source

    subroutine store_secondary_branch_state(i_out)
    implicit none
    integer, intent(in) :: i_out
    integer :: j, m_idx, u_idx, v_idx, e_idx, g_idx
    real(8) :: e_cum,g_cum,e_delta,g_delta,diag_total,gamma_m_total
    real(8) :: density_factor,branch_weight,n1,n_excess,n_pre,n4,e4,p4,p3,e3_sec,e_ad,comp,beta_s
    real(8) :: gamma_c_j,gamma43_j,n3,gamma_m_j,b_i,gam_e_max

        diag_total=zero; gamma_m_total=zero
        do j=1,active_density_jump_count
            m_idx=6+j
            u_idx=6+active_density_jump_count+j
            v_idx=6+2*active_density_jump_count+j
            e_idx=6+3*active_density_jump_count+j
            g_idx=6+4*active_density_jump_count+j
            Secondary_M3(j,i_out)=Y(m_idx)*para_m_ej
            Secondary_U3(j,i_out)=Y(u_idx)*para_m_ej*para_c**2
            Secondary_V3(j,i_out)=Y(v_idx)*V3_scale
            if (Secondary_V3(j,i_out) > zero) &
                Secondary_B3(j,i_out)=dsqrt(8d0*pi*b_r*Secondary_U3(j,i_out)/Secondary_V3(j,i_out))
            e_cum=Y(e_idx)*para_m_ej*para_c**2
            g_cum=Y(g_idx)*para_m_ej*para_c**2
            e_delta=e_cum-secondary_dissipated_prev(j)
            g_delta=g_cum-secondary_gammam_prev(j)
            if (e_delta < zero) error stop 'secondary RS dissipated energy must not decrease'
            Secondary_dissipated_energy(i_out)=Secondary_dissipated_energy(i_out)+e_delta
            Secondary_electron_injected_energy(i_out)=Secondary_electron_injected_energy(i_out)+e_r*e_delta
            if (e_delta > zero) then
                Secondary_branch_gamma_m(j,i_out)=g_delta/e_delta
                call secondary_reverse_density_branch_state(R(I_tobs),j,density_factor,branch_weight)
                if (branch_weight > zero .and. R_Gamma(i_out) > one) then
                    n1=density_factor
                    n_excess=branch_weight
                    n_pre=n1-n_excess
                    if (n_pre <= zero) error stop 'secondary branch output found non-positive pre-bump density'
                    n4=4d0*R_Gamma(i_out)*n_pre
                    e4=4d0*R_Gamma(i_out)*(R_Gamma(i_out)-one)*n_pre*Para_m_p*Para_c**2
                    p4=e4/3d0
                    if (j > 1) then
                        if (secondary_parent_upstream_available(j,Secondary_M3(j-1,i_out)/para_m_ej, &
                                                                Secondary_U3(j-1,i_out)/(para_m_ej*para_c**2), &
                                                                Secondary_V3(j-1,i_out)/V3_scale)) then
                            n4=Secondary_M3(j-1,i_out)/(Para_m_p*Secondary_V3(j-1,i_out))
                            e4=Secondary_U3(j-1,i_out)/Secondary_V3(j-1,i_out)
                            p4=e4/3d0
                        end if
                    end if
                    call secondary_reverse_contact_rh(R_Gamma(i_out),n1,n4,e4,p4,gamma_c_j,p3,gamma43_j,comp,beta_s)
                    if (comp > zero) then
                        e3_sec=3d0*p3
                        e_ad=e4*comp**(4d0/3d0)
                        n3=comp*n4
                        if (p_r <= two) error stop 'secondary RS injection requires p_r > 2'
                        b_i=dsqrt(8d0*pi*b_r*e3_sec)
                        gamma_m_j=one+e_r/f_e_r*(p_r-two)/(p_r-one)*(e3_sec-e_ad)/(n3*Para_m_e*Para_c**2)
                        gam_e_max=3d0*Para_m_energy/dsqrt(8d0*b_i*Para_e**3)
                        if (gamma_m_j >= gam_e_max) error stop 'secondary RS electron injection exceeds maximum'
                        diag_total=diag_total+e_delta; gamma_m_total=gamma_m_total+g_delta
                        Secondary_gamma_contact(i_out)=Secondary_gamma_contact(i_out)+e_delta*gamma_c_j
                        Secondary_pressure_3(i_out)=Secondary_pressure_3(i_out)+e_delta*p3
                        Secondary_gamma_43(i_out)=Secondary_gamma_43(i_out)+e_delta*(gamma43_j-one)
                        Secondary_beta_rs(i_out)=Secondary_beta_rs(i_out)+e_delta*beta_s
                        Secondary_u_diss(i_out)=Secondary_u_diss(i_out)+(e3_sec-e_ad)
                        Secondary_branch_gamma_contact(j,i_out)=gamma_c_j
                        Secondary_branch_gamma_43(j,i_out)=gamma43_j
                        Secondary_branch_compression(j,i_out)=comp
                        Secondary_branch_beta_rs(j,i_out)=beta_s
                        Secondary_branch_u_diss(j,i_out)=e3_sec-e_ad
                    end if
                end if
            end if
            secondary_dissipated_prev(j)=e_cum
            secondary_gammam_prev(j)=g_cum
            Secondary_M3_total(i_out)=Secondary_M3_total(i_out)+Secondary_M3(j,i_out)
            Secondary_U3_total(i_out)=Secondary_U3_total(i_out)+Secondary_U3(j,i_out)
            Secondary_V3_total(i_out)=Secondary_V3_total(i_out)+Secondary_V3(j,i_out)
        end do
        if (diag_total > zero) then
            Secondary_gamma_contact(i_out)=Secondary_gamma_contact(i_out)/diag_total
            Secondary_pressure_3(i_out)=Secondary_pressure_3(i_out)/diag_total
            Secondary_gamma_43(i_out)=one+Secondary_gamma_43(i_out)/diag_total
            Secondary_beta_rs(i_out)=Secondary_beta_rs(i_out)/diag_total
        end if
        if (Secondary_V3_total(i_out) > zero) &
            Secondary_B3_total(i_out)=dsqrt(8d0*pi*b_r*Secondary_U3_total(i_out)/Secondary_V3_total(i_out))
        if (Secondary_V3_total(i_out) > zero) then
            Secondary_pressure_total(i_out)=Secondary_U3_total(i_out)/(3d0*Secondary_V3_total(i_out))
            Secondary_enthalpy_density_total(i_out)=Secondary_M3_total(i_out)*para_c**2/Secondary_V3_total(i_out)+ &
                Secondary_U3_total(i_out)/Secondary_V3_total(i_out)+Secondary_pressure_total(i_out)
            if (diag_total > zero .and. Secondary_gamma_contact(i_out) > one) then
                Secondary_nu_m(i_out)=4.2d6*Secondary_B3_total(i_out)* &
                    (Secondary_gamma_contact(i_out)*(one-dsqrt(one-Secondary_gamma_contact(i_out)**(-2)))*(one+z))**(-1) * &
                    (gamma_m_total/diag_total)**2
            end if
            if (Secondary_B3_total(i_out) > zero .and. R_Gamma(i_out) > one .and. R_Tobs(i_out) > zero) then
                associate(gamma_cool => 7.7d8*(one+z)/(R_Gamma(i_out)*Secondary_B3_total(i_out)**2*R_Tobs(i_out)))
                    Secondary_nu_c(i_out)=4.2d6*Secondary_B3_total(i_out)*gamma_cool*gamma_cool/ &
                        (R_Gamma(i_out)*(one-dsqrt(one-R_Gamma(i_out)**(-2)))*(one+z))
                end associate
            end if
        end if
    end subroutine store_secondary_branch_state

    logical function secondary_parent_upstream_available(jump_index,parent_mass,parent_energy,parent_volume)
    implicit none
    integer, intent(in) :: jump_index
    real(8), intent(in) :: parent_mass,parent_energy,parent_volume
        secondary_parent_upstream_available=.false.
        if (jump_index <= 1) return
        if (parent_mass <= zero .or. parent_energy <= zero .or. parent_volume <= zero) return
        secondary_parent_upstream_available=.true.
    end function secondary_parent_upstream_available

    subroutine secondary_reverse_density_branch_state(radius,jump_index,density_factor,branch_weight)
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
    end subroutine secondary_reverse_density_branch_state
end subroutine dynamics_reverse

subroutine secondary_reverse_contact_rh(gamma4,n1,n4,e4,p4,gamma_c,p3,gamma43,n3_over_n4,beta_rs)
    use constants
    implicit none
    integer :: I
    real(8), intent(in) :: gamma4,n1,n4,e4,p4
    real(8), intent(out) :: gamma_c,p3,gamma43,n3_over_n4,beta_rs
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
    call solve_reverse_shock_compression(n3_over_n4,beta_rs)

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

    subroutine solve_reverse_shock_compression(comp,beta_rs_out)
    implicit none
    real(8), intent(out) :: comp,beta_rs_out
    integer :: K
    real(8) :: beta_s_lo,beta_s_hi,beta_s_mid,g_lo,g_hi,g_mid,eps_beta,beta_comp_one

        if (beta4 <= beta_c) error stop 'secondary_reverse_contact_rh reverse shock speed bracket is empty'
        beta_rs_out=zero
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
        beta_rs_out=0.5d0*(beta_s_lo+beta_s_hi)
        call shock_momentum_difference(beta_rs_out,g_mid,comp)
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

subroutine secondary_reverse_profile(Num_R,Num_jump,R,Tobs_axis,Gamma4,dNe_ISM,Jump_R,Jump_factor,Jump_width, &
                                     Epsilon_e,Epsilon_b,p_e,f_e,z, &
                                     gamma_contact,pressure_3,gamma_43,comp_ratio,beta_rs,u_diss,active_weight, &
                                     m3_reservoir,u3_reservoir,v3_reservoir,b3_reservoir,gamma_m_shell, &
                                     dissipated_energy,electron_injected_energy,pressure_total,enthalpy_density_total, &
                                     m3_branch,u3_branch,v3_branch,b3_branch,gamma_m_branch, &
                                     nu_m,nu_c,event_active, &
                                     start_radius,shock_end_radius,start_tobs_axis,shock_end_tobs_axis)
    use constants
    implicit none
    integer, intent(in) :: Num_R,Num_jump
    integer :: I,J,K
    real(8), intent(in) :: R(Num_R),Tobs_axis(Num_R),Gamma4(Num_R),dNe_ISM
    real(8), intent(in) :: Jump_R(Num_jump),Jump_factor(Num_jump),Jump_width(Num_jump)
    real(8), intent(in) :: Epsilon_e,Epsilon_b,p_e,f_e,z
    real(8), intent(out) :: gamma_contact(Num_R),pressure_3(Num_R),gamma_43(Num_R),comp_ratio(Num_R)
    real(8), intent(out) :: beta_rs(Num_R),u_diss(Num_R),active_weight(Num_R)
    real(8), intent(out) :: m3_reservoir(Num_R),u3_reservoir(Num_R),v3_reservoir(Num_R),b3_reservoir(Num_R)
    real(8), intent(out) :: gamma_m_shell(Num_R),dissipated_energy(Num_R),electron_injected_energy(Num_R)
    real(8), intent(out) :: pressure_total(Num_R),enthalpy_density_total(Num_R)
    real(8), intent(out) :: m3_branch(Num_jump,Num_R),u3_branch(Num_jump,Num_R),v3_branch(Num_jump,Num_R)
    real(8), intent(out) :: b3_branch(Num_jump,Num_R),gamma_m_branch(Num_jump,Num_R)
    real(8), intent(out) :: nu_m(Num_R),nu_c(Num_R)
    logical, intent(out) :: event_active(Num_jump)
    real(8), intent(out) :: start_radius(Num_jump),shock_end_radius(Num_jump)
    real(8), intent(out) :: start_tobs_axis(Num_jump),shock_end_tobs_axis(Num_jump)
    real(8) :: density_factor,branch_weight,n1,n_excess,n_pre,n4,e4,p4,p3,e3,e_ad,comp,beta_s
    real(8) :: gamma_c_j,gamma43_j
    real(8) :: beta4,beta_c,n3,d_radius,shell_mass,shell_volume,u_inj,gamma_m,gam_e_max,b_i,volume_k,energy_k
    real(8) :: diag_weight,diag_total
    real(8) :: doppler_den,gamma_cool
    real(8) :: width_cm,lower_bound,upper_bound,start_root,end_root

    gamma_contact=zero; pressure_3=zero; gamma_43=one; comp_ratio=zero; beta_rs=zero
    u_diss=zero; active_weight=zero; m3_reservoir=zero; u3_reservoir=zero; v3_reservoir=zero
    b3_reservoir=zero; gamma_m_shell=zero; dissipated_energy=zero; electron_injected_energy=zero
    pressure_total=zero; enthalpy_density_total=zero
    m3_branch=zero; u3_branch=zero; v3_branch=zero; b3_branch=zero; gamma_m_branch=zero
    nu_m=zero; nu_c=zero
    event_active=.false.; start_radius=zero; shock_end_radius=zero
    start_tobs_axis=zero; shock_end_tobs_axis=zero

    do I=2,Num_R
        diag_total=zero
        do J=1,Num_jump
            call secondary_reverse_density_branch(R(I),J,density_factor,branch_weight)
            active_weight(I)=active_weight(I)+branch_weight
            if (branch_weight <= zero .or. Gamma4(I) <= one) cycle
            n1=density_factor
            n_excess=branch_weight
            n_pre=n1-n_excess
            if (n_pre <= zero) error stop 'secondary_reverse_profile found non-positive pre-bump density'
            n4=4d0*Gamma4(I)*n_pre
            e4=4d0*Gamma4(I)*(Gamma4(I)-one)*n_pre*Para_m_p*Para_c**2
            p4=e4/3d0
            call secondary_reverse_contact_rh(Gamma4(I),n1,n4,e4,p4,gamma_c_j,p3,gamma43_j,comp,beta_s)
            if (comp <= zero) cycle
            e3=3d0*p3
            e_ad=e4*comp**(4d0/3d0)
            if (e3 <= e_ad) cycle
            beta4=dsqrt(one-Gamma4(I)**(-2))
            beta_c=dsqrt(one-gamma_c_j**(-2))
            n3=comp*n4
            d_radius=R(I)-R(I-1)
            shell_mass=4d0*pi*R(I)*R(I)*d_radius*n4*Para_m_p*Gamma4(I)*(beta4-beta_s)/beta_c
            shell_volume=shell_mass/(n3*Para_m_p)
            if (p_e <= two) error stop 'secondary_reverse_profile requires p_e > 2 for secondary RS injection'
            b_i=dsqrt(8d0*pi*Epsilon_b*3d0*p3)
            gamma_m=one+Epsilon_e/f_e*(p_e-two)/(p_e-one)*(e3-e_ad)/(n3*Para_m_e*Para_c**2)
            gam_e_max=3d0*Para_m_energy/dsqrt(8d0*b_i*Para_e**3)
            if (gamma_m >= gam_e_max) error stop 'secondary_reverse_profile electron injection exceeds synchrotron maximum'
            u_inj=(e3-e_ad)*shell_volume
            gamma_m_branch(J,I)=gamma_m
            diag_weight=u_inj
            diag_total=diag_total+diag_weight
            gamma_contact(I)=gamma_contact(I)+diag_weight*gamma_c_j
            pressure_3(I)=pressure_3(I)+diag_weight*p3
            comp_ratio(I)=comp_ratio(I)+diag_weight*comp
            beta_rs(I)=beta_rs(I)+diag_weight*beta_s
            gamma_43(I)=gamma_43(I)+diag_weight*(gamma43_j-one)
            u_diss(I)=u_diss(I)+(e3-e_ad)
            gamma_m_shell(I)=gamma_m_shell(I)+diag_weight*gamma_m
            dissipated_energy(I)=dissipated_energy(I)+u_inj
            electron_injected_energy(I)=electron_injected_energy(I)+Epsilon_e*u_inj
            do K=I,Num_R
                volume_k=shell_volume*(R(K)/R(I))**3*(gamma_c_j/Gamma4(K))
                energy_k=u_inj*(shell_volume/volume_k)**(one/3d0)
                m3_branch(J,K)=m3_branch(J,K)+shell_mass
                u3_branch(J,K)=u3_branch(J,K)+energy_k
                v3_branch(J,K)=v3_branch(J,K)+volume_k
                m3_reservoir(K)=m3_reservoir(K)+shell_mass
                u3_reservoir(K)=u3_reservoir(K)+energy_k
                v3_reservoir(K)=v3_reservoir(K)+volume_k
            end do
        end do
        if (diag_total > zero) then
            gamma_contact(I)=gamma_contact(I)/diag_total
            pressure_3(I)=pressure_3(I)/diag_total
            comp_ratio(I)=comp_ratio(I)/diag_total
            beta_rs(I)=beta_rs(I)/diag_total
            gamma_43(I)=one+gamma_43(I)/diag_total
            gamma_m_shell(I)=gamma_m_shell(I)/diag_total
        else
            gamma_43(I)=one
        end if
    end do

    do I=1,Num_R
        if (v3_reservoir(I) > zero) then
            pressure_total(I)=u3_reservoir(I)/(3d0*v3_reservoir(I))
            enthalpy_density_total(I)=m3_reservoir(I)*Para_c**2/v3_reservoir(I)+u3_reservoir(I)/v3_reservoir(I)+pressure_total(I)
            b3_reservoir(I)=dsqrt(8d0*pi*Epsilon_b*u3_reservoir(I)/v3_reservoir(I))
            if (gamma_m_shell(I) > one) then
                doppler_den=gamma_contact(I)*(one-dsqrt(one-gamma_contact(I)**(-2)))*(one+z)
                nu_m(I)=4.2d6*b3_reservoir(I)*gamma_m_shell(I)*gamma_m_shell(I)/doppler_den
            end if
            doppler_den=Gamma4(I)*(one-dsqrt(one-Gamma4(I)**(-2)))*(one+z)
            gamma_cool=7.7d8*(one+z)/Gamma4(I)/(b3_reservoir(I)*b3_reservoir(I))/Tobs_axis(I)
            nu_c(I)=4.2d6*b3_reservoir(I)*gamma_cool*gamma_cool/doppler_den
        end if
        do J=1,Num_jump
            if (v3_branch(J,I) > zero) b3_branch(J,I)=dsqrt(8d0*pi*Epsilon_b*u3_branch(J,I)/v3_branch(J,I))
        end do
    end do

    do J=1,Num_jump
        width_cm=Jump_width(J)*Jump_R(J)
        lower_bound=Jump_R(J)-4d0*width_cm
        upper_bound=nearest(Jump_R(J),lower_bound-Jump_R(J))
        call secondary_reverse_scan_event_window(J,lower_bound,upper_bound,start_root,end_root,event_active(J))
        if (event_active(J)) then
            start_radius(J)=start_root
            shock_end_radius(J)=end_root
            call secondary_reverse_interp_tobs(Num_R,R,Tobs_axis,start_root,start_tobs_axis(J))
            call secondary_reverse_interp_tobs(Num_R,R,Tobs_axis,end_root,shock_end_tobs_axis(J))
        end if
    end do

contains

    subroutine secondary_reverse_scan_event_window(jump_index,lower_bound,upper_bound,start_root,end_root,active)
    implicit none
    integer, intent(in) :: jump_index
    integer :: K
    logical, intent(out) :: active
    real(8), intent(in) :: lower_bound,upper_bound
    real(8), intent(out) :: start_root,end_root
    real(8), parameter :: scan_intervals=64d0
    real(8) :: r_lo,r_hi,g_lo,g_hi

        active=.false.; start_root=zero; end_root=zero
        r_lo=lower_bound
        call secondary_reverse_signed_source(r_lo,jump_index,g_lo)
        if (g_lo > zero) then
            active=.true.
            start_root=r_lo
        end if
        do K=1,64
            r_hi=lower_bound+(upper_bound-lower_bound)*dble(K)/scan_intervals
            call secondary_reverse_signed_source(r_hi,jump_index,g_hi)
            if (.not. active .and. g_lo <= zero .and. g_hi > zero) then
                active=.true.
                start_root=secondary_reverse_event_root_bracket(r_lo,r_hi,jump_index)
            end if
            if (active .and. g_lo > zero .and. g_hi <= zero) then
                end_root=secondary_reverse_event_root_bracket(r_lo,r_hi,jump_index)
                return
            end if
            r_lo=r_hi; g_lo=g_hi
        end do
        if (active) end_root=upper_bound
    end subroutine secondary_reverse_scan_event_window

    ! 对单个 density bump 取上升段权重；总密度仍包含其它 bump 的背景压力贡献。
    subroutine secondary_reverse_density_branch(radius,jump_index,density_factor,branch_weight)
    implicit none
    integer, intent(in) :: jump_index
    integer :: K
    real(8), intent(in) :: radius
    real(8), intent(out) :: density_factor,branch_weight
    real(8) :: x,width,profile

        density_factor=one; branch_weight=zero
        do K=1,Num_jump
            x=radius-Jump_R(K)
            width=Jump_width(K)*Jump_R(K)
            profile=(Jump_factor(K)-one)*dexp(-(x*x)/(2d0*width*width))
            density_factor=density_factor+profile
            if (K == jump_index .and. x >= -4d0*width .and. x < zero) branch_weight=profile
        end do
    end subroutine secondary_reverse_density_branch

    real(8) function secondary_reverse_event_root_bracket(lo_r_in,hi_r_in,jump_index)
    implicit none
    integer :: K
    integer, intent(in) :: jump_index
    real(8), intent(in) :: lo_r_in,hi_r_in
    real(8) :: lo_r,hi_r,mid_r,g_lo,g_hi,g_mid

        lo_r=lo_r_in; hi_r=hi_r_in
        call secondary_reverse_signed_source(lo_r,jump_index,g_lo)
        call secondary_reverse_signed_source(hi_r,jump_index,g_hi)
        if (g_lo == zero) then
            secondary_reverse_event_root_bracket=lo_r
            return
        end if
        if (g_hi == zero) then
            secondary_reverse_event_root_bracket=hi_r
            return
        end if
        if (g_lo*g_hi > zero) then
            error stop 'secondary_reverse_event_root_bracket requires a signed source bracket'
        end if
        do K=1,80
            mid_r=0.5d0*(lo_r+hi_r)
            call secondary_reverse_signed_source(mid_r,jump_index,g_mid)
            if (g_lo*g_mid <= zero) then
                hi_r=mid_r; g_hi=g_mid
            else
                lo_r=mid_r; g_lo=g_mid
            end if
        end do
        secondary_reverse_event_root_bracket=0.5d0*(lo_r+hi_r)
    end function secondary_reverse_event_root_bracket

    subroutine secondary_reverse_signed_source(radius,jump_index,source)
    implicit none
    integer, intent(in) :: jump_index
    real(8), intent(in) :: radius
    real(8), intent(out) :: source
    real(8) :: gamma_loc,density_factor_loc,local_weight_loc
    real(8) :: n1_loc,n_excess_loc,n_pre_loc,n4_loc,e4_loc,p4_loc,p3_loc
    real(8) :: gamma_contact_loc,gamma43_loc,comp_loc,beta_s_loc,e3_loc,e_ad_loc

        call secondary_reverse_interp_gamma(radius,gamma_loc)
        call secondary_reverse_density_branch(radius,jump_index,density_factor_loc,local_weight_loc)
        if (local_weight_loc <= zero .or. gamma_loc <= one) then
            source=-one
            return
        end if
        n1_loc=dNe_ISM*density_factor_loc
        n_excess_loc=dNe_ISM*local_weight_loc
        n_pre_loc=n1_loc-n_excess_loc
        if (n_pre_loc <= zero) error stop 'secondary_reverse_event root found non-positive pre-bump density'
        n4_loc=4d0*gamma_loc*n_pre_loc
        e4_loc=4d0*gamma_loc*(gamma_loc-one)*n_pre_loc*Para_m_p*Para_c**2
        p4_loc=e4_loc/3d0
        call secondary_reverse_contact_rh(gamma_loc,n1_loc,n4_loc,e4_loc,p4_loc, &
                                          gamma_contact_loc,p3_loc,gamma43_loc,comp_loc,beta_s_loc)
        if (comp_loc <= zero) then
            source=-one
            return
        end if
        e3_loc=3d0*p3_loc
        e_ad_loc=e4_loc*comp_loc**(4d0/3d0)
        source=e3_loc-e_ad_loc
    end subroutine secondary_reverse_signed_source

    subroutine secondary_reverse_interp_gamma(radius,gamma_root)
    implicit none
    integer :: K
    real(8), intent(in) :: radius
    real(8), intent(out) :: gamma_root
    real(8) :: ratio

        if (radius <= R(1)) then
            gamma_root=Gamma4(1)
            return
        end if
        do K=1,Num_R-1
            if (radius <= R(K+1)) then
                ratio=(radius-R(K))/(R(K+1)-R(K))
                gamma_root=Gamma4(K)+ratio*(Gamma4(K+1)-Gamma4(K))
                return
            end if
        end do
        gamma_root=Gamma4(Num_R)
    end subroutine secondary_reverse_interp_gamma

    subroutine secondary_reverse_interp_tobs(Num_R,R,Tobs_axis,root,tobs_root)
    implicit none
    integer, intent(in) :: Num_R
    integer :: K
    real(8), intent(in) :: R(Num_R),Tobs_axis(Num_R),root
    real(8), intent(out) :: tobs_root
    real(8) :: ratio

        if (root <= R(1)) then
            tobs_root=Tobs_axis(1)
            return
        end if
        do K=1,Num_R-1
            if (root <= R(K+1)) then
                ratio=(root-R(K))/(R(K+1)-R(K))
                tobs_root=Tobs_axis(K)+ratio*(Tobs_axis(K+1)-Tobs_axis(K))
                return
            end if
        end do
        tobs_root=Tobs_axis(Num_R)
    end subroutine secondary_reverse_interp_tobs
end subroutine secondary_reverse_profile

subroutine reverse_dynamics_rhs(dB3,T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross, &
             T,Y,D,M,para_m_ej,V3_scale,Delta_0,eta_0,A_star,dNe_ISM,R_tr,f_jump,f_wide,R0, &
             Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r,sigma_r)
    use constants
    use dynamics_common, only: dynamics_external_density_base, dynamics_external_density_profile, &
                               rs_mag_comp, rs_b4_up, reverse_rhs_phase, &
                               active_density_jump_count, active_density_jump_r, active_density_jump_factor, &
                               active_density_jump_width
    implicit none
    integer, intent(in) :: M
    real(8), intent(inout) :: dB3,T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross
    real(8), intent(in) :: T,para_m_ej,V3_scale,Delta_0,eta_0,A_star,dNe_ISM,R_tr,f_jump,f_wide,R0
    real(8), intent(in) :: Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r,sigma_r
    real(8), intent(in) :: Y(M)
    real(8), intent(out) :: D(M)
    real(8), parameter :: reverse_synch_b_coeff=0.39d0, reverse_gamma_c_precise_coeff=7.739d8
    real(8) :: gam2,RR,para_m2,para_m3,U3,V3,dNe,u2,u4,Delta,para1,para_n4,beta4,beta2,gam34,para_n3,betars
    real(8) :: dB2,gam_c2,gam_m2,eps2,e3,gam_c3,gam_m3,eps3,dgam2_1,dgam2_2,dgam2,dR,dm2,dm3
    real(8) :: thermal_gamma3,ad3,dV3_exp,dV3_shock,dU3_shock,dU3_ad,dU3,dV3
    real(8) :: secondary_m_total,secondary_u_total,secondary_v_total,secondary_p_total,secondary_inertia_mass
    real(8) :: comp_ratio,rho4,B4_ordered,B3_ordered
    logical :: pre_crossing

    call decode_reverse_state()
    D=zero
    pre_crossing=(reverse_rhs_phase == 1 .or. (reverse_rhs_phase == 0 .and. para_m_ej > para_m3))

    call compute_region2_radiative_efficiency()
    e3=U3/V3
    call compute_region3_field()
    call compute_region3_radiative_efficiency()
    call compute_secondary_inertia_mass()

    dgam2_1=-para1*((gam2*gam2-one)*dNe+(gam2*gam34-eta_0)*(beta4-betars)*eta_0*para_n4)
    dgam2_2=(para_m2+para_m3+(one-eps2)*(two*gam2-one)*para_m2+(one-eps3)*(gam34-one)*para_m3+ &
             (one-eps3)*gam2*para_m3*(eta_0*(one-beta4*beta2)-eta_0*beta4/(gam2**2*beta2)) + &
             secondary_inertia_mass)
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
        dgam2_1=-u2**2*dm2/dR
        dgam2_2=para_m_ej+secondary_inertia_mass+(eps2+two*(one-eps2)*gam2)*para_m2
        dgam2=dgam2_1/dgam2_2
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
    D(1:6)=[dgam2,dR,dm2,dm3/para_m_ej,dU3/(para_m_ej*para_c**2),dV3/V3_scale]
    call compute_secondary_branch_derivatives()

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

    ! Secondary reservoirs 随 contact 共动；其 comoving enthalpy 给 Gamma 方程增加有效惯性。
    subroutine compute_secondary_inertia_mass()
    implicit none
    integer :: j, m_idx, u_idx, v_idx

        secondary_m_total=zero; secondary_u_total=zero; secondary_v_total=zero
        secondary_p_total=zero; secondary_inertia_mass=zero
        do j=1,active_density_jump_count
            m_idx=6+j
            u_idx=6+active_density_jump_count+j
            v_idx=6+2*active_density_jump_count+j
            if (v_idx > M) error stop 'secondary RS branch state exceeds reverse dynamics vector'
            secondary_m_total=secondary_m_total+Y(m_idx)*para_m_ej
            secondary_u_total=secondary_u_total+Y(u_idx)*para_m_ej*para_c**2
            secondary_v_total=secondary_v_total+Y(v_idx)*V3_scale
        end do
        if (secondary_v_total > zero) then
            secondary_p_total=secondary_u_total/(3d0*secondary_v_total)
            secondary_inertia_mass=(secondary_m_total*para_c**2+secondary_u_total+ &
                                    secondary_p_total*secondary_v_total)/para_c**2
            secondary_inertia_mass=secondary_inertia_mass+secondary_p_total*secondary_v_total/(gam2*gam2*para_c**2)
        end if
    end subroutine compute_secondary_inertia_mass

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
