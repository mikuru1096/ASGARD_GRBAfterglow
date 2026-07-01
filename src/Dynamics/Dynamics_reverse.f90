! 反向激波动力学主驱动。
! 顺序: unpack ejecta/medium -> resolve RS start/crossing -> integrate main RS state
!       -> update density-jump secondary RS events -> compute hydro/MHD jump state
!       -> output region-3 thermal, magnetic, and secondary-branch records.
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
    use constants
    use dynamics_common, only: dynamics_rk4_reverse, &
                               dynamics_rk4_reverse_pre_m3, &
                               dynamics_reverse_rhs_iface, dynamics_external_density_profile, &
                               rs_vegas_ud, rs_vegas_comp, rs_mag_specific_internal, &
                               reverse_waiting_phase, reverse_pre_crossing_phase, &
                               dynamics_set_density_jump_profile, &
                               active_density_jump_count, density_jump_max, active_density_jump_r, &
                               active_density_jump_factor, active_density_jump_width
    use reverse_jump_conditions, only: secondary_reverse_contact_rh
    implicit none
    integer, intent(in) :: n,Num_R
    integer, parameter :: density_boundary_r0_index = 27
    integer :: I_tobs, Num_R1, Num_state, j_event
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
    real(8) :: Delta_0,para_m_ej,V3_scale,para_m2,para_m3,DM_0,R_dec,T00,t_dec,Grid_Tobs_bin,T_log10,dB3
    real(8) :: dT_grid,R_dec_ism,R_dec_wind
    real(8) :: T_state,T_target,T_pressure_event,dB3_wait_trial,dB3_pressure_event
    real(8) :: u2_init,u4_init,Delta_init,para_n4_init,gam34_init,para_n3_init,comp_init
    real(8) :: event_prev_radius,event_prev_gamma,event_prev_tobs,event_curr_radius,event_curr_gamma,event_curr_tobs
    real(8) :: event_prev_source(density_jump_max),event_curr_source(density_jump_max)
    real(8) :: secondary_dissipated_prev(density_jump_max),secondary_gammam_prev(density_jump_max)
    real(8),allocatable :: Y(:),event_prev_state(:),event_curr_state(:)
    real(8),allocatable :: wait_trial_state(:),pressure_event_state(:)
    logical :: Secondary_event_closed(density_jump_max), pressure_ready_trial

    Eta_0=Boundary(1); R(1)=Boundary(4); Epsilon_e=Boundary(5); Epsilon_b=Boundary(6); p_f=Boundary(7); z=Boundary(8)
    dNe_ISM=Boundary(11); A_star=Boundary(12); E_iso=Boundary(14); T_log10_duration=Boundary(15); f_e=Boundary(16)
    R_tr=Boundary(21); f_jump=Boundary(22); f_wide=Boundary(23)
    if (n >= density_boundary_r0_index) then
        R0=Boundary(density_boundary_r0_index)
    else
        R0=Boundary(n)
    end if
    call dynamics_set_density_jump_profile(Boundary,n)
    Num_state=6+5*active_density_jump_count
    allocate(Y(Num_state),event_prev_state(Num_state),event_curr_state(Num_state), &
             wait_trial_state(Num_state),pressure_event_state(Num_state))
    Y=zero
    Delta_0=Delta_t*para_c
    if (sigma_r <= zero) then
        para_m_ej=E_iso/eta_0/para_c**2
    else
        para_m_ej=E_iso/(one+sigma_r)/eta_0/para_c**2
    end if

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
    comp_init=rs_vegas_comp(gam34_init,sigma_r)
    para_n3_init=comp_init*para_n4_init
    V3_scale=para_m_ej/(para_n3_init*Para_m_p)
    Y(1:6)=[R_Gamma(1),R(1),para_m2,para_m3/para_m_ej,rs_mag_specific_internal(gam34_init,sigma_r)*para_m3/para_m_ej, &
            para_m3/(para_n3_init*Para_m_p)/V3_scale]

    DM_0=E_iso/((Eta_0-one)*4d0*pi*Para_m_p_E)
    R_dec_ism=(dNe_ISM*Eta_0/DM_0)**(-one/3d0)
    if (A_star > zero) then
        R_dec_wind=DM_0/(2.0d35*A_star*Eta_0)
        R_dec=min(R_dec_wind,R_dec_ism)
    else
        R_dec=R_dec_ism
    end if
    T_cross=-1d0; R_cross=zero; e3_cross=zero; gam20=one
    U3_cross=zero; V3_cross=zero; M3_cross=zero; gam_m_cross=zero
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
    event_prev_source=-one
    do j_event=1,active_density_jump_count
        call secondary_reverse_event_source(j_event,event_prev_radius,event_prev_gamma,event_prev_state,event_prev_source(j_event))
    end do

    do I_tobs=1,Num_R
        dT_grid=ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        T_target=T00+dT_grid
        call advance_reverse_state_to_target(T_target)
        R_Tobs(I_tobs)=T_target*(one+z); R_Gamma(I_tobs)=Y(1); R(I_tobs)=Y(2); M2(I_tobs)=Y(3)
        event_curr_radius=Y(2); event_curr_gamma=Y(1); event_curr_tobs=R_Tobs(I_tobs)
        event_curr_state=Y
        event_curr_source=-one
        do j_event=1,active_density_jump_count
            call secondary_reverse_event_source(j_event,event_curr_radius,event_curr_gamma, &
                                                event_curr_state,event_curr_source(j_event))
        end do
        call update_secondary_reverse_events()
        event_prev_radius=event_curr_radius; event_prev_gamma=event_curr_gamma
        event_prev_tobs=event_curr_tobs; event_prev_source=event_curr_source
        event_prev_state=event_curr_state
        M3(I_tobs)=Y(4)*para_m_ej; U3_th(I_tobs)=Y(5)*para_m_ej*para_c**2
        V3_comoving(I_tobs)=Y(6)*V3_scale; B3(I_tobs)=dB3
        call store_secondary_branch_state(I_tobs)
        Gamma34_inst(I_tobs)=(Y(1)*Y(1)+Eta_0*Eta_0-one)/(Eta_0*Y(1)+dsqrt(Y(1)*Y(1)-one)*u4_init)
    end do

    do j_event=1,active_density_jump_count
        if (Secondary_event_active(j_event) .and. .not. Secondary_event_closed(j_event)) then
            Secondary_event_closed(j_event)=.true.
            Secondary_end_radius(j_event)=event_prev_radius
            Secondary_end_tobs_axis(j_event)=event_prev_tobs
        end if
    end do
    deallocate(Y,event_prev_state,event_curr_state,wait_trial_state,pressure_event_state)

contains

    subroutine advance_reverse_state_to_target(T_target_in)
    implicit none
    real(8), intent(in) :: T_target_in
    real(8) :: T_local,H_local

        do while (T_state < T_target_in)
            if (T_cross < zero .and. Y(4) < one) then
                if (reverse_shock_pressure_ready_state(Y)) then
                    call dynamics_rk4_reverse_pre_m3(reverse_dynamics_rhs,dB3,T_cross,R_cross,e3_cross,gam20, &
                                                     U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross, &
                                                     T_state,T_target_in,Y,Num_state,para_m_ej,V3_scale,Delta_0,eta_0, &
                                                     A_star,dNe_ISM,R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e,p_f, &
                                                     f_e,e_r,b_r,p_r,f_e_r,sigma_r)
                else
                    call waiting_trial(T_target_in,wait_trial_state,dB3_wait_trial)
                    pressure_ready_trial=reverse_shock_pressure_ready_state(wait_trial_state)
                    if (.not. pressure_ready_trial) then
                        Y=wait_trial_state
                        dB3=dB3_wait_trial
                        T_state=T_target_in
                    else
                        call locate_waiting_event(T_target_in,T_pressure_event,pressure_event_state,dB3_pressure_event)
                        Y=pressure_event_state
                        dB3=dB3_pressure_event
                        T_state=T_pressure_event
                    end if
                end if
            else
                H_local=T_target_in-T_state
                T_local=T_state
                call dynamics_rk4_reverse(reverse_dynamics_rhs,dB3,T_cross,R_cross,e3_cross,gam20, &
                                          U3_cross,V3_cross,M3_cross, &
                                          gam_m_cross,B3_ordered_cross,T_local,H_local,Y,Num_state,para_m_ej,V3_scale,Delta_0, &
                                          eta_0,A_star,dNe_ISM,R_tr,f_jump,f_wide,R0, &
                                          Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r,sigma_r,reverse_pre_crossing_phase)
                T_state=T_target_in
            end if
        end do
    end subroutine advance_reverse_state_to_target

    logical function reverse_shock_pressure_ready_state(state)
    implicit none
    real(8), intent(in) :: state(Num_state)
    real(8) :: radius,gamma_cd,u_cd,u4,gamma43,delta_shell,n1,n4
    real(8) :: u_down,u_up,sigma_cr
    logical :: pressure_ready,fast_ready

        if (sigma_r <= zero) then
            reverse_shock_pressure_ready_state=.true.
            return
        end if
        radius=state(2)
        gamma_cd=state(1)
        u_cd=dsqrt(gamma_cd*gamma_cd-one)
        u4=dsqrt(Eta_0*Eta_0-one)
        gamma43=(gamma_cd*gamma_cd+Eta_0*Eta_0-one)/(Eta_0*gamma_cd+u_cd*u4)
        delta_shell=max(Delta_0,radius/(Eta_0*Eta_0))
        call dynamics_external_density_profile(A_star,dNe_ISM,radius,R0,1,R_tr,f_jump,f_wide,n1)
        n4=para_m_ej/(4d0*pi*Para_m_p*radius*radius*Eta_0*delta_shell)
        sigma_cr=two*(4d0*gamma_cd+3d0)*(gamma_cd-one)*n1/(3d0*n4)
        pressure_ready=(sigma_cr >= sigma_r)
        u_down=rs_vegas_ud(gamma43,sigma_r)
        u_up=dsqrt((one+u_down*u_down)*(gamma43-one)*(gamma43+one))+u_down*gamma43
        fast_ready=(u_up*u_up > sigma_r)
        reverse_shock_pressure_ready_state=pressure_ready .and. fast_ready
    end function reverse_shock_pressure_ready_state

    subroutine waiting_trial(T_end,state_out,dB3_out)
    implicit none
    real(8), intent(in) :: T_end
    real(8), intent(out) :: state_out(Num_state), dB3_out
    real(8) :: T_local,H_local,T_cross_trial,R_cross_trial,e3_cross_trial,gam20_trial
    real(8) :: U3_cross_trial,V3_cross_trial,M3_cross_trial,gam_m_cross_trial,B3_ordered_cross_trial

        state_out=Y
        dB3_out=dB3
        T_cross_trial=T_cross; R_cross_trial=R_cross; e3_cross_trial=e3_cross; gam20_trial=gam20
        U3_cross_trial=U3_cross; V3_cross_trial=V3_cross; M3_cross_trial=M3_cross
        gam_m_cross_trial=gam_m_cross; B3_ordered_cross_trial=B3_ordered_cross
        T_local=T_state
        H_local=T_end-T_state
        call dynamics_rk4_reverse(reverse_dynamics_rhs,dB3_out,T_cross_trial,R_cross_trial,e3_cross_trial, &
                                  gam20_trial,U3_cross_trial,V3_cross_trial,M3_cross_trial, &
                                  gam_m_cross_trial,B3_ordered_cross_trial,T_local,H_local,state_out, &
                                  Num_state,para_m_ej,V3_scale,Delta_0,eta_0,A_star,dNe_ISM, &
                                  R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e,p_f,f_e,e_r,b_r,p_r,f_e_r, &
                                  sigma_r,reverse_waiting_phase)
    end subroutine waiting_trial

    subroutine locate_waiting_event(T_end,event_T,event_state,dB3_event)
    implicit none
    integer :: iter
    real(8), intent(in) :: T_end
    real(8), intent(out) :: event_T,event_state(Num_state),dB3_event
    real(8) :: T_lo,T_hi,T_mid

        T_lo=T_state
        T_hi=T_end
        event_state=wait_trial_state
        dB3_event=dB3_wait_trial
        do iter=1,32
            T_mid=0.5d0*(T_lo+T_hi)
            call waiting_trial(T_mid,wait_trial_state,dB3_wait_trial)
            if (reverse_shock_pressure_ready_state(wait_trial_state)) then
                T_hi=T_mid
                event_state=wait_trial_state
                dB3_event=dB3_wait_trial
            else
                T_lo=T_mid
            end if
        end do
        event_T=T_hi
    end subroutine locate_waiting_event

    ! Secondary RS 事件诊断：用动力学积分接受态的 mechanical source 符号定位 start/end。
    subroutine update_secondary_reverse_events()
    implicit none
    integer :: j,k,n_scan
    real(8) :: radius_lo,radius_hi,gamma_lo,gamma_hi,tobs_lo,tobs_hi,source_lo,source_hi
    real(8) :: width,overlap_lo,overlap_hi,scan_length
    real(8) :: root_radius,root_tobs
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
                radius_lo=radius_hi; gamma_lo=gamma_hi; tobs_lo=tobs_hi; source_lo=source_hi; state_lo=state_hi
            end do
        end do
    end subroutine update_secondary_reverse_events

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
    integer :: j, m_idx, u_idx, v_idx, e_idx, g_idx, parent_idx
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
                        parent_idx = j-1
                        if (secondary_parent_upstream_available(j,Secondary_M3(parent_idx,i_out)/para_m_ej, &
                                                                Secondary_U3(parent_idx,i_out)/(para_m_ej*para_c**2), &
                                                                Secondary_V3(parent_idx,i_out)/V3_scale)) then
                            n4=Secondary_M3(parent_idx,i_out)/(Para_m_p*Secondary_V3(parent_idx,i_out))
                            e4=Secondary_U3(parent_idx,i_out)/Secondary_V3(parent_idx,i_out)
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
