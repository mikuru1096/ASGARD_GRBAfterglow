! Secondary reverse-shock profile diagnostics for external density-jump branches.
subroutine secondary_reverse_profile(Num_R,Num_jump,R,Tobs_axis,Gamma4,dNe_ISM,Jump_R,Jump_factor,Jump_width, &
                                     Epsilon_e,Epsilon_b,p_e,f_e,z, &
                                     gamma_contact,pressure_3,gamma_43,comp_ratio,beta_rs,u_diss,active_weight, &
                                     m3_reservoir,u3_reservoir,v3_reservoir,b3_reservoir,gamma_m_shell, &
                                     dissipated_energy,electron_injected_energy,pressure_total,enthalpy_density_total, &
                                     m3_branch,u3_branch,v3_branch,b3_branch,gamma_m_branch, &
                                     nu_m,nu_c,event_active, &
                                     start_radius,shock_end_radius,start_tobs_axis,shock_end_tobs_axis)
    use constants
    use reverse_jump_conditions, only: secondary_reverse_contact_rh
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
