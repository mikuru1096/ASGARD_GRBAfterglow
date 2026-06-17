module electron_reverse_kernel
    use constants
    use dynamics_common, only: dynamics_external_density_profile, dynamics_reverse_gamma_extrema
    use electron_injection_profiles, only: electron_exp_cutoff_factor, electron_profile_log_cell_edges
    use electron_transport_common, only: electron_fullhide_flux_split_step, electron_dnx_to_dndgamma_exp_centers
    use electron_radiation_kernel, only: get_syn_selected, get_nu_a
    use electron_cooling_kernel, only: prepare_forward_cooling_aux
    implicit none
contains

! 反向激波电子演化主驱动：注入→同步+IC冷却→隐式输运推进，IC模式遵循公开标准选择。
subroutine electron_reverse_evolve(Delta_0,e_r,b_r,p_r,f_e_r,eta_0,Epsilon_e,Epsilon_b,z,A_star,dNe_ISM,para_m_ej, &
                                   R_tr,f_jump,f_wide,R0, &
                                   T_cross,R_cross,U3_cross,M3_cross,R_Tobs,R_Gamma,R,B3,M3_shell,U3_shell,V_seed, &
                                   Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads,gam_e,dN_gam_e)
    implicit none
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer :: I_tobs,I_gam_e,L1,L
    real(8), intent(in) :: Delta_0,e_r,b_r,p_r,f_e_r,eta_0,Epsilon_e,Epsilon_b,z,A_star,dNe_ISM,para_m_ej
    real(8), intent(in) :: R_tr,f_jump,f_wide,R0
    real(8), intent(in) :: T_cross,R_cross,U3_cross,M3_cross
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),B3(Num_R),M3_shell(Num_R),U3_shell(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R)
    real(8), parameter :: reverse_gamma_c_coeff=7.7d8, reverse_synch_b_coeff=0.39d0, reverse_adv_coeff=1.35d-19
    real(8) :: factor2,dB,gamma34,Gam_e_max,Gam_e_m,Gam_e_c,dNe,DB_min,Gam_e_max_max,Gam_e_min_global,d_x,R_loc,R_Gamma_loc,Delta
    real(8) :: R_n4,beta4,beta2,u2,u4,f_r,dDR,dDD,Qshell,cooling_scale
    real(8) :: thermal_scale_lo,thermal_scale_hi,thermal_loss_rate,adiabatic_rate
    real(8) :: injection_rate,inj_hi,inj_width,mass_lo,mass_hi
    real(8), allocatable :: dEl(:),x(:),dF1(:),temp3(:),dN_x(:),x_edge(:)
    real(8), allocatable :: dB3_serial(:),P_syn(:),Seed_syn(:),cooling_aux(:),Compton(:)

    allocate(dEl(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e),temp3(Num_gam_e-1),dN_x(Num_gam_e),x_edge(Num_gam_e+1), &
             dB3_serial(Num_R),P_syn(Num_nu),Seed_syn(Num_nu), &
             cooling_aux(Num_gam_e),Compton(Num_gam_e))
    dB3_serial=B3

    factor2=(p_r-two)/(p_r-one)*e_r*Para_m_p_div_m_e
    if (p_r < 2.05d0) factor2=0.05d0/1.05d0*e_r*Para_m_p_div_m_e
    beta4=dsqrt(one-one/eta_0**2); u4=dsqrt(eta_0*eta_0-one)
    dB3_serial(1)=dB3_serial(min(2,Num_R))
    dB=dB3_serial(1); gamma34=1.001d0
    call dynamics_reverse_gamma_extrema(dB,gamma34,factor2,f_e_r,Gam_e_max,Gam_e_m)
    Gam_e_c=reverse_gamma_c_coeff/(one+dsqrt(e_r/b_r))/R_Gamma(1)/dB**2/(R_Tobs(1)/two)

    call dynamics_external_density_profile(A_star,dNe_ISM,R(1),R0,1,R_tr,f_jump,f_wide,dNe)
    DB_min=reverse_synch_b_coeff*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    Gam_e_min_global=one

    do I_tobs=2,Num_R
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        beta2=dsqrt(one-one/R_Gamma_loc**2)
        u2=dsqrt(R_Gamma_loc*R_Gamma_loc-one)
        gamma34=(R_Gamma_loc*R_Gamma_loc+eta_0*eta_0-one)/(eta_0*R_Gamma_loc+u2*u4)
        dB=(dB3_serial(I_tobs)+dB3_serial(I_tobs-1))/two
        call dynamics_reverse_gamma_extrema(dB,gamma34,factor2,f_e_r,Gam_e_max,Gam_e_m)
    end do
    if (Gam_e_max_max <= Gam_e_min_global) error stop "electron_reverse_evolve: reverse electron grid maximum must exceed minimum."

    do I_gam_e=1,Num_gam_e
        if (Num_gam_e == 1) then
            gam_e(I_gam_e)=Gam_e_min_global
        else
            gam_e(I_gam_e)=Gam_e_min_global*ten**(dlog10(Gam_e_max_max/Gam_e_min_global)*(I_gam_e-1)/(Num_gam_e-1))
        end if
        dN_gam_e(I_gam_e,1)=zero
    end do

    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
    d_x=dlog10(gam_e(2)/gam_e(1))
    call electron_profile_log_cell_edges(Num_gam_e,gam_e,x_edge)

    do I_tobs=2,Num_R
        call prepare_reverse_shell_state(I_tobs)
        call compute_reverse_thermal_loss(I_tobs)
        call compute_reverse_injection_rate(I_tobs)
        L1=max(100,min(1000,int(dDD/dDR)))
        dDR=dDD/L1
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)

        call compute_reverse_cooling_loss(I_tobs)
        call advance_reverse_transport_shell(I_tobs)
    end do

    deallocate(dEl,x,dF1,temp3,dN_x,x_edge,dB3_serial,P_syn,Seed_syn,cooling_aux,Compton)

contains

    subroutine prepare_reverse_shell_state(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs

        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        Delta=max(Delta_0,R_loc/Eta_0**2)
        R_n4=para_m_ej/(4d0*pi*Para_m_p*R_loc*R_loc*Eta_0*Delta)
        beta2=dsqrt(one-one/R_Gamma_loc**2)
        u2=dsqrt(R_Gamma_loc*R_Gamma_loc-one)
        gamma34=(R_Gamma_loc*R_Gamma_loc+eta_0*eta_0-one)/(eta_0*R_Gamma_loc+u2*u4)
        dB=(dB3_serial(I_tobs)+dB3_serial(I_tobs-1))/two
        call dynamics_reverse_gamma_extrema(dB,gamma34,factor2,f_e_r,Gam_e_max,Gam_e_m)
        Gam_e_c=reverse_gamma_c_coeff*(one+z)/R_Gamma_loc/dB**2/R_Tobs(I_tobs)
        f_r=reverse_adv_coeff/beta2/R_Gamma_loc*dB**2/pi
        dDR=0.7d0/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)

    end subroutine prepare_reverse_shell_state

    subroutine compute_reverse_thermal_loss(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs

        thermal_loss_rate=zero
        if (R(I_tobs) > R_cross) then
            if (R(I_tobs-1) < R_cross) then
                thermal_scale_lo=U3_cross/M3_cross
                thermal_scale_hi=U3_shell(I_tobs)/M3_shell(I_tobs)
            else
                thermal_scale_lo=U3_shell(I_tobs-1)/M3_shell(I_tobs-1)
                thermal_scale_hi=U3_shell(I_tobs)/M3_shell(I_tobs)
            end if
            if (thermal_scale_lo <= zero .or. thermal_scale_hi <= zero) &
                error stop "electron_reverse_evolve: post-crossing thermal scale must be positive."
            if (R(I_tobs-1) < R_cross) then
                thermal_loss_rate=-dlog(thermal_scale_hi/thermal_scale_lo)/(R(I_tobs)-R_cross)
            else
                thermal_loss_rate=-dlog(thermal_scale_hi/thermal_scale_lo)/dDD
            end if
            if (thermal_loss_rate <= zero) &
                error stop "electron_reverse_evolve: post-crossing thermal scale must decrease."
        end if
    end subroutine compute_reverse_thermal_loss

    subroutine compute_reverse_injection_rate(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs

        injection_rate=zero
        if (R(I_tobs-1) < R_cross) then
            inj_hi=min(R(I_tobs),R_cross)
            if (inj_hi > R(I_tobs-1)) then
                mass_lo=M3_shell(I_tobs-1)
                if (I_tobs == 2) mass_lo=zero
                mass_hi=M3_shell(I_tobs)
                if (R(I_tobs) > R_cross) mass_hi=M3_cross
                if (mass_hi < mass_lo) error stop "electron_reverse_evolve: reverse swept mass must not decrease."
                inj_width=inj_hi-R(I_tobs-1)
                injection_rate=f_e_r*(mass_hi-mass_lo)/(Para_m_p*inj_width)
            end if
        end if
    end subroutine compute_reverse_injection_rate

    subroutine compute_reverse_cooling_loss(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs

        call get_syn_selected(index_syn_intger,R(I_tobs-1),dB,Num_gam_e,Num_nu,n_threads, &
                              gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_syn,Seed_syn)
        call prepare_forward_cooling_aux(index_Y,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,Seed_syn,cooling_aux)
        select case(index_Y)
        case(0)
            dEl=f_r*gam_e
        case(1)
            cooling_scale=one/(beta2*R_Gamma_loc*Para_c)
            dEl=(f_r+cooling_aux*cooling_scale)*gam_e
        case(2)
            Qshell=4d0*pi*R(I_tobs-1)*R(I_tobs-1)*Para_c
            Compton=one+cooling_aux/Qshell/(4d0*R_Gamma_loc*R_Gamma_loc*R_n4*Para_m_p_E)
            Gam_e_max=Gam_e_max/dsqrt(Compton(Num_gam_e))
            dEl=f_r*Compton*gam_e
        case default
            error stop 'electron_reverse_evolve: index_Y must be 0, 1, or 2.'
        end select
    end subroutine compute_reverse_cooling_loss

    subroutine advance_reverse_transport_shell(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs

        do L=1,L1
            R_loc=R_loc+dDR
            if (R_cross >= R_loc) then
                call reverse_build_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m,Gam_e_max,injection_rate,p_r,dF1)
            else
                dF1=zero
            end if
            if (R_loc <= R_cross) then
                adiabatic_rate=one/R_loc
            else
                adiabatic_rate=thermal_loss_rate
            end if
            temp3=((dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two+adiabatic_rate)/dlog(ten)
            call electron_fullhide_flux_split_step(Num_gam_e,dDR,d_x,temp3,dF1,dN_x,x,.true.)
            dN_x=x
            if (L == L1) call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
        end do
    end subroutine advance_reverse_transport_shell
end subroutine electron_reverse_evolve

subroutine electron_secondary_reverse_evolve(e_r,b_r,p_r,f_e_r,z,R_Tobs,R_Gamma,R,B3,M3_shell,U3_shell,V3_shell,Gam_m_shell, &
                                             V_seed,Num_nu,Num_R,Num_gam_e,index_syn_intger,n_threads,gam_e,dN_gam_e)
    implicit none
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,index_syn_intger,n_threads
    integer :: I_tobs,I_gam_e,L1,L
    real(8), intent(in) :: e_r,b_r,p_r,f_e_r,z
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),B3(Num_R),M3_shell(Num_R),U3_shell(Num_R)
    real(8), intent(in) :: V3_shell(Num_R),Gam_m_shell(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R)
    real(8), parameter :: secondary_adv_coeff=1.35d-19
    real(8) :: dB,Gam_e_max,Gam_e_m,Gam_e_max_max,Gam_e_min_global,d_x,R_loc,R_Gamma_loc,beta2
    real(8) :: f_r,dDR,dDD,injection_rate,mass_lo,mass_hi,inj_width,adiabatic_rate
    real(8), allocatable :: dEl(:),x(:),dF1(:),temp3(:),dN_x(:),x_edge(:)

    allocate(dEl(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e),temp3(Num_gam_e-1),dN_x(Num_gam_e),x_edge(Num_gam_e+1))

    Gam_e_min_global=one
    Gam_e_max_max=zero
    do I_tobs=2,Num_R
        dB=(B3(I_tobs)+B3(I_tobs-1))/two
        if (dB > zero .and. Gam_m_shell(I_tobs) > one) then
            Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
            Gam_e_max_max=max(Gam_e_max_max,Gam_e_max)
        end if
    end do
    if (Gam_e_max_max <= Gam_e_min_global) error stop "electron_secondary_reverse_evolve: empty secondary electron grid."

    do I_gam_e=1,Num_gam_e
        if (Num_gam_e == 1) then
            gam_e(I_gam_e)=Gam_e_min_global
        else
            gam_e(I_gam_e)=Gam_e_min_global*ten**(dlog10(Gam_e_max_max/Gam_e_min_global)*(I_gam_e-1)/(Num_gam_e-1))
        end if
        dN_gam_e(I_gam_e,1)=zero
    end do

    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
    d_x=dlog10(gam_e(2)/gam_e(1))
    call electron_profile_log_cell_edges(Num_gam_e,gam_e,x_edge)

    do I_tobs=2,Num_R
        if (M3_shell(I_tobs) <= zero .and. M3_shell(I_tobs-1) <= zero) then
            dN_gam_e(:,I_tobs)=dN_gam_e(:,I_tobs-1)
            cycle
        end if
        call prepare_secondary_shell_state(I_tobs)
        call compute_secondary_injection_rate(I_tobs)
        call compute_secondary_adiabatic_rate(I_tobs)
        call advance_secondary_transport_shell(I_tobs)
    end do

    deallocate(dEl,x,dF1,temp3,dN_x,x_edge)

contains

    subroutine prepare_secondary_shell_state(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs

        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        beta2=dsqrt(one-one/R_Gamma_loc**2)
        dB=(B3(I_tobs)+B3(I_tobs-1))/two
        if (dB <= zero) error stop "electron_secondary_reverse_evolve: secondary reservoir requires B3 > 0."
        if (Gam_m_shell(I_tobs-1) > one) then
            Gam_e_m=(Gam_m_shell(I_tobs)+Gam_m_shell(I_tobs-1))/two
        else
            Gam_e_m=Gam_m_shell(I_tobs)
        end if
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        f_r=secondary_adv_coeff/beta2/R_Gamma_loc*dB**2/pi
        dDD=R(I_tobs)-R(I_tobs-1)
        dDR=0.7d0/(f_r*Gam_e_max+one/R(I_tobs-1))
        L1=max(100,min(1000,int(dDD/dDR)))
        dDR=dDD/L1
        dEl=f_r*gam_e
    end subroutine prepare_secondary_shell_state

    subroutine compute_secondary_injection_rate(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs

        mass_lo=M3_shell(I_tobs-1)
        mass_hi=M3_shell(I_tobs)
        if (mass_hi < mass_lo) error stop "electron_secondary_reverse_evolve: secondary swept mass must not decrease."
        inj_width=R(I_tobs)-R(I_tobs-1)
        injection_rate=f_e_r*(mass_hi-mass_lo)/(Para_m_p*inj_width)
    end subroutine compute_secondary_injection_rate

    subroutine compute_secondary_adiabatic_rate(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs

        if (V3_shell(I_tobs) <= zero .or. V3_shell(I_tobs-1) <= zero) then
            adiabatic_rate=zero
        else
            adiabatic_rate=dlog(V3_shell(I_tobs)/V3_shell(I_tobs-1))/(3d0*(R(I_tobs)-R(I_tobs-1)))
        end if
    end subroutine compute_secondary_adiabatic_rate

    subroutine advance_secondary_transport_shell(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs

        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)
        do L=1,L1
            R_loc=R_loc+dDR
            if (injection_rate > zero) then
                call reverse_build_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m,Gam_e_max,injection_rate,p_r,dF1)
            else
                dF1=zero
            end if
            temp3=((dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two+adiabatic_rate)/dlog(ten)
            call electron_fullhide_flux_split_step(Num_gam_e,dDR,d_x,temp3,dF1,dN_x,x,.true.)
            dN_x=x
            if (L == L1) call electron_dnx_to_dndgamma_exp_centers(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
        end do
    end subroutine advance_secondary_transport_shell
end subroutine electron_secondary_reverse_evolve

! Secondary RS 源项历史：逐半径壳层调用同步辐射和SSA核，返回给统一EATS投影使用。
subroutine electron_secondary_reverse_synchrotron(index_syn_intger,Num_nu,Num_R,Num_gam_e,n_threads,R,R_Gamma,B3,gam_e, &
                                                  dN_gam_e,V_seed,z,L_syn_spec,Seed_syn,Nu_a)
    implicit none
    integer, intent(in) :: index_syn_intger,Num_nu,Num_R,Num_gam_e,n_threads
    integer :: I_tobs
    real(8), intent(in) :: R(Num_R),R_Gamma(Num_R),B3(Num_R),gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R)
    real(8), intent(in) :: V_seed(Num_nu),z
    real(8), intent(out) :: L_syn_spec(Num_nu,Num_R),Seed_syn(Num_nu,Num_R),Nu_a(Num_R)
    real(8) :: doppler_den

    L_syn_spec=zero
    Seed_syn=zero
    Nu_a=zero
    do I_tobs=1,Num_R
        if (B3(I_tobs) <= zero) cycle
        call get_syn_selected(index_syn_intger,R(I_tobs),B3(I_tobs),Num_gam_e,Num_nu,n_threads, &
                              gam_e,dN_gam_e(:,I_tobs),V_seed,L_syn_spec(:,I_tobs),Seed_syn(:,I_tobs))
        doppler_den=R_Gamma(I_tobs)*(one-dsqrt(one-R_Gamma(I_tobs)**(-2)))*(one+z)
        call get_nu_a(R(I_tobs),B3(I_tobs),Num_gam_e,gam_e,dN_gam_e(:,I_tobs),Nu_a(I_tobs))
        Nu_a(I_tobs)=Nu_a(I_tobs)/doppler_den
    end do
end subroutine electron_secondary_reverse_synchrotron

! Secondary RS 分支辐射归并：每个 reservoir 独立输运电子谱，再在源项半径网格上求和。
subroutine electron_secondary_reverse_branch_synchrotron(e_r,b_r,p_r,f_e_r,z,R_Tobs,R_Gamma,R,B3_branch, &
                                                         M3_branch,U3_branch,V3_branch,Gam_m_branch,V_seed, &
                                                         Num_jump,Num_nu,Num_R,Num_gam_e,index_syn_intger,n_threads, &
                                                         Branch_L_syn_spec,L_syn_spec)
    implicit none
    integer, intent(in) :: Num_jump,Num_nu,Num_R,Num_gam_e,index_syn_intger,n_threads
    integer :: I_jump
    real(8), intent(in) :: e_r,b_r,p_r,f_e_r,z
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: B3_branch(Num_jump,Num_R),M3_branch(Num_jump,Num_R),U3_branch(Num_jump,Num_R)
    real(8), intent(in) :: V3_branch(Num_jump,Num_R),Gam_m_branch(Num_jump,Num_R)
    real(8), intent(out) :: Branch_L_syn_spec(Num_jump,Num_nu,Num_R),L_syn_spec(Num_nu,Num_R)
    real(8), allocatable :: gam_e_branch(:),dN_branch(:,:),seed_dummy(:,:),nu_a_dummy(:)

    allocate(gam_e_branch(Num_gam_e),dN_branch(Num_gam_e,Num_R),seed_dummy(Num_nu,Num_R),nu_a_dummy(Num_R))
    Branch_L_syn_spec=zero
    L_syn_spec=zero
    do I_jump=1,Num_jump
        if (.not. any(M3_branch(I_jump,:) > zero)) cycle
        call electron_secondary_reverse_evolve(e_r,b_r,p_r,f_e_r,z,R_Tobs,R_Gamma,R,B3_branch(I_jump,:), &
                                               M3_branch(I_jump,:),U3_branch(I_jump,:),V3_branch(I_jump,:), &
                                               Gam_m_branch(I_jump,:),V_seed,Num_nu,Num_R,Num_gam_e, &
                                               index_syn_intger,n_threads,gam_e_branch,dN_branch)
        call electron_secondary_reverse_synchrotron(index_syn_intger,Num_nu,Num_R,Num_gam_e,n_threads,R,R_Gamma, &
                                                    B3_branch(I_jump,:),gam_e_branch,dN_branch,V_seed,z, &
                                                    Branch_L_syn_spec(I_jump,:,:),seed_dummy,nu_a_dummy)
        L_syn_spec=L_syn_spec+Branch_L_syn_spec(I_jump,:,:)
    end do
    deallocate(gam_e_branch,dN_branch,seed_dummy,nu_a_dummy)
end subroutine electron_secondary_reverse_branch_synchrotron

subroutine electron_secondary_reverse_branch_reaccelerated(e_r,b_r,p_r,f_e_r,z,R_Tobs,R_Gamma,R,B3_branch, &
                                                           M3_branch,U3_branch,V3_branch,Gam_m_branch,Gamma43_branch, &
                                                           Comp_branch,Parent_branch,V_seed,Num_jump,Num_nu,Num_R,Num_gam_e, &
                                                           index_syn_intger,n_threads,gam_e,dN_total,Branch_L_syn_spec, &
                                                           L_syn_spec,Branch_seed_energy,Branch_reaccel_energy)
    implicit none
    integer, intent(in) :: Num_jump,Num_nu,Num_R,Num_gam_e,index_syn_intger,n_threads
    integer, intent(in) :: Parent_branch(Num_jump)
    integer :: I_tobs,I_jump,I_gam_e,L1,L,parent
    real(8), intent(in) :: e_r,b_r,p_r,f_e_r,z
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: B3_branch(Num_jump,Num_R),M3_branch(Num_jump,Num_R),U3_branch(Num_jump,Num_R)
    real(8), intent(in) :: V3_branch(Num_jump,Num_R),Gam_m_branch(Num_jump,Num_R),Gamma43_branch(Num_jump,Num_R)
    real(8), intent(in) :: Comp_branch(Num_jump,Num_R)
    real(8), intent(out) :: gam_e(Num_gam_e),dN_total(Num_gam_e,Num_R)
    real(8), intent(out) :: Branch_L_syn_spec(Num_jump,Num_nu,Num_R),L_syn_spec(Num_nu,Num_R)
    real(8), intent(out) :: Branch_seed_energy(Num_jump,Num_R),Branch_reaccel_energy(Num_jump,Num_R)
    real(8), parameter :: secondary_adv_coeff=1.35d-19
    real(8) :: dB,Gam_e_max,Gam_e_m,Gam_e_max_max,Gam_e_min_global,d_x,R_loc,R_Gamma_loc,beta2
    real(8) :: f_r,dDR,dDD,injection_rate,mass_lo,mass_hi,mass_delta,fresh_mass,adiabatic_rate
    real(8) :: parent_mass,transfer_mass,transfer_fraction,boost_factor,seed_energy,out_energy
    logical :: dissipative_shell
    real(8), allocatable :: dEl(:),x(:),dF1(:),temp3(:),x_edge(:),dN_x(:,:),dN_work(:,:)
    real(8), allocatable :: dN_seed(:),dN_boost(:),dN_reaccel(:),dN_gamma_branch(:,:,:),seed_dummy(:,:),nu_a_dummy(:)
    real(8), allocatable :: branch_mass_available(:),fresh_mass_branch(:)

    allocate(dEl(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e),temp3(Num_gam_e-1),x_edge(Num_gam_e+1), &
             dN_x(Num_jump,Num_gam_e),dN_work(Num_jump,Num_gam_e),dN_seed(Num_gam_e),dN_boost(Num_gam_e), &
             dN_reaccel(Num_gam_e),dN_gamma_branch(Num_jump,Num_gam_e,Num_R),seed_dummy(Num_nu,Num_R), &
             nu_a_dummy(Num_R),branch_mass_available(Num_jump),fresh_mass_branch(Num_jump))

    call build_secondary_reaccel_gamma_grid()
    call electron_profile_log_cell_edges(Num_gam_e,gam_e,x_edge)
    d_x=dlog10(gam_e(2)/gam_e(1))
    dN_gamma_branch=zero; dN_total=zero; dN_x=zero; branch_mass_available=zero
    fresh_mass_branch=zero
    Branch_L_syn_spec=zero; L_syn_spec=zero; Branch_seed_energy=zero; Branch_reaccel_energy=zero

    do I_tobs=2,Num_R
        dN_work=dN_gamma_branch(:,:,I_tobs-1)*spread(gam_e*dlog(ten),1,Num_jump)
        call transfer_reaccelerated_parent_electrons(I_tobs)
        do I_jump=1,Num_jump
            call advance_reaccelerated_branch_shell(I_tobs,I_jump)
        end do
    end do
    do I_jump=1,Num_jump
        dN_total=dN_total+dN_gamma_branch(I_jump,:,:)
        if (.not. any(dN_gamma_branch(I_jump,:,:) > zero)) cycle
        call electron_secondary_reverse_synchrotron(index_syn_intger,Num_nu,Num_R,Num_gam_e,n_threads,R,R_Gamma, &
                                                    B3_branch(I_jump,:),gam_e,dN_gamma_branch(I_jump,:,:),V_seed,z, &
                                                    Branch_L_syn_spec(I_jump,:,:),seed_dummy,nu_a_dummy)
        L_syn_spec=L_syn_spec+Branch_L_syn_spec(I_jump,:,:)
    end do
    deallocate(dEl,x,dF1,temp3,x_edge,dN_x,dN_work,dN_seed,dN_boost,dN_reaccel,dN_gamma_branch, &
               seed_dummy,nu_a_dummy,branch_mass_available,fresh_mass_branch)

contains

    subroutine build_secondary_reaccel_gamma_grid()
    implicit none

        Gam_e_min_global=one
        Gam_e_max_max=zero
        do I_tobs=2,Num_R
            do I_jump=1,Num_jump
                dB=(B3_branch(I_jump,I_tobs)+B3_branch(I_jump,I_tobs-1))/two
                if (dB > zero .and. Gam_m_branch(I_jump,I_tobs) > one) then
                    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
                    Gam_e_max_max=max(Gam_e_max_max,Gam_e_max)
                end if
            end do
        end do
        if (Gam_e_max_max <= Gam_e_min_global) &
            error stop "electron_secondary_reverse_branch_reaccelerated: empty electron grid."
        do I_gam_e=1,Num_gam_e
            if (Num_gam_e == 1) then
                gam_e(I_gam_e)=Gam_e_min_global
            else
                gam_e(I_gam_e)=Gam_e_min_global*ten**(dlog10(Gam_e_max_max/Gam_e_min_global)*(I_gam_e-1)/(Num_gam_e-1))
            end if
        end do
    end subroutine build_secondary_reaccel_gamma_grid

    subroutine transfer_reaccelerated_parent_electrons(i_shell)
    implicit none
    integer, intent(in) :: i_shell

        fresh_mass_branch=zero
        do I_jump=1,Num_jump
            mass_delta=M3_branch(I_jump,i_shell)-M3_branch(I_jump,i_shell-1)
            if (mass_delta < zero) error stop "electron_secondary_reverse_branch_reaccelerated: branch mass decreased."
            dissipative_shell=Gam_m_branch(I_jump,i_shell) > one .and. Gamma43_branch(I_jump,i_shell) > one
            fresh_mass=zero
            if (dissipative_shell) fresh_mass=mass_delta
            if (Parent_branch(I_jump) > 0 .and. dissipative_shell .and. mass_delta > zero) then
                parent=Parent_branch(I_jump)
                parent_mass=branch_mass_available(parent)
                if (parent_mass > zero) then
                    transfer_mass=min(mass_delta,parent_mass)
                    transfer_fraction=transfer_mass/parent_mass
                    dN_seed=dN_work(parent,:)*transfer_fraction
                    dN_work(parent,:)=dN_work(parent,:)-dN_seed
                    boost_factor=Gamma43_branch(I_jump,i_shell)
                    if (Comp_branch(I_jump,i_shell) > one) boost_factor=boost_factor*Comp_branch(I_jump,i_shell)**(one/3d0)
                    call boost_log_distribution(boost_factor,dN_seed,dN_boost)
                    call dsa_reaccelerate_distribution(p_r,dN_boost,dN_reaccel)
                    seed_energy=distribution_energy_from_log(dN_seed)
                    out_energy=distribution_energy_from_log(dN_reaccel)
                    Branch_seed_energy(I_jump,i_shell)=Branch_seed_energy(I_jump,i_shell)+seed_energy
                    Branch_reaccel_energy(I_jump,i_shell)=Branch_reaccel_energy(I_jump,i_shell)+out_energy
                    dN_work(I_jump,:)=dN_work(I_jump,:)+dN_reaccel
                    branch_mass_available(parent)=branch_mass_available(parent)-transfer_mass
                    branch_mass_available(I_jump)=branch_mass_available(I_jump)+transfer_mass
                    fresh_mass=mass_delta-transfer_mass
                end if
            end if
            fresh_mass_branch(I_jump)=fresh_mass
            branch_mass_available(I_jump)=branch_mass_available(I_jump)+fresh_mass
            dN_x(I_jump,:)=dN_work(I_jump,:)
        end do
    end subroutine transfer_reaccelerated_parent_electrons

    subroutine advance_reaccelerated_branch_shell(i_shell,jump_index)
    implicit none
    integer, intent(in) :: i_shell,jump_index

        if (M3_branch(jump_index,i_shell) <= zero .and. M3_branch(jump_index,i_shell-1) <= zero) then
            dN_gamma_branch(jump_index,:,i_shell)=dN_x(jump_index,:)/gam_e/dlog(ten)
            return
        end if
        call prepare_branch_shell(i_shell,jump_index)
        call compute_branch_injection(i_shell,jump_index)
        do L=1,L1
            R_loc=R_loc+dDR
            if (injection_rate > zero) then
                call reverse_build_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m,Gam_e_max,injection_rate,p_r,dF1)
            else
                dF1=zero
            end if
            temp3=((dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two+adiabatic_rate)/dlog(ten)
            call electron_fullhide_flux_split_step(Num_gam_e,dDR,d_x,temp3,dF1,dN_x(jump_index,:),x,.true.)
            dN_x(jump_index,:)=x
        end do
        dN_gamma_branch(jump_index,:,i_shell)=dN_x(jump_index,:)/gam_e/dlog(ten)
    end subroutine advance_reaccelerated_branch_shell

    subroutine prepare_branch_shell(i_shell,jump_index)
    implicit none
    integer, intent(in) :: i_shell,jump_index

        R_loc=R(i_shell-1)
        R_Gamma_loc=(R_Gamma(i_shell)+R_Gamma(i_shell-1))/two
        beta2=dsqrt(one-one/R_Gamma_loc**2)
        dB=(B3_branch(jump_index,i_shell)+B3_branch(jump_index,i_shell-1))/two
        if (dB <= zero) error stop "electron_secondary_reverse_branch_reaccelerated: active branch requires B3 > 0."
        if (Gam_m_branch(jump_index,i_shell-1) > one) then
            Gam_e_m=(Gam_m_branch(jump_index,i_shell)+Gam_m_branch(jump_index,i_shell-1))/two
        else
            Gam_e_m=Gam_m_branch(jump_index,i_shell)
        end if
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        f_r=secondary_adv_coeff/beta2/R_Gamma_loc*dB**2/pi
        dDD=R(i_shell)-R(i_shell-1)
        dDR=0.7d0/(f_r*Gam_e_max+one/R(i_shell-1))
        L1=max(100,min(1000,int(dDD/dDR)))
        dDR=dDD/L1
        dEl=f_r*gam_e
        if (V3_branch(jump_index,i_shell) <= zero .or. V3_branch(jump_index,i_shell-1) <= zero) then
            adiabatic_rate=zero
        else
            adiabatic_rate=dlog(V3_branch(jump_index,i_shell)/V3_branch(jump_index,i_shell-1))/(3d0*dDD)
        end if
    end subroutine prepare_branch_shell

    subroutine compute_branch_injection(i_shell,jump_index)
    implicit none
    integer, intent(in) :: i_shell,jump_index

        mass_lo=M3_branch(jump_index,i_shell-1)
        mass_hi=M3_branch(jump_index,i_shell)
        if (mass_hi < mass_lo) error stop "electron_secondary_reverse_branch_reaccelerated: swept mass must not decrease."
        injection_rate=f_e_r*fresh_mass_branch(jump_index)/(Para_m_p*dDD)
    end subroutine compute_branch_injection

    subroutine boost_log_distribution(boost,dN_in,dN_out)
    implicit none
    real(8), intent(in) :: boost,dN_in(Num_gam_e)
    real(8), intent(out) :: dN_out(Num_gam_e)
    integer :: i_src,i_hi
    real(8) :: x_src,pos,frac

        if (boost <= one) error stop "electron_secondary_reverse_branch_reaccelerated: boost must exceed unity."
        dN_out=zero
        do I_gam_e=1,Num_gam_e
            x_src=dlog10(gam_e(I_gam_e)/boost)
            if (x_src < dlog10(gam_e(1)) .or. x_src > dlog10(gam_e(Num_gam_e))) cycle
            pos=(x_src-dlog10(gam_e(1)))/d_x+one
            i_src=int(pos)
            if (i_src < 1) cycle
            if (i_src >= Num_gam_e) then
                dN_out(I_gam_e)=dN_in(Num_gam_e)
            else
                i_hi=i_src+1
                frac=pos-dble(i_src)
                dN_out(I_gam_e)=(one-frac)*dN_in(i_src)+frac*dN_in(i_hi)
            end if
        end do
    end subroutine boost_log_distribution

    subroutine dsa_reaccelerate_distribution(p,dN_seed_log,dN_out_log)
    implicit none
    real(8), intent(in) :: p,dN_seed_log(Num_gam_e)
    real(8), intent(out) :: dN_out_log(Num_gam_e)
    integer :: i
    real(8) :: integral,dN_dgamma

        if (p <= two) error stop "electron_secondary_reverse_branch_reaccelerated: DSA requires p > 2."
        integral=zero
        dN_out_log=zero
        do i=1,Num_gam_e
            dN_dgamma=dN_seed_log(i)/(gam_e(i)*dlog(ten))
            integral=integral+dN_dgamma*gam_e(i)**(p-one)*(gam_e(i)*dlog(ten)*d_x)
            dN_out_log(i)=(p-one)*gam_e(i)**(-p)*integral*gam_e(i)*dlog(ten)
        end do
    end subroutine dsa_reaccelerate_distribution

    real(8) function distribution_energy_from_log(dN_log)
    implicit none
    real(8), intent(in) :: dN_log(Num_gam_e)

        distribution_energy_from_log=sum((gam_e-one)*Para_m_e*Para_c**2*dN_log)*d_x
    end function distribution_energy_from_log
end subroutine electron_secondary_reverse_branch_reaccelerated

subroutine reverse_build_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m,Gam_e_max,injection_rate,p,dF1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e,I_q
    real(8), intent(in) :: x_edge(Num_gam_e+1),Gam_e_m,Gam_e_max,injection_rate,p
    real(8), intent(out) :: dF1(Num_gam_e)
    real(8), parameter :: xi(3)=(/-dsqrt(3d0/5d0),zero,dsqrt(3d0/5d0)/)
    real(8), parameter :: wi(3)=(/5d0/9d0,8d0/9d0,5d0/9d0/)
    real(8) :: cell_lo,cell_hi,dx_cell,half_dx,x_mid,x_eval,gam,cutoff_factor,cell_sum,shape_norm

    dF1=zero
    shape_norm=zero
    if (Gam_e_max <= zero .or. injection_rate <= zero) return

    do I_gam_e=1,Num_gam_e
        cell_lo=x_edge(I_gam_e)
        cell_hi=x_edge(I_gam_e+1)
        dx_cell=cell_hi-cell_lo
        if (dx_cell <= zero) cycle
        half_dx=0.5d0*dx_cell
        x_mid=0.5d0*(cell_lo+cell_hi)
        cell_sum=zero
        do I_q=1,3
            x_eval=x_mid+half_dx*xi(I_q)
            gam=ten**x_eval
            if (gam > Gam_e_m) then
                cutoff_factor=electron_exp_cutoff_factor(gam,Gam_e_max)
                cell_sum=cell_sum+wi(I_q)*gam*dlog(ten)*(gam-one)**(-p)*cutoff_factor
            end if
        end do
        cell_sum=half_dx*cell_sum
        dF1(I_gam_e)=cell_sum/dx_cell
        shape_norm=shape_norm+cell_sum
    end do
    dF1=injection_rate*dF1/shape_norm
end subroutine reverse_build_source_term_exp_cutoff_edges

end module electron_reverse_kernel
