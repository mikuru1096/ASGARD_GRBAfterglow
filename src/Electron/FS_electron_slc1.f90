subroutine fs_electron_slc1(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                            gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    !$ use omp_lib
    use constants
    use electron_common
    use get_Y
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R), &
                            Seed_syn(Num_nu,Num_R),V_m(Num_R),V_c(Num_R),V_a(Num_R)

    real(8), allocatable :: dEl(:),dEL_mean(:),dEL_mean_base(:),dN_x(:),dN_step(:),dF1(:)
    real(8), allocatable :: gam_e_rad(:),dN_gam_e_rad(:)
    logical :: is_uniform_density,anchor_gamma_m,anchor_gamma_c
    integer :: Num_gam_rad

    allocate(dEl(Num_gam_e),dEL_mean(Num_gam_e-1),dEL_mean_base(Num_gam_e-1),dN_x(Num_gam_e), &
             dN_step(Num_gam_e),dF1(Num_gam_e),gam_e_rad(Num_gam_e),dN_gam_e_rad(Num_gam_e))

    Eta_0=Boundary(1)
    R_ini=Boundary(4)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    f_e=Boundary(16)
    R_tr=Boundary(21)
    f_jump=Boundary(22)
    f_wide=Boundary(23)
    R0=Boundary(n)

    P_syn=zero
    Seed_syn=zero
    V_m=zero
    V_c=zero
    V_a=zero

    call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)
    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-one)
    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    call electron_build_gamma_grid(Num_gam_e,Gam_e_max_max,Gam_e)
    call electron_initial_powerlaw_exp_cutoff(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,dN_gam_e(:,1))

    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
    d_x=dlog10(gam_e(2)/gam_e(1))
    is_uniform_density=(A_star <= zero .and. f_jump == one)

    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        call electron_external_density(A_star,dNe_ISM,R_loc,R0,R_tr,f_jump,f_wide,1,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        dNe_shell=dNe

        beta_Gam=sqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=max(20,min(200,Int(dDD/dDR)))
        dDR=dDD/L1

        call electron_prepare_radiation_spectrum(Num_gam_e,gam_e,dN_gam_e(:,I_tobs-1), &
                                                 Num_gam_rad,gam_e_rad,dN_gam_e_rad)
        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        call get_nu_a(R_loc,DB,Num_gam_rad,gam_e_rad(1:Num_gam_rad),dN_gam_e_rad(1:Num_gam_rad),temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))

        call get_syn_selected(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                              gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_syn(:,I_tobs),Seed_syn(:,I_tobs))
        call get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                 R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                 P_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl)
        call electron_loss_mean(Num_gam_e,dEl,dEL_mean_base)

        do L=1,L1
            R_left=R_loc
            dNe_left=dNe
            R_right=R_left+dDR
            if (.not. is_uniform_density) then
                call electron_external_density(A_star,dNe_ISM,R_right,R0,R_tr,f_jump,f_wide,1,dNe_right)
                R_mid=0.5d0*(R_left+R_right)
                call electron_external_density(A_star,dNe_ISM,R_mid,R0,R_tr,f_jump,f_wide,1,dNe_mid)
            else
                dNe_right=dNe_left
                R_mid=0.5d0*(R_left+R_right)
                dNe_mid=dNe_left
            end if
            DB_step=0.39d0*dsqrt(Epsilon_b*dNe_mid*(R_Gamma_loc*(R_Gamma_loc-one)))
            Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
            temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
            call electron_gamma_m_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
            Gam_e_m_p_step=(one-p)/(Gam_e_max_step**(one-p)-Gam_e_m_step**(one-p))
            call electron_injection_prefactor(R_mid,dDR,dNe_mid,f_e,Gam_e_m_p_step,Q)
            call electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
            if (dNe_shell > zero) then
                dEL_mean=dEL_mean_base*(dNe_mid/dNe_shell)
            else
                dEL_mean=dEL_mean_base
            end if
            call electron_semi_lagrangian_step(Num_gam_e,dDR,d_x,dEL_mean+one/R_mid/dlog(ten),dF1,dN_x,dN_step)
            dN_x=dN_step
            R_loc=R_right
            dNe=dNe_right
        end do
        dN_gam_e(:,I_tobs)=dN_x/gam_e/dlog(ten)
    end do

    deallocate(dEl,dEL_mean,dEL_mean_base,dN_x,dN_step,dF1,gam_e_rad,dN_gam_e_rad)
end subroutine fs_electron_slc1

subroutine fs_electron_slc1_mmg2(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                                 gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a,work_x_edge_log10,work_dN_x)
    !$ use omp_lib
    use constants
    use electron_common
    use get_Y
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R), &
                            Seed_syn(Num_nu,Num_R),V_m(Num_R),V_c(Num_R),V_a(Num_R), &
                            work_x_edge_log10(Num_gam_e+1,Num_R),work_dN_x(Num_gam_e,Num_R)

    real(8), allocatable :: dEl(:),dEL_mean(:),dEL_mean_base(:),dN_x_work(:),dN_step(:),dN_public(:),dF1(:)
    real(8), allocatable :: gam_e_work(:),dN_gam_e_work(:)
    real(8), allocatable :: x_ref_edge(:),x_work_edge(:),x_target_edge(:),monitor(:)
    real(8) :: mmg_weights(5),mmg_window_dex(3),mmg_smooth,x_cut_low,x_cut_high
    real(8) :: x_lower_target,x_upper_target,x_lower_margin,x_upper_margin,x_char_low,x_char_m,x_char_c
    logical :: is_uniform_density,anchor_gamma_m,anchor_gamma_c

    allocate(dEl(Num_gam_e),dEL_mean(Num_gam_e-1),dEL_mean_base(Num_gam_e-1), &
             dN_x_work(Num_gam_e),dN_step(Num_gam_e),dN_public(Num_gam_e),dF1(Num_gam_e), &
             gam_e_work(Num_gam_e),dN_gam_e_work(Num_gam_e), &
             x_ref_edge(Num_gam_e+1),x_work_edge(Num_gam_e+1),x_target_edge(Num_gam_e+1),monitor(Num_gam_e))

    Eta_0=Boundary(1)
    R_ini=Boundary(4)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    f_e=Boundary(16)
    R_tr=Boundary(21)
    f_jump=Boundary(22)
    f_wide=Boundary(23)
    R0=Boundary(n)

    mmg_weights=(/6d0,4d0,8d0,2d0,1d0/)
    mmg_window_dex=(/0.18d0,0.18d0,0.10d0/)
    mmg_smooth=0.5d0

    P_syn=zero
    Seed_syn=zero
    V_m=zero
    V_c=zero
    V_a=zero
    work_x_edge_log10=zero
    work_dN_x=zero

    call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,Para_N_e_ini)
    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-one)
    call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    call electron_build_gamma_grid(Num_gam_e,Gam_e_max_max,Gam_e)
    call electron_log_cell_edges(Num_gam_e,gam_e,x_ref_edge)
    x_work_edge=x_ref_edge
    call electron_gamma_from_edges(Num_gam_e,x_work_edge,gam_e_work)
    call electron_initial_powerlaw_exp_cutoff_edges(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,x_work_edge,dN_x_work)
    dN_gam_e_work=dN_x_work/gam_e_work/dlog(ten)
    call electron_conservative_remap_nonuniform(Num_gam_e,x_work_edge,x_ref_edge,dN_x_work,dN_public)
    dN_gam_e(:,1)=dN_public/gam_e/dlog(ten)
    work_x_edge_log10(:,1)=x_work_edge
    work_dN_x(:,1)=dN_x_work

    is_uniform_density=(A_star <= zero .and. f_jump == one)

    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        call electron_external_density(A_star,dNe_ISM,R_loc,R0,R_tr,f_jump,f_wide,1,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
        call electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
        Gam_e_m_p=(one-p)/(Gam_e_max**(one-p)-Gam_e_m**(one-p))
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        dNe_shell=dNe

        call electron_find_low_energy_front(Num_gam_e,x_work_edge,dN_x_work,x_cut_low)
        call electron_find_high_energy_front(Num_gam_e,x_work_edge,dN_x_work,x_cut_high)
        x_char_m=dlog10(max(Gam_e_m,1d0))
        x_char_c=dlog10(max(Gam_e_c,1d0))
        x_char_low=min(x_char_m,x_char_c)
        x_lower_margin=max(2.5d0*mmg_window_dex(1),0.20d0)
        x_lower_target=min(x_cut_low,x_char_low)-x_lower_margin
        x_lower_target=max(x_ref_edge(1),min(x_lower_target,x_work_edge(Num_gam_e+1)-1d-4))
        if (x_lower_target < x_work_edge(1)-1d-6 .or. x_lower_target > x_work_edge(1)+1d-6) then
            call electron_rescale_lower_boundary(Num_gam_e,x_work_edge,x_lower_target,x_target_edge)
            call electron_conservative_remap_nonuniform(Num_gam_e,x_work_edge,x_target_edge,dN_x_work,dN_step)
            x_work_edge=x_target_edge
            dN_x_work=dN_step
            call electron_gamma_from_edges(Num_gam_e,x_work_edge,gam_e_work)
            dN_gam_e_work=dN_x_work/gam_e_work/dlog(ten)
        end if

        x_upper_margin=max(3d0*mmg_window_dex(3),0.25d0)
        x_upper_target=max(dlog10(max(Gam_e_max,1d0))+x_upper_margin,x_cut_high+x_upper_margin)
        x_upper_target=min(x_ref_edge(Num_gam_e+1),max(x_upper_target,x_ref_edge(1)+1d-4))
        if (x_upper_target < x_work_edge(Num_gam_e+1)-1d-6 .or. x_upper_target > x_work_edge(Num_gam_e+1)+1d-6) then
            call electron_rescale_upper_boundary(Num_gam_e,x_work_edge,x_upper_target,x_target_edge)
            call electron_conservative_remap_nonuniform(Num_gam_e,x_work_edge,x_target_edge,dN_x_work,dN_step)
            x_work_edge=x_target_edge
            dN_x_work=dN_step
            call electron_gamma_from_edges(Num_gam_e,x_work_edge,gam_e_work)
            dN_gam_e_work=dN_x_work/gam_e_work/dlog(ten)
        end if

        call electron_build_moving_monitor(Num_gam_e,x_work_edge,dN_x_work,Gam_e_m,Gam_e_c,Gam_e_max, &
                                           mmg_weights,mmg_window_dex,monitor,x_cut_low,x_cut_high)
        anchor_gamma_m=(x_char_m > x_cut_low+max(6d0*mmg_window_dex(1),0.8d0)) .and. &
     &                 (x_char_m < x_cut_high-max(2d0*mmg_window_dex(1),0.2d0))
        anchor_gamma_c=(x_char_c > x_cut_low+max(4d0*mmg_window_dex(2),0.4d0)) .and. &
     &                 (x_char_c < x_cut_high-max(2d0*mmg_window_dex(2),0.2d0))
        call electron_equidistribute_edges(Num_gam_e,x_work_edge,monitor,mmg_smooth,x_target_edge)
        call electron_anchor_low_energy_edges(Num_gam_e,x_work_edge,x_target_edge,x_cut_low,mmg_window_dex(1))
        if (anchor_gamma_m) call electron_anchor_characteristic_edges(Num_gam_e,x_target_edge,x_char_m,mmg_window_dex(1))
        if (anchor_gamma_c) call electron_anchor_characteristic_edges(Num_gam_e,x_target_edge,x_char_c,mmg_window_dex(2))
        call electron_anchor_high_energy_edges(Num_gam_e,x_work_edge,x_target_edge,x_cut_high,mmg_window_dex(3))
        if (anchor_gamma_m) call electron_anchor_characteristic_edges(Num_gam_e,x_target_edge,x_char_m,mmg_window_dex(1))
        if (anchor_gamma_c) call electron_anchor_characteristic_edges(Num_gam_e,x_target_edge,x_char_c,mmg_window_dex(2))
        call electron_conservative_remap_nonuniform(Num_gam_e,x_work_edge,x_target_edge,dN_x_work,dN_step)
        x_work_edge=x_target_edge
        dN_x_work=dN_step
        call electron_gamma_from_edges(Num_gam_e,x_work_edge,gam_e_work)
        dN_gam_e_work=dN_x_work/gam_e_work/dlog(ten)

        beta_Gam=sqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=max(10,min(120,Int(dDD/dDR)))
        dDR=dDD/L1

        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        call get_nu_a_nonuniform(R_loc,DB,Num_gam_e,x_work_edge,dN_x_work,temp)
        V_a(I_tobs-1)=temp/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))

        call electron_conservative_remap_nonuniform(Num_gam_e,x_work_edge,x_ref_edge,dN_x_work,dN_public)
        dN_gam_e(:,I_tobs)=dN_public/gam_e/dlog(ten)
        call get_syn_selected(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                              gam_e,dN_gam_e(:,I_tobs),V_seed,P_syn(:,I_tobs),Seed_syn(:,I_tobs))
        work_x_edge_log10(:,I_tobs)=x_work_edge
        work_dN_x(:,I_tobs)=dN_x_work
        call get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                 R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e_work,V_seed, &
                                 P_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl)
        call electron_loss_mean(Num_gam_e,dEl,dEL_mean_base)

        do L=1,L1
            R_left=R_loc
            dNe_left=dNe
            R_right=R_left+dDR
            if (.not. is_uniform_density) then
                call electron_external_density(A_star,dNe_ISM,R_right,R0,R_tr,f_jump,f_wide,1,dNe_right)
                R_mid=0.5d0*(R_left+R_right)
                call electron_external_density(A_star,dNe_ISM,R_mid,R0,R_tr,f_jump,f_wide,1,dNe_mid)
            else
                dNe_right=dNe_left
                R_mid=0.5d0*(R_left+R_right)
                dNe_mid=dNe_left
            end if
            DB_step=0.39d0*dsqrt(Epsilon_b*dNe_mid*(R_Gamma_loc*(R_Gamma_loc-one)))
            Gam_e_max_step=3d0*Para_m_energy/dsqrt(8d0*DB_step*Para_e**3)
            temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma_loc-one)
            call electron_gamma_m_exact(p,temp_gam,Gam_e_max_step,Gam_e_m_step)
            Gam_e_m_p_step=(one-p)/(Gam_e_max_step**(one-p)-Gam_e_m_step**(one-p))
            call electron_injection_prefactor(R_mid,dDR,dNe_mid,f_e,Gam_e_m_p_step,Q)
            call electron_build_source_term_exp_cutoff_edges(Num_gam_e,x_work_edge,Gam_e_m_step,Gam_e_max_step,Q,p,dF1)
            if (dNe_shell > zero) then
                dEL_mean=dEL_mean_base*(dNe_mid/dNe_shell)
            else
                dEL_mean=dEL_mean_base
            end if
            call electron_semi_lagrangian_step_nonuniform(Num_gam_e,dDR,x_work_edge, &
                                                          dEL_mean+one/R_mid/dlog(ten),dF1,dN_x_work,dN_step)
            dN_x_work=dN_step
            R_loc=R_right
            dNe=dNe_right
        end do
    end do

    deallocate(dEl,dEL_mean,dEL_mean_base,dN_x_work,dN_step,dN_public,dF1, &
               gam_e_work,dN_gam_e_work,x_ref_edge,x_work_edge,x_target_edge,monitor)
end subroutine fs_electron_slc1_mmg2
