!f2py: skip
module electron_common
    use constants
    use adaptive_resampling_mod, only: adaptive_resampling_log
    use electron_injection_profiles, only: electron_profile_log_cell_edges, &
                                     electron_initial_powerlaw_exp_cutoff, electron_initial_powerlaw_exp_cutoff_edges, &
                                     electron_add_thermal_population
    implicit none
    integer, parameter :: radiation_resample_threshold = 180
    integer, parameter :: radiation_resample_target = 160
    integer, parameter :: radiation_resample_smoothness = 4
    integer, parameter :: electron_initial_grid_gamma = 0
    integer, parameter :: electron_initial_grid_log_edges = 1
contains

! 解包公共 Boundary 数组字段。
subroutine electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                    E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: Boundary(n)
    real(8), intent(out) :: Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star
    real(8), intent(out) :: E_iso,T_log10_duration,f_e,R_tr,f_jump,f_wide,R0

    Eta_0=Boundary(1)
    R_ini=Boundary(4)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    E_iso=Boundary(14)
    T_log10_duration=Boundary(15)
    f_e=Boundary(16)
    R_tr=Boundary(21)
    f_jump=Boundary(22)
    f_wide=Boundary(23)
    ! Modern Boundary arrays reserve slot 27 for the external-density radius scale R0.
    if (n >= 27) then
        R0=Boundary(27)
    else
        R0=Boundary(n)
    end if
end subroutine electron_unpack_boundary

! 初始化电子能谱，并按需生成 gamma 或 log-gamma 边界网格。
subroutine electron_initialize_spectrum(Num_gam_e,Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                        grid_mode,gam_e,dN_init,x_edge,thermal_electrons,f_e,four_v)
    implicit none
    integer, intent(in) :: Num_gam_e,grid_mode
    integer :: I_gam_e
    real(8), intent(in) :: Gam_e_max_max,Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max
    real(8), intent(out) :: gam_e(Num_gam_e),dN_init(Num_gam_e)
    real(8), intent(out), optional :: x_edge(Num_gam_e+1)
    integer, intent(in), optional :: thermal_electrons
    real(8), intent(in), optional :: f_e,four_v

    do I_gam_e=1,Num_gam_e
        gam_e(I_gam_e)=3d0*ten**(dlog10(Gam_e_max_max)*(I_gam_e-1)/(Num_gam_e-1))
    end do
    select case (grid_mode)
    case (electron_initial_grid_gamma)
        call electron_initial_powerlaw_exp_cutoff(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,dN_init)
    case (electron_initial_grid_log_edges)
        call electron_profile_log_cell_edges(Num_gam_e,gam_e,x_edge)
        call electron_initial_powerlaw_exp_cutoff_edges(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,x_edge,dN_init)
    end select
    if (present(thermal_electrons)) then
        if (thermal_electrons /= 0) then
            call electron_add_thermal_population(Num_gam_e,gam_e,four_v,Para_N_e_ini*(one-f_e),dN_init)
        end if
    end if
end subroutine electron_initialize_spectrum

! 按完整 p 分支计算注入最小洛伦兹因子。
subroutine electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    implicit none
    real(8), intent(in) :: p,temp_gam,Gam_e_max
    real(8), intent(out) :: Gam_e_m
    real(8) :: eps,temp

    Gam_e_m=(p-two)/(p-one)*temp_gam+one
    if (p>two) then
        Gam_e_m=(p-two)/(p-one)*temp_gam+one
    else if (p<two .and. p>one) then
        Gam_e_m=((two-p)/(p-one)*temp_gam*Gam_e_max**(p-two))**(one/(p-one))+one
    else if (p==two) then
        eps=1d-5
        Gam_e_m=one
        temp=temp_gam/log(Gam_e_max/Gam_e_m)
        do while (abs(one-Gam_e_m/temp)>eps)
            temp=temp_gam/log(Gam_e_max/Gam_e_m)
            if (Gam_e_m>temp) then
                Gam_e_m=0.5d0*(Gam_e_m+temp)
            else
                Gam_e_m=0.5d0*(Gam_e_m+Gam_e_max)
            end if
        end do
    end if
end subroutine electron_gamma_m_exact

! 由面平均冷却率求冷却洛伦兹因子。
subroutine electron_gamma_c_from_loss_mean(Num_gam_e,gam_e,dEL_mean,R_loc,Gam_e_c)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e, I_cross
    real(8), intent(in) :: gam_e(Num_gam_e),dEL_mean(Num_gam_e-1),R_loc
    real(8), intent(out) :: Gam_e_c
    real(8) :: coeff_target,gam_mid(Num_gam_e-1)
    real(8) :: x0,x1,y0,y1,ytarget,xroot

    coeff_target=one/(R_loc*dlog(ten))
    do I_gam_e=1,Num_gam_e-1
        gam_mid(I_gam_e)=dsqrt(gam_e(I_gam_e)*gam_e(I_gam_e+1))
    end do
    if (dEL_mean(Num_gam_e-1) <= coeff_target) then
        Gam_e_c=gam_mid(Num_gam_e-1)
        return
    end if
    if (dEL_mean(1) >= coeff_target) then
        Gam_e_c=gam_mid(1)
        return
    end if
    I_cross=0
    do I_gam_e=Num_gam_e-2,1,-1
        if (dEL_mean(I_gam_e) <= coeff_target .and. dEL_mean(I_gam_e+1) > coeff_target) then
            I_cross=I_gam_e
            exit
        end if
    end do
    if (I_cross == 0) then
        Gam_e_c=gam_mid(Num_gam_e-1)
        return
    end if
    x0=dlog(gam_mid(I_cross))
    x1=dlog(gam_mid(I_cross+1))
    y0=dlog(max(dEL_mean(I_cross),tiny(one)))
    y1=dlog(max(dEL_mean(I_cross+1),tiny(one)))
    ytarget=dlog(coeff_target)
    if (y1 == y0) then
        xroot=0.5d0*(x0+x1)
    else
        xroot=x0+(ytarget-y0)*(x1-x0)/(y1-y0)
    end if
    Gam_e_c=dexp(xroot)
end subroutine electron_gamma_c_from_loss_mean

! 计算当前径向壳层注入源项的归一化系数。
subroutine electron_injection_prefactor(R_loc,dDR,dNe,f_e,Gam_e_m_p,Q)
    implicit none
    real(8), intent(in) :: R_loc,dDR,dNe,f_e,Gam_e_m_p
    real(8), intent(out) :: Q

    Q=4d0/3d0*pi*(3d0*R_loc**2+dDR*(3d0*R_loc+dDR))*dNe*f_e*Gam_e_m_p
end subroutine electron_injection_prefactor

! 在电子网格上定位注入区间覆盖的单元范围。
subroutine electron_source_bounds(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,src_lo,src_hi)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer, intent(out) :: src_lo,src_hi
    integer :: I_gam_e
    real(8), intent(in) :: gam_e(Num_gam_e),Gam_e_m,Gam_e_max
    real(8) :: x_edge(Num_gam_e+1),x_lo,x_hi

    call electron_profile_log_cell_edges(Num_gam_e,gam_e,x_edge)
    x_lo=dlog10(Gam_e_m)
    x_hi=dlog10(Gam_e_max)
    src_lo=Num_gam_e+1
    do I_gam_e=1,Num_gam_e
        if (x_edge(I_gam_e+1) > x_lo) then
            src_lo=I_gam_e
            exit
        end if
    end do
    src_hi=0
    do I_gam_e=Num_gam_e,1,-1
        if (x_edge(I_gam_e) < x_hi) then
            src_hi=I_gam_e
            exit
        end if
    end do
    if (src_hi < src_lo) then
        src_lo=1
        src_hi=0
    end if
end subroutine electron_source_bounds

! 为辐射积分压缩活跃电子谱网格。
subroutine electron_prepare_radiation_spectrum(Num_gam_e,gam_e,dN_gam_e,Num_gam_rad,gam_e_rad,dN_gam_e_rad)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer, intent(out) :: Num_gam_rad
    integer :: I_gam_e,left_pos,right_pos,first_pos,last_pos,n_active,m_target,n_resampled,info
    integer :: idx(Num_gam_e),out_idx,src_idx,last_added_idx
    real(8), intent(in) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e)
    real(8), intent(out) :: gam_e_rad(Num_gam_e),dN_gam_e_rad(Num_gam_e)

    first_pos = 0
    last_pos = 0
    left_pos = 1
    right_pos = Num_gam_e
    do while (left_pos <= right_pos .and. (first_pos == 0 .or. last_pos == 0))
        if (first_pos == 0) then
            if (dN_gam_e(left_pos) > zero) then
                first_pos = left_pos
            else
                left_pos = left_pos + 1
            end if
        end if
        if (last_pos == 0 .and. right_pos >= left_pos) then
            if (dN_gam_e(right_pos) > zero) then
                last_pos = right_pos
            else
                right_pos = right_pos - 1
            end if
        end if
    end do
    gam_e_rad = gam_e
    dN_gam_e_rad = dN_gam_e
    Num_gam_rad = Num_gam_e
    if (first_pos == 0) return
    n_active = last_pos - first_pos + 1
    if (n_active <= radiation_resample_threshold) return
    m_target = min(n_active, radiation_resample_target)
    if (m_target >= n_active) return
    call adaptive_resampling_log(gam_e(first_pos:last_pos), dN_gam_e(first_pos:last_pos), n_active, &
                                 m_target, radiation_resample_smoothness, idx(1:m_target), n_resampled, info)
    if (info /= 0 .or. n_resampled <= 0) return
    gam_e_rad = zero
    dN_gam_e_rad = zero
    out_idx = 0
    if (first_pos > 1) then
        out_idx = out_idx + 1
        gam_e_rad(out_idx) = gam_e(first_pos-1)
        dN_gam_e_rad(out_idx) = dN_gam_e(first_pos-1)
    end if
    out_idx = out_idx + 1
    gam_e_rad(out_idx) = gam_e(first_pos)
    dN_gam_e_rad(out_idx) = dN_gam_e(first_pos)
    last_added_idx = first_pos
    do I_gam_e = 1, n_resampled
        src_idx = first_pos + idx(I_gam_e) - 1
        if (src_idx > first_pos .and. src_idx < last_pos) then
            if (src_idx /= last_added_idx) then
                out_idx = out_idx + 1
                gam_e_rad(out_idx) = gam_e(src_idx)
                dN_gam_e_rad(out_idx) = dN_gam_e(src_idx)
                last_added_idx = src_idx
            end if
        end if
    end do
    if (last_pos > first_pos) then
        out_idx = out_idx + 1
        gam_e_rad(out_idx) = gam_e(last_pos)
        dN_gam_e_rad(out_idx) = dN_gam_e(last_pos)
        last_added_idx = last_pos
    end if
    if (last_pos < Num_gam_e) then
        out_idx = out_idx + 1
        gam_e_rad(out_idx) = gam_e(last_pos+1)
        dN_gam_e_rad(out_idx) = dN_gam_e(last_pos+1)
    end if
    Num_gam_rad = out_idx
end subroutine electron_prepare_radiation_spectrum

! 计算两个数组之间的最大相对误差。
subroutine electron_max_relative_error(Num_gam_e,x_ref,x_trial,error_max)
    implicit none
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: x_ref(Num_gam_e),x_trial(Num_gam_e)
    real(8), intent(out) :: error_max
    integer :: I_gam_e
    real(8) :: denom

    error_max=zero
    do I_gam_e=1,Num_gam_e
        denom=max(abs(x_ref(I_gam_e)),1d-99)
        error_max=max(error_max,abs(x_trial(I_gam_e)-x_ref(I_gam_e))/denom)
    end do
end subroutine electron_max_relative_error

! 按风介质或均匀介质条件初始化外部电子密度和累计电子数。
subroutine electron_initial_density(A_star,dNe_ISM,R_ini,R_start,R0,dNe,Para_N_e_ini)
    implicit none
    real(8), intent(in) :: A_star,dNe_ISM,R_ini,R_start,R0
    real(8), intent(out) :: dNe,Para_N_e_ini
    real(8) :: dNe_wind

    if (A_star > zero) then
        dNe_wind=A_star*3.0d35/R_start**2
        Para_N_e_ini=4d0*pi*R_ini*A_star*3.0d35
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
        Para_N_e_ini=4d0/3d0*pi*R_ini**3*dNe_ISM
    end if
    if (R_start < R0) then
        if (A_star > zero) then
            dNe=A_star*3.0d35/R0**2*4d0
        else
            dNe=dNe_ISM
        end if
        if (A_star > zero) then
            Para_N_e_ini=4d0*pi*R_ini*A_star*3.0d35
        else
            Para_N_e_ini=4d0/3d0*pi*R_ini**3*dNe_ISM
        end if
    end if
end subroutine electron_initial_density

end module electron_common
