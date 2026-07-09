!f2py: skip
module electron_common
    use constants
    use dynamics_density_profile, only: set_density_profile
    use electron_injection_profiles, only: log_edges, &
                                     init_powerlaw, init_edges, &
                                     thermal_pop
    implicit none
    public :: imodeg, imodelog, tail_factor
    integer, parameter :: imodeg = 0, imodelog = 1
    real(8), parameter :: tail_factor = 30d0
contains

! 解包公共 Boundary 数组字段。
! Unpack public Boundary-array fields.
subroutine electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                    E_iso,tdur_log,f_e,R_tr,f_jump,f_wide,R0)
    implicit none
    integer, intent(in) :: n
    integer, parameter :: r0_slot = 27
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(out) :: Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star
    real(8), intent(out) :: E_iso,tdur_log,f_e,R_tr,f_jump,f_wide,R0

    Eta_0=Boundary(1)
    R_ini=Boundary(4)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    E_iso=Boundary(14)
    tdur_log=Boundary(15)
    f_e=Boundary(16)
    R_tr=Boundary(21)
    f_jump=Boundary(22)
    f_wide=Boundary(23)
    if (n >= r0_slot) then
        R0=Boundary(r0_slot)
    else
        R0=Boundary(n)
    end if
    call set_density_profile(Boundary,n)
end subroutine electron_unpack_boundary

! 初始化电子能谱，并按需生成 gamma 或 log-gamma 边界网格。
! Initialize the electron spectrum and optionally build gamma or log-gamma cell edges.
subroutine electron_initialize_spectrum(ng,gemax0,ninit,p,gm,gc,gemax, &
                                        grid_mode,gam_e,dN_init,x_edge,thermal_electrons,f_e,four_v)
    implicit none
    integer, intent(in) :: ng,grid_mode
    integer :: ig
    real(8), intent(in) :: gemax0,ninit,p,gm,gc,gemax
    real(8), intent(out), dimension(ng) :: gam_e,dN_init
    real(8), intent(out), optional, dimension(ng+1) :: x_edge
    integer, intent(in), optional :: thermal_electrons
    real(8), intent(in), optional :: f_e,four_v

    do ig=1,ng
        gam_e(ig)=3d0*dexp(dlog(tail_factor*gemax0/3d0)* &
                                 (ig-1)/(ng-1))
    end do
    select case (grid_mode)
    case (imodeg)
        call init_powerlaw(ninit,p,gm,gc,gemax,ng,gam_e,dN_init)
    case (imodelog)
        call log_edges(ng,gam_e,x_edge)
        call init_edges(ninit,p,gm,gc,gemax,ng,x_edge,dN_init)
    end select
    if (present(thermal_electrons)) then
        if (thermal_electrons /= 0) then
            call thermal_pop(ng,gam_e,four_v,ninit*(1d0-f_e),dN_init)
        end if
    end if
end subroutine electron_initialize_spectrum

! 按完整 p 分支计算注入最小洛伦兹因子。
! Compute the injection minimum Lorentz factor with the full p-branch formula.
subroutine electron_gm_exact(p,temp_gam,gemax,gm)
    implicit none
    real(8), intent(in) :: p,temp_gam,gemax
    real(8), intent(out) :: gm
    real(8) :: eps,temp

    if (p>2d0) then
        gm=(p-2d0)/(p-1d0)*temp_gam+1d0
    else if (p<2d0 .and. p>1d0) then
        gm=((2d0-p)/(p-1d0)*temp_gam*gemax**(p-2d0))**(1d0/(p-1d0))+1d0
    else if (p==2d0) then
        eps=1d-5
        gm=1d0
        temp=temp_gam/log(gemax/gm)
        do while (abs(1d0-gm/temp)>eps)
            temp=temp_gam/log(gemax/gm)
            if (gm>temp) then
                gm=0.5d0*(gm+temp)
            else
                gm=0.5d0*(gm+gemax)
            end if
        end do
    end if
end subroutine electron_gm_exact

! 由面平均冷却率求冷却洛伦兹因子。
! Compute the cooling Lorentz factor from face-averaged cooling rates.
subroutine electron_gc_loss(ng,gam_e,dEL_mean,R_loc,gc)
    implicit none
    integer, intent(in) :: ng
    integer :: ig, I_cross
    real(8), intent(in), dimension(ng) :: gam_e
    real(8), intent(in), dimension(ng-1) :: dEL_mean
    real(8), intent(in) :: R_loc
    real(8), intent(out) :: gc
    real(8), dimension(ng-1) :: gam_mid
    real(8) :: coeff_target, x0, x1, y0, y1, ytarget, xroot

    coeff_target=1d0/(R_loc)
    do ig=1,ng-1
        gam_mid(ig)=dsqrt(gam_e(ig)*gam_e(ig+1))
    end do
    if (dEL_mean(ng-1) <= coeff_target) then
        gc=gam_mid(ng-1)
        return
    end if
    if (dEL_mean(1) >= coeff_target) then
        gc=gam_mid(1)
        return
    end if
    I_cross=0
    do ig=ng-2,1,-1
        if (dEL_mean(ig) <= coeff_target .and. dEL_mean(ig+1) > coeff_target) then
            I_cross=ig
            exit
        end if
    end do
    if (I_cross == 0) then
        gc=gam_mid(ng-1)
        return
    end if
    x0=dlog(gam_mid(I_cross))
    x1=dlog(gam_mid(I_cross+1))
    y0=dlog(max(dEL_mean(I_cross),tiny(1d0)))
    y1=dlog(max(dEL_mean(I_cross+1),tiny(1d0)))
    ytarget=dlog(coeff_target)
    if (y1 == y0) then
        xroot=0.5d0*(x0+x1)
    else
        xroot=x0+(ytarget-y0)*(x1-x0)/(y1-y0)
    end if
    gc=dexp(xroot)
end subroutine electron_gc_loss

! 计算从 R_left 开始的径向区间内平均单位半径注入源归一化。
! Compute the per-unit-radius injection normalization averaged over the interval starting at R_left.
subroutine electron_injection_prefactor(R_left,dDR,dNe,f_e,gmfac,Q)
    implicit none
    real(8), intent(in) :: R_left,dDR,dNe,f_e,gmfac
    real(8), intent(out) :: Q

    Q=4d0/3d0*pi*(3d0*R_left**2+dDR*(3d0*R_left+dDR))*dNe*f_e*gmfac
end subroutine electron_injection_prefactor

! 在电子网格上定位注入区间覆盖的单元范围。
! Locate the electron-grid cell range covered by the injection interval.
subroutine electron_source_bounds(ng,gam_e,gm,gemax,src_lo,src_hi)
    implicit none
    integer, intent(in) :: ng
    integer, intent(out) :: src_lo,src_hi
    integer :: ig
    real(8), intent(in), dimension(ng) :: gam_e
    real(8), intent(in) :: gm,gemax
    real(8), dimension(ng+1) :: x_edge
    real(8) :: x_lo,x_hi

    call log_edges(ng,gam_e,x_edge)
    x_lo=dlog(gm)
    x_hi=dlog(tail_factor*gemax)
    src_lo=ng+1
    do ig=1,ng
        if (x_edge(ig+1) > x_lo) then
            src_lo=ig
            exit
        end if
    end do
    src_hi=0
    do ig=ng,1,-1
        if (x_edge(ig) < x_hi) then
            src_hi=ig
            exit
        end if
    end do
    if (src_hi < src_lo) then
        src_lo=1
        src_hi=0
    end if
end subroutine electron_source_bounds

! 计算两个数组之间的最大相对误差。
! Compute the maximum relative error between 2 arrays.
subroutine electron_relerr_max(ng,x_ref,x_trial,error_max)
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: x_ref,x_trial
    real(8), intent(out) :: error_max
    integer :: ig
    real(8) :: denom

    error_max=0d0
    do ig=1,ng
        denom=max(abs(x_ref(ig)),1d-99)
        error_max=max(error_max,abs(x_trial(ig)-x_ref(ig))/denom)
    end do
end subroutine electron_relerr_max

! 按风介质或均匀介质条件初始化外部电子密度和累计电子数。
! Initialize external electron density and cumulative electron count for wind or ISM media.
subroutine electron_initial_density(A_star,dNe_ISM,R_ini,R_start,R0,dNe,ninit)
    implicit none
    real(8), intent(in) :: A_star,dNe_ISM,R_ini,R_start,R0
    real(8), intent(out) :: dNe,ninit
    real(8) :: dNe_wind,wind_norm,r_floor,r_cap,r_join,r_base

    if (A_star > 0d0) then
        wind_norm=A_star*3.0d35
        r_floor=dsqrt(4d0*wind_norm/dNe_ISM)
        if (R0 > 0d0) then
            r_cap=min(R_ini,R0)
            ninit=4d0/3d0*pi*r_cap**3*wind_norm/R0**2
            if (R_ini > R0) then
                r_join=min(R_ini,r_floor)
                if (r_join > R0) ninit=ninit+4d0*pi*wind_norm*(r_join-R0)
                r_base=max(R0,r_floor)
                if (R_ini > r_base) ninit=ninit+4d0/3d0*pi*dNe_ISM*(R_ini**3-r_base**3)
            end if
        else
            r_join=min(R_ini,r_floor)
            ninit=4d0*pi*wind_norm*r_join
            if (R_ini > r_join) ninit=ninit+4d0/3d0*pi*dNe_ISM*(R_ini**3-r_join**3)
        end if
        if (R0 > 0d0 .and. R_start < R0) then
            dNe=wind_norm/R0**2
        else
            dNe_wind=wind_norm/R_start**2
            if (dNe_wind <= dNe_ISM/4d0) then
                dNe=dNe_ISM
            else
                dNe=dNe_wind
            end if
        end if
    else
        dNe=dNe_ISM
        ninit=4d0/3d0*pi*R_ini**3*dNe_ISM
    end if
end subroutine electron_initial_density

end module electron_common
