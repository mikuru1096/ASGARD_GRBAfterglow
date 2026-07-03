! 计算同步自康普顿(SSC)辐射谱。
! Compute the synchrotron self-Compton spectrum.
!
! 均匀电子网格路径: 对电子谱和种子光子场做双重 Simpson 积分。
! Uniform electron grid path: nested Simpson integration over electrons and seed photons.
subroutine ssc_spec(R,gam_e,dN_gam_e,V_seed,seed,Num_nu,Num_R,Num_gam_e,n_threads, P_SSC_spec,seed_SSC)
    use constants
    use rad_common
    !$ use omp_lib
    implicit none
    !***********************************************************
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,n_threads
    real(8), dimension(Num_R), intent(in) :: R
    real(8), dimension(Num_gam_e), intent(in) :: gam_e
    real(8), dimension(Num_gam_e,Num_R), intent(in) :: dN_gam_e
    real(8), dimension(Num_nu), intent(in) :: V_seed
    real(8), dimension(Num_nu,Num_R), intent(in) :: seed
    real(8), dimension(Num_nu,Num_R), intent(out) :: P_SSC_spec,seed_SSC

    real(8), dimension(:), allocatable :: simpson_weights,V_weights,E_seed,inv_gam,inv_gam2,radius_inv2
    real(8), dimension(:), allocatable :: vinv,logvinv
    real(8), dimension(:,:), allocatable :: q_pref,kn_pref,logq_pref,dnwg,weighted_seed
    real(8), dimension(:,:), allocatable :: tail_dn,tail_inv2
    integer, dimension(:), allocatable :: gamma_start
    integer, dimension(:,:), allocatable :: gamma_low,gamma_high
    integer :: i_r,i_nu
    real(8) :: para_hEme,h_nu,h_gam,dnu3,dgam3,tmp,tmp2

    allocate (simpson_weights(Num_gam_e), V_weights(Num_nu))
    allocate (E_seed(Num_nu), inv_gam(Num_gam_e), inv_gam2(Num_gam_e))
    allocate (radius_inv2(Num_R), vinv(Num_nu), logvinv(Num_nu))
    allocate (q_pref(Num_gam_e,Num_nu), kn_pref(Num_gam_e,Num_nu), logq_pref(Num_gam_e,Num_nu))
    allocate (dnwg(Num_gam_e,Num_R), weighted_seed(Num_nu,Num_R))
    allocate (tail_dn(Num_gam_e+1,Num_R), tail_inv2(Num_gam_e+1,Num_R))
    allocate (gamma_start(Num_nu), gamma_low(Num_nu,Num_nu), gamma_high(Num_nu,Num_nu))

    para_hEme = Para_h/para_m_energy

    h_nu = log(V_seed(2))-log(V_seed(1))
    h_gam = log(gam_e(2))-log(gam_e(1))
    E_seed = V_seed*para_hEme
    inv_gam = 1d0/gam_e
    inv_gam2 = inv_gam*inv_gam
    vinv = 1d0/V_seed
    logvinv = log(vinv)
    radius_inv2 = 1d0/(R*R)
    q_pref = 0d0
    kn_pref = 0d0
    logq_pref = 0d0

    call uniform_weights()
    call kn_tables()
    call uniform_bounds()

    P_SSC_spec=0d0
    seed_SSC=0d0

    dnu3 = h_nu/3d0
    dgam3 = h_gam/3d0

    !$OMP PARALLEL num_threads(n_threads), private(i_r, i_nu)
    !$OMP DO COLLAPSE(2) SCHEDULE(GUIDED,4)
    do i_r=1,Num_R
        do i_nu=1,Num_nu
            call uniform_point(i_r, i_nu)
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    tmp=0.75d0*Para_c*Para_h*Para_SigmaT
    P_SSC_spec=P_SSC_spec*tmp

    tmp2=4d0*pi*Para_c*Para_h
    seed_SSC=seed_SSC/tmp2*tmp

    deallocate(simpson_weights, V_weights, E_seed, inv_gam, inv_gam2)
    deallocate(radius_inv2, vinv, logvinv)
    deallocate(q_pref, kn_pref, logq_pref)
    deallocate(dnwg, weighted_seed, tail_dn, tail_inv2)
    deallocate(gamma_start, gamma_low, gamma_high)

    return

contains

integer function first_gt(gfloor, lower)
    implicit none
    real(8), intent(in) :: gfloor
    integer, intent(in) :: lower
    integer :: i_low, i_high, i_mid

    i_low=lower
    i_high=Num_gam_e+1
    do while (i_low < i_high)
        i_mid=(i_low+i_high)/2
        if (gam_e(i_mid) <= gfloor) then
            i_low=i_mid+1
        else
            i_high=i_mid
        end if
    end do
    first_gt=i_low
end function first_gt

integer function first_ge(gfloor, lower)
    implicit none
    real(8), intent(in) :: gfloor
    integer, intent(in) :: lower
    integer :: i_low, i_high, i_mid

    i_low=lower
    i_high=Num_gam_e+1
    do while (i_low < i_high)
        i_mid=(i_low+i_high)/2
        if (gam_e(i_mid) < gfloor) then
            i_low=i_mid+1
        else
            i_high=i_mid
        end if
    end do
    first_ge=i_low
end function first_ge

subroutine uniform_weights()
    implicit none
    integer :: i_shell, i_gamma

    call compute_simpson_weights(simpson_weights, Num_gam_e)
    call compute_simpson_weights(V_weights, Num_nu)
    do i_shell=1,Num_R
        dnwg(:,i_shell)=dN_gam_e(:,i_shell)*simpson_weights*inv_gam
        weighted_seed(:,i_shell)=seed(:,i_shell)*V_weights
        tail_dn(Num_gam_e+1,i_shell)=0d0
        tail_inv2(Num_gam_e+1,i_shell)=0d0
        do i_gamma=Num_gam_e,1,-1
            tail_dn(i_gamma,i_shell)=tail_dn(i_gamma+1,i_shell)+dnwg(i_gamma,i_shell)
            tail_inv2(i_gamma,i_shell)=tail_inv2(i_gamma+1,i_shell)+dnwg(i_gamma,i_shell)*inv_gam2(i_gamma)
        end do
    end do
end subroutine uniform_weights

subroutine kn_tables()
    implicit none
    integer :: i_seed, i_gamma
    real(8) :: temp_norm, q_gamma

    do i_seed=1,Num_nu
        gamma_start(i_seed)=first_gt(E_seed(i_seed), 1)
        do i_gamma=gamma_start(i_seed),Num_gam_e
            temp_norm=gam_e(i_gamma)-E_seed(i_seed)
            q_pref(i_gamma,i_seed)=V_seed(i_seed)/(4d0*gam_e(i_gamma)*temp_norm)
            q_gamma=E_seed(i_seed)/temp_norm
            kn_pref(i_gamma,i_seed)=q_gamma*q_gamma/(2d0*(1d0+q_gamma))
            logq_pref(i_gamma,i_seed)=log(q_pref(i_gamma,i_seed))
        end do
    end do
end subroutine kn_tables

subroutine uniform_bounds()
    implicit none
    integer :: i_obs, i_seed
    real(8) :: gamma_floor

    do i_obs=1,Num_nu
        do i_seed=1,i_obs-1
            gamma_floor=0.5d0*(E_seed(i_obs)+sqrt(E_seed(i_obs)*E_seed(i_obs)+V_seed(i_obs)*vinv(i_seed)))
            gamma_high(i_seed,i_obs)=first_gt(gamma_floor, gamma_start(i_obs))
        end do
        do i_seed=i_obs,Num_nu
            gamma_floor=0.5d0*sqrt(V_seed(i_seed)/V_seed(i_obs))
            gamma_low(i_seed,i_obs)=first_ge(gamma_floor, 1)
        end do
    end do
end subroutine uniform_bounds

real(8) function low_sum(i_shell, i_obs)
    implicit none
    integer, intent(in) :: i_shell, i_obs
    integer :: i_seed, i_gamma
    real(8) :: q_coeff, q, logq, fssc, gsum, emiss

    low_sum=0d0
    do i_seed=1,i_obs-1
        gsum=0d0
        !$OMP SIMD REDUCTION(+:gsum)
        do i_gamma=gamma_high(i_seed,i_obs),Num_gam_e
            q_coeff=q_pref(i_gamma,i_obs)
            q=q_coeff*vinv(i_seed)
            if (q >= 1d0) cycle
            logq=logq_pref(i_gamma,i_obs)+logvinv(i_seed)
            fssc=2d0*q*(logq-q)+1d0+q+kn_pref(i_gamma,i_obs)*(1d0-q)
            gsum=gsum+dnwg(i_gamma,i_shell)*fssc
        end do
        emiss=dgam3*gsum
        low_sum=low_sum+weighted_seed(i_seed,i_shell)*emiss
    end do
end function low_sum

real(8) function high_tail(i_shell, i_obs)
    implicit none
    integer, intent(in) :: i_shell, i_obs
    integer :: i_seed, i_gamma
    real(8) :: ratio_v, gsum, emiss

    high_tail=0d0
    do i_seed=i_obs,Num_nu
        ratio_v=V_seed(i_obs)*vinv(i_seed)
        i_gamma=gamma_low(i_seed,i_obs)
        if (i_gamma <= Num_gam_e) then
            gsum=ratio_v*tail_dn(i_gamma,i_shell)-0.25d0*tail_inv2(i_gamma,i_shell)
        else
            gsum=0d0
        end if
        emiss=dgam3*gsum
        high_tail=high_tail+weighted_seed(i_seed,i_shell)*emiss
    end do
end function high_tail

subroutine uniform_point(i_shell, i_obs)
    implicit none
    integer, intent(in) :: i_shell, i_obs
    real(8) :: integ, power, flux

    if (gamma_start(i_obs) > Num_gam_e) return

    integ=dnu3*(low_sum(i_shell, i_obs)+high_tail(i_shell, i_obs))
    power=integ*V_seed(i_obs)
    P_SSC_spec(i_obs,i_shell)=P_SSC_spec(i_obs,i_shell)+power
    flux=integ*radius_inv2(i_shell)
    seed_SSC(i_obs,i_shell)=seed_SSC(i_obs,i_shell)+flux
end subroutine uniform_point

end subroutine ssc_spec

! 计算非均匀网格 SSC 辐射谱。
! Compute the SSC spectrum on a non-uniform electron grid.
!
! 电子谱用分段线性重构，种子光子积分用二点 Gauss-Legendre，斜率用 minmod 限制。
! The electron spectrum is linearly reconstructed; seed photons use 2-point
! Gauss-Legendre integration; slopes are minmod-limited.
subroutine ssc_spec_nonuniform(R,x_edge_log10,dN_x,V_seed,seed,Num_nu,Num_R,Num_gam_e,n_threads, &
                               P_SSC_spec,seed_SSC)
    use constants
    use rad_common
    !$ use omp_lib
    implicit none
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,n_threads
    real(8), dimension(Num_R), intent(in) :: R
    real(8), dimension(Num_gam_e+1,Num_R), intent(in) :: x_edge_log10
    real(8), dimension(Num_gam_e,Num_R), intent(in) :: dN_x
    real(8), dimension(Num_nu), intent(in) :: V_seed
    real(8), dimension(Num_nu,Num_R), intent(in) :: seed
    real(8), dimension(Num_nu,Num_R), intent(out) :: P_SSC_spec,seed_SSC

    real(8), dimension(:), allocatable :: x_seed,radius_inv2
    real(8), dimension(:,:), allocatable :: x_center,slope_q,tail_gamma,tail_inv2
    real(8) :: para_hEme,w2,tmp,tmp2
    integer :: i_r,i_nu

    allocate(x_seed(Num_nu),radius_inv2(Num_R),x_center(Num_gam_e,Num_R),slope_q(Num_gam_e,Num_R))
    allocate(tail_gamma(Num_gam_e+1,Num_R),tail_inv2(Num_gam_e+1,Num_R))

    para_hEme=Para_h/para_m_energy
    x_seed=dlog(V_seed)
    radius_inv2=1d0/(R*R)
    w2=1d0/dsqrt(3d0)

    call build_state()

    P_SSC_spec=0d0
    seed_SSC=0d0

    !$OMP PARALLEL num_threads(n_threads), private(i_r,i_nu)
    !$OMP DO COLLAPSE(2) SCHEDULE(GUIDED,4)
    do i_r=1,Num_R
        do i_nu=1,Num_nu
            call nonuniform_point(i_r, i_nu)
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    tmp=0.75d0*Para_c*Para_h*Para_SigmaT
    P_SSC_spec=P_SSC_spec*tmp

    tmp2=4d0*pi*Para_c*Para_h
    seed_SSC=seed_SSC/tmp2*tmp

    deallocate(x_seed,radius_inv2,x_center,slope_q)
    deallocate(tail_gamma,tail_inv2)

contains

subroutine build_state()
    implicit none
    integer :: i_shell, i_gamma

    do i_shell=1,Num_R
        call build_shell(i_shell)
        tail_gamma(Num_gam_e+1,i_shell)=0d0
        tail_inv2(Num_gam_e+1,i_shell)=0d0
        do i_gamma=Num_gam_e,1,-1
            call add_tail(i_shell, i_gamma)
        end do
    end do
end subroutine build_state

subroutine build_shell(i_shell)
    implicit none
    integer, intent(in) :: i_shell
    integer :: i_gamma

    do i_gamma=1,Num_gam_e
        if (x_edge_log10(i_gamma+1,i_shell) <= x_edge_log10(i_gamma,i_shell)) then
            error stop "ssc_spec_nonuniform: electron edge grid must be strictly increasing in every shell."
        end if
    end do
    do i_gamma=1,Num_gam_e
        x_center(i_gamma,i_shell)=0.5d0*(x_edge_log10(i_gamma,i_shell)+x_edge_log10(i_gamma+1,i_shell))
    end do
    do i_gamma=1,Num_gam_e
        call limit_slope(i_shell, i_gamma)
    end do
end subroutine build_shell

subroutine limit_slope(i_shell, i_gamma)
    implicit none
    integer, intent(in) :: i_shell, i_gamma
    real(8) :: dx_log10, left_slope, right_slope, slope_lim

    dx_log10=x_edge_log10(i_gamma+1,i_shell)-x_edge_log10(i_gamma,i_shell)
    if (Num_gam_e == 1) then
        slope_q(i_gamma,i_shell)=0d0
    else if (i_gamma == 1) then
        right_slope=(dN_x(2,i_shell)-dN_x(1,i_shell))/(x_center(2,i_shell)-x_center(1,i_shell))
        slope_q(i_gamma,i_shell)=right_slope
    else if (i_gamma == Num_gam_e) then
        left_slope=(dN_x(Num_gam_e,i_shell)-dN_x(Num_gam_e-1,i_shell))/ &
                   (x_center(Num_gam_e,i_shell)-x_center(Num_gam_e-1,i_shell))
        slope_q(i_gamma,i_shell)=left_slope
    else
        left_slope=(dN_x(i_gamma,i_shell)-dN_x(i_gamma-1,i_shell))/ &
                   (x_center(i_gamma,i_shell)-x_center(i_gamma-1,i_shell))
        right_slope=(dN_x(i_gamma+1,i_shell)-dN_x(i_gamma,i_shell))/ &
                    (x_center(i_gamma+1,i_shell)-x_center(i_gamma,i_shell))
        slope_q(i_gamma,i_shell)=ssc_minmod(left_slope,right_slope)
    end if
    slope_lim=2d0*dN_x(i_gamma,i_shell)/dx_log10
    if (abs(slope_q(i_gamma,i_shell)) > slope_lim) then
        slope_q(i_gamma,i_shell)=sign(slope_lim,slope_q(i_gamma,i_shell))
    end if
end subroutine limit_slope

subroutine add_tail(i_shell, i_gamma)
    implicit none
    integer, intent(in) :: i_shell, i_gamma
    real(8) :: x0_cell, x1_cell

    x0_cell=x_edge_log10(i_gamma,i_shell)
    x1_cell=x_edge_log10(i_gamma+1,i_shell)
    tail_gamma(i_gamma,i_shell)=tail_gamma(i_gamma+1,i_shell)+ &
        gamma_moment(x0_cell,x1_cell,dN_x(i_gamma,i_shell),slope_q(i_gamma,i_shell),x_center(i_gamma,i_shell),2d0)
    tail_inv2(i_gamma,i_shell)=tail_inv2(i_gamma+1,i_shell)+ &
        gamma_moment(x0_cell,x1_cell,dN_x(i_gamma,i_shell),slope_q(i_gamma,i_shell),x_center(i_gamma,i_shell),4d0)
end subroutine add_tail

real(8) function gauss_point(xm, dx, i_gauss)
    implicit none
    real(8), intent(in) :: xm, dx
    integer, intent(in) :: i_gauss

    if (i_gauss == 1) then
        gauss_point=xm-dx*w2
    else
        gauss_point=xm+dx*w2
    end if
end function gauss_point

subroutine nonuniform_point(i_shell, i_obs)
    implicit none
    integer, intent(in) :: i_shell, i_obs
    integer :: i_seed, i_gauss
    real(8) :: ratio_v, gfloor, integ, emiss, power, flux, vobs, eobs
    real(8) :: x0, x1, xm, dx, xloc, vseed
    real(8) :: sloc, swgt

    vobs=V_seed(i_obs)
    eobs=vobs*para_hEme
    integ=0d0
    do i_seed=1,Num_nu-1
        x0=x_seed(i_seed)
        x1=x_seed(i_seed+1)
        xm=0.5d0*(x0+x1)
        dx=0.5d0*(x1-x0)
        do i_gauss=1,2
            xloc=gauss_point(xm, dx, i_gauss)
            vseed=dexp(xloc)
            sloc=rad_interp(V_seed(i_seed),V_seed(i_seed+1),seed(i_seed,i_shell),seed(i_seed+1,i_shell),vseed)
            if (sloc <= 0d0) cycle
            swgt=dx*sloc
            if (vseed < vobs) then
                gfloor=0.5d0*(eobs+dsqrt(eobs*eobs+vobs/vseed))
                emiss=low_integral(i_shell,dlog10(gfloor),vobs,vseed,eobs)
            else
                ratio_v=vobs/vseed
                gfloor=0.5d0*dsqrt(vseed/vobs)
                emiss=high_tail(i_shell,dlog10(gfloor),ratio_v)
            end if
            integ=integ+swgt*emiss
        end do
    end do
    power=integ*vobs
    P_SSC_spec(i_obs,i_shell)=power
    flux=integ*radius_inv2(i_shell)
    seed_SSC(i_obs,i_shell)=flux
end subroutine nonuniform_point

! 二分查找第一个右边界超过 x_floor 的电子单元。
! Find the first electron cell whose right edge exceeds x_floor.
integer function first_cell(x_edge_col,n,x_floor)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n+1), intent(in) :: x_edge_col
    real(8), intent(in) :: x_floor
    integer :: left,right,mid

    if (x_floor <= x_edge_col(1)) then
        first_cell=1
        return
    end if
    if (x_floor >= x_edge_col(n+1)) then
        first_cell=n+1
        return
    end if

    left=1
    right=n
    do while (left < right)
        mid=(left+right)/2
        if (x_edge_col(mid+1) > x_floor) then
            right=mid
        else
            left=mid+1
        end if
    end do
    first_cell=left
end function first_cell

! minmod 斜率限制器: 同号取绝对值最小者，异号返回零。
! Minmod slope limiter: keep the smaller same-sign slope, otherwise zero.
real(8) function ssc_minmod(a,b)
    implicit none
    real(8), intent(in) :: a,b

    if (a*b <= 0d0) then
        ssc_minmod=0d0
    else
        ssc_minmod=sign(min(abs(a),abs(b)),a)
    end if
end function ssc_minmod

! 线性重构剖面 q(x) = qbar + slope*(x - xc)。
! Linear reconstructed profile q(x) = qbar + slope*(x - xc).
real(8) function profile_value(x,qbar,slope,xc)
    implicit none
    real(8), intent(in) :: x,qbar,slope,xc
    profile_value=qbar+slope*(x-xc)
    if (profile_value < 0d0) then
        error stop "ssc_spec_nonuniform: linear reconstruction became negative."
    end if
end function profile_value

! 线性重构剖面在 [x0,x1] 上的 gamma^(-power) 矩积分。
! Analytic gamma^(-power) moment of the linear reconstructed profile on [x0,x1].
real(8) function gamma_moment(x0,x1,qbar,slope,xc,power)
    implicit none
    real(8), intent(in) :: x0,x1,qbar,slope,xc,power
    real(8) :: alpha,exp0,exp1,i0,i1

    if (x1 <= x0) then
        gamma_moment=0d0
        return
    end if

    alpha=power*dlog(1d1)
    exp0=1d1**(-power*x0)
    exp1=1d1**(-power*x1)
    i0=(exp0-exp1)/alpha
    i1=((x0/alpha)+1d0/(alpha*alpha))*exp0-((x1/alpha)+1d0/(alpha*alpha))*exp1
    gamma_moment=qbar*i0+slope*(i1-xc*i0)
end function gamma_moment

! 低能种子区(nu_seed < nu_obs)单个电子单元的 SSC 散射积分。
! SSC scattering integral per electron cell in the low-seed branch.
real(8) function low_cell(x0,x1,qbar,slope,xc,vobs,vseed,eobs)
    implicit none
    real(8), intent(in) :: x0,x1,qbar,slope,xc,vobs,vseed,eobs
    real(8) :: xm,dx,w2_loc,xg,gam,temp,q,q_gamma,kn_coeff,fssc,val,q_loc
    integer :: i_pt

    if (x1 <= x0 .or. qbar == 0d0) then
        low_cell=0d0
        return
    end if

    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    w2_loc=1d0/dsqrt(3d0)
    low_cell=0d0
    do i_pt=1,2
        if (i_pt == 1) then
            xg=xm-dx*w2_loc
        else
            xg=xm+dx*w2_loc
        end if
        q_loc=profile_value(xg,qbar,slope,xc)
        if (q_loc <= 0d0) cycle
        gam=1d1**xg
        temp=gam-eobs
        if (temp <= 0d0) cycle
        q=vobs/(4d0*gam*temp*vseed)
        if (q <= 0d0 .or. q >= 1d0) cycle
        q_gamma=eobs/temp
        kn_coeff=q_gamma*q_gamma/(2d0*(1d0+q_gamma))
        fssc=2d0*q*(dlog(q)-q)+1d0+q+kn_coeff*(1d0-q)
        val=q_loc*fssc/(gam*gam)
        low_cell=low_cell+val
    end do
    low_cell=low_cell*dx
end function low_cell

! 低能种子区(nu_seed < nu_obs)的完整 SSC 散射积分。
! Full low-seed SSC integral over the electron grid.
real(8) function low_integral(i_r,x_floor,vobs,vseed,eobs)
    implicit none
    integer, intent(in) :: i_r
    real(8), intent(in) :: x_floor,vobs,vseed,eobs
    integer :: i_start,i_gam
    real(8) :: x0,x1

    i_start=first_cell(x_edge_log10(:,i_r),Num_gam_e,x_floor)
    if (i_start > Num_gam_e) then
        low_integral=0d0
        return
    end if

    low_integral=0d0
    do i_gam=i_start,Num_gam_e
        x0=x_edge_log10(i_gam,i_r)
        if (x_floor > x0) x0=x_floor
        x1=x_edge_log10(i_gam+1,i_r)
        low_integral=low_integral+low_cell(x0,x1,dN_x(i_gam,i_r),slope_q(i_gam,i_r),x_center(i_gam,i_r), &
                                           vobs,vseed,eobs)
    end do
end function low_integral

! 高能种子区(nu_seed >= nu_obs)SSC 尾部积分，使用预计算 gamma 矩。
! High-seed SSC tail integral using precomputed gamma moments.
real(8) function high_tail(i_r,x_floor,ratio_v)
    implicit none
    integer, intent(in) :: i_r
    real(8), intent(in) :: x_floor,ratio_v
    integer :: i_start
    real(8) :: x0,x1,part2,part4

    i_start=first_cell(x_edge_log10(:,i_r),Num_gam_e,x_floor)
    if (i_start > Num_gam_e) then
        high_tail=0d0
        return
    end if

    x0=x_edge_log10(i_start,i_r)
    if (x_floor > x0) x0=x_floor
    x1=x_edge_log10(i_start+1,i_r)
    part2=gamma_moment(x0,x1,dN_x(i_start,i_r),slope_q(i_start,i_r),x_center(i_start,i_r),2d0)
    part4=gamma_moment(x0,x1,dN_x(i_start,i_r),slope_q(i_start,i_r),x_center(i_start,i_r),4d0)
    high_tail=ratio_v*(part2+tail_gamma(i_start+1,i_r))-0.25d0*(part4+tail_inv2(i_start+1,i_r))
end function high_tail

end subroutine ssc_spec_nonuniform
