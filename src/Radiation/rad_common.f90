!f2py: skip
module rad_common
    use constants
    implicit none
    private

    public :: sampled_weights, rad_interp, transfer_factor, &
              syn_kernel, pair_grid, &
              pair_sigma, pair_delta, pair_tau, &
              syn_seed_chi, syn_flux_chi

contains

! 实际节点上的 Simpson 权重；偶数样点用 Cartwright 末区间校正。
! Simpson weights on sampled nodes, with Cartwright's last-interval correction for even n.
pure subroutine sampled_weights(x, w, n)
    integer, intent(in) :: n
    integer :: i,last
    real(8), dimension(n), intent(in) :: x
    real(8), dimension(n), intent(out) :: w
    real(8) :: h0,h1,hsum,alpha,beta,eta

    w=0d0
    if (n == 2) then
        h0=0.5d0*(x(2)-x(1))
        w=h0
        return
    end if

    last=n
    if (mod(n,2) == 0) last=n-1
    do i=1,last-2,2
        h0=x(i+1)-x(i)
        h1=x(i+2)-x(i+1)
        hsum=h0+h1
        w(i)=w(i)+hsum*(2d0-h1/h0)/6d0
        w(i+1)=w(i+1)+hsum**3/(6d0*h0*h1)
        w(i+2)=w(i+2)+hsum*(2d0-h0/h1)/6d0
    end do

    if (last < n) then
        h0=x(n-1)-x(n-2)
        h1=x(n)-x(n-1)
        hsum=h0+h1
        alpha=(2d0*h1*h1+3d0*h0*h1)/(6d0*hsum)
        beta=(h1*h1+3d0*h0*h1)/(6d0*h0)
        eta=h1**3/(6d0*h0*hsum)
        w(n)=w(n)+alpha
        w(n-1)=w(n-1)+beta
        w(n-2)=w(n-2)-eta
    end if
end subroutine sampled_weights

! 幂律插值: 两端为正时用 log-log，否则用线性插值。
! Power-law interpolation: log-log for positive endpoints, otherwise linear.
real(8) function rad_interp(v0,v1,y0,y1,v)
    implicit none
    real(8), intent(in) :: v0,v1,y0,y1,v
    real(8) :: slope

    if (v1 <= v0) then
        error stop "rad_interp requires a strictly increasing frequency interval."
    end if

    if (v <= v0) then
        rad_interp=y0
        return
    end if
    if (v >= v1) then
        rad_interp=y1
        return
    end if

    if (y0 > 0d0 .and. y1 > 0d0) then
        slope=dlog(y1/y0)/dlog(v1/v0)
        rad_interp=y0*(v/v0)**slope
    else
        rad_interp=y0+(y1-y0)*(v-v0)/(v1-v0)
    end if
end function rad_interp

! 辐射转移因子: (1-exp(-tau))/tau；五阶级数是同一解析函数在零点的展开。
! Transfer factor: the fifth-order branch is the analytic Taylor expansion at tau=0.
subroutine transfer_factor(Tau, factor)
    implicit none
    real(8), intent(in) :: Tau
    real(8), intent(out) :: factor
    if (dabs(Tau) < 1d-3) then
        factor=1d0+Tau*(-0.5d0+Tau*(1d0/6d0+Tau*(-1d0/24d0+Tau*(1d0/120d0-Tau/720d0))))
    else
        factor=(1d0-dexp(-Tau))/Tau
    end if
end subroutine transfer_factor

elemental real(8) function syn_kernel(x, ratio_v_pow, factor) result(Fx)
    implicit none
    real(8), intent(in) :: x, ratio_v_pow, factor

    ! 同一近似同步核的低频渐近、指数段和高频尾分段。
    ! Same approximate synchrotron kernel split into low, main, and high-frequency branches.
    if (x > 30d0) then
        Fx = 0d0
    else if (x < 3d-6) then
        Fx = 1.81d0/dsqrt(ratio_v_pow)
    else
        Fx = 1.81d0*dexp(-x)/dsqrt(ratio_v_pow+factor)
    end if
end function syn_kernel

! 准备湮灭计算网格: 归一化光子能量、中心频率和体积元。
! Prepare pair-production grid: normalized photon energy, midpoint frequency, and measure.
subroutine pair_grid(V_seed, Num_nu, ep1, ep2, dVloc, V_mid)
    implicit none
    integer, intent(in) :: Num_nu
    real(8), dimension(Num_nu), intent(in) :: V_seed
    real(8), dimension(1,Num_nu), intent(out) :: ep1
    real(8), dimension(Num_nu-1,1), intent(out) :: ep2
    real(8), dimension(Num_nu-1), intent(out) :: dVloc,V_mid
    real(8), dimension(Num_nu) :: x_seed
    real(8) :: para_hEme

    para_hEme=Para_h/para_m_energy
    ep1(1,:)=para_hEme*V_seed
    x_seed=dlog(V_seed)
    V_mid=dexp(0.5d0*(x_seed(1:Num_nu-1)+x_seed(2:Num_nu)))
    ep2(:,1)=para_hEme*V_mid
    dVloc=V_mid*(x_seed(2:Num_nu)-x_seed(1:Num_nu-1))
end subroutine pair_grid

! 光子-光子对产生截面。
! Photon-photon pair-production cross section.
elemental real(8) function pair_sigma(s_center) result(sigma_pair)
    implicit none
    real(8), intent(in) :: s_center

    if (s_center <= 1d0) then
        sigma_pair = 0d0
        return
    end if

    sigma_pair=pair_delta(s_center-1d0)
end function pair_sigma

! Evaluate the Breit-Wheeler cross section from s-1, retaining beta^2 at the nearest threshold.
elemental real(8) function pair_delta(delta) result(sigma_pair)
    implicit none
    real(8), intent(in) :: delta
    real(8) :: inv_s,beta_sq,beta_loc,log_s,log_term

    inv_s=1d0/(1d0+delta)
    if (delta < 1d0) then
        beta_sq=delta*inv_s
        log_s=2d0*atanh(delta/(2d0+delta))
    else
        beta_sq=1d0-inv_s
        log_s=dlog(1d0+delta)
    end if
    beta_loc=dsqrt(beta_sq)
    log_term=log_s+4d0*atanh(beta_loc/(2d0+beta_loc))
    sigma_pair=3d0*Para_sigmaT*inv_s*((3d0-beta_sq*beta_sq)*log_term- &
               2d0*beta_loc*(2d0-beta_sq))/16d0
end function pair_delta

! 对头碰撞近似下计算光子-光子对产生光深。
! Pair-production optical depth in the head-on approximation.
subroutine pair_tau(V_seed,Num_nu,Seed_target,dx_cm,Tau_pair)
    implicit none
    integer, intent(in) :: Num_nu
    real(8), dimension(Num_nu), intent(in) :: V_seed,Seed_target
    real(8), intent(in) :: dx_cm
    real(8), dimension(Num_nu), intent(out) :: Tau_pair
    integer :: i_nu,i_seg
    real(8) :: para_hEme,ep1,v_mid,dv_mid,seed_mid,s_center
    real(8), dimension(Num_nu) :: x_seed

    Tau_pair = 0d0
    if (dx_cm <= 0d0) return

    para_hEme = Para_h/para_m_energy
    x_seed = dlog(V_seed)

    do i_nu = 1, Num_nu
        ep1 = para_hEme*V_seed(i_nu)
        do i_seg = 1, Num_nu-1
            v_mid = dexp(0.5d0*(x_seed(i_seg)+x_seed(i_seg+1)))
            dv_mid = v_mid*(x_seed(i_seg+1)-x_seed(i_seg))
            seed_mid = rad_interp(V_seed(i_seg),V_seed(i_seg+1),Seed_target(i_seg),Seed_target(i_seg+1),v_mid)
            s_center = ep1*para_hEme*v_mid
            Tau_pair(i_nu) = Tau_pair(i_nu) + pair_sigma(s_center)*seed_mid*dv_mid
        end do
        Tau_pair(i_nu) = dx_cm*Tau_pair(i_nu)
    end do
end subroutine pair_tau

! chi 批量同步辐射核心: 同一半径且 chi 上磁场相同时复用 F(nu/nu_c) 核。
! Chi-batch synchrotron core: reuse F(nu/nu_c) when B is shared across chi.
! 中文：Q_weight 只缩放每个 chi cell 的径向 SSA 光深；发射功率仍保持 cell-local 电子谱定义。
! English: Q_weight scales only the radial SSA optical depth of each chi cell;
! emitted power remains defined by the cell-local electron spectrum.
subroutine syn_seed_chi(R_loc,Num_gam_e,Num_nu,Num_chi,gam_e,DNe_chi,V_seed,DB_chi,Q_weight,ssa_prefactor, &
                                             P_emit,P_syn,Seed_syn,Tau_syn)
    implicit none
    integer, intent(in) :: Num_gam_e,Num_nu,Num_chi
    real(8), dimension(Num_gam_e), intent(in) :: gam_e
    real(8), dimension(Num_gam_e,Num_chi), intent(in) :: DNe_chi
    real(8), dimension(Num_nu), intent(in) :: V_seed
    real(8), dimension(Num_chi), intent(in) :: DB_chi
    real(8), dimension(Num_chi), intent(in) :: Q_weight
    real(8), intent(in) :: R_loc,ssa_prefactor
    real(8), dimension(Num_nu,Num_chi), intent(out) :: P_emit,P_syn,Seed_syn,Tau_syn
    real(8) :: factor,temp,r2,norm,dbref,db,sumv,tau,power,trans,vc,fx,x
    real(8) :: vcinv,vcpow,vcm13
    real(8), dimension(Num_nu) :: intg,taug
    real(8) :: emit_w,tau_w
    real(8), dimension(Num_gam_e-1) :: gmean,dgam,vcarr,vpow23,ewgt,twgt
    real(8), dimension(Num_gam_e) :: ginv2
    real(8), dimension(Num_nu) :: vpm23,vp13
    real(8), dimension(Num_gam_e-1,Num_nu) :: fxgrid
    integer :: i_chi,i_gam,i_nu,ilow,itail
    logical :: uniform

    dbref = DB_chi(1)
    uniform = .true.
    do i_chi = 2, Num_chi
        if (DB_chi(i_chi) /= dbref) uniform = .false.
    end do

    factor = (3.62d0/pi)**2
    temp = dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    r2 = R_loc*R_loc
    norm = 4d0*pi*Para_c*Para_h
    do i_gam = 1, Num_gam_e
        ginv2(i_gam) = 1d0/(gam_e(i_gam)*gam_e(i_gam))
    end do
    do i_gam = 1, Num_gam_e-1
        gmean(i_gam) = (gam_e(i_gam)+gam_e(i_gam+1))**2/4d0
        dgam(i_gam) = (gam_e(i_gam+1)-gam_e(i_gam))/2d0
        vc = 4.2d6*gmean(i_gam)*dbref
        vcarr(i_gam) = 1d0/vc
        vpow23(i_gam) = vc**(2d0/3d0)
    end do
    do i_nu = 1, Num_nu
        vpm23(i_nu) = V_seed(i_nu)**(-2d0/3d0)
        vp13(i_nu) = V_seed(i_nu)**(1d0/3d0)
    end do
    if (.not. uniform) then
        do i_chi = 1, Num_chi
            db = DB_chi(i_chi)
            intg = 0d0
            taug = 0d0
            ilow = 0
            itail = 0
            do i_gam = 1, Num_gam_e-1
                emit_w = (DNe_chi(i_gam,i_chi)+DNe_chi(i_gam+1,i_chi))*dgam(i_gam)
                tau_w = gmean(i_gam)*(DNe_chi(i_gam,i_chi)*ginv2(i_gam) - DNe_chi(i_gam+1,i_chi)*ginv2(i_gam+1))
                vc = 4.2d6*gmean(i_gam)*db
                vcinv = 1d0/vc
                vcpow = vc**(2d0/3d0)
                vcm13 = 1d0/(vc**(1d0/3d0))
                do while (ilow < Num_nu .and. V_seed(ilow+1) < 3d-6*vc)
                    ilow = ilow+1
                end do
                if (itail < ilow) itail = ilow
                do while (itail < Num_nu .and. V_seed(itail+1) <= 30d0*vc)
                    itail = itail+1
                end do
                do i_nu = 1, ilow
                    fx = 1.81d0*vp13(i_nu)*vcm13
                    intg(i_nu) = intg(i_nu) + emit_w*fx
                    taug(i_nu) = taug(i_nu) + tau_w*fx
                end do
                do i_nu = ilow+1, itail
                    x = V_seed(i_nu)*vcinv
                    fx = syn_kernel(x,vcpow*vpm23(i_nu),factor)
                    intg(i_nu) = intg(i_nu) + emit_w*fx
                    taug(i_nu) = taug(i_nu) + tau_w*fx
                end do
            end do
            do i_nu = 1, Num_nu
                power = temp*db*intg(i_nu)
                tau = ssa_prefactor*taug(i_nu)*db*Q_weight(i_chi)/(4d0*pi*r2*V_seed(i_nu)*V_seed(i_nu))
                P_emit(i_nu,i_chi) = power
                Tau_syn(i_nu,i_chi) = tau
                call transfer_factor(tau,trans)
                P_syn(i_nu,i_chi) = power*trans
                Seed_syn(i_nu,i_chi) = P_syn(i_nu,i_chi)/(r2*V_seed(i_nu)*norm)
            end do
        end do
        return
    end if
    do i_nu = 1, Num_nu
        do i_gam = 1, Num_gam_e-1
            x = V_seed(i_nu)*vcarr(i_gam)
            fxgrid(i_gam,i_nu) = syn_kernel(x,vpow23(i_gam)*vpm23(i_nu),factor)
        end do
    end do

    do i_chi = 1, Num_chi
        do i_gam = 1, Num_gam_e-1
            ewgt(i_gam) = (DNe_chi(i_gam,i_chi)+DNe_chi(i_gam+1,i_chi))*dgam(i_gam)
            twgt(i_gam) = gmean(i_gam)*(DNe_chi(i_gam,i_chi)*ginv2(i_gam) - DNe_chi(i_gam+1,i_chi)*ginv2(i_gam+1))
        end do
        do i_nu = 1, Num_nu
            sumv = 0d0
            tau = 0d0
            do i_gam = 1, Num_gam_e-1
                sumv = sumv + ewgt(i_gam)*fxgrid(i_gam,i_nu)
                tau = tau + twgt(i_gam)*fxgrid(i_gam,i_nu)
            end do
            power = temp*dbref*sumv
            tau = ssa_prefactor*tau*dbref*Q_weight(i_chi)/(4d0*pi*r2*V_seed(i_nu)*V_seed(i_nu))
            P_emit(i_nu,i_chi) = power
            Tau_syn(i_nu,i_chi) = tau
            call transfer_factor(tau,trans)
            P_syn(i_nu,i_chi) = power*trans
            Seed_syn(i_nu,i_chi) = P_syn(i_nu,i_chi)/(r2*V_seed(i_nu)*norm)
        end do
    end do
end subroutine syn_seed_chi

subroutine syn_flux_chi(R_loc,Num_gam_e,Num_nu,Num_chi,gam_e,DNe_chi,V_seed,DB_chi,Q_weight, &
                                                 ssa_prefactor,P_syn,Tau_syn)
    implicit none
    integer, intent(in) :: Num_gam_e,Num_nu,Num_chi
    real(8), dimension(Num_gam_e), intent(in) :: gam_e
    real(8), dimension(Num_gam_e,Num_chi), intent(in) :: DNe_chi
    real(8), dimension(Num_nu), intent(in) :: V_seed
    real(8), dimension(Num_chi), intent(in) :: DB_chi
    real(8), dimension(Num_chi), intent(in) :: Q_weight
    real(8), intent(in) :: R_loc,ssa_prefactor
    real(8), dimension(Num_nu,Num_chi), intent(out) :: P_syn,Tau_syn
    real(8) :: factor,temp,r2,dbref,db,sumv,tau,power,vc,fx,x
    real(8), dimension(Num_nu) :: intg,taug
    real(8) :: emit_w,tau_w,vcinv,vcpow,vcm13
    real(8), dimension(Num_gam_e-1) :: gmean,dgam,vcarr,vpow23,ewgt,twgt
    real(8), dimension(Num_gam_e) :: ginv2
    real(8), dimension(Num_nu) :: vpm23,vp13
    real(8), dimension(Num_gam_e-1,Num_nu) :: fxgrid
    integer :: i_chi,i_gam,i_nu,ilow,itail
    logical :: uniform

    dbref = DB_chi(1)
    uniform = .true.
    do i_chi = 2, Num_chi
        if (DB_chi(i_chi) /= dbref) uniform = .false.
    end do
    factor = (3.62d0/pi)**2
    temp = dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    r2 = R_loc*R_loc
    do i_gam = 1, Num_gam_e
        ginv2(i_gam) = 1d0/(gam_e(i_gam)*gam_e(i_gam))
    end do
    do i_gam = 1, Num_gam_e-1
        gmean(i_gam) = (gam_e(i_gam)+gam_e(i_gam+1))**2/4d0
        dgam(i_gam) = (gam_e(i_gam+1)-gam_e(i_gam))/2d0
    end do
    do i_nu = 1, Num_nu
        vpm23(i_nu) = V_seed(i_nu)**(-2d0/3d0)
        vp13(i_nu) = V_seed(i_nu)**(1d0/3d0)
    end do

    if (uniform) then
        do i_gam = 1, Num_gam_e-1
            vc = 4.2d6*gmean(i_gam)*dbref
            vcarr(i_gam) = 1d0/vc
            vpow23(i_gam) = vc**(2d0/3d0)
        end do
        do i_nu = 1, Num_nu
            do i_gam = 1, Num_gam_e-1
                x = V_seed(i_nu)*vcarr(i_gam)
                fxgrid(i_gam,i_nu) = syn_kernel(x,vpow23(i_gam)*vpm23(i_nu),factor)
            end do
        end do
    end if

    if (uniform) then
        do i_chi = 1, Num_chi
            do i_gam = 1, Num_gam_e-1
                ewgt(i_gam) = (DNe_chi(i_gam,i_chi)+DNe_chi(i_gam+1,i_chi))*dgam(i_gam)
                twgt(i_gam) = gmean(i_gam)*(DNe_chi(i_gam,i_chi)*ginv2(i_gam) - DNe_chi(i_gam+1,i_chi)*ginv2(i_gam+1))
            end do
            do i_nu = 1, Num_nu
                sumv = 0d0
                tau = 0d0
                do i_gam = 1, Num_gam_e-1
                    sumv = sumv + ewgt(i_gam)*fxgrid(i_gam,i_nu)
                    tau = tau + twgt(i_gam)*fxgrid(i_gam,i_nu)
                end do
                power = temp*dbref*sumv
                tau = ssa_prefactor*tau*dbref*Q_weight(i_chi)/(4d0*pi*r2*V_seed(i_nu)*V_seed(i_nu))
                Tau_syn(i_nu,i_chi) = tau
                P_syn(i_nu,i_chi) = power
            end do
        end do
    else
        do i_chi = 1, Num_chi
            db = DB_chi(i_chi)
            intg = 0d0
            taug = 0d0
            ilow = 0
            itail = 0
            do i_gam = 1, Num_gam_e-1
                emit_w = (DNe_chi(i_gam,i_chi)+DNe_chi(i_gam+1,i_chi))*dgam(i_gam)
                tau_w = gmean(i_gam)*(DNe_chi(i_gam,i_chi)*ginv2(i_gam) - DNe_chi(i_gam+1,i_chi)*ginv2(i_gam+1))
                vc = 4.2d6*gmean(i_gam)*db
                vcinv = 1d0/vc
                vcpow = vc**(2d0/3d0)
                vcm13 = 1d0/(vc**(1d0/3d0))
                do while (ilow < Num_nu .and. V_seed(ilow+1) < 3d-6*vc)
                    ilow = ilow+1
                end do
                if (itail < ilow) itail = ilow
                do while (itail < Num_nu .and. V_seed(itail+1) <= 30d0*vc)
                    itail = itail+1
                end do
                do i_nu = 1, ilow
                    fx = 1.81d0*vp13(i_nu)*vcm13
                    intg(i_nu) = intg(i_nu) + emit_w*fx
                    taug(i_nu) = taug(i_nu) + tau_w*fx
                end do
                do i_nu = ilow+1, itail
                    x = V_seed(i_nu)*vcinv
                    fx = syn_kernel(x,vcpow*vpm23(i_nu),factor)
                    intg(i_nu) = intg(i_nu) + emit_w*fx
                    taug(i_nu) = taug(i_nu) + tau_w*fx
                end do
            end do
            do i_nu = 1, Num_nu
                power = temp*db*intg(i_nu)
                tau = ssa_prefactor*taug(i_nu)*db*Q_weight(i_chi)/(4d0*pi*r2*V_seed(i_nu)*V_seed(i_nu))
                Tau_syn(i_nu,i_chi) = tau
                P_syn(i_nu,i_chi) = power
            end do
        end do
    end if
end subroutine syn_flux_chi

end module rad_common
