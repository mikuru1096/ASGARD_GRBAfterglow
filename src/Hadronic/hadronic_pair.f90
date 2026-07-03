!f2py: skip
module hadronic_pair
    use constants
    use hadronic_base, only: electron_m, check_grid
    implicit none
    private

    real(8), parameter :: reg = 1.0d-50
    real(8), parameter :: tssc = dsqrt(Para_SigmaT)/Para_c
    real(8), parameter :: nssc = Para_SigmaT**(-1.5d0)
    real(8), parameter :: ratessc = nssc/tssc

    public :: pair_operator

contains

! 光子-光子对产生主算子：吸收率、pair 注入谱和能量闭合。
! Main gamma-gamma pair operator: absorption, pair injection, and energy closure.
subroutine pair_operator(nph,eph,phden,ne,ee,ncom,phloss,qpair,qtot,pabs,pinj)
    integer, intent(in) :: nph,ne,ncom
    real(8), intent(in), dimension(nph) :: eph,phden
    real(8), intent(in), dimension(ne) :: ee
    real(8), intent(out), dimension(nph) :: phloss
    real(8), intent(out), dimension(ne) :: qpair,qtot
    real(8), intent(out) :: pabs,pinj
    real(8) :: dlnph,dlne,dln
    real(8), dimension(nph) :: xph,phcgs,phssc,ppng
    real(8), dimension(ne) :: ge,epspair
    real(8), dimension(nph) :: afpair
    integer :: iphmin

    call check_grid(nph,eph,"pair photon grid",dlnph)
    call check_grid(ne,ee,"pair electron grid",dlne)
    dln = dlnph

    if (dabs(dlnph-dlne) > dmax1(1.0d-14,1.0d-10*dabs(dlnph))) then
        error stop "pair production requires photon/electron grids with the same logarithmic spacing."
    end if

    iphmin = grid_offset(ee(1),eph(1),dln)
    if (iphmin < 0) then
        error stop "pair production requires electron grid minimum >= photon grid minimum."
    end if

    xph = eph/electron_m
    ge = ee/electron_m
    phcgs = eph*phden
    phssc = phcgs/nssc

    call photon_loss(nph,xph,dln,phssc,afpair)

    ppng = phssc/(xph*xph)
    call pair_injection(ne,ge,nph,xph,ppng,dln,iphmin,ncom,afpair,phssc,epspair)

    phloss = afpair/tssc
    qpair = (epspair*ratessc)/ee
    qtot = 2d0*qpair
    pabs = sum(xph(iphmin+1:nph)*afpair(iphmin+1:nph)*phssc(iphmin+1:nph))*ratessc*electron_m
    pinj = sum(2d0*ge*epspair)*ratessc*electron_m
end subroutine pair_operator

! 注入谱积分：双重卷积得到 pair source，再按吸收功率重标定。
! Pair-source integral: double convolution followed by absorbed-power closure.
subroutine pair_injection(ne,ge,nph,xph,ppng,dln,iphmin,ncom,afpair,phssc,epspair)
    integer, intent(in) :: ne,nph,iphmin,ncom
    real(8), intent(in), dimension(ne) :: ge
    real(8), intent(in), dimension(nph) :: xph,ppng,afpair,phssc
    real(8), intent(in) :: dln
    real(8), intent(out), dimension(ne) :: epspair
    integer :: i

    epspair = 0d0
    !$OMP PARALLEL DO if((ne-2)*(nph-iphmin) >= 8192) schedule(dynamic,1) private(i)
    do i=2,ne-1
        epspair(i) = pair_bin(ge(i))
    end do
    !$OMP END PARALLEL DO

    call energy_closure
    epspair(1) = 0d0

contains

    ! 单个电子 bin 的 γγ 注入积分。
    ! Gamma-gamma injection integral for a single electron bin.
    real(8) function pair_bin(gm)
        real(8), intent(in) :: gm
        integer :: j,k,outer,inner,kmax
        real(8) :: accum

        outer = outer_pp(gm,dln,iphmin,nph)
        accum = 0d0
        do j=outer+1,nph
            call energy_window(gm,j,inner,kmax)
            do k=inner+1,kmax
                accum = accum + rggd1(gm,xph(j),xph(k))*ppng(k)*ppng(j)
            end do
        end do
        pair_bin = accum*dln*dln*0.75d0*gm
    end function pair_bin

    ! 当前电子能量对应的目标光子积分窗。
    ! Target-photon integration window for the current electron energy.
    subroutine energy_window(gm,j,inner,kmax)
        real(8), intent(in) :: gm
        integer, intent(in) :: j
        integer, intent(out) :: inner,kmax

        inner = min(inner_pp(xph(j),gm,dln,iphmin,nph),nph-1)
        kmax = min(-j + 1 + 2*(iphmin+1) + ncom,nph)
    end subroutine energy_window

    ! 把未归一的 pair source 调到吸收功率的一半；两种电荷相加后守恒。
    ! Renormalize a single-species pair source to half the absorbed power.
    subroutine energy_closure
        real(8) :: ralpha,reps

        ralpha = sum( &
            xph(iphmin+1:nph) * afpair(iphmin+1:nph) * &
            phssc(iphmin+1:nph) &
        )
        reps = sum(ge*epspair)
        if (reps > 1.0d-100) then
            epspair = epspair*(0.5d0*ralpha/reps)
        end if
    end subroutine energy_closure
end subroutine pair_injection

! 光子损失核：对光子场卷积 phibar，得到每个光子能量的吸收率。
! Photon-loss kernel: convolve phibar over the photon field.
subroutine photon_loss(nph,xph,dln,phssc,afpair)
    integer, intent(in) :: nph
    real(8), intent(in), dimension(nph) :: xph,phssc
    real(8), intent(in) :: dln
    real(8), intent(out), dimension(nph) :: afpair
    integer :: i,j
    real(8) :: accum

    !$OMP PARALLEL DO if(nph*nph >= 8192) schedule(static) private(i,j,accum)
    do i=1,nph
        accum = 0d0
        do j=1,nph
            accum = accum + 0.375d0*phibar(xph(i),xph(j))*phssc(j)/(xph(i)*xph(i)*xph(j)*xph(j))
        end do
        afpair(i) = dln*accum
    end do
    !$OMP END PARALLEL DO
end subroutine photon_loss

! φbar(s)：各向同性光子场中 pair-production 截面的角度平均核，s = a*b。
! phibar(s): angle-averaged pair-production cross-section kernel.
real(8) function phibar(a,b)
    real(8), intent(in) :: a,b
    real(8) :: s,w
    s = a*b
    if (s <= 1d0) then
        phibar = 0d0
    else if (s <= 1.1d0) then
        phibar = phibar1(s)
    else if (s < 5.0d0) then
        w = -1d0 + 2d0*s*(1d0 + dsqrt(1d0 - 1d0/s))
        phibar = phibar2(s,w)
    else
        phibar = phibar3(s)
    end if
end function phibar

! 近阈区：Taylor 展开至 (s-1)^(7/2)。
! Near threshold: Taylor expansion through (s-1)^(7/2).
real(8) function phibar1(s)
    real(8), intent(in) :: s
    real(8) :: s1,s2
    s1 = s - 1d0
    s2 = s1*dsqrt(s1)
    phibar1 = s2*(1.333333d0 + 1.2d0*s1 - 253.0d0*s1*s1/70.0d0)
end function phibar1

! 中间区：使用 w 参数化的解析表达式。
! Intermediate branch: analytic w-parametrized expression.
real(8) function phibar2(s,w)
    real(8), intent(in) :: s,w
    real(8) :: v,u
    v = dlog(w)
    u = 1d0 - 1d0/s
    phibar2 = ( &
        (2d0 - 4.0d0*s)*dsqrt(u) + v*(4.0d0*dlog(1d0+w) - 3.0d0*v + s*(1d0 + u*u)) - &
        3.289868d0 + phisum(w) &
    )
end function phibar2

! 高能区：使用 s >= 5 的渐近展开。
! High-energy branch: asymptotic expansion for s >= 5.
real(8) function phibar3(s)
    real(8), intent(in) :: s
    real(8) :: w
    w = dlog(4.0d0*s)
    phibar3 = (2d0*s + w)*(w - 2d0) + (w + 1.125d0)/s - 0.289868d0
end function phibar3

! phibar2 的辅助级数 Σ(-1)^i/(i^2 w^i)。
! Auxiliary series used by phibar2.
real(8) function phisum(w)
    real(8), intent(in) :: w
    integer :: i
    phisum = 0d0
    do i=1,14
        phisum = phisum + dble(sign_int(i))/(w**i*dble(i*i))
    end do
    phisum = -4.0d0*phisum
end function phisum

! 内层积分下界：给定目标光子 x 和 pair 电子 gamma，求伴随光子起点。
! Inner integration bound: lower target-photon index for the given x and gamma.
integer function inner_pp(x,gm,dln,iphmin,nph)
    real(8), intent(in) :: x,gm,dln
    integer, intent(in) :: iphmin,nph
    real(8) :: fval,arg

    if (x <= 0.5d0) then
        fval = fppm(x,gm)
        if (dabs(fval) <= 1.0d-300) then
            inner_pp = 0
            return
        end if
        arg = gm - x + 0.5d0*(fval + 1d0/fval)
        if (arg <= 0d0) then
            inner_pp = 0
            return
        end if
        inner_pp = min(max(int(dlog(arg)/dln + dble(iphmin) + 1d0),0),nph)
        return
    end if
    if (x < 1d0 .and. gm < gmb(x)) then
        fval = fppm(x,gm)
        if (dabs(fval) <= 1.0d-300) then
            inner_pp = 0
            return
        end if
        arg = gm - x + 0.5d0*(fval + 1d0/fval)
        if (arg <= 0d0) then
            inner_pp = 0
            return
        end if
        inner_pp = min(max(int(dlog(arg)/dln + dble(iphmin) + 1d0),0),nph)
        return
    end if
    if (x > 1d0 .and. gm < gmb(x)) then
        fval = fppp(x,gm)
        if (dabs(fval) <= 1.0d-300) then
            inner_pp = 0
            return
        end if
        arg = gm - x + 0.5d0*(fval + 1d0/fval)
        if (arg <= 0d0) then
            inner_pp = 0
            return
        end if
        inner_pp = min(max(int(dlog(arg)/dln + dble(iphmin) + 1d0),0),nph)
        return
    end if
    inner_pp = 0
end function inner_pp

! 外层积分下界：由 pair 电子 gamma 决定第一个可产生 pair 的光子 bin。
! Outer integration bound: first photon bin kinematically allowed for this gamma.
integer function outer_pp(gm,dln,iphmin,nph)
    real(8), intent(in) :: gm,dln
    integer, intent(in) :: iphmin,nph
    outer_pp = min(max(int(dlog(xlow(gm))/dln + dble(iphmin)),0),nph-1)
end function outer_pp

! 最小光子能量 x_- = (gamma - sqrt(gamma^2-1))/2。
! Minimum photon energy on the pair-production kinematic boundary.
real(8) function xlow(gm)
    real(8), intent(in) :: gm
    xlow = 0.5d0*gm*(1d0 - betag(gm))
end function xlow

! 分支边界 gamma_b(x)：区分 f_- 和 f_+ 控制的积分域。
! Branch boundary gamma_b(x): selects the f_minus or f_plus integration branch.
real(8) function gmb(x)
    real(8), intent(in) :: x
    gmb = x - (x - 1d0)/(2d0*x - 1d0)
end function gmb

! f_- 边界函数：用于低能侧 x1 积分下界。
! f_minus boundary function for the low-energy side of the x1 integral.
real(8) function fppm(x,gp)
    real(8), intent(in) :: x,gp
    fppm = 2d0*x - gp*(1d0 - betag(gp))
end function fppm

! f_+ 边界函数：用于高能侧 x1 积分下界。
! f_plus boundary function for the high-energy side of the x1 integral.
real(8) function fppp(x,gp)
    real(8), intent(in) :: x,gp
    fppp = 2d0*x - gp*(1d0 + betag(gp))
end function fppp

! 微分截面核 R_gg(gamma, x, x1)：Aharonian+1983 公式和运动学截断。
! Differential kernel R_gg(gamma, x, x1): Aharonian+1983 with kinematic cuts.
real(8) function rggd1(gm,x,x1)
    real(8), intent(in) :: gm,x,x1
    real(8) :: gp,kl,ku,tval,xup,xlo,u,v,hnup,hnlo,hiup,hilo
    real(8) :: td2a,td2b,td2c,td2d,denom

    gp = x + x1 - gm
    if (gp <= 1d0 .or. gm <= 1d0) then
        rggd1 = 0d0
        return
    end if

    if (gp > 1d1 .and. gm > 1d1) then
        kl = 0.25d0/gp + 0.09375d0/gp**3 + 0.25d0/gm + 0.09375d0/gm**3
        ku = 0.5d0*(gp + dsqrt(gp*gp - 1d0) + gm + dsqrt(gm*gm - 1d0))
    else if (gp > 1d1) then
        kl = 0.5d0*(0.5d0/gp + gm - dsqrt(gm*gm - 1d0))
        ku = gp + 0.5d0*(gm + dsqrt(gm*gm - 1d0))
    else if (gm > 1d1) then
        kl = 0.5d0*(gp - dsqrt(gp*gp - 1d0) + 0.5d0/gm)
        ku = gm + 0.5d0*(gp + dsqrt(gp*gp - 1d0))
    else
        kl = 0.5d0*(gp + gm - dsqrt(gp*gp - 1d0) - dsqrt(gm*gm - 1d0))
        ku = 0.5d0*(gp + gm + dsqrt(gp*gp - 1d0) + dsqrt(gm*gm - 1d0))
    end if

    if (x < kl .or. x1 < kl .or. x > ku .or. x1 > ku) then
        rggd1 = 0d0
        return
    end if

    tval = x + x1
    xup = xcmu(gm,x,x1)
    xlo = xcml(gm,x,x1)
    u = 4.0d0*xup*xup/(tval*tval)
    v = 4.0d0*xlo*xlo/(tval*tval)

    hnup = ((gm - x)*(gm - x) - 1d0)*xup*xup/(x*x1)
    hnlo = ((gm - x)*(gm - x) - 1d0)*xlo*xlo/(x*x1)
    hiup = ((gm - x1)*(gm - x1) - 1d0)*xup*xup/(x*x1)
    hilo = ((gm - x1)*(gm - x1) - 1d0)*xlo*xlo/(x*x1)

    td2a = td2(gm,x,x1,xup)
    td2b = td2(gm,x,x1,xlo)
    td2c = td2(gm,x1,x,xup)
    td2d = td2(gm,x1,x,xlo)

    if (hnup > 1000.0d0 .and. hnlo > 1000.0d0) then
        denom = dsqrt((gm - x)*(gm - x) - 1d0)
        td2a = -(xup*xup + 2d0*dlog(xup))/denom
        td2b = -(xlo*xlo + 2d0*dlog(xlo))/denom
    end if
    if (hiup > 1000.0d0 .and. hilo > 1000.0d0) then
        denom = dsqrt((gm - x1)*(gm - x1) - 1d0)
        td2c = -(xup*xup + 2d0*dlog(xup))/denom
        td2d = -(xlo*xlo + 2d0*dlog(xlo))/denom
    end if

    if (u < 0.01d0 .and. v < 0.01d0) then
        rggd1 = -0.25d0*(0.5d0*tval*(u - v) + td2a - td2b + td2c - td2d)
    else
        rggd1 = -0.25d0*(tval*(dsqrt(1d0 - v) - dsqrt(1d0 - u + 1.0d-9)) + td2a - td2b + td2c - td2d)
    end if
end function rggd1

! 质心系能量上界：两个不变量给出的较小值。
! Upper center-of-momentum energy from the smaller invariant bound.
real(8) function xcmu(gm,x,x1)
    real(8), intent(in) :: gm,x,x1
    real(8) :: gp
    gp = x + x1 - gm
    if (gp > 15.0d0 .and. gm > 15.0d0) then
        xcmu = min(dsqrt(x*x1),dsqrt(gp*gm))
    else
        xcmu = min(dsqrt(x*x1),dsqrt(0.5d0*(gm*gp + 1d0 + dsqrt(gm*gm - 1d0)*dsqrt(gp*gp - 1d0))))
    end if
end function xcmu

! 质心系能量下界：按 gamma 大小选择数值稳定的渐近分支。
! Lower center-of-momentum energy with gamma-dependent asymptotic branches.
real(8) function xcml(gm,x,x1)
    real(8), intent(in) :: gm,x,x1
    real(8) :: gp
    gp = x + x1 - gm
    if (gp > 1d1) then
        if (gm > 1d1) then
            xcml = 0.5d0*dsqrt(2d0 + gp/gm + gm/gp)
            return
        end if
        if (gm < 1.001d0) then
            xcml = 1.001d0
            return
        end if
        xcml = dsqrt(0.5d0*((gm - dsqrt(gm*gm - 1d0 + reg))*(gp - 0.5d0/gp) + 1d0))
        return
    end if
    if (gp < 1.001d0) then
        xcml = 1.001d0
        return
    end if
    if (gm > 1d1) then
        xcml = dsqrt(0.5d0*((gp - dsqrt(gp*gp - 1d0 + reg))*(gm - 0.5d0/gm) + 1d0))
        return
    end if
    if (gm < 1.001d0) then
        xcml = 1.001d0
        return
    end if
    xcml = dsqrt(0.5d0*(gp*gm + 1d0 - dsqrt((gp*gp - 1d0)*(gm*gm - 1d0))))
end function xcml

! 微分截面中的 T(gamma,x,x1,xcm) 辅助函数。
! Auxiliary T(gamma,x,x1,xcm) term in the differential kernel.
real(8) function td2(gm,x,x1,xcm)
    real(8), intent(in) :: gm,x,x1,xcm
    real(8) :: y,y2,q,h,ah
    y = dsqrt(x*x1)
    y2 = y*y
    q = xcm/y
    h = ((gm - x)*(gm - x) - 1d0)*q*q
    if (h > 1000.0d0) then
        td2 = -(xcm*xcm + h/(xcm*xcm) + 0.886294d0 - 0.5d0*(gm*(x1 - x) + x*x)/y2 + dlog(h))*q/dsqrt(h)
        return
    end if
    if (h < -0.999d0) then
        td2 = q*(q*q*(y2 - 1d0)*(-1.5707963d0) + 0.5d0*((x1 - x)/y2 - 6.28319d0))
        return
    end if
    ah = dsqrt(1d0 + h)
    td2 = q*q*q*(y2 - 1d0)*a0ratio(h) - ah/(xcm*y) + &
                   0.5d0*q*((1d0 + (gm*(x1 - x) + x*x)/y2 - 2d0*q*q)/ah - 4.0d0*a0(h))
end function td2

! (a0(h) - sqrt(1+h))/h：大 h 使用渐近式。
! Stable ratio (a0(h) - sqrt(1+h))/h; large h uses asymptotics.
real(8) function a0ratio(h)
    real(8), intent(in) :: h
    if (h <= 20.0d0) then
        a0ratio = a0ratiof(h)
    else
        a0ratio = -(1d0 - (0.693147d0 + 0.5d0*dlog(h))/h)/dsqrt(h)
    end if
end function a0ratio

! a0(h)：h>0 用 asinh，h<0 用 asin。
! a0(h): asinh branch for h>0, asin branch for h<0.
real(8) function a0(h)
    real(8), intent(in) :: h
    if (h < 20.0d0) then
        a0 = a0f(h)
    else
        a0 = (0.693147d0 + 0.5d0*dlog(h))/dsqrt(h)
    end if
end function a0

! a0 基础实现：小 |h| 用 Taylor 展开。
! Base a0 implementation; small |h| uses a Taylor expansion.
real(8) function a0f(h)
    real(8), intent(in) :: h
    if (h > 1.0d-6) then
        a0f = dlog(dsqrt(h) + dsqrt(1d0 + h))/dsqrt(h)
    else if (h < -1.0d-6) then
        a0f = dasin(dsqrt(-h) - 1.0d-8)/dsqrt(-h)
    else
        a0f = 1d0 - h/6.0d0
    end if
end function a0f

! a0ratio 基础实现：小 |h| 用 Taylor 展开。
! Base a0ratio implementation; small |h| uses a Taylor expansion.
real(8) function a0ratiof(h)
    real(8), intent(in) :: h
    if (h >= -0.1d0 .and. h <= 0.1d0) then
        a0ratiof = -0.666667d0 + (0.2d0 - (0.107143d0 - 0.0694444d0*h)*h)*h
    else
        a0ratiof = (a0f(h) - dsqrt(1d0 + h + 1.0d-8))/h
    end if
end function a0ratiof

! 相对论 beta(gamma)：低 gamma 直接算，高 gamma 用渐近展开。
! Relativistic beta(gamma): exact low-gamma branch and high-gamma asymptotics.
real(8) function betag(gamma)
    real(8), intent(in) :: gamma
    if (gamma < 1d1) then
        betag = dsqrt(1d0 - 1d0/(gamma*gamma) + reg)
    else
        betag = 1d0 - 0.5d0/(gamma*gamma)*(1d0 + 0.25d0/(gamma*gamma))
    end if
end function betag

! 级数符号：偶数为 +1，奇数为 -1。
! Alternating-series sign: even is +1 and odd is -1.
integer function sign_int(i)
    integer, intent(in) :: i
    if (mod(i,2) == 0) then
        sign_int = 1
    else
        sign_int = -1
    end if
end function sign_int

! 两套对数网格最小能量的整数偏移。
! Integer offset between the logarithmic grid minima.
integer function grid_offset(emin,phmin,dln)
    real(8), intent(in) :: emin,phmin,dln
    real(8) :: ratio
    integer :: rounded
    ratio = dlog(emin/phmin)/dln
    rounded = nint(ratio)
    if (dabs(ratio - dble(rounded)) > 1.0d-9) then
        error stop "pair production grid minima are not aligned on the shared logarithmic lattice."
    end if
    grid_offset = rounded
end function grid_offset

end module hadronic_pair
