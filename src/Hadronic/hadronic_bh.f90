!f2py: skip
module hadronic_bh
    use constants
    use hadronic_base, only: proton_m, electron_m, check_grid
    implicit none
    private

    real(8), parameter :: alpha = 1d0/137d0
    real(8), parameter :: omega_max = 6d2
    real(8), parameter :: qnode(6) = [ &
        0.12523340851146891547d0, 0.36783149899818019375d0, &
        0.58731795428661744730d0, 0.76990267419430468704d0, &
        0.90411725637047485668d0, 0.98156063424671925069d0]
    real(8), parameter :: qweight(6) = [ &
        0.24914704581340278500d0, 0.23349253653835480876d0, &
        0.20316742672306592175d0, 0.16007832854334622633d0, &
        0.10693932599531843096d0, 0.04717533638651182719d0]

    public :: bh_calc, bh_pair_edge

contains

! Bethe-Heitler 主算子：计算 pair 注入、质子损失和光子 sink。
! Main Bethe-Heitler operator: one-charge pair injection, fractional proton loss,
! and photon sink. The transport owner converts the fractional loss to dgamma/dt.
subroutine bh_calc(np,ep,pden,nph,eph,phden,ne,ee,qpair,ploss,phloss)
    integer, intent(in) :: np,nph,ne
    real(8), intent(in), dimension(np) :: ep,pden
    real(8), intent(in), dimension(nph) :: eph,phden
    real(8), intent(in), dimension(ne) :: ee
    real(8), intent(out), dimension(ne) :: qpair
    real(8), intent(out), dimension(np) :: ploss
    real(8), intent(out), dimension(nph) :: phloss
    integer :: ip,iph,ie
    real(8) :: dep,dph,de,kinj,kern,emoment
    real(8), dimension(np) :: plog,gp
    real(8), dimension(nph) :: phlog,lphlog,xph
    real(8), dimension(ne) :: qlog,xe

    call check_grid(np,ep,"ep",dep)
    call check_grid(nph,eph,"eph",dph)
    call check_grid(ne,ee,"ee",de)

    plog = ep*pden
    phlog = eph*phden
    gp = ep/proton_m
    xph = eph/electron_m
    xe = ee/electron_m

    kinj = Para_c*3d0*Para_SigmaT*alpha/(16d0*pi)

    qlog = 0d0
    phloss = 0d0
    lphlog = 0d0

    do ip=1,np
        emoment = 0d0
        do iph=1,nph
            if (2d0*gp(ip)*xph(iph) <= 2d0) cycle

            ! 同一 kernel 同时累加 pair source 和 photon sink。
            ! The same kernel contributes to the pair source and photon sink.
            do ie=1,ne
                kern = bh_pair(xe(ie),gp(ip),xph(iph))
                qlog(ie) = qlog(ie) + kern*plog(ip)*phlog(iph)
                lphlog(iph) = lphlog(iph) + kern*plog(ip)*phlog(iph)
                emoment = emoment + kern*ee(ie)*phlog(iph)
            end do
        end do
        ! Fractional proton loss is the two-charge energy moment of this same
        ! differential pair kernel, including its omega <= 600 support.
        ploss(ip) = -2d0*kinj*dph*de*emoment/ep(ip)
    end do
    qlog = kinj*dep*dph*qlog
    qpair = qlog/ee
    lphlog = kinj*dep*de*lphlog
    where (phlog > 0d0)
        phloss = lphlog/phlog
    elsewhere
        phloss = 0d0
    end where
end subroutine bh_calc

! Bethe-Heitler pair 产生核：给定 p/ph 能量，返回 e± 微分产额。
! Bethe-Heitler pair kernel for the differential e+/e- yield.
real(8) function bh_pair(ee,gam,ephoton)
    real(8), intent(in) :: ee,gam,ephoton
    real(8) :: upper,lower

    upper = min(2d0*gam*ephoton,omega_max)
    lower = ((gam + ee)*(gam + ee))/(2d0*gam*ee)
    if (upper <= lower) then
        bh_pair = 0d0
        return
    end if

    bh_pair = bh_quad3(bh_outer,lower,upper,gam,ee)*ee/(2d0*gam*gam*gam*ephoton*ephoton)
end function bh_pair

! BH 运动学电子上边界：由 omega 上限的正根给出最大 gamma_e/gamma_p。
! BH kinematic electron upper edge from the positive gamma_e/gamma_p root at the omega cap.
pure real(8) function bh_pair_edge(gpmax) result(edge)
    real(8), intent(in) :: gpmax

    edge=gpmax*(omega_max-1d0+dsqrt(omega_max*(omega_max-2d0)))
end function bh_pair_edge

! Bethe-Heitler 外层积分核：对 omega 积分。
! Outer Bethe-Heitler integral over omega.
real(8) function bh_outer(gam,ee,omega)
    real(8), intent(in) :: gam,ee,omega
    real(8) :: lower,upper

    lower = (gam*gam + ee*ee)/(2d0*gam*ee)
    upper = omega - 1d0

    bh_outer = bh_quad4(bh_inner,lower,upper,gam,ee,omega)
end function bh_outer

! Bethe-Heitler 内层积分核：对 ebar 积分。
! Inner Bethe-Heitler integral over ebar.
real(8) function bh_inner(gam,ee,omega,ebar)
    real(8), intent(in) :: gam,ee,omega,ebar
    real(8) :: pbar,ksi

    pbar = dsqrt(ebar*ebar - 1d0)
    ksi = (gam*ebar - ee)/(gam*pbar)
    bh_inner = omega*bh_sigma(omega,ebar,ksi)/pbar
end function bh_inner

! Bethe-Heitler 微分截面 sigma(omega, ebar, ksi)，来自 Blumenthal 1970。
! Analytic Bethe-Heitler differential cross-section from Blumenthal 1970.
real(8) function bh_sigma(omega,ebar,ksi)
    real(8), intent(in) :: omega,ebar,ksi
    real(8) :: k,g,u,p,g1,p1,dloc,tloc,d1t,y1,yloc
    real(8) :: a0,a1,a2,a3,a4,a50,a51,a52,a53,a6,a7

    k = omega
    g = ebar
    u = ksi

    p = dsqrt(g*g - 1d0)
    g1 = k - g
    p1 = dsqrt(g1*g1 - 1d0)
    dloc = g - u*p
    tloc = dsqrt(k*k + p*p - 2d0*k*p*u)
    d1t = dlog((tloc + p1)/(tloc - p1))
    y1 = dlog((g1 + p1)/(g1 - p1))/p1
    yloc = 2d0*dlog((g*g1 + p*p1 + 1d0)/k)/(p*p)

    a0 = p*p1/(k**3)
    a1 = -4d0*(1d0 - u*u)*(2d0*g*g + 1d0)/(p*p*dloc**4)
    a2 = (5d0*g*g - 2d0*g*g1 + 3d0)/(p*p*dloc*dloc)
    a3 = (p*p - k*k)/(tloc*tloc*dloc*dloc)
    a4 = 2d0*g1/(p*p*dloc)
    a50 = yloc/(p*p1)
    a51 = 2d0*g*(1d0 - u*u)*(3d0*k + p*p*g1)/(dloc**4)
    a52 = (2d0*g*g*(g*g + g1*g1) - 7d0*g*g - 3d0*g*g1 - g1*g1 + 1d0)/(dloc*dloc)
    a53 = k*(g*g - g*g1 - 1d0)/dloc
    a6 = -(d1t/(p1*tloc))*(2d0/(dloc*dloc) - 3d0*k/dloc - k*(p*p - k*k)/(tloc*tloc*dloc))
    a7 = -2d0*y1/dloc
    bh_sigma = a0*(a1 + a2 + a3 + a4 + a50*(a51 + a52 + a53) + a6 + a7)
end function bh_sigma

! BH 外层核使用的十二点 Gauss--Legendre 积分。
! Twelve-point Gauss-Legendre quadrature used by the outer BH kernel.
real(8) function bh_quad3(func,lower,upper,var0,var1)
    interface
        function func(arg0,arg1,arg2) result(value)
            real(8), intent(in) :: arg0,arg1,arg2
            real(8) :: value
        end function func
    end interface
    real(8), intent(in) :: lower,upper,var0,var1
    integer :: i
    real(8) :: half,mid,delta

    half = 0.5d0*(upper - lower)
    mid = 0.5d0*(upper + lower)
    bh_quad3 = 0d0
    do i=1,6
        delta = half*qnode(i)
        bh_quad3 = bh_quad3 + qweight(i)*( &
            func(var0,var1,mid-delta) + func(var0,var1,mid+delta))
    end do
    bh_quad3 = half*bh_quad3
end function bh_quad3

! BH 内层核使用的十二点 Gauss--Legendre 积分。
! Twelve-point Gauss-Legendre quadrature used by the inner BH kernel.
real(8) function bh_quad4(func,lower,upper,var0,var1,var2)
    interface
        function func(arg0,arg1,arg2,arg3) result(value)
            real(8), intent(in) :: arg0,arg1,arg2,arg3
            real(8) :: value
        end function func
    end interface
    real(8), intent(in) :: lower,upper,var0,var1,var2
    integer :: i
    real(8) :: half,mid,delta

    half = 0.5d0*(upper - lower)
    mid = 0.5d0*(upper + lower)
    bh_quad4 = 0d0
    do i=1,6
        delta = half*qnode(i)
        bh_quad4 = bh_quad4 + qweight(i)*( &
            func(var0,var1,var2,mid-delta) + func(var0,var1,var2,mid+delta))
    end do
    bh_quad4 = half*bh_quad4
end function bh_quad4

end module hadronic_bh
