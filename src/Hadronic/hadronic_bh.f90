!f2py: skip
module hadronic_bh
    use constants
    use hadronic_base, only: proton_m, electron_m, check_grid
    implicit none
    private

    real(8), parameter :: alpha = 1d0/137d0
    integer, parameter :: nout = 5
    integer, parameter :: nin = 5

    public :: bh_calc

contains

! Bethe-Heitler 主算子：计算 pair 注入、质子损失和光子 sink。
! Main Bethe-Heitler operator for pair injection, proton loss, and photon sink.
subroutine bh_calc(np,ep,pden,nph,eph,phden,ne,ee,qpair,ploss,phloss)
    integer, intent(in) :: np,nph,ne
    real(8), intent(in), dimension(np) :: ep,pden
    real(8), intent(in), dimension(nph) :: eph,phden
    real(8), intent(in), dimension(ne) :: ee
    real(8), intent(out), dimension(ne) :: qpair
    real(8), intent(out), dimension(np) :: ploss
    real(8), intent(out), dimension(nph) :: phloss
    integer :: ip,iph,ie
    real(8) :: dep,dph,de,kinj,kloss,kern
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
    kloss = (3d0/(8d0*pi))*alpha*Para_SigmaT*Para_c*(electron_m/proton_m)

    qlog = 0d0
    ploss = 0d0
    phloss = 0d0
    lphlog = 0d0

    do ip=1,np
        do iph=1,nph
            ploss(ip) = ploss(ip) + loss_point(gp(ip),xph(iph),phlog(iph))
            if (2d0*gp(ip)*xph(iph) <= 2d0) cycle

            ! 同一 kernel 同时累加 pair source 和 photon sink。
            ! The same kernel contributes to the pair source and photon sink.
            do ie=1,ne
                kern = bh_pair(xe(ie),gp(ip),xph(iph))
                qlog(ie) = qlog(ie) + kern*plog(ip)*phlog(iph)
                lphlog(iph) = lphlog(iph) + kern*plog(ip)*phlog(iph)
            end do
        end do
    end do

    qlog = kinj*dep*dph*qlog
    qpair = qlog/ee
    lphlog = kinj*dep*de*lphlog
    where (phlog > 0d0)
        phloss = lphlog/phlog
    elsewhere
        phloss = 0d0
    end where

contains

    real(8) function loss_point(gam,ephoton,phlog1)
        real(8), intent(in) :: gam,ephoton,phlog1

        loss_point = -kloss*dph*proton_loss(gam,ephoton)*phlog1
    end function loss_point
end subroutine bh_calc

! Bethe-Heitler pair 产生核：给定 p/ph 能量，返回 e± 微分产额。
! Bethe-Heitler pair kernel for the differential e+/e- yield.
real(8) function bh_pair(ee,gam,ephoton)
    real(8), intent(in) :: ee,gam,ephoton
    real(8) :: upper,lower,diff

    upper = 2d0*gam*ephoton
    lower = ((gam + ee)*(gam + ee))/(2d0*gam*ee)
    diff = (upper - lower)/lower
    if (diff <= 1d-8) then
        bh_pair = 0d0
        return
    end if

    bh_pair = bh_rk3(bh_outer,lower,upper,gam,ee,nout)*ee/(2d0*gam*gam*gam*ephoton*ephoton)
end function bh_pair

! Bethe-Heitler 外层积分核：对 omega 积分。
! Outer Bethe-Heitler integral over omega.
real(8) function bh_outer(gam,ee,omega)
    real(8), intent(in) :: gam,ee,omega
    real(8), parameter :: eps = 1d-5
    real(8) :: lower,upper,diff

    lower = (gam*gam + ee*ee)/(2d0*gam*ee) + eps
    upper = omega - 1d0 - eps
    diff = (upper - lower)/lower
    if (diff <= -1d-4) then
        error stop "Bethe-Heitler outer integral received invalid limits."
    end if
    if (diff < 1d-4) then
        bh_outer = 0d0
        return
    end if

    bh_outer = bh_rk4(bh_inner,lower,upper,gam,ee,omega,nin)
end function bh_outer

! Bethe-Heitler 内层积分核：对 ebar 积分。
! Inner Bethe-Heitler integral over ebar.
real(8) function bh_inner(gam,ee,omega,ebar)
    real(8), intent(in) :: gam,ee,omega,ebar
    real(8), parameter :: eps = 1d-20
    real(8) :: pbar,ksi

    pbar = dsqrt(ebar*ebar - 1d0 + eps)
    ksi = (gam*ebar - ee)/(gam*pbar)
    if (ksi < -1d0 .or. ksi > 1d0 .or. ebar < 1d0 .or. ebar > omega - 1d0) then
        bh_inner = 0d0
        return
    end if

    bh_inner = omega*bh_sigma(omega,ebar,ksi)/pbar
end function bh_inner

! Bethe-Heitler 微分截面 sigma(omega, ebar, ksi)，来自 Blumenthal 1970。
! Analytic Bethe-Heitler differential cross-section from Blumenthal 1970.
real(8) function bh_sigma(omega,ebar,ksi)
    real(8), intent(in) :: omega,ebar,ksi
    real(8), parameter :: eps = 1d-30
    real(8) :: k,g,u,p,g1,p1,dloc,tloc,d1t,y1,yloc
    real(8) :: a0,a1,a2,a3,a4,a50,a51,a52,a53,a6,a7

    if (omega < 2.0001d0 .or. omega > 600d0) then
        bh_sigma = 0d0
        return
    end if

    k = omega
    g = ebar
    u = ksi

    p = dsqrt(g*g - 1d0 + eps)
    g1 = k - g
    p1 = dsqrt(g1*g1 - 1d0 + eps)
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

! Bethe-Heitler 质子能量损失核。
! Proton energy-loss kernel for Bethe-Heitler interactions.
real(8) function proton_loss(gam,ephoton)
    real(8), intent(in) :: gam,ephoton

    proton_loss = 0.5d0*bh_phi(2d0*gam*ephoton)/(gam*gam*ephoton*ephoton)
end function proton_loss

! Bethe-Heitler 能量损失核 phi(x)：低能多项式，高能对数近似。
! Piecewise approximation for the Bethe-Heitler energy-loss kernel phi(x).
real(8) function bh_phi(xval)
    real(8), intent(in) :: xval
    real(8) :: zval,zlog

    if (xval < 2d0 + 1d-6) then
        bh_phi = 0d0
        return
    end if
    if (xval < 25d0) then
        zval = xval - 2d0
        bh_phi = (pi/12d0)*zval**4 / &
                 (1d0 + 0.8048d0*zval + 0.1459d0*zval*zval + 1.137d-3*zval**3 - 3.879d-6*zval**4)
        return
    end if

    zlog = dlog(xval)
    bh_phi = xval*(-86.07d0 + 50.96d0*zlog - 14.45d0*zlog*zlog + (8d0/3d0)*zlog**3) / &
             (1d0 - 2.910d0/xval - 78.35d0/(xval*xval) - 1837d0/(xval**3))
end function bh_phi

! 三参数 RK 3/8 积分器；用于 BH 外层核。
! Three-argument RK 3/8 integrator used by the outer BH kernel.
real(8) function bh_rk3(func,lower,upper,var0,var1,nbin)
    interface
        function func(arg0,arg1,arg2) result(value)
            real(8), intent(in) :: arg0,arg1,arg2
            real(8) :: value
        end function func
    end interface
    integer, intent(in) :: nbin
    real(8), intent(in) :: lower,upper,var0,var1
    integer :: i
    real(8) :: dx,xval,k1,k2,k4

    dx = (upper - lower)/dble(nbin)
    xval = lower
    bh_rk3 = 0d0
    do i=1,nbin
        k1 = dx*func(var0,var1,xval)
        k2 = dx*func(var0,var1,xval + 0.5d0*dx)
        k4 = dx*func(var0,var1,xval + dx)
        bh_rk3 = bh_rk3 + (k1 + k4)/6d0 + (2d0/3d0)*k2
        xval = xval + dx
    end do
end function bh_rk3

! 四参数 RK 3/8 积分器；用于 BH 内层核。
! Four-argument RK 3/8 integrator used by the inner BH kernel.
real(8) function bh_rk4(func,lower,upper,var0,var1,var2,nbin)
    interface
        function func(arg0,arg1,arg2,arg3) result(value)
            real(8), intent(in) :: arg0,arg1,arg2,arg3
            real(8) :: value
        end function func
    end interface
    integer, intent(in) :: nbin
    real(8), intent(in) :: lower,upper,var0,var1,var2
    integer :: i
    real(8) :: dx,xval,k1,k2,k4

    dx = (upper - lower)/dble(nbin)
    xval = lower
    bh_rk4 = 0d0
    do i=1,nbin
        k1 = dx*func(var0,var1,var2,xval)
        k2 = dx*func(var0,var1,var2,xval + 0.5d0*dx)
        k4 = dx*func(var0,var1,var2,xval + dx)
        bh_rk4 = bh_rk4 + (k1 + k4)/6d0 + (2d0/3d0)*k2
        xval = xval + dx
    end do
end function bh_rk4

end module hadronic_bh
