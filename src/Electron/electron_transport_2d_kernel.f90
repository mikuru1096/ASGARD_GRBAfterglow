!f2py: skip
module electron_transport_2d
  use constants
  use electron_transport_common, only: prepare_implicit_coeffs, backward_sweep, &
                             prepare_remap, ppm_eval_cell, &
                             char_update, cooling_affine, cooling_piecewise
  use electron_injection_profiles, only: log_edges
  implicit none
  private
  real(8), parameter :: shock_comp = 4d0

  ! 有限 q 质量坐标的 2D 电子输运核：几何、q 输运、能量输运分开暴露。
  ! Finite-q mass-coordinate 2D electron transport kernel: geometry, q transport, and energy transport are separate.
    public :: q_geometry, q_cell_geometry, shock_state, downstream_grid, q_divergence, q_step_limit
    public :: advance_q_implicit, advance_energy_chi, advance_pwncr_q, advance_pwncr_energy
    public :: advance_stoch_chi, advance_q_charint, advance_q_diffusion, advance_energy_charint

contains

! 构建有限主动壳层质量坐标q：q=0为激波前沿，外边界停在强激波压缩给出的接触面量级。
! Build the finite active shell mass coordinate: q=0 at the shock front, outer face near the contact scale.
subroutine q_geometry(nchi, dq, q_face, q_grid)
    integer, intent(in) :: nchi
    real(8), intent(out), dimension(0:nchi) :: q_face
    real(8), intent(out), dimension(nchi) :: q_grid
    real(8), intent(out) :: dq
    real(8) :: qactive
    integer :: ichi

    qactive = 1d0-(1d0-1d0/shock_comp)**shock_comp
    dq = qactive/dble(nchi)
    do ichi = 0, nchi
        q_face(ichi) = dble(ichi)*dq
    end do
    do ichi = 1, nchi
        q_grid(ichi) = (dble(ichi)-0.5d0)*dq
    end do
end subroutine q_geometry

! 单个 q 点的桥接几何：BM 超相对论极限和平面强激波极限按四速度权重混合。
! Geometry at a single q point: bridge the BM ultra-relativistic limit and planar strong-shock limit by four-velocity weight.
subroutine q_point(kmed,rshell,gf,qloc,rcell,gcell,bcell)
    integer, intent(in) :: kmed
    real(8), intent(in) :: rshell,gf,qloc
    real(8), intent(out) :: rcell,gcell,bcell
    real(8) :: qtail, alpha, uf, bf, w, chibm, logxbm, logxst, xrat, ubm, bst, ust, ycell, ucell

    qtail = 1d0-qloc
    alpha = dble(4-kmed)/dble(3-kmed)
    uf = dsqrt(gf*gf-1d0)
    bf = uf/gf
    w = uf*uf/(1d0+uf*uf)
    chibm = qtail**(-alpha)
    logxbm = -(chibm-1d0)/(4d0*dble(4-kmed)*gf*gf)
    logxst = dlog(qtail)/(shock_comp*dble(3-kmed))
    xrat = dexp(w*logxbm+(1d0-w)*logxst)
    ubm = uf/dsqrt(chibm)
    bst = bf*dexp(logxst)
    ust = bst/dsqrt(1d0-bst*bst)
    ycell = w*dasinh(ubm)+(1d0-w)*dasinh(ust)
    ucell = dsinh(ycell)
    gcell = dsqrt(1d0+ucell*ucell)
    bcell = ucell/gcell
    rcell = rshell*xrat
end subroutine q_point

! 计算 q 单元中心的半径、Lorentz 因子、速度以及相对激波速度。
! Compute radius, Lorentz factor, velocity, and shock-relative velocity at q-cell centers.
subroutine q_cell_geometry(nchi,kmed,rshell,gf,q_grid,rcell,gcell,bcell,brelsh)
    integer, intent(in) :: nchi,kmed
    real(8), intent(in), dimension(nchi) :: q_grid
    real(8), intent(in) :: rshell,gf
    real(8), intent(out), dimension(nchi) :: rcell,gcell,bcell,brelsh
    real(8) :: bshock, uf, w, ushbm, ushst, ysh, ush
    integer :: ichi

    ! 激波速度使用同一个 BM/ST 桥接式。
    ! Shock speed uses the same BM/ST bridge.
    uf = dsqrt(gf*gf-1d0)
    w = uf*uf/(1d0+uf*uf)
    ushbm = dsqrt(2d0*gf*gf-1d0)
    ushst = shock_comp*uf/(shock_comp-1d0)
    ysh = w*dasinh(ushbm)+(1d0-w)*dasinh(ushst)
    ush = dsinh(ysh)
    bshock = ush/dsqrt(1d0+ush*ush)
    do ichi = 1, nchi
        call q_point(kmed,rshell,gf,q_grid(ichi), &
                              rcell(ichi),gcell(ichi),bcell(ichi))
        brelsh(ichi) = (bshock-bcell(ichi))/(1d0-bshock*bcell(ichi))
    end do
end subroutine q_cell_geometry

! 将q网格转换为物理深度和局域下游共动系空间网格。
! Convert the q grid to lab-frame depth and local downstream comoving spatial grid.
subroutine downstream_grid(nchi,kmed,rshell,gf,q_face,q_grid, &
                                            xface,xcf,xcmid,dxco, &
                                            rcell,gcell,bcell,brelsh)
    integer, intent(in) :: nchi,kmed
    real(8), intent(in), dimension(0:nchi) :: q_face
    real(8), intent(in), dimension(nchi) :: q_grid
    real(8), intent(in) :: rshell,gf
    real(8), intent(out), dimension(0:nchi) :: xface,xcf
    real(8), intent(out), dimension(nchi) :: xcmid, dxco, rcell, gcell, bcell, brelsh
    real(8), dimension(0:nchi) :: rface
    real(8) :: gface,bface,dx_lab
    integer :: ichi

    do ichi = 0, nchi
        call q_point(kmed,rshell,gf,q_face(ichi),rface(ichi),gface,bface)
        xface(ichi) = rshell-rface(ichi)
    end do
    call q_cell_geometry(nchi,kmed,rshell,gf,q_grid,rcell,gcell,bcell,brelsh)
    xcf(0) = 0d0
    do ichi = 1, nchi
        dx_lab = xface(ichi)-xface(ichi-1)
        dxco(ichi) = gcell(ichi)*dx_lab
        xcf(ichi) = xcf(ichi-1) + dxco(ichi)
        xcmid(ichi) = 0.5d0*(xcf(ichi-1)+xcf(ichi))
    end do
end subroutine downstream_grid

! 计算激波输运状态量：激波速度β_sh，下游速度β_2及其相对激波的速度β_2_sh。
! Compute shock transport state: shock beta, downstream beta, and downstream velocity relative to the shock.
subroutine shock_state(gsh, bsh, b2, b2sh)
    real(8), intent(in) :: gsh
    real(8), intent(out) :: bsh, b2, b2sh
    real(8) :: g2

    bsh = dsqrt(1d0-1d0/gsh**2)
    g2 = gsh/dsqrt(2d0)
    if (g2 > 1d0) then
        b2 = dsqrt(1d0-1d0/g2**2)
    else
        b2 = 0d0
    end if
    b2sh = (bsh-b2)/(1d0-bsh*b2)
end subroutine shock_state

! q下游局部散度：用桥接半径/速度场计算(∇·v)/(3β_f c)，返回ln γ方程系数。
! Downstream q-local divergence: compute the ln-gamma adiabatic coefficient from bridged radius/velocity fields.
subroutine q_divergence(nchi,kmed,rloc,gf,bf,q_grid,adlog)
    integer, intent(in) :: nchi,kmed
    real(8), intent(in), dimension(nchi) :: q_grid
    real(8), intent(in) :: rloc,gf,bf
    real(8), intent(out), dimension(nchi) :: adlog
    real(8), dimension(nchi) :: rcell,gcell,bcell,brelsh
    real(8) :: dbdr,divc
    integer :: ichi

    call q_cell_geometry(nchi,kmed,rloc,gf,q_grid,rcell,gcell,bcell,brelsh)
    dbdr = (bcell(2)-bcell(1))/(rcell(2)-rcell(1))
    divc = 2d0*bcell(1)/rcell(1) + dbdr
    adlog(1) = divc/(3d0*bf)
    do ichi = 2, nchi-1
        dbdr = (bcell(ichi+1)-bcell(ichi-1))/(rcell(ichi+1)-rcell(ichi-1))
        divc = 2d0*bcell(ichi)/rcell(ichi) + dbdr
        adlog(ichi) = divc/(3d0*bf)
    end do
    dbdr = (bcell(nchi)-bcell(nchi-1))/(rcell(nchi)-rcell(nchi-1))
    divc = 2d0*bcell(nchi)/rcell(nchi) + dbdr
    adlog(nchi) = divc/(3d0*bf)
end subroutine q_divergence

! 估计q方向对流子步长限制。
! Estimate the q-advection substep limit.
real(8) function q_step_limit(nchi,kmed,rloc,dq,q_face,cfl)
    integer, intent(in) :: nchi,kmed
    real(8), intent(in), dimension(0:nchi) :: q_face
    real(8), intent(in) :: rloc,dq,cfl
    real(8) :: aqface,maxq,qactive
    integer :: ichi

    maxq = 0d0
    qactive = 1d0-(1d0-1d0/shock_comp)**shock_comp
    do ichi = 1, nchi
        aqface = dble(3-kmed)*(qactive-q_face(ichi))/rloc
        maxq = max(maxq,dabs(aqface))
    end do
    q_step_limit = huge(1d0)
    if (maxq > 0d0) q_step_limit = cfl*dq/maxq
end function q_step_limit

! 线性速度场内回溯到单元边界所需的路径长度。
! Path length required to trace back to a cell boundary in a linear velocity field.
real(8) function eta_hit_time(etastart, etabound, acell, bcell)
    real(8), intent(in) :: etastart, etabound, acell, bcell
    real(8) :: shift, num, den

    if (dabs(bcell) <= 1d-30) then
        eta_hit_time = (etabound-etastart)/(-acell)
    else
        shift = acell/bcell
        num = etastart+shift
        den = etabound+shift
        eta_hit_time = dlog(num/den)/bcell
    end if
end function eta_hit_time

! 回溯 q 面的特征线位置，用于保守 PPM 重映射。
! Trace q-face characteristic positions for conservative PPM remapping.
subroutine trace_faces(nchi, dR_step, etaface, aetaface, etaback)
    integer, intent(in) :: nchi
    real(8), intent(in), dimension(0:nchi) :: etaface,aetaface
    real(8), intent(in) :: dR_step
    real(8), intent(out), dimension(0:nchi) :: etaback
    integer :: iface, icell
    real(8) :: etacur, etatrial, etabound, srem, hit, acell, bcell, acur, deta

    etaback(0) = etaface(0)
    do iface = 1, nchi
        etacur = etaface(iface)
        acur = aetaface(iface)
        if (acur == 0d0) then
            etaback(iface) = etacur
            cycle
        end if
        if (acur < 0d0 .and. iface == nchi) then
            etaback(iface) = etaface(nchi)
            cycle
        end if
        if (acur > 0d0) then
            icell = iface
        else
            icell = iface + 1
        end if
        srem = dR_step
        do while (srem > 0d0)
            deta = etaface(icell)-etaface(icell-1)
            bcell = (aetaface(icell)-aetaface(icell-1))/deta
            acell = aetaface(icell-1)-bcell*etaface(icell-1)
            acur = acell+bcell*etacur
            if (acur == 0d0) exit
            ! 线性 A(eta) 场的解析回溯。
            ! Analytic traceback in a linear A(eta) field.
            if (dabs(bcell) <= 1d-30) then
                etatrial = etacur - acell*srem
            else
                etatrial = (etacur + acell/bcell)*dexp(-bcell*srem) - acell/bcell
            end if
            if (acur > 0d0) then
                etabound = etaface(icell-1)
                if (etatrial > etabound) then
                    etacur = etatrial
                    exit
                end if
                hit = eta_hit_time(etacur, etabound, acell, bcell)
                etacur = etabound
                srem = srem-hit
                icell = icell-1
                if (icell < 1) exit
            else
                etabound = etaface(icell)
                if (etatrial < etabound) then
                    etacur = etatrial
                    exit
                end if
                hit = eta_hit_time(etacur, etabound, acell, bcell)
                etacur = etabound
                srem = srem-hit
                icell = icell+1
                if (icell > nchi) exit
            end if
        end do
        etaback(iface) = etacur
    end do
end subroutine trace_faces

! 把 q 面速度拆成左右迎风通量系数。
! Split q-face speeds into left/right upwind flux coefficients.
subroutine split_q_faces(nchi,aqface,advleft,advrgt)
    implicit none
    integer, intent(in) :: nchi
    integer :: iface
    real(8), intent(in), dimension(1:nchi) :: aqface
    real(8), intent(out), dimension(1:nchi) :: advleft,advrgt

    advleft=0d0
    advrgt=0d0
    do iface=1,nchi-1
        if (aqface(iface) >= 0d0) then
            advleft(iface)=aqface(iface)
        else
            advrgt(iface)=aqface(iface)
        end if
    end do
    if (aqface(nchi) > 0d0) advleft(nchi)=aqface(nchi)
end subroutine split_q_faces

! q 深度度量的倒数 dq/ds，用于把空间扩散写成 q 坐标扩散。
! Inverse q-depth metric dq/ds, used to express spatial diffusion in q coordinates.
real(8) function q_inv_metric(kmed,rshell,gf,qloc)
    integer, intent(in) :: kmed
    real(8), intent(in) :: rshell,gf,qloc
    real(8) :: qtail, alpha, uf, w, chibm, logxbm, logxst, xrat, dchidq, dlogbm, dlogst, dlogdq, dxdq

    qtail = 1d0-qloc
    alpha = dble(4-kmed)/dble(3-kmed)
    uf = dsqrt(gf*gf-1d0)
    w = uf*uf/(1d0+uf*uf)
    chibm = qtail**(-alpha)
    logxbm = -(chibm-1d0)/(4d0*dble(4-kmed)*gf*gf)
    logxst = dlog(qtail)/(shock_comp*dble(3-kmed))
    xrat = dexp(w*logxbm+(1d0-w)*logxst)
    dchidq = alpha*chibm/qtail
    dlogbm = -dchidq/(4d0*dble(4-kmed)*gf*gf)
    dlogst = -1d0/(shock_comp*dble(3-kmed)*qtail)
    dlogdq = w*dlogbm+(1d0-w)*dlogst
    dxdq = xrat*dlogdq
    q_inv_metric = 1d0/(-rshell*dxdq)
end function q_inv_metric

! 构建 q 面扩散基系数；是否包含外边界逃逸由 outer 控制。
! Build q-face diffusion base coefficients; outer controls whether the outer escape face is included.
subroutine q_diff_faces(nchi,dq,q_face,kmed,rloc,gf,bsh,outer, &
                                   dlbase,drbase)
    implicit none
    integer, intent(in) :: nchi,kmed
    integer :: iface
    logical, intent(in) :: outer
    real(8), intent(in), dimension(0:nchi) :: q_face
    real(8), intent(in) :: dq,rloc,gf,bsh
    real(8), intent(out), dimension(1:nchi) :: dlbase,drbase
    real(8) :: dpref,dqds

    dlbase=0d0
    drbase=0d0
    do iface=1,nchi-1
        dqds = q_inv_metric(kmed,rloc,gf,q_face(iface))
        dpref = dqds*dqds/(bsh*para_c)
        dlbase(iface) = dpref/dq
        drbase(iface) = -dpref/dq
    end do
    if (outer) then
        dqds = q_inv_metric(kmed,rloc,gf,q_face(nchi))
        dpref = dqds*dqds/(bsh*para_c)
        dlbase(nchi) = dpref/dq
    end if
end subroutine q_diff_faces

! 构建 q 方向迎风对流三对角基矩阵。
! Build the q-direction upwind-advection tridiagonal base matrix.
subroutine build_q_adv(nchi, lamq, aqface, lbase, dbase, ubase)
    integer, intent(in) :: nchi
    real(8), intent(in), dimension(1:nchi) :: aqface
    real(8), intent(in) :: lamq
    real(8), intent(out), dimension(nchi) :: lbase,dbase,ubase
    real(8), dimension(1:nchi) :: advleft,advrgt
    integer :: ichi

    call split_q_faces(nchi,aqface,advleft,advrgt)
    lbase = 0d0
    dbase = 1d0
    ubase = 0d0
    do ichi = 1, nchi
        dbase(ichi) = dbase(ichi) + lamq*advleft(ichi)
        if (ichi < nchi) ubase(ichi) = ubase(ichi) + lamq*advrgt(ichi)
    end do
    do ichi = 2, nchi
        lbase(ichi) = lbase(ichi) - lamq*advleft(ichi-1)
        dbase(ichi) = dbase(ichi) - lamq*advrgt(ichi-1)
    end do
end subroutine build_q_adv

! 把 q 扩散项加入已有三对角矩阵。
! Add the q-diffusion term to an existing tridiagonal matrix.
subroutine add_q_diff(nchi, lamq, kappa_row, dlbase, &
                                     drbase, escape, lower, diag, upper)
    integer, intent(in) :: nchi
    logical, intent(in) :: escape
    real(8), intent(in), dimension(nchi) :: kappa_row
    real(8), intent(in) :: lamq
    real(8), intent(in), dimension(1:nchi) :: dlbase,drbase
    real(8), intent(inout), dimension(nchi) :: lower,diag,upper
    real(8) :: kface
    integer :: ichi

    do ichi = 1, nchi
        if (ichi < nchi) then
            kface = 0.5d0*(kappa_row(ichi)+kappa_row(ichi+1))
            diag(ichi) = diag(ichi) + lamq*kface*dlbase(ichi)
            upper(ichi) = upper(ichi) + lamq*kface*drbase(ichi)
        else if (escape) then
            diag(ichi) = diag(ichi) + lamq*kappa_row(ichi)*dlbase(ichi)
        end if
    end do
    do ichi = 2, nchi
        kface = 0.5d0*(kappa_row(ichi-1)+kappa_row(ichi))
        lower(ichi) = lower(ichi) - lamq*kface*dlbase(ichi-1)
        diag(ichi) = diag(ichi) - lamq*kface*drbase(ichi-1)
    end do
end subroutine add_q_diff

! 构建 q 输运共享基矩阵：对流基矩阵加扩散面系数。
! Build the shared q-transport base matrix: advection base plus diffusion face coefficients.
subroutine build_q_base(nchi,dq,q_face,kmed,rloc,gf,bsh,dR_step, &
                                         escape,lamq,dlbase,drbase, &
                                         lbase,dbase,ubase)
    integer, intent(in) :: nchi,kmed
    logical, intent(in) :: escape
    real(8), intent(in), dimension(0:nchi) :: q_face
    real(8), intent(in) :: dq,rloc,gf,bsh,dR_step
    real(8), intent(out), dimension(1:nchi) :: dlbase,drbase
    real(8), intent(out) :: lamq
    real(8), intent(out), dimension(nchi) :: lbase,dbase,ubase
    real(8), dimension(1:nchi) :: aqface
    real(8) :: qactive
    integer :: iface

    lamq = dR_step/dq
    qactive = 1d0-(1d0-1d0/shock_comp)**shock_comp
    do iface = 1, nchi
        aqface(iface) = dble(3-kmed)*(qactive-q_face(iface))/rloc
    end do
    call q_diff_faces(nchi,dq,q_face,kmed,rloc,gf,bsh,escape, &
                                 dlbase,drbase)
    call build_q_adv(nchi,lamq,aqface,lbase,dbase,ubase)
end subroutine build_q_base

! Thomas算法求解三对角线性方程组。
! Solve a tridiagonal linear system with the Thomas algorithm.
subroutine solve_tridiagonal(ncell, lower, diag, upper, rhs, sol)
    integer, intent(in) :: ncell
    real(8), intent(in), dimension(ncell) :: lower,diag,upper,rhs
    real(8), intent(out), dimension(ncell) :: sol
    real(8), dimension(ncell) :: cprime,dprime
    real(8) :: denom
    integer :: icell

    denom = diag(1)
    cprime(1) = upper(1)/denom
    dprime(1) = rhs(1)/denom

    do icell = 2, ncell
        denom = diag(icell) - lower(icell)*cprime(icell-1)
        if (icell < ncell) cprime(icell) = upper(icell)/denom
        dprime(icell) = (rhs(icell)-lower(icell)*dprime(icell-1))/denom
    end do

    sol(ncell) = dprime(ncell)
    do icell = ncell-1, 1, -1
        sol(icell) = dprime(icell) - cprime(icell)*sol(icell+1)
    end do
end subroutine solve_tridiagonal

! q方向纯对流推进（特征线积分）：PPM重构+特征线回溯。
! Pure q-advection step by characteristic integration: PPM reconstruction plus face traceback.
subroutine advance_q_charint(ulog, ng, nchi, active_hi, dq, q_face, &
                                       kmed, rloc, srcq, dR_step, n_threads)
    integer, intent(in) :: ng, nchi, active_hi, kmed, n_threads
    real(8), intent(inout), dimension(ng,nchi) :: ulog
    real(8), intent(in), dimension(0:nchi) :: q_face
    real(8), intent(in), dimension(ng) :: srcq
    real(8), intent(in) :: dq,rloc,dR_step

    real(8), dimension(0:nchi) :: aqface,q_back
    real(8) :: qactive
    real(8), dimension(nchi) :: qin, qout, ppm_left, ppm_right
    real(8), dimension(0:nchi) :: prefix
    integer, dimension(0:nchi) :: qcell
    integer :: ig, iface, left, right, mid

    qactive = 1d0-(1d0-1d0/shock_comp)**shock_comp
    do iface = 1, nchi
        aqface(iface) = dble(3-kmed)*(qactive-q_face(iface))/rloc
    end do
    aqface(0) = dble(3-kmed)*(qactive-q_face(0))/rloc
    call trace_faces(nchi, dR_step, q_face, aqface, q_back)
    do iface = 0, nchi
        if (q_back(iface) <= q_face(0)) then
            qcell(iface) = 1
        else if (q_back(iface) >= q_face(nchi)) then
            qcell(iface) = nchi
        else
            left = 1
            right = nchi
            do while (left < right)
                mid = (left+right)/2
                if (q_face(mid) >= q_back(iface)) then
                    right = mid
                else
                    left = mid+1
                end if
            end do
            qcell(iface) = left
        end if
    end do

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. active_hi*nchi >= 512) schedule(static) &
    !$OMP& private(ig,iface,qin,qout,ppm_left,ppm_right,prefix)
    do ig = 1, active_hi
        qin = ulog(ig, :)
        qin(1) = qin(1) + dR_step*srcq(ig)
        call prepare_remap(nchi, q_face, qin, ppm_left, ppm_right, prefix)
        do iface = 1, nchi
            qout(iface) = &
                (ppm_eval_cell(nchi, q_face, qin, &
                                                   ppm_left, ppm_right, prefix, q_back(iface), qcell(iface)) - &
                 ppm_eval_cell(nchi, q_face, qin, &
                                                   ppm_left, ppm_right, prefix, q_back(iface-1), qcell(iface-1))) / dq
        end do
        ulog(ig, :) = qout
    end do
    !$OMP END PARALLEL DO
end subroutine advance_q_charint

! q方向纯扩散推进（隐式）：三对角求解扩散项 ∂_q(κ∂_q U)。
! Pure q-diffusion implicit step: solve the tridiagonal discretization of d_q(kappa d_q U).
subroutine advance_q_diffusion(ulog, ng, nchi, active_hi, dq, q_face, kmed, &
                                        rloc, gf, bsh, kappa2_chi, dR_step, n_threads)
    integer, intent(in) :: ng, nchi, active_hi, kmed, n_threads
    real(8), intent(inout), dimension(ng,nchi) :: ulog
    real(8), intent(in), dimension(0:nchi) :: q_face
    real(8), intent(in) :: dq,rloc,gf
    real(8), intent(in), dimension(ng,nchi) :: kappa2_chi
    real(8), intent(in) :: bsh,dR_step

    real(8), dimension(1:nchi) :: dlbase,drbase
    real(8), dimension(nchi) :: lower,diag,upper,rhs,sol
    real(8) :: lamq
    integer :: ig

    lamq = dR_step/dq
    call q_diff_faces(nchi,dq,q_face,kmed,rloc,gf,bsh,.false., &
                                 dlbase,drbase)

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. active_hi*nchi >= 512) schedule(static) &
    !$OMP& private(ig,lower,diag,upper,rhs,sol)
    do ig = 1, active_hi
        lower = 0d0
        diag = 1d0
        upper = 0d0
        call add_q_diff(nchi,lamq,kappa2_chi(ig,:), &
                                       dlbase,drbase,.false.,lower,diag,upper)
        rhs = ulog(ig, :)
        call solve_tridiagonal(nchi, lower, diag, upper, rhs, sol)
        ulog(ig, :) = max(0d0, sol)
    end do
    !$OMP END PARALLEL DO
end subroutine advance_q_diffusion

! q方向对流+扩散联合隐式推进：迎风对流+中心扩散，三对角求解。
! Joint implicit q transport: upwind advection plus centered diffusion in a tridiagonal solve.
subroutine advance_q_implicit(ulog, ng, nchi, active_hi, dq, q_face, kmed, rloc, gf, &
                              bsh, kappa2_chi, srcq, dR_step, n_threads)
    integer, intent(in) :: ng, nchi, active_hi, kmed, n_threads
    real(8), intent(inout), dimension(ng,nchi) :: ulog
    real(8), intent(in), dimension(0:nchi) :: q_face
    real(8), intent(in) :: dq,rloc,gf
    real(8), intent(in), dimension(ng,nchi) :: kappa2_chi
    real(8), intent(in), dimension(ng) :: srcq
    real(8), intent(in) :: bsh, dR_step

    real(8), dimension(1:nchi) :: dlbase,drbase
    real(8), dimension(nchi) :: lbase, dbase, ubase, lower, diag, upper, rhs, sol
    real(8) :: lamq
    integer :: ig

    call build_q_base(nchi,dq,q_face,kmed,rloc,gf,bsh,dR_step,.false., &
                                       lamq,dlbase,drbase, &
                                       lbase,dbase,ubase)

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. active_hi*nchi >= 512) schedule(static) &
    !$OMP& private(ig,lower,diag,upper,rhs,sol)
    do ig = 1, active_hi
        lower = lbase
        diag = dbase
        upper = ubase
        call add_q_diff(nchi,lamq,kappa2_chi(ig,:), &
                                       dlbase,drbase,.false.,lower,diag,upper)
        rhs = ulog(ig, :)
        rhs(1) = rhs(1) + dR_step*srcq(ig)

        call solve_tridiagonal(nchi, lower, diag, upper, rhs, sol)
        ulog(ig, :) = max(0d0, sol)
    end do
    !$OMP END PARALLEL DO
end subroutine advance_q_implicit

! PWN/CR q方向隐式输运：保守扩散，q=1为最深扫掠质量边界。
! PWN/CR implicit q transport: conservative diffusion with q=1 as the deepest swept-mass boundary.
subroutine advance_pwncr_q(ulog, ng, nchi, active_hi, dq, q_face, kmed, &
                                    rloc, gf, bsh, kappa2_chi, srcq, &
                                    dR_step, escape, n_threads)
    integer, intent(in) :: ng, nchi, active_hi, kmed, n_threads
    logical, intent(in) :: escape
    real(8), intent(inout), dimension(ng,nchi) :: ulog
    real(8), intent(in), dimension(0:nchi) :: q_face
    real(8), intent(in) :: dq,rloc,gf,bsh
    real(8), intent(in), dimension(ng,nchi) :: kappa2_chi
    real(8), intent(in), dimension(ng) :: srcq
    real(8), intent(in) :: dR_step
    real(8), dimension(1:nchi) :: dlbase,drbase
    real(8), dimension(nchi) :: lbase, dbase, ubase, lower, diag, upper, rhs, sol
    real(8) :: lamq
    integer :: ig

    call build_q_base(nchi,dq,q_face,kmed,rloc,gf,bsh,dR_step,escape, &
                                       lamq,dlbase,drbase, &
                                       lbase,dbase,ubase)

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. active_hi*nchi >= 512) schedule(static) &
    !$OMP& private(ig,lower,diag,upper,rhs,sol)
    do ig = 1, active_hi
        lower = lbase
        diag = dbase
        upper = ubase
        call add_q_diff(nchi,lamq,kappa2_chi(ig,:), &
                                       dlbase,drbase,escape,lower,diag,upper)
        rhs = ulog(ig, :)
        rhs(1) = rhs(1) + dR_step*srcq(ig)
        call solve_tridiagonal(nchi, lower, diag, upper, rhs, sol)
        ulog(ig, :) = max(0d0, sol)
    end do
    !$OMP END PARALLEL DO
end subroutine advance_pwncr_q

! 能量维（log γ）冷却推进：隐式迎风格式，对每个q柱独立求解三对角。
! Energy-dimension cooling step in log gamma: independent implicit upwind tridiagonal solve per q column.
subroutine advance_energy_chi(ulog, ng, nchi, dloss, rloc, dxg, dR_step, n_threads)
    integer, intent(in) :: ng, nchi, n_threads
    real(8), intent(inout), dimension(ng,nchi) :: ulog
    real(8), intent(in), dimension(ng-1,nchi) :: dloss
    real(8), intent(in) :: rloc,dxg,dR_step

    real(8), dimension(ng-1) :: xicoef,up
    real(8), dimension(ng) :: principal
    real(8), dimension(ng-1) :: temp1
    real(8), dimension(ng) :: rhs,sol
    real(8) :: CFL,adcoef
    integer :: ichi

    CFL = dR_step/dxg
    adcoef = 1d0/(rloc)

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. nchi*ng >= 512) schedule(static) &
    !$OMP& private(ichi,xicoef,up,principal,temp1,rhs,sol)
    do ichi = 1, nchi
        xicoef = dloss(:, ichi) + adcoef
        up = -CFL*xicoef
        call prepare_implicit_coeffs(ng, 1d0, up, principal, temp1)
        rhs = ulog(:, ichi) / principal
        call backward_sweep(ng, temp1, rhs, sol)
        ulog(:, ichi) = max(0d0, sol)
    end do
    !$OMP END PARALLEL DO
end subroutine advance_energy_chi

! PWN/CR能量维冷却：辐射损失加q局部散度给出的柱依赖绝热项。
! PWN/CR energy cooling: radiative loss plus q-column-dependent adiabatic term from local divergence.
subroutine advance_pwncr_energy(ulog, ng, nchi, dloss, adlog, &
                                             dxg, dR_step, n_threads)
    integer, intent(in) :: ng, nchi, n_threads
    real(8), intent(inout), dimension(ng,nchi) :: ulog
    real(8), intent(in), dimension(ng-1,nchi) :: dloss
    real(8), intent(in), dimension(nchi) :: adlog
    real(8), intent(in) :: dxg, dR_step
    real(8), dimension(ng-1) :: xicoef,up
    real(8), dimension(ng) :: principal
    real(8), dimension(ng-1) :: temp1
    real(8), dimension(ng) :: rhs,sol
    real(8) :: CFL
    integer :: ichi

    CFL = dR_step/dxg
    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. nchi*ng >= 512) schedule(static) &
    !$OMP& private(ichi,xicoef,up,principal,temp1,rhs,sol)
    do ichi = 1, nchi
        xicoef = dloss(:, ichi) + adlog(ichi)
        up = -CFL*xicoef
        call prepare_implicit_coeffs(ng, 1d0, up, principal, temp1)
        rhs = ulog(:, ichi) / principal
        call backward_sweep(ng, temp1, rhs, sol)
        ulog(:, ichi) = max(0d0, sol)
    end do
    !$OMP END PARALLEL DO
end subroutine advance_pwncr_energy

! 随机再加速：log γ空间零通量保守扩散步，用于Strang分裂的半步或整步。
! Stochastic re-acceleration: 0d0-flux conservative diffusion in log gamma for Strang half/full steps.
subroutine advance_stoch_chi(ulog, ng, nchi, stoch, &
                                                  rloc, dxg, dR_step, n_threads)
    integer, intent(in) :: ng, nchi, n_threads
    real(8), intent(inout), dimension(ng,nchi) :: ulog
    real(8), intent(in) :: stoch, rloc, dxg, dR_step
    real(8), dimension(ng) :: lower,diag,upper,rhs,sol
    real(8) :: lamg
    integer :: ichi, ig

    lamg = dR_step*stoch/(rloc*dxg*dxg)
    lower = 0d0
    diag = 1d0
    upper = 0d0
    do ig = 1, ng
        if (ig > 1) lower(ig) = -lamg
        if (ig < ng) upper(ig) = -lamg
        if (ig == 1 .or. ig == ng) then
            diag(ig) = 1d0 + lamg
        else
            diag(ig) = 1d0 + 2d0*lamg
        end if
    end do

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. nchi*ng >= 512) schedule(static) &
    !$OMP& private(ichi,rhs,sol)
    do ichi = 1, nchi
        rhs = ulog(:, ichi)
        call solve_tridiagonal(ng, lower, diag, upper, rhs, sol)
        ulog(:, ichi) = max(0d0, sol)
    end do
    !$OMP END PARALLEL DO
end subroutine advance_stoch_chi

! 能量维冷却推进（特征线积分版）：对每个χ柱独立做电子冷却特征线更新。
! Characteristic energy-cooling step: update electron cooling independently in every chi column.
subroutine advance_energy_charint(ulog, ng, nchi, gam_e, DB_chi, dEl_chi, rloc, &
                                               gsh, bsh, index_Y, dR_step, chihi, n_threads, &
                                               srcq)
    integer, intent(in) :: ng, nchi, index_Y, n_threads
    integer, intent(in), optional :: chihi
    real(8), intent(inout), dimension(ng,nchi) :: ulog
    real(8), intent(in), dimension(ng) :: gam_e
    real(8), intent(in), dimension(nchi) :: DB_chi
    real(8), intent(in), dimension(ng,nchi) :: dEl_chi
    real(8), intent(in) :: rloc, gsh, bsh, dR_step
    real(8), intent(in), optional, dimension(ng) :: srcq

    real(8), dimension(ng+1) :: x_edge
    real(8), dimension(ng) :: uin,uout,srczero
    real(8) :: arad, adcoef, ushock
    integer :: ichi, chi_hi, coolmode
    logical :: srccol

    call log_edges(ng, gam_e, x_edge)
    srczero = 0d0

    chi_hi = nchi
    if (present(chihi)) chi_hi = max(1, min(nchi, chihi))
    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. chi_hi*ng >= 512) schedule(static) &
    !$OMP& private(ichi,uin,uout,srccol,coolmode,arad,adcoef,ushock)
    do ichi = 1, chi_hi
        uin = ulog(:, ichi)
        srccol = present(srcq) .and. ichi == 1
        coolmode = cooling_piecewise
        arad = 0d0
        adcoef = 0d0
        if (index_Y == 0) then
            ushock = bsh*gsh
            if (ushock <= 0d0) error stop "advance_energy_charint requires bsh*gsh > 0"
            arad = 1.35d-19*DB_chi(ichi)**2/(ushock*pi)
            adcoef = 1d0/rloc
            coolmode = cooling_affine
        end if
        if (srccol) then
            call char_update(ng, dR_step, x_edge, coolmode, &
                                                 arad, adcoef, gam_e, dEl_chi(:,ichi), rloc, &
                                                 1d0, srcq, uin, uout)
        else
            call char_update(ng, dR_step, x_edge, coolmode, &
                                                arad, adcoef, gam_e, dEl_chi(:,ichi), rloc, &
                                                0d0, srczero, uin, uout)
        end if
        ulog(:, ichi) = uout
    end do
    !$OMP END PARALLEL DO
end subroutine advance_energy_charint

end module electron_transport_2d
