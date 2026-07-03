!f2py: skip
module hadronic_base
    use constants
    use quantum_synch, only: quantum_chi_parameter, quantum_cooling
    implicit none
    private

    ! 粒子静质量 [GeV/c^2]，来自 constants。
    ! Rest masses [GeV/c^2] imported from the shared constants module.
    real(8), parameter, public :: electron_m = Para_m_e_GeV
    real(8), parameter, public :: proton_m = Para_m_p_GeV
    real(8), parameter, public :: neutron_m = Para_m_n_GeV
    real(8), parameter, public :: pion_m = Para_m_pi_charged_GeV
    real(8), parameter, public :: pion0_m = Para_m_pi0_GeV
    real(8), parameter, public :: muon_m = Para_m_mu_GeV

    public :: build_grid, source_bounds, build_edges
    public :: shell_dt, dyn_time, proton_limit
    public :: check_grid, quant_factor

contains

! 构建质子对数均匀 Lorentz factor 网格。
! Build a logarithmically uniform proton Lorentz-factor grid.
subroutine build_grid(ng,gmin,gmax,gamma)
    integer, intent(in) :: ng
    real(8), intent(in) :: gmin,gmax
    real(8), dimension(ng), intent(out) :: gamma
    integer :: ig
    real(8) :: xmin,xmax,dx

    if (ng < 1) error stop "hadronic gamma grid must contain at least 1 point."
    if (gmin <= 1d0) error stop "hadronic gamma grid minimum must exceed 1."
    if (gmax <= gmin) error stop "hadronic gamma grid maximum must exceed minimum."
    if (ng == 1) then
        gamma(1) = gmin
        return
    end if

    xmin = dlog10(gmin)
    xmax = dlog10(gmax)
    dx = (xmax - xmin)/dble(ng - 1)
    do ig=1,ng
        gamma(ig) = 1d1**(xmin + dx*dble(ig - 1))
    end do
end subroutine build_grid

! 查找注入范围 [gmin, gmax] 在 Lorentz factor 网格上的索引边界。
! Find the grid-index bounds for the injection interval [gmin, gmax].
subroutine source_bounds(ng,gamma,gmin,gmax,ilo,ihi)
    integer, intent(in) :: ng
    real(8), dimension(ng), intent(in) :: gamma
    real(8), intent(in) :: gmin,gmax
    integer, intent(out) :: ilo,ihi
    integer :: ig

    ilo = 1
    ihi = ng
    do ig=1,ng
        if (gamma(ig) >= gmin) then
            ilo = ig
            exit
        end if
    end do
    do ig=ng,1,-1
        if (gamma(ig) <= gmax) then
            ihi = ig
            exit
        end if
    end do
end subroutine source_bounds

! 由网格中心构造 cell edge；内部边界取相邻中心的几何平均。
! Build cell edges from centers with geometric means between adjacent centers.
subroutine build_edges(ng,gamma,edge)
    integer, intent(in) :: ng
    real(8), dimension(ng), intent(in) :: gamma
    real(8), dimension(ng+1), intent(out) :: edge
    integer :: ig

    if (ng == 1) then
        if (gamma(1) <= 1d0) error stop "hadronic gamma center must exceed 1."
        edge(1) = 0.5d0*gamma(1)
        edge(2) = 2d0*gamma(1)
        return
    end if

    if (gamma(1) < 1d0) error stop "hadronic gamma grid minimum must be at least 1."
    edge(1) = gamma(1)*dsqrt(gamma(1)/gamma(2))
    do ig=2,ng
        edge(ig) = dsqrt(gamma(ig-1)*gamma(ig))
    end do
    edge(ng+1) = gamma(ng)*dsqrt(gamma(ng)/gamma(ng-1))
end subroutine build_edges

! 计算壳层观测者时间步长 dt = t_obs(i) - t_obs(i-1)。
! Return the observer-time shell step dt = t_obs(i) - t_obs(i-1).
real(8) function shell_dt(tobs,ishell)
    integer, intent(in) :: ishell
    real(8), dimension(*), intent(in) :: tobs

    if (ishell <= 1) then
        shell_dt = tobs(1)
    else
        shell_dt = tobs(ishell) - tobs(ishell-1)
    end if
    if (shell_dt <= 0d0) error stop "hadronic shell dt must be positive."
end function shell_dt

! 计算共动动力学时标 t_dyn = R / (Gamma * c)。
! Compute the comoving dynamical timescale t_dyn = R / (Gamma * c).
real(8) function dyn_time(radius,gamma)
    real(8), intent(in) :: radius,gamma

    if (radius <= 0d0) error stop "hadronic dynamical time requires positive radius."
    if (gamma < 1d0) error stop "hadronic dynamical time requires Gamma >= 1."
    dyn_time = radius/(gamma*Para_c)
end function dyn_time

! 由动力学限制和同步冷却限制估计质子最大 Lorentz factor。
! Estimate the proton maximum Lorentz factor from dynamical and synchrotron limits.
real(8) function proton_limit(bfield,tdyn,eta)
    real(8), intent(in) :: bfield,tdyn,eta
    real(8) :: gdyn,gsyn

    if (bfield <= 0d0) error stop "proton_limit requires bfield > 0."
    if (tdyn <= 0d0) error stop "proton_limit requires tdyn > 0."
    if (eta <= 0d0) error stop "proton_limit requires eta > 0."
    gdyn = Para_e*bfield*tdyn/(eta*Para_m_p*Para_c)
    gsyn = dsqrt(6d0*pi*Para_e/(eta*Para_sigmaT*bfield))*Para_m_p_DIV_m_e
    proton_limit = min(gdyn,gsyn)
end function proton_limit

! 验证正值、严格递增、对数均匀的能量网格。
! Validate a positive, strictly increasing, logarithmically uniform energy grid.
subroutine check_grid(ngrid,grid,name,dln)
    integer, intent(in) :: ngrid
    real(8), dimension(ngrid), intent(in) :: grid
    character(*), intent(in) :: name
    real(8), intent(out), optional :: dln
    integer :: i
    real(8) :: step,di

    if (ngrid < 2) error stop trim(name)//" must contain at least 2 points."
    if (grid(1) <= 0d0) error stop trim(name)//" must be strictly positive."
    do i=2,ngrid
        if (grid(i) <= 0d0) error stop trim(name)//" must be strictly positive."
        if (grid(i) <= grid(i-1)) error stop trim(name)//" must be strictly increasing."
    end do
    step = dlog(grid(2)/grid(1))
    do i=3,ngrid
        di = dlog(grid(i)/grid(i-1))
        if (dabs(di-step) > dmax1(1d-12,1d-6*dabs(step))) &
            error stop trim(name)//" must be logarithmically uniform."
    end do
    if (present(dln)) dln = step
end subroutine check_grid

! Landau 量子同步压制因子；chi 由量子同步模块计算。
! Landau quantum-synchrotron suppression factor with chi from the radiation kernel.
real(8) function quant_factor(gamma,bfield,mass)
    real(8), intent(in) :: gamma,bfield,mass
    real(8) :: chi

    if (bfield <= 0d0 .or. gamma <= 1d0) then
        quant_factor = 1d0
        return
    end if
    chi = quantum_chi_parameter(gamma,bfield,mass)
    quant_factor = quantum_cooling(chi)
end function quant_factor

end module hadronic_base
