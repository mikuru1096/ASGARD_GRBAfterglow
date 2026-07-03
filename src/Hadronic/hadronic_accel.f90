!f2py: skip
module hadronic_accel
    use constants
    use hadronic_base, only: proton_m, neutron_m, pion_m, muon_m
    implicit none
    private

    real(8), parameter :: gevgram = 1d9 * Para_eV2erg / (Para_c * Para_c)

    public :: species_info
    public :: accel_time
    public :: syn_time
    public :: ext_time
    public :: inject_rate
    public :: gamma_limit
    public :: accel_calc

contains

! 查询 hadronic species 的质量和电荷。
! Return mass and charge for each hadronic species.
subroutine species_info(species,mgev,charge,mg,qabs)
    character(len=*), intent(in) :: species
    real(8), intent(out) :: mgev,mg,qabs
    integer, intent(out) :: charge

    select case (trim(species))
    case ("proton")
        mgev = proton_m
        charge = 1
    case ("neutron")
        mgev = neutron_m
        charge = 0
    case ("pion_plus")
        mgev = pion_m
        charge = 1
    case ("pion_minus")
        mgev = pion_m
        charge = -1
    case ("muon_plus")
        mgev = muon_m
        charge = 1
    case ("muon_minus")
        mgev = muon_m
        charge = -1
    case default
        error stop "species_info: unsupported species."
    end select

    mg = mgev*gevgram
    qabs = dabs(dble(charge))*Para_e
end subroutine species_info

! 费米加速时标 t_acc = eta * gamma * m * c / (|q| * B)。
! Fermi acceleration time in a magnetic field.
subroutine accel_time(ng,species,gamma,bfield,eta,tacc)
    integer, intent(in) :: ng
    character(len=*), intent(in) :: species
    real(8), intent(in), dimension(ng) :: gamma
    real(8), intent(in) :: bfield,eta
    real(8), intent(out), dimension(ng) :: tacc
    integer :: ig,charge
    real(8) :: mgev,mg,qabs

    call species_info(species,mgev,charge,mg,qabs)
    if (bfield <= 0d0) error stop "accel_time: bfield must be > 0."
    if (eta <= 0d0) error stop "accel_time: eta must be > 0."
    if (qabs <= 0d0) error stop "accel_time: species is neutral."

    do ig=1,ng
        if (gamma(ig) <= 0d0) error stop "accel_time: gamma must be > 0."
        tacc(ig) = eta*gamma(ig)*mg*Para_c/(qabs*bfield)
    end do
end subroutine accel_time

! 同步辐射冷却时标。
! Synchrotron cooling time for the requested species.
subroutine syn_time(ng,species,gamma,bfield,tsyn)
    integer, intent(in) :: ng
    character(len=*), intent(in) :: species
    real(8), intent(in), dimension(ng) :: gamma
    real(8), intent(in) :: bfield
    real(8), intent(out), dimension(ng) :: tsyn
    integer :: ig,charge
    real(8) :: mgev,mg,qabs

    call species_info(species,mgev,charge,mg,qabs)
    if (bfield <= 0d0) error stop "syn_time: bfield must be > 0."

    do ig=1,ng
        if (gamma(ig) <= 0d0) error stop "syn_time: gamma must be > 0."
        tsyn(ig) = 6d0*pi*(mg**3)*Para_c / &
                   (Para_SigmaT*(Para_m_e**2)*(bfield**2)*gamma(ig))
    end do
end subroutine syn_time

! 外部光子场冷却时标 t_ext = gamma / rate。
! External-field cooling time from a tabulated cooling rate.
subroutine ext_time(ng,gamma,cool,text)
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: gamma,cool
    real(8), intent(out), dimension(ng) :: text
    integer :: ig

    do ig=1,ng
        if (gamma(ig) <= 0d0) error stop "ext_time: gamma must be > 0."
        if (cool(ig) <= 0d0) error stop "ext_time: cool must be > 0."
        text(ig) = gamma(ig)/cool(ig)
    end do
end subroutine ext_time

! 注入谱：幂律形状、可选指数截断，并用能量光度归一化。
! Injection spectrum: power law, optional cutoff, energy-luminosity normalization.
subroutine inject_rate(ng,gamma,species,lum,pidx,gmin,gmax,gcut,cut,qinj)
    integer, intent(in) :: ng
    character(len=*), intent(in) :: species
    real(8), intent(in), dimension(ng) :: gamma
    real(8), intent(in) :: lum,pidx,gmin,gmax,gcut
    logical, intent(in) :: cut
    real(8), intent(out), dimension(ng) :: qinj
    integer :: ig,charge
    real(8) :: mgev,mg,qabs,mener,norm,q0
    real(8), dimension(ng) :: shape,integrand

    call species_info(species,mgev,charge,mg,qabs)
    if (gmin <= 0d0 .or. gmax <= 0d0 .or. gmin >= gmax) then
        error stop "inject_rate: require 0 < gmin < gmax."
    end if
    if (lum <= 0d0) error stop "inject_rate: lum must be > 0."
    if (cut .and. gcut <= 0d0) error stop "inject_rate: gcut must be > 0."

    do ig=1,ng
        if (gamma(ig) <= 0d0) error stop "inject_rate: gamma must be > 0."
        shape(ig) = 0d0
        if (gamma(ig) >= gmin .and. gamma(ig) <= gmax) then
            shape(ig) = gamma(ig)**(-pidx)
            if (cut) shape(ig) = shape(ig)*dexp(-gamma(ig)/gcut)
        end if
        integrand(ig) = shape(ig)*gamma(ig)
    end do

    norm = trapz(ng,gamma,integrand)
    if (norm <= 0d0) error stop "inject_rate: normalization integral must be > 0."

    mener = mgev*Para_GeV2erg
    q0 = lum/(mener*norm)
    qinj = q0*shape
end subroutine inject_rate

! 最大洛伦兹因子：取动力学、同步冷却和外部冷却限制的最小值。
! Maximum Lorentz factor from dynamical, synchrotron, and external cooling limits.
subroutine gamma_limit(species,bfield,radius,gbulk,eta,nscan,gscan,xrate,use_ext, &
                       gmax,gdyn,gsyn,gext,has_xlim)
    character(len=*), intent(in) :: species
    real(8), intent(in) :: bfield,radius,gbulk,eta
    integer, intent(in) :: nscan
    real(8), intent(in), dimension(nscan) :: gscan,xrate
    logical, intent(in) :: use_ext
    real(8), intent(out) :: gmax,gdyn,gsyn,gext
    logical, intent(out) :: has_xlim
    integer :: charge,icross
    real(8) :: mgev,mg,qabs
    real(8), dimension(nscan) :: tacc,text
    real(8) :: rlo,rhi,x0,x1,y0,y1,xroot

    call species_info(species,mgev,charge,mg,qabs)
    if (qabs <= 0d0) error stop "gamma_limit: species is neutral."
    if (bfield <= 0d0 .or. radius <= 0d0 .or. gbulk <= 0d0 .or. eta <= 0d0) then
        error stop "gamma_limit: physical inputs must be > 0."
    end if

    call init_limit

    if (.not. use_ext) return

    call accel_time(nscan,species,gscan,bfield,eta,tacc)
    call ext_time(nscan,gscan,xrate,text)
    icross = find_cross()
    if (icross == 0) return

    call apply_xlimit(icross)

contains

    ! 无外部冷却时的动力学和同步辐射限制。
    ! Dynamical and synchrotron limits before external cooling is applied.
    subroutine init_limit
        real(8) :: tdyn

        tdyn = radius/(gbulk*Para_c)
        gdyn = qabs*bfield*tdyn/(eta*mg*Para_c)
        gsyn = dsqrt(6d0*pi*qabs*(mg**2)/(eta*Para_SigmaT*(Para_m_e**2)*bfield))

        has_xlim = .false.
        gext = 0d0
        gmax = dmin1(gdyn,gsyn)
    end subroutine init_limit

    ! 查找 t_acc 与 t_ext 的交叉区间。
    ! Find the bracket where acceleration and external cooling times cross.
    integer function find_cross()
        integer :: ig

        find_cross = 0
        do ig=1,nscan-1
            rlo = tacc(ig)-text(ig)
            rhi = tacc(ig+1)-text(ig+1)
            if (rlo*rhi <= 0d0) then
                find_cross = ig
                exit
            end if
        end do
    end function find_cross

    ! 用对数插值给出外部冷却限制。
    ! Apply the external cooling limit through log-space interpolation.
    subroutine apply_xlimit(icross)
        integer, intent(in) :: icross

        x0 = dlog(gscan(icross))
        x1 = dlog(gscan(icross+1))
        y0 = dlog(tacc(icross)/text(icross))
        y1 = dlog(tacc(icross+1)/text(icross+1))
        xroot = x0-y0*(x1-x0)/(y1-y0)
        gext = dexp(xroot)

        has_xlim = .true.
        gmax = dmin1(gmax,gext)
    end subroutine apply_xlimit
end subroutine gamma_limit

! 一次性计算加速、冷却、注入和最大能量限制。
! Compute acceleration, cooling, injection, and maximum-energy limits together.
subroutine accel_calc(ng,species,gamma,bfield,eta,lum,pidx,gmin,ginj,gcut,cut, &
                      radius,gbulk,nscan,gscan,xrate,use_ext,tacc,tsyn,qinj, &
                      gmax,gdyn,gsyn,gext,has_xlim)
    integer, intent(in) :: ng,nscan
    character(len=*), intent(in) :: species
    real(8), intent(in), dimension(ng) :: gamma
    real(8), intent(in) :: bfield,eta,lum,pidx,gmin,ginj,gcut,radius,gbulk
    real(8), intent(in), dimension(nscan) :: gscan,xrate
    logical, intent(in) :: cut,use_ext
    real(8), intent(out), dimension(ng) :: tacc,tsyn,qinj
    real(8), intent(out) :: gmax,gdyn,gsyn,gext
    logical, intent(out) :: has_xlim

    call accel_time(ng,species,gamma,bfield,eta,tacc)
    call syn_time(ng,species,gamma,bfield,tsyn)
    call inject_rate(ng,gamma,species,lum,pidx,gmin,ginj,gcut,cut,qinj)
    call gamma_limit(species,bfield,radius,gbulk,eta,nscan,gscan,xrate,use_ext, &
                     gmax,gdyn,gsyn,gext,has_xlim)
end subroutine accel_calc

! 梯形法则积分。
! Trapezoidal integration for a tabulated function.
real(8) function trapz(nx,x,y)
    integer, intent(in) :: nx
    real(8), intent(in), dimension(nx) :: x,y
    integer :: ix

    trapz = 0d0
    if (nx <= 1) return
    do ix=1,nx-1
        trapz = trapz + 5d-1*(y(ix)+y(ix+1))*(x(ix+1)-x(ix))
    end do
end function trapz

end module hadronic_accel
