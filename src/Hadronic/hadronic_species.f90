!f2py: skip
module hadronic_species
    use constants
    implicit none
    private

    real(8), parameter :: nmass = Para_m_n_GeV
    real(8), parameter :: pimass = Para_m_pi_charged_GeV
    real(8), parameter :: mumass = Para_m_mu_GeV
    real(8), parameter :: ntaudec = 879.4d0
    real(8), parameter :: pitaudec = 2.6033d-8
    real(8), parameter :: mutaudec = 2.1969811d-6
    real(8), parameter :: gev2g = Para_GeV2erg / (Para_c * Para_c)
    real(8), parameter :: dlntol = 1d-3

    public :: div_rate
    public :: synch_loss
    public :: ad_loss
    public :: species_advance

contains

! 球对称膨胀发散率：3*v_exp/r，用于绝热冷却。
! Spherical expansion divergence: 3*v_exp/r for adiabatic cooling.
real(8) function div_rate(rad,vexp)
    real(8), intent(in) :: rad, vexp

    if (rad <= 0d0) then
        error stop "hadronic species transport: rad must be positive."
    end if
    if (vexp < 0d0) then
        error stop "hadronic species transport: vexp must be non-negative."
    end if
    div_rate = 3d0*vexp/rad
end function div_rate

! 同步冷却率：charged species 使用 Thomson 缩放，中性粒子返回 0。
! Synchrotron cooling rate: charged species use Thomson scaling; neutral species return 0.
subroutine synch_loss(ng,gamma,bfield,mass,charge,dgdt)
    integer, intent(in) :: ng
    integer, intent(in) :: charge
    real(8), intent(in), dimension(ng) :: gamma
    real(8), intent(in) :: bfield,mass
    real(8), intent(out), dimension(ng) :: dgdt
    integer :: ig
    real(8) :: massg, ub, sigt

    call validate_gamma(ng,gamma)
    if (bfield < 0d0) then
        error stop "hadronic species transport: bfield must be non-negative."
    end if
    if (mass <= 0d0) then
        error stop "hadronic species transport: mass must be positive."
    end if
    if (charge == 0) then
        dgdt = 0d0
        return
    end if

    massg = mass*gev2g
    ub = bfield*bfield/(8d0*pi)
    sigt = Para_sigmaT*(charge**4)*(Para_m_e/massg)**2
    do ig=1,ng
        dgdt(ig) = -(4d0/3d0)*sigt*ub*gamma(ig)**2/(massg*Para_c)
    end do
end subroutine synch_loss

! 绝热冷却率 dgamma/dt = -(div v / 3) * gamma。
! Adiabatic cooling rate dgamma/dt = -(div v / 3) * gamma.
subroutine ad_loss(ng,gamma,divr,dgdt)
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: gamma
    real(8), intent(in) :: divr
    real(8), intent(out), dimension(ng) :: dgdt
    integer :: ig

    call validate_gamma(ng,gamma)
    if (divr < 0d0) then
        error stop "hadronic species transport: divr must be non-negative."
    end if

    do ig=1,ng
        dgdt(ig) = -(divr/3d0)*gamma(ig)
    end do
end subroutine ad_loss

! 七分量 species 输运：neutron、pi+/pi-、mu- L/R、mu+ L/R。
! Seven-component species transport: neutron, pi+/pi-, mu- L/R, and mu+ L/R.
subroutine species_advance(ng,gamma,dt,bfield,divr, &
                           n0,pip0,pim0,muml0,mumr0,mupl0,mupr0, &
                           nsrc,pipsrc,pimsrc,mumlsrc,mumrsrc,muplsrc,muprsrc, &
                           n1,pip1,pim1,muml1,mumr1,mupl1,mupr1)
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: gamma
    real(8), intent(in) :: dt,bfield,divr
    ! 上一步粒子谱；previous spectra for the seven species.
    real(8), intent(in), dimension(ng) :: n0,pip0,pim0,muml0,mumr0,mupl0,mupr0
    ! 源项；source terms for the seven species.
    real(8), intent(in), dimension(ng) :: nsrc,pipsrc,pimsrc,mumlsrc,mumrsrc,muplsrc,muprsrc
    ! 下一步粒子谱；advanced spectra for the seven species.
    real(8), intent(out), dimension(ng) :: n1,pip1,pim1,muml1,mumr1,mupl1,mupr1
    real(8), dimension(ng) :: dgpion,dgmuon,dgad,dgtot

    if (dt <= 0d0) then
        error stop "hadronic species transport: dt must be positive."
    end if

    call validate_inputs
    call synch_loss(ng,gamma,bfield,nmass,0,dgtot)
    call synch_loss(ng,gamma,bfield,pimass,1,dgpion)
    call synch_loss(ng,gamma,bfield,mumass,1,dgmuon)
    call ad_loss(ng,gamma,divr,dgad)

    dgtot = dgad
    call advance_one(ng,gamma,n0,nsrc,dt,ntaudec,dgtot,n1)

    dgtot = dgpion + dgad
    call advance_one(ng,gamma,pip0,pipsrc,dt,pitaudec,dgtot,pip1)
    call advance_one(ng,gamma,pim0,pimsrc,dt,pitaudec,dgtot,pim1)

    dgtot = dgmuon + dgad
    call advance_muons

contains

    subroutine validate_inputs
        call validate_gamma(ng,gamma)
        call validate_positive(ng,n0,"n0")
        call validate_positive(ng,pip0,"pip0")
        call validate_positive(ng,pim0,"pim0")
        call validate_positive(ng,muml0,"muml0")
        call validate_positive(ng,mumr0,"mumr0")
        call validate_positive(ng,mupl0,"mupl0")
        call validate_positive(ng,mupr0,"mupr0")
        call validate_positive(ng,nsrc,"nsrc")
        call validate_positive(ng,pipsrc,"pipsrc")
        call validate_positive(ng,pimsrc,"pimsrc")
        call validate_positive(ng,mumlsrc,"mumlsrc")
        call validate_positive(ng,mumrsrc,"mumrsrc")
        call validate_positive(ng,muplsrc,"muplsrc")
        call validate_positive(ng,muprsrc,"muprsrc")
    end subroutine validate_inputs

    subroutine advance_muons
        call advance_one(ng,gamma,muml0,mumlsrc,dt,mutaudec,dgtot,muml1)
        call advance_one(ng,gamma,mumr0,mumrsrc,dt,mutaudec,dgtot,mumr1)
        call advance_one(ng,gamma,mupl0,muplsrc,dt,mutaudec,dgtot,mupl1)
        call advance_one(ng,gamma,mupr0,muprsrc,dt,mutaudec,dgtot,mupr1)
    end subroutine advance_muons
end subroutine species_advance

! 单一 species 推进：half-decay -> upwind transport -> half-decay。
! Single-species advance: half-decay, upwind transport, then half-decay.
subroutine advance_one(ng,gamma,y0,qsrc,dt,tau0,dgdt,y1)
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: gamma,y0,qsrc
    real(8), intent(in) :: dt
    real(8), intent(in), dimension(ng) :: dgdt
    real(8), intent(in) :: tau0
    real(8), intent(out), dimension(ng) :: y1
    integer :: ig, iface, isub, nsub
    real(8) :: dlng, cmax, dtsub
    real(8), dimension(ng) :: taudec,nxi,qxi,ac
    real(8), dimension(ng+1) :: aface,flux
    real(8), dimension(ng+2) :: nedge
    real(8), dimension(ng) :: divf
    real(8), dimension(ng) :: halfdecay

    call validate_positive(ng,y0,"y0")
    call validate_positive(ng,qsrc,"qsrc")
    dlng = log_spacing(ng,gamma)

    do ig=1,ng
        taudec(ig) = tau0*gamma(ig)
        nxi(ig) = gamma(ig)*y0(ig)
        qxi(ig) = gamma(ig)*qsrc(ig)
        ac(ig) = dgdt(ig)/gamma(ig)
    end do

    aface(1) = ac(1)
    aface(ng+1) = ac(ng)
    do iface=2,ng
        aface(iface) = 0.5d0*(ac(iface-1)+ac(iface))
    end do

    cmax = 0d0
    do iface=1,ng+1
        cmax = max(cmax,abs(aface(iface))*dt/dlng)
    end do
    nsub = max(1,ceiling(cmax))
    dtsub = dt/dble(nsub)

    do ig=1,ng
        halfdecay(ig) = dexp(-0.5d0*dtsub/taudec(ig))
    end do

    do isub=1,nsub
        nxi = nxi*halfdecay
        nedge(1) = 0d0
        nedge(ng+2) = 0d0
        nedge(2:ng+1) = nxi

        do iface=1,ng+1
            if (aface(iface) >= 0d0) then
                flux(iface) = aface(iface)*nedge(iface)
            else
                flux(iface) = aface(iface)*nedge(iface+1)
            end if
        end do

        do ig=1,ng
            divf(ig) = (flux(ig+1)-flux(ig))/dlng
            nxi(ig) = nxi(ig) + dtsub*(-divf(ig)+qxi(ig))
            if (nxi(ig) < 0d0) nxi(ig) = 0d0
        end do
        nxi = nxi*halfdecay
    end do

    do ig=1,ng
        y1(ig) = nxi(ig)/gamma(ig)
    end do
end subroutine advance_one

! gamma 网格检查：至少两个点、全正、严格递增。
! Gamma-grid check: at least 2 points, positive values, strictly increasing.
subroutine validate_gamma(ng,gamma)
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: gamma
    integer :: ig

    if (ng < 2) then
        error stop "hadronic species transport: gamma grid must have at least 2 points."
    end if
    do ig=1,ng
        if (gamma(ig) <= 0d0) then
            error stop "hadronic species transport: gamma grid must be positive."
        end if
    end do
    do ig=2,ng
        if (gamma(ig) <= gamma(ig-1)) then
            error stop "hadronic species transport: gamma grid must be strictly increasing."
        end if
    end do
end subroutine validate_gamma

! 非负数组检查；内部 transport/source 谱不允许负值。
! Non-negative array check for internal transport and source spectra.
subroutine validate_positive(ng,values,name)
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: values
    character(len=*), intent(in) :: name
    integer :: ig

    do ig=1,ng
        if (values(ig) < 0d0) then
            error stop "hadronic species transport: "//trim(name)//" must be non-negative."
        end if
    end do
end subroutine validate_positive

! gamma 网格平均对数间距检查。
! Average logarithmic gamma-grid spacing check.
real(8) function log_spacing(ng,gamma)
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: gamma
    integer :: ig
    real(8) :: dlog_i, dlog_ref

    dlog_ref = 0d0
    do ig=1,ng-1
        dlog_ref = dlog_ref + dlog(gamma(ig+1)/gamma(ig))
    end do
    dlog_ref = dlog_ref/dble(ng-1)

    do ig=1,ng-1
        dlog_i = dlog(gamma(ig+1)/gamma(ig))
        if (abs(dlog_i-dlog_ref) > dlntol*abs(dlog_ref)) then
            error stop "hadronic species transport: gamma grid must be approximately log-spaced."
        end if
    end do
    log_spacing = dlog_ref
end function log_spacing

end module hadronic_species
