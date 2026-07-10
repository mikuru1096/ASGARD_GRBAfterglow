!f2py: skip
module hadronic_pp
    use constants
    use hadronic_base, only: check_grid
    implicit none
    private

    real(8), parameter :: pmass = Para_m_p_GeV
    real(8), parameter :: p0mass = Para_m_pi0_GeV
    real(8), parameter :: mbarn = 1.0d-27

    public :: pp_source
    public :: pp_sigma
    public :: pp_threshold

contains

! pp delta 近似：Kelner+2006 非弹性截面，输出 gamma、neutrino 和 e+/e- source。
! pp delta approximation using Kelner+2006 inelastic cross section.
subroutine pp_source(np,ep,pden,ntarget,ng,egam,nnu,enu,npair,epair, &
                    qgam,qnu,qpair,ploss,kappa,fpion,fneutral)
    integer, intent(in) :: np,ng,nnu,npair
    real(8), intent(in), dimension(np) :: ep,pden
    real(8), intent(in) :: ntarget
    real(8), intent(in), dimension(ng) :: egam
    real(8), intent(in), dimension(nnu) :: enu
    real(8), intent(in), dimension(npair) :: epair
    real(8), intent(out), dimension(ng) :: qgam
    real(8), intent(out), dimension(nnu) :: qnu
    real(8), intent(out), dimension(npair) :: qpair
    real(8), intent(out), dimension(np) :: ploss
    real(8), intent(in), optional :: kappa,fpion
    real(8), intent(in), optional :: fneutral
    real(8), dimension(np) :: siginel,coll,prate
    real(8) :: kloc,fpiloc,fneut,fch
    real(8) :: xgam,xnu,xpair

    call validate_inputs
    call set_options
    call pp_sigma(np,ep,siginel)
    coll = ntarget*Para_c*siginel
    prate = coll*pden
    ! Continuous proton cooling is a single-particle dgamma/dt; only secondary
    ! production depends on the proton distribution normalization.
    ploss = -kloc*coll*(ep/pmass - 1d0)

    call emit_secondaries

contains

    subroutine validate_inputs
        call check_grid(np,ep,"ep")
        call check_grid(ng,egam,"egam")
        call check_grid(nnu,enu,"enu")
        call check_grid(npair,epair,"epair")

        if (ntarget < 0d0) then
            error stop "pp_source: ntarget must be non-negative."
        end if
    end subroutine validate_inputs

    subroutine set_options
        kloc = 0.5d0
        if (present(kappa)) kloc = kappa
        if (kloc <= 0d0 .or. kloc > 1d0) then
            error stop "pp_source: kappa must be in (0, 1]."
        end if

        fpiloc = 0.17d0
        if (present(fpion)) fpiloc = fpion
        if (fpiloc <= 0d0 .or. fpiloc >= 1d0) then
            error stop "pp_source: fpion must be in (0, 1)."
        end if

        fneut = 1d0/3d0
        if (present(fneutral)) fneut = fneutral
        if (fneut < 0d0 .or. fneut > 1d0) then
            error stop "pp_source: fneutral must be in [0, 1]."
        end if

        fch = 1d0 - fneut
        xgam = 0.5d0*fpiloc
        xnu = 0.25d0*fpiloc
        xpair = xnu
    end subroutine set_options

    subroutine emit_secondaries
        call secondary_source(ng,egam,np,ep,prate,xgam,2d0*fneut,qgam)
        call secondary_source(nnu,enu,np,ep,prate,xnu,3d0*fch,qnu)
        call secondary_source(npair,epair,np,ep,prate,xpair,fch,qpair)
    end subroutine emit_secondaries
end subroutine pp_source

! Kelner+2006 pp 非弹性截面参数化。
! Kelner+2006 inelastic pp cross-section parametrization.
subroutine pp_sigma(np,ep,sigma)
    integer, intent(in) :: np
    real(8), intent(in), dimension(np) :: ep
    real(8), intent(out), dimension(np) :: sigma
    integer :: ip
    real(8) :: eth,lterm,cutoff

    call check_grid(np,ep,"ep")
    if (minval(ep) <= pmass) then
        error stop "pp_sigma: ep must exceed proton rest energy."
    end if

    eth = pp_threshold()
    sigma = 0d0
    do ip=1,np
        if (ep(ip) <= eth) cycle
        lterm = dlog(ep(ip)/1d3)
        cutoff = 1d0 - (eth/ep(ip))**4
        sigma(ip) = (34.3d0 + 1.88d0*lterm + 0.25d0*lterm*lterm) * &
                              cutoff*cutoff*mbarn
    end do
end subroutine pp_sigma

! pp 反应阈值总能量：E_th = m_p + 2 m_pi0 + m_pi0^2/(2 m_p)。
! Total proton energy at the pp pion-production threshold.
real(8) function pp_threshold()
    pp_threshold = pmass + 2d0*p0mass + p0mass*p0mass/(2d0*pmass)
end function pp_threshold

! delta 近似 source：secondary energy = frac * parent kinetic energy，再乘多重数。
! Delta-approximation source with secondary energy = frac * parent kinetic energy.
subroutine secondary_source(ns,es,npar,epar,rpar,frac,mult,qsec)
    integer, intent(in) :: ns,npar
    real(8), intent(in), dimension(ns) :: es
    real(8), intent(in), dimension(npar) :: epar
    real(8), intent(in), dimension(npar) :: rpar
    real(8), intent(in) :: frac,mult
    real(8), intent(out), dimension(ns) :: qsec
    real(8), dimension(ns) :: epeval,psamp

    epeval = pmass + es/frac
    call pos_interp(npar,epar,rpar,ns,epeval,psamp)
    qsec = (mult/frac)*psamp
end subroutine secondary_source

! 正值 log-log 插值：只使用 y>0 的 bracket。
! Positive log-log interpolation using only positive-y brackets.
subroutine pos_interp(nx,x,y,nnew,xnew,ynew)
    integer, intent(in) :: nx,nnew
    real(8), intent(in), dimension(nx) :: x,y
    real(8), intent(in), dimension(nnew) :: xnew
    real(8), intent(out), dimension(nnew) :: ynew
    integer, dimension(nx) :: posidx
    integer :: npos
    real(8), dimension(nx) :: xp,yp
    integer :: i,ipos

    call check_grid(nx,x,"interpolation x")

    npos = 0
    do i=1,nx
        if (y(i) > 0d0) then
            npos = npos + 1
            posidx(npos) = i
        end if
    end do

    ynew = 0d0
    if (npos < 2) return

    do i=1,npos
        xp(i) = x(posidx(i))
        yp(i) = y(posidx(i))
    end do

    do i=1,nnew
        if (xnew(i) < xp(1) .or. xnew(i) > xp(npos)) cycle
        ipos = upper_bracket(npos,xp,xnew(i))
        ynew(i) = dexp(log_eval(xp(ipos),xp(ipos+1),yp(ipos),yp(ipos+1),xnew(i)))
    end do
end subroutine pos_interp

! 单调数组 bracket 查找。
! Bracket search for a monotonic array.
integer function upper_bracket(nx,x,xeval)
    integer, intent(in) :: nx
    real(8), intent(in), dimension(nx) :: x
    real(8), intent(in) :: xeval
    integer :: ilo,ihi,imid

    if (xeval <= x(1)) then
        upper_bracket = 1
        return
    end if
    if (xeval >= x(nx)) then
        upper_bracket = nx - 1
        return
    end if

    ilo = 1
    ihi = nx
    do while (ihi - ilo > 1)
        imid = (ilo + ihi)/2
        if (x(imid) <= xeval) then
            ilo = imid
        else
            ihi = imid
        end if
    end do
    upper_bracket = ilo
end function upper_bracket

! 单次 log-log 线性插值。
! Single log-log linear interpolation.
real(8) function log_eval(x0,x1,y0,y1,xeval)
    real(8), intent(in) :: x0,x1,y0,y1,xeval
    real(8) :: lx0,lx1,ly0,ly1,lxeval,frac

    lx0 = dlog(x0)
    lx1 = dlog(x1)
    ly0 = dlog(y0)
    ly1 = dlog(y1)
    lxeval = dlog(xeval)
    frac = (lxeval - lx0)/(lx1 - lx0)
    log_eval = ly0 + frac*(ly1 - ly0)
end function log_eval

end module hadronic_pp
