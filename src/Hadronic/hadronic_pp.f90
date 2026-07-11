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
                    model,qgam,qnu,qpair,ploss,kappa,fpion,fneutral)
    integer, intent(in) :: np,ng,nnu,npair,model
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
    if (model >= 0) call pp_gamma(np,ep,pden,ntarget,ng,egam,model,qgam)

contains

    subroutine validate_inputs
        call check_grid(np,ep,"ep")
        call check_grid(ng,egam,"egam")
        call check_grid(nnu,enu,"enu")
        call check_grid(npair,epair,"epair")

        if (ntarget < 0d0) then
            error stop "pp_source: ntarget must be non-negative."
        end if
        if (model < -1 .or. model > 3) error stop "pp_source: invalid gamma model."
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

! Kafexhiu+2014 / AM3 v1.2.0 neutral-pion gamma spectrum. Proton density is dN/(dE dV).
! SIBYLL/QGSJET switch at the paper/AM3-normalization 100 GeV; AM3 F uses an inconsistent 50 GeV.
subroutine pp_gamma(np,ep,pden,ntarget,ng,egam,model,qgam)
    integer, intent(in) :: np,ng,model
    real(8), intent(in), dimension(np) :: ep,pden
    real(8), intent(in), dimension(ng) :: egam
    real(8), intent(in) :: ntarget
    real(8), intent(out), dimension(ng) :: qgam
    integer :: ip,ig
    real(8) :: dln,tp,amp,weight,emax,y0

    dln=dlog(ep(2)/ep(1))
    qgam=0d0
    do ip=1,np
        tp=ep(ip)-pmass
        if (tp <= pp_tkin()) cycle
        weight=ntarget*Para_c*pden(ip)*ep(ip)*dln
        amp=mbarn*pp_amax(tp,model)
        emax=pp_egmax(tp)
        y0=emax+p0mass*p0mass/(4d0*emax)
        do ig=1,ng
            qgam(ig)=qgam(ig)+weight*amp*pp_shape(tp,egam(ig),model,emax,y0)
        end do
    end do
end subroutine pp_gamma

! Peak differential cross section in mbarn/GeV.
real(8) function pp_amax(tp,model)
    real(8), intent(in) :: tp
    integer, intent(in) :: model
    real(8) :: theta,ltheta,b1,b2,b3

    if (tp < 1d0) then
        pp_amax=5.9d0*pp_sigpi(tp,model)/pp_epimax(tp)
        return
    end if
    if (model == 0 .and. tp > 100d0) then
        b1=10.77d0; b2=-0.412d0; b3=0.01264d0
    else if (model == 1 .and. tp > 100d0) then
        b1=13.16d0; b2=-0.4419d0; b3=0.01439d0
    else if (model == 3 .and. tp > 50d0) then
        b1=9.06d0; b2=-0.3795d0; b3=0.01105d0
    else if (tp < 5d0) then
        b1=9.53d0; b2=-0.52d0; b3=0.054d0
    else
        b1=9.13d0; b2=-0.35d0; b3=0.0097d0
    end if
    theta=tp/pmass
    ltheta=dlog(theta)
    pp_amax=b1*theta**b2*dexp(b3*ltheta*ltheta)*pp_sigpi(tp,model)/pmass
end function pp_amax

! Dimensionless spectral shape. SIBYLL/QGSJET switch at the published 100 GeV.
real(8) function pp_shape(tp,egam,model,emax,y0)
    real(8), intent(in) :: tp,egam,emax,y0
    integer, intent(in) :: model
    real(8) :: y,x,c

    pp_shape=0d0
    if (egam <= 0d0 .or. egam >= emax) return
    y=egam+p0mass*p0mass/(4d0*egam)
    x=(y-p0mass)/(y0-p0mass)
    if (x < 0d0 .or. x >= 1d0) return
    if (model == 0 .and. tp > 100d0) then
        c=3.55d0*p0mass/y0
        pp_shape=(1d0-dsqrt(x))**3.6d0/(1d0+x/c)
    else if (model == 1 .and. tp > 100d0) then
        c=3.55d0*p0mass/y0
        pp_shape=(1d0-dsqrt(x))**4.5d0/(1d0+x/c)
    else if (model == 3 .and. tp > 50d0) then
        c=3.5d0*p0mass/y0
        pp_shape=(1d0-dsqrt(x))**4d0/(1d0+x/c)
    else
        pp_shape=pp_geant(tp,x,y0)
    end if
end function pp_shape

real(8) function pp_geant(tp,x,y0)
    real(8), intent(in) :: tp,x,y0
    real(8) :: c,q,mu,kappa

    c=3d0*p0mass/y0
    if (tp < 1d0) then
        kappa=3.29d0-0.2d0*(tp/pmass)**(-1.5d0)
        pp_geant=(1d0-x)**kappa
    else if (tp <= 4d0) then
        q=(tp-1d0)/pmass
        mu=1.25d0*q**1.25d0*dexp(-1.25d0*q)
        pp_geant=(1d0-x)**(mu+2.45d0)/(1d0+x/c)**(mu+1.45d0)
    else if (tp <= 20d0) then
        q=(tp-1d0)/pmass
        mu=1.25d0*q**1.25d0*dexp(-1.25d0*q)
        pp_geant=(1d0-x)**(1.5d0*mu+4.95d0)/(1d0+x/c)**(mu+1.5d0)
    else if (tp <= 100d0) then
        pp_geant=(1d0-dsqrt(x))**4.2d0/(1d0+x/c)
    else
        pp_geant=(1d0-dsqrt(x))**4.9d0/(1d0+x/c)
    end if
end function pp_geant

! Inclusive pi0 production cross section in mbarn.
real(8) function pp_sigpi(tp,model)
    real(8), intent(in) :: tp
    integer, intent(in) :: model
    pp_sigpi=pp_sig1(tp)+pp_sig2(tp)+pp_siginel(tp)*pp_mult(tp,model)
end function pp_sigpi

real(8) function pp_sig1(tp)
    real(8), intent(in) :: tp
    real(8), parameter :: mres=1.1883d0, gres0=0.2264d0
    real(8) :: s,x,eta,gres,kres,fbw

    pp_sig1=0d0
    if (tp <= pp_tkin() .or. tp > 2d0) return
    s=2d0*pmass*(tp+2d0*pmass)
    x=dsqrt(s)-pmass
    eta=dsqrt((s-p0mass*p0mass-4d0*pmass*pmass)**2- &
              16d0*p0mass*p0mass*pmass*pmass)/(2d0*p0mass*dsqrt(s))
    gres=dsqrt(mres*mres*(mres*mres+gres0*gres0))
    kres=dsqrt(8d0)*mres*gres0*gres/(dacos(-1d0)*dsqrt(mres*mres+gres))
    fbw=pmass*kres/((x*x-mres*mres)**2+(mres*gres0)**2)
    pp_sig1=7.66d-3*eta**1.95d0*(1d0+eta+eta**5)*fbw**1.86d0
end function pp_sig1

real(8) function pp_sig2(tp)
    real(8), intent(in) :: tp
    pp_sig2=0d0
    if (tp >= 0.56d0 .and. tp <= 2d0) pp_sig2=5.7d0/(1d0+dexp(-9.3d0*(tp-1.4d0)))
end function pp_sig2

real(8) function pp_siginel(tp)
    real(8), intent(in) :: tp
    real(8) :: lterm,cutoff

    lterm=dlog(tp/pp_tkin())
    cutoff=1d0-(pp_tkin()/tp)**1.9d0
    pp_siginel=(30.7d0-0.96d0*lterm+0.18d0*lterm*lterm)*cutoff**3
end function pp_siginel

real(8) function pp_mult(tp,model)
    real(8), intent(in) :: tp
    integer, intent(in) :: model
    real(8) :: xi,a1,a2,a3,a4,a5,q

    pp_mult=0d0
    if (tp <= 2d0) return
    if (tp < 5d0) then
        q=(tp-pp_tkin())/pmass
        pp_mult=-6d-3+0.237d0*q-0.023d0*q*q
        return
    end if
    if (model == 0 .and. tp > 100d0) then
        a1=5.436d0; a2=0.254d0; a3=0.072d0; a4=0.075d0; a5=0.166d0
    else if (model == 1 .and. tp > 100d0) then
        a1=0.908d0; a2=0.0009d0; a3=6.089d0; a4=0.176d0; a5=0.448d0
    else if (model == 3 .and. tp > 50d0) then
        a1=0.652d0; a2=0.0016d0; a3=0.488d0; a4=0.1928d0; a5=0.483d0
    else
        a1=0.728d0; a2=0.596d0; a3=0.491d0; a4=0.2503d0; a5=0.117d0
    end if
    xi=(tp-3d0)/pmass
    pp_mult=a1*xi**a4*(1d0+dexp(-a2*xi**a5))*(1d0-dexp(-a3*xi**0.25d0))
end function pp_mult

real(8) function pp_tkin()
    pp_tkin=2d0*p0mass+p0mass*p0mass/(2d0*pmass)
end function pp_tkin

real(8) function pp_epimax(tp)
    real(8), intent(in) :: tp
    real(8) :: s,gcm,ecm,pcm,bcm

    s=2d0*pmass*(tp+2d0*pmass)
    gcm=(tp+2d0*pmass)/dsqrt(s)
    ecm=(s-4d0*pmass*pmass+p0mass*p0mass)/(2d0*dsqrt(s))
    pcm=dsqrt(ecm*ecm-p0mass*p0mass)
    bcm=dsqrt(1d0-gcm**(-2))
    pp_epimax=gcm*(ecm+pcm*bcm)
end function pp_epimax

real(8) function pp_egmax(tp)
    real(8), intent(in) :: tp
    real(8) :: gpi

    gpi=pp_epimax(tp)/p0mass
    pp_egmax=0.5d0*p0mass*gpi*(1d0+dsqrt(1d0-gpi**(-2)))
end function pp_egmax

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
