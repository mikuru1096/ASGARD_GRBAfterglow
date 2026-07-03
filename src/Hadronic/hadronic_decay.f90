!f2py: skip
module hadronic_decay
    use constants
    use hadronic_base, only: check_grid
    implicit none
    private

    real(8), parameter :: pion_m = Para_m_pi_charged_GeV
    real(8), parameter :: muon_m = Para_m_mu_GeV
    real(8), parameter :: pion_tau = 2.6033d-8
    real(8), parameter :: muon_tau = 2.1969811d-6

    public :: pi0_gamma, pion_decay, muon_decay, decay_hummer

contains

! pi0 衰变为两个光子；把 pi0 source 谱投到 gamma 谱。
! Convert the neutral-pion source spectrum into a double-photon gamma source.
subroutine pi0_gamma(npi,epi,qpi0,ngam,egam,qgam)
    integer, intent(in) :: npi,ngam
    real(8), dimension(npi), intent(in) :: epi,qpi0
    real(8), dimension(ngam), intent(in) :: egam
    real(8), dimension(ngam), intent(out) :: qgam
    integer :: ip,ig
    real(8) :: dln,eout,epion,r,sumk

    qgam = 0d0
    dln = log_spacing(npi,epi)
    if (dln <= 0d0) return

    do ig=1,ngam
        eout = egam(ig)
        sumk = 0d0
        do ip=1,npi
            epion = epi(ip)
            r = 2d0*eout/epion
            if (r > 0d0 .and. r <= 1.00001d0) sumk = sumk + 2d0*r*epion*qpi0(ip)
        end do
        qgam(ig) = dln*sumk/eout
    end do
end subroutine pi0_gamma

! charged pion 衰变；输出 helicity-resolved muon 和 prompt neutrino source。
! Charged-pion decay into helicity-resolved muons and prompt neutrinos.
subroutine pion_decay(npi,epi,qpip,qpim,nmu,emu,qmupr,qmupl,qmuml,qmumr,nnu,enu,qnumu,qnumub)
    integer, intent(in) :: npi,nmu,nnu
    real(8), dimension(npi), intent(in) :: epi,qpip,qpim
    real(8), dimension(nmu), intent(in) :: emu
    real(8), dimension(nnu), intent(in) :: enu
    real(8), dimension(nmu), intent(out) :: qmupr,qmupl,qmuml,qmumr
    real(8), dimension(nnu), intent(out) :: qnumu,qnumub
    integer :: im,inu
    real(8) :: dln,rpi,em,en
    real(8) :: smupr,smupl,smuml,smumr,snup,snum
    real(8), dimension(npi) :: qpip_log,qpim_log,tdec

    qmupr = 0d0
    qmupl = 0d0
    qmuml = 0d0
    qmumr = 0d0
    qnumu = 0d0
    qnumub = 0d0

    dln = log_spacing(npi,epi)
    if (dln <= 0d0) return
    rpi = (muon_m/pion_m)**2

    call build_pion

    do im=1,nmu
        em = emu(im)
        call pion_madd(em,smupr,smupl,smuml,smumr)
        qmupr(im) = dln*smupr/em
        qmupl(im) = dln*smupl/em
        qmuml(im) = dln*smuml/em
        qmumr(im) = dln*smumr/em
    end do

    do inu=1,nnu
        en = enu(inu)
        call pion_nadd(en,snup,snum)
        qnumu(inu) = dln*snup/en
        qnumub(inu) = dln*snum/en
    end do

contains

    subroutine build_pion
        integer :: ip

        do ip=1,npi
            tdec(ip) = (epi(ip)/pion_m)*pion_tau
            qpip_log(ip) = epi(ip)*qpip(ip)/tdec(ip)
            qpim_log(ip) = epi(ip)*qpim(ip)/tdec(ip)
        end do
    end subroutine build_pion

    subroutine pion_madd(em,smupr,smupl,smuml,smumr)
        real(8), intent(in) :: em
        real(8), intent(out) :: smupr,smupl,smuml,smumr
        integer :: ip
        real(8) :: epion,r,f1,f2

        smupr = 0d0
        smupl = 0d0
        smuml = 0d0
        smumr = 0d0
        do ip=1,npi
            epion = epi(ip)
            r = em/epion
            if (r < 1d0 .and. r >= rpi) then
                f1 = rpi*(1d0-r)/((1d0-rpi)**2*r)
                f2 = (r-rpi)/((1d0-rpi)**2*r)
                smupr = smupr + r*f1*qpip_log(ip)
                smupl = smupl + r*f2*qpip_log(ip)
                smuml = smuml + r*f1*qpim_log(ip)
                smumr = smumr + r*f2*qpim_log(ip)
            end if
        end do
    end subroutine pion_madd

    subroutine pion_nadd(en,snup,snum)
        real(8), intent(in) :: en
        real(8), intent(out) :: snup,snum
        integer :: ip
        real(8) :: epion,r,fnu

        snup = 0d0
        snum = 0d0
        do ip=1,npi
            epion = epi(ip)
            r = en/epion
            if (r >= 0d0 .and. r <= 1d0-rpi) then
                fnu = 1d0/(1d0-rpi)
                snup = snup + r*fnu*qpip_log(ip)
                snum = snum + r*fnu*qpim_log(ip)
            end if
        end do
    end subroutine pion_nadd
end subroutine pion_decay

! muon 衰变；输出 electron/positron 与 flavor-resolved neutrino source。
! Muon decay into electron/positron and flavor-resolved neutrino sources.
subroutine muon_decay(nmu,emu,qmupr,qmupl,qmuml,qmumr,nnu,enu,qnue,qnueb,qnumu,qnumub,ne,ee,qem,qep)
    integer, intent(in) :: nmu,nnu,ne
    real(8), dimension(nmu), intent(in) :: emu,qmupr,qmupl,qmuml,qmumr
    real(8), dimension(nnu), intent(in) :: enu
    real(8), dimension(ne), intent(in) :: ee
    real(8), dimension(nnu), intent(out) :: qnue,qnueb,qnumu,qnumub
    real(8), dimension(ne), intent(out) :: qem,qep
    integer :: im,inu,ie
    real(8) :: dln,en,elec,snue,snueb,snumu,snumub,sem,sep
    real(8), dimension(nmu) :: qpl_log,qpr_log,qml_log,qmr_log,tdec

    qnue = 0d0
    qnueb = 0d0
    qnumu = 0d0
    qnumub = 0d0
    qem = 0d0
    qep = 0d0

    dln = log_spacing(nmu,emu)
    if (dln <= 0d0) return

    call build_muon

    do inu=1,nnu
        en = enu(inu)
        call muon_nadd(en,snue,snueb,snumu,snumub)
        qnue(inu) = dln*snue/en
        qnueb(inu) = dln*snueb/en
        qnumu(inu) = dln*snumu/en
        qnumub(inu) = dln*snumub/en
    end do

    do ie=1,ne
        elec = ee(ie)
        call muon_eadd(elec,sem,sep)
        qep(ie) = dln*sep/elec
        qem(ie) = dln*sem/elec
    end do

contains

    subroutine build_muon
        do im=1,nmu
            tdec(im) = (emu(im)/muon_m)*muon_tau
            qpl_log(im) = emu(im)*qmupl(im)/tdec(im)
            qpr_log(im) = emu(im)*qmupr(im)/tdec(im)
            qml_log(im) = emu(im)*qmuml(im)/tdec(im)
            qmr_log(im) = emu(im)*qmumr(im)/tdec(im)
        end do
    end subroutine build_muon

    subroutine muon_nadd(en,snue,snueb,snumu,snumub)
        real(8), intent(in) :: en
        real(8), intent(out) :: snue,snueb,snumu,snumub
        real(8) :: em,x,f1p,f1m,f2p,f2m

        snue = 0d0
        snueb = 0d0
        snumu = 0d0
        snumub = 0d0
        do im=1,nmu
            em = emu(im)
            x = en/em
            f1p = nu1_decay(x,1d0)
            f1m = nu1_decay(x,-1d0)
            f2p = nu2_decay(x,1d0)
            f2m = nu2_decay(x,-1d0)
            snueb = snueb + x*(f2p*qml_log(im) + f2m*qmr_log(im))
            snue = snue + x*(f2m*qpl_log(im) + f2p*qpr_log(im))
            snumub = snumub + x*(f1m*qpl_log(im) + f1p*qpr_log(im))
            snumu = snumu + x*(f1p*qml_log(im) + f1m*qmr_log(im))
        end do
    end subroutine muon_nadd

    subroutine muon_eadd(elec,sem,sep)
        real(8), intent(in) :: elec
        real(8), intent(out) :: sem,sep
        real(8) :: em,x

        sem = 0d0
        sep = 0d0
        do im=1,nmu
            em = emu(im)
            x = elec/em
            if (x > 0d0 .and. x <= 1.00001d0) then
                sep = sep + x*(4d0/3d0)*(1d0-x*x*x)*(qpl_log(im) + qpr_log(im))
                sem = sem + x*(4d0/3d0)*(1d0-x*x*x)*(qml_log(im) + qmr_log(im))
            end if
        end do
    end subroutine muon_eadd
end subroutine muon_decay

! Hummer 2010 衰变聚合：统一输出 gamma、prompt/muon neutrino 与 e± source。
! Hummer 2010 decay aggregate for gamma, prompt/muon neutrino, and e± sources.
subroutine decay_hummer(ngrid,ehad,qpi0,qpip,qpim,ngam,egam,nnu,enu,nproc,eproc, &
                        qgam,qproc,qmupr,qmupl,qmuml,qmumr,qprompt,qmunu,qmue,qnu)
    integer, intent(in) :: ngrid,ngam,nnu,nproc
    real(8), dimension(ngrid), intent(in) :: ehad,qpi0,qpip,qpim
    real(8), dimension(ngam), intent(in) :: egam
    real(8), dimension(nnu), intent(in) :: enu
    real(8), dimension(nproc), intent(in) :: eproc
    real(8), dimension(ngam), intent(out) :: qgam
    real(8), dimension(3,nproc), intent(out) :: qproc
    real(8), dimension(ngrid), intent(out) :: qmupr,qmupl,qmuml,qmumr
    real(8), dimension(nnu), intent(out) :: qprompt,qmunu,qnu
    real(8), dimension(nproc), intent(out) :: qmue
    integer :: iproc
    real(8), dimension(nproc) :: qpgam,qem,qep
    real(8), dimension(nnu) :: qnumu,qnumub,qnue,qnueb,qmnumu,qmnumub

    call pi0_gamma(ngrid,ehad,qpi0,ngam,egam,qgam)
    call pi0_gamma(ngrid,ehad,qpi0,nproc,eproc,qpgam)
    call pion_decay(ngrid,ehad,qpip,qpim,ngrid,ehad,qmupr,qmupl,qmuml,qmumr,nnu,enu,qnumu,qnumub)
    call muon_decay(ngrid,ehad,qmupr,qmupl,qmuml,qmumr,nnu,enu,qnue,qnueb,qmnumu,qmnumub,nproc,eproc,qem,qep)

    qprompt = qnumu + qnumub
    qmunu = qnue + qnueb + qmnumu + qmnumub
    qmue = qem + qep
    qnu = qprompt + qmunu

    qproc = 0d0
    qproc(1,:) = qpgam
    do iproc=1,nproc
        qproc(2,iproc) = log_interp(nnu,enu,qprompt,eproc(iproc))
        qproc(3,iproc) = log_interp(nnu,enu,qmunu,eproc(iproc)) + qmue(iproc)
    end do
end subroutine decay_hummer

! 返回对数网格间距；网格点不足时返回 0d0。
! Return the logarithmic grid spacing, or 0d0 when the grid is too short.
real(8) function log_spacing(ngrid,grid)
    integer, intent(in) :: ngrid
    real(8), dimension(ngrid), intent(in) :: grid
    real(8) :: dln

    if (ngrid <= 1) then
        log_spacing = 0d0
    else
        call check_grid(ngrid,grid,"energy_gev",dln)
        log_spacing = dln
    end if
end function log_spacing

! muon decay neutrino 谱函数 f_nu1(x,h)。
! Muon-decay neutrino spectrum function f_nu1(x,h).
real(8) function nu1_decay(x,h)
    real(8), intent(in) :: x,h

    if (x >= 0d0 .and. x <= 1d0) then
        nu1_decay = (5d0/3d0 - 3d0*x*x + 4d0*x*x*x/3d0) + &
                    h*(-1d0/3d0 + 3d0*x*x - 8d0*x*x*x/3d0)
    else
        nu1_decay = 0d0
    end if
end function nu1_decay

! muon decay neutrino 谱函数 f_nu2(x,h)。
! Muon-decay neutrino spectrum function f_nu2(x,h).
real(8) function nu2_decay(x,h)
    real(8), intent(in) :: x,h

    if (x >= 0d0 .and. x <= 1d0) then
        nu2_decay = (2d0 - 6d0*x*x + 4d0*x*x*x) + &
                    h*(2d0 - 12d0*x + 18d0*x*x - 8d0*x*x*x)
    else
        nu2_decay = 0d0
    end if
end function nu2_decay

! 对数能量网格上的线性插值；范围外返回 0d0。
! Linear interpolation on a logarithmic energy grid, returning 0d0 outside range.
real(8) function log_interp(ngrid,grid,rate,eval)
    integer, intent(in) :: ngrid
    real(8), dimension(ngrid), intent(in) :: grid,rate
    real(8), intent(in) :: eval
    integer :: ilo,ihi,imid
    real(8) :: w

    if (ngrid <= 0) then
        log_interp = 0d0
        return
    end if
    if (eval < grid(1) .or. eval > grid(ngrid)) then
        log_interp = 0d0
        return
    end if
    if (ngrid == 1) then
        log_interp = rate(1)
        return
    end if

    ilo = 1
    ihi = ngrid
    do while (ihi - ilo > 1)
        imid = (ilo + ihi)/2
        if (grid(imid) <= eval) then
            ilo = imid
        else
            ihi = imid
        end if
    end do

    if (grid(ihi) <= grid(ilo)) then
        log_interp = rate(ilo)
        return
    end if

    w = dlog(eval/grid(ilo))/dlog(grid(ihi)/grid(ilo))
    log_interp = (1d0-w)*rate(ilo) + w*rate(ihi)
end function log_interp

end module hadronic_decay
