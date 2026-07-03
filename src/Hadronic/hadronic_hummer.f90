!f2py: skip
module hadronic_hummer
    use constants
    use hadronic_base
    use hadronic_pg, only: pg_hummer
    use hadronic_decay, only: decay_hummer
    use hadronic_species, only: species_advance
    implicit none
    private

    public :: hummer_shell

contains

! 旧 Hummer p-gamma 壳层路径：先估计 photon survival，再推进 proton/neutron/pi/mu。
! Legacy Hummer p-gamma shell path: estimate photon survival, then advance proton/neutron/pi/mu.
subroutine hummer_shell(ngp,nph,nnu,dt,rad,gbulk,bfield,vol, &
                        gam,dn0,nuph,seed,nunu, &
                        incpg,incnu,n0,pip0,pim0,muml0,mumr0, &
                        mupl0,mupr0,dn1,n1,pip1,pim1,muml1,mumr1, &
                        mupl1,mupr1,pgsurv,pglum,nulum)
    implicit none
    integer, intent(in) :: ngp,nph,nnu,incpg,incnu
    real(8), intent(in) :: dt,rad,gbulk,bfield,vol
    real(8), intent(in), dimension(ngp) :: gam,dn0
    real(8), intent(in), dimension(nph) :: nuph,seed
    real(8), intent(in), dimension(nnu) :: nunu
    real(8), intent(in), dimension(ngp) :: n0,pip0,pim0
    real(8), intent(in), dimension(ngp) :: muml0,mumr0,mupl0,mupr0
    real(8), intent(out), dimension(ngp) :: dn1,n1,pip1,pim1
    real(8), intent(out), dimension(ngp) :: muml1,mumr1,mupl1,mupr1
    real(8), intent(out), dimension(nph) :: pgsurv,pglum
    real(8), intent(out), dimension(nnu) :: nulum

    dn1=dn0
    n1=n0; pip1=pip0; pim1=pim0
    muml1=muml0; mumr1=mumr0; mupl1=mupl0; mupr1=mupr0
    pgsurv=1d0
    pglum=0d0
    nulum=0d0

    if (incpg == 0 .and. incnu == 0) return

    block
        real(8), dimension(ngp) :: ehad,pden,nden
        real(8), dimension(nph) :: eph,phtau,phatt
        real(8), dimension(nnu) :: enu
        real(8), dimension(ngp) :: qpi0,qpip,qpim,preinj,nreinj
        real(8), dimension(ngp) :: ploss,nloss
        real(8), dimension(nph) :: phloss,taupg,grate
        real(8), dimension(3,nph) :: proc
        real(8), dimension(ngp) :: mupr,mupl,mumnl,mumnr
        real(8), dimension(nnu) :: nuprompt,numuon
        real(8), dimension(nph) :: emuon
        real(8) :: tpath,rdt,tau
        real(8), dimension(ngp) :: nsrc,pipsrc,pimsrc
        real(8), dimension(ngp) :: mumlsrc,mumrsrc,muplsrc,muprsrc
        real(8) :: divr
        integer :: i

        ehad=gam*proton_m
        pden=dn0/(vol*proton_m)
        nden=n0/(vol*neutron_m)
        eph=nuph*Para_h_GeV
        phtau=seed/Para_h_GeV
        enu=nunu*Para_h_GeV

        ! 第一次 p-gamma 卷积只用来估计 photon escape/survival。
        ! First p-gamma convolution estimates photon escape/survival.
        call pg_hummer(ngp,nph,ehad,pden,eph,phtau,nden,qpi0,qpip,qpim,preinj,nreinj, &
                       ploss,nloss,phloss)
        tpath=rad/(12d0*gbulk*Para_c)
        taupg=phloss*tpath
        do i=1,nph
            tau=taupg(i)
            if (tau > 1d-6) then
                pgsurv(i)=(1d0-exp(-tau))/tau
            else if (tau > 0d0) then
                pgsurv(i)=1d0-tau/2d0+tau*tau/6d0
            end if
        end do

        phatt=phtau*pgsurv
        ! 第二次 p-gamma 卷积使用 attenuated photon field，输出真实源项。
        ! Second p-gamma convolution uses the attenuated photon field for sources.
        call pg_hummer(ngp,nph,ehad,pden,eph,phatt,nden,qpi0,qpip,qpim,preinj,nreinj, &
                       ploss,nloss,phloss)

        do i=1,ngp
            rdt=ploss(i)*dt
            if (ploss(i) > 0d0) then
                dn1(i)=dn0(i)*exp(-rdt)+(1d0-exp(-rdt))/ploss(i) &
                           *vol*proton_m*preinj(i)
            else
                dn1(i)=dn0(i)+dt*vol*proton_m*preinj(i)
            end if
        end do

        call decay_hummer(ngp,ehad,qpi0,qpip,qpim,nph,eph,nnu,enu,nph,eph,grate,proc, &
                          mupr,mupl,mumnl,mumnr,nuprompt,numuon,emuon,nulum)
        call map_secondaries(ngp,gam,ehad,vol,nreinj,qpip,qpim, &
                             mumnl,mumnr,mupl,mupr,nsrc,pipsrc,pimsrc, &
                             mumlsrc,mumrsrc,muplsrc,muprsrc)
        divr=3d0*gbulk*Para_c/rad
        call species_advance(ngp,gam,dt,bfield,divr, &
                             n0,pip0,pim0,muml0,mumr0,mupl0,mupr0, &
                             nsrc,pipsrc,pimsrc,mumlsrc,mumrsrc,muplsrc,muprsrc, &
                             n1,pip1,pim1,muml1,mumr1,mupl1,mupr1)
        call neutron_loss(ngp,gam,dt,ehad,nloss,n1)
        if (incpg /= 0) pglum=vol*grate*eph*Para_h_GeV*Para_GeV2erg*pgsurv
        if (incnu /= 0) nulum=vol*nulum*enu*Para_h_GeV*Para_GeV2erg
    end block
end subroutine hummer_shell

! p-gamma secondary source 从 proton gamma 网格映射到各 species 网格。
! Map p-gamma secondary sources from the proton-gamma grid to each species grid.
subroutine map_secondaries(ngp,gam,ehad,vol,nreinj,qpip,qpim, &
                           mumnl,mumnr,mupl,mupr,nsrc,pipsrc,pimsrc, &
                           mumlsrc,mumrsrc,muplsrc,muprsrc)
    implicit none
    integer, intent(in) :: ngp
    real(8), intent(in), dimension(ngp) :: gam,ehad
    real(8), intent(in) :: vol
    real(8), intent(in), dimension(ngp) :: nreinj,qpip,qpim
    real(8), intent(in), dimension(ngp) :: mumnl,mumnr,mupl,mupr
    real(8), intent(out), dimension(ngp) :: nsrc,pipsrc,pimsrc
    real(8), intent(out), dimension(ngp) :: mumlsrc,mumrsrc,muplsrc,muprsrc
    real(8), dimension(ngp) :: en,epi,emu
    real(8), dimension(ngp) :: srcgev

    en=gam*neutron_m
    epi=gam*pion_m
    emu=gam*muon_m

    call pos_interp(ngp,ehad,nreinj,en,srcgev)
    nsrc=vol*srcgev*neutron_m
    call pos_interp(ngp,ehad,qpip,epi,srcgev)
    pipsrc=vol*srcgev*pion_m
    call pos_interp(ngp,ehad,qpim,epi,srcgev)
    pimsrc=vol*srcgev*pion_m
    call pos_interp(ngp,ehad,mumnl,emu,srcgev)
    mumlsrc=vol*srcgev*muon_m
    call pos_interp(ngp,ehad,mumnr,emu,srcgev)
    mumrsrc=vol*srcgev*muon_m
    call pos_interp(ngp,ehad,mupl,emu,srcgev)
    muplsrc=vol*srcgev*muon_m
    call pos_interp(ngp,ehad,mupr,emu,srcgev)
    muprsrc=vol*srcgev*muon_m
end subroutine map_secondaries

! neutron 的 p-gamma loss 按同一 gamma 网格插值后显式扣除。
! Apply neutron p-gamma loss after interpolation to the shared gamma grid.
subroutine neutron_loss(ngp,gam,dt,ehad,nloss,n1)
    implicit none
    integer, intent(in) :: ngp
    real(8), intent(in), dimension(ngp) :: gam,ehad,nloss
    real(8), intent(in) :: dt
    real(8), intent(inout), dimension(ngp) :: n1
    real(8), dimension(ngp) :: en,nlossi
    integer :: i

    en=gam*neutron_m
    call pos_interp(ngp,ehad,nloss,en,nlossi)
    do i=1,ngp
        if (1d0-dt*nlossi(i) < 0d0) error stop "hadronic p-gamma neutron survival became negative."
        n1(i)=n1(i)*(1d0-dt*nlossi(i))
    end do
end subroutine neutron_loss

! 正值 log-log 插值：只在输入两端都为正时插值。
! Positive log-log interpolation; interpolate only when both bracketing values are positive.
subroutine pos_interp(ngrid,xsrc,ysrc,xdst,ydst)
    implicit none
    integer, intent(in) :: ngrid
    real(8), intent(in), dimension(ngrid) :: xsrc,ysrc,xdst
    real(8), intent(out), dimension(ngrid) :: ydst
    integer :: i,j
    real(8) :: lx,frac

    ydst=0d0
    do i=1,ngrid
        if (xdst(i) < xsrc(1) .or. xdst(i) > xsrc(ngrid)) cycle
        do j=1,ngrid-1
            if (xdst(i) >= xsrc(j) .and. xdst(i) <= xsrc(j+1)) then
                if (ysrc(j) > 0d0 .and. ysrc(j+1) > 0d0) then
                    lx=log(xdst(i)/xsrc(j))/log(xsrc(j+1)/xsrc(j))
                    frac=log(ysrc(j+1)/ysrc(j))*lx
                    ydst(i)=ysrc(j)*exp(frac)
                end if
                exit
            end if
        end do
    end do
end subroutine pos_interp

end module hadronic_hummer
