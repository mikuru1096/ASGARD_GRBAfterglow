module hadronic_formal
    use hadronic_shell
    implicit none
    private
    public :: formal_transport
contains

! Formal 1D hadronic 壳层序列：Python 只传数组，Fortran 负责物理推进顺序。
! Formal 1D hadronic shell sequence: Python passes arrays; Fortran owns the shell transport order.
subroutine formal_transport(tobs,gbulk,radius,bfield,nuseed,seed,ge, &
        einj,ppdens,pidx,eta,synidx,do_psyn, &
        do_pg,do_nu,do_bh,do_hic,do_pp, &
        ppmodel,qsyn,nth,nnu,nr,ne,ngp,nout,gp,gsec,dnp, &
        psyn,ssyn,pgam,nuout,pnu,pbh,sbh,dnebh, &
        qesrc,taubh,bhloss_ph,phic,dnn,dnpip, &
        dnpim,dnmuml,dnmumr,dnmupl, &
        dnmupr,ppionsyn,pmuonsyn,ppionic,pmuonic, &
        taupg,pgsurv,procpow)
    use constants
    use hadronic_base, only: build_grid, dyn_time
    use hadronic_bh, only: bh_calc
    use hadronic_decay, only: decay_hummer
    use hadronic_pg, only: pg_hummer
    use hadronic_rad, only: proton_syn
    use hadronic_transport, only: proton_inject
    implicit none
    integer, intent(in) :: synidx,do_psyn,do_pg,do_nu
    integer, intent(in) :: do_bh,do_hic,do_pp,ppmodel,qsyn,nth
    integer, intent(in) :: nnu,nr,ne,ngp,nout
    real(8), intent(in), dimension(nr) :: tobs,gbulk,radius,bfield
    real(8), intent(in), dimension(nnu) :: nuseed
    real(8), intent(in), dimension(nnu,nr) :: seed
    real(8), intent(in), dimension(ne) :: ge
    real(8), intent(in), dimension(nr) :: einj
    real(8), intent(in), dimension(nr) :: ppdens
    real(8), intent(in) :: pidx,eta
    real(8), intent(out), dimension(ngp) :: gp,gsec
    real(8), intent(out), dimension(ngp,nr) :: dnp
    real(8), intent(out), dimension(nnu,nr) :: psyn,ssyn,pgam
    real(8), intent(out), dimension(nout) :: nuout
    real(8), intent(out), dimension(nout,nr) :: pnu
    real(8), intent(out), dimension(nnu,nr) :: pbh
    real(8), intent(out), dimension(nnu,nr) :: sbh
    real(8), intent(out), dimension(ne,nr) :: dnebh
    real(8), intent(out), dimension(ne,nr) :: qesrc
    real(8), intent(out), dimension(nnu,nr) :: taubh
    real(8), intent(out), dimension(nnu,nr) :: bhloss_ph,phic
    real(8), intent(out), dimension(ngp,nr) :: dnn,dnpip
    real(8), intent(out), dimension(ngp,nr) :: dnpim,dnmuml
    real(8), intent(out), dimension(ngp,nr) :: dnmumr,dnmupl
    real(8), intent(out), dimension(ngp,nr) :: dnmupr
    real(8), intent(out), dimension(nnu,nr) :: ppionsyn
    real(8), intent(out), dimension(nnu,nr) :: pmuonsyn,ppionic
    real(8), intent(out), dimension(nnu,nr) :: pmuonic,taupg
    real(8), intent(out), dimension(nnu,nr) :: pgsurv
    real(8), intent(out), dimension(3,ngp,nr) :: procpow
    integer :: ir,inu,nalign
    real(8) :: dt,tdyn,dr,gmin,gmax,vol,divrate,dln
    real(8), dimension(nr) :: vols
    real(8), dimension(ngp) :: ehad,zrate
    real(8), dimension(nnu) :: eph
    real(8), dimension(ne) :: ee
    real(8), dimension(nout) :: enu
    real(8), dimension(ngp) :: dnprev,dnnext,dntry,qinj
    real(8), dimension(ngp) :: pdens,ndens
    real(8), dimension(nnu) :: phtau
    real(8), dimension(nnu) :: phden
    real(8), dimension(ngp) :: qpi0,qpip,qpim
    real(8), dimension(ngp) :: preinj,nreinj,pgloss,nloss
    real(8), dimension(nnu) :: pgloss_ph,qgam
    real(8), dimension(3,nnu) :: qproc
    real(8), dimension(nout) :: qnu
    real(8), dimension(ngp) :: qmupr,qmupl,qmuml,qmumr
    real(8), dimension(nout) :: qprompt,qmunu
    real(8), dimension(nnu) :: qmue
    real(8), dimension(ne) :: qbh
    real(8), dimension(ngp) :: bhfrac,bhloss
    real(8), dimension(nnu) :: qppg
    real(8), dimension(nout) :: qppnu
    real(8), dimension(ne) :: qpp
    real(8), dimension(ngp) :: pploss
    real(8), dimension(ngp) :: ppploss
    real(8), dimension(nnu) :: lppg
    real(8), dimension(nout) :: lppnu
    real(8), dimension(ne,nr) :: qse
    real(8), dimension(ngp) :: nprev,nnext,pipp,pipn
    real(8), dimension(ngp) :: pimp,pimn,mumlp,mumln
    real(8), dimension(ngp) :: mumrp,mumrn,muplp,mupln
    real(8), dimension(ngp) :: muprp,muprn

    ! 初始化质子/二级粒子/光子能格和输出缓存。
    ! Initialize proton, secondary-particle, photon grids, and output buffers.
    call scan_pmax(nr,radius,gbulk,bfield,eta,gmax)
    call build_grid(ngp,1d0+1d-3,gmax,gp)
    gsec=gp
    ehad=gp*Para_m_p_GeV
    ee=ge*Para_m_e_GeV
    eph=nuseed*Para_h_GeV
    zrate=0d0
    dnprev=0d0
    nprev=0d0; pipp=0d0; pimp=0d0
    mumlp=0d0; mumrp=0d0; muplp=0d0; muprp=0d0
    dnp=0d0; psyn=0d0; ssyn=0d0; pgam=0d0; pnu=0d0
    pbh=0d0; sbh=0d0; dnebh=0d0; qesrc=0d0
    taubh=0d0; bhloss_ph=0d0; phic=0d0; qse=0d0
    dnn=0d0; dnpip=0d0; dnpim=0d0
    dnmuml=0d0; dnmumr=0d0; dnmupl=0d0
    dnmupr=0d0; ppionsyn=0d0; pmuonsyn=0d0
    ppionic=0d0; pmuonic=0d0; taupg=0d0; pgsurv=1d0
    procpow=0d0

    if (nout > 1) then
        do inu=1,nout
            nuout(inu)=1d1**(dlog10(1d-3*Para_GeV2Hz)+dble(inu-1)* &
                           (dlog10(1d8*Para_GeV2Hz)-dlog10(1d-3*Para_GeV2Hz))/dble(nout-1))
        end do
    else
        nuout(1)=1d-3*Para_GeV2Hz
    end if
    enu=Para_h_GeV*nuout

    dln=dlog(ehad(2)/ehad(1))
    nalign=int(ceiling((dlog(eph(nnu))-dlog(eph(1)))/dln))+1
    call shell_volumes(nr,radius,vols)

    ! 按壳层推进：质子注入 -> pγ/BH/pp -> 衰变/二级辐射 -> 光子存活。
    ! Shell loop: proton injection -> p-gamma/BH/pp -> decay/secondary radiation -> photon survival.
    do ir=1,nr
        call shell_geom(nr,radius,gbulk,ir,dr,dt)
        tdyn=dyn_time(radius(ir),gbulk(ir))
        gmin=max(gp(1),gbulk(ir))
        vol=vols(ir)
        call proton_inject(ngp,gp,pidx,einj(ir),gmin,gp(ngp),qinj)
        call proton_step(ngp,gp,dnprev,qinj,bfield(ir), &
                                               tdyn,Para_m_p_GeV,qsyn,zrate,zrate,zrate, &
                                               zrate,vol,dt,dntry)
        call shell_density(ngp,dntry,Para_m_p_GeV,vol,pdens)
        call shell_density(ngp,nprev,Para_m_n_GeV,vol,ndens)

        ! pγ 先估计光子吸收，再用 survival 后的 photon field 重算正式源项。
        ! p-gamma is evaluated before and after photon survival to align loss and source terms.
        eph(1:nnu)=Para_h_GeV*nuseed(1:nnu)
        phtau(1:nnu)=seed(1:nnu,ir)/Para_h_GeV
        call pg_hummer(ngp,nnu,ehad,pdens, &
                                             eph,phtau,ndens,qpi0, &
                                             qpip,qpim,preinj,nreinj,pgloss, &
                                             nloss,pgloss_ph)
        call decay_hummer(ngp,ehad,qpi0,qpip,qpim, &
                                                nnu,eph,nout,enu,nnu, &
                                                eph,qgam,qproc,qmupr,qmupl, &
                                                qmuml,qmumr,qprompt,qmunu, &
                                                qmue,qnu)
        call photon_loss(nnu,nr,radius,gbulk,ir,pgloss_ph, &
                                             taupg(:,ir),pgsurv(:,ir))
        phden=phtau*pgsurv(:,ir)
        call pg_hummer(ngp,nnu,ehad,pdens, &
                                             eph,phden,ndens,qpi0, &
                                             qpip,qpim,preinj,nreinj,pgloss, &
                                             nloss,pgloss_ph)
        call decay_hummer(ngp,ehad,qpi0,qpip,qpim, &
                                                nnu,eph,nout,enu,nnu, &
                                                eph,qgam,qproc,qmupr,qmupl, &
                                                qmuml,qmumr,qprompt,qmunu, &
                                                qmue,qnu)

        ! BH 和 pp 只在开关打开时加入 proton loss 与二级 e±/gamma/neutrino 源项。
        ! BH and pp add proton losses plus secondary e±/gamma/neutrino sources only when enabled.
        bhloss=0d0; qbh=0d0; bhloss_ph(:,ir)=0d0
        if (do_bh /= 0) then
            call bh_calc(ngp,ehad,pdens,nnu, &
                                                 eph,phden,ne,ee, &
                                                 qbh,bhfrac,bhloss_ph(:,ir))
            if (any(bhfrac > 0d0)) error stop "Bethe-Heitler fractional loss rate must be non-positive."
            bhloss=-gp*bhfrac
            call photon_loss(nnu,nr,radius,gbulk,ir,bhloss_ph(:,ir), &
                                                 taubh(:,ir),phtau)
        end if

        ppploss=0d0; qpp=0d0; lppg=0d0; lppnu=0d0
        if (do_pp /= 0) then
            call pp_delta(ngp,ehad,pdens,ppdens(ir), &
                                            nnu,eph,nout,enu,ne,ee, &
                                            0.5d0,0.17d0,1d0/3d0,ppmodel,qppg,qppnu,qpp, &
                                            pploss)
            if (any(pploss > 0d0)) error stop "pp proton loss rate must be non-positive."
            ppploss=-pploss
            call rate_lum(nnu,eph,qppg,vol,lppg)
            call rate_lum(nout,enu,qppnu,vol,lppnu)
        end if

        call proton_step(ngp,gp,dnprev,qinj,bfield(ir), &
                                               tdyn,Para_m_p_GeV,qsyn,bhloss,ppploss,pgloss, &
                                               preinj,vol,dt,dnnext)
        dnp(:,ir)=dnnext

        if (do_psyn /= 0) then
            call proton_syn(radius(ir),bfield(ir),ngp,nnu,gp,dnnext,nuseed, &
                                               psyn(:,ir),ssyn(:,ir))
        end if

        ! secondary species 与 pion/muon 辐射使用同一个壳层 photon field。
        ! Secondary species and pion/muon radiation use the same shell photon field.
        divrate=3d0/tdyn
        call species_step(ngp,ngp,gp,ehad,nreinj, &
                                                qpip,qpim,qmuml,qmumr,qmupl, &
                                                qmupr,nloss,dt,bfield(ir),divrate, &
                                                vol,nprev,pipp,pimp,mumlp,mumrp, &
                                                muplp,muprp,nnext,pipn,pimn,mumln, &
                                                mumrn,mupln,muprn)
        dnn(:,ir)=nnext; dnpip(:,ir)=pipn; dnpim(:,ir)=pimn
        dnmuml(:,ir)=mumln; dnmumr(:,ir)=mumrn
        dnmupl(:,ir)=mupln; dnmupr(:,ir)=muprn

        call secondary_project(ngp,nnu,nalign,ehad, &
                                                       eph,phden,pipn,pimn,mumln, &
                                                       mumrn,mupln,muprn,vol,bfield(ir), &
                                                       ppionsyn(:,ir),pmuonsyn(:,ir), &
                                                       ppionic(:,ir),pmuonic(:,ir))
        if (do_pg /= 0) then
            call rate_lum(nnu,eph,qgam,vol, &
                                                         pgam(:,ir))
            call proc_power(ngp,nnu,3,ehad,dnnext,eph, &
                                           qproc,vol,procpow(:,:,ir))
        end if
        if (do_pp /= 0) pgam(:,ir)=pgam(:,ir)+lppg
        if (do_nu /= 0) then
            call rate_lum(nout,enu,qnu,vol,pnu(:,ir))
            if (do_pp /= 0) pnu(:,ir)=pnu(:,ir)+lppnu
        end if
        if (do_hic /= 0) then
            call hic_project(ngp,nnu,nalign,ehad,eph, &
                                           phden,pdens,vol,phic(:,ir))
        end if
        if (do_bh /= 0 .or. do_pp /= 0) then
            call pair_content(ne,qpp,qbh,do_pp, &
                                                 do_bh,vol,qse(:,ir))
        end if

        dnprev=dnnext
        nprev=nnext; pipp=pipn; pimp=pimn
        mumlp=mumln; mumrp=mumrn; muplp=mupln; muprp=muprn
    end do

    ! BH/pp 产生的二级电子最后统一进入 electron sequence。
    ! Secondary electrons from BH/pp are advanced through the electron sequence at the end.
    if (do_bh /= 0 .or. do_pp /= 0) then
        call electron_seq(ne,nnu,nr,ge,radius,gbulk,bfield, &
                                                     nuseed,qse,synidx,nth,qsyn, &
                                                     dnebh,pbh,sbh,qesrc)
    end if
    psyn=psyn*pgsurv
    ssyn=ssyn*pgsurv
    pgam=pgam*pgsurv
    if (do_bh /= 0) then
        pbh=pbh*pgsurv
        sbh=sbh*pgsurv
    end if
    if (do_hic /= 0) phic=phic*pgsurv
    ppionsyn=ppionsyn*pgsurv
    pmuonsyn=pmuonsyn*pgsurv
    ppionic=ppionic*pgsurv
    pmuonic=pmuonic*pgsurv
end subroutine formal_transport

end module hadronic_formal
