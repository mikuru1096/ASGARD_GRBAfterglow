! 电磁对级联核心：gamma-gamma 产对、pair cooling、同步辐射和壳层光子闭合。
! Electromagnetic pair-cascade core: gamma-gamma pair production, pair cooling, synchrotron emission, and shell photon closure.
module hadronic_cascade
    use constants
    use hadronic_base, only: electron_m, check_grid
    use hadronic_pair, only: pair_operator
    use electron_radiation_kernel, only: syn_state
    implicit none
    private

    public :: cascade_step, cascade_seq

contains

! 单步级联：gamma-gamma 产对 -> 同步冷却 -> pair synch photon source。
! Single cascade step: gamma-gamma pair production -> synchrotron cooling -> pair-synch photon source.
subroutine cascade_step(nph,eph,phden, &
                                  ne,ee,bfield,ptime, &
                                  cspec,phloss,pabs)
    integer, intent(in) :: nph,ne
    real(8), intent(in), dimension(nph) :: eph,phden
    real(8), intent(in), dimension(ne) :: ee
    real(8), intent(in) :: bfield,ptime
    real(8), intent(out), dimension(nph) :: cspec,phloss
    real(8), intent(out) :: pabs

    integer :: comfac,nlep,np
    real(8), dimension(ne) :: ge
    real(8) :: dln,csyn,ginj,gcool
    real(8), dimension(ne) :: qpair,qtot
    real(8), dimension(nph) :: ploss
    real(8), dimension(ne) :: pairs
    real(8) :: pinj,nuc,x,fx,temp

    nlep = ne
    np = nph
    comfac = 138  ! 与 AM3 一致。 / Match AM3 center-of-momentum sampling.
    call validate_step

    call produce_pairs

    call build_egamma
    if (bfield <= 0d0) then
        cspec(1:np) = 0d0
        return
    end if

    call cool_pairs

    call emit_syn

contains

    subroutine validate_step
        if (bfield < 0d0) error stop "pair cascade requires bfield >= 0."
        if (ptime < 0d0) error stop "pair cascade requires ptime >= 0."
        if (any(phden(1:np) < 0d0)) error stop "pair cascade requires non-negative photon density."
        call check_grid(np,eph,"pair cascade eph")
        call check_grid(nlep,ee,"pair cascade ee",dln)
    end subroutine validate_step

    subroutine produce_pairs
        call pair_operator(np,eph,phden, &
                                                nlep,ee, &
                                                comfac,ploss, &
                                                qpair, &
                                                qtot,pabs, &
                                                pinj)
        phloss(1:np) = ploss(1:np)
    end subroutine produce_pairs

    subroutine build_egamma
        ge(1:nlep) = ee(1:nlep) / electron_m
        if (ge(1) < 1d0) error stop "pair cascade electron grid must start at gamma >= 1."
    end subroutine build_egamma

    subroutine cool_pairs
        integer :: ie

        pairs(1:nlep) = 0d0
        csyn = Para_SigmaT * bfield**2 / (6d0*pi*Para_m_e*Para_c)

        do ie=2,nlep-1
            if (qtot(ie) <= 0d0) cycle
            ginj = ge(ie)
            ! 冷却后 Lorentz 因子: 1/gamma_c = 1/gamma_inj + csyn * ptime。
            ! Cooled Lorentz factor after synchrotron loss over the path time.
            gcool = ginj / (1d0 + csyn * ginj * ptime)
            call cool_deposit(nlep,ge,dln,ginj,gcool, &
                                          qtot(ie),csyn,pairs)
        end do
    end subroutine cool_pairs

    subroutine emit_syn
        integer :: ie,iph

        cspec(1:np) = 0d0
        temp = dsqrt(3d0)*Para_e**3 / Para_m_p_e

        do iph=1,np
            do ie=2,nlep-1
                if (pairs(ie) <= 0d0) cycle
                nuc = 4.2d6 * bfield * ge(ie)*ge(ie)
                x = eph(iph)*Para_GeV2Hz / nuc
                if (x > 1d2) cycle
                fx = 1.80842d0 * x**(1d0/3d0) * dexp(-x)
                cspec(iph) = cspec(iph) + temp*bfield*pairs(ie)*fx*ge(ie)*dln
            end do
        end do
    end subroutine emit_syn
end subroutine cascade_step

! 壳层序列级联：把 pair population 和 photon density 沿半径序列推进。
! Shell-sequence cascade: advect pair population and photon density through the radial shell sequence.
subroutine cascade_seq(nph,ne,nshell,eph,phprim, &
                                     ee,freq,radius,gbulk,tobs, &
                                     bfield,nth,synidx,nsub, &
                                     phloss,tau_pair,pden,plum, &
                                     pseed,cphden,pabs,pinjout)
    integer, intent(in) :: nph,ne,nshell,nth,synidx,nsub
    real(8), intent(in), dimension(nph) :: eph
    real(8), intent(in), dimension(nph,nshell) :: phprim
    real(8), intent(in), dimension(ne) :: ee
    real(8), intent(in), dimension(nph) :: freq
    real(8), intent(in), dimension(nshell) :: radius,gbulk,tobs,bfield
    real(8), intent(out), dimension(nph,nshell) :: phloss,tau_pair
    real(8), intent(out), dimension(ne,nshell) :: pden
    real(8), intent(out), dimension(nph,nshell) :: plum
    real(8), intent(out), dimension(nph,nshell) :: pseed,cphden
    real(8), intent(out), dimension(nshell) :: pabs,pinjout

    integer :: ir,isub,comfac
    real(8), dimension(ne) :: ge,prev,cur
    real(8), dimension(ne+1) :: edge
    real(8), dimension(nph) :: phden,src,loss,qsyn
    real(8), dimension(ne) :: qpair,qtot,psrc,eloss
    real(8), dimension(ne) :: next
    real(8), dimension(nph) :: phnext
    real(8) :: dt,dtsub,tesc,vol
    real(8) :: absloc,injloc

    comfac=138
    call validate_seq
    call build_grids

    phloss=0d0; tau_pair=0d0; pden=0d0
    plum=0d0; pseed=0d0; cphden=0d0
    pabs=0d0; pinjout=0d0; prev=0d0

    do ir=1,nshell
        call shell_geom(ir,dt,tesc,vol)
        dtsub=dt/dble(nsub)
        phden(1:nph)=phprim(1:nph,ir)
        cur(1:ne)=prev(1:ne)

        do isub=1,nsub
            call pair_operator(nph,eph,phden, &
                                                   ne,ee,comfac,loss, &
                                                   qpair,qtot, &
                                                   absloc,injloc)
            psrc(1:ne)=vol*qtot(1:ne)*electron_m*dtsub
            call electron_loss(ne,ge,bfield(ir), &
                                     radius(ir)/(gbulk(ir)*Para_c),eloss)
            call advance_energy(ne,ge,edge,cur,psrc,eloss,dtsub,next)
            cur(1:ne)=next(1:ne)

            call pair_syn(ir,cur,plum(:,ir),pseed(:,ir))
            qsyn(1:nph)=pseed(1:nph,ir)/Para_h_GeV/tesc
            src(1:nph)=phprim(1:nph,ir)/tesc+qsyn(1:nph)
            call advance_photon(nph,phden,src,loss+1d0/tesc,dtsub,phnext)
            phden(1:nph)=phnext(1:nph)
            pabs(ir)=absloc
            pinjout(ir)=injloc
        end do

        call pair_operator(nph,eph,phden, &
                                               ne,ee,comfac,loss, &
                                               qpair,qtot,absloc,injloc)
        call pair_syn(ir,cur,plum(:,ir),pseed(:,ir))
        phloss(1:nph,ir)=loss(1:nph)
        tau_pair(1:nph,ir)=loss(1:nph)*tesc
        pden(1:ne,ir)=cur(1:ne)
        cphden(1:nph,ir)=pseed(1:nph,ir)/Para_h_GeV
        prev(1:ne)=cur(1:ne)
    end do

contains

    subroutine validate_seq
        if (nph < 2 .or. ne < 2 .or. nshell < 1) error stop "pair cascade sequence grid is too small."
        if (nsub < 1) error stop "pair cascade sequence requires positive substeps."
        if (any(phprim < 0d0)) error stop "pair cascade sequence requires non-negative photons."
        if (any(gbulk < 1d0)) error stop "pair cascade sequence requires gbulk >= 1."
        if (any(radius <= 0d0) .or. any(bfield < 0d0)) error stop "pair cascade sequence received bad shell state."
        call check_grid(nph,eph,"pair cascade eph")
        call check_grid(ne,ee,"pair cascade ee")
        call check_grid(nph,freq,"pair cascade freq")
    end subroutine validate_seq

    subroutine build_grids
        integer :: ie

        ge(1:ne)=ee(1:ne)/electron_m
        if (ge(1) < 1d0) error stop "pair cascade sequence electron gamma grid must start at gamma >= 1."
        edge(1)=ge(1)*dsqrt(ge(1)/ge(2))
        do ie=2,ne
            edge(ie)=dsqrt(ge(ie-1)*ge(ie))
        end do
        edge(ne+1)=ge(ne)*dsqrt(ge(ne)/ge(ne-1))
    end subroutine build_grids

    subroutine shell_geom(ir0,dtout,tescout,volout)
        integer, intent(in) :: ir0
        real(8), intent(out) :: dtout,tescout,volout

        if (ir0 == 1) then
            dtout=tobs(1)
            volout=(4d0/3d0)*pi*radius(1)**3
        else
            dtout=tobs(ir0)-tobs(ir0-1)
            volout=(4d0/3d0)*pi*(radius(ir0)**3-radius(ir0-1)**3)
        end if
        if (dtout <= 0d0) error stop "pair cascade sequence times must be strictly increasing."
        tescout=radius(ir0)/(12d0*gbulk(ir0)*Para_c)
    end subroutine shell_geom

    subroutine pair_syn(ir0,pstate,psyn,seed)
        integer, intent(in) :: ir0
        real(8), intent(in), dimension(ne) :: pstate
        real(8), intent(out), dimension(nph) :: psyn,seed
        real(8), dimension(nph) :: pemit,tausyn

        if (bfield(ir0) <= 0d0) then
            psyn=0d0; seed=0d0
            return
        end if
        call syn_state(synidx,radius(ir0),bfield(ir0), &
                                    ne,nph,nth,ge,pstate,freq,pemit,psyn, &
                                    seed,tausyn)
    end subroutine pair_syn
end subroutine cascade_seq

! ------------------------------------------------------------
! 将注入 pair 按同步冷却轨迹沉积到 gamma grid。
! Deposit injected pairs along the synchrotron-cooling track on the gamma grid.
! 物理权重: N(gamma) dgamma = injected number * cooling time / dln_gamma。
! Physical weight: N(gamma) dgamma = injected number * cooling time / dln_gamma.
subroutine cool_deposit(nlep,ge,dln,ginj,gcool,inj,csyn,pairs)
    integer, intent(in) :: nlep
    real(8), intent(in), dimension(nlep) :: ge
    real(8), intent(in) :: dln,ginj,gcool,inj,csyn
    real(8), intent(inout), dimension(nlep) :: pairs
    integer :: k
    real(8) :: weight

    do k=1,nlep
        if (ge(k) > gcool .and. ge(k) <= ginj) then
            weight = 1d0 / (csyn * ge(k)*ge(k))  ! t_cool(gamma)
            pairs(k) = pairs(k) + inj * weight * dln
        end if
    end do
    ! 不需要额外归一化，weight 已经直接是 dN/dgamma 的贡献。
    ! No extra normalization is applied: weight is already the dN/dgamma contribution.

end subroutine cool_deposit

! pair 能量损失率：同步冷却加动力学逃逸。
! Pair energy-loss rate: synchrotron cooling plus dynamical escape.
subroutine electron_loss(ne,ge,bfield,tdyn,loss)
    integer, intent(in) :: ne
    real(8), intent(in), dimension(ne) :: ge
    real(8), intent(in) :: bfield,tdyn
    real(8), intent(out), dimension(ne) :: loss
    real(8) :: csyn

    csyn=Para_SigmaT*bfield*bfield/(6d0*pi*Para_m_e*Para_c)
    loss(1:ne)=csyn*ge(1:ne)*ge(1:ne)+ge(1:ne)/tdyn
end subroutine electron_loss

! 沿冷却轨迹推进 pair content，再加本步注入。
! Advance pair content along cooling trajectories, then add the step source.
subroutine advance_energy(ne,ge,edge,prev,source,loss,dt,next)
    integer, intent(in) :: ne
    real(8), intent(in), dimension(ne) :: ge,prev,source
    real(8), intent(in), dimension(ne+1) :: edge
    real(8), intent(in), dimension(ne) :: loss
    real(8), intent(in) :: dt
    real(8), intent(out), dimension(ne) :: next
    integer :: ie,target
    real(8), dimension(ne) :: content,nextcon
    real(8) :: gcool

    content(1:ne)=prev(1:ne)*(edge(2:ne+1)-edge(1:ne))
    nextcon=0d0
    do ie=1,ne
        gcool=ge(ie)-loss(ie)*dt
        target=gamma_bin(ne,edge,gcool)
        if (target >= 1 .and. target <= ne) nextcon(target)=nextcon(target)+content(ie)
    end do
    next(1:ne)=nextcon(1:ne)/(edge(2:ne+1)-edge(1:ne))+source(1:ne)
end subroutine advance_energy

! 定位冷却后的 Lorentz factor 所在 cell。
! Locate the cell containing the cooled Lorentz factor.
integer function gamma_bin(ne,edge,gval)
    integer, intent(in) :: ne
    real(8), intent(in), dimension(ne+1) :: edge
    real(8), intent(in) :: gval
    integer :: ie

    gamma_bin=0
    do ie=1,ne
        if (gval >= edge(ie) .and. gval < edge(ie+1)) then
            gamma_bin=ie
            return
        end if
    end do
end function gamma_bin

! 光子场解析推进：源项常数，吸收率在子步内常数。
! Analytic photon update with constant source and loss over a substep.
subroutine advance_photon(nph,prev,src,loss,dt,next)
    integer, intent(in) :: nph
    real(8), intent(in), dimension(nph) :: prev,src,loss
    real(8), intent(in) :: dt
    real(8), intent(out), dimension(nph) :: next
    integer :: iph
    real(8) :: atten

    do iph=1,nph
        if (loss(iph) > 0d0) then
            atten=dexp(-loss(iph)*dt)
            next(iph)=prev(iph)*atten+src(iph)*(1d0-atten)/loss(iph)
        else
            next(iph)=prev(iph)+src(iph)*dt
        end if
    end do
end subroutine advance_photon

end module hadronic_cascade
