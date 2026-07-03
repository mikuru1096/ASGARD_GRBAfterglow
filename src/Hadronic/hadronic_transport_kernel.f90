!f2py: skip
module hadronic_transport
    use constants
    use hadronic_base
    implicit none

contains

! 质子幂律注入：qinj ∝ γ^(-p)，按给定能量预算归一化。
subroutine proton_inject(ng,gp,pidx,ebudget,gmin,gmax,qinj)
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: gp
    real(8), intent(in) :: pidx,ebudget,gmin,gmax
    real(8), intent(out), dimension(ng) :: qinj
    integer :: ig,ilo,ihi
    real(8), dimension(ng+1) :: edge
    real(8) :: norm,moment,weight,gmid,dgam

    call build_edges(ng,gp,edge)
    call source_bounds(ng,gp,gmin,gmax,ilo,ihi)
    qinj=0d0
    if (ebudget <= 0d0 .or. ihi < ilo) return

    ! 先积分能量矩作为归一化分母。 / First integrate the energy moment for normalization.
    norm=0d0
    do ig=ilo,ihi
        gmid=gp(ig)
        dgam=edge(ig+1)-edge(ig)
        weight=gmid**(-pidx)
        moment=(gmid-1d0)*Para_m_p_E
        norm=norm+weight*moment*dgam
    end do
    if (norm <= 0d0) return
    ! 再写入每个 gamma bin 的注入率。 / Then write the injection rate per gamma bin.
    do ig=ilo,ihi
        gmid=gp(ig)
        qinj(ig)=ebudget*gmid**(-pidx)/norm
    end do
end subroutine proton_inject

! 计算质子能量损失率：绝热冷却 dγ/dt = γ/t_dyn 和同步冷却 dγ/dt ∝ γ²。
subroutine proton_loss(ng,gp,bfield,tdyn,lad,lsyn,ltot)
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: gp
    real(8), intent(in) :: bfield,tdyn
    real(8), intent(out), dimension(ng) :: lad,lsyn,ltot
    integer :: ig
    real(8) :: coeff_syn

    if (bfield < 0d0) error stop "hadronic proton loss rates require bfield >= 0."
    if (tdyn <= 0d0) error stop "hadronic proton loss rates require tdyn > 0."
    coeff_syn=Para_sigmaT*bfield*bfield/(6d0*pi*Para_m_e*Para_c) / (Para_m_p_DIV_m_e**3)
    do ig=1,ng
        lad(ig)=gp(ig)/tdyn
        lsyn(ig)=coeff_syn*gp(ig)*gp(ig)
        ltot(ig)=lad(ig)+lsyn(ig)
    end do
end subroutine proton_loss

! 对数gamma空间迎风输运推进：用边心公式计算损失通量，更新粒子谱。
subroutine advance_loggamma(ng,gp,prev,qinj,ltot,dt,next)
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: gp,prev,qinj,ltot
    real(8), intent(in) :: dt
    real(8), intent(out), dimension(ng) :: next
    integer :: ig
    real(8), dimension(ng+1) :: edge,flux
    real(8) :: dgam
    real(8) :: ledge,src

    call build_edges(ng,gp,edge)
    if (dt <= 0d0) error stop "hadronic energy advance requires dt > 0."
    ! 边界通量使用相邻 cell 的损失率平均。 / Boundary flux uses averaged neighboring loss rates.
    flux=0d0
    do ig=2,ng
        ledge=0.5d0*(ltot(ig-1)+ltot(ig))
        flux(ig)=ledge*prev(ig)
    end do
    flux(1)=ltot(1)*prev(1)
    flux(ng+1)=0d0
    ! 通量散度和本步注入一起更新。 / Flux divergence and current injection advance together.
    do ig=1,ng
        dgam=edge(ig+1)-edge(ig)
        if (dgam <= 0d0) error stop "hadronic energy advance requires positive gamma cell width."
        src=qinj(ig)
        next(ig)=prev(ig)+dt*(flux(ig+1)-flux(ig))/dgam+src
        if (next(ig) < 0d0) error stop "hadronic energy advance produced negative particle density."
    end do
end subroutine advance_loggamma

end module hadronic_transport
