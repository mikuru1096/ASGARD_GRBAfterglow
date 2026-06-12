!f2py: skip
module hadronic_transport_kernel
    use constants
    use hadronic_common
    implicit none

contains

! 质子幂律注入：Q_inj ∝ γ^(-p)，按给定能量预算归一化。
subroutine hadronic_proton_injection_powerlaw(Num_gam_p,gam_p,p_p,energy_budget_erg,gam_p_min,gam_p_max,Q_inj)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p
    real(8), intent(in) :: gam_p(Num_gam_p),p_p,energy_budget_erg,gam_p_min,gam_p_max
    real(8), intent(out) :: Q_inj(Num_gam_p)
    integer :: I_gam,i_lo,i_hi
    real(8) :: gam_edge(Num_gam_p+1),norm,moment,weight,gam_mid,dgam

    call hadronic_build_gamma_edges(Num_gam_p,gam_p,gam_edge)
    call hadronic_source_bounds(Num_gam_p,gam_p,gam_p_min,gam_p_max,i_lo,i_hi)
    Q_inj=zero
    if (energy_budget_erg <= zero .or. i_hi < i_lo) return

    call accumulate_powerlaw_energy_moment()
    if (norm <= zero) return
    call write_powerlaw_injection()

contains

    subroutine accumulate_powerlaw_energy_moment()
    implicit none

    norm=zero
    do I_gam=i_lo,i_hi
        gam_mid=gam_p(I_gam)
        dgam=gam_edge(I_gam+1)-gam_edge(I_gam)
        weight=gam_mid**(-p_p)
        moment=(gam_mid-one)*Para_m_p_E
        norm=norm+weight*moment*dgam
    end do
    end subroutine accumulate_powerlaw_energy_moment

    subroutine write_powerlaw_injection()
    implicit none

    do I_gam=i_lo,i_hi
        gam_mid=gam_p(I_gam)
        Q_inj(I_gam)=energy_budget_erg*gam_mid**(-p_p)/norm
    end do
    end subroutine write_powerlaw_injection
end subroutine hadronic_proton_injection_powerlaw

! 计算质子能量损失率：绝热冷却 dγ/dt = γ/t_dyn 和同步冷却 dγ/dt ∝ γ²。
subroutine hadronic_proton_loss_rates(Num_gam_p,gam_p,B_field_g,t_dyn_s,loss_ad,loss_syn,loss_total)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p
    real(8), intent(in) :: gam_p(Num_gam_p),B_field_g,t_dyn_s
    real(8), intent(out) :: loss_ad(Num_gam_p),loss_syn(Num_gam_p),loss_total(Num_gam_p)
    integer :: I_gam
    real(8) :: coeff_syn

    if (B_field_g < zero) error stop "hadronic proton loss rates require B_field_g >= 0."
    if (t_dyn_s <= zero) error stop "hadronic proton loss rates require t_dyn_s > 0."
    coeff_syn=Para_sigmaT*B_field_g*B_field_g/(6d0*pi*Para_m_e*Para_c) / (Para_m_p_div_m_e**3)
    do I_gam=1,Num_gam_p
        loss_ad(I_gam)=gam_p(I_gam)/t_dyn_s
        loss_syn(I_gam)=coeff_syn*gam_p(I_gam)*gam_p(I_gam)
        loss_total(I_gam)=loss_ad(I_gam)+loss_syn(I_gam)
    end do
end subroutine hadronic_proton_loss_rates

! 对数gamma空间迎风输运推进：用边心公式计算损失通量，更新粒子谱。
subroutine hadronic_advance_energy_loggamma(Num_gam_p,gam_p,dN_prev,Q_inj,loss_total,dt_s,dN_next)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p
    real(8), intent(in) :: gam_p(Num_gam_p),dN_prev(Num_gam_p),Q_inj(Num_gam_p),loss_total(Num_gam_p),dt_s
    real(8), intent(out) :: dN_next(Num_gam_p)
    integer :: I_gam
    real(8) :: gam_edge(Num_gam_p+1),flux_edge(Num_gam_p+1),dgam_cell
    real(8) :: loss_edge,source_loc

    call hadronic_build_gamma_edges(Num_gam_p,gam_p,gam_edge)
    if (dt_s <= zero) error stop "hadronic energy advance requires dt_s > 0."
    call build_loss_flux_edges()
    call apply_flux_divergence_with_injection()

contains

    subroutine build_loss_flux_edges()
    implicit none

    flux_edge=zero
    do I_gam=2,Num_gam_p
        loss_edge=0.5d0*(loss_total(I_gam-1)+loss_total(I_gam))
        flux_edge(I_gam)=loss_edge*dN_prev(I_gam)
    end do
    flux_edge(1)=loss_total(1)*dN_prev(1)
    flux_edge(Num_gam_p+1)=zero
    end subroutine build_loss_flux_edges

    subroutine apply_flux_divergence_with_injection()
    implicit none

    do I_gam=1,Num_gam_p
        dgam_cell=gam_edge(I_gam+1)-gam_edge(I_gam)
        if (dgam_cell <= zero) error stop "hadronic energy advance requires positive gamma cell width."
        ! Q_inj is normalized to the shell energy increment, so it is already
        ! the injected particle content of this step rather than a time rate.
        source_loc=Q_inj(I_gam)
        dN_next(I_gam)=dN_prev(I_gam)+dt_s*(flux_edge(I_gam+1)-flux_edge(I_gam))/dgam_cell+source_loc
        if (dN_next(I_gam) < zero) error stop "hadronic energy advance produced negative particle density."
    end do
    end subroutine apply_flux_divergence_with_injection
end subroutine hadronic_advance_energy_loggamma

end module hadronic_transport_kernel
