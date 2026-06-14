!f2py: skip
module hadronic_common
    use constants
    use quantum_synchrotron_kernel, only: quantum_chi_parameter, quantum_syn_cooling_factor
    implicit none
    ! 粒子静质量 [GeV/c^2]，引用 constants 模块定义。
    real(8), parameter :: hadronic_electron_mass_gev = Para_m_e_GeV
    real(8), parameter :: hadronic_proton_mass_gev = Para_m_p_GeV
    real(8), parameter :: hadronic_neutron_mass_gev = Para_m_n_GeV
    real(8), parameter :: hadronic_pion_charged_mass_gev = Para_m_pi_charged_GeV
    real(8), parameter :: hadronic_pion0_mass_gev = Para_m_pi0_GeV
    real(8), parameter :: hadronic_muon_mass_gev = Para_m_mu_GeV

contains

! 构建质子的对数均匀洛伦兹因子网格。
subroutine hadronic_build_gamma_p_grid(Num_gam_p,gam_p_min,gam_p_max,gam_p)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p
    real(8), intent(in) :: gam_p_min,gam_p_max
    real(8), intent(out) :: gam_p(Num_gam_p)
    integer :: I_gam
    real(8) :: x_min,x_max,dx

    if (Num_gam_p < 1) error stop "hadronic gamma grid must contain at least one point."
    if (gam_p_min <= one) error stop "hadronic gamma grid minimum must exceed 1."
    if (gam_p_max <= gam_p_min) error stop "hadronic gamma grid maximum must exceed minimum."
    if (Num_gam_p == 1) then
        gam_p(1)=gam_p_min
        return
    end if

    x_min=dlog10(gam_p_min)
    x_max=dlog10(gam_p_max)
    dx=(x_max-x_min)/dble(Num_gam_p-1)
    do I_gam=1,Num_gam_p
        gam_p(I_gam)=ten**(x_min+dx*dble(I_gam-1))
    end do
end subroutine hadronic_build_gamma_p_grid

! 在洛伦兹因子网格中查找源项范围 [gam_p_min, gam_p_max] 对应的索引上下界。
subroutine hadronic_source_bounds(Num_gam_p,gam_p,gam_p_min,gam_p_max,i_lo,i_hi)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p
    real(8), intent(in) :: gam_p(Num_gam_p),gam_p_min,gam_p_max
    integer, intent(out) :: i_lo,i_hi
    integer :: I_gam

    i_lo=1
    i_hi=Num_gam_p
    do I_gam=1,Num_gam_p
        if (gam_p(I_gam) >= gam_p_min) then
            i_lo=I_gam
            exit
        end if
    end do
    do I_gam=Num_gam_p,1,-1
        if (gam_p(I_gam) <= gam_p_max) then
            i_hi=I_gam
            exit
        end if
    end do
end subroutine hadronic_source_bounds

! 由网格中心值构建网格单元边界，边界取相邻中心值的几何平均。
subroutine hadronic_build_gamma_edges(Num_gam_p,gam_p,gam_edge)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p
    real(8), intent(in) :: gam_p(Num_gam_p)
    real(8), intent(out) :: gam_edge(Num_gam_p+1)
    integer :: I_gam

    if (Num_gam_p == 1) then
        if (gam_p(1) <= one) error stop "hadronic gamma center must exceed 1."
        gam_edge(1)=0.5d0*gam_p(1)
        gam_edge(2)=2d0*gam_p(1)
        return
    end if

    if (gam_p(1) < one) error stop "hadronic gamma grid minimum must be at least 1."
    gam_edge(1)=gam_p(1)*dsqrt(gam_p(1)/gam_p(2))
    do I_gam=2,Num_gam_p
        gam_edge(I_gam)=dsqrt(gam_p(I_gam-1)*gam_p(I_gam))
    end do
    gam_edge(Num_gam_p+1)=gam_p(Num_gam_p)*dsqrt(gam_p(Num_gam_p)/gam_p(Num_gam_p-1))
end subroutine hadronic_build_gamma_edges

! 计算壳层在观测者时间中的时间步长 dt = R_tobs(i) - R_tobs(i-1)。
real(8) function hadronic_shell_dt(R_tobs,i_shell)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: i_shell
    real(8), intent(in) :: R_tobs(*)

    if (i_shell <= 1) then
        hadronic_shell_dt=R_tobs(1)
    else
        hadronic_shell_dt=R_tobs(i_shell)-R_tobs(i_shell-1)
    end if
    if (hadronic_shell_dt <= zero) error stop "hadronic shell dt must be positive."
end function hadronic_shell_dt

! 计算动力学时标 t_dyn = R / (Gamma_bulk * c)。
real(8) function hadronic_dynamical_time(R_loc,Gamma_bulk)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: R_loc,Gamma_bulk

    if (R_loc <= zero) error stop "hadronic dynamical time requires positive radius."
    if (Gamma_bulk < one) error stop "hadronic dynamical time requires Gamma_bulk >= 1."
    hadronic_dynamical_time=R_loc/(Gamma_bulk*Para_c)
end function hadronic_dynamical_time

! 由加速-冷却平衡估计质子最大洛伦兹因子，取动力学限制和同步冷却限制的较小值。
real(8) function hadronic_gamma_p_max(B_field_g,t_dyn_s,eta_acc)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: B_field_g,t_dyn_s,eta_acc
    real(8) :: gam_dyn,gam_syn

    if (B_field_g <= zero) error stop "hadronic_gamma_p_max requires B_field_g > 0."
    if (t_dyn_s <= zero) error stop "hadronic_gamma_p_max requires t_dyn_s > 0."
    if (eta_acc <= zero) error stop "hadronic_gamma_p_max requires eta_acc > 0."
    gam_dyn=Para_e*B_field_g*t_dyn_s/(eta_acc*Para_m_p*Para_c)
    gam_syn=dsqrt(6d0*pi*Para_e/(eta_acc*Para_sigmaT*B_field_g)) * Para_m_p_div_m_e
    hadronic_gamma_p_max=min(gam_dyn,gam_syn)
end function hadronic_gamma_p_max

! 验证能量网格为严格递增、正值且对数均匀分布，可选返回对数间距。
subroutine hadronic_validate_log_grid(num_grid,energy_grid,name,dln_out)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: num_grid
    real(8), intent(in) :: energy_grid(num_grid)
    character(*), intent(in) :: name
    real(8), intent(out), optional :: dln_out
    integer :: i
    real(8) :: dln_local,dln_i

    if (num_grid < 2) error stop trim(name)//" must contain at least two points."
    if (energy_grid(1) <= zero) error stop trim(name)//" must be strictly positive."
    do i=2,num_grid
        if (energy_grid(i) <= zero) error stop trim(name)//" must be strictly positive."
        if (energy_grid(i) <= energy_grid(i-1)) error stop trim(name)//" must be strictly increasing."
    end do
    dln_local = dlog(energy_grid(2)/energy_grid(1))
    do i=3,num_grid
        dln_i = dlog(energy_grid(i)/energy_grid(i-1))
        if (dabs(dln_i-dln_local) > dmax1(1d-12,1d-6*dabs(dln_local))) &
            error stop trim(name)//" must be logarithmically uniform."
    end do
    if (present(dln_out)) dln_out = dln_local
end subroutine hadronic_validate_log_grid

! Landau 量子同步压制因子: f(χ) = 1/(1+√2·χ^(2/3))²
! 适用于质子同步冷却率修正。
! χ = γ * B / B_crit * m_e/m_particle
real(8) function hadronic_quantum_syn_cooling_factor(gamma,b_field_g,mass_gev)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: gamma,b_field_g,mass_gev
    real(8) :: chi

    if (b_field_g <= 0d0 .or. gamma <= 1d0) then
        hadronic_quantum_syn_cooling_factor = 1d0
        return
    end if
    chi = quantum_chi_parameter(gamma,b_field_g,mass_gev)
    hadronic_quantum_syn_cooling_factor = quantum_syn_cooling_factor(chi)
end function hadronic_quantum_syn_cooling_factor

end module hadronic_common
