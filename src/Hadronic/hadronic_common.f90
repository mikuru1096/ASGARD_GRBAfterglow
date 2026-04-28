!f2py: skip
module hadronic_common
    use constants
    implicit none
    real(8), parameter :: hadronic_eta_acc_floor = 1d-12
    real(8), parameter :: hadronic_bfield_floor = 1d-30
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

    if (Num_gam_p <= 1) then
        gam_p(1)=max(gam_p_min,one)
        return
    end if

    x_min=dlog10(max(gam_p_min,one))
    x_max=dlog10(max(gam_p_max,ten))
    dx=(x_max-x_min)/dble(Num_gam_p-1)
    do I_gam=1,Num_gam_p
        gam_p(I_gam)=ten**(x_min+dx*dble(I_gam-1))
    end do
end subroutine hadronic_build_gamma_p_grid

! 初始化粒子数密度为零。
subroutine hadronic_initial_density(Num_gam_p,dN_gam_p)
    integer, intent(in) :: Num_gam_p
    real(8), intent(out) :: dN_gam_p(Num_gam_p)

    dN_gam_p=zero
end subroutine hadronic_initial_density

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
        gam_edge(1)=max(one,0.5d0*gam_p(1))
        gam_edge(2)=2d0*gam_p(1)
        return
    end if

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
        hadronic_shell_dt=max(R_tobs(1),one)
    else
        hadronic_shell_dt=max(R_tobs(i_shell)-R_tobs(i_shell-1),one)
    end if
end function hadronic_shell_dt

! 计算动力学时标 t_dyn = R / (Gamma_bulk * c)。
real(8) function hadronic_dynamical_time(R_loc,Gamma_bulk)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: R_loc,Gamma_bulk

    hadronic_dynamical_time=max(R_loc/(max(Gamma_bulk,one)*Para_c),one)
end function hadronic_dynamical_time

! 由加速-冷却平衡估计质子最大洛伦兹因子，取动力学限制和同步冷却限制的较小值。
real(8) function hadronic_gamma_p_max(B_field_g,t_dyn_s,eta_acc)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: B_field_g,t_dyn_s,eta_acc
    real(8) :: gam_dyn,gam_syn

    gam_dyn=Para_e*B_field_g*t_dyn_s/(max(eta_acc,hadronic_eta_acc_floor)*Para_m_p*Para_c)
    gam_syn=dsqrt(6d0*pi*Para_e/(max(eta_acc,hadronic_eta_acc_floor)*Para_sigmaT* &
            max(B_field_g,hadronic_bfield_floor))) * Para_m_p_div_m_e
    hadronic_gamma_p_max=max(ten,min(gam_dyn,gam_syn))
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
    real(8), parameter :: b_crit=4.414d13
    real(8) :: chi,chi23

    if (b_field_g <= 0d0 .or. gamma <= 1d0) then
        hadronic_quantum_syn_cooling_factor = 1d0
        return
    end if
    chi = gamma * b_field_g / b_crit * (hadronic_electron_mass_gev / max(mass_gev,1d-30))
    if (chi <= 1d-6) then
        hadronic_quantum_syn_cooling_factor = 1d0
        return
    end if
    chi23 = chi**(2d0/3d0)
    hadronic_quantum_syn_cooling_factor = 1d0 / (1d0 + dsqrt(2d0)*chi23)**2
end function hadronic_quantum_syn_cooling_factor

end module hadronic_common
