! 量子同步辐射修正：Landau 压制因子 (AM3 Klinger+23 eq.A9)
module quantum_synchrotron_kernel
    use constants
    implicit none
    private

    public :: quantum_syn_cooling_factor
    public :: quantum_chi_parameter

    ! 电子临界磁场 B_c = m_e² c³ / (e ħ) [G]
    real(8), parameter :: b_crit_electron = 4.414d13

contains

! 量子参数 χ = γ * B / (B_crit * m_e/m_particle)
! 对于电子: χ_e = γ * B / B_crit_electron
! 对于质子: χ_p = χ_e * m_e/m_p
pure real(8) function quantum_chi_parameter(gamma,b_field_g,mass_gev)
    real(8), intent(in) :: gamma,b_field_g,mass_gev
    real(8) :: mass_ratio

    mass_ratio = Para_m_e_GeV / mass_gev
    quantum_chi_parameter = gamma * b_field_g / b_crit_electron * mass_ratio
end function quantum_chi_parameter

! Landau 量子压制因子: f(χ) = 1 / (1 + √2 * χ^(2/3))²
! 物理意义: 量子电动力学修正下同步辐射冷却率与经典值之比
! 当 χ ≪ 1: f → 1 (经典极限)
! 当 χ ≫ 1: f → 1/(2 * χ^(4/3)) (强压制)
pure real(8) function quantum_syn_cooling_factor(chi)
    real(8), intent(in) :: chi
    real(8) :: chi23

    if (chi <= 1d-6) then
        quantum_syn_cooling_factor = 1d0
        return
    end if

    chi23 = chi**(2d0/3d0)
    quantum_syn_cooling_factor = 1d0 / (1d0 + dsqrt(2d0) * chi23)**2
end function quantum_syn_cooling_factor

end module quantum_synchrotron_kernel
