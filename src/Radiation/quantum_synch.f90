! 量子同步辐射修正: Landau 压制因子(AM3 Klinger+23 eq.A9)。
! Quantum synchrotron correction: Landau suppression factor.
module quantum_synch
    use constants
    implicit none
    private

    public :: quantum_cooling
    public :: quantum_chi_parameter

    ! 电子临界磁场 B_c = m_e^2 c^3 / (e hbar) [G]。
    ! Electron critical magnetic field B_c = m_e^2 c^3 / (e hbar) [G].
    real(8), parameter :: bcrit = 4.414d13

contains

! 量子参数 chi = gamma * B / B_c * (m_e / m_particle)。
! Quantum parameter chi = gamma * B / B_c * (m_e / m_particle).
pure real(8) function quantum_chi_parameter(gamma,b_field_g,mass_gev)
    real(8), intent(in) :: gamma,b_field_g,mass_gev
    real(8) :: mass_ratio

    mass_ratio = Para_m_e_GeV / mass_gev
    quantum_chi_parameter = gamma * b_field_g / bcrit * mass_ratio
end function quantum_chi_parameter

! Landau 压制因子: f(chi) = 1 / (1 + sqrt(2) * chi^(2/3))^2。
! Landau suppression factor: ratio of QED-corrected cooling to classical cooling.
pure real(8) function quantum_cooling(chi)
    real(8), intent(in) :: chi
    real(8) :: chi23

    if (chi <= 1d-6) then
        quantum_cooling = 1d0
        return
    end if

    chi23 = chi**(2d0/3d0)
    quantum_cooling = 1d0 / (1d0 + dsqrt(2d0) * chi23)**2
end function quantum_cooling

end module quantum_synch
