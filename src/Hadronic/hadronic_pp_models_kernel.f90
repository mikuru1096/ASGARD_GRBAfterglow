! pp 多模型谱形核 (SIBYLL/QGSJET/Geant4/Pythia8)
! 从 AM3 ProtonProton.cc (Klinger+23) 移植参数化公式
module hadronic_pp_models_kernel
    use constants
    use hadronic_common, only: hadronic_proton_mass_gev, hadronic_pion0_mass_gev
    implicit none
    private

    public :: hadronic_pp_spectral_shape
    public :: hadronic_pp_pi0_source_spectrum
    public :: hadronic_pp_threshold_kinetic_energy_gev

    integer, parameter :: MODEL_SIBYLL=0, MODEL_QGSJET=1, MODEL_GEANT4=2, MODEL_PYTHIA8=3
    real(8), parameter :: mbarn_to_cm2 = 1d-27
    real(8), parameter :: mpi0 = 0.1349768d0
    real(8), parameter :: mp = 0.9382720813d0

contains

! pp 产生 π⁰ 的动能阈值: Tp_th = 2*m_pi0 + m_pi0²/(2*m_p)
real(8) function hadronic_pp_threshold_kinetic_energy_gev()
    hadronic_pp_threshold_kinetic_energy_gev = 2d0*mpi0 + mpi0*mpi0/(2d0*mp)
end function

! ------------------------------------------------------------
! 谱形函数 F(Tp, Eγ): 归一化无量纲形状
! 输入: Tp 质子动能 [GeV], Egam 光子能量 [GeV], model 模型索引
real(8) function hadronic_pp_spectral_shape(Tp,Egam,model)
    real(8), intent(in) :: Tp,Egam
    integer, intent(in) :: model
    real(8) :: Y,Y0,X,C
    real(8) :: Egam_max_val

    Egam_max_val = Egam_max(Tp)
    if (Egam <= 0d0 .or. Egam >= Egam_max_val .or. Tp <= 0d0) then
        hadronic_pp_spectral_shape = 0d0; return
    end if

    Y = Egam + mpi0*mpi0/(4d0*Egam)
    Y0 = Egam_max_val + mpi0*mpi0/(4d0*Egam_max_val)
    X = (Y - mpi0) / (Y0 - mpi0)
    C = 3.55d0 * mpi0 / Y0  ! SIBYLL/QGSJET default; Geant4 overrides

    select case (model)
    case (MODEL_SIBYLL)
        hadronic_pp_spectral_shape = high_energy_or_geant4_shape(3.6d0)

    case (MODEL_QGSJET)
        hadronic_pp_spectral_shape = high_energy_or_geant4_shape(4.5d0)

    case (MODEL_GEANT4)
        hadronic_pp_spectral_shape = F_geant4_local(Tp,Egam)

    case (MODEL_PYTHIA8)
        hadronic_pp_spectral_shape = high_energy_or_geant4_shape(4.2d0)

    case default
        hadronic_pp_spectral_shape = F_geant4_local(Tp,Egam)
    end select

contains

    real(8) function high_energy_or_geant4_shape(exponent)
        real(8), intent(in) :: exponent

        if (Tp > 50d0 .and. X >= 0d0 .and. X < 1d0) then
            high_energy_or_geant4_shape = (1d0 - dsqrt(X))**exponent / (1d0 + X/C)
        else
            high_energy_or_geant4_shape = F_geant4_local(Tp,Egam)
        end if
    end function high_energy_or_geant4_shape
end function

! ------------------------------------------------------------
! Geant4 谱形 (所有模型的低能/中能 fallback)
real(8) function F_geant4_local(Tp,Egam)
    real(8), intent(in) :: Tp,Egam
    real(8) :: Y,Y0,X,C,theta,kappa,q,mu,beta_f,gamma_f

    Y = Egam + mpi0*mpi0/(4d0*Egam)
    Y0 = Egam_max(Tp) + mpi0*mpi0/(4d0*Egam_max(Tp))
    X = (Y - mpi0) / (Y0 - mpi0)
    C = 3d0 * mpi0 / Y0
    theta = Tp / mp
    q = (Tp - 1d0) / mp

    if (X < 0d0 .or. X >= 1d0) then
        F_geant4_local = 0d0; return
    end if

    if (Tp < 1d0) then
        kappa = 3.29d0 - 0.2d0 * theta**(-1.5d0)
        F_geant4_local = (1d0 - X)**kappa
    else if (Tp <= 4d0) then
        mu = 1.25d0 * q**1.25d0 * dexp(-1.25d0*q)
        beta_f = mu + 2.45d0; gamma_f = mu + 1.45d0
        F_geant4_local = (1d0 - X)**beta_f / ((1d0 + X/C)**gamma_f)
    else if (Tp <= 20d0) then
        mu = 1.25d0 * q**1.25d0 * dexp(-1.25d0*q)
        beta_f = 1.5d0*mu + 4.95d0; gamma_f = mu + 1.5d0
        F_geant4_local = (1d0 - X)**beta_f / ((1d0 + X/C)**gamma_f)
    else if (Tp <= 100d0) then
        F_geant4_local = (1d0 - dsqrt(X))**4.2d0 / (1d0 + X/C)
    else
        F_geant4_local = (1d0 - dsqrt(X))**4.9d0 / (1d0 + X/C)
    end if
end function

! ------------------------------------------------------------
! π⁰ 衰变光子最大能量
real(8) function Egam_max(Tp)
    real(8), intent(in) :: Tp
    Egam_max = 0.5d0 * (Tp + dsqrt(Tp*(Tp + 2d0*mp)))
end function

! ------------------------------------------------------------
! π⁰ 源谱: dN_γ/(dE_γ dt dV) [cm⁻³ s⁻¹ GeV⁻¹], 调用谱形函数
! n_H: 靶质子数密度 [cm⁻³], n_p(Tp): 入射质子谱 [cm⁻³ GeV⁻¹]
subroutine hadronic_pp_pi0_source_spectrum(num_p,Tp_grid,n_p_per_gev, &
    num_g,Egam_grid,n_H,model,pi0_gamma_rate)
    integer, intent(in) :: num_p,num_g,model
    real(8), intent(in) :: Tp_grid(num_p),n_p_per_gev(num_p)
    real(8), intent(in) :: Egam_grid(num_g),n_H
    real(8), intent(out) :: pi0_gamma_rate(num_g)
    integer :: i_p,i_g
    real(8) :: Tp,dln_Tp,Amax_val,sigma_pi0_val,dsde,c0

    pi0_gamma_rate(1:num_g) = 0d0
    if (num_p < 2) return
    dln_Tp = dlog(Tp_grid(2)/Tp_grid(1))
    c0 = Para_c * mbarn_to_cm2 * n_H

    do i_p=1,num_p
        Tp = Tp_grid(i_p)
        if (Tp <= hadronic_pp_threshold_kinetic_energy_gev() .or. n_p_per_gev(i_p) <= 0d0) cycle
        sigma_pi0_val = sigma_pi0_model(Tp,model)
        Amax_val = Amax_model(Tp,model)

        do i_g=1,num_g
            if (Egam_grid(i_g) <= 0d0 .or. Egam_grid(i_g) >= Egam_max(Tp)) cycle
            dsde = mbarn_to_cm2 * Amax_val * hadronic_pp_spectral_shape(Tp,Egam_grid(i_g),model)
            pi0_gamma_rate(i_g) = pi0_gamma_rate(i_g) + &
                c0 * n_p_per_gev(i_p) * dsde * dln_Tp
        end do
    end do
end subroutine

! ------------------------------------------------------------
! 总 π⁰ 产生截面 [mbarn]
real(8) function sigma_pi0_model(Tp,model)
    real(8), intent(in) :: Tp
    integer, intent(in) :: model
    real(8) :: sigma_1pi,sigma_2pi,s_inel,mult

    sigma_1pi = sigma_1_pi_local(Tp)
    sigma_2pi = sigma_2_pi_local(Tp)
    s_inel = sigma_inel_local(Tp)

    select case (model)
    case (MODEL_SIBYLL);  mult = multip_pi0_SIBYLL(Tp)
    case (MODEL_QGSJET);  mult = multip_pi0_QGSJET(Tp)
    case (MODEL_GEANT4);  mult = multip_pi0_Geant4(Tp)
    case (MODEL_PYTHIA8); mult = multip_pi0_Pythia8(Tp)
    case default;         mult = multip_pi0_SIBYLL(Tp)
    end select
    sigma_pi0_model = sigma_1pi + sigma_2pi + s_inel * mult
end function

! ------------------------------------------------------------
! Amax 归一化因子 [GeV⁻¹]
real(8) function Amax_model(Tp,model)
    real(8), intent(in) :: Tp
    integer, intent(in) :: model
    real(8) :: theta_p,Ltheta,b1,b2,b3

    theta_p = Tp / mp
    if (theta_p <= 1d-6) then; Amax_model = 0d0; return; end if
    Ltheta = dlog(theta_p)

    select case (model)
    case (MODEL_SIBYLL)
        b1=5.9d0; b2=0.9d0; b3=0.16d0
    case (MODEL_QGSJET)
        b1=5.7d0; b2=0.8d0; b3=0.18d0
    case (MODEL_PYTHIA8)
        b1=6.2d0; b2=0.85d0; b3=0.14d0
    case default
        b1=5.9d0; b2=0.9d0; b3=0.15d0
    end select
    Amax_model = b1 * theta_p**b2 * dexp(b3 * Ltheta*Ltheta) * &
                 sigma_pi0_model(Tp,model) / mp
end function

! ---- 低能截面 (Kafexhiu+2014) ----
real(8) function sigma_1_pi_local(Tp)
    real(8), intent(in) :: Tp
    real(8) :: eta
    if (Tp <= hadronic_pp_threshold_kinetic_energy_gev()) then
        sigma_1_pi_local = 0d0; return
    end if
    eta = (Tp - hadronic_pp_threshold_kinetic_energy_gev()) / 0.4d0
    sigma_1_pi_local = 7.66d-3 * eta**2.95d0 * dexp(-1.15d0*eta)  ! mbarn
end function

real(8) function sigma_2_pi_local(Tp)
    real(8), intent(in) :: Tp
    real(8) :: thr2,eta
    thr2 = 2d0*hadronic_pp_threshold_kinetic_energy_gev()
    if (Tp <= thr2) then; sigma_2_pi_local = 0d0; return; end if
    eta = (Tp - thr2) / 0.4d0
    sigma_2_pi_local = 5.7d-3 * eta**3.2d0 * dexp(-1.35d0*eta)
end function

real(8) function sigma_inel_local(Tp)
    real(8), intent(in) :: Tp
    if (Tp <= 0d0) then; sigma_inel_local = 0d0; return; end if
    sigma_inel_local = 30.7d0 - 0.96d0*dlog(Tp) + 0.36d0*(dlog(Tp))**2  ! Kafexhiu+2014 fit
end function

! ---- π⁰ 多重数 ----
real(8) function multip_pi0_SIBYLL(Tp)
    real(8), intent(in) :: Tp
    real(8) :: xi,a1,a2,a4,a5
    if (Tp < 50d0) then; multip_pi0_SIBYLL = multip_pi0_Geant4(Tp); return; end if
    xi = Tp / mp
    a1=0.48d0; a2=0.13d0; a4=0.30d0; a5=1.45d0
    multip_pi0_SIBYLL = a1 * xi**a4 * (1d0 + dexp(-a2*xi**a5))
end function

real(8) function multip_pi0_QGSJET(Tp)
    real(8), intent(in) :: Tp
    real(8) :: xi,a1,a2,a4,a5
    if (Tp < 50d0) then; multip_pi0_QGSJET = multip_pi0_Geant4(Tp); return; end if
    xi = Tp / mp
    a1=0.52d0; a2=0.11d0; a4=0.28d0; a5=1.52d0
    multip_pi0_QGSJET = a1 * xi**a4 * (1d0 + dexp(-a2*xi**a5))
end function

real(8) function multip_pi0_Geant4(Tp)
    real(8), intent(in) :: Tp
    real(8) :: xi
    xi = Tp / mp
    multip_pi0_Geant4 = 0.45d0 * xi**0.25d0 * (1d0 + 0.1d0*xi**(-0.5d0))
end function

real(8) function multip_pi0_Pythia8(Tp)
    real(8), intent(in) :: Tp
    real(8) :: xi,a1,a2,a4,a5
    if (Tp < 50d0) then; multip_pi0_Pythia8 = multip_pi0_Geant4(Tp); return; end if
    xi = Tp / mp
    a1=0.55d0; a2=0.12d0; a4=0.27d0; a5=1.48d0
    multip_pi0_Pythia8 = a1 * xi**a4 * (1d0 + dexp(-a2*xi**a5))
end function

end module hadronic_pp_models_kernel
