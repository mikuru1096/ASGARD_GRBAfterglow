!f2py: skip
! pp 多模型谱形核：SIBYLL/QGSJET/Geant4/Pythia8。
! Multi-model pp spectral kernel from the AM3/Klinger+23 parametrization.
module pp_models
    use constants
    implicit none
    private

    public :: pp_shape
    public :: pi0_source
    public :: pp_threshold

    integer, parameter :: MODEL_SIBYLL=0, MODEL_QGSJET=1, MODEL_GEANT4=2, MODEL_PYTHIA8=3
    real(8), parameter :: mbarn = 1d-27
    real(8), parameter :: mpi0 = 0.1349768d0
    real(8), parameter :: mp = 0.9382720813d0

contains

! pp -> pi0 的动能阈值。
! Kinetic threshold for pp -> pi0 production.
real(8) function pp_threshold()
    pp_threshold = 2d0*mpi0 + mpi0*mpi0/(2d0*mp)
end function

! 谱形函数 F(tp, egam)：归一化无量纲形状。
! Dimensionless spectral shape F(tp, egam).
real(8) function pp_shape(tp,egam,model)
    real(8), intent(in) :: tp,egam
    integer, intent(in) :: model
    real(8) :: y,y0,x,c
    real(8) :: emax

    emax = egam_max(tp)
    if (egam <= 0d0 .or. egam >= emax .or. tp <= 0d0) then
        pp_shape = 0d0; return
    end if

    y = egam + mpi0*mpi0/(4d0*egam)
    y0 = emax + mpi0*mpi0/(4d0*emax)
    x = (y - mpi0) / (y0 - mpi0)
    c = 3.55d0 * mpi0 / y0  ! SIBYLL/QGSJET 默认；Geant4 分支会覆盖。 / Default except Geant4.

    select case (model)
    case (MODEL_SIBYLL)
        pp_shape = high_shape(3.6d0)

    case (MODEL_QGSJET)
        pp_shape = high_shape(4.5d0)

    case (MODEL_GEANT4)
        pp_shape = geant4_shape(tp,egam)

    case (MODEL_PYTHIA8)
        pp_shape = high_shape(4.2d0)

    case default
        error stop "pp_shape: unsupported pp model."
    end select

contains

    real(8) function high_shape(exponent)
        real(8), intent(in) :: exponent

        if (tp > 50d0 .and. x >= 0d0 .and. x < 1d0) then
            high_shape = (1d0 - dsqrt(x))**exponent / (1d0 + x/c)
        else
            high_shape = geant4_shape(tp,egam)
        end if
    end function high_shape
end function

! Geant4 谱形；高能模型在低/中能区也回到该分支。
! Geant4 branch, reused by high-energy models in low/intermediate energies.
real(8) function geant4_shape(tp,egam)
    real(8), intent(in) :: tp,egam
    real(8) :: y,y0,x,c,theta,kappa,q,mu,betaf,gamf

    y = egam + mpi0*mpi0/(4d0*egam)
    y0 = egam_max(tp) + mpi0*mpi0/(4d0*egam_max(tp))
    x = (y - mpi0) / (y0 - mpi0)
    c = 3d0 * mpi0 / y0
    theta = tp / mp
    q = (tp - 1d0) / mp

    if (x < 0d0 .or. x >= 1d0) then
        geant4_shape = 0d0; return
    end if

    if (tp < 1d0) then
        kappa = 3.29d0 - 0.2d0 * theta**(-1.5d0)
        geant4_shape = (1d0 - x)**kappa
    else if (tp <= 4d0) then
        mu = 1.25d0 * q**1.25d0 * dexp(-1.25d0*q)
        betaf = mu + 2.45d0; gamf = mu + 1.45d0
        geant4_shape = (1d0 - x)**betaf / ((1d0 + x/c)**gamf)
    else if (tp <= 20d0) then
        mu = 1.25d0 * q**1.25d0 * dexp(-1.25d0*q)
        betaf = 1.5d0*mu + 4.95d0; gamf = mu + 1.5d0
        geant4_shape = (1d0 - x)**betaf / ((1d0 + x/c)**gamf)
    else if (tp <= 100d0) then
        geant4_shape = (1d0 - dsqrt(x))**4.2d0 / (1d0 + x/c)
    else
        geant4_shape = (1d0 - dsqrt(x))**4.9d0 / (1d0 + x/c)
    end if
end function

! pi0 衰变光子的最大能量。
! Maximum gamma-ray energy from pi0 decay.
real(8) function egam_max(tp)
    real(8), intent(in) :: tp
    egam_max = 0.5d0 * (tp + dsqrt(tp*(tp + 2d0*mp)))
end function

! pi0 源谱：dN_gamma/(dE_gamma dt dV)，调用谱形函数。
! pi0 source spectrum using the selected pp spectral shape.
subroutine pi0_source(np,tpgrid,npgev, &
    ng,egrid,nh,model,qpi0)
    integer, intent(in) :: np,ng,model
    real(8), intent(in), dimension(np) :: tpgrid,npgev
    real(8), intent(in), dimension(ng) :: egrid
    real(8), intent(in) :: nh
    real(8), intent(out), dimension(ng) :: qpi0
    integer :: ip,ig
    real(8) :: tp,dlntp,aval,sigval,dsde,c0

    qpi0(1:ng) = 0d0
    if (np < 2) return
    dlntp = dlog(tpgrid(2)/tpgrid(1))
    c0 = Para_c * mbarn * nh

    do ip=1,np
        tp = tpgrid(ip)
        if (tp <= pp_threshold() .or. npgev(ip) <= 0d0) cycle
        sigval = sigpi0(tp,model)
        aval = amax_model(tp,model)

        do ig=1,ng
            if (egrid(ig) <= 0d0 .or. egrid(ig) >= egam_max(tp)) cycle
            dsde = mbarn * aval * pp_shape(tp,egrid(ig),model)
            qpi0(ig) = qpi0(ig) + &
                c0 * npgev(ip) * dsde * dlntp
        end do
    end do
end subroutine

! 总 pi0 产生截面，单位 mbarn。
! Total pi0 production cross section in mbarn.
real(8) function sigpi0(tp,model)
    real(8), intent(in) :: tp
    integer, intent(in) :: model
    real(8) :: sig1,sig2,sinel,mult

    sig1 = sig1pi(tp)
    sig2 = sig2pi(tp)
    sinel = siginel(tp)

    select case (model)
    case (MODEL_SIBYLL);  mult = mult_sibyll(tp)
    case (MODEL_QGSJET);  mult = mult_qgsjet(tp)
    case (MODEL_GEANT4);  mult = mult_geant4(tp)
    case (MODEL_PYTHIA8); mult = mult_pythia(tp)
    case default
        error stop "sigpi0: unsupported pp model."
    end select
    sigpi0 = sig1 + sig2 + sinel * mult
end function

! amax 归一化因子，单位 GeV^-1。
! amax normalization factor in GeV^-1.
real(8) function amax_model(tp,model)
    real(8), intent(in) :: tp
    integer, intent(in) :: model
    real(8) :: thetap,ltheta,b1,b2,b3

    thetap = tp / mp
    if (thetap <= 1d-6) then; amax_model = 0d0; return; end if
    ltheta = dlog(thetap)

    select case (model)
    case (MODEL_SIBYLL)
        b1=5.9d0; b2=0.9d0; b3=0.16d0
    case (MODEL_QGSJET)
        b1=5.7d0; b2=0.8d0; b3=0.18d0
    case (MODEL_PYTHIA8)
        b1=6.2d0; b2=0.85d0; b3=0.14d0
    case default
        error stop "amax_model: unsupported pp model."
    end select
    amax_model = b1 * thetap**b2 * dexp(b3 * ltheta*ltheta) * &
                 sigpi0(tp,model) / mp
end function

! 低能截面分支，来自 Kafexhiu+2014。
! Low-energy cross-section branches from Kafexhiu+2014.
real(8) function sig1pi(tp)
    real(8), intent(in) :: tp
    real(8) :: eta
    if (tp <= pp_threshold()) then
        sig1pi = 0d0; return
    end if
    eta = (tp - pp_threshold()) / 0.4d0
    sig1pi = 7.66d-3 * eta**2.95d0 * dexp(-1.15d0*eta)  ! mbarn
end function

real(8) function sig2pi(tp)
    real(8), intent(in) :: tp
    real(8) :: thr2,eta
    thr2 = 2d0*pp_threshold()
    if (tp <= thr2) then; sig2pi = 0d0; return; end if
    eta = (tp - thr2) / 0.4d0
    sig2pi = 5.7d-3 * eta**3.2d0 * dexp(-1.35d0*eta)
end function

real(8) function siginel(tp)
    real(8), intent(in) :: tp
    if (tp <= 0d0) then; siginel = 0d0; return; end if
    siginel = 30.7d0 - 0.96d0*dlog(tp) + 0.36d0*(dlog(tp))**2  ! Kafexhiu+2014 fit
end function

! pi0 多重数模型。
! pi0 multiplicity models.
real(8) function mult_sibyll(tp)
    real(8), intent(in) :: tp
    real(8) :: xi,a1,a2,a4,a5
    if (tp < 50d0) then; mult_sibyll = mult_geant4(tp); return; end if
    xi = tp / mp
    a1=0.48d0; a2=0.13d0; a4=0.30d0; a5=1.45d0
    mult_sibyll = a1 * xi**a4 * (1d0 + dexp(-a2*xi**a5))
end function

real(8) function mult_qgsjet(tp)
    real(8), intent(in) :: tp
    real(8) :: xi,a1,a2,a4,a5
    if (tp < 50d0) then; mult_qgsjet = mult_geant4(tp); return; end if
    xi = tp / mp
    a1=0.52d0; a2=0.11d0; a4=0.28d0; a5=1.52d0
    mult_qgsjet = a1 * xi**a4 * (1d0 + dexp(-a2*xi**a5))
end function

real(8) function mult_geant4(tp)
    real(8), intent(in) :: tp
    real(8) :: xi
    xi = tp / mp
    mult_geant4 = 0.45d0 * xi**0.25d0 * (1d0 + 0.1d0*xi**(-0.5d0))
end function

real(8) function mult_pythia(tp)
    real(8), intent(in) :: tp
    real(8) :: xi,a1,a2,a4,a5
    if (tp < 50d0) then; mult_pythia = mult_geant4(tp); return; end if
    xi = tp / mp
    a1=0.55d0; a2=0.12d0; a4=0.27d0; a5=1.48d0
    mult_pythia = a1 * xi**a4 * (1d0 + dexp(-a2*xi**a5))
end function

end module pp_models
