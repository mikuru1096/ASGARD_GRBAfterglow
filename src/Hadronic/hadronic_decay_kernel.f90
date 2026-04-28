!f2py: skip
module hadronic_decay_kernel
    use constants
    use hadronic_common, only: hadronic_validate_log_grid
    implicit none
    private
    real(8), parameter :: pi_plus_mass_gev = Para_m_pi_charged_GeV
    real(8), parameter :: muon_mass_gev = Para_m_mu_GeV
    real(8), parameter :: charged_pion_decay_s = 2.6033d-8
    real(8), parameter :: muon_decay_s = 2.1969811d-6

    public :: hadronic_pi0_to_gamma_operator
    public :: hadronic_pion_decay_operator
    public :: hadronic_muon_decay_operator
    public :: hadronic_hummer2010_decay_operator

contains

! π0衰变为两个光子的算子：将π0注入谱转换为光子产生率。
subroutine hadronic_pi0_to_gamma_operator(Num_pion,pion_energy_gev,pion0_source_rate_per_gev, &
                                          Num_gamma,gamma_energy_gev,gamma_rate_per_gev)
    integer, intent(in) :: Num_pion,Num_gamma
    real(8), intent(in) :: pion_energy_gev(Num_pion),pion0_source_rate_per_gev(Num_pion)
    real(8), intent(in) :: gamma_energy_gev(Num_gamma)
    real(8), intent(out) :: gamma_rate_per_gev(Num_gamma)
    integer :: i_gamma,i_pion
    real(8) :: dln_pi,eout,epi,r,sumk

    gamma_rate_per_gev = zero
    dln_pi = hadronic_log_spacing(Num_pion,pion_energy_gev)
    if (dln_pi <= zero) return

    do i_gamma=1,Num_gamma
        eout = gamma_energy_gev(i_gamma)
        sumk = zero
        do i_pion=1,Num_pion
            epi = pion_energy_gev(i_pion)
            r = two*eout/epi
            if (r > zero .and. r <= 1.00001d0) then
                sumk = sumk + two*r*epi*pion0_source_rate_per_gev(i_pion)
            end if
        end do
        gamma_rate_per_gev(i_gamma) = dln_pi*sumk/eout
    end do
end subroutine hadronic_pi0_to_gamma_operator

! 带电π介子衰变算子：计算μ子（左右手征）和中微子的产生率。
subroutine hadronic_pion_decay_operator(Num_pion,pion_energy_gev,pion_plus_source_rate_per_gev, &
                                        pion_minus_source_rate_per_gev,Num_mu,muon_energy_gev, &
                                        muon_plus_right_source_rate_per_gev,muon_plus_left_source_rate_per_gev, &
                                        muon_minus_left_source_rate_per_gev,muon_minus_right_source_rate_per_gev, &
                                        Num_nu,neutrino_energy_gev,prompt_numu_rate_per_gev, &
                                        prompt_numubar_rate_per_gev)
    integer, intent(in) :: Num_pion,Num_mu,Num_nu
    real(8), intent(in) :: pion_energy_gev(Num_pion), pion_plus_source_rate_per_gev(Num_pion), &
        pion_minus_source_rate_per_gev(Num_pion), muon_energy_gev(Num_mu), neutrino_energy_gev(Num_nu)
    real(8), intent(out) :: muon_plus_right_source_rate_per_gev(Num_mu), muon_plus_left_source_rate_per_gev(Num_mu), &
        muon_minus_left_source_rate_per_gev(Num_mu), muon_minus_right_source_rate_per_gev(Num_mu)
    real(8), intent(out) :: prompt_numu_rate_per_gev(Num_nu), prompt_numubar_rate_per_gev(Num_nu)
    integer :: i_mu, i_nu, i_pion
    real(8) :: dln_pi, rpi, emu, enu, epi, r, fmu1, fmu2, fpi2nu
    real(8) :: rate_pip_log(Num_pion), rate_pim_log(Num_pion), tdec_pi(Num_pion)
    real(8) :: sum_mu_pr, sum_mu_pl, sum_mu_ml, sum_mu_mr, sum_nu_p, sum_nu_m

    muon_plus_right_source_rate_per_gev = zero
    muon_plus_left_source_rate_per_gev = zero
    muon_minus_left_source_rate_per_gev = zero
    muon_minus_right_source_rate_per_gev = zero
    prompt_numu_rate_per_gev = zero
    prompt_numubar_rate_per_gev = zero

    dln_pi = hadronic_log_spacing(Num_pion,pion_energy_gev)
    if (dln_pi <= zero) return
    rpi = (muon_mass_gev/pi_plus_mass_gev)**2

    do i_pion=1,Num_pion
        tdec_pi(i_pion) = (pion_energy_gev(i_pion)/pi_plus_mass_gev)*charged_pion_decay_s
        rate_pip_log(i_pion) = pion_energy_gev(i_pion)*pion_plus_source_rate_per_gev(i_pion)/tdec_pi(i_pion)
        rate_pim_log(i_pion) = pion_energy_gev(i_pion)*pion_minus_source_rate_per_gev(i_pion)/tdec_pi(i_pion)
    end do

    do i_mu=1,Num_mu
        emu = muon_energy_gev(i_mu)
        sum_mu_pr = zero
        sum_mu_pl = zero
        sum_mu_ml = zero
        sum_mu_mr = zero
        do i_pion=1,Num_pion
            epi = pion_energy_gev(i_pion)
            r = emu/epi
            if (r < one .and. r >= rpi) then
                fmu1 = rpi*(one-r)/((one-rpi)**2*r)
                fmu2 = (r-rpi)/((one-rpi)**2*r)
                sum_mu_pr = sum_mu_pr + r*fmu1*rate_pip_log(i_pion)
                sum_mu_pl = sum_mu_pl + r*fmu2*rate_pip_log(i_pion)
                sum_mu_ml = sum_mu_ml + r*fmu1*rate_pim_log(i_pion)
                sum_mu_mr = sum_mu_mr + r*fmu2*rate_pim_log(i_pion)
            end if
        end do
        muon_plus_right_source_rate_per_gev(i_mu) = dln_pi*sum_mu_pr/emu
        muon_plus_left_source_rate_per_gev(i_mu) = dln_pi*sum_mu_pl/emu
        muon_minus_left_source_rate_per_gev(i_mu) = dln_pi*sum_mu_ml/emu
        muon_minus_right_source_rate_per_gev(i_mu) = dln_pi*sum_mu_mr/emu
    end do

    do i_nu=1,Num_nu
        enu = neutrino_energy_gev(i_nu)
        sum_nu_p = zero
        sum_nu_m = zero
        do i_pion=1,Num_pion
            epi = pion_energy_gev(i_pion)
            r = enu/epi
            if (r >= zero .and. r <= one-rpi) then
                fpi2nu = one/(one-rpi)
                sum_nu_p = sum_nu_p + r*fpi2nu*rate_pip_log(i_pion)
                sum_nu_m = sum_nu_m + r*fpi2nu*rate_pim_log(i_pion)
            end if
        end do
        prompt_numu_rate_per_gev(i_nu) = dln_pi*sum_nu_p/enu
        prompt_numubar_rate_per_gev(i_nu) = dln_pi*sum_nu_m/enu
    end do
end subroutine hadronic_pion_decay_operator

! μ子衰变算子：计算电子、正电子和各类中微子的产生率。
subroutine hadronic_muon_decay_operator(Num_mu,muon_energy_gev,muon_plus_right_source_rate_per_gev, &
                                        muon_plus_left_source_rate_per_gev,muon_minus_left_source_rate_per_gev, &
                                        muon_minus_right_source_rate_per_gev,Num_nu,neutrino_energy_gev, &
                                        nu_e_rate_per_gev,nu_ebar_rate_per_gev,nu_mu_rate_per_gev, &
                                        nu_mubar_rate_per_gev,Num_e,electron_energy_gev,electron_minus_rate_per_gev, &
                                        electron_plus_rate_per_gev)
    integer, intent(in) :: Num_mu,Num_nu,Num_e
    real(8), intent(in) :: muon_energy_gev(Num_mu), muon_plus_right_source_rate_per_gev(Num_mu), &
        muon_plus_left_source_rate_per_gev(Num_mu), muon_minus_left_source_rate_per_gev(Num_mu), &
        muon_minus_right_source_rate_per_gev(Num_mu), neutrino_energy_gev(Num_nu), electron_energy_gev(Num_e)
    real(8), intent(out) :: nu_e_rate_per_gev(Num_nu), nu_ebar_rate_per_gev(Num_nu), &
        nu_mu_rate_per_gev(Num_nu), nu_mubar_rate_per_gev(Num_nu), &
        electron_minus_rate_per_gev(Num_e), electron_plus_rate_per_gev(Num_e)
    integer :: i_mu,i_nu,i_e
    real(8) :: dln_mu, enu, ee, emu, x, fnu1_p, fnu1_m, fnu2_p, fnu2_m
    real(8) :: mu_pl_log(Num_mu), mu_pr_log(Num_mu), mu_ml_log(Num_mu), mu_mr_log(Num_mu), tdec_mu(Num_mu)
    real(8) :: sum_nue, sum_nueb, sum_numu, sum_numub, sum_em, sum_ep

    nu_e_rate_per_gev = zero
    nu_ebar_rate_per_gev = zero
    nu_mu_rate_per_gev = zero
    nu_mubar_rate_per_gev = zero
    electron_minus_rate_per_gev = zero
    electron_plus_rate_per_gev = zero

    dln_mu = hadronic_log_spacing(Num_mu,muon_energy_gev)
    if (dln_mu <= zero) return

    do i_mu=1,Num_mu
        tdec_mu(i_mu) = (muon_energy_gev(i_mu)/muon_mass_gev)*muon_decay_s
        mu_pl_log(i_mu) = muon_energy_gev(i_mu)*muon_plus_left_source_rate_per_gev(i_mu)/tdec_mu(i_mu)
        mu_pr_log(i_mu) = muon_energy_gev(i_mu)*muon_plus_right_source_rate_per_gev(i_mu)/tdec_mu(i_mu)
        mu_ml_log(i_mu) = muon_energy_gev(i_mu)*muon_minus_left_source_rate_per_gev(i_mu)/tdec_mu(i_mu)
        mu_mr_log(i_mu) = muon_energy_gev(i_mu)*muon_minus_right_source_rate_per_gev(i_mu)/tdec_mu(i_mu)
    end do

    do i_nu=1,Num_nu
        enu = neutrino_energy_gev(i_nu)
        sum_nue = zero
        sum_nueb = zero
        sum_numu = zero
        sum_numub = zero
        do i_mu=1,Num_mu
            emu = muon_energy_gev(i_mu)
            x = enu/emu
            fnu1_p = hadronic_fnu1_decay(x,one)
            fnu1_m = hadronic_fnu1_decay(x,-one)
            fnu2_p = hadronic_fnu2_decay(x,one)
            fnu2_m = hadronic_fnu2_decay(x,-one)
            sum_nueb = sum_nueb + x*(fnu2_p*mu_ml_log(i_mu) + fnu2_m*mu_mr_log(i_mu))
            sum_nue = sum_nue + x*(fnu2_m*mu_pl_log(i_mu) + fnu2_p*mu_pr_log(i_mu))
            sum_numub = sum_numub + x*(fnu1_m*mu_pl_log(i_mu) + fnu1_p*mu_pr_log(i_mu))
            sum_numu = sum_numu + x*(fnu1_p*mu_ml_log(i_mu) + fnu1_m*mu_mr_log(i_mu))
        end do
        nu_e_rate_per_gev(i_nu) = dln_mu*sum_nue/enu
        nu_ebar_rate_per_gev(i_nu) = dln_mu*sum_nueb/enu
        nu_mu_rate_per_gev(i_nu) = dln_mu*sum_numu/enu
        nu_mubar_rate_per_gev(i_nu) = dln_mu*sum_numub/enu
    end do

    do i_e=1,Num_e
        ee = electron_energy_gev(i_e)
        sum_em = zero
        sum_ep = zero
        do i_mu=1,Num_mu
            emu = muon_energy_gev(i_mu)
            x = ee/emu
            if (x > zero .and. x <= 1.00001d0) then
                sum_ep = sum_ep + x*(4d0/3d0)*(one-x*x*x)*(mu_pl_log(i_mu) + mu_pr_log(i_mu))
                sum_em = sum_em + x*(4d0/3d0)*(one-x*x*x)*(mu_ml_log(i_mu) + mu_mr_log(i_mu))
            end if
        end do
        electron_plus_rate_per_gev(i_e) = dln_mu*sum_ep/ee
        electron_minus_rate_per_gev(i_e) = dln_mu*sum_em/ee
    end do
end subroutine hadronic_muon_decay_operator

! Hummer2010统一衰变算子：整合π0、π±和μ子衰变，输出光子、中微子和带电轻子谱。
subroutine hadronic_hummer2010_decay_operator( &
    Num_gam_p,hadron_energy_gev,pion0_source_rate_per_gev,pion_plus_source_rate_per_gev, &
    pion_minus_source_rate_per_gev,Num_gamma,gamma_energy_gev,Num_nu,neutrino_energy_gev, &
    Num_proc,process_energy_gev,gamma_rate_per_gev,process_rate_per_gev, &
    muon_plus_right_source_rate_per_gev,muon_plus_left_source_rate_per_gev, &
    muon_minus_left_source_rate_per_gev,muon_minus_right_source_rate_per_gev, &
    prompt_pion_neutrino_rate_per_gev,muon_neutrino_rate_per_gev, &
    muon_electron_rate_per_gev,neutrino_rate_per_gev &
)
    integer, intent(in) :: Num_gam_p,Num_gamma,Num_nu,Num_proc
    real(8), intent(in) :: hadron_energy_gev(Num_gam_p), pion0_source_rate_per_gev(Num_gam_p), &
        pion_plus_source_rate_per_gev(Num_gam_p), pion_minus_source_rate_per_gev(Num_gam_p), &
        gamma_energy_gev(Num_gamma), neutrino_energy_gev(Num_nu), process_energy_gev(Num_proc)
    real(8), intent(out) :: gamma_rate_per_gev(Num_gamma), process_rate_per_gev(3,Num_proc)
    real(8), intent(out) :: muon_plus_right_source_rate_per_gev(Num_gam_p), muon_plus_left_source_rate_per_gev(Num_gam_p), &
        muon_minus_left_source_rate_per_gev(Num_gam_p), muon_minus_right_source_rate_per_gev(Num_gam_p)
    real(8), intent(out) :: prompt_pion_neutrino_rate_per_gev(Num_nu), muon_neutrino_rate_per_gev(Num_nu), &
        muon_electron_rate_per_gev(Num_proc), neutrino_rate_per_gev(Num_nu)
    real(8) :: process_gamma_rate_per_gev(Num_proc), prompt_numu_rate_per_gev(Num_nu), &
        prompt_numubar_rate_per_gev(Num_nu), nu_e_rate_per_gev(Num_nu), &
        nu_ebar_rate_per_gev(Num_nu), nu_mu_rate_per_gev(Num_nu), nu_mubar_rate_per_gev(Num_nu), &
        electron_minus_rate_per_gev(Num_proc), electron_plus_rate_per_gev(Num_proc)
    integer :: i_proc

    call hadronic_pi0_to_gamma_operator(Num_gam_p,hadron_energy_gev,pion0_source_rate_per_gev, &
                                        Num_gamma,gamma_energy_gev,gamma_rate_per_gev)
    call hadronic_pi0_to_gamma_operator(Num_gam_p,hadron_energy_gev,pion0_source_rate_per_gev, &
                                        Num_proc,process_energy_gev,process_gamma_rate_per_gev)
    call hadronic_pion_decay_operator( &
        Num_gam_p,hadron_energy_gev,pion_plus_source_rate_per_gev,pion_minus_source_rate_per_gev, &
        Num_gam_p,hadron_energy_gev,muon_plus_right_source_rate_per_gev,muon_plus_left_source_rate_per_gev, &
        muon_minus_left_source_rate_per_gev,muon_minus_right_source_rate_per_gev, &
        Num_nu,neutrino_energy_gev,prompt_numu_rate_per_gev,prompt_numubar_rate_per_gev &
    )
    call hadronic_muon_decay_operator( &
        Num_gam_p,hadron_energy_gev,muon_plus_right_source_rate_per_gev,muon_plus_left_source_rate_per_gev, &
        muon_minus_left_source_rate_per_gev,muon_minus_right_source_rate_per_gev, &
        Num_nu,neutrino_energy_gev,nu_e_rate_per_gev,nu_ebar_rate_per_gev,nu_mu_rate_per_gev, &
        nu_mubar_rate_per_gev,Num_proc,process_energy_gev,electron_minus_rate_per_gev,electron_plus_rate_per_gev &
    )

    prompt_pion_neutrino_rate_per_gev = prompt_numu_rate_per_gev + prompt_numubar_rate_per_gev
    muon_neutrino_rate_per_gev = nu_e_rate_per_gev + nu_ebar_rate_per_gev + nu_mu_rate_per_gev + &
                                 nu_mubar_rate_per_gev
    muon_electron_rate_per_gev = electron_minus_rate_per_gev + electron_plus_rate_per_gev
    neutrino_rate_per_gev = prompt_pion_neutrino_rate_per_gev + muon_neutrino_rate_per_gev

    process_rate_per_gev = zero
    process_rate_per_gev(1,:) = process_gamma_rate_per_gev
    do i_proc=1,Num_proc
        process_rate_per_gev(2,i_proc) = hadronic_log_interpolate(Num_nu,neutrino_energy_gev,prompt_pion_neutrino_rate_per_gev, &
                                                                  process_energy_gev(i_proc))
        process_rate_per_gev(3,i_proc) = hadronic_log_interpolate(Num_nu,neutrino_energy_gev,muon_neutrino_rate_per_gev, &
                                                                  process_energy_gev(i_proc)) + muon_electron_rate_per_gev(i_proc)
    end do
end subroutine hadronic_hummer2010_decay_operator

! 获取能量网格的对数间距（安全返回零当网格点数不足）。
real(8) function hadronic_log_spacing(num_energy,energy_gev)
    integer, intent(in) :: num_energy
    real(8), intent(in) :: energy_gev(num_energy)
    real(8) :: dln_local
    if (num_energy <= 1) then
        hadronic_log_spacing = zero
    else
        call hadronic_validate_log_grid(num_energy,energy_gev,"energy_gev",dln_local)
        hadronic_log_spacing = dln_local
    end if
end function hadronic_log_spacing

! 中微子衰变谱函数 f_ν1(x, h)，x = E_ν/E_μ，h = ±1 为螺旋度。
real(8) function hadronic_fnu1_decay(x,h)
    real(8), intent(in) :: x,h
    if (x >= zero .and. x <= one) then
        hadronic_fnu1_decay = (5d0/3d0 - 3d0*x*x + 4d0*x*x*x/3d0) + h*(-1d0/3d0 + 3d0*x*x - 8d0*x*x*x/3d0)
    else
        hadronic_fnu1_decay = zero
    end if
end function hadronic_fnu1_decay

! 中微子衰变谱函数 f_ν2(x, h)，x = E_ν/E_μ，h = ±1 为螺旋度。
real(8) function hadronic_fnu2_decay(x,h)
    real(8), intent(in) :: x,h
    if (x >= zero .and. x <= one) then
        hadronic_fnu2_decay = (2d0 - 6d0*x*x + 4d0*x*x*x) + h*(2d0 - 12d0*x + 18d0*x*x - 8d0*x*x*x)
    else
        hadronic_fnu2_decay = zero
    end if
end function hadronic_fnu2_decay

! 在对数空间中线性插值：二分查找定位 + 对数权重。
real(8) function hadronic_log_interpolate(num_grid,energy_grid,rate_grid,energy_value)
    integer, intent(in) :: num_grid
    real(8), intent(in) :: energy_grid(num_grid),rate_grid(num_grid),energy_value
    integer :: i_lo,i_hi,i_mid
    real(8) :: weight

    if (num_grid <= 0) then
        hadronic_log_interpolate = zero
        return
    end if
    if (energy_value < energy_grid(1) .or. energy_value > energy_grid(num_grid)) then
        hadronic_log_interpolate = zero
        return
    end if
    if (num_grid == 1) then
        hadronic_log_interpolate = rate_grid(1)
        return
    end if

    i_lo = 1
    i_hi = num_grid
    do while (i_hi - i_lo > 1)
        i_mid = (i_lo + i_hi)/2
        if (energy_grid(i_mid) <= energy_value) then
            i_lo = i_mid
        else
            i_hi = i_mid
        end if
    end do

    if (energy_grid(i_hi) <= energy_grid(i_lo)) then
        hadronic_log_interpolate = rate_grid(i_lo)
        return
    end if

    weight = dlog(energy_value/energy_grid(i_lo))/dlog(energy_grid(i_hi)/energy_grid(i_lo))
    hadronic_log_interpolate = (one-weight)*rate_grid(i_lo) + weight*rate_grid(i_hi)
end function hadronic_log_interpolate

end module hadronic_decay_kernel
