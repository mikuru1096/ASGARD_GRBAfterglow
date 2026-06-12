!f2py: skip
module hadronic_interaction_kernel
    use constants
    use hadronic_common, only: hadronic_proton_mass_gev
    implicit none
    real(8), parameter :: microbarn_to_cm2 = 1.0d-30

    real(8), parameter :: it1_er(28) = (/ &
        0.34d0, 0.64d0, 0.64d0, 0.64d0, 0.75d0, 0.75d0, 0.75d0, 0.75d0, 0.77d0, 1.03d0, 1.03d0, 1.03d0, &
        1.04d0, 1.04d0, 1.04d0, 1.04d0, 1.05d0, 1.05d0, 1.05d0, 1.05d0, 1.45d0, 1.45d0, 1.45d0, 1.45d0, &
        1.56d0, 1.56d0, 1.56d0, 1.56d0 /)
    real(8), parameter :: it1_bout(28) = (/ &
        1.00d0, 0.67d0, 0.33d0, 0.33d0, 0.52d0, 0.42d0, 0.42d0, 0.06d0, 0.45d0, 0.75d0, 0.14d0, 0.14d0, &
        0.64d0, 0.22d0, 0.22d0, 0.14d0, 0.14d0, 0.55d0, 0.55d0, 0.31d0, 0.14d0, 0.13d0, 0.13d0, 0.73d0, &
        0.37d0, 0.39d0, 0.39d0, 0.24d0 /)
    real(8), parameter :: it1_g(28) = (/ &
        66.69d0, 4.260d0, 4.260d0, 4.260d0, 27.01d0, 27.01d0, 27.01d0, 27.01d0, 6.590d0, 3.190d0, 3.190d0, 3.190d0, &
        15.72d0, 15.72d0, 15.72d0, 15.72d0, 21.68d0, 21.68d0, 21.68d0, 21.68d0, 3.030d0, 3.030d0, 3.030d0, 3.030d0, &
        16.46d0, 16.46d0, 16.46d0, 16.46d0 /)
    real(8), parameter :: it1_chi(28) = (/ &
        0.22d0, 0.29d0, 0.14d0, 0.19d0, 0.31d0, 0.17d0, 0.18d0, 0.22d0, 0.32d0, 0.35d0, 0.23d0, 0.17d0, &
        0.35d0, 0.24d0, 0.17d0, 0.23d0, 0.35d0, 0.24d0, 0.16d0, 0.23d0, 0.38d0, 0.29d0, 0.15d0, 0.23d0, &
        0.39d0, 0.30d0, 0.15d0, 0.23d0 /)
    real(8), parameter :: it1_chi_pn(28) = (/ &
        0.78d0, 0.71d0, 1.00d0, 0.67d0, 0.69d0, 1.00d0, 0.65d0, 0.56d0, 0.68d0, 0.65d0, 1.00d0, 0.60d0, &
        0.65d0, 1.00d0, 0.59d0, 0.54d0, 0.65d0, 1.00d0, 0.60d0, 0.54d0, 0.62d0, 1.00d0, 0.56d0, 0.54d0, &
        0.61d0, 1.00d0, 0.55d0, 0.54d0 /)
    real(8), parameter :: it1_mp_pi0(28) = (/ &
        0.667d0, 0.333d0, 0.333d0, 0.333d0, 0.333d0, 0.333d0, 0.333d0, 0.667d0, 0.333d0, 0.333d0, 0.333d0, 0.333d0, &
        0.333d0, 0.333d0, 0.333d0, 0.667d0, 0.667d0, 0.067d0, 0.400d0, 0.333d0, 0.667d0, 0.067d0, 0.400d0, 0.333d0, &
        0.667d0, 0.067d0, 0.400d0, 0.333d0 /)
    real(8), parameter :: it1_mp_pip(28) = (/ &
        0.333d0, 0.667d0, 0.167d0, 0.611d0, 0.667d0, 0.167d0, 0.611d0, 1.000d0, 0.667d0, 0.667d0, 0.167d0, 0.611d0, &
        0.667d0, 0.167d0, 0.611d0, 1.000d0, 0.333d0, 0.533d0, 0.422d0, 1.000d0, 0.333d0, 0.533d0, 0.422d0, 1.000d0, &
        0.333d0, 0.533d0, 0.422d0, 1.000d0 /)
    real(8), parameter :: it1_mp_pim(28) = (/ &
        0.000d0, 0.000d0, 0.500d0, 0.056d0, 0.000d0, 0.500d0, 0.056d0, 0.333d0, 0.000d0, 0.000d0, 0.500d0, 0.056d0, &
        0.000d0, 0.500d0, 0.056d0, 0.333d0, 0.000d0, 0.400d0, 0.178d0, 0.667d0, 0.000d0, 0.400d0, 0.178d0, 0.667d0, &
        0.000d0, 0.400d0, 0.178d0, 0.667d0 /)
    real(8), parameter :: it1_mp_pro(28) = (/ &
        0.667d0, 0.333d0, 0.000d0, 0.778d0, 0.333d0, 0.000d0, 0.778d0, 0.333d0, 0.333d0, 0.333d0, 0.000d0, 0.778d0, &
        0.333d0, 0.000d0, 0.778d0, 0.333d0, 0.667d0, 0.000d0, 0.622d0, 0.667d0, 0.667d0, 0.000d0, 0.622d0, 0.667d0, &
        0.667d0, 0.000d0, 0.622d0, 0.667d0 /)
    real(8), parameter :: it1_mp_ntr(28) = (/ &
        0.333d0, 0.667d0, 0.000d0, 0.222d0, 0.667d0, 0.000d0, 0.222d0, 0.667d0, 0.667d0, 0.667d0, 0.000d0, 0.222d0, &
        0.667d0, 0.000d0, 0.222d0, 0.667d0, 0.333d0, 0.000d0, 0.378d0, 0.333d0, 0.333d0, 0.000d0, 0.378d0, 0.333d0, &
        0.333d0, 0.000d0, 0.378d0, 0.333d0 /)
    real(8), parameter :: it2_em(9) = (/ 0.17d0, 0.56d0, 10.0d0, 0.4d0, 1.58d0, 10.0d0, 0.4d0, 1.58d0, 10.0d0 /)
    real(8), parameter :: it2_ex(9) = (/ 0.56d0, 10.0d0, 1.0d6, 1.58d0, 10.0d0, 1.0d6, 1.58d0, 10.0d0, 1.0d6 /)
    real(8), parameter :: it2_chi(9) = (/ 0.13d0, 0.05d0, 0.001d0, 0.08d0, 0.02d0, 0.001d0, 0.2d0, 0.2d0, 0.2d0 /)
    real(8), parameter :: it2_chi_pn(9) = (/ 0.87d0, 0.95d0, 0.999d0, 0.72d0, 0.78d0, 0.799d0, 1.00d0, 1.00d0, 1.00d0 /)
    real(8), parameter :: it2_mp_pi0(9) = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.1667d0, 0.1667d0, 0.1667d0 /)
    real(8), parameter :: it2_mp_pip(9) = (/ 1.0d0, 1.0d0, 1.0d0, 0.25d0, 0.25d0, 0.25d0, 0.75d0, 0.75d0, 0.75d0 /)
    real(8), parameter :: it2_mp_pim(9) = (/ 0.0d0, 0.0d0, 0.0d0, 0.75d0, 0.75d0, 0.75d0, 0.0833d0, 0.0833d0, 0.0833d0 /)
    real(8), parameter :: it2_mp_pro(9) = (/ 0.000d0, 0.000d0, 0.000d0, 0.833d0, 0.833d0, 0.833d0, 0.0d0, 0.0d0, 0.0d0 /)
    real(8), parameter :: it2_mp_ntr(9) = (/ 1.000d0, 1.000d0, 1.000d0, 0.167d0, 0.167d0, 0.167d0, 0.0d0, 0.0d0, 0.0d0 /)
    real(8), parameter :: it3_sig(14) = (/ &
        60.0d0, 60.0d0, 85.0d0, 85.0d0, 120.0d0, 120.0d0, 120.0d0, 120.0d0, &
        120.0d0, 120.0d0, 120.0d0, 120.0d0, 120.0d0, 120.0d0 /)
    real(8), parameter :: it3_em(14) = (/ &
        0.5d0, 0.5d0, 0.9d0, 0.9d0, 1.5d0, 1.5d0, 5.0d0, 5.0d0, 50.0d0, 50.0d0, &
        500.0d0, 500.0d0, 5000.0d0, 5000.0d0 /)
    real(8), parameter :: it3_ex(14) = (/ &
        0.9d0, 0.9d0, 1.5d0, 1.5d0, 5.0d0, 5.0d0, 50.0d0, 50.0d0, 500.0d0, 500.0d0, &
        5000.0d0, 5000.0d0, 5.0d6, 5.0d6 /)
    real(8), parameter :: it3_chi(14) = (/ &
        0.1d0, 0.4d0, 0.15d0, 0.35d0, 0.15d0, 0.35d0, 0.07d0, 0.35d0, 0.02d0, 0.5d0, &
        0.007d0, 0.5d0, 0.002d0, 0.6d0 /)
    real(8), parameter :: it3_chi_pn(14) = (/ &
        0.73d0, 1.00d0, 0.66d0, 1.00d0, 0.61d0, 1.00d0, 0.51d0, 1.00d0, 0.55d0, 1.00d0, &
        0.56d0, 1.00d0, 0.56d0, 1.00d0 /)
    real(8), parameter :: it3_mp_pi0(14) = (/ &
        0.32d0, 0.17d0, 0.42d0, 0.19d0, 0.59d0, 0.16d0, 1.38d0, 0.16d0, 3.01d0, 0.20d0, &
        5.13d0, 0.27d0, 7.59d0, 0.26d0 /)
    real(8), parameter :: it3_mp_pip(14) = (/ &
        0.34d0, 0.29d0, 0.31d0, 0.35d0, 0.57d0, 0.21d0, 1.37d0, 0.25d0, 2.86d0, 0.21d0, &
        4.68d0, 0.29d0, 6.80d0, 0.27d0 /)
    real(8), parameter :: it3_mp_pim(14) = (/ &
        0.04d0, 0.05d0, 0.07d0, 0.08d0, 0.30d0, 0.13d0, 1.11d0, 0.23d0, 2.64d0, 0.14d0, &
        4.57d0, 0.12d0, 6.65d0, 0.13d0 /)
    real(8), parameter :: it3_mp_pro(14) = (/ &
        0.46d0, 0.00d0, 0.49d0, 0.00d0, 0.65d0, 0.00d0, 0.72d0, 0.00d0, 0.71d0, 0.00d0, &
        0.72d0, 0.00d0, 0.71d0, 0.00d0 /)
    real(8), parameter :: it3_mp_ntr(14) = (/ &
        0.54d0, 0.00d0, 0.51d0, 0.00d0, 0.35d0, 0.00d0, 0.28d0, 0.00d0, 0.29d0, 0.00d0, &
        0.28d0, 0.00d0, 0.29d0, 0.00d0 /)

contains

! Hummer2010光介子产生主算子：计算pγ相互作用产生的π介子、重子再注入及粒子损失率。
subroutine hadronic_pg_hummer2010_operator(Num_gam_p,Num_nu,hadron_energy_gev,hadron_density_per_gev,photon_energy_gev, &
                                           photon_density_per_gev,neutron_density_per_gev,pion0_source_rate_per_gev, &
                                           pion_plus_source_rate_per_gev,pion_minus_source_rate_per_gev, &
                                           proton_reinjection_rate_per_gev,neutron_reinjection_rate_per_gev, &
                                           proton_loss_rate,neutron_loss_rate,photon_loss_rate)
    integer, intent(in) :: Num_gam_p,Num_nu
    real(8), intent(in) :: hadron_energy_gev(Num_gam_p), hadron_density_per_gev(Num_gam_p), &
        photon_energy_gev(Num_nu), photon_density_per_gev(Num_nu), neutron_density_per_gev(Num_gam_p)
    real(8), intent(out) :: pion0_source_rate_per_gev(Num_gam_p), pion_plus_source_rate_per_gev(Num_gam_p), &
        pion_minus_source_rate_per_gev(Num_gam_p), proton_reinjection_rate_per_gev(Num_gam_p), &
        neutron_reinjection_rate_per_gev(Num_gam_p), proton_loss_rate(Num_gam_p), &
        neutron_loss_rate(Num_gam_p), photon_loss_rate(Num_nu)
    integer :: i_gam
    real(8) :: dln_ep, dln_eg, gamma_p, src_p, src_n, ep
    real(8) :: proton_log_density(Num_gam_p), neutron_log_density(Num_gam_p), photon_log_density(Num_nu)
    real(8) :: ph_conv_res(28), ph_conv_dir(9), ph_conv_mul(14)

    dln_ep = dlog(hadron_energy_gev(2)/hadron_energy_gev(1))
    dln_eg = dlog(photon_energy_gev(2)/photon_energy_gev(1))
    proton_log_density = hadron_energy_gev * hadron_density_per_gev
    neutron_log_density = hadron_energy_gev * neutron_density_per_gev
    photon_log_density = photon_energy_gev * photon_density_per_gev

    pion0_source_rate_per_gev = zero
    pion_plus_source_rate_per_gev = zero
    pion_minus_source_rate_per_gev = zero
    proton_reinjection_rate_per_gev = zero
    neutron_reinjection_rate_per_gev = zero
    proton_loss_rate = zero
    neutron_loss_rate = zero
    photon_loss_rate = zero

    do i_gam=1,Num_gam_p
        src_p = proton_log_density(i_gam)
        src_n = neutron_log_density(i_gam)
        if (src_p <= zero .and. src_n <= zero) cycle
        ep = hadron_energy_gev(i_gam)
        gamma_p = ep / hadronic_proton_mass_gev
        call hadronic_pg_family_rates_res(gamma_p,Num_nu,photon_energy_gev,photon_log_density,dln_eg,ph_conv_res)
        call hadronic_pg_family_rates_dir(gamma_p,Num_nu,photon_energy_gev,photon_log_density,dln_eg,ph_conv_dir)
        call hadronic_pg_family_rates_mul(gamma_p,Num_nu,photon_energy_gev,photon_log_density,dln_eg,ph_conv_mul)

        call hadronic_pg_deposit_family(28,i_gam,dln_ep,src_p,src_n,it1_chi,it1_mp_pi0,it1_mp_pip,it1_mp_pim, &
                                        it1_mp_pi0,it1_mp_pim,it1_mp_pip,ph_conv_res,pion0_source_rate_per_gev, &
                                        pion_plus_source_rate_per_gev,pion_minus_source_rate_per_gev)
        call hadronic_pg_deposit_family(9,i_gam,dln_ep,src_p,src_n,it2_chi,it2_mp_pi0,it2_mp_pip,it2_mp_pim, &
                                        it2_mp_pi0,it2_mp_pim,it2_mp_pip,ph_conv_dir,pion0_source_rate_per_gev, &
                                        pion_plus_source_rate_per_gev,pion_minus_source_rate_per_gev)
        call hadronic_pg_deposit_family(14,i_gam,dln_ep,src_p,src_n,it3_chi,it3_mp_pi0,it3_mp_pip,it3_mp_pim, &
                                        it3_mp_pi0,it3_mp_pim,it3_mp_pip,ph_conv_mul,pion0_source_rate_per_gev, &
                                        pion_plus_source_rate_per_gev,pion_minus_source_rate_per_gev)

        call hadronic_pg_deposit_baryons(28,i_gam,dln_ep,src_p,src_n,it1_chi_pn,it1_mp_pro,it1_mp_ntr,it1_mp_ntr,it1_mp_pro, &
                                         ph_conv_res,proton_reinjection_rate_per_gev,neutron_reinjection_rate_per_gev)
        call hadronic_pg_deposit_baryons(9,i_gam,dln_ep,src_p,src_n,it2_chi_pn,it2_mp_pro,it2_mp_ntr,it2_mp_ntr,it2_mp_pro, &
                                         ph_conv_dir,proton_reinjection_rate_per_gev,neutron_reinjection_rate_per_gev)
        call hadronic_pg_deposit_baryons(14,i_gam,dln_ep,src_p,src_n,it3_chi_pn,it3_mp_pro,it3_mp_ntr,it3_mp_ntr,it3_mp_pro, &
                                         ph_conv_mul,proton_reinjection_rate_per_gev,neutron_reinjection_rate_per_gev)

        proton_loss_rate(i_gam) = hadronic_pg_loss_coeff_res(ph_conv_res) + hadronic_pg_loss_coeff_dir(ph_conv_dir) + &
                                  hadronic_pg_loss_coeff_mul(ph_conv_mul)
        neutron_loss_rate(i_gam) = proton_loss_rate(i_gam)

        call accumulate_pg_photon_loss(gamma_p,src_p+src_n)
    end do

    pion0_source_rate_per_gev = pion0_source_rate_per_gev / hadron_energy_gev
    pion_plus_source_rate_per_gev = pion_plus_source_rate_per_gev / hadron_energy_gev
    pion_minus_source_rate_per_gev = pion_minus_source_rate_per_gev / hadron_energy_gev
    proton_reinjection_rate_per_gev = proton_reinjection_rate_per_gev / hadron_energy_gev
    neutron_reinjection_rate_per_gev = neutron_reinjection_rate_per_gev / hadron_energy_gev

contains

    subroutine accumulate_pg_photon_loss(gamma_p,source_log_density)
        real(8), intent(in) :: gamma_p,source_log_density
        integer :: i_nu,it
        real(8) :: y

        do i_nu=1,Num_nu
            y = gamma_p * photon_energy_gev(i_nu)
            do it=1,28
                if (it1_mp_pro(it) > 1d-4 .or. it1_mp_ntr(it) > 1d-4) then
                    photon_loss_rate(i_nu) = photon_loss_rate(i_nu) + &
                                             source_log_density * hadronic_pg_family_photon_loss_res(it,y,dln_ep)
                end if
            end do
            do it=1,9
                if (it2_mp_pro(it) > 1d-4 .or. it2_mp_ntr(it) > 1d-4) then
                    photon_loss_rate(i_nu) = photon_loss_rate(i_nu) + &
                                             source_log_density * hadronic_pg_family_photon_loss_dir(it,y,dln_ep)
                end if
            end do
            do it=1,14
                if (it3_mp_pro(it) > 1d-4 .or. it3_mp_ntr(it) > 1d-4) then
                    photon_loss_rate(i_nu) = photon_loss_rate(i_nu) + &
                                             source_log_density * hadronic_pg_family_photon_loss_mul(it,y,dln_ep)
                end if
            end do
        end do
    end subroutine accumulate_pg_photon_loss
end subroutine hadronic_pg_hummer2010_operator

! 将pγ相互作用产生的π介子多重数按能量偏移χ分配到网格上。
subroutine hadronic_pg_deposit_family(nfam,i_gam,dln_ep,src_p,src_n,chi_arr,mult_pi0_p,mult_pip_p,mult_pim_p, &
                                      mult_pi0_n,mult_pip_n,mult_pim_n,ph_conv,qpi0,qpip,qpim)
    integer, intent(in) :: nfam,i_gam
    real(8), intent(in) :: dln_ep,src_p,src_n,chi_arr(nfam),mult_pi0_p(nfam),mult_pip_p(nfam),mult_pim_p(nfam)
    real(8), intent(in) :: mult_pi0_n(nfam),mult_pip_n(nfam),mult_pim_n(nfam),ph_conv(nfam)
    real(8), intent(inout) :: qpi0(:),qpip(:),qpim(:)
    integer :: it
    do it=1,nfam
        if (ph_conv(it) <= zero) cycle
        call hadronic_pg_deposit_shifted(qpi0,i_gam,chi_arr(it),dln_ep, &
                                         ph_conv(it)*(src_p*mult_pi0_p(it)+src_n*mult_pi0_n(it)))
        call hadronic_pg_deposit_shifted(qpip,i_gam,chi_arr(it),dln_ep, &
                                         ph_conv(it)*(src_p*mult_pip_p(it)+src_n*mult_pip_n(it)))
        call hadronic_pg_deposit_shifted(qpim,i_gam,chi_arr(it),dln_ep, &
                                         ph_conv(it)*(src_p*mult_pim_p(it)+src_n*mult_pim_n(it)))
    end do
end subroutine hadronic_pg_deposit_family

! 将pγ相互作用产生的重子（质子/中子）按能量偏移χ分配到网格上。
subroutine hadronic_pg_deposit_baryons(nfam,i_gam,dln_ep,src_p,src_n,chi_arr,mult_pro_p,mult_ntr_p,mult_pro_n,mult_ntr_n, &
                                       ph_conv,qpro,qntr)
    integer, intent(in) :: nfam,i_gam
    real(8), intent(in) :: dln_ep,src_p,src_n,chi_arr(nfam),mult_pro_p(nfam),mult_ntr_p(nfam)
    real(8), intent(in) :: mult_pro_n(nfam),mult_ntr_n(nfam),ph_conv(nfam)
    real(8), intent(inout) :: qpro(:),qntr(:)
    integer :: it
    do it=1,nfam
        if (ph_conv(it) <= zero) cycle
        call hadronic_pg_deposit_shifted(qpro,i_gam,chi_arr(it),dln_ep, &
                                         ph_conv(it)*(src_p*mult_pro_p(it)+src_n*mult_pro_n(it)))
        call hadronic_pg_deposit_shifted(qntr,i_gam,chi_arr(it),dln_ep, &
                                         ph_conv(it)*(src_p*mult_ntr_p(it)+src_n*mult_ntr_n(it)))
    end do
end subroutine hadronic_pg_deposit_baryons

! 将权重点按对数能量偏移分配到目标数组的相邻网格点上（线性插值沉积）。
subroutine hadronic_pg_deposit_shifted(target,i_gam,chi,dln_ep,weight)
    real(8), intent(inout) :: target(:)
    integer, intent(in) :: i_gam
    real(8), intent(in) :: chi,dln_ep,weight
    integer :: ileft,iright
    real(8) :: pos,frac
    if (weight == zero .or. chi <= zero) return
    pos = dble(i_gam) + dlog(chi)/dln_ep
    ileft = floor(pos)
    frac = pos - dble(ileft)
    iright = ileft + 1
    if (ileft >= 1 .and. ileft <= size(target)) target(ileft) = target(ileft) + (one-frac)*weight
    if (iright >= 1 .and. iright <= size(target)) target(iright) = target(iright) + frac*weight
end subroutine hadronic_pg_deposit_shifted

! 计算共振区（Δ共振）pγ相互作用率，对光子场卷积。
subroutine hadronic_pg_family_rates_res(gamma_p,Num_nu,photon_energy_gev,photon_log_density,dln_eg,ph_conv)
    integer, intent(in) :: Num_nu
    real(8), intent(in) :: gamma_p,photon_energy_gev(Num_nu),photon_log_density(Num_nu),dln_eg
    real(8), intent(out) :: ph_conv(28)
    integer :: it,i_nu
    real(8) :: y,sumk
    if (gamma_p <= zero) error stop "hadronic_pg_family_rates_res: gamma_p must be positive."
    do it=1,28
        sumk = zero
        do i_nu=1,Num_nu
            y = gamma_p * photon_energy_gev(i_nu)
            if (photon_energy_gev(i_nu) <= zero) error stop "hadronic_pg_family_rates_res: photon energy grid must be positive."
            if (two*y >= it1_er(it)) sumk = sumk + photon_log_density(i_nu) / (two*y*y)
        end do
        ph_conv(it) = Para_c * microbarn_to_cm2 * dln_eg * (it1_er(it)*it1_bout(it)*it1_g(it)) * sumk
    end do
end subroutine hadronic_pg_family_rates_res

! 计算直接多π产生区pγ相互作用率，对光子场卷积。
subroutine hadronic_pg_family_rates_dir(gamma_p,Num_nu,photon_energy_gev,photon_log_density,dln_eg,ph_conv)
    integer, intent(in) :: Num_nu
    real(8), intent(in) :: gamma_p,photon_energy_gev(Num_nu),photon_log_density(Num_nu),dln_eg
    real(8), intent(out) :: ph_conv(9)
    integer :: it,i_nu
    real(8) :: y,sumk
    do it=1,9
        sumk = zero
        do i_nu=1,Num_nu
            y = gamma_p * photon_energy_gev(i_nu)
            sumk = sumk + photon_log_density(i_nu) * hadronic_pg_kernel_dir(it,y)
        end do
        ph_conv(it) = Para_c * microbarn_to_cm2 * dln_eg * sumk
    end do
end subroutine hadronic_pg_family_rates_dir

! 计算多π产生区（高能）pγ相互作用率，对光子场卷积。
subroutine hadronic_pg_family_rates_mul(gamma_p,Num_nu,photon_energy_gev,photon_log_density,dln_eg,ph_conv)
    integer, intent(in) :: Num_nu
    real(8), intent(in) :: gamma_p,photon_energy_gev(Num_nu),photon_log_density(Num_nu),dln_eg
    real(8), intent(out) :: ph_conv(14)
    integer :: it,i_nu
    real(8) :: y,sumk
    do it=1,14
        sumk = zero
        do i_nu=1,Num_nu
            y = gamma_p * photon_energy_gev(i_nu)
            sumk = sumk + photon_log_density(i_nu) * hadronic_pg_kernel_mul(it,y)
        end do
        ph_conv(it) = Para_c * microbarn_to_cm2 * dln_eg * it3_sig(it) * sumk
    end do
end subroutine hadronic_pg_family_rates_mul

! 共振区质子/中子能量损失系数（求和有效通道）。
real(8) function hadronic_pg_loss_coeff_res(ph_conv)
    real(8), intent(in) :: ph_conv(28)
    integer :: it
    hadronic_pg_loss_coeff_res = zero
    do it=1,28
        if (it1_mp_pro(it) > 1d-4 .or. it1_mp_ntr(it) > 1d-4) hadronic_pg_loss_coeff_res = hadronic_pg_loss_coeff_res + ph_conv(it)
    end do
end function hadronic_pg_loss_coeff_res

! 直接多π区质子/中子能量损失系数。
real(8) function hadronic_pg_loss_coeff_dir(ph_conv)
    real(8), intent(in) :: ph_conv(9)
    integer :: it
    hadronic_pg_loss_coeff_dir = zero
    do it=1,9
        if (it2_mp_pro(it) > 1d-4 .or. it2_mp_ntr(it) > 1d-4) hadronic_pg_loss_coeff_dir = hadronic_pg_loss_coeff_dir + ph_conv(it)
    end do
end function hadronic_pg_loss_coeff_dir

! 高能多π区质子/中子能量损失系数。
real(8) function hadronic_pg_loss_coeff_mul(ph_conv)
    real(8), intent(in) :: ph_conv(14)
    integer :: it
    hadronic_pg_loss_coeff_mul = zero
    do it=1,14
        if (it3_mp_pro(it) > 1d-4 .or. it3_mp_ntr(it) > 1d-4) then
            hadronic_pg_loss_coeff_mul = hadronic_pg_loss_coeff_mul + ph_conv(it)
        end if
    end do
end function hadronic_pg_loss_coeff_mul

! 共振区光子损失率核函数。
real(8) function hadronic_pg_family_photon_loss_res(it,y,dln_ep)
    integer, intent(in) :: it
    real(8), intent(in) :: y,dln_ep
    hadronic_pg_family_photon_loss_res = Para_c * microbarn_to_cm2 * dln_ep * &
                                         (it1_er(it)*it1_bout(it)*it1_g(it)) * hadronic_pg_kernel_res(it,y)
end function hadronic_pg_family_photon_loss_res

! 直接多π区光子损失率核函数。
real(8) function hadronic_pg_family_photon_loss_dir(it,y,dln_ep)
    integer, intent(in) :: it
    real(8), intent(in) :: y,dln_ep
    hadronic_pg_family_photon_loss_dir = Para_c * microbarn_to_cm2 * dln_ep * hadronic_pg_kernel_dir(it,y)
end function hadronic_pg_family_photon_loss_dir

! 高能多π区光子损失率核函数。
real(8) function hadronic_pg_family_photon_loss_mul(it,y,dln_ep)
    integer, intent(in) :: it
    real(8), intent(in) :: y,dln_ep
    hadronic_pg_family_photon_loss_mul = Para_c * microbarn_to_cm2 * dln_ep * &
                                         it3_sig(it) * hadronic_pg_kernel_mul(it,y)
end function hadronic_pg_family_photon_loss_mul

! 共振区pγ截面核：阶跃函数形式 σ ∝ 1/y^2，阈值为 ε_r/2。
real(8) function hadronic_pg_kernel_res(it,y)
    integer, intent(in) :: it
    real(8), intent(in) :: y
    if (y <= zero) error stop "hadronic_pg_kernel_res: y must be positive."
    if (two*y >= it1_er(it)) then
        hadronic_pg_kernel_res = one / (two * y * y)
    else
        hadronic_pg_kernel_res = zero
    end if
end function hadronic_pg_kernel_res

! 高能多π区pγ截面核：分段函数形式，含常数平台和1/y^2尾部。
real(8) function hadronic_pg_kernel_mul(it,y)
    integer, intent(in) :: it
    real(8), intent(in) :: y
    if (y <= zero) error stop "hadronic_pg_kernel_mul: y must be positive."
    if (y < 0.5d0*it3_em(it)) then
        hadronic_pg_kernel_mul = zero
    else if (y < 0.5d0*it3_ex(it)) then
        hadronic_pg_kernel_mul = one - it3_em(it)*it3_em(it)/(4d0*y*y)
    else
        hadronic_pg_kernel_mul = (it3_ex(it)*it3_ex(it)-it3_em(it)*it3_em(it))/(4d0*y*y)
    end if
end function hadronic_pg_kernel_mul

! 直接多π区pγ截面核：通过分段单调函数 I(z) 的差分计算。
real(8) function hadronic_pg_kernel_dir(it,y)
    integer, intent(in) :: it
    real(8), intent(in) :: y
    real(8) :: zz
    if (y <= zero) error stop "hadronic_pg_kernel_dir: y must be positive."
    if (y < 0.5d0*it2_em(it)) then
        hadronic_pg_kernel_dir = zero
    else if (y < 0.5d0*it2_ex(it)) then
        zz = two*y
        hadronic_pg_kernel_dir = (hadronic_pg_idir(it,zz)-hadronic_pg_idir(it,it2_em(it))) / (two*y*y)
    else
        hadronic_pg_kernel_dir = (hadronic_pg_idir(it,it2_ex(it))-hadronic_pg_idir(it,it2_em(it))) / (two*y*y)
    end if
end function hadronic_pg_kernel_dir

! 直接多π区截面积分函数 I(z)：低z用四次多项式，高z用解析近似。
real(8) function hadronic_pg_idir(it,z)
    integer, intent(in) :: it
    real(8), intent(in) :: z
    real(8) :: x
    hadronic_pg_idir = zero
    if (z <= zero) error stop "hadronic_pg_idir: z must be positive."
    x = dlog10(0.5d0*z)
    if (it <= 3) then
        if (z >= 0.17d0 .and. z < 0.96d0) then
            hadronic_pg_idir = 35.95d0 + 84.08d0*x + 110.76d0*x*x + 102.73d0*x*x*x + 40.47d0*x*x*x*x
        else if (z >= 0.96d0) then
            hadronic_pg_idir = 30.20d0 + 40.55d0*x + 2.031d0*x*x - 0.3879d0*x*x*x + 0.02504d0*x*x*x*x
        end if
    else
        if (z >= 0.4d0) hadronic_pg_idir = -3.4083d0 + 16.2864d0/z + 40.7160d0*dlog(z)
    end if
end function hadronic_pg_idir

end module hadronic_interaction_kernel
