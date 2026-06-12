!f2py: skip
module hadronic_bethe_heitler_kernel
    use constants
    use hadronic_common, only: hadronic_proton_mass_gev, hadronic_electron_mass_gev, &
                               hadronic_validate_log_grid
    implicit none
    private

    real(8), parameter :: fine_structure_alpha = 1.0d0/137.0d0
    integer, parameter :: bh_outer_bins = 5
    integer, parameter :: bh_inner_bins = 5

    public :: hadronic_bethe_heitler_operator

contains

! Bethe-Heitler过程主算子：计算质子-光子碰撞产生的正负电子对注入率和质子能量损失率。
subroutine hadronic_bethe_heitler_operator(num_p,proton_energy_gev,proton_density_per_gev,num_ph,photon_energy_gev, &
                                           photon_density_per_gev,num_e,electron_energy_gev,pair_rate_per_gev, &
                                           proton_loss_rate,photon_loss_rate)
    integer, intent(in) :: num_p,num_ph,num_e
    real(8), intent(in) :: proton_energy_gev(num_p),proton_density_per_gev(num_p)
    real(8), intent(in) :: photon_energy_gev(num_ph),photon_density_per_gev(num_ph)
    real(8), intent(in) :: electron_energy_gev(num_e)
    real(8), intent(out) :: pair_rate_per_gev(num_e),proton_loss_rate(num_p),photon_loss_rate(num_ph)
    real(8) :: dln_ep,dln_eph,dln_ee
    real(8) :: proton_log_density(num_p),photon_log_density(num_ph)
    real(8) :: pair_log_source(num_e),photon_log_loss(num_ph)
    real(8) :: gp_arr(num_p),eph_dimless(num_ph),ee_dimless(num_e)
    real(8) :: k_inj_rate,k_loss_rate
    integer :: j_p,k_ph

    call hadronic_validate_log_grid(num_p,proton_energy_gev,"proton_energy_gev")
    call hadronic_validate_log_grid(num_ph,photon_energy_gev,"photon_energy_gev")
    call hadronic_validate_log_grid(num_e,electron_energy_gev,"electron_energy_gev")

    dln_ep = hadronic_uniform_log_spacing(num_p,proton_energy_gev,"proton_energy_gev")
    dln_eph = hadronic_uniform_log_spacing(num_ph,photon_energy_gev,"photon_energy_gev")
    dln_ee = hadronic_uniform_log_spacing(num_e,electron_energy_gev,"electron_energy_gev")

    proton_log_density = proton_energy_gev*proton_density_per_gev
    photon_log_density = photon_energy_gev*photon_density_per_gev
    gp_arr = proton_energy_gev/hadronic_proton_mass_gev
    eph_dimless = photon_energy_gev/hadronic_electron_mass_gev
    ee_dimless = electron_energy_gev/hadronic_electron_mass_gev

    k_inj_rate = para_c*3.0d0*Para_SigmaT*fine_structure_alpha/(16.0d0*pi)
    k_loss_rate = (3.0d0/(8.0d0*pi))*fine_structure_alpha*Para_SigmaT*para_c &
        *(hadronic_electron_mass_gev/hadronic_proton_mass_gev)

    pair_log_source = zero
    proton_loss_rate = zero
    photon_loss_rate = zero
    photon_log_loss = zero

    do j_p=1,num_p
        do k_ph=1,num_ph
            proton_loss_rate(j_p) = proton_loss_rate(j_p) + &
                                    bh_proton_loss_point(gp_arr(j_p),eph_dimless(k_ph),photon_log_density(k_ph))
            if (2.0d0*gp_arr(j_p)*eph_dimless(k_ph) <= 2.0d0) cycle
            call accumulate_bh_pair_source(k_ph,gp_arr(j_p),proton_log_density(j_p), &
                                           eph_dimless(k_ph),photon_log_density(k_ph))
        end do
    end do

    pair_log_source = k_inj_rate*dln_ep*dln_eph*pair_log_source
    pair_rate_per_gev = pair_log_source/electron_energy_gev
    photon_log_loss = k_inj_rate*dln_ep*dln_ee*photon_log_loss
    where (photon_log_density > zero)
        photon_loss_rate = photon_log_loss/photon_log_density
    elsewhere
        photon_loss_rate = zero
    end where

contains

    real(8) function bh_proton_loss_point(gp,eph,photon_log_value)
        real(8), intent(in) :: gp,eph,photon_log_value

        bh_proton_loss_point = -k_loss_rate*dln_eph*hadronic_bh_kernel_proton_loss(gp,eph)*photon_log_value
    end function bh_proton_loss_point

    subroutine accumulate_bh_pair_source(i_ph,gp,proton_log_value,eph,photon_log_value)
        integer, intent(in) :: i_ph
        real(8), intent(in) :: gp,proton_log_value,eph,photon_log_value
        integer :: i_e
        real(8) :: kernel_value

        do i_e=1,num_e
            kernel_value = hadronic_bh_kernel_electron_generation(ee_dimless(i_e),gp,eph)
            pair_log_source(i_e) = pair_log_source(i_e) + kernel_value*proton_log_value*photon_log_value
            photon_log_loss(i_ph) = photon_log_loss(i_ph) + kernel_value*proton_log_value*photon_log_value
        end do
    end subroutine accumulate_bh_pair_source
end subroutine hadronic_bethe_heitler_operator

! 检查网格为对数均匀并返回对数间距。
real(8) function hadronic_uniform_log_spacing(num_values,values,name)
    integer, intent(in) :: num_values
    real(8), intent(in) :: values(num_values)
    character(*), intent(in) :: name
    real(8) :: dln_local
    call hadronic_validate_log_grid(num_values,values,name,dln_local)
    hadronic_uniform_log_spacing = dln_local
end function hadronic_uniform_log_spacing

! Bethe-Heitler电子产生核：计算给定质子/光子能量下产生能量为ee的电子的微分谱。
real(8) function hadronic_bh_kernel_electron_generation(ee,gp,eph)
    real(8), intent(in) :: ee,gp,eph
    real(8) :: upper,lower,diff

    upper = 2.0d0*gp*eph
    lower = ((gp + ee)*(gp + ee))/(2.0d0*gp*ee)
    diff = (upper - lower)/lower
    if (diff <= 1.0d-8) then
        hadronic_bh_kernel_electron_generation = zero
        return
    end if

    hadronic_bh_kernel_electron_generation = hadronic_rk4_3(hadronic_bh_outer,lower,upper,gp,ee,bh_outer_bins) * &
                                             ee/(2.0d0*gp*gp*gp*eph*eph)
end function hadronic_bh_kernel_electron_generation

! Bethe-Heitler外层积分核（对omega积分）。
real(8) function hadronic_bh_outer(gp,ee,omega)
    real(8), intent(in) :: gp,ee,omega
    real(8), parameter :: reg_outer = 1.0d-5
    real(8) :: lower,upper,diff

    lower = (gp*gp + ee*ee)/(2.0d0*gp*ee) + reg_outer
    upper = omega - 1.0d0 - reg_outer
    diff = (upper - lower)/lower
    if (diff <= -1.0d-4) then
        error stop "Bethe-Heitler outer integral received invalid limits."
    end if
    if (diff < 1.0d-4) then
        hadronic_bh_outer = zero
        return
    end if

    hadronic_bh_outer = hadronic_rk4_4(hadronic_bh_inner,lower,upper,gp,ee,omega,bh_inner_bins)
end function hadronic_bh_outer

! Bethe-Heitler内层积分核（对ebar积分）。
real(8) function hadronic_bh_inner(gp,ee,omega,ebar)
    real(8), intent(in) :: gp,ee,omega,ebar
    real(8), parameter :: reg_inner = 1.0d-20
    real(8) :: pbar,ksi

    pbar = dsqrt(ebar*ebar - 1.0d0 + reg_inner)
    ksi = (gp*ebar - ee)/(gp*pbar)
    if (ksi < -1.0d0 .or. ksi > 1.0d0 .or. ebar < 1.0d0 .or. ebar > omega - 1.0d0) then
        hadronic_bh_inner = zero
        return
    end if

    hadronic_bh_inner = omega*hadronic_bh_sigma_w(omega,ebar,ksi)/pbar
end function hadronic_bh_inner

! Bethe-Heitler微分截面 σ(ω, ε̄, ξ) 的解析表达式（Blumenthal 1970）。
real(8) function hadronic_bh_sigma_w(omega,ebar,ksi)
    real(8), intent(in) :: omega,ebar,ksi
    real(8), parameter :: reg_sigma = 1.0d-30
    real(8) :: k,g,u,p,g1,p1,dloc,tloc,d1t,y1,yloc
    real(8) :: a0,a1,a2,a3,a4,a50,a51,a52,a53,a6,a7

    if (omega < 2.0001d0 .or. omega > 600.0d0) then
        hadronic_bh_sigma_w = zero
        return
    end if

    k = omega
    g = ebar
    u = ksi

    p = dsqrt(g*g - 1.0d0 + reg_sigma)
    g1 = k - g
    p1 = dsqrt(g1*g1 - 1.0d0 + reg_sigma)
    dloc = g - u*p
    tloc = dsqrt(k*k + p*p - 2.0d0*k*p*u)
    d1t = dlog((tloc + p1)/(tloc - p1))
    y1 = dlog((g1 + p1)/(g1 - p1))/p1
    yloc = 2.0d0*dlog((g*g1 + p*p1 + 1.0d0)/k)/(p*p)

    a0 = p*p1/(k**3)
    a1 = -4.0d0*(1.0d0 - u*u)*(2.0d0*g*g + 1.0d0)/(p*p*dloc**4)
    a2 = (5.0d0*g*g - 2.0d0*g*g1 + 3.0d0)/(p*p*dloc*dloc)
    a3 = (p*p - k*k)/(tloc*tloc*dloc*dloc)
    a4 = 2.0d0*g1/(p*p*dloc)
    a50 = yloc/(p*p1)
    a51 = 2.0d0*g*(1.0d0 - u*u)*(3.0d0*k + p*p*g1)/(dloc**4)
    a52 = (2.0d0*g*g*(g*g + g1*g1) - 7.0d0*g*g - 3.0d0*g*g1 - g1*g1 + 1.0d0)/(dloc*dloc)
    a53 = k*(g*g - g*g1 - 1.0d0)/dloc
    a6 = -(d1t/(p1*tloc))*(2.0d0/(dloc*dloc) - 3.0d0*k/dloc - k*(p*p - k*k)/(tloc*tloc*dloc))
    a7 = -2.0d0*y1/dloc
    hadronic_bh_sigma_w = a0*(a1 + a2 + a3 + a4 + a50*(a51 + a52 + a53) + a6 + a7)
end function hadronic_bh_sigma_w

! Bethe-Heitler质子能量损失核函数。
real(8) function hadronic_bh_kernel_proton_loss(gp,eph)
    real(8), intent(in) :: gp,eph
    hadronic_bh_kernel_proton_loss = 0.5d0*hadronic_bh_eloss_kernel_phi(2.0d0*gp*eph)/(gp*gp*eph*eph)
end function hadronic_bh_kernel_proton_loss

! Bethe-Heitler能量损失核 φ(x) 的分段近似（低能多项式/高能对数）。
real(8) function hadronic_bh_eloss_kernel_phi(xval)
    real(8), intent(in) :: xval
    real(8) :: zval,zlog

    if (xval < 2.0d0 + 1.0d-6) then
        hadronic_bh_eloss_kernel_phi = zero
        return
    end if
    if (xval < 25.0d0) then
        zval = xval - 2.0d0
        hadronic_bh_eloss_kernel_phi = (pi/12.0d0)*zval**4 / &
                                       (1.0d0 + 0.8048d0*zval + 0.1459d0*zval*zval + 1.137d-3*zval**3 - 3.879d-6*zval**4)
        return
    end if

    zlog = dlog(xval)
    hadronic_bh_eloss_kernel_phi = xval*(-86.07d0 + 50.96d0*zlog - 14.45d0*zlog*zlog + (8.0d0/3.0d0)*zlog**3) / &
                                   (1.0d0 - 2.910d0/xval - 78.35d0/(xval*xval) - 1837.0d0/(xval**3))
end function hadronic_bh_eloss_kernel_phi

! 三参数版本的三阶Runge-Kutta（RK4-3/8规则）数值积分。
real(8) function hadronic_rk4_3(func,lower,upper,var0,var1,n_bins)
    interface
        function func(arg0,arg1,arg2) result(value)
            real(8), intent(in) :: arg0,arg1,arg2
            real(8) :: value
        end function func
    end interface
    integer, intent(in) :: n_bins
    real(8), intent(in) :: lower,upper,var0,var1
    integer :: i
    real(8) :: dx,xval,k1,k2,k4

    dx = (upper - lower)/dble(n_bins)
    xval = lower
    hadronic_rk4_3 = zero
    do i=1,n_bins
        k1 = dx*func(var0,var1,xval)
        k2 = dx*func(var0,var1,xval + 0.5d0*dx)
        k4 = dx*func(var0,var1,xval + dx)
        hadronic_rk4_3 = hadronic_rk4_3 + (k1 + k4)/6.0d0 + (2.0d0/3.0d0)*k2
        xval = xval + dx
    end do
end function hadronic_rk4_3

! 四参数版本的三阶Runge-Kutta（RK4-3/8规则）数值积分。
real(8) function hadronic_rk4_4(func,lower,upper,var0,var1,var2,n_bins)
    interface
        function func(arg0,arg1,arg2,arg3) result(value)
            real(8), intent(in) :: arg0,arg1,arg2,arg3
            real(8) :: value
        end function func
    end interface
    integer, intent(in) :: n_bins
    real(8), intent(in) :: lower,upper,var0,var1,var2
    integer :: i
    real(8) :: dx,xval,k1,k2,k4

    dx = (upper - lower)/dble(n_bins)
    xval = lower
    hadronic_rk4_4 = zero
    do i=1,n_bins
        k1 = dx*func(var0,var1,var2,xval)
        k2 = dx*func(var0,var1,var2,xval + 0.5d0*dx)
        k4 = dx*func(var0,var1,var2,xval + dx)
        hadronic_rk4_4 = hadronic_rk4_4 + (k1 + k4)/6.0d0 + (2.0d0/3.0d0)*k2
        xval = xval + dx
    end do
end function hadronic_rk4_4

end module hadronic_bethe_heitler_kernel
