!f2py: skip
module hadronic_pair_production_kernel
    use constants
    use hadronic_common, only: hadronic_electron_mass_gev, hadronic_validate_log_grid
    implicit none
    private

    real(8), parameter :: sigma_t_cgs = 6.6524587158d-25
    real(8), parameter :: reg = 1.0d-50
    real(8), parameter :: t_ssc_cgs = dsqrt(sigma_t_cgs)/para_c
    real(8), parameter :: n_ssc_cgs = sigma_t_cgs**(-1.5d0)
    real(8), parameter :: rate_ssc_cgs = n_ssc_cgs/t_ssc_cgs

    public :: hadronic_pair_production_operator

contains

! 光子-光子对产生主算子：计算光子吸收率、正负电子对注入谱和能量守恒校验。
subroutine hadronic_pair_production_operator(num_ph,photon_energy_gev,photon_density_per_gev,num_e,electron_energy_gev, &
                                             max_com_energy_factor,photon_loss_rate,pair_injection_rate_per_gev_per_species, &
                                             pair_injection_rate_per_gev_total,absorbed_power_gev_per_cm3_s, &
                                             injected_power_gev_per_cm3_s)
    integer, intent(in) :: num_ph,num_e,max_com_energy_factor
    real(8), intent(in) :: photon_energy_gev(num_ph),photon_density_per_gev(num_ph),electron_energy_gev(num_e)
    real(8), intent(out) :: photon_loss_rate(num_ph),pair_injection_rate_per_gev_per_species(num_e)
    real(8), intent(out) :: pair_injection_rate_per_gev_total(num_e),absorbed_power_gev_per_cm3_s
    real(8), intent(out) :: injected_power_gev_per_cm3_s
    real(8) :: dln_ph,dln_e,dln
    real(8) :: x_ph(num_ph),gm_e(num_e),photons_log_cgs(num_ph),photons_log_ssc(num_ph),pp_ng(num_ph)
    real(8) :: afpair_ssc(num_ph),epspair_ssc(num_e)
    integer :: ind_min_energy_pho

    dln_ph = hadronic_uniform_log_spacing(num_ph,photon_energy_gev)
    dln_e = hadronic_uniform_log_spacing(num_e,electron_energy_gev)
    dln = dln_ph

    if (dabs(dln_ph-dln_e) > dmax1(1.0d-14,1.0d-10*dabs(dln_ph))) then
        error stop "pair production requires photon/electron grids with the same logarithmic spacing."
    end if

    ind_min_energy_pho = hadronic_grid_index_offset(electron_energy_gev(1),photon_energy_gev(1),dln)
    if (ind_min_energy_pho < 0) then
        error stop "pair production requires electron grid minimum >= photon grid minimum."
    end if

    x_ph = photon_energy_gev/hadronic_electron_mass_gev
    gm_e = electron_energy_gev/hadronic_electron_mass_gev
    photons_log_cgs = photon_energy_gev*photon_density_per_gev
    photons_log_ssc = photons_log_cgs/n_ssc_cgs

    call hadronic_build_photon_loss_kernel(num_ph,x_ph,dln,photons_log_ssc,afpair_ssc)

    pp_ng = photons_log_ssc/(x_ph*x_ph)
    call hadronic_calc_pair_injection_full(num_e,gm_e,num_ph,x_ph,pp_ng,dln,ind_min_energy_pho,max_com_energy_factor, &
                                           afpair_ssc,photons_log_ssc,epspair_ssc)

    photon_loss_rate = afpair_ssc/t_ssc_cgs
    pair_injection_rate_per_gev_per_species = (epspair_ssc*rate_ssc_cgs)/electron_energy_gev
    pair_injection_rate_per_gev_total = two*pair_injection_rate_per_gev_per_species
    absorbed_power_gev_per_cm3_s = sum(x_ph(ind_min_energy_pho+1:num_ph)*afpair_ssc(ind_min_energy_pho+1:num_ph)* &
                                       photons_log_ssc(ind_min_energy_pho+1:num_ph))*rate_ssc_cgs*hadronic_electron_mass_gev
    injected_power_gev_per_cm3_s = sum(two*gm_e*epspair_ssc)*rate_ssc_cgs*hadronic_electron_mass_gev
end subroutine hadronic_pair_production_operator

! 计算完整光子-光子对注入谱：双重卷积+能量守恒重标定。
subroutine hadronic_calc_pair_injection_full(num_e,gm_e,num_ph,x_ph,pp_ng,dln,ind_min_energy_pho,max_com_energy_factor, &
                                             afpair_ssc,photons_log_ssc,epspair_ssc)
    integer, intent(in) :: num_e,num_ph,ind_min_energy_pho,max_com_energy_factor
    real(8), intent(in) :: gm_e(num_e),x_ph(num_ph),pp_ng(num_ph),dln,afpair_ssc(num_ph),photons_log_ssc(num_ph)
    real(8), intent(out) :: epspair_ssc(num_e)
    integer :: i,j,k,outer,inner,kmax
    real(8) :: accum,r_alpha_phot,r_eps_raw

    epspair_ssc = zero
    !$OMP PARALLEL DO if((num_e-2)*(num_ph-ind_min_energy_pho) >= 8192) schedule(dynamic,1) &
    !$OMP& private(i,j,k,outer,inner,kmax,accum)
    do i=2,num_e-1
        outer = hadronic_outer_pp(gm_e(i),dln,ind_min_energy_pho,num_ph)
        accum = zero
        do j=outer+1,num_ph
            inner = min(hadronic_inner_pp(x_ph(j),gm_e(i),dln,ind_min_energy_pho,num_ph),num_ph-1)
            kmax = min(-j + 1 + 2*(ind_min_energy_pho+1) + max_com_energy_factor,num_ph)
            do k=inner+1,kmax
                accum = accum + hadronic_rgg_d1(gm_e(i),x_ph(j),x_ph(k))*pp_ng(k)*pp_ng(j)
            end do
        end do
        epspair_ssc(i) = accum*dln*dln*0.75d0*gm_e(i)
    end do
    !$OMP END PARALLEL DO

    r_alpha_phot = sum( &
        x_ph(ind_min_energy_pho+1:num_ph) * afpair_ssc(ind_min_energy_pho+1:num_ph) * &
        photons_log_ssc(ind_min_energy_pho+1:num_ph) &
    )
    r_eps_raw = sum(gm_e*epspair_ssc)
    if (r_eps_raw > 1.0d-100) then
        epspair_ssc = epspair_ssc*(0.5d0*r_alpha_phot/r_eps_raw)
    end if
    epspair_ssc(1) = zero
end subroutine hadronic_calc_pair_injection_full

! 构建光子损失核：对光子场卷积phibar函数得到各能量光子的吸收率。
subroutine hadronic_build_photon_loss_kernel(num_ph,x_ph,dln,photons_log_ssc,afpair_ssc)
    integer, intent(in) :: num_ph
    real(8), intent(in) :: x_ph(num_ph),dln,photons_log_ssc(num_ph)
    real(8), intent(out) :: afpair_ssc(num_ph)
    integer :: i,j
    real(8) :: accum

    !$OMP PARALLEL DO if(num_ph*num_ph >= 8192) schedule(static) private(i,j,accum)
    do i=1,num_ph
        accum = zero
        do j=1,num_ph
            accum = accum + 0.375d0*hadronic_phibar(x_ph(i),x_ph(j))*photons_log_ssc(j)/(x_ph(i)*x_ph(i)*x_ph(j)*x_ph(j))
        end do
        afpair_ssc(i) = dln*accum
    end do
    !$OMP END PARALLEL DO
end subroutine hadronic_build_photon_loss_kernel

! φ̄(s) 函数：各向同性光子场中对产生截面的角度平均核，s = a*b。
real(8) function hadronic_phibar(a,b)
    real(8), intent(in) :: a,b
    real(8) :: s,w
    s = a*b
    if (s <= one) then
        hadronic_phibar = zero
    else if (s <= 1.1d0) then
        hadronic_phibar = hadronic_phibar1(s)
    else if (s < 5.0d0) then
        w = -one + two*s*(one + dsqrt(one - one/s))
        hadronic_phibar = hadronic_phibar2(s,w)
    else
        hadronic_phibar = hadronic_phibar3(s)
    end if
end function hadronic_phibar

! φ̄ 近阈区近似 (s-1 ≪ 1)：泰勒展开至 (s-1)^(7/2)。
real(8) function hadronic_phibar1(s)
    real(8), intent(in) :: s
    real(8) :: s1,s2
    s1 = s - one
    s2 = s1*dsqrt(s1)
    hadronic_phibar1 = s2*(1.333333d0 + 1.2d0*s1 - 253.0d0*s1*s1/70.0d0)
end function hadronic_phibar1

! φ̄ 中间区解析表达式 (1.1 ≤ s < 5)：使用变量 w 的参数化形式。
real(8) function hadronic_phibar2(s,w)
    real(8), intent(in) :: s,w
    real(8) :: v,u
    v = dlog(w)
    u = one - one/s
    hadronic_phibar2 = ( &
        (two - 4.0d0*s)*dsqrt(u) + v*(4.0d0*dlog(one+w) - 3.0d0*v + s*(one + u*u)) - &
        3.289868d0 + hadronic_phisum(w) &
    )
end function hadronic_phibar2

! φ̄ 高能渐近展开 (s ≥ 5)：主导项 (2s+ln 4s)(ln 4s-2)。
real(8) function hadronic_phibar3(s)
    real(8), intent(in) :: s
    real(8) :: w
    w = dlog(4.0d0*s)
    hadronic_phibar3 = (two*s + w)*(w - two) + (w + 1.125d0)/s - 0.289868d0
end function hadronic_phibar3

! φ̄2中用到的辅助级数求和 Σ(-1)^i/(i^2 w^i)。
real(8) function hadronic_phisum(w)
    real(8), intent(in) :: w
    integer :: i
    hadronic_phisum = zero
    do i=1,14
        hadronic_phisum = hadronic_phisum + dble(hadronic_sign_int(i))/(w**i*dble(i*i))
    end do
    hadronic_phisum = -4.0d0*hadronic_phisum
end function hadronic_phisum

! 对产生积分内边界索引：由光子能量x和电子gamma确定积分下限。
integer function hadronic_inner_pp(x,gm,dln,ind_min_energy_pho,num_ph)
    real(8), intent(in) :: x,gm,dln
    integer, intent(in) :: ind_min_energy_pho,num_ph
    real(8) :: fval,arg

    if (x <= 0.5d0) then
        fval = hadronic_fpp_m(x,gm)
        if (dabs(fval) <= 1.0d-300) then
            hadronic_inner_pp = 0
            return
        end if
        arg = gm - x + 0.5d0*(fval + one/fval)
        if (arg <= zero) then
            hadronic_inner_pp = 0
            return
        end if
        hadronic_inner_pp = min(max(int(dlog(arg)/dln + dble(ind_min_energy_pho) + one),0),num_ph)
        return
    end if
    if (x < one .and. gm < hadronic_gm_b(x)) then
        fval = hadronic_fpp_m(x,gm)
        if (dabs(fval) <= 1.0d-300) then
            hadronic_inner_pp = 0
            return
        end if
        arg = gm - x + 0.5d0*(fval + one/fval)
        if (arg <= zero) then
            hadronic_inner_pp = 0
            return
        end if
        hadronic_inner_pp = min(max(int(dlog(arg)/dln + dble(ind_min_energy_pho) + one),0),num_ph)
        return
    end if
    if (x > one .and. gm < hadronic_gm_b(x)) then
        fval = hadronic_fpp_p(x,gm)
        if (dabs(fval) <= 1.0d-300) then
            hadronic_inner_pp = 0
            return
        end if
        arg = gm - x + 0.5d0*(fval + one/fval)
        if (arg <= zero) then
            hadronic_inner_pp = 0
            return
        end if
        hadronic_inner_pp = min(max(int(dlog(arg)/dln + dble(ind_min_energy_pho) + one),0),num_ph)
        return
    end if
    hadronic_inner_pp = 0
end function hadronic_inner_pp

! 对产生积分外边界索引：由电子gamma确定积分上限。
integer function hadronic_outer_pp(gm,dln,ind_min_energy_pho,num_ph)
    real(8), intent(in) :: gm,dln
    integer, intent(in) :: ind_min_energy_pho,num_ph
    hadronic_outer_pp = min(max(int(dlog(hadronic_x_l(gm))/dln + dble(ind_min_energy_pho)),0),num_ph-1)
end function hadronic_outer_pp

! 最小光子能量 x_- = (γ - √(γ²-1))/2。
real(8) function hadronic_x_l(gm)
    real(8), intent(in) :: gm
    hadronic_x_l = 0.5d0*gm*(one - hadronic_beta(gm))
end function hadronic_x_l

! 分支边界函数 γ_b(x) = x - (x-1)/(2x-1)，区分积分域。
real(8) function hadronic_gm_b(x)
    real(8), intent(in) :: x
    hadronic_gm_b = x - (x - one)/(two*x - one)
end function hadronic_gm_b

! 辅助函数 f_- = 2x - γ(1-β)，用于确定x_1积分边界。
real(8) function hadronic_fpp_m(x,gp)
    real(8), intent(in) :: x,gp
    hadronic_fpp_m = two*x - gp*(one - hadronic_beta(gp))
end function hadronic_fpp_m

! 辅助函数 f_+ = 2x - γ(1+β)，用于确定x_1积分边界。
real(8) function hadronic_fpp_p(x,gp)
    real(8), intent(in) :: x,gp
    hadronic_fpp_p = two*x - gp*(one + hadronic_beta(gp))
end function hadronic_fpp_p

! 对产生微分截面核 R_γγ(γ, x, x₁)：Aharonian+1983公式，含运动学截断。
real(8) function hadronic_rgg_d1(gm,x,x1)
    real(8), intent(in) :: gm,x,x1
    real(8) :: gp,kl,ku,tval,x_upper,x_lower,u,v,h_nh_u,h_nh_l,h_ih_u,h_ih_l,td2a,td2b,td2c,td2d,denom

    gp = x + x1 - gm
    if (gp <= one .or. gm <= one) then
        hadronic_rgg_d1 = zero
        return
    end if

    if (gp > ten .and. gm > ten) then
        kl = 0.25d0/gp + 0.09375d0/gp**3 + 0.25d0/gm + 0.09375d0/gm**3
        ku = 0.5d0*(gp + dsqrt(gp*gp - one) + gm + dsqrt(gm*gm - one))
    else if (gp > ten) then
        kl = 0.5d0*(0.5d0/gp + gm - dsqrt(gm*gm - one))
        ku = gp + 0.5d0*(gm + dsqrt(gm*gm - one))
    else if (gm > ten) then
        kl = 0.5d0*(gp - dsqrt(gp*gp - one) + 0.5d0/gm)
        ku = gm + 0.5d0*(gp + dsqrt(gp*gp - one))
    else
        kl = 0.5d0*(gp + gm - dsqrt(gp*gp - one) - dsqrt(gm*gm - one))
        ku = 0.5d0*(gp + gm + dsqrt(gp*gp - one) + dsqrt(gm*gm - one))
    end if

    if (x < kl .or. x1 < kl .or. x > ku .or. x1 > ku) then
        hadronic_rgg_d1 = zero
        return
    end if

    tval = x + x1
    x_upper = hadronic_xcm_u(gm,x,x1)
    x_lower = hadronic_xcm_l(gm,x,x1)
    u = 4.0d0*x_upper*x_upper/(tval*tval)
    v = 4.0d0*x_lower*x_lower/(tval*tval)

    h_nh_u = ((gm - x)*(gm - x) - one)*x_upper*x_upper/(x*x1)
    h_nh_l = ((gm - x)*(gm - x) - one)*x_lower*x_lower/(x*x1)
    h_ih_u = ((gm - x1)*(gm - x1) - one)*x_upper*x_upper/(x*x1)
    h_ih_l = ((gm - x1)*(gm - x1) - one)*x_lower*x_lower/(x*x1)

    td2a = hadronic_td2(gm,x,x1,x_upper)
    td2b = hadronic_td2(gm,x,x1,x_lower)
    td2c = hadronic_td2(gm,x1,x,x_upper)
    td2d = hadronic_td2(gm,x1,x,x_lower)

    if (h_nh_u > 1000.0d0 .and. h_nh_l > 1000.0d0) then
        denom = dsqrt((gm - x)*(gm - x) - one)
        td2a = -(x_upper*x_upper + two*dlog(x_upper))/denom
        td2b = -(x_lower*x_lower + two*dlog(x_lower))/denom
    end if
    if (h_ih_u > 1000.0d0 .and. h_ih_l > 1000.0d0) then
        denom = dsqrt((gm - x1)*(gm - x1) - one)
        td2c = -(x_upper*x_upper + two*dlog(x_upper))/denom
        td2d = -(x_lower*x_lower + two*dlog(x_lower))/denom
    end if

    if (u < 0.01d0 .and. v < 0.01d0) then
        hadronic_rgg_d1 = -0.25d0*(0.5d0*tval*(u - v) + td2a - td2b + td2c - td2d)
    else
        hadronic_rgg_d1 = -0.25d0*(tval*(dsqrt(one - v) - dsqrt(one - u + 1.0d-9)) + td2a - td2b + td2c - td2d)
    end if
end function hadronic_rgg_d1

! 质心系能量上界 x_cm,upper = min(√(xx₁), √(γγ_p))。
real(8) function hadronic_xcm_u(gm,x,x1)
    real(8), intent(in) :: gm,x,x1
    real(8) :: gp
    gp = x + x1 - gm
    if (gp > 15.0d0 .and. gm > 15.0d0) then
        hadronic_xcm_u = min(dsqrt(x*x1),dsqrt(gp*gm))
    else
        hadronic_xcm_u = min(dsqrt(x*x1),dsqrt(0.5d0*(gm*gp + one + dsqrt(gm*gm - one)*dsqrt(gp*gp - one))))
    end if
end function hadronic_xcm_u

! 质心系能量下界：含相对论/非相对论/极端相对论三种渐近分支。
real(8) function hadronic_xcm_l(gm,x,x1)
    real(8), intent(in) :: gm,x,x1
    real(8) :: gp
    gp = x + x1 - gm
    if (gp > ten) then
        if (gm > ten) then
            hadronic_xcm_l = 0.5d0*dsqrt(two + gp/gm + gm/gp)
            return
        end if
        if (gm < 1.001d0) then
            hadronic_xcm_l = 1.001d0
            return
        end if
        hadronic_xcm_l = dsqrt(0.5d0*((gm - dsqrt(gm*gm - one + reg))*(gp - 0.5d0/gp) + one))
        return
    end if
    if (gp < 1.001d0) then
        hadronic_xcm_l = 1.001d0
        return
    end if
    if (gm > ten) then
        hadronic_xcm_l = dsqrt(0.5d0*((gp - dsqrt(gp*gp - one + reg))*(gm - 0.5d0/gm) + one))
        return
    end if
    if (gm < 1.001d0) then
        hadronic_xcm_l = 1.001d0
        return
    end if
    hadronic_xcm_l = dsqrt(0.5d0*(gp*gm + one - dsqrt((gp*gp - one)*(gm*gm - one))))
end function hadronic_xcm_l

! 微分截面中的 T(γ,x,x₁,x_cm) 辅助函数。
real(8) function hadronic_td2(gm,x,x1,xcm)
    real(8), intent(in) :: gm,x,x1,xcm
    real(8) :: y,y2,q,h,ah
    y = dsqrt(x*x1)
    y2 = y*y
    q = xcm/y
    h = ((gm - x)*(gm - x) - one)*q*q
    if (h > 1000.0d0) then
        hadronic_td2 = -(xcm*xcm + h/(xcm*xcm) + 0.886294d0 - 0.5d0*(gm*(x1 - x) + x*x)/y2 + dlog(h))*q/dsqrt(h)
        return
    end if
    if (h < -0.999d0) then
        hadronic_td2 = q*(q*q*(y2 - one)*(-1.5707963d0) + 0.5d0*((x1 - x)/y2 - 6.28319d0))
        return
    end if
    ah = dsqrt(one + h)
    hadronic_td2 = q*q*q*(y2 - one)*hadronic_a_a0_h(h) - ah/(xcm*y) + &
                   0.5d0*q*((one + (gm*(x1 - x) + x*x)/y2 - two*q*q)/ah - 4.0d0*hadronic_a0(h))
end function hadronic_td2

! (a₀(h) - √(1+h))/h 的稳定计算。
real(8) function hadronic_a_a0_h(h)
    real(8), intent(in) :: h
    if (h <= 20.0d0) then
        hadronic_a_a0_h = hadronic_a_a0_hf(h)
    else
        hadronic_a_a0_h = -(one - (0.693147d0 + 0.5d0*dlog(h))/h)/dsqrt(h)
    end if
end function hadronic_a_a0_h

! a₀(h) = arcsinh(√h)/√h (h>0) 或 arcsin(√|h|)/√|h| (h<0)。
real(8) function hadronic_a0(h)
    real(8), intent(in) :: h
    if (h < 20.0d0) then
        hadronic_a0 = hadronic_a0f(h)
    else
        hadronic_a0 = (0.693147d0 + 0.5d0*dlog(h))/dsqrt(h)
    end if
end function hadronic_a0

! a₀(h) 基础实现：h>0用asinh, h<0用asin, |h|→0用泰勒展开。
real(8) function hadronic_a0f(h)
    real(8), intent(in) :: h
    if (h > 1.0d-6) then
        hadronic_a0f = dlog(dsqrt(h) + dsqrt(one + h))/dsqrt(h)
    else if (h < -1.0d-6) then
        hadronic_a0f = dasin(dsqrt(-h) - 1.0d-8)/dsqrt(-h)
    else
        hadronic_a0f = one - h/6.0d0
    end if
end function hadronic_a0f

! (a₀(h) - √(1+h))/h 基础实现：|h|≤0.1用泰勒展开提高精度。
real(8) function hadronic_a_a0_hf(h)
    real(8), intent(in) :: h
    if (h >= -0.1d0 .and. h <= 0.1d0) then
        hadronic_a_a0_hf = -0.666667d0 + (0.2d0 - (0.107143d0 - 0.0694444d0*h)*h)*h
    else
        hadronic_a_a0_hf = (hadronic_a0f(h) - dsqrt(one + h + 1.0d-8))/h
    end if
end function hadronic_a_a0_hf

! 相对论β因子：γ<10用精确公式，γ≥10用渐近展开。
real(8) function hadronic_beta(gamma)
    real(8), intent(in) :: gamma
    if (gamma < ten) then
        hadronic_beta = dsqrt(one - one/(gamma*gamma) + reg)
    else
        hadronic_beta = one - 0.5d0/(gamma*gamma)*(one + 0.25d0/(gamma*gamma))
    end if
end function hadronic_beta

! 符号函数：偶数为+1，奇数为-1。
integer function hadronic_sign_int(i)
    integer, intent(in) :: i
    if (mod(i,2) == 0) then
        hadronic_sign_int = 1
    else
        hadronic_sign_int = -1
    end if
end function hadronic_sign_int

! 验证并返回网格的对数均匀间距。
real(8) function hadronic_uniform_log_spacing(num_grid,grid)
    integer, intent(in) :: num_grid
    real(8), intent(in) :: grid(num_grid)
    real(8) :: dln_local
    call hadronic_validate_log_grid(num_grid,grid,"pair_production_grid",dln_local)
    hadronic_uniform_log_spacing = dln_local
end function hadronic_uniform_log_spacing

! 计算两网格最小能量值在对数空间中的整数索引偏移量。
integer function hadronic_grid_index_offset(e_min,ph_min,dln)
    real(8), intent(in) :: e_min,ph_min,dln
    real(8) :: ratio
    integer :: rounded
    ratio = dlog(e_min/ph_min)/dln
    rounded = nint(ratio)
    if (dabs(ratio - dble(rounded)) > 1.0d-9) then
        error stop "pair production grid minima are not aligned on the shared logarithmic lattice."
    end if
    hadronic_grid_index_offset = rounded
end function hadronic_grid_index_offset

end module hadronic_pair_production_kernel
