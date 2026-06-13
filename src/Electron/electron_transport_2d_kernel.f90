!f2py: skip
module electron_transport_2d_kernel
  use constants
  use electron_transport_common, only: electron_prepare_implicit_coeffs_common, electron_backward_sweep_common, &
                             electron_prepare_conservative_remap_nonuniform, electron_ppm_prefix_eval_nonuniform, &
                             electron_characteristic_update, electron_cooling_affine, electron_cooling_piecewise
  use electron_injection_profiles, only: electron_profile_log_cell_edges
  implicit real(8)(a-h,o-z)
  private

  public :: compute_log_chi_geometry, get_shock_transport_state
  public :: compute_downstream_comoving_grid
  public :: bm_beta2_lab, bm_beta2_shock
  public :: compute_bm_divergence_chi
  public :: compute_logchi_eta_step_limit
  public :: advance_eta_logchi_implicit, advance_energy_loggamma_chi
  public :: advance_eta_logchi_pwncr_implicit, advance_energy_loggamma_chi_pwncr
  public :: advance_energy_stochastic_loggamma_chi
  public :: advance_eta_logchi_advection_charint, advance_eta_logchi_diffusion_implicit
  public :: advance_energy_loggamma_chi_charint

contains

! 构建对数χ几何网格：计算加速度项a=8Γ²/R及其导数，生成η=log₁₀(χ)均匀网格。
subroutine compute_log_chi_geometry(Num_R, Num_chi, R, R_Gamma, a_arr, dln_a_dR_arr, &
                                    chi_max_global, deta, eta_face, eta_grid, chi_face, chi_grid)
    integer, intent(in) :: Num_R, Num_chi
    real(8), intent(in) :: R(Num_R), R_Gamma(Num_R)
    real(8), intent(out) :: a_arr(Num_R), dln_a_dR_arr(Num_R)
    real(8), intent(out) :: chi_max_global, deta
    real(8), intent(out) :: eta_face(0:Num_chi), eta_grid(Num_chi)
    real(8), intent(out) :: chi_face(0:Num_chi), chi_grid(Num_chi)

    integer :: I_R, I_chi

    do I_R = 1, Num_R
        a_arr(I_R) = 8d0*R_Gamma(I_R)*R_Gamma(I_R)/R(I_R)
    end do

    if (Num_R > 1) then
        dln_a_dR_arr(1) = (dlog(a_arr(2))-dlog(a_arr(1))) / (R(2)-R(1))
        dln_a_dR_arr(Num_R) = (dlog(a_arr(Num_R))-dlog(a_arr(Num_R-1))) / (R(Num_R)-R(Num_R-1))
        do I_R = 2, Num_R-1
            dln_a_dR_arr(I_R) = (dlog(a_arr(I_R+1))-dlog(a_arr(I_R-1))) / (R(I_R+1)-R(I_R-1))
        end do
    else
        dln_a_dR_arr(1) = -one/R(1)
    end if

    chi_max_global = one + 8d0*maxval(R_Gamma*R_Gamma)
    deta = dlog10(chi_max_global) / dble(Num_chi)

    do I_chi = 0, Num_chi
        eta_face(I_chi) = dble(I_chi)*deta
        chi_face(I_chi) = ten**eta_face(I_chi)
    end do

    do I_chi = 1, Num_chi
        eta_grid(I_chi) = (dble(I_chi)-0.5d0)*deta
        chi_grid(I_chi) = ten**eta_grid(I_chi)
    end do
end subroutine compute_log_chi_geometry

! 将激波系χ网格转换为下游共动系空间网格，含相对论长度收缩修正。
subroutine compute_downstream_comoving_grid(Num_chi,R_shell,Gamma_sh,chi_face,chi_grid,x_face,x_comov_face,x_comov_center,dx_comov)
    integer, intent(in) :: Num_chi
    integer :: I_chi
    real(8), intent(in) :: R_shell,Gamma_sh,chi_face(0:Num_chi),chi_grid(Num_chi)
    real(8), intent(out) :: x_face(0:Num_chi),x_comov_face(0:Num_chi),x_comov_center(Num_chi),dx_comov(Num_chi)
    real(8) :: dx_shock,beta_rel,gamma_rel,chi_mid

    x_face = (chi_face-one)*R_shell/(8d0*Gamma_sh*Gamma_sh)
    x_comov_face(0) = zero
    do I_chi = 1, Num_chi
        dx_shock = x_face(I_chi)-x_face(I_chi-1)
        chi_mid = chi_grid(I_chi)
        beta_rel = bm_beta2_shock(Gamma_sh, chi_mid)
        gamma_rel = one/dsqrt(max(one-beta_rel*beta_rel, tiny(one)))
        dx_comov(I_chi) = gamma_rel*dx_shock
        x_comov_face(I_chi) = x_comov_face(I_chi-1) + dx_comov(I_chi)
        x_comov_center(I_chi) = 0.5d0*(x_comov_face(I_chi-1)+x_comov_face(I_chi))
    end do
end subroutine compute_downstream_comoving_grid

! 计算激波输运状态量：激波速度β_sh，下游速度β_2及其相对激波的速度β_2_sh。
subroutine get_shock_transport_state(Gamma_sh, beta_sh, beta_2, beta_2_sh)
    real(8), intent(in) :: Gamma_sh
    real(8), intent(out) :: beta_sh, beta_2, beta_2_sh
    real(8) :: Gamma_2

    beta_sh = dsqrt(one-one/Gamma_sh**2)
    Gamma_2 = Gamma_sh/dsqrt(two)
    if (Gamma_2 > one) then
        beta_2 = dsqrt(one-one/Gamma_2**2)
    else
        beta_2 = zero
    end if
    beta_2_sh = (beta_sh-beta_2)/(one-beta_sh*beta_2)
end subroutine get_shock_transport_state

! BM解中实验室系下游速度：Γ_2 = Γ_sh/√(2χ)，β_2 = √(1-1/Γ_2²)。
real(8) function bm_beta2_lab(Gamma_sh, chi_loc)
    real(8), intent(in) :: Gamma_sh, chi_loc
    real(8) :: Gamma_2

    Gamma_2 = Gamma_sh/dsqrt(two*chi_loc)
    if (Gamma_2 > one) then
        bm_beta2_lab = dsqrt(one-one/Gamma_2**2)
    else
        bm_beta2_lab = zero
    end if
end function bm_beta2_lab

! BM解中下游相对激波的速度：β_2_sh = (β_sh-β_2)/(1-β_sh·β_2)。
real(8) function bm_beta2_shock(Gamma_sh, chi_loc)
    real(8), intent(in) :: Gamma_sh, chi_loc
    real(8) :: beta_sh, beta_2

    beta_sh = dsqrt(one-one/Gamma_sh**2)
    beta_2 = bm_beta2_lab(Gamma_sh, chi_loc)
    bm_beta2_shock = (beta_sh-beta_2)/(one-beta_sh*beta_2)
end function bm_beta2_shock

! BM下游局部散度：用实验室系径向速度场β(χ)计算(∇·v)/(3β_sh c)，返回log10 γ方程系数。
subroutine compute_bm_divergence_chi(Num_chi, R_loc, Gamma_sh, beta_sh, chi_grid, adiabatic_log_coeff)
    integer, intent(in) :: Num_chi
    real(8), intent(in) :: R_loc, Gamma_sh, beta_sh, chi_grid(Num_chi)
    real(8), intent(out) :: adiabatic_log_coeff(Num_chi)
    real(8) :: Gamma_2, beta_2, radius_ratio, div_over_c, uniform_coeff
    integer :: I_chi

    uniform_coeff = one/(R_loc*dlog(ten))
    do I_chi = 1, Num_chi
        Gamma_2 = Gamma_sh/dsqrt(two*chi_grid(I_chi))
        radius_ratio = one - (chi_grid(I_chi)-one)/(8d0*Gamma_sh*Gamma_sh)
        if (Gamma_2 >= two .and. radius_ratio > zero) then
            beta_2 = bm_beta2_lab(Gamma_sh, chi_grid(I_chi))
            div_over_c = (two*beta_2/radius_ratio + 8d0/beta_2)/R_loc
            adiabatic_log_coeff(I_chi) = div_over_c/(3d0*beta_sh*dlog(ten))
        else
            adiabatic_log_coeff(I_chi) = uniform_coeff
        end if
    end do
end subroutine compute_bm_divergence_chi

! 计算η网格面上的输运系数A(η) = [a·β₂_sh/(χ·β_sh) + (χ-1)/χ·dln a/dR] / ln(10)。
real(8) function eta_face_transport_coeff(Gamma_sh, chi_loc, a_loc, dln_a_dR_loc, beta_sh)
    real(8), intent(in) :: Gamma_sh, chi_loc, a_loc, dln_a_dR_loc, beta_sh
    real(8) :: beta_2_sh

    beta_2_sh = bm_beta2_shock(Gamma_sh, chi_loc)
    eta_face_transport_coeff = (a_loc*beta_2_sh/(chi_loc*beta_sh) + &
                               ((chi_loc-one)/chi_loc)*dln_a_dR_loc) / dlog(ten)
end function eta_face_transport_coeff

subroutine eta_face_transport_coeffs(Gamma_sh, Num_chi, chi_face, a_loc, dln_a_dR_loc, beta_sh, A_eta_face)
    integer, intent(in) :: Num_chi
    real(8), intent(in) :: Gamma_sh, chi_face(0:Num_chi), a_loc, dln_a_dR_loc, beta_sh
    real(8), intent(out) :: A_eta_face(1:Num_chi)
    integer :: I_face

    do I_face = 1, Num_chi
        A_eta_face(I_face) = eta_face_transport_coeff(Gamma_sh, chi_face(I_face), a_loc, dln_a_dR_loc, beta_sh)
    end do
end subroutine eta_face_transport_coeffs

subroutine eta_face_transport_coeffs_all(Gamma_sh, Num_chi, chi_face, a_loc, dln_a_dR_loc, beta_sh, A_eta_face)
    integer, intent(in) :: Num_chi
    real(8), intent(in) :: Gamma_sh, chi_face(0:Num_chi), a_loc, dln_a_dR_loc, beta_sh
    real(8), intent(out) :: A_eta_face(0:Num_chi)
    integer :: I_face

    do I_face = 0, Num_chi
        A_eta_face(I_face) = eta_face_transport_coeff(Gamma_sh, chi_face(I_face), a_loc, dln_a_dR_loc, beta_sh)
    end do
end subroutine eta_face_transport_coeffs_all

! 估计χ方向对流子步长限制。
real(8) function compute_logchi_eta_step_limit(Num_chi,R_loc,R_Gamma_loc,beta_sh, &
                                               dln_a_dR_loc,deta,chi_face,cfl_factor)
    implicit none
    integer, intent(in) :: Num_chi
    integer :: I_chi
    real(8), intent(in) :: R_loc,R_Gamma_loc,beta_sh,dln_a_dR_loc,deta,cfl_factor
    real(8), intent(in) :: chi_face(0:Num_chi)
    real(8) :: A_eta_face,max_eta_coeff

    max_eta_coeff=zero
    do I_chi=1,Num_chi
        A_eta_face=eta_face_transport_coeff(R_Gamma_loc,chi_face(I_chi),8d0*R_Gamma_loc*R_Gamma_loc/R_loc, &
                                            dln_a_dR_loc,beta_sh)
        max_eta_coeff=max(max_eta_coeff,dabs(A_eta_face))
    end do
    compute_logchi_eta_step_limit=huge(one)
    if (max_eta_coeff > zero) compute_logchi_eta_step_limit=cfl_factor*deta/max_eta_coeff
end function compute_logchi_eta_step_limit

real(8) function eta_linear_back_position(eta_start, a_cell, b_cell, lag)
    real(8), intent(in) :: eta_start, a_cell, b_cell, lag
    real(8) :: shift

    if (dabs(b_cell) <= 1d-30) then
        eta_linear_back_position = eta_start - a_cell*lag
    else
        shift = a_cell/b_cell
        eta_linear_back_position = (eta_start+shift)*dexp(-b_cell*lag) - shift
    end if
end function eta_linear_back_position

real(8) function eta_linear_hit_time(eta_start, eta_bound, a_cell, b_cell)
    real(8), intent(in) :: eta_start, eta_bound, a_cell, b_cell
    real(8) :: shift, num, den

    if (dabs(b_cell) <= 1d-30) then
        eta_linear_hit_time = (eta_bound-eta_start)/(-a_cell)
    else
        shift = a_cell/b_cell
        num = eta_start+shift
        den = eta_bound+shift
        eta_linear_hit_time = dlog(num/den)/b_cell
    end if
end function eta_linear_hit_time

subroutine eta_trace_back_faces_piecewise(Num_chi, dR_step, eta_face, A_eta_face, eta_back)
    integer, intent(in) :: Num_chi
    real(8), intent(in) :: dR_step, eta_face(0:Num_chi), A_eta_face(0:Num_chi)
    real(8), intent(out) :: eta_back(0:Num_chi)
    integer :: I_face, I_cell
    real(8) :: eta_cur, eta_trial, eta_bound, s_rem, s_hit, a_cell, b_cell, A_cur, d_eta

    eta_back(0) = eta_face(0)
    do I_face = 1, Num_chi
        eta_cur = eta_face(I_face)
        A_cur = A_eta_face(I_face)
        if (A_cur == zero) then
            eta_back(I_face) = eta_cur
            cycle
        end if
        if (A_cur < zero .and. I_face == Num_chi) then
            eta_back(I_face) = eta_face(Num_chi)
            cycle
        end if
        if (A_cur > zero) then
            I_cell = I_face
        else
            I_cell = I_face + 1
        end if
        s_rem = dR_step
        do while (s_rem > zero)
            d_eta = eta_face(I_cell)-eta_face(I_cell-1)
            b_cell = (A_eta_face(I_cell)-A_eta_face(I_cell-1))/d_eta
            a_cell = A_eta_face(I_cell-1)-b_cell*eta_face(I_cell-1)
            A_cur = a_cell+b_cell*eta_cur
            if (A_cur == zero) exit
            eta_trial = eta_linear_back_position(eta_cur, a_cell, b_cell, s_rem)
            if (A_cur > zero) then
                eta_bound = eta_face(I_cell-1)
                if (eta_trial > eta_bound) then
                    eta_cur = eta_trial
                    exit
                end if
                s_hit = eta_linear_hit_time(eta_cur, eta_bound, a_cell, b_cell)
                eta_cur = eta_bound
                s_rem = s_rem-s_hit
                I_cell = I_cell-1
                if (I_cell < 1) exit
            else
                eta_bound = eta_face(I_cell)
                if (eta_trial < eta_bound) then
                    eta_cur = eta_trial
                    exit
                end if
                s_hit = eta_linear_hit_time(eta_cur, eta_bound, a_cell, b_cell)
                eta_cur = eta_bound
                s_rem = s_rem-s_hit
                I_cell = I_cell+1
                if (I_cell > Num_chi) exit
            end if
        end do
        eta_back(I_face) = eta_cur
    end do
end subroutine eta_trace_back_faces_piecewise

subroutine eta_split_advection_faces(Num_chi,A_eta_face,adv_face_left,adv_face_right)
    implicit none
    integer, intent(in) :: Num_chi
    integer :: I_face
    real(8), intent(in) :: A_eta_face(1:Num_chi)
    real(8), intent(out) :: adv_face_left(1:Num_chi),adv_face_right(1:Num_chi)

    adv_face_left=zero
    adv_face_right=zero
    do I_face=1,Num_chi-1
        if (A_eta_face(I_face) >= zero) then
            adv_face_left(I_face)=A_eta_face(I_face)
        else
            adv_face_right(I_face)=A_eta_face(I_face)
        end if
    end do
    if (A_eta_face(Num_chi) > zero) adv_face_left(Num_chi)=A_eta_face(Num_chi)
end subroutine eta_split_advection_faces

subroutine eta_diffusion_face_coeffs(Num_chi,deta,chi_face,a_loc,beta_sh,include_outer_face, &
                                     diff_face_left_base,diff_face_right_base)
    implicit none
    integer, intent(in) :: Num_chi
    integer :: I_face
    logical, intent(in) :: include_outer_face
    real(8), intent(in) :: deta,chi_face(0:Num_chi),a_loc,beta_sh
    real(8), intent(out) :: diff_face_left_base(1:Num_chi),diff_face_right_base(1:Num_chi)
    real(8) :: diff_prefactor,ln10

    ln10=dlog(ten)
    diff_face_left_base=zero
    diff_face_right_base=zero
    do I_face=1,Num_chi-1
        diff_prefactor=a_loc*a_loc/(chi_face(I_face)*chi_face(I_face)*beta_sh*para_c*ln10*ln10)
        diff_face_left_base(I_face)=diff_prefactor/deta+0.5d0*ln10*diff_prefactor
        diff_face_right_base(I_face)=-diff_prefactor/deta+0.5d0*ln10*diff_prefactor
    end do
    if (include_outer_face) then
        diff_prefactor=a_loc*a_loc/(chi_face(Num_chi)*chi_face(Num_chi)*beta_sh*para_c*ln10*ln10)
        diff_face_left_base(Num_chi)=diff_prefactor/deta+0.5d0*ln10*diff_prefactor
    end if
end subroutine eta_diffusion_face_coeffs

! Thomas算法求解三对角线性方程组。
subroutine solve_tridiagonal(Num_cell, lower, diag, upper, rhs, sol)
    integer, intent(in) :: Num_cell
    real(8), intent(in) :: lower(Num_cell), diag(Num_cell), upper(Num_cell), rhs(Num_cell)
    real(8), intent(out) :: sol(Num_cell)
    real(8) :: cprime(Num_cell), dprime(Num_cell), denom
    integer :: I_cell

    denom = diag(1)
    cprime(1) = upper(1)/denom
    dprime(1) = rhs(1)/denom

    do I_cell = 2, Num_cell
        denom = diag(I_cell) - lower(I_cell)*cprime(I_cell-1)
        if (I_cell < Num_cell) cprime(I_cell) = upper(I_cell)/denom
        dprime(I_cell) = (rhs(I_cell)-lower(I_cell)*dprime(I_cell-1))/denom
    end do

    sol(Num_cell) = dprime(Num_cell)
    do I_cell = Num_cell-1, 1, -1
        sol(I_cell) = dprime(I_cell) - cprime(I_cell)*sol(I_cell+1)
    end do
end subroutine solve_tridiagonal

! η方向纯对流推进（特征线积分）：PPM重构+特征线回溯。
subroutine advance_eta_logchi_advection_charint(U_log, Num_gam_e, Num_chi, active_hi, deta, eta_face, chi_face, &
                                                Gamma_sh, a_loc, dln_a_dR_loc, beta_sh, source_eta1, dR_step)
    integer, intent(in) :: Num_gam_e, Num_chi, active_hi
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: deta, eta_face(0:Num_chi), chi_face(0:Num_chi)
    real(8), intent(in) :: Gamma_sh, a_loc, dln_a_dR_loc, beta_sh, source_eta1(Num_gam_e), dR_step

    real(8) :: A_eta_face(0:Num_chi), eta_back(0:Num_chi)
    real(8) :: U_eta_in(Num_chi), U_eta_out(Num_chi)
    real(8) :: q_left(Num_chi), q_right(Num_chi), prefix(0:Num_chi)
    integer :: I_gam_e

    call eta_face_transport_coeffs_all(Gamma_sh, Num_chi, chi_face, a_loc, dln_a_dR_loc, beta_sh, A_eta_face)
    call eta_trace_back_faces_piecewise(Num_chi, dR_step, eta_face, A_eta_face, eta_back)

    do I_gam_e = 1, active_hi
        U_eta_in = U_log(I_gam_e, :)
        U_eta_in(1) = U_eta_in(1) + dR_step*source_eta1(I_gam_e)
        call electron_prepare_conservative_remap_nonuniform(Num_chi, eta_face, U_eta_in, q_left, q_right, prefix)
        do I_face = 1, Num_chi
            U_eta_out(I_face) = &
                (electron_ppm_prefix_eval_nonuniform(Num_chi, eta_face, U_eta_in, &
                                                     q_left, q_right, prefix, eta_back(I_face)) - &
                 electron_ppm_prefix_eval_nonuniform(Num_chi, eta_face, U_eta_in, &
                                                     q_left, q_right, prefix, eta_back(I_face-1))) / deta
        end do
        U_log(I_gam_e, :) = max(zero, U_eta_out)
    end do
end subroutine advance_eta_logchi_advection_charint

! η方向纯扩散推进（隐式）：三对角求解扩散项 ∂_η(κ∂_η U)。
subroutine advance_eta_logchi_diffusion_implicit(U_log, Num_gam_e, Num_chi, active_hi, deta, chi_face, Gamma_sh, a_loc, &
                                                 dln_a_dR_loc, beta_sh, kappa2_chi, dR_step, n_threads)
    integer, intent(in) :: Num_gam_e, Num_chi, active_hi, n_threads
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: deta, chi_face(0:Num_chi), Gamma_sh, a_loc, dln_a_dR_loc
    real(8), intent(in) :: beta_sh, kappa2_chi(Num_gam_e,Num_chi), dR_step

    real(8) :: diff_face_left_base(1:Num_chi), diff_face_right_base(1:Num_chi)
    real(8) :: lower(Num_chi), diag(Num_chi), upper(Num_chi), rhs(Num_chi), sol(Num_chi)
    real(8) :: lambda_eta, kappa_face
    integer :: I_gam_e, I_chi

    lambda_eta = dR_step/deta
    call eta_diffusion_face_coeffs(Num_chi,deta,chi_face,a_loc,beta_sh,.false., &
                                   diff_face_left_base,diff_face_right_base)

    do I_gam_e = 1, active_hi
        lower = zero
        diag = one
        upper = zero
        do I_chi = 1, Num_chi
            if (I_chi < Num_chi) then
                kappa_face = 0.5d0*(kappa2_chi(I_gam_e,I_chi)+kappa2_chi(I_gam_e,I_chi+1))
                diag(I_chi) = diag(I_chi) + lambda_eta*kappa_face*diff_face_left_base(I_chi)
                upper(I_chi) = upper(I_chi) + lambda_eta*kappa_face*diff_face_right_base(I_chi)
            end if
            if (I_chi > 1) then
                kappa_face = 0.5d0*(kappa2_chi(I_gam_e,I_chi-1)+kappa2_chi(I_gam_e,I_chi))
                lower(I_chi) = lower(I_chi) - lambda_eta*kappa_face*diff_face_left_base(I_chi-1)
                diag(I_chi) = diag(I_chi) - lambda_eta*kappa_face*diff_face_right_base(I_chi-1)
            end if
        end do
        rhs = U_log(I_gam_e, :)
        call solve_tridiagonal(Num_chi, lower, diag, upper, rhs, sol)
        U_log(I_gam_e, :) = max(zero, sol)
    end do
end subroutine advance_eta_logchi_diffusion_implicit

! η方向对流+扩散联合隐式推进：迎风对流+中心扩散，三对角求解。
subroutine advance_eta_logchi_implicit(U_log, Num_gam_e, Num_chi, active_hi, deta, chi_face, Gamma_sh, a_loc, dln_a_dR_loc, &
                                       beta_sh, kappa2_chi, source_eta1, dR_step, n_threads)
    integer, intent(in) :: Num_gam_e, Num_chi, active_hi, n_threads
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: deta, chi_face(0:Num_chi), Gamma_sh, a_loc, dln_a_dR_loc
      real(8), intent(in) :: beta_sh, kappa2_chi(Num_gam_e,Num_chi), source_eta1(Num_gam_e)
    real(8), intent(in) :: dR_step

    real(8) :: A_eta_face(1:Num_chi), adv_face_left(1:Num_chi), adv_face_right(1:Num_chi)
    real(8) :: diff_face_left_base(1:Num_chi), diff_face_right_base(1:Num_chi)
      real(8) :: lower_base(Num_chi), diag_base(Num_chi), upper_base(Num_chi)
      real(8) :: lower(Num_chi), diag(Num_chi), upper(Num_chi), rhs(Num_chi), sol(Num_chi)
      real(8) :: lambda_eta, kappa_face
    integer :: I_gam_e, I_chi

    lambda_eta = dR_step/deta
    call eta_face_transport_coeffs(Gamma_sh, Num_chi, chi_face, a_loc, dln_a_dR_loc, beta_sh, A_eta_face)

    call eta_split_advection_faces(Num_chi,A_eta_face,adv_face_left,adv_face_right)
    call eta_diffusion_face_coeffs(Num_chi,deta,chi_face,a_loc,beta_sh,.false., &
                                   diff_face_left_base,diff_face_right_base)

    lower_base = zero
    diag_base = one
    upper_base = zero
    do I_chi = 1, Num_chi
        diag_base(I_chi) = diag_base(I_chi) + lambda_eta*adv_face_left(I_chi)
        if (I_chi < Num_chi) upper_base(I_chi) = upper_base(I_chi) + lambda_eta*adv_face_right(I_chi)
        if (I_chi > 1) then
            lower_base(I_chi) = lower_base(I_chi) - lambda_eta*adv_face_left(I_chi-1)
            diag_base(I_chi) = diag_base(I_chi) - lambda_eta*adv_face_right(I_chi-1)
        end if
    end do

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. active_hi*Num_chi >= 512) schedule(static) &
    !$OMP& private(I_gam_e,I_chi,kappa_face,lower,diag,upper,rhs,sol)
    do I_gam_e = 1, active_hi
        lower = lower_base
        diag = diag_base
        upper = upper_base
        do I_chi = 1, Num_chi
            if (I_chi < Num_chi) then
                kappa_face = 0.5d0*(kappa2_chi(I_gam_e,I_chi)+kappa2_chi(I_gam_e,I_chi+1))
                diag(I_chi) = diag(I_chi) + lambda_eta*kappa_face*diff_face_left_base(I_chi)
                upper(I_chi) = upper(I_chi) + lambda_eta*kappa_face*diff_face_right_base(I_chi)
            end if
            if (I_chi > 1) then
                kappa_face = 0.5d0*(kappa2_chi(I_gam_e,I_chi-1)+kappa2_chi(I_gam_e,I_chi))
                lower(I_chi) = lower(I_chi) - lambda_eta*kappa_face*diff_face_left_base(I_chi-1)
                diag(I_chi) = diag(I_chi) - lambda_eta*kappa_face*diff_face_right_base(I_chi-1)
            end if
        end do
        rhs = U_log(I_gam_e, :)
        rhs(1) = rhs(1) + dR_step*source_eta1(I_gam_e)

        call solve_tridiagonal(Num_chi, lower, diag, upper, rhs, sol)
        U_log(I_gam_e, :) = max(zero, sol)
    end do
    !$OMP END PARALLEL DO
end subroutine advance_eta_logchi_implicit

! PWN/CR η方向隐式输运：保守扩散，外边界可选零通量或自由逃逸。
subroutine advance_eta_logchi_pwncr_implicit(U_log, Num_gam_e, Num_chi, active_hi, deta, chi_face, Gamma_sh, &
                                             a_loc, dln_a_dR_loc, beta_sh, kappa2_chi, source_eta1, &
                                             dR_step, free_outer_escape, n_threads)
    integer, intent(in) :: Num_gam_e, Num_chi, active_hi, n_threads
    logical, intent(in) :: free_outer_escape
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: deta, chi_face(0:Num_chi), Gamma_sh, a_loc, dln_a_dR_loc, beta_sh
    real(8), intent(in) :: kappa2_chi(Num_gam_e,Num_chi), source_eta1(Num_gam_e), dR_step
    real(8) :: A_eta_face(1:Num_chi), adv_face_left(1:Num_chi), adv_face_right(1:Num_chi)
    real(8) :: diff_face_left_base(1:Num_chi), diff_face_right_base(1:Num_chi)
    real(8) :: lower_base(Num_chi), diag_base(Num_chi), upper_base(Num_chi)
    real(8) :: lower(Num_chi), diag(Num_chi), upper(Num_chi), rhs(Num_chi), sol(Num_chi)
    real(8) :: lambda_eta, kappa_face
    integer :: I_gam_e, I_chi

    lambda_eta = dR_step/deta
    call eta_face_transport_coeffs(Gamma_sh, Num_chi, chi_face, a_loc, dln_a_dR_loc, beta_sh, A_eta_face)

    call eta_split_advection_faces(Num_chi,A_eta_face,adv_face_left,adv_face_right)
    call eta_diffusion_face_coeffs(Num_chi,deta,chi_face,a_loc,beta_sh,free_outer_escape, &
                                   diff_face_left_base,diff_face_right_base)

    lower_base = zero
    diag_base = one
    upper_base = zero
    do I_chi = 1, Num_chi
        diag_base(I_chi) = diag_base(I_chi) + lambda_eta*adv_face_left(I_chi)
        if (I_chi < Num_chi) upper_base(I_chi) = upper_base(I_chi) + lambda_eta*adv_face_right(I_chi)
        if (I_chi > 1) then
            lower_base(I_chi) = lower_base(I_chi) - lambda_eta*adv_face_left(I_chi-1)
            diag_base(I_chi) = diag_base(I_chi) - lambda_eta*adv_face_right(I_chi-1)
        end if
    end do

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. active_hi*Num_chi >= 512) schedule(static) &
    !$OMP& private(I_gam_e,I_chi,kappa_face,lower,diag,upper,rhs,sol)
    do I_gam_e = 1, active_hi
        lower = lower_base
        diag = diag_base
        upper = upper_base
        do I_chi = 1, Num_chi
            if (I_chi < Num_chi) then
                kappa_face = 0.5d0*(kappa2_chi(I_gam_e,I_chi)+kappa2_chi(I_gam_e,I_chi+1))
                diag(I_chi) = diag(I_chi) + lambda_eta*kappa_face*diff_face_left_base(I_chi)
                upper(I_chi) = upper(I_chi) + lambda_eta*kappa_face*diff_face_right_base(I_chi)
            else if (free_outer_escape) then
                diag(I_chi) = diag(I_chi) + lambda_eta*kappa2_chi(I_gam_e,I_chi)*diff_face_left_base(I_chi)
            end if
            if (I_chi > 1) then
                kappa_face = 0.5d0*(kappa2_chi(I_gam_e,I_chi-1)+kappa2_chi(I_gam_e,I_chi))
                lower(I_chi) = lower(I_chi) - lambda_eta*kappa_face*diff_face_left_base(I_chi-1)
                diag(I_chi) = diag(I_chi) - lambda_eta*kappa_face*diff_face_right_base(I_chi-1)
            end if
        end do
        rhs = U_log(I_gam_e, :)
        rhs(1) = rhs(1) + dR_step*source_eta1(I_gam_e)
        call solve_tridiagonal(Num_chi, lower, diag, upper, rhs, sol)
        U_log(I_gam_e, :) = max(zero, sol)
    end do
    !$OMP END PARALLEL DO
end subroutine advance_eta_logchi_pwncr_implicit

! 能量维（log γ）冷却推进：隐式迎风格式，对每个χ柱独立求解三对角。
subroutine advance_energy_loggamma_chi(U_log, Num_gam_e, Num_chi, dEL_mean_chi, R_loc, d_x_E, dR_step, n_threads)
    integer, intent(in) :: Num_gam_e, Num_chi, n_threads
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: dEL_mean_chi(Num_gam_e-1, Num_chi), R_loc, d_x_E, dR_step

    real(8) :: coeff_xi(Num_gam_e-1), up(Num_gam_e-1), principal(Num_gam_e)
    real(8) :: temp1(Num_gam_e-1), rhs(Num_gam_e), sol(Num_gam_e), CFL, ad_coeff
    integer :: I_chi

    CFL = dR_step/d_x_E
    ad_coeff = one/(R_loc*dlog(ten))

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. Num_chi*Num_gam_e >= 512) schedule(static) &
    !$OMP& private(I_chi,coeff_xi,up,principal,temp1,rhs,sol)
    do I_chi = 1, Num_chi
        coeff_xi = dEL_mean_chi(:, I_chi) + ad_coeff
        up = -CFL*coeff_xi
        call electron_prepare_implicit_coeffs_common(Num_gam_e, one, up, principal, temp1)
        rhs = U_log(:, I_chi) / principal
        call electron_backward_sweep_common(Num_gam_e, temp1, rhs, sol)
        U_log(:, I_chi) = max(zero, sol)
        U_log(Num_gam_e, I_chi) = zero
    end do
    !$OMP END PARALLEL DO
end subroutine advance_energy_loggamma_chi

! PWN/CR能量维冷却：辐射损失加BM局部散度给出的χ依赖绝热项。
subroutine advance_energy_loggamma_chi_pwncr(U_log, Num_gam_e, Num_chi, dEL_mean_chi, adiabatic_log_coeff, &
                                             d_x_E, dR_step, n_threads)
    integer, intent(in) :: Num_gam_e, Num_chi, n_threads
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: dEL_mean_chi(Num_gam_e-1, Num_chi), adiabatic_log_coeff(Num_chi)
    real(8), intent(in) :: d_x_E, dR_step
    real(8) :: coeff_xi(Num_gam_e-1), up(Num_gam_e-1), principal(Num_gam_e)
    real(8) :: temp1(Num_gam_e-1), rhs(Num_gam_e), sol(Num_gam_e), CFL
    integer :: I_chi

    CFL = dR_step/d_x_E
    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. Num_chi*Num_gam_e >= 512) schedule(static) &
    !$OMP& private(I_chi,coeff_xi,up,principal,temp1,rhs,sol)
    do I_chi = 1, Num_chi
        coeff_xi = dEL_mean_chi(:, I_chi) + adiabatic_log_coeff(I_chi)
        up = -CFL*coeff_xi
        call electron_prepare_implicit_coeffs_common(Num_gam_e, one, up, principal, temp1)
        rhs = U_log(:, I_chi) / principal
        call electron_backward_sweep_common(Num_gam_e, temp1, rhs, sol)
        U_log(:, I_chi) = max(zero, sol)
        U_log(Num_gam_e, I_chi) = zero
    end do
    !$OMP END PARALLEL DO
end subroutine advance_energy_loggamma_chi_pwncr

! 随机再加速：log γ空间零通量保守扩散步，用于Strang分裂的半步或整步。
subroutine advance_energy_stochastic_loggamma_chi(U_log, Num_gam_e, Num_chi, stochastic_accel_norm, &
                                                  R_loc, d_x_E, dR_step, n_threads)
    integer, intent(in) :: Num_gam_e, Num_chi, n_threads
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: stochastic_accel_norm, R_loc, d_x_E, dR_step
    real(8) :: lower(Num_gam_e), diag(Num_gam_e), upper(Num_gam_e), rhs(Num_gam_e), sol(Num_gam_e)
    real(8) :: lambda_gamma
    integer :: I_chi, I_gam_e

    lambda_gamma = dR_step*stochastic_accel_norm/(R_loc*d_x_E*d_x_E)
    lower = zero
    diag = one
    upper = zero
    do I_gam_e = 1, Num_gam_e
        if (I_gam_e > 1) lower(I_gam_e) = -lambda_gamma
        if (I_gam_e < Num_gam_e) upper(I_gam_e) = -lambda_gamma
        if (I_gam_e == 1 .or. I_gam_e == Num_gam_e) then
            diag(I_gam_e) = one + lambda_gamma
        else
            diag(I_gam_e) = one + two*lambda_gamma
        end if
    end do

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. Num_chi*Num_gam_e >= 512) schedule(static) &
    !$OMP& private(I_chi,rhs,sol)
    do I_chi = 1, Num_chi
        rhs = U_log(:, I_chi)
        call solve_tridiagonal(Num_gam_e, lower, diag, upper, rhs, sol)
        U_log(:, I_chi) = max(zero, sol)
    end do
    !$OMP END PARALLEL DO
end subroutine advance_energy_stochastic_loggamma_chi

! 能量维冷却推进（特征线积分版）：对每个χ柱独立做电子冷却特征线更新。
subroutine advance_energy_loggamma_chi_charint(U_log, Num_gam_e, Num_chi, gam_e, DB_chi, dEl_chi, R_loc, &
                                               Gamma_sh, beta_sh, index_Y, dR_step, active_chi_hi, n_threads)
    integer, intent(in) :: Num_gam_e, Num_chi, index_Y, n_threads
    integer, intent(in), optional :: active_chi_hi
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: gam_e(Num_gam_e), DB_chi(Num_chi), dEl_chi(Num_gam_e, Num_chi)
    real(8), intent(in) :: R_loc, Gamma_sh, beta_sh, dR_step

    real(8) :: x_edge(Num_gam_e+1), U_in(Num_gam_e), U_out(Num_gam_e)
    real(8) :: a_rad, b_ad
    integer :: I_chi, chi_hi

    call electron_profile_log_cell_edges(Num_gam_e, gam_e, x_edge)

    chi_hi = Num_chi
    if (present(active_chi_hi)) chi_hi = max(1, min(Num_chi, active_chi_hi))
    do I_chi = 1, chi_hi
        U_in = U_log(:, I_chi)
        if (index_Y == 0) then
            a_rad = 1.35d-19*DB_chi(I_chi)**2/(max(beta_sh*Gamma_sh, tiny(one))*pi)
            b_ad = one/R_loc
            call electron_characteristic_update(Num_gam_e, dR_step, x_edge, electron_cooling_affine, &
                                                 a_rad, b_ad, gam_e, dEl_chi(:,I_chi), R_loc, zero, U_in, U_in, U_out)
        else
            call electron_characteristic_update(Num_gam_e, dR_step, x_edge, electron_cooling_piecewise, &
                                                 zero, zero, gam_e, dEl_chi(:,I_chi), R_loc, zero, U_in, U_in, U_out)
        end if
        U_log(:, I_chi) = U_out
    end do
end subroutine advance_energy_loggamma_chi_charint

end module electron_transport_2d_kernel
