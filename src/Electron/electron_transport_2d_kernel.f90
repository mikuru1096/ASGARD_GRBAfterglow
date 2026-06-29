!f2py: skip
module electron_transport_2d_kernel
  use constants
  use electron_transport_common, only: electron_prepare_implicit_coeffs_common, electron_backward_sweep_common, &
                             electron_prepare_conservative_remap_nonuniform, electron_ppm_prefix_eval_nonuniform, &
                             electron_characteristic_update, electron_cooling_affine, electron_cooling_piecewise
  use electron_injection_profiles, only: electron_profile_log_cell_edges
  implicit none
  private
  real(8), parameter :: sigma_strong_shock = 4d0

  public :: compute_q_geometry, compute_q_cell_geometry, get_shock_transport_state
  public :: compute_downstream_comoving_grid
  public :: compute_q_divergence
  public :: compute_q_step_limit
  public :: advance_q_implicit, advance_energy_loggamma_chi
  public :: advance_q_pwncr_implicit, advance_energy_loggamma_chi_pwncr
  public :: advance_energy_stochastic_loggamma_chi
  public :: advance_q_advection_charint, advance_q_diffusion_implicit
  public :: advance_energy_loggamma_chi_charint

contains

! 构建有限主动壳层质量坐标q：q=0为激波前沿，外边界停在强激波压缩给出的接触面量级。
subroutine compute_q_geometry(Num_chi, dq, q_face, q_grid)
    integer, intent(in) :: Num_chi
    real(8), intent(out) :: dq, q_face(0:Num_chi), q_grid(Num_chi)
    real(8) :: q_active
    integer :: I_chi

    q_active = one-(one-one/sigma_strong_shock)**sigma_strong_shock
    dq = q_active/dble(Num_chi)
    do I_chi = 0, Num_chi
        q_face(I_chi) = dble(I_chi)*dq
    end do
    do I_chi = 1, Num_chi
        q_grid(I_chi) = (dble(I_chi)-0.5d0)*dq
    end do
end subroutine compute_q_geometry

subroutine q_geometry_point(k_medium,R_shell,Gamma_f,q_loc,radius_cell,gamma_cell,beta_cell)
    integer, intent(in) :: k_medium
    real(8), intent(in) :: R_shell,Gamma_f,q_loc
    real(8), intent(out) :: radius_cell,gamma_cell,beta_cell
    real(8) :: q_tail,alpha,u_f,beta_f,w,chi_bm,log_x_bm,log_x_st
    real(8) :: x_ratio,u_bm,beta_st,u_st,y_cell,u_cell

    q_tail = one-q_loc
    alpha = dble(4-k_medium)/dble(3-k_medium)
    u_f = dsqrt(Gamma_f*Gamma_f-one)
    beta_f = u_f/Gamma_f
    w = u_f*u_f/(one+u_f*u_f)
    chi_bm = q_tail**(-alpha)
    log_x_bm = -(chi_bm-one)/(4d0*dble(4-k_medium)*Gamma_f*Gamma_f)
    log_x_st = dlog(q_tail)/(sigma_strong_shock*dble(3-k_medium))
    x_ratio = dexp(w*log_x_bm+(one-w)*log_x_st)
    u_bm = u_f/dsqrt(chi_bm)
    beta_st = beta_f*dexp(log_x_st)
    u_st = beta_st/dsqrt(one-beta_st*beta_st)
    y_cell = w*dasinh(u_bm)+(one-w)*dasinh(u_st)
    u_cell = dsinh(y_cell)
    gamma_cell = dsqrt(one+u_cell*u_cell)
    beta_cell = u_cell/gamma_cell
    radius_cell = R_shell*x_ratio
end subroutine q_geometry_point

subroutine compute_q_cell_geometry(Num_chi,k_medium,R_shell,Gamma_f,q_grid,radius_cell,gamma_cell,beta_cell,beta_rel_sh)
    integer, intent(in) :: Num_chi,k_medium
    real(8), intent(in) :: R_shell,Gamma_f,q_grid(Num_chi)
    real(8), intent(out) :: radius_cell(Num_chi),gamma_cell(Num_chi),beta_cell(Num_chi),beta_rel_sh(Num_chi)
    real(8) :: beta_shock, u_f, w, u_sh_bm, u_sh_st, y_sh, u_sh
    integer :: I_chi

    ! Inlined from shock_beta_from_front_gamma
    u_f = dsqrt(Gamma_f*Gamma_f-one)
    w = u_f*u_f/(one+u_f*u_f)
    u_sh_bm = dsqrt(two*Gamma_f*Gamma_f-one)
    u_sh_st = sigma_strong_shock*u_f/(sigma_strong_shock-one)
    y_sh = w*dasinh(u_sh_bm)+(one-w)*dasinh(u_sh_st)
    u_sh = dsinh(y_sh)
    beta_shock = u_sh/dsqrt(one+u_sh*u_sh)
    do I_chi = 1, Num_chi
        call q_geometry_point(k_medium,R_shell,Gamma_f,q_grid(I_chi), &
                              radius_cell(I_chi),gamma_cell(I_chi),beta_cell(I_chi))
        beta_rel_sh(I_chi) = (beta_shock-beta_cell(I_chi))/(one-beta_shock*beta_cell(I_chi))
    end do
end subroutine compute_q_cell_geometry

! 将q网格转换为物理深度和局域下游共动系空间网格。
subroutine compute_downstream_comoving_grid(Num_chi,k_medium,R_shell,Gamma_f,q_face,q_grid, &
                                            x_face,x_comov_face,x_comov_center,dx_comov, &
                                            radius_cell,gamma_cell,beta_cell,beta_rel_sh)
    integer, intent(in) :: Num_chi,k_medium
    real(8), intent(in) :: R_shell,Gamma_f,q_face(0:Num_chi),q_grid(Num_chi)
    real(8), intent(out) :: x_face(0:Num_chi),x_comov_face(0:Num_chi),x_comov_center(Num_chi),dx_comov(Num_chi)
    real(8), intent(out) :: radius_cell(Num_chi),gamma_cell(Num_chi),beta_cell(Num_chi),beta_rel_sh(Num_chi)
    real(8) :: radius_face(0:Num_chi),gamma_face,beta_face,dx_lab
    integer :: I_chi

    do I_chi = 0, Num_chi
        call q_geometry_point(k_medium,R_shell,Gamma_f,q_face(I_chi),radius_face(I_chi),gamma_face,beta_face)
        x_face(I_chi) = R_shell-radius_face(I_chi)
    end do
    call compute_q_cell_geometry(Num_chi,k_medium,R_shell,Gamma_f,q_grid,radius_cell,gamma_cell,beta_cell,beta_rel_sh)
    x_comov_face(0) = zero
    do I_chi = 1, Num_chi
        dx_lab = x_face(I_chi)-x_face(I_chi-1)
        dx_comov(I_chi) = gamma_cell(I_chi)*dx_lab
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

! q下游局部散度：用桥接半径/速度场计算(∇·v)/(3β_f c)，返回log10 γ方程系数。
subroutine compute_q_divergence(Num_chi,k_medium,R_loc,Gamma_f,beta_f,q_grid,adiabatic_log_coeff)
    integer, intent(in) :: Num_chi,k_medium
    real(8), intent(in) :: R_loc,Gamma_f,beta_f,q_grid(Num_chi)
    real(8), intent(out) :: adiabatic_log_coeff(Num_chi)
    real(8) :: radius_cell(Num_chi),gamma_cell(Num_chi),beta_cell(Num_chi),beta_rel_sh(Num_chi)
    real(8) :: dbeta_dr,div_over_c
    integer :: I_chi

    call compute_q_cell_geometry(Num_chi,k_medium,R_loc,Gamma_f,q_grid,radius_cell,gamma_cell,beta_cell,beta_rel_sh)
    dbeta_dr = (beta_cell(2)-beta_cell(1))/(radius_cell(2)-radius_cell(1))
    div_over_c = two*beta_cell(1)/radius_cell(1) + dbeta_dr
    adiabatic_log_coeff(1) = div_over_c/(3d0*beta_f*dlog(ten))
    do I_chi = 2, Num_chi-1
        dbeta_dr = (beta_cell(I_chi+1)-beta_cell(I_chi-1))/(radius_cell(I_chi+1)-radius_cell(I_chi-1))
        div_over_c = two*beta_cell(I_chi)/radius_cell(I_chi) + dbeta_dr
        adiabatic_log_coeff(I_chi) = div_over_c/(3d0*beta_f*dlog(ten))
    end do
    dbeta_dr = (beta_cell(Num_chi)-beta_cell(Num_chi-1))/(radius_cell(Num_chi)-radius_cell(Num_chi-1))
    div_over_c = two*beta_cell(Num_chi)/radius_cell(Num_chi) + dbeta_dr
    adiabatic_log_coeff(Num_chi) = div_over_c/(3d0*beta_f*dlog(ten))
end subroutine compute_q_divergence

real(8) function q_face_transport_coeff(k_medium,R_loc,q_loc)
    integer, intent(in) :: k_medium
    real(8), intent(in) :: R_loc,q_loc
    real(8) :: q_active

    q_active = one-(one-one/sigma_strong_shock)**sigma_strong_shock
    q_face_transport_coeff = dble(3-k_medium)*(q_active-q_loc)/R_loc
end function q_face_transport_coeff

subroutine q_face_transport_coeffs(k_medium,R_loc,Num_chi,q_face,A_q_face)
    integer, intent(in) :: k_medium,Num_chi
    real(8), intent(in) :: R_loc,q_face(0:Num_chi)
    real(8), intent(out) :: A_q_face(1:Num_chi)
    integer :: I_face

    do I_face = 1, Num_chi
        A_q_face(I_face) = q_face_transport_coeff(k_medium,R_loc,q_face(I_face))
    end do
end subroutine q_face_transport_coeffs

! 估计q方向对流子步长限制。
real(8) function compute_q_step_limit(Num_chi,k_medium,R_loc,dq,q_face,cfl_factor)
    integer, intent(in) :: Num_chi,k_medium
    real(8), intent(in) :: R_loc,dq,q_face(0:Num_chi),cfl_factor
    real(8) :: A_q_face,max_q_coeff
    integer :: I_chi

    max_q_coeff = zero
    do I_chi = 1, Num_chi
        A_q_face = q_face_transport_coeff(k_medium,R_loc,q_face(I_chi))
        max_q_coeff = max(max_q_coeff,dabs(A_q_face))
    end do
    compute_q_step_limit = huge(one)
    if (max_q_coeff > zero) compute_q_step_limit = cfl_factor*dq/max_q_coeff
end function compute_q_step_limit

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
            ! Inlined from eta_linear_back_position
            if (dabs(b_cell) <= 1d-30) then
                eta_trial = eta_cur - a_cell*s_rem
            else
                eta_trial = (eta_cur + a_cell/b_cell)*dexp(-b_cell*s_rem) - a_cell/b_cell
            end if
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

subroutine q_split_advection_faces(Num_chi,A_q_face,adv_face_left,adv_face_right)
    implicit none
    integer, intent(in) :: Num_chi
    integer :: I_face
    real(8), intent(in) :: A_q_face(1:Num_chi)
    real(8), intent(out) :: adv_face_left(1:Num_chi),adv_face_right(1:Num_chi)

    adv_face_left=zero
    adv_face_right=zero
    do I_face=1,Num_chi-1
        if (A_q_face(I_face) >= zero) then
            adv_face_left(I_face)=A_q_face(I_face)
        else
            adv_face_right(I_face)=A_q_face(I_face)
        end if
    end do
    if (A_q_face(Num_chi) > zero) adv_face_left(Num_chi)=A_q_face(Num_chi)
end subroutine q_split_advection_faces

real(8) function q_depth_inverse_metric(k_medium,R_shell,Gamma_f,q_loc)
    integer, intent(in) :: k_medium
    real(8), intent(in) :: R_shell,Gamma_f,q_loc
    real(8) :: q_tail,alpha,u_f,w,chi_bm,log_x_bm,log_x_st
    real(8) :: x_ratio,dchi_dq,dlog_x_bm_dq,dlog_x_st_dq,dlog_x_dq,dx_dq

    q_tail = one-q_loc
    alpha = dble(4-k_medium)/dble(3-k_medium)
    u_f = dsqrt(Gamma_f*Gamma_f-one)
    w = u_f*u_f/(one+u_f*u_f)
    chi_bm = q_tail**(-alpha)
    log_x_bm = -(chi_bm-one)/(4d0*dble(4-k_medium)*Gamma_f*Gamma_f)
    log_x_st = dlog(q_tail)/(sigma_strong_shock*dble(3-k_medium))
    x_ratio = dexp(w*log_x_bm+(one-w)*log_x_st)
    dchi_dq = alpha*chi_bm/q_tail
    dlog_x_bm_dq = -dchi_dq/(4d0*dble(4-k_medium)*Gamma_f*Gamma_f)
    dlog_x_st_dq = -one/(sigma_strong_shock*dble(3-k_medium)*q_tail)
    dlog_x_dq = w*dlog_x_bm_dq+(one-w)*dlog_x_st_dq
    dx_dq = x_ratio*dlog_x_dq
    q_depth_inverse_metric = one/(-R_shell*dx_dq)
end function q_depth_inverse_metric

subroutine q_diffusion_face_coeffs(Num_chi,dq,q_face,k_medium,R_loc,Gamma_f,beta_sh,include_outer_face, &
                                   diff_face_left_base,diff_face_right_base)
    implicit none
    integer, intent(in) :: Num_chi,k_medium
    integer :: I_face
    logical, intent(in) :: include_outer_face
    real(8), intent(in) :: dq,q_face(0:Num_chi),R_loc,Gamma_f,beta_sh
    real(8), intent(out) :: diff_face_left_base(1:Num_chi),diff_face_right_base(1:Num_chi)
    real(8) :: diff_prefactor,dqds

    diff_face_left_base=zero
    diff_face_right_base=zero
    do I_face=1,Num_chi-1
        dqds = q_depth_inverse_metric(k_medium,R_loc,Gamma_f,q_face(I_face))
        diff_prefactor = dqds*dqds/(beta_sh*para_c)
        diff_face_left_base(I_face) = diff_prefactor/dq
        diff_face_right_base(I_face) = -diff_prefactor/dq
    end do
    if (include_outer_face) then
        dqds = q_depth_inverse_metric(k_medium,R_loc,Gamma_f,q_face(Num_chi))
        diff_prefactor = dqds*dqds/(beta_sh*para_c)
        diff_face_left_base(Num_chi) = diff_prefactor/dq
    end if
end subroutine q_diffusion_face_coeffs

subroutine build_q_advection_base_matrix(Num_chi, lambda_q, A_q_face, lower_base, diag_base, upper_base)
    integer, intent(in) :: Num_chi
    real(8), intent(in) :: lambda_q, A_q_face(1:Num_chi)
    real(8), intent(out) :: lower_base(Num_chi), diag_base(Num_chi), upper_base(Num_chi)
    real(8) :: adv_face_left(1:Num_chi), adv_face_right(1:Num_chi)
    integer :: I_chi

    call q_split_advection_faces(Num_chi,A_q_face,adv_face_left,adv_face_right)
    lower_base = zero
    diag_base = one
    upper_base = zero
    do I_chi = 1, Num_chi
        diag_base(I_chi) = diag_base(I_chi) + lambda_q*adv_face_left(I_chi)
        if (I_chi < Num_chi) upper_base(I_chi) = upper_base(I_chi) + lambda_q*adv_face_right(I_chi)
    end do
    do I_chi = 2, Num_chi
        lower_base(I_chi) = lower_base(I_chi) - lambda_q*adv_face_left(I_chi-1)
        diag_base(I_chi) = diag_base(I_chi) - lambda_q*adv_face_right(I_chi-1)
    end do
end subroutine build_q_advection_base_matrix

subroutine add_q_diffusion_to_matrix(Num_chi, lambda_q, kappa_row, diff_face_left_base, &
                                     diff_face_right_base, free_outer_escape, lower, diag, upper)
    integer, intent(in) :: Num_chi
    logical, intent(in) :: free_outer_escape
    real(8), intent(in) :: lambda_q, kappa_row(Num_chi)
    real(8), intent(in) :: diff_face_left_base(1:Num_chi), diff_face_right_base(1:Num_chi)
    real(8), intent(inout) :: lower(Num_chi), diag(Num_chi), upper(Num_chi)
    real(8) :: kappa_face
    integer :: I_chi

    do I_chi = 1, Num_chi
        if (I_chi < Num_chi) then
            kappa_face = 0.5d0*(kappa_row(I_chi)+kappa_row(I_chi+1))
            diag(I_chi) = diag(I_chi) + lambda_q*kappa_face*diff_face_left_base(I_chi)
            upper(I_chi) = upper(I_chi) + lambda_q*kappa_face*diff_face_right_base(I_chi)
        else if (free_outer_escape) then
            diag(I_chi) = diag(I_chi) + lambda_q*kappa_row(I_chi)*diff_face_left_base(I_chi)
        end if
    end do
    do I_chi = 2, Num_chi
        kappa_face = 0.5d0*(kappa_row(I_chi-1)+kappa_row(I_chi))
        lower(I_chi) = lower(I_chi) - lambda_q*kappa_face*diff_face_left_base(I_chi-1)
        diag(I_chi) = diag(I_chi) - lambda_q*kappa_face*diff_face_right_base(I_chi-1)
    end do
end subroutine add_q_diffusion_to_matrix

subroutine build_q_transport_base_matrix(Num_chi,dq,q_face,k_medium,R_loc,Gamma_f,beta_sh,dR_step, &
                                         free_outer_escape,lambda_q,diff_face_left_base,diff_face_right_base, &
                                         lower_base,diag_base,upper_base)
    integer, intent(in) :: Num_chi,k_medium
    logical, intent(in) :: free_outer_escape
    real(8), intent(in) :: dq,q_face(0:Num_chi),R_loc,Gamma_f,beta_sh,dR_step
    real(8), intent(out) :: lambda_q,diff_face_left_base(1:Num_chi),diff_face_right_base(1:Num_chi)
    real(8), intent(out) :: lower_base(Num_chi),diag_base(Num_chi),upper_base(Num_chi)
    real(8) :: A_q_face(1:Num_chi)

    lambda_q = dR_step/dq
    call q_face_transport_coeffs(k_medium, R_loc, Num_chi, q_face, A_q_face)
    call q_diffusion_face_coeffs(Num_chi,dq,q_face,k_medium,R_loc,Gamma_f,beta_sh,free_outer_escape, &
                                 diff_face_left_base,diff_face_right_base)
    call build_q_advection_base_matrix(Num_chi,lambda_q,A_q_face,lower_base,diag_base,upper_base)
end subroutine build_q_transport_base_matrix

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

! q方向纯对流推进（特征线积分）：PPM重构+特征线回溯。
subroutine advance_q_advection_charint(U_log, Num_gam_e, Num_chi, active_hi, dq, q_face, &
                                       k_medium, R_loc, source_q1, dR_step)
    integer, intent(in) :: Num_gam_e, Num_chi, active_hi, k_medium
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: dq, q_face(0:Num_chi), R_loc, source_q1(Num_gam_e), dR_step

    real(8) :: A_q_face(0:Num_chi), q_back(0:Num_chi)
    real(8) :: U_q_in(Num_chi), U_q_out(Num_chi)
    real(8) :: ppm_left(Num_chi), ppm_right(Num_chi), prefix(0:Num_chi)
    integer :: I_gam_e, I_face

    call q_face_transport_coeffs(k_medium, R_loc, Num_chi, q_face, A_q_face(1:))
    A_q_face(0) = q_face_transport_coeff(k_medium, R_loc, q_face(0))
    call eta_trace_back_faces_piecewise(Num_chi, dR_step, q_face, A_q_face, q_back)

    do I_gam_e = 1, active_hi
        U_q_in = U_log(I_gam_e, :)
        U_q_in(1) = U_q_in(1) + dR_step*source_q1(I_gam_e)
        call electron_prepare_conservative_remap_nonuniform(Num_chi, q_face, U_q_in, ppm_left, ppm_right, prefix)
        do I_face = 1, Num_chi
            U_q_out(I_face) = &
                (electron_ppm_prefix_eval_nonuniform(Num_chi, q_face, U_q_in, &
                                                     ppm_left, ppm_right, prefix, q_back(I_face)) - &
                 electron_ppm_prefix_eval_nonuniform(Num_chi, q_face, U_q_in, &
                                                     ppm_left, ppm_right, prefix, q_back(I_face-1))) / dq
        end do
        U_log(I_gam_e, :) = U_q_out
    end do
end subroutine advance_q_advection_charint

! q方向纯扩散推进（隐式）：三对角求解扩散项 ∂_q(κ∂_q U)。
subroutine advance_q_diffusion_implicit(U_log, Num_gam_e, Num_chi, active_hi, dq, q_face, k_medium, &
                                        R_loc, Gamma_f, beta_sh, kappa2_chi, dR_step, n_threads)
    integer, intent(in) :: Num_gam_e, Num_chi, active_hi, k_medium, n_threads
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: dq, q_face(0:Num_chi), R_loc, Gamma_f
    real(8), intent(in) :: beta_sh, kappa2_chi(Num_gam_e,Num_chi), dR_step

    real(8) :: diff_face_left_base(1:Num_chi), diff_face_right_base(1:Num_chi)
    real(8) :: lower(Num_chi), diag(Num_chi), upper(Num_chi), rhs(Num_chi), sol(Num_chi)
    real(8) :: lambda_q
    integer :: I_gam_e

    lambda_q = dR_step/dq
    call q_diffusion_face_coeffs(Num_chi,dq,q_face,k_medium,R_loc,Gamma_f,beta_sh,.false., &
                                 diff_face_left_base,diff_face_right_base)

    do I_gam_e = 1, active_hi
        lower = zero
        diag = one
        upper = zero
        call add_q_diffusion_to_matrix(Num_chi,lambda_q,kappa2_chi(I_gam_e,:), &
                                       diff_face_left_base,diff_face_right_base,.false.,lower,diag,upper)
        rhs = U_log(I_gam_e, :)
        call solve_tridiagonal(Num_chi, lower, diag, upper, rhs, sol)
        U_log(I_gam_e, :) = max(zero, sol)
    end do
end subroutine advance_q_diffusion_implicit

! q方向对流+扩散联合隐式推进：迎风对流+中心扩散，三对角求解。
subroutine advance_q_implicit(U_log, Num_gam_e, Num_chi, active_hi, dq, q_face, k_medium, R_loc, Gamma_f, &
                              beta_sh, kappa2_chi, source_q1, dR_step, n_threads)
    integer, intent(in) :: Num_gam_e, Num_chi, active_hi, k_medium, n_threads
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: dq, q_face(0:Num_chi), R_loc, Gamma_f
    real(8), intent(in) :: beta_sh, kappa2_chi(Num_gam_e,Num_chi), source_q1(Num_gam_e)
    real(8), intent(in) :: dR_step

    real(8) :: diff_face_left_base(1:Num_chi), diff_face_right_base(1:Num_chi)
    real(8) :: lower_base(Num_chi), diag_base(Num_chi), upper_base(Num_chi)
    real(8) :: lower(Num_chi), diag(Num_chi), upper(Num_chi), rhs(Num_chi), sol(Num_chi)
    real(8) :: lambda_q
    integer :: I_gam_e

    call build_q_transport_base_matrix(Num_chi,dq,q_face,k_medium,R_loc,Gamma_f,beta_sh,dR_step,.false., &
                                       lambda_q,diff_face_left_base,diff_face_right_base, &
                                       lower_base,diag_base,upper_base)

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. active_hi*Num_chi >= 512) schedule(static) &
    !$OMP& private(I_gam_e,lower,diag,upper,rhs,sol)
    do I_gam_e = 1, active_hi
        lower = lower_base
        diag = diag_base
        upper = upper_base
        call add_q_diffusion_to_matrix(Num_chi,lambda_q,kappa2_chi(I_gam_e,:), &
                                       diff_face_left_base,diff_face_right_base,.false.,lower,diag,upper)
        rhs = U_log(I_gam_e, :)
        rhs(1) = rhs(1) + dR_step*source_q1(I_gam_e)

        call solve_tridiagonal(Num_chi, lower, diag, upper, rhs, sol)
        U_log(I_gam_e, :) = max(zero, sol)
    end do
    !$OMP END PARALLEL DO
end subroutine advance_q_implicit

! PWN/CR q方向隐式输运：保守扩散，q=1为最深扫掠质量边界。
subroutine advance_q_pwncr_implicit(U_log, Num_gam_e, Num_chi, active_hi, dq, q_face, k_medium, &
                                    R_loc, Gamma_f, beta_sh, kappa2_chi, source_q1, &
                                    dR_step, free_outer_escape, n_threads)
    integer, intent(in) :: Num_gam_e, Num_chi, active_hi, k_medium, n_threads
    logical, intent(in) :: free_outer_escape
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: dq, q_face(0:Num_chi), R_loc, Gamma_f, beta_sh
    real(8), intent(in) :: kappa2_chi(Num_gam_e,Num_chi), source_q1(Num_gam_e), dR_step
    real(8) :: diff_face_left_base(1:Num_chi), diff_face_right_base(1:Num_chi)
    real(8) :: lower_base(Num_chi), diag_base(Num_chi), upper_base(Num_chi)
    real(8) :: lower(Num_chi), diag(Num_chi), upper(Num_chi), rhs(Num_chi), sol(Num_chi)
    real(8) :: lambda_q
    integer :: I_gam_e

    call build_q_transport_base_matrix(Num_chi,dq,q_face,k_medium,R_loc,Gamma_f,beta_sh,dR_step,free_outer_escape, &
                                       lambda_q,diff_face_left_base,diff_face_right_base, &
                                       lower_base,diag_base,upper_base)

    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. active_hi*Num_chi >= 512) schedule(static) &
    !$OMP& private(I_gam_e,lower,diag,upper,rhs,sol)
    do I_gam_e = 1, active_hi
        lower = lower_base
        diag = diag_base
        upper = upper_base
        call add_q_diffusion_to_matrix(Num_chi,lambda_q,kappa2_chi(I_gam_e,:), &
                                       diff_face_left_base,diff_face_right_base,free_outer_escape,lower,diag,upper)
        rhs = U_log(I_gam_e, :)
        rhs(1) = rhs(1) + dR_step*source_q1(I_gam_e)
        call solve_tridiagonal(Num_chi, lower, diag, upper, rhs, sol)
        U_log(I_gam_e, :) = max(zero, sol)
    end do
    !$OMP END PARALLEL DO
end subroutine advance_q_pwncr_implicit

! 能量维（log γ）冷却推进：隐式迎风格式，对每个q柱独立求解三对角。
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
    end do
    !$OMP END PARALLEL DO
end subroutine advance_energy_loggamma_chi

! PWN/CR能量维冷却：辐射损失加q局部散度给出的柱依赖绝热项。
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
                                               Gamma_sh, beta_sh, index_Y, dR_step, active_chi_hi, n_threads, &
                                               source_q1)
    integer, intent(in) :: Num_gam_e, Num_chi, index_Y, n_threads
    integer, intent(in), optional :: active_chi_hi
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: gam_e(Num_gam_e), DB_chi(Num_chi), dEl_chi(Num_gam_e, Num_chi)
    real(8), intent(in) :: R_loc, Gamma_sh, beta_sh, dR_step
    real(8), intent(in), optional :: source_q1(Num_gam_e)

    real(8) :: x_edge(Num_gam_e+1), U_in(Num_gam_e), U_out(Num_gam_e), source_zero(Num_gam_e)
    real(8) :: a_rad, b_ad, shock_four_velocity
    integer :: I_chi, chi_hi, cooling_mode
    logical :: source_column

    call electron_profile_log_cell_edges(Num_gam_e, gam_e, x_edge)
    source_zero = zero

    chi_hi = Num_chi
    if (present(active_chi_hi)) chi_hi = max(1, min(Num_chi, active_chi_hi))
    do I_chi = 1, chi_hi
        U_in = U_log(:, I_chi)
        source_column = present(source_q1) .and. I_chi == 1
        cooling_mode = electron_cooling_piecewise
        a_rad = zero
        b_ad = zero
        if (index_Y == 0) then
            shock_four_velocity = beta_sh*Gamma_sh
            if (shock_four_velocity <= zero) error stop "advance_energy_loggamma_chi_charint requires beta_sh*Gamma_sh > 0"
            a_rad = 1.35d-19*DB_chi(I_chi)**2/(shock_four_velocity*pi)
            b_ad = one/R_loc
            cooling_mode = electron_cooling_affine
        end if
        if (source_column) then
            call electron_characteristic_update(Num_gam_e, dR_step, x_edge, cooling_mode, &
                                                 a_rad, b_ad, gam_e, dEl_chi(:,I_chi), R_loc, &
                                                 one, source_q1, U_in, U_out)
        else
            call electron_characteristic_update(Num_gam_e, dR_step, x_edge, cooling_mode, &
                                                a_rad, b_ad, gam_e, dEl_chi(:,I_chi), R_loc, &
                                                zero, source_zero, U_in, U_out)
        end if
        U_log(:, I_chi) = U_out
    end do
end subroutine advance_energy_loggamma_chi_charint

end module electron_transport_2d_kernel
