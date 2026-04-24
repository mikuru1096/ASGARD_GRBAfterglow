!f2py: skip
module electron_transport_2d_kernel
  use constants
  use electron_common, only: electron_prepare_implicit_coeffs, electron_backward_sweep, &
                             electron_prepare_conservative_remap_nonuniform, &
                             electron_ppm_prefix_eval_nonuniform, &
                             electron_characteristic_step_affine_u, electron_characteristic_step_piecewise_u, &
                             electron_characteristic_transport_affine_u, &
                             electron_characteristic_transport_piecewise_u
  use electron_injection_profiles, only: electron_profile_log_cell_edges
  implicit real(8)(a-h,o-z)
  private

  public :: compute_log_chi_geometry, get_shock_transport_state
  public :: compute_downstream_comoving_grid
  public :: bm_beta2_lab, bm_beta2_shock
  public :: advance_eta_logchi_implicit, advance_energy_loggamma_chi
  public :: advance_eta_logchi_advection_charint, advance_eta_logchi_diffusion_implicit
  public :: advance_energy_loggamma_chi_charint

contains

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

real(8) function bm_beta2_shock(Gamma_sh, chi_loc)
    real(8), intent(in) :: Gamma_sh, chi_loc
    real(8) :: beta_sh, beta_2

    beta_sh = dsqrt(one-one/Gamma_sh**2)
    beta_2 = bm_beta2_lab(Gamma_sh, chi_loc)
    bm_beta2_shock = (beta_sh-beta_2)/(one-beta_sh*beta_2)
end function bm_beta2_shock

subroutine eta_face_transport_coeffs(Gamma_sh, Num_chi, chi_face, a_loc, dln_a_dR_loc, beta_sh, A_eta_face)
    integer, intent(in) :: Num_chi
    real(8), intent(in) :: Gamma_sh, chi_face(0:Num_chi), a_loc, dln_a_dR_loc, beta_sh
    real(8), intent(out) :: A_eta_face(1:Num_chi)
    real(8) :: beta_2_sh
    integer :: I_face

    do I_face = 1, Num_chi
        beta_2_sh = bm_beta2_shock(Gamma_sh, chi_face(I_face))
        A_eta_face(I_face) = (a_loc*beta_2_sh/(chi_face(I_face)*beta_sh) + &
                             ((chi_face(I_face)-one)/chi_face(I_face))*dln_a_dR_loc) / dlog(ten)
    end do
end subroutine eta_face_transport_coeffs

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

subroutine advance_eta_logchi_advection_charint(U_log, Num_gam_e, Num_chi, active_hi, deta, eta_face, chi_face, &
                                                Gamma_sh, a_loc, dln_a_dR_loc, beta_sh, source_eta1, dR_step)
    integer, intent(in) :: Num_gam_e, Num_chi, active_hi
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: deta, eta_face(0:Num_chi), chi_face(0:Num_chi)
    real(8), intent(in) :: Gamma_sh, a_loc, dln_a_dR_loc, beta_sh, source_eta1(Num_gam_e), dR_step

    real(8) :: A_eta_face(1:Num_chi), eta_back(0:Num_chi)
    real(8) :: U_eta_in(Num_chi), U_eta_out(Num_chi)
    real(8) :: q_left(Num_chi), q_right(Num_chi), prefix(0:Num_chi)
    integer :: I_gam_e, I_face

    call eta_face_transport_coeffs(Gamma_sh, Num_chi, chi_face, a_loc, dln_a_dR_loc, beta_sh, A_eta_face)

    eta_back(0) = eta_face(0)
    do I_face = 1, Num_chi
        eta_back(I_face) = eta_face(I_face) - dR_step*A_eta_face(I_face)
    end do

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

subroutine advance_eta_logchi_diffusion_implicit(U_log, Num_gam_e, Num_chi, active_hi, deta, chi_face, Gamma_sh, a_loc, &
                                                 dln_a_dR_loc, beta_sh, kappa2_chi, dR_step)
    integer, intent(in) :: Num_gam_e, Num_chi, active_hi
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: deta, chi_face(0:Num_chi), Gamma_sh, a_loc, dln_a_dR_loc
    real(8), intent(in) :: beta_sh, kappa2_chi(Num_gam_e,Num_chi), dR_step

    real(8) :: diff_face_left_base(1:Num_chi), diff_face_right_base(1:Num_chi)
    real(8) :: lower(Num_chi), diag(Num_chi), upper(Num_chi), rhs(Num_chi), sol(Num_chi)
    real(8) :: lambda_eta, diff_prefactor, coeff_left, coeff_right, ln10, kappa_face
    integer :: I_gam_e, I_face, I_chi

    ln10 = dlog(ten)
    lambda_eta = dR_step/deta

    diff_face_left_base = zero
    diff_face_right_base = zero
    do I_face = 1, Num_chi-1
        diff_prefactor = a_loc*a_loc/(chi_face(I_face)*chi_face(I_face)*beta_sh*para_c*ln10*ln10)
        coeff_left = diff_prefactor/deta + 0.5d0*ln10*diff_prefactor
        coeff_right = -diff_prefactor/deta + 0.5d0*ln10*diff_prefactor
        diff_face_left_base(I_face) = coeff_left
        diff_face_right_base(I_face) = coeff_right
    end do

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

subroutine advance_eta_logchi_implicit(U_log, Num_gam_e, Num_chi, active_hi, deta, chi_face, Gamma_sh, a_loc, dln_a_dR_loc, &
                                       beta_sh, kappa2_chi, source_eta1, dR_step)
    integer, intent(in) :: Num_gam_e, Num_chi, active_hi
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: deta, chi_face(0:Num_chi), Gamma_sh, a_loc, dln_a_dR_loc
      real(8), intent(in) :: beta_sh, kappa2_chi(Num_gam_e,Num_chi), source_eta1(Num_gam_e)
    real(8), intent(in) :: dR_step

    real(8) :: A_eta_face(1:Num_chi), adv_face_left(1:Num_chi), adv_face_right(1:Num_chi)
    real(8) :: diff_face_left_base(1:Num_chi), diff_face_right_base(1:Num_chi)
      real(8) :: lower_base(Num_chi), diag_base(Num_chi), upper_base(Num_chi)
      real(8) :: lower(Num_chi), diag(Num_chi), upper(Num_chi), rhs(Num_chi), sol(Num_chi)
      real(8) :: lambda_eta, diff_prefactor, coeff_left, coeff_right, ln10, kappa_face
    integer :: I_gam_e, I_face, I_chi

    ln10 = dlog(ten)
    lambda_eta = dR_step/deta
    call eta_face_transport_coeffs(Gamma_sh, Num_chi, chi_face, a_loc, dln_a_dR_loc, beta_sh, A_eta_face)

    adv_face_left = zero
    adv_face_right = zero
    diff_face_left_base = zero
    diff_face_right_base = zero

    do I_face = 1, Num_chi-1
        if (A_eta_face(I_face) >= zero) then
            adv_face_left(I_face) = A_eta_face(I_face)
        else
            adv_face_right(I_face) = A_eta_face(I_face)
        end if
        diff_prefactor = a_loc*a_loc/(chi_face(I_face)*chi_face(I_face)*beta_sh*para_c*ln10*ln10)
        coeff_left = diff_prefactor/deta + 0.5d0*ln10*diff_prefactor
        coeff_right = -diff_prefactor/deta + 0.5d0*ln10*diff_prefactor
        diff_face_left_base(I_face) = coeff_left
        diff_face_right_base(I_face) = coeff_right
    end do

    if (A_eta_face(Num_chi) > zero) adv_face_left(Num_chi) = A_eta_face(Num_chi)

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
end subroutine advance_eta_logchi_implicit

subroutine advance_energy_loggamma_chi(U_log, Num_gam_e, Num_chi, dEL_mean_chi, R_loc, d_x_E, dR_step)
    integer, intent(in) :: Num_gam_e, Num_chi
    real(8), intent(inout) :: U_log(Num_gam_e, Num_chi)
    real(8), intent(in) :: dEL_mean_chi(Num_gam_e-1, Num_chi), R_loc, d_x_E, dR_step

    real(8) :: coeff_xi(Num_gam_e-1), up(Num_gam_e-1), principal(Num_gam_e)
    real(8) :: temp1(Num_gam_e-1), rhs(Num_gam_e), sol(Num_gam_e), CFL, ad_coeff
    integer :: I_chi

    CFL = dR_step/d_x_E
    ad_coeff = one/(R_loc*dlog(ten))

    do I_chi = 1, Num_chi
        coeff_xi = dEL_mean_chi(:, I_chi) + ad_coeff
        up = -CFL*coeff_xi
        call electron_prepare_implicit_coeffs(Num_gam_e, one, up, principal, temp1)
        rhs = U_log(:, I_chi) / principal
        call electron_backward_sweep(Num_gam_e, temp1, rhs, sol)
        U_log(:, I_chi) = max(zero, sol)
        U_log(Num_gam_e, I_chi) = zero
    end do
end subroutine advance_energy_loggamma_chi

subroutine advance_energy_loggamma_chi_charint(U_log, Num_gam_e, Num_chi, gam_e, DB_chi, dEl_chi, R_loc, &
                                               Gamma_sh, beta_sh, index_Y, dR_step, active_chi_hi)
    integer, intent(in) :: Num_gam_e, Num_chi, index_Y
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
            call electron_characteristic_transport_affine_u(Num_gam_e, dR_step, x_edge, a_rad, b_ad, U_in, U_out)
        else
            call electron_characteristic_transport_piecewise_u(Num_gam_e, dR_step, x_edge, gam_e, &
                                                               dEl_chi(:,I_chi), R_loc, U_in, U_out)
        end if
        U_log(:, I_chi) = U_out
    end do
end subroutine advance_energy_loggamma_chi_charint

end module electron_transport_2d_kernel
