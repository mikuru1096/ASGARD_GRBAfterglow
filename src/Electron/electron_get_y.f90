subroutine get_nu_a(R_loc, DB, gam_e, dN_gam_e, V_a, Num_gam_e)
    use electron_radiation_kernel, only: get_nu_a_impl => get_nu_a
    implicit none
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: R_loc, DB
    real(8), intent(in) :: gam_e(Num_gam_e), dN_gam_e(Num_gam_e)
    real(8), intent(out) :: V_a

    call get_nu_a_impl(R_loc, DB, Num_gam_e, gam_e, dN_gam_e, V_a)
end subroutine get_nu_a

subroutine get_syn_selected(index_syn_intger, R_loc, DB, n_threads, gam_e, dN_gam_e, V_seed, P_syn, Seed_syn, Num_gam_e, Num_nu)
    use electron_radiation_kernel, only: get_syn_selected_impl => get_syn_selected
    implicit none
    integer, intent(in) :: index_syn_intger, n_threads, Num_gam_e, Num_nu
    real(8), intent(in) :: R_loc, DB
    real(8), intent(in) :: gam_e(Num_gam_e), dN_gam_e(Num_gam_e), V_seed(Num_nu)
    real(8), intent(out) :: P_syn(Num_nu), Seed_syn(Num_nu)

    call get_syn_selected_impl(index_syn_intger, R_loc, DB, Num_gam_e, Num_nu, n_threads, gam_e, dN_gam_e, V_seed, P_syn, Seed_syn)
end subroutine get_syn_selected

subroutine get_syn_transfer(R_loc, DB, n_threads, gam_e, dN_gam_e, V_seed, Transfer_syn, Num_gam_e, Num_nu)
    use electron_radiation_kernel, only: get_syn_state
    implicit none
    integer, intent(in) :: n_threads, Num_gam_e, Num_nu
    real(8), intent(in) :: R_loc, DB
    real(8), intent(in) :: gam_e(Num_gam_e), dN_gam_e(Num_gam_e), V_seed(Num_nu)
    real(8), intent(out) :: Transfer_syn(Num_nu)
    real(8) :: P_emit(Num_nu), P_syn(Num_nu), Seed_syn(Num_nu), Tau_syn(Num_nu)
    integer :: I_nu

    call get_syn_state(R_loc, DB, Num_gam_e, Num_nu, n_threads, gam_e, dN_gam_e, V_seed, &
                       P_emit, P_syn, Seed_syn, Tau_syn)
    do I_nu = 1, Num_nu
        if (P_emit(I_nu) > 0.0d0 .and. P_syn(I_nu) >= 0.0d0) then
            Transfer_syn(I_nu) = P_syn(I_nu) / P_emit(I_nu)
        else
            Transfer_syn(I_nu) = 1.0d0
        endif
    enddo
end subroutine get_syn_transfer
