module electron_forward_kernel
    use constants
    use electron_cooling_kernel, only: get_forward_cooling
    implicit none
contains

subroutine electron_forward_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star,E_iso, &
                                            T_log10_duration,f_e,R_tr,f_jump,f_wide,R0)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: Boundary(n)
    real(8), intent(out) :: Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star,E_iso,T_log10_duration,f_e
    real(8), intent(out) :: R_tr,f_jump,f_wide,R0

    Eta_0=Boundary(1); R_ini=Boundary(4); Epsilon_e=Boundary(5); Epsilon_b=Boundary(6); p=Boundary(7); z=Boundary(8)
    dNe_ISM=Boundary(11); A_star=Boundary(12); E_iso=Boundary(14); T_log10_duration=Boundary(15); f_e=Boundary(16)
    R_tr=Boundary(21); f_jump=Boundary(22); f_wide=Boundary(23); R0=Boundary(n)
end subroutine electron_forward_unpack_boundary

subroutine electron_forward_cooling_step(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                                         R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                         P_syn,Seed_syn,dEl)
    implicit none
    integer, intent(in) :: index_Y,Num_gam_e,Num_nu,n_threads
    real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
    real(8), intent(inout) :: Gam_e_max
    real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu)
    real(8), intent(out) :: dEl(Num_gam_e)

    call get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc, &
                             R_Gamma_loc,beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                             P_syn,Seed_syn,dEl)
end subroutine electron_forward_cooling_step

subroutine electron_forward_implicit_source_step(Num_gam_e,dDR,dF1,dN_x,gam_e,dlog_ten,dN_gam_col)
    implicit none
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: dDR,dlog_ten,dF1(Num_gam_e),gam_e(Num_gam_e)
    real(8), intent(inout) :: dN_x(Num_gam_e)
    real(8), intent(out) :: dN_gam_col(Num_gam_e)

    dN_x=max(zero,dN_x+dDR*dF1)
    dN_gam_col=dN_x/gam_e/dlog_ten
end subroutine electron_forward_implicit_source_step

end module electron_forward_kernel
