subroutine fs_electron_charint_2d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                  Num_chi,index_Y,index_syn_intger,n_threads, &
                                  gam_e,dN_gam_e,dN_gam_e_total,P_syn,Seed_syn,V_m,V_c,V_a)
    implicit real(8)(a-h,o-z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,Num_chi,index_Y,index_syn_intger,n_threads
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e, Num_chi, Num_R)
    real(8), intent(out) :: dN_gam_e_total(Num_gam_e, Num_R)
    real(8), intent(out) :: gam_e(Num_gam_e)
    real(8), intent(out) :: P_syn(Num_nu, Num_R)
    real(8), intent(out) :: Seed_syn(Num_nu, Num_R)
    real(8), intent(out) :: V_m(Num_R), V_c(Num_R), V_a(Num_R)

    call fs_electron_transport_2d_core(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                       Num_chi,index_Y,index_syn_intger,n_threads, &
                                       gam_e,dN_gam_e,dN_gam_e_total,P_syn,Seed_syn,V_m,V_c,V_a, &
                                       .true., 'charint_2d')
end subroutine fs_electron_charint_2d
