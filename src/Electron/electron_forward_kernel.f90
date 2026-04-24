module electron_forward_kernel
    use constants
    use electron_cooling_kernel, only: get_SSA_numerical, get_SSA_numerical_batch, get_IC_numerical
    use electron_y_kernel, only: get_Y_Nakar, get_Y_Fan
    implicit none
    private

    public :: get_forward_cooling, prepare_forward_cooling_aux, prepare_forward_cooling_aux_batch
    public :: assemble_forward_cooling_split, assemble_forward_cooling_split_batch
contains

subroutine prepare_forward_cooling_aux(index_Y,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,Seed_syn,cooling_aux)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: cooling_aux(Num_gam_e)

    cooling_aux=zero
    select case(index_Y)
    case(1)
        call get_IC_numerical(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn,cooling_aux)
    case(2)
        call get_Y_Nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,cooling_aux)
    case default
    end select
end subroutine prepare_forward_cooling_aux

subroutine prepare_forward_cooling_aux_batch(index_Y,Num_gam_e,Num_nu,Num_chi,n_threads,gam_e,V_seed,P_syn,Seed_syn,cooling_aux)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,Num_chi,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu,Num_chi),Seed_syn(Num_nu,Num_chi)
real(8), intent(out) :: cooling_aux(Num_gam_e,Num_chi)
integer :: I_chi

    cooling_aux=zero
    select case(index_Y)
    case(1)
        do I_chi=1,Num_chi
            call get_IC_numerical(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn(:,I_chi),cooling_aux(:,I_chi))
        end do
    case(2)
        do I_chi=1,Num_chi
            call get_Y_Nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn(:,I_chi),cooling_aux(:,I_chi))
        end do
    case default
    end select
end subroutine prepare_forward_cooling_aux_batch

subroutine assemble_forward_cooling_split(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                          beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                          P_syn_ic,Seed_syn_ic,Seed_syn_ssa,cooling_aux,dEl)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
real(8), intent(inout) :: Gam_e_max
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn_ic(Num_nu),Seed_syn_ic(Num_nu), &
                       Seed_syn_ssa(Num_nu),cooling_aux(Num_gam_e)
real(8), intent(out) :: dEl(Num_gam_e)
real(8) :: Compton(Num_gam_e),dot_gam_e_SSA(Num_gam_e)

    call get_SSA_numerical(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn_ssa,dot_gam_e_SSA)
    call assemble_forward_cooling_from_terms(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                             beta_Gam,dNe,Num_gam_e,gam_e,Compton,dot_gam_e_SSA,cooling_aux,dEl)
end subroutine assemble_forward_cooling_split

subroutine assemble_forward_cooling_split_batch(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                                beta_Gam,dNe,Num_gam_e,Num_nu,Num_chi,n_threads,gam_e,V_seed, &
                                                P_syn_ic,Seed_syn_ic,Seed_syn_ssa,cooling_aux,dEl)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,Num_chi,n_threads
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
real(8), intent(inout) :: Gam_e_max
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn_ic(Num_nu,Num_chi),Seed_syn_ic(Num_nu,Num_chi), &
                       Seed_syn_ssa(Num_nu,Num_chi),cooling_aux(Num_gam_e,Num_chi)
real(8), intent(out) :: dEl(Num_gam_e,Num_chi)
real(8) :: Compton(Num_gam_e),Gam_e_max_cell,dot_gam_e_SSA(Num_gam_e,Num_chi)
integer :: I_chi

    call get_SSA_numerical_batch(DB,Num_gam_e,Num_nu,Num_chi,n_threads,gam_e,V_seed,Seed_syn_ssa,dot_gam_e_SSA)
    do I_chi=1,Num_chi
       Gam_e_max_cell=Gam_e_max
       call assemble_forward_cooling_from_terms(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max_cell,R_loc,R_Gamma_loc, &
                                                beta_Gam,dNe,Num_gam_e,gam_e,Compton,dot_gam_e_SSA(:,I_chi), &
                                                cooling_aux(:,I_chi),dEl(:,I_chi))
    end do
end subroutine assemble_forward_cooling_split_batch

subroutine assemble_forward_cooling_from_terms(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                               beta_Gam,dNe,Num_gam_e,gam_e,Compton,dot_gam_e_SSA,cooling_aux,dEl)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
real(8), intent(inout) :: Gam_e_max
real(8), intent(in) :: gam_e(Num_gam_e)
real(8), intent(inout) :: Compton(Num_gam_e)
real(8), intent(in) :: dot_gam_e_SSA(Num_gam_e),cooling_aux(Num_gam_e)
real(8), intent(out) :: dEl(Num_gam_e)

    cooling_scale=one/(beta_Gam*R_Gamma_loc)
    ssa_scale=cooling_scale/para_c
    f_r=1.35d-19*DB**2*cooling_scale/pi

    select case(index_Y)
    case(0)
        dEl=(f_r-dot_gam_e_SSA*ssa_scale)*gam_e
    case(1)
        dEl=(f_r+(cooling_aux-dot_gam_e_SSA)*ssa_scale)*gam_e
    case(2)
        Q=4d0*pi*R_loc*R_loc*para_c
        Compton=one+cooling_aux/Q/(4d0*R_Gamma_loc*R_Gamma_loc*dNe*Para_m_p_E)
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA*ssa_scale)*gam_e
    case(3)
        call get_Y_Fan(Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,Compton)
        Compton=one+Compton
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA*ssa_scale)*gam_e
    case default
        print*, 'invalid Compton case, check your chosen model!'
        stop
    end select
end subroutine assemble_forward_cooling_from_terms

subroutine get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                               beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,Seed_syn,dEl)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
real(8), intent(inout) :: Gam_e_max
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dEl(Num_gam_e)
real(8) :: cooling_aux(Num_gam_e)
real(8) :: Compton(Num_gam_e),dot_gam_e_SSA(Num_gam_e)

    call prepare_forward_cooling_aux(index_Y,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,Seed_syn,cooling_aux)
    call get_SSA_numerical(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn,dot_gam_e_SSA)
    call assemble_forward_cooling_from_terms(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                             beta_Gam,dNe,Num_gam_e,gam_e,Compton,dot_gam_e_SSA,cooling_aux,dEl)
end subroutine get_forward_cooling

end module electron_forward_kernel
