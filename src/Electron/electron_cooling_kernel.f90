!f2py: skip
module electron_cooling_kernel
  use constants
  use electron_cooling_ssa_kernel, only: electron_cooling_ssa_loss, electron_cooling_ssa_loss_batch
  use electron_cooling_ic_kernel, only: electron_cooling_ic_loss
  use electron_cooling_y_kernel, only: electron_cooling_y_nakar, electron_cooling_y_fan
  private

  public :: get_forward_cooling, prepare_forward_cooling_aux, prepare_forward_cooling_aux_batch
  public :: assemble_forward_cooling_split, assemble_forward_cooling_split_batch

contains
! 根据index_Y准备正向激波冷却辅助量：IC数值积分或Nakar Y参数。
subroutine prepare_forward_cooling_aux(index_Y,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,Seed_syn,cooling_aux)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: cooling_aux(Num_gam_e)
real(8) :: P_syn_column(Num_nu,1),Seed_syn_column(Num_nu,1),cooling_aux_column(Num_gam_e,1)

    P_syn_column(:,1)=P_syn
    Seed_syn_column(:,1)=Seed_syn
    call prepare_forward_cooling_aux_batch(index_Y,Num_gam_e,Num_nu,1,n_threads,gam_e,V_seed, &
                                           P_syn_column,Seed_syn_column,cooling_aux_column)
    cooling_aux=cooling_aux_column(:,1)
end subroutine prepare_forward_cooling_aux

! 批量版prepare_forward_cooling_aux：对多个χ列分别计算冷却辅助量。
subroutine prepare_forward_cooling_aux_batch(index_Y,Num_gam_e,Num_nu,Num_chi,n_threads,gam_e,V_seed,P_syn,Seed_syn,cooling_aux)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,Num_chi,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu,Num_chi),Seed_syn(Num_nu,Num_chi)
real(8), intent(out) :: cooling_aux(Num_gam_e,Num_chi)
integer :: I_chi

    cooling_aux=zero
    select case(index_Y)
    case(0)
    case(1)
        do I_chi=1,Num_chi
            call electron_cooling_ic_loss(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn(:,I_chi),cooling_aux(:,I_chi))
        end do
    case(2)
        do I_chi=1,Num_chi
            call electron_cooling_y_nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn(:,I_chi),cooling_aux(:,I_chi))
        end do
    case default
        error stop 'prepare_forward_cooling_aux_batch: index_Y must be 0, 1, or 2.'
    end select
end subroutine prepare_forward_cooling_aux_batch

! 组装正向激波冷却率 dγ/dR：index_Y=0为纯同步冷却，IC/Compton分支加入SSA回热。
subroutine assemble_forward_cooling_split(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                          beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                          Seed_syn_ssa,cooling_aux,dEl)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
real(8), intent(inout) :: Gam_e_max
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn_ssa(Num_nu),cooling_aux(Num_gam_e)
real(8), intent(out) :: dEl(Num_gam_e)
real(8) :: Compton(Num_gam_e),dot_gam_e_SSA(Num_gam_e)

    if (index_Y == 0) then
        dot_gam_e_SSA=zero
    else
        call electron_cooling_ssa_loss(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn_ssa,dot_gam_e_SSA)
    end if
    call assemble_forward_cooling_from_terms(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                             beta_Gam,dNe,Num_gam_e,gam_e,Compton,dot_gam_e_SSA,cooling_aux,dEl)
end subroutine assemble_forward_cooling_split

! 批量版assemble_forward_cooling_split：对多个χ列分别组装冷却率。
subroutine assemble_forward_cooling_split_batch(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                                beta_Gam,dNe,Num_gam_e,Num_nu,Num_chi,n_threads,gam_e,V_seed, &
                                                Seed_syn_ssa,cooling_aux,dEl)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,Num_chi,n_threads
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
real(8), intent(inout) :: Gam_e_max
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn_ssa(Num_nu,Num_chi),cooling_aux(Num_gam_e,Num_chi)
real(8), intent(out) :: dEl(Num_gam_e,Num_chi)
real(8) :: Compton(Num_gam_e),Gam_e_max_cell,dot_gam_e_SSA(Num_gam_e,Num_chi)
integer :: I_chi

    if (index_Y == 0) then
        dot_gam_e_SSA=zero
    else
        call electron_cooling_ssa_loss_batch(DB,Num_gam_e,Num_nu,Num_chi,n_threads,gam_e,V_seed,Seed_syn_ssa,dot_gam_e_SSA)
    end if
    do I_chi=1,Num_chi
       Gam_e_max_cell=Gam_e_max
       call assemble_forward_cooling_from_terms(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max_cell,R_loc,R_Gamma_loc, &
                                                beta_Gam,dNe,Num_gam_e,gam_e,Compton,dot_gam_e_SSA(:,I_chi), &
                                                cooling_aux(:,I_chi),dEl(:,I_chi))
    end do
end subroutine assemble_forward_cooling_split_batch

! 由各项组装正向激波冷却率：dγ/dt = (f_r*Y - SSA)*γ，支持4种Compton Y方案。
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
        dEl=f_r*gam_e
    case(1)
        dEl=(f_r+(cooling_aux-dot_gam_e_SSA)*ssa_scale)*gam_e
    case(2)
        Q=4d0*pi*R_loc*R_loc*para_c
        Compton=one+cooling_aux/Q/(4d0*R_Gamma_loc*R_Gamma_loc*dNe*Para_m_p_E)
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA*ssa_scale)*gam_e
    case(3)
        call electron_cooling_y_fan(Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,Compton)
        Compton=one+Compton
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA*ssa_scale)*gam_e
    case default
        error stop 'assemble_forward_cooling_from_terms: index_Y must be 0, 1, 2, or 3.'
    end select
end subroutine assemble_forward_cooling_from_terms

! 正向激波冷却主入口：准备IC辅助量，按index_Y选择是否计算SSA回热，再组装冷却率。
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
    if (index_Y == 0) then
        dot_gam_e_SSA=zero
    else
        call electron_cooling_ssa_loss(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn,dot_gam_e_SSA)
    end if
    call assemble_forward_cooling_from_terms(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                             beta_Gam,dNe,Num_gam_e,gam_e,Compton,dot_gam_e_SSA,cooling_aux,dEl)
end subroutine get_forward_cooling


end module electron_cooling_kernel
