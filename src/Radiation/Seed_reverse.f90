! 计算反向激波区域的同步辐射谱及其种子光子场，包含外介质密度剖面。
subroutine seed_reverse(T_cross,R_cross,e3_cross,gam20, Delta_t,b_r, &
                        Boundary,R_Tobs,R_gamma,R,gam_e,dN_gam_e,V_seed,n,Num_nu,Num_R,Num_gam_e,n_threads, P_syn_spec,seed_syn)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    use radiation_common
    IMPLICIT REAL(8)(A-H,O-Z)
    !***********************************************************
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,n_threads
    integer :: I_R
    real(8), intent(in) :: T_cross,R_cross,e3_cross,gam20, Delta_t,b_r
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_gamma(Num_R),R(Num_R)
    real(8), intent(in) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R),V_seed(Num_nu)
    real(8), intent(out) :: P_syn_spec(Num_nu,Num_R),seed_syn(Num_nu,Num_R)

    allocatable :: P_emit_shell(:), Tau_shell(:)
    allocate(P_emit_shell(Num_nu),Tau_shell(Num_nu))

    P_syn_spec=zero
    seed_syn=zero
    
    E_b = Boundary(6)
    dNe_ISM = Boundary(11)
    A_star = Boundary(12)
    f_e = Boundary(16)
    R0 = Boundary(n)
    Delta_0=Delta_t*para_c
    
    if (gam20 < one) then
        goto 100
    end if
    
    do I_R=1,Num_R
        Gam0=R_gamma(I_R)
        call dynamics_external_density_profile(A_star,dNe_ISM,R(I_R),R0,0,one,one,one,dNe)
        
        e2=4d0*Gam0*Gam0*dNe*Para_m_p*para_c*para_c
        if (R(I_R) < R_cross) then
            e3=e2
        else
            e3=e3_cross*(R_cross/R(I_R))**3*Gam0/gam20
        end if
        DB=dsqrt(8d0*pi*b_r*e3)
        call radiation_syn_seed_core(R(I_R),DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e(:,I_R),V_seed,1.025d4, &
                                     P_emit_shell,P_syn_spec(:,I_R),seed_syn(:,I_R),Tau_shell)
    end do

100 continue

    deallocate(P_emit_shell,Tau_shell)
    
    
    return
end subroutine seed_reverse
