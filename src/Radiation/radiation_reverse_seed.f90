! 计算反向激波区域的同步辐射谱及其种子光子场，包含外介质密度剖面。
subroutine seed_reverse(T_cross,R_cross,e3_cross,gam20, Delta_t,b_r, &
                        Boundary,R_Tobs,R_gamma,R,gam_e,dN_gam_e,V_seed,n,Num_nu,Num_R,Num_gam_e,n_threads, P_syn_spec,seed_syn)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_boundary_r0, dynamics_external_density_profile, dynamics_set_density_jump_profile
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

    Eta_0 = Boundary(1)
    E_b = Boundary(6)
    dNe_ISM = Boundary(11)
    A_star = Boundary(12)
    E_iso = Boundary(14)
    f_e = Boundary(16)
    call dynamics_boundary_r0(Boundary,n,R0)
    call dynamics_set_density_jump_profile(Boundary,n)
    Delta_0=Delta_t*para_c
    para_m_ej=E_iso/Eta_0/para_c**2

    if (gam20 >= one) then
        do I_R=1,Num_R
            Gam0=R_gamma(I_R)
            call dynamics_external_density_profile(A_star,dNe_ISM,R(I_R),R0,0,one,one,one,dNe)
            e3=reverse_region3_energy_density(I_R,Gam0)
            DB=dsqrt(8d0*pi*b_r*e3)
            call radiation_syn_seed_core(R(I_R),DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e(:,I_R),V_seed,1.025d4, &
                                         P_emit_shell,P_syn_spec(:,I_R),seed_syn(:,I_R),Tau_shell)
        end do
    end if

    deallocate(P_emit_shell,Tau_shell)

    return

contains

    real(8) function reverse_region3_energy_density(I_R,Gam0)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_R
    real(8), intent(in) :: Gam0

        if (R(I_R) < R_cross) then
            Delta=max(Delta_0,R(I_R)/Eta_0**2)
            u2=dsqrt(Gam0*Gam0-one)
            u4=dsqrt(Eta_0*Eta_0-one)
            gam34=(Gam0*Gam0+Eta_0*Eta_0-one)/(Eta_0*Gam0+u2*u4)
            para_n4=para_m_ej/(4d0*pi*Para_m_p*R(I_R)*R(I_R)*Eta_0*Delta)
            para_n3=(4d0*gam34+3d0)*para_n4
            reverse_region3_energy_density=(gam34-one)*para_n3*Para_m_p_E
        else
            reverse_region3_energy_density=e3_cross*(R_cross/R(I_R))**3*Gam0/gam20
        end if
    end function reverse_region3_energy_density
end subroutine seed_reverse
