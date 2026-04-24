!f2py: skip
module radiation_common
    use constants
    use dynamics_common, only: dynamics_external_density_profile
    implicit none
    real(8), parameter :: tau_transfer_floor = 1.0d-4
    real(8), parameter :: tau_transfer_guard = 1.0d-5

contains

subroutine compute_simpson_weights(weights, n)
    integer, intent(in) :: n
    integer :: i
    real(8), intent(out) :: weights(n)

    weights = 1.0d0
    if (n >= 3) then
        do i = 2, n - 1
            if (mod(i, 2) == 0) then
                weights(i) = 4.0d0
            else
                weights(i) = 2.0d0
            end if
        end do
    end if
end subroutine compute_simpson_weights

real(8) function radiation_powerlaw_interp(v0,v1,y0,y1,v)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: v0,v1,y0,y1,v
    real(8) :: slope

    if (v <= v0) then
        radiation_powerlaw_interp=y0
        return
    end if
    if (v >= v1) then
        radiation_powerlaw_interp=y1
        return
    end if

    if (v1 <= v0) then
        radiation_powerlaw_interp=0.5d0*(y0+y1)
    else if (y0 > zero .and. y1 > zero) then
        slope=dlog(y1/y0)/dlog(v1/v0)
        radiation_powerlaw_interp=y0*(v/v0)**slope
    else
        radiation_powerlaw_interp=y0+(y1-y0)*(v-v0)/(v1-v0)
    end if
end function radiation_powerlaw_interp

subroutine radiation_transfer_factor(Tau, factor)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: Tau
    real(8), intent(out) :: factor
    real(8) :: Tau_eff

    Tau_eff = Tau
    if ((Tau_eff - tau_transfer_floor) < tau_transfer_guard) Tau_eff = tau_transfer_floor
    factor = (one - dexp(-Tau_eff)) / Tau_eff
end subroutine radiation_transfer_factor

subroutine radiation_prepare_annihilation_grid(V_seed, Num_nu, ep1, ep2, dVloc, V_mid)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_nu
    real(8), intent(in) :: V_seed(Num_nu)
    real(8), intent(out) :: ep1(1,Num_nu),ep2(Num_nu-1,1),dVloc(Num_nu-1),V_mid(Num_nu-1)
    real(8) :: x_seed(Num_nu)

    para_hEme=Para_h/para_m_energy
    ep1(1,:)=para_hEme*V_seed
    x_seed=dlog(V_seed)
    V_mid=dexp(0.5d0*(x_seed(1:Num_nu-1)+x_seed(2:Num_nu)))
    ep2(:,1)=para_hEme*V_mid
    dVloc=V_mid*(x_seed(2:Num_nu)-x_seed(1:Num_nu-1))
end subroutine radiation_prepare_annihilation_grid

elemental real(8) function radiation_pair_cross_section(s_center) result(sigma_pair)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: s_center
    real(8) :: beta_sq,beta_loc,log_term

    if (s_center <= one) then
        sigma_pair = zero
        return
    end if

    beta_sq = one-one/s_center
    beta_loc = dsqrt(beta_sq)
    log_term = dlog((one+beta_loc)/(one-beta_loc))
    sigma_pair = (3.0d0/16.0d0)*Para_sigmaT*(one-beta_sq) * &
                 ((3.0d0-beta_sq*beta_sq)*log_term-two*beta_loc*(two-beta_sq))
end function radiation_pair_cross_section

subroutine radiation_pair_tau_headon_segment(V_seed,Num_nu,Seed_target,dx_cm,Tau_pair)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_nu
    real(8), intent(in) :: V_seed(Num_nu),Seed_target(Num_nu),dx_cm
    real(8), intent(out) :: Tau_pair(Num_nu)
    integer :: I_nu,I_seg
    real(8) :: para_hEme,ep1,v_mid,dv_mid,seed_mid,s_center
    real(8) :: x_seed(Num_nu)

    Tau_pair = zero
    if (dx_cm <= zero) return

    para_hEme = Para_h/para_m_energy
    x_seed = dlog(V_seed)

    do I_nu = 1, Num_nu
        ep1 = para_hEme*V_seed(I_nu)
        do I_seg = 1, Num_nu-1
            v_mid = dexp(0.5d0*(x_seed(I_seg)+x_seed(I_seg+1)))
            dv_mid = v_mid*(x_seed(I_seg+1)-x_seed(I_seg))
            seed_mid = radiation_powerlaw_interp(V_seed(I_seg),V_seed(I_seg+1), &
                                                 Seed_target(I_seg),Seed_target(I_seg+1),v_mid)
            s_center = ep1*para_hEme*v_mid
            Tau_pair(I_nu) = Tau_pair(I_nu) + radiation_pair_cross_section(s_center)*seed_mid*dv_mid
        end do
        Tau_pair(I_nu) = dx_cm*Tau_pair(I_nu)
    end do
end subroutine radiation_pair_tau_headon_segment

subroutine radiation_syn_seed_core(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,ssa_prefactor, &
                                   P_emit,P_syn,Seed_syn,Tau_syn)
    !$ use omp_lib
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e,Num_nu,n_threads
    real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu),ssa_prefactor
    real(8), intent(out) :: P_emit(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu),Tau_syn(Num_nu)
    integer :: I_nu,I_gam_e,work_items
    real(8) :: factor,Temp_syn,Rariv2,temp_para
    real(8) :: V_cal,dInteg,Tau,Vc,x,ratio_v,Fx,P_v,gam_e_mean2,dN,dgam_e,ddN,temp_abs

    factor=(3.62d0/pi)**2
    Temp_syn=dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    Rariv2=R_loc*R_loc
    work_items=Num_nu*(Num_gam_e-1)
    if (n_threads <= 1 .or. work_items < 8192) then
        do I_nu=1,Num_nu
            V_cal=V_seed(I_nu)
            dInteg=zero
            Tau=zero
            do I_gam_e=1,Num_gam_e-1
                gam_e_mean2=(gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0
                Vc=(4.2d6)*gam_e_mean2*DB
                x=V_cal/Vc
                ratio_v=Vc/V_cal
                Fx=1.81d0*dexp(-x)/dsqrt(ratio_v**(2d0/3d0)+factor)
                dN=(dN_gam_e(I_gam_e)+dN_gam_e(I_gam_e+1))/two
                dgam_e=gam_e(I_gam_e+1)-gam_e(I_gam_e)
                dInteg=dInteg+dN*Fx*dgam_e
                ddN=dN_gam_e(I_gam_e)/(gam_e(I_gam_e)*gam_e(I_gam_e)) - &
                    dN_gam_e(I_gam_e+1)/(gam_e(I_gam_e+1)*gam_e(I_gam_e+1))
                Tau=Tau+gam_e_mean2*ddN*Fx
            end do
            P_v=Temp_syn*DB*dInteg
            Tau=ssa_prefactor*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal)
            P_emit(I_nu)=P_v
            Tau_syn(I_nu)=Tau
            call radiation_transfer_factor(Tau,temp_abs)
            P_syn(I_nu)=P_v*temp_abs
            Seed_syn(I_nu)=P_syn(I_nu)/(Rariv2*V_cal)
        end do
    else
        !$OMP PARALLEL num_threads(n_threads), private(I_nu,I_gam_e,V_cal,dInteg,Tau,Vc,x,ratio_v,Fx,P_v, &
        !$OMP& gam_e_mean2,dN,dgam_e,ddN,temp_abs)
        !$OMP DO SCHEDULE(STATIC)
        do I_nu=1,Num_nu
            V_cal=V_seed(I_nu)
            dInteg=zero
            Tau=zero
            do I_gam_e=1,Num_gam_e-1
                gam_e_mean2=(gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0
                Vc=(4.2d6)*gam_e_mean2*DB
                x=V_cal/Vc
                ratio_v=Vc/V_cal
                Fx=1.81d0*dexp(-x)/dsqrt(ratio_v**(2d0/3d0)+factor)
                dN=(dN_gam_e(I_gam_e)+dN_gam_e(I_gam_e+1))/two
                dgam_e=gam_e(I_gam_e+1)-gam_e(I_gam_e)
                dInteg=dInteg+dN*Fx*dgam_e
                ddN=dN_gam_e(I_gam_e)/(gam_e(I_gam_e)*gam_e(I_gam_e)) - &
                    dN_gam_e(I_gam_e+1)/(gam_e(I_gam_e+1)*gam_e(I_gam_e+1))
                Tau=Tau+gam_e_mean2*ddN*Fx
            end do
            P_v=Temp_syn*DB*dInteg
            Tau=ssa_prefactor*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal)
            P_emit(I_nu)=P_v
            Tau_syn(I_nu)=Tau
            call radiation_transfer_factor(Tau,temp_abs)
            P_syn(I_nu)=P_v*temp_abs
            Seed_syn(I_nu)=P_syn(I_nu)/(Rariv2*V_cal)
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    temp_para=4d0*pi*Para_c*Para_h
    Seed_syn=Seed_syn/temp_para
end subroutine radiation_syn_seed_core

subroutine radiation_external_density(A_star,dNe_ISM,R_loc,R0,dNe)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: A_star,dNe_ISM,R_loc,R0
    real(8), intent(out) :: dNe

    call dynamics_external_density_profile(A_star,dNe_ISM,R_loc,R0,0,one,one,one,dNe)
end subroutine radiation_external_density

end module radiation_common
