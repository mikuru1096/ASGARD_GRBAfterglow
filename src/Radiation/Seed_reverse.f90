subroutine seed_reverse(T_cross,R_cross,e3_cross,gam20, Delta_t,b_r, &
                        Boundary,R_Tobs,R_gamma,R,gam_e,dN_gam_e,V_seed,n,Num_nu,Num_R,Num_gam_e,n_threads, P_syn_spec,seed_syn)
    !$ use omp_lib
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    !***********************************************************
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,n_threads
    real(8), intent(in) :: T_cross,R_cross,e3_cross,gam20, Delta_t,b_r
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_gamma(Num_R),R(Num_R)
    real(8), intent(in) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R),V_seed(Num_nu)
    real(8), intent(out) :: P_syn_spec(Num_nu,Num_R),seed_syn(Num_nu,Num_R)

    allocatable :: gam_e_pow(:)
    allocate(gam_e_pow(Num_gam_e))

    P_syn_spec=zero
    seed_syn=zero
    
    E_b = Boundary(6)
    dNe_ISM = Boundary(11)
    A_star = Boundary(12)
    f_e = Boundary(16)
    R0 = Boundary(n)
    gam_e_pow=gam_e*gam_e
    
    factor=(3.62d0/pi)**2
    Temp_syn=dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    
    beta4=sqrt(one-G00**(-2))
    Delta_0=Delta_t*para_c
    
    if (gam20 < one) then
        goto 100
    end if
    
    !$ call omp_set_dynamic(.true.)
    !$OMP PARALLEL num_threads(n_threads)
    !$OMP DO SIMD
    do I_R=1,Num_R
        Gam0=R_gamma(I_R)
        Beta=dsqrt(one-Gam0**(-2))
        Gth=Gam0-one
        Rariv2=R(I_R)**2
            
        if (A_star >= zero) then
            dNe_wind=A_star*3.0d35/Rariv2
            if (dNe_wind <= dNe_ISM/4.0d0) then
                dNe=dNe_ISM
            else
                dNe=dNe_wind
            end if
        else
            dNe=dNe_ISM
        end if
        
        if (R(I_R)<R0) then
            dNe=A_star*3.0d35/R0**2
        end if
        
        e2=4d0*Gam0*Gam0*dNe*Para_m_p*para_c*para_c
        if (R(I_R) < R_cross) then
            e3=e2
        else
            e3=e3_cross*(R_cross/R(I_R))**3*Gam0/gam20
        end if
        DB=dsqrt(8d0*pi*b_r*e3)
        do I_nu=1,Num_nu
            dInteg=zero
            Tau=zero
            do I_gam_e=1,Num_gam_e-1
                gam_e_mean2=(gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0
                Vc=(4.2d6)*gam_e_mean2*DB !Which is $\nu_c$
                x=V_seed(I_nu)/Vc !Which is ($\nu/\nu_c$)
                Fx=1.81d0*dexp(-x)/dsqrt(x**(-2d0/3d0)+factor) !Approximate function of synchrotron radiation spectrum
     !          Fx=2.149d0*x**(one/3.0d0)*dexp(-x) !!Another approximate function
                dN=(dN_gam_e(I_gam_e,I_R)+dN_gam_e(I_gam_e+1,I_R))/two
                dgam_e=gam_e(I_gam_e+1)-gam_e(I_gam_e)
                dInteg=dInteg+dN*Fx*dgam_e
                !====================  [SSA]  ======================
                ddN=dN_gam_e(I_gam_e,I_R)/gam_e_pow(I_gam_e)-dN_gam_e(I_gam_e+1,I_R)/gam_e_pow(I_gam_e+1)
                Tau=Tau+gam_e_mean2*ddN*Fx
                !Synchrotron self absorption effect
                !===================================================
            end do
            P_v=Temp_syn*DB*dInteg ! with units in erg/Hz/cm^2/s
            Tau=1.025d4*Tau/Rariv2/4d0/pi*DB/V_seed(I_nu)**2 !Synchrotron self absorption effect
            if ((Tau-1d-4) < 1d-5) Tau=1d-4
            P_v=P_v*(one-dexp(-Tau))/Tau !Radiation transfer equation for the emission-absortion plasma

            seed_syn(I_nu,I_R)=seed_syn(I_nu,I_R)+P_v/(Rariv2*V_seed(I_nu))
            P_syn_spec(I_nu,I_R)=P_syn_spec(I_nu,I_R)+P_v
            !Power of synchrotron radiation that be emitted, intrinsic, send to main program
        end do
    end do
    !$OMP END DO SIMD
    !$OMP END PARALLEL

100 continue
    
    temp_para=4d0*pi*Para_c*Para_h
    seed_syn=seed_syn/temp_para
    !seed photons to form SSC radiation, intrinsic, send to module 'SSC_spec' and 'annihilation'
    
    deallocate(gam_e_pow)
    
    
    return
end subroutine seed_reverse
