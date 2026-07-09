! 电子 1D 三层隐式格式主驱动：BDF2 时间推进，前两步用单步隐式启动。
! 1D 3-level implicit electron driver: BDF2 time stepping, bootstrapped by single-step implicit updates.
subroutine fs_t2g1_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                            gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    use constants
    use dynamics_density_profile, only: density_profile
    use electron_common
    use electron_injection_profiles, only: source_edges
    use electron_radiation_kernel, only: syn_state, nua_fromtau
    use electron_cooling_kernel, only: forward_cooling
    use electron_transport_common, only: prepare_implicit_coeffs, backward_sweep, &
                                         dnx_dgamma
    implicit none
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer :: I_tobs,L,L1
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(out), dimension(Num_gam_e,Num_R) :: dN_gam_e
    real(8), intent(out), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_nu,Num_R) :: P_syn,Seed_syn
    real(8), intent(out), dimension(Num_R) :: V_m,V_c,V_a
    real(8) :: Eta_0,Epsilon_e,Epsilon_b,p,z,rho_ism,A_star,E_iso,tdur_log,f_e
    real(8) :: R_tr,f_jump,f_wide,R0,dNe,einit,DB,gemax,DB_min,gemax0
    real(8) :: temp_gam,gm,gc,d_x,rloc,gloc,gmp,rho_shell
    real(8) :: beta,f_r,dDR,dDD,CFL,temp,dbstep,gmax_step,gm_step,gmp_step,Q
    real(8), dimension(Num_nu) :: pemit,tau

    real(8), allocatable, dimension(:) :: dEl,delmean,delbase,principal,x,dF1,up,spec,spec_prev,xedge,temp1 &
        & ,temp2,temp3
    allocate (dEl(Num_gam_e),delmean(Num_gam_e-1),delbase(Num_gam_e-1),principal(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e), &
              up(Num_gam_e-1),spec(Num_gam_e),spec_prev(Num_gam_e),xedge(Num_gam_e+1), &
              temp1(Num_gam_e-1),temp2(Num_gam_e),temp3(Num_gam_e-1))

    call electron_unpack_boundary(Boundary,n,Eta_0,Epsilon_e,Epsilon_b,p,z,rho_ism,A_star, &
                                  E_iso,tdur_log,f_e,R_tr,f_jump,f_wide,R0)

    P_syn=0d0
    Seed_syn=0d0
    V_m=0d0
    V_c=0d0
    V_a=0d0

    call electron_initial_density(A_star,rho_ism,R(1),R0,R_tr,f_jump,f_wide,dNe,einit)

    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-1d0)))
    gemax=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-1d0)))
    gemax0=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-1d0)
    call electron_gm_exact(p,temp_gam,gemax,gm)
    gc=7.7d8/(1d0+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/2d0)
    call electron_initialize_spectrum(Num_gam_e,gemax0,f_e*einit,p,gm,gc,gemax, &
                                      imodelog,gam_e,spec,xedge)
    call dnx_dgamma(Num_gam_e,xedge,gam_e,spec,dN_gam_e(:,1))
    spec_prev = spec
    d_x=dlog(gam_e(2)/gam_e(1))

    do I_tobs=2,Num_R
        call prepare_t2g1_shell(I_tobs)
        if (I_tobs == 2) then
            spec=dN_gam_e(:,1)*gam_e
            spec_prev = spec
        end if
        call write_t2g1_cooling(I_tobs)

        do L=1,L1
            call advance_t2g1_substep(I_tobs,L)
        end do
    end do

    deallocate(dEl,delmean,delbase,principal,x,dF1,up,spec,spec_prev,xedge,temp1,temp2,temp3)

contains

    ! 准备当前输出壳层的动力学和注入尺度。
    ! Prepare the current output shell dynamics and injection scales.
    subroutine prepare_t2g1_shell(I_tobs)
        implicit none
        integer, intent(in) :: I_tobs

        rloc=R(I_tobs-1)
        gloc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/2d0
        call density_profile(A_star,rho_ism,rloc,R0,0,R_tr,f_jump,f_wide,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(gloc*(gloc-1d0)))
        gemax=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(gloc-1d0)
        call electron_gm_exact(p,temp_gam,gemax,gm)
        gmp=(1d0-p)/(gemax**(1d0-p)-gm**(1d0-p))
        gc=7.7d8*(1d0+z)/gloc/DB**2/R_Tobs(I_tobs)
        rho_shell=dNe

        beta=sqrt(1d0-1d0/gloc**2)
        f_r=(1.35d-19)/beta/gloc*DB**2/pi
        dDR=0.1/(f_r*gemax+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=max(100,min(1000,Int(dDD/dDR)))
        dDR=dDD/L1
        CFL=dDR/d_x
    end subroutine prepare_t2g1_shell

    ! 写同步辐射诊断，同时冻结本壳层冷却系数。
    ! Write synchrotron diagnostics and freeze the cooling coefficients for this shell.
    subroutine write_t2g1_cooling(I_tobs)
        implicit none
        integer, intent(in) :: I_tobs

        V_m(I_tobs-1)=4.2d6*DB*gm*gm/(gloc*(1d0-beta)*(1d0+z))
        V_c(I_tobs-1)=4.2d6*DB*gc*gc/(gloc*(1d0-beta)*(1d0+z))

        call syn_state(index_syn_intger,rloc,DB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,I_tobs-1),V_seed,pemit, &
                                    P_syn(:,I_tobs),Seed_syn(:,I_tobs),tau)
        call nua_fromtau(Num_nu,V_seed,tau,temp)
        V_a(I_tobs-1)=temp/(gloc*(1d0-beta)*(1d0+z))

        call forward_cooling(1,index_Y,Epsilon_e,Epsilon_b,p,DB,gm,gc,gemax,rloc, &
                             gloc,beta,dNe,Num_gam_e,Num_nu,1,n_threads,gam_e,V_seed, &
                             P_syn(:,I_tobs),Seed_syn(:,I_tobs),Seed_syn(:,I_tobs),temp2,dEl)

        delmean=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/2d0
        delbase=delmean
    end subroutine write_t2g1_cooling

    ! 一个 T2G1 子步：更新源项，组装隐式系数，然后后扫求解。
    ! One T2G1 substep: update the source, assemble implicit coefficients, then solve by backward sweep.
    subroutine advance_t2g1_substep(I_tobs,L)
        implicit none
        integer, intent(in) :: I_tobs,L

        rloc=rloc+dDR

        call density_profile(A_star,rho_ism,rloc,R0,1,R_tr,f_jump,f_wide,dNe)
        dbstep=0.39d0*dsqrt(Epsilon_b*dNe*(gloc*(gloc-1d0)))
        gmax_step=3d0*Para_m_energy/dsqrt(8d0*dbstep*Para_e**3)
        temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(gloc-1d0)
        call electron_gm_exact(p,temp_gam,gmax_step,gm_step)
        gmp_step=(1d0-p)/(gmax_step**(1d0-p)-gm_step**(1d0-p))

        call electron_injection_prefactor(rloc-dDR,dDR,dNe,f_e,gmp_step,Q)
        call source_edges(Num_gam_e,xedge,gm_step,gmax_step,Q,p,dF1)
        if (rho_shell > 0d0) then
            delmean=delbase*(dNe/rho_shell)
        else
            delmean=delbase
        end if

        temp3=delmean+1d0/rloc
        up=-CFL*temp3

        if (I_tobs == 2 .and. L <= 2) then
            call prepare_implicit_coeffs(Num_gam_e,1d0,up,principal,temp1)
            temp2 = (spec + dDR * dF1) / principal
        else
            call prepare_implicit_coeffs(Num_gam_e,1.5d0,up,principal,temp1)
            temp2 = ( (2d0)*spec - 0.5d0*spec_prev + dF1 * dDR ) / principal
        end if

        call backward_sweep(Num_gam_e,temp1,temp2,x)

        spec_prev = spec
        spec = x

        if (L1 == L) then
            call dnx_dgamma(Num_gam_e,xedge,gam_e,spec,dN_gam_e(:,I_tobs))
        end if
    end subroutine advance_t2g1_substep
end subroutine fs_t2g1_1d
