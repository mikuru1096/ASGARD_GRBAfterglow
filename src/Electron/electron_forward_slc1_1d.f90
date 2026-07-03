! 电子 1D 半拉格朗日格式主驱动：初始化电子谱，再逐壳层更新辐射、冷却和谱输运。
! 1D semi-Lagrangian electron driver: initialize the spectrum, then update radiation, cooling, and spectral transport shell by shell.
subroutine fs_slc1_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                            gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    use constants
    use dynamics_density_profile, only: density_profile
    use electron_common
    use electron_transport_common, only: semi_lagrangian_step, dnx_dgamma
    use electron_injection_profiles, only: source_edges
    use electron_radiation_kernel, only: syn_state, nua_fromtau
    use electron_cooling_kernel, only: get_forward_cooling
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(out), dimension(Num_gam_e,Num_R) :: dN_gam_e
    real(8), intent(out), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_nu,Num_R) :: P_syn,Seed_syn
    real(8), intent(out), dimension(Num_R) :: V_m,V_c,V_a

    real(8), allocatable, dimension(:) :: dEl,delmean,delbase,spec,spec_step,dF1,xedge
    real(8), dimension(Num_nu) :: pemit,tau
    logical :: uniform

    allocate(dEl(Num_gam_e),delmean(Num_gam_e-1),delbase(Num_gam_e-1),spec(Num_gam_e), &
             spec_step(Num_gam_e),dF1(Num_gam_e),xedge(Num_gam_e+1))

    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,rho_ism,A_star, &
                                  E_iso,tdur_log,f_e,R_tr,f_jump,f_wide,R0)

    P_syn=0d0
    Seed_syn=0d0
    V_m=0d0
    V_c=0d0
    V_a=0d0

    call electron_initial_density(A_star,rho_ism,R_ini,R(1),R0,dNe,einit)
    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-1d0)))
    gemax=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-1d0)))
    gemax0=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(R_Gamma(1)-1d0)
    call electron_gm_exact(p,temp_gam,gemax,gm)
    gc=7.7d8/(1d0+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/2d0)
    call electron_initialize_spectrum(Num_gam_e,gemax0,einit,p,gm,gc,gemax, &
                                      imodelog,gam_e,spec,xedge)
    call dnx_dgamma(Num_gam_e,xedge,gam_e,spec,dN_gam_e(:,1))
    d_x=dlog10(gam_e(2)/gam_e(1))
    uniform=(A_star <= 0d0 .and. f_jump == 1d0)

    ! 每个输出壳层先固定辐射诊断，再在半拉格朗日子步中推进电子谱。
    ! Each output shell freezes diagnostics first, then advances the electron spectrum with semi-Lagrangian substeps.
    do I_tobs=2,Num_R
        rloc=R(I_tobs-1)
        gloc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/2d0
        call density_profile(A_star,rho_ism,rloc,R0,1,R_tr,f_jump,f_wide,dNe)

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

        V_m(I_tobs-1)=4.2d6*DB*gm*gm/(gloc*(1d0-beta)*(1d0+z))
        V_c(I_tobs-1)=4.2d6*DB*gc*gc/(gloc*(1d0-beta)*(1d0+z))

        call syn_state(index_syn_intger,rloc,DB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,I_tobs-1),V_seed,pemit,P_syn(:,I_tobs), &
                                    Seed_syn(:,I_tobs),tau)
        call nua_fromtau(Num_nu,V_seed,tau,temp)
        V_a(I_tobs-1)=temp/(gloc*(1d0-beta)*(1d0+z))
        call get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,gm,gc,gemax,rloc, &
                                 gloc,beta,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                 P_syn(:,I_tobs),Seed_syn(:,I_tobs),dEl)
        delbase=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/2d0/dlog(1d1)

        ! 子步中非均匀介质使用中点密度计算注入，端点密度推进到下一步。
        ! For nonuniform media, each substep uses midpoint density for injection and endpoint density for the next state.
        do L=1,L1
            rleft=rloc
            rho_left=dNe
            rright=rleft+dDR
            if (.not. uniform) then
                call density_profile(A_star,rho_ism,rright,R0,1,R_tr,f_jump,f_wide,rho_right)
                rmid=0.5d0*(rleft+rright)
                call density_profile(A_star,rho_ism,rmid,R0,1,R_tr,f_jump,f_wide,rho_mid)
            else
                rho_right=rho_left
                rmid=0.5d0*(rleft+rright)
                rho_mid=rho_left
            end if
            dbstep=0.39d0*dsqrt(Epsilon_b*rho_mid*(gloc*(gloc-1d0)))
            gmax_step=3d0*Para_m_energy/dsqrt(8d0*dbstep*Para_e**3)
            temp_gam=Epsilon_e/f_e*para_m_p/para_m_e*(gloc-1d0)
            call electron_gm_exact(p,temp_gam,gmax_step,gm_step)
            gmp_step=(1d0-p)/(gmax_step**(1d0-p)-gm_step**(1d0-p))
            call electron_injection_prefactor(rmid,dDR,rho_mid,f_e,gmp_step,Q)
            call source_edges(Num_gam_e,xedge,gm_step,gmax_step,Q,p,dF1)
            if (rho_shell > 0d0) then
                delmean=delbase*(rho_mid/rho_shell)
            else
                delmean=delbase
            end if
            call semi_lagrangian_step(Num_gam_e,dDR,d_x,delmean+1d0/rmid/dlog(1d1),dF1,spec,spec_step)
            spec=spec_step
            rloc=rright
            dNe=rho_right
        end do
        call dnx_dgamma(Num_gam_e,xedge,gam_e,spec,dN_gam_e(:,I_tobs))
    end do

    deallocate(dEl,delmean,delbase,spec,spec_step,dF1,xedge)
end subroutine fs_slc1_1d
