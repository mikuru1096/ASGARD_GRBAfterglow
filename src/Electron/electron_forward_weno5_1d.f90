! 电子1D WENO5格式主驱动：三阶TVD Runge-Kutta + WENO5通量重构。
! 1D electron WENO5 driver: third-order TVD Runge-Kutta plus WENO5 flux reconstruction.
subroutine fs_weno5_1d(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,n_threads, &
                                gam_e,dN_gam_e,P_syn,Seed_syn)
    use constants
    use dynamics_density_profile, only: density_profile
    use electron_common
    use electron_injection_profiles, only: source_edges
    use rad_common, only: syn_seed_chi
    use electron_cooling_kernel, only: forward_cooling
    use electron_transport_common, only: dnx_dgamma
    implicit none
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,n_threads
    integer :: I_tobs,L,L1
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(out), dimension(Num_gam_e,Num_R) :: dN_gam_e
    real(8), intent(out), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_nu,Num_R) :: P_syn,Seed_syn
    real(8) :: Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star,E_iso,tdur_log,f_e
    real(8) :: R_tr,f_jump,f_wide,R0,dNe,einit,DB,DB_min,gmax,gmax0
    real(8) :: temp_gam,gm,gc,dxlog,cfl_target,rloc,gloc,gmp
    real(8) :: beta,f_r,dDR,dDD,maxspeed,CFL,Q

    real(8),allocatable,dimension (:) :: dEl,deladv,spec,xedge,dF1
    real(8),allocatable,dimension (:,:) :: stage
    real(8),allocatable,dimension (:) :: specext,fpext,flux,delext
    allocate (dEl(Num_gam_e),deladv(Num_gam_e),spec(Num_gam_e),xedge(Num_gam_e+1),dF1(Num_gam_e))
    allocate(specext(1-3:Num_gam_e+3),stage(3, 1-3:Num_gam_e+3),&
             fpext(1-3:Num_gam_e+3),flux(0:Num_gam_e),delext(1-3:Num_gam_e+3))

    ! 解包公开边界参数，内部只用短名参与公式。
    ! Unpack the public boundary vector once; formulas below use short internal names.
    call electron_unpack_boundary(Boundary,n,Eta_0,R_ini,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star, &
                                  E_iso,tdur_log,f_e,R_tr,f_jump,f_wide,R0)

    P_syn=0d0
    Seed_syn=0d0

    ! 初始 shock-front 电子谱使用解析注入近似，随后进入 WENO 输运。
    ! The initial shock-front spectrum uses the analytic injection approximation before WENO transport.
    call electron_initial_density(A_star,dNe_ISM,R_ini,R(1),R0,dNe,einit)

    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-1d0)))
    gmax=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-1d0)))
    gmax0=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    temp_gam=Epsilon_e/f_e*Para_m_p_DIV_m_e*(R_Gamma(1)-1d0)
    call electron_gm_exact(p,temp_gam,gmax,gm)
    gc=7.7d8/(1d0+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/2d0)
    call electron_initialize_spectrum(Num_gam_e,gmax0,einit,p,gm,gc,gmax,imodelog,gam_e,spec,xedge)
    call dnx_dgamma(Num_gam_e,xedge,gam_e,spec,dN_gam_e(:,1))

    dxlog=dlog(gam_e(2)/gam_e(1))
    cfl_target=0.8d0

    ! 每个半径壳层先冻结冷却速度，再用 RK-WENO 子步推进电子数。
    ! Each radial shell freezes the cooling speed, then advances electron counts with RK-WENO substeps.
    do I_tobs=2,Num_R
        call prepare_weno_shell(I_tobs)
        call write_weno_cooling(I_tobs)

        do L=1,L1
            call advance_weno_substep(L)
        end do
    end do
    call write_finaldiag()

    deallocate (dEl,deladv,spec,xedge,dF1,specext,stage,fpext,flux,delext)

contains

    ! 准备当前输出壳层的 shock 尺度和 WENO 子步长度。
    ! Prepare the current output shell shock scales and WENO substep length.
    subroutine prepare_weno_shell(I_tobs)
        implicit none
        integer, intent(in) :: I_tobs

        rloc=R(I_tobs-1)
        gloc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/2d0
        call density_profile(A_star,dNe_ISM,rloc,R0,0,R_tr,f_jump,f_wide,dNe)

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(gloc*(gloc-1d0)))
        gmax=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        temp_gam=Epsilon_e/f_e*Para_m_p_DIV_m_e*(gloc-1d0)
        call electron_gm_exact(p,temp_gam,gmax,gm)
        gmp=(1d0-p)/(gmax**(1d0-p)-gm**(1d0-p))
        gc=7.7d8*(1d0+z)/gloc/DB**2/R_Tobs(I_tobs)

        beta=dsqrt(1d0-1d0/gloc**2)
        f_r=(1.35d-19)/beta/gloc*DB**2/pi
        dDR=0.1d0/(f_r*gmax+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        spec=dN_gam_e(:,I_tobs-1)*gam_e
    end subroutine prepare_weno_shell

    ! 写壳层同步辐射输出，并把冷却速度转成对数能量坐标速度。
    ! Write shell synchrotron output and convert cooling speed to log-energy coordinate speed.
    subroutine write_weno_cooling(I_tobs)
        implicit none
        integer, intent(in) :: I_tobs
        real(8), dimension(1) :: DB_chi
        real(8), dimension(1) :: Weight_chi
        real(8), dimension(Num_gam_e,1) :: specchi
        real(8), dimension(Num_nu,1) :: pemit,pshell,seedshell,taushell

        DB_chi(1)=DB
        Weight_chi(1)=1d0
        specchi(:,1)=dN_gam_e(:,I_tobs-1)
        call syn_seed_chi(rloc,Num_gam_e,Num_nu,1,gam_e,specchi,V_seed,DB_chi,Weight_chi,1.046d4, &
                                               pemit,pshell,seedshell,taushell)
        P_syn(:,I_tobs)=pshell(:,1)
        Seed_syn(:,I_tobs)=seedshell(:,1)

        call forward_cooling(1,index_Y,Epsilon_e,Epsilon_b,p,DB,gm,gc,gmax,rloc, &
                             gloc,beta,dNe,Num_gam_e,Num_nu,1,n_threads,gam_e,V_seed, &
                             P_syn(:,I_tobs),Seed_syn(:,I_tobs),Seed_syn(:,I_tobs),deladv,dEl)

        dEl(1:Num_gam_e-1)=(dEl(1:Num_gam_e-1)+dEl(2:Num_gam_e))*0.5d0
        dEl(Num_gam_e)=dEl(Num_gam_e-1)*0.5d0
        deladv=(dEl+1d0/rloc)

        L1=Int(dDD/dDR)+10
        maxspeed=maxval(abs(deladv))
        L1=max(L1,ceiling(dDD*maxspeed/(cfl_target*dxlog)))
        dDR=dDD/L1
        CFL=dDR/dxlog
    end subroutine write_weno_cooling

    ! 最后一个输出点没有后续推进，需与最终电子谱同步刷新辐射数组。
    ! The final output point has no following advance, so refresh radiation from the final electron spectrum.
    subroutine write_finaldiag()
        implicit none
        real(8), dimension(1) :: DB_chi,Weight_chi
        real(8), dimension(Num_gam_e,1) :: specchi
        real(8), dimension(Num_nu,1) :: pemit,pshell,seedshell,taushell

        rloc=R(Num_R)
        gloc=R_Gamma(Num_R)
        call density_profile(A_star,dNe_ISM,rloc,R0,1,R_tr,f_jump,f_wide,dNe)
        DB=0.39d0*dsqrt(Epsilon_b*dNe*(gloc*(gloc-1d0)))

        DB_chi(1)=DB
        Weight_chi(1)=1d0
        specchi(:,1)=dN_gam_e(:,Num_R)
        call syn_seed_chi(rloc,Num_gam_e,Num_nu,1,gam_e,specchi,V_seed,DB_chi,Weight_chi,1.046d4, &
                                               pemit,pshell,seedshell,taushell)
        P_syn(:,Num_R)=pshell(:,1)
        Seed_syn(:,Num_R)=seedshell(:,1)
    end subroutine write_finaldiag

    ! 一个 RK-WENO 子步：更新注入项，重构通量，最后投影到非负电子数。
    ! One RK-WENO substep: update injection, reconstruct fluxes, then project electron counts nonnegative.
    subroutine advance_weno_substep(L)
        implicit none
        integer, intent(in) :: L
        integer :: j

        rloc=rloc+dDR

        call density_profile(A_star,dNe_ISM,rloc,R0,1,R_tr,f_jump,f_wide,dNe)

        deladv=(dEl+1d0/rloc)
        call electron_injection_prefactor(rloc,dDR,dNe,f_e,gmp,Q)
        call source_edges(Num_gam_e,xedge,gm,gmax,Q,p,dF1)

        call load_weno_state()

        do j=1,3
            call compute_weno_fluxes()
            call rk_weno_stage(j)
        end do

        spec = specext(1:Num_gam_e)
        where(spec < 0.0d0) spec = 0.0d0

        if (L1 == L) then
            call dnx_dgamma(Num_gam_e,xedge,gam_e,spec,dN_gam_e(:,I_tobs))
        end if
    end subroutine advance_weno_substep

    ! ghost cell 只服务 WENO stencil，不参与物理输出。
    ! Ghost cells only serve the WENO stencil; they are not physical output cells.
    subroutine load_weno_state()
        implicit none

        specext(1-3:0) = spec(1)
        specext(1:Num_gam_e) = spec
        specext(Num_gam_e+1:Num_gam_e+3) = spec(Num_gam_e)
        stage(1,:) = specext

        delext(1-3:0) = deladv(1)
        delext(1:Num_gam_e) = deladv
        delext(Num_gam_e+1:Num_gam_e+3) = deladv(Num_gam_e)
    end subroutine load_weno_state

    ! 按冷却速度方向选择左右偏 WENO stencil。
    ! Choose the left- or right-biased WENO stencil from the cooling-speed sign.
    subroutine compute_weno_fluxes()
        implicit none
        integer :: ig

        specext(1-3:0) = specext(1)
        specext(Num_gam_e+1:Num_gam_e+3) = specext(Num_gam_e)

        fpext = delext * specext

        do ig = 0, Num_gam_e
            if (delext(ig) <= 0.0d0) then
                flux(ig) = weno5_positive_flux(fpext(ig-2:ig+2))
            else
                flux(ig) = weno5_negative_flux(fpext(ig-1:ig+3))
            end if
        end do
    end subroutine compute_weno_fluxes

    ! 三阶 TVD RK stage，源项按当前子步注入量显式加入。
    ! Third-order TVD RK stage with explicit source injection for the current substep.
    subroutine rk_weno_stage(j)
        implicit none
        integer, intent(in) :: j
        integer :: i

        if(j==1) then
            do i = 1, Num_gam_e
                specext(i) = stage(1,i) + CFL*(flux(i)-flux(i-1)) + dF1(i)*dDR
            end do
            stage(2,:) = specext
        else if(j==2) then
            do i = 1, Num_gam_e
                specext(i) = 0.75d0*stage(1,i) + 0.25d0*(stage(2,i) + &
                             CFL*(flux(i)-flux(i-1))) + 0.25d0*dF1(i)*dDR
            end do
            stage(3,:) = specext
        else if(j==3) then
            do i = 1, Num_gam_e
                specext(i) = (stage(1,i) + 2.0d0*(stage(3,i) + &
                             CFL*(flux(i)-flux(i-1))))/3.0d0 + 2.0d0/3.0d0*dF1(i)*dDR
            end do
        end if
    end subroutine rk_weno_stage

! WENO5 正通量重构 f+(f_{i-2},...,f_{i+2})。
! WENO5 positive-flux reconstruction f+(f_{i-2},...,f_{i+2}).
function weno5_positive_flux(fps)
    implicit none
    real(8), intent(in), dimension(-2:2) :: fps
    real(8) :: weno5_positive_flux
    real(8), dimension(3) :: omega,fu,beta
    real(8), dimension(3) :: alpha,fomega
    real(8) :: tao5,totalpha,eps

    omega(1) = 0.1d0;   omega(2) = 0.6d0;   omega(3) = 0.3d0
    eps=1d-30

    fu(1) =  1.0d0/3.0d0*fps(-2) - 7.0d0/6.0d0*fps(-1) + 11.0d0/6.0d0*fps(0)
    fu(2) = -1.0d0/6.0d0*fps(-1) + 5.0d0/6.0d0*fps(0)  + 1.0d0/3.0d0*fps(1)
    fu(3) =  1.0d0/3.0d0*fps(0)  + 5.0d0/6.0d0*fps(1)  - 1.0d0/6.0d0*fps(2)

    beta(1) = 13.0d0/12.0d0*( fps(-2) - 2.0d0*fps(-1) + fps(0) )**2 &
            + 0.25d0*( fps(-2) - 4.0d0*fps(-1) + 3.0d0*fps(0) )**2
    beta(2) = 13.0d0/12.0d0*( fps(-1) - 2.0d0*fps(0) + fps(1) )**2 &
            + 0.25d0*( fps(1) - fps(-1) )**2
    beta(3) = 13.0d0/12.0d0*( fps(0) - 2.0d0*fps(1) + fps(2) )**2 &
            + 0.25d0*( 3.0d0*fps(0) - 4.0d0*fps(1) + fps(2) )**2

    tao5 = abs(beta(1) - beta(3))

    alpha(:) = omega(:)*( 1.0d0 + (tao5/(beta(:)+eps))**2 )
    totalpha = alpha(1) + alpha(2) + alpha(3)

    if(totalpha < eps) then
        weno5_positive_flux = fu(2)
    else
        fomega(:) = alpha(:)/totalpha
        weno5_positive_flux = fu(1)*fomega(1) + fu(2)*fomega(2) + fu(3)*fomega(3)
    end if

end function weno5_positive_flux

! WENO5 负通量重构 f-(f_{i-1},...,f_{i+3})。
! WENO5 negative-flux reconstruction f-(f_{i-1},...,f_{i+3}).
function weno5_negative_flux(fms)
    implicit none
    real(8), intent(in), dimension(-1:3) :: fms
    real(8) :: weno5_negative_flux
    real(8), dimension(3) :: omega,fu,beta
    real(8), dimension(3) :: alpha,fomega
    real(8) :: tao5,totalpha,eps

    omega(1) = 0.1d0;   omega(2) = 0.6d0;   omega(3) = 0.3d0
    eps=1d-30

    fu(1) =  1.0d0/3.0d0*fms(3) - 7.0d0/6.0d0*fms(2) + 11.0d0/6.0d0*fms(1)
    fu(2) = -1.0d0/6.0d0*fms(2) + 5.0d0/6.0d0*fms(1) + 1.0d0/3.0d0*fms(0)
    fu(3) =  1.0d0/3.0d0*fms(1) + 5.0d0/6.0d0*fms(0) - 1.0d0/6.0d0*fms(-1)

    beta(1) = 13.0d0/12.0d0*( fms(3) - 2.0d0*fms(2) + fms(1) )**2 &
            + 0.25d0*( fms(3) - 4.0d0*fms(2) + 3.0d0*fms(1) )**2
    beta(2) = 13.0d0/12.0d0*( fms(2) - 2.0d0*fms(1) + fms(0) )**2 &
            + 0.25d0*( fms(2) - fms(0) )**2
    beta(3) = 13.0d0/12.0d0*( fms(1) - 2.0d0*fms(0) + fms(-1) )**2 &
            + 0.25d0*( 3.0d0*fms(1) - 4.0d0*fms(0) + fms(-1) )**2

    tao5 = abs(beta(1) - beta(3))

    alpha(:) = omega(:)*( 1.0d0 + (tao5/(beta(:)+eps))**2 )
    totalpha = alpha(1) + alpha(2) + alpha(3)

    if(totalpha < eps) then
        weno5_negative_flux = fu(2)
    else
        fomega(:) = alpha(:)/totalpha
        weno5_negative_flux = fu(1)*fomega(1) + fu(2)*fomega(2) + fu(3)*fomega(3)
    end if

end function weno5_negative_flux
end subroutine fs_weno5_1d
