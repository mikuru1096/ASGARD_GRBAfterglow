! Created by rj on 2021/2/4.

!##################################################################################################
!When an intrinsic SED 'F_tot' observed in the comoving frame send in,
!to produce the observed SED 'F_tot_obs' after consider the EATS and Doppler boosting effect.
!##################################################################################################

! Shell-level top-hat projection.
! 顺序: angular EATS geometry -> Doppler/redshift factor -> log-frequency interpolation
!       -> observer-time accumulation over theta/phi cells.
subroutine sed_interpolation(Boundary,R_Tobs1,R_gamma,R,F_tot,V_seed,V_obs,Tobs, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads, F_tot_obs)
    !$ use omp_lib
    use constants
    use interpolation_common
    implicit none
    !##############################################################################################
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_Tobs) :: Tobs
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu_obs) :: V_obs
    real(8), intent(in), dimension(Num_R) :: R_Tobs1,R_gamma,R
    real(8), intent(in), dimension(Num_nu,Num_R) :: F_tot
    real(8), intent(out), dimension(Num_nu_obs,Num_Tobs) :: F_tot_obs


    real(8), allocatable, dimension(:,:) :: F_tot_obs_temp
    real(8), allocatable, dimension(:) :: V_obs_log,V_seed_log
    real(8), dimension(Num_R) :: R_Tobs_theta
    real(8), dimension(Num_nu) :: F_tot_theta,F_tot_log_theta,V_seed_log_theta
    real(8), dimension(Num_Tobs) :: T_sorted
    integer, dimension(Num_Tobs) :: T_order
    real(8) :: z,OpeningAngle_jet,Tv,dPhi,phi_scale,dtheta,Taa_lower,Taa_boundary,Taa_center,domega
    real(8) :: Phi_center,DMu,Ratio,DG,Beta,doppler,lgamlo,lgamhi,ldomega,ldopred
    integer :: I_Theta,i_Phi,Iobs,K2,II,last_k2
    logical :: time_ordered
    allocate (F_tot_obs_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))

    F_tot_obs=0d0
    F_tot_obs_temp=0d0

    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)

    dPhi = pi / Num_Phi
    phi_scale = 2d0
    if (Num_Phi == 1) then
        dPhi = pi / 1440d0
        phi_scale = 2d0 * 1440d0
    end if
    dtheta=OpeningAngle_jet/Num_Theta

    V_obs_log = log(V_obs)
    V_seed_log = log(V_seed)
    call time_order(Tobs,Num_Tobs,T_order,T_sorted,time_ordered)
    !$OMP PARALLEL num_threads(n_threads), reduction(+:F_tot_obs_temp), private(I_Theta, Taa_lower, Taa_boundary, &
    !$OMP& Taa_center, domega, i_Phi, Phi_center, DMu, Iobs, K2, II, Ratio, DG, Beta, doppler, &
    !$OMP& R_Tobs_theta, F_tot_theta, F_tot_log_theta, V_seed_log_theta, &
    !$OMP& last_k2, lgamlo, lgamhi, ldomega, ldopred)
    !$OMP DO SCHEDULE(GUIDED,4)
    do I_Theta=1,Num_Theta
       Taa_lower=dtheta*(I_Theta-1)
       Taa_boundary=dtheta*I_Theta
       Taa_center=dtheta*(I_Theta-0.5)
       domega=(dcos(Taa_lower)-dcos(Taa_boundary))*dPhi
       ldomega=log(domega)-log(4.0d0*pi)
       do i_Phi=1,Num_Phi
          Phi_center=(i_Phi-0.5)*dPhi
          DMu=dcos(Tv)*dcos(Taa_center)+dsin(Tv)*dsin(Taa_center)*dcos(Phi_center)
          R_Tobs_theta=R_Tobs1+R*(1d0-DMu)*(1d0+z)/Para_c
          II=1
          last_k2=0
          do Iobs=1,Num_Tobs
             if (T_sorted(Iobs) < R_Tobs_theta(1) .or. T_sorted(Iobs) >= R_Tobs_theta(Num_R)) cycle
             do while (II < Num_R-1 .and. T_sorted(Iobs) >= R_Tobs_theta(II+1))
                II=II+1
             end do
             K2=II
             if (K2 /= last_k2) then
                 lgamlo=log(R_gamma(K2))
                 lgamhi=log(R_gamma(K2+1))
                 last_k2=K2
             end if
             Ratio=(T_sorted(Iobs)-R_Tobs_theta(K2))/(R_Tobs_theta(K2+1)-R_Tobs_theta(K2))
             DG=exp(lgamlo+Ratio*(lgamhi-lgamlo))
             ! 时间方向用非负通量重构；频率方向仍在 interpolation_common 中按 log SED 插值。
             F_tot_theta=(1d0-Ratio)*F_tot(:,K2)+Ratio*F_tot(:,K2+1)
             F_tot_log_theta=-huge(1d0)
             where(F_tot_theta > 0d0) F_tot_log_theta=log(F_tot_theta)
             Beta=dsqrt(1d0-DG**(-2))

             doppler=DG*(1d0-Beta*DMu) !Doppler factor, changed with R
             ldopred=log(doppler)+log(1d0+z)
             F_tot_log_theta=F_tot_log_theta+ldomega-3d0*log(doppler)
             V_seed_log_theta=V_seed_log-ldopred
             call accum_logsed(V_seed_log_theta,F_tot_log_theta, &
                                                   Num_nu,V_obs_log,Num_nu_obs,F_tot_obs_temp(:,Iobs))
         end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    if (time_ordered) then
        F_tot_obs=F_tot_obs_temp*phi_scale
    else
        do Iobs=1,Num_Tobs
            F_tot_obs(:,T_order(Iobs))=F_tot_obs_temp(:,Iobs)*phi_scale
        end do
    end if

    deallocate (F_tot_obs_temp,V_obs_log,V_seed_log)


    return
end subroutine sed_interpolation

! 壳层级EATS自适应角向积分：每个基础theta cell用一阶/二阶中点规则估计角向积分误差并递归细分。
subroutine sed_adaptive_theta(Boundary,R_Tobs1,R_gamma,R,F_tot,V_seed,V_obs,Tobs, &
                             adaptive_rtol,addepthmax, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads,F_tot_obs)
    !$ use omp_lib
    use constants
    use interpolation_common
    implicit none
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads
    integer, intent(in) :: addepthmax
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_Tobs) :: Tobs
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu_obs) :: V_obs
    real(8), intent(in), dimension(Num_R) :: R_Tobs1,R_gamma,R
    real(8), intent(in), dimension(Num_nu,Num_R) :: F_tot
    real(8), intent(in) :: adaptive_rtol
    real(8), intent(out), dimension(Num_nu_obs,Num_Tobs) :: F_tot_obs
    real(8), allocatable, dimension(:,:) :: F_tot_obs_temp
    real(8), allocatable, dimension(:) :: V_obs_log,V_seed_log
    real(8), dimension(Num_nu_obs,Num_Tobs) :: cell_obs
    real(8), dimension(Num_Tobs) :: T_sorted
    integer, dimension(Num_Tobs) :: T_order
    real(8) :: z,OpeningAngle_jet,Tv,dPhi,phi_scale,dtheta,Taa_lower,Taa_boundary,Phi_center
    integer :: I_Theta,i_Phi,Iobs
    logical :: time_ordered

    if (addepthmax == 0) then
        call sed_interpolation(Boundary,R_Tobs1,R_gamma,R,F_tot,V_seed,V_obs,Tobs, &
                               n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads,F_tot_obs)
        return
    end if

    allocate(F_tot_obs_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))

    F_tot_obs = 0d0
    F_tot_obs_temp = 0d0
    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)
    dPhi = pi/Num_Phi
    phi_scale = 2d0
    if (Num_Phi == 1) then
        dPhi = pi/1440d0
        phi_scale = 2d0*1440d0
    end if
    dtheta = OpeningAngle_jet/Num_Theta
    V_obs_log = dlog(V_obs)
    V_seed_log = dlog(V_seed)
    call time_order(Tobs,Num_Tobs,T_order,T_sorted,time_ordered)

    !$OMP PARALLEL num_threads(n_threads), reduction(+:F_tot_obs_temp), private(I_Theta, &
    !$OMP& i_Phi,Phi_center,Taa_lower,Taa_boundary,cell_obs)
    !$OMP DO SCHEDULE(GUIDED,4)
    do I_Theta = 1, Num_Theta
        Taa_lower = dtheta*(I_Theta-1)
        Taa_boundary = dtheta*I_Theta
        do i_Phi = 1, Num_Phi
            Phi_center = (i_Phi-0.5d0)*dPhi
            cell_obs = 0d0
            call integrate_theta_cell(Taa_lower,Taa_boundary,Phi_center,0,cell_obs)
            F_tot_obs_temp = F_tot_obs_temp + cell_obs
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    if (time_ordered) then
        F_tot_obs=F_tot_obs_temp*phi_scale
    else
        do Iobs=1,Num_Tobs
            F_tot_obs(:,T_order(Iobs))=F_tot_obs_temp(:,Iobs)*phi_scale
        end do
    end if
    deallocate(F_tot_obs_temp,V_obs_log,V_seed_log)
    return

contains

recursive subroutine integrate_theta_cell(theta_lo,theta_hi,phi_center,depth,accum_obs)
    implicit none
    integer, intent(in) :: depth
    real(8), intent(in) :: theta_lo,theta_hi,phi_center
    real(8), intent(inout), dimension(Num_nu_obs,Num_Tobs) :: accum_obs
    real(8) :: theta_mid,theta_lmid,theta_rmid,err_norm,flux_norm
    real(8), dimension(Num_nu_obs,Num_Tobs) :: coarse_obs,fine_obs

    theta_mid = 0.5d0*(theta_lo+theta_hi)
    coarse_obs = 0d0
    call project_theta_sample(theta_lo,theta_hi,theta_mid,phi_center,coarse_obs)
    if (depth >= addepthmax) then
        accum_obs = accum_obs + coarse_obs
        return
    end if
    theta_lmid = 0.5d0*(theta_lo+theta_mid)
    theta_rmid = 0.5d0*(theta_mid+theta_hi)
    fine_obs = 0d0
    call project_theta_sample(theta_lo,theta_mid,theta_lmid,phi_center,fine_obs)
    call project_theta_sample(theta_mid,theta_hi,theta_rmid,phi_center,fine_obs)
    flux_norm = sum(abs(fine_obs))
    err_norm = sum(abs(fine_obs-coarse_obs))
    if (flux_norm == 0d0 .or. err_norm <= adaptive_rtol*flux_norm) then
        accum_obs = accum_obs + fine_obs
    else
        call integrate_theta_cell(theta_lo,theta_mid,phi_center,depth+1,accum_obs)
        call integrate_theta_cell(theta_mid,theta_hi,phi_center,depth+1,accum_obs)
    end if
end subroutine integrate_theta_cell

subroutine project_theta_sample(theta_lo,theta_hi,theta_center,phi_center,local_obs)
    implicit none
    real(8), intent(in) :: theta_lo,theta_hi,theta_center,phi_center
    real(8), intent(inout), dimension(Num_nu_obs,Num_Tobs) :: local_obs
    real(8), dimension(Num_R) :: R_Tobs_theta
    real(8) :: lgamlo,lgamhi,ldomega
    real(8) :: domega,DMu,Ratio
    integer :: Iobs,K2,II
    integer :: last_k2

    domega = (dcos(theta_lo)-dcos(theta_hi))*dPhi
    ldomega = dlog(domega)-dlog(4.0d0*pi)
    DMu = dcos(Tv)*dcos(theta_center)+dsin(Tv)*dsin(theta_center)*dcos(phi_center)
    R_Tobs_theta = R_Tobs1+R*(1d0-DMu)*(1d0+z)/Para_c
    II = 1
    last_k2 = 0
    do Iobs = 1, Num_Tobs
        if (T_sorted(Iobs) < R_Tobs_theta(1) .or. T_sorted(Iobs) >= R_Tobs_theta(Num_R)) cycle
        do while (II < Num_R-1 .and. T_sorted(Iobs) >= R_Tobs_theta(II+1))
            II = II + 1
        end do
        K2 = II
        if (K2 /= last_k2) then
            lgamlo = dlog(R_gamma(K2))
            lgamhi = dlog(R_gamma(K2+1))
            last_k2 = K2
        end if
        Ratio = (T_sorted(Iobs)-R_Tobs_theta(K2))/(R_Tobs_theta(K2+1)-R_Tobs_theta(K2))
        call projshellseg(Iobs,K2,Ratio,DMu,ldomega, &
                                   lgamlo,lgamhi,local_obs)
    end do
end subroutine project_theta_sample

subroutine projshellseg(K1,K2,Ratio,DMu,ldomega,lgamlo,lgamhi,local_obs)
    implicit none
    integer, intent(in) :: K1,K2
    real(8), intent(in) :: Ratio,DMu,ldomega,lgamlo,lgamhi
    real(8), intent(inout), dimension(Num_nu_obs,Num_Tobs) :: local_obs
    real(8), dimension(Num_nu) :: F_tot_theta,F_tot_log_theta,V_seed_log_theta
    real(8) :: DG,Beta,doppler,ldopred

    DG = dexp(lgamlo+Ratio*(lgamhi-lgamlo))
    F_tot_theta = (1d0-Ratio)*F_tot(:,K2)+Ratio*F_tot(:,K2+1)
    F_tot_log_theta = -huge(1d0)
    where(F_tot_theta > 0d0) F_tot_log_theta = dlog(F_tot_theta)
    Beta = dsqrt(1d0-DG**(-2))
    doppler = DG*(1d0-Beta*DMu)
    ldopred = dlog(doppler)+dlog(1d0+z)
    F_tot_log_theta = F_tot_log_theta+ldomega-3d0*dlog(doppler)
    V_seed_log_theta = V_seed_log-ldopred
    call accum_logsed(V_seed_log_theta,F_tot_log_theta, &
                                          Num_nu,V_obs_log,Num_nu_obs,local_obs(:,K1))
end subroutine projshellseg
end subroutine sed_adaptive_theta

! 将χ分辨有限厚壳层积分到完整top-hat角网格：θ/φ EATS + χ厚度EATS + 局域Doppler。
subroutine sed_interpolation_chi(Boundary,R_Tobs1,R_front,F_chi,Tau_chi,R_chi,Gamma_chi,Chi_weight,V_seed,V_obs,Tobs, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_chi,Num_R,n_threads,F_tot_obs)
    !$ use omp_lib
    use constants
    use interpolation_common
    implicit none
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_chi,Num_R,n_threads
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs1,R_front
    real(8), intent(in), dimension(Num_Tobs) :: Tobs
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu_obs) :: V_obs
    real(8), intent(in), dimension(Num_nu,Num_chi,Num_R) :: F_chi,Tau_chi
    real(8), intent(in), dimension(Num_chi,Num_R) :: R_chi,Gamma_chi,Chi_weight
    real(8), intent(out), dimension(Num_nu_obs,Num_Tobs) :: F_tot_obs
    real(8), allocatable, dimension(:,:) :: F_temp
    real(8), allocatable, dimension(:) :: V_obs_log,V_seed_log
    real(8), allocatable, dimension(:,:,:) :: Tau_prefix
    real(8), allocatable, dimension(:) :: Angle_dmu,Angle_logw
    real(8), dimension(Num_R) :: R_Tobs_chi
    real(8), dimension(Num_Tobs) :: T_sorted
    integer, dimension(Num_Tobs) :: T_order
    real(8) :: ldomega,z,OpeningAngle_jet,Tv,dPhi,phi_scale,dtheta
    real(8) :: Taa_lower,Taa_boundary,Taa_center,domega,Phi_center,DMu,Ratio
    real(8) :: lgamlo,lgamhi,segment_lo,segment_hi,cos_tv,sin_tv
    integer :: I_chi,I_Theta,i_Phi,Iobs,K2,II,I_ang,last_k2,k_start,lowerreal8
    logical :: monotonic_chi,time_ordered

    allocate(F_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))
    allocate(Tau_prefix(Num_nu,0:Num_chi,Num_R))
    allocate(Angle_dmu(Num_Theta*Num_Phi),Angle_logw(Num_Theta*Num_Phi))
    F_tot_obs = 0d0
    F_temp = 0d0
    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)
    dPhi = pi/Num_Phi
    phi_scale = 2d0
    if (Num_Phi == 1) then
        dPhi = pi/1440d0
        phi_scale = 2d0*1440d0
    end if
    dtheta = OpeningAngle_jet/Num_Theta
    V_obs_log = dlog(V_obs)
    V_seed_log = dlog(V_seed)
    call time_order(Tobs,Num_Tobs,T_order,T_sorted,time_ordered)
    Tau_prefix(:,0,:) = 0d0
    do I_chi = 1, Num_chi
        Tau_prefix(:,I_chi,:) = Tau_prefix(:,I_chi-1,:) + Tau_chi(:,I_chi,:)
    end do
    cos_tv = dcos(Tv)
    sin_tv = dsin(Tv)
    do I_ang = 1, Num_Theta*Num_Phi
        I_Theta = (I_ang-1)/Num_Phi + 1
        i_Phi = I_ang - (I_Theta-1)*Num_Phi
        Taa_lower = dtheta*(I_Theta-1)
        Taa_boundary = dtheta*I_Theta
        Taa_center = dtheta*(I_Theta-0.5d0)
        Phi_center = (i_Phi-0.5d0)*dPhi
        domega = (dcos(Taa_lower)-dcos(Taa_boundary))*dPhi
        Angle_logw(I_ang) = dlog(domega)-dlog(4d0*pi)
        Angle_dmu(I_ang) = cos_tv*dcos(Taa_center)+sin_tv*dsin(Taa_center)*dcos(Phi_center)
    end do

    !$OMP PARALLEL num_threads(n_threads), reduction(+:F_temp), private(I_ang,ldomega,DMu, &
    !$OMP& I_chi,R_Tobs_chi,monotonic_chi,last_k2,Iobs,K2,II,Ratio, &
    !$OMP& lgamlo,lgamhi,segment_lo,segment_hi,k_start)
    !$OMP DO SCHEDULE(STATIC)
    do I_ang = 1, Num_Theta*Num_Phi
        ldomega = Angle_logw(I_ang)
        DMu = Angle_dmu(I_ang)
        do I_chi = 1, Num_chi
                R_Tobs_chi = R_Tobs1 + (1d0+z)*(R_front - R_chi(I_chi,:)*DMu)/Para_c
                monotonic_chi = all(R_Tobs_chi(2:Num_R) > R_Tobs_chi(1:Num_R-1))
                if (monotonic_chi) then
                    II = 1
                    last_k2 = 0
                    do Iobs = 1, Num_Tobs
                        if (T_sorted(Iobs) < R_Tobs_chi(1) .or. T_sorted(Iobs) >= R_Tobs_chi(Num_R)) cycle
                        do while (II < Num_R-1 .and. T_sorted(Iobs) >= R_Tobs_chi(II+1))
                            II = II + 1
                        end do
                        K2 = II
                        if (K2 /= last_k2) then
                            lgamlo = dlog(Gamma_chi(I_chi,K2))
                            lgamhi = dlog(Gamma_chi(I_chi,K2+1))
                            last_k2 = K2
                        end if
                        Ratio = (T_sorted(Iobs)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                        call project_chi(I_chi,K2,Ratio,DMu,ldomega, &
                                                      lgamlo,lgamhi,F_temp(:,Iobs))
                    end do
                else
                    do K2 = 1, Num_R-1
                        segment_lo = min(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                        segment_hi = max(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                        k_start = lowerreal8(T_sorted,Num_Tobs,segment_lo)
                        lgamlo = dlog(Gamma_chi(I_chi,K2))
                        lgamhi = dlog(Gamma_chi(I_chi,K2+1))
                        do Iobs = k_start, Num_Tobs
                            if (T_sorted(Iobs) >= segment_hi) exit
                            Ratio = (T_sorted(Iobs)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                            call project_chi(I_chi,K2,Ratio,DMu,ldomega, &
                                                          lgamlo,lgamhi,F_temp(:,Iobs))
                        end do
                    end do
                end if
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    if (time_ordered) then
        F_tot_obs=F_temp*phi_scale
    else
        do Iobs=1,Num_Tobs
            F_tot_obs(:,T_order(Iobs))=F_temp(:,Iobs)*phi_scale
        end do
    end if
    deallocate(F_temp,V_obs_log,V_seed_log,Tau_prefix,Angle_dmu,Angle_logw)
    return

contains

subroutine project_chi(I_chi,K2,Ratio,DMu,ldomega,lgamlo,lgamhi,flux_accum)
    implicit none
    integer, intent(in) :: I_chi,K2
    real(8), intent(in) :: Ratio,DMu,ldomega,lgamlo,lgamhi
    real(8), intent(inout), dimension(Num_nu_obs) :: flux_accum
    real(8), dimension(Num_nu) :: F_theta
    real(8) :: ldopred,lfluxw,doppler

    call chi_state(I_chi,K2,Ratio,DMu,lgamlo,lgamhi, &
                                   ldopred,doppler,F_theta)
    lfluxw = ldomega - 3d0*dlog(doppler)
    call accum_shifted(V_seed_log,F_theta,Num_nu,V_obs_log,Num_nu_obs, &
                                                     ldopred,lfluxw,flux_accum)
end subroutine project_chi

subroutine chi_state(I_chi,K2,Ratio,DMu,lgamlo,lgamhi, &
                                     ldopred,doppler,F_theta)
    implicit none
    integer, intent(in) :: I_chi,K2
    real(8), intent(in) :: Ratio,DMu,lgamlo,lgamhi
    real(8), intent(out), dimension(Num_nu) :: F_theta
    real(8), intent(out) :: ldopred,doppler
    real(8) :: DG,Beta
    real(8), external :: chi_escape
    integer :: I_nu_inl

    DG = dexp(lgamlo+Ratio*(lgamhi-lgamlo))
    Beta = dsqrt(1d0-DG**(-2))
    doppler = DG*(1d0-Beta*DMu)
    ldopred = dlog(doppler)+dlog(1d0+z)
    ! Inlined from accumulate_chi_cell_source
    F_theta = 0d0
    do I_nu_inl = 1, Num_nu
        F_theta(I_nu_inl) = (1d0-Ratio) * F_chi(I_nu_inl,I_chi,K2)*Chi_weight(I_chi,K2) &
            * chi_escape(Tau_prefix(I_nu_inl,I_chi-1,K2), Tau_chi(I_nu_inl,I_chi,K2)) &
          + Ratio * F_chi(I_nu_inl,I_chi,K2+1)*Chi_weight(I_chi,K2+1) &
            * chi_escape(Tau_prefix(I_nu_inl,I_chi-1,K2+1), Tau_chi(I_nu_inl,I_chi,K2+1))
    end do
end subroutine chi_state
end subroutine sed_interpolation_chi

! Direct-electron chi projection: build chi-local synchrotron spectra on V_seed, then use top-hat chi EATS.
subroutine sed_chi_electron(Boundary,R_Tobs1,R_front,DNe_chi,B_chi,R_chi,Gamma_chi,Chi_weight, &
                             gam_e,V_seed,V_obs,Tobs,n,Num_gam_e,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi, &
                             Num_chi,Num_R,n_threads,F_tot_obs)
    use constants
    use rad_common, only: syn_flux_chi
    implicit none
    integer, intent(in) :: n,Num_gam_e,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_chi,Num_R,n_threads
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs1,R_front
    real(8), intent(in), dimension(Num_Tobs) :: Tobs
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu_obs) :: V_obs
    real(8), intent(in), dimension(Num_gam_e) :: gam_e
    real(8), intent(in), dimension(Num_gam_e,Num_chi,Num_R) :: DNe_chi
    real(8), intent(in), dimension(Num_chi,Num_R) :: B_chi
    real(8), intent(in), dimension(Num_chi,Num_R) :: R_chi,Gamma_chi,Chi_weight
    real(8), intent(out), dimension(Num_nu_obs,Num_Tobs) :: F_tot_obs
    real(8), allocatable, dimension(:,:,:) :: F_chi,Tau_chi
    real(8), allocatable, dimension(:,:) :: P_syn,Tau_syn
    real(8) :: DL,z,flux_prefactor
    integer :: K2

    DL = Boundary(13)
    z = Boundary(8)
    flux_prefactor = (1d0+z)/(4d0*pi*DL*DL)
    allocate(F_chi(Num_nu,Num_chi,Num_R),Tau_chi(Num_nu,Num_chi,Num_R), &
             P_syn(Num_nu,Num_chi),Tau_syn(Num_nu,Num_chi))
    do K2 = 1, Num_R
        call syn_flux_chi(R_front(K2),Num_gam_e,Num_nu,Num_chi,gam_e,DNe_chi(:,:,K2), &
                                                   V_seed,B_chi(:,K2),Chi_weight(:,K2),1.046d4, &
                                                   P_syn,Tau_syn)
        F_chi(:,:,K2) = P_syn*flux_prefactor
        Tau_chi(:,:,K2) = Tau_syn
    end do
    call sed_interpolation_chi(Boundary,R_Tobs1,R_front,F_chi,Tau_chi,R_chi,Gamma_chi,Chi_weight,V_seed,V_obs,Tobs, &
                               n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_chi,Num_R,n_threads,F_tot_obs)
    deallocate(F_chi,Tau_chi,P_syn,Tau_syn)
end subroutine sed_chi_electron

! Structured precomputed-ring chi projection: Python supplies 1d0 theta-ring chi spectra and SSA grids.
subroutine sed_chi_ring(Boundary,R_Tobs1,R_front, &
                             F_ring,Tau_ring,R_chi,Gamma_chi,Chi_weight,V_seed,V_obs,Tobs,theta_lo,theta_hi, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R,F_tot_obs)
    implicit none
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Tobs1,R_front
    real(8), intent(in), dimension(Num_Tobs) :: Tobs
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu_obs) :: V_obs
    real(8), intent(in), dimension(Num_nu,Num_chi,Num_R) :: F_ring,Tau_ring
    real(8), intent(in), dimension(Num_chi,Num_R) :: R_chi,Gamma_chi,Chi_weight
    real(8), intent(in) :: theta_lo,theta_hi
    real(8), intent(out), dimension(Num_nu_obs,Num_Tobs) :: F_tot_obs

    call sed_chiring_core(Boundary,R_Tobs1,R_front,F_ring,Tau_ring,R_chi,Gamma_chi, &
                          Chi_weight,V_seed,V_obs,Tobs,theta_lo,theta_hi,1d0, &
                          n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R,F_tot_obs)
end subroutine sed_chi_ring

subroutine sed_chiring_batchlum(Boundary,R_Tobs1,R_front, &
                             F_ring,Tau_ring,R_chi,Gamma_chi,Chi_weight,V_seed,V_obs,Tobs,theta_lo,theta_hi, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R,Num_ring,F_tot_obs)
    use constants
    implicit none
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R,Num_ring
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_R,Num_ring) :: R_Tobs1,R_front
    real(8), intent(in), dimension(Num_Tobs) :: Tobs
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu_obs) :: V_obs
    real(8), intent(in) :: F_ring(Num_nu,Num_chi,Num_R,Num_ring),Tau_ring(Num_nu,Num_chi,Num_R,Num_ring)
    real(8), intent(in) :: R_chi(Num_chi,Num_R,Num_ring),Gamma_chi(Num_chi,Num_R,Num_ring)
    real(8), intent(in) :: Chi_weight(Num_chi,Num_R,Num_ring)
    real(8), intent(in), dimension(Num_ring) :: theta_lo,theta_hi
    real(8), intent(out), dimension(Num_nu_obs,Num_Tobs) :: F_tot_obs
    real(8), allocatable, dimension(:,:) :: Fring
    real(8) :: z,DL,fluxscale
    integer :: Iring

    z = Boundary(8)
    DL = Boundary(13)
    fluxscale = (1d0+z)/(4d0*pi*DL*DL)
    allocate(Fring(Num_nu_obs,Num_Tobs))
    F_tot_obs = 0d0
    do Iring = 1, Num_ring
        call sed_chiring_core(Boundary,R_Tobs1(:,Iring),R_front(:,Iring), &
                              F_ring(:,:,:,Iring),Tau_ring(:,:,:,Iring), &
                              R_chi(:,:,Iring),Gamma_chi(:,:,Iring),Chi_weight(:,:,Iring), &
                              V_seed,V_obs,Tobs,theta_lo(Iring),theta_hi(Iring),fluxscale, &
                              n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R,Fring)
        F_tot_obs = F_tot_obs + Fring
    end do
    deallocate(Fring)
end subroutine sed_chiring_batchlum

subroutine sed_chiring_batchlum_ray(Boundary,R_Tobs1,R_front, &
                             F_ring,Tau_ring,R_chi,Gamma_chi,Chi_weight,V_seed,V_obs,Tobs,theta_lo,theta_hi, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R,Num_ring,F_tot_obs)
    ! 中文：每个发射 patch 自身作为光线起点；前景 SSA 光深由同一观测时刻的投影遮挡和共动系弦长积分给出。
    ! English: Each emitting patch is the ray origin; foreground SSA optical depth comes from
    ! same-observer-time projected occultation and comoving chord integration.
    use constants
    implicit none
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R,Num_ring
    real(8), intent(in) :: Boundary(n), R_Tobs1(Num_R,Num_ring), R_front(Num_R,Num_ring), Tobs(Num_Tobs)
    real(8), intent(in) :: V_seed(Num_nu), V_obs(Num_nu_obs), theta_lo(Num_ring), theta_hi(Num_ring)
    real(8), intent(in) :: F_ring(Num_nu,Num_chi,Num_R,Num_ring), Tau_ring(Num_nu,Num_chi,Num_R,Num_ring)
    real(8), intent(in) :: R_chi(Num_chi,Num_R,Num_ring), Gamma_chi(Num_chi,Num_R,Num_ring), Chi_weight(Num_chi,Num_R,Num_ring)
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    real(8), allocatable, dimension(:) :: evsrc,evtau,evtau0,evdepth,evratio,evldop,evpath,evlogw,evamp,evx,evy
    real(8), allocatable, dimension(:) :: evax,evay,evbx,evby,evdetinv,evxmin,evxmax,evymin,evymax,evrpos,evrin,evrout,evgam,evdrcom
    real(8), allocatable, dimension(:) :: evdmu,evnth,evnph,evhalfth,evhalfph,hitpath,evnr,evntp,evnpp,evrlo,evrhi,evslo,evshi
    real(8), allocatable, dimension(:) :: Vobslog,Vseedlog,leaflo,leafhi,leafc,leafwt,phiwleaf,cphleaf,sphleaf,cthring,sthring
    real(8), allocatable :: logwgrid(:,:),dmugrid(:,:),nthgrid(:,:),nphgrid(:,:)
    integer, allocatable, dimension(:) :: eviring,evichi,evk2,hitfg,hitstart
    real(8) :: Tbase(Num_R,Num_ring), z,DL,fluxscale,Tv,costv,sintv,rootphi,phiw,phi,logzp1,delayfac,dvseed
    real(8) :: theta,cth,sth,cph,sph,dmu,domega,logw,ratio,ntheta,nphi,logg,dg,beta,dop,target,src,tau,rpos,xsky,ysky,ta,tb
    integer :: maxev,nleaf,nev,nhit,hitcap,iring,ichi,inu,kobs,k2,ii,ileaf,ihit
    logical, allocatable :: tmono(:,:,:)

    z = Boundary(8)
    DL = Boundary(13)
    Tv = Boundary(10)
    costv = dcos(Tv)
    sintv = dsin(Tv)
    logzp1 = dlog(1d0+z)
    delayfac = (1d0+z)/Para_c
    fluxscale = (1d0+z)/(4d0*pi*DL*DL)
    rootphi = 2d0*pi/dble(Num_phi_patch)
    allocate(leaflo(Num_phi_patch),leafhi(Num_phi_patch),leafc(Num_phi_patch),leafwt(Num_phi_patch))
    allocate(Vobslog(Num_nu_obs),Vseedlog(Num_nu))
    Vobslog = dlog(V_obs)
    Vseedlog = dlog(V_seed)
    dvseed = Vseedlog(2) - Vseedlog(1)
    F_tot_obs = 0d0

    call init_leaves()
    allocate(phiwleaf(nleaf),cphleaf(nleaf),sphleaf(nleaf),cthring(Num_ring),sthring(Num_ring))
    allocate(logwgrid(Num_ring,nleaf),dmugrid(Num_ring,nleaf),nthgrid(Num_ring,nleaf),nphgrid(Num_ring,nleaf))
    allocate(tmono(Num_chi,Num_ring,nleaf))
    call init_geom()
    call init_tbase()
    call init_tmono()
    maxev = count(tmono) + (Num_R-1)*count(.not. tmono)
    allocate(evsrc(maxev),evtau(maxev),evtau0(maxev),evdepth(maxev))
    allocate(evratio(maxev),evldop(maxev),evpath(maxev),evlogw(maxev),evamp(maxev))
    allocate(evx(maxev),evy(maxev),evax(maxev),evay(maxev),evbx(maxev),evby(maxev),evdetinv(maxev))
    allocate(evxmin(maxev),evxmax(maxev),evymin(maxev),evymax(maxev))
    allocate(evrpos(maxev),evrin(maxev),evrout(maxev),evgam(maxev),evdrcom(maxev),evdmu(maxev),evnth(maxev),evnph(maxev))
    allocate(evnr(maxev),evntp(maxev),evnpp(maxev),evrlo(maxev),evrhi(maxev),evslo(maxev),evshi(maxev))
    hitcap = max(1,maxev)
    allocate(evhalfth(maxev),evhalfph(maxev),hitpath(hitcap))
    allocate(eviring(maxev),evichi(maxev),evk2(maxev),hitfg(hitcap),hitstart(maxev+1))
    do kobs = 1, Num_Tobs
        call leaves_flux(F_tot_obs(:,kobs))
    end do
    deallocate(evsrc,evtau,evtau0,evdepth,Vobslog,Vseedlog)
    deallocate(evratio,evldop,evpath,evlogw,evamp,evx,evy,evax,evay,evbx,evby,evdetinv)
    deallocate(evxmin,evxmax,evymin,evymax,eviring,evichi,evk2,hitfg,hitstart)
    deallocate(evrpos,evrin,evrout,evgam,evdrcom,evdmu,evnth,evnph,evhalfth,evhalfph,hitpath)
    deallocate(evnr,evntp,evnpp,evrlo,evrhi,evslo,evshi)
    deallocate(leaflo,leafhi,leafc,leafwt)
    deallocate(phiwleaf,cphleaf,sphleaf,cthring,sthring,logwgrid,dmugrid,nthgrid,nphgrid)
    deallocate(tmono)
    return

contains

subroutine init_leaves()
    implicit none
    integer :: i

    if (mod(Num_phi_patch,2) == 0) then
        nleaf = Num_phi_patch/2
        do i = 1, nleaf
            leaflo(i) = (dble(i)-1d0)*rootphi
            leafhi(i) = dble(i)*rootphi
            leafc(i) = (dble(i)-0.5d0)*rootphi
            leafwt(i) = 2d0
        end do
    else
        nleaf = Num_phi_patch
        do i = 1, nleaf
            leaflo(i) = (dble(i)-1d0)*rootphi
            leafhi(i) = dble(i)*rootphi
            leafc(i) = (dble(i)-0.5d0)*rootphi
            leafwt(i) = 1d0
        end do
    end if
end subroutine init_leaves

subroutine init_geom()
    implicit none
    integer :: jleaf,jring

    do jring = 1, Num_ring
        theta = 0.5d0*(theta_lo(jring)+theta_hi(jring))
        cthring(jring) = dcos(theta)
        sthring(jring) = dsin(theta)
    end do
    do jleaf = 1, nleaf
        phiwleaf(jleaf) = leafhi(jleaf)-leaflo(jleaf)
        cphleaf(jleaf) = dcos(leafc(jleaf))
        sphleaf(jleaf) = dsin(leafc(jleaf))
        do jring = 1, Num_ring
            domega = (dcos(theta_lo(jring))-dcos(theta_hi(jring)))*phiwleaf(jleaf)*leafwt(jleaf)
            logwgrid(jring,jleaf) = dlog(domega)-dlog(4d0*pi)
            dmugrid(jring,jleaf) = costv*cthring(jring) + sintv*sthring(jring)*cphleaf(jleaf)
            nthgrid(jring,jleaf) = sintv*cthring(jring)*cphleaf(jleaf) - costv*sthring(jring)
            nphgrid(jring,jleaf) = -sintv*sphleaf(jleaf)
        end do
    end do
end subroutine init_geom

subroutine init_tbase()
    implicit none

    do iring = 1, Num_ring
        Tbase(:,iring) = R_Tobs1(:,iring) + delayfac*R_front(:,iring)
    end do
end subroutine init_tbase

subroutine init_tmono()
    implicit none
    integer :: jleaf,jring,jrad
    real(8) :: mudelay,tprev(Num_chi),tnext(Num_chi)

    tmono = .true.
    do jleaf = 1, nleaf
        do jring = 1, Num_ring
            mudelay = delayfac*dmugrid(jring,jleaf)
            tprev = Tbase(1,jring) - mudelay*R_chi(:,1,jring)
            do jrad = 2, Num_R
                tnext = Tbase(jrad,jring) - mudelay*R_chi(:,jrad,jring)
                tmono(:,jring,jleaf) = tmono(:,jring,jleaf) .and. tnext > tprev
                tprev = tnext
            end do
        end do
    end do
end subroutine init_tmono

subroutine leaves_flux(fluxvec)
    implicit none
    real(8), intent(out), dimension(Num_nu_obs) :: fluxvec

    nev = 0
    do ileaf = 1, nleaf
        call sample_phi(ileaf)
    end do
    if (Num_chi > 1) call build_hits()
    do inu = 1, Num_nu_obs
        call fill_radiation()
        fluxvec(inu) = ray_flux()
    end do
end subroutine leaves_flux

subroutine sample_phi(leafid)
    implicit none
    integer, intent(in) :: leafid
    integer :: mid
    real(8) :: mudelay,tmid,tview

    tview = Tobs(kobs)
    phiw = phiwleaf(leafid)
    phi = leafc(leafid)
    cph = cphleaf(leafid)
    sph = sphleaf(leafid)
    do iring = 1, Num_ring
        cth = cthring(iring)
        sth = sthring(iring)
        logw = logwgrid(iring,leafid)
        dmu = dmugrid(iring,leafid)
        mudelay = delayfac*dmu
        ntheta = nthgrid(iring,leafid)
        nphi = nphgrid(iring,leafid)
        do ichi = 1, Num_chi
            if (tmono(ichi,iring,leafid)) then
                ta = Tbase(1,iring) - mudelay*R_chi(ichi,1,iring)
                tb = Tbase(Num_R,iring) - mudelay*R_chi(ichi,Num_R,iring)
                if (.not.(tview >= ta .and. tview < tb)) cycle
                ii = 1
                k2 = Num_R
                do while (ii+1 < k2)
                    mid = (ii+k2)/2
                    tmid = Tbase(mid,iring) - mudelay*R_chi(ichi,mid,iring)
                    if (tview >= tmid) then
                        ii = mid
                        ta = tmid
                    else
                        k2 = mid
                        tb = tmid
                    end if
                end do
                k2 = ii
                ratio = (tview-ta)/(tb-ta)
                call sample_event(iring,ichi,k2,ratio)
            else
                ta = Tbase(1,iring) - mudelay*R_chi(ichi,1,iring)
                do k2 = 1, Num_R-1
                    tb = Tbase(k2+1,iring) - mudelay*R_chi(ichi,k2+1,iring)
                    if (tview >= min(ta,tb) .and. tview < max(ta,tb)) then
                        ratio = (tview-ta)/(tb-ta)
                        call sample_event(iring,ichi,k2,ratio)
                    end if
                    ta = tb
                end do
            end if
        end do
    end do
end subroutine sample_phi

subroutine sample_event(iring,ichi,k2,ratio)
    implicit none
    integer, intent(in) :: iring,ichi,k2
    real(8), intent(in) :: ratio
    real(8) :: routv,rinv

    logg = dlog(Gamma_chi(ichi,k2,iring)) + &
           ratio*dlog(Gamma_chi(ichi,k2+1,iring)/Gamma_chi(ichi,k2,iring))
    dg = dexp(logg)
    beta = dsqrt(1d0-dg**(-2))
    dop = dg*(1d0-beta*dmu)
    rpos = (1d0-ratio)*R_chi(ichi,k2,iring) + ratio*R_chi(ichi,k2+1,iring)
    routv = (1d0-ratio)*outer_face(iring,ichi,k2) + ratio*outer_face(iring,ichi,k2+1)
    rinv = (1d0-ratio)*inner_face(iring,ichi,k2) + ratio*inner_face(iring,ichi,k2+1)
    xsky = rpos*(sth*cph*costv - cth*sintv)
    ysky = rpos*sth*sph
    nev = nev + 1
    eviring(nev) = iring
    evichi(nev) = ichi
    evk2(nev) = k2
    evratio(nev) = ratio
    evldop(nev) = dlog(dop)
    evlogw(nev) = logw
    evamp(nev) = fluxscale*dexp(logw-3d0*evldop(nev))
    evdepth(nev) = rpos*dmu
    evx(nev) = xsky
    evy(nev) = ysky
    evrpos(nev) = rpos
    evrin(nev) = rinv
    evrout(nev) = routv
    evgam(nev) = dg
    evdrcom(nev) = dg*(routv-rinv)
    evdmu(nev) = dmu
    evnth(nev) = ntheta
    evnph(nev) = nphi
    ! 中文：观测方向先做相对局部流体的光行差变换，再用共动 slab 弦长计算 SSA 路径。
    ! English: The observer direction is aberrated into the local comoving frame before the SSA path
    ! is measured as a chord through the comoving slab.
    evnr(nev) = (dmu-beta)/(dop/dg)
    evntp(nev) = ntheta/dop
    evnpp(nev) = nphi/dop
    evrlo(nev) = dg*(rinv-rpos)
    evrhi(nev) = dg*(routv-rpos)
    call radial_range(nev)
    evhalfth(nev) = 0.5d0*rpos*(theta_hi(iring)-theta_lo(iring))
    evhalfph(nev) = 0.5d0*rpos*sth*phiw
    evpath(nev) = pathfactor(nev,0d0,0d0)
    evax(nev) = 0.5d0*rpos*(theta_hi(iring)-theta_lo(iring))*(cth*cph*costv + sth*sintv)
    evay(nev) = 0.5d0*rpos*(theta_hi(iring)-theta_lo(iring))*cth*sph
    evbx(nev) = -0.5d0*rpos*phiw*sth*sph*costv
    evby(nev) = 0.5d0*rpos*phiw*sth*cph
    evdetinv(nev) = 1d0/(evax(nev)*evby(nev) - evay(nev)*evbx(nev))
    evxmin(nev) = min(xsky+evax(nev)+evbx(nev),xsky+evax(nev)-evbx(nev), &
                      xsky-evax(nev)+evbx(nev),xsky-evax(nev)-evbx(nev))
    evxmax(nev) = max(xsky+evax(nev)+evbx(nev),xsky+evax(nev)-evbx(nev), &
                      xsky-evax(nev)+evbx(nev),xsky-evax(nev)-evbx(nev))
    evymin(nev) = min(ysky+evay(nev)+evby(nev),ysky+evay(nev)-evby(nev), &
                      ysky-evay(nev)+evby(nev),ysky-evay(nev)-evby(nev))
    evymax(nev) = max(ysky+evay(nev)+evby(nev),ysky+evay(nev)-evby(nev), &
                      ysky-evay(nev)+evby(nev),ysky-evay(nev)-evby(nev))
end subroutine sample_event

subroutine build_hits()
    implicit none
    integer, allocatable, dimension(:) :: bincount,binstart,bincur,binitem
    real(8) :: gridx0,gridy0,dxbin,dybin
    integer :: isrc,jfg,nbin,ncell,nitem,ix,iy,ixlo,ixhi,iylo,iyhi,icell,ipos
    real(8) :: path

    nhit = 0
    if (nev == 0) return
    nbin = max(1,int(dsqrt(dble(nev))))
    ncell = nbin*nbin
    allocate(bincount(ncell),binstart(ncell+1),bincur(ncell))
    bincount = 0
    gridx0 = minval(evxmin(1:nev))
    gridy0 = minval(evymin(1:nev))
    dxbin = (maxval(evxmax(1:nev))-gridx0)/dble(nbin)
    dybin = (maxval(evymax(1:nev))-gridy0)/dble(nbin)
    do jfg = 1, nev
        ixlo = max(1,min(nbin,1+int((evxmin(jfg)-gridx0)/dxbin)))
        ixhi = max(1,min(nbin,1+int((evxmax(jfg)-gridx0)/dxbin)))
        iylo = max(1,min(nbin,1+int((evymin(jfg)-gridy0)/dybin)))
        iyhi = max(1,min(nbin,1+int((evymax(jfg)-gridy0)/dybin)))
        do iy = iylo, iyhi
            do ix = ixlo, ixhi
                bincount(ix+(iy-1)*nbin) = bincount(ix+(iy-1)*nbin) + 1
            end do
        end do
    end do
    binstart(1) = 1
    do icell = 1, ncell
        binstart(icell+1) = binstart(icell) + bincount(icell)
    end do
    nitem = binstart(ncell+1) - 1
    allocate(binitem(nitem))
    bincur = binstart(1:ncell)
    do jfg = 1, nev
        ixlo = max(1,min(nbin,1+int((evxmin(jfg)-gridx0)/dxbin)))
        ixhi = max(1,min(nbin,1+int((evxmax(jfg)-gridx0)/dxbin)))
        iylo = max(1,min(nbin,1+int((evymin(jfg)-gridy0)/dybin)))
        iyhi = max(1,min(nbin,1+int((evymax(jfg)-gridy0)/dybin)))
        do iy = iylo, iyhi
            do ix = ixlo, ixhi
                icell = ix + (iy-1)*nbin
                binitem(bincur(icell)) = jfg
                bincur(icell) = bincur(icell) + 1
            end do
        end do
    end do
    do isrc = 1, nev
        hitstart(isrc) = nhit + 1
        ix = max(1,min(nbin,1+int((evx(isrc)-gridx0)/dxbin)))
        iy = max(1,min(nbin,1+int((evy(isrc)-gridy0)/dybin)))
        icell = ix + (iy-1)*nbin
        do ipos = binstart(icell), binstart(icell+1)-1
            jfg = binitem(ipos)
            if (evdepth(jfg) <= evdepth(isrc)) cycle
            if (hitpatch(isrc,jfg,path)) call append_hit(jfg,path)
        end do
    end do
    hitstart(nev+1) = nhit + 1
    deallocate(bincount,binstart,bincur,binitem)
end subroutine build_hits

subroutine append_hit(jfg,path)
    implicit none
    integer, intent(in) :: jfg
    real(8), intent(in) :: path
    integer :: newcap
    integer, allocatable, dimension(:) :: fgnew
    real(8), allocatable, dimension(:) :: pathnew

    if (nhit == hitcap) then
        newcap = 2*hitcap
        allocate(fgnew(newcap))
        allocate(pathnew(newcap))
        fgnew(1:nhit) = hitfg(1:nhit)
        pathnew(1:nhit) = hitpath(1:nhit)
        call move_alloc(fgnew,hitfg)
        call move_alloc(pathnew,hitpath)
        hitcap = newcap
    end if
    nhit = nhit + 1
    hitfg(nhit) = jfg
    hitpath(nhit) = path
end subroutine append_hit

subroutine fill_radiation()
    implicit none
    integer :: iev,jring,jchi,jk,jseed
    real(8) :: rfac,frac,taulo,tauhi,srclo,srchi

    do iev = 1, nev
        jring = eviring(iev)
        jchi = evichi(iev)
        jk = evk2(iev)
        rfac = evratio(iev)
        target = Vobslog(inu) + evldop(iev) + logzp1
        if (target <= Vseedlog(1) .or. target > Vseedlog(Num_nu)) then
            evsrc(iev) = 0d0
            evtau(iev) = 0d0
            evtau0(iev) = 0d0
            cycle
        end if
        jseed = min(Num_nu-1,1+int((target-Vseedlog(1))/dvseed))
        frac = (target - Vseedlog(jseed))/(Vseedlog(jseed+1)-Vseedlog(jseed))
        if (Num_chi == 1) then
            srclo = interpescidx(jring,jchi,jk,jseed,frac,evpath(iev)) * Chi_weight(jchi,jk,jring)
            srchi = interpescidx(jring,jchi,jk+1,jseed,frac,evpath(iev)) * Chi_weight(jchi,jk+1,jring)
            src = ((1d0-rfac)*srclo + rfac*srchi)*evamp(iev)
            tau = 0d0
        else
            taulo = interplogidx(Tau_ring(:,jchi,jk,jring),jseed,frac)
            tauhi = interplogidx(Tau_ring(:,jchi,jk+1,jring),jseed,frac)
            srclo = interplogidx(F_ring(:,jchi,jk,jring),jseed,frac) * Chi_weight(jchi,jk,jring)
            srchi = interplogidx(F_ring(:,jchi,jk+1,jring),jseed,frac) * Chi_weight(jchi,jk+1,jring)
            tau = (1d0-rfac)*taulo + rfac*tauhi
            src = ((1d0-rfac)*srclo + rfac*srchi)*evamp(iev)
        end if
        if (src <= 0d0) src = 0d0
        evsrc(iev) = src
        evtau0(iev) = tau
        evtau(iev) = tau*evpath(iev)
    end do
end subroutine fill_radiation

real(8) function interpescidx(jring,jchi,jk,jseed,frac,path)
    implicit none
    integer, intent(in) :: jring,jchi,jk,jseed
    real(8), intent(in) :: frac,path
    real(8) :: y0,y1,t0,t1

    t0 = Tau_ring(jseed,jchi,jk,jring)*path
    t1 = Tau_ring(jseed+1,jchi,jk,jring)*path
    y0 = F_ring(jseed,jchi,jk,jring)*local_escape(t0)
    y1 = F_ring(jseed+1,jchi,jk,jring)*local_escape(t1)
    if (y0 > 0d0 .and. y1 > 0d0) then
        interpescidx = dexp(dlog(y0)+frac*dlog(y1/y0))
    else
        interpescidx = (1d0-frac)*y0 + frac*y1
    end if
end function interpescidx

real(8) function local_escape(taucell)
    implicit none
    real(8), intent(in) :: taucell

    if (taucell > 1d-6) then
        local_escape = (1d0-dexp(-taucell))/taucell
    else if (taucell > 0d0) then
        local_escape = 1d0-0.5d0*taucell+taucell*taucell/6d0
    else
        local_escape = 1d0
    end if
end function local_escape

subroutine radial_range(iev)
    implicit none
    integer, intent(in) :: iev
    real(8) :: tmp

    if (evnr(iev) == 0d0) then
        evslo(iev) = -1d300
        evshi(iev) = 1d300
        return
    end if
    evslo(iev) = evrlo(iev)/evnr(iev)
    evshi(iev) = evrhi(iev)/evnr(iev)
    if (evslo(iev) > evshi(iev)) then
        tmp = evslo(iev)
        evslo(iev) = evshi(iev)
        evshi(iev) = tmp
    end if
end subroutine radial_range

real(8) function pathfactor(iev,u,v)
    implicit none
    integer, intent(in) :: iev
    real(8), intent(in) :: u,v
    real(8) :: slo,shi,chord
    logical :: hit

    slo = evslo(iev)
    shi = evshi(iev)
    hit = .true.
    call cut_slab(u*evhalfth(iev),evntp(iev),-evhalfth(iev),evhalfth(iev),slo,shi,hit)
    call cut_slab(v*evhalfph(iev),evnpp(iev),-evhalfph(iev),evhalfph(iev),slo,shi,hit)
    if (hit .and. shi > slo) then
        chord = shi - slo
        pathfactor = chord/evdrcom(iev)
    else
        pathfactor = 0d0
    end if
end function pathfactor

subroutine cut_slab(x0,ndir,xlo,xhi,slo,shi,hit)
    implicit none
    real(8), intent(in) :: x0,ndir,xlo,xhi
    real(8), intent(inout) :: slo,shi
    logical, intent(inout) :: hit
    real(8) :: sa,sb,tmp

    if (.not. hit) return
    if (ndir == 0d0) then
        if (x0 < xlo .or. x0 > xhi) hit = .false.
        return
    end if
    sa = (xlo-x0)/ndir
    sb = (xhi-x0)/ndir
    if (sa > sb) then
        tmp = sa
        sa = sb
        sb = tmp
    end if
    if (sa > slo) slo = sa
    if (sb < shi) shi = sb
    if (shi <= slo) hit = .false.
end subroutine cut_slab

real(8) function outer_face(iring,ichi,k2)
    implicit none
    integer, intent(in) :: iring,ichi,k2

    if (ichi == 1) then
        outer_face = R_front(k2,iring)
    else
        outer_face = 0.5d0*(R_chi(ichi-1,k2,iring)+R_chi(ichi,k2,iring))
    end if
end function outer_face

real(8) function inner_face(iring,ichi,k2)
    implicit none
    integer, intent(in) :: iring,ichi,k2

    if (Num_chi == 1) then
        inner_face = R_chi(ichi,k2,iring) - (R_front(k2,iring)-R_chi(ichi,k2,iring))
    else if (ichi < Num_chi) then
        inner_face = 0.5d0*(R_chi(ichi,k2,iring)+R_chi(ichi+1,k2,iring))
    else
        inner_face = R_chi(ichi,k2,iring) - 0.5d0*(R_chi(ichi-1,k2,iring)-R_chi(ichi,k2,iring))
    end if
end function inner_face

real(8) function interplogidx(values,jseed,frac)
    implicit none
    real(8), intent(in), dimension(Num_nu) :: values
    integer, intent(in) :: jseed
    real(8), intent(in) :: frac
    real(8) :: y0,y1

    y0 = values(jseed)
    y1 = values(jseed+1)
    if (y0 > 0d0 .and. y1 > 0d0) then
        interplogidx = dexp(dlog(y0)+frac*dlog(y1/y0))
    else
        interplogidx = (1d0-frac)*y0 + frac*y1
    end if
end function interplogidx

real(8) function ray_flux()
    implicit none
    integer :: iev
    real(8) :: taufront,esc

    ray_flux = 0d0
    if (Num_chi == 1) then
        do iev = 1, nev
            if (evsrc(iev) <= 0d0 .and. evtau(iev) <= 0d0) cycle
            if (evtau(iev) > 1d-6) then
                esc = (1d0-dexp(-evtau(iev)))/evtau(iev)
            else if (evtau(iev) > 0d0) then
                esc = 1d0-0.5d0*evtau(iev)+evtau(iev)*evtau(iev)/6d0
            else
                esc = 1d0
            end if
            ray_flux = ray_flux + evsrc(iev)*esc
        end do
        return
    end if
    do iev = 1, nev
        if (evsrc(iev) <= 0d0 .and. evtau(iev) <= 0d0) cycle
        taufront = 0d0
        do ihit = hitstart(iev), hitstart(iev+1)-1
            taufront = taufront + evtau0(hitfg(ihit))*hitpath(ihit)
        end do
        if (evtau(iev) > 1d-6) then
            esc = dexp(-taufront)*(1d0-dexp(-evtau(iev)))/evtau(iev)
        else if (evtau(iev) > 0d0) then
            esc = dexp(-taufront)*(1d0-0.5d0*evtau(iev)+evtau(iev)*evtau(iev)/6d0)
        else
            esc = dexp(-taufront)
        end if
        ray_flux = ray_flux + evsrc(iev)*esc
    end do
end function ray_flux

logical function hitpatch(isrc,jfg,path)
    implicit none
    integer, intent(in) :: isrc,jfg
    real(8), intent(out) :: path
    real(8) :: dx,dy,u,v

    path = 0d0
    if (evx(isrc) < evxmin(jfg) .or. evx(isrc) > evxmax(jfg)) then
        hitpatch = .false.
        return
    end if
    if (evy(isrc) < evymin(jfg) .or. evy(isrc) > evymax(jfg)) then
        hitpatch = .false.
        return
    end if
    dx = evx(isrc) - evx(jfg)
    dy = evy(isrc) - evy(jfg)
    u = (dx*evby(jfg) - dy*evbx(jfg))*evdetinv(jfg)
    v = (evax(jfg)*dy - evay(jfg)*dx)*evdetinv(jfg)
    if (dabs(u) <= 1d0 .and. dabs(v) <= 1d0) then
        path = pathfactor(jfg,u,v)
        hitpatch = path > 0d0
    else
        hitpatch = .false.
    end if
end function hitpatch
end subroutine sed_chiring_batchlum_ray

subroutine sed_chiring_core(Boundary,R_Tobs1,R_front, &
                             F_ring,Tau_ring,R_chi,Gamma_chi,Chi_weight,V_seed,V_obs,Tobs,theta_lo,theta_hi, &
                             fluxscale, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R,F_tot_obs)
    use constants
    use interpolation_common
    implicit none
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R
    real(8), intent(in) :: Boundary(n), R_Tobs1(Num_R), R_front(Num_R), Tobs(Num_Tobs), V_seed(Num_nu), V_obs(Num_nu_obs)
    real(8), intent(in) :: F_ring(Num_nu,Num_chi,Num_R), Tau_ring(Num_nu,Num_chi,Num_R), R_chi(Num_chi,Num_R)
    real(8), intent(in) :: Gamma_chi(Num_chi,Num_R), Chi_weight(Num_chi,Num_R), theta_lo, theta_hi, fluxscale
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    real(8), allocatable :: F_temp(:,:), V_obs_log(:), V_seed_log(:), Tau_prefix(:,:,:)
    real(8) :: R_Tobs_chi(Num_R), ldomega,lgamlo,lgamhi,segment_lo,segment_hi
    real(8) :: T_sorted(Num_Tobs)
    integer :: T_order(Num_Tobs)
    real(8) :: z,Tv,cos_tv,sin_tv,theta_center,domega,dPhi,DMu,Ratio,phi_lo,phi_hi
    integer :: I_chi,i_Phi,Iobs,K2,II,last_k2,k_start
    logical :: monotonic_chi,time_ordered

    allocate(F_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))
    allocate(Tau_prefix(Num_nu,0:Num_chi,Num_R))
    F_tot_obs = 0d0
    F_temp = 0d0
    z = Boundary(8)
    Tv = Boundary(10)
    dPhi = 2d0*pi/Num_phi_patch
    V_obs_log = dlog(V_obs)
    V_seed_log = dlog(V_seed)
    call time_order(Tobs,Num_Tobs,T_order,T_sorted,time_ordered)
    cos_tv = dcos(Tv)
    sin_tv = dsin(Tv)

    Tau_prefix(:,0,:) = 0d0
    do I_chi = 1, Num_chi
        Tau_prefix(:,I_chi,:) = Tau_prefix(:,I_chi-1,:) + Tau_ring(:,I_chi,:)
    end do
    theta_center = 0.5d0*(theta_lo+theta_hi)
    do i_Phi = 1, Num_phi_patch
        phi_lo = (i_Phi-1)*dPhi
        phi_hi = i_Phi*dPhi
        call projphisamp(phi_lo,phi_hi,0.5d0*(phi_lo+phi_hi),F_temp)
    end do
    if (time_ordered) then
        F_tot_obs=F_temp
    else
        do Iobs=1,Num_Tobs
            F_tot_obs(:,T_order(Iobs))=F_temp(:,Iobs)
        end do
    end if
    deallocate(F_temp,V_obs_log,V_seed_log,Tau_prefix)
    return

contains

subroutine projphisamp(phi_lo,phi_hi,phi_center,local_obs)
    implicit none
    real(8), intent(in) :: phi_lo,phi_hi,phi_center
    real(8), intent(inout), dimension(Num_nu_obs,Num_Tobs) :: local_obs
    integer :: lowerreal8

    domega = (dcos(theta_lo)-dcos(theta_hi))*(phi_hi-phi_lo)
    ldomega = dlog(domega)-dlog(4d0*pi)
    DMu = cos_tv*dcos(theta_center)+sin_tv*dsin(theta_center)*dcos(phi_center)
    do I_chi = 1, Num_chi
        R_Tobs_chi = R_Tobs1 + (1d0+z)*(R_front-R_chi(I_chi,:)*DMu)/Para_c
        monotonic_chi = all(R_Tobs_chi(2:Num_R) > R_Tobs_chi(1:Num_R-1))
        if (monotonic_chi) then
            II = 1
            last_k2 = 0
            do Iobs = 1, Num_Tobs
                if (T_sorted(Iobs) < R_Tobs_chi(1) .or. T_sorted(Iobs) >= R_Tobs_chi(Num_R)) cycle
                do while (II < Num_R-1 .and. T_sorted(Iobs) >= R_Tobs_chi(II+1))
                    II = II + 1
                end do
                K2 = II
                if (K2 /= last_k2) then
                    lgamlo = dlog(Gamma_chi(I_chi,K2))
                    lgamhi = dlog(Gamma_chi(I_chi,K2+1))
                    last_k2 = K2
                end if
                Ratio = (T_sorted(Iobs)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                call project_ring(I_chi,K2,Iobs,Ratio,DMu,ldomega, &
                                  lgamlo,lgamhi,local_obs)
            end do
        else
            do K2 = 1, Num_R-1
                segment_lo = min(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                segment_hi = max(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                k_start = lowerreal8(T_sorted,Num_Tobs,segment_lo)
                lgamlo = dlog(Gamma_chi(I_chi,K2))
                lgamhi = dlog(Gamma_chi(I_chi,K2+1))
                do Iobs = k_start, Num_Tobs
                    if (T_sorted(Iobs) >= segment_hi) exit
                    Ratio = (T_sorted(Iobs)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                    call project_ring(I_chi,K2,Iobs,Ratio,DMu,ldomega, &
                                      lgamlo,lgamhi,local_obs)
                end do
            end do
        end if
    end do
end subroutine projphisamp

subroutine project_ring(I_chi,K2,K1,Ratio,DMu,ldomega,lgamlo,lgamhi,local_obs)
    implicit none
    integer, intent(in) :: I_chi,K2,K1
    real(8), intent(in) :: Ratio,DMu,ldomega,lgamlo,lgamhi
    real(8), intent(inout), dimension(Num_nu_obs,Num_Tobs) :: local_obs
    real(8), dimension(Num_nu) :: Ftheta
    real(8) :: ldopred,lfluxw,doppler,DG,Beta,targetlo,targethi
    integer :: lo,hi,nwin,lowerreal8

    DG = dexp(lgamlo+Ratio*(lgamhi-lgamlo))
    Beta = dsqrt(1d0-DG**(-2))
    doppler = DG*(1d0-Beta*DMu)
    ldopred = dlog(doppler)+dlog(1d0+z)
    targetlo = V_obs_log(1) + ldopred
    targethi = V_obs_log(Num_nu_obs) + ldopred
    lo = max(1,lowerreal8(V_seed_log,Num_nu,targetlo)-1)
    hi = min(Num_nu,lowerreal8(V_seed_log,Num_nu,targethi))
    nwin = hi - lo + 1
    if (nwin < 2) return
    call ring_source(I_chi,K2,Ratio,lo,hi,Ftheta(lo:hi))
    lfluxw = ldomega - 3d0*dlog(doppler)
    call accum_shifted(V_seed_log(lo:hi),Ftheta(lo:hi),nwin,V_obs_log,Num_nu_obs, &
                                                     ldopred,lfluxw,local_obs(:,K1))
end subroutine project_ring

subroutine ring_source(I_chi,K2,Ratio,lo,hi,Ftheta)
    implicit none
    integer, intent(in) :: I_chi,K2,lo,hi
    real(8), intent(in) :: Ratio
    real(8), intent(out), dimension(hi-lo+1) :: Ftheta
    real(8) :: source_lo,source_hi,tau_front_lo,tau_front_hi,tau_cell_lo,tau_cell_hi
    real(8), external :: chi_escape
    integer :: I_nu,j

    do I_nu = lo, hi
        j = I_nu - lo + 1
        tau_front_lo = Tau_prefix(I_nu,I_chi-1,K2)
        tau_front_hi = Tau_prefix(I_nu,I_chi-1,K2+1)
        tau_cell_lo = Tau_ring(I_nu,I_chi,K2)
        tau_cell_hi = Tau_ring(I_nu,I_chi,K2+1)
        source_lo = F_ring(I_nu,I_chi,K2)*fluxscale*Chi_weight(I_chi,K2) &
            * chi_escape(tau_front_lo,tau_cell_lo)
        source_hi = F_ring(I_nu,I_chi,K2+1)*fluxscale*Chi_weight(I_chi,K2+1) &
            * chi_escape(tau_front_hi,tau_cell_hi)
        Ftheta(j) = (1d0-Ratio)*source_lo + Ratio*source_hi
    end do
end subroutine ring_source
end subroutine sed_chiring_core

! Shared projection helpers.
real(8) function chi_escape(tau_front,tau_cell)
    use constants
    implicit none
    real(8), intent(in) :: tau_front,tau_cell
    if (tau_cell > 1d-6) then
        chi_escape = dexp(-tau_front)*(1d0-dexp(-tau_cell))/tau_cell
    else if (tau_cell > 0d0) then
        ! Taylor expansion of the removable tau_cell -> 0 singularity.
        chi_escape = dexp(-tau_front)*(1d0 - 0.5d0*tau_cell + tau_cell*tau_cell/6d0)
    else
        chi_escape = dexp(-tau_front)
    end if
end function chi_escape

integer function lowerreal8(values,n,x)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in), dimension(n) :: values
    real(8), intent(in) :: x
    integer :: lo,hi,mid
    lo = 1
    hi = n + 1
    do while (lo < hi)
        mid = (lo + hi) / 2
        if (values(mid) < x) then
            lo = mid + 1
        else
            hi = mid
        end if
    end do
    lowerreal8 = lo
end function lowerreal8
