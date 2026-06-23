! Created by rj on 2021/2/4.

!##################################################################################################
!When an intrinsic SED 'F_tot' observed in the comoving frame send in,
!to produce the observed SED 'F_tot_obs' after consider the EATS and Doppler boosting effect.
!##################################################################################################

! 将共动系SED插值到观测系：EATS时间修正 + 多普勒增亮 + 红移，均匀角网格。
subroutine sed_interpolation(Boundary,R_Tobs1,R_gamma,R,F_tot,V_seed,V_obs,Tobs, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads, F_tot_obs)
    !$ use omp_lib
    use constants
    use interpolation_common
    IMPLICIT REAL(8)(A-H,O-Z)
    !##############################################################################################
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads
    real(8), intent(in) :: Boundary(n)
    real(8), intent(in) :: Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: R_Tobs1(Num_R),R_gamma(Num_R),R(Num_R),F_tot(Num_nu,Num_R)
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    
    
    real(8), allocatable :: F_tot_obs_temp(:,:),V_obs_log(:),V_seed_log(:)
    real(8) :: R_Tobs_theta(Num_R),F_tot_theta(Num_nu),F_tot_log_theta(Num_nu), &
               V_seed_log_theta(Num_nu),log_gamma_lo,log_gamma_hi,log_domega_4pi, &
               log_doppler_redshift
    integer :: last_k2
    allocate (F_tot_obs_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))

    F_tot_obs=zero
    F_tot_obs_temp=zero
    
    G00 = Boundary(1)
    R00 = Boundary(4)
    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)
    
    dPhi = pi / Num_Phi
    phi_scale = two
    if (Num_Phi == 1) then
        dPhi = pi / 1440d0
        phi_scale = two * 1440d0
    end if
    dtheta=OpeningAngle_jet/Num_Theta

    V_obs_log = log(V_obs)
    V_seed_log = log(V_seed)
    !$OMP PARALLEL num_threads(n_threads), reduction(+:F_tot_obs_temp), private(I_Theta, Taa_lower, Taa_boundary, &
    !$OMP& Taa_center, domega, i_Phi, Phi_center, DMu, K1, II, K2, Ratio, DG, Beta, doppler, &
    !$OMP& R_Tobs_theta, F_tot_theta, F_tot_log_theta, V_seed_log_theta, &
    !$OMP& last_k2, log_gamma_lo, log_gamma_hi, log_domega_4pi, log_doppler_redshift)
    !$OMP DO SCHEDULE(GUIDED,4)
    do I_Theta=1,Num_Theta
       Taa_lower=dtheta*(I_Theta-1)
       Taa_boundary=dtheta*I_Theta
       Taa_center=dtheta*(I_Theta-0.5)
       domega=(dcos(Taa_lower)-dcos(Taa_boundary))*dPhi
       log_domega_4pi=log(domega)-log(4.0d0*pi)
       do i_Phi=1,Num_Phi
          Phi_center=(i_Phi-0.5)*dPhi
          DMu=dcos(Tv)*dcos(Taa_center)+dsin(Tv)*dsin(Taa_center)*dcos(Phi_center)
          R_Tobs_theta=R_Tobs1+R*(one-DMu)*(one+z)/Para_c 
          II=1
          last_k2=0
          do K1=1,Num_Tobs
             if (Tobs(K1) < R_Tobs_theta(1).or.Tobs(K1) > R_Tobs_theta(Num_R)) cycle
             do while (II < Num_R-1 .and. Tobs(K1) >= R_Tobs_theta(II+1))
                II=II+1
             end do
             K2=II
             if ((Tobs(K1) >= R_Tobs_theta(K2)).and.(Tobs(K1) < R_Tobs_theta(K2+1))) then
                 if (K2 /= last_k2) then
                     log_gamma_lo=log(R_gamma(K2))
                     log_gamma_hi=log(R_gamma(K2+1))
                     last_k2=K2
                 end if
                 Ratio=(Tobs(K1)-R_Tobs_theta(K2))/(R_Tobs_theta(K2+1)-R_Tobs_theta(K2))
                 DG=exp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
                 ! 时间方向用非负通量重构；频率方向仍在 interpolation_common 中按 log SED 插值。
                 F_tot_theta=(one-Ratio)*F_tot(:,K2)+Ratio*F_tot(:,K2+1)
                 F_tot_log_theta=-huge(one)
                 where(F_tot_theta > zero) F_tot_log_theta=log(F_tot_theta)
                 Beta=dsqrt(one-DG**(-2))

                 doppler=DG*(one-Beta*DMu) !Doppler factor, changed with R
                 log_doppler_redshift=log(doppler)+log(one+z)
                 F_tot_log_theta=F_tot_log_theta+log_domega_4pi-3d0*log(doppler)
                 V_seed_log_theta=V_seed_log-log_doppler_redshift
                 call interpolation_accumulate_log_sed(V_seed_log_theta,F_tot_log_theta, &
                                                       Num_nu,V_obs_log,Num_nu_obs,F_tot_obs_temp(:,K1))
             end if
         end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    F_tot_obs=F_tot_obs_temp*phi_scale
    
    deallocate (F_tot_obs_temp,V_obs_log,V_seed_log)
    
    
    return
end subroutine sed_interpolation

! 壳层级EATS自适应角向积分：每个基础theta cell用一阶/二阶中点规则估计角向积分误差并递归细分。
subroutine sed_interpolation_adaptive_theta(Boundary,R_Tobs1,R_gamma,R,F_tot,V_seed,V_obs,Tobs, &
                             adaptive_rtol,adaptive_max_depth, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads,F_tot_obs)
    !$ use omp_lib
    use constants
    use interpolation_common
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads
    integer, intent(in) :: adaptive_max_depth
    real(8), intent(in) :: Boundary(n),Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: R_Tobs1(Num_R),R_gamma(Num_R),R(Num_R),F_tot(Num_nu,Num_R)
    real(8), intent(in) :: adaptive_rtol
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    real(8), allocatable :: F_tot_obs_temp(:,:),V_obs_log(:),V_seed_log(:)
    real(8) :: cell_obs(Num_nu_obs,Num_Tobs)

    if (adaptive_max_depth == 0) then
        call sed_interpolation(Boundary,R_Tobs1,R_gamma,R,F_tot,V_seed,V_obs,Tobs, &
                               n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads,F_tot_obs)
        return
    end if

    allocate(F_tot_obs_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))

    F_tot_obs = zero
    F_tot_obs_temp = zero
    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)
    dPhi = pi/Num_Phi
    phi_scale = two
    if (Num_Phi == 1) then
        dPhi = pi/1440d0
        phi_scale = two*1440d0
    end if
    dtheta = OpeningAngle_jet/Num_Theta
    V_obs_log = dlog(V_obs)
    V_seed_log = dlog(V_seed)

    !$OMP PARALLEL num_threads(n_threads), reduction(+:F_tot_obs_temp), private(I_Theta, &
    !$OMP& i_Phi,Phi_center,Taa_lower,Taa_boundary,cell_obs)
    !$OMP DO SCHEDULE(GUIDED,4)
    do I_Theta = 1, Num_Theta
        Taa_lower = dtheta*(I_Theta-1)
        Taa_boundary = dtheta*I_Theta
        do i_Phi = 1, Num_Phi
            Phi_center = (i_Phi-0.5d0)*dPhi
            cell_obs = zero
            call integrate_theta_cell(Taa_lower,Taa_boundary,Phi_center,0,cell_obs)
            F_tot_obs_temp = F_tot_obs_temp + cell_obs
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    F_tot_obs = F_tot_obs_temp*phi_scale
    deallocate(F_tot_obs_temp,V_obs_log,V_seed_log)
    return

contains

recursive subroutine integrate_theta_cell(theta_lo,theta_hi,phi_center,depth,accum_obs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: depth
    real(8), intent(in) :: theta_lo,theta_hi,phi_center
    real(8), intent(inout) :: accum_obs(Num_nu_obs,Num_Tobs)
    real(8) :: theta_mid,theta_lmid,theta_rmid,err_norm,flux_norm
    real(8) :: coarse_obs(Num_nu_obs,Num_Tobs),fine_obs(Num_nu_obs,Num_Tobs)

    theta_mid = 0.5d0*(theta_lo+theta_hi)
    coarse_obs = zero
    call project_theta_sample(theta_lo,theta_hi,theta_mid,phi_center,coarse_obs)
    if (depth >= adaptive_max_depth) then
        accum_obs = accum_obs + coarse_obs
        return
    end if
    theta_lmid = 0.5d0*(theta_lo+theta_mid)
    theta_rmid = 0.5d0*(theta_mid+theta_hi)
    fine_obs = zero
    call project_theta_sample(theta_lo,theta_mid,theta_lmid,phi_center,fine_obs)
    call project_theta_sample(theta_mid,theta_hi,theta_rmid,phi_center,fine_obs)
    flux_norm = sum(abs(fine_obs))
    err_norm = sum(abs(fine_obs-coarse_obs))
    if (flux_norm == zero .or. err_norm <= adaptive_rtol*flux_norm) then
        accum_obs = accum_obs + fine_obs
    else
        call integrate_theta_cell(theta_lo,theta_mid,phi_center,depth+1,accum_obs)
        call integrate_theta_cell(theta_mid,theta_hi,phi_center,depth+1,accum_obs)
    end if
end subroutine integrate_theta_cell

subroutine project_theta_sample(theta_lo,theta_hi,theta_center,phi_center,local_obs)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: theta_lo,theta_hi,theta_center,phi_center
    real(8), intent(inout) :: local_obs(Num_nu_obs,Num_Tobs)
    real(8) :: R_Tobs_theta(Num_R),log_gamma_lo,log_gamma_hi,log_domega_4pi
    integer :: last_k2

    domega = (dcos(theta_lo)-dcos(theta_hi))*dPhi
    log_domega_4pi = dlog(domega)-dlog(4.0d0*pi)
    DMu = dcos(Tv)*dcos(theta_center)+dsin(Tv)*dsin(theta_center)*dcos(phi_center)
    R_Tobs_theta = R_Tobs1+R*(one-DMu)*(one+z)/Para_c
    II = 1
    last_k2 = 0
    do K1 = 1, Num_Tobs
        if (Tobs(K1) < R_Tobs_theta(1) .or. Tobs(K1) > R_Tobs_theta(Num_R)) cycle
        do while (II < Num_R-1 .and. Tobs(K1) >= R_Tobs_theta(II+1))
            II = II + 1
        end do
        K2 = II
        if (Tobs(K1) >= R_Tobs_theta(K2) .and. Tobs(K1) < R_Tobs_theta(K2+1)) then
            if (K2 /= last_k2) then
                log_gamma_lo = dlog(R_gamma(K2))
                log_gamma_hi = dlog(R_gamma(K2+1))
                last_k2 = K2
            end if
            Ratio = (Tobs(K1)-R_Tobs_theta(K2))/(R_Tobs_theta(K2+1)-R_Tobs_theta(K2))
            call project_shell_segment(K1,K2,Ratio,DMu,log_domega_4pi, &
                                       log_gamma_lo,log_gamma_hi,local_obs)
        end if
    end do
end subroutine project_theta_sample

subroutine project_shell_segment(K1,K2,Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi,local_obs)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: K1,K2
    real(8), intent(in) :: Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi
    real(8), intent(inout) :: local_obs(Num_nu_obs,Num_Tobs)
    real(8) :: F_tot_theta(Num_nu),F_tot_log_theta(Num_nu),V_seed_log_theta(Num_nu)
    real(8) :: log_doppler_redshift

    DG = dexp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
    F_tot_theta = (one-Ratio)*F_tot(:,K2)+Ratio*F_tot(:,K2+1)
    F_tot_log_theta = -huge(one)
    where(F_tot_theta > zero) F_tot_log_theta = dlog(F_tot_theta)
    Beta = dsqrt(one-DG**(-2))
    doppler = DG*(one-Beta*DMu)
    log_doppler_redshift = dlog(doppler)+dlog(one+z)
    F_tot_log_theta = F_tot_log_theta+log_domega_4pi-3d0*dlog(doppler)
    V_seed_log_theta = V_seed_log-log_doppler_redshift
    call interpolation_accumulate_log_sed(V_seed_log_theta,F_tot_log_theta, &
                                          Num_nu,V_obs_log,Num_nu_obs,local_obs(:,K1))
end subroutine project_shell_segment
end subroutine sed_interpolation_adaptive_theta

real(8) function chi_ssa_cell_escape(tau_front,tau_cell)
    use constants
    implicit real(8)(a-h,o-z)
    real(8), intent(in) :: tau_front,tau_cell
    if (tau_cell > 1d-6) then
        chi_ssa_cell_escape = dexp(-tau_front)*(one-dexp(-tau_cell))/tau_cell
    else if (tau_cell > zero) then
        ! Taylor expansion of the removable tau_cell -> 0 singularity.
        chi_ssa_cell_escape = dexp(-tau_front)*(one - 0.5d0*tau_cell + tau_cell*tau_cell/6d0)
    else
        chi_ssa_cell_escape = dexp(-tau_front)
    end if
end function chi_ssa_cell_escape

integer function lower_bound_real8(values,n,x)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: values(n),x
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
    lower_bound_real8 = lo
end function lower_bound_real8

! 将χ分辨有限厚壳层积分到完整top-hat角网格：θ/φ EATS + χ厚度EATS + 局域Doppler。
subroutine sed_interpolation_chi(Boundary,R_Tobs1,R_front,F_chi,Tau_chi,R_chi,Gamma_chi,Chi_weight,V_seed,V_obs,Tobs, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_chi,Num_R,n_threads,F_tot_obs)
    !$ use omp_lib
    use constants
    use interpolation_common
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_chi,Num_R,n_threads
    real(8), intent(in) :: Boundary(n),R_Tobs1(Num_R),R_front(Num_R),Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: F_chi(Num_nu,Num_chi,Num_R),Tau_chi(Num_nu,Num_chi,Num_R)
    real(8), intent(in) :: R_chi(Num_chi,Num_R),Gamma_chi(Num_chi,Num_R),Chi_weight(Num_chi,Num_R)
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    real(8), allocatable :: F_temp(:,:),V_obs_log(:),V_seed_log(:),Tau_prefix(:,:,:)
    real(8), allocatable :: Angle_dmu(:),Angle_log_domega(:)
    real(8) :: R_Tobs_chi(Num_R),log_domega_4pi
    real(8) :: log_gamma_lo,log_gamma_hi,segment_lo,segment_hi,cos_tv,sin_tv
    integer :: last_k2,k_start,lower_bound_real8,I_ang
    logical :: monotonic_chi

    allocate(F_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))
    allocate(Tau_prefix(Num_nu,0:Num_chi,Num_R))
    allocate(Angle_dmu(Num_Theta*Num_Phi),Angle_log_domega(Num_Theta*Num_Phi))
    F_tot_obs = zero
    F_temp = zero
    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)
    dPhi = pi/Num_Phi
    phi_scale = two
    if (Num_Phi == 1) then
        dPhi = pi/1440d0
        phi_scale = two*1440d0
    end if
    dtheta = OpeningAngle_jet/Num_Theta
    V_obs_log = dlog(V_obs)
    V_seed_log = dlog(V_seed)
    Tau_prefix(:,0,:) = zero
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
        Angle_log_domega(I_ang) = dlog(domega)-dlog(4d0*pi)
        Angle_dmu(I_ang) = cos_tv*dcos(Taa_center)+sin_tv*dsin(Taa_center)*dcos(Phi_center)
    end do

    !$OMP PARALLEL num_threads(n_threads), reduction(+:F_temp), private(I_ang,log_domega_4pi,DMu, &
    !$OMP& I_chi,R_Tobs_chi,monotonic_chi,II,last_k2,K1,K2,Ratio, &
    !$OMP& log_gamma_lo,log_gamma_hi,segment_lo,segment_hi,k_start)
    !$OMP DO SCHEDULE(STATIC)
    do I_ang = 1, Num_Theta*Num_Phi
        log_domega_4pi = Angle_log_domega(I_ang)
        DMu = Angle_dmu(I_ang)
        do I_chi = 1, Num_chi
                R_Tobs_chi = R_Tobs1 + (one+z)*(R_front - R_chi(I_chi,:)*DMu)/Para_c
                monotonic_chi = all(R_Tobs_chi(2:Num_R) > R_Tobs_chi(1:Num_R-1))
                if (monotonic_chi) then
                    II = 1
                    last_k2 = 0
                    do K1 = 1, Num_Tobs
                        if (Tobs(K1) < R_Tobs_chi(1) .or. Tobs(K1) > R_Tobs_chi(Num_R)) cycle
                        do while (II < Num_R-1 .and. Tobs(K1) >= R_Tobs_chi(II+1))
                            II = II + 1
                        end do
                        K2 = II
                        if (Tobs(K1) >= R_Tobs_chi(K2) .and. Tobs(K1) < R_Tobs_chi(K2+1)) then
                            if (K2 /= last_k2) then
                                log_gamma_lo = dlog(Gamma_chi(I_chi,K2))
                                log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1))
                                last_k2 = K2
                            end if
                            Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                            call project_chi_segment_flux(I_chi,K2,K1,Ratio,DMu,log_domega_4pi, &
                                                          log_gamma_lo,log_gamma_hi)
                        end if
                    end do
                else
                    do K2 = 1, Num_R-1
                        segment_lo = min(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                        segment_hi = max(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                        k_start = lower_bound_real8(Tobs,Num_Tobs,segment_lo)
                        log_gamma_lo = dlog(Gamma_chi(I_chi,K2))
                        log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1))
                        do K1 = k_start, Num_Tobs
                            if (Tobs(K1) >= segment_hi) exit
                            Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                            call project_chi_segment_flux(I_chi,K2,K1,Ratio,DMu,log_domega_4pi, &
                                                          log_gamma_lo,log_gamma_hi)
                        end do
                    end do
                end if
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    F_tot_obs = F_temp*phi_scale
    deallocate(F_temp,V_obs_log,V_seed_log,Tau_prefix,Angle_dmu,Angle_log_domega)
    return

contains

subroutine project_chi_segment_flux(I_chi,K2,K1,Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_chi,K2,K1
    real(8), intent(in) :: Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi
    real(8) :: F_theta(Num_nu)
    real(8) :: log_doppler_redshift,log_flux_weight,doppler

    call compute_chi_segment_state(I_chi,K2,Ratio,DMu,log_gamma_lo,log_gamma_hi, &
                                   log_doppler_redshift,doppler,F_theta)
    log_flux_weight = log_domega_4pi - 3d0*dlog(doppler)
    call interpolation_accumulate_shifted_linear_sed(V_seed_log,F_theta,Num_nu,V_obs_log,Num_nu_obs, &
                                                     log_doppler_redshift,log_flux_weight,F_temp(:,K1))
end subroutine project_chi_segment_flux

subroutine compute_chi_segment_state(I_chi,K2,Ratio,DMu,log_gamma_lo,log_gamma_hi, &
                                     log_doppler_redshift,doppler,F_theta)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_chi,K2
    real(8), intent(in) :: Ratio,DMu,log_gamma_lo,log_gamma_hi
    real(8), intent(out) :: log_doppler_redshift,doppler,F_theta(Num_nu)
    real(8) :: DG,Beta

    DG = dexp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
    Beta = dsqrt(one-DG**(-2))
    doppler = DG*(one-Beta*DMu)
    log_doppler_redshift = dlog(doppler)+dlog(one+z)
    call accumulate_chi_cell_source(I_chi,K2,Ratio,F_theta)
end subroutine compute_chi_segment_state

subroutine accumulate_chi_cell_source(I_chi,K2,Ratio,F_theta)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_chi,K2
    real(8), intent(in) :: Ratio
    real(8), intent(out) :: F_theta(Num_nu)
    real(8) :: source_lo,source_hi,tau_front_lo,tau_front_hi,tau_cell_lo,tau_cell_hi
    integer :: I_nu

    F_theta = zero
    do I_nu = 1, Num_nu
        tau_front_lo = Tau_prefix(I_nu,I_chi-1,K2)
        tau_front_hi = Tau_prefix(I_nu,I_chi-1,K2+1)
        tau_cell_lo = Tau_chi(I_nu,I_chi,K2)
        tau_cell_hi = Tau_chi(I_nu,I_chi,K2+1)
        source_lo = F_chi(I_nu,I_chi,K2)*Chi_weight(I_chi,K2)*chi_ssa_cell_escape(tau_front_lo,tau_cell_lo)
        source_hi = F_chi(I_nu,I_chi,K2+1)*Chi_weight(I_chi,K2+1)*chi_ssa_cell_escape(tau_front_hi,tau_cell_hi)
        F_theta(I_nu) = (one-Ratio)*source_lo + Ratio*source_hi
    end do
end subroutine accumulate_chi_cell_source
end subroutine sed_interpolation_chi

! 将轴对称结构化喷流的χ分辨有限厚壳层一次性积分到观测系。
subroutine sed_interpolation_chi_structured_axisym(Boundary,R_Tobs1,R_front,F_chi,Tau_chi,R_chi,Gamma_chi,Chi_weight, &
                             V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs,Num_Tobs,Num_theta_patch,Num_phi_patch,Num_chi,Num_R, &
                             n_threads,F_tot_obs)
    !$ use omp_lib
    use constants
    use interpolation_common
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_theta_patch,Num_phi_patch,Num_chi,Num_R,n_threads
    real(8), intent(in) :: Boundary(n),R_Tobs1(Num_R,Num_theta_patch),R_front(Num_R,Num_theta_patch)
    real(8), intent(in) :: Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: F_chi(Num_nu,Num_chi,Num_R,Num_theta_patch)
    real(8), intent(in) :: Tau_chi(Num_nu,Num_chi,Num_R,Num_theta_patch)
    real(8), intent(in) :: R_chi(Num_chi,Num_R,Num_theta_patch),Gamma_chi(Num_chi,Num_R,Num_theta_patch)
    real(8), intent(in) :: Chi_weight(Num_chi,Num_R,Num_theta_patch)
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    real(8), allocatable :: F_temp(:,:),V_obs_log(:),V_seed_log(:),Tau_prefix(:,:,:,:)
    real(8) :: R_Tobs_chi(Num_R),log_domega_4pi,log_gamma_lo,log_gamma_hi,segment_lo,segment_hi
    real(8) :: cos_tv,sin_tv,theta_lo,theta_hi,theta_center,phi_center,domega
    integer :: last_k2,k_start,lower_bound_real8
    logical :: monotonic_chi

    allocate(F_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))
    allocate(Tau_prefix(Num_nu,0:Num_chi,Num_R,Num_theta_patch))
    F_tot_obs = zero
    F_temp = zero
    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)
    dtheta = OpeningAngle_jet/Num_theta_patch
    dPhi = two*pi/Num_phi_patch
    V_obs_log = dlog(V_obs)
    V_seed_log = dlog(V_seed)
    Tau_prefix(:,0,:,:) = zero
    do I_Theta = 1, Num_theta_patch
        do I_chi = 1, Num_chi
            Tau_prefix(:,I_chi,:,I_Theta) = Tau_prefix(:,I_chi-1,:,I_Theta) + Tau_chi(:,I_chi,:,I_Theta)
        end do
    end do
    cos_tv = dcos(Tv)
    sin_tv = dsin(Tv)

    !$OMP PARALLEL num_threads(n_threads), reduction(+:F_temp), private(I_Theta,i_Phi,theta_lo,theta_hi, &
    !$OMP& theta_center,phi_center,domega,log_domega_4pi,DMu,I_chi,R_Tobs_chi,monotonic_chi, &
    !$OMP& II,last_k2,K1,K2,Ratio,log_gamma_lo,log_gamma_hi,segment_lo,segment_hi,k_start)
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
    do I_Theta = 1, Num_theta_patch
        do i_Phi = 1, Num_phi_patch
            theta_lo = dtheta*(I_Theta-1)
            theta_hi = dtheta*I_Theta
            theta_center = dtheta*(I_Theta-0.5d0)
            phi_center = (i_Phi-0.5d0)*dPhi
            domega = (dcos(theta_lo)-dcos(theta_hi))*dPhi
            log_domega_4pi = dlog(domega)-dlog(4d0*pi)
            DMu = cos_tv*dcos(theta_center)+sin_tv*dsin(theta_center)*dcos(phi_center)
            do I_chi = 1, Num_chi
                R_Tobs_chi = R_Tobs1(:,I_Theta) + (one+z)*(R_front(:,I_Theta)-R_chi(I_chi,:,I_Theta)*DMu)/Para_c
                monotonic_chi = all(R_Tobs_chi(2:Num_R) > R_Tobs_chi(1:Num_R-1))
                if (monotonic_chi) then
                    II = 1
                    last_k2 = 0
                    do K1 = 1, Num_Tobs
                        if (Tobs(K1) < R_Tobs_chi(1) .or. Tobs(K1) > R_Tobs_chi(Num_R)) cycle
                        do while (II < Num_R-1 .and. Tobs(K1) >= R_Tobs_chi(II+1))
                            II = II + 1
                        end do
                        K2 = II
                        if (Tobs(K1) >= R_Tobs_chi(K2) .and. Tobs(K1) < R_Tobs_chi(K2+1)) then
                            if (K2 /= last_k2) then
                                log_gamma_lo = dlog(Gamma_chi(I_chi,K2,I_Theta))
                                log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1,I_Theta))
                                last_k2 = K2
                            end if
                            Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                            call project_structured_chi_segment(I_Theta,I_chi,K2,K1,Ratio,DMu,log_domega_4pi, &
                                                               log_gamma_lo,log_gamma_hi)
                        end if
                    end do
                else
                    do K2 = 1, Num_R-1
                        segment_lo = min(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                        segment_hi = max(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                        k_start = lower_bound_real8(Tobs,Num_Tobs,segment_lo)
                        log_gamma_lo = dlog(Gamma_chi(I_chi,K2,I_Theta))
                        log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1,I_Theta))
                        do K1 = k_start, Num_Tobs
                            if (Tobs(K1) >= segment_hi) exit
                            Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                            call project_structured_chi_segment(I_Theta,I_chi,K2,K1,Ratio,DMu,log_domega_4pi, &
                                                               log_gamma_lo,log_gamma_hi)
                        end do
                    end do
                end if
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    F_tot_obs = F_temp
    deallocate(F_temp,V_obs_log,V_seed_log,Tau_prefix)
    return

contains

subroutine project_structured_chi_segment(I_Theta,I_chi,K2,K1,Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_Theta,I_chi,K2,K1
    real(8), intent(in) :: Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi
    real(8) :: F_theta(Num_nu),log_doppler_redshift,log_flux_weight,doppler

    call compute_structured_chi_state(I_Theta,I_chi,K2,Ratio,DMu,log_gamma_lo,log_gamma_hi, &
                                      log_doppler_redshift,doppler,F_theta)
    log_flux_weight = log_domega_4pi - 3d0*dlog(doppler)
    call interpolation_accumulate_shifted_linear_sed(V_seed_log,F_theta,Num_nu,V_obs_log,Num_nu_obs, &
                                                     log_doppler_redshift,log_flux_weight,F_temp(:,K1))
end subroutine project_structured_chi_segment

subroutine compute_structured_chi_state(I_Theta,I_chi,K2,Ratio,DMu,log_gamma_lo,log_gamma_hi, &
                                        log_doppler_redshift,doppler,F_theta)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_Theta,I_chi,K2
    real(8), intent(in) :: Ratio,DMu,log_gamma_lo,log_gamma_hi
    real(8), intent(out) :: log_doppler_redshift,doppler,F_theta(Num_nu)
    real(8) :: DG,Beta

    DG = dexp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
    Beta = dsqrt(one-DG**(-2))
    doppler = DG*(one-Beta*DMu)
    log_doppler_redshift = dlog(doppler)+dlog(one+z)
    call accumulate_structured_chi_source(I_Theta,I_chi,K2,Ratio,F_theta)
end subroutine compute_structured_chi_state

subroutine accumulate_structured_chi_source(I_Theta,I_chi,K2,Ratio,F_theta)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_Theta,I_chi,K2
    real(8), intent(in) :: Ratio
    real(8), intent(out) :: F_theta(Num_nu)
    real(8) :: source_lo,source_hi,tau_front_lo,tau_front_hi,tau_cell_lo,tau_cell_hi
    integer :: I_nu

    F_theta = zero
    do I_nu = 1, Num_nu
        tau_front_lo = Tau_prefix(I_nu,I_chi-1,K2,I_Theta)
        tau_front_hi = Tau_prefix(I_nu,I_chi-1,K2+1,I_Theta)
        tau_cell_lo = Tau_chi(I_nu,I_chi,K2,I_Theta)
        tau_cell_hi = Tau_chi(I_nu,I_chi,K2+1,I_Theta)
        source_lo = F_chi(I_nu,I_chi,K2,I_Theta)*Chi_weight(I_chi,K2,I_Theta) &
                    *chi_ssa_cell_escape(tau_front_lo,tau_cell_lo)
        source_hi = F_chi(I_nu,I_chi,K2+1,I_Theta)*Chi_weight(I_chi,K2+1,I_Theta) &
                    *chi_ssa_cell_escape(tau_front_hi,tau_cell_hi)
        F_theta(I_nu) = (one-Ratio)*source_lo + Ratio*source_hi
    end do
end subroutine accumulate_structured_chi_source
end subroutine sed_interpolation_chi_structured_axisym

! direct-electronχ投影的批量同步谱版本：先按V_seed生成χ谱，再进入快速EATS累加。
subroutine sed_interpolation_chi_electron_cached(Boundary,R_Tobs1,R_front,DNe_chi,B_chi,R_chi,Gamma_chi,Chi_weight, &
                             gam_e,V_seed,V_obs,Tobs,n,Num_gam_e,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi, &
                             Num_chi,Num_R,n_threads,F_tot_obs)
    use constants
    use radiation_common, only: radiation_syn_flux_tau_chi_batch_core
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_gam_e,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_chi,Num_R,n_threads
    real(8), intent(in) :: Boundary(n),R_Tobs1(Num_R),R_front(Num_R),Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: gam_e(Num_gam_e),DNe_chi(Num_gam_e,Num_chi,Num_R),B_chi(Num_chi,Num_R)
    real(8), intent(in) :: R_chi(Num_chi,Num_R),Gamma_chi(Num_chi,Num_R),Chi_weight(Num_chi,Num_R)
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    real(8), allocatable :: F_chi(:,:,:),Tau_chi(:,:,:),P_syn(:,:),Tau_syn(:,:)

    DL = Boundary(13)
    z = Boundary(8)
    flux_prefactor = (one+z)/(4d0*pi*DL*DL)
    allocate(F_chi(Num_nu,Num_chi,Num_R),Tau_chi(Num_nu,Num_chi,Num_R), &
             P_syn(Num_nu,Num_chi),Tau_syn(Num_nu,Num_chi))
    do K2 = 1, Num_R
        call radiation_syn_flux_tau_chi_batch_core(R_front(K2),Num_gam_e,Num_nu,Num_chi,gam_e,DNe_chi(:,:,K2), &
                                                   V_seed,B_chi(:,K2),1.046d4,P_syn,Tau_syn)
        F_chi(:,:,K2) = P_syn*flux_prefactor
        Tau_chi(:,:,K2) = Tau_syn
    end do
    call sed_interpolation_chi(Boundary,R_Tobs1,R_front,F_chi,Tau_chi,R_chi,Gamma_chi,Chi_weight,V_seed,V_obs,Tobs, &
                               n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_chi,Num_R,n_threads,F_tot_obs)
    deallocate(F_chi,Tau_chi,P_syn,Tau_syn)
end subroutine sed_interpolation_chi_electron_cached

! 结构化direct-electronχ投影的批量同步谱版本：单次调用内逐theta ring流式生成χ谱并累加。
subroutine sed_interpolation_chi_structured_axisym_electron_cached(Boundary,R_Tobs1,R_front,DNe_chi,B_chi,R_chi, &
                             Gamma_chi,Chi_weight,gam_e,V_seed,V_obs,Tobs,n,Num_gam_e,Num_nu,Num_nu_obs, &
                             Num_Tobs,Num_theta_patch,Num_phi_patch,Num_chi,Num_R,n_threads,F_tot_obs)
    use constants
    use interpolation_common
    use radiation_common, only: radiation_syn_flux_tau_chi_batch_core
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_gam_e,Num_nu,Num_nu_obs,Num_Tobs,Num_theta_patch,Num_phi_patch,Num_chi,Num_R,n_threads
    real(8), intent(in) :: Boundary(n),R_Tobs1(Num_R,Num_theta_patch),R_front(Num_R,Num_theta_patch)
    real(8), intent(in) :: Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs),gam_e(Num_gam_e)
    real(8), intent(in) :: DNe_chi(Num_gam_e,Num_chi,Num_R,Num_theta_patch),B_chi(Num_chi,Num_R,Num_theta_patch)
    real(8), intent(in) :: R_chi(Num_chi,Num_R,Num_theta_patch),Gamma_chi(Num_chi,Num_R,Num_theta_patch)
    real(8), intent(in) :: Chi_weight(Num_chi,Num_R,Num_theta_patch)
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    real(8), allocatable :: F_temp(:,:),V_obs_log(:),V_seed_log(:),F_ring(:,:,:),Tau_ring(:,:,:),Tau_prefix(:,:,:)
    real(8), allocatable :: P_syn(:,:),Tau_syn(:,:)
    real(8) :: R_Tobs_chi(Num_R),log_domega_4pi,log_gamma_lo,log_gamma_hi,segment_lo,segment_hi
    real(8) :: cos_tv,sin_tv,theta_lo,theta_hi,theta_center,phi_center,domega
    integer :: last_k2,k_start,lower_bound_real8
    logical :: monotonic_chi

    allocate(F_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))
    allocate(F_ring(Num_nu,Num_chi,Num_R),Tau_ring(Num_nu,Num_chi,Num_R),Tau_prefix(Num_nu,0:Num_chi,Num_R))
    allocate(P_syn(Num_nu,Num_chi),Tau_syn(Num_nu,Num_chi))
    F_tot_obs = zero
    F_temp = zero
    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)
    DL = Boundary(13)
    flux_prefactor = (one+z)/(4d0*pi*DL*DL)
    dtheta = OpeningAngle_jet/Num_theta_patch
    dPhi = two*pi/Num_phi_patch
    V_obs_log = dlog(V_obs)
    V_seed_log = dlog(V_seed)
    cos_tv = dcos(Tv)
    sin_tv = dsin(Tv)

    do I_Theta = 1, Num_theta_patch
        do K2 = 1, Num_R
            call radiation_syn_flux_tau_chi_batch_core(R_front(K2,I_Theta),Num_gam_e,Num_nu,Num_chi,gam_e, &
                                                       DNe_chi(:,:,K2,I_Theta),V_seed,B_chi(:,K2,I_Theta), &
                                                       1.046d4,P_syn,Tau_syn)
            F_ring(:,:,K2) = P_syn*flux_prefactor
            Tau_ring(:,:,K2) = Tau_syn
        end do
        Tau_prefix(:,0,:) = zero
        do I_chi = 1, Num_chi
            Tau_prefix(:,I_chi,:) = Tau_prefix(:,I_chi-1,:) + Tau_ring(:,I_chi,:)
        end do
        theta_lo = dtheta*(I_Theta-1)
        theta_hi = dtheta*I_Theta
        theta_center = dtheta*(I_Theta-0.5d0)
        do i_Phi = 1, Num_phi_patch
            phi_center = (i_Phi-0.5d0)*dPhi
            domega = (dcos(theta_lo)-dcos(theta_hi))*dPhi
            log_domega_4pi = dlog(domega)-dlog(4d0*pi)
            DMu = cos_tv*dcos(theta_center)+sin_tv*dsin(theta_center)*dcos(phi_center)
            do I_chi = 1, Num_chi
                R_Tobs_chi = R_Tobs1(:,I_Theta) + (one+z)*(R_front(:,I_Theta)-R_chi(I_chi,:,I_Theta)*DMu)/Para_c
                monotonic_chi = all(R_Tobs_chi(2:Num_R) > R_Tobs_chi(1:Num_R-1))
                if (monotonic_chi) then
                    II = 1
                    last_k2 = 0
                    do K1 = 1, Num_Tobs
                        if (Tobs(K1) < R_Tobs_chi(1) .or. Tobs(K1) > R_Tobs_chi(Num_R)) cycle
                        do while (II < Num_R-1 .and. Tobs(K1) >= R_Tobs_chi(II+1))
                            II = II + 1
                        end do
                        K2 = II
                        if (Tobs(K1) >= R_Tobs_chi(K2) .and. Tobs(K1) < R_Tobs_chi(K2+1)) then
                            if (K2 /= last_k2) then
                                log_gamma_lo = dlog(Gamma_chi(I_chi,K2,I_Theta))
                                log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1,I_Theta))
                                last_k2 = K2
                            end if
                            Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                            call project_cached_structured_segment(I_chi,K2,K1,Ratio,DMu,log_domega_4pi, &
                                                                   log_gamma_lo,log_gamma_hi)
                        end if
                    end do
                else
                    do K2 = 1, Num_R-1
                        segment_lo = min(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                        segment_hi = max(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                        k_start = lower_bound_real8(Tobs,Num_Tobs,segment_lo)
                        log_gamma_lo = dlog(Gamma_chi(I_chi,K2,I_Theta))
                        log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1,I_Theta))
                        do K1 = k_start, Num_Tobs
                            if (Tobs(K1) >= segment_hi) exit
                            Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                            call project_cached_structured_segment(I_chi,K2,K1,Ratio,DMu,log_domega_4pi, &
                                                                   log_gamma_lo,log_gamma_hi)
                        end do
                    end do
                end if
            end do
        end do
    end do
    F_tot_obs = F_temp
    deallocate(F_temp,V_obs_log,V_seed_log,F_ring,Tau_ring,Tau_prefix,P_syn,Tau_syn)
    return

contains

subroutine project_cached_structured_segment(I_chi,K2,K1,Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_chi,K2,K1
    real(8), intent(in) :: Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi
    real(8) :: F_theta(Num_nu),log_doppler_redshift,log_flux_weight,doppler

    call compute_cached_structured_state(I_chi,K2,Ratio,DMu,log_gamma_lo,log_gamma_hi, &
                                         log_doppler_redshift,doppler,F_theta)
    log_flux_weight = log_domega_4pi - 3d0*dlog(doppler)
    call interpolation_accumulate_shifted_linear_sed(V_seed_log,F_theta,Num_nu,V_obs_log,Num_nu_obs, &
                                                     log_doppler_redshift,log_flux_weight,F_temp(:,K1))
end subroutine project_cached_structured_segment

subroutine compute_cached_structured_state(I_chi,K2,Ratio,DMu,log_gamma_lo,log_gamma_hi, &
                                           log_doppler_redshift,doppler,F_theta)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_chi,K2
    real(8), intent(in) :: Ratio,DMu,log_gamma_lo,log_gamma_hi
    real(8), intent(out) :: log_doppler_redshift,doppler,F_theta(Num_nu)
    real(8) :: DG,Beta

    DG = dexp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
    Beta = dsqrt(one-DG**(-2))
    doppler = DG*(one-Beta*DMu)
    log_doppler_redshift = dlog(doppler)+dlog(one+z)
    call accumulate_cached_structured_source(I_chi,K2,Ratio,F_theta)
end subroutine compute_cached_structured_state

subroutine accumulate_cached_structured_source(I_chi,K2,Ratio,F_theta)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_chi,K2
    real(8), intent(in) :: Ratio
    real(8), intent(out) :: F_theta(Num_nu)
    real(8) :: source_lo,source_hi,tau_front_lo,tau_front_hi,tau_cell_lo,tau_cell_hi
    integer :: I_nu

    F_theta = zero
    do I_nu = 1, Num_nu
        tau_front_lo = Tau_prefix(I_nu,I_chi-1,K2)
        tau_front_hi = Tau_prefix(I_nu,I_chi-1,K2+1)
        tau_cell_lo = Tau_ring(I_nu,I_chi,K2)
        tau_cell_hi = Tau_ring(I_nu,I_chi,K2+1)
        source_lo = F_ring(I_nu,I_chi,K2)*Chi_weight(I_chi,K2,I_Theta) &
                    *chi_ssa_cell_escape(tau_front_lo,tau_cell_lo)
        source_hi = F_ring(I_nu,I_chi,K2+1)*Chi_weight(I_chi,K2+1,I_Theta) &
                    *chi_ssa_cell_escape(tau_front_hi,tau_cell_hi)
        F_theta(I_nu) = (one-Ratio)*source_lo + Ratio*source_hi
    end do
end subroutine accumulate_cached_structured_source
end subroutine sed_interpolation_chi_structured_axisym_electron_cached

! 结构化direct-electronχ投影的单θ-ring版本：用于Python流式solve/project/accumulate。
subroutine sed_interpolation_chi_structured_axisym_electron_cached_ring(Boundary,R_Tobs1,R_front,DNe_chi, &
                             B_chi,R_chi,Gamma_chi,Chi_weight,gam_e,V_seed,V_obs,Tobs,theta_lo,theta_hi, &
                             n,Num_gam_e,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R,F_tot_obs)
    use constants
    use interpolation_common
    use radiation_common, only: radiation_syn_flux_tau_chi_batch_core
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_gam_e,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R
    real(8), intent(in) :: Boundary(n),R_Tobs1(Num_R),R_front(Num_R),Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: gam_e(Num_gam_e),DNe_chi(Num_gam_e,Num_chi,Num_R),B_chi(Num_chi,Num_R)
    real(8), intent(in) :: R_chi(Num_chi,Num_R),Gamma_chi(Num_chi,Num_R),Chi_weight(Num_chi,Num_R)
    real(8), intent(in) :: theta_lo,theta_hi
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    real(8), allocatable :: F_temp(:,:),V_obs_log(:),V_seed_log(:),F_ring(:,:,:),Tau_ring(:,:,:),Tau_prefix(:,:,:)
    real(8), allocatable :: P_syn(:,:),Tau_syn(:,:)
    real(8) :: R_Tobs_chi(Num_R),log_domega_4pi,log_gamma_lo,log_gamma_hi,segment_lo,segment_hi
    real(8) :: cos_tv,sin_tv,theta_center,phi_center,domega,dPhi
    integer :: last_k2,k_start,lower_bound_real8
    logical :: monotonic_chi

    allocate(F_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))
    allocate(F_ring(Num_nu,Num_chi,Num_R),Tau_ring(Num_nu,Num_chi,Num_R),Tau_prefix(Num_nu,0:Num_chi,Num_R))
    allocate(P_syn(Num_nu,Num_chi),Tau_syn(Num_nu,Num_chi))
    F_tot_obs = zero
    F_temp = zero
    z = Boundary(8)
    Tv = Boundary(10)
    DL = Boundary(13)
    flux_prefactor = (one+z)/(4d0*pi*DL*DL)
    dPhi = two*pi/Num_phi_patch
    V_obs_log = dlog(V_obs)
    V_seed_log = dlog(V_seed)
    cos_tv = dcos(Tv)
    sin_tv = dsin(Tv)

    do K2 = 1, Num_R
        call radiation_syn_flux_tau_chi_batch_core(R_front(K2),Num_gam_e,Num_nu,Num_chi,gam_e, &
                                                   DNe_chi(:,:,K2),V_seed,B_chi(:,K2),1.046d4,P_syn,Tau_syn)
        F_ring(:,:,K2) = P_syn*flux_prefactor
        Tau_ring(:,:,K2) = Tau_syn
    end do
    Tau_prefix(:,0,:) = zero
    do I_chi = 1, Num_chi
        Tau_prefix(:,I_chi,:) = Tau_prefix(:,I_chi-1,:) + Tau_ring(:,I_chi,:)
    end do
    theta_center = 0.5d0*(theta_lo+theta_hi)
    do i_Phi = 1, Num_phi_patch
        phi_center = (i_Phi-0.5d0)*dPhi
        domega = (dcos(theta_lo)-dcos(theta_hi))*dPhi
        log_domega_4pi = dlog(domega)-dlog(4d0*pi)
        DMu = cos_tv*dcos(theta_center)+sin_tv*dsin(theta_center)*dcos(phi_center)
        do I_chi = 1, Num_chi
            R_Tobs_chi = R_Tobs1 + (one+z)*(R_front-R_chi(I_chi,:)*DMu)/Para_c
            monotonic_chi = all(R_Tobs_chi(2:Num_R) > R_Tobs_chi(1:Num_R-1))
            if (monotonic_chi) then
                II = 1
                last_k2 = 0
                do K1 = 1, Num_Tobs
                    if (Tobs(K1) < R_Tobs_chi(1) .or. Tobs(K1) > R_Tobs_chi(Num_R)) cycle
                    do while (II < Num_R-1 .and. Tobs(K1) >= R_Tobs_chi(II+1))
                        II = II + 1
                    end do
                    K2 = II
                    if (Tobs(K1) >= R_Tobs_chi(K2) .and. Tobs(K1) < R_Tobs_chi(K2+1)) then
                        if (K2 /= last_k2) then
                            log_gamma_lo = dlog(Gamma_chi(I_chi,K2))
                            log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1))
                            last_k2 = K2
                        end if
                        Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                        call project_cached_ring_segment(I_chi,K2,K1,Ratio,DMu,log_domega_4pi, &
                                                         log_gamma_lo,log_gamma_hi)
                    end if
                end do
            else
                do K2 = 1, Num_R-1
                    segment_lo = min(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                    segment_hi = max(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                    k_start = lower_bound_real8(Tobs,Num_Tobs,segment_lo)
                    log_gamma_lo = dlog(Gamma_chi(I_chi,K2))
                    log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1))
                    do K1 = k_start, Num_Tobs
                        if (Tobs(K1) >= segment_hi) exit
                        Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                        call project_cached_ring_segment(I_chi,K2,K1,Ratio,DMu,log_domega_4pi, &
                                                         log_gamma_lo,log_gamma_hi)
                    end do
                end do
            end if
        end do
    end do
    F_tot_obs = F_temp
    deallocate(F_temp,V_obs_log,V_seed_log,F_ring,Tau_ring,Tau_prefix,P_syn,Tau_syn)
    return

contains

subroutine project_cached_ring_segment(I_chi,K2,K1,Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_chi,K2,K1
    real(8), intent(in) :: Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi
    real(8) :: F_theta(Num_nu),log_doppler_redshift,log_flux_weight,doppler

    call compute_cached_ring_state(I_chi,K2,Ratio,DMu,log_gamma_lo,log_gamma_hi, &
                                   log_doppler_redshift,doppler,F_theta)
    log_flux_weight = log_domega_4pi - 3d0*dlog(doppler)
    call interpolation_accumulate_shifted_linear_sed(V_seed_log,F_theta,Num_nu,V_obs_log,Num_nu_obs, &
                                                     log_doppler_redshift,log_flux_weight,F_temp(:,K1))
end subroutine project_cached_ring_segment

subroutine compute_cached_ring_state(I_chi,K2,Ratio,DMu,log_gamma_lo,log_gamma_hi, &
                                     log_doppler_redshift,doppler,F_theta)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_chi,K2
    real(8), intent(in) :: Ratio,DMu,log_gamma_lo,log_gamma_hi
    real(8), intent(out) :: log_doppler_redshift,doppler,F_theta(Num_nu)
    real(8) :: DG,Beta

    DG = dexp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
    Beta = dsqrt(one-DG**(-2))
    doppler = DG*(one-Beta*DMu)
    log_doppler_redshift = dlog(doppler)+dlog(one+z)
    call accumulate_cached_ring_source(I_chi,K2,Ratio,F_theta)
end subroutine compute_cached_ring_state

subroutine accumulate_cached_ring_source(I_chi,K2,Ratio,F_theta)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_chi,K2
    real(8), intent(in) :: Ratio
    real(8), intent(out) :: F_theta(Num_nu)
    real(8) :: source_lo,source_hi,tau_front_lo,tau_front_hi,tau_cell_lo,tau_cell_hi
    integer :: I_nu

    F_theta = zero
    do I_nu = 1, Num_nu
        tau_front_lo = Tau_prefix(I_nu,I_chi-1,K2)
        tau_front_hi = Tau_prefix(I_nu,I_chi-1,K2+1)
        tau_cell_lo = Tau_ring(I_nu,I_chi,K2)
        tau_cell_hi = Tau_ring(I_nu,I_chi,K2+1)
        source_lo = F_ring(I_nu,I_chi,K2)*Chi_weight(I_chi,K2)*chi_ssa_cell_escape(tau_front_lo,tau_cell_lo)
        source_hi = F_ring(I_nu,I_chi,K2+1)*Chi_weight(I_chi,K2+1)*chi_ssa_cell_escape(tau_front_hi,tau_cell_hi)
        F_theta(I_nu) = (one-Ratio)*source_lo + Ratio*source_hi
    end do
end subroutine accumulate_cached_ring_source
end subroutine sed_interpolation_chi_structured_axisym_electron_cached_ring

! 预计算光谱版：F_ring=P_syn*flux_prefactor 和 Tau_ring 已在Python侧从transport输出准备好，
! 跳过 radiation_syn_flux_tau_chi_batch_core 重算。
subroutine sed_interpolation_chi_structured_axisym_ring_precomputed(Boundary,R_Tobs1,R_front, &
                             F_ring,Tau_ring,R_chi,Gamma_chi,Chi_weight,V_seed,V_obs,Tobs,theta_lo,theta_hi, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R,F_tot_obs)
    use constants
    use interpolation_common
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_phi_patch,Num_chi,Num_R
    real(8), intent(in) :: Boundary(n),R_Tobs1(Num_R),R_front(Num_R),Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: F_ring(Num_nu,Num_chi,Num_R),Tau_ring(Num_nu,Num_chi,Num_R)
    real(8), intent(in) :: R_chi(Num_chi,Num_R),Gamma_chi(Num_chi,Num_R),Chi_weight(Num_chi,Num_R)
    real(8), intent(in) :: theta_lo,theta_hi
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    real(8), allocatable :: F_temp(:,:),V_obs_log(:),V_seed_log(:),Tau_prefix(:,:,:)
    real(8) :: R_Tobs_chi(Num_R),log_domega_4pi,log_gamma_lo,log_gamma_hi,segment_lo,segment_hi
    real(8) :: cos_tv,sin_tv,theta_center,phi_center,domega,dPhi
    integer :: last_k2,k_start,lower_bound_real8
    logical :: monotonic_chi

    allocate(F_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))
    allocate(Tau_prefix(Num_nu,0:Num_chi,Num_R))
    F_tot_obs = zero
    F_temp = zero
    z = Boundary(8)
    Tv = Boundary(10)
    dPhi = two*pi/Num_phi_patch
    V_obs_log = dlog(V_obs)
    V_seed_log = dlog(V_seed)
    cos_tv = dcos(Tv)
    sin_tv = dsin(Tv)

    Tau_prefix(:,0,:) = zero
    do I_chi = 1, Num_chi
        Tau_prefix(:,I_chi,:) = Tau_prefix(:,I_chi-1,:) + Tau_ring(:,I_chi,:)
    end do
    theta_center = 0.5d0*(theta_lo+theta_hi)
    do i_Phi = 1, Num_phi_patch
        phi_center = (i_Phi-0.5d0)*dPhi
        domega = (dcos(theta_lo)-dcos(theta_hi))*dPhi
        log_domega_4pi = dlog(domega)-dlog(4d0*pi)
        DMu = cos_tv*dcos(theta_center)+sin_tv*dsin(theta_center)*dcos(phi_center)
        do I_chi = 1, Num_chi
            R_Tobs_chi = R_Tobs1 + (one+z)*(R_front-R_chi(I_chi,:)*DMu)/Para_c
            monotonic_chi = all(R_Tobs_chi(2:Num_R) > R_Tobs_chi(1:Num_R-1))
            if (monotonic_chi) then
                II = 1
                last_k2 = 0
                do K1 = 1, Num_Tobs
                    if (Tobs(K1) < R_Tobs_chi(1) .or. Tobs(K1) > R_Tobs_chi(Num_R)) cycle
                    do while (II < Num_R-1 .and. Tobs(K1) >= R_Tobs_chi(II+1))
                        II = II + 1
                    end do
                    K2 = II
                    if (Tobs(K1) >= R_Tobs_chi(K2) .and. Tobs(K1) < R_Tobs_chi(K2+1)) then
                        if (K2 /= last_k2) then
                            log_gamma_lo = dlog(Gamma_chi(I_chi,K2))
                            log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1))
                            last_k2 = K2
                        end if
                        Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                        call project_precomputed_ring_segment(I_chi,K2,K1,Ratio,DMu,log_domega_4pi, &
                                                              log_gamma_lo,log_gamma_hi)
                    end if
                end do
            else
                do K2 = 1, Num_R-1
                    segment_lo = min(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                    segment_hi = max(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                    k_start = lower_bound_real8(Tobs,Num_Tobs,segment_lo)
                    log_gamma_lo = dlog(Gamma_chi(I_chi,K2))
                    log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1))
                    do K1 = k_start, Num_Tobs
                        if (Tobs(K1) >= segment_hi) exit
                        Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                        call project_precomputed_ring_segment(I_chi,K2,K1,Ratio,DMu,log_domega_4pi, &
                                                              log_gamma_lo,log_gamma_hi)
                    end do
                end do
            end if
        end do
    end do
    F_tot_obs = F_temp
    deallocate(F_temp,V_obs_log,V_seed_log,Tau_prefix)
    return

contains

subroutine project_precomputed_ring_segment(I_chi,K2,K1,Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_chi,K2,K1
    real(8), intent(in) :: Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi
    real(8) :: F_theta(Num_nu),log_doppler_redshift,log_flux_weight,doppler,DG,Beta

    DG = dexp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
    Beta = dsqrt(one-DG**(-2))
    doppler = DG*(one-Beta*DMu)
    log_doppler_redshift = dlog(doppler)+dlog(one+z)
    call accumulate_precomputed_ring_source(I_chi,K2,Ratio,F_theta)
    log_flux_weight = log_domega_4pi - 3d0*dlog(doppler)
    call interpolation_accumulate_shifted_linear_sed(V_seed_log,F_theta,Num_nu,V_obs_log,Num_nu_obs, &
                                                     log_doppler_redshift,log_flux_weight,F_temp(:,K1))
end subroutine project_precomputed_ring_segment

subroutine accumulate_precomputed_ring_source(I_chi,K2,Ratio,F_theta)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: I_chi,K2
    real(8), intent(in) :: Ratio
    real(8), intent(out) :: F_theta(Num_nu)
    real(8) :: source_lo,source_hi,tau_front_lo,tau_front_hi,tau_cell_lo,tau_cell_hi
    integer :: I_nu

    F_theta = zero
    do I_nu = 1, Num_nu
        tau_front_lo = Tau_prefix(I_nu,I_chi-1,K2)
        tau_front_hi = Tau_prefix(I_nu,I_chi-1,K2+1)
        tau_cell_lo = Tau_ring(I_nu,I_chi,K2)
        tau_cell_hi = Tau_ring(I_nu,I_chi,K2+1)
        source_lo = F_ring(I_nu,I_chi,K2)*Chi_weight(I_chi,K2)*chi_ssa_cell_escape(tau_front_lo,tau_cell_lo)
        source_hi = F_ring(I_nu,I_chi,K2+1)*Chi_weight(I_chi,K2+1)*chi_ssa_cell_escape(tau_front_hi,tau_cell_hi)
        F_theta(I_nu) = (one-Ratio)*source_lo + Ratio*source_hi
    end do
end subroutine accumulate_precomputed_ring_source
end subroutine sed_interpolation_chi_structured_axisym_ring_precomputed

! 单频同步辐射点计算：与radiation_syn_seed_core同一核，用于direct chi投影。
subroutine chi_synch_point(R_loc,DB,Num_gam_e,gam_e,dN_gam_e,V_cal,P_syn,Tau_syn)
    use constants
    use radiation_common, only: radiation_syn_kernel_value, radiation_transfer_factor
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_cal
    real(8), intent(out) :: P_syn,Tau_syn
    real(8) :: factor,Temp_syn,Rariv2,dInteg,Tau,x,ratio_v_pow,Fx,P_v,transfer

    factor=(3.62d0/pi)**2
    Temp_syn=dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    Rariv2=R_loc*R_loc
    dInteg=zero
    Tau=zero
    do I_gam_e=1,Num_gam_e-1
        gam_mid2=(gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0
        dN_seg=(dN_gam_e(I_gam_e)+dN_gam_e(I_gam_e+1))*(gam_e(I_gam_e+1)-gam_e(I_gam_e))/two
        ddN=dN_gam_e(I_gam_e)/(gam_e(I_gam_e)*gam_e(I_gam_e)) - &
            dN_gam_e(I_gam_e+1)/(gam_e(I_gam_e+1)*gam_e(I_gam_e+1))
        Vc=(4.2d6)*gam_mid2*DB
        x=V_cal/Vc
        ratio_v_pow=(Vc/V_cal)**(2d0/3d0)
        Fx=radiation_syn_kernel_value(x,ratio_v_pow,factor)
        dInteg=dInteg+dN_seg*Fx
        Tau=Tau+gam_mid2*ddN*Fx
    end do
    P_v=Temp_syn*DB*dInteg
    Tau=1.046d4*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal)
    Tau_syn=Tau
    call radiation_transfer_factor(Tau,transfer)
    P_syn=P_v*transfer
end subroutine chi_synch_point
