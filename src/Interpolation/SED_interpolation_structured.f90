! Created by rj on 2021/2/4.

!##################################################################################################
!When an intrinsic SED 'P_tot' observed in the comoving frame send in,
!to produce the observed SED 'P_tot_obs' after consider the EATS and Dopper boosting effect.
!##################################################################################################

! SEDEATS++
subroutine sed_interpolation_structured(Boundary, angle_narrow_jet, R_Tobs1,R_gamma,R,P_tot,V_seed,V_obs,Tobs, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads, P_tot_obs)
    !$ use omp_lib
    use constants
    use interpolation_common
    implicit none
    !##############################################################################################
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in) :: angle_narrow_jet
    real(8), intent(in), dimension(Num_Tobs) :: Tobs
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu_obs) :: V_obs
    real(8), intent(in), dimension(Num_R,Num_theta) :: R_Tobs1,R_gamma
    real(8), intent(in), dimension(Num_R,Num_theta) :: R
    real(8), intent(in), dimension(Num_nu,Num_R,Num_theta) :: P_tot
    real(8), intent(out), dimension(Num_nu_obs,Num_Tobs) :: P_tot_obs


    real(8), allocatable, dimension(:,:) :: P_tot_obs_temp
    real(8), allocatable, dimension(:) :: V_obs_log,V_seed_log
    real(8), dimension(Num_R) :: R_Tobs_theta
    real(8), dimension(Num_nu) :: specLo,specHi
    real(8), dimension(Num_Tobs) :: T_sorted
    integer, dimension(Num_Tobs) :: T_order
    real(8) :: z,OpeningAngle_jet,Tv,dPhi,phi_scale,dtheta,Taa_lower,Taa_boundary,Taa_center,domega
    real(8) :: Phi_center,DMu,Ratio,log_gamma_lo,log_gamma_hi,log_domega_4pi
    integer :: I_Theta,i_Phi,Iobs,K2,II,last_k2
    logical :: time_ordered
    allocate (P_tot_obs_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))


    P_tot_obs_temp=0d0
    P_tot_obs=0d0

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

    !$OMP PARALLEL num_threads(n_threads), reduction(+:P_tot_obs_temp), private(I_Theta, &
    !$OMP& Taa_lower, Taa_boundary, Taa_center, domega, i_Phi, Phi_center, DMu, &
    !$OMP& Iobs, K2, II, Ratio, R_Tobs_theta, &
    !$OMP& specLo, specHi, last_k2, log_gamma_lo, log_gamma_hi, log_domega_4pi)
    !$OMP DO SCHEDULE(GUIDED,4)
    do I_Theta=1,Num_Theta
       Taa_lower=dtheta*(I_Theta-1)
       Taa_boundary=dtheta*I_Theta
       Taa_center=dtheta*(I_Theta-0.5d0)
       if (angle_narrow_jet>Taa_center) cycle
       domega=(dcos(Taa_lower)-dcos(Taa_boundary))*dPhi
       log_domega_4pi=log(domega)-log(4.0d0*pi)
       do i_Phi=1,Num_Phi
          Phi_center=(i_Phi-0.5)*dPhi
          DMu=dcos(Tv)*dcos(Taa_center)+dsin(Tv)*dsin(Taa_center)*dcos(Phi_center)
          R_Tobs_theta=R_Tobs1(:,I_Theta)+R(:,I_Theta)*(1d0-DMu)*(1d0+z)/Para_c
          II=1
          last_k2=0
          do Iobs=1,Num_Tobs
             if (T_sorted(Iobs) < R_Tobs_theta(1) .or. T_sorted(Iobs) >= R_Tobs_theta(Num_R)) cycle
             do while (II < Num_R-1 .and. T_sorted(Iobs) >= R_Tobs_theta(II+1))
                II=II+1
             end do
             K2=II
             if (K2 /= last_k2) then
                 log_gamma_lo=log(R_gamma(K2,I_Theta))
                 log_gamma_hi=log(R_gamma(K2+1,I_Theta))
                 specLo=P_tot(:,K2,I_Theta)
                 specHi=P_tot(:,K2+1,I_Theta)
                 last_k2=K2
             end if
             Ratio=(T_sorted(Iobs)-R_Tobs_theta(K2))/(R_Tobs_theta(K2+1)-R_Tobs_theta(K2))
             call project_structured_segment(Ratio,DMu,log_domega_4pi, &
                                             log_gamma_lo,log_gamma_hi,specLo,specHi,P_tot_obs_temp(:,Iobs))
         end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    if (time_ordered) then
        P_tot_obs=P_tot_obs_temp*phi_scale
    else
        do Iobs=1,Num_Tobs
            P_tot_obs(:,T_order(Iobs))=P_tot_obs_temp(:,Iobs)*phi_scale
        end do
    end if

    deallocate (P_tot_obs_temp,V_obs_log,V_seed_log)


    return

contains

subroutine project_structured_segment(Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi, &
                                      specLo,specHi,flux_accum)
    implicit none
    real(8), intent(in) :: Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi
    real(8), intent(in), dimension(Num_nu) :: specLo,specHi
    real(8), intent(inout), dimension(Num_nu_obs) :: flux_accum
    real(8), dimension(Num_nu) :: specTheta,logTheta,vShift
    real(8) :: DG,Beta,doppler,log_doppler_redshift

    DG=exp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
    specTheta=(1d0-Ratio)*specLo+Ratio*specHi
    logTheta=-huge(1d0)
    where(specTheta > 0d0) logTheta=log(specTheta)
    Beta=dsqrt(1d0-DG**(-2))
    doppler=DG*(1d0-Beta*DMu)
    log_doppler_redshift=log(doppler)+log(1d0+z)
    logTheta=logTheta+log_domega_4pi-3d0*log(doppler)
    vShift=V_seed_log-log_doppler_redshift
    call accum_logsed(vShift,logTheta, &
                                          Num_nu,V_obs_log,Num_nu_obs,flux_accum)
end subroutine project_structured_segment
end subroutine sed_interpolation_structured

!  theta-phi
subroutine sed_structured_phi(Boundary,R_Tobs1,R_gamma,R,P_tot,V_seed,V_obs,Tobs, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_R,n_threads, P_tot_obs)
    !$ use omp_lib
    use constants
    use interpolation_common
    implicit none
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_R,n_threads
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_Tobs) :: Tobs
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu_obs) :: V_obs
    real(8), intent(in), dimension(Num_R,Num_Theta,Num_Phi) :: R_Tobs1,R_gamma
    real(8), intent(in), dimension(Num_R,Num_Theta,Num_Phi) :: R
    real(8), intent(in), dimension(Num_nu,Num_R,Num_Theta,Num_Phi) :: P_tot
    real(8), intent(out), dimension(Num_nu_obs,Num_Tobs) :: P_tot_obs

    real(8), allocatable, dimension(:,:) :: P_tot_obs_temp
    real(8), allocatable, dimension(:) :: V_obs_log,V_seed_log
    real(8), dimension(Num_R) :: R_Tobs_theta
    real(8), dimension(Num_nu) :: specLo,specHi
    real(8), dimension(Num_Tobs) :: T_sorted
    integer, dimension(Num_Tobs) :: T_order
    real(8) :: z,OpeningAngle_jet,Tv,dtheta,dPhi,Taa_lower,Taa_boundary,Taa_center,domega
    real(8) :: Phi_center,DMu,Ratio,log_gamma_lo,log_gamma_hi,log_domega_4pi
    integer :: I_Theta,i_Phi,Iobs,K2,II,last_k2
    logical :: time_ordered
    allocate (P_tot_obs_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))

    P_tot_obs_temp=0d0
    P_tot_obs=0d0
    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)
    dtheta=OpeningAngle_jet/Num_Theta
    dPhi=2d0*pi/Num_Phi
    V_obs_log = log(V_obs)
    V_seed_log = log(V_seed)
    call time_order(Tobs,Num_Tobs,T_order,T_sorted,time_ordered)

    !$OMP PARALLEL num_threads(n_threads), reduction(+:P_tot_obs_temp), private(I_Theta, &
    !$OMP& Taa_lower,Taa_boundary,Taa_center,domega,i_Phi,Phi_center,DMu, &
    !$OMP& Iobs,K2,II,Ratio,R_Tobs_theta,specLo,specHi, &
    !$OMP& last_k2,log_gamma_lo,log_gamma_hi,log_domega_4pi)
    !$OMP DO COLLAPSE(2) SCHEDULE(GUIDED,4)
    do I_Theta=1,Num_Theta
       do i_Phi=1,Num_Phi
          Taa_lower=dtheta*(I_Theta-1)
          Taa_boundary=dtheta*I_Theta
          Taa_center=dtheta*(I_Theta-0.5d0)
          Phi_center=(i_Phi-0.5d0)*dPhi
          domega=(dcos(Taa_lower)-dcos(Taa_boundary))*dPhi
          log_domega_4pi=log(domega)-log(4.0d0*pi)
          DMu=dcos(Tv)*dcos(Taa_center)+dsin(Tv)*dsin(Taa_center)*dcos(Phi_center)
          R_Tobs_theta=R_Tobs1(:,I_Theta,i_Phi)+R(:,I_Theta,i_Phi)*(1d0-DMu)*(1d0+z)/Para_c
          II=1
          last_k2=0
          do Iobs=1,Num_Tobs
             if (T_sorted(Iobs) < R_Tobs_theta(1) .or. T_sorted(Iobs) >= R_Tobs_theta(Num_R)) cycle
             do while (II < Num_R-1 .and. T_sorted(Iobs) >= R_Tobs_theta(II+1))
                II=II+1
             end do
             K2=II
             if (K2 /= last_k2) then
                 log_gamma_lo=log(R_gamma(K2,I_Theta,i_Phi))
                 log_gamma_hi=log(R_gamma(K2+1,I_Theta,i_Phi))
                 specLo=P_tot(:,K2,I_Theta,i_Phi)
                 specHi=P_tot(:,K2+1,I_Theta,i_Phi)
                 last_k2=K2
             end if
             Ratio=(T_sorted(Iobs)-R_Tobs_theta(K2))/(R_Tobs_theta(K2+1)-R_Tobs_theta(K2))
             call project_phi(Ratio,DMu,log_domega_4pi, &
                                                 log_gamma_lo,log_gamma_hi,specLo,specHi,P_tot_obs_temp(:,Iobs))
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    if (time_ordered) then
        P_tot_obs=P_tot_obs_temp
    else
        do Iobs=1,Num_Tobs
            P_tot_obs(:,T_order(Iobs))=P_tot_obs_temp(:,Iobs)
        end do
    end if
    deallocate (P_tot_obs_temp,V_obs_log,V_seed_log)
    return

contains

subroutine project_phi(Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi, &
                                          specLo,specHi,flux_accum)
    implicit none
    real(8), intent(in) :: Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi
    real(8), intent(in), dimension(Num_nu) :: specLo,specHi
    real(8), intent(inout), dimension(Num_nu_obs) :: flux_accum
    real(8), dimension(Num_nu) :: specTheta,logTheta,vShift
    real(8) :: DG,Beta,doppler,log_doppler_redshift

    DG=exp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
    specTheta=(1d0-Ratio)*specLo+Ratio*specHi
    logTheta=-huge(1d0)
    where(specTheta > 0d0) logTheta=log(specTheta)
    Beta=dsqrt(1d0-DG**(-2))
    doppler=DG*(1d0-Beta*DMu)
    log_doppler_redshift=log(doppler)+log(1d0+z)
    logTheta=logTheta+log_domega_4pi-3d0*log(doppler)
    vShift=V_seed_log-log_doppler_redshift
    call accum_logsed(vShift,logTheta, &
                                          Num_nu,V_obs_log,Num_nu_obs,flux_accum)
end subroutine project_phi
end subroutine sed_structured_phi
