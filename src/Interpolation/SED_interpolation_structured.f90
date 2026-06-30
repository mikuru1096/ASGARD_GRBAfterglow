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
    real(8), intent(in) :: Boundary(n),angle_narrow_jet
    real(8), intent(in) :: Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: R_Tobs1(Num_R,Num_theta),R_gamma(Num_R,Num_theta)
    real(8), intent(in) :: R(Num_R,Num_theta),P_tot(Num_nu,Num_R,Num_theta)
    real(8), intent(out) :: P_tot_obs(Num_nu_obs,Num_Tobs)
    
    
    real(8), allocatable :: P_tot_obs_temp(:,:),V_obs_log(:),V_seed_log(:)
    real(8) :: R_Tobs_theta(Num_R),P_tot_log_lo(Num_nu),P_tot_log_hi(Num_nu)
    real(8) :: z,OpeningAngle_jet,Tv,dPhi,phi_scale,dtheta,Taa_lower,Taa_boundary,Taa_center,domega
    real(8) :: Phi_center,DMu,Ratio,log_gamma_lo,log_gamma_hi,log_domega_4pi
    integer :: I_Theta,i_Phi,K1,K2,II,last_k2
    allocate (P_tot_obs_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))
    
    
    P_tot_obs_temp=zero
    P_tot_obs=zero
    
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
    
    !$OMP PARALLEL num_threads(n_threads), reduction(+:P_tot_obs_temp), private(I_Theta, &
    !$OMP& Taa_lower, Taa_boundary, Taa_center, domega, i_Phi, Phi_center, DMu, &
    !$OMP& K1, II, K2, Ratio, R_Tobs_theta, &
    !$OMP& P_tot_log_lo, P_tot_log_hi, last_k2, log_gamma_lo, log_gamma_hi, log_domega_4pi)
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
          R_Tobs_theta=R_Tobs1(:,I_Theta)+R(:,I_Theta)*(one-DMu)*(one+z)/Para_c
          II=1
          last_k2=0
          do K1=1,Num_Tobs
             if (Tobs(K1) < R_Tobs_theta(1) .or. Tobs(K1) > R_Tobs_theta(Num_R)) cycle
             do while (II < Num_R-1 .and. Tobs(K1) >= R_Tobs_theta(II+1))
                II=II+1
             end do
             K2=II
             if ((Tobs(K1) >= R_Tobs_theta(K2)).and.(Tobs(K1) < R_Tobs_theta(K2+1))) then
                 if (K2 /= last_k2) then
                     log_gamma_lo=log(R_gamma(K2,I_Theta))
                     log_gamma_hi=log(R_gamma(K2+1,I_Theta))
                     P_tot_log_lo=-huge(one)
                     P_tot_log_hi=-huge(one)
                     where(P_tot(:,K2,I_Theta) > zero) P_tot_log_lo=log(P_tot(:,K2,I_Theta))
                     where(P_tot(:,K2+1,I_Theta) > zero) P_tot_log_hi=log(P_tot(:,K2+1,I_Theta))
                     last_k2=K2
                 end if
                 Ratio=(Tobs(K1)-R_Tobs_theta(K2))/(R_Tobs_theta(K2+1)-R_Tobs_theta(K2))
                 call project_structured_segment(K1,Ratio,DMu,log_domega_4pi, &
                                                 log_gamma_lo,log_gamma_hi,P_tot_log_lo,P_tot_log_hi)
             end if
         end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    P_tot_obs=P_tot_obs_temp*phi_scale
    
    deallocate (P_tot_obs_temp,V_obs_log,V_seed_log)
    
    
    return

contains

subroutine project_structured_segment(K1,Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi, &
                                      P_tot_log_lo,P_tot_log_hi)
    implicit none
    integer, intent(in) :: K1
    real(8), intent(in) :: Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi
    real(8), intent(in) :: P_tot_log_lo(Num_nu),P_tot_log_hi(Num_nu)
    real(8) :: DP_theta(Num_nu),P_tot_log_theta(Num_nu),V_seed_log_theta(Num_nu)
    real(8) :: DG,Beta,doppler,log_doppler_redshift

    DG=exp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
    DP_theta=P_tot_log_lo+Ratio*(P_tot_log_hi-P_tot_log_lo)
    Beta=dsqrt(one-DG**(-2))
    doppler=DG*(one-Beta*DMu)
    log_doppler_redshift=log(doppler)+log(one+z)
    P_tot_log_theta=max(-199d0,DP_theta+log_domega_4pi-3d0*log(doppler))
    V_seed_log_theta=V_seed_log-log_doppler_redshift
    call interpolation_accumulate_log_sed(V_seed_log_theta,P_tot_log_theta, &
                                          Num_nu,V_obs_log,Num_nu_obs,P_tot_obs_temp(:,K1))
end subroutine project_structured_segment
end subroutine sed_interpolation_structured

!  theta-phi 
subroutine sed_interpolation_structured_phi(Boundary,R_Tobs1,R_gamma,R,P_tot,V_seed,V_obs,Tobs, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_R,n_threads, P_tot_obs)
    !$ use omp_lib
    use constants
    use interpolation_common
    implicit none
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_Phi,Num_R,n_threads
    real(8), intent(in) :: Boundary(n),Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: R_Tobs1(Num_R,Num_Theta,Num_Phi),R_gamma(Num_R,Num_Theta,Num_Phi)
    real(8), intent(in) :: R(Num_R,Num_Theta,Num_Phi),P_tot(Num_nu,Num_R,Num_Theta,Num_Phi)
    real(8), intent(out) :: P_tot_obs(Num_nu_obs,Num_Tobs)

    real(8), allocatable :: P_tot_obs_temp(:,:),V_obs_log(:),V_seed_log(:)
    real(8) :: R_Tobs_theta(Num_R),P_tot_log_lo(Num_nu),P_tot_log_hi(Num_nu)
    real(8) :: z,OpeningAngle_jet,Tv,dtheta,dPhi,Taa_lower,Taa_boundary,Taa_center,domega
    real(8) :: Phi_center,DMu,Ratio,log_gamma_lo,log_gamma_hi,log_domega_4pi
    integer :: I_Theta,i_Phi,K1,K2,II,last_k2
    allocate (P_tot_obs_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))

    P_tot_obs_temp=zero
    P_tot_obs=zero
    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)
    dtheta=OpeningAngle_jet/Num_Theta
    dPhi=two*pi/Num_Phi
    V_obs_log = log(V_obs)
    V_seed_log = log(V_seed)

    !$OMP PARALLEL num_threads(n_threads), reduction(+:P_tot_obs_temp), private(I_Theta, &
    !$OMP& Taa_lower,Taa_boundary,Taa_center,domega,i_Phi,Phi_center,DMu, &
    !$OMP& K1,II,K2,Ratio,R_Tobs_theta,P_tot_log_lo,P_tot_log_hi, &
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
          R_Tobs_theta=R_Tobs1(:,I_Theta,i_Phi)+R(:,I_Theta,i_Phi)*(one-DMu)*(one+z)/Para_c
          II=1
          last_k2=0
          do K1=1,Num_Tobs
             if (Tobs(K1) < R_Tobs_theta(1) .or. Tobs(K1) > R_Tobs_theta(Num_R)) cycle
             do while (II < Num_R-1 .and. Tobs(K1) >= R_Tobs_theta(II+1))
                II=II+1
             end do
             K2=II
             if ((Tobs(K1) >= R_Tobs_theta(K2)).and.(Tobs(K1) < R_Tobs_theta(K2+1))) then
                 if (K2 /= last_k2) then
                     log_gamma_lo=log(R_gamma(K2,I_Theta,i_Phi))
                     log_gamma_hi=log(R_gamma(K2+1,I_Theta,i_Phi))
                     P_tot_log_lo=-huge(one)
                     P_tot_log_hi=-huge(one)
                     where(P_tot(:,K2,I_Theta,i_Phi) > zero) P_tot_log_lo=log(P_tot(:,K2,I_Theta,i_Phi))
                     where(P_tot(:,K2+1,I_Theta,i_Phi) > zero) P_tot_log_hi=log(P_tot(:,K2+1,I_Theta,i_Phi))
                     last_k2=K2
                 end if
                 Ratio=(Tobs(K1)-R_Tobs_theta(K2))/(R_Tobs_theta(K2+1)-R_Tobs_theta(K2))
                 call project_structured_phi_segment(K1,Ratio,DMu,log_domega_4pi, &
                                                     log_gamma_lo,log_gamma_hi,P_tot_log_lo,P_tot_log_hi)
             end if
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    P_tot_obs=P_tot_obs_temp
    deallocate (P_tot_obs_temp,V_obs_log,V_seed_log)
    return

contains

subroutine project_structured_phi_segment(K1,Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi, &
                                          P_tot_log_lo,P_tot_log_hi)
    implicit none
    integer, intent(in) :: K1
    real(8), intent(in) :: Ratio,DMu,log_domega_4pi,log_gamma_lo,log_gamma_hi
    real(8), intent(in) :: P_tot_log_lo(Num_nu),P_tot_log_hi(Num_nu)
    real(8) :: DP_theta(Num_nu),P_tot_log_theta(Num_nu),V_seed_log_theta(Num_nu)
    real(8) :: DG,Beta,doppler,log_doppler_redshift

    DG=exp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
    DP_theta=P_tot_log_lo+Ratio*(P_tot_log_hi-P_tot_log_lo)
    Beta=dsqrt(one-DG**(-2))
    doppler=DG*(one-Beta*DMu)
    log_doppler_redshift=log(doppler)+log(one+z)
    P_tot_log_theta=max(-199d0,DP_theta+log_domega_4pi-3d0*log(doppler))
    V_seed_log_theta=V_seed_log-log_doppler_redshift
    call interpolation_accumulate_log_sed(V_seed_log_theta,P_tot_log_theta, &
                                          Num_nu,V_obs_log,Num_nu_obs,P_tot_obs_temp(:,K1))
end subroutine project_structured_phi_segment
end subroutine sed_interpolation_structured_phi
