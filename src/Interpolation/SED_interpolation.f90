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
    real(8) :: R_Tobs_theta(Num_R),DP_theta(Num_nu),F_tot_log_theta(Num_nu), &
               V_seed_log_theta(Num_nu), &
               F_tot_log_lo(Num_nu),F_tot_log_hi(Num_nu),log_gamma_lo,log_gamma_hi,log_domega_4pi, &
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
    !$OMP PARALLEL num_threads(n_threads), reduction(+:F_tot_obs_temp), private(I_Theta, Taa_boundary, &
    !$OMP& Taa_center, domega, i_Phi, Phi_center, DMu, K1, II, K2, Ratio, DG, Beta, doppler, &
    !$OMP& R_Tobs_theta, DP_theta, F_tot_log_theta, V_seed_log_theta, F_tot_log_lo, F_tot_log_hi, &
    !$OMP& last_k2, log_gamma_lo, log_gamma_hi, log_domega_4pi, log_doppler_redshift)
    !$OMP DO SCHEDULE(GUIDED,4)
    do I_Theta=1,Num_Theta
       Taa_boundary=dtheta*I_Theta
       Taa_center=dtheta*(I_Theta-0.5)
       domega=dsin(Taa_boundary)*dtheta*dPhi
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
                     F_tot_log_lo=log(F_tot(:,K2))
                     F_tot_log_hi=log(F_tot(:,K2+1))
                     last_k2=K2
                 end if
                 Ratio=(Tobs(K1)-R_Tobs_theta(K2))/(R_Tobs_theta(K2+1)-R_Tobs_theta(K2))
                 DG=exp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
                 DP_theta=F_tot_log_lo+Ratio*(F_tot_log_hi-F_tot_log_lo)
                 !logarithm interpolation to get the intrinsic SED at EATS.
                 Beta=dsqrt(one-DG**(-2))

                 doppler=DG*(one-Beta*DMu) !Doppler factor, changed with R
                 log_doppler_redshift=log(doppler)+log(one+z)
                 F_tot_log_theta=max(-199d0,DP_theta+log_domega_4pi-3d0*log(doppler))
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

! 输出二维随机冲击面磁场的角元级Stokes：同一EATS/Doppler角元同时累加I和局域Q。
subroutine sed_interpolation_shock_random_stokes(Boundary,R_Tobs1,R_gamma,R,F_tot,P_tot,V_seed,V_obs,Tobs, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads, I_obs,Q_obs,U_obs)
    !$ use omp_lib
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads
    real(8), intent(in) :: Boundary(n),Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: R_Tobs1(Num_R),R_gamma(Num_R),R(Num_R),F_tot(Num_nu,Num_R),P_tot(Num_nu,Num_R)
    real(8), intent(out) :: I_obs(Num_nu_obs,Num_Tobs),Q_obs(Num_nu_obs,Num_Tobs),U_obs(Num_nu_obs,Num_Tobs)

    real(8), allocatable :: I_obs_temp(:,:),Q_obs_temp(:,:),U_obs_temp(:,:),V_obs_log(:),V_seed_log(:)
    real(8) :: R_Tobs_theta(Num_R),F_log_theta(Num_nu),P_log_theta(Num_nu),V_seed_log_theta(Num_nu), &
               F_log_lo(Num_nu),F_log_hi(Num_nu),P_log_lo(Num_nu),P_log_hi(Num_nu), &
               log_gamma_lo,log_gamma_hi,log_domega_4pi,log_doppler_redshift
    real(8) :: beta,mu_prime,anisotropy,proj_x,proj_y,proj_norm2,q_angle,stokes_weight,ratio_freq
    integer :: last_k2,i_nu_obs,i_nu_src
    allocate (I_obs_temp(Num_nu_obs,Num_Tobs),Q_obs_temp(Num_nu_obs,Num_Tobs), &
              U_obs_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu_obs),V_seed_log(Num_nu))

    I_obs=zero; Q_obs=zero; U_obs=zero
    I_obs_temp=zero; Q_obs_temp=zero; U_obs_temp=zero
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
    !$OMP PARALLEL num_threads(n_threads), reduction(+:I_obs_temp,Q_obs_temp,U_obs_temp), private(I_Theta, &
    !$OMP& Taa_boundary, Taa_center, domega, i_Phi, Phi_center, DMu, K1, II, K2, Ratio, DG, beta, doppler, &
    !$OMP& R_Tobs_theta, F_log_theta, P_log_theta, V_seed_log_theta, F_log_lo, F_log_hi, P_log_lo, P_log_hi, &
    !$OMP& last_k2, log_gamma_lo, log_gamma_hi, log_domega_4pi, log_doppler_redshift, mu_prime, anisotropy, &
    !$OMP& proj_x, proj_y, proj_norm2, q_angle, stokes_weight, i_nu_obs, i_nu_src, ratio_freq)
    !$OMP DO SCHEDULE(GUIDED,4)
    do I_Theta=1,Num_Theta
       Taa_boundary=dtheta*I_Theta
       Taa_center=dtheta*(I_Theta-0.5)
       domega=dsin(Taa_boundary)*dtheta*dPhi
       log_domega_4pi=log(domega)-log(4.0d0*pi)
       do i_Phi=1,Num_Phi
          Phi_center=(i_Phi-0.5)*dPhi
          DMu=dcos(Tv)*dcos(Taa_center)+dsin(Tv)*dsin(Taa_center)*dcos(Phi_center)
          proj_x=dsin(Tv)*dcos(Taa_center)-dcos(Tv)*dsin(Taa_center)*dcos(Phi_center)
          proj_y=dsin(Taa_center)*dsin(Phi_center)
          proj_norm2=one-DMu*DMu
          if (proj_norm2 <= zero) then
              q_angle=zero
          else
              q_angle=(proj_x*proj_x-proj_y*proj_y)/proj_norm2
          end if
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
                     log_gamma_lo=log(R_gamma(K2)); log_gamma_hi=log(R_gamma(K2+1))
                     F_log_lo=log(F_tot(:,K2)); F_log_hi=log(F_tot(:,K2+1))
                     P_log_lo=log(P_tot(:,K2)); P_log_hi=log(P_tot(:,K2+1))
                     last_k2=K2
                 end if
                 Ratio=(Tobs(K1)-R_Tobs_theta(K2))/(R_Tobs_theta(K2+1)-R_Tobs_theta(K2))
                 DG=exp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
                 F_log_theta=F_log_lo+Ratio*(F_log_hi-F_log_lo)
                 P_log_theta=P_log_lo+Ratio*(P_log_hi-P_log_lo)
                 beta=dsqrt(one-DG**(-2))
                 doppler=DG*(one-beta*DMu)
                 mu_prime=(DMu-beta)/(one-beta*DMu)
                 anisotropy=(one-mu_prime*mu_prime)/(one+mu_prime*mu_prime)
                 stokes_weight=anisotropy*q_angle
                 log_doppler_redshift=log(doppler)+log(one+z)
                 F_log_theta=max(-199d0,F_log_theta+log_domega_4pi-3d0*log(doppler))
                 P_log_theta=max(-199d0,P_log_theta+log_domega_4pi-3d0*log(doppler))
                 V_seed_log_theta=V_seed_log-log_doppler_redshift
                 i_nu_src=1
                 do i_nu_obs=1,Num_nu_obs
                    if (V_obs_log(i_nu_obs) <= V_seed_log_theta(1)) cycle
                    if (V_obs_log(i_nu_obs) > V_seed_log_theta(Num_nu)) exit
                    do while (i_nu_src < Num_nu-1 .and. V_obs_log(i_nu_obs) > V_seed_log_theta(i_nu_src+1))
                       i_nu_src=i_nu_src+1
                    end do
                    if (V_obs_log(i_nu_obs) > V_seed_log_theta(i_nu_src) .and. &
                        V_obs_log(i_nu_obs) <= V_seed_log_theta(i_nu_src+1)) then
                        ratio_freq=(V_obs_log(i_nu_obs)-V_seed_log_theta(i_nu_src)) / &
                                   (V_seed_log_theta(i_nu_src+1)-V_seed_log_theta(i_nu_src))
                        I_obs_temp(i_nu_obs,K1)=I_obs_temp(i_nu_obs,K1)+ &
                            exp(F_log_theta(i_nu_src)+ratio_freq*(F_log_theta(i_nu_src+1)-F_log_theta(i_nu_src)))
                        Q_obs_temp(i_nu_obs,K1)=Q_obs_temp(i_nu_obs,K1)+stokes_weight * &
                            exp(P_log_theta(i_nu_src)+ratio_freq*(P_log_theta(i_nu_src+1)-P_log_theta(i_nu_src)))
                    end if
                 end do
             end if
         end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    I_obs=I_obs_temp*phi_scale
    Q_obs=Q_obs_temp*phi_scale
    U_obs=U_obs_temp*phi_scale
    deallocate (I_obs_temp,Q_obs_temp,U_obs_temp,V_obs_log,V_seed_log)
    return
end subroutine sed_interpolation_shock_random_stokes
