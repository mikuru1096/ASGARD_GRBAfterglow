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
    !$OMP PARALLEL num_threads(n_threads), reduction(+:F_tot_obs_temp), private(I_Theta, Taa_boundary, &
    !$OMP& Taa_center, domega, i_Phi, Phi_center, DMu, K1, II, K2, Ratio, DG, Beta, doppler, &
    !$OMP& R_Tobs_theta, F_tot_theta, F_tot_log_theta, V_seed_log_theta, &
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

! 将单个真实球面面元的共动系SED投影到观测系：固定视角DMu，用显式dOmega权重累加EATS+Doppler。
subroutine sed_interpolation_surface_element(Boundary,R_Tobs1,R_gamma,R,F_tot,V_seed,V_obs,Tobs,Domega, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_R, F_tot_obs)
    use constants
    use interpolation_common
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_R
    real(8), intent(in) :: Boundary(n),Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs),Domega
    real(8), intent(in) :: R_Tobs1(Num_R),R_gamma(Num_R),R(Num_R),F_tot(Num_nu,Num_R)
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    real(8), allocatable :: V_obs_log(:),V_seed_log(:)
    real(8) :: R_Tobs_theta(Num_R),F_tot_theta(Num_nu),F_tot_log_theta(Num_nu), &
               V_seed_log_theta(Num_nu),log_gamma_lo,log_gamma_hi,log_domega_4pi,log_doppler_redshift
    integer :: last_k2
    allocate (V_obs_log(Num_nu_obs),V_seed_log(Num_nu))

    F_tot_obs=zero
    z = Boundary(8)
    Tv = Boundary(10)
    DMu = dcos(Tv)
    log_domega_4pi=log(Domega)-log(4.0d0*pi)
    V_obs_log = log(V_obs)
    V_seed_log = log(V_seed)
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
               last_k2=K2
           end if
           Ratio=(Tobs(K1)-R_Tobs_theta(K2))/(R_Tobs_theta(K2+1)-R_Tobs_theta(K2))
           DG=exp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
           F_tot_theta=(one-Ratio)*F_tot(:,K2)+Ratio*F_tot(:,K2+1)
           F_tot_log_theta=-huge(one)
           where(F_tot_theta > zero) F_tot_log_theta=log(F_tot_theta)
           Beta=dsqrt(one-DG**(-2))
           doppler=DG*(one-Beta*DMu)
           log_doppler_redshift=log(doppler)+log(one+z)
           F_tot_log_theta=F_tot_log_theta+log_domega_4pi-3d0*log(doppler)
           V_seed_log_theta=V_seed_log-log_doppler_redshift
           call interpolation_accumulate_log_sed(V_seed_log_theta,F_tot_log_theta, &
                                                 Num_nu,V_obs_log,Num_nu_obs,F_tot_obs(:,K1))
       end if
    end do
    deallocate (V_obs_log,V_seed_log)
    return
end subroutine sed_interpolation_surface_element
