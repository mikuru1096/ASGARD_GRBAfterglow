! Created by rj on 2021/2/4.

!##################################################################################################
!When an intrinsic SED 'F_tot' observed in the comoving frame send in,
!to produce the observed SED 'F_tot_obs' after consider the EATS and Doppler boosting effect.
!##################################################################################################

subroutine sed_interpolation(Boundary,R_Tobs1,R_gamma,R,F_tot,V_seed,V_obs,Tobs, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads, F_tot_obs)
    !$ use omp_lib
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    !##############################################################################################
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads
    real(8), intent(in) :: Boundary(n)
    real(8), intent(in) :: Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: R_Tobs1(Num_R),R_gamma(Num_R),R(Num_R),F_tot(Num_nu,Num_R)
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    
    
    allocatable :: R_Tobs(:,:),DP(:,:),V_seed_temp(:,:),F_tot_temp(:,:),F_tot_obs_temp(:,:),V_obs_log(:)
    allocate (R_Tobs(Num_R,Num_Theta),DP(Num_nu,Num_Theta),V_seed_temp(Num_nu,Num_Theta), &
              F_tot_temp(Num_nu,Num_Theta),F_tot_obs_temp(Num_nu_obs,Num_Tobs),V_obs_log(Num_nu))

    F_tot_obs=zero
    F_tot_obs_temp=zero
    
    G00 = Boundary(1)
    R00 = Boundary(4)
    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)
    
    dPhi=pi/Num_Phi
    if (Num_Phi==1) dPhi=pi/1440d0
    dtheta=OpeningAngle_jet/Num_Theta

    V_obs_log = log(V_obs)
    !$OMP PARALLEL num_threads(n_threads),reduction(+:F_tot_obs_temp)
    !$OMP DO SIMD
    do I_Theta=1,Num_Theta
       Taa_boundary=dtheta*I_Theta
       Taa_center=dtheta*(I_Theta-0.5)
       domega=dsin(Taa_boundary)*dtheta*dPhi
       do i_Phi=1,Num_Phi
          Phi_center=(i_Phi-0.5)*dPhi
          DMu=dcos(Tv)*dcos(Taa_center)+dsin(Tv)*dsin(Taa_center)*dcos(Phi_center)
          R_Tobs(:,I_Theta)=R_Tobs1+R*(one-DMu)*(one+z)/Para_c 
          do K1=1,Num_Tobs
             II=1
             if (Tobs(K1) < R_Tobs(1,I_Theta).or.Tobs(K1) > R_Tobs(Num_R,I_Theta)) cycle
             do K2=II,Num_R-1
               if ((Tobs(K1) >= R_Tobs(K2,I_Theta)).and.(Tobs(K1) < R_Tobs(K2+1,I_Theta))) then
                   Ratio=(Tobs(K1)-R_Tobs(K2,I_Theta))/(R_Tobs(K2+1,I_Theta)-R_Tobs(K2,I_Theta))
                   DG=exp(log(R_gamma(K2))+Ratio*(log(R_gamma(K2+1))-log(R_gamma(K2))))
!                   DR=exp(log(R(K2))+Ratio*(log(R(K2+1))-log(R(K2))))
                   DP(:,I_Theta)=log(F_tot(:,K2))+Ratio*(log(F_tot(:,K2+1))-log(F_tot(:,K2)))
                   !logarithm interpolation to get the intrinsic SED at EATS.
                   Beta=dsqrt(one-DG**(-2))

                   doppler=DG*(one-Beta*DMu) !Doppler factor, changed with R
                   V_seed_temp(:,I_Theta)=log(V_seed/(doppler*(one+z))) !For frequency that has decayed with D, and shifted with (1+z)
                   F_tot_temp(:,I_Theta)=max(-399d0,DP(:,I_Theta)+log(domega)-log(4.0*pi)-3d0*log(doppler)) !For flux that has decayed with D^3

                   do i_nu1=1,Num_nu_obs
                      do i_nu2=1,Num_nu-1
                         if (V_obs_log(i_nu1) > V_seed_temp(i_nu2,I_Theta) .and. &
                             V_obs_log(i_nu1) <= V_seed_temp(i_nu2+1,I_Theta)) then
                             Ratio=(V_obs_log(i_nu1)-V_seed_temp(i_nu2,I_Theta))/ &
                                   (V_seed_temp(i_nu2+1,I_Theta)-V_seed_temp(i_nu2,I_Theta))
                             F_tot_obs_temp(i_nu1,K1)=F_tot_obs_temp(i_nu1,K1)+exp(F_tot_temp(i_nu2,I_Theta)+ &
                                                   Ratio*(F_tot_temp(i_nu2+1,I_Theta)-F_tot_temp(i_nu2,I_Theta)))
                             !another logarithm interpolation from the Doppler boosted freqency to the observed frequency
                         end if
                      end do
                   end do
                   II=K2
                   exit
               end if
            end do   
         end do
       end do
    end do
    !$OMP END DO SIMD
    !$OMP END PARALLEL
    
    F_tot_obs=F_tot_obs_temp*two!*Num_Phi
    if (Num_phi==1) then
        F_tot_obs=F_tot_obs_temp*two*1440d0
    end if
    
    deallocate (R_Tobs,DP,V_seed_temp,F_tot_temp,F_tot_obs_temp,V_obs_log)
    
    
    return
end subroutine sed_interpolation
