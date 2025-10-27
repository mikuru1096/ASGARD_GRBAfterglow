! Created by rj on 2021/2/4.

!##################################################################################################
!When an intrinsic SED 'P_tot' observed in the comoving frame send in,
!to produce the observed SED 'P_tot_obs' after consider the EATS and Dopper boosting effect.
!##################################################################################################

subroutine sed_interpolation_structured(Boundary, angle_narrow_jet, R_Tobs1,R_gamma,R,P_tot,V_seed,V_obs,Tobs, &
                             n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads, P_tot_obs)
    !$ use omp_lib
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    !##############################################################################################
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_Theta,Num_R,Num_Phi,n_threads
    real(8), intent(in) :: Boundary(n),angle_narrow_jet
    real(8), intent(in) :: Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(in) :: R_Tobs1(Num_R,Num_theta),R_gamma(Num_R,Num_theta)
    real(8), intent(in) :: R(Num_R,Num_theta),P_tot(Num_nu,Num_R,Num_theta)
    real(8), intent(out) :: P_tot_obs(Num_nu_obs,Num_Tobs)
    
    
    allocatable :: R_Tobs(:,:),DP(:,:),V_seed_temp(:,:),P_tot_temp(:,:),P_tot_obs_temp(:,:)
    allocate (R_Tobs(Num_R,Num_Theta),DP(Num_nu,Num_Theta),V_seed_temp(Num_nu,Num_Theta), &
              P_tot_temp(Num_nu,Num_Theta),P_tot_obs_temp(Num_nu_obs,Num_Tobs))
    
    
    V_seed_temp=zero
    DP=zero
    R_Tobs=zero
    P_tot_temp=zero
    P_tot_obs_temp=zero
    P_tot_obs=zero
    
    G00 = Boundary(1)
    R00 = Boundary(4)
    z = Boundary(8)
    OpeningAngle_jet = Boundary(9)
    Tv = Boundary(10)
    
    dPhi=pi/Num_Phi !Num_Phi
    if (Num_Phi==1) then
        dPhi=pi/1440d0
    end if
    dtheta=OpeningAngle_jet/Num_Theta

    !$OMP PARALLEL num_threads(n_threads),reduction(+:P_tot_obs_temp)
    !$OMP DO SIMD
    do I_Theta=1,Num_Theta
       Taa=dtheta*(I_Theta+0.5)
       if (angle_narrow_jet>Taa) cycle
       domega=dsin(Taa)*dtheta*dPhi
       do i_Phi=1,Num_Phi
          Phi=i_Phi*dPhi
          DMu=dcos(Tv)*dcos(Taa)+dsin(Tv)*dsin(Taa)*dcos(Phi)
          R_Tobs(:,I_Theta)=R_Tobs1(:,I_Theta)+R(:,I_Theta)*(one-DMu)*(one+z)/Para_c       
          do K1=1,Num_Tobs
             II=1
             if (Tobs(K1) < R_Tobs(1,I_Theta) .or. Tobs(K1) > R_Tobs(Num_R,I_Theta)) cycle
             do K2=II,Num_R-1
                if ((Tobs(K1) >= R_Tobs(K2,I_Theta)).and.(Tobs(K1) < R_Tobs(K2+1,I_Theta))) then
                   Ratio=(Tobs(K1)-R_Tobs(K2,I_Theta))/(R_Tobs(K2+1,I_Theta)-R_Tobs(K2,I_Theta))
                   DG=R_gamma(K2,I_Theta)+Ratio*(R_gamma(K2+1,I_Theta)-R_gamma(K2,I_Theta))
                   DR=R(K2,I_Theta)+Ratio*(R(K2+1,I_Theta)-R(K2,I_Theta))
                   DP(:,I_Theta)=P_tot(:,K2,I_Theta)+Ratio*(P_tot(:,K2+1,I_Theta)-P_tot(:,K2,I_Theta))
                   !linear interpolation to get the intrinsic SED at EATS.
                   Beta=dsqrt(one-DG**(-2))

                   doppler=DG*(one-Beta*DMu) !Doppler factor, changed with R
                   doppler3=doppler*doppler*doppler
                   V_seed_temp(:,I_Theta)=V_seed/(doppler*(one+z)) !For frequency that has decayed with D, and shifted with (1+z)
                   P_tot_temp(:,I_Theta)=DP(:,I_Theta)*domega/(4.0*pi*doppler3)
                   !For flux that has decayed with D^3

                   do i_nu1=1,Num_nu_obs
                      do i_nu2=1,Num_nu-1
                         if (V_obs(i_nu1) > V_seed_temp(i_nu2,I_Theta) .and. &
                             V_obs(i_nu1) <= V_seed_temp(i_nu2+1,I_Theta)) then
                             Ratio=(V_obs(i_nu1)-V_seed_temp(i_nu2,I_Theta))/ &
                                   (V_seed_temp(i_nu2+1,I_Theta)-V_seed_temp(i_nu2,I_Theta))
                             P_tot_obs_temp(i_nu1,K1)=P_tot_obs_temp(i_nu1,K1)+P_tot_temp(i_nu2,I_Theta)+ &
                                                   Ratio*(P_tot_temp(i_nu2+1,I_Theta)-P_tot_temp(i_nu2,I_Theta))
                             !another linear interpolation from the Doppler boosted freqency to the observed frequency
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

    P_tot_obs=P_tot_obs_temp*two!*Num_Phi
    if (Num_phi==1) then
        P_tot_obs=P_tot_obs_temp*two*1440d0
    end if
    
    deallocate (R_Tobs,DP,V_seed_temp,P_tot_temp,P_tot_obs_temp)
    
    
    return
end subroutine sed_interpolation_structured
