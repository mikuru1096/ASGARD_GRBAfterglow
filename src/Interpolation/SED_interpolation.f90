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
    real(8) :: R_Tobs_chi(Num_R),F_theta(Num_nu)
    real(8) :: log_domega_4pi,log_doppler_redshift,log_gamma_lo,log_gamma_hi
    real(8) :: source_lo,source_hi,tau_front_lo,tau_front_hi,tau_cell_lo,tau_cell_hi
    real(8) :: segment_lo,segment_hi,cos_tv,sin_tv,log_flux_weight
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
        Taa_boundary = dtheta*I_Theta
        Taa_center = dtheta*(I_Theta-0.5d0)
        Phi_center = (i_Phi-0.5d0)*dPhi
        domega = dsin(Taa_boundary)*dtheta*dPhi
        Angle_log_domega(I_ang) = dlog(domega)-dlog(4d0*pi)
        Angle_dmu(I_ang) = cos_tv*dcos(Taa_center)+sin_tv*dsin(Taa_center)*dcos(Phi_center)
    end do

    !$OMP PARALLEL num_threads(n_threads), reduction(+:F_temp), private(I_ang,log_domega_4pi,DMu, &
    !$OMP& I_chi,R_Tobs_chi,monotonic_chi,II,last_k2,K1,K2, &
    !$OMP& log_gamma_lo,log_gamma_hi,Ratio,DG,Beta,doppler,log_doppler_redshift,F_theta,log_flux_weight, &
    !$OMP& I_nu,tau_front_lo,tau_front_hi,tau_cell_lo,tau_cell_hi, &
    !$OMP& source_lo,source_hi,temp,segment_lo,segment_hi,k_start)
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
                            DG = dexp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
                            Beta = dsqrt(one-DG**(-2))
                            doppler = DG*(one-Beta*DMu)
                            log_doppler_redshift = dlog(doppler)+dlog(one+z)
                            F_theta = zero
                            do I_nu = 1, Num_nu
                                tau_front_lo = Tau_prefix(I_nu,I_chi-1,K2)
                                tau_front_hi = Tau_prefix(I_nu,I_chi-1,K2+1)
                                tau_cell_lo = Tau_chi(I_nu,I_chi,K2)
                                tau_cell_hi = Tau_chi(I_nu,I_chi,K2+1)
                                source_lo = F_chi(I_nu,I_chi,K2)*Chi_weight(I_chi,K2)* &
                                            chi_ssa_cell_escape(tau_front_lo,tau_cell_lo)
                                source_hi = F_chi(I_nu,I_chi,K2+1)*Chi_weight(I_chi,K2+1)* &
                                            chi_ssa_cell_escape(tau_front_hi,tau_cell_hi)
                                F_theta(I_nu) = (one-Ratio)*source_lo + Ratio*source_hi
                            end do
                            log_flux_weight = log_domega_4pi - 3d0*dlog(doppler)
                            call interpolation_accumulate_shifted_linear_sed(V_seed_log,F_theta,Num_nu, &
                                                                             V_obs_log,Num_nu_obs, &
                                                                             log_doppler_redshift,log_flux_weight,F_temp(:,K1))
                        end if
                    end do
                else
                    do K2 = 1, Num_R-1
                        segment_lo = min(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                        segment_hi = max(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                        k_start = lower_bound_real8(Tobs,Num_Tobs,segment_lo)
                        do K1 = k_start, Num_Tobs
                            if (Tobs(K1) >= segment_hi) exit
                            Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                            log_gamma_lo = dlog(Gamma_chi(I_chi,K2))
                            log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1))
                            DG = dexp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
                            Beta = dsqrt(one-DG**(-2))
                            doppler = DG*(one-Beta*DMu)
                            log_doppler_redshift = dlog(doppler)+dlog(one+z)
                            F_theta = zero
                            do I_nu = 1, Num_nu
                                tau_front_lo = Tau_prefix(I_nu,I_chi-1,K2)
                                tau_front_hi = Tau_prefix(I_nu,I_chi-1,K2+1)
                                tau_cell_lo = Tau_chi(I_nu,I_chi,K2)
                                tau_cell_hi = Tau_chi(I_nu,I_chi,K2+1)
                                source_lo = F_chi(I_nu,I_chi,K2)*Chi_weight(I_chi,K2)* &
                                            chi_ssa_cell_escape(tau_front_lo,tau_cell_lo)
                                source_hi = F_chi(I_nu,I_chi,K2+1)*Chi_weight(I_chi,K2+1)* &
                                            chi_ssa_cell_escape(tau_front_hi,tau_cell_hi)
                                F_theta(I_nu) = (one-Ratio)*source_lo + Ratio*source_hi
                            end do
                            log_flux_weight = log_domega_4pi - 3d0*dlog(doppler)
                            call interpolation_accumulate_shifted_linear_sed(V_seed_log,F_theta,Num_nu, &
                                                                             V_obs_log,Num_nu_obs, &
                                                                             log_doppler_redshift,log_flux_weight,F_temp(:,K1))
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
end subroutine sed_interpolation_chi

! 将χ分辨有限厚壳层面元投影到观测系：每个χ体元使用局域半径、局域Γ和向外SSA存活率。
subroutine sed_interpolation_chi_surface_element(Boundary,R_Tobs1,R_front,F_chi,Tau_chi,R_chi,Gamma_chi,Chi_weight, &
                             V_seed,V_obs,Tobs,Domega,n,Num_nu,Num_nu_obs,Num_Tobs,Num_chi,Num_R,F_tot_obs)
    use constants
    use interpolation_common
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_nu_obs,Num_Tobs,Num_chi,Num_R
    real(8), intent(in) :: Boundary(n),R_Tobs1(Num_R),R_front(Num_R),Domega
    real(8), intent(in) :: F_chi(Num_nu,Num_chi,Num_R),Tau_chi(Num_nu,Num_chi,Num_R)
    real(8), intent(in) :: R_chi(Num_chi,Num_R),Gamma_chi(Num_chi,Num_R),Chi_weight(Num_chi,Num_R)
    real(8), intent(in) :: Tobs(Num_Tobs),V_seed(Num_nu),V_obs(Num_nu_obs)
    real(8), intent(out) :: F_tot_obs(Num_nu_obs,Num_Tobs)
    real(8), allocatable :: V_obs_log(:),V_seed_log(:)
    real(8) :: R_Tobs_chi(Num_R),F_log_theta(Num_nu),V_seed_log_theta(Num_nu)
    real(8) :: log_domega_4pi,log_doppler_redshift,log_gamma_lo,log_gamma_hi
    real(8) :: source_lo,source_hi,tau_front_lo,tau_front_hi,tau_cell_lo,tau_cell_hi
    real(8) :: segment_lo,segment_hi
    integer :: last_k2,k_start,lower_bound_real8
    logical :: monotonic_chi

    allocate(V_obs_log(Num_nu_obs),V_seed_log(Num_nu))
    F_tot_obs = zero
    z = Boundary(8)
    Tv = Boundary(10)
    DMu = dcos(Tv)
    log_domega_4pi = dlog(Domega)-dlog(4d0*pi)
    V_obs_log = dlog(V_obs)
    V_seed_log = dlog(V_seed)

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
                    DG = dexp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
                    Beta = dsqrt(one-DG**(-2))
                    doppler = DG*(one-Beta*DMu)
                    log_doppler_redshift = dlog(doppler)+dlog(one+z)
                    F_log_theta = -huge(one)
                    do I_nu = 1, Num_nu
                        tau_front_lo = zero
                        tau_front_hi = zero
                        do J_chi = 1, I_chi-1
                            tau_front_lo = tau_front_lo + Tau_chi(I_nu,J_chi,K2)
                            tau_front_hi = tau_front_hi + Tau_chi(I_nu,J_chi,K2+1)
                        end do
                        tau_cell_lo = Tau_chi(I_nu,I_chi,K2)
                        tau_cell_hi = Tau_chi(I_nu,I_chi,K2+1)
                        source_lo = F_chi(I_nu,I_chi,K2)*Chi_weight(I_chi,K2)* &
                                    chi_ssa_cell_escape(tau_front_lo,tau_cell_lo)
                        source_hi = F_chi(I_nu,I_chi,K2+1)*Chi_weight(I_chi,K2+1)* &
                                    chi_ssa_cell_escape(tau_front_hi,tau_cell_hi)
                        temp = (one-Ratio)*source_lo + Ratio*source_hi
                        if (temp > zero) F_log_theta(I_nu) = dlog(temp)
                    end do
                    F_log_theta = F_log_theta + log_domega_4pi - 3d0*dlog(doppler)
                    V_seed_log_theta = V_seed_log - log_doppler_redshift
                    call interpolation_accumulate_log_sed(V_seed_log_theta,F_log_theta,Num_nu,V_obs_log,Num_nu_obs,F_tot_obs(:,K1))
                end if
            end do
        else
            do K2 = 1, Num_R-1
                segment_lo = min(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                segment_hi = max(R_Tobs_chi(K2),R_Tobs_chi(K2+1))
                k_start = lower_bound_real8(Tobs,Num_Tobs,segment_lo)
                do K1 = k_start, Num_Tobs
                    if (Tobs(K1) >= segment_hi) exit
                    Ratio = (Tobs(K1)-R_Tobs_chi(K2))/(R_Tobs_chi(K2+1)-R_Tobs_chi(K2))
                    log_gamma_lo = dlog(Gamma_chi(I_chi,K2))
                    log_gamma_hi = dlog(Gamma_chi(I_chi,K2+1))
                    DG = dexp(log_gamma_lo+Ratio*(log_gamma_hi-log_gamma_lo))
                    Beta = dsqrt(one-DG**(-2))
                    doppler = DG*(one-Beta*DMu)
                    log_doppler_redshift = dlog(doppler)+dlog(one+z)
                    F_log_theta = -huge(one)
                    do I_nu = 1, Num_nu
                        tau_front_lo = zero
                        tau_front_hi = zero
                        do J_chi = 1, I_chi-1
                            tau_front_lo = tau_front_lo + Tau_chi(I_nu,J_chi,K2)
                            tau_front_hi = tau_front_hi + Tau_chi(I_nu,J_chi,K2+1)
                        end do
                        tau_cell_lo = Tau_chi(I_nu,I_chi,K2)
                        tau_cell_hi = Tau_chi(I_nu,I_chi,K2+1)
                        source_lo = F_chi(I_nu,I_chi,K2)*Chi_weight(I_chi,K2)* &
                                    chi_ssa_cell_escape(tau_front_lo,tau_cell_lo)
                        source_hi = F_chi(I_nu,I_chi,K2+1)*Chi_weight(I_chi,K2+1)* &
                                    chi_ssa_cell_escape(tau_front_hi,tau_cell_hi)
                        temp = (one-Ratio)*source_lo + Ratio*source_hi
                        if (temp > zero) F_log_theta(I_nu) = dlog(temp)
                    end do
                    F_log_theta = F_log_theta + log_domega_4pi - 3d0*dlog(doppler)
                    V_seed_log_theta = V_seed_log - log_doppler_redshift
                    call interpolation_accumulate_log_sed(V_seed_log_theta,F_log_theta,Num_nu,V_obs_log,Num_nu_obs,F_tot_obs(:,K1))
                end do
            end do
        end if
    end do
    deallocate(V_obs_log,V_seed_log)
    return
end subroutine sed_interpolation_chi_surface_element
