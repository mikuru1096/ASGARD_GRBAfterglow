subroutine ssc_spec(R,gam_e,dN_gam_e,V_seed,seed,Num_nu,Num_R,Num_gam_e,n_threads, P_SSC_spec,seed_SSC)
    use constants
    use radiation_common
    !$ use omp_lib
    IMPLICIT REAL(8)(A-H,O-Z)
    !***********************************************************
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,n_threads
    real(8), intent(in) :: R(Num_R),gam_e(Num_gam_e)
    real(8), intent(in) :: dN_gam_e(Num_gam_e,Num_R),V_seed(Num_nu),seed(Num_nu,Num_R)
    real(8), intent(out) :: P_SSC_spec(Num_nu,Num_R),seed_SSC(Num_nu,Num_R)
        
    real(8), allocatable :: simpson_weights(:), V_weights(:)
    real(8), allocatable :: E_seed(:), inv_gam(:), inv_gam2(:), radius_inv2(:)
    real(8), allocatable :: inv_v_seed(:), log_inv_v_seed(:)
    real(8), allocatable :: q_prefactor(:,:), q_gamma_coeff(:,:), kn_prefactor(:,:)
    real(8), allocatable :: log_q_prefactor(:,:)
    real(8), allocatable :: weighted_dN_over_gam(:,:), weighted_seed(:,:)
    real(8), allocatable :: tail_weighted_dN(:,:), tail_weighted_dN_inv2(:,:)
    integer, allocatable :: gamma_start(:), gamma_low(:,:), gamma_high(:,:)
    integer :: i_low, i_high, i_mid
    real(8) :: gamma_floor, temp_norm, h_nu_third, h_gam_third, q_coeff, kn_coeff

    allocate (simpson_weights(Num_gam_e), V_weights(Num_nu))
    allocate (E_seed(Num_nu), inv_gam(Num_gam_e), inv_gam2(Num_gam_e))
    allocate (radius_inv2(Num_R), inv_v_seed(Num_nu), log_inv_v_seed(Num_nu))
    allocate (q_prefactor(Num_gam_e,Num_nu), q_gamma_coeff(Num_gam_e,Num_nu))
    allocate (kn_prefactor(Num_gam_e,Num_nu), log_q_prefactor(Num_gam_e,Num_nu))
    allocate (weighted_dN_over_gam(Num_gam_e,Num_R), weighted_seed(Num_nu,Num_R))
    allocate (tail_weighted_dN(Num_gam_e+1,Num_R), tail_weighted_dN_inv2(Num_gam_e+1,Num_R))
    allocate (gamma_start(Num_nu), gamma_low(Num_nu,Num_nu), gamma_high(Num_nu,Num_nu))
    
!    call system_clock(int1)
    
    para_hEme = Para_h/para_m_energy

    h_nu = log(V_seed(2))-log(V_seed(1))
    h_gam = log(gam_e(2))-log(gam_e(1))
    E_seed = V_seed*para_hEme
    inv_gam = one/gam_e
    inv_gam2 = inv_gam*inv_gam
    inv_v_seed = one/V_seed
    log_inv_v_seed = log(inv_v_seed)
    radius_inv2 = one/(R*R)
    q_prefactor = zero
    q_gamma_coeff = zero
    kn_prefactor = zero
    log_q_prefactor = zero

    call compute_simpson_weights(simpson_weights, Num_gam_e)
    call compute_simpson_weights(V_weights, Num_nu)
    do I_R=1,Num_R
        weighted_dN_over_gam(:,I_R)=dN_gam_e(:,I_R)*simpson_weights*inv_gam
        weighted_seed(:,I_R)=seed(:,I_R)*V_weights
        tail_weighted_dN(Num_gam_e+1,I_R)=zero
        tail_weighted_dN_inv2(Num_gam_e+1,I_R)=zero
        do i_game=Num_gam_e,1,-1
            tail_weighted_dN(i_game,I_R)=tail_weighted_dN(i_game+1,I_R)+weighted_dN_over_gam(i_game,I_R)
            tail_weighted_dN_inv2(i_game,I_R)=tail_weighted_dN_inv2(i_game+1,I_R)+ &
                                              weighted_dN_over_gam(i_game,I_R)*inv_gam2(i_game)
        end do
    end do

    do I_nu=1,Num_nu
        i_low=1
        i_high=Num_gam_e+1
        do while (i_low < i_high)
            i_mid=(i_low+i_high)/2
            if (gam_e(i_mid) <= E_seed(I_nu)) then
                i_low=i_mid+1
            else
                i_high=i_mid
            end if
        end do
        gamma_start(I_nu)=i_low
        do i_game=i_low,Num_gam_e
            temp_norm = gam_e(i_game)-E_seed(I_nu)
            q_prefactor(i_game,I_nu)=V_seed(I_nu)/(4.0d0*gam_e(i_game)*temp_norm)
            q_gamma_coeff(i_game,I_nu)=E_seed(I_nu)/temp_norm
            kn_prefactor(i_game,I_nu)=q_gamma_coeff(i_game,I_nu)*q_gamma_coeff(i_game,I_nu)/(two*(one+q_gamma_coeff(i_game,I_nu)))
            log_q_prefactor(i_game,I_nu)=log(q_prefactor(i_game,I_nu))
        end do
    end do

    do I_nu=1,Num_nu
        do Nu_s=1,I_nu-1
            gamma_floor=0.5d0*(E_seed(I_nu)+sqrt(E_seed(I_nu)*E_seed(I_nu)+V_seed(I_nu)*inv_v_seed(Nu_s)))
            i_low=gamma_start(I_nu)
            i_high=Num_gam_e+1
            do while (i_low < i_high)
                i_mid=(i_low+i_high)/2
                if (gam_e(i_mid) <= gamma_floor) then
                    i_low=i_mid+1
                else
                    i_high=i_mid
                end if
            end do
            gamma_high(Nu_s,I_nu)=i_low
        end do
        do Nu_s=I_nu,Num_nu
            gamma_floor=0.5d0*sqrt(V_seed(Nu_s)/V_seed(I_nu))
            i_low=1
            i_high=Num_gam_e+1
            do while (i_low < i_high)
                i_mid=(i_low+i_high)/2
                if (gam_e(i_mid) < gamma_floor) then
                    i_low=i_mid+1
                else
                    i_high=i_mid
                end if
            end do
            gamma_low(Nu_s,I_nu)=i_low
        end do
    end do

    P_SSC_spec=zero
    seed_SSC=zero
    
    !$ call omp_set_dynamic(.true.)
    h_nu_third = h_nu/3.0d0
    h_gam_third = h_gam/3.0d0

    !$OMP PARALLEL num_threads(n_threads), private(I_R, I_nu, Nu_s, i_game, i, II, Vloc, ratio_v, q_coeff, &
    !$OMP& Ephoton2eV, dInteg, simpson_sum_nu, val1, val2, emission_int2, &
    !$OMP& simpson_sum_gam, P_v, F1, fssc, temp, q, q_gamma, log_q)
    !$OMP DO COLLAPSE(2) SCHEDULE(GUIDED,4)
    do I_R=1,Num_R
        do I_nu=1,Num_nu
            Vloc=V_seed(I_nu)
            Ephoton2eV=E_seed(I_nu)
            II=gamma_start(I_nu)
            if (II > Num_gam_e) cycle
            
            dInteg=zero
            simpson_sum_nu = zero
            do Nu_s=1,I_nu-1
               simpson_sum_gam = zero
               !$OMP SIMD REDUCTION(+:simpson_sum_gam)
               do i_game = gamma_high(Nu_s,I_nu),Num_gam_e
                  q_coeff = q_prefactor(i_game,I_nu)
                  q = q_coeff*inv_v_seed(Nu_s)
                  if (q >= one) cycle
                  log_q = log_q_prefactor(i_game,I_nu) + log_inv_v_seed(Nu_s)
                  fssc = two*q*(log_q-q)+one+q+kn_prefactor(i_game,I_nu)*(one-q)
                  val2 = weighted_dN_over_gam(i_game, I_R) * fssc
                  simpson_sum_gam = simpson_sum_gam + val2
               end do
               emission_int2 = h_gam_third * simpson_sum_gam
               simpson_sum_nu = simpson_sum_nu + weighted_seed(Nu_s, I_R) * emission_int2
            end do
            do Nu_s=I_nu,Num_nu
               ratio_v = Vloc*inv_v_seed(Nu_s)
               i_game=gamma_low(Nu_s,I_nu)
               if (i_game <= Num_gam_e) then
                  simpson_sum_gam = ratio_v*tail_weighted_dN(i_game,I_R) - &
                                    0.25d0*tail_weighted_dN_inv2(i_game,I_R)
               else
                  simpson_sum_gam = zero
               end if
               emission_int2 = h_gam_third * simpson_sum_gam
               simpson_sum_nu = simpson_sum_nu + weighted_seed(Nu_s, I_R) * emission_int2
            end do
            dInteg = h_nu_third * simpson_sum_nu
            
            P_v=dInteg*Vloc
            P_SSC_spec(I_nu,I_R)=P_SSC_spec(I_nu,I_R)+P_v
            F1=dInteg*radius_inv2(I_R)
            seed_SSC(I_nu,I_R)=seed_SSC(I_nu,I_R)+F1
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    Temp_para=0.75d0*Para_c*Para_h*Para_SigmaT
    P_SSC_spec=P_SSC_spec*Temp_para
    
    Temp_para2=4.0d0*pi*Para_c*Para_h
    seed_SSC=seed_SSC/Temp_para2*Temp_para
    
!    call system_clock(int2)
!    print*, 'time=', (int2-int1)/1000.0
    
    deallocate(simpson_weights, V_weights, E_seed, inv_gam, inv_gam2)
    deallocate(radius_inv2, inv_v_seed, log_inv_v_seed)
    deallocate(q_prefactor, q_gamma_coeff, kn_prefactor, log_q_prefactor)
    deallocate(weighted_dN_over_gam, weighted_seed, tail_weighted_dN, tail_weighted_dN_inv2)
    deallocate(gamma_start, gamma_low, gamma_high)

    return
end subroutine ssc_spec

subroutine ssc_spec_nonuniform(R,x_edge_log10,dN_x,V_seed,seed,Num_nu,Num_R,Num_gam_e,n_threads, &
                               P_SSC_spec,seed_SSC)
    use constants
    !$ use omp_lib
    implicit REAL(8)(A-H,O-Z)
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,n_threads
    real(8), intent(in) :: R(Num_R),x_edge_log10(Num_gam_e+1,Num_R),dN_x(Num_gam_e,Num_R)
    real(8), intent(in) :: V_seed(Num_nu),seed(Num_nu,Num_R)
    real(8), intent(out) :: P_SSC_spec(Num_nu,Num_R),seed_SSC(Num_nu,Num_R)

    real(8), allocatable :: x_seed(:),radius_inv2(:)
    real(8), allocatable :: x_center(:,:),slope_q(:,:),tail_gamma(:,:),tail_gamma_inv2(:,:)
    real(8) :: ratio_v,gamma_floor,dInteg,emission_int2,P_v,F1,Vloc,Ephoton2eV
    real(8) :: x_seed0,x_seed1,xm_seed,dx_seed,x_seed_loc,V_seed_loc,E_seed_loc
    real(8) :: seed_loc,seed_weight,dx_log10,x0_cell,x1_cell,w2
    real(8) :: left_slope,right_slope,slope_lim
    integer :: I_gauss

    allocate(x_seed(Num_nu),radius_inv2(Num_R),x_center(Num_gam_e,Num_R),slope_q(Num_gam_e,Num_R))
    allocate(tail_gamma(Num_gam_e+1,Num_R),tail_gamma_inv2(Num_gam_e+1,Num_R))

    para_hEme=Para_h/para_m_energy
    x_seed=dlog(V_seed)
    radius_inv2=one/(R*R)
    w2=one/dsqrt(3d0)

    do I_R=1,Num_R
        do i_game=1,Num_gam_e
            x_center(i_game,I_R)=0.5d0*(x_edge_log10(i_game,I_R)+x_edge_log10(i_game+1,I_R))
        end do
        do i_game=1,Num_gam_e
            dx_log10=max(x_edge_log10(i_game+1,I_R)-x_edge_log10(i_game,I_R),1d-30)
            if (Num_gam_e == 1) then
                slope_q(i_game,I_R)=zero
            else if (i_game == 1) then
                right_slope=(dN_x(2,I_R)-dN_x(1,I_R))/max(x_center(2,I_R)-x_center(1,I_R),1d-30)
                slope_q(i_game,I_R)=right_slope
            else if (i_game == Num_gam_e) then
                left_slope=(dN_x(Num_gam_e,I_R)-dN_x(Num_gam_e-1,I_R))/ &
                           max(x_center(Num_gam_e,I_R)-x_center(Num_gam_e-1,I_R),1d-30)
                slope_q(i_game,I_R)=left_slope
            else
                left_slope=(dN_x(i_game,I_R)-dN_x(i_game-1,I_R))/ &
                           max(x_center(i_game,I_R)-x_center(i_game-1,I_R),1d-30)
                right_slope=(dN_x(i_game+1,I_R)-dN_x(i_game,I_R))/ &
                            max(x_center(i_game+1,I_R)-x_center(i_game,I_R),1d-30)
                slope_q(i_game,I_R)=ssc_minmod(left_slope,right_slope)
            end if
            slope_lim=two*dN_x(i_game,I_R)/dx_log10
            if (abs(slope_q(i_game,I_R)) > slope_lim) then
                slope_q(i_game,I_R)=sign(slope_lim,slope_q(i_game,I_R))
            end if
        end do
        tail_gamma(Num_gam_e+1,I_R)=zero
        tail_gamma_inv2(Num_gam_e+1,I_R)=zero
        do i_game=Num_gam_e,1,-1
            x0_cell=x_edge_log10(i_game,I_R)
            x1_cell=x_edge_log10(i_game+1,I_R)
            tail_gamma(i_game,I_R)=tail_gamma(i_game+1,I_R)+ &
                                   linear_gamma_moment(x0_cell,x1_cell,dN_x(i_game,I_R),slope_q(i_game,I_R), &
                                                       x_center(i_game,I_R),2d0)
            tail_gamma_inv2(i_game,I_R)=tail_gamma_inv2(i_game+1,I_R)+ &
                                        linear_gamma_moment(x0_cell,x1_cell,dN_x(i_game,I_R),slope_q(i_game,I_R), &
                                                            x_center(i_game,I_R),4d0)
        end do
    end do

    P_SSC_spec=zero
    seed_SSC=zero

    !$OMP PARALLEL num_threads(n_threads), private(I_R,I_nu,Nu_s,I_gauss,Vloc,Ephoton2eV,ratio_v,gamma_floor, &
    !$OMP& dInteg,emission_int2,P_v,F1,x_seed0,x_seed1,xm_seed,dx_seed,x_seed_loc,V_seed_loc,E_seed_loc, &
    !$OMP& seed_loc,seed_weight)
    !$OMP DO COLLAPSE(2) SCHEDULE(GUIDED,4)
    do I_R=1,Num_R
        do I_nu=1,Num_nu
            Vloc=V_seed(I_nu)
            Ephoton2eV=Vloc*para_hEme
            dInteg=zero
            do Nu_s=1,Num_nu-1
                x_seed0=x_seed(Nu_s)
                x_seed1=x_seed(Nu_s+1)
                if (x_seed1 <= x_seed0) cycle
                xm_seed=0.5d0*(x_seed0+x_seed1)
                dx_seed=0.5d0*(x_seed1-x_seed0)
                do I_gauss=1,2
                    if (I_gauss == 1) then
                        x_seed_loc=xm_seed-dx_seed*w2
                    else
                        x_seed_loc=xm_seed+dx_seed*w2
                    end if
                    V_seed_loc=dexp(x_seed_loc)
                    E_seed_loc=V_seed_loc*para_hEme
                    seed_loc=seed_log_interp(V_seed(Nu_s),V_seed(Nu_s+1),seed(Nu_s,I_R), &
                                             seed(Nu_s+1,I_R),V_seed_loc)
                    if (seed_loc <= zero) cycle
                    seed_weight=dx_seed*seed_loc
                    if (V_seed_loc < Vloc) then
                        gamma_floor=0.5d0*(Ephoton2eV+dsqrt(Ephoton2eV*Ephoton2eV+Vloc/V_seed_loc))
                        emission_int2=ssc_low_gamma_integral(I_R,dlog10(gamma_floor),Vloc,V_seed_loc,Ephoton2eV)
                    else
                        ratio_v=Vloc/V_seed_loc
                        gamma_floor=0.5d0*dsqrt(V_seed_loc/Vloc)
                        emission_int2=ssc_high_gamma_tail(I_R,dlog10(gamma_floor),ratio_v)
                    end if
                    dInteg=dInteg+seed_weight*emission_int2
                end do
            end do
            P_v=dInteg*Vloc
            P_SSC_spec(I_nu,I_R)=P_v
            F1=dInteg*radius_inv2(I_R)
            seed_SSC(I_nu,I_R)=F1
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    Temp_para=0.75d0*Para_c*Para_h*Para_SigmaT
    P_SSC_spec=P_SSC_spec*Temp_para

    Temp_para2=4.0d0*pi*Para_c*Para_h
    seed_SSC=seed_SSC/Temp_para2*Temp_para

    deallocate(x_seed,radius_inv2,x_center,slope_q)
    deallocate(tail_gamma,tail_gamma_inv2)

contains

subroutine first_greater_monotonic(arr,n,target,idx)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: arr(n),target
    integer, intent(out) :: idx
    integer :: left,right,mid

    if (arr(1) > target) then
        idx=1
        return
    end if
    if (arr(n) <= target) then
        idx=n+1
        return
    end if

    left=1
    right=n
    do while (left < right)
        mid=(left+right)/2
        if (arr(mid) > target) then
            right=mid
        else
            left=mid+1
        end if
    end do
    idx=left
end subroutine first_greater_monotonic

integer function first_cell_above_edge(x_edge_col,n,x_floor)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x_edge_col(n+1),x_floor
    integer :: left,right,mid

    if (x_floor <= x_edge_col(1)) then
        first_cell_above_edge=1
        return
    end if
    if (x_floor >= x_edge_col(n+1)) then
        first_cell_above_edge=n+1
        return
    end if

    left=1
    right=n
    do while (left < right)
        mid=(left+right)/2
        if (x_edge_col(mid+1) > x_floor) then
            right=mid
        else
            left=mid+1
        end if
    end do
    first_cell_above_edge=left
end function first_cell_above_edge

real(8) function inv_gamma_power_integral(x0,x1,power)
    implicit REAL(8)(A-H,O-Z)
    real(8), intent(in) :: x0,x1,power

    if (x1 <= x0) then
        inv_gamma_power_integral=zero
    else
        inv_gamma_power_integral=(ten**(-power*x0)-ten**(-power*x1))/(power*dlog(ten))
    end if
end function inv_gamma_power_integral

real(8) function ssc_minmod(a,b)
    implicit REAL(8)(A-H,O-Z)
    real(8), intent(in) :: a,b

    if (a*b <= zero) then
        ssc_minmod=zero
    else
        ssc_minmod=sign(min(abs(a),abs(b)),a)
    end if
end function ssc_minmod

real(8) function linear_profile_value(x,qbar,slope,xc)
    implicit REAL(8)(A-H,O-Z)
    real(8), intent(in) :: x,qbar,slope,xc

    linear_profile_value=max(zero,qbar+slope*(x-xc))
end function linear_profile_value

real(8) function linear_gamma_moment(x0,x1,qbar,slope,xc,power)
    implicit REAL(8)(A-H,O-Z)
    real(8), intent(in) :: x0,x1,qbar,slope,xc,power
    real(8) :: alpha,exp0,exp1,i0,i1

    if (x1 <= x0) then
        linear_gamma_moment=zero
        return
    end if

    alpha=power*dlog(ten)
    exp0=ten**(-power*x0)
    exp1=ten**(-power*x1)
    i0=(exp0-exp1)/alpha
    i1=((x0/alpha)+one/(alpha*alpha))*exp0-((x1/alpha)+one/(alpha*alpha))*exp1
    linear_gamma_moment=qbar*i0+slope*(i1-xc*i0)
end function linear_gamma_moment

real(8) function ssc_low_gamma_cell(x0,x1,qbar,slope,xc,Vloc,V_seed_loc,Ephoton2eV)
    implicit REAL(8)(A-H,O-Z)
    real(8), intent(in) :: x0,x1,qbar,slope,xc,Vloc,V_seed_loc,Ephoton2eV
    real(8) :: xm,dx,w2_loc,xg,gam,temp,q,q_gamma,kn_coeff,fssc,val,q_loc
    integer :: I_pt

    if (x1 <= x0 .or. qbar == zero) then
        ssc_low_gamma_cell=zero
        return
    end if

    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    w2_loc=one/dsqrt(3d0)
    ssc_low_gamma_cell=zero
    do I_pt=1,2
        if (I_pt == 1) then
            xg=xm-dx*w2_loc
        else
            xg=xm+dx*w2_loc
        end if
        q_loc=linear_profile_value(xg,qbar,slope,xc)
        if (q_loc <= zero) cycle
        gam=ten**xg
        temp=gam-Ephoton2eV
        if (temp <= zero) cycle
        q=Vloc/(4.0d0*gam*temp*V_seed_loc)
        if (q <= zero .or. q >= one) cycle
        q_gamma=Ephoton2eV/temp
        kn_coeff=q_gamma*q_gamma/(two*(one+q_gamma))
        fssc=two*q*(dlog(q)-q)+one+q+kn_coeff*(one-q)
        val=q_loc*fssc/(gam*gam)
        ssc_low_gamma_cell=ssc_low_gamma_cell+val
    end do
    ssc_low_gamma_cell=ssc_low_gamma_cell*dx
end function ssc_low_gamma_cell

real(8) function ssc_low_gamma_integral(I_R,x_floor,Vloc,V_seed_loc,Ephoton2eV)
    implicit REAL(8)(A-H,O-Z)
    integer, intent(in) :: I_R
    real(8), intent(in) :: x_floor,Vloc,V_seed_loc,Ephoton2eV
    integer :: i_start,i_game
    real(8) :: x0,x1

    i_start=first_cell_above_edge(x_edge_log10(:,I_R),Num_gam_e,x_floor)
    if (i_start > Num_gam_e) then
        ssc_low_gamma_integral=zero
        return
    end if

    ssc_low_gamma_integral=zero
    do i_game=i_start,Num_gam_e
        x0=max(x_edge_log10(i_game,I_R),x_floor)
        x1=x_edge_log10(i_game+1,I_R)
        ssc_low_gamma_integral=ssc_low_gamma_integral+ &
            ssc_low_gamma_cell(x0,x1,dN_x(i_game,I_R),slope_q(i_game,I_R),x_center(i_game,I_R), &
                               Vloc,V_seed_loc,Ephoton2eV)
    end do
end function ssc_low_gamma_integral

real(8) function ssc_high_gamma_tail(I_R,x_floor,ratio_v)
    implicit REAL(8)(A-H,O-Z)
    integer, intent(in) :: I_R
    real(8), intent(in) :: x_floor,ratio_v
    integer :: i_start
    real(8) :: x0,x1,part2,part4

    i_start=first_cell_above_edge(x_edge_log10(:,I_R),Num_gam_e,x_floor)
    if (i_start > Num_gam_e) then
        ssc_high_gamma_tail=zero
        return
    end if

    x0=max(x_edge_log10(i_start,I_R),x_floor)
    x1=x_edge_log10(i_start+1,I_R)
    part2=linear_gamma_moment(x0,x1,dN_x(i_start,I_R),slope_q(i_start,I_R),x_center(i_start,I_R),2d0)
    part4=linear_gamma_moment(x0,x1,dN_x(i_start,I_R),slope_q(i_start,I_R),x_center(i_start,I_R),4d0)
    ssc_high_gamma_tail=ratio_v*(part2+tail_gamma(i_start+1,I_R))- &
                        0.25d0*(part4+tail_gamma_inv2(i_start+1,I_R))
end function ssc_high_gamma_tail

real(8) function seed_log_interp(v0,v1,y0,y1,v)
    implicit REAL(8)(A-H,O-Z)
    real(8), intent(in) :: v0,v1,y0,y1,v
    real(8) :: slope

    if (v <= v0) then
        seed_log_interp=y0
        return
    end if
    if (v >= v1) then
        seed_log_interp=y1
        return
    end if

    if (v1 <= v0) then
        seed_log_interp=0.5d0*(y0+y1)
    else if (y0 > zero .and. y1 > zero) then
        slope=dlog(y1/y0)/dlog(v1/v0)
        seed_log_interp=y0*(v/v0)**slope
    else
        seed_log_interp=y0+(y1-y0)*(v-v0)/(v1-v0)
    end if
end function seed_log_interp

end subroutine ssc_spec_nonuniform
