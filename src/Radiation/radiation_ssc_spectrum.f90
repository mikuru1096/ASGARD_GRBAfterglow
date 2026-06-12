! 计算同步自康普顿（SSC）辐射谱：均匀对数网格上对电子谱和种子光子双重Simpson积分。
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
    real(8), allocatable :: q_prefactor(:,:), kn_prefactor(:,:)
    real(8), allocatable :: log_q_prefactor(:,:)
    real(8), allocatable :: weighted_dN_over_gam(:,:), weighted_seed(:,:)
    real(8), allocatable :: tail_weighted_dN(:,:), tail_weighted_dN_inv2(:,:)
    integer, allocatable :: gamma_start(:), gamma_low(:,:), gamma_high(:,:)
    integer :: I_R, I_nu
    real(8) :: h_nu_third, h_gam_third

    allocate (simpson_weights(Num_gam_e), V_weights(Num_nu))
    allocate (E_seed(Num_nu), inv_gam(Num_gam_e), inv_gam2(Num_gam_e))
    allocate (radius_inv2(Num_R), inv_v_seed(Num_nu), log_inv_v_seed(Num_nu))
    allocate (q_prefactor(Num_gam_e,Num_nu), kn_prefactor(Num_gam_e,Num_nu))
    allocate (log_q_prefactor(Num_gam_e,Num_nu))
    allocate (weighted_dN_over_gam(Num_gam_e,Num_R), weighted_seed(Num_nu,Num_R))
    allocate (tail_weighted_dN(Num_gam_e+1,Num_R), tail_weighted_dN_inv2(Num_gam_e+1,Num_R))
    allocate (gamma_start(Num_nu), gamma_low(Num_nu,Num_nu), gamma_high(Num_nu,Num_nu))

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
    kn_prefactor = zero
    log_q_prefactor = zero

    call prepare_uniform_weights_and_tails()
    call prepare_kn_tables()
    call prepare_uniform_gamma_bounds()

    P_SSC_spec=zero
    seed_SSC=zero
    
    !$ call omp_set_dynamic(.true.)
    h_nu_third = h_nu/3.0d0
    h_gam_third = h_gam/3.0d0

    !$OMP PARALLEL num_threads(n_threads), private(I_R, I_nu)
    !$OMP DO COLLAPSE(2) SCHEDULE(GUIDED,4)
    do I_R=1,Num_R
        do I_nu=1,Num_nu
            call accumulate_uniform_point(I_R, I_nu)
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    Temp_para=0.75d0*Para_c*Para_h*Para_SigmaT
    P_SSC_spec=P_SSC_spec*Temp_para
    
    Temp_para2=4.0d0*pi*Para_c*Para_h
    seed_SSC=seed_SSC/Temp_para2*Temp_para
    
    deallocate(simpson_weights, V_weights, E_seed, inv_gam, inv_gam2)
    deallocate(radius_inv2, inv_v_seed, log_inv_v_seed)
    deallocate(q_prefactor, kn_prefactor, log_q_prefactor)
    deallocate(weighted_dN_over_gam, weighted_seed, tail_weighted_dN, tail_weighted_dN_inv2)
    deallocate(gamma_start, gamma_low, gamma_high)

    return

contains

integer function first_gamma_gt(gamma_floor, lower_bound)
    implicit none
    real(8), intent(in) :: gamma_floor
    integer, intent(in) :: lower_bound
    integer :: i_low, i_high, i_mid

    i_low=lower_bound
    i_high=Num_gam_e+1
    do while (i_low < i_high)
        i_mid=(i_low+i_high)/2
        if (gam_e(i_mid) <= gamma_floor) then
            i_low=i_mid+1
        else
            i_high=i_mid
        end if
    end do
    first_gamma_gt=i_low
end function first_gamma_gt

integer function first_gamma_ge(gamma_floor, lower_bound)
    implicit none
    real(8), intent(in) :: gamma_floor
    integer, intent(in) :: lower_bound
    integer :: i_low, i_high, i_mid

    i_low=lower_bound
    i_high=Num_gam_e+1
    do while (i_low < i_high)
        i_mid=(i_low+i_high)/2
        if (gam_e(i_mid) < gamma_floor) then
            i_low=i_mid+1
        else
            i_high=i_mid
        end if
    end do
    first_gamma_ge=i_low
end function first_gamma_ge

subroutine prepare_uniform_weights_and_tails()
    implicit none
    integer :: i_shell, i_gamma

    call compute_simpson_weights(simpson_weights, Num_gam_e)
    call compute_simpson_weights(V_weights, Num_nu)
    do i_shell=1,Num_R
        weighted_dN_over_gam(:,i_shell)=dN_gam_e(:,i_shell)*simpson_weights*inv_gam
        weighted_seed(:,i_shell)=seed(:,i_shell)*V_weights
        tail_weighted_dN(Num_gam_e+1,i_shell)=zero
        tail_weighted_dN_inv2(Num_gam_e+1,i_shell)=zero
        do i_gamma=Num_gam_e,1,-1
            tail_weighted_dN(i_gamma,i_shell)=tail_weighted_dN(i_gamma+1,i_shell)+ &
                                              weighted_dN_over_gam(i_gamma,i_shell)
            tail_weighted_dN_inv2(i_gamma,i_shell)=tail_weighted_dN_inv2(i_gamma+1,i_shell)+ &
                                                   weighted_dN_over_gam(i_gamma,i_shell)*inv_gam2(i_gamma)
        end do
    end do
end subroutine prepare_uniform_weights_and_tails

subroutine prepare_kn_tables()
    implicit none
    integer :: i_seed, i_gamma
    real(8) :: temp_norm, q_gamma

    do i_seed=1,Num_nu
        gamma_start(i_seed)=first_gamma_gt(E_seed(i_seed), 1)
        do i_gamma=gamma_start(i_seed),Num_gam_e
            temp_norm=gam_e(i_gamma)-E_seed(i_seed)
            q_prefactor(i_gamma,i_seed)=V_seed(i_seed)/(4.0d0*gam_e(i_gamma)*temp_norm)
            q_gamma=E_seed(i_seed)/temp_norm
            kn_prefactor(i_gamma,i_seed)=q_gamma*q_gamma/(two*(one+q_gamma))
            log_q_prefactor(i_gamma,i_seed)=log(q_prefactor(i_gamma,i_seed))
        end do
    end do
end subroutine prepare_kn_tables

subroutine prepare_uniform_gamma_bounds()
    implicit none
    integer :: i_obs, i_seed
    real(8) :: gamma_floor

    do i_obs=1,Num_nu
        do i_seed=1,i_obs-1
            gamma_floor=0.5d0*(E_seed(i_obs)+sqrt(E_seed(i_obs)*E_seed(i_obs)+ &
                         V_seed(i_obs)*inv_v_seed(i_seed)))
            gamma_high(i_seed,i_obs)=first_gamma_gt(gamma_floor, gamma_start(i_obs))
        end do
        do i_seed=i_obs,Num_nu
            gamma_floor=0.5d0*sqrt(V_seed(i_seed)/V_seed(i_obs))
            gamma_low(i_seed,i_obs)=first_gamma_ge(gamma_floor, 1)
        end do
    end do
end subroutine prepare_uniform_gamma_bounds

real(8) function low_seed_sum(i_shell, i_obs)
    implicit none
    integer, intent(in) :: i_shell, i_obs
    integer :: i_seed, i_gamma
    real(8) :: q_coeff, q, log_q, fssc, simpson_sum_gam, emission_int2

    low_seed_sum=zero
    do i_seed=1,i_obs-1
        simpson_sum_gam=zero
        !$OMP SIMD REDUCTION(+:simpson_sum_gam)
        do i_gamma=gamma_high(i_seed,i_obs),Num_gam_e
            q_coeff=q_prefactor(i_gamma,i_obs)
            q=q_coeff*inv_v_seed(i_seed)
            if (q >= one) cycle
            log_q=log_q_prefactor(i_gamma,i_obs)+log_inv_v_seed(i_seed)
            fssc=two*q*(log_q-q)+one+q+kn_prefactor(i_gamma,i_obs)*(one-q)
            simpson_sum_gam=simpson_sum_gam+weighted_dN_over_gam(i_gamma,i_shell)*fssc
        end do
        emission_int2=h_gam_third*simpson_sum_gam
        low_seed_sum=low_seed_sum+weighted_seed(i_seed,i_shell)*emission_int2
    end do
end function low_seed_sum

real(8) function high_seed_tail_sum(i_shell, i_obs)
    implicit none
    integer, intent(in) :: i_shell, i_obs
    integer :: i_seed, i_gamma
    real(8) :: ratio_v, simpson_sum_gam, emission_int2

    high_seed_tail_sum=zero
    do i_seed=i_obs,Num_nu
        ratio_v=V_seed(i_obs)*inv_v_seed(i_seed)
        i_gamma=gamma_low(i_seed,i_obs)
        if (i_gamma <= Num_gam_e) then
            simpson_sum_gam=ratio_v*tail_weighted_dN(i_gamma,i_shell) - &
                            0.25d0*tail_weighted_dN_inv2(i_gamma,i_shell)
        else
            simpson_sum_gam=zero
        end if
        emission_int2=h_gam_third*simpson_sum_gam
        high_seed_tail_sum=high_seed_tail_sum+weighted_seed(i_seed,i_shell)*emission_int2
    end do
end function high_seed_tail_sum

subroutine accumulate_uniform_point(i_shell, i_obs)
    implicit none
    integer, intent(in) :: i_shell, i_obs
    real(8) :: dInteg, P_v, F1

    if (gamma_start(i_obs) > Num_gam_e) return

    dInteg=h_nu_third*(low_seed_sum(i_shell, i_obs)+high_seed_tail_sum(i_shell, i_obs))
    P_v=dInteg*V_seed(i_obs)
    P_SSC_spec(i_obs,i_shell)=P_SSC_spec(i_obs,i_shell)+P_v
    F1=dInteg*radius_inv2(i_shell)
    seed_SSC(i_obs,i_shell)=seed_SSC(i_obs,i_shell)+F1
end subroutine accumulate_uniform_point

end subroutine ssc_spec

! 计算非均匀网格SSC辐射谱：PPM重构电子谱 + Gauss-Legendre种子积分 + minmod斜率限制器。
subroutine ssc_spec_nonuniform(R,x_edge_log10,dN_x,V_seed,seed,Num_nu,Num_R,Num_gam_e,n_threads, &
                               P_SSC_spec,seed_SSC)
    use constants
    use radiation_common
    !$ use omp_lib
    implicit REAL(8)(A-H,O-Z)
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,n_threads
    real(8), intent(in) :: R(Num_R),x_edge_log10(Num_gam_e+1,Num_R),dN_x(Num_gam_e,Num_R)
    real(8), intent(in) :: V_seed(Num_nu),seed(Num_nu,Num_R)
    real(8), intent(out) :: P_SSC_spec(Num_nu,Num_R),seed_SSC(Num_nu,Num_R)

    real(8), allocatable :: x_seed(:),radius_inv2(:)
    real(8), allocatable :: x_center(:,:),slope_q(:,:),tail_gamma(:,:),tail_gamma_inv2(:,:)
    real(8) :: w2
    integer :: I_R,I_nu

    allocate(x_seed(Num_nu),radius_inv2(Num_R),x_center(Num_gam_e,Num_R),slope_q(Num_gam_e,Num_R))
    allocate(tail_gamma(Num_gam_e+1,Num_R),tail_gamma_inv2(Num_gam_e+1,Num_R))

    call validate_nonuniform_inputs()

    para_hEme=Para_h/para_m_energy
    x_seed=dlog(V_seed)
    radius_inv2=one/(R*R)
    w2=one/dsqrt(3d0)

    call build_nonuniform_reconstruction_and_tails()

    P_SSC_spec=zero
    seed_SSC=zero

    !$OMP PARALLEL num_threads(n_threads), private(I_R,I_nu)
    !$OMP DO COLLAPSE(2) SCHEDULE(GUIDED,4)
    do I_R=1,Num_R
        do I_nu=1,Num_nu
            call accumulate_nonuniform_point(I_R, I_nu)
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

subroutine validate_nonuniform_inputs()
    implicit none
    integer :: i_seed

    if (any(R <= zero)) error stop "ssc_spec_nonuniform: radius must be positive."
    if (any(V_seed <= zero)) error stop "ssc_spec_nonuniform: seed frequency grid must be positive."
    if (any(seed < zero)) error stop "ssc_spec_nonuniform: seed photon field must be non-negative."
    if (any(dN_x < zero)) error stop "ssc_spec_nonuniform: electron spectrum must be non-negative."
    do i_seed=1,Num_nu-1
        if (V_seed(i_seed+1) <= V_seed(i_seed)) then
            error stop "ssc_spec_nonuniform: seed frequency grid must be strictly increasing."
        end if
    end do
end subroutine validate_nonuniform_inputs

subroutine build_nonuniform_reconstruction_and_tails()
    implicit none
    integer :: i_shell, i_gamma

    do i_shell=1,Num_R
        call build_shell_reconstruction(i_shell)
        tail_gamma(Num_gam_e+1,i_shell)=zero
        tail_gamma_inv2(Num_gam_e+1,i_shell)=zero
        do i_gamma=Num_gam_e,1,-1
            call add_shell_tail_cell(i_shell, i_gamma)
        end do
    end do
end subroutine build_nonuniform_reconstruction_and_tails

subroutine build_shell_reconstruction(i_shell)
    implicit none
    integer, intent(in) :: i_shell
    integer :: i_gamma

    do i_gamma=1,Num_gam_e
        if (x_edge_log10(i_gamma+1,i_shell) <= x_edge_log10(i_gamma,i_shell)) then
            error stop "ssc_spec_nonuniform: electron edge grid must be strictly increasing in every shell."
        end if
    end do
    do i_gamma=1,Num_gam_e
        x_center(i_gamma,i_shell)=0.5d0*(x_edge_log10(i_gamma,i_shell)+x_edge_log10(i_gamma+1,i_shell))
    end do
    do i_gamma=1,Num_gam_e
        call set_limited_cell_slope(i_shell, i_gamma)
    end do
end subroutine build_shell_reconstruction

subroutine set_limited_cell_slope(i_shell, i_gamma)
    implicit none
    integer, intent(in) :: i_shell, i_gamma
    real(8) :: dx_log10, left_slope, right_slope, slope_lim

    dx_log10=x_edge_log10(i_gamma+1,i_shell)-x_edge_log10(i_gamma,i_shell)
    if (Num_gam_e == 1) then
        slope_q(i_gamma,i_shell)=zero
    else if (i_gamma == 1) then
        right_slope=(dN_x(2,i_shell)-dN_x(1,i_shell))/(x_center(2,i_shell)-x_center(1,i_shell))
        slope_q(i_gamma,i_shell)=right_slope
    else if (i_gamma == Num_gam_e) then
        left_slope=(dN_x(Num_gam_e,i_shell)-dN_x(Num_gam_e-1,i_shell))/ &
                   (x_center(Num_gam_e,i_shell)-x_center(Num_gam_e-1,i_shell))
        slope_q(i_gamma,i_shell)=left_slope
    else
        left_slope=(dN_x(i_gamma,i_shell)-dN_x(i_gamma-1,i_shell))/ &
                   (x_center(i_gamma,i_shell)-x_center(i_gamma-1,i_shell))
        right_slope=(dN_x(i_gamma+1,i_shell)-dN_x(i_gamma,i_shell))/ &
                    (x_center(i_gamma+1,i_shell)-x_center(i_gamma,i_shell))
        slope_q(i_gamma,i_shell)=ssc_minmod(left_slope,right_slope)
    end if
    slope_lim=two*dN_x(i_gamma,i_shell)/dx_log10
    if (abs(slope_q(i_gamma,i_shell)) > slope_lim) then
        slope_q(i_gamma,i_shell)=sign(slope_lim,slope_q(i_gamma,i_shell))
    end if
end subroutine set_limited_cell_slope

subroutine add_shell_tail_cell(i_shell, i_gamma)
    implicit none
    integer, intent(in) :: i_shell, i_gamma
    real(8) :: x0_cell, x1_cell

    x0_cell=x_edge_log10(i_gamma,i_shell)
    x1_cell=x_edge_log10(i_gamma+1,i_shell)
    tail_gamma(i_gamma,i_shell)=tail_gamma(i_gamma+1,i_shell)+ &
        linear_gamma_moment(x0_cell,x1_cell,dN_x(i_gamma,i_shell),slope_q(i_gamma,i_shell), &
                            x_center(i_gamma,i_shell),2d0)
    tail_gamma_inv2(i_gamma,i_shell)=tail_gamma_inv2(i_gamma+1,i_shell)+ &
        linear_gamma_moment(x0_cell,x1_cell,dN_x(i_gamma,i_shell),slope_q(i_gamma,i_shell), &
                            x_center(i_gamma,i_shell),4d0)
end subroutine add_shell_tail_cell

real(8) function gauss2_abscissa(xm, dx, i_gauss)
    implicit none
    real(8), intent(in) :: xm, dx
    integer, intent(in) :: i_gauss

    if (i_gauss == 1) then
        gauss2_abscissa=xm-dx*w2
    else
        gauss2_abscissa=xm+dx*w2
    end if
end function gauss2_abscissa

subroutine accumulate_nonuniform_point(i_shell, i_obs)
    implicit none
    integer, intent(in) :: i_shell, i_obs
    integer :: i_seed, i_gauss
    real(8) :: ratio_v, gamma_floor, dInteg, emission_int2, P_v, F1, Vloc, Ephoton2eV
    real(8) :: x_seed0, x_seed1, xm_seed, dx_seed, x_seed_loc, V_seed_loc
    real(8) :: seed_loc, seed_weight

    Vloc=V_seed(i_obs)
    Ephoton2eV=Vloc*para_hEme
    dInteg=zero
    do i_seed=1,Num_nu-1
        x_seed0=x_seed(i_seed)
        x_seed1=x_seed(i_seed+1)
        xm_seed=0.5d0*(x_seed0+x_seed1)
        dx_seed=0.5d0*(x_seed1-x_seed0)
        do i_gauss=1,2
            x_seed_loc=gauss2_abscissa(xm_seed, dx_seed, i_gauss)
            V_seed_loc=dexp(x_seed_loc)
            seed_loc=radiation_powerlaw_interp(V_seed(i_seed),V_seed(i_seed+1),seed(i_seed,i_shell), &
                                                seed(i_seed+1,i_shell),V_seed_loc)
            if (seed_loc <= zero) cycle
            seed_weight=dx_seed*seed_loc
            if (V_seed_loc < Vloc) then
                gamma_floor=0.5d0*(Ephoton2eV+dsqrt(Ephoton2eV*Ephoton2eV+Vloc/V_seed_loc))
                emission_int2=ssc_low_gamma_integral(i_shell,dlog10(gamma_floor),Vloc,V_seed_loc,Ephoton2eV)
            else
                ratio_v=Vloc/V_seed_loc
                gamma_floor=0.5d0*dsqrt(V_seed_loc/Vloc)
                emission_int2=ssc_high_gamma_tail(i_shell,dlog10(gamma_floor),ratio_v)
            end if
            dInteg=dInteg+seed_weight*emission_int2
        end do
    end do
    P_v=dInteg*Vloc
    P_SSC_spec(i_obs,i_shell)=P_v
    F1=dInteg*radius_inv2(i_shell)
    seed_SSC(i_obs,i_shell)=F1
end subroutine accumulate_nonuniform_point

! 二分查找第一个网格单元，其右边界 > x_floor。
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

! minmod斜率限制器：同号取绝对值最小者，异号返回零。
real(8) function ssc_minmod(a,b)
    implicit REAL(8)(A-H,O-Z)
    real(8), intent(in) :: a,b

    if (a*b <= zero) then
        ssc_minmod=zero
    else
        ssc_minmod=sign(min(abs(a),abs(b)),a)
    end if
end function ssc_minmod

! 计算线性重构剖面在x点的值：q(x) = q̄ + slope*(x - xc)，截断负值为零。
real(8) function linear_profile_value(x,qbar,slope,xc)
    implicit REAL(8)(A-H,O-Z)
    real(8), intent(in) :: x,qbar,slope,xc
    linear_profile_value=qbar+slope*(x-xc)
    if (linear_profile_value < zero) then
        error stop "ssc_spec_nonuniform: linear reconstruction became negative."
    end if
end function linear_profile_value

! 计算线性重构剖面在[x0,x1]上的γ^(-power)矩积分（解析公式）。
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

! 低能种子区（ν_seed < ν_obs）单个网格单元的SSC散射积分。
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

! 低能种子区（ν_seed < ν_obs）的完整SSC散射积分，遍历电子能格。
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
        x0=x_edge_log10(i_game,I_R)
        if (x_floor > x0) x0=x_floor
        x1=x_edge_log10(i_game+1,I_R)
        ssc_low_gamma_integral=ssc_low_gamma_integral+ &
            ssc_low_gamma_cell(x0,x1,dN_x(i_game,I_R),slope_q(i_game,I_R),x_center(i_game,I_R), &
                               Vloc,V_seed_loc,Ephoton2eV)
    end do
end function ssc_low_gamma_integral

! 高能种子区（ν_seed ≥ ν_obs）SSC尾部积分，利用预计算的gamma矩加速。
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

    x0=x_edge_log10(i_start,I_R)
    if (x_floor > x0) x0=x_floor
    x1=x_edge_log10(i_start+1,I_R)
    part2=linear_gamma_moment(x0,x1,dN_x(i_start,I_R),slope_q(i_start,I_R),x_center(i_start,I_R),2d0)
    part4=linear_gamma_moment(x0,x1,dN_x(i_start,I_R),slope_q(i_start,I_R),x_center(i_start,I_R),4d0)
    ssc_high_gamma_tail=ratio_v*(part2+tail_gamma(i_start+1,I_R))- &
                        0.25d0*(part4+tail_gamma_inv2(i_start+1,I_R))
end function ssc_high_gamma_tail

end subroutine ssc_spec_nonuniform
