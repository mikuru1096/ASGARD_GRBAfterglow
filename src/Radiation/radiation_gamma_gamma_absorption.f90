! Created by  on 2021/1/31.
! 计算γ-γ湮灭吸收因子：对各向异性光子场积分对产生截面，含多普勒角度权重。
subroutine annihilation(R_gamma,R,V_seed,seed_syn,seed_ssc,tau_extra,Num_nu,Num_R,n_threads, absorption)
    !$ use omp_lib
    use constants
    use radiation_common
    implicit none
    integer, intent(in) :: Num_R,Num_nu,n_threads
    real(8), intent(in) :: R_gamma(Num_R),R(Num_R),V_seed(Num_nu),seed_syn(Num_nu,Num_R),seed_ssc(Num_nu,Num_R), &
                           tau_extra(Num_nu,Num_R)
    real(8), intent(out) :: absorption(Num_nu,Num_R)

    real(8), allocatable :: seed_tot(:,:),seed_tot_mid(:,:),ep1(:,:),ep2(:,:),ep2ep1(:,:),dVloc(:),V_mid(:), &
                            dRariv_Sigma(:),beta(:),seed_tot_mid_dVloc(:,:),dcos_grid(:),z_grid(:), &
                            cos_weight(:,:),sigma_kernel(:,:,:)
    integer, allocatable :: nu2_start(:,:),nu2_stop(:,:)
    integer :: Num_cos,I_R,Nu_s1,i_cos
    real(8) :: dcos_bin,Cross_Area,Tau,Tau1,temp_abs
                
    allocate (seed_tot(Num_nu,Num_R),dRariv_Sigma(Num_R),ep1(1,Num_nu),ep2(Num_nu-1,1),ep2ep1(Num_nu-1,Num_nu), &
            dVloc(Num_nu-1),V_mid(Num_nu-1),seed_tot_mid(Num_nu-1,Num_R),beta(Num_R), &
            seed_tot_mid_dVloc(Num_nu-1,Num_R))

    absorption=zero

    Num_cos=50
    dcos_bin=two/Num_cos
    Cross_Area=3.0d0/16.0d0*Para_sigmaT
    allocate (dcos_grid(Num_cos+1),z_grid(Num_cos+1),cos_weight(Num_cos+1,Num_R))
    allocate (sigma_kernel(Num_nu-1,Num_nu,Num_cos+1))
    allocate (nu2_start(Num_cos+1,Num_nu),nu2_stop(Num_cos+1,Num_nu))
    absorption=zero
    seed_tot=zero
    sigma_kernel=zero

    seed_tot=seed_syn+seed_ssc
    
    dRariv_Sigma=Cross_Area*R/(12.0d0*R_gamma)
    beta=dsqrt(one-R_gamma**(-2))

    call radiation_prepare_annihilation_grid(V_seed,Num_nu,ep1,ep2,dVloc,V_mid)
    call compute_mid_seed_field()
    ep2ep1=matmul(ep2,ep1)
    call prepare_pair_angle_weights()
    call build_pair_sigma_kernel()
    
    !$OMP PARALLEL num_threads(n_threads), private(I_R, Nu_s1, Tau, i_cos, Tau1, temp_abs)
    !$OMP DO
    do i_R=1,Num_R
        seed_tot_mid_dVloc(:,i_R)=seed_tot_mid(:,i_R)*dVloc
    end do
    !$OMP END DO
    
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
    do I_R=1,Num_R
        do Nu_s1=1,Num_nu
            Tau=zero
            do i_cos=1,Num_cos+1
                if (nu2_stop(i_cos,Nu_s1) <= 0) cycle
                Tau1=dot_product( &
                    sigma_kernel(nu2_start(i_cos,Nu_s1):nu2_stop(i_cos,Nu_s1),Nu_s1,i_cos), &
                    seed_tot_mid_dVloc(nu2_start(i_cos,Nu_s1):nu2_stop(i_cos,Nu_s1),I_R) &
                )
                Tau=Tau+Tau1*cos_weight(i_cos,I_R)
            end do
            Tau=Tau*dRariv_Sigma(I_R)/two + tau_extra(Nu_s1,I_R)
            call radiation_transfer_factor(Tau,temp_abs)
            absorption(Nu_s1,I_R)=absorption(Nu_s1,I_R)+temp_abs
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(seed_tot,dRariv_Sigma,ep1,ep2,ep2ep1,dVloc,V_mid,seed_tot_mid,beta,seed_tot_mid_dVloc)
    deallocate(dcos_grid,z_grid,cos_weight,sigma_kernel,nu2_start,nu2_stop)

    return

contains

subroutine compute_mid_seed_field()
    implicit none
    integer :: Nu_s2,i_R

    do Nu_s2=1,Num_nu-1
        do i_R=1,Num_R
            seed_tot_mid(Nu_s2,i_R)=radiation_powerlaw_interp(V_seed(Nu_s2),V_seed(Nu_s2+1), &
                                                              seed_tot(Nu_s2,i_R),seed_tot(Nu_s2+1,i_R), &
                                                              V_mid(Nu_s2))
        end do
    end do
end subroutine compute_mid_seed_field

subroutine prepare_pair_angle_weights()
    implicit none
    integer :: i_cos,i_R

    do i_cos=1,Num_cos+1
        dcos_grid(i_cos)=-one+dcos_bin*(i_cos-one)
        z_grid(i_cos)=(one-dcos_grid(i_cos))/two
    end do
    do i_R=1,Num_R
        cos_weight(:,i_R)=dcos_bin*(one-beta(i_R)*dcos_grid)
    end do
end subroutine prepare_pair_angle_weights

subroutine build_pair_sigma_kernel()
    implicit none
    integer :: Nu_s1,i_cos,Nu_s2
    real(8) :: s_center

    do Nu_s1=1,Num_nu
        do i_cos=1,Num_cos+1
            call set_pair_energy_window(i_cos,Nu_s1)
            if (nu2_stop(i_cos,Nu_s1) <= 0) cycle
            do Nu_s2=nu2_start(i_cos,Nu_s1),nu2_stop(i_cos,Nu_s1)
                s_center=ep2ep1(Nu_s2,Nu_s1)*z_grid(i_cos)
                if (s_center <= one .or. s_center >= 1.0d12) cycle
                sigma_kernel(Nu_s2,Nu_s1,i_cos)=radiation_pair_cross_section(s_center)/Cross_Area
            end do
        end do
    end do
end subroutine build_pair_sigma_kernel

subroutine set_pair_energy_window(i_cos,Nu_s1)
    implicit none
    integer, intent(in) :: i_cos,Nu_s1
    real(8) :: lo_target,hi_target

    if (z_grid(i_cos) <= zero) then
        call clear_pair_energy_window(i_cos,Nu_s1)
        return
    end if

    lo_target=one/z_grid(i_cos)
    hi_target=1.0d12/z_grid(i_cos)
    if (ep2ep1(Num_nu-1,Nu_s1) <= lo_target .or. ep2ep1(1,Nu_s1) >= hi_target) then
        call clear_pair_energy_window(i_cos,Nu_s1)
        return
    end if

    nu2_start(i_cos,Nu_s1)=first_pair_index_gt(Nu_s1,lo_target)
    nu2_stop(i_cos,Nu_s1)=last_pair_index_lt(Nu_s1,hi_target)
    if (nu2_start(i_cos,Nu_s1) > nu2_stop(i_cos,Nu_s1)) nu2_stop(i_cos,Nu_s1)=0
end subroutine set_pair_energy_window

subroutine clear_pair_energy_window(i_cos,Nu_s1)
    implicit none
    integer, intent(in) :: i_cos,Nu_s1

    nu2_start(i_cos,Nu_s1)=1
    nu2_stop(i_cos,Nu_s1)=0
end subroutine clear_pair_energy_window

integer function first_pair_index_gt(Nu_s1,target)
    implicit none
    integer, intent(in) :: Nu_s1
    real(8), intent(in) :: target
    integer :: i_low,i_high,i_mid

    i_low=1
    i_high=Num_nu-1
    do while (i_low < i_high)
        i_mid=(i_low+i_high)/2
        if (ep2ep1(i_mid,Nu_s1) <= target) then
            i_low=i_mid+1
        else
            i_high=i_mid
        end if
    end do
    first_pair_index_gt=i_low
end function first_pair_index_gt

integer function last_pair_index_lt(Nu_s1,target)
    implicit none
    integer, intent(in) :: Nu_s1
    real(8), intent(in) :: target
    integer :: i_low,i_high,i_mid

    i_low=1
    i_high=Num_nu-1
    do while (i_low < i_high)
        i_mid=(i_low+i_high+1)/2
        if (ep2ep1(i_mid,Nu_s1) < target) then
            i_low=i_mid
        else
            i_high=i_mid-1
        end if
    end do
    last_pair_index_lt=i_low
end function last_pair_index_lt
end subroutine annihilation
