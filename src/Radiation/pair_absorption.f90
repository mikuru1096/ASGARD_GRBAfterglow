! 计算 gamma-gamma 湮灭吸收因子。
! Compute gamma-gamma annihilation absorption.
!
! 对各向异性光子场积分对产生截面，并包含多普勒角度权重。
! Integrate the pair-production cross section over an anisotropic photon field
! with Doppler angular weights.
subroutine annihilation(R_gamma,R,V_seed,seed_syn,seed_ssc,tau_extra,Num_nu,Num_R,n_threads, absorption)
    !$ use omp_lib
    use constants
    use rad_common
    implicit none
    integer, intent(in) :: Num_R,Num_nu,n_threads
    real(8), dimension(Num_R), intent(in) :: R_gamma,R
    real(8), dimension(Num_nu), intent(in) :: V_seed
    real(8), dimension(Num_nu,Num_R), intent(in) :: seed_syn,seed_ssc,tau_extra
    real(8), dimension(Num_nu,Num_R), intent(out) :: absorption

    real(8), dimension(:), allocatable :: dVloc,V_mid,geom,beta,mu,z
    real(8), dimension(:,:), allocatable :: seedsum,seedmid,ep1,ep2,ep2ep1,seedw,angw
    real(8), dimension(:,:,:), allocatable :: sigker
    integer, dimension(:,:), allocatable :: nu_start,nu_stop
    integer :: ncos,i_r,n1,i_cos
    real(8) :: dmu,area,tau,tau1,trans

    allocate(seedsum(Num_nu,Num_R),geom(Num_R),ep1(1,Num_nu),ep2(Num_nu-1,1),ep2ep1(Num_nu-1,Num_nu), &
             dVloc(Num_nu-1),V_mid(Num_nu-1),seedmid(Num_nu-1,Num_R),beta(Num_R),seedw(Num_nu-1,Num_R))

    absorption=0d0

    ncos=50
    dmu=2d0/ncos
    area=3d0/16d0*Para_sigmaT
    allocate(mu(ncos+1),z(ncos+1),angw(ncos+1,Num_R))
    allocate(sigker(Num_nu-1,Num_nu,ncos+1))
    allocate(nu_start(ncos+1,Num_nu),nu_stop(ncos+1,Num_nu))
    absorption=0d0
    seedsum=0d0
    sigker=0d0

    seedsum=seed_syn+seed_ssc

    geom=area*R/(12d0*R_gamma)
    beta=dsqrt(1d0-R_gamma**(-2))

    call pair_grid(V_seed,Num_nu,ep1,ep2,dVloc,V_mid)
    call mid_seed()
    ep2ep1=matmul(ep2,ep1)
    call angle_weights()
    call build_sigma()

    !$OMP PARALLEL num_threads(n_threads), private(i_r, n1, tau, i_cos, tau1, trans)
    !$OMP DO
    do i_r=1,Num_R
        seedw(:,i_r)=seedmid(:,i_r)*dVloc
    end do
    !$OMP END DO

    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
    do i_r=1,Num_R
        do n1=1,Num_nu
            tau=0d0
            do i_cos=1,ncos+1
                if (nu_stop(i_cos,n1) <= 0) cycle
                tau1=dot_product( &
                    sigker(nu_start(i_cos,n1):nu_stop(i_cos,n1),n1,i_cos), &
                    seedw(nu_start(i_cos,n1):nu_stop(i_cos,n1),i_r) &
                )
                tau=tau+tau1*angw(i_cos,i_r)
            end do
            tau=tau*geom(i_r)/2d0 + tau_extra(n1,i_r)
            call transfer_factor(tau,trans)
            absorption(n1,i_r)=absorption(n1,i_r)+trans
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(seedsum,geom,ep1,ep2,ep2ep1,dVloc,V_mid,seedmid,beta,seedw)
    deallocate(mu,z,angw,sigker,nu_start,nu_stop)

    return

contains

subroutine mid_seed()
    implicit none
    integer :: n2,i_r

    do n2=1,Num_nu-1
        do i_r=1,Num_R
            seedmid(n2,i_r)=rad_interp(V_seed(n2),V_seed(n2+1),seedsum(n2,i_r),seedsum(n2+1,i_r),V_mid(n2))
        end do
    end do
end subroutine mid_seed

subroutine angle_weights()
    implicit none
    integer :: i_cos,i_r

    do i_cos=1,ncos+1
        mu(i_cos)=-1d0+dmu*(i_cos-1d0)
        z(i_cos)=(1d0-mu(i_cos))/2d0
    end do
    do i_r=1,Num_R
        angw(:,i_r)=dmu*(1d0-beta(i_r)*mu)
    end do
end subroutine angle_weights

subroutine build_sigma()
    implicit none
    integer :: n1,i_cos,n2
    real(8) :: s

    do n1=1,Num_nu
        do i_cos=1,ncos+1
            call set_window(i_cos,n1)
            if (nu_stop(i_cos,n1) <= 0) cycle
            do n2=nu_start(i_cos,n1),nu_stop(i_cos,n1)
                s=ep2ep1(n2,n1)*z(i_cos)
                if (s <= 1d0 .or. s >= 1d12) cycle
                sigker(n2,n1,i_cos)=pair_sigma(s)/area
            end do
        end do
    end do
end subroutine build_sigma

subroutine set_window(i_cos,n1)
    implicit none
    integer, intent(in) :: i_cos,n1
    real(8) :: lo,hi

    if (z(i_cos) <= 0d0) then
        call clear_window(i_cos,n1)
        return
    end if

    lo=1d0/z(i_cos)
    hi=1d12/z(i_cos)
    if (ep2ep1(Num_nu-1,n1) <= lo .or. ep2ep1(1,n1) >= hi) then
        call clear_window(i_cos,n1)
        return
    end if

    nu_start(i_cos,n1)=first_gt(n1,lo)
    nu_stop(i_cos,n1)=last_lt(n1,hi)
    if (nu_start(i_cos,n1) > nu_stop(i_cos,n1)) nu_stop(i_cos,n1)=0
end subroutine set_window

subroutine clear_window(i_cos,n1)
    implicit none
    integer, intent(in) :: i_cos,n1

    nu_start(i_cos,n1)=1
    nu_stop(i_cos,n1)=0
end subroutine clear_window

integer function first_gt(n1,target)
    implicit none
    integer, intent(in) :: n1
    real(8), intent(in) :: target
    integer :: i_low,i_high,i_mid

    i_low=1
    i_high=Num_nu-1
    do while (i_low < i_high)
        i_mid=(i_low+i_high)/2
        if (ep2ep1(i_mid,n1) <= target) then
            i_low=i_mid+1
        else
            i_high=i_mid
        end if
    end do
    first_gt=i_low
end function first_gt

integer function last_lt(n1,target)
    implicit none
    integer, intent(in) :: n1
    real(8), intent(in) :: target
    integer :: i_low,i_high,i_mid

    i_low=1
    i_high=Num_nu-1
    do while (i_low < i_high)
        i_mid=(i_low+i_high+1)/2
        if (ep2ep1(i_mid,n1) < target) then
            i_low=i_mid
        else
            i_high=i_mid-1
        end if
    end do
    last_lt=i_low
end function last_lt
end subroutine annihilation
