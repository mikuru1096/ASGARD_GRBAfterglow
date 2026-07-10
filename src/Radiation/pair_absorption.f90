! 计算 gamma-gamma 湮灭吸收因子。
! Compute gamma-gamma annihilation absorption.
!
! 对局域共动各向同性光子场积分 Breit-Wheeler 截面和 Gould-Schreder 相对速度权重。
! Integrate the Breit-Wheeler cross section over a local comoving isotropic photon field
! with the Gould-Schreder relative-speed weight.
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

    integer, parameter :: nq=16
    real(8), parameter :: tq(nq)=[ &
        0.005299532504175031d0,0.027712488463383700d0,0.067184398806084122d0,0.122297795822498501d0, &
        0.191061877798678115d0,0.270991611171386315d0,0.359198224610370542d0,0.452493745081181287d0, &
        0.547506254918818769d0,0.640801775389629458d0,0.729008388828613629d0,0.808938122201321885d0, &
        0.877702204177501555d0,0.932815601193915930d0,0.972287511536616300d0,0.994700467495824969d0]
    real(8), parameter :: wq(nq)=[ &
        0.013576229705877088d0,0.031126761969323728d0,0.047579255841246303d0,0.062314485627767036d0, &
        0.074797994408288354d0,0.084578259697501323d0,0.091301707522461820d0,0.094725305227534320d0, &
        0.094725305227534320d0,0.091301707522461820d0,0.084578259697501323d0,0.074797994408288354d0, &
        0.062314485627767036d0,0.047579255841246303d0,0.031126761969323728d0,0.013576229705877088d0]
    real(8), dimension(:), allocatable :: dVloc,V_mid,geom
    real(8), dimension(:,:), allocatable :: ep1,ep2,seedw,kerrel
    integer, dimension(:), allocatable :: nu_start
    integer :: i_r,n1
    real(8) :: tau,trans

    ! 零局域靶场不产生额外 gamma-gamma 光深，只保留 tau_extra 壳层转移。
    ! A zero local target adds no pair opacity, leaving only the tau_extra shell transfer.
    if (.not. any(seed_syn /= 0d0) .and. .not. any(seed_ssc /= 0d0)) then
        !$OMP PARALLEL DO COLLAPSE(2) num_threads(n_threads) SCHEDULE(STATIC) private(i_r,n1)
        do i_r=1,Num_R
            do n1=1,Num_nu
                call transfer_factor(tau_extra(n1,i_r),absorption(n1,i_r))
            end do
        end do
        !$OMP END PARALLEL DO
        return
    end if

    allocate(geom(Num_R),ep1(1,Num_nu),ep2(Num_nu-1,1),dVloc(Num_nu-1),V_mid(Num_nu-1), &
             seedw(Num_nu-1,Num_R))

    allocate(kerrel(Num_nu-1,Num_nu))
    allocate(nu_start(Num_nu))

    geom=R/(12d0*R_gamma)

    call pair_grid(V_seed,Num_nu,ep1,ep2,dVloc,V_mid)
    call weight_seed()
    call build_moment()

    !$OMP PARALLEL DO COLLAPSE(2) num_threads(n_threads) SCHEDULE(STATIC) &
    !$OMP& private(i_r,n1,tau,trans)
    do i_r=1,Num_R
        do n1=1,Num_nu
            tau=0d0
            if (nu_start(n1) < Num_nu) then
                tau=dot_product(kerrel(nu_start(n1):Num_nu-1,n1), &
                                seedw(nu_start(n1):Num_nu-1,i_r))*geom(i_r)
            end if
            tau=tau+tau_extra(n1,i_r)
            call transfer_factor(tau,trans)
            absorption(n1,i_r)=trans
        end do
    end do
    !$OMP END PARALLEL DO

contains

subroutine weight_seed()
    implicit none
    integer :: n2,i_r

    do n2=1,Num_nu-1
        do i_r=1,Num_R
            seedw(n2,i_r)=dVloc(n2)*rad_interp(V_seed(n2),V_seed(n2+1), &
                seed_syn(n2,i_r)+seed_ssc(n2,i_r),seed_syn(n2+1,i_r)+seed_ssc(n2+1,i_r),V_mid(n2))
        end do
    end do
end subroutine weight_seed

! 各向同性靶场的 (1-mu)/2 角平均只需一个正定二维截面矩。
! The map t^2=(s-1)/(x-1) removes the threshold square-root endpoint and integrates exact support.
subroutine build_moment()
    implicit none
    integer :: n1,n2,nlo,iq
    real(8) :: xprod,delta,ds,acc

    kerrel=0d0
    nu_start=Num_nu
    !$OMP PARALLEL DO num_threads(n_threads) SCHEDULE(STATIC) &
    !$OMP& private(n1,n2,nlo,iq,xprod,delta,ds,acc)
    do n1=1,Num_nu
        if (ep2(Num_nu-1,1)*ep1(1,n1) <= 1d0) cycle
        nlo=first_gt(n1,1d0)
        nu_start(n1)=nlo
        do n2=nlo,Num_nu-1
            xprod=ep2(n2,1)*ep1(1,n1)
            delta=xprod-1d0
            acc=0d0
            do iq=1,nq
                ds=delta*tq(iq)*tq(iq)
                acc=acc+wq(iq)*tq(iq)*(1d0+ds)*pair_delta(ds)
            end do
            kerrel(n2,n1)=4d0*delta/xprod/xprod*acc
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine build_moment

integer function first_gt(n1,target)
    implicit none
    integer, intent(in) :: n1
    real(8), intent(in) :: target
    integer :: i_low,i_high,i_mid

    i_low=1
    i_high=Num_nu-1
    do while (i_low < i_high)
        i_mid=(i_low+i_high)/2
        if (ep2(i_mid,1)*ep1(1,n1) <= target) then
            i_low=i_mid+1
        else
            i_high=i_mid
        end if
    end do
    first_gt=i_low
end function first_gt
end subroutine annihilation
