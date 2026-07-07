module electron_radiation_kernel
  use constants
  use rad_common, only: syn_seed_chi, syn_kernel, &
                              transfer_factor
  private

  integer, parameter :: index_syn_fixed_grid_legacy=1, index_syn_fixed_grid=2, index_syn_cyclotron=4

  public :: greater_window
  public :: besselk, syn_state
  public :: syn_transfer, nua_solve
  public :: nua_path
  public :: reduce_grid, project_syn
  public :: nua_fromtau
  public :: pl_interp, log_gauss2, pl_integral

contains

! 在单调递增数组arr的[left,right]区间内二分查找第一个 > target 的索引。
subroutine greater_from(arr,left,right,target,idx)
implicit none
integer, intent(in) :: left,right
real(8), intent(in), dimension(right) :: arr
real(8), intent(in) :: target
integer, intent(out) :: idx
integer :: lo,hi,mid

    lo=left
    hi=right
    do while (lo < hi)
        mid=(lo+hi)/2
        if (arr(mid) > target) then
            hi=mid
        else
            lo=mid+1
        end if
    end do
    idx=lo
end subroutine greater_from

! 在全数组arr(1:n)中二分查找第一个 > target 的索引，含边界快速返回。
! 从start_idx开始在arr(1:n)中查找第一个 > target 的索引（带游标加速）。
subroutine greater_window(arr,n,start_idx,target,idx)
implicit none
integer, intent(in) :: n,start_idx
real(8), intent(in), dimension(n) :: arr
real(8), intent(in) :: target
integer, intent(out) :: idx

    if (arr(start_idx) > target) then
        idx=start_idx
        return
    end if
    if (arr(n) <= target) then
        idx=n+1
        return
    end if

    call greater_from(arr,start_idx,n,target,idx)
end subroutine greater_window

!****************************************************************************************
!**************************** Besselk function interpolation ****************************
!****************************************************************************************
function besselk(var)
    REAL(8) :: val,var,besselk
    integer :: i
    REAL(8),dimension(100),parameter :: minus_theta= &
                 [1.0000000e-04,1.1233240e-04,1.2618569e-04,1.4174742e-04,1.5922828e-04,1.7886495e-04,2.0092330e-04,&
                 2.2570197e-04,2.5353645e-04,2.8480359e-04,3.1992671e-04,3.5938137e-04,4.0370173e-04,4.5348785e-04,&
                 5.0941380e-04,5.7223677e-04,6.4280731e-04,7.2208090e-04,8.1113083e-04,9.1116276e-04,1.0235310e-03,&
                 1.1497570e-03,1.2915497e-03,1.4508288e-03,1.6297508e-03,1.8307383e-03,2.0565123e-03,2.3101297e-03,&
                 2.5950242e-03,2.9150531e-03,3.2745492e-03,3.6783798e-03,4.1320124e-03,4.6415888e-03,5.2140083e-03,&
                 5.8570208e-03,6.5793322e-03,7.3907220e-03,8.3021757e-03,9.3260335e-03,1.0476158e-02,1.1768120e-02,&
                 1.3219411e-02,1.4849683e-02,1.6681005e-02,1.8738174e-02,2.1049041e-02,2.3644894e-02,2.6560878e-02,&
                 2.9836472e-02,3.3516027e-02,3.7649358e-02,4.2292429e-02,4.7508102e-02,5.3366992e-02,5.9948425e-02,&
                 6.7341507e-02,7.5646333e-02,8.4975344e-02,9.5454846e-02,1.0722672e-01,1.2045035e-01,1.3530478e-01,&
                 1.5199111e-01,1.7073526e-01,1.9179103e-01,2.1544347e-01,2.4201283e-01,2.7185882e-01,3.0538555e-01,&
                 3.4304693e-01,3.8535286e-01,4.3287613e-01,4.8626016e-01,5.4622772e-01,6.1359073e-01,6.8926121e-01,&
                 7.7426368e-01,8.6974900e-01,9.7700996e-01,1.0974988e+00,1.2328467e+00,1.3848864e+00,1.5556761e+00,&
                 1.7475284e+00,1.9630407e+00,2.2051307e+00,2.4770764e+00,2.7825594e+00,3.1257158e+00,3.5111917e+00,&
                 3.9442061e+00,4.4306215e+00,4.9770236e+00,5.5908102e+00,6.2802914e+00,7.0548023e+00,7.9248290e+00,&
                 8.9021509e+00,1.0000000e+01]
    REAL(8),dimension(100),parameter :: besselk2= &
                 [2.0000000e+08,1.5849658e+08,1.2560583e+08,9.9540471e+07,7.8884121e+07,6.2514316e+07,4.9541527e+07,&
                 3.9260813e+07,3.1113522e+07,2.4656934e+07,1.9540199e+07,1.5485273e+07,1.2271814e+07,9.7252027e+06,&
                 7.7070567e+06,6.1077105e+06,4.8402560e+06,3.8358200e+06,3.0398217e+06,2.4090066e+06,1.9090964e+06,&
                 1.5129262e+06,1.1989680e+06,9.5016153e+05,7.5298666e+05,5.9672895e+05,4.7289738e+05,3.7476298e+05,&
                 2.9699315e+05,2.3536189e+05,1.8652017e+05,1.4781394e+05,1.1713992e+05,9.2831277e+04,7.3567095e+04,&
                 5.8300561e+04,4.6202094e+04,3.6614266e+04,2.9016076e+04,2.2994640e+04,1.8222755e+04,1.4441118e+04,&
                 1.1444235e+04,9.0692572e+03,7.1871275e+03,5.6955719e+03,4.5135397e+03,3.5767994e+03,2.8344487e+03,&
                 2.2461486e+03,1.7799308e+03,1.4104612e+03,1.1176629e+03,8.8562540e+02,7.0173970e+02,5.5601353e+02,&
                 4.4052817e+02,3.4900815e+02,2.7648028e+02,2.1900342e+02,1.7345426e+02,1.3735766e+02,1.0875212e+02,&
                 8.6083185e+01,6.8119011e+01,5.3883384e+01,4.2602693e+01,3.3663886e+01,2.6581150e+01,2.0969514e+01,&
                 1.6523922e+01,1.3002650e+01,1.0214166e+01,8.0067134e+00,6.2600579e+00,4.8789400e+00,3.7878895e+00,&
                 2.9271104e+00,2.2492185e+00,1.7166521e+00,1.2996171e+00,9.7445621e-01,7.2235470e-01,5.2831390e-01,&
                 3.8033947e-01,2.6880181e-01,1.8593568e-01,1.2545279e-01,8.2245762e-02,5.2165420e-02,3.1855005e-02,&
                 1.8626561e-02,1.0365680e-02,5.4526367e-03,2.6905368e-03,1.2347515e-03,5.2199962e-04,2.0112029e-04,&
                 6.9778918e-05,2.1509817e-05]

   if (var <= 0d0) error stop 'besselk requires var > 0'
   if (var < minus_theta(1)) then
       besselk=2d0/var**2
       return
   end if
   if (var > minus_theta(100)) then
       besselk=dsqrt(pi/(2d0*var))*dexp(-var)*(1d0+15d0/(8d0*var)+105d0/(128d0*var**2))
       return
   end if

   besselk=0d0

   do i=1,99
       if (var >= minus_theta(i) .and. var <= minus_theta(i+1)) then
           val = (var-minus_theta(i))/(minus_theta(i+1)-minus_theta(i))
           besselk = besselk2(i)+val*(besselk2(i+1)-besselk2(i))
           exit
       end if
   end do
end function besselk


subroutine syn_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                  P_emit,P_syn,Seed_syn,Tau_syn)
implicit none
integer, intent(in) :: index_syn_intger,Num_gam_e,Num_nu,n_threads
real(8), intent(in), dimension(Num_gam_e) :: gam_e,dN_gam_e
real(8), intent(in), dimension(Num_nu) :: V_seed
real(8), intent(in) :: R_loc,DB
real(8), intent(out), dimension(Num_nu) :: P_emit,P_syn,Seed_syn,Tau_syn

    select case(index_syn_intger)
    case(index_syn_fixed_grid_legacy,index_syn_fixed_grid)
        call syn_fixed(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                      P_emit,P_syn,Seed_syn,Tau_syn)
    case(index_syn_cyclotron)
        call syn_cyclotron(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                     P_emit,P_syn,Seed_syn,Tau_syn)
    case default
        error stop 'syn_state: index_syn_intger must be 1, 2, or 4.'
    end select
end subroutine syn_state

! 固定频率网格同步辐射：Simpson 积分给发射率，有限体积端点差给 SSA 光深。
! Fixed-grid synchrotron: Simpson integration for emissivity and finite-volume endpoint drops for SSA depth.
subroutine syn_fixed(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                    P_emit,P_syn,Seed_syn,Tau_syn)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in), dimension(Num_gam_e) :: gam_e,dN_gam_e
real(8), intent(in), dimension(Num_nu) :: V_seed
real(8), intent(in) :: R_loc,DB
real(8), intent(out), dimension(Num_nu) :: P_emit,P_syn,Seed_syn,Tau_syn
real(8), dimension(Num_gam_e) :: dN1,wsim,wemit
real(8), dimension(Num_gam_e-1) :: ddN, gmean2, wtau
real(8), dimension(Num_nu) :: vpow
real(8), dimension(Num_gam_e) :: vcinv,vcpow
real(8), dimension(Num_gam_e-1) :: tauinv,taupow
integer :: I_nu

    factor=(3.62d0/pi)**2
    Temp_syn=dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    r2=R_loc*R_loc
    dN1=dN_gam_e/(gam_e*gam_e)
    ddN=dN1(1:Num_gam_e-1)-dN1(2:Num_gam_e)
    h=dlog(gam_e(2))-dlog(gam_e(1))

    ! 预计算只依赖电子网格的权重，频率循环内只做 kernel 卷积。
    ! Precompute electron-grid weights so the frequency loop only evaluates the radiation kernel.
    do I_gam_e=1,Num_gam_e
        if (I_gam_e == 1 .or. I_gam_e == Num_gam_e) then
            wsim(I_gam_e)=1d0
        else if (mod(I_gam_e,2) == 0) then
            wsim(I_gam_e)=4d0
        else
            wsim(I_gam_e)=2d0
        end if
        Vc=4.2d6*gam_e(I_gam_e)*gam_e(I_gam_e)*DB
        vcinv(I_gam_e)=1d0/Vc
        vcpow(I_gam_e)=Vc**(2d0/3d0)
        wemit(I_gam_e)=wsim(I_gam_e)*dN_gam_e(I_gam_e)*gam_e(I_gam_e)
    end do
    do I_gam_e=1,Num_gam_e-1
        gmean2(I_gam_e)=(gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0
        Vc=4.2d6*gmean2(I_gam_e)*DB
        tauinv(I_gam_e)=1d0/Vc
        taupow(I_gam_e)=Vc**(2d0/3d0)
        wtau(I_gam_e)=gmean2(I_gam_e)*ddN(I_gam_e)
    end do
    do I_nu=1,Num_nu
        vpow(I_nu)=V_seed(I_nu)**(-2d0/3d0)
    end do

    !$OMP PARALLEL num_threads(n_threads), private(I_nu)
    !$OMP DO SCHEDULE(STATIC)
    do I_nu=1,Num_nu
        call store_syn(I_nu)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    hfac=4d0*pi*Para_c*Para_h
    Seed_syn=Seed_syn/hfac

    return

contains

real(8) function emit_simpson(V_cal,V_powm23)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: V_cal,V_powm23
real(8) :: vsum,x,Fx
integer :: I_gam_e

    vsum=0d0
    do I_gam_e=1,Num_gam_e
        x=V_cal*vcinv(I_gam_e)
        Fx=syn_kernel(x,vcpow(I_gam_e)*V_powm23,factor)
        vsum=vsum+wemit(I_gam_e)*Fx
    end do
    emit_simpson=h*vsum/3d0
end function emit_simpson

real(8) function tau_simpson(V_cal,V_powm23)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: V_cal,V_powm23
real(8) :: x,Fx
integer :: I_gam_e

    tau_simpson=0d0
    do I_gam_e=1,Num_gam_e-1
        x=V_cal*tauinv(I_gam_e)
        Fx=syn_kernel(x,taupow(I_gam_e)*V_powm23,factor)
        tau_simpson=tau_simpson+wtau(I_gam_e)*Fx
    end do
end function tau_simpson

subroutine store_syn(I_nu)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: I_nu
real(8) :: V_cal,integ,Tau,pv

    V_cal=V_seed(I_nu)
    integ=emit_simpson(V_cal,vpow(I_nu))
    Tau=tau_simpson(V_cal,vpow(I_nu))
    pv=Temp_syn*DB*integ
    Tau=1.046d4*Tau*DB/(4d0*pi*r2*V_cal*V_cal)
    if ((Tau-1d-4) < 1d-5) Tau=1d-4
    P_emit(I_nu)=pv
    Tau_syn(I_nu)=Tau
    P_syn(I_nu)=pv*(1d0-dexp(-Tau))/Tau
    Seed_syn(I_nu)=P_syn(I_nu)/(r2*V_cal)
end subroutine store_syn
end subroutine syn_fixed

! 构建约化对数频率网格：在85%高频处加密采样，用于冷却计算的轻量级频率表。
subroutine reduce_grid(Num_nu_in,V_in,Num_nu_out,V_out)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu_in,Num_nu_out
real(8), intent(in), dimension(Num_nu_in) :: V_in
real(8), intent(out), dimension(Num_nu_out) :: V_out
integer :: I_nu
real(8) :: x0,x1,x_tail_start

    if (Num_nu_out <= 1) then
        V_out(1)=V_in(1)
        return
    end if

    x0=dlog(V_in(1))
    x1=dlog(V_in(Num_nu_in))
    if (Num_nu_out >= 5) then
        x_tail_start=x0+0.85d0*(x1-x0)
        do I_nu=1,Num_nu_out-1
            V_out(I_nu)=dexp(x0+(x_tail_start-x0)*dble(I_nu-1)/dble(Num_nu_out-2))
        end do
        V_out(Num_nu_out)=dexp(x1)
    else
        do I_nu=1,Num_nu_out
            V_out(I_nu)=dexp(x0+(x1-x0)*dble(I_nu-1)/dble(Num_nu_out-1))
        end do
    end if
end subroutine reduce_grid

! 将同步辐射状态(P,Seed,Tau)从原频率网格幂律插值投影到约化频率网格。
subroutine project_syn(Num_nu_in,V_in,P_in,Seed_in,Tau_in,Num_nu_out,V_out,P_out,Seed_out,Tau_out)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu_in,Num_nu_out
integer :: I_out, idx_lo, idx_hi
real(8), intent(in), dimension(Num_nu_in) :: V_in,P_in,Seed_in,Tau_in
real(8), intent(in), dimension(Num_nu_out) :: V_out
real(8), intent(out), dimension(Num_nu_out) :: P_out,Seed_out,Tau_out
real(8) :: V_tar

    idx_lo=1
    do I_out=1,Num_nu_out
        V_tar=min(max(V_out(I_out),V_in(1)),V_in(Num_nu_in))
        if (V_tar <= V_in(1)) then
            P_out(I_out)=max(P_in(1),0d0)
            Seed_out(I_out)=max(Seed_in(1),0d0)
            Tau_out(I_out)=max(Tau_in(1),1d-4)
            cycle
        end if
        if (V_tar >= V_in(Num_nu_in)) then
            P_out(I_out)=max(P_in(Num_nu_in),0d0)
            Seed_out(I_out)=max(Seed_in(Num_nu_in),0d0)
            Tau_out(I_out)=max(Tau_in(Num_nu_in),1d-4)
            cycle
        end if

        call greater_window(V_in,Num_nu_in,idx_lo,V_tar,idx_hi)
        idx_lo=max(1,min(idx_hi-1,Num_nu_in-1))
        P_out(I_out)=max(pl_interp(V_in(idx_lo),V_in(idx_hi),P_in(idx_lo),P_in(idx_hi),V_tar),0d0)
        Seed_out(I_out)=max(pl_interp(V_in(idx_lo),V_in(idx_hi),Seed_in(idx_lo),Seed_in(idx_hi),V_tar),0d0)
        Tau_out(I_out)=max(pl_interp(V_in(idx_lo),V_in(idx_hi),Tau_in(idx_lo),Tau_in(idx_hi),V_tar),1d-4)
    end do
end subroutine project_syn

! 同步+非相对论回旋发射核：γ<2 的电子使用基频回旋发射，γ>=2 仍走标准同步核。
subroutine syn_cyclotron(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                   P_emit,P_syn,Seed_syn,Tau_syn)
implicit none
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in), dimension(Num_gam_e) :: gam_e,dN_gam_e
real(8), intent(in), dimension(Num_nu) :: V_seed
real(8), intent(in) :: R_loc,DB
real(8), intent(out), dimension(Num_nu) :: P_emit,P_syn,Seed_syn,Tau_syn
real(8), dimension(Num_gam_e) :: dN_syn
integer :: I_gam_e

    dN_syn=0d0
    do I_gam_e=1,Num_gam_e
        if (gam_e(I_gam_e) >= 2d0) dN_syn(I_gam_e)=dN_gam_e(I_gam_e)
    end do
    call syn_fixed(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_syn,V_seed, &
                                  P_emit,P_syn,Seed_syn,Tau_syn)
    call add_cyclotron(R_loc,DB,Num_gam_e,Num_nu,gam_e,dN_gam_e,V_seed,P_emit,P_syn,Seed_syn,Tau_syn)
end subroutine syn_cyclotron

subroutine add_cyclotron(R_loc,DB,Num_gam_e,Num_nu,gam_e,dN_gam_e,V_seed,P_emit,P_syn,Seed_syn,Tau_syn)
implicit none
integer, intent(in) :: Num_gam_e,Num_nu
integer :: I_gam_e,I_nu
real(8), intent(in), dimension(Num_gam_e) :: gam_e,dN_gam_e
real(8), intent(in), dimension(Num_nu) :: V_seed,Tau_syn
real(8), intent(in) :: R_loc,DB
real(8), intent(inout), dimension(Num_nu) :: P_emit,P_syn,Seed_syn
real(8), dimension(Num_nu+1) :: nu_edge
real(8) :: nu_b, nu0, gamma_mid, beta2, n_e_seg, p_total, p_nu, transfer, r2, hfac

    if (DB <= 0d0) return
    nu_edge(1)=V_seed(1)*dsqrt(V_seed(1)/V_seed(2))
    do I_nu=2,Num_nu
        nu_edge(I_nu)=dsqrt(V_seed(I_nu-1)*V_seed(I_nu))
    end do
    nu_edge(Num_nu+1)=V_seed(Num_nu)*dsqrt(V_seed(Num_nu)/V_seed(Num_nu-1))
    nu_b=para_e*DB/(2d0*pi*para_m_e*para_c)
    r2=R_loc*R_loc
    hfac=4d0*pi*para_c*para_h
    do I_gam_e=1,Num_gam_e-1
        gamma_mid=0.5d0*(gam_e(I_gam_e)+gam_e(I_gam_e+1))
        if (gamma_mid >= 2d0) cycle
        n_e_seg=0.5d0*(dN_gam_e(I_gam_e)+dN_gam_e(I_gam_e+1))*(gam_e(I_gam_e+1)-gam_e(I_gam_e))
        if (n_e_seg <= 0d0) cycle
        beta2=1d0-1d0/(gamma_mid*gamma_mid)
        nu0=nu_b/gamma_mid
        if (nu0 < nu_edge(1) .or. nu0 >= nu_edge(Num_nu+1)) cycle
        call greater_window(nu_edge,Num_nu+1,1,nu0,I_nu)
        I_nu=max(1,min(Num_nu,I_nu-1))
        p_total=n_e_seg*(4d0/3d0)*para_sigmat*para_c*(DB*DB/(8d0*pi))*gamma_mid*gamma_mid*beta2
        p_nu=p_total/(nu_edge(I_nu+1)-nu_edge(I_nu))
        call transfer_factor(Tau_syn(I_nu),transfer)
        P_emit(I_nu)=P_emit(I_nu)+p_nu
        P_syn(I_nu)=P_syn(I_nu)+p_nu*transfer
        Seed_syn(I_nu)=Seed_syn(I_nu)+p_nu*transfer/(r2*V_seed(I_nu)*hfac)
    end do
end subroutine add_cyclotron

! 同步辐射F(x)核：F(x) = 1.81 e^(-x)/√(x^(-2/3)+factor)，x=ν/ν_c。
real(8) function syn_fx(gam,V_cal,DB,factor)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: gam,V_cal,DB,factor
real(8) :: Vc,x,ratio_v

    Vc=(4.2d6)*gam*gam*DB
    x=V_cal/Vc
    ratio_v=Vc/V_cal
    syn_fx=syn_kernel(x,ratio_v**(2d0/3d0),factor)
end function syn_fx

! 线性插值：x∈[x0,x1]→y，等距时取平均。
real(8) function linear_interp(x0,x1,y0,y1,x)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,y0,y1,x

    if (x1 <= x0) then
        linear_interp=0.5d0*(y0+y1)
    else
        linear_interp=y0+(y1-y0)*(x-x0)/(x1-x0)
    end if
end function linear_interp

! 同步辐射被积函数：dN/dx * F(x) * γ，用于发射功率积分。
real(8) function syn_integrand(x,x0,x1,dN0,dN1,V_cal,DB,factor)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x,x0,x1,dN0,dN1,V_cal,DB,factor
real(8) :: gam,dN

    gam=dexp(x)
    dN=max(0d0,linear_interp(x0,x1,dN0,dN1,x))
    syn_integrand=dN*syn_fx(gam,V_cal,DB,factor)*gam
end function syn_integrand

! 幂律插值：两端为正时做log-log线性插值，否则退化为线性插值。
real(8) function pl_interp(v0,v1,y0,y1,v)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: v0,v1,y0,y1,v
real(8) :: slope

    if (v <= v0) then
        pl_interp=y0
        return
    end if
    if (v >= v1) then
        pl_interp=y1
        return
    end if

    if (v1 <= v0) then
        pl_interp=0.5d0*(y0+y1)
    else if (y0 > 0d0 .and. y1 > 0d0) then
        slope=dlog(y1/y0)/dlog(v1/v0)
        pl_interp=y0*(v/v0)**slope
    else
        pl_interp=y0+(y1-y0)*(v-v0)/(v1-v0)
    end if
end function pl_interp

! 对数空间2点Gauss-Legendre节点和权重：(v_g1, v_g2)为频率节点，(w_g1, w_g2)含dv权重。
subroutine log_gauss2(v0,v1,v_g1,v_g2,w_g1,w_g2)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: v0,v1
real(8), intent(out) :: v_g1,v_g2,w_g1,w_g2
real(8) :: x0,x1,xm,dx,w2

    x0=dlog(v0)
    x1=dlog(v1)
    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    w2=1d0/dsqrt(3d0)
    v_g1=dexp(xm-dx*w2)
    v_g2=dexp(xm+dx*w2)
    w_g1=dx*v_g1
    w_g2=dx*v_g2
end subroutine log_gauss2

! 对幂律插值函数在[v0,v1]上做2点Gauss-Legendre积分。
real(8) function pl_integral(v0,v1,y0,y1)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: v0,v1,y0,y1
real(8) :: vg1,vg2,wg1,wg2

    if (v1 <= v0) then
        pl_integral=0d0
        return
    end if

    call log_gauss2(v0,v1,vg1,vg2,wg1,wg2)
    pl_integral= &
        wg1*pl_interp(v0,v1,y0,y1,vg1)+ &
        wg2*pl_interp(v0,v1,y0,y1,vg2)
end function pl_integral

! SSA光深核函数：γ²F(ν/νc)，电子谱导数由有限体积端点差单独给出。
real(8) function tau_kernel(x,V_cal,DB,factor)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x,V_cal,DB,factor
real(8) :: gam

    gam=dexp(x)
    tau_kernel=gam*gam*syn_fx(gam,V_cal,DB,factor)
end function tau_kernel

! 单网格单元同步发射功率的2点+3点Gauss积分（用于自适应误差估计）。
subroutine syn_gauss(x0,x1,dN0,dN1,V_cal,DB,factor,p2,p3)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,dN0,dN1,V_cal,DB,factor
real(8), intent(out) :: p2,p3
real(8) :: xm, dx, w2, w3a, x2a, x2b, x3a, x3b

    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    w2=1d0/dsqrt(3d0)
    x2a=xm-dx*w2
    x2b=xm+dx*w2
    p2=dx*(syn_integrand(x2a,x0,x1,dN0,dN1,V_cal,DB,factor)+ &
           syn_integrand(x2b,x0,x1,dN0,dN1,V_cal,DB,factor))

    w3a=dsqrt(3d0/5d0)
    x3a=xm-dx*w3a
    x3b=xm+dx*w3a
    p3=dx*((5d0/9d0)*syn_integrand(x3a,x0,x1,dN0,dN1,V_cal,DB,factor)+ &
           (8d0/9d0)*syn_integrand(xm ,x0,x1,dN0,dN1,V_cal,DB,factor)+ &
           (5d0/9d0)*syn_integrand(x3b,x0,x1,dN0,dN1,V_cal,DB,factor))
end subroutine syn_gauss

! 单网格单元SSA光深的2点+3点Gauss积分：有限体积端点差给出 -dq，Gauss积分给核函数平均。
subroutine tau_gauss(x0,x1,dN10,dN11,V_cal,DB,factor,t2,t3)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,dN10,dN11,V_cal,DB,factor
real(8), intent(out) :: t2,t3
real(8) :: xm, dx, width, w2, w3a, q_drop, x2a, x2b, x3a, x3b

    if (x1 <= x0) then
        t2=0d0
        t3=0d0
        return
    end if

    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    width=x1-x0
    q_drop=dN10-dN11
    w2=1d0/dsqrt(3d0)
    x2a=xm-dx*w2
    x2b=xm+dx*w2
    t2=q_drop/width*dx*(tau_kernel(x2a,V_cal,DB,factor)+ &
                        tau_kernel(x2b,V_cal,DB,factor))

    w3a=dsqrt(3d0/5d0)
    x3a=xm-dx*w3a
    x3b=xm+dx*w3a
    t3=q_drop/width*dx*((5d0/9d0)*tau_kernel(x3a,V_cal,DB,factor)+ &
                        (8d0/9d0)*tau_kernel(xm ,V_cal,DB,factor)+ &
                        (5d0/9d0)*tau_kernel(x3b,V_cal,DB,factor))
end subroutine tau_gauss

! 自适应同步辐射单元积分：发射率用线性重构，SSA光深用保守端点差重构。
subroutine syn_adapt(x0,x1,dN0,dN1,dN10,dN11, &
                                      V_cal,DB,factor,rel_tol,p_int,tau_int)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,dN0,dN1,dN10,dN11,V_cal,DB,factor,rel_tol
real(8), intent(out) :: p_int,tau_int
real(8) :: p2, p3, t2, t3, xm, dNm, dN1m, err_p, err_t, ref_p, ref_t, p3_l, p3_r, t3_l, t3_r

    call syn_gauss(x0,x1,dN0,dN1,V_cal,DB,factor,p2,p3)
    call tau_gauss(x0,x1,dN10,dN11,V_cal,DB,factor,t2,t3)
    ref_p=max(abs(p3),1d-30)
    ref_t=max(abs(t3),1d-30)
    err_p=abs(p3-p2)/ref_p
    err_t=abs(t3-t2)/ref_t

    if (max(err_p,err_t) <= rel_tol) then
        p_int=p3
        tau_int=t3
    else
        xm=0.5d0*(x0+x1)
        dNm=0.5d0*(dN0+dN1)
        dN1m=0.5d0*(dN10+dN11)
        call syn_gauss(x0,xm,dN0,dNm,V_cal,DB,factor,p2,p3_l)
        call syn_gauss(xm,x1,dNm,dN1,V_cal,DB,factor,p2,p3_r)
        call tau_gauss(x0,xm,dN10,dN1m,V_cal,DB,factor,t2,t3_l)
        call tau_gauss(xm,x1,dN1m,dN11,V_cal,DB,factor,t2,t3_r)
        p_int=p3_l+p3_r
        tau_int=t3_l+t3_r
    end if
end subroutine syn_adapt

! 计算同步辐射转移函数：Transfer = P_absorbed / P_emit，即(1-e^(-τ))/τ。
subroutine syn_transfer(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,Transfer_syn)
!$ use omp_lib
implicit none
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in), dimension(Num_gam_e) :: gam_e,dN_gam_e
real(8), intent(in), dimension(Num_nu) :: V_seed
real(8), intent(in) :: R_loc,DB
real(8), intent(out), dimension(Num_nu) :: Transfer_syn
real(8), dimension(1) :: DB_chi, Weight_chi
real(8), dimension(Num_gam_e,1) :: DNe_chi
real(8), dimension(Num_nu,1) :: P_emit,P_syn,Seed_syn,Tau_syn
integer :: I_nu

    DB_chi(1)=DB
    Weight_chi(1)=1d0
    DNe_chi(:,1)=dN_gam_e
    call syn_seed_chi(R_loc,Num_gam_e,Num_nu,1,gam_e,DNe_chi,V_seed,DB_chi,Weight_chi,1.046d4, &
                                           P_emit,P_syn,Seed_syn,Tau_syn)
    !$OMP PARALLEL DO SCHEDULE(STATIC) num_threads(n_threads) private(I_nu)
    do I_nu=1,Num_nu
        if (P_emit(I_nu,1) > 0d0 .and. P_syn(I_nu,1) >= 0d0) then
            Transfer_syn(I_nu)=P_syn(I_nu,1)/P_emit(I_nu,1)
        else
            Transfer_syn(I_nu)=1d0
        end if
    end do
    !$OMP END PARALLEL DO
end subroutine syn_transfer

! 计算SSA频率ν_a：二分搜索+对数割线法求解τ(ν_a)=1。
subroutine nua_solve(R_loc,DB,Num_gam_e,gam_e,dN_gam_e, &
                       V_a)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e
real(8), intent(in), dimension(Num_gam_e) :: gam_e,dN_gam_e
real(8), intent(in) :: R_loc,DB
real(8), intent(out) :: V_a
real(8), parameter :: rel_tol=5d-4

real(8),allocatable,dimension (:) :: dN1,x_gam
allocate (dN1(Num_gam_e),x_gam(Num_gam_e))

    factor=(3.62d0/pi)**2
    r2=R_loc*R_loc
    dN1=dN_gam_e/(gam_e*gam_e)
    x_gam=dlog(gam_e)
    V_a_floor=1d1**4d0
    V_a_min=1d1**(-20d0)
    V_a_cap=1d0
    do I_gam_e=1,Num_gam_e-1
       Vc=4.2d6*((gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0)*DB
       if (Vc > V_a_cap) V_a_cap=Vc
    end do
    V_a_cap=max(1d1**14d0, min(1d1**30d0, 1d1*V_a_cap))

    V_low=V_a_floor
    call evaluate_tau(V_low,Tau_low)
    if (Tau_low <= 1d0) then
        do I_nu=1,26
           if (V_low <= V_a_min) exit
           V_high=V_low
           Tau_high=Tau_low
           V_low=max(V_a_min,V_low/1d1)
           call evaluate_tau(V_low,Tau_low)
           if (Tau_low > 1d0) exit
        end do

        if (Tau_low > 1d0) then
           call refine_tau(V_low,Tau_low,V_high,Tau_high,V_a)
        else
           V_a=V_low
        end if
    else
        V_high=V_low
        Tau_high=Tau_low
       do I_nu=1,26
          V_low=V_high
          Tau_low=Tau_high
          if (V_high >= V_a_cap) exit
          V_high=min(V_a_cap, 1d1*V_high)
          call evaluate_tau(V_high,Tau_high)
          if (Tau_high <= 1d0) exit
       end do

       if (Tau_high > 1d0) then
          error stop 'nua_solve: tau remains above unity at the adaptive upper bound.'
       else
          if (Tau_low <= 1d0 .or. Tau_low == Tau_high) then
             V_a=V_high
          else
             call refine_tau(V_low,Tau_low,V_high,Tau_high,V_a)
          end if
       end if
    end if

    deallocate (dN1,x_gam)

return

contains

! 计算给定频率ν处的SSA光深（遍历电子能格自适应积分）。
subroutine evaluate_tau(V_cal,Tau)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: V_cal
real(8), intent(out) :: Tau
real(8) :: tau_cell,discard

    Tau=0d0
    do I_gam_e=1,Num_gam_e-1
       call syn_adapt(x_gam(I_gam_e),x_gam(I_gam_e+1), &
                                       dN_gam_e(I_gam_e),dN_gam_e(I_gam_e+1), &
                                       dN1(I_gam_e),dN1(I_gam_e+1), &
                                       V_cal,DB,factor,rel_tol,discard,tau_cell)
       Tau=Tau+tau_cell
    end do
    Tau=1.046d4*Tau*DB/(4d0*pi*r2*V_cal*V_cal)
end subroutine evaluate_tau

! 对数空间割线法迭代精化ν_a：在[V_low,V_high]内求解τ(ν)=1，最多10次迭代。
subroutine refine_tau(V_low,Tau_low,V_high,Tau_high,V_root)
implicit REAL(8)(A-H,O-Z)
real(8), intent(inout) :: V_low,Tau_low,V_high,Tau_high
real(8), intent(out) :: V_root
real(8) :: log_v_low, log_v_high, log_v_mid, log_tau_low, log_tau_high, V_mid, Tau_mid
integer :: I_iter

    log_v_low=dlog(V_low)
    log_v_high=dlog(V_high)

    do I_iter=1,10
        log_tau_low=dlog(Tau_low)
        log_tau_high=dlog(Tau_high)
        if (log_tau_low == log_tau_high) then
            log_v_mid=0.5d0*(log_v_low+log_v_high)
        else
            log_v_mid=log_v_low-log_tau_low*(log_v_high-log_v_low)/(log_tau_high-log_tau_low)
        end if
        log_v_mid=max(log_v_low,min(log_v_high,log_v_mid))
        V_mid=dexp(log_v_mid)
        call evaluate_tau(V_mid,Tau_mid)

        if (Tau_mid > 1d0) then
            V_low=V_mid
            Tau_low=Tau_mid
            log_v_low=log_v_mid
        else
            V_high=V_mid
            Tau_high=Tau_mid
            log_v_high=log_v_mid
        end if

        if (abs(log_v_high-log_v_low) <= 5d-3) exit
    end do

    log_tau_low=dlog(Tau_low)
    log_tau_high=dlog(Tau_high)
    if (log_tau_low == log_tau_high) then
        V_root=dexp(0.5d0*(log_v_low+log_v_high))
    else
        V_root=dexp(log_v_low-log_tau_low*(log_v_high-log_v_low)/(log_tau_high-log_tau_low))
    end if
end subroutine refine_tau
end subroutine nua_solve

! 2D逐列ν_a：对每个χ列累加光深后求解τ_cumul(ν_a)=1，输出各列ν_a(χ)。
subroutine nua_path(Num_nu,Num_chi,V_seed,Tau_chi,V_a_chi)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu,Num_chi
integer :: I_chi
real(8), intent(in), dimension(Num_nu) :: V_seed
real(8), intent(in), dimension(Num_nu,Num_chi) :: Tau_chi
real(8), intent(out), dimension(Num_chi) :: V_a_chi
real(8), dimension(Num_nu) :: Tau_path

    Tau_path=0d0
    do I_chi=1,Num_chi
        Tau_path=Tau_path+Tau_chi(:,I_chi)
        call nua_fromtau(Num_nu,V_seed,Tau_path,V_a_chi(I_chi))
    end do
end subroutine nua_path

! 从光深网格求解ν_a：在τ穿越1的位置对数插值。
subroutine nua_fromtau(Num_nu,V_seed,Tau_grid,V_a)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu
integer :: I_nu,I_cross
real(8), intent(in), dimension(Num_nu) :: V_seed,Tau_grid
real(8), intent(out) :: V_a
real(8) :: Tau_left,Tau_right

    if (Tau_grid(1) <= 1d0) then
        V_a=V_seed(1)
        return
    end if

    if (Tau_grid(Num_nu) >= 1d0) then
        V_a=V_seed(Num_nu)
        return
    end if

    I_cross=0
    do I_nu=2,Num_nu
        if (Tau_grid(I_nu) <= 1d0) then
            I_cross=I_nu
            exit
        end if
    end do

    if (I_cross == 0) then
        if (Tau_grid(Num_nu-1) <= 0d0 .or. Tau_grid(Num_nu) <= 0d0) &
            error stop 'nua_fromtau requires positive tau values for log interpolation'
        call tau_root(V_seed(Num_nu-1),Tau_grid(Num_nu-1), &
                                      V_seed(Num_nu),Tau_grid(Num_nu),V_a)
        if (V_a < V_seed(Num_nu)) V_a=V_seed(Num_nu)
        return
    end if

    Tau_left=Tau_grid(I_cross-1)
    Tau_right=Tau_grid(I_cross)
    if (Tau_left <= 0d0 .or. Tau_right <= 0d0) &
        error stop 'nua_fromtau requires positive tau values for log interpolation'
    call tau_root(V_seed(I_cross-1),Tau_left,V_seed(I_cross),Tau_right,V_a)
end subroutine nua_fromtau

! 对数-对数空间线性插值求解τ(ν)=1的根ν_a。
subroutine tau_root(V_left,Tau_left,V_right,Tau_right,V_root)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: V_left,Tau_left,V_right,Tau_right
real(8), intent(out) :: V_root
real(8) :: log_v_left,log_v_right,log_tau_left,log_tau_right,log_v_root

    log_v_left=dlog(V_left)
    log_v_right=dlog(V_right)
    log_tau_left=dlog(Tau_left)
    log_tau_right=dlog(Tau_right)

    if (log_tau_left == log_tau_right) then
        log_v_root=0.5d0*(log_v_left+log_v_right)
    else
        log_v_root=log_v_left-log_tau_left*(log_v_right-log_v_left)/(log_tau_right-log_tau_left)
    end if
    V_root=dexp(log_v_root)
end subroutine tau_root

end module electron_radiation_kernel
