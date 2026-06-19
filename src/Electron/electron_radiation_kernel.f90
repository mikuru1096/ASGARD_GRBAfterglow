module electron_radiation_kernel
  use constants
  use radiation_common, only: radiation_syn_seed_core, radiation_transfer_factor
  use synchrotron_polarization_kernel, only: synchrotron_polarized_components
  private

    public :: first_greater_monotonic, first_greater_monotonic_window
    public :: besselk, get_syn, get_syn_state, get_syn_cyclotron_state, get_syn_selected, get_syn_selected_state
    public :: get_syn_transfer, get_syn_polarization_selected, get_nu_a
    public :: get_nu_a_2d_path, get_nu_a_2d_cell_path, reduce_syn_shell_from_chi
    public :: build_reduced_log_grid, project_syn_state_logbands
    public :: get_syn_adaptive, get_syn_adaptive_state, get_nu_a_from_tau_grid
    public :: electron_powerlaw_interp, electron_log_gauss2_interval, electron_integrate_powerlaw_segment, electron_ssa_segment

contains

! 在单调递增数组arr的[left,right]区间内二分查找第一个 > target 的索引。
subroutine first_greater_monotonic_from(arr,left,right,target,idx)
implicit none
integer, intent(in) :: left,right
real(8), intent(in) :: arr(right),target
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
end subroutine first_greater_monotonic_from

! 在全数组arr(1:n)中二分查找第一个 > target 的索引，含边界快速返回。
subroutine first_greater_monotonic(arr,n,target,idx)
implicit none
integer, intent(in) :: n
real(8), intent(in) :: arr(n),target
integer, intent(out) :: idx

    if (arr(1) > target) then
        idx=1
        return
    end if
    if (arr(n) <= target) then
        idx=n+1
        return
    end if

    call first_greater_monotonic_from(arr,1,n,target,idx)
end subroutine first_greater_monotonic

! 从start_idx开始在arr(1:n)中查找第一个 > target 的索引（带游标加速）。
subroutine first_greater_monotonic_window(arr,n,start_idx,target,idx)
implicit none
integer, intent(in) :: n,start_idx
real(8), intent(in) :: arr(n),target
integer, intent(out) :: idx

    if (arr(start_idx) > target) then
        idx=start_idx
        return
    end if
    if (arr(n) <= target) then
        idx=n+1
        return
    end if

    call first_greater_monotonic_from(arr,start_idx,n,target,idx)
end subroutine first_greater_monotonic_window

!****************************************************************************************
!**************************** Besselk function interpolation ****************************
!****************************************************************************************
function besselk(var)
    REAL(8) :: val,var,besselk
    integer :: i
    REAL(8),dimension(100),parameter :: minus_theta= &
                 (/1.0000000e-04,1.1233240e-04,1.2618569e-04,1.4174742e-04,1.5922828e-04,1.7886495e-04,2.0092330e-04,& 
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
                 8.9021509e+00,1.0000000e+01/)
    REAL(8),dimension(100),parameter :: besselk2= &
                 (/2.0000000e+08,1.5849658e+08,1.2560583e+08,9.9540471e+07,7.8884121e+07,6.2514316e+07,4.9541527e+07,&
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
                 6.9778918e-05,2.1509817e-05/)

   if (var <= zero) error stop 'besselk requires var > 0'
   if (var < minus_theta(1)) then
       besselk=two/var**2
       return
   end if
   if (var > minus_theta(100)) then
       besselk=dsqrt(pi/(two*var))*dexp(-var)*(one+15d0/(8d0*var)+105d0/(128d0*var**2))
       return
   end if

   besselk=zero
   
   do i=1,99
       if (var >= minus_theta(i) .and. var <= minus_theta(i+1)) then
           val = (var-minus_theta(i))/(minus_theta(i+1)-minus_theta(i))
           besselk = besselk2(i)+val*(besselk2(i+1)-besselk2(i))
           exit
       end if
   end do

   return
end function besselk


!****************************************************************************************
!************************** get syn power and number density ****************************
!****************************************************************************************
subroutine get_syn(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                       P_syn,Seed_syn)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) ::P_syn(Num_nu),Seed_syn(Num_nu)
real(8) :: P_emit(Num_nu),Tau_syn(Num_nu)

    call get_syn_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,P_emit,P_syn,Seed_syn,Tau_syn)
end subroutine get_syn

! GitHub基线的复合Simpson电子能量积分；case(2)不能退化成区间中点核，否则高频尾会断崖式归零。
subroutine get_syn_simpson(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                           P_syn,Seed_syn)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: P_syn(Num_nu),Seed_syn(Num_nu)
real(8) :: P_emit(Num_nu),Tau_syn(Num_nu)

    call get_syn_simpson_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                               P_emit,P_syn,Seed_syn,Tau_syn)
end subroutine get_syn_simpson

subroutine get_syn_simpson_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                 P_emit,P_syn,Seed_syn,Tau_syn)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: P_emit(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu),Tau_syn(Num_nu)
real(8) :: dN1(Num_gam_e),ddN(Num_gam_e-1),simpson_weight(Num_gam_e),emit_weight(Num_gam_e)
real(8) :: gam_e_mean2_arr(Num_gam_e-1),tau_weight(Num_gam_e-1),V_seed_powm23(Num_nu)
real(8) :: Vc_emit_inv(Num_gam_e),Vc_emit_pow23(Num_gam_e),Vc_tau_inv(Num_gam_e-1),Vc_tau_pow23(Num_gam_e-1)
integer :: I_nu

    factor=(3.62d0/pi)**2
    Temp_syn=dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    Rariv2=R_loc*R_loc
    dN1=dN_gam_e/(gam_e*gam_e)
    ddN=dN1(1:Num_gam_e-1)-dN1(2:Num_gam_e)
    h=dlog(gam_e(2))-dlog(gam_e(1))
    do I_gam_e=1,Num_gam_e
        if (I_gam_e == 1 .or. I_gam_e == Num_gam_e) then
            simpson_weight(I_gam_e)=one
        else if (mod(I_gam_e,2) == 0) then
            simpson_weight(I_gam_e)=4d0
        else
            simpson_weight(I_gam_e)=2d0
        end if
        Vc=4.2d6*gam_e(I_gam_e)*gam_e(I_gam_e)*DB
        Vc_emit_inv(I_gam_e)=one/Vc
        Vc_emit_pow23(I_gam_e)=Vc**(2d0/3d0)
        emit_weight(I_gam_e)=simpson_weight(I_gam_e)*dN_gam_e(I_gam_e)*gam_e(I_gam_e)
    end do
    do I_gam_e=1,Num_gam_e-1
        gam_e_mean2_arr(I_gam_e)=(gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0
        Vc=4.2d6*gam_e_mean2_arr(I_gam_e)*DB
        Vc_tau_inv(I_gam_e)=one/Vc
        Vc_tau_pow23(I_gam_e)=Vc**(2d0/3d0)
        tau_weight(I_gam_e)=gam_e_mean2_arr(I_gam_e)*ddN(I_gam_e)
    end do
    do I_nu=1,Num_nu
        V_seed_powm23(I_nu)=V_seed(I_nu)**(-2d0/3d0)
    end do

    !$OMP PARALLEL num_threads(n_threads), private(I_nu)
    !$OMP DO SCHEDULE(STATIC)
    do I_nu=1,Num_nu
        call accumulate_simpson_syn_point(I_nu)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    temp_para=4d0*pi*Para_c*Para_h
    Seed_syn=Seed_syn/temp_para

    return

contains

real(8) function simpson_emission_integral(V_cal,V_powm23)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: V_cal,V_powm23
real(8) :: simpson_sum,x,Fx
integer :: I_gam_e

    simpson_sum=zero
    do I_gam_e=1,Num_gam_e
        x=V_cal*Vc_emit_inv(I_gam_e)
        Fx=1.81d0*dexp(-x)/dsqrt(Vc_emit_pow23(I_gam_e)*V_powm23+factor)
        simpson_sum=simpson_sum+emit_weight(I_gam_e)*Fx
    end do
    simpson_emission_integral=h*simpson_sum/3d0
end function simpson_emission_integral

real(8) function simpson_ssa_tau_integral(V_cal,V_powm23)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: V_cal,V_powm23
real(8) :: x,Fx
integer :: I_gam_e

    simpson_ssa_tau_integral=zero
    do I_gam_e=1,Num_gam_e-1
        x=V_cal*Vc_tau_inv(I_gam_e)
        Fx=1.81d0*dexp(-x)/dsqrt(Vc_tau_pow23(I_gam_e)*V_powm23+factor)
        simpson_ssa_tau_integral=simpson_ssa_tau_integral+tau_weight(I_gam_e)*Fx
    end do
end function simpson_ssa_tau_integral

subroutine accumulate_simpson_syn_point(I_nu)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: I_nu
real(8) :: V_cal,dInteg,Tau,P_v

    V_cal=V_seed(I_nu)
    dInteg=simpson_emission_integral(V_cal,V_seed_powm23(I_nu))
    Tau=simpson_ssa_tau_integral(V_cal,V_seed_powm23(I_nu))
    P_v=Temp_syn*DB*dInteg
    Tau=1.046d4*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal)
    if ((Tau-1d-4) < 1d-5) Tau=1d-4
    P_emit(I_nu)=P_v
    Tau_syn(I_nu)=Tau
    P_syn(I_nu)=P_v*(one-dexp(-Tau))/Tau
    Seed_syn(I_nu)=P_syn(I_nu)/(Rariv2*V_cal)
end subroutine accumulate_simpson_syn_point
end subroutine get_syn_simpson_state

! 构建约化对数频率网格：在85%高频处加密采样，用于冷却计算的轻量级频率表。
subroutine build_reduced_log_grid(Num_nu_in,V_in,Num_nu_out,V_out)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu_in,Num_nu_out
real(8), intent(in) :: V_in(Num_nu_in)
real(8), intent(out) :: V_out(Num_nu_out)
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
end subroutine build_reduced_log_grid

! 将同步辐射状态(P,Seed,Tau)从原频率网格幂律插值投影到约化频率网格。
subroutine project_syn_state_logbands(Num_nu_in,V_in,P_in,Seed_in,Tau_in,Num_nu_out,V_out,P_out,Seed_out,Tau_out)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu_in,Num_nu_out
integer :: I_out, idx_lo, idx_hi
real(8), intent(in) :: V_in(Num_nu_in),P_in(Num_nu_in),Seed_in(Num_nu_in),Tau_in(Num_nu_in),V_out(Num_nu_out)
real(8), intent(out) :: P_out(Num_nu_out),Seed_out(Num_nu_out),Tau_out(Num_nu_out)
real(8) :: V_tar

    idx_lo=1
    do I_out=1,Num_nu_out
        V_tar=min(max(V_out(I_out),V_in(1)),V_in(Num_nu_in))
        if (V_tar <= V_in(1)) then
            P_out(I_out)=max(P_in(1),zero)
            Seed_out(I_out)=max(Seed_in(1),zero)
            Tau_out(I_out)=max(Tau_in(1),1d-4)
            cycle
        end if
        if (V_tar >= V_in(Num_nu_in)) then
            P_out(I_out)=max(P_in(Num_nu_in),zero)
            Seed_out(I_out)=max(Seed_in(Num_nu_in),zero)
            Tau_out(I_out)=max(Tau_in(Num_nu_in),1d-4)
            cycle
        end if

        call first_greater_monotonic_window(V_in,Num_nu_in,idx_lo,V_tar,idx_hi)
        idx_lo=max(1,min(idx_hi-1,Num_nu_in-1))
        P_out(I_out)=max(electron_powerlaw_interp(V_in(idx_lo),V_in(idx_hi),P_in(idx_lo),P_in(idx_hi),V_tar),zero)
        Seed_out(I_out)=max(electron_powerlaw_interp(V_in(idx_lo),V_in(idx_hi),Seed_in(idx_lo),Seed_in(idx_hi),V_tar),zero)
        Tau_out(I_out)=max(electron_powerlaw_interp(V_in(idx_lo),V_in(idx_hi),Tau_in(idx_lo),Tau_in(idx_hi),V_tar),1d-4)
    end do
end subroutine project_syn_state_logbands

! 调用radiation_syn_seed_core计算同步辐射发射功率、SSA光深、转移后谱和种子光子场。
subroutine get_syn_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                         P_emit,P_syn,Seed_syn,Tau_syn)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) ::P_emit(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu),Tau_syn(Num_nu)
    call radiation_syn_seed_core(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,1.046d4, &
                                 P_emit,P_syn,Seed_syn,Tau_syn)
end subroutine get_syn_state

! 同步+非相对论回旋发射核：γ<2 的电子使用基频回旋发射，γ>=2 仍走标准同步核。
subroutine get_syn_cyclotron_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                   P_emit,P_syn,Seed_syn,Tau_syn)
implicit none
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: P_emit(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu),Tau_syn(Num_nu)
real(8) :: dN_syn(Num_gam_e)
integer :: i

    dN_syn=zero
    do i=1,Num_gam_e
        if (gam_e(i) >= two) dN_syn(i)=dN_gam_e(i)
    end do
    call get_syn_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_syn,V_seed,P_emit,P_syn,Seed_syn,Tau_syn)
    call add_cyclotron_fundamental(R_loc,DB,Num_gam_e,Num_nu,gam_e,dN_gam_e,V_seed,P_emit,P_syn,Seed_syn,Tau_syn)
end subroutine get_syn_cyclotron_state

subroutine get_syn_cyclotron(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,P_syn,Seed_syn)
implicit none
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: P_syn(Num_nu),Seed_syn(Num_nu)
real(8) :: P_emit(Num_nu),Tau_syn(Num_nu)

    call get_syn_cyclotron_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                 P_emit,P_syn,Seed_syn,Tau_syn)
end subroutine get_syn_cyclotron

subroutine add_cyclotron_fundamental(R_loc,DB,Num_gam_e,Num_nu,gam_e,dN_gam_e,V_seed,P_emit,P_syn,Seed_syn,Tau_syn)
implicit none
integer, intent(in) :: Num_gam_e,Num_nu
integer :: i,j
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu),Tau_syn(Num_nu)
real(8), intent(inout) :: P_emit(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu)
real(8) :: nu_edge(Num_nu+1),nu_b,nu0,gamma_mid,beta2,n_e_seg,p_total,p_nu,transfer
real(8) :: r2,temp_para

    if (DB <= zero) return
    call build_log_frequency_edges(Num_nu,V_seed,nu_edge)
    nu_b=Para_e*DB/(two*pi*Para_m_e*Para_c)
    r2=R_loc*R_loc
    temp_para=4d0*pi*Para_c*Para_h
    do i=1,Num_gam_e-1
        gamma_mid=0.5d0*(gam_e(i)+gam_e(i+1))
        if (gamma_mid >= two) cycle
        n_e_seg=0.5d0*(dN_gam_e(i)+dN_gam_e(i+1))*(gam_e(i+1)-gam_e(i))
        if (n_e_seg <= zero) cycle
        beta2=one-one/(gamma_mid*gamma_mid)
        nu0=nu_b/gamma_mid
        if (nu0 < nu_edge(1) .or. nu0 >= nu_edge(Num_nu+1)) cycle
        call first_greater_monotonic(nu_edge,Num_nu+1,nu0,j)
        j=max(1,min(Num_nu,j-1))
        p_total=n_e_seg*(4d0/3d0)*Para_SigmaT*Para_c*(DB*DB/(8d0*pi))*gamma_mid*gamma_mid*beta2
        p_nu=p_total/(nu_edge(j+1)-nu_edge(j))
        call radiation_transfer_factor(Tau_syn(j),transfer)
        P_emit(j)=P_emit(j)+p_nu
        P_syn(j)=P_syn(j)+p_nu*transfer
        Seed_syn(j)=Seed_syn(j)+p_nu*transfer/(r2*V_seed(j)*temp_para)
    end do
end subroutine add_cyclotron_fundamental

subroutine build_log_frequency_edges(Num_nu,V_seed,nu_edge)
implicit none
integer, intent(in) :: Num_nu
integer :: i
real(8), intent(in) :: V_seed(Num_nu)
real(8), intent(out) :: nu_edge(Num_nu+1)

    nu_edge(1)=V_seed(1)*dsqrt(V_seed(1)/V_seed(2))
    do i=2,Num_nu
        nu_edge(i)=dsqrt(V_seed(i-1)*V_seed(i))
    end do
    nu_edge(Num_nu+1)=V_seed(Num_nu)*dsqrt(V_seed(Num_nu)/V_seed(Num_nu-1))
end subroutine build_log_frequency_edges

! 同步辐射F(x)核：F(x) = 1.81 e^(-x)/√(x^(-2/3)+factor)，x=ν/ν_c。
real(8) function electron_syn_fx(gam,V_cal,DB,factor)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: gam,V_cal,DB,factor
real(8) :: Vc,x,ratio_v

    Vc=(4.2d6)*gam*gam*DB
    x=V_cal/Vc
    ratio_v=Vc/V_cal
    electron_syn_fx=1.81d0*dexp(-x)/dsqrt(ratio_v**(2d0/3d0)+factor)
end function electron_syn_fx

! 线性插值：x∈[x0,x1]→y，等距时取平均。
real(8) function electron_linear_interp(x0,x1,y0,y1,x)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,y0,y1,x

    if (x1 <= x0) then
        electron_linear_interp=0.5d0*(y0+y1)
    else
        electron_linear_interp=y0+(y1-y0)*(x-x0)/(x1-x0)
    end if
end function electron_linear_interp

! 同步辐射被积函数：dN/dx * F(x) * γ，用于发射功率积分。
real(8) function electron_syn_integrand_x(x,x0,x1,dN0,dN1,V_cal,DB,factor)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x,x0,x1,dN0,dN1,V_cal,DB,factor
real(8) :: gam,dN

    gam=dexp(x)
    dN=max(zero,electron_linear_interp(x0,x1,dN0,dN1,x))
    electron_syn_integrand_x=dN*electron_syn_fx(gam,V_cal,DB,factor)*gam
end function electron_syn_integrand_x

! 幂律插值：两端为正时做log-log线性插值，否则退化为线性插值。
real(8) function electron_powerlaw_interp(v0,v1,y0,y1,v)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: v0,v1,y0,y1,v
real(8) :: slope

    if (v <= v0) then
        electron_powerlaw_interp=y0
        return
    end if
    if (v >= v1) then
        electron_powerlaw_interp=y1
        return
    end if

    if (v1 <= v0) then
        electron_powerlaw_interp=0.5d0*(y0+y1)
    else if (y0 > zero .and. y1 > zero) then
        slope=dlog(y1/y0)/dlog(v1/v0)
        electron_powerlaw_interp=y0*(v/v0)**slope
    else
        electron_powerlaw_interp=y0+(y1-y0)*(v-v0)/(v1-v0)
    end if
end function electron_powerlaw_interp

! 对数空间2点Gauss-Legendre节点和权重：(v_g1, v_g2)为频率节点，(w_g1, w_g2)含dv权重。
subroutine electron_log_gauss2_interval(v0,v1,v_g1,v_g2,w_g1,w_g2)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: v0,v1
real(8), intent(out) :: v_g1,v_g2,w_g1,w_g2
real(8) :: x0,x1,xm,dx,w2

    x0=dlog(v0)
    x1=dlog(v1)
    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    w2=one/dsqrt(3d0)
    v_g1=dexp(xm-dx*w2)
    v_g2=dexp(xm+dx*w2)
    w_g1=dx*v_g1
    w_g2=dx*v_g2
end subroutine electron_log_gauss2_interval

! 对幂律插值函数在[v0,v1]上做2点Gauss-Legendre积分。
real(8) function electron_integrate_powerlaw_segment(v0,v1,y0,y1)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: v0,v1,y0,y1
real(8) :: vg1,vg2,wg1,wg2

    if (v1 <= v0) then
        electron_integrate_powerlaw_segment=zero
        return
    end if

    call electron_log_gauss2_interval(v0,v1,vg1,vg2,wg1,wg2)
    electron_integrate_powerlaw_segment= &
        wg1*electron_powerlaw_interp(v0,v1,y0,y1,vg1)+ &
        wg2*electron_powerlaw_interp(v0,v1,y0,y1,vg2)
end function electron_integrate_powerlaw_segment

! SSA冷却率段积分：mode=1低频Σ∝ν^(-5/3)，mode=2高频Σ∝(ν_c/ν)e^(-ν/ν_uplim)。
real(8) function electron_ssa_segment(v0,v1,seed0,seed1,sigma_prefactor,mode,Cyclotron_nu,V_uplim)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: mode
real(8), intent(in) :: v0,v1,seed0,seed1,sigma_prefactor,Cyclotron_nu,V_uplim
real(8) :: vg1,vg2,wg1,wg2,seed_loc,sigma_loc

    if (v1 <= v0) then
        electron_ssa_segment=zero
        return
    end if

    call electron_log_gauss2_interval(v0,v1,vg1,vg2,wg1,wg2)
    electron_ssa_segment=zero

    seed_loc=electron_powerlaw_interp(v0,v1,seed0,seed1,vg1)
    if (mode == 1) then
        sigma_loc=sigma_prefactor*vg1**(-5d0/3d0)
    else
        sigma_loc=sigma_prefactor*(Cyclotron_nu/vg1)*dexp(-vg1/V_uplim)
    end if
    electron_ssa_segment=electron_ssa_segment+wg1*sigma_loc*seed_loc*para_h*vg1*para_c

    seed_loc=electron_powerlaw_interp(v0,v1,seed0,seed1,vg2)
    if (mode == 1) then
        sigma_loc=sigma_prefactor*vg2**(-5d0/3d0)
    else
        sigma_loc=sigma_prefactor*(Cyclotron_nu/vg2)*dexp(-vg2/V_uplim)
    end if
    electron_ssa_segment=electron_ssa_segment+wg2*sigma_loc*seed_loc*para_h*vg2*para_c
end function electron_ssa_segment

! SSA光深核函数：γ²F(ν/νc)，电子谱导数由有限体积端点差单独给出。
real(8) function electron_tau_kernel_x(x,V_cal,DB,factor)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x,V_cal,DB,factor
real(8) :: gam

    gam=dexp(x)
    electron_tau_kernel_x=gam*gam*electron_syn_fx(gam,V_cal,DB,factor)
end function electron_tau_kernel_x

! 单网格单元同步发射功率的2点+3点Gauss积分（用于自适应误差估计）。
subroutine electron_syn_gauss_cell(x0,x1,dN0,dN1,V_cal,DB,factor,p2,p3)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,dN0,dN1,V_cal,DB,factor
real(8), intent(out) :: p2,p3
real(8) :: xm,dx,w2,w3a
real(8) :: x2a,x2b,x3a,x3b

    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    w2=one/dsqrt(3d0)
    x2a=xm-dx*w2
    x2b=xm+dx*w2
    p2=dx*(electron_syn_integrand_x(x2a,x0,x1,dN0,dN1,V_cal,DB,factor)+ &
           electron_syn_integrand_x(x2b,x0,x1,dN0,dN1,V_cal,DB,factor))

    w3a=dsqrt(3d0/5d0)
    x3a=xm-dx*w3a
    x3b=xm+dx*w3a
    p3=dx*((5d0/9d0)*electron_syn_integrand_x(x3a,x0,x1,dN0,dN1,V_cal,DB,factor)+ &
           (8d0/9d0)*electron_syn_integrand_x(xm ,x0,x1,dN0,dN1,V_cal,DB,factor)+ &
           (5d0/9d0)*electron_syn_integrand_x(x3b,x0,x1,dN0,dN1,V_cal,DB,factor))
end subroutine electron_syn_gauss_cell

! 单网格单元SSA光深的2点+3点Gauss积分：有限体积端点差给出 -dq，Gauss积分给核函数平均。
subroutine electron_tau_gauss_cell(x0,x1,dN10,dN11,V_cal,DB,factor,t2,t3)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,dN10,dN11,V_cal,DB,factor
real(8), intent(out) :: t2,t3
real(8) :: xm,dx,width,w2,w3a,q_drop
real(8) :: x2a,x2b,x3a,x3b

    if (x1 <= x0) then
        t2=zero
        t3=zero
        return
    end if

    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    width=x1-x0
    q_drop=dN10-dN11
    w2=one/dsqrt(3d0)
    x2a=xm-dx*w2
    x2b=xm+dx*w2
    t2=q_drop/width*dx*(electron_tau_kernel_x(x2a,V_cal,DB,factor)+ &
                        electron_tau_kernel_x(x2b,V_cal,DB,factor))

    w3a=dsqrt(3d0/5d0)
    x3a=xm-dx*w3a
    x3b=xm+dx*w3a
    t3=q_drop/width*dx*((5d0/9d0)*electron_tau_kernel_x(x3a,V_cal,DB,factor)+ &
                        (8d0/9d0)*electron_tau_kernel_x(xm ,V_cal,DB,factor)+ &
                        (5d0/9d0)*electron_tau_kernel_x(x3b,V_cal,DB,factor))
end subroutine electron_tau_gauss_cell

! 自适应同步辐射单元积分：发射率用线性重构，SSA光深用保守端点差重构。
subroutine electron_syn_cell_adaptive(x0,x1,dN0,dN1,dN10,dN11, &
                                      V_cal,DB,factor,rel_tol,p_int,tau_int)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,dN0,dN1,dN10,dN11,V_cal,DB,factor,rel_tol
real(8), intent(out) :: p_int,tau_int
real(8) :: p2,p3,t2,t3,xm,dNm,dN1m,err_p,err_t,ref_p,ref_t
real(8) :: p3_l,p3_r,t3_l,t3_r

    call electron_syn_gauss_cell(x0,x1,dN0,dN1,V_cal,DB,factor,p2,p3)
    call electron_tau_gauss_cell(x0,x1,dN10,dN11,V_cal,DB,factor,t2,t3)
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
        call electron_syn_gauss_cell(x0,xm,dN0,dNm,V_cal,DB,factor,p2,p3_l)
        call electron_syn_gauss_cell(xm,x1,dNm,dN1,V_cal,DB,factor,p2,p3_r)
        call electron_tau_gauss_cell(x0,xm,dN10,dN1m,V_cal,DB,factor,t2,t3_l)
        call electron_tau_gauss_cell(xm,x1,dN1m,dN11,V_cal,DB,factor,t2,t3_r)
        p_int=p3_l+p3_r
        tau_int=t3_l+t3_r
    end if
end subroutine electron_syn_cell_adaptive

! 自适应同步辐射计算：发射率自适应积分，SSA光深按有限体积端点差守恒积分。
subroutine get_syn_adaptive_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                  P_emit,P_syn,Seed_syn,Tau_syn)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: P_emit(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu),Tau_syn(Num_nu)
real(8), parameter :: rel_tol=5d-4
real(8), allocatable :: dN1(:),x_gam(:)
real(8) :: factor,Temp_syn,Rariv2,temp_para
integer :: I_nu

    allocate(dN1(Num_gam_e),x_gam(Num_gam_e))

    factor=(3.62d0/pi)**2
    Temp_syn=dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    Rariv2=R_loc*R_loc
    dN1=dN_gam_e/(gam_e*gam_e)
    x_gam=dlog(gam_e)

    !$OMP PARALLEL num_threads(n_threads), private(I_nu)
    !$OMP DO SCHEDULE(STATIC)
    do I_nu=1,Num_nu
       call accumulate_adaptive_syn_point(I_nu)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    temp_para=4d0*pi*Para_c*Para_h
    Seed_syn=Seed_syn/temp_para

    deallocate(dN1,x_gam)

    return

contains

subroutine adaptive_syn_integrals(V_cal,dInteg,Tau)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: V_cal
real(8), intent(out) :: dInteg,Tau
real(8) :: cell_int,tau_cell
integer :: I_gam_e

    dInteg=zero
    Tau=zero
    do I_gam_e=1,Num_gam_e-1
       call electron_syn_cell_adaptive(x_gam(I_gam_e),x_gam(I_gam_e+1), &
                                       dN_gam_e(I_gam_e),dN_gam_e(I_gam_e+1), &
                                       dN1(I_gam_e),dN1(I_gam_e+1), &
                                       V_cal,DB,factor,rel_tol,cell_int,tau_cell)
       dInteg=dInteg+cell_int
       Tau=Tau+tau_cell
    end do
end subroutine adaptive_syn_integrals

subroutine accumulate_adaptive_syn_point(I_nu)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: I_nu
real(8) :: V_cal,dInteg,Tau,P_v

    V_cal=V_seed(I_nu)
    call adaptive_syn_integrals(V_cal,dInteg,Tau)
    P_v=Temp_syn*DB*dInteg
    Tau=1.046d4*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal)
    if ((Tau-1d-4) < 1d-5) Tau=1d-4
    P_emit(I_nu)=P_v
    Tau_syn(I_nu)=Tau
    P_syn(I_nu)=P_v*(one-dexp(-Tau))/Tau
    Seed_syn(I_nu)=P_syn(I_nu)/(Rariv2*V_cal)
end subroutine accumulate_adaptive_syn_point
end subroutine get_syn_adaptive_state

subroutine get_syn_adaptive(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                            P_syn,Seed_syn)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: P_syn(Num_nu),Seed_syn(Num_nu)
real(8) :: P_emit(Num_nu),Tau_syn(Num_nu)

    call get_syn_adaptive_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                P_emit,P_syn,Seed_syn,Tau_syn)
end subroutine get_syn_adaptive

    ! 同步辐射计算选择器：index=4 在标准同步核之外加入非相对论回旋基频发射。
subroutine get_syn_selected(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                            P_syn,Seed_syn)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_syn_intger,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: P_syn(Num_nu),Seed_syn(Num_nu)
real(8) :: P_emit(Num_nu),Tau_syn(Num_nu)

    call get_syn_selected_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                P_emit,P_syn,Seed_syn,Tau_syn)
end subroutine get_syn_selected

subroutine get_syn_selected_state(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                  P_emit,P_syn,Seed_syn,Tau_syn)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_syn_intger,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: P_emit(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu),Tau_syn(Num_nu)
real(8) :: h_ref,h_loc

    select case(index_syn_intger)
    case(1)
        if (Num_gam_e > 2) then
            h_ref=dlog(gam_e(2))-dlog(gam_e(1))
            do I_gam_e=3,Num_gam_e
                h_loc=dlog(gam_e(I_gam_e))-dlog(gam_e(I_gam_e-1))
                if (abs(h_loc-h_ref) > 1d-6*max(abs(h_ref),1d-30)) then
                    call get_syn_adaptive_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                                P_emit,P_syn,Seed_syn,Tau_syn)
                    return
                end if
            end do
        end if
        call get_syn_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,P_emit,P_syn,Seed_syn,Tau_syn)
    case(2)
        call get_syn_simpson_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                   P_emit,P_syn,Seed_syn,Tau_syn)
    case(3)
        call get_syn_adaptive_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                    P_emit,P_syn,Seed_syn,Tau_syn)
    case(4)
        call get_syn_cyclotron_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                                     P_emit,P_syn,Seed_syn,Tau_syn)
    case default
        print*, 'invalid synchrotron integral case, check your chosen model!'
        stop
    end select
end subroutine get_syn_selected_state

! 同步辐射偏振核：现有总谱给强度，F/G偏振核直接积分给频率依赖Pi。
subroutine get_syn_polarization_selected(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads, &
                                         gam_e,dN_gam_e,V_seed,p_index,P_perp,P_parallel,Pi_nu)
implicit none
integer, intent(in) :: index_syn_intger,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu),p_index
real(8), intent(out) :: P_perp(Num_nu),P_parallel(Num_nu),Pi_nu(Num_nu)
real(8) :: P_syn(Num_nu),Seed_syn(Num_nu),Pi_emit(Num_nu)

    if (p_index <= zero) error stop "get_syn_polarization_selected requires p_index > 0."
    call get_syn_selected(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                          P_syn,Seed_syn)
    call get_syn_polarization_fraction(DB,Num_gam_e,Num_nu,gam_e,dN_gam_e,V_seed,Pi_emit)
    Pi_nu=Pi_emit
    P_perp=0.5d0*(one+Pi_nu)*P_syn
    P_parallel=0.5d0*(one-Pi_nu)*P_syn
end subroutine get_syn_polarization_selected

! 输出专用局域偏振率核：对电子谱直接卷积(F+G)/2和(F-G)/2两个偏振发射核。
subroutine get_syn_polarization_fraction(DB,Num_gam_e,Num_nu,gam_e,dN_gam_e,V_seed,Pi_nu)
implicit none
integer, intent(in) :: Num_gam_e,Num_nu
real(8), intent(in) :: DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: Pi_nu(Num_nu)
integer :: I_nu,I_gam_e
real(8) :: V_cal,dPerp,dParallel,total_pol,gam_e_mean2,Vc,x,dN,dgam_e,perp_kernel,parallel_kernel

    if (DB <= zero) error stop "get_syn_polarization_fraction requires DB > 0."
    do I_nu=1,Num_nu
        V_cal=V_seed(I_nu)
        if (V_cal <= zero) error stop "get_syn_polarization_fraction requires positive frequency."
        dPerp=zero
        dParallel=zero
        do I_gam_e=1,Num_gam_e-1
            dN=(dN_gam_e(I_gam_e)+dN_gam_e(I_gam_e+1))/two
            if (dN <= zero) cycle
            gam_e_mean2=(gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0
            Vc=(4.2d6)*gam_e_mean2*DB
            x=V_cal/Vc
            dgam_e=gam_e(I_gam_e+1)-gam_e(I_gam_e)
            call synchrotron_polarized_components(x,perp_kernel,parallel_kernel)
            dPerp=dPerp+dN*perp_kernel*dgam_e
            dParallel=dParallel+dN*parallel_kernel*dgam_e
        end do
        total_pol=dPerp+dParallel
        if (total_pol > zero) then
            Pi_nu(I_nu)=(dPerp-dParallel)/total_pol
        else
            Pi_nu(I_nu)=zero
        end if
    end do
end subroutine get_syn_polarization_fraction

! 计算同步辐射转移函数：Transfer = P_absorbed / P_emit，即(1-e^(-τ))/τ。
subroutine get_syn_transfer(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,Transfer_syn)
implicit none
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: Transfer_syn(Num_nu)
real(8) :: P_emit(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu),Tau_syn(Num_nu)
integer :: I_nu

    call get_syn_state(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,P_emit,P_syn,Seed_syn,Tau_syn)
    do I_nu=1,Num_nu
        if (P_emit(I_nu) > 0d0 .and. P_syn(I_nu) >= 0d0) then
            Transfer_syn(I_nu)=P_syn(I_nu)/P_emit(I_nu)
        else
            Transfer_syn(I_nu)=one
        end if
    end do
end subroutine get_syn_transfer

! 计算SSA频率ν_a：二分搜索+对数割线法求解τ(ν_a)=1。
subroutine get_nu_a(R_loc,DB,Num_gam_e,gam_e,dN_gam_e, &
                       V_a)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e)
real(8), intent(out) :: V_a
real(8), parameter :: rel_tol=5d-4

real(8),allocatable,dimension (:) :: dN1,x_gam
allocate (dN1(Num_gam_e),x_gam(Num_gam_e))

    factor=(3.62d0/pi)**2
    Rariv2=R_loc*R_loc
    dN1=dN_gam_e/(gam_e*gam_e)
    x_gam=dlog(gam_e)
    V_a_floor=ten**4d0
    V_a_min=ten**(-20d0)
    V_a_cap=one
    do I_gam_e=1,Num_gam_e-1
       Vc=4.2d6*((gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0)*DB
       if (Vc > V_a_cap) V_a_cap=Vc
    end do
    V_a_cap=max(ten**14d0, min(ten**30d0, ten*V_a_cap))

    V_low=V_a_floor
    call evaluate_tau(V_low,Tau_low)
    if (Tau_low <= one) then
        do I_nu=1,26
           if (V_low <= V_a_min) exit
           V_high=V_low
           Tau_high=Tau_low
           V_low=max(V_a_min,V_low/ten)
           call evaluate_tau(V_low,Tau_low)
           if (Tau_low > one) exit
        end do

        if (Tau_low > one) then
           call refine_nu_a_bracket(V_low,Tau_low,V_high,Tau_high,V_a)
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
          V_high=min(V_a_cap, ten*V_high)
          call evaluate_tau(V_high,Tau_high)
          if (Tau_high <= one) exit
       end do

       if (Tau_high > one) then
          V_a=V_high
          print*, 'nu_a_comoving larger than adaptive upper bound!', V_a_cap
       else
          if (Tau_low <= one .or. Tau_low == Tau_high) then
             V_a=V_high
          else
             call refine_nu_a_bracket(V_low,Tau_low,V_high,Tau_high,V_a)
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

    Tau=zero
    do I_gam_e=1,Num_gam_e-1
       call electron_syn_cell_adaptive(x_gam(I_gam_e),x_gam(I_gam_e+1), &
                                       dN_gam_e(I_gam_e),dN_gam_e(I_gam_e+1), &
                                       dN1(I_gam_e),dN1(I_gam_e+1), &
                                       V_cal,DB,factor,rel_tol,discard,tau_cell)
       Tau=Tau+tau_cell
    end do
    Tau=1.046d4*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal)
end subroutine evaluate_tau

! 对数空间割线法迭代精化ν_a：在[V_low,V_high]内求解τ(ν)=1，最多10次迭代。
subroutine refine_nu_a_bracket(V_low,Tau_low,V_high,Tau_high,V_root)
implicit REAL(8)(A-H,O-Z)
real(8), intent(inout) :: V_low,Tau_low,V_high,Tau_high
real(8), intent(out) :: V_root
real(8) :: log_v_low,log_v_high,log_v_mid,log_tau_low,log_tau_high
real(8) :: V_mid,Tau_mid
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

        if (Tau_mid > one) then
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
end subroutine refine_nu_a_bracket
end subroutine get_nu_a

! 2D路径积分ν_a：累加所有χ列的光深后求解τ_total(ν_a)=1。
subroutine get_nu_a_2d_path(Num_nu,Num_chi,V_seed,Tau_chi,V_a)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu,Num_chi
integer :: I_nu
real(8), intent(in) :: V_seed(Num_nu),Tau_chi(Num_nu,Num_chi)
real(8), intent(out) :: V_a
real(8) :: Tau_path(Num_nu)

    Tau_path=zero
    do I_nu=1,Num_nu
        Tau_path(I_nu)=sum(Tau_chi(I_nu,:))
    end do

    call get_nu_a_from_tau_grid(Num_nu,V_seed,Tau_path,V_a)
end subroutine get_nu_a_2d_path

! 2D逐列ν_a：对每个χ列累加光深后求解τ_cumul(ν_a)=1，输出各列ν_a(χ)。
subroutine get_nu_a_2d_cell_path(Num_nu,Num_chi,V_seed,Tau_chi,V_a_chi)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu,Num_chi
integer :: I_chi
real(8), intent(in) :: V_seed(Num_nu),Tau_chi(Num_nu,Num_chi)
real(8), intent(out) :: V_a_chi(Num_chi)
real(8) :: Tau_path(Num_nu)

    Tau_path=zero
    do I_chi=1,Num_chi
        Tau_path=Tau_path+Tau_chi(:,I_chi)
        call get_nu_a_from_tau_grid(Num_nu,V_seed,Tau_path,V_a_chi(I_chi))
    end do
end subroutine get_nu_a_2d_cell_path

! 将χ分辨的同步辐射谱加权求和为壳层平均谱：Σ χ dη ln(10) * Q(χ)。
subroutine reduce_syn_shell_from_chi(Num_nu,Num_chi,deta,chi_grid,P_chi,Seed_chi,P_shell,Seed_shell)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu,Num_chi
integer :: I_chi
real(8), intent(in) :: deta,chi_grid(Num_chi),P_chi(Num_nu,Num_chi),Seed_chi(Num_nu,Num_chi)
real(8), intent(out) :: P_shell(Num_nu),Seed_shell(Num_nu)
real(8) :: weight,ln10

    ln10=dlog(ten)
    P_shell=zero
    Seed_shell=zero
    do I_chi=1,Num_chi
        weight=chi_grid(I_chi)*deta*ln10
        P_shell=P_shell+weight*P_chi(:,I_chi)
        Seed_shell=Seed_shell+weight*Seed_chi(:,I_chi)
    end do
end subroutine reduce_syn_shell_from_chi

! 从光深网格求解ν_a：在τ穿越1的位置对数插值。
subroutine get_nu_a_from_tau_grid(Num_nu,V_seed,Tau_grid,V_a)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu
integer :: I_nu,I_cross
real(8), intent(in) :: V_seed(Num_nu),Tau_grid(Num_nu)
real(8), intent(out) :: V_a
real(8) :: Tau_left,Tau_right

    if (Tau_grid(1) <= one) then
        V_a=V_seed(1)
        return
    end if

    if (Tau_grid(Num_nu) >= one) then
        V_a=V_seed(Num_nu)
        return
    end if

    I_cross=0
    do I_nu=2,Num_nu
        if (Tau_grid(I_nu) <= one) then
            I_cross=I_nu
            exit
        end if
    end do

    if (I_cross == 0) then
        call interpolate_log_tau_root(V_seed(Num_nu-1),max(Tau_grid(Num_nu-1),tiny(one)), &
                                      V_seed(Num_nu),max(Tau_grid(Num_nu),tiny(one)),V_a)
        if (V_a < V_seed(Num_nu)) V_a=V_seed(Num_nu)
        return
    end if

    Tau_left=max(Tau_grid(I_cross-1),tiny(one))
    Tau_right=max(Tau_grid(I_cross),tiny(one))
    call interpolate_log_tau_root(V_seed(I_cross-1),Tau_left,V_seed(I_cross),Tau_right,V_a)
end subroutine get_nu_a_from_tau_grid

! 对数-对数空间线性插值求解τ(ν)=1的根ν_a。
subroutine interpolate_log_tau_root(V_left,Tau_left,V_right,Tau_right,V_root)
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
end subroutine interpolate_log_tau_root

end module electron_radiation_kernel
