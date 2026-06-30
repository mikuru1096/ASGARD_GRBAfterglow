!f2py: skip
module electron_cooling_ssa_kernel
  use constants
  use electron_radiation_kernel, only: first_greater_monotonic_window, electron_powerlaw_interp, &
                                       electron_log_gauss2_interval
  private

  public :: electron_cooling_ssa_loss_batch

  integer, parameter :: ssa_segment_low = 1
  integer, parameter :: ssa_segment_high = 2
  integer, save :: ssa_seed_num_nu_cache=0
  logical, save :: ssa_seed_cache_ready=.false.
  real(8), allocatable, save :: ssa_seed_v_cache(:), ssa_seed_v_low(:), ssa_seed_v_high(:), &
                                ssa_seed_vg1(:), ssa_seed_vg2(:), ssa_seed_wg1(:), ssa_seed_wg2(:)
  !$omp threadprivate(ssa_seed_num_nu_cache,ssa_seed_cache_ready,ssa_seed_v_cache,ssa_seed_v_low,ssa_seed_v_high)
  !$omp threadprivate(ssa_seed_vg1,ssa_seed_vg2,ssa_seed_wg1,ssa_seed_wg2)

contains
! 确保SSA种子频率缓存已按当前种子网格分配。
subroutine ensure_ssa_seed_cache(Num_nu,V_seed)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu
real(8), intent(in) :: V_seed(Num_nu)
logical :: cache_match

    cache_match=.false.
    if (ssa_seed_cache_ready) then
        if (ssa_seed_num_nu_cache == Num_nu) then
            if (allocated(ssa_seed_v_cache)) then
                cache_match=all(ssa_seed_v_cache == V_seed)
            end if
        end if
    end if

    if (cache_match) return

    if (allocated(ssa_seed_v_cache)) deallocate(ssa_seed_v_cache, ssa_seed_v_low, ssa_seed_v_high, &
                                                ssa_seed_vg1, ssa_seed_vg2, ssa_seed_wg1, ssa_seed_wg2)
    allocate(ssa_seed_v_cache(Num_nu), ssa_seed_v_low(Num_nu-1), ssa_seed_v_high(Num_nu-1), &
             ssa_seed_vg1(Num_nu-1), ssa_seed_vg2(Num_nu-1), ssa_seed_wg1(Num_nu-1), ssa_seed_wg2(Num_nu-1))
    ssa_seed_v_cache=V_seed
    ssa_seed_v_low=V_seed(1:Num_nu-1)
    ssa_seed_v_high=V_seed(2:Num_nu)
    do I_nu=1,Num_nu-1
        call electron_log_gauss2_interval(ssa_seed_v_low(I_nu),ssa_seed_v_high(I_nu), &
                                          ssa_seed_vg1(I_nu),ssa_seed_vg2(I_nu), &
                                          ssa_seed_wg1(I_nu),ssa_seed_wg2(I_nu))
    end do
    ssa_seed_num_nu_cache=Num_nu
    ssa_seed_cache_ready=.true.
end subroutine ensure_ssa_seed_cache

! 推进SSA种子频率游标至第一个 > V_lowlim 的位置。
subroutine advance_ssa_seed_cursor(Num_nu,V_lowlim,low_idx)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_nu
real(8), intent(in) :: V_lowlim
integer, intent(inout) :: low_idx

    if (low_idx < 1) low_idx=1
    if (low_idx > Num_nu+1) low_idx=Num_nu+1

    do while (low_idx > 1)
        if (ssa_seed_v_cache(low_idx-1) <= V_lowlim) exit
        low_idx=low_idx-1
    end do

    do while (low_idx <= Num_nu)
        if (ssa_seed_v_cache(low_idx) > V_lowlim) exit
        low_idx=low_idx+1
    end do
end subroutine advance_ssa_seed_cursor

! 构建SSA几何映射：对每个电子γ，计算种子频率的低频/高频索引范围和截面前因子。
subroutine build_ssa_geometry(DB,Num_gam_e,Num_nu,gam_e,V_low_idx,V_low_first,V_low_last,V_high_first, &
                              sigma_prefactor_low,sigma_prefactor_high,Cyclotron_nu)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu
real(8), intent(in) :: DB,gam_e(Num_gam_e)
integer, intent(out) :: V_low_idx(Num_gam_e),V_low_first(Num_gam_e),V_low_last(Num_gam_e),V_high_first(Num_gam_e)
real(8), intent(out) :: sigma_prefactor_low(Num_gam_e),sigma_prefactor_high(Num_gam_e),Cyclotron_nu
integer :: low_idx,upper_idx

    B_cr=4.4d13
    Temp1=2.5042d-22*B_cr/DB
    Temp2=7.787d-22*B_cr/DB
    Cyclotron_nu=para_e*DB/(two*pi*para_m_e*para_c)

    low_idx=1
    do I_gam_e=1,Num_gam_e
       gam=gam_e(I_gam_e)
       gam2=gam*gam
       gam3=gam2*gam
       V_lowlim=Cyclotron_nu/gam
       V_uplim=1.5d0*gam2*Cyclotron_nu

       call advance_ssa_seed_cursor(Num_nu,V_lowlim,low_idx)
       V_low_idx(I_gam_e)=low_idx

       if (low_idx <= Num_nu) then
          call first_greater_monotonic_window(ssa_seed_v_cache,Num_nu,low_idx,V_uplim,upper_idx)
          sigma_prefactor_low(I_gam_e)=Temp1*(3d0*V_lowlim)**(5d0/3d0)
          sigma_prefactor_high(I_gam_e)=Temp2/gam3
          V_low_first(I_gam_e)=max(1,low_idx-1)
          V_low_last(I_gam_e)=min(Num_nu-1,upper_idx-1)
          V_high_first(I_gam_e)=max(1,upper_idx-1)
       else
          sigma_prefactor_low(I_gam_e)=zero
          sigma_prefactor_high(I_gam_e)=zero
          V_low_first(I_gam_e)=1
          V_low_last(I_gam_e)=0
          V_high_first(I_gam_e)=Num_nu
       end if
    end do
end subroutine build_ssa_geometry

! SSA冷却率：对多个χ列的种子光子场同时计算；单列调用传 Num_chi=1。
subroutine electron_cooling_ssa_loss_batch(DB,Num_gam_e,Num_nu,Num_chi,n_threads,gam_e,V_seed,Seed_syn_batch,dot_gam_e_batch)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,Num_chi,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn_batch(Num_nu,Num_chi)
real(8), intent(out) :: dot_gam_e_batch(Num_gam_e,Num_chi)
integer, parameter :: parallel_work_threshold=512
integer :: V_low_idx(Num_gam_e),V_low_first(Num_gam_e),V_low_last(Num_gam_e),V_high_first(Num_gam_e)
integer :: I_nu,I_chi,I_gam_e,work_items,thread_count
logical :: use_parallel
real(8) :: V_seed_low(Num_nu-1),V_seed_high(Num_nu-1),V_seed_g1(Num_nu-1),V_seed_g2(Num_nu-1)
real(8) :: V_seed_w1(Num_nu-1),V_seed_w2(Num_nu-1),sigma_low(Num_gam_e),sigma_high(Num_gam_e)
real(8) :: Cyclotron_nu,Seed_g1(Num_nu-1,Num_chi),Seed_g2(Num_nu-1,Num_chi)
real(8) :: Low_prefix(0:Num_nu-1,Num_chi),High_amp1(Num_nu-1,Num_chi),High_amp2(Num_nu-1,Num_chi)

    call ensure_ssa_seed_cache(Num_nu,V_seed)
    V_seed_low=ssa_seed_v_low
    V_seed_high=ssa_seed_v_high
    V_seed_g1=ssa_seed_vg1
    V_seed_g2=ssa_seed_vg2
    V_seed_w1=ssa_seed_wg1
    V_seed_w2=ssa_seed_wg2
    Low_prefix(0,:)=zero
    do I_chi=1,Num_chi
        do I_nu=1,Num_nu-1
            Seed_g1(I_nu,I_chi)=electron_powerlaw_interp(V_seed_low(I_nu),V_seed_high(I_nu), &
                                                         Seed_syn_batch(I_nu,I_chi),Seed_syn_batch(I_nu+1,I_chi), &
                                                         V_seed_g1(I_nu))
            Seed_g2(I_nu,I_chi)=electron_powerlaw_interp(V_seed_low(I_nu),V_seed_high(I_nu), &
                                                         Seed_syn_batch(I_nu,I_chi),Seed_syn_batch(I_nu+1,I_chi), &
                                                         V_seed_g2(I_nu))
            Low_prefix(I_nu,I_chi)=Low_prefix(I_nu-1,I_chi)+para_h*para_c*(V_seed_w1(I_nu)* &
                                    Seed_g1(I_nu,I_chi)*V_seed_g1(I_nu)**(-2d0/3d0)+ &
                                    V_seed_w2(I_nu)*Seed_g2(I_nu,I_chi)*V_seed_g2(I_nu)**(-2d0/3d0))
            High_amp1(I_nu,I_chi)=V_seed_w1(I_nu)*Seed_g1(I_nu,I_chi)
            High_amp2(I_nu,I_chi)=V_seed_w2(I_nu)*Seed_g2(I_nu,I_chi)
        end do
    end do
    call build_ssa_geometry(DB,Num_gam_e,Num_nu,gam_e,V_low_idx,V_low_first,V_low_last,V_high_first, &
                            sigma_low,sigma_high,Cyclotron_nu)

    dot_gam_e_batch=zero
    work_items=Num_gam_e*Num_chi*Num_nu
    thread_count=max(1,n_threads)
    use_parallel=(n_threads > 1 .and. work_items >= parallel_work_threshold)
    !$OMP PARALLEL DO collapse(2) if(use_parallel) num_threads(thread_count) schedule(static) &
    !$OMP& private(I_chi,I_gam_e)
    do I_chi=1,Num_chi
        do I_gam_e=1,Num_gam_e
            call accumulate_ssa_batch_gamma(I_gam_e,I_chi,dot_gam_e_batch(I_gam_e,I_chi))
        end do
    end do
    !$OMP END PARALLEL DO

contains

subroutine accumulate_ssa_batch_gamma(I_gam_e,I_chi,dot_val)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: I_gam_e,I_chi
real(8), intent(out) :: dot_val
integer :: I_nu,low_full_first,low_full_last,high_full_first
real(8) :: gam,gam2,V_lowlim,V_uplim,inv_uplim,high_prefactor,ssa_sum,cell_low,cell_high

    dot_val=zero
    if (V_low_idx(I_gam_e) > Num_nu) return

    gam=gam_e(I_gam_e)
    gam2=gam*gam
    V_lowlim=Cyclotron_nu/gam
    V_uplim=1.5d0*gam2*Cyclotron_nu
    inv_uplim=one/V_uplim
    high_prefactor=sigma_high(I_gam_e)*Cyclotron_nu*para_h*para_c
    ssa_sum=zero

    low_full_first=V_low_first(I_gam_e)
    low_full_last=V_low_last(I_gam_e)
    if (low_full_first <= low_full_last) then
        cell_low=max(V_seed_low(low_full_first),V_lowlim)
        cell_high=min(V_seed_high(low_full_first),V_uplim)
        if (cell_high > cell_low) then
            if (cell_low /= V_seed_low(low_full_first) .or. cell_high /= V_seed_high(low_full_first)) then
                ssa_sum=ssa_sum+clipped_ssa_batch_segment(cell_low,cell_high,low_full_first,I_chi, &
                                                          sigma_low(I_gam_e),ssa_segment_low,Cyclotron_nu,V_uplim)
                low_full_first=low_full_first+1
            end if
        else
            low_full_first=low_full_first+1
        end if
    end if
    if (low_full_last >= low_full_first) then
        cell_low=max(V_seed_low(low_full_last),V_lowlim)
        cell_high=min(V_seed_high(low_full_last),V_uplim)
        if (cell_high > cell_low) then
            if (cell_low /= V_seed_low(low_full_last) .or. cell_high /= V_seed_high(low_full_last)) then
                ssa_sum=ssa_sum+clipped_ssa_batch_segment(cell_low,cell_high,low_full_last,I_chi, &
                                                          sigma_low(I_gam_e),ssa_segment_low,Cyclotron_nu,V_uplim)
                low_full_last=low_full_last-1
            end if
        else
            low_full_last=low_full_last-1
        end if
    end if
    if (low_full_last >= low_full_first) then
        ssa_sum=ssa_sum+sigma_low(I_gam_e)*(Low_prefix(low_full_last,I_chi)-Low_prefix(low_full_first-1,I_chi))
    end if

    high_full_first=V_high_first(I_gam_e)
    if (high_full_first <= Num_nu-1) then
        cell_low=max(V_seed_low(high_full_first),V_uplim)
        cell_high=V_seed_high(high_full_first)
        if (cell_high > cell_low) then
            if (cell_low /= V_seed_low(high_full_first)) then
                ssa_sum=ssa_sum+clipped_ssa_batch_segment(cell_low,cell_high,high_full_first,I_chi, &
                                                          sigma_high(I_gam_e),ssa_segment_high,Cyclotron_nu,V_uplim)
                high_full_first=high_full_first+1
            end if
        end if
    end if
    do I_nu=high_full_first,Num_nu-1
        ssa_sum=ssa_sum+high_prefactor*(High_amp1(I_nu,I_chi)*dexp(-V_seed_g1(I_nu)*inv_uplim)+ &
                                        High_amp2(I_nu,I_chi)*dexp(-V_seed_g2(I_nu)*inv_uplim))
    end do

    dot_val=ssa_sum
end subroutine accumulate_ssa_batch_gamma

real(8) function clipped_ssa_batch_segment(cell_low,cell_high,I_nu,I_chi,sigma_prefactor,mode,Cyclotron_nu,V_uplim)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: I_nu,I_chi,mode
real(8), intent(in) :: cell_low,cell_high,sigma_prefactor,Cyclotron_nu,V_uplim
integer :: I_quad
real(8) :: seed_loc,sigma_loc,vg(2),wg(2)

    call electron_log_gauss2_interval(cell_low,cell_high,vg(1),vg(2),wg(1),wg(2))
    clipped_ssa_batch_segment=zero

    do I_quad=1,2
        seed_loc=electron_powerlaw_interp(V_seed_low(I_nu),V_seed_high(I_nu), &
                                          Seed_syn_batch(I_nu,I_chi),Seed_syn_batch(I_nu+1,I_chi),vg(I_quad))
        if (mode == ssa_segment_low) then
            sigma_loc=sigma_prefactor*vg(I_quad)**(-5d0/3d0)
        else
            sigma_loc=sigma_prefactor*(Cyclotron_nu/vg(I_quad))*dexp(-vg(I_quad)/V_uplim)
        end if
        clipped_ssa_batch_segment=clipped_ssa_batch_segment+ &
                                  wg(I_quad)*sigma_loc*seed_loc*para_h*vg(I_quad)*para_c
    end do
end function clipped_ssa_batch_segment
end subroutine electron_cooling_ssa_loss_batch

end module electron_cooling_ssa_kernel
