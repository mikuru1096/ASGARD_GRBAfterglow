!f2py: skip
module electron_cooling_kernel
  use constants
  use radiation_common, only: compute_simpson_weights
  use electron_radiation_kernel, only: first_greater_monotonic, first_greater_monotonic_window, &
                                       electron_powerlaw_interp, electron_integrate_powerlaw_segment, &
                                       electron_log_gauss2_interval, electron_ssa_segment
  private

  public :: electron_cooling_ssa_loss, electron_cooling_ssa_loss_batch, electron_cooling_ic_loss
  public :: electron_cooling_ic_loss_emissivity_budget
  public :: electron_cooling_y_nakar, electron_cooling_y_fan
  public :: get_forward_cooling, prepare_forward_cooling_aux, prepare_forward_cooling_aux_batch
  public :: assemble_forward_cooling_split, assemble_forward_cooling_split_batch

  integer, save :: ssa_seed_num_nu_cache=0
  logical, save :: ssa_seed_cache_ready=.false.
  real(8), allocatable, save :: ssa_seed_v_cache(:), ssa_seed_v_low(:), ssa_seed_v_high(:), &
                                ssa_seed_vg1(:), ssa_seed_vg2(:), ssa_seed_wg1(:), ssa_seed_wg2(:)
  integer, save :: ssa_geom_num_gam_cache=0
  integer, save :: ssa_geom_num_chi_cache=0
  integer, allocatable, save :: ssa_low_idx_cache(:), ssa_low_first_cache(:), ssa_low_last_cache(:), ssa_high_first_cache(:)
  real(8), allocatable, save :: ssa_prefactor_low_cache(:), ssa_prefactor_high_cache(:), ssa_dot_batch_cache(:,:)
  integer, save :: ic_num_gam_cache=0, ic_num_nu_cache=0
  logical, save :: ic_grid_cache_ready=.false.
  real(8), allocatable, save :: ic_d_nu_cache(:), ic_gam_e_mean_cache(:), &
                                ic_e_seed_cache(:), ic_x_seed_cache(:), ic_v_seed_mid_cache(:)
  integer, save :: y_nakar_num_gam_cache=0, y_nakar_num_nu_cache=0
  integer, allocatable, save :: y_nakar_idx_cache(:)
  real(8), allocatable, save :: y_nakar_hat_nu_cache(:), y_nakar_prefix_cache(:), y_nakar_v_cache(:), y_nakar_vloc_cache(:), &
                                y_nakar_vg1(:), y_nakar_vg2(:), y_nakar_wg1(:), y_nakar_wg2(:)

contains

! 确保SSA几何工作数组已按当前网格大小分配（含缓存检查）。
subroutine ensure_ssa_geometry_workspace(Num_gam_e,Num_chi)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_chi

    if (.not. allocated(ssa_low_idx_cache) .or. ssa_geom_num_gam_cache /= Num_gam_e) then
        if (allocated(ssa_low_idx_cache)) then
            deallocate(ssa_low_idx_cache,ssa_low_first_cache,ssa_low_last_cache,ssa_high_first_cache, &
                       ssa_prefactor_low_cache,ssa_prefactor_high_cache)
        end if
        allocate(ssa_low_idx_cache(Num_gam_e),ssa_low_first_cache(Num_gam_e),ssa_low_last_cache(Num_gam_e), &
                 ssa_high_first_cache(Num_gam_e),ssa_prefactor_low_cache(Num_gam_e),ssa_prefactor_high_cache(Num_gam_e))
        ssa_geom_num_gam_cache=Num_gam_e
    end if

    if (Num_chi > 0) then
        if (.not. allocated(ssa_dot_batch_cache) .or. ssa_geom_num_gam_cache /= Num_gam_e .or. &
            ssa_geom_num_chi_cache /= Num_chi) then
            if (allocated(ssa_dot_batch_cache)) deallocate(ssa_dot_batch_cache)
            allocate(ssa_dot_batch_cache(Num_gam_e,Num_chi))
            ssa_geom_num_chi_cache=Num_chi
        end if
    end if
end subroutine ensure_ssa_geometry_workspace

! 确保IC网格缓存已计算（种子频率中点值、间距、电子能量中点值等）。
subroutine ensure_ic_grid_cache(Num_gam_e,Num_nu,gam_e,V_seed)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu)
logical :: rebuild

    rebuild=.not. ic_grid_cache_ready
    if (.not. rebuild) then
        rebuild = (ic_num_gam_cache /= Num_gam_e) .or. (ic_num_nu_cache /= Num_nu)
    end if
    if (.not. rebuild) then
        if (allocated(ic_x_seed_cache)) rebuild = any(ic_x_seed_cache /= dlog(V_seed))
        if (.not. rebuild .and. allocated(ic_gam_e_mean_cache)) rebuild = &
            any(ic_gam_e_mean_cache /= (gam_e(1:Num_gam_e-1)+gam_e(2:Num_gam_e))/two)
    end if
    if (.not. rebuild) return

    if (allocated(ic_d_nu_cache)) deallocate(ic_d_nu_cache,ic_gam_e_mean_cache,ic_e_seed_cache,ic_x_seed_cache,ic_v_seed_mid_cache)
    allocate(ic_d_nu_cache(Num_nu-1),ic_gam_e_mean_cache(Num_gam_e-1),ic_e_seed_cache(Num_nu-1), &
             ic_x_seed_cache(Num_nu),ic_v_seed_mid_cache(Num_nu-1))

    para_hEme = Para_h/para_m_energy
    ic_x_seed_cache=dlog(V_seed)
    ic_v_seed_mid_cache=dexp(0.5d0*(ic_x_seed_cache(1:Num_nu-1)+ic_x_seed_cache(2:Num_nu)))
    ic_d_nu_cache=ic_v_seed_mid_cache*(ic_x_seed_cache(2:Num_nu)-ic_x_seed_cache(1:Num_nu-1))
    ic_gam_e_mean_cache=(gam_e(1:Num_gam_e-1)+gam_e(2:Num_gam_e))/two
    ic_e_seed_cache=ic_v_seed_mid_cache*para_hEme
    ic_num_gam_cache=Num_gam_e
    ic_num_nu_cache=Num_nu
    ic_grid_cache_ready=.true.
end subroutine ensure_ic_grid_cache

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

! 用缓存的完整频率单元log-Gauss节点计算SSA段积分；裁剪段仍使用通用electron_ssa_segment。
real(8) function ssa_cached_cell_segment(I_nu,seed0,seed1,sigma_prefactor,mode,Cyclotron_nu,V_uplim)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: I_nu,mode
real(8), intent(in) :: seed0,seed1,sigma_prefactor,Cyclotron_nu,V_uplim
real(8) :: seed_loc,sigma_loc,vg,wg

    vg=ssa_seed_vg1(I_nu)
    wg=ssa_seed_wg1(I_nu)
    seed_loc=electron_powerlaw_interp(ssa_seed_v_low(I_nu),ssa_seed_v_high(I_nu),seed0,seed1,vg)
    if (mode == 1) then
        sigma_loc=sigma_prefactor*vg**(-5d0/3d0)
    else
        sigma_loc=sigma_prefactor*(Cyclotron_nu/vg)*dexp(-vg/V_uplim)
    end if
    ssa_cached_cell_segment=wg*sigma_loc*seed_loc*para_h*vg*para_c

    vg=ssa_seed_vg2(I_nu)
    wg=ssa_seed_wg2(I_nu)
    seed_loc=electron_powerlaw_interp(ssa_seed_v_low(I_nu),ssa_seed_v_high(I_nu),seed0,seed1,vg)
    if (mode == 1) then
        sigma_loc=sigma_prefactor*vg**(-5d0/3d0)
    else
        sigma_loc=sigma_prefactor*(Cyclotron_nu/vg)*dexp(-vg/V_uplim)
    end if
    ssa_cached_cell_segment=ssa_cached_cell_segment+wg*sigma_loc*seed_loc*para_h*vg*para_c
end function ssa_cached_cell_segment

! 对给定种子光子场累加SSA冷却率：低频Σ ∝ ν^(5/3)，高频Σ ∝ 1/γ³。
subroutine accumulate_ssa_for_seed(Num_gam_e,Num_nu,n_threads,gam_e,Seed_syn,V_low_idx,V_low_first,V_low_last, &
                                   V_high_first,sigma_prefactor_low,sigma_prefactor_high,Cyclotron_nu,dot_gam_e)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),Seed_syn(Num_nu)
integer, intent(in) :: V_low_idx(Num_gam_e),V_low_first(Num_gam_e),V_low_last(Num_gam_e),V_high_first(Num_gam_e)
real(8), intent(in) :: sigma_prefactor_low(Num_gam_e),sigma_prefactor_high(Num_gam_e),Cyclotron_nu
real(8), intent(out) :: dot_gam_e(Num_gam_e)
integer, parameter :: parallel_work_threshold=512
integer :: work_items

    dot_gam_e=zero

    work_items=Num_gam_e*Num_nu
    if (n_threads <= 1 .or. work_items < parallel_work_threshold) then
       do I_gam_e=1,Num_gam_e
          call accumulate_ssa_single_gamma(I_gam_e,dot_gam_e(I_gam_e))
       end do
    else
       !$OMP PARALLEL DO num_threads(n_threads) schedule(static) &
       !$OMP& private(I_gam_e)
       do I_gam_e=1,Num_gam_e
          call accumulate_ssa_single_gamma(I_gam_e,dot_gam_e(I_gam_e))
       end do
       !$OMP END PARALLEL DO
    end if

contains

subroutine accumulate_ssa_single_gamma(I_gam_e,dot_val)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: I_gam_e
real(8), intent(out) :: dot_val
integer :: I_nu
real(8) :: gam,gam2,V_lowlim,V_uplim,ssa_sum,cell_low,cell_high

    dot_val=zero
    if (V_low_idx(I_gam_e) > Num_nu) return

    gam=gam_e(I_gam_e)
    gam2=gam*gam
    V_lowlim=Cyclotron_nu/gam
    V_uplim=1.5d0*gam2*Cyclotron_nu
    ssa_sum=zero

    do I_nu=V_low_first(I_gam_e),V_low_last(I_gam_e)
       cell_low=max(ssa_seed_v_low(I_nu),V_lowlim)
       cell_high=min(ssa_seed_v_high(I_nu),V_uplim)
       if (cell_high > cell_low) then
          if (cell_low == ssa_seed_v_low(I_nu) .and. cell_high == ssa_seed_v_high(I_nu)) then
              ssa_sum=ssa_sum+ssa_cached_cell_segment(I_nu,Seed_syn(I_nu),Seed_syn(I_nu+1), &
                                                      sigma_prefactor_low(I_gam_e),1,Cyclotron_nu,V_uplim)
          else
              ssa_sum=ssa_sum+electron_ssa_segment(cell_low,cell_high,Seed_syn(I_nu),Seed_syn(I_nu+1), &
                                                   sigma_prefactor_low(I_gam_e),1,Cyclotron_nu,V_uplim)
          end if
       end if
    end do

    do I_nu=V_high_first(I_gam_e),Num_nu-1
       cell_low=max(ssa_seed_v_low(I_nu),V_uplim)
       cell_high=ssa_seed_v_high(I_nu)
       if (cell_high > cell_low) then
          if (cell_low == ssa_seed_v_low(I_nu)) then
              ssa_sum=ssa_sum+ssa_cached_cell_segment(I_nu,Seed_syn(I_nu),Seed_syn(I_nu+1), &
                                                      sigma_prefactor_high(I_gam_e),2,Cyclotron_nu,V_uplim)
          else
              ssa_sum=ssa_sum+electron_ssa_segment(cell_low,cell_high,Seed_syn(I_nu),Seed_syn(I_nu+1), &
                                                   sigma_prefactor_high(I_gam_e),2,Cyclotron_nu,V_uplim)
          end if
       end if
    end do

    dot_val=ssa_sum
end subroutine accumulate_ssa_single_gamma
end subroutine accumulate_ssa_for_seed

! 计算同步自吸收（SSA）冷却率 dγ/dt：构建几何→累加低频+高频贡献。
subroutine electron_cooling_ssa_loss(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn, dot_gam_e)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dot_gam_e(Num_gam_e)
real(8) :: Cyclotron_nu

    call ensure_ssa_seed_cache(Num_nu,V_seed)
    call ensure_ssa_geometry_workspace(Num_gam_e,0)
    call build_ssa_geometry(DB,Num_gam_e,Num_nu,gam_e,ssa_low_idx_cache,ssa_low_first_cache,ssa_low_last_cache, &
                            ssa_high_first_cache, &
                            ssa_prefactor_low_cache,ssa_prefactor_high_cache,Cyclotron_nu)
    call accumulate_ssa_for_seed(Num_gam_e,Num_nu,n_threads,gam_e,Seed_syn,ssa_low_idx_cache,ssa_low_first_cache, &
                                 ssa_low_last_cache,ssa_high_first_cache, &
                                 ssa_prefactor_low_cache,ssa_prefactor_high_cache,Cyclotron_nu,dot_gam_e)

end subroutine electron_cooling_ssa_loss

! 批量模式下计算单个(γ, χ)单元的SSA冷却率。
subroutine accumulate_ssa_batch_cell(Num_gam_e,Num_nu,Num_chi,I_gam_e,I_chi,gam_e,Seed_syn_batch,Cyclotron_nu,dot_gam_e_val)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,Num_chi,I_gam_e,I_chi
real(8), intent(in) :: gam_e(Num_gam_e),Seed_syn_batch(Num_nu,Num_chi)
real(8), intent(in) :: Cyclotron_nu
real(8), intent(out) :: dot_gam_e_val
integer :: I_nu
real(8) :: gam,gam2,V_lowlim,V_uplim,ssa_sum,cell_low,cell_high

    dot_gam_e_val=zero
    if (ssa_low_idx_cache(I_gam_e) > Num_nu) return

    gam=gam_e(I_gam_e)
    gam2=gam*gam
    V_lowlim=Cyclotron_nu/gam
    V_uplim=1.5d0*gam2*Cyclotron_nu
    ssa_sum=zero

    do I_nu=ssa_low_first_cache(I_gam_e),ssa_low_last_cache(I_gam_e)
        cell_low=max(ssa_seed_v_low(I_nu),V_lowlim)
        cell_high=min(ssa_seed_v_high(I_nu),V_uplim)
        if (cell_high > cell_low) then
            if (cell_low == ssa_seed_v_low(I_nu) .and. cell_high == ssa_seed_v_high(I_nu)) then
                ssa_sum=ssa_sum+ssa_cached_cell_segment(I_nu,Seed_syn_batch(I_nu,I_chi), &
                                                        Seed_syn_batch(I_nu+1,I_chi), &
                                                        ssa_prefactor_low_cache(I_gam_e),1,Cyclotron_nu,V_uplim)
            else
                ssa_sum=ssa_sum+electron_ssa_segment(cell_low,cell_high,Seed_syn_batch(I_nu,I_chi), &
                                                     Seed_syn_batch(I_nu+1,I_chi), &
                                                     ssa_prefactor_low_cache(I_gam_e),1,Cyclotron_nu,V_uplim)
            end if
        end if
    end do

    do I_nu=ssa_high_first_cache(I_gam_e),Num_nu-1
        cell_low=max(ssa_seed_v_low(I_nu),V_uplim)
        cell_high=ssa_seed_v_high(I_nu)
        if (cell_high > cell_low) then
            if (cell_low == ssa_seed_v_low(I_nu)) then
                ssa_sum=ssa_sum+ssa_cached_cell_segment(I_nu,Seed_syn_batch(I_nu,I_chi), &
                                                        Seed_syn_batch(I_nu+1,I_chi), &
                                                        ssa_prefactor_high_cache(I_gam_e),2,Cyclotron_nu,V_uplim)
            else
                ssa_sum=ssa_sum+electron_ssa_segment(cell_low,cell_high,Seed_syn_batch(I_nu,I_chi), &
                                                     Seed_syn_batch(I_nu+1,I_chi), &
                                                     ssa_prefactor_high_cache(I_gam_e),2,Cyclotron_nu,V_uplim)
            end if
        end if
    end do

    dot_gam_e_val=ssa_sum
end subroutine accumulate_ssa_batch_cell

! 批量SSA冷却率：对多个χ列的种子光子场同时计算SSA，可OpenMP并行。
subroutine electron_cooling_ssa_loss_batch(DB,Num_gam_e,Num_nu,Num_chi,n_threads,gam_e,V_seed,Seed_syn_batch,dot_gam_e_batch)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,Num_chi,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn_batch(Num_nu,Num_chi)
real(8), intent(out) :: dot_gam_e_batch(Num_gam_e,Num_chi)
integer :: I_chi,I_gam_e
real(8) :: Cyclotron_nu,ssa_sum_cell
integer :: batch_work

    call ensure_ssa_seed_cache(Num_nu,V_seed)
    call ensure_ssa_geometry_workspace(Num_gam_e,Num_chi)
    call build_ssa_geometry(DB,Num_gam_e,Num_nu,gam_e,ssa_low_idx_cache,ssa_low_first_cache,ssa_low_last_cache, &
                            ssa_high_first_cache, &
                            ssa_prefactor_low_cache,ssa_prefactor_high_cache,Cyclotron_nu)

    dot_gam_e_batch=zero
    batch_work = Num_gam_e*Num_chi
    if (n_threads <= 1 .or. batch_work < 512) then
       do I_chi=1,Num_chi
          do I_gam_e=1,Num_gam_e
             call accumulate_ssa_batch_cell(Num_gam_e,Num_nu,Num_chi,I_gam_e,I_chi,gam_e,Seed_syn_batch,Cyclotron_nu,ssa_sum_cell)
             dot_gam_e_batch(I_gam_e,I_chi)=ssa_sum_cell
          end do
       end do
    else
       !$OMP PARALLEL DO num_threads(n_threads) collapse(2) schedule(static) &
       !$OMP& private(I_chi,I_gam_e,ssa_sum_cell)
       do I_chi=1,Num_chi
          do I_gam_e=1,Num_gam_e
             call accumulate_ssa_batch_cell(Num_gam_e,Num_nu,Num_chi,I_gam_e,I_chi,gam_e,Seed_syn_batch,Cyclotron_nu,ssa_sum_cell)
             dot_gam_e_batch(I_gam_e,I_chi)=ssa_sum_cell
          end do
       end do
       !$OMP END PARALLEL DO
    end if
end subroutine electron_cooling_ssa_loss_batch

! 数值计算逆康普顿（IC）冷却率：双重积分（种子光子×散射截面），含Jones/Blumenthal核。
subroutine electron_cooling_ic_loss(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn, dot_gam_e)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dot_gam_e(Num_gam_e)

real(8),allocatable,dimension (:) :: photon_number
integer :: work_items

    call ensure_ic_grid_cache(Num_gam_e,Num_nu,gam_e,V_seed)
    allocate (photon_number(Num_nu-1))
    
    dot_gam_e=zero

    do I_nu=1,Num_nu-1
       photon_number(I_nu)=electron_powerlaw_interp(V_seed(I_nu),V_seed(I_nu+1), &
                                                    Seed_syn(I_nu),Seed_syn(I_nu+1),ic_v_seed_mid_cache(I_nu))
    end do

    work_items=(Num_gam_e-1)*(Num_nu-1)*(Num_nu-1)
    if (n_threads <= 1 .or. work_items < 8192) then
       do i_gam_e=1,Num_gam_e-1
          call accumulate_ic_gamma_loss(i_gam_e,dot_gam_e(i_gam_e))
       end do
    else
       !$OMP PARALLEL num_threads(n_threads), &
       !$OMP& private(i_gam_e)
       !$OMP DO SCHEDULE(STATIC)
       do i_gam_e=1,Num_gam_e-1
          call accumulate_ic_gamma_loss(i_gam_e,dot_gam_e(i_gam_e))
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    end if

    dot_gam_e=dot_gam_e/gam_e/gam_e*para_h*Para_h*Para_SigmaT/para_m_energy
    dot_gam_e(Num_gam_e)=0.99*dot_gam_e(Num_gam_e-1)

deallocate (photon_number)

contains

subroutine accumulate_ic_gamma_loss(i_gam_e,dot_val)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: i_gam_e
real(8), intent(out) :: dot_val
integer :: I_nu,Nu_s
real(8) :: game,game_pow,var,dInteg2,V_t,E_t2eV,Vloc,E_Vloc2eV,uplim,temp,q,fssc,kn_factor

    dot_val=zero
    game=ic_gam_e_mean_cache(i_gam_e)
    game_pow=game*game
    var=0.25d0/game_pow
    do I_nu=1,Num_nu-1
       dInteg2=zero
       V_t=ic_v_seed_mid_cache(I_nu)
       E_t2eV=ic_e_seed_cache(I_nu)
       kn_factor=4d0*game*E_t2eV
       uplim=(4d0*game_pow*E_t2eV)/(one+kn_factor)
       do Nu_s=1,Num_nu-1
          fssc=zero
          Vloc=ic_v_seed_mid_cache(Nu_s)
          E_Vloc2eV=ic_e_seed_cache(Nu_s)
          if (Vloc > var*V_t .and. Vloc <= V_t) then
             fssc=Vloc/V_t-var
          else
             if (E_Vloc2eV > uplim) exit
             temp=game-E_Vloc2eV
             if (temp <= zero) exit
             q=E_Vloc2eV/(kn_factor*temp)
             if (q <= zero) cycle
             if (q >= one) exit
             fssc=two*q*(log(q)-q)+one+q+ &
                  0.5d0*(one-q)*(4d0*game*E_t2eV*q)**2/(1+4d0*game*q*E_t2eV)
          end if
          dInteg2=dInteg2+Vloc*fssc*ic_d_nu_cache(Nu_s)
       end do
       dot_val=dot_val+photon_number(I_nu)/V_t*ic_d_nu_cache(I_nu)*dInteg2
    end do
end subroutine accumulate_ic_gamma_loss
end subroutine electron_cooling_ic_loss

! IC cooling from the same Jones/KN emissivity kernel used by radiation_ssc_spectrum.
subroutine electron_cooling_ic_loss_emissivity_budget(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn,dot_gam_over_gam)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dot_gam_over_gam(Num_gam_e)
real(8), allocatable :: seed_weights(:),obs_weights(:),e_seed(:),inv_v_seed(:)
integer :: i_gam,work_items

    allocate(seed_weights(Num_nu),obs_weights(Num_nu),e_seed(Num_nu),inv_v_seed(Num_nu))
    call compute_simpson_weights(seed_weights,Num_nu)
    call compute_simpson_weights(obs_weights,Num_nu)
    h_nu=dlog(V_seed(2))-dlog(V_seed(1))
    h_nu_third=h_nu/3.0d0
    para_hEme=Para_h/para_m_energy
    e_seed=V_seed*para_hEme
    inv_v_seed=one/V_seed
    temp_norm_ic=0.75d0*Para_c*Para_h*Para_SigmaT/Para_m_energy
    dot_gam_over_gam=zero

    work_items=Num_gam_e*Num_nu*Num_nu
    if (n_threads <= 1 .or. work_items < 8192) then
        do i_gam=1,Num_gam_e
            call accumulate_budget_gamma(i_gam,dot_gam_over_gam(i_gam))
        end do
    else
        !$OMP PARALLEL num_threads(n_threads), private(i_gam)
        !$OMP DO SCHEDULE(STATIC)
        do i_gam=1,Num_gam_e
            call accumulate_budget_gamma(i_gam,dot_gam_over_gam(i_gam))
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end if

    deallocate(seed_weights,obs_weights,e_seed,inv_v_seed)

contains

subroutine accumulate_budget_gamma(i_gam,loss_val)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: i_gam
real(8), intent(out) :: loss_val
integer :: i_obs,i_seed
real(8) :: gam,gam2,seed_sum,power_log

    gam=gam_e(i_gam)
    gam2=gam*gam
    loss_val=zero
    do i_obs=1,Num_nu
        if (gam_e(Num_gam_e) <= e_seed(i_obs)) cycle
        seed_sum=zero
        do i_seed=1,Num_nu
            if (Seed_syn(i_seed) <= zero) cycle
            if (i_seed < i_obs) then
                seed_sum=seed_sum+seed_weights(i_seed)*Seed_syn(i_seed)*low_seed_kernel(gam,i_obs,i_seed)/gam2
            else
                seed_sum=seed_sum+seed_weights(i_seed)*Seed_syn(i_seed)*high_seed_kernel(gam,i_obs,i_seed)/gam2
            end if
        end do
        power_log=temp_norm_ic*V_seed(i_obs)*h_nu_third*seed_sum
        loss_val=loss_val+obs_weights(i_obs)*power_log
    end do
    loss_val=h_nu_third*loss_val/gam
end subroutine accumulate_budget_gamma

real(8) function low_seed_kernel(gam,i_obs,i_seed)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: gam
integer, intent(in) :: i_obs,i_seed
real(8) :: temp,q,log_q,q_gamma,kn_coeff

    low_seed_kernel=zero
    temp=gam-e_seed(i_obs)
    if (temp <= zero) return
    q=V_seed(i_obs)/(4.0d0*gam*temp*V_seed(i_seed))
    if (q <= zero .or. q >= one) return
    log_q=dlog(q)
    q_gamma=e_seed(i_obs)/temp
    kn_coeff=q_gamma*q_gamma/(two*(one+q_gamma))
    low_seed_kernel=two*q*(log_q-q)+one+q+kn_coeff*(one-q)
end function low_seed_kernel

real(8) function high_seed_kernel(gam,i_obs,i_seed)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: gam
integer, intent(in) :: i_obs,i_seed

    high_seed_kernel=max(zero,V_seed(i_obs)*inv_v_seed(i_seed)-0.25d0/(gam*gam))
end function high_seed_kernel
end subroutine electron_cooling_ic_loss_emissivity_budget

! 确保Nakar Y参数工作数组已分配（缓存hat_nu、频率段Gauss节点和查找区间）。
subroutine ensure_y_nakar_workspace(Num_gam_e,Num_nu,gam_e,V_seed)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu)
logical :: rebuild

    rebuild=.not. allocated(y_nakar_hat_nu_cache)
    if (.not. rebuild) rebuild = (y_nakar_num_gam_cache /= Num_gam_e)
    if (.not. rebuild) rebuild = (y_nakar_num_nu_cache /= Num_nu)
    if (.not. rebuild) rebuild = any(y_nakar_hat_nu_cache /= Para_m_energy/Para_h/gam_e)
    if (.not. rebuild) rebuild = any(y_nakar_v_cache /= V_seed)
    if (.not. rebuild) return

    if (allocated(y_nakar_hat_nu_cache)) deallocate(y_nakar_hat_nu_cache,y_nakar_prefix_cache,y_nakar_v_cache, &
                                                    y_nakar_vloc_cache,y_nakar_idx_cache,y_nakar_vg1,y_nakar_vg2, &
                                                    y_nakar_wg1,y_nakar_wg2)
    allocate(y_nakar_hat_nu_cache(Num_gam_e),y_nakar_prefix_cache(Num_nu),y_nakar_v_cache(Num_nu), &
             y_nakar_vloc_cache(Num_gam_e),y_nakar_idx_cache(Num_gam_e),y_nakar_vg1(Num_nu-1), &
             y_nakar_vg2(Num_nu-1),y_nakar_wg1(Num_nu-1),y_nakar_wg2(Num_nu-1))
    y_nakar_hat_nu_cache=Para_m_energy/Para_h/gam_e
    y_nakar_v_cache=V_seed
    do I_nu=1,Num_nu-1
        call electron_log_gauss2_interval(V_seed(I_nu),V_seed(I_nu+1), &
                                          y_nakar_vg1(I_nu),y_nakar_vg2(I_nu), &
                                          y_nakar_wg1(I_nu),y_nakar_wg2(I_nu))
    end do
    do I_Compton=1,Num_gam_e
        if (y_nakar_hat_nu_cache(I_Compton) <= V_seed(1)) then
            y_nakar_idx_cache(I_Compton)=0
            y_nakar_vloc_cache(I_Compton)=V_seed(1)
        else
            call first_greater_monotonic_window(V_seed,Num_nu,2,y_nakar_hat_nu_cache(I_Compton),I_nu)
            y_nakar_idx_cache(I_Compton)=I_nu
            if (I_nu <= Num_nu) y_nakar_vloc_cache(I_Compton)=min(y_nakar_hat_nu_cache(I_Compton),V_seed(I_nu))
        end if
    end do
    y_nakar_num_gam_cache=Num_gam_e
    y_nakar_num_nu_cache=Num_nu
end subroutine ensure_y_nakar_workspace

! Nakar+2009 Compton Y参数：Y(γ) = ∫_{ν̂(γ)}^{ν_max} P_syn(ν) dν，谱形依赖。
subroutine electron_cooling_y_nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn, Compton)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu)
real(8), intent(out) :: Compton(Num_gam_e)
    integer, parameter :: parallel_work_threshold=512
    integer :: I_Compton,I_nu,work_items

    call ensure_y_nakar_workspace(Num_gam_e,Num_nu,gam_e,V_seed)

    Compton=zero
    var_Compensation=zero
    y_nakar_prefix_cache(1)=zero
    do I_nu=2,Num_nu
       y_nakar_prefix_cache(I_nu)=y_nakar_prefix_cache(I_nu-1)+ &
            y_nakar_wg1(I_nu-1)*electron_powerlaw_interp(V_seed(I_nu-1),V_seed(I_nu), &
                                                          P_syn(I_nu-1),P_syn(I_nu),y_nakar_vg1(I_nu-1))+ &
            y_nakar_wg2(I_nu-1)*electron_powerlaw_interp(V_seed(I_nu-1),V_seed(I_nu), &
                                                          P_syn(I_nu-1),P_syn(I_nu),y_nakar_vg2(I_nu-1))
    end do

    work_items=Num_gam_e*Num_nu
    if (n_threads <= 1 .or. work_items < parallel_work_threshold) then
       do I_Compton=1,Num_gam_e
          call accumulate_y_nakar_point(I_Compton,Compton(I_Compton))
       end do
    else
       !$OMP PARALLEL num_threads(n_threads), private(I_Compton)
       !$OMP DO SCHEDULE(STATIC)
       do I_Compton=1,Num_gam_e
          call accumulate_y_nakar_point(I_Compton,Compton(I_Compton))
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    end if

contains

subroutine accumulate_y_nakar_point(I_Compton,Compton_val)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: I_Compton
real(8), intent(out) :: Compton_val
integer :: I_nu
real(8) :: V_loc

    Compton_val=zero
    I_nu=y_nakar_idx_cache(I_Compton)
    if (I_nu == 0) return
    if (I_nu <= Num_nu) then
       V_loc=y_nakar_vloc_cache(I_Compton)
       Compton_val=y_nakar_prefix_cache(I_nu-1)+ &
                   electron_integrate_powerlaw_segment(V_seed(I_nu-1),V_loc, &
                       P_syn(I_nu-1), &
                       electron_powerlaw_interp(V_seed(I_nu-1),V_seed(I_nu), &
                                                P_syn(I_nu-1),P_syn(I_nu),V_loc))
    else
       Compton_val=y_nakar_prefix_cache(Num_nu)
    end if
end subroutine accumulate_y_nakar_point
end subroutine electron_cooling_y_nakar

! Fan+2008 Compton Y参数：解析分段η_NK(γ) × η_rad，含快/慢冷却和谱指数分支。
subroutine electron_cooling_y_fan(Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e, Compton)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,gam_e(Num_gam_e)
real(8), intent(out) :: Compton(Num_gam_e)
integer :: i_gam_e

    eta=(Gam_e_m/Gam_e_c)**(p-two)
    if (eta-one > 0.001) eta=one

    do i_gam_e=1,Num_gam_e
        if (Num_gam_e > i_gam_e) then
            hat_gam=5.4246D6/sqrt(DB*gam_e(i_gam_e+1))
            if (Gam_e_m > Gam_e_c) then
                if (hat_gam < Gam_e_c) then
                    eta_NK=zero
                else if (hat_gam < Gam_e_m) then
                    if (p>2) then
                        Step1=(p-1)/(p-2)*Gam_e_m-Gam_e_c
                        eta_NK=(hat_gam-Gam_e_c)/Step1
                    else
                        Step1=Gam_e_m**(p-1)*Gam_e_max**(2-p)-(p-1)*Gam_e_m-(2-p)*Gam_e_c
                        eta_NK=(2-p)*(hat_gam-Gam_e_c)/Step1
                    end if
                else
                    if (p>2) then
                        Step2=Gam_e_m**(p-1)*hat_gam**(2-p)
                        Step3=(p-1)*Gam_e_m-(p-2)*Gam_e_c
                        eta_NK=1-Step2/Step3
                    else
                        Step2=Gam_e_m**(p-1)*Gam_e_max**(2-p)-(p-1)*Gam_e_m-(2-p)*Gam_e_c
                        Step3=Gam_e_m**(p-1)*(Gam_e_max**(2-p)-hat_gam**(2-p))
                        eta_NK=1-Step2/Step3
                    end if
                end if
            else if (hat_gam < Gam_e_m) then
                eta_NK=zero
            else if (hat_gam < Gam_e_c) then
                if (p>2) then
                    Step4=Gam_e_c**(3-p)/(p-2.0)-Gam_e_m**(3-p)
                    eta_NK=(hat_gam**(3-p)-Gam_e_m**(3-p))/Step4
                else
                    Step4=(3-p)*Gam_e_c*Gam_e_max**(2-p)-Gam_e_c**(3-p)-(2-p)*Gam_e_m**(3-p)
                    eta_NK=(2-p)*(hat_gam**(3-p)-Gam_e_m**(3-p))/Step4
                end if
            else
                if (p>2) then
                    Step5=(3-p)*Gam_e_c*hat_gam**(2-p)
                    Step6=Gam_e_c**(3.0-p)-(p-2)*Gam_e_m**(3.0-p)
                    eta_NK=1-Step5/Step6
                else
                    Step5=(3-p)*Gam_e_c*(Gam_e_max**(2-p)-hat_gam**(2-p))
                    Step6=(3-p)*Gam_e_c*Gam_e_max**(2-p)-Gam_e_c**(3-p)-(2-p)*Gam_e_m**(3-p)
                    eta_NK=1-Step5/Step6
                end if
            end if
            Compton(i_gam_e)=0.5d0*(-1.0+sqrt(1.0+4.0*eta*eta_NK*Epsilon_e/Epsilon_b))
        else
            Compton(i_gam_e)=0.99*Compton(i_gam_e-1)
        end if
    end do
end subroutine electron_cooling_y_fan

! 根据index_Y准备正激波冷却辅助量：IC数值积分或Nakar Y参数。
subroutine prepare_forward_cooling_aux(index_Y,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,Seed_syn,cooling_aux)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: cooling_aux(Num_gam_e)

    cooling_aux=zero
    select case(index_Y)
    case(1)
        call electron_cooling_ic_loss(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn,cooling_aux)
    case(2)
        call electron_cooling_y_nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,cooling_aux)
    case default
    end select
end subroutine prepare_forward_cooling_aux

! 批量版prepare_forward_cooling_aux：对多个χ列分别计算冷却辅助量。
subroutine prepare_forward_cooling_aux_batch(index_Y,Num_gam_e,Num_nu,Num_chi,n_threads,gam_e,V_seed,P_syn,Seed_syn,cooling_aux)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,Num_chi,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu,Num_chi),Seed_syn(Num_nu,Num_chi)
real(8), intent(out) :: cooling_aux(Num_gam_e,Num_chi)
integer :: I_chi

    cooling_aux=zero
    select case(index_Y)
    case(1)
        do I_chi=1,Num_chi
            call electron_cooling_ic_loss(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn(:,I_chi),cooling_aux(:,I_chi))
        end do
    case(2)
        do I_chi=1,Num_chi
            call electron_cooling_y_nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn(:,I_chi),cooling_aux(:,I_chi))
        end do
    case default
    end select
end subroutine prepare_forward_cooling_aux_batch

! 组装正激波冷却率 dγ/dR：SSA + 同步 + IC/Compton Y，分离SSA和IC计算。
subroutine assemble_forward_cooling_split(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                          beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed, &
                                          Seed_syn_ssa,cooling_aux,dEl)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
real(8), intent(inout) :: Gam_e_max
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn_ssa(Num_nu),cooling_aux(Num_gam_e)
real(8), intent(out) :: dEl(Num_gam_e)
real(8) :: Compton(Num_gam_e),dot_gam_e_SSA(Num_gam_e)

    call electron_cooling_ssa_loss(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn_ssa,dot_gam_e_SSA)
    call assemble_forward_cooling_from_terms(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                             beta_Gam,dNe,Num_gam_e,gam_e,Compton,dot_gam_e_SSA,cooling_aux,dEl)
end subroutine assemble_forward_cooling_split

! 批量版assemble_forward_cooling_split：对多个χ列分别组装冷却率。
subroutine assemble_forward_cooling_split_batch(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                                beta_Gam,dNe,Num_gam_e,Num_nu,Num_chi,n_threads,gam_e,V_seed, &
                                                Seed_syn_ssa,cooling_aux,dEl)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,Num_chi,n_threads
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
real(8), intent(inout) :: Gam_e_max
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn_ssa(Num_nu,Num_chi),cooling_aux(Num_gam_e,Num_chi)
real(8), intent(out) :: dEl(Num_gam_e,Num_chi)
real(8) :: Compton(Num_gam_e),Gam_e_max_cell,dot_gam_e_SSA(Num_gam_e,Num_chi)
integer :: I_chi

    call electron_cooling_ssa_loss_batch(DB,Num_gam_e,Num_nu,Num_chi,n_threads,gam_e,V_seed,Seed_syn_ssa,dot_gam_e_SSA)
    do I_chi=1,Num_chi
       Gam_e_max_cell=Gam_e_max
       call assemble_forward_cooling_from_terms(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max_cell,R_loc,R_Gamma_loc, &
                                                beta_Gam,dNe,Num_gam_e,gam_e,Compton,dot_gam_e_SSA(:,I_chi), &
                                                cooling_aux(:,I_chi),dEl(:,I_chi))
    end do
end subroutine assemble_forward_cooling_split_batch

! 由各项组装正激波冷却率：dγ/dt = (f_r*Y - SSA)*γ，支持4种Compton Y方案。
subroutine assemble_forward_cooling_from_terms(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                               beta_Gam,dNe,Num_gam_e,gam_e,Compton,dot_gam_e_SSA,cooling_aux,dEl)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
real(8), intent(inout) :: Gam_e_max
real(8), intent(in) :: gam_e(Num_gam_e)
real(8), intent(inout) :: Compton(Num_gam_e)
real(8), intent(in) :: dot_gam_e_SSA(Num_gam_e),cooling_aux(Num_gam_e)
real(8), intent(out) :: dEl(Num_gam_e)

    cooling_scale=one/(beta_Gam*R_Gamma_loc)
    ssa_scale=cooling_scale/para_c
    f_r=1.35d-19*DB**2*cooling_scale/pi

    select case(index_Y)
    case(0)
        dEl=(f_r-dot_gam_e_SSA*ssa_scale)*gam_e
    case(1)
        dEl=(f_r+(cooling_aux-dot_gam_e_SSA)*ssa_scale)*gam_e
    case(2)
        Q=4d0*pi*R_loc*R_loc*para_c
        Compton=one+cooling_aux/Q/(4d0*R_Gamma_loc*R_Gamma_loc*dNe*Para_m_p_E)
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA*ssa_scale)*gam_e
    case(3)
        call electron_cooling_y_fan(Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,Compton)
        Compton=one+Compton
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA*ssa_scale)*gam_e
    case default
        print*, 'invalid Compton case, check your chosen model!'
        stop
    end select
end subroutine assemble_forward_cooling_from_terms

! 正激波冷却主入口：准备IC辅助量→计算SSA→组装冷却率 dγ/dt。
subroutine get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                               beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,Seed_syn,dEl)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
real(8), intent(inout) :: Gam_e_max
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dEl(Num_gam_e)
real(8) :: cooling_aux(Num_gam_e)
real(8) :: Compton(Num_gam_e),dot_gam_e_SSA(Num_gam_e)

    call prepare_forward_cooling_aux(index_Y,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,Seed_syn,cooling_aux)
    call electron_cooling_ssa_loss(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn,dot_gam_e_SSA)
    call assemble_forward_cooling_from_terms(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                                             beta_Gam,dNe,Num_gam_e,gam_e,Compton,dot_gam_e_SSA,cooling_aux,dEl)
end subroutine get_forward_cooling

end module electron_cooling_kernel
