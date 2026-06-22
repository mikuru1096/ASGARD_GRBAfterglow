!f2py: skip
module electron_seed_history_kernel
  use constants
  implicit none
  private

  public :: integrate_downstream_proper_time
  public :: accumulate_comoving_history_fields
  public :: history_transfer_weight

  integer, save :: hist_cache_num_shell=0, hist_cache_num_chi=0, hist_cache_num_nu=0, hist_cache_built_shells=0
  real(8), allocatable, save :: hist_inv_dx_ws(:,:), hist_log_transfer_prefix_ws(:,:,:)
  logical, allocatable, save :: hist_valid_map_ws(:)
  integer, allocatable, save :: hist_idx_lo_map_ws(:)
  real(8), allocatable, save :: hist_log_frac_map_ws(:)
  integer, save :: hist_seed_num_nu_cache=0
  logical, save :: hist_seed_cache_ready=.false.
  real(8), allocatable, save :: hist_v_seed_cache(:), hist_log_v_seed_cache(:), hist_inv_dlog_v_cell_cache(:)
  !$omp threadprivate(hist_cache_num_shell,hist_cache_num_chi,hist_cache_num_nu,hist_cache_built_shells)
  !$omp threadprivate(hist_inv_dx_ws,hist_log_transfer_prefix_ws,hist_valid_map_ws,hist_idx_lo_map_ws,hist_log_frac_map_ws)
  !$omp threadprivate(hist_seed_num_nu_cache,hist_seed_cache_ready,hist_v_seed_cache,hist_log_v_seed_cache)
  !$omp threadprivate(hist_inv_dlog_v_cell_cache)

contains

! 确保历史场工作数组、频率映射和种子频率缓存已按当前网格分配。
subroutine ensure_history_cache(Num_shell,Num_chi,Num_nu,V_seed)
implicit none
integer, intent(in) :: Num_shell,Num_chi,Num_nu
real(8), intent(in) :: V_seed(Num_nu)
logical :: rebuild_main, rebuild_map, cache_match

    rebuild_main=.not. allocated(hist_inv_dx_ws)
    if (.not. rebuild_main) then
        rebuild_main = (hist_cache_num_shell /= Num_shell) .or. (hist_cache_num_chi /= Num_chi) .or. (hist_cache_num_nu /= Num_nu)
    end if
    if (rebuild_main) then
        if (allocated(hist_inv_dx_ws)) deallocate(hist_inv_dx_ws)
        if (allocated(hist_log_transfer_prefix_ws)) deallocate(hist_log_transfer_prefix_ws)
        allocate(hist_inv_dx_ws(Num_chi,Num_shell),hist_log_transfer_prefix_ws(Num_nu,0:Num_chi,Num_shell))
        hist_cache_num_shell=Num_shell
        hist_cache_num_chi=Num_chi
        hist_cache_num_nu=Num_nu
        hist_cache_built_shells=0
    end if

    rebuild_map=.not. allocated(hist_valid_map_ws)
    if (.not. rebuild_map) rebuild_map = (size(hist_valid_map_ws) /= Num_nu)
    if (rebuild_map) then
        if (allocated(hist_valid_map_ws)) deallocate(hist_valid_map_ws,hist_idx_lo_map_ws,hist_log_frac_map_ws)
        allocate(hist_valid_map_ws(Num_nu),hist_idx_lo_map_ws(Num_nu),hist_log_frac_map_ws(Num_nu))
    end if

    cache_match = .false.
    if (hist_seed_cache_ready) then
        if (hist_seed_num_nu_cache == Num_nu) then
            if (allocated(hist_v_seed_cache)) cache_match=all(hist_v_seed_cache == V_seed)
        end if
    end if
    if (cache_match) return

    if (allocated(hist_v_seed_cache)) deallocate(hist_v_seed_cache,hist_log_v_seed_cache,hist_inv_dlog_v_cell_cache)
    allocate(hist_v_seed_cache(Num_nu),hist_log_v_seed_cache(Num_nu),hist_inv_dlog_v_cell_cache(Num_nu-1))
    hist_v_seed_cache=V_seed
    hist_log_v_seed_cache=dlog(V_seed)
    hist_inv_dlog_v_cell_cache=one/(hist_log_v_seed_cache(2:Num_nu)-hist_log_v_seed_cache(1:Num_nu-1))
    hist_seed_num_nu_cache=Num_nu
    hist_seed_cache_ready=.true.
end subroutine ensure_history_cache

! 沿半径积分下游共动固有时间。
subroutine integrate_downstream_proper_time(Num_shell,R_cm,Gamma_bulk,proper_time_s)
implicit none
integer, intent(in) :: Num_shell
integer :: I_shell
real(8), intent(in) :: R_cm(Num_shell),Gamma_bulk(Num_shell)
real(8), intent(out) :: proper_time_s(Num_shell)
real(8) :: gamma_mean,beta_mean,dR

    proper_time_s=zero
    do I_shell=2,Num_shell
        gamma_mean=0.5d0*(Gamma_bulk(I_shell-1)+Gamma_bulk(I_shell))
        beta_mean=sqrt(one-one/gamma_mean**2)
        dR=R_cm(I_shell)-R_cm(I_shell-1)
        proper_time_s(I_shell)=proper_time_s(I_shell-1)+dR/(beta_mean*gamma_mean*Para_c)
    end do
end subroutine integrate_downstream_proper_time

! 把可由光行时连接的历史壳层辐射叠加到当前共动光子场。
subroutine accumulate_comoving_history_fields(target_t,Num_shell,Num_chi,Num_nu,proper_time_s,V_seed, &
                                              x_face_hist,x_center_hist,dx_hist,beta_hist,tau_hist,P_hist,Seed_hist, &
                                              P_eff,Seed_eff,n_threads)
implicit none
integer, intent(in) :: target_t,Num_shell,Num_chi,Num_nu,n_threads
integer :: I_src_t,I_src_chi,I_tgt_chi
real(8), intent(in) :: proper_time_s(Num_shell),V_seed(Num_nu)
real(8), intent(in) :: x_face_hist(0:Num_chi,Num_shell),x_center_hist(Num_chi,Num_shell),dx_hist(Num_chi,Num_shell)
real(8), intent(in) :: beta_hist(Num_chi,Num_shell)
real(8), intent(in) :: tau_hist(Num_nu,Num_chi,Num_shell),P_hist(Num_nu,Num_chi,Num_shell),Seed_hist(Num_nu,Num_chi,Num_shell)
real(8), intent(out) :: P_eff(Num_nu,Num_chi),Seed_eff(Num_nu,Num_chi)
real(8) :: delta_tau_total,dtau_src

    P_eff = P_hist(:,:,target_t)
    Seed_eff = Seed_hist(:,:,target_t)
    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. target_t*Num_chi*Num_chi*Num_nu >= 512) &
    !$OMP& schedule(static) private(I_tgt_chi,I_src_t,I_src_chi,delta_tau_total,dtau_src)
    do I_tgt_chi = 1, Num_chi
        call ensure_history_cache(Num_shell,Num_chi,Num_nu,V_seed)
        call build_shell_transfer_cache(target_t,Num_chi,Num_nu,dx_hist,tau_hist, &
                                        hist_inv_dx_ws,hist_log_transfer_prefix_ws)
        do I_src_t = 1, target_t-1
            delta_tau_total = proper_time_s(target_t)-proper_time_s(I_src_t)
            if (delta_tau_total <= zero) cycle
            if (I_src_t == 1) then
                dtau_src = proper_time_s(2)-proper_time_s(1)
            else
                dtau_src = proper_time_s(I_src_t)-proper_time_s(I_src_t-1)
            end if
            do I_src_chi = 1, Num_chi
                call accumulate_history_source_cell(I_src_t,I_src_chi,I_tgt_chi,delta_tau_total,dtau_src)
            end do
        end do
    end do
    !$OMP END PARALLEL DO

contains

    subroutine accumulate_history_source_cell(src_t,src_chi,tgt_chi,delta_tau_total_src,dtau_src)
    implicit none
    integer, intent(in) :: src_t,src_chi,tgt_chi
    real(8), intent(in) :: delta_tau_total_src,dtau_src
    real(8) :: source_weight,x_src

        x_src = x_center_hist(src_chi,src_t)
        source_weight = min(one, Para_c*dtau_src/max(dx_hist(src_chi,src_t), tiny(one)))
        call accumulate_history_target_cell(src_t,src_chi,tgt_chi,delta_tau_total_src,source_weight,x_src)
    end subroutine accumulate_history_source_cell

    subroutine accumulate_history_target_cell(src_t,src_chi,tgt_chi,delta_tau_total_src,source_weight,x_src)
    implicit none
    integer, intent(in) :: src_t,src_chi,tgt_chi
    integer :: I_seg,I_nu,I_lo
    real(8), intent(in) :: delta_tau_total_src,source_weight,x_src
    real(8) :: attenuation(Num_nu),amp_p,amp_seed,doppler_rel,dt_seg,x_curr,x_prev,x_tgt

        x_tgt = x_center_hist(tgt_chi,target_t)
        if (x_src < x_tgt) return
        if (Para_c*delta_tau_total_src < x_src-x_tgt) return
        doppler_rel = relative_doppler_backward(beta_hist(src_chi,src_t),beta_hist(tgt_chi,target_t))
        call build_doppler_map(Num_nu,V_seed,doppler_rel,hist_valid_map_ws,hist_idx_lo_map_ws,hist_log_frac_map_ws)

        attenuation = one
        x_prev = x_src
        do I_seg = src_t+1, target_t
            dt_seg = proper_time_s(I_seg)-proper_time_s(I_seg-1)
            if (dt_seg <= zero) cycle
            x_curr = max(x_tgt, x_prev-Para_c*dt_seg)
            call apply_shell_path_attenuation(Num_chi,Num_nu,x_prev,x_curr,x_face_hist(:,I_seg), &
                                              hist_inv_dx_ws(:,I_seg),tau_hist(:,:,I_seg), &
                                              hist_log_transfer_prefix_ws(:,:,I_seg),attenuation)
            x_prev = x_curr
            if (x_prev <= x_tgt) exit
        end do
        if (doppler_rel <= zero .or. source_weight <= zero) return
        amp_p = source_weight*doppler_rel**3
        amp_seed = source_weight*doppler_rel**2
        do I_nu = 1, Num_nu
            if (.not. hist_valid_map_ws(I_nu)) cycle
            I_lo = hist_idx_lo_map_ws(I_nu)
            P_eff(I_nu,tgt_chi) = P_eff(I_nu,tgt_chi) + amp_p*attenuation(I_nu) * &
                loglog_interp_mapped(P_hist(I_lo,src_chi,src_t), &
                                     P_hist(I_lo+1,src_chi,src_t),hist_log_frac_map_ws(I_nu))
            Seed_eff(I_nu,tgt_chi) = Seed_eff(I_nu,tgt_chi) + amp_seed*attenuation(I_nu) * &
                loglog_interp_mapped(Seed_hist(I_lo,src_chi,src_t), &
                                     Seed_hist(I_lo+1,src_chi,src_t),hist_log_frac_map_ws(I_nu))
        end do
    end subroutine accumulate_history_target_cell
end subroutine accumulate_comoving_history_fields

! 构造当前相对多普勒因子下目标频率到源频率网格的映射。
subroutine build_doppler_map(Num_nu,V_seed,doppler_rel,valid_map,idx_lo_map,log_frac_map)
implicit none
integer, intent(in) :: Num_nu
logical, intent(out) :: valid_map(Num_nu)
integer, intent(out) :: idx_lo_map(Num_nu)
real(8), intent(in) :: V_seed(Num_nu),doppler_rel
real(8), intent(out) :: log_frac_map(Num_nu)
integer :: I_nu,I_lo
real(8) :: nu_src

    valid_map = .false.
    idx_lo_map = 1
    log_frac_map = zero
    if (doppler_rel <= zero) return

    I_lo = 1
    do I_nu = 1, Num_nu
        nu_src = V_seed(I_nu)/doppler_rel
        if (nu_src < V_seed(1) .or. nu_src > V_seed(Num_nu)) cycle
        do while (I_lo < Num_nu-1)
            if (V_seed(I_lo+1) > nu_src) exit
            I_lo = I_lo + 1
        end do
        idx_lo_map(I_nu) = I_lo
        if (hist_inv_dlog_v_cell_cache(I_lo) <= zero) cycle
        log_frac_map(I_nu) = (dlog(nu_src)-hist_log_v_seed_cache(I_lo))*hist_inv_dlog_v_cell_cache(I_lo)
        valid_map(I_nu) = .true.
    end do
end subroutine build_doppler_map

! 按预计算对数分数做 log-log 插值。
! 按预计算对数分数做 log-log 插值：y = y0 * exp(log_frac * log(y1/y0))。
real(8) function loglog_interp_mapped(y0,y1,log_frac)
implicit none
real(8), intent(in) :: y0,y1,log_frac

    if (y0 <= zero .or. y1 <= zero) then
        loglog_interp_mapped = zero
    else
        loglog_interp_mapped = y0*dexp(log_frac*dlog(y1/y0))
    end if
end function loglog_interp_mapped

! 计算从历史源区到当前目标区的相对多普勒因子。
! 计算从历史源区到当前目标区的相对多普勒因子：D = γ_rel(1+β_rel)。
real(8) function relative_doppler_backward(beta_src,beta_tgt)
implicit none
real(8), intent(in) :: beta_src,beta_tgt
real(8) :: beta_rel,gamma_rel

    beta_rel = (beta_tgt-beta_src)/(one-beta_tgt*beta_src)
    gamma_rel = one/dsqrt(max(one-beta_rel*beta_rel, tiny(one)))
    relative_doppler_backward = gamma_rel*(one+beta_rel)
end function relative_doppler_backward

! 缓存每个壳层 chi 单元宽度倒数，供历史光线路径吸收计算复用。
subroutine build_shell_transfer_cache(Num_shell,Num_chi,Num_nu,dx_hist,tau_hist,inv_dx_hist,log_transfer_prefix)
implicit none
integer, intent(in) :: Num_shell,Num_chi,Num_nu
integer :: I_shell,I_chi,I_nu
real(8), intent(in) :: dx_hist(Num_chi,Num_shell)
real(8), intent(in) :: tau_hist(Num_nu,Num_chi,Num_shell)
real(8), intent(out) :: inv_dx_hist(Num_chi,Num_shell)
real(8), intent(out) :: log_transfer_prefix(Num_nu,0:Num_chi,Num_shell)

    if (Num_shell < hist_cache_built_shells) hist_cache_built_shells=0
    do I_shell = hist_cache_built_shells+1, Num_shell
        log_transfer_prefix(:,0,I_shell) = zero
        do I_chi = 1, Num_chi
            inv_dx_hist(I_chi,I_shell) = one/max(dx_hist(I_chi,I_shell), tiny(one))
            do I_nu = 1, Num_nu
                log_transfer_prefix(I_nu,I_chi,I_shell) = log_transfer_prefix(I_nu,I_chi-1,I_shell) + &
                    dlog(history_transfer_weight(tau_hist(I_nu,I_chi,I_shell)))
            end do
        end do
    end do
    hist_cache_built_shells=max(hist_cache_built_shells,Num_shell)
end subroutine build_shell_transfer_cache

! 用壳层内路径段的吸收权重对 attenuation 做累乘。
subroutine apply_shell_path_attenuation(Num_chi,Num_nu,x_start,x_stop,x_face,inv_dx_cell,tau_cell, &
                                        log_transfer_prefix,attenuation)
implicit none
integer, intent(in) :: Num_chi,Num_nu
integer :: I_start,I_stop,I_nu
real(8), intent(in) :: x_start,x_stop,x_face(0:Num_chi),inv_dx_cell(Num_chi),tau_cell(Num_nu,Num_chi)
real(8), intent(in) :: log_transfer_prefix(Num_nu,0:Num_chi)
real(8), intent(inout) :: attenuation(Num_nu)
real(8) :: frac_start,frac_stop

    if (x_start <= x_stop) return
    call locate_path_cell(Num_chi,x_face,x_start,.true.,I_start)
    call locate_path_cell(Num_chi,x_face,x_stop,.false.,I_stop)
    if (I_start == I_stop) then
        if (x_start > x_stop) then
            do I_nu = 1, Num_nu
                attenuation(I_nu) = attenuation(I_nu) * &
                    history_transfer_weight((x_start-x_stop)*inv_dx_cell(I_start)*tau_cell(I_nu,I_start))
            end do
        end if
        return
    end if

    frac_start = (x_start-x_face(I_start-1))*inv_dx_cell(I_start)
    frac_stop = (x_face(I_stop)-x_stop)*inv_dx_cell(I_stop)
    if (frac_start > zero) then
        do I_nu = 1, Num_nu
            attenuation(I_nu) = attenuation(I_nu) * history_transfer_weight(frac_start*tau_cell(I_nu,I_start))
        end do
    end if
    if (I_stop < I_start-1) then
        do I_nu = 1, Num_nu
            attenuation(I_nu) = attenuation(I_nu) * &
                dexp(log_transfer_prefix(I_nu,I_start-1)-log_transfer_prefix(I_nu,I_stop))
        end do
    end if
    if (frac_stop > zero) then
        do I_nu = 1, Num_nu
            attenuation(I_nu) = attenuation(I_nu) * history_transfer_weight(frac_stop*tau_cell(I_nu,I_stop))
        end do
    end if
end subroutine apply_shell_path_attenuation

! 在下游面网格中定位路径端点所在单元。
subroutine locate_path_cell(Num_chi,x_face,x_pos,use_upper_face,I_cell)
implicit none
integer, intent(in) :: Num_chi
integer, intent(out) :: I_cell
integer :: left,right,mid
real(8), intent(in) :: x_face(0:Num_chi),x_pos
logical, intent(in) :: use_upper_face

    if (x_pos >= x_face(Num_chi)) then
        I_cell = Num_chi
        return
    end if
    if (use_upper_face) then
        if (x_pos <= x_face(1)) then
            I_cell = 1
            return
        end if
    else
        if (x_pos <= x_face(0)) then
            I_cell = 1
            return
        end if
    end if

    left = 1
    right = Num_chi
    do while (left < right)
        if (use_upper_face) then
            mid = (left+right+1)/2
            if (x_face(mid) < x_pos) then
                left = mid
            else
                right = mid-1
            end if
        else
            mid = (left+right)/2
            if (x_face(mid) > x_pos) then
                right = mid
            else
                left = mid + 1
            end if
        end if
    end do
    if (use_upper_face) then
        I_cell = left + 1
    else
        I_cell = left
    end if
end subroutine locate_path_cell

! 把光深转换为均匀源函数逃逸/传输权重。
elemental real(8) function history_transfer_weight(tau_segment)
implicit none
real(8), intent(in) :: tau_segment
real(8) :: tau_loc

    tau_loc=max(zero,tau_segment)
    if (tau_loc < 1d-10) then
        history_transfer_weight=one-0.5d0*tau_loc+tau_loc*tau_loc/6d0
    else
        history_transfer_weight=(one-dexp(-tau_loc))/tau_loc
    end if
end function history_transfer_weight

end module electron_seed_history_kernel
