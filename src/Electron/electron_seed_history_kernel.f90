!f2py: skip
module electron_seed_history_kernel
  use constants
  implicit none
  private

  ! 历史种子光子场：按下游固有时把旧壳层辐射输运到当前 chi 网格。
  ! Seed-photon history field: transport older shell emission to the current chi grid in downstream proper time.
  public :: integrate_proper_time
  public :: accumulate_history_fields
  public :: advance_history_stream
  public :: history_transfer_weight

  integer, save :: cache_shells=0, cache_chi=0, cache_nu=0, built_shells=0
  real(8), allocatable, save, dimension(:,:) :: inv_dx
  real(8), allocatable, save, dimension(:,:,:) :: tau_prefix
  logical, allocatable, save, dimension(:) :: map_valid
  integer, allocatable, save, dimension(:) :: map_idx
  real(8), allocatable, save, dimension(:) :: map_frac
  integer, save :: seed_nu=0
  logical, save :: seed_ready=.false.
  real(8), allocatable, save, dimension(:) :: seed_v,seed_logv,seed_invdlog
  !$omp threadprivate(cache_shells,cache_chi,cache_nu,built_shells)
  !$omp threadprivate(inv_dx,tau_prefix,map_valid,map_idx,map_frac)
  !$omp threadprivate(seed_nu,seed_ready,seed_v,seed_logv)
  !$omp threadprivate(seed_invdlog)

contains

! 确保历史场工作数组、频率映射和种子频率缓存已按当前网格分配。
! Ensure history work arrays, frequency maps, and seed-frequency cache match the current grid.
subroutine ensure_cache(Num_shell,Num_chi,Num_nu,V_seed)
implicit none
integer, intent(in) :: Num_shell,Num_chi,Num_nu
real(8), intent(in), dimension(Num_nu) :: V_seed
logical :: rebuild_main, rebuild_map, cache_match

    rebuild_main=.not. allocated(inv_dx)
    if (.not. rebuild_main) then
        rebuild_main = (cache_shells /= Num_shell) .or. (cache_chi /= Num_chi) .or. (cache_nu /= Num_nu)
    end if
    if (rebuild_main) then
        if (allocated(inv_dx)) deallocate(inv_dx)
        if (allocated(tau_prefix)) deallocate(tau_prefix)
        allocate(inv_dx(Num_chi,Num_shell),tau_prefix(Num_nu,0:Num_chi,Num_shell))
        cache_shells=Num_shell
        cache_chi=Num_chi
        cache_nu=Num_nu
        built_shells=0
    end if

    rebuild_map=.not. allocated(map_valid)
    if (.not. rebuild_map) rebuild_map = (size(map_valid) /= Num_nu)
    if (rebuild_map) then
        if (allocated(map_valid)) deallocate(map_valid,map_idx,map_frac)
        allocate(map_valid(Num_nu),map_idx(Num_nu),map_frac(Num_nu))
    end if

    cache_match = .false.
    if (seed_ready) then
        if (seed_nu == Num_nu) then
            if (allocated(seed_v)) cache_match=all(seed_v == V_seed)
        end if
    end if
    if (cache_match) return

    if (allocated(seed_v)) deallocate(seed_v,seed_logv,seed_invdlog)
    allocate(seed_v(Num_nu),seed_logv(Num_nu),seed_invdlog(Num_nu-1))
    seed_v=V_seed
    seed_logv=dlog(V_seed)
    seed_invdlog=1d0/(seed_logv(2:Num_nu)-seed_logv(1:Num_nu-1))
    seed_nu=Num_nu
    seed_ready=.true.
end subroutine ensure_cache

! 沿半径积分下游共动固有时间。
! Integrate downstream comoving proper time along the radius grid.
subroutine integrate_proper_time(Num_shell,r_cm,gamma_bulk,tprop)
implicit none
integer, intent(in) :: Num_shell
integer :: ishell
real(8), intent(in), dimension(Num_shell) :: r_cm,gamma_bulk
real(8), intent(out), dimension(Num_shell) :: tprop
real(8) :: gamma_mean,beta_mean,dR

    tprop=0d0
    do ishell=2,Num_shell
        gamma_mean=0.5d0*(gamma_bulk(ishell-1)+gamma_bulk(ishell))
        beta_mean=sqrt(1d0-1d0/gamma_mean**2)
        dR=r_cm(ishell)-r_cm(ishell-1)
        tprop(ishell)=tprop(ishell-1)+dR/(beta_mean*gamma_mean*Para_c)
    end do
end subroutine integrate_proper_time

! 把可由光行时连接的历史壳层辐射叠加到当前共动光子场。
! Accumulate light-travel-connected older shell emission into the current comoving photon field.
subroutine accumulate_history_fields(target_t,Num_shell,Num_chi,Num_nu,tprop,V_seed, &
                                     xface,xmid,dxcell,beta,tau,pemit,seed, &
                                     peff,seeff,n_threads)
implicit none
integer, intent(in) :: target_t,Num_shell,Num_chi,Num_nu,n_threads
integer :: isrct,isrc,itgt
real(8), intent(in), dimension(Num_shell) :: tprop
real(8), intent(in), dimension(Num_nu) :: V_seed
real(8), intent(in), dimension(0:Num_chi,Num_shell) :: xface
real(8), intent(in), dimension(Num_chi,Num_shell) :: xmid,dxcell
real(8), intent(in), dimension(Num_chi,Num_shell) :: beta
real(8), intent(in), dimension(Num_nu,Num_chi,Num_shell) :: tau,pemit,seed
real(8), intent(out), dimension(Num_nu,Num_chi) :: peff,seeff
real(8) :: dtot,dtsrc

    peff = pemit(:,:,target_t)
    seeff = seed(:,:,target_t)
    !$OMP PARALLEL DO num_threads(n_threads) if(n_threads > 1 .and. target_t*Num_chi*Num_chi*Num_nu >= 512) &
    !$OMP& schedule(static) private(itgt,isrct,isrc,dtot,dtsrc)
    do itgt = 1, Num_chi
        call ensure_cache(Num_shell,Num_chi,Num_nu,V_seed)
        call build_transfer_cache(target_t,Num_chi,Num_nu,dxcell,tau, &
                                  inv_dx,tau_prefix)
        if (target_t > 1) then
            isrct = 1
            dtot = tprop(target_t)-tprop(isrct)
            if (dtot > 0d0) then
                dtsrc = tprop(2)-tprop(1)
                do isrc = 1, Num_chi
                    call add_source_cell(isrct,isrc,itgt,dtot,dtsrc)
                end do
            end if
        end if
        do isrct = 2, target_t-1
            dtot = tprop(target_t)-tprop(isrct)
            if (dtot <= 0d0) cycle
            dtsrc = tprop(isrct)-tprop(isrct-1)
            do isrc = 1, Num_chi
                call add_source_cell(isrct,isrc,itgt,dtot,dtsrc)
            end do
        end do
    end do
    !$OMP END PARALLEL DO

contains

    ! 源单元给出一个可达时间片和有效发射长度。
    ! The source cell supplies a reachable time slice and effective emitting length.
    subroutine add_source_cell(src_t,src_chi,tgt_chi,dtot,dtsrc)
    implicit none
    integer, intent(in) :: src_t,src_chi,tgt_chi
    real(8), intent(in) :: dtot,dtsrc
    real(8) :: swgt,xsrc

        xsrc = xmid(src_chi,src_t)
        if (dxcell(src_chi,src_t) <= 0d0) error stop 'history source cell requires positive dxcell'
        swgt = min(1d0, Para_c*dtsrc/dxcell(src_chi,src_t))
        call add_target_cell(src_t,src_chi,tgt_chi,dtot,swgt,xsrc)
    end subroutine add_source_cell

    ! 目标单元沿光线路径累乘吸收，并把源频率映射到目标频率。
    ! The target cell multiplies path absorption and maps source frequencies to target frequencies.
    subroutine add_target_cell(src_t,src_chi,tgt_chi,dtot,swgt,xsrc)
    implicit none
    integer, intent(in) :: src_t,src_chi,tgt_chi
    integer :: iseg,inu,ilo
    real(8), intent(in) :: dtot,swgt,xsrc
    real(8), dimension(Num_nu) :: attenuation
    real(8) :: amp_p,amp_seed,doprel,dt_seg,xcurr,xprev,xtgt

        xtgt = xmid(tgt_chi,target_t)
        if (xsrc < xtgt) return
        if (Para_c*dtot < xsrc-xtgt) return
        doprel = relative_doppler(beta(src_chi,src_t),beta(tgt_chi,target_t))
        call build_map(Num_nu,V_seed,doprel,map_valid,map_idx,map_frac)

        attenuation = 1d0
        xprev = xsrc
        do iseg = src_t+1, target_t
            dt_seg = tprop(iseg)-tprop(iseg-1)
            if (dt_seg <= 0d0) cycle
            xcurr = max(xtgt, xprev-Para_c*dt_seg)
            call apply_path_tau(Num_chi,Num_nu,xprev,xcurr,xface(:,iseg), &
                                inv_dx(:,iseg),tau(:,:,iseg), &
                                tau_prefix(:,:,iseg),attenuation)
            xprev = xcurr
            if (xprev <= xtgt) exit
        end do
        if (doprel <= 0d0 .or. swgt <= 0d0) return
        amp_p = swgt*doprel**3
        amp_seed = swgt*doprel**2
        do inu = 1, Num_nu
            if (.not. map_valid(inu)) cycle
            ilo = map_idx(inu)
            peff(inu,tgt_chi) = peff(inu,tgt_chi) + amp_p*attenuation(inu) * &
                log_interp(pemit(ilo,src_chi,src_t), &
                                     pemit(ilo+1,src_chi,src_t),map_frac(inu))
            seeff(inu,tgt_chi) = seeff(inu,tgt_chi) + amp_seed*attenuation(inu) * &
                log_interp(seed(ilo,src_chi,src_t), &
                                     seed(ilo+1,src_chi,src_t),map_frac(inu))
        end do
    end subroutine add_target_cell
end subroutine accumulate_history_fields

! 以一阶特征线递推下游历史光子场，避免每个目标壳层回扫所有过去壳层。
! Advance the downstream photon history with first-order characteristics instead of rescanning all past shells.
subroutine advance_history_stream(prev_t,target_t,Num_shell,Num_chi,Num_nu,tprop,V_seed, &
                                  xface,xmid,dxcell,beta,tau,pemit,seed, &
                                  pstream,sstream)
implicit none
integer, intent(in) :: prev_t,target_t,Num_shell,Num_chi,Num_nu
integer :: isrc,itgt,inu,ilo
real(8), intent(in), dimension(Num_shell) :: tprop
real(8), intent(in), dimension(Num_nu) :: V_seed
real(8), intent(in), dimension(0:Num_chi,Num_shell) :: xface
real(8), intent(in), dimension(Num_chi,Num_shell) :: xmid,dxcell
real(8), intent(in), dimension(Num_chi,Num_shell) :: beta
real(8), intent(in), dimension(Num_nu,Num_chi,Num_shell) :: tau,pemit,seed
real(8), intent(inout), dimension(Num_nu,Num_chi) :: pstream,sstream
real(8), dimension(Num_nu,Num_chi) :: pnext,snext
real(8), dimension(Num_nu) :: attenuation
real(8) :: dtau,path_hi,path_lo,seg_lo,seg_hi,xsrc,xtgt,doprel,weight
real(8) :: amp_p,amp_seed

    call ensure_cache(Num_shell,Num_chi,Num_nu,V_seed)
    call build_transfer_cache(target_t,Num_chi,Num_nu,dxcell,tau, &
                              inv_dx,tau_prefix)
    dtau = tprop(target_t)-tprop(prev_t)
    pnext = 0d0
    snext = 0d0
    do itgt = 1, Num_chi
        xtgt = xmid(itgt,target_t)
        path_lo = xtgt
        path_hi = xtgt + Para_c*dtau
        if (path_hi >= xface(0,prev_t) .and. path_hi <= xface(Num_chi,prev_t)) then
            call locate_path_cell(Num_chi,xface(:,prev_t),path_hi,.false.,isrc)
            call add_mapped_cell(isrc,path_hi,itgt,1d0,pstream,sstream)
        end if
        do isrc = 1, Num_chi
            seg_lo = max(xface(isrc-1,prev_t),path_lo)
            seg_hi = min(xface(isrc,prev_t),path_hi)
            if (seg_hi <= seg_lo) cycle
            if (dxcell(isrc,prev_t) <= 0d0) error stop 'history stream cell requires positive dxcell'
            weight = (seg_hi-seg_lo)/dxcell(isrc,prev_t)
            xsrc = 0.5d0*(seg_lo+seg_hi)
            call add_mapped_cell(isrc,xsrc,itgt,weight,pemit(:,:,prev_t),seed(:,:,prev_t))
        end do
    end do
    pstream = pnext
    sstream = snext

contains

    ! 把上一壳层流式历史场或本壳层发射映射到当前目标单元。
    ! Map either the previous streamed history field or current-shell emission into a single target cell.
    subroutine add_mapped_cell(src_chi,xsrcpos,tgt_chi,swgt,P_src,Seed_src)
    implicit none
    integer, intent(in) :: src_chi,tgt_chi
    real(8), intent(in), dimension(Num_nu,Num_chi) :: P_src,Seed_src
    real(8), intent(in) :: xsrcpos,swgt

        xtgt = xmid(tgt_chi,target_t)
        if (xsrcpos < xtgt) return
        if (swgt <= 0d0) return
        xsrc = xsrcpos
        doprel = relative_doppler(beta(src_chi,prev_t),beta(tgt_chi,target_t))
        if (doprel <= 0d0) return
        call build_map(Num_nu,V_seed,doprel,map_valid,map_idx,map_frac)
        attenuation = 1d0
        call apply_path_tau(Num_chi,Num_nu,xsrc,xtgt,xface(:,target_t), &
                            inv_dx(:,target_t),tau(:,:,target_t), &
                            tau_prefix(:,:,target_t),attenuation)
        amp_p = swgt*doprel**3
        amp_seed = swgt*doprel**2
        do inu = 1, Num_nu
            if (.not. map_valid(inu)) cycle
            ilo = map_idx(inu)
            pnext(inu,tgt_chi) = pnext(inu,tgt_chi) + amp_p*attenuation(inu) * &
                log_interp(P_src(ilo,src_chi),P_src(ilo+1,src_chi),map_frac(inu))
            snext(inu,tgt_chi) = snext(inu,tgt_chi) + amp_seed*attenuation(inu) * &
                log_interp(Seed_src(ilo,src_chi),Seed_src(ilo+1,src_chi),map_frac(inu))
        end do
    end subroutine add_mapped_cell
end subroutine advance_history_stream

! 构造当前相对多普勒因子下目标频率到源频率网格的映射。
! Build the target-to-source frequency map for the current relative Doppler factor.
subroutine build_map(Num_nu,V_seed,doprel,validmap,idxmap,fracmap)
implicit none
integer, intent(in) :: Num_nu
logical, intent(out), dimension(Num_nu) :: validmap
integer, intent(out), dimension(Num_nu) :: idxmap
real(8), intent(in), dimension(Num_nu) :: V_seed
real(8), intent(in) :: doprel
real(8), intent(out), dimension(Num_nu) :: fracmap
integer :: inu,ilo
real(8) :: nusrc

    validmap = .false.
    idxmap = 1
    fracmap = 0d0
    if (doprel <= 0d0) return

    ilo = 1
    do inu = 1, Num_nu
        nusrc = V_seed(inu)/doprel
        if (nusrc < V_seed(1) .or. nusrc > V_seed(Num_nu)) cycle
        do while (ilo < Num_nu-1)
            if (V_seed(ilo+1) > nusrc) exit
            ilo = ilo + 1
        end do
        idxmap(inu) = ilo
        if (seed_invdlog(ilo) <= 0d0) cycle
        fracmap(inu) = (dlog(nusrc)-seed_logv(ilo))*seed_invdlog(ilo)
        validmap(inu) = .true.
    end do
end subroutine build_map

! 按预计算对数分数做 log-log 插值：y = y0 * exp(log_frac * log(y1/y0))。
! Use the precomputed logarithmic fraction for log-log interpolation.
real(8) function log_interp(y0,y1,log_frac)
implicit none
real(8), intent(in) :: y0,y1,log_frac

    if (y0 <= 0d0 .or. y1 <= 0d0) then
        log_interp = 0d0
    else
        log_interp = y0*dexp(log_frac*dlog(y1/y0))
    end if
end function log_interp

! 计算从历史源区到当前目标区的相对多普勒因子：D = γ_rel(1+β_rel)。
! Compute the relative Doppler factor from a history source cell to a current target cell.
real(8) function relative_doppler(bsrc,btgt)
implicit none
real(8), intent(in) :: bsrc,btgt
real(8) :: brel,grel

    brel = (btgt-bsrc)/(1d0-btgt*bsrc)
    if (dabs(brel) >= 1d0) error stop 'relative_doppler requires subluminal relative beta'
    grel = 1d0/dsqrt(1d0-brel*brel)
    relative_doppler = grel*(1d0+brel)
end function relative_doppler

! 缓存每个壳层 chi 单元宽度倒数，供历史光线路径吸收计算复用。
! Cache inverse chi-cell widths and cumulative log transfer for path absorption reuse.
subroutine build_transfer_cache(Num_shell,Num_chi,Num_nu,dxcell,tau,inv_dxcell,logprefix)
implicit none
integer, intent(in) :: Num_shell,Num_chi,Num_nu
integer :: ishell,ichi,inu
real(8), intent(in), dimension(Num_chi,Num_shell) :: dxcell
real(8), intent(in), dimension(Num_nu,Num_chi,Num_shell) :: tau
real(8), intent(out), dimension(Num_chi,Num_shell) :: inv_dxcell
real(8), intent(out), dimension(Num_nu,0:Num_chi,Num_shell) :: logprefix

    if (Num_shell < built_shells) built_shells=0
    do ishell = built_shells+1, Num_shell
        logprefix(:,0,ishell) = 0d0
        do ichi = 1, Num_chi
            if (dxcell(ichi,ishell) <= 0d0) error stop 'history transfer cache requires positive dxcell'
            inv_dxcell(ichi,ishell) = 1d0/dxcell(ichi,ishell)
            do inu = 1, Num_nu
                logprefix(inu,ichi,ishell) = logprefix(inu,ichi-1,ishell) + &
                    dlog(history_transfer_weight(tau(inu,ichi,ishell)))
            end do
        end do
    end do
    built_shells=max(built_shells,Num_shell)
end subroutine build_transfer_cache

! 用壳层内路径段的吸收权重对 attenuation 做累乘。
! Multiply attenuation by the absorption weight of a path segment inside a single shell.
subroutine apply_path_tau(Num_chi,Num_nu,xstart,xstop,xface,invcell,tau_cell, &
                          logprefix,attenuation)
implicit none
integer, intent(in) :: Num_chi,Num_nu
integer :: istart,istop,inu
real(8), intent(in), dimension(0:Num_chi) :: xface
real(8), intent(in), dimension(Num_chi) :: invcell
real(8), intent(in), dimension(Num_nu,Num_chi) :: tau_cell
real(8), intent(in) :: xstart,xstop
real(8), intent(in), dimension(Num_nu,0:Num_chi) :: logprefix
real(8), intent(inout), dimension(Num_nu) :: attenuation
real(8) :: fstart,fstop

    if (xstart <= xstop) return
    call locate_path_cell(Num_chi,xface,xstart,.true.,istart)
    call locate_path_cell(Num_chi,xface,xstop,.false.,istop)
    if (istart == istop) then
        if (xstart > xstop) then
            do inu = 1, Num_nu
                attenuation(inu) = attenuation(inu) * &
                    history_transfer_weight((xstart-xstop)*invcell(istart)*tau_cell(inu,istart))
            end do
        end if
        return
    end if

    fstart = (xstart-xface(istart-1))*invcell(istart)
    fstop = (xface(istop)-xstop)*invcell(istop)
    if (fstart > 0d0) then
        do inu = 1, Num_nu
            attenuation(inu) = attenuation(inu) * history_transfer_weight(fstart*tau_cell(inu,istart))
        end do
    end if
    if (istop < istart-1) then
        do inu = 1, Num_nu
            attenuation(inu) = attenuation(inu) * &
                dexp(logprefix(inu,istart-1)-logprefix(inu,istop))
        end do
    end if
    if (fstop > 0d0) then
        do inu = 1, Num_nu
            attenuation(inu) = attenuation(inu) * history_transfer_weight(fstop*tau_cell(inu,istop))
        end do
    end if
end subroutine apply_path_tau

! 在下游面网格中定位路径端点所在单元。
! Locate the downstream face-grid cell containing a path endpoint.
subroutine locate_path_cell(Num_chi,xface,xpos,use_upper,icell)
implicit none
integer, intent(in) :: Num_chi
integer, intent(out) :: icell
integer :: left,right,mid
real(8), intent(in), dimension(0:Num_chi) :: xface
real(8), intent(in) :: xpos
logical, intent(in) :: use_upper

    if (xpos >= xface(Num_chi)) then
        icell = Num_chi
        return
    end if
    if (use_upper) then
        if (xpos <= xface(1)) then
            icell = 1
            return
        end if
    else
        if (xpos <= xface(0)) then
            icell = 1
            return
        end if
    end if

    left = 1
    right = Num_chi
    do while (left < right)
        if (use_upper) then
            mid = (left+right+1)/2
            if (xface(mid) < xpos) then
                left = mid
            else
                right = mid-1
            end if
        else
            mid = (left+right)/2
            if (xface(mid) > xpos) then
                right = mid
            else
                left = mid + 1
            end if
        end if
    end do
    if (use_upper) then
        icell = left + 1
    else
        icell = left
    end if
end subroutine locate_path_cell

! 把光深转换为均匀源函数逃逸/传输权重。
! Convert optical depth into the escape/transfer weight for a uniform source function.
elemental real(8) function history_transfer_weight(tauseg)
implicit none
real(8), intent(in) :: tauseg
real(8) :: tauloc

    tauloc=max(0d0,tauseg)
    if (tauloc < 1d-10) then
        history_transfer_weight=1d0-0.5d0*tauloc+tauloc*tauloc/6d0
    else
        history_transfer_weight=(1d0-dexp(-tauloc))/tauloc
    end if
end function history_transfer_weight

end module electron_seed_history_kernel
