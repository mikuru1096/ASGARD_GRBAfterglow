!f2py: skip
module electron_seed_history_kernel
  use constants
  implicit none
  private

  ! 历史种子光子场：按下游固有时把旧壳层辐射输运到当前 chi 网格。
  ! Seed-photon history field: transport older shell emission to the current chi grid in downstream proper time.
  public :: integrate_proper_time, advance_history_stream, history_transfer_weight

contains

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

! 以一阶特征线递推下游历史光子场，避免每个目标壳层回扫所有过去壳层。
! Advance the downstream photon history with first-order characteristics instead of rescanning all past shells.
subroutine advance_history_stream(prev_t,target_t,Num_shell,Num_chi,Num_nu,tprop,V_seed, &
                                  xface,xmid,dxcell,beta,tau,pemit,seed, &
                                  pstream,sstream,hist_inv,hist_prefix)
implicit none
integer, intent(in) :: prev_t,target_t,Num_shell,Num_chi,Num_nu
integer :: isrc,itgt
real(8), intent(in), dimension(Num_shell) :: tprop
real(8), intent(in), dimension(Num_nu) :: V_seed
real(8), intent(in), dimension(0:Num_chi,Num_shell) :: xface
real(8), intent(in), dimension(Num_chi,Num_shell) :: xmid, dxcell, beta
real(8), intent(in), dimension(Num_nu,Num_chi,Num_shell) :: tau,pemit,seed
real(8), intent(inout), dimension(Num_nu,Num_chi) :: pstream,sstream
real(8), intent(out), dimension(Num_chi) :: hist_inv
real(8), intent(out), dimension(Num_nu,0:Num_chi) :: hist_prefix
real(8), dimension(Num_nu,Num_chi) :: pnext,snext
real(8), dimension(Num_nu) :: attenuation,seed_logv
real(8), dimension(Num_nu-1) :: seed_invdlog
real(8) :: dtau, path_hi, path_lo, seg_lo, seg_hi, xsrc, xtgt, doprel, weight, amp_p, amp_seed

    seed_logv=dlog(V_seed)
    seed_invdlog=1d0/(seed_logv(2:Num_nu)-seed_logv(1:Num_nu-1))
    call build_transfer(Num_chi,Num_nu,dxcell(:,target_t),tau(:,:,target_t),hist_inv,hist_prefix)
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
    integer :: inu,ilo
    real(8), intent(in), dimension(Num_nu,Num_chi) :: P_src,Seed_src
    real(8), intent(in) :: xsrcpos,swgt
    real(8) :: nusrc,frac

        xtgt = xmid(tgt_chi,target_t)
        if (xsrcpos < xtgt) return
        if (swgt <= 0d0) return
        xsrc = xsrcpos
        doprel = relative_doppler(beta(src_chi,prev_t),beta(tgt_chi,target_t))
        if (doprel <= 0d0) return
        attenuation = 1d0
        call apply_path_tau(Num_chi,Num_nu,xsrc,xtgt,xface(:,target_t), &
                            hist_inv,tau(:,:,target_t),hist_prefix,attenuation)
        amp_p = swgt*doprel**3
        amp_seed = swgt*doprel**2
        ilo = 1
        do inu = 1, Num_nu
            nusrc = V_seed(inu)/doprel
            if (nusrc < V_seed(1) .or. nusrc > V_seed(Num_nu)) cycle
            do while (ilo < Num_nu-1)
                if (V_seed(ilo+1) > nusrc) exit
                ilo = ilo + 1
            end do
            if (seed_invdlog(ilo) <= 0d0) cycle
            frac = (dlog(nusrc)-seed_logv(ilo))*seed_invdlog(ilo)
            pnext(inu,tgt_chi) = pnext(inu,tgt_chi) + amp_p*attenuation(inu) * &
                log_interp(P_src(ilo,src_chi),P_src(ilo+1,src_chi),frac)
            snext(inu,tgt_chi) = snext(inu,tgt_chi) + amp_seed*attenuation(inu) * &
                log_interp(Seed_src(ilo,src_chi),Seed_src(ilo+1,src_chi),frac)
        end do
    end subroutine add_mapped_cell
end subroutine advance_history_stream

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

! 构造当前目标壳层的 chi 单元宽度倒数和累积对数传输。
! Build inverse chi-cell widths and cumulative log transfer for the current target shell.
subroutine build_transfer(Num_chi,Num_nu,dxcell,tau,invcell,prefix)
implicit none
integer, intent(in) :: Num_chi,Num_nu
integer :: ichi,inu
real(8), intent(in), dimension(Num_chi) :: dxcell
real(8), intent(in), dimension(Num_nu,Num_chi) :: tau
real(8), intent(out), dimension(Num_chi) :: invcell
real(8), intent(out), dimension(Num_nu,0:Num_chi) :: prefix

    prefix(:,0) = 0d0
    do ichi = 1, Num_chi
        if (dxcell(ichi) <= 0d0) error stop 'history transfer requires positive dxcell'
        invcell(ichi) = 1d0/dxcell(ichi)
        do inu = 1, Num_nu
            prefix(inu,ichi) = prefix(inu,ichi-1) + &
                dlog(history_transfer_weight(tau(inu,ichi)))
        end do
    end do
end subroutine build_transfer

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
