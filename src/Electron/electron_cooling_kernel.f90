!f2py: skip
module electron_cooling_kernel
  use constants
  use electron_ssa_kernel, only: electron_ssa_loss, invalidate_ssa_cache
  use electron_ic_kernel, only: electron_ic_loss, electron_ic_loss_batch, invalidate_ic_cache
  use electron_y_kernel, only: electron_y_nakar, electron_y_nakar_batch, electron_y_fan, invalidate_y_cache
  private

  public :: forward_cooling, cooling_reset

contains

! 清空当前 OpenMP worker 的电子冷却积分缓存，供外层 structured worker 启动时调用。
! Clear the current OpenMP worker's cooling quadrature caches when a structured worker starts.
subroutine cooling_reset()
implicit none

    call invalidate_ic_cache()
    call invalidate_ssa_cache()
    call invalidate_y_cache()
end subroutine cooling_reset

! 正向激波电子冷却统一入口：mode=0 只准备 Compton auxiliary，mode=1 由 seed 计算 auxiliary 并组装冷却，
! mode=2 复用调用方传入的 auxiliary 组装冷却。
! Forward-shock cooling entry: mode=0 prepares Compton auxiliary only, mode=1 computes auxiliary
! from seeds and assembles losses, and mode=2 reuses caller-supplied auxiliary for loss assembly.
subroutine forward_cooling(mode,iy,ee,eb,p,db,gm,gc,gmax,rloc,rg,beta,dens,ng,nnu,nchi,nthr, &
                           gam,vseed,psyn,seed,seedssa,aux,del)
implicit none
integer, intent(in) :: mode,iy,ng,nnu,nchi,nthr
real(8), intent(in) :: ee,eb,p,db,gm,gc,gmax,rloc,rg,beta,dens
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in), dimension(nnu) :: vseed
real(8), intent(in), dimension(*) :: psyn,seed,seedssa
real(8), intent(inout), dimension(*) :: aux
real(8), intent(out), dimension(*) :: del
real(8), dimension(ng) :: comp
integer :: ic,g0
real(8) :: cscale,ssascale,fr,qvol

    select case(iy)
    case(0)
        if (mode /= 2) then
            do ic=1,nchi
                g0=(ic-1)*ng
                aux(g0+1:g0+ng)=0d0
            end do
        end if
    case(1)
        if (mode /= 2) then
            if (nchi == 1) then
                call electron_ic_loss(ng,nnu,nthr,gam,vseed,seed,aux)
            else
                call electron_ic_loss_batch(ng,nnu,nchi,nthr,gam,vseed,seed,aux)
            end if
        end if
    case(2)
        if (mode /= 2) then
            if (nchi == 1) then
                call electron_y_nakar(ng,nnu,nthr,gam,vseed,psyn,aux)
            else
                call electron_y_nakar_batch(ng,nnu,nchi,nthr,gam,vseed,psyn,aux)
            end if
        end if
    case(3)
        if (mode /= 2) then
            do ic=1,nchi
                g0=(ic-1)*ng
                aux(g0+1:g0+ng)=0d0
            end do
        end if
    case default
        error stop 'forward_cooling: index_Y must be 0, 1, 2, or 3.'
    end select
    if (mode == 0) return

    cscale=1d0/(beta*rg)
    ssascale=cscale/para_c
    fr=1.35d-19*db**2*cscale/pi
    if (iy == 0) then
        do ic=1,nchi
            g0=(ic-1)*ng
            del(g0+1:g0+ng)=fr*gam
        end do
        return
    end if

    call electron_ssa_loss(db,ng,nnu,nchi,nthr,gam,vseed,seedssa,del)
    if (iy == 1) then
        do ic=1,nchi
            g0=(ic-1)*ng
            del(g0+1:g0+ng)=(fr+(aux(g0+1:g0+ng)-del(g0+1:g0+ng))*ssascale)*gam
        end do
    else if (iy == 2) then
        qvol=4d0*pi*rloc*rloc*para_c
        do ic=1,nchi
            g0=(ic-1)*ng
            comp=1d0+aux(g0+1:g0+ng)/qvol/(db*db/(8d0*pi))
            del(g0+1:g0+ng)=(fr*comp-del(g0+1:g0+ng)*ssascale)*gam
        end do
    else
        call electron_y_fan(ee,eb,p,db,gm,gc,gmax,ng,gam,comp)
        do ic=1,nchi
            g0=(ic-1)*ng
            del(g0+1:g0+ng)=(fr*(1d0+comp)-del(g0+1:g0+ng)*ssascale)*gam
        end do
    end if
end subroutine forward_cooling

end module electron_cooling_kernel
