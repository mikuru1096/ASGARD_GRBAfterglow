!f2py: skip
module electron_cooling_kernel
  use constants
  use electron_ssa_kernel, only: electron_ssa_loss
  use electron_ic_kernel, only: electron_ic_loss
  use electron_y_kernel, only: electron_y_nakar, electron_y_fan
  private

  public :: forward_cooling

contains

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
real(8), dimension(ng) :: comp,dotssa
integer :: ic,ig,nu0,g0
real(8) :: cscale,ssascale,fr,qvol

    if (mode /= 2) then
        select case(iy)
        case(0,3)
            do ic=1,nchi
                g0=(ic-1)*ng
                aux(g0+1:g0+ng)=0d0
            end do
        case(1)
            do ic=1,nchi
                nu0=(ic-1)*nnu
                g0=(ic-1)*ng
                call electron_ic_loss(ng,nnu,nthr,gam,vseed,seed(nu0+1),aux(g0+1))
            end do
        case(2)
            do ic=1,nchi
                nu0=(ic-1)*nnu
                g0=(ic-1)*ng
                call electron_y_nakar(ng,nnu,nthr,gam,vseed,psyn(nu0+1),aux(g0+1))
            end do
        case default
            error stop 'forward_cooling: index_Y must be 0, 1, 2, or 3.'
        end select
    end if
    if (mode == 0) return

    cscale=1d0/(beta*rg)
    ssascale=cscale/para_c
    fr=1.35d-19*db**2*cscale/pi
    if (iy /= 0) call electron_ssa_loss(db,ng,nnu,nchi,nthr,gam,vseed,seedssa,del)

    do ic=1,nchi
        g0=(ic-1)*ng
        if (iy == 0) then
            dotssa=0d0
        else
            do ig=1,ng
                dotssa(ig)=del(g0+ig)
            end do
        end if
        select case(iy)
        case(0)
            do ig=1,ng
                del(g0+ig)=fr*gam(ig)
            end do
        case(1)
            do ig=1,ng
                del(g0+ig)=(fr+(aux(g0+ig)-dotssa(ig))*ssascale)*gam(ig)
            end do
        case(2)
            qvol=4d0*pi*rloc*rloc*para_c
            do ig=1,ng
                comp(ig)=1d0+aux(g0+ig)/qvol/(4d0*rg*rg*dens*Para_m_p_E)
                del(g0+ig)=(fr*comp(ig)-dotssa(ig)*ssascale)*gam(ig)
            end do
        case(3)
            call electron_y_fan(ee,eb,p,db,gm,gc,gmax,ng,gam,comp)
            do ig=1,ng
                del(g0+ig)=(fr*(1d0+comp(ig))-dotssa(ig)*ssascale)*gam(ig)
            end do
        case default
            error stop 'forward_cooling: index_Y must be 0, 1, 2, or 3.'
        end select
    end do
end subroutine forward_cooling

end module electron_cooling_kernel
