!f2py: skip
module electron_cooling_kernel
  use constants
  use electron_ssa_kernel, only: electron_ssa_loss
  use electron_ic_kernel, only: electron_ic_loss
  use electron_y_kernel, only: electron_y_nakar, electron_y_fan
  private

  public :: get_forward_cooling, forward_cooling_aux
  public :: forward_cooling_batch

contains
! 为每个 chi 列准备 Compton 辅助量；index_Y 决定 IC/Nakar/Fan 路径。
! Prepare the Compton auxiliary field per chi column; index_Y selects IC, Nakar, or Fan cooling.
subroutine forward_cooling_aux(iy,ng,nnu,nchi,nthr,gam,vseed,psyn,seed,aux)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: iy,ng,nnu,nchi,nthr
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in), dimension(nnu) :: vseed
real(8), intent(in), dimension(nnu,nchi) :: psyn,seed
real(8), intent(out), dimension(ng,nchi) :: aux
integer :: ic

    aux=0d0
    select case(iy)
    case(0,3)
    case(1)
        do ic=1,nchi
            call electron_ic_loss(ng,nnu,nthr,gam,vseed,seed(:,ic),aux(:,ic))
        end do
    case(2)
        do ic=1,nchi
            call electron_y_nakar(ng,nnu,nthr,gam,vseed,psyn(:,ic),aux(:,ic))
        end do
    case default
        error stop 'forward_cooling_aux: index_Y must be 0, 1, 2, or 3.'
    end select
end subroutine forward_cooling_aux

! 把 synchrotron、SSA 和 Compton 项合成每个 chi 列的电子冷却率。
! Assemble synchrotron, SSA, and Compton terms into the electron cooling rate per chi column.
subroutine forward_cooling_batch(iy,ee,eb,p,db,gm,gc,gmax,rloc,rg,beta,dens,ng,nnu,nchi,nthr,gam,vseed,seedssa,aux,del)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: iy,ng,nnu,nchi,nthr
real(8), intent(in) :: ee,eb,p,db,gm,gc,rloc,rg,beta,dens
real(8), intent(inout) :: gmax
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in), dimension(nnu) :: vseed
real(8), intent(in), dimension(nnu,nchi) :: seedssa
real(8), intent(in), dimension(ng,nchi) :: aux
real(8), intent(out), dimension(ng,nchi) :: del
real(8), dimension(ng) :: comp
real(8), dimension(ng,nchi) :: dotssa
real(8) :: gcell
integer :: ic

    if (iy == 0) then
        dotssa=0d0
    else
        call electron_ssa_loss(db,ng,nnu,nchi,nthr,gam,vseed,seedssa,dotssa)
    end if
    do ic=1,nchi
       gcell=gmax
       call forward_cooling_terms(iy,ee,eb,p,db,gm,gc,gcell,rloc,rg,beta,dens,ng,gam,comp,dotssa(:,ic),aux(:,ic),del(:,ic))
    end do
end subroutine forward_cooling_batch

! 单列冷却率公式：d gamma/dt = (synch*Y - SSA)*gamma。
! Single-column cooling formula: d gamma/dt = (synch*Y - SSA)*gamma.
subroutine forward_cooling_terms(iy,ee,eb,p,db,gm,gc,gmax,rloc,rg,beta,dens,ng,gam,comp,dotssa,aux,del)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: iy,ng
real(8), intent(in) :: ee,eb,p,db,gm,gc,rloc,rg,beta,dens
real(8), intent(inout) :: gmax
real(8), intent(in), dimension(ng) :: gam
real(8), intent(inout), dimension(ng) :: comp
real(8), intent(in), dimension(ng) :: dotssa,aux
real(8), intent(out), dimension(ng) :: del

    cscale=1d0/(beta*rg)
    ssascale=cscale/para_c
    fr=1.35d-19*db**2*cscale/pi

    select case(iy)
    case(0)
        del=fr*gam
    case(1)
        del=(fr+(aux-dotssa)*ssascale)*gam
    case(2)
        qvol=4d0*pi*rloc*rloc*para_c
        comp=1d0+aux/qvol/(4d0*rg*rg*dens*Para_m_p_E)
        gmax=gmax/sqrt(comp(ng))
        del=(fr*comp-dotssa*ssascale)*gam
    case(3)
        call electron_y_fan(ee,eb,p,db,gm,gc,gmax,ng,gam,comp)
        comp=1d0+comp
        gmax=gmax/sqrt(comp(ng))
        del=(fr*comp-dotssa*ssascale)*gam
    case default
        error stop 'forward_cooling_terms: index_Y must be 0, 1, 2, or 3.'
    end select
end subroutine forward_cooling_terms

! 单壳层主入口：包装成一列 chi，再复用批量冷却装配。
! Single-shell entry: wrap the state as 1 chi column and reuse the batch assembler.
subroutine get_forward_cooling(iy,ee,eb,p,db,gm,gc,gmax,rloc,rg,beta,dens,ng,nnu,nthr,gam,vseed,psyn,seed,del)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: iy,ng,nnu,nthr
real(8), intent(in) :: ee,eb,p,db,gm,gc,rloc,rg,beta,dens
real(8), intent(inout) :: gmax
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in), dimension(nnu) :: vseed,psyn,seed
real(8), intent(out), dimension(ng) :: del
real(8), dimension(nnu,1) :: pscol,secol
real(8), dimension(ng,1) :: auxcol,delcol

    pscol(:,1)=psyn
    secol(:,1)=seed
    call forward_cooling_aux(iy,ng,nnu,1,nthr,gam,vseed,pscol,secol,auxcol)
    call forward_cooling_batch(iy,ee,eb,p,db,gm,gc,gmax,rloc,rg,beta,dens,ng,nnu,1,nthr,gam,vseed,secol,auxcol,delcol)
    del=delcol(:,1)
end subroutine get_forward_cooling


end module electron_cooling_kernel
