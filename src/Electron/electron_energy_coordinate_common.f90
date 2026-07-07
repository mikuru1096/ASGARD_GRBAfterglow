!f2py: skip
module electron_coord_common
    use constants
    implicit none
    private

    integer, parameter, public :: coord_loggamma=1,coord_fourvel=2
    real(8), parameter, public :: fourvel_scale=2d0

  public :: coord_from_xg, xg_from_coord, gamma_from_coord, dxg_dcoord
  public :: build_fourvel_grid

contains

! 在 log-gamma 与 log four-velocity 坐标之间转换网格坐标。
! Convert grid coordinates between log-gamma and log four-velocity forms.
pure real(8) function coord_from_xg(kind,scale,xg) result(coord)
    integer, intent(in) :: kind
    real(8), intent(in) :: scale,xg
    real(8) :: gamma,fourv

    select case (kind)
    case (coord_loggamma)
        coord=xg
    case (coord_fourvel)
        gamma=dexp(xg)
        fourv=gamma*gamma-1d0
        coord=dlog(1d0+fourv/scale)
    end select
end function coord_from_xg

pure real(8) function xg_from_coord(kind,scale,coord) result(xg)
    integer, intent(in) :: kind
    real(8), intent(in) :: scale,coord

    select case (kind)
    case (coord_loggamma)
        xg=coord
    case (coord_fourvel)
        xg=dlog(dsqrt(1d0+scale*(dexp(coord)-1d0)))
    end select
end function xg_from_coord

pure real(8) function gamma_from_coord(kind,scale,coord) result(gamma)
    integer, intent(in) :: kind
    real(8), intent(in) :: scale,coord

    select case (kind)
    case (coord_loggamma)
        gamma=dexp(coord)
    case (coord_fourvel)
        gamma=dsqrt(1d0+scale*(dexp(coord)-1d0))
    end select
end function gamma_from_coord

pure real(8) function dxg_dcoord(kind,scale,coord) result(dxdy)
    integer, intent(in) :: kind
    real(8), intent(in) :: scale,coord
    real(8) :: gamma

    select case (kind)
    case (coord_loggamma)
        dxdy=1d0
    case (coord_fourvel)
        gamma=dsqrt(1d0+scale*(dexp(coord)-1d0))
        dxdy=scale*dexp(coord)/(2d0*gamma*gamma)
    end select
end function dxg_dcoord

! 构造 log four-velocity 网格，同时返回对应 log-gamma cell edge。
! Build the log four-velocity grid and return matching log-gamma cell edges.
subroutine build_fourvel_grid(ng,gmin,gmax,gscale,gam,coord_edge,xedge)
    integer, intent(in) :: ng
    integer :: i
    real(8), intent(in) :: gmin,gmax,gscale
    real(8), intent(out), dimension(ng) :: gam
    real(8), intent(out), dimension(ng+1) :: coord_edge,xedge
    real(8) :: cmin,cmax,cmid,scale

    if (gscale <= 1d0) error stop 'build_fourvel_grid requires gscale > 1.'
    if (gmax <= gmin) error stop 'build_fourvel_grid requires gmax > gmin.'
    scale=gscale*gscale-1d0
    cmin=dlog(1d0+(gmin*gmin-1d0)/scale)
    cmax=dlog(1d0+(gmax*gmax-1d0)/scale)
    do i=1,ng+1
        coord_edge(i)=cmin+(cmax-cmin)*dble(i-1)/dble(ng)
        xedge(i)=dlog(dsqrt(1d0+scale*(dexp(coord_edge(i))-1d0)))
    enddo
    do i=1,ng
        cmid=0.5d0*(coord_edge(i)+coord_edge(i+1))
        gam(i)=dsqrt(1d0+scale*(dexp(cmid)-1d0))
    enddo
end subroutine build_fourvel_grid

end module electron_coord_common
