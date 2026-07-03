!f2py: skip
module electron_shell_transport
    use constants
    use electron_coord_common, only: coord_fourvel, dxg_dcoord
    use electron_transport_common, only: flux_split_nonuniform, &
                                         prepare_exp_source
    implicit none
    private

    ! 1D shell 电子输运公共入口：选择后端、推进坐标空间通量、投影回 dN/dgamma。
    ! Shared 1D shell electron transport entry points: select backend, advance coordinate flux, project to dN/dgamma.
    integer, parameter, public :: solver_fullhide = 1
    integer, parameter, public :: solver_dg = 2

    public :: resolve_solver
    public :: shell_coord_step, coord_to_dgamma

contains

! 解析可选 solver id；缺省走 fullhide 1D。
! Resolve an optional solver id; fullhide 1D is the default.
integer function resolve_solver(solver_id) result(resolved)
    integer, intent(in), optional :: solver_id

    resolved = solver_fullhide
    if (present(solver_id)) resolved = solver_id
    if (resolved /= solver_fullhide .and. resolved /= solver_dg) &
        error stop 'resolve_solver: unsupported electron solver id.'
end function resolve_solver

! 在 four-velocity 坐标中推进单个 shell 的 conservative flux-split step。
! Advance a single shell with the conservative flux-split step in four-velocity coordinates.
subroutine shell_coord_step(Num_gam_e,dDR,coord_edge,coord_scale,dEl,adiabatic_rate, &
                            dF1,n_in,n_out)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in), dimension(Num_gam_e+1) :: coord_edge
    real(8), intent(in), dimension(Num_gam_e) :: dEl
    real(8), intent(in) :: dDR,coord_scale,adiabatic_rate
    real(8), intent(in), dimension(Num_gam_e) :: dF1,n_in
    real(8), intent(out), dimension(Num_gam_e) :: n_out
    real(8), dimension(Num_gam_e-1) :: vface
    real(8) :: yface,jface
    integer :: i

    do i = 1, Num_gam_e - 1
        yface = coord_edge(i + 1)
        jface = dlog(1d1)*dxg_dcoord(coord_fourvel, coord_scale, yface)
        vface(i) = ((dEl(i) + dEl(i + 1))/2d0 + adiabatic_rate)/jface
    enddo
    call flux_split_nonuniform(Num_gam_e,dDR,coord_edge,vface,dF1, &
                                                      n_in,n_out,.true.)
end subroutine shell_coord_step

! 把坐标空间单元含量投影回 gamma 中心的 dN/dgamma 诊断量。
! Project coordinate-cell content back to center-sampled dN/dgamma diagnostics.
subroutine coord_to_dgamma(Num_gam_e,coord_edge,coord_scale,gam_e,ncoord,ndg)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in), dimension(Num_gam_e+1) :: coord_edge
    real(8), intent(in), dimension(Num_gam_e) :: gam_e,ncoord
    real(8), intent(in) :: coord_scale
    real(8), intent(out), dimension(Num_gam_e) :: ndg
    real(8), dimension(Num_gam_e) :: exp_slope,center_density
    real(8), dimension(0:Num_gam_e) :: exp_prefix
    real(8) :: ymid,dxdy
    integer :: i

    call prepare_exp_source(Num_gam_e,coord_edge,ncoord, &
                                                   exp_slope,center_density,exp_prefix)
    do i = 1, Num_gam_e
        ymid = 0.5d0*(coord_edge(i) + coord_edge(i + 1))
        dxdy = dxg_dcoord(coord_fourvel, coord_scale, ymid)
        ndg(i) = center_density(i)/(gam_e(i)*dlog(1d1)*dxdy)
    enddo
end subroutine coord_to_dgamma

end module electron_shell_transport
