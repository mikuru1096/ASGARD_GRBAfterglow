!f2py: skip
module electron_shell_transport_common
    use constants
    use electron_energy_coordinate_common, only: electron_coord_log_four_velocity_sq, electron_dxgamma_dcoord
    use electron_transport_common, only: electron_fullhide_flux_split_step, &
                                         electron_fullhide_flux_split_step_nonuniform, &
                                         electron_prepare_exponential_source_remap, &
                                         transport_fullhide_spacetime_sequence => electron_fullhide_spacetime_sequence, &
                                         transport_fullhide_step => electron_fullhide_step
    implicit none
    private

    integer, parameter, public :: electron_solver_fullhide_1d = 1
    integer, parameter, public :: electron_solver_dg_1d = 2

    public :: electron_resolve_1d_solver_id, electron_shell_fullhide_step, electron_shell_flux_split_step
    public :: electron_shell_flux_split_coord_step, electron_shell_dcoord_to_dndgamma_exp_centers
    public :: electron_shell_fullhide_spacetime_sequence

contains

integer function electron_resolve_1d_solver_id(solver_id) result(resolved)
    integer, intent(in), optional :: solver_id

    resolved = electron_solver_fullhide_1d
    if (present(solver_id)) resolved = solver_id
    if (resolved /= electron_solver_fullhide_1d .and. resolved /= electron_solver_dg_1d) &
        error stop 'electron_resolve_1d_solver_id: unsupported electron solver id.'
end function electron_resolve_1d_solver_id

subroutine electron_shell_fullhide_step(Num_gam_e,R_loc,dDR,d_x,dEL_mean,dF1,dN_x_in,dN_x_out)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: R_loc,dDR,d_x,dEL_mean(Num_gam_e-1),dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)

    call transport_fullhide_step(Num_gam_e,R_loc,dDR,d_x,dEL_mean,dF1,dN_x_in,dN_x_out)
end subroutine electron_shell_fullhide_step

subroutine electron_shell_fullhide_spacetime_sequence(Num_gam_e,Num_step,face_coupling,source_step,dN_x_in,dN_x_out)
    integer, intent(in) :: Num_gam_e,Num_step
    real(8), intent(in) :: face_coupling(Num_gam_e-1,Num_step),source_step(Num_gam_e,Num_step),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)

    call transport_fullhide_spacetime_sequence(Num_gam_e,Num_step,face_coupling,source_step,dN_x_in,dN_x_out)
end subroutine electron_shell_fullhide_spacetime_sequence

subroutine electron_shell_flux_split_step(Num_gam_e,dDR,d_x,dEl,adiabatic_rate,dF1,dN_x_in,dN_x_out)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: dDR,d_x,dEl(Num_gam_e),adiabatic_rate,dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    real(8) :: face_speed(Num_gam_e-1)

    face_speed=((dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two+adiabatic_rate)/dlog(ten)
    call electron_fullhide_flux_split_step(Num_gam_e,dDR,d_x,face_speed,dF1,dN_x_in,dN_x_out,.true.)
end subroutine electron_shell_flux_split_step

subroutine electron_shell_flux_split_coord_step(Num_gam_e,dDR,coord_edge,coord_scale,dEl,adiabatic_rate, &
                                                dF1,dN_coord_in,dN_coord_out)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: dDR,coord_edge(Num_gam_e+1),coord_scale,dEl(Num_gam_e),adiabatic_rate
    real(8), intent(in) :: dF1(Num_gam_e),dN_coord_in(Num_gam_e)
    real(8), intent(out) :: dN_coord_out(Num_gam_e)
    real(8) :: face_speed(Num_gam_e-1),face_coord,face_jac
    integer :: i

    do i = 1, Num_gam_e - 1
        face_coord = coord_edge(i + 1)
        face_jac = dlog(ten)*electron_dxgamma_dcoord(electron_coord_log_four_velocity_sq, coord_scale, face_coord)
        face_speed(i) = ((dEl(i) + dEl(i + 1))/two + adiabatic_rate)/face_jac
    enddo
    call electron_fullhide_flux_split_step_nonuniform(Num_gam_e,dDR,coord_edge,face_speed,dF1, &
                                                      dN_coord_in,dN_coord_out,.true.)
end subroutine electron_shell_flux_split_coord_step

subroutine electron_shell_dcoord_to_dndgamma_exp_centers(Num_gam_e,coord_edge,coord_scale,gam_e,dN_coord,dN_gam_e)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: coord_edge(Num_gam_e+1),coord_scale,gam_e(Num_gam_e),dN_coord(Num_gam_e)
    real(8), intent(out) :: dN_gam_e(Num_gam_e)
    real(8) :: slope(Num_gam_e),center(Num_gam_e),prefix(0:Num_gam_e),coord_mid,dxdy
    integer :: i

    call electron_prepare_exponential_source_remap(Num_gam_e,coord_edge,dN_coord,slope,center,prefix)
    do i = 1, Num_gam_e
        coord_mid = 0.5d0*(coord_edge(i) + coord_edge(i + 1))
        dxdy = electron_dxgamma_dcoord(electron_coord_log_four_velocity_sq, coord_scale, coord_mid)
        dN_gam_e(i) = center(i)/(gam_e(i)*dlog(ten)*dxdy)
    enddo
end subroutine electron_shell_dcoord_to_dndgamma_exp_centers

end module electron_shell_transport_common
