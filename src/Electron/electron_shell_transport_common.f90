!f2py: skip
module electron_shell_transport_common
    use constants
    use electron_transport_common, only: electron_fullhide_flux_split_step, &
                                         transport_fullhide_spacetime_sequence => electron_fullhide_spacetime_sequence, &
                                         transport_fullhide_step => electron_fullhide_step
    implicit none
    private

    integer, parameter, public :: electron_solver_fullhide_1d = 1
    integer, parameter, public :: electron_solver_dg_1d = 2

    public :: electron_resolve_1d_solver_id, electron_shell_fullhide_step, electron_shell_flux_split_step
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

end module electron_shell_transport_common
