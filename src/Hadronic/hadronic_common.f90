!f2py: skip
module hadronic_common
    use constants
    implicit none
    real(8), parameter :: hadronic_eta_acc_floor = 1d-12
    real(8), parameter :: hadronic_bfield_floor = 1d-30

contains

subroutine hadronic_build_gamma_p_grid(Num_gam_p,gam_p_min,gam_p_max,gam_p)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p
    real(8), intent(in) :: gam_p_min,gam_p_max
    real(8), intent(out) :: gam_p(Num_gam_p)
    integer :: I_gam
    real(8) :: x_min,x_max,dx

    if (Num_gam_p <= 1) then
        gam_p(1)=max(gam_p_min,one)
        return
    end if

    x_min=dlog10(max(gam_p_min,one))
    x_max=dlog10(max(gam_p_max,ten))
    dx=(x_max-x_min)/dble(Num_gam_p-1)
    do I_gam=1,Num_gam_p
        gam_p(I_gam)=ten**(x_min+dx*dble(I_gam-1))
    end do
end subroutine hadronic_build_gamma_p_grid

subroutine hadronic_initial_density(Num_gam_p,dN_gam_p)
    integer, intent(in) :: Num_gam_p
    real(8), intent(out) :: dN_gam_p(Num_gam_p)

    dN_gam_p=zero
end subroutine hadronic_initial_density

subroutine hadronic_source_bounds(Num_gam_p,gam_p,gam_p_min,gam_p_max,i_lo,i_hi)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p
    real(8), intent(in) :: gam_p(Num_gam_p),gam_p_min,gam_p_max
    integer, intent(out) :: i_lo,i_hi
    integer :: I_gam

    i_lo=1
    i_hi=Num_gam_p
    do I_gam=1,Num_gam_p
        if (gam_p(I_gam) >= gam_p_min) then
            i_lo=I_gam
            exit
        end if
    end do
    do I_gam=Num_gam_p,1,-1
        if (gam_p(I_gam) <= gam_p_max) then
            i_hi=I_gam
            exit
        end if
    end do
end subroutine hadronic_source_bounds

subroutine hadronic_build_gamma_edges(Num_gam_p,gam_p,gam_edge)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p
    real(8), intent(in) :: gam_p(Num_gam_p)
    real(8), intent(out) :: gam_edge(Num_gam_p+1)
    integer :: I_gam

    if (Num_gam_p == 1) then
        gam_edge(1)=max(one,0.5d0*gam_p(1))
        gam_edge(2)=2d0*gam_p(1)
        return
    end if

    gam_edge(1)=gam_p(1)*dsqrt(gam_p(1)/gam_p(2))
    do I_gam=2,Num_gam_p
        gam_edge(I_gam)=dsqrt(gam_p(I_gam-1)*gam_p(I_gam))
    end do
    gam_edge(Num_gam_p+1)=gam_p(Num_gam_p)*dsqrt(gam_p(Num_gam_p)/gam_p(Num_gam_p-1))
end subroutine hadronic_build_gamma_edges

real(8) function hadronic_shell_dt(R_tobs,i_shell)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: i_shell
    real(8), intent(in) :: R_tobs(*)

    if (i_shell <= 1) then
        hadronic_shell_dt=max(R_tobs(1),one)
    else
        hadronic_shell_dt=max(R_tobs(i_shell)-R_tobs(i_shell-1),one)
    end if
end function hadronic_shell_dt

real(8) function hadronic_dynamical_time(R_loc,Gamma_bulk)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: R_loc,Gamma_bulk

    hadronic_dynamical_time=max(R_loc/(max(Gamma_bulk,one)*Para_c),one)
end function hadronic_dynamical_time

real(8) function hadronic_gamma_p_max(B_field_g,t_dyn_s,eta_acc)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: B_field_g,t_dyn_s,eta_acc
    real(8) :: gam_dyn,gam_syn

    gam_dyn=Para_e*B_field_g*t_dyn_s/(max(eta_acc,hadronic_eta_acc_floor)*Para_m_p*Para_c)
    gam_syn=dsqrt(6d0*pi*Para_e/(max(eta_acc,hadronic_eta_acc_floor)*Para_sigmaT* &
            max(B_field_g,hadronic_bfield_floor))) * Para_m_p_div_m_e
    hadronic_gamma_p_max=max(ten,min(gam_dyn,gam_syn))
end function hadronic_gamma_p_max

end module hadronic_common
