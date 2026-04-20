module dynamics_common
    use constants
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none

contains

subroutine dynamics_deceleration_radius(A_star,dNe_ISM,Eta_0,DM_0,R_dec)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: A_star,dNe_ISM,Eta_0,DM_0
    real(8), intent(out) :: R_dec

    R_dec_ISM=(dNe_ISM*Eta_0/DM_0)**(-one/3d0)
    if (A_star > zero) then
        R_dec_wind=DM_0/(2.0d35*A_star*Eta_0)
        R_dec=min(R_dec_wind,R_dec_ISM)
    else
        R_dec=R_dec_ISM
    end if
end subroutine dynamics_deceleration_radius

subroutine dynamics_external_density_base(A_star,dNe_ISM,RR,dNe)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: A_star,dNe_ISM,RR
    real(8), intent(out) :: dNe

    if (A_star > zero) then
        dNe_wind=A_star*3.0d35/RR**2
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
    end if
end subroutine dynamics_external_density_base

subroutine dynamics_rk4_coefficients(HH,A)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: HH
    real(8), intent(out) :: A(4)

    A(1)=0.5d0*HH
    A(2)=A(1)
    A(3)=HH
    A(4)=HH
end subroutine dynamics_rk4_coefficients

subroutine dynamics_rk4_error_n(Y,G,M,P)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: M
    integer :: I
    real(8), intent(in) :: Y(M),G(M)
    real(8), intent(out) :: P

    P=0.0d0
    do I=1,M
        if (.not. ieee_is_finite(Y(I)) .or. .not. ieee_is_finite(G(I))) then
            P=huge(one)
            return
        end if
        if (Y(I)+G(I) <= zero) then
            P=huge(one)
            return
        end if
        Q=2.0d0*abs(Y(I)-G(I))/(Y(I)+G(I))
        if (Q > P) P=Q
    end do
end subroutine dynamics_rk4_error_n

subroutine dynamics_rk4_error4(Y,G,P)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: Y(4),G(4)
    real(8), intent(out) :: P

    call dynamics_rk4_error_n(Y,G,4,P)
end subroutine dynamics_rk4_error4

end module dynamics_common
