module electron_transport_common
    use constants
    implicit none
contains

subroutine electron_prepare_implicit_coeffs_common(Num_gam_e,diag_base,up,principal,temp1)
    implicit none
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: diag_base,up(Num_gam_e-1)
    real(8), intent(out) :: principal(Num_gam_e),temp1(Num_gam_e-1)

    principal(2:Num_gam_e)=diag_base-up
    principal(1)=principal(2)
    temp1=up/(principal(2:Num_gam_e)+principal(1:Num_gam_e-1))*two
end subroutine electron_prepare_implicit_coeffs_common

subroutine electron_backward_sweep_common(Num_gam_e,temp1,temp2,x)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: temp1(Num_gam_e-1),temp2(Num_gam_e)
    real(8), intent(out) :: x(Num_gam_e)

    x(Num_gam_e)=temp2(Num_gam_e)
    do I_gam_e=Num_gam_e-1,1,-1
        x(I_gam_e)=max(zero,temp2(I_gam_e)-temp1(I_gam_e)*x(I_gam_e+1))
    end do
end subroutine electron_backward_sweep_common

end module electron_transport_common
