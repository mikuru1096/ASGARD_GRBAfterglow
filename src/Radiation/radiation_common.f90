module radiation_common
    use constants
    implicit none

contains

subroutine compute_simpson_weights(weights, n)
    integer, intent(in) :: n
    integer :: i
    real(8), intent(out) :: weights(n)

    weights = 1.0d0
    if (n >= 3) then
        do i = 2, n - 1
            if (mod(i, 2) == 0) then
                weights(i) = 4.0d0
            else
                weights(i) = 2.0d0
            end if
        end do
    end if
end subroutine compute_simpson_weights

real(8) function radiation_powerlaw_interp(v0,v1,y0,y1,v)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: v0,v1,y0,y1,v
    real(8) :: slope

    if (v <= v0) then
        radiation_powerlaw_interp=y0
        return
    end if
    if (v >= v1) then
        radiation_powerlaw_interp=y1
        return
    end if

    if (v1 <= v0) then
        radiation_powerlaw_interp=0.5d0*(y0+y1)
    else if (y0 > zero .and. y1 > zero) then
        slope=dlog(y1/y0)/dlog(v1/v0)
        radiation_powerlaw_interp=y0*(v/v0)**slope
    else
        radiation_powerlaw_interp=y0+(y1-y0)*(v-v0)/(v1-v0)
    end if
end function radiation_powerlaw_interp

subroutine radiation_transfer_factor(Tau, factor)
    implicit real(8)(A-H,O-Z)
    real(8), intent(inout) :: Tau
    real(8), intent(out) :: factor

    if ((Tau - 1.0d-4) < 1.0d-5) Tau = 1.0d-4
    factor = (one - dexp(-Tau)) / Tau
end subroutine radiation_transfer_factor

subroutine radiation_prepare_annihilation_grid(V_seed, Num_nu, ep1, ep2, dVloc, V_mid)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_nu
    real(8), intent(in) :: V_seed(Num_nu)
    real(8), intent(out) :: ep1(1,Num_nu),ep2(Num_nu-1,1),dVloc(Num_nu-1),V_mid(Num_nu-1)
    real(8) :: x_seed(Num_nu)

    para_hEme=Para_h/para_m_energy
    ep1(1,:)=para_hEme*V_seed
    x_seed=dlog(V_seed)
    V_mid=dexp(0.5d0*(x_seed(1:Num_nu-1)+x_seed(2:Num_nu)))
    ep2(:,1)=para_hEme*V_mid
    dVloc=V_mid*(x_seed(2:Num_nu)-x_seed(1:Num_nu-1))
end subroutine radiation_prepare_annihilation_grid

subroutine radiation_external_density(A_star,dNe_ISM,R_loc,R0,dNe)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: A_star,dNe_ISM,R_loc,R0
    real(8), intent(out) :: dNe

    if (A_star >= zero) then
        dNe_wind=A_star*3.0d35/R_loc**2
        if (dNe_wind <= dNe_ISM/4.0d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
    end if

    if (R_loc < R0) then
        dNe=A_star*3.0d35/R0**2
    end if
end subroutine radiation_external_density

end module radiation_common
