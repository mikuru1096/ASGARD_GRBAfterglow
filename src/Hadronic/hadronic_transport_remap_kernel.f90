!f2py: skip
module hadronic_transport_remap_kernel
    use constants
    use hadronic_common, only: hadronic_build_gamma_edges
    implicit none

contains

subroutine hadronic_advance_energy_loggamma_remap(num_gamma,gamma,dn_prev,q_inj,loss_total,dt_s,dn_next)
    implicit none
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: gamma(num_gamma),dn_prev(num_gamma),q_inj(num_gamma),loss_total(num_gamma),dt_s
    real(8), intent(out) :: dn_next(num_gamma)
    real(8) :: gamma_edge(num_gamma+1),content,cooled_gamma,dgamma
    integer :: i,target

    call hadronic_build_gamma_edges(num_gamma,gamma,gamma_edge)
    dn_next=zero
    call deposit_cooled_bin_content()
    call restore_density_units()

contains

    subroutine deposit_cooled_bin_content()
    implicit none

    do i=1,num_gamma
        dgamma=gamma_edge(i+1)-gamma_edge(i)
        if (dgamma <= zero) error stop "hadronic remap transport requires positive gamma cell width."
        content=dn_prev(i)*dgamma
        cooled_gamma=gamma(i)-loss_total(i)*dt_s
        target=hadronic_remap_target(num_gamma,gamma_edge,cooled_gamma)
        if (target >= 1 .and. target <= num_gamma) dn_next(target)=dn_next(target)+content
    end do
    end subroutine deposit_cooled_bin_content

    subroutine restore_density_units()
    implicit none

    do i=1,num_gamma
        dgamma=gamma_edge(i+1)-gamma_edge(i)
        dn_next(i)=dn_next(i)/dgamma+q_inj(i)
    end do
    end subroutine restore_density_units
end subroutine hadronic_advance_energy_loggamma_remap

integer function hadronic_remap_target(num_gamma,gamma_edge,value)
    implicit none
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: gamma_edge(num_gamma+1),value
    integer :: lo,hi,mid

    if (value < gamma_edge(1) .or. value >= gamma_edge(num_gamma+1)) then
        hadronic_remap_target=0
        return
    end if
    lo=1
    hi=num_gamma+1
    do while (hi-lo > 1)
        mid=(lo+hi)/2
        if (gamma_edge(mid) <= value) then
            lo=mid
        else
            hi=mid
        end if
    end do
    hadronic_remap_target=lo
end function hadronic_remap_target

end module hadronic_transport_remap_kernel
