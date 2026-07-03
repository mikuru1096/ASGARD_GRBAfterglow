!f2py: skip
module hadronic_remap
    use constants
    use hadronic_base, only: build_edges
    implicit none

contains

! 冷却后把 cell content 保守重映射到同一 gamma grid。
! Remap cooled cell content conservatively back onto the same gamma grid.
subroutine remap_loggamma(ng,gamma,prev,qinj,ltot,dt,next)
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng) :: gamma,prev,qinj,ltot
    real(8), intent(in) :: dt
    real(8), intent(out), dimension(ng) :: next
    real(8), dimension(ng+1) :: edge
    real(8) :: content,gcool,dgamma
    integer :: i,target

    call build_edges(ng,gamma,edge)
    next=0d0
    ! 先按 cell 宽度保存粒子数 content。 / First store particle content using cell widths.
    do i=1,ng
        dgamma=edge(i+1)-edge(i)
        if (dgamma <= 0d0) error stop "hadronic remap transport requires positive gamma cell width."
        content=prev(i)*dgamma
        gcool=gamma(i)-ltot(i)*dt
        target=remap_target(ng,edge,gcool)
        if (target >= 1 .and. target <= ng) next(target)=next(target)+content
    end do
    ! 再除回 cell 宽度并加入本步源项。 / Then restore density units and add current sources.
    do i=1,ng
        dgamma=edge(i+1)-edge(i)
        next(i)=next(i)/dgamma+qinj(i)
    end do
end subroutine remap_loggamma

! 返回冷却后 gamma 所在的目标 cell，下界外返回 0。
! Return the destination cell for cooled gamma; values below range return 0.
integer function remap_target(ng,edge,value)
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng+1) :: edge
    real(8), intent(in) :: value
    integer :: lo,hi,mid

    if (value < edge(1) .or. value >= edge(ng+1)) then
        remap_target=0
        return
    end if
    lo=1
    hi=ng+1
    do while (hi-lo > 1)
        mid=(lo+hi)/2
        if (edge(mid) <= value) then
            lo=mid
        else
            hi=mid
        end if
    end do
    remap_target=lo
end function remap_target

end module hadronic_remap
