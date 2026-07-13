module dynamics_density_profile
    use iso_fortran_env, only: real128
    use constants
    implicit none
    private
    public :: jump_max, jump_count, profile_count
    public :: jump_radius, jump_factor, jump_width
    public :: density_profile, density_moment, set_density_profile, uniform_density
    public :: secondary_event_count, secondary_event_window, secondary_branch_density, secondary_knot

    integer, parameter :: jump_max = 8
    integer, parameter :: jump_slot = 28
    integer, parameter :: profile_slot = jump_slot+1+3*jump_max
    integer :: jump_count = 0, profile_count = 0, profile_step = 0
    integer :: profile_event_count = 0
    integer :: knot_count = 0
    logical :: density_ready = .false.
    real(8) :: density_rr=0d0,density_r0=0d0,density_value=0d0
    real(8), dimension(jump_max) :: profile_event_start=0d0, profile_event_end=0d0
    real(8), dimension(jump_max) :: profile_event_base=0d0
    real(8), dimension(jump_max) :: jump_radius= 0d0
    real(8), dimension(jump_max) :: jump_factor= 1d0
    real(8), dimension(jump_max) :: jump_width= 1d0
    real(8), allocatable, dimension(:) :: profile_radius
    real(8), allocatable, dimension(:) :: profile_knots
    real(8), allocatable, dimension(:) :: profile_logr, profile_logn, profile_logw, profile_power
    !$omp threadprivate(jump_count,profile_count,profile_step,jump_radius,jump_factor,jump_width, &
    !$omp& profile_radius,profile_logr,profile_logn,profile_logw,profile_power, &
    !$omp& profile_event_count,profile_event_start,profile_event_end,profile_event_base, &
    !$omp& knot_count,profile_knots,density_ready,density_rr,density_r0,density_value)

contains

pure logical function uniform_density(A_star, f_jump)
    implicit none
    real(8), intent(in) :: A_star, f_jump

    uniform_density = A_star <= 0d0 .and. f_jump == 1d0 .and. &
                      jump_count == 0 .and. profile_count == 0
end function uniform_density

subroutine density_profile(A_star, dNe_ISM, RR, R0, apply_jump, R_tr, f_jump, f_wide, dNe)
    implicit none
    integer, intent(in) :: apply_jump
    integer :: i
    real(8), intent(in) :: A_star, dNe_ISM, RR, R0, R_tr, f_jump, f_wide
    real(8), intent(out) :: dNe
    real(8) :: dNe_base, dNe_wind, enhancement, width_cm

    if (profile_count > 0 .and. density_ready) then
        if (RR == density_rr .and. R0 == density_r0) then
            dNe = density_value
            return
        end if
    end if

    select case (profile_count)
    case (0)
        if (A_star > 0d0) then
            dNe_wind = A_star*3.0d35/RR**2
            if (dNe_wind <= dNe_ISM/4d0) then
                dNe = dNe_ISM
            else
                dNe = dNe_wind
            end if
        else
            dNe = dNe_ISM
        end if
    case default
        call tab_density(RR, dNe)
    end select

    if (apply_jump /= 0 .and. profile_count == 0) then
        select case (jump_count)
        case (0)
            if (f_jump /= 1d0) dNe=dNe*(1d0+(f_jump-1d0)* &
                exp(-((RR-R_tr)*(RR-R_tr))/(2d0*(f_wide*R_tr)*(f_wide*R_tr))))
        case default
            dNe_base = dNe
            enhancement = 1d0
            do i = 1, jump_count
                if (jump_factor(i) == 1d0) cycle
                width_cm = jump_width(i)*jump_radius(i)
                enhancement = enhancement+(jump_factor(i)-1d0)* &
                              exp(-((RR-jump_radius(i))*(RR-jump_radius(i)))/ &
                              (2d0*width_cm*width_cm))
            end do
            dNe = dNe_base*enhancement
        end select
    end if

    if (RR < R0) then
        select case (profile_count)
        case (0)
            if (A_star > 0d0) then
                dNe = A_star*3.0d35/R0**2
            else
                dNe = dNe_ISM
            end if
        case default
            call tab_density(R0, dNe)
        end select
    end if
    if (profile_count > 0) then
        density_rr = RR
        density_r0 = R0
        density_value = dNe
        density_ready = .true.
    end if
end subroutine density_profile

! 计算真实外介质的径向扫掠矩 integral_0^R r^2 n(r)dr。
! Compute the radial swept-density moment integral_0^R r^2 n(r)dr.
subroutine density_moment(A_star,dNe_ISM,R,R0,R_tr,f_jump,f_wide,moment)
    implicit none
    integer :: j,njump
    real(8), intent(in) :: A_star,dNe_ISM,R,R0,R_tr,f_jump,f_wide
    real(8), intent(out) :: moment
    real(8) :: a,b,wind_norm,r_floor,center,factor,width,g0,g2
    logical :: wind_piece

    if (profile_count > 0) then
        call tab_moment(R,R0,moment)
        return
    end if

    if (A_star <= 0d0) then
        moment=dNe_ISM*R**3/3d0
        a=max(R0,0d0)
        if (R > a) then
            wind_piece=.false.
            call add_jumps(a,R,dNe_ISM)
        end if
        return
    end if

    wind_norm=A_star*3d35
    r_floor=sqrt(4d0*wind_norm/dNe_ISM)
    moment=0d0
    a=0d0
    if (R0 > 0d0) then
        b=min(R,R0)
        moment=wind_norm/R0**2*b**3/3d0
        if (R <= R0) return
        a=R0
    end if

    b=min(R,r_floor)
    if (b > a) then
        moment=moment+wind_norm*(b-a)
        wind_piece=.true.
        call add_jumps(a,b,wind_norm)
    end if
    a=max(a,r_floor)
    if (R > a) then
        moment=moment+dNe_ISM*(R**3-a**3)/3d0
        wind_piece=.false.
        call add_jumps(a,R,dNe_ISM)
    end if

contains

subroutine add_jumps(r_left,r_right,norm)
    implicit none
    real(8), intent(in) :: r_left,r_right,norm

    njump=max(1,jump_count)
    do j=1,njump
        if (jump_count == 0) then
            center=R_tr
            factor=f_jump
            width=f_wide*R_tr
        else
            center=jump_radius(j)
            factor=jump_factor(j)
            width=jump_width(j)*center
        end if
        if (factor == 1d0) cycle
        call gauss_moments(r_left,r_right,center,width,g0,g2)
        if (wind_piece) then
            moment=moment+norm*(factor-1d0)*g0
        else
            moment=moment+norm*(factor-1d0)*g2
        end if
    end do
end subroutine add_jumps

end subroutine density_moment

subroutine set_density_profile(Boundary, n)
    implicit none
    integer, intent(in) :: n
    integer :: i, radius_index, factor_index, width_index
    real(8), intent(in), dimension(n) :: Boundary
    real(8), allocatable, dimension(:) :: profile_density
    real(real128) :: lr0,lr1,ln0,ln1,ratn,ratr

    jump_count = 0
    profile_count = 0
    profile_step = 0
    profile_event_count = 0
    knot_count = 0
    density_ready = .false.
    profile_event_start = 0d0
    profile_event_end = 0d0
    profile_event_base = 0d0
    if (allocated(profile_radius)) &
        deallocate(profile_radius,profile_logr,profile_logn,profile_logw,profile_power)
    if (allocated(profile_knots)) deallocate(profile_knots)
    if (n < jump_slot) return
    jump_count = nint(Boundary(jump_slot))
    if (jump_count < 0 .or. jump_count > jump_max) &
        error stop 'density jump count outside supported range'
    radius_index = jump_slot+1
    factor_index = radius_index+jump_max
    width_index = factor_index+jump_max
    if (jump_count > 0 .and. n < width_index+jump_max-1) &
        error stop 'boundary is missing density jump arrays'
    do i = 1, jump_count
        jump_radius(i) = Boundary(radius_index+i-1)
        jump_factor(i) = Boundary(factor_index+i-1)
        jump_width(i) = Boundary(width_index+i-1)
        if (jump_radius(i) <= 0d0 .or. jump_factor(i) <= 0d0 .or. &
            jump_width(i) <= 0d0) &
            error stop 'density jump radii, factors, and widths must be positive'
    end do
    if (n < profile_slot) return
    profile_count = nint(Boundary(profile_slot))
    if (profile_count < 0) &
        error stop 'density profile point count outside supported range'
    if (profile_count > 0 .and. (jump_count > 0 .or. Boundary(22) /= 1d0)) &
        error stop 'density profile and density jumps are mutually exclusive'
    if (profile_count == 1) &
        error stop 'density profile requires at least 2 points'
    radius_index = profile_slot+1
    factor_index = radius_index+profile_count
    if (profile_count > 0 .and. n < factor_index+profile_count-1) &
        error stop 'boundary is missing density profile arrays'
    if (profile_count == 0) return
    profile_step = 1
    do while (2*profile_step < profile_count)
        profile_step = 2*profile_step
    end do
    allocate(profile_radius(profile_count),profile_logr(profile_count),profile_logn(profile_count), &
             profile_logw(profile_count),profile_power(profile_count-1),profile_density(profile_count))
    do i = 1, profile_count
        profile_radius(i) = Boundary(radius_index+i-1)
        profile_density(i) = Boundary(factor_index+i-1)
        if (profile_radius(i) <= 0d0 .or. profile_density(i) <= 0d0) &
            error stop 'density profile radii and densities must be positive'
        lr0=log(real(profile_radius(i),real128))
        ln0=log(real(profile_density(i),real128))
        profile_logr(i)=real(lr0,kind(1d0))
        profile_logn(i)=real(ln0,kind(1d0))
        profile_logw(i)=real(ln0+3*lr0,kind(1d0))
    end do
    do i = 2, profile_count
        if (profile_radius(i) <= profile_radius(i-1)) &
            error stop 'density profile radii must be strictly increasing'
        lr0=log(real(profile_radius(i-1),real128))
        lr1=log(real(profile_radius(i),real128))
        ln0=log(real(profile_density(i-1),real128))
        ln1=log(real(profile_density(i),real128))
        ratn=real(profile_density(i),real128)/real(profile_density(i-1),real128)
        ratr=real(profile_radius(i-1),real128)/real(profile_radius(i),real128)
        if (ratn == ratr**3) then
            profile_power(i-1)=0d0
        else
            profile_power(i-1)=real(((ln1+3*lr1)-(ln0+3*lr0))/(lr1-lr0),kind(1d0))
        end if
    end do
    call build_profile_events()
end subroutine set_density_profile

subroutine build_profile_events()
    implicit none
    integer :: i, i_start, i_end
    logical :: rising

    ! q=d ln(r^3 n)/d ln r; q>3 is exactly d n/dR>0 for the log-log table.
    ! Consecutive rising intervals form one finite-width compression source.
    profile_event_count = 0
    if (profile_count < 2) return
    allocate(profile_knots(profile_count))
    knot_count = 0
    do i = 1, profile_count
        rising = .false.
        if (i < profile_count) rising = profile_power(i) > 3d0
        if (i > 1) rising = rising .or. profile_power(i-1) > 3d0
        if (rising) then
            knot_count = knot_count+1
            profile_knots(knot_count) = profile_radius(i)
        end if
    end do
    i = 1
    do while (i < profile_count)
        if (profile_power(i) <= 3d0) then
            i = i + 1
            cycle
        end if
        i_start = i
        do while (i < profile_count .and. profile_power(i) > 3d0)
            i = i + 1
        end do
        i_end = i
        profile_event_count = profile_event_count + 1
        if (profile_event_count > jump_max) &
            error stop 'tabulated density profile has more than 8 compression intervals'
        profile_event_start(profile_event_count) = profile_radius(i_start)
        profile_event_end(profile_event_count) = profile_radius(i_end)
        profile_event_base(profile_event_count) = exp(profile_logn(i_start))
    end do
end subroutine build_profile_events

integer function secondary_event_count()
    implicit none

    secondary_event_count = jump_count + profile_event_count
end function secondary_event_count

subroutine secondary_event_window(j, r_left, r_right, width, center)
    implicit none
    integer, intent(in) :: j
    real(8), intent(out) :: r_left, r_right, width, center
    integer :: jp

    if (j <= jump_count) then
        center = jump_radius(j)
        width = jump_width(j)*center
        r_left = center-4d0*width
        r_right = center
    else
        jp = j-jump_count
        r_left = profile_event_start(jp)
        r_right = profile_event_end(jp)
        center = sqrt(r_left*r_right)
        width = 0.5d0*(r_right-r_left)
    end if
end subroutine secondary_event_window

! 返回两个半径之间最先遇到的 tabulated 上升段 knot。
! Return the first tabulated rising-segment knot between two radii.
subroutine secondary_knot(r_left,r_right,knot,found)
    implicit none
    integer :: lo,hi,mid
    real(8), intent(in) :: r_left,r_right
    real(8), intent(out) :: knot
    logical, intent(out) :: found
    real(8) :: lower,upper

    knot=0d0
    found=.false.
    if (knot_count == 0 .or. r_right <= r_left) return
    lower=r_left*(1d0+1d-12)
    upper=r_right*(1d0-1d-12)
    lo=0
    hi=knot_count+1
    do while (hi-lo > 1)
        mid=(lo+hi)/2
        if (profile_knots(mid) <= lower) then
            lo=mid
        else
            hi=mid
        end if
    end do
    if (hi > knot_count .or. profile_knots(hi) >= upper) return
    knot=profile_knots(hi)
    found=.true.
end subroutine secondary_knot

subroutine secondary_branch_density(A_star,dNe_ISM,RR,R0,apply_jump,R_tr,f_jump,f_wide,j, &
                                    dens_all,dens_bump,dens_known)
    implicit none
    integer, intent(in) :: apply_jump, j
    real(8), intent(in) :: A_star,dNe_ISM,RR,R0,R_tr,f_jump,f_wide
    real(8), intent(in), optional :: dens_known
    real(8), intent(out) :: dens_all,dens_bump
    integer :: k,jp
    real(8) :: x,width,prof,enh,nbase

    ! The tabulated branch uses the exact local profile minus its pre-rise state;
    ! no Gaussian replacement or post-processing is applied.
    if (present(dens_known)) then
        dens_all=dens_known
    else
        call density_profile(A_star,dNe_ISM,RR,R0,apply_jump,R_tr,f_jump,f_wide,dens_all)
    end if
    dens_bump = 0d0
    if (j > jump_count) then
        jp = j-jump_count
        if (RR >= profile_event_start(jp) .and. RR < profile_event_end(jp)) &
            dens_bump = max(0d0,dens_all-profile_event_base(jp))
        return
    end if

    enh = 1d0
    do k=1,jump_count
        x=RR-jump_radius(k)
        width=jump_width(k)*jump_radius(k)
        prof=(jump_factor(k)-1d0)*exp(-(x*x)/(2d0*width*width))
        enh=enh+prof
    end do
    nbase=dens_all/enh
    x=RR-jump_radius(j)
    width=jump_width(j)*jump_radius(j)
    prof=(jump_factor(j)-1d0)*exp(-(x*x)/(2d0*width*width))
    if (x >= -4d0*width .and. x < 0d0) dens_bump=nbase*prof
end subroutine secondary_branch_density

subroutine tab_moment(R,R0,moment)
    implicit none
    integer :: i
    real(8), intent(in) :: R,R0
    real(8), intent(out) :: moment
    real(8) :: a,b,lo,hi,q,span,x,numer,tval,delta,logw

    moment=0d0
    a=0d0
    if (R0 > 0d0) then
        b=min(R,R0)
        moment=exp(tab_logdensity(R0)+3d0*log(b)-log(3d0))
        if (R <= R0) return
        a=R0
    else if (profile_power(1) <= 0d0) then
        error stop 'density profile swept moment diverges at the origin'
    end if

    do i=1,profile_count-1
        lo=a
        if (i > 1) lo=max(lo,profile_radius(i))
        hi=R
        if (i < profile_count-1) hi=min(hi,profile_radius(i+1))
        if (hi <= lo) cycle
        q=profile_power(i)
        if (lo == 0d0) then
            delta=log(hi)-profile_logr(i)
            logw=profile_logw(i)+q*delta
            moment=moment+exp(logw-log(q))
        else
            span=log(hi/lo)
            if (q == 0d0) then
                logw=profile_logw(i)
                moment=moment+exp(logw+log(span))
            else
                x=q*span
                tval=tanh(x/2d0)
                if (x < 0d0) then
                    delta=log(lo)-profile_logr(i)
                    numer=2d0*tval/(1d0-tval)
                else
                    delta=log(hi)-profile_logr(i)
                    numer=2d0*tval/(1d0+tval)
                end if
                logw=profile_logw(i)+q*delta
                moment=moment+exp(logw+log(abs(numer))-log(abs(q)))
            end if
        end if
    end do
end subroutine tab_moment

pure subroutine gauss_moments(a,b,center,width,g0,g2)
    implicit none
    real(8), intent(in) :: a,b,center,width
    real(8), intent(out) :: g0,g2
    real(8) :: ua,ub,ea,eb,derf

    ua=(a-center)/(sqrt(2d0)*width)
    ub=(b-center)/(sqrt(2d0)*width)
    ea=exp(-ua**2)
    eb=exp(-ub**2)
    if (ub <= 0d0) then
        derf=erfc(-ub)-erfc(-ua)
    else if (ua >= 0d0) then
        derf=erfc(ua)-erfc(ub)
    else
        derf=erf(ub)-erf(ua)
    end if
    g0=width*sqrt(pi/2d0)*derf
    g2=(center**2+width**2)*g0-width**2*((b+center)*eb-(a+center)*ea)
end subroutine gauss_moments

function tab_logdensity(RR) result(logn)
    implicit none
    integer :: lo, mid, step
    real(8), intent(in) :: RR
    real(8) :: logn,delta

    lo=1
    step=profile_step
    do while (step > 0)
        mid=min(lo+step,profile_count-1)
        if (profile_radius(mid) <= RR) lo=mid
        step=step/2
    end do
    delta=log(RR)-profile_logr(lo)
    logn=profile_logn(lo)+profile_power(lo)*delta-3d0*delta
end function tab_logdensity

subroutine tab_density(RR, dNe)
    implicit none
    real(8), intent(in) :: RR
    real(8), intent(out) :: dNe

    dNe=exp(tab_logdensity(RR))
end subroutine tab_density

end module dynamics_density_profile
