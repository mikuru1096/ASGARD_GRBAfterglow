module dynamics_density_profile
    use iso_fortran_env, only: real128
    use constants
    implicit none
    private
    public :: jump_max, jump_count, profile_count
    public :: jump_radius, jump_factor, jump_width
    public :: density_profile, density_moment, set_density_profile, uniform_density

    integer, parameter :: jump_max = 8, profile_max = 96
    integer, parameter :: jump_slot = 28
    integer, parameter :: profile_slot = jump_slot+1+3*jump_max
    integer :: jump_count = 0, profile_count = 0
    real(8), dimension(jump_max) :: jump_radius= 0d0
    real(8), dimension(jump_max) :: jump_factor= 1d0
    real(8), dimension(jump_max) :: jump_width= 1d0
    real(8), dimension(profile_max) :: profile_radius= 0d0
    real(8), dimension(profile_max) :: profile_density= 0d0
    real(8), dimension(profile_max) :: profile_logr= 0d0, profile_logn= 0d0, profile_logw= 0d0
    real(8), dimension(profile_max-1) :: profile_power= 0d0

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
    real(real128) :: lr0,lr1,ln0,ln1,ratn,ratr

    jump_count = 0
    profile_count = 0
    jump_radius = 0d0
    jump_factor = 1d0
    jump_width = 1d0
    profile_radius = 0d0
    profile_density = 0d0
    profile_logr = 0d0
    profile_logn = 0d0
    profile_logw = 0d0
    profile_power = 0d0
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
    if (profile_count < 0 .or. profile_count > profile_max) &
        error stop 'density profile point count outside supported range'
    if (profile_count > 0 .and. (jump_count > 0 .or. Boundary(22) /= 1d0)) &
        error stop 'density profile and density jumps are mutually exclusive'
    if (profile_count == 1) &
        error stop 'density profile requires at least 2d0 points'
    radius_index = profile_slot+1
    factor_index = radius_index+profile_max
    if (profile_count > 0 .and. n < factor_index+profile_max-1) &
        error stop 'boundary is missing density profile arrays'
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
end subroutine set_density_profile

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
    integer :: lo, hi, mid
    real(8), intent(in) :: RR
    real(8) :: logn,delta

    if (RR <= profile_radius(1)) then
        lo=1
    else if (RR >= profile_radius(profile_count)) then
        lo=profile_count-1
    else
        lo = 1; hi = profile_count
        do while (hi-lo > 1)
            mid = (lo+hi)/2
            if (profile_radius(mid) <= RR) then
                lo = mid
            else
                hi = mid
            end if
        end do
    end if
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
