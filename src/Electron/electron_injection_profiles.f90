!f2py: skip
! 电子注入谱工具：统一处理快/慢冷却初始谱、FS源项、RS动能源项和热电子注入。
! Electron injection-profile utilities: fast/slow-cooling spectra, FS sources, RS kinetic sources,
! and thermal injection.
module electron_injection_profiles
    use constants
    use electron_coord_common, only: coord_fourvel, coord_from_xg, &
                                                 xg_from_coord, gamma_from_coord, &
                                                 dxg_dcoord
    use electron_radiation_kernel, only: besselk
    implicit none
    private

    public :: exp_cutoff, dnx_cutoff, pl_params
    public :: log_edges
    public :: init_powerlaw, init_edges
    public :: init_coord
    public :: source_edges, source_coord
    public :: kinetic_edges, kinetic_coord
    public :: add_thermal, thermal_pop

contains

! 指数截断因子：γ > gmax 时返回 exp(1-γ/gmax)，否则返回 1。
! Exponential cutoff factor: return exp(1-gamma/gmax) above gmax, otherwise 1.
pure real(8) function exp_cutoff(gam,Gam_e_max)
    implicit none
    real(8), intent(in) :: gam,Gam_e_max

    exp_cutoff=1d0
    if (gam > Gam_e_max) exp_cutoff=dexp(1d0-gam/Gam_e_max)
end function exp_cutoff

! 根据快冷却/慢冷却分支返回幂律注入参数：系数 coeff 和谱指数 slope。
! Return power-law injection parameters, coefficient and slope, for fast/slow cooling branches.
subroutine pl_params(Para_N_e_ini,p,Gam_e_m,Gam_e_c,gam,active,coeff,slope)
    implicit none
    logical, intent(out) :: active
    real(8), intent(in) :: Para_N_e_ini,p,Gam_e_m,Gam_e_c,gam
    real(8), intent(out) :: coeff,slope

    active=.false.
    coeff=0d0
    slope=0d0
    if (Gam_e_m > Gam_e_c) then
        if (gam < Gam_e_c) return
        active=.true.
        coeff=Para_N_e_ini*Gam_e_c
        slope=merge(2d0,p+1d0,gam < Gam_e_m)
        if (gam >= Gam_e_m) coeff=coeff*Gam_e_m**(p-1d0)
    else
        if (gam < Gam_e_m) return
        active=.true.
        coeff=Para_N_e_ini*(p-1d0)*Gam_e_m**(p-1d0)
        slope=merge(p,p+1d0,gam < Gam_e_c)
        if (gam >= Gam_e_c) coeff=coeff*Gam_e_c
    end if
end subroutine pl_params

! 幂律+指数截断的 dN/dx 值：dN/dx = coeff * ln(10) * gamma^(1-slope) * cutoff(gamma)。
! Power law with exponential cutoff in dN/dx form.
real(8) function dnx_cutoff(x,coeff,slope,Gam_e_max)
    implicit none
    real(8), intent(in) :: x,coeff,slope,Gam_e_max
    real(8) :: gam,cutoff_factor

    gam=dexp(x)
    if (gam <= 0d0 .or. coeff <= 0d0) then
        dnx_cutoff=0d0
        return
    end if

    cutoff_factor=1d0
    if (Gam_e_max > 0d0 .and. gam > Gam_e_max) cutoff_factor=dexp(1d0-gam/Gam_e_max)
    dnx_cutoff=coeff*gam**(1d0-slope)*cutoff_factor
end function dnx_cutoff

! 3点 Gauss-Legendre 积分：在 [xlo, xhi] 上对 dN/dx 积分。
! Three-point Gauss-Legendre integral of dN/dx over [xlo, xhi].
real(8) function dnx_gauss3(coeff,slope,Gam_e_max,x_lo,x_hi)
    implicit none
    integer :: I_q
    real(8), intent(in) :: coeff,slope,Gam_e_max,x_lo,x_hi
    real(8), parameter, dimension(3) :: xi=[-dsqrt(3d0/5d0),0d0,dsqrt(3d0/5d0)]
    real(8), parameter, dimension(3) :: wi=[5d0/9d0,8d0/9d0,5d0/9d0]
    real(8) :: half_dx,x_mid,x_eval,quad

    if (x_hi <= x_lo) then
        dnx_gauss3=0d0
        return
    end if

    half_dx=0.5d0*(x_hi-x_lo)
    x_mid=0.5d0*(x_hi+x_lo)
    quad=0d0
    do I_q=1,3
        x_eval=x_mid+half_dx*xi(I_q)
        quad=quad+wi(I_q)*dnx_cutoff(x_eval,coeff,slope,Gam_e_max)
    end do
    dnx_gauss3=half_dx*quad
end function dnx_gauss3

! 3点 Gauss-Legendre 积分：在四速度坐标 y 上对 dN/dy=(dN/dx)(dx/dy) 积分。
! Three-point Gauss-Legendre integral in four-velocity coordinate y.
real(8) function dny_gauss3(coord_scale,coeff,slope,Gam_e_max,y_lo,y_hi)
    implicit none
    integer :: I_q
    real(8), intent(in) :: coord_scale,coeff,slope,Gam_e_max,y_lo,y_hi
    real(8), parameter, dimension(3) :: xi=[-dsqrt(3d0/5d0),0d0,dsqrt(3d0/5d0)]
    real(8), parameter, dimension(3) :: wi=[5d0/9d0,8d0/9d0,5d0/9d0]
    real(8) :: half_dy,y_mid,y_eval,x_eval,quad

    if (y_hi <= y_lo) then
        dny_gauss3=0d0
        return
    end if

    half_dy=0.5d0*(y_hi-y_lo)
    y_mid=0.5d0*(y_hi+y_lo)
    quad=0d0
    do I_q=1,3
        y_eval=y_mid+half_dy*xi(I_q)
        x_eval=xg_from_coord(coord_fourvel,coord_scale,y_eval)
        quad=quad+wi(I_q)*dnx_cutoff(x_eval,coeff,slope,Gam_e_max) &
             *dxg_dcoord(coord_fourvel,coord_scale,y_eval)
    end do
    dny_gauss3=half_dy*quad
end function dny_gauss3

! 将激活区间上的 dN/dx 积分累加到 acc，在 gmax 处自动分段处理截断。
! Accumulate active dN/dx support into acc, split at gmax for the cutoff kink.
subroutine dnx_segment(cell_lo,cell_hi,active_lo,active_hi,coeff,slope,Gam_e_max,acc)
    implicit none
    real(8), intent(in) :: cell_lo,cell_hi,active_lo,active_hi,coeff,slope,Gam_e_max
    real(8), intent(inout) :: acc
    real(8) :: x_lo,x_hi,x_cut

    if (coeff <= 0d0 .or. cell_hi <= cell_lo .or. active_hi <= active_lo) return
    x_lo=max(cell_lo,active_lo)
    x_hi=min(cell_hi,active_hi)
    if (x_hi <= x_lo) return

    x_cut=dlog(Gam_e_max)
    if (x_lo < x_cut .and. x_hi > x_cut) then
        acc=acc+dnx_gauss3(coeff,slope,Gam_e_max,x_lo,x_cut)
        acc=acc+dnx_gauss3(coeff,slope,Gam_e_max,x_cut,x_hi)
    else
        acc=acc+dnx_gauss3(coeff,slope,Gam_e_max,x_lo,x_hi)
    end if
end subroutine dnx_segment

! 将激活区间上的 dN/dy 积分累加到 acc，在 gmax 处自动分段处理截断。
! Accumulate active dN/dy support into acc, split at gmax for the cutoff kink.
subroutine dny_segment(cell_lo,cell_hi,active_lo,active_hi,coord_scale,coeff,slope,Gam_e_max,acc)
    implicit none
    real(8), intent(in) :: cell_lo,cell_hi,active_lo,active_hi,coord_scale,coeff,slope,Gam_e_max
    real(8), intent(inout) :: acc
    real(8) :: y_lo,y_hi,y_cut

    if (coeff <= 0d0 .or. cell_hi <= cell_lo .or. active_hi <= active_lo) return
    y_lo=max(cell_lo,active_lo)
    y_hi=min(cell_hi,active_hi)
    if (y_hi <= y_lo) return

    y_cut=coord_from_xg(coord_fourvel,coord_scale,dlog(Gam_e_max))
    if (y_lo < y_cut .and. y_hi > y_cut) then
        acc=acc+dny_gauss3(coord_scale,coeff,slope,Gam_e_max,y_lo,y_cut)
        acc=acc+dny_gauss3(coord_scale,coeff,slope,Gam_e_max,y_cut,y_hi)
    else
        acc=acc+dny_gauss3(coord_scale,coeff,slope,Gam_e_max,y_lo,y_hi)
    end if
end subroutine dny_segment

! 由网格中心值推导 ln(gamma) 单元边界。
! Build ln(gamma) cell edges from grid-center values.
subroutine log_edges(Num_gam_e,gam_e,x_edge)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_gam_e+1) :: x_edge

    x_edge(1)=dlog(gam_e(1))-0.5d0*(dlog(gam_e(2))-dlog(gam_e(1)))
    do I_gam_e=2,Num_gam_e
        x_edge(I_gam_e)=0.5d0*(dlog(gam_e(I_gam_e-1))+dlog(gam_e(I_gam_e)))
    end do
    x_edge(Num_gam_e+1)=dlog(gam_e(Num_gam_e))+0.5d0*(dlog(gam_e(Num_gam_e))-dlog(gam_e(Num_gam_e-1)))
end subroutine log_edges

! 生成快/慢冷却幂律+指数截断的初始电子谱 dN/dgamma（网格中心值）。
! Build the fast/slow-cooling initial electron spectrum dN/dgamma at grid centers.
subroutine init_powerlaw(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,dN_gam_e_1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in), dimension(Num_gam_e) :: gam_e
    real(8), intent(in) :: Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max
    real(8), intent(out), dimension(Num_gam_e) :: dN_gam_e_1
    logical :: active
    real(8) :: coeff,slope

    dN_gam_e_1=0d0
    if (Gam_e_max <= 0d0) return

    do I_gam_e=1,Num_gam_e
        call pl_params(Para_N_e_ini,p,Gam_e_m,Gam_e_c,gam_e(I_gam_e),active,coeff,slope)
        if (active) dN_gam_e_1(I_gam_e)=coeff*gam_e(I_gam_e)**(-slope)*exp_cutoff(gam_e(I_gam_e),Gam_e_max)
    end do
end subroutine init_powerlaw

! 生成快/慢冷却幂律+指数截断的初始电子谱 dN/dx（网格单元平均，保正/守恒）。
! Build the fast/slow-cooling initial dN/dx as conservative positive cell averages.
subroutine init_edges(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,x_edge,dN_x_1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in), dimension(Num_gam_e+1) :: x_edge
    real(8), intent(in) :: Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max
    real(8), intent(out), dimension(Num_gam_e) :: dN_x_1
    real(8) :: cell_lo,cell_hi,dx_cell,seg_sum,x_m,x_c,coeff_lo,coeff_hi,huge_x

    dN_x_1=0d0
    if (Gam_e_max <= 0d0) return

    x_m=dlog(Gam_e_m)
    x_c=dlog(Gam_e_c)
    huge_x=1d300

    do I_gam_e=1,Num_gam_e
        cell_lo=x_edge(I_gam_e)
        cell_hi=x_edge(I_gam_e+1)
        dx_cell=cell_hi-cell_lo
        if (dx_cell <= 0d0) cycle

        seg_sum=0d0
        if (Gam_e_m > Gam_e_c) then
            coeff_lo=Para_N_e_ini*Gam_e_c
            coeff_hi=Para_N_e_ini*Gam_e_c*Gam_e_m**(p-1d0)
            call dnx_segment(cell_lo,cell_hi,x_c,x_m,coeff_lo,2d0,Gam_e_max,seg_sum)
            call dnx_segment(cell_lo,cell_hi,x_m,huge_x,coeff_hi,p+1d0,Gam_e_max,seg_sum)
        else
            coeff_lo=Para_N_e_ini*(p-1d0)*Gam_e_m**(p-1d0)
            coeff_hi=coeff_lo*Gam_e_c
            call dnx_segment(cell_lo,cell_hi,x_m,x_c,coeff_lo,p,Gam_e_max,seg_sum)
            call dnx_segment(cell_lo,cell_hi,x_c,huge_x,coeff_hi,p+1d0,Gam_e_max,seg_sum)
        end if
        dN_x_1(I_gam_e)=seg_sum/dx_cell
    end do
end subroutine init_edges

! 生成快/慢冷却幂律+指数截断的初始电子谱 dN/dy，y 为四速度坐标。
! Build the initial dN/dy in four-velocity coordinate y.
subroutine init_coord(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                                            Num_gam_e,coord_edge,coord_scale,dN_y_1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in), dimension(Num_gam_e+1) :: coord_edge
    real(8), intent(in) :: Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,coord_scale
    real(8), intent(out), dimension(Num_gam_e) :: dN_y_1
    real(8) :: cell_lo,cell_hi,dy_cell,seg_sum,y_m,y_c,coeff_lo,coeff_hi,huge_y

    dN_y_1=0d0
    if (Gam_e_max <= 0d0) return

    y_m=coord_from_xg(coord_fourvel,coord_scale,dlog(Gam_e_m))
    y_c=coord_from_xg(coord_fourvel,coord_scale,dlog(Gam_e_c))
    huge_y=1d300

    do I_gam_e=1,Num_gam_e
        cell_lo=coord_edge(I_gam_e)
        cell_hi=coord_edge(I_gam_e+1)
        dy_cell=cell_hi-cell_lo
        if (dy_cell <= 0d0) cycle

        seg_sum=0d0
        if (Gam_e_m > Gam_e_c) then
            coeff_lo=Para_N_e_ini*Gam_e_c
            coeff_hi=Para_N_e_ini*Gam_e_c*Gam_e_m**(p-1d0)
            call dny_segment(cell_lo,cell_hi,y_c,y_m,coord_scale,coeff_lo,2d0,Gam_e_max,seg_sum)
            call dny_segment(cell_lo,cell_hi,y_m,huge_y,coord_scale,coeff_hi,p+1d0,Gam_e_max,seg_sum)
        else
            coeff_lo=Para_N_e_ini*(p-1d0)*Gam_e_m**(p-1d0)
            coeff_hi=coeff_lo*Gam_e_c
            call dny_segment(cell_lo,cell_hi,y_m,y_c,coord_scale,coeff_lo,p,Gam_e_max,seg_sum)
            call dny_segment(cell_lo,cell_hi,y_c,huge_y,coord_scale,coeff_hi,p+1d0,Gam_e_max,seg_sum)
        end if
        dN_y_1(I_gam_e)=seg_sum/dy_cell
    end do
end subroutine init_coord

! 构建幂律+指数截断源项 dF/dx（网格单元平均，保正/守恒）。
! Build the power-law cutoff source dF/dx as conservative positive cell averages.
subroutine source_edges(Num_gam_e,x_edge,Gam_e_m,Gam_e_max,Q,p,dF1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in), dimension(Num_gam_e+1) :: x_edge
    real(8), intent(in) :: Gam_e_m,Gam_e_max,Q,p
    real(8), intent(out), dimension(Num_gam_e) :: dF1
    real(8) :: cell_lo,cell_hi,dx_cell,seg_sum,x_m,huge_x

    dF1=0d0
    if (Gam_e_max <= 0d0 .or. Q <= 0d0) return

    x_m=dlog(Gam_e_m)
    huge_x=1d300
    do I_gam_e=1,Num_gam_e
        cell_lo=x_edge(I_gam_e)
        cell_hi=x_edge(I_gam_e+1)
        dx_cell=cell_hi-cell_lo
        if (dx_cell <= 0d0) cycle
        seg_sum=0d0
        call dnx_segment(cell_lo,cell_hi,x_m,huge_x,Q,p,Gam_e_max,seg_sum)
        dF1(I_gam_e)=seg_sum/dx_cell
    end do
end subroutine source_edges

! 构建 FS 幂律+指数截断源项 dF/dy，y 为四速度坐标，物理谱仍为 gamma^(-p)。
! Build the FS source dF/dy in four-velocity coordinate y; the physical spectrum remains gamma^(-p).
subroutine source_coord(Num_gam_e,coord_edge,coord_scale,Gam_e_m,Gam_e_max,Q,p,dF1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in), dimension(Num_gam_e+1) :: coord_edge
    real(8), intent(in) :: coord_scale,Gam_e_m,Gam_e_max,Q,p
    real(8), intent(out), dimension(Num_gam_e) :: dF1
    real(8) :: cell_lo,cell_hi,dy_cell,seg_sum,y_m,huge_y

    dF1=0d0
    if (Gam_e_max <= 0d0 .or. Q <= 0d0) return

    y_m=coord_from_xg(coord_fourvel,coord_scale,dlog(Gam_e_m))
    huge_y=1d300
    do I_gam_e=1,Num_gam_e
        cell_lo=coord_edge(I_gam_e)
        cell_hi=coord_edge(I_gam_e+1)
        dy_cell=cell_hi-cell_lo
        if (dy_cell <= 0d0) cycle
        seg_sum=0d0
        call dny_segment(cell_lo,cell_hi,y_m,huge_y,coord_scale,Q,p,Gam_e_max,seg_sum)
        dF1(I_gam_e)=seg_sum/dy_cell
    end do
end subroutine source_coord

! 构建以动能幂律归一的反向激波源项：dN/dx 正比 gamma*(gamma-1)^(-p)*cutoff。
! Build the reverse-shock kinetic source normalized from gamma*(gamma-1)^(-p)*cutoff.
subroutine kinetic_edges(Num_gam_e,x_edge,Gam_e_m,Gam_e_max,Q,p,dF1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e,I_q
    real(8), intent(in), dimension(Num_gam_e+1) :: x_edge
    real(8), intent(in) :: Gam_e_m,Gam_e_max,Q,p
    real(8), intent(out), dimension(Num_gam_e) :: dF1
    real(8), parameter, dimension(3) :: xi=[-dsqrt(3d0/5d0),0d0,dsqrt(3d0/5d0)]
    real(8), parameter, dimension(3) :: wi=[5d0/9d0,8d0/9d0,5d0/9d0]
    real(8) :: cell_lo,cell_hi,dx_cell,half_dx,x_mid,x_eval,gam,cutoff_factor,cell_sum,shape_norm

    dF1=0d0
    shape_norm=0d0
    if (Gam_e_max <= 0d0 .or. Q <= 0d0) return

    do I_gam_e=1,Num_gam_e
        cell_lo=x_edge(I_gam_e)
        cell_hi=x_edge(I_gam_e+1)
        dx_cell=cell_hi-cell_lo
        if (dx_cell <= 0d0) cycle
        half_dx=0.5d0*dx_cell
        x_mid=0.5d0*(cell_lo+cell_hi)
        cell_sum=0d0
        do I_q=1,3
            x_eval=x_mid+half_dx*xi(I_q)
            gam=dexp(x_eval)
            if (gam > Gam_e_m) then
                cutoff_factor=exp_cutoff(gam,Gam_e_max)
                cell_sum=cell_sum+wi(I_q)*gam*(gam-1d0)**(-p)*cutoff_factor
            end if
        end do
        cell_sum=half_dx*cell_sum
        dF1(I_gam_e)=cell_sum/dx_cell
        shape_norm=shape_norm+cell_sum
    end do
    if (shape_norm <= 0d0) error stop 'kinetic electron source has empty active support'
    dF1=Q*dF1/shape_norm
end subroutine kinetic_edges

! 构建以动能幂律归一的源项 dN/dy；y 为 log four-velocity-squared 坐标。
! Build the kinetic-normalized source dN/dy in log four-velocity-squared coordinate.
subroutine kinetic_coord(Num_gam_e,coord_edge,coord_scale,Gam_e_m,Gam_e_max,Q,p,dF1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e,I_q
    real(8), intent(in), dimension(Num_gam_e+1) :: coord_edge
    real(8), intent(in) :: coord_scale,Gam_e_m,Gam_e_max,Q,p
    real(8), intent(out), dimension(Num_gam_e) :: dF1
    real(8), parameter, dimension(3) :: xi=[-dsqrt(3d0/5d0),0d0,dsqrt(3d0/5d0)]
    real(8), parameter, dimension(3) :: wi=[5d0/9d0,8d0/9d0,5d0/9d0]
    real(8) :: cell_lo,cell_hi,dy_cell,half_dy,y_mid,y_eval,gam,dxdy,cell_sum,shape_norm

    dF1=0d0
    shape_norm=0d0
    if (Gam_e_max <= 0d0 .or. Q <= 0d0) return

    do I_gam_e=1,Num_gam_e
        cell_lo=coord_edge(I_gam_e)
        cell_hi=coord_edge(I_gam_e+1)
        dy_cell=cell_hi-cell_lo
        half_dy=0.5d0*dy_cell
        y_mid=0.5d0*(cell_lo+cell_hi)
        cell_sum=0d0
        do I_q=1,3
            y_eval=y_mid+half_dy*xi(I_q)
            gam=gamma_from_coord(coord_fourvel,coord_scale,y_eval)
            if (gam > Gam_e_m) then
                dxdy=dxg_dcoord(coord_fourvel,coord_scale,y_eval)
                cell_sum=cell_sum+wi(I_q)*gam*dxdy*(gam-1d0)**(-p)*exp_cutoff(gam,Gam_e_max)
            end if
        end do
        cell_sum=half_dy*cell_sum
        dF1(I_gam_e)=cell_sum/dy_cell
        shape_norm=shape_norm+cell_sum
    end do
    if (shape_norm <= 0d0) error stop 'kinetic electron coordinate source has empty active support'
    dF1=Q*dF1/shape_norm
end subroutine kinetic_coord

! 由 shock-front 四速度估算 Maxwell-Juttner 温度 theta。
! Estimate Maxwell-Juttner temperature theta from shock-front four velocity.
pure real(8) function thermal_theta(four_v)
    implicit none
    real(8), intent(in) :: four_v

    if (four_v <= 0d0) error stop 'thermal_theta requires four_v > 0'
    thermal_theta=four_v*(four_v+1.07d0*four_v*four_v)/(3d0*(1d0+four_v+1.07d0*four_v*four_v))
end function thermal_theta

! 构建单位归一的热电子 dN/dx 形状。
! Build a unit-normalized thermal-electron dN/dx shape.
subroutine thermal_shape(Num_gam_e,gam_e,theta,shape_dnx)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in), dimension(Num_gam_e) :: gam_e
    real(8), intent(in) :: theta
    real(8), intent(out), dimension(Num_gam_e) :: shape_dnx
    real(8), dimension(Num_gam_e) :: shape_dgam
    real(8) :: k2_theta,norm_dgam

    if (theta <= 0d0) error stop 'thermal_shape requires theta > 0'
    if (any(gam_e <= 1d0)) error stop 'thermal electron grid requires gamma > 1'

    k2_theta=besselk(1d0/theta)
    if (k2_theta <= 0d0) error stop 'besselk normalization vanished in thermal electron source'

    do I_gam_e=1,Num_gam_e
        shape_dgam(I_gam_e)=gam_e(I_gam_e)**2*dsqrt(1d0-1d0/gam_e(I_gam_e)**2) &
                          *dexp(-gam_e(I_gam_e)/theta)/(theta*k2_theta)
    end do
    norm_dgam=sum((shape_dgam(2:Num_gam_e)+shape_dgam(1:Num_gam_e-1)) &
            *(gam_e(2:Num_gam_e)-gam_e(1:Num_gam_e-1)))/2d0
    if (norm_dgam <= 0d0) error stop 'thermal electron distribution normalization is non-positive'

    shape_dnx=shape_dgam/norm_dgam*gam_e
end subroutine thermal_shape

! 将热电子源项加入已有 dF1。
! Add thermal-electron source counts into the existing dF1 source array.
subroutine add_thermal(Num_gam_e,gam_e,four_v,total_count,dF1)
    implicit none
    integer, intent(in) :: Num_gam_e
    real(8), intent(in), dimension(Num_gam_e) :: gam_e
    real(8), intent(in) :: four_v,total_count
    real(8), intent(inout), dimension(Num_gam_e) :: dF1
    real(8), dimension(Num_gam_e) :: shape_dnx
    real(8) :: theta

    if (total_count <= 0d0) return
    theta=thermal_theta(four_v)
    call thermal_shape(Num_gam_e,gam_e,theta,shape_dnx)
    dF1=dF1+total_count*shape_dnx
end subroutine add_thermal

! 将热电子种群直接加入 dN/dgamma 谱。
! Add a thermal-electron population directly into the dN/dgamma spectrum.
subroutine thermal_pop(Num_gam_e,gam_e,four_v,total_count,dN_gam_e)
    implicit none
    integer, intent(in) :: Num_gam_e
    real(8), intent(in), dimension(Num_gam_e) :: gam_e
    real(8), intent(in) :: four_v,total_count
    real(8), intent(inout), dimension(Num_gam_e) :: dN_gam_e
    real(8), dimension(Num_gam_e) :: shape_dnx
    real(8) :: theta

    if (total_count <= 0d0) return
    theta=thermal_theta(four_v)
    call thermal_shape(Num_gam_e,gam_e,theta,shape_dnx)
    dN_gam_e=dN_gam_e+total_count*shape_dnx/(gam_e)
end subroutine thermal_pop

end module electron_injection_profiles
