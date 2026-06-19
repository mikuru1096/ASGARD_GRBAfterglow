!f2py: skip
module electron_injection_profiles
    use constants
    use electron_energy_coordinate_common, only: electron_coord_log_four_velocity_sq, electron_coord_from_xgamma, &
                                                 electron_xgamma_from_coord, electron_gamma_from_coord, &
                                                 electron_dxgamma_dcoord
    use electron_radiation_kernel, only: besselk
    implicit none

contains

! 指数截断因子：γ > Gam_e_max 时返回 exp(1-γ/Gam_e_max)，否则返回1。
pure real(8) function electron_exp_cutoff_factor(gam,Gam_e_max)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: gam,Gam_e_max

    electron_exp_cutoff_factor=one
    if (gam > Gam_e_max) electron_exp_cutoff_factor=dexp(one-gam/Gam_e_max)
end function electron_exp_cutoff_factor

! 根据快冷却/慢冷却分支返回幂律注入参数：系数coeff和谱指数slope。
subroutine electron_initial_powerlaw_params(Para_N_e_ini,p,Gam_e_m,Gam_e_c,gam,active,coeff,slope)
    implicit real(8)(A-H,O-Z)
    logical, intent(out) :: active
    real(8), intent(in) :: Para_N_e_ini,p,Gam_e_m,Gam_e_c,gam
    real(8), intent(out) :: coeff,slope

    active=.false.
    coeff=zero
    slope=zero
    if (Gam_e_m > Gam_e_c) then
        if (gam < Gam_e_c) return
        active=.true.
        coeff=Para_N_e_ini*Gam_e_c
        slope=merge(2d0,p+one,gam < Gam_e_m)
        if (gam >= Gam_e_m) coeff=coeff*Gam_e_m**(p-one)
    else
        if (gam < Gam_e_m) return
        active=.true.
        coeff=Para_N_e_ini*(p-one)*Gam_e_m**(p-one)
        slope=merge(p,p+one,gam < Gam_e_c)
        if (gam >= Gam_e_c) coeff=coeff*Gam_e_c
    end if
end subroutine electron_initial_powerlaw_params

! 幂律+指数截断的 dN/dx 值：dN/dx = coeff * ln(10) * γ^(1-slope) * cutoff(γ)。
real(8) function electron_dnx_powerlaw_cutoff_value(x,coeff,slope,Gam_e_max)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: x,coeff,slope,Gam_e_max
    real(8) :: gam,cutoff_factor

    gam=ten**x
    if (gam <= zero .or. coeff <= zero) then
        electron_dnx_powerlaw_cutoff_value=zero
        return
    end if

    cutoff_factor=one
    if (Gam_e_max > zero .and. gam > Gam_e_max) cutoff_factor=dexp(one-gam/Gam_e_max)
    electron_dnx_powerlaw_cutoff_value=coeff*dlog(ten)*gam**(one-slope)*cutoff_factor
end function electron_dnx_powerlaw_cutoff_value

! 3点Gauss-Legendre积分：在[x_lo, x_hi]上对dN/dx进行数值积分。
real(8) function electron_dnx_gauss3_integral(coeff,slope,Gam_e_max,x_lo,x_hi)
    implicit real(8)(A-H,O-Z)
    integer :: I_q
    real(8), intent(in) :: coeff,slope,Gam_e_max,x_lo,x_hi
    real(8), parameter :: xi(3)=(/-dsqrt(3d0/5d0),zero,dsqrt(3d0/5d0)/)
    real(8), parameter :: wi(3)=(/5d0/9d0,8d0/9d0,5d0/9d0/)
    real(8) :: half_dx,x_mid,x_eval,quad

    if (x_hi <= x_lo) then
        electron_dnx_gauss3_integral=zero
        return
    end if

    half_dx=0.5d0*(x_hi-x_lo)
    x_mid=0.5d0*(x_hi+x_lo)
    quad=zero
    do I_q=1,3
        x_eval=x_mid+half_dx*xi(I_q)
        quad=quad+wi(I_q)*electron_dnx_powerlaw_cutoff_value(x_eval,coeff,slope,Gam_e_max)
    end do
    electron_dnx_gauss3_integral=half_dx*quad
end function electron_dnx_gauss3_integral

! 3点Gauss-Legendre积分：在四速度坐标y上对dN/dy=(dN/dx)(dx/dy)积分。
real(8) function electron_dny_gauss3_integral(coord_scale,coeff,slope,Gam_e_max,y_lo,y_hi)
    implicit none
    integer :: I_q
    real(8), intent(in) :: coord_scale,coeff,slope,Gam_e_max,y_lo,y_hi
    real(8), parameter :: xi(3)=(/-dsqrt(3d0/5d0),zero,dsqrt(3d0/5d0)/)
    real(8), parameter :: wi(3)=(/5d0/9d0,8d0/9d0,5d0/9d0/)
    real(8) :: half_dy,y_mid,y_eval,x_eval,quad

    if (y_hi <= y_lo) then
        electron_dny_gauss3_integral=zero
        return
    end if

    half_dy=0.5d0*(y_hi-y_lo)
    y_mid=0.5d0*(y_hi+y_lo)
    quad=zero
    do I_q=1,3
        y_eval=y_mid+half_dy*xi(I_q)
        x_eval=electron_xgamma_from_coord(electron_coord_log_four_velocity_sq,coord_scale,y_eval)
        quad=quad+wi(I_q)*electron_dnx_powerlaw_cutoff_value(x_eval,coeff,slope,Gam_e_max) &
             *electron_dxgamma_dcoord(electron_coord_log_four_velocity_sq,coord_scale,y_eval)
    end do
    electron_dny_gauss3_integral=half_dy*quad
end function electron_dny_gauss3_integral

! 将激活区间上的dN/dx积分累加到acc，在Gam_e_max处自动分段处理截断。
subroutine electron_add_dnx_segment(cell_lo,cell_hi,active_lo,active_hi,coeff,slope,Gam_e_max,acc)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: cell_lo,cell_hi,active_lo,active_hi,coeff,slope,Gam_e_max
    real(8), intent(inout) :: acc
    real(8) :: x_lo,x_hi,x_cut

    if (coeff <= zero .or. cell_hi <= cell_lo .or. active_hi <= active_lo) return
    x_lo=max(cell_lo,active_lo)
    x_hi=min(cell_hi,active_hi)
    if (x_hi <= x_lo) return

    x_cut=dlog10(max(Gam_e_max,1d-300))
    if (x_lo < x_cut .and. x_hi > x_cut) then
        acc=acc+electron_dnx_gauss3_integral(coeff,slope,Gam_e_max,x_lo,x_cut)
        acc=acc+electron_dnx_gauss3_integral(coeff,slope,Gam_e_max,x_cut,x_hi)
    else
        acc=acc+electron_dnx_gauss3_integral(coeff,slope,Gam_e_max,x_lo,x_hi)
    end if
end subroutine electron_add_dnx_segment

! 将激活区间上的dN/dy积分累加到acc，在Gam_e_max处自动分段处理截断。
subroutine electron_add_dny_segment(cell_lo,cell_hi,active_lo,active_hi,coord_scale,coeff,slope,Gam_e_max,acc)
    implicit none
    real(8), intent(in) :: cell_lo,cell_hi,active_lo,active_hi,coord_scale,coeff,slope,Gam_e_max
    real(8), intent(inout) :: acc
    real(8) :: y_lo,y_hi,y_cut

    if (coeff <= zero .or. cell_hi <= cell_lo .or. active_hi <= active_lo) return
    y_lo=max(cell_lo,active_lo)
    y_hi=min(cell_hi,active_hi)
    if (y_hi <= y_lo) return

    y_cut=electron_coord_from_xgamma(electron_coord_log_four_velocity_sq,coord_scale,dlog10(Gam_e_max))
    if (y_lo < y_cut .and. y_hi > y_cut) then
        acc=acc+electron_dny_gauss3_integral(coord_scale,coeff,slope,Gam_e_max,y_lo,y_cut)
        acc=acc+electron_dny_gauss3_integral(coord_scale,coeff,slope,Gam_e_max,y_cut,y_hi)
    else
        acc=acc+electron_dny_gauss3_integral(coord_scale,coeff,slope,Gam_e_max,y_lo,y_hi)
    end if
end subroutine electron_add_dny_segment

! 由网格中心值推导log10(gamma)的单元边界。
subroutine electron_profile_log_cell_edges(Num_gam_e,gam_e,x_edge)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: gam_e(Num_gam_e)
    real(8), intent(out) :: x_edge(Num_gam_e+1)

    x_edge(1)=dlog10(gam_e(1))-0.5d0*(dlog10(gam_e(2))-dlog10(gam_e(1)))
    do I_gam_e=2,Num_gam_e
        x_edge(I_gam_e)=0.5d0*(dlog10(gam_e(I_gam_e-1))+dlog10(gam_e(I_gam_e)))
    end do
    x_edge(Num_gam_e+1)=dlog10(gam_e(Num_gam_e))+0.5d0*(dlog10(gam_e(Num_gam_e))-dlog10(gam_e(Num_gam_e-1)))
end subroutine electron_profile_log_cell_edges

! 生成快/慢冷却幂律+指数截断的初始电子谱 dN/dγ（网格中心值）。
subroutine electron_initial_powerlaw_exp_cutoff(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,dN_gam_e_1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,gam_e(Num_gam_e)
    real(8), intent(out) :: dN_gam_e_1(Num_gam_e)
    logical :: active
    real(8) :: coeff,slope

    dN_gam_e_1=zero
    if (Gam_e_max <= zero) return

    do I_gam_e=1,Num_gam_e
        call electron_initial_powerlaw_params(Para_N_e_ini,p,Gam_e_m,Gam_e_c,gam_e(I_gam_e),active,coeff,slope)
        if (active) dN_gam_e_1(I_gam_e)=coeff*gam_e(I_gam_e)**(-slope)*electron_exp_cutoff_factor(gam_e(I_gam_e),Gam_e_max)
    end do
end subroutine electron_initial_powerlaw_exp_cutoff

! 生成快/慢冷却幂律+指数截断的初始电子谱 dN/dx（网格单元平均，保正/守恒）。
subroutine electron_initial_powerlaw_exp_cutoff_edges(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,x_edge,dN_x_1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,x_edge(Num_gam_e+1)
    real(8), intent(out) :: dN_x_1(Num_gam_e)
    real(8) :: cell_lo,cell_hi,dx_cell,seg_sum,x_m,x_c,coeff_lo,coeff_hi,huge_x

    dN_x_1=zero
    if (Gam_e_max <= zero) return

    x_m=dlog10(max(Gam_e_m,1d-300))
    x_c=dlog10(max(Gam_e_c,1d-300))
    huge_x=1d300

    do I_gam_e=1,Num_gam_e
        cell_lo=x_edge(I_gam_e)
        cell_hi=x_edge(I_gam_e+1)
        dx_cell=cell_hi-cell_lo
        if (dx_cell <= zero) cycle

        seg_sum=zero
        if (Gam_e_m > Gam_e_c) then
            coeff_lo=Para_N_e_ini*Gam_e_c
            coeff_hi=Para_N_e_ini*Gam_e_c*Gam_e_m**(p-one)
            call electron_add_dnx_segment(cell_lo,cell_hi,x_c,x_m,coeff_lo,2d0,Gam_e_max,seg_sum)
            call electron_add_dnx_segment(cell_lo,cell_hi,x_m,huge_x,coeff_hi,p+one,Gam_e_max,seg_sum)
        else
            coeff_lo=Para_N_e_ini*(p-one)*Gam_e_m**(p-one)
            coeff_hi=coeff_lo*Gam_e_c
            call electron_add_dnx_segment(cell_lo,cell_hi,x_m,x_c,coeff_lo,p,Gam_e_max,seg_sum)
            call electron_add_dnx_segment(cell_lo,cell_hi,x_c,huge_x,coeff_hi,p+one,Gam_e_max,seg_sum)
        end if
        dN_x_1(I_gam_e)=seg_sum/dx_cell
    end do
end subroutine electron_initial_powerlaw_exp_cutoff_edges

! 生成快/慢冷却幂律+指数截断的初始电子谱 dN/dy，y为四速度坐标。
subroutine electron_initial_powerlaw_exp_cutoff_coord_edges(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max, &
                                                            Num_gam_e,coord_edge,coord_scale,dN_y_1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,coord_edge(Num_gam_e+1),coord_scale
    real(8), intent(out) :: dN_y_1(Num_gam_e)
    real(8) :: cell_lo,cell_hi,dy_cell,seg_sum,y_m,y_c,coeff_lo,coeff_hi,huge_y

    dN_y_1=zero
    if (Gam_e_max <= zero) return

    y_m=electron_coord_from_xgamma(electron_coord_log_four_velocity_sq,coord_scale,dlog10(Gam_e_m))
    y_c=electron_coord_from_xgamma(electron_coord_log_four_velocity_sq,coord_scale,dlog10(Gam_e_c))
    huge_y=1d300

    do I_gam_e=1,Num_gam_e
        cell_lo=coord_edge(I_gam_e)
        cell_hi=coord_edge(I_gam_e+1)
        dy_cell=cell_hi-cell_lo
        if (dy_cell <= zero) cycle

        seg_sum=zero
        if (Gam_e_m > Gam_e_c) then
            coeff_lo=Para_N_e_ini*Gam_e_c
            coeff_hi=Para_N_e_ini*Gam_e_c*Gam_e_m**(p-one)
            call electron_add_dny_segment(cell_lo,cell_hi,y_c,y_m,coord_scale,coeff_lo,2d0,Gam_e_max,seg_sum)
            call electron_add_dny_segment(cell_lo,cell_hi,y_m,huge_y,coord_scale,coeff_hi,p+one,Gam_e_max,seg_sum)
        else
            coeff_lo=Para_N_e_ini*(p-one)*Gam_e_m**(p-one)
            coeff_hi=coeff_lo*Gam_e_c
            call electron_add_dny_segment(cell_lo,cell_hi,y_m,y_c,coord_scale,coeff_lo,p,Gam_e_max,seg_sum)
            call electron_add_dny_segment(cell_lo,cell_hi,y_c,huge_y,coord_scale,coeff_hi,p+one,Gam_e_max,seg_sum)
        end if
        dN_y_1(I_gam_e)=seg_sum/dy_cell
    end do
end subroutine electron_initial_powerlaw_exp_cutoff_coord_edges

! 构建幂律+指数截断源项 dF/dx（网格中心值）。
subroutine electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,Q,p,dF1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: gam_e(Num_gam_e),Gam_e_m,Gam_e_max,Q,p
    real(8), intent(out) :: dF1(Num_gam_e)

    dF1=zero
    if (Gam_e_max <= zero) return

    do I_gam_e=1,Num_gam_e
        if (gam_e(I_gam_e) <= Gam_e_m) cycle
        dF1(I_gam_e)=Q*gam_e(I_gam_e)**(one-p)*dlog(ten)*electron_exp_cutoff_factor(gam_e(I_gam_e),Gam_e_max)
    end do
end subroutine electron_build_source_term_exp_cutoff

! 构建幂律+指数截断源项 dF/dx（网格单元平均，保正/守恒）。
subroutine electron_build_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m,Gam_e_max,Q,p,dF1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1),Gam_e_m,Gam_e_max,Q,p
    real(8), intent(out) :: dF1(Num_gam_e)
    real(8) :: cell_lo,cell_hi,dx_cell,seg_sum,x_m,huge_x

    dF1=zero
    if (Gam_e_max <= zero .or. Q <= zero) return

    x_m=dlog10(max(Gam_e_m,1d-300))
    huge_x=1d300
    do I_gam_e=1,Num_gam_e
        cell_lo=x_edge(I_gam_e)
        cell_hi=x_edge(I_gam_e+1)
        dx_cell=cell_hi-cell_lo
        if (dx_cell <= zero) cycle
        seg_sum=zero
        call electron_add_dnx_segment(cell_lo,cell_hi,x_m,huge_x,Q,p,Gam_e_max,seg_sum)
        dF1(I_gam_e)=seg_sum/dx_cell
    end do
end subroutine electron_build_source_term_exp_cutoff_edges

! 构建FS幂律+指数截断源项 dF/dy，y为四速度坐标，物理谱仍为γ^{-p}。
subroutine electron_build_source_term_exp_cutoff_coord_edges(Num_gam_e,coord_edge,coord_scale,Gam_e_m,Gam_e_max,Q,p,dF1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: coord_edge(Num_gam_e+1),coord_scale,Gam_e_m,Gam_e_max,Q,p
    real(8), intent(out) :: dF1(Num_gam_e)
    real(8) :: cell_lo,cell_hi,dy_cell,seg_sum,y_m,huge_y

    dF1=zero
    if (Gam_e_max <= zero .or. Q <= zero) return

    y_m=electron_coord_from_xgamma(electron_coord_log_four_velocity_sq,coord_scale,dlog10(Gam_e_m))
    huge_y=1d300
    do I_gam_e=1,Num_gam_e
        cell_lo=coord_edge(I_gam_e)
        cell_hi=coord_edge(I_gam_e+1)
        dy_cell=cell_hi-cell_lo
        if (dy_cell <= zero) cycle
        seg_sum=zero
        call electron_add_dny_segment(cell_lo,cell_hi,y_m,huge_y,coord_scale,Q,p,Gam_e_max,seg_sum)
        dF1(I_gam_e)=seg_sum/dy_cell
    end do
end subroutine electron_build_source_term_exp_cutoff_coord_edges

! 构建以动能幂律归一的反向激波源项：dN/dx ∝ γ (γ-1)^(-p) exp cutoff。
subroutine electron_build_kinetic_source_term_exp_cutoff_edges(Num_gam_e,x_edge,Gam_e_m,Gam_e_max,Q,p,dF1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e,I_q
    real(8), intent(in) :: x_edge(Num_gam_e+1),Gam_e_m,Gam_e_max,Q,p
    real(8), intent(out) :: dF1(Num_gam_e)
    real(8), parameter :: xi(3)=(/-dsqrt(3d0/5d0),zero,dsqrt(3d0/5d0)/)
    real(8), parameter :: wi(3)=(/5d0/9d0,8d0/9d0,5d0/9d0/)
    real(8) :: cell_lo,cell_hi,dx_cell,half_dx,x_mid,x_eval,gam,cutoff_factor,cell_sum,shape_norm

    dF1=zero
    shape_norm=zero
    if (Gam_e_max <= zero .or. Q <= zero) return

    do I_gam_e=1,Num_gam_e
        cell_lo=x_edge(I_gam_e)
        cell_hi=x_edge(I_gam_e+1)
        dx_cell=cell_hi-cell_lo
        if (dx_cell <= zero) cycle
        half_dx=0.5d0*dx_cell
        x_mid=0.5d0*(cell_lo+cell_hi)
        cell_sum=zero
        do I_q=1,3
            x_eval=x_mid+half_dx*xi(I_q)
            gam=ten**x_eval
            if (gam > Gam_e_m) then
                cutoff_factor=electron_exp_cutoff_factor(gam,Gam_e_max)
                cell_sum=cell_sum+wi(I_q)*gam*dlog(ten)*(gam-one)**(-p)*cutoff_factor
            end if
        end do
        cell_sum=half_dx*cell_sum
        dF1(I_gam_e)=cell_sum/dx_cell
        shape_norm=shape_norm+cell_sum
    end do
    dF1=Q*dF1/shape_norm
end subroutine electron_build_kinetic_source_term_exp_cutoff_edges

! 构建以动能幂律归一的源项 dN/dy；y 为 log four-velocity-squared 坐标。
subroutine electron_build_kinetic_source_term_exp_cutoff_coord_edges(Num_gam_e,coord_edge,coord_scale,Gam_e_m,Gam_e_max,Q,p,dF1)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e,I_q
    real(8), intent(in) :: coord_edge(Num_gam_e+1),coord_scale,Gam_e_m,Gam_e_max,Q,p
    real(8), intent(out) :: dF1(Num_gam_e)
    real(8), parameter :: xi(3)=(/-dsqrt(3d0/5d0),zero,dsqrt(3d0/5d0)/)
    real(8), parameter :: wi(3)=(/5d0/9d0,8d0/9d0,5d0/9d0/)
    real(8) :: cell_lo,cell_hi,dy_cell,half_dy,y_mid,y_eval,gam,dxdy,cell_sum,shape_norm

    dF1=zero
    shape_norm=zero
    if (Gam_e_max <= zero .or. Q <= zero) return

    do I_gam_e=1,Num_gam_e
        cell_lo=coord_edge(I_gam_e)
        cell_hi=coord_edge(I_gam_e+1)
        dy_cell=cell_hi-cell_lo
        half_dy=0.5d0*dy_cell
        y_mid=0.5d0*(cell_lo+cell_hi)
        cell_sum=zero
        do I_q=1,3
            y_eval=y_mid+half_dy*xi(I_q)
            gam=electron_gamma_from_coord(electron_coord_log_four_velocity_sq,coord_scale,y_eval)
            if (gam > Gam_e_m) then
                dxdy=electron_dxgamma_dcoord(electron_coord_log_four_velocity_sq,coord_scale,y_eval)
                cell_sum=cell_sum+wi(I_q)*gam*dlog(ten)*dxdy*(gam-one)**(-p)*electron_exp_cutoff_factor(gam,Gam_e_max)
            end if
        end do
        cell_sum=half_dy*cell_sum
        dF1(I_gam_e)=cell_sum/dy_cell
        shape_norm=shape_norm+cell_sum
    end do
    dF1=Q*dF1/shape_norm
end subroutine electron_build_kinetic_source_term_exp_cutoff_coord_edges

pure real(8) function electron_thermal_theta(four_v)
    implicit none
    real(8), intent(in) :: four_v

    if (four_v <= zero) error stop 'electron_thermal_theta requires four_v > 0'
    electron_thermal_theta=four_v*(four_v+1.07d0*four_v*four_v)/(3d0*(one+four_v+1.07d0*four_v*four_v))
end function electron_thermal_theta

subroutine electron_build_thermal_shape_dnx(Num_gam_e,gam_e,theta,shape_dnx)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: gam_e(Num_gam_e),theta
    real(8), intent(out) :: shape_dnx(Num_gam_e)
    real(8) :: shape_dgam(Num_gam_e),k2_theta,norm_dgam

    if (theta <= zero) error stop 'electron_build_thermal_shape_dnx requires theta > 0'
    if (any(gam_e <= one)) error stop 'thermal electron grid requires gamma > 1'

    k2_theta=besselk(one/theta)
    if (k2_theta <= zero) error stop 'besselk normalization vanished in thermal electron source'

    do I_gam_e=1,Num_gam_e
        shape_dgam(I_gam_e)=gam_e(I_gam_e)**2*dsqrt(one-one/gam_e(I_gam_e)**2) &
                          *dexp(-gam_e(I_gam_e)/theta)/(theta*k2_theta)
    end do
    norm_dgam=sum((shape_dgam(2:Num_gam_e)+shape_dgam(1:Num_gam_e-1)) &
            *(gam_e(2:Num_gam_e)-gam_e(1:Num_gam_e-1)))/two
    if (norm_dgam <= zero) error stop 'thermal electron distribution normalization is non-positive'

    shape_dnx=shape_dgam/norm_dgam*gam_e*dlog(ten)
end subroutine electron_build_thermal_shape_dnx

subroutine electron_add_thermal_source_term(Num_gam_e,gam_e,four_v,total_count,dF1)
    implicit none
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: gam_e(Num_gam_e),four_v,total_count
    real(8), intent(inout) :: dF1(Num_gam_e)
    real(8) :: shape_dnx(Num_gam_e),theta

    if (total_count <= zero) return
    theta=electron_thermal_theta(four_v)
    call electron_build_thermal_shape_dnx(Num_gam_e,gam_e,theta,shape_dnx)
    dF1=dF1+total_count*shape_dnx
end subroutine electron_add_thermal_source_term

subroutine electron_add_thermal_population(Num_gam_e,gam_e,four_v,total_count,dN_gam_e)
    implicit none
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: gam_e(Num_gam_e),four_v,total_count
    real(8), intent(inout) :: dN_gam_e(Num_gam_e)
    real(8) :: shape_dnx(Num_gam_e),theta

    if (total_count <= zero) return
    theta=electron_thermal_theta(four_v)
    call electron_build_thermal_shape_dnx(Num_gam_e,gam_e,theta,shape_dnx)
    dN_gam_e=dN_gam_e+total_count*shape_dnx/(gam_e*dlog(ten))
end subroutine electron_add_thermal_population

end module electron_injection_profiles
