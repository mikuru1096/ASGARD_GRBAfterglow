module electron_transport_common
    use constants
    implicit none
    integer, parameter :: charint_quad_order = 4
    real(8), parameter :: charint_cfl_relax = 32d0
    real(8), parameter :: charint_substep_rtol_relax = 8d0
    integer, parameter :: electron_cooling_affine = 0
    integer, parameter :: electron_cooling_piecewise = 1
    real(8), parameter :: inv_log_ten = 4.3429448190325182765d-1
    real(8), parameter :: tiny_u_char = 1d-300
    real(8), parameter :: charint_quad_nodes(charint_quad_order) = &
        (/6.9431844202973714d-2, 3.3000947820757187d-1, 6.6999052179242813d-1, 9.3056815579702629d-1/)
    real(8), parameter :: charint_quad_weights(charint_quad_order) = &
        (/1.7392742256872693d-1, 3.2607257743127307d-1, 3.2607257743127307d-1, 1.7392742256872693d-1/)
contains

! 准备隐式迎风输运系数：计算主对角元principal和上三角系数temp1。
subroutine electron_prepare_implicit_coeffs_common(Num_gam_e,diag_base,up,principal,temp1)
    implicit none
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: diag_base,up(Num_gam_e-1)
    real(8), intent(out) :: principal(Num_gam_e),temp1(Num_gam_e-1)

    principal(2:Num_gam_e)=diag_base-up
    principal(1)=principal(2)
    temp1=up/(principal(2:Num_gam_e)+principal(1:Num_gam_e-1))*two
end subroutine electron_prepare_implicit_coeffs_common

! 隐式迎风输运后向回代求解：从高能向低能端回代，解截断为非负。
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

! 用三个点做二次插值。
real(8) function electron_quadratic_interp3(x0,y0,x1,y1,x2,y2,x)
    implicit none
    real(8), intent(in) :: x0,y0,x1,y1,x2,y2,x
    real(8) :: l0,l1,l2
    l0=((x-x1)*(x-x2))/max((x0-x1)*(x0-x2),1d-30)
    l1=((x-x0)*(x-x2))/max((x1-x0)*(x1-x2),1d-30)
    l2=((x-x0)*(x-x1))/max((x2-x0)*(x2-x1),1d-30)
    electron_quadratic_interp3=y0*l0+y1*l1+y2*l2
end function electron_quadratic_interp3



! 构造非均匀网格 PPM 单元左右界面值并做单调性限制。
subroutine electron_ppm_interfaces_nonuniform(Num_gam_e,x_edge,q,q_left,q_right)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1),q(Num_gam_e)
    real(8), intent(out) :: q_left(Num_gam_e),q_right(Num_gam_e)
    real(8) :: x_center(Num_gam_e),dx_cell(Num_gam_e)
    real(8) :: q_min,q_max,q_bar,dq
    do I_gam_e=1,Num_gam_e
        x_center(I_gam_e)=0.5d0*(x_edge(I_gam_e)+x_edge(I_gam_e+1))
        dx_cell(I_gam_e)=x_edge(I_gam_e+1)-x_edge(I_gam_e)
    end do
    q_left=q
    q_right=q
    if (Num_gam_e >= 3) then
        q_right(1)=electron_quadratic_interp3(x_center(1),q(1),x_center(2),q(2),x_center(3),q(3),x_edge(2))
        q_left(Num_gam_e)=electron_quadratic_interp3(x_center(Num_gam_e-2),q(Num_gam_e-2), &
                                                     x_center(Num_gam_e-1),q(Num_gam_e-1), &
                                                     x_center(Num_gam_e),q(Num_gam_e),x_edge(Num_gam_e))
    end if
    do I_gam_e=2,Num_gam_e-1
        q_left(I_gam_e)=electron_quadratic_interp3(x_center(I_gam_e-1),q(I_gam_e-1), &
                                                   x_center(I_gam_e),q(I_gam_e), &
                                                   x_center(I_gam_e+1),q(I_gam_e+1),x_edge(I_gam_e))
        q_right(I_gam_e)=electron_quadratic_interp3(x_center(I_gam_e-1),q(I_gam_e-1), &
                                                    x_center(I_gam_e),q(I_gam_e), &
                                                    x_center(I_gam_e+1),q(I_gam_e+1),x_edge(I_gam_e+1))
        q_min=min(q(I_gam_e-1),min(q(I_gam_e),q(I_gam_e+1)))
        q_max=max(q(I_gam_e-1),max(q(I_gam_e),q(I_gam_e+1)))
        q_left(I_gam_e)=max(q_min,min(q_max,q_left(I_gam_e)))
        q_right(I_gam_e)=max(q_min,min(q_max,q_right(I_gam_e)))
        if ((q(I_gam_e+1)-q(I_gam_e))*(q(I_gam_e)-q(I_gam_e-1)) <= zero) then
            q_left(I_gam_e)=q(I_gam_e)
            q_right(I_gam_e)=q(I_gam_e)
        else
            q_bar=0.5d0*(q_left(I_gam_e)+q_right(I_gam_e))
            dq=q_right(I_gam_e)-q_left(I_gam_e)
            if (dq*(q(I_gam_e)-q_bar) > dq*dq/6d0) then
                q_left(I_gam_e)=3d0*q(I_gam_e)-2d0*q_right(I_gam_e)
            else if (dq*(q(I_gam_e)-q_bar) < -dq*dq/6d0) then
                q_right(I_gam_e)=3d0*q(I_gam_e)-2d0*q_left(I_gam_e)
            end if
        end if
    end do
    if (Num_gam_e >= 2) then
        q_min=min(q(1),q(2))
        q_max=max(q(1),q(2))
        q_left(1)=max(q_min,min(q_max,q_left(1)))
        q_right(1)=max(q_min,min(q_max,q_right(1)))
        q_min=min(q(Num_gam_e-1),q(Num_gam_e))
        q_max=max(q(Num_gam_e-1),q(Num_gam_e))
        q_left(Num_gam_e)=max(q_min,min(q_max,q_left(Num_gam_e)))
        q_right(Num_gam_e)=max(q_min,min(q_max,q_right(Num_gam_e)))
    end if
end subroutine electron_ppm_interfaces_nonuniform



! 计算非均匀网格 PPM 分段抛物线的前缀积分。
subroutine electron_ppm_prefix_nonuniform(Num_gam_e,x_edge,q,q_left,q_right,prefix)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1),q(Num_gam_e),q_left(Num_gam_e),q_right(Num_gam_e)
    real(8), intent(out) :: prefix(0:Num_gam_e)
    real(8) :: dx_cell
    prefix(0)=zero
    do I_gam_e=1,Num_gam_e
        dx_cell=x_edge(I_gam_e+1)-x_edge(I_gam_e)
        prefix(I_gam_e)=prefix(I_gam_e-1)+ &
            electron_ppm_cell_int(q(I_gam_e),q_left(I_gam_e),q_right(I_gam_e), &
                                  x_edge(I_gam_e),dx_cell,x_edge(I_gam_e),x_edge(I_gam_e+1))
    end do
end subroutine electron_ppm_prefix_nonuniform



! 在任意位置求非均匀网格 PPM 前缀积分值。
real(8) function electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge,q,q_left,q_right,prefix,x_eval)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e,left,right,mid
    real(8), intent(in) :: x_edge(Num_gam_e+1),q(Num_gam_e),q_left(Num_gam_e),q_right(Num_gam_e),prefix(0:Num_gam_e),x_eval
    real(8) :: xa,dx_cell
    xa=max(x_edge(1),min(x_eval,x_edge(Num_gam_e+1)))
    if (xa <= x_edge(1)) then
        electron_ppm_prefix_eval_nonuniform=zero
        return
    end if
    if (xa >= x_edge(Num_gam_e+1)) then
        electron_ppm_prefix_eval_nonuniform=prefix(Num_gam_e)
        return
    end if
    left=1
    right=Num_gam_e
    do while (left < right)
        mid=(left+right)/2
        if (x_edge(mid+1) >= xa) then
            right=mid
        else
            left=mid+1
        end if
    end do
    I_gam_e=left
    dx_cell=x_edge(I_gam_e+1)-x_edge(I_gam_e)
    electron_ppm_prefix_eval_nonuniform=prefix(I_gam_e-1)+ &
        electron_ppm_cell_int(q(I_gam_e),q_left(I_gam_e),q_right(I_gam_e), &
                              x_edge(I_gam_e),dx_cell,x_edge(I_gam_e),xa)
end function electron_ppm_prefix_eval_nonuniform



! 预处理非均匀网格守恒重映射所需的 PPM 界面值和前缀积分。
subroutine electron_prepare_conservative_remap_nonuniform(Num_gam_e,x_edge,q,q_left,q_right,prefix)
    implicit none
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1),q(Num_gam_e)
    real(8), intent(out) :: q_left(Num_gam_e),q_right(Num_gam_e),prefix(0:Num_gam_e)
    call electron_ppm_interfaces_nonuniform(Num_gam_e,x_edge,q,q_left,q_right)
    call electron_ppm_prefix_nonuniform(Num_gam_e,x_edge,q,q_left,q_right,prefix)
end subroutine electron_prepare_conservative_remap_nonuniform



! 将 x=log10(gamma) 边界转换为 u=1/gamma 边界。
subroutine electron_u_edges_from_x(Num_gam_e,x_edge,u_edge)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1)
    real(8), intent(out) :: u_edge(Num_gam_e+1)
    do I_gam_e=1,Num_gam_e+1
        u_edge(I_gam_e)=ten**(-x_edge(I_gam_e))
    end do
end subroutine electron_u_edges_from_x



! 将 u=1/gamma 转换回 x=log10(gamma)。
real(8) function electron_x_from_u(u)
    implicit none
    real(8), intent(in) :: u
    electron_x_from_u=-dlog(max(u,tiny_u_char))*inv_log_ten
end function electron_x_from_u



! 批量回溯多个滞后时间下的全局仿射 u 特征线边界。
subroutine electron_trace_affine_u_edges_batch(Num_gam_e,Num_lag,u_now_edge,lag_arr,a_u,b_u,x_back_batch)
    implicit none
    integer, intent(in) :: Num_gam_e,Num_lag
    integer :: I_gam_e,I_lag
    real(8), intent(in) :: u_now_edge(Num_gam_e+1),lag_arr(Num_lag),a_u,b_u
    real(8), intent(out) :: x_back_batch(Num_gam_e+1,Num_lag)
    real(8) :: u_prev,shift,exp_fac(Num_lag)
    if (abs(b_u) <= 1d-30) then
        do I_lag=1,Num_lag
            do I_gam_e=1,Num_gam_e+1
                u_prev=u_now_edge(I_gam_e)-a_u*lag_arr(I_lag)
                x_back_batch(I_gam_e,I_lag)=electron_x_from_u(u_prev)
            end do
        end do
        return
    end if
    shift=a_u/b_u
    do I_lag=1,Num_lag
        exp_fac(I_lag)=dexp(-b_u*lag_arr(I_lag))
    end do
    do I_lag=1,Num_lag
        do I_gam_e=1,Num_gam_e+1
            u_prev=(u_now_edge(I_gam_e)+shift)*exp_fac(I_lag)-shift
            x_back_batch(I_gam_e,I_lag)=electron_x_from_u(u_prev)
        end do
    end do
end subroutine electron_trace_affine_u_edges_batch



! 由冷却率在 u 空间构造分段仿射特征线系数。
subroutine electron_build_piecewise_affine_u(Num_gam_e,x_edge,gam_e,dEl,R_loc,u_edge,a_cell,b_cell)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1),gam_e(Num_gam_e),dEl(Num_gam_e),R_loc
    real(8), intent(out) :: u_edge(Num_gam_e+1),a_cell(Num_gam_e),b_cell(Num_gam_e)
    real(8) :: u_center(Num_gam_e),w_center(Num_gam_e),w_edge(Num_gam_e+1)
    real(8) :: frac,du,slope
    call electron_u_edges_from_x(Num_gam_e,x_edge,u_edge)
    do I_gam_e=1,Num_gam_e
        u_center(I_gam_e)=one/gam_e(I_gam_e)
        w_center(I_gam_e)=u_center(I_gam_e)*dEl(I_gam_e)
    end do
    w_edge(1)=w_center(1)
    do I_gam_e=2,Num_gam_e
        if (abs(u_center(I_gam_e)-u_center(I_gam_e-1)) <= 1d-40) then
            w_edge(I_gam_e)=0.5d0*(w_center(I_gam_e-1)+w_center(I_gam_e))
        else
            frac=(u_edge(I_gam_e)-u_center(I_gam_e-1))/(u_center(I_gam_e)-u_center(I_gam_e-1))
            w_edge(I_gam_e)=w_center(I_gam_e-1)+frac*(w_center(I_gam_e)-w_center(I_gam_e-1))
        end if
    end do
    w_edge(Num_gam_e+1)=w_center(Num_gam_e)
    do I_gam_e=1,Num_gam_e
        du=u_edge(I_gam_e+1)-u_edge(I_gam_e)
        if (abs(du) <= 1d-40) then
            slope=zero
        else
            slope=(w_edge(I_gam_e+1)-w_edge(I_gam_e))/du
        end if
        a_cell(I_gam_e)=w_edge(I_gam_e)-slope*u_edge(I_gam_e)
        b_cell(I_gam_e)=slope+one/R_loc
    end do
end subroutine electron_build_piecewise_affine_u



! 借助上次单元索引在递减 u 网格中定位给定 u。
subroutine electron_find_u_cell_desc_hint(Num_gam_e,u_edge,u_val,I_cell)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer, intent(inout) :: I_cell
    real(8), intent(in) :: u_edge(Num_gam_e+1),u_val
    if (u_val >= u_edge(1)) then
        I_cell=1
        return
    end if
    if (u_val <= u_edge(Num_gam_e+1)) then
        I_cell=Num_gam_e
        return
    end if
    I_cell=max(1,min(Num_gam_e,I_cell))
    do while (I_cell < Num_gam_e .and. u_val <= u_edge(I_cell+1))
        I_cell=I_cell+1
    end do
    do while (I_cell > 1 .and. u_val > u_edge(I_cell))
        I_cell=I_cell-1
    end do
end subroutine electron_find_u_cell_desc_hint



! 从指定单元开始沿分段仿射 u 特征线回溯单个边界。
subroutine electron_trace_piecewise_affine_u_edge_from_cell(Num_gam_e,u_edge,a_cell,b_cell,lag,u_now,I_cell_start,u_back)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer, intent(in) :: I_cell_start
    integer :: I_cell
    real(8), intent(in) :: u_edge(Num_gam_e+1),a_cell(Num_gam_e),b_cell(Num_gam_e),lag,u_now
    real(8), intent(out) :: u_back
    real(8) :: s_rem,u_cur,g_cur,u_bound,s_hit,a_cur,b_cur,shift,num,den,vel_back
    if (lag <= zero) then
        u_back=u_now
        return
    end if
    I_cell=max(1,min(Num_gam_e,I_cell_start))
    u_cur=u_now
    s_rem=lag
    do while (s_rem > 1d-14*max(lag,one))
        a_cur=a_cell(I_cell)
        b_cur=b_cell(I_cell)
        g_cur=a_cur+b_cur*u_cur
        if (abs(g_cur) <= 1d-40) then
            u_back=u_cur
            return
        end if
        if (g_cur > zero) then
            if (I_cell == Num_gam_e) then
                if (abs(b_cur) <= 1d-30) then
                    u_back=u_cur-a_cur*s_rem
                else
                    shift=a_cur/b_cur
                    u_back=(u_cur+shift)*dexp(-b_cur*s_rem)-shift
                end if
                return
            end if
            u_bound=u_edge(I_cell+1)
            s_hit=1d300
            if (abs(u_bound-u_cur) <= 1d-30*max(one,abs(u_cur),abs(u_bound))) then
                s_hit=zero
            else if (abs(b_cur) <= 1d-30) then
                vel_back=-a_cur
                if (abs(vel_back) > 1d-40) s_hit=(u_bound-u_cur)/vel_back
            else
                shift=a_cur/b_cur
                num=u_cur+shift
                den=u_bound+shift
                if (num > zero .and. den > zero) s_hit=dlog(num/den)/b_cur
            end if
            if (s_hit < zero) s_hit=1d300
            if (s_hit >= s_rem .or. s_hit >= 1d290) then
                if (abs(b_cur) <= 1d-30) then
                    u_back=u_cur-a_cur*s_rem
                else
                    if (abs(b_cur) > 1d-30) shift=a_cur/b_cur
                    u_back=(u_cur+shift)*dexp(-b_cur*s_rem)-shift
                end if
                return
            end if
            s_rem=s_rem-s_hit
            I_cell=min(Num_gam_e,I_cell+1)
            u_cur=max(u_edge(Num_gam_e+1),u_bound*(one-1d-12))
        else
            if (I_cell == 1) then
                u_bound=u_edge(1)
                s_hit=1d300
                if (abs(u_bound-u_cur) <= 1d-30*max(one,abs(u_cur),abs(u_bound))) then
                    s_hit=zero
                else if (abs(b_cur) <= 1d-30) then
                    vel_back=-a_cur
                    if (abs(vel_back) > 1d-40) s_hit=(u_bound-u_cur)/vel_back
                else
                    shift=a_cur/b_cur
                    num=u_cur+shift
                    den=u_bound+shift
                    if (num > zero .and. den > zero) s_hit=dlog(num/den)/b_cur
                end if
                if (s_hit < zero) s_hit=1d300
                if (s_hit >= s_rem .or. s_hit >= 1d290) then
                    if (abs(b_cur) <= 1d-30) then
                        u_back=u_cur-a_cur*s_rem
                    else
                        if (abs(b_cur) > 1d-30) shift=a_cur/b_cur
                        u_back=(u_cur+shift)*dexp(-b_cur*s_rem)-shift
                    end if
                else
                    u_back=u_edge(1)
                end if
                return
            end if
            u_bound=u_edge(I_cell)
            s_hit=1d300
            if (abs(u_bound-u_cur) <= 1d-30*max(one,abs(u_cur),abs(u_bound))) then
                s_hit=zero
            else if (abs(b_cur) <= 1d-30) then
                vel_back=-a_cur
                if (abs(vel_back) > 1d-40) s_hit=(u_bound-u_cur)/vel_back
            else
                shift=a_cur/b_cur
                num=u_cur+shift
                den=u_bound+shift
                if (num > zero .and. den > zero) s_hit=dlog(num/den)/b_cur
            end if
            if (s_hit < zero) s_hit=1d300
            if (s_hit >= s_rem .or. s_hit >= 1d290) then
                if (abs(b_cur) <= 1d-30) then
                    u_back=u_cur-a_cur*s_rem
                else
                    if (abs(b_cur) > 1d-30) shift=a_cur/b_cur
                    u_back=(u_cur+shift)*dexp(-b_cur*s_rem)-shift
                end if
                return
            end if
            s_rem=s_rem-s_hit
            I_cell=max(1,I_cell-1)
            u_cur=min(u_edge(1),u_bound*(one+1d-12))
        end if
    end do
    u_back=u_cur
end subroutine electron_trace_piecewise_affine_u_edge_from_cell



! 沿分段仿射 u 特征线回溯整组网格边界。
subroutine electron_trace_piecewise_affine_u_edges(Num_gam_e,u_now_edge,u_edge,a_cell,b_cell,lag,x_back)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e,I_cell
    real(8), intent(in) :: u_now_edge(Num_gam_e+1),u_edge(Num_gam_e+1),a_cell(Num_gam_e),b_cell(Num_gam_e),lag
    real(8), intent(out) :: x_back(Num_gam_e+1)
    real(8) :: u_prev
    I_cell=1
    do I_gam_e=1,Num_gam_e+1
        call electron_find_u_cell_desc_hint(Num_gam_e,u_edge,u_now_edge(I_gam_e),I_cell)
        call electron_trace_piecewise_affine_u_edge_from_cell(Num_gam_e,u_edge,a_cell,b_cell,lag,u_now_edge(I_gam_e),I_cell,u_prev)
        x_back(I_gam_e)=electron_x_from_u(u_prev)
    end do
end subroutine electron_trace_piecewise_affine_u_edges



! 批量回溯多个滞后时间下的分段仿射 u 特征线边界。
subroutine electron_trace_piecewise_affine_u_edges_batch(Num_gam_e,Num_lag,u_now_edge,u_edge,a_cell,b_cell,lag_arr,x_back_batch)
    implicit none
    integer, intent(in) :: Num_gam_e,Num_lag
    integer :: I_lag
    real(8), intent(in) :: u_now_edge(Num_gam_e+1),u_edge(Num_gam_e+1),a_cell(Num_gam_e),b_cell(Num_gam_e),lag_arr(Num_lag)
    real(8), intent(out) :: x_back_batch(Num_gam_e+1,Num_lag)
    do I_lag=1,Num_lag
        call electron_trace_piecewise_affine_u_edges(Num_gam_e,u_now_edge,u_edge,a_cell,b_cell,lag_arr(I_lag),x_back_batch(:,I_lag))
    end do
end subroutine electron_trace_piecewise_affine_u_edges_batch



! 用已回溯边界和源项 PPM 积分完成特征线输运核心更新。
subroutine electron_characteristic_core(Num_gam_e,dDR,x_edge,x_back_batch,source_scale,dF1,dN_x_in,dN_x_out)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_quad,I_gam_e
    real(8), intent(in) :: dDR,x_edge(Num_gam_e+1),source_scale,dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(in) :: x_back_batch(Num_gam_e+1,5)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    real(8) :: ql(Num_gam_e),qr(Num_gam_e),prefix(0:Num_gam_e),x_back(Num_gam_e+1),dx_cur
    real(8) :: dN_source(Num_gam_e),dN_quad(Num_gam_e),qls(Num_gam_e),qrs(Num_gam_e),prefixs(0:Num_gam_e)
    call electron_prepare_conservative_remap_nonuniform(Num_gam_e,x_edge,dN_x_in,ql,qr,prefix)
    x_back=x_back_batch(:,1)
    do I_gam_e=1,Num_gam_e
        dx_cur=max(x_edge(I_gam_e+1)-x_edge(I_gam_e),1d-30)
        dN_x_out(I_gam_e)=(electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge,dN_x_in,ql,qr,prefix,x_back(I_gam_e+1))- &
                           electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge,dN_x_in,ql,qr,prefix,x_back(I_gam_e)))/dx_cur
    end do
    dN_x_out=max(zero,dN_x_out)
    if (source_scale == zero) return
    call electron_prepare_conservative_remap_nonuniform(Num_gam_e,x_edge,dF1,qls,qrs,prefixs)
    dN_source=zero
    do I_quad=1,charint_quad_order
        x_back=x_back_batch(:,I_quad+1)
        do I_gam_e=1,Num_gam_e
            dx_cur=max(x_edge(I_gam_e+1)-x_edge(I_gam_e),1d-30)
            dN_quad(I_gam_e)=(electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge,dF1,qls,qrs,prefixs,x_back(I_gam_e+1))- &
                              electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge,dF1,qls,qrs,prefixs,x_back(I_gam_e)))/dx_cur
        end do
        dN_source=dN_source+charint_quad_weights(I_quad)*max(zero,dN_quad)
    end do
    dN_x_out=max(zero,dN_x_out+dDR*source_scale*dN_source)
end subroutine electron_characteristic_core



! 统一执行仿射或分段仿射 u 冷却下的特征线输运更新。
subroutine electron_characteristic_update(Num_gam_e,dDR,x_edge,cooling_mode,a_u,b_u,gam_e,dEl,R_loc, &
                                            source_scale,dF1,dN_x_in,dN_x_out)
    implicit none
    integer, intent(in) :: Num_gam_e,cooling_mode
    real(8), intent(in) :: dDR,x_edge(Num_gam_e+1),a_u,b_u,gam_e(Num_gam_e),dEl(Num_gam_e)
    real(8), intent(in) :: R_loc,source_scale,dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    real(8) :: lag_arr(5),x_back_batch(Num_gam_e+1,5),u_edge(Num_gam_e+1)
    real(8) :: a_cell(Num_gam_e),b_cell(Num_gam_e)
    integer :: Num_lag,I_quad
    lag_arr(1)=dDR
    do I_quad=1,charint_quad_order
        lag_arr(I_quad+1)=charint_quad_nodes(I_quad)*dDR
    end do
    Num_lag=charint_quad_order+1
    if (source_scale == zero) Num_lag=1
    select case (cooling_mode)
    case (electron_cooling_affine)
        call electron_u_edges_from_x(Num_gam_e,x_edge,u_edge)
        call electron_trace_affine_u_edges_batch(Num_gam_e,Num_lag,u_edge,lag_arr,a_u,b_u,x_back_batch)
    case (electron_cooling_piecewise)
        call electron_build_piecewise_affine_u(Num_gam_e,x_edge,gam_e,dEl,R_loc,u_edge,a_cell,b_cell)
        call electron_trace_piecewise_affine_u_edges_batch(Num_gam_e,Num_lag,u_edge,u_edge,a_cell,b_cell,lag_arr,x_back_batch)
    case default
        error stop 'unknown electron_characteristic_update cooling_mode'
    end select
    call electron_characteristic_core(Num_gam_e,dDR,x_edge,x_back_batch,source_scale,dF1, &
                                           dN_x_in,dN_x_out)
end subroutine electron_characteristic_update




! 积分单个 PPM 抛物线单元在给定区间内的质量。
real(8) function electron_ppm_cell_int(qc,q_left,q_right,cell_lo,d_x,xl,xr)
    implicit none
    real(8), intent(in) :: qc,q_left,q_right,cell_lo,d_x,xl,xr
    real(8) :: xi_l,xi_r,q6,bcoef
    xi_l=(xl-cell_lo)/d_x
    xi_r=(xr-cell_lo)/d_x
    q6=6d0*qc-3d0*(q_left+q_right)
    bcoef=q_right-q_left+q6
    electron_ppm_cell_int=d_x*( &
        q_left*(xi_r-xi_l) + &
        0.5d0*bcoef*(xi_r*xi_r-xi_l*xi_l) - &
        (q6/3d0)*(xi_r*xi_r*xi_r-xi_l*xi_l*xi_l))
end function electron_ppm_cell_int



! 均匀 x 网格源项半步-半拉格朗日输运-源项半步更新。
subroutine electron_semi_lagrangian_step(Num_gam_e,dDR,d_x,face_coeff,dF1,dN_x_in,dN_x_out)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: dDR,d_x,face_coeff(Num_gam_e-1),dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    real(8) :: x_edge(Num_gam_e+1),dN_half(Num_gam_e),ql(Num_gam_e),qr(Num_gam_e),prefix(0:Num_gam_e),dep_l,dep_r,dx_cell
    do I_gam_e=1,Num_gam_e+1
        x_edge(I_gam_e)=dble(I_gam_e-1)*d_x
    end do
    dN_half=max(zero,dN_x_in+0.5d0*dDR*dF1)
    call electron_prepare_conservative_remap_nonuniform(Num_gam_e,x_edge,dN_half,ql,qr,prefix)
    do I_gam_e=1,Num_gam_e
        dx_cell=d_x
        dep_l=x_edge(I_gam_e)+dDR*face_coeff(max(1,I_gam_e-1))
        dep_r=x_edge(I_gam_e+1)+dDR*face_coeff(min(Num_gam_e-1,I_gam_e))
        dep_l=max(x_edge(1),min(dep_l,x_edge(Num_gam_e+1)))
        dep_r=max(x_edge(1),min(dep_r,x_edge(Num_gam_e+1)))
        dN_x_out(I_gam_e)=(electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge,dN_half,ql,qr,prefix,dep_r)- &
                           electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge,dN_half,ql,qr,prefix,dep_l))/dx_cell
    end do
    dN_x_out=max(zero,dN_x_out+0.5d0*dDR*dF1)
end subroutine electron_semi_lagrangian_step



! 用face速度正负分裂执行隐式守恒输运，适用于冷却和SSA加热共存的符号切换。
subroutine electron_fullhide_flux_split_step(Num_gam_e,dDR,d_x,face_speed,dF1,dN_x_in,dN_x_out,close_low_boundary)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: dDR,d_x,face_speed(Num_gam_e-1),dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    logical, intent(in), optional :: close_low_boundary
    logical :: low_boundary_closed
    real(8) :: lambda,denom,a_plus(Num_gam_e-1),a_minus(Num_gam_e-1)
    real(8) :: lower(Num_gam_e),diag(Num_gam_e),upper(Num_gam_e),rhs(Num_gam_e),cprime(Num_gam_e),dprime(Num_gam_e)

    lambda=dDR/d_x
    a_plus=max(face_speed,zero)
    a_minus=min(face_speed,zero)
    lower=zero; diag=one; upper=zero
    rhs=dN_x_in+dDR*dF1
    low_boundary_closed=.false.
    if (present(close_low_boundary)) low_boundary_closed=close_low_boundary

    do I_gam_e=1,Num_gam_e
        if (I_gam_e == 1) then
            if (.not. low_boundary_closed) diag(I_gam_e)=diag(I_gam_e)+lambda*a_plus(1)
        else
            diag(I_gam_e)=diag(I_gam_e)+lambda*a_plus(I_gam_e-1)
            lower(I_gam_e)=-lambda*(-a_minus(I_gam_e-1))
        end if
        if (I_gam_e == Num_gam_e) then
            diag(I_gam_e)=diag(I_gam_e)+lambda*(-a_minus(Num_gam_e-1))
        else
            diag(I_gam_e)=diag(I_gam_e)+lambda*(-a_minus(I_gam_e))
            upper(I_gam_e)=-lambda*a_plus(I_gam_e)
        end if
    end do

    cprime(1)=upper(1)/diag(1)
    dprime(1)=rhs(1)/diag(1)
    do I_gam_e=2,Num_gam_e
        denom=diag(I_gam_e)-lower(I_gam_e)*cprime(I_gam_e-1)
        if (I_gam_e < Num_gam_e) cprime(I_gam_e)=upper(I_gam_e)/denom
        dprime(I_gam_e)=(rhs(I_gam_e)-lower(I_gam_e)*dprime(I_gam_e-1))/denom
    end do
    dN_x_out(Num_gam_e)=dprime(Num_gam_e)
    do I_gam_e=Num_gam_e-1,1,-1
        dN_x_out(I_gam_e)=dprime(I_gam_e)-cprime(I_gam_e)*dN_x_out(I_gam_e+1)
    end do
    dN_x_out=max(zero,dN_x_out)
end subroutine electron_fullhide_flux_split_step

! 执行 GitHub 基线 fullhide 电子输运的上三角隐式冷却和注入更新。
subroutine electron_fullhide_step(Num_gam_e,R_loc,dDR,d_x,dEL_mean,dF1,dN_x_in,dN_x_out)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: R_loc,dDR,d_x
    real(8), intent(in) :: dEL_mean(Num_gam_e-1),dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    real(8) :: face_speed(Num_gam_e-1),up(Num_gam_e-1),principal(Num_gam_e),coupling(Num_gam_e-1),rhs(Num_gam_e)

    face_speed=dEL_mean+one/R_loc/dlog(ten)
    up=-(dDR/d_x)*face_speed
    principal(2:Num_gam_e)=one-up
    principal(1)=principal(2)
    coupling=two*up/(principal(2:Num_gam_e)+principal(1:Num_gam_e-1))
    rhs=(dN_x_in+dDR*dF1)/principal

    dN_x_out(Num_gam_e)=rhs(Num_gam_e)
    do I_gam_e=Num_gam_e-1,1,-1
        dN_x_out(I_gam_e)=max(zero,rhs(I_gam_e)-coupling(I_gam_e)*dN_x_out(I_gam_e+1))
    end do
end subroutine electron_fullhide_step

! CPU版space-time fullhide推进：使用GPU同源离散，但按step-major顺序回代以避免反对角同步开销。
subroutine electron_fullhide_spacetime_sequence(Num_gam_e,Num_step,face_coupling,source_step, &
                                                dN_x_in,dN_x_out)
    !$ use omp_lib
    implicit none
    integer, intent(in) :: Num_gam_e,Num_step
    integer :: I_gam_e,I_step
    real(8), intent(in) :: face_coupling(Num_gam_e-1,Num_step),source_step(Num_gam_e,Num_step)
    real(8), intent(in) :: dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    real(8) :: diag(Num_gam_e,Num_step),upper(Num_gam_e,Num_step),rhs(Num_gam_e,Num_step)
    real(8) :: face_left,face_right,previous_step,upper_neighbor

    diag=one
    upper=zero
    rhs=source_step
    rhs(:,1)=rhs(:,1)+dN_x_in

    do I_step=1,Num_step
        do I_gam_e=1,Num_gam_e
            face_left=zero
            face_right=zero
            if (I_gam_e > 1) face_left=face_coupling(I_gam_e-1,I_step)
            if (I_gam_e < Num_gam_e) face_right=face_coupling(I_gam_e,I_step)

            if (I_gam_e == 1) then
                diag(I_gam_e,I_step)=one+face_right
                upper(I_gam_e,I_step)=-face_right
            else if (I_gam_e == Num_gam_e) then
                diag(I_gam_e,I_step)=one+face_left
            else
                diag(I_gam_e,I_step)=one+face_left
                upper(I_gam_e,I_step)=-face_right
            end if
        end do
    end do

    do I_step=1,Num_step
        do I_gam_e=Num_gam_e,1,-1
            previous_step=zero
            if (I_step > 1) previous_step=rhs(I_gam_e,I_step-1)
            upper_neighbor=zero
            if (I_gam_e < Num_gam_e) upper_neighbor=rhs(I_gam_e+1,I_step)
            rhs(I_gam_e,I_step)=max(zero,(rhs(I_gam_e,I_step)+previous_step &
                                  -upper(I_gam_e,I_step)*upper_neighbor)/diag(I_gam_e,I_step))
        end do
    end do

    dN_x_out=rhs(:,Num_step)
end subroutine electron_fullhide_spacetime_sequence



end module electron_transport_common
