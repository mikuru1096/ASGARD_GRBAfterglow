module electron_transport_common
    use constants
    implicit none
    integer, parameter :: char_order = 4
    real(8), parameter :: cfl_relax = 32d0, rtol_relax = 8d0
    integer, parameter :: cooling_affine = 0, cooling_piecewise = 1
    real(8), parameter :: tiny_u = 1d-300
    real(8), parameter, dimension(char_order) :: qnodes0= [&
        6.9431844202973714d-2, 3.3000947820757187d-1, &
        6.6999052179242813d-1, 9.3056815579702629d-1]
    real(8), parameter, dimension(char_order) :: qweights0= [&
        1.7392742256872693d-1, 3.2607257743127307d-1, &
        3.2607257743127307d-1, 1.7392742256872693d-1]
    real(8), parameter, dimension(char_order) :: qnodes= qnodes0**4, qweights= 4d0*qweights0*qnodes0**3
contains

! 准备隐式迎风输运系数：计算主对角元principal和上三角系数temp1。
! Prepare implicit upwind coefficients: diagonal terms and backward-sweep coupling.
subroutine prepare_implicit_coeffs(ng,diag_base,up,principal,temp1)
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng-1) :: up
    real(8), intent(in) :: diag_base
    real(8), intent(out), dimension(ng) :: principal
    real(8), intent(out), dimension(ng-1) :: temp1

    principal(2:ng)=diag_base-up
    principal(1)=principal(2)
    temp1=up/(principal(2:ng)+principal(1:ng-1))*2d0
end subroutine prepare_implicit_coeffs

! 隐式迎风输运后向回代求解：从高能向低能端回代，解截断为非负。
! Solve the implicit upwind system from high to low energy with a nonnegative projection.
subroutine backward_sweep(ng,temp1,temp2,x)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in), dimension(ng-1) :: temp1
    real(8), intent(in), dimension(ng) :: temp2
    real(8), intent(out), dimension(ng) :: x

    x(ng)=temp2(ng)
    do ig=ng-1,1,-1
        x(ig)=max(0d0,temp2(ig)-temp1(ig)*x(ig+1))
    end do
end subroutine backward_sweep

! 用三个点做二次插值。
! Quadratic interpolation over 3 neighboring samples.
real(8) function quad_interp3(x0,y0,x1,y1,x2,y2,x)
    implicit none
    real(8), intent(in) :: x0,y0,x1,y1,x2,y2,x
    real(8) :: l0,l1,l2
    l0=((x-x1)*(x-x2))/max((x0-x1)*(x0-x2),1d-30)
    l1=((x-x0)*(x-x2))/max((x1-x0)*(x1-x2),1d-30)
    l2=((x-x0)*(x-x1))/max((x2-x0)*(x2-x1),1d-30)
    quad_interp3=y0*l0+y1*l1+y2*l2
end function quad_interp3



! 构造非均匀网格 PPM 单元左右界面值并做单调性限制。
! Build left/right PPM interface states on a nonuniform grid with monotonic limiting.
subroutine ppm_interfaces(ng,xedge,q,ql,qr)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(in), dimension(ng) :: q
    real(8), intent(out), dimension(ng) :: ql,qr
    real(8), dimension(ng) :: xmid,dx
    real(8) :: qmin,qmax,qbar,dq
    do ig=1,ng
        xmid(ig)=0.5d0*(xedge(ig)+xedge(ig+1))
        dx(ig)=xedge(ig+1)-xedge(ig)
    end do
    ql=q
    qr=q
    if (ng >= 3) then
        qr(1)=quad_interp3(xmid(1),q(1),xmid(2),q(2),xmid(3),q(3),xedge(2))
        ql(ng)=quad_interp3(xmid(ng-2),q(ng-2), &
                                                     xmid(ng-1),q(ng-1), &
                                                     xmid(ng),q(ng),xedge(ng))
    end if
    do ig=2,ng-1
        ql(ig)=quad_interp3(xmid(ig-1),q(ig-1), &
                                                   xmid(ig),q(ig), &
                                                   xmid(ig+1),q(ig+1),xedge(ig))
        qr(ig)=quad_interp3(xmid(ig-1),q(ig-1), &
                                                    xmid(ig),q(ig), &
                                                    xmid(ig+1),q(ig+1),xedge(ig+1))
        qmin=min(q(ig-1),min(q(ig),q(ig+1)))
        qmax=max(q(ig-1),max(q(ig),q(ig+1)))
        ql(ig)=max(qmin,min(qmax,ql(ig)))
        qr(ig)=max(qmin,min(qmax,qr(ig)))
        if ((q(ig+1)-q(ig))*(q(ig)-q(ig-1)) <= 0d0) then
            ql(ig)=q(ig)
            qr(ig)=q(ig)
        else
            qbar=0.5d0*(ql(ig)+qr(ig))
            dq=qr(ig)-ql(ig)
            if (dq*(q(ig)-qbar) > dq*dq/6d0) then
                ql(ig)=3d0*q(ig)-2d0*qr(ig)
            else if (dq*(q(ig)-qbar) < -dq*dq/6d0) then
                qr(ig)=3d0*q(ig)-2d0*ql(ig)
            end if
        end if
    end do
    if (ng >= 2) then
        qmin=min(q(1),q(2))
        qmax=max(q(1),q(2))
        ql(1)=max(qmin,min(qmax,ql(1)))
        qr(1)=max(qmin,min(qmax,qr(1)))
        qmin=min(q(ng-1),q(ng))
        qmax=max(q(ng-1),q(ng))
        ql(ng)=max(qmin,min(qmax,ql(ng)))
        qr(ng)=max(qmin,min(qmax,qr(ng)))
    end if
end subroutine ppm_interfaces



! 对 PPM 单元多项式做守恒正性缩放；单元平均值不变，光滑正区间保持原三阶重构。
! Positivity-scale a PPM cell polynomial while preserving the cell average.
subroutine ppm_positive_cell(qc,ql,qr)
    implicit none
    real(8), intent(in) :: qc
    real(8), intent(inout) :: ql,qr
    real(8) :: q6,bcoef,xi_vertex,p_min,theta

    p_min=min(ql,qr)
    q6=6d0*qc-3d0*(ql+qr)
    bcoef=qr-ql+q6
    if (q6 < 0d0) then
        xi_vertex=bcoef/(2d0*q6)
        if (xi_vertex > 0d0 .and. xi_vertex < 1d0) then
            p_min=min(p_min,ql+bcoef*xi_vertex-q6*xi_vertex*xi_vertex)
        end if
    end if
    if (p_min < 0d0) then
        if (qc == 0d0) then
            ql=0d0
            qr=0d0
        else
            theta=qc/(qc-p_min)
            ql=qc+theta*(ql-qc)
            qr=qc+theta*(qr-qc)
        end if
    end if
end subroutine ppm_positive_cell



! 计算非均匀网格 PPM 分段抛物线的前缀积分。
! Accumulate prefix mass for the nonuniform-grid PPM parabolas.
subroutine ppm_prefix(ng,xedge,q,ql,qr,prefix)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(in), dimension(ng) :: q,ql,qr
    real(8), intent(out), dimension(0:ng) :: prefix
    real(8) :: dx
    prefix(0)=0d0
    do ig=1,ng
        dx=xedge(ig+1)-xedge(ig)
        prefix(ig)=prefix(ig-1)+ &
            ppm_cell_int(q(ig),ql(ig),qr(ig), &
                                  xedge(ig),dx,xedge(ig),xedge(ig+1))
    end do
end subroutine ppm_prefix



! 在任意位置求非均匀网格 PPM 前缀积分值。
! Evaluate the PPM prefix mass at an arbitrary coordinate.
real(8) function ppm_eval_prefix(ng,xedge,q,ql,qr,prefix,x_eval)
    implicit none
    integer, intent(in) :: ng
    integer :: ig,left,right,mid
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(in), dimension(ng) :: q,ql,qr
    real(8), intent(in), dimension(0:ng) :: prefix
    real(8), intent(in) :: x_eval
    real(8) :: xa,dx
    xa=max(xedge(1),min(x_eval,xedge(ng+1)))
    if (xa <= xedge(1)) then
        ppm_eval_prefix=0d0
        return
    end if
    if (xa >= xedge(ng+1)) then
        ppm_eval_prefix=prefix(ng)
        return
    end if
    left=1
    right=ng
    do while (left < right)
        mid=(left+right)/2
        if (xedge(mid+1) >= xa) then
            right=mid
        else
            left=mid+1
        end if
    end do
    ig=left
    dx=xedge(ig+1)-xedge(ig)
    ppm_eval_prefix=prefix(ig-1)+ &
        ppm_cell_int(q(ig),ql(ig),qr(ig), &
                              xedge(ig),dx,xedge(ig),xa)
end function ppm_eval_prefix



real(8) function ppm_eval_cell(ng,xedge,q,ql,qr,prefix,x_eval,icell)
    implicit none
    integer, intent(in) :: ng,icell
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(in), dimension(ng) :: q,ql,qr
    real(8), intent(in), dimension(0:ng) :: prefix
    real(8), intent(in) :: x_eval
    real(8) :: xa,dx

    xa=max(xedge(1),min(x_eval,xedge(ng+1)))
    if (xa <= xedge(1)) then
        ppm_eval_cell=0d0
        return
    end if
    if (xa >= xedge(ng+1)) then
        ppm_eval_cell=prefix(ng)
        return
    end if
    dx=xedge(icell+1)-xedge(icell)
    ppm_eval_cell=prefix(icell-1)+ &
        ppm_cell_int(q(icell),ql(icell),qr(icell), &
                                  xedge(icell),dx,xedge(icell),xa)
end function ppm_eval_cell



! 预处理非均匀网格守恒重映射所需的 PPM 界面值和前缀积分。
! Prepare PPM interface states and prefix mass for conservative remapping.
subroutine prepare_remap(ng,xedge,q,ql,qr,prefix)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(in), dimension(ng) :: q
    real(8), intent(out), dimension(ng) :: ql,qr
    real(8), intent(out), dimension(0:ng) :: prefix

    call ppm_interfaces(ng,xedge,q,ql,qr)
    do ig=1,ng
        call ppm_positive_cell(q(ig),ql(ig),qr(ig))
    end do
    call ppm_prefix(ng,xedge,q,ql,qr,prefix)
end subroutine prepare_remap



! 对正的解析源项做守恒指数重建：q(x)=A exp[s(x-xc)]，A由单元平均值确定。
! Reconstruct a positive analytic source as a conservative exponential profile.
subroutine prepare_exp_source(ng,xedge,q,sslope,samp,sprefix)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(in), dimension(ng) :: q
    real(8), intent(out), dimension(ng) :: sslope,samp
    real(8), intent(out), dimension(0:ng) :: sprefix
    real(8), dimension(ng) :: xmid
    real(8) :: dx,slope_l,slope_r,z

    do ig=1,ng
        xmid(ig)=0.5d0*(xedge(ig)+xedge(ig+1))
    end do
    sslope=0d0
    do ig=2,ng-1
        if (q(ig-1) > 0d0 .and. q(ig) > 0d0 .and. q(ig+1) > 0d0) then
            slope_l=(dlog(q(ig))-dlog(q(ig-1)))/(xmid(ig)-xmid(ig-1))
            slope_r=(dlog(q(ig+1))-dlog(q(ig)))/(xmid(ig+1)-xmid(ig))
            if (slope_l*slope_r > 0d0) sslope(ig)=sign(min(abs(slope_l),abs(slope_r)),slope_l)
        end if
    end do
    sprefix(0)=0d0
    do ig=1,ng
        dx=xedge(ig+1)-xedge(ig)
        z=0.5d0*sslope(ig)*dx
        if (abs(z) <= 1d-8) then
            samp(ig)=q(ig)
        else
            samp(ig)=q(ig)*z/dsinh(z)
        end if
        sprefix(ig)=sprefix(ig-1)+q(ig)*dx
    end do
end subroutine prepare_exp_source



real(8) function exp_source_int(cell_lo,cell_hi,amp,slope,xl,xr)
    implicit none
    real(8), intent(in) :: cell_lo,cell_hi,amp,slope,xl,xr
    real(8) :: xmid,x_left,x_right

    xmid=0.5d0*(cell_lo+cell_hi)
    x_left=max(cell_lo,min(xl,cell_hi))
    x_right=max(cell_lo,min(xr,cell_hi))
    if (x_right <= x_left) then
        exp_source_int=0d0
    else if (abs(slope*(cell_hi-cell_lo)) <= 1d-8) then
        exp_source_int=amp*(x_right-x_left)
    else
        exp_source_int=amp*(dexp(slope*(x_right-xmid))- &
                                          dexp(slope*(x_left-xmid)))/slope
    end if
end function exp_source_int



real(8) function exp_source_prefix(ng,xedge,sslope,samp,sprefix,x_eval)
    implicit none
    integer, intent(in) :: ng
    integer :: ig,left,right,mid
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(in), dimension(ng) :: sslope,samp
    real(8), intent(in), dimension(0:ng) :: sprefix
    real(8), intent(in) :: x_eval
    real(8) :: xa

    xa=max(xedge(1),min(x_eval,xedge(ng+1)))
    if (xa <= xedge(1)) then
        exp_source_prefix=0d0
        return
    end if
    if (xa >= xedge(ng+1)) then
        exp_source_prefix=sprefix(ng)
        return
    end if
    left=1
    right=ng
    do while (left < right)
        mid=(left+right)/2
        if (xedge(mid+1) >= xa) then
            right=mid
        else
            left=mid+1
        end if
    end do
    ig=left
    exp_source_prefix=sprefix(ig-1)+ &
        exp_source_int(xedge(ig),xedge(ig+1),samp(ig), &
                                     sslope(ig),xedge(ig),xa)
end function exp_source_prefix



! 将 dN/dx 的单元平均值转换为中心 dN/dgamma，用局部指数形状匹配幂律和指数尾巴。
! Convert cell-averaged dN/dx to centered dN/dgamma with a local exponential shape.
subroutine dnx_dgamma(ng,xedge,gam_e,nx,ndg)
    implicit none
    integer, intent(in) :: ng
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(in), dimension(ng) :: gam_e,nx
    real(8), intent(out), dimension(ng) :: ndg
    real(8), dimension(ng) :: eslope,nxmid
    real(8), dimension(0:ng) :: eprefix

    call prepare_exp_source(ng,xedge,nx,eslope,nxmid,eprefix)
    ndg=nxmid/(gam_e)
end subroutine dnx_dgamma



! 将 x=ln(gamma) 边界转换为 u=1/gamma 边界。
! Convert x=ln(gamma) edges to u=1/gamma edges.
subroutine u_edges(ng,xedge,uedge)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(out), dimension(ng+1) :: uedge
    do ig=1,ng+1
        uedge(ig)=dexp(-xedge(ig))
    end do
end subroutine u_edges



! 将 u=1/gamma 转换回 x=ln(gamma)。
! Convert u=1/gamma back to x=ln(gamma).
real(8) function x_from_u(u)
    implicit none
    real(8), intent(in) :: u
    x_from_u=-dlog(max(u,tiny_u))
end function x_from_u



! 批量回溯多个滞后时间下的全局仿射 u 特征线边界。
! Trace global affine-u characteristic edges for several lag times.
subroutine trace_affine_u(ng,nlag,unow,lag_arr,a_u,b_u,xbatch)
    implicit none
    integer, intent(in) :: ng,nlag
    integer :: ig,ilag
    real(8), intent(in), dimension(ng+1) :: unow
    real(8), intent(in), dimension(nlag) :: lag_arr
    real(8), intent(in) :: a_u,b_u
    real(8), intent(out), dimension(ng+1,nlag) :: xbatch
    real(8), dimension(nlag) :: exp_fac
    real(8) :: uprev,shift
    if (abs(b_u) <= 1d-30) then
        do ilag=1,nlag
            do ig=1,ng+1
                uprev=unow(ig)-a_u*lag_arr(ilag)
                xbatch(ig,ilag)=x_from_u(uprev)
            end do
        end do
        return
    end if
    shift=a_u/b_u
    do ilag=1,nlag
        exp_fac(ilag)=dexp(-b_u*lag_arr(ilag))
    end do
    do ilag=1,nlag
        do ig=1,ng+1
            uprev=(unow(ig)+shift)*exp_fac(ilag)-shift
            xbatch(ig,ilag)=x_from_u(uprev)
        end do
    end do
end subroutine trace_affine_u



! 由冷却率在 u 空间构造分段仿射特征线系数。
! Build piecewise affine-u characteristic coefficients from the cooling rate.
subroutine build_piece_u(ng,xedge,gam_e,dEl,adrate,uedge,acell,bcell)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(in), dimension(ng) :: gam_e,dEl
    real(8), intent(in) :: adrate
    real(8), intent(out), dimension(ng+1) :: uedge
    real(8), intent(out), dimension(ng) :: acell,bcell
    real(8), dimension(ng) :: u_center,w_center
    real(8), dimension(ng+1) :: w_edge
    real(8) :: frac,du,slope
    call u_edges(ng,xedge,uedge)
    do ig=1,ng
        u_center(ig)=1d0/gam_e(ig)
        w_center(ig)=u_center(ig)*dEl(ig)
    end do
    w_edge(1)=w_center(1)
    do ig=2,ng
        if (abs(u_center(ig)-u_center(ig-1)) <= 1d-40) then
            w_edge(ig)=0.5d0*(w_center(ig-1)+w_center(ig))
        else
            frac=(uedge(ig)-u_center(ig-1))/(u_center(ig)-u_center(ig-1))
            w_edge(ig)=w_center(ig-1)+frac*(w_center(ig)-w_center(ig-1))
        end if
    end do
    w_edge(ng+1)=w_center(ng)
    do ig=1,ng
        du=uedge(ig+1)-uedge(ig)
        if (abs(du) <= 1d-40) then
            slope=0d0
        else
            slope=(w_edge(ig+1)-w_edge(ig))/du
        end if
        acell(ig)=w_edge(ig)-slope*uedge(ig)
        bcell(ig)=slope+adrate
    end do
end subroutine build_piece_u



! 借助上次单元索引在递减 u 网格中定位给定 u。
! Locate a u value on the decreasing grid using the previous cell index.
subroutine find_u_cell(ng,uedge,uval,I_cell)
    implicit none
    integer, intent(in) :: ng
    integer, intent(inout) :: I_cell
    real(8), intent(in), dimension(ng+1) :: uedge
    real(8), intent(in) :: uval
    if (uval >= uedge(1)) then
        I_cell=1
        return
    end if
    if (uval <= uedge(ng+1)) then
        I_cell=ng
        return
    end if
    I_cell=max(1,min(ng,I_cell))
    do while (I_cell < ng .and. uval <= uedge(I_cell+1))
        I_cell=I_cell+1
    end do
    do while (I_cell > 1 .and. uval > uedge(I_cell))
        I_cell=I_cell-1
    end do
end subroutine find_u_cell



! 从指定单元开始沿分段仿射 u 特征线回溯单个边界。
! Trace a single edge backward along the piecewise affine-u characteristic.
subroutine trace_piece_edge(ng,uedge,acell,bcell,lag,unow,icell,uback)
    implicit none
    integer, intent(in) :: ng, icell
    integer :: I_cell
    real(8), intent(in), dimension(ng+1) :: uedge
    real(8), intent(in), dimension(ng) :: acell,bcell
    real(8), intent(in) :: lag,unow
    real(8), intent(out) :: uback
    real(8) :: s_rem,ucur,g_cur,ubound,s_hit,acur,bcur,shift,num,den,vel_back
    if (lag <= 0d0) then
        uback=unow
        return
    end if
    I_cell=max(1,min(ng,icell))
    ucur=unow
    s_rem=lag
    do while (s_rem > 1d-14*max(lag,1d0))
        acur=acell(I_cell)
        bcur=bcell(I_cell)
        g_cur=acur+bcur*ucur
        if (abs(g_cur) <= 1d-40) then
            uback=ucur
            return
        end if
        if (g_cur > 0d0) then
            if (I_cell == ng) then
                if (abs(bcur) <= 1d-30) then
                    uback=ucur-acur*s_rem
                else
                    shift=acur/bcur
                    uback=(ucur+shift)*dexp(-bcur*s_rem)-shift
                end if
                return
            end if
            ubound=uedge(I_cell+1)
            s_hit=1d300
            if (abs(ubound-ucur) <= 1d-30*max(1d0,abs(ucur),abs(ubound))) then
                s_hit=0d0
            else if (abs(bcur) <= 1d-30) then
                vel_back=-acur
                if (abs(vel_back) > 1d-40) s_hit=(ubound-ucur)/vel_back
            else
                shift=acur/bcur
                num=ucur+shift
                den=ubound+shift
                if (num > 0d0 .and. den > 0d0) s_hit=dlog(num/den)/bcur
            end if
            if (s_hit < 0d0) s_hit=1d300
            if (s_hit >= s_rem .or. s_hit >= 1d290) then
                if (abs(bcur) <= 1d-30) then
                    uback=ucur-acur*s_rem
                else
                    if (abs(bcur) > 1d-30) shift=acur/bcur
                    uback=(ucur+shift)*dexp(-bcur*s_rem)-shift
                end if
                return
            end if
            s_rem=s_rem-s_hit
            I_cell=min(ng,I_cell+1)
            ucur=max(uedge(ng+1),ubound*(1d0-1d-12))
        else
            if (I_cell == 1) then
                ubound=uedge(1)
                s_hit=1d300
                if (abs(ubound-ucur) <= 1d-30*max(1d0,abs(ucur),abs(ubound))) then
                    s_hit=0d0
                else if (abs(bcur) <= 1d-30) then
                    vel_back=-acur
                    if (abs(vel_back) > 1d-40) s_hit=(ubound-ucur)/vel_back
                else
                    shift=acur/bcur
                    num=ucur+shift
                    den=ubound+shift
                    if (num > 0d0 .and. den > 0d0) s_hit=dlog(num/den)/bcur
                end if
                if (s_hit < 0d0) s_hit=1d300
                if (s_hit >= s_rem .or. s_hit >= 1d290) then
                    if (abs(bcur) <= 1d-30) then
                        uback=ucur-acur*s_rem
                    else
                        if (abs(bcur) > 1d-30) shift=acur/bcur
                        uback=(ucur+shift)*dexp(-bcur*s_rem)-shift
                    end if
                else
                    uback=uedge(1)
                end if
                return
            end if
            ubound=uedge(I_cell)
            s_hit=1d300
            if (abs(ubound-ucur) <= 1d-30*max(1d0,abs(ucur),abs(ubound))) then
                s_hit=0d0
            else if (abs(bcur) <= 1d-30) then
                vel_back=-acur
                if (abs(vel_back) > 1d-40) s_hit=(ubound-ucur)/vel_back
            else
                shift=acur/bcur
                num=ucur+shift
                den=ubound+shift
                if (num > 0d0 .and. den > 0d0) s_hit=dlog(num/den)/bcur
            end if
            if (s_hit < 0d0) s_hit=1d300
            if (s_hit >= s_rem .or. s_hit >= 1d290) then
                if (abs(bcur) <= 1d-30) then
                    uback=ucur-acur*s_rem
                else
                    if (abs(bcur) > 1d-30) shift=acur/bcur
                    uback=(ucur+shift)*dexp(-bcur*s_rem)-shift
                end if
                return
            end if
            s_rem=s_rem-s_hit
            I_cell=max(1,I_cell-1)
            ucur=min(uedge(1),ubound*(1d0+1d-12))
        end if
    end do
    uback=ucur
end subroutine trace_piece_edge



! 批量回溯多个滞后时间下的分段仿射 u 特征线边界。
! Trace piecewise affine-u characteristic edges for several lag times.
subroutine trace_piece_u(ng,nlag,unow,uedge,acell,bcell,lag_arr,xbatch)
    implicit none
    integer, intent(in) :: ng,nlag
    integer :: ig,ilag,I_cell
    real(8), intent(in), dimension(ng+1) :: unow,uedge
    real(8), intent(in), dimension(ng) :: acell,bcell
    real(8), intent(in), dimension(nlag) :: lag_arr
    real(8), intent(out), dimension(ng+1,nlag) :: xbatch
    real(8) :: uprev
    do ilag=1,nlag
        I_cell=1
        do ig=1,ng+1
            call find_u_cell(ng,uedge,unow(ig),I_cell)
            call trace_piece_edge(ng,uedge,acell,bcell,lag_arr(ilag), &
                                                                  unow(ig),I_cell,uprev)
            xbatch(ig,ilag)=x_from_u(uprev)
        end do
    end do
end subroutine trace_piece_u



! 用已回溯边界和源项守恒指数积分完成特征线输运核心更新。
! Apply the characteristic transport update from traced edges and source integrals.
subroutine char_core(ng,dDR,xedge,xbatch,srcscale,dF1,nxin,nxout)
    implicit none
    integer, intent(in) :: ng
    integer :: iq,ig
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(in), dimension(ng) :: dF1,nxin
    real(8), intent(in) :: dDR,srcscale
    real(8), intent(in), dimension(ng+1,char_order+1) :: xbatch
    real(8), intent(out), dimension(ng) :: nxout
    real(8), dimension(ng) :: ql,qr
    real(8), dimension(0:ng) :: prefix
    real(8), dimension(ng+1) :: xback
    real(8) :: dx_cur
    real(8), dimension(ng) :: nsrc, nquad, sslope, samp
    real(8), dimension(0:ng) :: sprefix
    call prepare_remap(ng,xedge,nxin,ql,qr,prefix)
    xback=xbatch(:,1)
    do ig=1,ng
        dx_cur=max(xedge(ig+1)-xedge(ig),1d-30)
        nxout(ig)=(ppm_eval_prefix(ng,xedge,nxin,ql,qr,prefix,xback(ig+1))- &
                           ppm_eval_prefix(ng,xedge,nxin,ql,qr,prefix,xback(ig)))/dx_cur
    end do
    if (srcscale == 0d0) return
    call prepare_exp_source(ng,xedge,dF1,sslope,samp,sprefix)
    nsrc=0d0
    do iq=1,char_order
        xback=xbatch(:,iq+1)
        do ig=1,ng
            dx_cur=max(xedge(ig+1)-xedge(ig),1d-30)
            nquad(ig)=(exp_source_prefix(ng,xedge,sslope,samp, &
                                  sprefix,xback(ig+1))- &
                              exp_source_prefix(ng,xedge,sslope,samp, &
                                  sprefix,xback(ig)))/dx_cur
        end do
        nsrc=nsrc+qweights(iq)*nquad
    end do
    nxout=nxout+dDR*srcscale*nsrc
end subroutine char_core



! 用指定回溯边界直接做无源保守重映射。
! Perform source-free conservative remapping from supplied traced edges.
subroutine remap_edges(ng,xedge,xback,nxin,nxout,close_low)
    implicit none
    integer, intent(in) :: ng
    integer :: iq
    real(8), intent(in), dimension(ng+1) :: xedge,xback
    real(8), intent(in), dimension(ng) :: nxin
    real(8), intent(out), dimension(ng) :: nxout
    logical, intent(in), optional :: close_low
    real(8), dimension(ng+1) :: xwork
    real(8), dimension(ng+1,char_order+1) :: xbatch
    real(8), dimension(ng) :: dF_zero

    xwork=xback
    if (present(close_low)) then
        if (close_low) xwork(1)=xedge(1)
    end if
    xbatch(:,1)=xwork
    do iq=2,char_order+1
        xbatch(:,iq)=xwork
    end do
    dF_zero=0d0
    call char_core(ng,0d0,xedge,xbatch,0d0,dF_zero,nxin,nxout)
end subroutine remap_edges



! 统一执行仿射或分段仿射 u 冷却下的特征线输运更新。
! Dispatch the characteristic update for affine or piecewise affine-u cooling.
subroutine char_update(ng,dDR,xedge,coolmode,a_u,b_u,gam_e,dEl,R_loc, &
                                            srcscale,dF1,nxin,nxout,adrate)
    implicit none
    integer, intent(in) :: ng,coolmode
    real(8), intent(in), dimension(ng+1) :: xedge
    real(8), intent(in), dimension(ng) :: gam_e,dEl
    real(8), intent(in) :: dDR,a_u,b_u
    real(8), intent(in), dimension(ng) :: dF1,nxin
    real(8), intent(in) :: R_loc,srcscale
    real(8), intent(in), optional :: adrate
    real(8), intent(out), dimension(ng) :: nxout
    real(8), dimension(char_order+1) :: lag_arr
    real(8), dimension(ng+1,char_order+1) :: xbatch
    real(8), dimension(ng+1) :: uedge
    real(8) :: bpiece
    real(8), dimension(ng) :: acell,bcell
    integer :: nlag,iq
    lag_arr(1)=dDR
    do iq=1,char_order
        lag_arr(iq+1)=qnodes(iq)*dDR
    end do
    nlag=char_order+1
    if (srcscale == 0d0) nlag=1
    select case (coolmode)
    case (cooling_affine)
        call u_edges(ng,xedge,uedge)
        call trace_affine_u(ng,nlag,uedge,lag_arr,a_u,b_u,xbatch)
    case (cooling_piecewise)
        bpiece=1d0/R_loc
        if (present(adrate)) bpiece=adrate
        call build_piece_u(ng,xedge,gam_e,dEl,bpiece,uedge,acell,bcell)
        call trace_piece_u(ng,nlag,uedge,uedge,acell,bcell,lag_arr,xbatch)
    end select
    call char_core(ng,dDR,xedge,xbatch,srcscale,dF1, &
                                           nxin,nxout)
end subroutine char_update



! 积分单个 PPM 抛物线单元在给定区间内的质量。
! Integrate the mass of a PPM cell polynomial over the requested interval.
real(8) function ppm_cell_int(qc,ql,qr,cell_lo,d_x,xl,xr)
    implicit none
    real(8), intent(in) :: qc,ql,qr,cell_lo,d_x,xl,xr
    real(8) :: xi_l,xi_r,q6,bcoef
    xi_l=(xl-cell_lo)/d_x
    xi_r=(xr-cell_lo)/d_x
    q6=6d0*qc-3d0*(ql+qr)
    bcoef=qr-ql+q6
    ppm_cell_int=d_x*( &
        ql*(xi_r-xi_l) + &
        0.5d0*bcoef*(xi_r*xi_r-xi_l*xi_l) - &
        (q6/3d0)*(xi_r*xi_r*xi_r-xi_l*xi_l*xi_l))
end function ppm_cell_int



! 均匀 x 网格源项半步-半拉格朗日输运-源项半步更新。
! Advance a uniform-x cell average with source half-steps and semi-Lagrangian remap.
subroutine semi_lagrangian_step(ng,dDR,d_x,fcoef,dF1,nxin,nxout)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in), dimension(ng-1) :: fcoef
    real(8), intent(in), dimension(ng) :: dF1,nxin
    real(8), intent(in) :: dDR,d_x
    real(8), intent(out), dimension(ng) :: nxout
    real(8), dimension(ng+1) :: xedge
    real(8), dimension(ng) :: nhalf,ql,qr
    real(8), dimension(0:ng) :: prefix
    real(8) :: dep_l,dep_r,dx
    do ig=1,ng+1
        xedge(ig)=dble(ig-1)*d_x
    end do
    nhalf=max(0d0,nxin+0.5d0*dDR*dF1)
    call prepare_remap(ng,xedge,nhalf,ql,qr,prefix)
    do ig=1,ng
        dx=d_x
        dep_l=xedge(ig)+dDR*fcoef(max(1,ig-1))
        dep_r=xedge(ig+1)+dDR*fcoef(min(ng-1,ig))
        dep_l=max(xedge(1),min(dep_l,xedge(ng+1)))
        dep_r=max(xedge(1),min(dep_r,xedge(ng+1)))
        nxout(ig)=(ppm_eval_prefix(ng,xedge,nhalf,ql,qr,prefix,dep_r)- &
                           ppm_eval_prefix(ng,xedge,nhalf,ql,qr,prefix,dep_l))/dx
    end do
    nxout=max(0d0,nxout+0.5d0*dDR*dF1)
end subroutine semi_lagrangian_step



! 用face速度正负分裂执行隐式守恒输运，适用于冷却和SSA加热共存的符号切换。
! Split face speeds by sign for conservative implicit transport with cooling/SSA heating.
subroutine flux_split_step(ng,dDR,d_x,vface,dF1,nxin,nxout,close_low)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in), dimension(ng-1) :: vface
    real(8), intent(in), dimension(ng) :: dF1,nxin
    real(8), intent(in) :: dDR,d_x
    real(8), intent(out), dimension(ng) :: nxout
    logical, intent(in), optional :: close_low
    logical :: low_closed
    real(8), dimension(ng-1) :: a_plus,a_minus
    real(8) :: lambda,denom
    real(8), dimension(ng) :: lower,diag,upper,rhs,cprime,dprime

    lambda=dDR/d_x
    a_plus=max(vface,0d0)
    a_minus=min(vface,0d0)
    lower=0d0; diag=1d0; upper=0d0
    rhs=nxin+dDR*dF1
    low_closed=.false.
    if (present(close_low)) low_closed=close_low

    do ig=1,ng-1
        diag(ig)=diag(ig)+lambda*(-a_minus(ig))
        upper(ig)=-lambda*a_plus(ig)
        diag(ig+1)=diag(ig+1)+lambda*a_plus(ig)
        lower(ig+1)=lambda*a_minus(ig)
    end do
    if (.not. low_closed) diag(1)=diag(1)+lambda*a_plus(1)
    diag(ng)=diag(ng)+lambda*(-a_minus(ng-1))

    cprime(1)=upper(1)/diag(1)
    dprime(1)=rhs(1)/diag(1)
    do ig=2,ng
        denom=diag(ig)-lower(ig)*cprime(ig-1)
        if (ig < ng) cprime(ig)=upper(ig)/denom
        dprime(ig)=(rhs(ig)-lower(ig)*dprime(ig-1))/denom
    end do
    nxout(ng)=dprime(ng)
    do ig=ng-1,1,-1
        nxout(ig)=dprime(ig)-cprime(ig)*nxout(ig+1)
    end do
    nxout=max(0d0,nxout)
end subroutine flux_split_step

! 非均匀能量坐标上的隐式守恒输运；vface>0 表示向低能端移动。
! Conservative implicit transport on a nonuniform energy coordinate; vface>0 moves down.
subroutine flux_split_nonuniform(ng,dDR,coord_edge,vface,dF1,dN_in,dN_out,close_low)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in), dimension(ng+1) :: coord_edge
    real(8), intent(in), dimension(ng-1) :: vface
    real(8), intent(in), dimension(ng) :: dF1,dN_in
    real(8), intent(in) :: dDR
    real(8), intent(out), dimension(ng) :: dN_out
    logical, intent(in), optional :: close_low
    logical :: low_closed
    real(8), dimension(ng-1) :: a_plus,a_minus
    real(8) :: lambda,denom,dx
    real(8), dimension(ng) :: lower,diag,upper,rhs,cprime,dprime

    a_plus=max(vface,0d0)
    a_minus=min(vface,0d0)
    lower=0d0; diag=1d0; upper=0d0
    rhs=dN_in+dDR*dF1
    low_closed=.false.
    if (present(close_low)) low_closed=close_low

    do ig=1,ng-1
        dx=coord_edge(ig+1)-coord_edge(ig)
        lambda=dDR/dx
        diag(ig)=diag(ig)+lambda*(-a_minus(ig))
        upper(ig)=-lambda*a_plus(ig)

        dx=coord_edge(ig+2)-coord_edge(ig+1)
        lambda=dDR/dx
        diag(ig+1)=diag(ig+1)+lambda*a_plus(ig)
        lower(ig+1)=lambda*a_minus(ig)
    end do
    dx=coord_edge(2)-coord_edge(1)
    if (.not. low_closed) diag(1)=diag(1)+(dDR/dx)*a_plus(1)
    dx=coord_edge(ng+1)-coord_edge(ng)
    diag(ng)=diag(ng)+(dDR/dx)*(-a_minus(ng-1))

    cprime(1)=upper(1)/diag(1)
    dprime(1)=rhs(1)/diag(1)
    do ig=2,ng
        denom=diag(ig)-lower(ig)*cprime(ig-1)
        if (ig < ng) cprime(ig)=upper(ig)/denom
        dprime(ig)=(rhs(ig)-lower(ig)*dprime(ig-1))/denom
    end do
    dN_out(ng)=dprime(ng)
    do ig=ng-1,1,-1
        dN_out(ig)=dprime(ig)-cprime(ig)*dN_out(ig+1)
    end do
    dN_out=max(0d0,dN_out)
end subroutine flux_split_nonuniform

! 非均匀能量坐标上的 integrated 隐式守恒输运；fdisp 是单个合并步内的坐标位移积分。
! Integrated implicit transport on a nonuniform energy coordinate using face displacements.
subroutine flux_seq_nonuniform(ng,coord_edge,fdisp,srcstep, &
                                                            dN_in,dN_out,close_low)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in), dimension(ng+1) :: coord_edge
    real(8), intent(in), dimension(ng-1) :: fdisp
    real(8), intent(in), dimension(ng) :: srcstep,dN_in
    real(8), intent(out), dimension(ng) :: dN_out
    logical, intent(in), optional :: close_low
    logical :: low_closed
    real(8), dimension(ng-1) :: a_plus,a_minus
    real(8) :: denom,dx
    real(8), dimension(ng) :: lower,diag,upper,rhs,cprime,dprime

    a_plus=max(fdisp,0d0)
    a_minus=min(fdisp,0d0)
    lower=0d0; diag=1d0; upper=0d0
    rhs=dN_in+srcstep
    low_closed=.false.
    if (present(close_low)) low_closed=close_low

    do ig=1,ng-1
        dx=coord_edge(ig+1)-coord_edge(ig)
        diag(ig)=diag(ig)+(-a_minus(ig))/dx
        upper(ig)=-a_plus(ig)/dx
    end do
    do ig=2,ng
        dx=coord_edge(ig+1)-coord_edge(ig)
        diag(ig)=diag(ig)+a_plus(ig-1)/dx
        lower(ig)=-(-a_minus(ig-1))/dx
    end do
    dx=coord_edge(2)-coord_edge(1)
    if (.not. low_closed) diag(1)=diag(1)+a_plus(1)/dx
    dx=coord_edge(ng+1)-coord_edge(ng)
    diag(ng)=diag(ng)+(-a_minus(ng-1))/dx

    cprime(1)=upper(1)/diag(1)
    dprime(1)=rhs(1)/diag(1)
    do ig=2,ng
        denom=diag(ig)-lower(ig)*cprime(ig-1)
        if (ig < ng) cprime(ig)=upper(ig)/denom
        dprime(ig)=(rhs(ig)-lower(ig)*dprime(ig-1))/denom
    end do
    dN_out(ng)=dprime(ng)
    do ig=ng-1,1,-1
        dN_out(ig)=dprime(ig)-cprime(ig)*dN_out(ig+1)
    end do
    dN_out=max(0d0,dN_out)
end subroutine flux_seq_nonuniform

! 执行 GitHub 基线 fullhide 电子输运的上三角隐式冷却和注入更新。
! Apply the baseline fullhide upper-triangular cooling and injection update.
subroutine fullhide_step(ng,R_loc,dDR,d_x,dEL_mean,dF1,nxin,nxout)
    implicit none
    integer, intent(in) :: ng
    integer :: ig
    real(8), intent(in) :: R_loc,dDR,d_x
    real(8), intent(in), dimension(ng-1) :: dEL_mean
    real(8), intent(in), dimension(ng) :: dF1,nxin
    real(8), intent(out), dimension(ng) :: nxout
    real(8), dimension(ng-1) :: vface,up,coupling
    real(8), dimension(ng) :: principal,rhs

    vface=dEL_mean+1d0/R_loc
    up=-(dDR/d_x)*vface
    principal(2:ng)=1d0-up
    principal(1)=principal(2)
    coupling=2d0*up/(principal(2:ng)+principal(1:ng-1))
    rhs=(nxin+dDR*dF1)/principal

    nxout(ng)=rhs(ng)
    do ig=ng-1,1,-1
        nxout(ig)=max(0d0,rhs(ig)-coupling(ig)*nxout(ig+1))
    end do
end subroutine fullhide_step

! 由离散谱峰附近三个点估计 log-parabola 峰值频率。
! Estimate the log-parabola peak frequency from the 3 samples around the discrete peak.
real(8) function logparabola_peak(nnu,nu,pnu)
    implicit none
    integer, intent(in) :: nnu
    integer :: inu
    real(8), intent(in), dimension(nnu) :: nu,pnu
    real(8) :: x_l,x_c,x_r,y_l,y_c,y_r,x_peak,denom_peak

    inu=maxloc(pnu, dim=1)
    logparabola_peak=max(nu(inu),tiny(1d0))
    if (inu > 1 .and. inu < nnu) then
        if (pnu(inu-1) > 0d0 .and. pnu(inu) > 0d0 .and. pnu(inu+1) > 0d0) then
            x_l=dlog(nu(inu-1))
            x_c=dlog(nu(inu))
            x_r=dlog(nu(inu+1))
            y_l=dlog(pnu(inu-1))
            y_c=dlog(pnu(inu))
            y_r=dlog(pnu(inu+1))
            denom_peak=y_l-2d0*y_c+y_r
            if (dabs(denom_peak) > tiny(1d0)) then
                x_peak=x_c+0.5d0*(y_l-y_r)*(x_c-x_l)/denom_peak
                x_peak=min(max(x_peak,x_l),x_r)
                logparabola_peak=dexp(x_peak)
            end if
        end if
    end if
end function logparabola_peak

! 合并注入源项与已有人口支撑区，确定能量网格活跃上界。
! Combine injection support and existing population support to choose the active energy top.
integer function gamma_active_hi(ng,dF1,shellpop,src_lo,src_hi,pop_peak)
    implicit none
    integer, intent(in) :: ng,src_lo,src_hi
    integer :: ig,source_hi
    real(8), intent(in), dimension(ng) :: dF1,shellpop
    real(8), intent(in) :: pop_peak
    real(8) :: src_peak,sh_peak,support

    source_hi=src_hi
    src_peak=maxval(dF1)
    if (src_peak > 0d0) then
        source_hi=max(2,src_lo)
        do ig=ng,1,-1
            if (dF1(ig) > 1d-12*src_peak) then
                source_hi=max(source_hi,min(ng,ig+1))
                exit
            end if
        end do
    end if

    gamma_active_hi=max(2,source_hi)
    sh_peak=max(pop_peak,src_peak)
    support=1d-12*sh_peak
    if (sh_peak > 0d0) then
        do ig=ng,1,-1
            if (shellpop(ig) > support) then
                gamma_active_hi=max(gamma_active_hi,min(ng,ig+1))
                exit
            end if
        end do
    end if
end function gamma_active_hi

! 由 χ 方向人口支撑区确定 χ 网格活跃上界。
! Choose the active chi-grid top from population support along chi.
integer function chi_active_hi(nchi,chipop,chipeak)
    implicit none
    integer, intent(in) :: nchi
    integer :: ichi
    real(8), intent(in), dimension(nchi) :: chipop
    real(8), intent(in) :: chipeak

    chi_active_hi=nchi
    if (chipeak > 0d0) then
        chi_active_hi=1
        do ichi=nchi,1,-1
            if (chipop(ichi) > 1d-10*chipeak) then
                chi_active_hi=min(nchi,ichi+1)
                exit
            end if
        end do
    end if
end function chi_active_hi

! 扫描 χ-resolved 能量输运系数，给子步长估计提供最大速度。
! Scan chi-resolved energy transport coefficients for the substep speed estimate.
real(8) function max_xi_chi(ng,nchi,dloss, &
                                           adlogchi,chipop,chipeak,active_hi)
    implicit none
    integer, intent(in) :: ng,nchi,active_hi
    integer :: ichi
    real(8), intent(in), dimension(ng-1,nchi) :: dloss
    real(8), intent(in), dimension(nchi) :: adlogchi, chipop
    real(8), intent(in) :: chipeak

    max_xi_chi=0d0
    if (active_hi > 1) then
        do ichi=1,nchi
            if (chipeak > 0d0) then
                if (chipop(ichi) <= 1d-10*chipeak) cycle
            end if
            max_xi_chi=max(max_xi_chi, &
                maxval(dabs(dloss(1:active_hi-1,ichi)+adlogchi(ichi))))
        end do
        if (max_xi_chi <= 0d0) then
            do ichi=1,nchi
                max_xi_chi=max(max_xi_chi, &
                    maxval(dabs(dloss(1:active_hi-1,ichi)+adlogchi(ichi))))
            end do
        end if
    end if
end function max_xi_chi

end module electron_transport_common
