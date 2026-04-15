module electron_common
    use constants
    use adaptive_resampling_mod, only: adaptive_resampling_log
    implicit none

    integer, parameter :: radiation_resample_threshold = 180
    integer, parameter :: radiation_resample_target = 160
    integer, parameter :: radiation_resample_smoothness = 4

contains

subroutine electron_build_gamma_grid(Num_gam_e,Gam_e_max_max,gam_e)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: Gam_e_max_max
    real(8), intent(out) :: gam_e(Num_gam_e)

    do I_gam_e=1,Num_gam_e
        gam_e(I_gam_e)=3d0*ten**(dlog10(Gam_e_max_max)*(I_gam_e-1)/(Num_gam_e-1))
    end do
end subroutine electron_build_gamma_grid

subroutine electron_log_cell_edges(Num_gam_e,gam_e,x_edge)
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
end subroutine electron_log_cell_edges

subroutine electron_gamma_from_edges(Num_gam_e,x_edge,gam_e)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1)
    real(8), intent(out) :: gam_e(Num_gam_e)

    do I_gam_e=1,Num_gam_e
        gam_e(I_gam_e)=ten**(0.5d0*(x_edge(I_gam_e)+x_edge(I_gam_e+1)))
    end do
end subroutine electron_gamma_from_edges

subroutine electron_mc_slopes_nonuniform(Num_gam_e,x_center,q,slope)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_center(Num_gam_e),q(Num_gam_e)
    real(8), intent(out) :: slope(Num_gam_e)
    real(8) :: dl,dr,dc

    slope=zero
    do I_gam_e=2,Num_gam_e-1
        dl=(q(I_gam_e)-q(I_gam_e-1))/max(x_center(I_gam_e)-x_center(I_gam_e-1),1d-30)
        dr=(q(I_gam_e+1)-q(I_gam_e))/max(x_center(I_gam_e+1)-x_center(I_gam_e),1d-30)
        dc=(q(I_gam_e+1)-q(I_gam_e-1))/max(x_center(I_gam_e+1)-x_center(I_gam_e-1),1d-30)
        slope(I_gam_e)=electron_minmod3(two*dl,dc,two*dr)
    end do
end subroutine electron_mc_slopes_nonuniform

subroutine electron_ppm_interfaces_nonuniform(Num_gam_e,x_edge,q,q_left,q_right)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1),q(Num_gam_e)
    real(8), intent(out) :: q_left(Num_gam_e),q_right(Num_gam_e)
    real(8) :: x_center(Num_gam_e),dx_cell(Num_gam_e),slope(Num_gam_e)
    real(8) :: q_min,q_max,q_bar,dq

    call electron_cell_geometry(Num_gam_e,x_edge,x_center,dx_cell)
    call electron_mc_slopes_nonuniform(Num_gam_e,x_center,q,slope)
    q_left=q
    q_right=q

    do I_gam_e=2,Num_gam_e-1
        q_left(I_gam_e)=q(I_gam_e)-0.5d0*slope(I_gam_e)*dx_cell(I_gam_e)
        q_right(I_gam_e)=q(I_gam_e)+0.5d0*slope(I_gam_e)*dx_cell(I_gam_e)
        q_min=min(q(I_gam_e-1),min(q(I_gam_e),q(I_gam_e+1)))
        q_max=max(q(I_gam_e-1),max(q(I_gam_e),q(I_gam_e+1)))
        q_left(I_gam_e)=max(q_min,min(q_max,q_left(I_gam_e)))
        q_right(I_gam_e)=max(q_min,min(q_max,q_right(I_gam_e)))
        if ((q_right(I_gam_e)-q(I_gam_e))*(q(I_gam_e)-q_left(I_gam_e)) <= zero) then
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
end subroutine electron_ppm_interfaces_nonuniform

subroutine electron_ppm_prefix_nonuniform(Num_gam_e,x_edge,q,q_left,q_right,prefix)
    implicit real(8)(A-H,O-Z)
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

real(8) function electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge,q,q_left,q_right,prefix,x_eval)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: left,right,mid,I_gam_e
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

subroutine electron_conservative_remap_nonuniform(Num_gam_e,x_edge_old,x_edge_new,q_old,q_new)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge_old(Num_gam_e+1),x_edge_new(Num_gam_e+1),q_old(Num_gam_e)
    real(8), intent(out) :: q_new(Num_gam_e)
    real(8) :: q_left(Num_gam_e),q_right(Num_gam_e),prefix(0:Num_gam_e),dx_cell

    call electron_ppm_interfaces_nonuniform(Num_gam_e,x_edge_old,q_old,q_left,q_right)
    call electron_ppm_prefix_nonuniform(Num_gam_e,x_edge_old,q_old,q_left,q_right,prefix)
    do I_gam_e=1,Num_gam_e
        dx_cell=max(x_edge_new(I_gam_e+1)-x_edge_new(I_gam_e),1d-30)
        q_new(I_gam_e)= &
            (electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge_old,q_old,q_left,q_right,prefix,x_edge_new(I_gam_e+1)) - &
             electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge_old,q_old,q_left,q_right,prefix,x_edge_new(I_gam_e)))/dx_cell
    end do
    q_new=max(zero,q_new)
end subroutine electron_conservative_remap_nonuniform

real(8) function electron_loglinear_cell_int(qlog_center,slope,x_center,xa,xb)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: qlog_center,slope,x_center,xa,xb
    real(8) :: alpha,beta

    if (xb <= xa) then
        electron_loglinear_cell_int=zero
        return
    end if

    if (abs(slope) <= 1d-14) then
        electron_loglinear_cell_int=(ten**qlog_center-one)*(xb-xa)
        return
    end if

    alpha=slope
    beta=qlog_center-alpha*x_center
    electron_loglinear_cell_int=(ten**beta)*(ten**(alpha*xb)-ten**(alpha*xa))/(alpha*dlog(ten))-(xb-xa)
end function electron_loglinear_cell_int

subroutine electron_log_prefix_nonuniform(Num_gam_e,x_edge,q,prefix,qlog,slope,x_center)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1),q(Num_gam_e)
    real(8), intent(out) :: prefix(0:Num_gam_e),qlog(Num_gam_e),slope(Num_gam_e),x_center(Num_gam_e)
    real(8) :: dx_cell(Num_gam_e)

    call electron_cell_geometry(Num_gam_e,x_edge,x_center,dx_cell)
    qlog=dlog10(one+max(q,zero))
    call electron_mc_slopes_nonuniform(Num_gam_e,x_center,qlog,slope)

    prefix(0)=zero
    do I_gam_e=1,Num_gam_e
        prefix(I_gam_e)=prefix(I_gam_e-1)+electron_loglinear_cell_int(qlog(I_gam_e),slope(I_gam_e),x_center(I_gam_e), &
                                                                      x_edge(I_gam_e),x_edge(I_gam_e+1))
    end do
end subroutine electron_log_prefix_nonuniform

real(8) function electron_log_prefix_eval_nonuniform(Num_gam_e,x_edge,prefix,qlog,slope,x_center,x_eval)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: left,right,mid,I_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1),prefix(0:Num_gam_e),qlog(Num_gam_e),slope(Num_gam_e),x_center(Num_gam_e),x_eval
    real(8) :: xa

    xa=max(x_edge(1),min(x_eval,x_edge(Num_gam_e+1)))
    if (xa <= x_edge(1)) then
        electron_log_prefix_eval_nonuniform=zero
        return
    end if
    if (xa >= x_edge(Num_gam_e+1)) then
        electron_log_prefix_eval_nonuniform=prefix(Num_gam_e)
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
    electron_log_prefix_eval_nonuniform=prefix(I_gam_e-1)+ &
        electron_loglinear_cell_int(qlog(I_gam_e),slope(I_gam_e),x_center(I_gam_e),x_edge(I_gam_e),xa)
end function electron_log_prefix_eval_nonuniform

subroutine electron_conservative_remap_log_nonuniform(Num_gam_e,x_edge_old,x_edge_new,q_old,q_new)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge_old(Num_gam_e+1),x_edge_new(Num_gam_e+1),q_old(Num_gam_e)
    real(8), intent(out) :: q_new(Num_gam_e)
    real(8) :: prefix(0:Num_gam_e),qlog(Num_gam_e),slope(Num_gam_e),x_center(Num_gam_e),dx_cell

    call electron_log_prefix_nonuniform(Num_gam_e,x_edge_old,q_old,prefix,qlog,slope,x_center)
    do I_gam_e=1,Num_gam_e
        dx_cell=max(x_edge_new(I_gam_e+1)-x_edge_new(I_gam_e),1d-30)
        q_new(I_gam_e)= &
            (electron_log_prefix_eval_nonuniform(Num_gam_e,x_edge_old,prefix,qlog,slope,x_center,x_edge_new(I_gam_e+1)) - &
             electron_log_prefix_eval_nonuniform(Num_gam_e,x_edge_old,prefix,qlog,slope,x_center,x_edge_new(I_gam_e)))/dx_cell
    end do
    q_new=max(zero,q_new)
end subroutine electron_conservative_remap_log_nonuniform

subroutine electron_semi_lagrangian_transport_nonuniform(Num_gam_e,dDR,x_edge,face_coeff,dN_x_in,dN_x_out)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: dDR,x_edge(Num_gam_e+1),face_coeff(Num_gam_e-1),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    real(8) :: vel_face(Num_gam_e+1),dep_l,dep_r,dx_cell,prefix(0:Num_gam_e)
    real(8) :: q_left(Num_gam_e),q_right(Num_gam_e)

    call electron_ppm_interfaces_nonuniform(Num_gam_e,x_edge,dN_x_in,q_left,q_right)
    call electron_ppm_prefix_nonuniform(Num_gam_e,x_edge,dN_x_in,q_left,q_right,prefix)

    vel_face(1)=-face_coeff(1)
    vel_face(2:Num_gam_e)=-face_coeff
    vel_face(Num_gam_e+1)=vel_face(Num_gam_e)

    do I_gam_e=1,Num_gam_e
        dep_l=x_edge(I_gam_e)-vel_face(I_gam_e)*dDR
        dep_r=x_edge(I_gam_e+1)-vel_face(I_gam_e+1)*dDR
        dx_cell=max(x_edge(I_gam_e+1)-x_edge(I_gam_e),1d-30)
        dN_x_out(I_gam_e)= &
            (electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge,dN_x_in,q_left,q_right,prefix,dep_r) - &
             electron_ppm_prefix_eval_nonuniform(Num_gam_e,x_edge,dN_x_in,q_left,q_right,prefix,dep_l))/dx_cell
    end do
    dN_x_out=max(zero,dN_x_out)
end subroutine electron_semi_lagrangian_transport_nonuniform

subroutine electron_semi_lagrangian_step_nonuniform(Num_gam_e,dDR,x_edge,face_coeff,dF1,dN_x_in,dN_x_out)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: dDR,x_edge(Num_gam_e+1),face_coeff(Num_gam_e-1),dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    real(8) :: dN_half(Num_gam_e)

    call electron_apply_source_step(Num_gam_e,0.5d0*dDR,dF1,dN_x_in,dN_half)
    call electron_semi_lagrangian_transport_nonuniform(Num_gam_e,dDR,x_edge,face_coeff,dN_half,dN_x_out)
    call electron_apply_source_step(Num_gam_e,0.5d0*dDR,dF1,dN_x_out,dN_x_out)
end subroutine electron_semi_lagrangian_step_nonuniform

subroutine electron_cell_geometry(Num_gam_e,x_edge,x_center,dx_cell)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1)
    real(8), intent(out) :: x_center(Num_gam_e),dx_cell(Num_gam_e)

    do I_gam_e=1,Num_gam_e
        x_center(I_gam_e)=0.5d0*(x_edge(I_gam_e)+x_edge(I_gam_e+1))
        dx_cell(I_gam_e)=x_edge(I_gam_e+1)-x_edge(I_gam_e)
    end do
end subroutine electron_cell_geometry

subroutine electron_find_high_energy_front(Num_gam_e,x_edge,dN_x,x_cut_high)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e,last_active
    real(8), intent(in) :: x_edge(Num_gam_e+1),dN_x(Num_gam_e)
    real(8), intent(out) :: x_cut_high
    real(8) :: x_center(Num_gam_e),dx_cell(Num_gam_e),peak_val,threshold,q0,q1,qth,frac

    call electron_cell_geometry(Num_gam_e,x_edge,x_center,dx_cell)
    peak_val=maxval(dN_x)
    if (peak_val <= zero) then
        x_cut_high=x_edge(Num_gam_e+1)
        return
    end if

    threshold=max(1d-30,1d-8*peak_val)
    last_active=0
    do I_gam_e=1,Num_gam_e
        if (dN_x(I_gam_e) > threshold) last_active=I_gam_e
    end do

    if (last_active <= 0) then
        x_cut_high=x_edge(Num_gam_e+1)
    else if (last_active >= Num_gam_e) then
        x_cut_high=x_edge(Num_gam_e+1)
    else
        q0=dlog10(one+max(dN_x(last_active),zero))
        q1=dlog10(one+max(dN_x(last_active+1),zero))
        qth=dlog10(one+threshold)
        if (abs(q1-q0) <= 1d-30) then
            x_cut_high=0.5d0*(x_center(last_active)+x_center(last_active+1))
        else
            frac=(qth-q0)/(q1-q0)
            frac=max(zero,min(one,frac))
            x_cut_high=x_center(last_active)+frac*(x_center(last_active+1)-x_center(last_active))
        end if
    end if
end subroutine electron_find_high_energy_front

subroutine electron_find_low_energy_front(Num_gam_e,x_edge,dN_x,x_cut_low)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e,first_active
    real(8), intent(in) :: x_edge(Num_gam_e+1),dN_x(Num_gam_e)
    real(8), intent(out) :: x_cut_low
    real(8) :: x_center(Num_gam_e),dx_cell(Num_gam_e),peak_val,threshold,q0,q1,qth,frac

    call electron_cell_geometry(Num_gam_e,x_edge,x_center,dx_cell)
    peak_val=maxval(dN_x)
    if (peak_val <= zero) then
        x_cut_low=x_edge(1)
        return
    end if

    threshold=max(1d-30,1d-15*peak_val)
    first_active=Num_gam_e+1
    do I_gam_e=1,Num_gam_e
        if (dN_x(I_gam_e) > threshold) then
            first_active=I_gam_e
            exit
        end if
    end do

    if (first_active > Num_gam_e) then
        x_cut_low=x_edge(1)
    else if (first_active <= 1) then
        x_cut_low=x_edge(1)
    else
        q0=dlog10(one+max(dN_x(first_active-1),zero))
        q1=dlog10(one+max(dN_x(first_active),zero))
        qth=dlog10(one+threshold)
        if (abs(q1-q0) <= 1d-30) then
            x_cut_low=0.5d0*(x_center(first_active-1)+x_center(first_active))
        else
            frac=(qth-q0)/(q1-q0)
            frac=max(zero,min(one,frac))
            x_cut_low=x_center(first_active-1)+frac*(x_center(first_active)-x_center(first_active-1))
        end if
    end if
end subroutine electron_find_low_energy_front

subroutine electron_build_moving_monitor(Num_gam_e,x_edge,dN_x,Gam_e_m,Gam_e_c,Gam_e_max, &
                                         weights,window_dex,monitor,x_cut_low,x_cut_high)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: x_edge(Num_gam_e+1),dN_x(Num_gam_e),Gam_e_m,Gam_e_c,Gam_e_max
    real(8), intent(in) :: weights(5),window_dex(3)
    real(8), intent(out) :: monitor(Num_gam_e)
    real(8), intent(out) :: x_cut_low,x_cut_high
    real(8) :: x_center(Num_gam_e),dx_cell(Num_gam_e),qx(Num_gam_e),qxx(Num_gam_e),qlog(Num_gam_e)
    real(8) :: max_qx,max_qxx,x_m,x_c,x_max,w_m,w_c,w_max,w_cut,w_low,arg

    call electron_cell_geometry(Num_gam_e,x_edge,x_center,dx_cell)
    do I_gam_e=1,Num_gam_e
        qlog(I_gam_e)=dlog10(one+max(dN_x(I_gam_e),zero))
    end do
    call electron_mc_slopes_nonuniform(Num_gam_e,x_center,qlog,qx)
    qxx=zero
    do I_gam_e=2,Num_gam_e-1
        qxx(I_gam_e)=two*( &
            (qlog(I_gam_e+1)-qlog(I_gam_e))/max(x_center(I_gam_e+1)-x_center(I_gam_e),1d-30) - &
            (qlog(I_gam_e)-qlog(I_gam_e-1))/max(x_center(I_gam_e)-x_center(I_gam_e-1),1d-30))/ &
            max(x_center(I_gam_e+1)-x_center(I_gam_e-1),1d-30)
    end do

    max_qx=maxval(abs(qx))
    max_qxx=maxval(abs(qxx))
    if (max_qx <= zero) max_qx=one
    if (max_qxx <= zero) max_qxx=one

    x_m=dlog10(max(Gam_e_m,ten**x_edge(1)))
    x_c=dlog10(max(Gam_e_c,ten**x_edge(1)))
    x_max=dlog10(max(Gam_e_max,ten**x_edge(1)))
    call electron_find_low_energy_front(Num_gam_e,x_edge,dN_x,x_cut_low)
    call electron_find_high_energy_front(Num_gam_e,x_edge,dN_x,x_cut_high)
    w_m=max(window_dex(1),1d-4)
    w_c=max(window_dex(2),1d-4)
    w_max=max(window_dex(3),1d-4)
    w_cut=max(0.5d0*w_max,1d-4)
    w_low=max(0.5d0*w_m,1d-4)

    do I_gam_e=1,Num_gam_e
        monitor(I_gam_e)=one
        arg=(x_center(I_gam_e)-x_cut_low)/w_low
        monitor(I_gam_e)=monitor(I_gam_e)+1.5d0*weights(1)*dexp(-0.5d0*arg*arg)
        arg=(x_center(I_gam_e)-x_m)/w_m
        monitor(I_gam_e)=monitor(I_gam_e)+weights(1)*dexp(-0.5d0*arg*arg)
        arg=(x_center(I_gam_e)-x_c)/w_c
        monitor(I_gam_e)=monitor(I_gam_e)+weights(2)*dexp(-0.5d0*arg*arg)
        arg=(x_center(I_gam_e)-x_max)/w_max
        monitor(I_gam_e)=monitor(I_gam_e)+weights(3)*dexp(-0.5d0*arg*arg)
        arg=(x_center(I_gam_e)-x_cut_high)/w_cut
        monitor(I_gam_e)=monitor(I_gam_e)+1.5d0*weights(3)*dexp(-0.5d0*arg*arg)
        monitor(I_gam_e)=monitor(I_gam_e)+weights(4)*abs(qx(I_gam_e))/max_qx
        monitor(I_gam_e)=monitor(I_gam_e)+weights(5)*abs(qxx(I_gam_e))/max_qxx
    end do
end subroutine electron_build_moving_monitor

subroutine electron_anchor_low_energy_edges(Num_gam_e,x_edge_old,x_edge_new,x_cut_low,window_dex)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: head_cells,I_edge,anchor_idx
    real(8), intent(in) :: x_edge_old(Num_gam_e+1),x_cut_low,window_dex
    real(8), intent(inout) :: x_edge_new(Num_gam_e+1)
    real(8) :: x_lo,x_next,x_hi,x_cut,eps_spacing,s,min_head_span

    head_cells=min(9,max(6,Num_gam_e/4))
    anchor_idx=head_cells+1
    if (anchor_idx >= Num_gam_e+1) return

    eps_spacing=1d-8
    x_lo=x_edge_old(1)
    min_head_span=max(0.8d0*max(window_dex,1d-4),2d-2)
    if (x_cut_low <= x_lo+min_head_span) return

    x_next=x_edge_new(anchor_idx+1)
    x_cut=max(min(x_cut_low,x_next-4d0*eps_spacing),x_lo+min_head_span)
    x_hi=min(x_next-2d0*eps_spacing,max(x_cut+1.8d0*max(window_dex,1d-4),x_lo+5d0*eps_spacing))

    if (x_cut >= x_hi-2d0*eps_spacing) return

    x_edge_new(anchor_idx)=x_hi
    x_edge_new(anchor_idx-1)=x_cut
    do I_edge=2,anchor_idx-2
        s=real(I_edge-1,8)/real(anchor_idx-2,8)
        x_edge_new(I_edge)=x_lo+(x_cut-x_lo)*(two*s-s*s)
    end do

    do I_edge=2,anchor_idx
        if (x_edge_new(I_edge) <= x_edge_new(I_edge-1)+eps_spacing) then
            x_edge_new(I_edge)=x_edge_new(I_edge-1)+eps_spacing
        end if
    end do
    if (x_edge_new(anchor_idx) >= x_edge_new(anchor_idx+1)-eps_spacing) then
        x_edge_new(anchor_idx)=x_edge_new(anchor_idx+1)-eps_spacing
    end if
end subroutine electron_anchor_low_energy_edges

subroutine electron_anchor_characteristic_edges(Num_gam_e,x_edge_new,x_break,window_dex)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_edge,anchor_idx
    real(8), intent(in) :: x_break,window_dex
    real(8), intent(inout) :: x_edge_new(Num_gam_e+1)
    real(8) :: eps_spacing,span,x_prev,x_next

    eps_spacing=1d-8
    if (x_break <= x_edge_new(2) .or. x_break >= x_edge_new(Num_gam_e)) return

    anchor_idx=2
    do I_edge=2,Num_gam_e
        if (x_edge_new(I_edge) <= x_break .and. x_break <= x_edge_new(I_edge+1)) then
            anchor_idx=I_edge
            exit
        end if
    end do
    anchor_idx=max(2,min(Num_gam_e-1,anchor_idx))

    span=max(0.45d0*window_dex,5d-3)
    x_prev=x_edge_new(anchor_idx-1)
    x_next=x_edge_new(anchor_idx+2)

    x_edge_new(anchor_idx)=max(x_prev+eps_spacing,min(x_break-span,x_next-2d0*eps_spacing))
    x_edge_new(anchor_idx+1)=max(x_edge_new(anchor_idx)+eps_spacing,min(x_break+span,x_next-eps_spacing))

    do I_edge=2,Num_gam_e+1
        if (x_edge_new(I_edge) <= x_edge_new(I_edge-1)+eps_spacing) then
            x_edge_new(I_edge)=x_edge_new(I_edge-1)+eps_spacing
        end if
    end do
end subroutine electron_anchor_characteristic_edges

subroutine electron_anchor_high_energy_edges(Num_gam_e,x_edge_old,x_edge_new,x_cut_high,window_dex)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: tail_cells,I_edge,anchor_idx
    real(8), intent(in) :: x_edge_old(Num_gam_e+1),x_cut_high,window_dex
    real(8), intent(inout) :: x_edge_new(Num_gam_e+1)
    real(8) :: x_hi,x_prev,x_lo,x_cut,eps_spacing,s

    tail_cells=min(4,max(3,Num_gam_e/8))
    anchor_idx=Num_gam_e+1-tail_cells
    if (anchor_idx <= 1) return

    eps_spacing=1d-8
    x_hi=x_edge_old(Num_gam_e+1)
    x_prev=x_edge_new(anchor_idx-1)
    x_cut=min(max(x_cut_high,x_prev+4d0*eps_spacing),x_hi-3d0*eps_spacing)
    x_lo=max(x_prev+2d0*eps_spacing,min(x_cut-1.2d0*max(window_dex,1d-4),x_hi-5d0*eps_spacing))

    if (x_cut <= x_lo+2d0*eps_spacing) return

    x_edge_new(anchor_idx)=x_lo
    x_edge_new(anchor_idx+1)=x_cut
    do I_edge=anchor_idx+2,Num_gam_e+1
        s=real(I_edge-(anchor_idx+1),8)/real(Num_gam_e-anchor_idx,8)
        x_edge_new(I_edge)=x_cut+(x_hi-x_cut)*s*s
    end do

    do I_edge=anchor_idx+1,Num_gam_e+1
        if (x_edge_new(I_edge) <= x_edge_new(I_edge-1)+eps_spacing) then
            x_edge_new(I_edge)=x_edge_new(I_edge-1)+eps_spacing
        end if
    end do
    if (x_edge_new(Num_gam_e+1) > x_hi) x_edge_new(Num_gam_e+1)=x_hi
end subroutine electron_anchor_high_energy_edges

subroutine electron_equidistribute_edges(Num_gam_e,x_edge_old,monitor,smooth,x_edge_new)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_edge,I_gam_e
    real(8), intent(in) :: x_edge_old(Num_gam_e+1),monitor(Num_gam_e),smooth
    real(8), intent(out) :: x_edge_new(Num_gam_e+1)
    real(8) :: accum(0:Num_gam_e),target_mass,mass_total,frac,x_target,eps_spacing

    accum(0)=zero
    do I_gam_e=1,Num_gam_e
        accum(I_gam_e)=accum(I_gam_e-1)+monitor(I_gam_e)*(x_edge_old(I_gam_e+1)-x_edge_old(I_gam_e))
    end do
    mass_total=accum(Num_gam_e)

    x_edge_new(1)=x_edge_old(1)
    x_edge_new(Num_gam_e+1)=x_edge_old(Num_gam_e+1)
    if (mass_total <= zero) then
        x_edge_new=x_edge_old
        return
    end if

    do I_edge=2,Num_gam_e
        target_mass=mass_total*real(I_edge-1,8)/real(Num_gam_e,8)
        x_target=x_edge_old(I_edge)
        do I_gam_e=1,Num_gam_e
            if (accum(I_gam_e) >= target_mass) then
                frac=(target_mass-accum(I_gam_e-1))/max(accum(I_gam_e)-accum(I_gam_e-1),1d-30)
                x_target=x_edge_old(I_gam_e)+frac*(x_edge_old(I_gam_e+1)-x_edge_old(I_gam_e))
                exit
            end if
        end do
        x_edge_new(I_edge)=smooth*x_edge_old(I_edge)+(one-smooth)*x_target
    end do

    eps_spacing=1d-8
    do I_edge=2,Num_gam_e+1
        if (x_edge_new(I_edge) <= x_edge_new(I_edge-1)+eps_spacing) then
            x_edge_new(I_edge)=x_edge_new(I_edge-1)+eps_spacing
        end if
    end do
    if (x_edge_new(Num_gam_e+1) > x_edge_old(Num_gam_e+1)) x_edge_new(Num_gam_e+1)=x_edge_old(Num_gam_e+1)
end subroutine electron_equidistribute_edges

subroutine electron_rescale_upper_boundary(Num_gam_e,x_edge_old,x_upper_new,x_edge_new)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_edge
    real(8), intent(in) :: x_edge_old(Num_gam_e+1),x_upper_new
    real(8), intent(out) :: x_edge_new(Num_gam_e+1)
    real(8) :: x_lo,x_hi_old,scale,eps_spacing

    x_lo=x_edge_old(1)
    x_hi_old=x_edge_old(Num_gam_e+1)
    x_edge_new(1)=x_lo
    x_edge_new(Num_gam_e+1)=x_upper_new
    if (x_hi_old <= x_lo .or. x_upper_new <= x_lo) then
        x_edge_new=x_edge_old
        return
    end if

    scale=(x_upper_new-x_lo)/(x_hi_old-x_lo)
    do I_edge=2,Num_gam_e
        x_edge_new(I_edge)=x_lo+(x_edge_old(I_edge)-x_lo)*scale
    end do

    eps_spacing=1d-8
    do I_edge=2,Num_gam_e+1
        if (x_edge_new(I_edge) <= x_edge_new(I_edge-1)+eps_spacing) then
            x_edge_new(I_edge)=x_edge_new(I_edge-1)+eps_spacing
        end if
    end do
end subroutine electron_rescale_upper_boundary

subroutine electron_rescale_lower_boundary(Num_gam_e,x_edge_old,x_lower_new,x_edge_new)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_edge
    real(8), intent(in) :: x_edge_old(Num_gam_e+1),x_lower_new
    real(8), intent(out) :: x_edge_new(Num_gam_e+1)
    real(8) :: x_lo_old,x_hi,scale,eps_spacing

    x_lo_old=x_edge_old(1)
    x_hi=x_edge_old(Num_gam_e+1)
    x_edge_new(1)=x_lower_new
    x_edge_new(Num_gam_e+1)=x_hi
    if (x_hi <= x_lo_old .or. x_hi <= x_lower_new) then
        x_edge_new=x_edge_old
        return
    end if

    scale=(x_hi-x_lower_new)/(x_hi-x_lo_old)
    do I_edge=2,Num_gam_e
        x_edge_new(I_edge)=x_hi-(x_hi-x_edge_old(I_edge))*scale
    end do

    eps_spacing=1d-8
    do I_edge=2,Num_gam_e+1
        if (x_edge_new(I_edge) <= x_edge_new(I_edge-1)+eps_spacing) then
            x_edge_new(I_edge)=x_edge_new(I_edge-1)+eps_spacing
        end if
    end do
end subroutine electron_rescale_lower_boundary

real(8) function electron_dnx_powerlaw_integral(A,slope,x_lo,x_hi)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: A,slope,x_lo,x_hi
    real(8) :: expo

    if (x_hi <= x_lo) then
        electron_dnx_powerlaw_integral=zero
        return
    end if

    expo=one-slope
    if (abs(expo) < 1d-12) then
        electron_dnx_powerlaw_integral=A*dlog(ten)*(x_hi-x_lo)
    else
        electron_dnx_powerlaw_integral=A*(ten**(expo*x_hi)-ten**(expo*x_lo))/expo
    end if
end function electron_dnx_powerlaw_integral

subroutine electron_fill_powerlaw_interval(Num_gam_e,gam_e,x_edge,A,slope,gam_lo,gam_hi,dNx)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: gam_e(Num_gam_e),x_edge(Num_gam_e+1),A,slope,gam_lo,gam_hi
    real(8), intent(inout) :: dNx(Num_gam_e)
    real(8) :: x_lo,x_hi,cell_lo,cell_hi,dx_cell

    if (gam_hi <= gam_lo) return
    x_lo=dlog10(gam_lo)
    x_hi=dlog10(gam_hi)

    do I_gam_e=1,Num_gam_e
        cell_lo=max(x_edge(I_gam_e),x_lo)
        cell_hi=min(x_edge(I_gam_e+1),x_hi)
        if (cell_hi > cell_lo) then
            dx_cell=x_edge(I_gam_e+1)-x_edge(I_gam_e)
            dNx(I_gam_e)=dNx(I_gam_e)+electron_dnx_powerlaw_integral(A,slope,cell_lo,cell_hi)/dx_cell
        end if
    end do
end subroutine electron_fill_powerlaw_interval

real(8) function electron_dnx_powerlaw_cutoff_value(x,coeff,slope,Gam_e_max)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: x,coeff,slope,Gam_e_max
    real(8) :: gam,cutoff_factor

    if (coeff <= zero .or. Gam_e_max <= zero) then
        electron_dnx_powerlaw_cutoff_value=zero
        return
    end if

    gam=ten**x
    cutoff_factor=one
    if (gam > Gam_e_max) cutoff_factor=dexp(one-gam/Gam_e_max)
    electron_dnx_powerlaw_cutoff_value=coeff*dlog(ten)*gam**(one-slope)*cutoff_factor
end function electron_dnx_powerlaw_cutoff_value

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

subroutine electron_initial_powerlaw(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,dN_gam_e_1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,gam_e(Num_gam_e)
    real(8), intent(out) :: dN_gam_e_1(Num_gam_e)
    do I_gam_e=1,Num_gam_e
        if (Gam_e_m > Gam_e_c) then
            if (Gam_e_c > Gam_e(I_gam_e) .or. Gam_e_max < Gam_e(I_gam_e)) then
                dN_gam_e_1(I_gam_e)=zero
            else
                Q1=Para_N_e_ini*Gam_e_c
                if (Gam_e_m > Gam_e(I_gam_e)) then
                    dN_gam_e_1(I_gam_e)=Q1*Gam_e(I_gam_e)**(-2)
                else
                    dN_gam_e_1(I_gam_e)=Q1*Gam_e_m**(p-one)*Gam_e(I_gam_e)**(-(p+one))
                end if
            end if
        else
            if (Gam_e_m > Gam_e(I_gam_e) .or. Gam_e_max < Gam_e(I_gam_e)) then
                dN_gam_e_1(I_gam_e)=zero
            else
                Q1=Para_N_e_ini*Gam_e_m**(p-one)
                if (Gam_e_c > Gam_e(I_gam_e)) then
                    dN_gam_e_1(I_gam_e)=Q1*Gam_e(I_gam_e)**(-p)
                else
                    dN_gam_e_1(I_gam_e)=Q1*Gam_e_c*Gam_e(I_gam_e)**(-(p+one))
                end if
            end if
        end if
    end do
end subroutine electron_initial_powerlaw

subroutine electron_initial_powerlaw_exp_cutoff(Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,dN_gam_e_1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: Para_N_e_ini,p,Gam_e_m,Gam_e_c,Gam_e_max,gam_e(Num_gam_e)
    real(8), intent(out) :: dN_gam_e_1(Num_gam_e)
    real(8) :: cutoff_factor

    dN_gam_e_1=zero
    if (Gam_e_max <= zero) return

    do I_gam_e=1,Num_gam_e
        if (Gam_e_m > Gam_e_c) then
            if (Gam_e_c > Gam_e(I_gam_e)) then
                dN_gam_e_1(I_gam_e)=zero
            else if (Gam_e_m > Gam_e(I_gam_e)) then
                cutoff_factor=one
                if (Gam_e(I_gam_e) > Gam_e_max) cutoff_factor=dexp(one-Gam_e(I_gam_e)/Gam_e_max)
                dN_gam_e_1(I_gam_e)=Para_N_e_ini*Gam_e_c*Gam_e(I_gam_e)**(-2)*cutoff_factor
            else
                cutoff_factor=one
                if (Gam_e(I_gam_e) > Gam_e_max) cutoff_factor=dexp(one-Gam_e(I_gam_e)/Gam_e_max)
                dN_gam_e_1(I_gam_e)=Para_N_e_ini*Gam_e_c*Gam_e_m**(p-one)*Gam_e(I_gam_e)**(-(p+one))*cutoff_factor
            end if
        else
            if (Gam_e_m > Gam_e(I_gam_e)) then
                dN_gam_e_1(I_gam_e)=zero
            else if (Gam_e_c > Gam_e(I_gam_e)) then
                cutoff_factor=one
                if (Gam_e(I_gam_e) > Gam_e_max) cutoff_factor=dexp(one-Gam_e(I_gam_e)/Gam_e_max)
                dN_gam_e_1(I_gam_e)=Para_N_e_ini*Gam_e_m**(p-one)*Gam_e(I_gam_e)**(-p)*cutoff_factor
            else
                cutoff_factor=one
                if (Gam_e(I_gam_e) > Gam_e_max) cutoff_factor=dexp(one-Gam_e(I_gam_e)/Gam_e_max)
                dN_gam_e_1(I_gam_e)=Para_N_e_ini*Gam_e_m**(p-one)*Gam_e_c*Gam_e(I_gam_e)**(-(p+one))*cutoff_factor
            end if
        end if
    end do
end subroutine electron_initial_powerlaw_exp_cutoff

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
            coeff_lo=Para_N_e_ini*Gam_e_m**(p-one)
            coeff_hi=coeff_lo*Gam_e_c
            call electron_add_dnx_segment(cell_lo,cell_hi,x_m,x_c,coeff_lo,p,Gam_e_max,seg_sum)
            call electron_add_dnx_segment(cell_lo,cell_hi,x_c,huge_x,coeff_hi,p+one,Gam_e_max,seg_sum)
        end if
        dN_x_1(I_gam_e)=seg_sum/dx_cell
    end do
end subroutine electron_initial_powerlaw_exp_cutoff_edges

subroutine electron_gamma_m_near_two(p,threshold,coeff,temp_gam,Gam_e_max,Gam_e_m)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: p,threshold,coeff,temp_gam,Gam_e_max
    real(8), intent(out) :: Gam_e_m

    Gam_e_m=(p-two)/(p-one)*temp_gam+one
    if (p<threshold .and. p>=two) then
        Gam_e_m=coeff/(one+coeff)*temp_gam+one
    else if (p<two .and. p>one) then
        Gam_e_m=((two-p)/(p-one)*temp_gam*Gam_e_max**(p-two))**(one/(p-one))+one
    end if
end subroutine electron_gamma_m_near_two

subroutine electron_gamma_m_exact(p,temp_gam,Gam_e_max,Gam_e_m)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: p,temp_gam,Gam_e_max
    real(8), intent(out) :: Gam_e_m

    Gam_e_m=(p-two)/(p-one)*temp_gam+one
    if (p>two) then
        Gam_e_m=(p-two)/(p-one)*temp_gam+one
    else if (p<two .and. p>one) then
        Gam_e_m=((two-p)/(p-one)*temp_gam*Gam_e_max**(p-two))**(one/(p-one))+one
    else if (p==two) then
        eps=1d-5
        Gam_e_m=one
        temp=temp_gam/log(Gam_e_max/Gam_e_m)
        do while (abs(one-Gam_e_m/temp)>eps)
            temp=temp_gam/log(Gam_e_max/Gam_e_m)
            if (Gam_e_m>temp) then
                Gam_e_m=0.5d0*(Gam_e_m+temp)
            else
                Gam_e_m=0.5d0*(Gam_e_m+Gam_e_max)
            end if
        end do
    end if
end subroutine electron_gamma_m_exact

subroutine electron_loss_mean(Num_gam_e,dEl,dEL_mean)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: dEl(Num_gam_e)
    real(8), intent(out) :: dEL_mean(Num_gam_e-1)

    dEL_mean=(dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two/dlog(ten)
end subroutine electron_loss_mean

subroutine electron_injection_prefactor(R_loc,dDR,dNe,f_e,Gam_e_m_p,Q)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: R_loc,dDR,dNe,f_e,Gam_e_m_p
    real(8), intent(out) :: Q

    Q=4d0/3d0*pi*(3d0*R_loc**2+dDR*(3d0*R_loc+dDR))*dNe*f_e*Gam_e_m_p
end subroutine electron_injection_prefactor

subroutine electron_build_source_term(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,Q,p,dF1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: gam_e(Num_gam_e),Gam_e_m,Gam_e_max,Q,p
    real(8), intent(out) :: dF1(Num_gam_e)

    dF1=zero
    where(gam_e<Gam_e_max .and. gam_e>Gam_e_m) dF1=Q*gam_e**(-p)*gam_e*dlog(ten)
end subroutine electron_build_source_term

subroutine electron_build_source_term_exp_cutoff(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,Q,p,dF1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: gam_e(Num_gam_e),Gam_e_m,Gam_e_max,Q,p
    real(8), intent(out) :: dF1(Num_gam_e)
    real(8) :: cutoff_factor

    dF1=zero
    if (Gam_e_max <= zero) return

    do I_gam_e=1,Num_gam_e
        if (gam_e(I_gam_e) <= Gam_e_m) cycle
        cutoff_factor=one
        if (gam_e(I_gam_e) > Gam_e_max) cutoff_factor=dexp(one-gam_e(I_gam_e)/Gam_e_max)
        dF1(I_gam_e)=Q*gam_e(I_gam_e)**(-p)*cutoff_factor*gam_e(I_gam_e)*dlog(ten)
    end do
end subroutine electron_build_source_term_exp_cutoff

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

subroutine electron_build_source_term_profile(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,Q,profile,dF1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: gam_e(Num_gam_e),Gam_e_m,Gam_e_max,Q,profile(Num_gam_e)
    real(8), intent(out) :: dF1(Num_gam_e)
    real(8) :: x_edge(Num_gam_e+1),x_lo,x_hi,cell_lo,cell_hi,dx_cell

    dF1=zero
    call electron_log_cell_edges(Num_gam_e,gam_e,x_edge)
    x_lo=dlog10(Gam_e_m)
    x_hi=dlog10(Gam_e_max)

    do I_gam_e=1,Num_gam_e
        cell_lo=max(x_edge(I_gam_e),x_lo)
        cell_hi=min(x_edge(I_gam_e+1),x_hi)
        if (cell_hi > cell_lo) then
            dx_cell=x_edge(I_gam_e+1)-x_edge(I_gam_e)
            dF1(I_gam_e)=Q*profile(I_gam_e)*(cell_hi-cell_lo)/dx_cell
        end if
    end do
end subroutine electron_build_source_term_profile

subroutine electron_source_bounds(Num_gam_e,gam_e,Gam_e_m,Gam_e_max,src_lo,src_hi)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer, intent(out) :: src_lo,src_hi
    integer :: I_gam_e
    real(8), intent(in) :: gam_e(Num_gam_e),Gam_e_m,Gam_e_max
    real(8) :: x_edge(Num_gam_e+1),x_lo,x_hi

    call electron_log_cell_edges(Num_gam_e,gam_e,x_edge)
    x_lo=dlog10(Gam_e_m)
    x_hi=dlog10(Gam_e_max)
    src_lo=Num_gam_e+1
    do I_gam_e=1,Num_gam_e
        if (x_edge(I_gam_e+1) > x_lo) then
            src_lo=I_gam_e
            exit
        end if
    end do

    src_hi=0
    do I_gam_e=Num_gam_e,1,-1
        if (x_edge(I_gam_e) < x_hi) then
            src_hi=I_gam_e
            exit
        end if
    end do

    if (src_hi < src_lo) then
        src_lo=1
        src_hi=0
    end if
end subroutine electron_source_bounds

subroutine electron_prepare_radiation_spectrum(Num_gam_e,gam_e,dN_gam_e,Num_gam_rad,gam_e_rad,dN_gam_e_rad)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer, intent(out) :: Num_gam_rad
    integer :: I_gam_e,first_pos,last_pos,n_active,m_target,n_resampled,info
    integer :: idx(Num_gam_e),out_idx,src_idx,last_added_idx
    real(8), intent(in) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e)
    real(8), intent(out) :: gam_e_rad(Num_gam_e),dN_gam_e_rad(Num_gam_e)

    gam_e_rad = gam_e
    dN_gam_e_rad = dN_gam_e
    Num_gam_rad = Num_gam_e

    first_pos = 0
    do I_gam_e = 1, Num_gam_e
        if (dN_gam_e(I_gam_e) > zero) then
            first_pos = I_gam_e
            exit
        end if
    end do
    if (first_pos == 0) return

    last_pos = first_pos
    do I_gam_e = Num_gam_e, first_pos, -1
        if (dN_gam_e(I_gam_e) > zero) then
            last_pos = I_gam_e
            exit
        end if
    end do

    n_active = last_pos - first_pos + 1
    if (n_active <= radiation_resample_threshold) return

    m_target = min(n_active, radiation_resample_target)
    if (m_target >= n_active) return

    call adaptive_resampling_log(gam_e(first_pos:last_pos), dN_gam_e(first_pos:last_pos), n_active, &
                                 m_target, radiation_resample_smoothness, idx(1:m_target), n_resampled, info)
    if (info /= 0 .or. n_resampled <= 0) return

    out_idx = 0
    if (first_pos > 1) then
        out_idx = out_idx + 1
        gam_e_rad(out_idx) = gam_e(first_pos-1)
        dN_gam_e_rad(out_idx) = dN_gam_e(first_pos-1)
    end if

    out_idx = out_idx + 1
    gam_e_rad(out_idx) = gam_e(first_pos)
    dN_gam_e_rad(out_idx) = dN_gam_e(first_pos)
    last_added_idx = first_pos

    do I_gam_e = 1, n_resampled
        src_idx = first_pos + idx(I_gam_e) - 1
        if (src_idx > first_pos .and. src_idx < last_pos) then
            if (src_idx /= last_added_idx) then
                out_idx = out_idx + 1
                gam_e_rad(out_idx) = gam_e(src_idx)
                dN_gam_e_rad(out_idx) = dN_gam_e(src_idx)
                last_added_idx = src_idx
            end if
        end if
    end do

    if (last_pos > first_pos) then
        out_idx = out_idx + 1
        gam_e_rad(out_idx) = gam_e(last_pos)
        dN_gam_e_rad(out_idx) = dN_gam_e(last_pos)
        last_added_idx = last_pos
    end if

    if (last_pos < Num_gam_e) then
        out_idx = out_idx + 1
        gam_e_rad(out_idx) = gam_e(last_pos+1)
        dN_gam_e_rad(out_idx) = dN_gam_e(last_pos+1)
    end if

    Num_gam_rad = out_idx
end subroutine electron_prepare_radiation_spectrum

subroutine electron_prepare_implicit_coeffs(Num_gam_e,diag_base,up,principal,temp1)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: diag_base,up(Num_gam_e-1)
    real(8), intent(out) :: principal(Num_gam_e),temp1(Num_gam_e-1)

    principal(2:Num_gam_e)=diag_base-up
    principal(1)=principal(2)
    temp1=up/(principal(2:Num_gam_e)+principal(1:Num_gam_e-1))*two
end subroutine electron_prepare_implicit_coeffs

subroutine electron_backward_sweep(Num_gam_e,temp1,temp2,x)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: temp1(Num_gam_e-1),temp2(Num_gam_e)
    real(8), intent(out) :: x(Num_gam_e)

    x(Num_gam_e)=temp2(Num_gam_e)
    do I_gam_e=Num_gam_e-1,1,-1
        x(I_gam_e)=max(zero,temp2(I_gam_e)-temp1(I_gam_e)*x(I_gam_e+1))
    end do
end subroutine electron_backward_sweep

real(8) function electron_minmod(a,b)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: a,b

    if (a*b <= zero) then
        electron_minmod=zero
    else if (abs(a) <= abs(b)) then
        electron_minmod=a
    else
        electron_minmod=b
    end if
end function electron_minmod

real(8) function electron_minmod3(a,b,c)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: a,b,c

    electron_minmod3=electron_minmod(a,electron_minmod(b,c))
end function electron_minmod3

subroutine electron_mc_slopes(Num_gam_e,d_x,q,slope)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: d_x,q(Num_gam_e)
    real(8), intent(out) :: slope(Num_gam_e)
    real(8) :: dl,dr,dc

    slope=zero
    do I_gam_e=2,Num_gam_e-1
        dl=(q(I_gam_e)-q(I_gam_e-1))/d_x
        dr=(q(I_gam_e+1)-q(I_gam_e))/d_x
        dc=0.5d0*(q(I_gam_e+1)-q(I_gam_e-1))/d_x
        slope(I_gam_e)=electron_minmod3(two*dl,dc,two*dr)
    end do
end subroutine electron_mc_slopes

real(8) function electron_cell_linear_integral(qc,slope,xc,xl,xr)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: qc,slope,xc,xl,xr

    electron_cell_linear_integral=qc*(xr-xl)+0.5d0*slope*((xr-xc)**2-(xl-xc)**2)
end function electron_cell_linear_integral

subroutine electron_ppm_interfaces(Num_gam_e,d_x,q,q_left,q_right)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: d_x,q(Num_gam_e)
    real(8), intent(out) :: q_left(Num_gam_e),q_right(Num_gam_e)
    real(8) :: slope(Num_gam_e),q_min,q_max,q_bar,dq

    call electron_mc_slopes(Num_gam_e,d_x,q,slope)
    q_left=q
    q_right=q

    do I_gam_e=2,Num_gam_e-1
        q_left(I_gam_e)=q(I_gam_e)-0.5d0*slope(I_gam_e)*d_x
        q_right(I_gam_e)=q(I_gam_e)+0.5d0*slope(I_gam_e)*d_x
        q_min=min(q(I_gam_e-1),min(q(I_gam_e),q(I_gam_e+1)))
        q_max=max(q(I_gam_e-1),max(q(I_gam_e),q(I_gam_e+1)))
        q_left(I_gam_e)=max(q_min,min(q_max,q_left(I_gam_e)))
        q_right(I_gam_e)=max(q_min,min(q_max,q_right(I_gam_e)))
        if ((q_right(I_gam_e)-q(I_gam_e))*(q(I_gam_e)-q_left(I_gam_e)) <= zero) then
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
end subroutine electron_ppm_interfaces

real(8) function electron_ppm_cell_int(qc,q_left,q_right,cell_lo,d_x,xl,xr)
    implicit real(8)(A-H,O-Z)
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

subroutine electron_linear_prefix_integral(Num_gam_e,d_x,q,slope,prefix)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: d_x,q(Num_gam_e),slope(Num_gam_e)
    real(8), intent(out) :: prefix(0:Num_gam_e)
    real(8) :: cell_lo,cell_hi,xc

    prefix(0)=zero
    do I_gam_e=1,Num_gam_e
        cell_lo=(I_gam_e-1)*d_x
        cell_hi=I_gam_e*d_x
        xc=0.5d0*(cell_lo+cell_hi)
        prefix(I_gam_e)=prefix(I_gam_e-1)+electron_cell_linear_integral(q(I_gam_e),slope(I_gam_e),xc,cell_lo,cell_hi)
    end do
end subroutine electron_linear_prefix_integral

subroutine electron_ppm_prefix(Num_gam_e,d_x,q,q_left,q_right,prefix)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: d_x,q(Num_gam_e),q_left(Num_gam_e),q_right(Num_gam_e)
    real(8), intent(out) :: prefix(0:Num_gam_e)
    real(8) :: cell_lo,cell_hi

    prefix(0)=zero
    do I_gam_e=1,Num_gam_e
        cell_lo=(I_gam_e-1)*d_x
        cell_hi=I_gam_e*d_x
        prefix(I_gam_e)=prefix(I_gam_e-1)+ &
            electron_ppm_cell_int(q(I_gam_e),q_left(I_gam_e),q_right(I_gam_e),cell_lo,d_x,cell_lo,cell_hi)
    end do
end subroutine electron_ppm_prefix

real(8) function electron_prefix_integral_eval(Num_gam_e,d_x,q,slope,prefix,x_eval)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: d_x,q(Num_gam_e),slope(Num_gam_e),prefix(0:Num_gam_e),x_eval
    real(8) :: xa,cell_lo,xc

    xa=max(zero,min(x_eval,Num_gam_e*d_x))
    if (xa <= zero) then
        electron_prefix_integral_eval=zero
        return
    end if
    if (xa >= Num_gam_e*d_x) then
        electron_prefix_integral_eval=prefix(Num_gam_e)
        return
    end if

    I_gam_e=max(1,min(Num_gam_e,int(xa/d_x)+1))
    cell_lo=(I_gam_e-1)*d_x
    xc=cell_lo+0.5d0*d_x
    electron_prefix_integral_eval=prefix(I_gam_e-1)+ &
        electron_cell_linear_integral(q(I_gam_e),slope(I_gam_e),xc,cell_lo,xa)
end function electron_prefix_integral_eval

real(8) function electron_ppm_prefix_eval(Num_gam_e,d_x,q,q_left,q_right,prefix,x_eval)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: d_x,q(Num_gam_e),q_left(Num_gam_e),q_right(Num_gam_e),prefix(0:Num_gam_e),x_eval
    real(8) :: xa,cell_lo

    xa=max(zero,min(x_eval,Num_gam_e*d_x))
    if (xa <= zero) then
        electron_ppm_prefix_eval=zero
        return
    end if
    if (xa >= Num_gam_e*d_x) then
        electron_ppm_prefix_eval=prefix(Num_gam_e)
        return
    end if

    I_gam_e=max(1,min(Num_gam_e,int(xa/d_x)+1))
    cell_lo=(I_gam_e-1)*d_x
    electron_ppm_prefix_eval=prefix(I_gam_e-1)+ &
        electron_ppm_cell_int(q(I_gam_e),q_left(I_gam_e),q_right(I_gam_e),cell_lo,d_x,cell_lo,xa)
end function electron_ppm_prefix_eval

subroutine electron_semi_lagrangian_transport(Num_gam_e,dDR,d_x,face_coeff,dN_x_in,dN_x_out)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    integer :: I_gam_e
    real(8), intent(in) :: dDR,d_x,face_coeff(Num_gam_e-1),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    real(8) :: vel_face(0:Num_gam_e),x_edge(0:Num_gam_e),dep_l,dep_r,prefix(0:Num_gam_e)
    real(8) :: q_left(Num_gam_e),q_right(Num_gam_e)

    call electron_ppm_interfaces(Num_gam_e,d_x,dN_x_in,q_left,q_right)
    call electron_ppm_prefix(Num_gam_e,d_x,dN_x_in,q_left,q_right,prefix)
    do I_gam_e=0,Num_gam_e
        x_edge(I_gam_e)=I_gam_e*d_x
    end do

    vel_face(0)=-face_coeff(1)
    vel_face(1:Num_gam_e-1)=-face_coeff
    vel_face(Num_gam_e)=vel_face(Num_gam_e-1)

    do I_gam_e=1,Num_gam_e
        dep_l=x_edge(I_gam_e-1)-vel_face(I_gam_e-1)*dDR
        dep_r=x_edge(I_gam_e)-vel_face(I_gam_e)*dDR
        dN_x_out(I_gam_e)= &
            (electron_ppm_prefix_eval(Num_gam_e,d_x,dN_x_in,q_left,q_right,prefix,dep_r)- &
             electron_ppm_prefix_eval(Num_gam_e,d_x,dN_x_in,q_left,q_right,prefix,dep_l))/d_x
    end do
    dN_x_out=max(zero,dN_x_out)
end subroutine electron_semi_lagrangian_transport

subroutine electron_apply_source_step(Num_gam_e,dDR,dF1,dN_x_in,dN_x_out)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: dDR,dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)

    dN_x_out=max(zero,dN_x_in+dDR*dF1)
end subroutine electron_apply_source_step

subroutine electron_semi_lagrangian_step(Num_gam_e,dDR,d_x,face_coeff,dF1,dN_x_in,dN_x_out)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: dDR,d_x,face_coeff(Num_gam_e-1),dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    real(8) :: dN_half(Num_gam_e)

    call electron_apply_source_step(Num_gam_e,0.5d0*dDR,dF1,dN_x_in,dN_half)
    call electron_semi_lagrangian_transport(Num_gam_e,dDR,d_x,face_coeff,dN_half,dN_x_out)
    call electron_apply_source_step(Num_gam_e,0.5d0*dDR,dF1,dN_x_out,dN_x_out)
end subroutine electron_semi_lagrangian_step

subroutine electron_fullhide_step(Num_gam_e,R_loc,dDR,d_x,dEL_mean,dF1,dN_x_in,dN_x_out)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: R_loc,dDR,d_x
    real(8), intent(in) :: dEL_mean(Num_gam_e-1),dF1(Num_gam_e),dN_x_in(Num_gam_e)
    real(8), intent(out) :: dN_x_out(Num_gam_e)
    real(8) :: principal(Num_gam_e),up(Num_gam_e-1),temp1(Num_gam_e-1),temp2(Num_gam_e),temp3(Num_gam_e-1)

    temp3=dEL_mean+one/R_loc/dlog(ten)
    up=-(dDR/d_x)*temp3
    call electron_prepare_implicit_coeffs(Num_gam_e,one,up,principal,temp1)
    temp2=dN_x_in/principal
    temp2=temp2+dDR*dF1/principal
    call electron_backward_sweep(Num_gam_e,temp1,temp2,dN_x_out)
end subroutine electron_fullhide_step

subroutine electron_max_relative_error(Num_gam_e,x_ref,x_trial,error_max)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_e
    real(8), intent(in) :: x_ref(Num_gam_e),x_trial(Num_gam_e)
    real(8), intent(out) :: error_max
    integer :: I_gam_e

    error_max=zero
    do I_gam_e=1,Num_gam_e
        denom=max(abs(x_ref(I_gam_e)),1d-99)
        error_max=max(error_max,abs(x_trial(I_gam_e)-x_ref(I_gam_e))/denom)
    end do
end subroutine electron_max_relative_error

subroutine electron_initial_density(A_star,dNe_ISM,R_ini,R_start,R0,dNe,Para_N_e_ini)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: A_star,dNe_ISM,R_ini,R_start,R0
    real(8), intent(out) :: dNe,Para_N_e_ini

    if (A_star > zero) then
        dNe_wind=A_star*3.0d35/R_start**2
        Para_N_e_ini=4d0*pi*R_ini*A_star*3.0d35
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
        Para_N_e_ini=4d0/3d0*pi*R_ini**3*dNe_ISM
    end if

    if (R_start < R0) then
        dNe=A_star*3.0d35/R0**2*4d0
        Para_N_e_ini=4d0/3d0*pi*R_ini**3*dNe_ISM
    end if
end subroutine electron_initial_density

subroutine electron_external_density(A_star,dNe_ISM,R_loc,R0,R_tr,f_jump,f_wide,apply_jump,dNe)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: apply_jump
    real(8), intent(in) :: A_star,dNe_ISM,R_loc,R0,R_tr,f_jump,f_wide
    real(8), intent(out) :: dNe

    if (A_star > zero) then
        dNe_wind=A_star*3.0d35/R_loc**2
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
        if (apply_jump /= 0) then
            dNe=dNe_ISM*(1d0+(f_jump-1d0)*exp(-(log10(R_loc)-log10(R_tr))**2/(2d0*f_wide*f_wide)))
        end if
    end if

    if (R_loc < R0) then
        dNe=A_star*3.0d35/R0**2
    end if
end subroutine electron_external_density


subroutine electron_ppm_point_values_nonuniform(Num_src,x_src_edge, &
     q_src,Num_tgt,x_tgt,q_tgt)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_src, Num_tgt
    real(8), intent(in) :: x_src_edge(Num_src+1), q_src(Num_src)
    real(8), intent(in) :: x_tgt(Num_tgt)
    real(8), intent(out) :: q_tgt(Num_tgt)
    real(8) :: q_left(Num_src), q_right(Num_src)
    real(8) :: xi, coeff_c
    integer :: i_tgt, i_src, lo, hi, mid

    call electron_ppm_interfaces_nonuniform(Num_src,x_src_edge, &
         q_src,q_left,q_right)

    do i_tgt = 1, Num_tgt
        q_tgt(i_tgt) = zero
    end do

    do i_tgt = 1, Num_tgt
        if (x_tgt(i_tgt) < x_src_edge(1)) cycle
        if (x_tgt(i_tgt) > x_src_edge(Num_src+1)) cycle
        ! binary search for source cell
        lo = 1
        hi = Num_src
        do while (lo < hi)
            mid = (lo + hi) / 2
            if (x_tgt(i_tgt) > x_src_edge(mid+1)) then
                lo = mid + 1
            else
                hi = mid
            end if
        end do
        i_src = lo
        ! PPM parabola evaluation in cell i_src
        xi = (x_tgt(i_tgt) - x_src_edge(i_src)) &
             / (x_src_edge(i_src+1) - x_src_edge(i_src))
        coeff_c = q_src(i_src) &
             - 0.5d0*(q_left(i_src) + q_right(i_src))
        q_tgt(i_tgt) = q_left(i_src) &
             + xi*(q_right(i_src) - q_left(i_src) + 6d0*coeff_c) &
             + 6d0*coeff_c*xi*(1d0 - xi)
    end do
end subroutine electron_ppm_point_values_nonuniform
end module electron_common
