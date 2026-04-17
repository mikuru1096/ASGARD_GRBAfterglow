module get_Y
  use constants
  use electron_common, only: electron_ppm_interfaces_nonuniform
  private

    public :: besselk, get_syn, get_syn_simpson, get_syn_selected, get_forward_cooling, get_nu_a
    public :: get_nu_a_nonuniform
    public :: get_syn_adaptive
    public :: get_SSA_numerical, get_IC_numerical, get_Y_Nakar, get_Y_Fan

contains

subroutine first_greater_monotonic(arr,n,target,idx)
implicit none
integer, intent(in) :: n
real(8), intent(in) :: arr(n),target
integer, intent(out) :: idx
integer :: left,right,mid

    if (arr(1) > target) then
        idx=1
        return
    end if
    if (arr(n) <= target) then
        idx=n+1
        return
    end if

    left=1
    right=n
    do while (left < right)
        mid=(left+right)/2
        if (arr(mid) > target) then
            right=mid
        else
            left=mid+1
        end if
    end do
    idx=left
end subroutine first_greater_monotonic

subroutine first_greater_monotonic_window(arr,n,start_idx,target,idx)
implicit none
integer, intent(in) :: n,start_idx
real(8), intent(in) :: arr(n),target
integer, intent(out) :: idx
integer :: left,right,mid

    if (arr(start_idx) > target) then
        idx=start_idx
        return
    end if
    if (arr(n) <= target) then
        idx=n+1
        return
    end if

    left=start_idx
    right=n
    do while (left < right)
        mid=(left+right)/2
        if (arr(mid) > target) then
            right=mid
        else
            left=mid+1
        end if
    end do
    idx=left
end subroutine first_greater_monotonic_window

!****************************************************************************************
!**************************** Besselk function interpolation ****************************
!****************************************************************************************
function besselk(var)
    REAL(8) :: val,var,besselk
    integer :: i
    REAL(8),dimension(100),parameter :: minus_theta= &
                 (/1.0000000e-04,1.1233240e-04,1.2618569e-04,1.4174742e-04,1.5922828e-04,1.7886495e-04,2.0092330e-04,& 
                 2.2570197e-04,2.5353645e-04,2.8480359e-04,3.1992671e-04,3.5938137e-04,4.0370173e-04,4.5348785e-04,&
                 5.0941380e-04,5.7223677e-04,6.4280731e-04,7.2208090e-04,8.1113083e-04,9.1116276e-04,1.0235310e-03,&
                 1.1497570e-03,1.2915497e-03,1.4508288e-03,1.6297508e-03,1.8307383e-03,2.0565123e-03,2.3101297e-03,&
                 2.5950242e-03,2.9150531e-03,3.2745492e-03,3.6783798e-03,4.1320124e-03,4.6415888e-03,5.2140083e-03,&
                 5.8570208e-03,6.5793322e-03,7.3907220e-03,8.3021757e-03,9.3260335e-03,1.0476158e-02,1.1768120e-02,&
                 1.3219411e-02,1.4849683e-02,1.6681005e-02,1.8738174e-02,2.1049041e-02,2.3644894e-02,2.6560878e-02,&
                 2.9836472e-02,3.3516027e-02,3.7649358e-02,4.2292429e-02,4.7508102e-02,5.3366992e-02,5.9948425e-02,&
                 6.7341507e-02,7.5646333e-02,8.4975344e-02,9.5454846e-02,1.0722672e-01,1.2045035e-01,1.3530478e-01,&
                 1.5199111e-01,1.7073526e-01,1.9179103e-01,2.1544347e-01,2.4201283e-01,2.7185882e-01,3.0538555e-01,&
                 3.4304693e-01,3.8535286e-01,4.3287613e-01,4.8626016e-01,5.4622772e-01,6.1359073e-01,6.8926121e-01,&
                 7.7426368e-01,8.6974900e-01,9.7700996e-01,1.0974988e+00,1.2328467e+00,1.3848864e+00,1.5556761e+00,&
                 1.7475284e+00,1.9630407e+00,2.2051307e+00,2.4770764e+00,2.7825594e+00,3.1257158e+00,3.5111917e+00,&
                 3.9442061e+00,4.4306215e+00,4.9770236e+00,5.5908102e+00,6.2802914e+00,7.0548023e+00,7.9248290e+00,&
                 8.9021509e+00,1.0000000e+01/)
    REAL(8),dimension(100),parameter :: besselk2= &
                 (/2.0000000e+08,1.5849658e+08,1.2560583e+08,9.9540471e+07,7.8884121e+07,6.2514316e+07,4.9541527e+07,&
                 3.9260813e+07,3.1113522e+07,2.4656934e+07,1.9540199e+07,1.5485273e+07,1.2271814e+07,9.7252027e+06,&
                 7.7070567e+06,6.1077105e+06,4.8402560e+06,3.8358200e+06,3.0398217e+06,2.4090066e+06,1.9090964e+06,&
                 1.5129262e+06,1.1989680e+06,9.5016153e+05,7.5298666e+05,5.9672895e+05,4.7289738e+05,3.7476298e+05,&
                 2.9699315e+05,2.3536189e+05,1.8652017e+05,1.4781394e+05,1.1713992e+05,9.2831277e+04,7.3567095e+04,&
                 5.8300561e+04,4.6202094e+04,3.6614266e+04,2.9016076e+04,2.2994640e+04,1.8222755e+04,1.4441118e+04,&
                 1.1444235e+04,9.0692572e+03,7.1871275e+03,5.6955719e+03,4.5135397e+03,3.5767994e+03,2.8344487e+03,&
                 2.2461486e+03,1.7799308e+03,1.4104612e+03,1.1176629e+03,8.8562540e+02,7.0173970e+02,5.5601353e+02,&
                 4.4052817e+02,3.4900815e+02,2.7648028e+02,2.1900342e+02,1.7345426e+02,1.3735766e+02,1.0875212e+02,&
                 8.6083185e+01,6.8119011e+01,5.3883384e+01,4.2602693e+01,3.3663886e+01,2.6581150e+01,2.0969514e+01,&
                 1.6523922e+01,1.3002650e+01,1.0214166e+01,8.0067134e+00,6.2600579e+00,4.8789400e+00,3.7878895e+00,&
                 2.9271104e+00,2.2492185e+00,1.7166521e+00,1.2996171e+00,9.7445621e-01,7.2235470e-01,5.2831390e-01,&
                 3.8033947e-01,2.6880181e-01,1.8593568e-01,1.2545279e-01,8.2245762e-02,5.2165420e-02,3.1855005e-02,&
                 1.8626561e-02,1.0365680e-02,5.4526367e-03,2.6905368e-03,1.2347515e-03,5.2199962e-04,2.0112029e-04,&
                 6.9778918e-05,2.1509817e-05/)

   besselk=0.0d0
   
   do i=1,99
       if (var > minus_theta(i) .and. var <= minus_theta(i+1)) then
           val = (var-minus_theta(i))/(minus_theta(i+1)-minus_theta(i))
           besselk = besselk2(i)+val*(besselk2(i+1)-besselk2(i))
           exit
       end if
   end do

   return
end function besselk


!****************************************************************************************
!************************** get syn power and number density ****************************
!****************************************************************************************
subroutine get_syn(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                       P_syn,Seed_syn)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) ::P_syn(Num_nu),Seed_syn(Num_nu)

real(8),allocatable,dimension (:) :: dN1,ddN
allocate (dN1(Num_gam_e),ddN(Num_gam_e-1))

    factor=(3.62d0/pi)**2
    Temp_syn=dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    Rariv2=R_loc*R_loc
    dN1=dN_gam_e/(gam_e*gam_e)
    ddN=dN1(1:Num_gam_e-1)-dN1(2:Num_gam_e)
    
    !$ call omp_set_dynamic(.true.)
    !$OMP PARALLEL num_threads(n_threads), private(I_nu,I_gam_e,V_cal,dInteg,Tau,Vc,x,Fx,P_v, &
    !$OMP& gam_e_mean2,dN,dgam_e,ratio_v)
    !$OMP DO SCHEDULE(STATIC)
    do I_nu=1,Num_nu
       V_cal=V_seed(I_nu)
       dInteg=zero
       Tau=zero
       do I_gam_e=1,Num_gam_e-1
          gam_e_mean2=(gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0
          Vc=(4.2d6)*gam_e_mean2*DB !Which is $\nu_c$
          x=V_cal/Vc !Which is ($\nu/\nu_c$)
          ratio_v=Vc/V_cal
          Fx=1.81d0*dexp(-x)/dsqrt(ratio_v**(2d0/3d0)+factor) !Approximate function of synchrotron radiation spectrum
!         Fx=2.149d0*x**(one/3.0d0)*dexp(-x) !!Another approximate function
          dN=(dN_gam_e(I_gam_e)+dN_gam_e(I_gam_e+1))/two
          dgam_e=gam_e(I_gam_e+1)-gam_e(I_gam_e)
          dInteg=dInteg+dN*Fx*dgam_e
          !====================  [SSA]  ======================
          Tau=Tau+gam_e_mean2*ddN(I_gam_e)*Fx
          !Synchrotron self absorption effect
          !===================================================
       end do
       P_v=Temp_syn*DB*dInteg ! with units in erg/Hz/s
       Tau=1.046d4*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal) !Synchrotron self absorption effect
       if ((Tau-1d-4) < 1d-5) Tau=1d-4
       P_syn(I_nu)=P_v*(one-dexp(-Tau))/Tau !Radiation transfer equation for the emission-absorption plasma
       Seed_syn(I_nu)=P_syn(I_nu)/(Rariv2*V_cal)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    temp_para=4d0*pi*Para_c*Para_h
    Seed_syn=Seed_syn/temp_para

deallocate (dN1,ddN)
end subroutine get_syn

real(8) function electron_syn_fx(gam,V_cal,DB,factor)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: gam,V_cal,DB,factor
real(8) :: Vc,x,ratio_v

    Vc=(4.2d6)*gam*gam*DB
    x=V_cal/Vc
    ratio_v=Vc/V_cal
    electron_syn_fx=1.81d0*dexp(-x)/dsqrt(ratio_v**(2d0/3d0)+factor)
end function electron_syn_fx

real(8) function electron_linear_interp(x0,x1,y0,y1,x)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,y0,y1,x

    if (x1 <= x0) then
        electron_linear_interp=0.5d0*(y0+y1)
    else
        electron_linear_interp=y0+(y1-y0)*(x-x0)/(x1-x0)
    end if
end function electron_linear_interp

real(8) function electron_syn_integrand_x(x,x0,x1,dN0,dN1,V_cal,DB,factor)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x,x0,x1,dN0,dN1,V_cal,DB,factor
real(8) :: gam,dN

    gam=dexp(x)
    dN=max(zero,electron_linear_interp(x0,x1,dN0,dN1,x))
    electron_syn_integrand_x=dN*electron_syn_fx(gam,V_cal,DB,factor)*gam
end function electron_syn_integrand_x

real(8) function electron_powerlaw_interp(v0,v1,y0,y1,v)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: v0,v1,y0,y1,v
real(8) :: slope

    if (v <= v0) then
        electron_powerlaw_interp=y0
        return
    end if
    if (v >= v1) then
        electron_powerlaw_interp=y1
        return
    end if

    if (v1 <= v0) then
        electron_powerlaw_interp=0.5d0*(y0+y1)
    else if (y0 > zero .and. y1 > zero) then
        slope=dlog(y1/y0)/dlog(v1/v0)
        electron_powerlaw_interp=y0*(v/v0)**slope
    else
        electron_powerlaw_interp=y0+(y1-y0)*(v-v0)/(v1-v0)
    end if
end function electron_powerlaw_interp

subroutine electron_log_gauss2_interval(v0,v1,v_g1,v_g2,w_g1,w_g2)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: v0,v1
real(8), intent(out) :: v_g1,v_g2,w_g1,w_g2
real(8) :: x0,x1,xm,dx,w2

    x0=dlog(v0)
    x1=dlog(v1)
    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    w2=one/dsqrt(3d0)
    v_g1=dexp(xm-dx*w2)
    v_g2=dexp(xm+dx*w2)
    w_g1=dx*v_g1
    w_g2=dx*v_g2
end subroutine electron_log_gauss2_interval

real(8) function electron_integrate_powerlaw_segment(v0,v1,y0,y1)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: v0,v1,y0,y1
real(8) :: vg1,vg2,wg1,wg2

    if (v1 <= v0) then
        electron_integrate_powerlaw_segment=zero
        return
    end if

    call electron_log_gauss2_interval(v0,v1,vg1,vg2,wg1,wg2)
    electron_integrate_powerlaw_segment= &
        wg1*electron_powerlaw_interp(v0,v1,y0,y1,vg1)+ &
        wg2*electron_powerlaw_interp(v0,v1,y0,y1,vg2)
end function electron_integrate_powerlaw_segment

real(8) function electron_ssa_segment(v0,v1,seed0,seed1,sigma_prefactor,mode,Cyclotron_nu,V_uplim)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: mode
real(8), intent(in) :: v0,v1,seed0,seed1,sigma_prefactor,Cyclotron_nu,V_uplim
real(8) :: vg1,vg2,wg1,wg2,seed_loc,sigma_loc

    if (v1 <= v0) then
        electron_ssa_segment=zero
        return
    end if

    call electron_log_gauss2_interval(v0,v1,vg1,vg2,wg1,wg2)
    electron_ssa_segment=zero

    seed_loc=electron_powerlaw_interp(v0,v1,seed0,seed1,vg1)
    if (mode == 1) then
        sigma_loc=sigma_prefactor*vg1**(-5d0/3d0)
    else
        sigma_loc=sigma_prefactor*(Cyclotron_nu/vg1)*dexp(-vg1/V_uplim)
    end if
    electron_ssa_segment=electron_ssa_segment+wg1*sigma_loc*seed_loc*para_h*vg1*para_c

    seed_loc=electron_powerlaw_interp(v0,v1,seed0,seed1,vg2)
    if (mode == 1) then
        sigma_loc=sigma_prefactor*vg2**(-5d0/3d0)
    else
        sigma_loc=sigma_prefactor*(Cyclotron_nu/vg2)*dexp(-vg2/V_uplim)
    end if
    electron_ssa_segment=electron_ssa_segment+wg2*sigma_loc*seed_loc*para_h*vg2*para_c
end function electron_ssa_segment

real(8) function electron_quadratic_derivative_x(x,xl,xc,xr,yl,yc,yr)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x,xl,xc,xr,yl,yc,yr

    electron_quadratic_derivative_x= &
        yl*(two*x-xc-xr)/((xl-xc)*(xl-xr))+ &
        yc*(two*x-xl-xr)/((xc-xl)*(xc-xr))+ &
        yr*(two*x-xl-xc)/((xr-xl)*(xr-xc))
end function electron_quadratic_derivative_x

subroutine electron_fill_quadratic_slopes(x_arr,y_arr,dy_arr,n)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: n
real(8), intent(in) :: x_arr(n),y_arr(n)
real(8), intent(out) :: dy_arr(n)
integer :: i

    dy_arr(1)=electron_quadratic_derivative_x(x_arr(1),x_arr(1),x_arr(2),x_arr(3), &
                                              y_arr(1),y_arr(2),y_arr(3))
    do i=2,n-1
        dy_arr(i)=electron_quadratic_derivative_x(x_arr(i),x_arr(i-1),x_arr(i),x_arr(i+1), &
                                                  y_arr(i-1),y_arr(i),y_arr(i+1))
    enddo
    dy_arr(n)=electron_quadratic_derivative_x(x_arr(n),x_arr(n-2),x_arr(n-1),x_arr(n), &
                                              y_arr(n-2),y_arr(n-1),y_arr(n))
end subroutine electron_fill_quadratic_slopes

real(8) function electron_hermite_interp_x(x,x0,x1,y0,y1,dy0,dy1)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x,x0,x1,y0,y1,dy0,dy1
real(8) :: h,s

    h=x1-x0
    s=(x-x0)/h
    electron_hermite_interp_x=(two*s*s*s-three*s*s+one)*y0+ &
                               (s*s*s-two*s*s+s)*h*dy0+ &
                               (-two*s*s*s+three*s*s)*y1+ &
                               (s*s*s-s*s)*h*dy1
end function electron_hermite_interp_x

real(8) function electron_hermite_derivative_x(x,x0,x1,y0,y1,dy0,dy1)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x,x0,x1,y0,y1,dy0,dy1
real(8) :: h,s

    h=x1-x0
    s=(x-x0)/h
    electron_hermite_derivative_x=((6d0*s*s-6d0*s)/h)*y0+ &
                                   (3d0*s*s-4d0*s+one)*dy0+ &
                                   ((-6d0*s*s+6d0*s)/h)*y1+ &
                                   (3d0*s*s-2d0*s)*dy1
end function electron_hermite_derivative_x

real(8) function electron_tau_integrand_x(x,x0,x1,dN10,dN11,ddN10,ddN11,V_cal,DB,factor)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x,x0,x1,dN10,dN11,ddN10,ddN11,V_cal,DB,factor
real(8) :: gam,d_dN1_dx

    gam=dexp(x)
    d_dN1_dx=electron_hermite_derivative_x(x,x0,x1,dN10,dN11,ddN10,ddN11)
    electron_tau_integrand_x=-d_dN1_dx*gam*gam*electron_syn_fx(gam,V_cal,DB,factor)
end function electron_tau_integrand_x

subroutine electron_syn_gauss_cell(x0,x1,dN0,dN1,V_cal,DB,factor,p2,p3)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,dN0,dN1,V_cal,DB,factor
real(8), intent(out) :: p2,p3
real(8) :: xm,dx,w2,w3a,w3b
real(8) :: x2a,x2b,x3a,x3b

    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    w2=one/dsqrt(3d0)
    x2a=xm-dx*w2
    x2b=xm+dx*w2
    p2=dx*(electron_syn_integrand_x(x2a,x0,x1,dN0,dN1,V_cal,DB,factor)+ &
           electron_syn_integrand_x(x2b,x0,x1,dN0,dN1,V_cal,DB,factor))

    w3a=dsqrt(3d0/5d0)
    x3a=xm-dx*w3a
    x3b=xm+dx*w3a
    p3=dx*((5d0/9d0)*electron_syn_integrand_x(x3a,x0,x1,dN0,dN1,V_cal,DB,factor)+ &
           (8d0/9d0)*electron_syn_integrand_x(xm ,x0,x1,dN0,dN1,V_cal,DB,factor)+ &
           (5d0/9d0)*electron_syn_integrand_x(x3b,x0,x1,dN0,dN1,V_cal,DB,factor))
end subroutine electron_syn_gauss_cell

subroutine electron_tau_gauss_cell(x0,x1,dN10,dN11,ddN10,ddN11,V_cal,DB,factor,t2,t3)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,dN10,dN11,ddN10,ddN11,V_cal,DB,factor
real(8), intent(out) :: t2,t3
real(8) :: xm,dx,w2,w3a
real(8) :: x2a,x2b,x3a,x3b

    if (x1 <= x0) then
        t2=zero
        t3=zero
        return
    end if

    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    w2=one/dsqrt(3d0)
    x2a=xm-dx*w2
    x2b=xm+dx*w2
    t2=dx*(electron_tau_integrand_x(x2a,x0,x1,dN10,dN11,ddN10,ddN11,V_cal,DB,factor)+ &
           electron_tau_integrand_x(x2b,x0,x1,dN10,dN11,ddN10,ddN11,V_cal,DB,factor))

    w3a=dsqrt(3d0/5d0)
    x3a=xm-dx*w3a
    x3b=xm+dx*w3a
    t3=dx*((5d0/9d0)*electron_tau_integrand_x(x3a,x0,x1,dN10,dN11,ddN10,ddN11,V_cal,DB,factor)+ &
           (8d0/9d0)*electron_tau_integrand_x(xm ,x0,x1,dN10,dN11,ddN10,ddN11,V_cal,DB,factor)+ &
           (5d0/9d0)*electron_tau_integrand_x(x3b,x0,x1,dN10,dN11,ddN10,ddN11,V_cal,DB,factor))
end subroutine electron_tau_gauss_cell

subroutine electron_syn_cell_adaptive(x0,x1,dN0,dN1,dN10,dN11,ddN10,ddN11, &
                                      V_cal,DB,factor,rel_tol,p_int,tau_int)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,dN0,dN1,dN10,dN11,ddN10,ddN11,V_cal,DB,factor,rel_tol
real(8), intent(out) :: p_int,tau_int
real(8) :: p2,p3,t2,t3,xm,dNm,dN1m,ddN1m,err_p,err_t,ref_p,ref_t
real(8) :: p3_l,p3_r,t3_l,t3_r

    call electron_syn_gauss_cell(x0,x1,dN0,dN1,V_cal,DB,factor,p2,p3)
    call electron_tau_gauss_cell(x0,x1,dN10,dN11,ddN10,ddN11,V_cal,DB,factor,t2,t3)
    ref_p=max(abs(p3),1d-30)
    ref_t=max(abs(t3),1d-30)
    err_p=abs(p3-p2)/ref_p
    err_t=abs(t3-t2)/ref_t

    if (max(err_p,err_t) <= rel_tol) then
        p_int=p3
        tau_int=t3
    else
        xm=0.5d0*(x0+x1)
        dNm=0.5d0*(dN0+dN1)
        dN1m=electron_hermite_interp_x(xm,x0,x1,dN10,dN11,ddN10,ddN11)
        ddN1m=electron_hermite_derivative_x(xm,x0,x1,dN10,dN11,ddN10,ddN11)
        call electron_syn_gauss_cell(x0,xm,dN0,dNm,V_cal,DB,factor,p2,p3_l)
        call electron_syn_gauss_cell(xm,x1,dNm,dN1,V_cal,DB,factor,p2,p3_r)
        call electron_tau_gauss_cell(x0,xm,dN10,dN1m,ddN10,ddN1m,V_cal,DB,factor,t2,t3_l)
        call electron_tau_gauss_cell(xm,x1,dN1m,dN11,ddN1m,ddN11,V_cal,DB,factor,t2,t3_r)
        p_int=p3_l+p3_r
        tau_int=t3_l+t3_r
    end if
end subroutine electron_syn_cell_adaptive

subroutine get_syn_adaptive(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                            P_syn,Seed_syn)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: P_syn(Num_nu),Seed_syn(Num_nu)
real(8), parameter :: rel_tol=5d-4
real(8), allocatable :: dN1(:),ddN1(:),x_gam(:)
real(8) :: factor,Temp_syn,Rariv2,temp_para

    allocate(dN1(Num_gam_e),ddN1(Num_gam_e),x_gam(Num_gam_e))

    factor=(3.62d0/pi)**2
    Temp_syn=dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    Rariv2=R_loc*R_loc
    dN1=dN_gam_e/(gam_e*gam_e)
    x_gam=dlog(gam_e)
    call electron_fill_quadratic_slopes(x_gam,dN1,ddN1,Num_gam_e)

    !$OMP PARALLEL num_threads(n_threads), private(I_nu,I_gam_e,V_cal,dInteg,Tau,P_v, &
    !$OMP& cell_int,tau_cell)
    !$OMP DO SCHEDULE(STATIC)
    do I_nu=1,Num_nu
       V_cal=V_seed(I_nu)
       dInteg=zero
       Tau=zero
       do I_gam_e=1,Num_gam_e-1
          call electron_syn_cell_adaptive(x_gam(I_gam_e),x_gam(I_gam_e+1), &
                                          dN_gam_e(I_gam_e),dN_gam_e(I_gam_e+1), &
                                          dN1(I_gam_e),dN1(I_gam_e+1), &
                                          ddN1(I_gam_e),ddN1(I_gam_e+1), &
                                          V_cal,DB,factor,rel_tol,cell_int,tau_cell)
          dInteg=dInteg+cell_int
          Tau=Tau+tau_cell
       end do
       P_v=Temp_syn*DB*dInteg
       Tau=1.046d4*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal)
       if ((Tau-1d-4) < 1d-5) Tau=1d-4
       P_syn(I_nu)=P_v*(one-dexp(-Tau))/Tau
       Seed_syn(I_nu)=P_syn(I_nu)/(Rariv2*V_cal)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    temp_para=4d0*pi*Para_c*Para_h
    Seed_syn=Seed_syn/temp_para

    deallocate(dN1,ddN1,x_gam)
end subroutine get_syn_adaptive

subroutine get_syn_selected(index_syn_intger,R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                            P_syn,Seed_syn)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_syn_intger,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) :: P_syn(Num_nu),Seed_syn(Num_nu)
real(8) :: h_ref,h_loc

    select case(index_syn_intger)
    case(1)
        if (Num_gam_e > 2) then
            h_ref=dlog(gam_e(2))-dlog(gam_e(1))
            do I_gam_e=3,Num_gam_e
                h_loc=dlog(gam_e(I_gam_e))-dlog(gam_e(I_gam_e-1))
                if (abs(h_loc-h_ref) > 1d-6*max(abs(h_ref),1d-30)) then
                    call get_syn_adaptive(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,P_syn,Seed_syn)
                    return
                end if
            end do
        end if
        call get_syn(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,P_syn,Seed_syn)
    case(2)
        call get_syn_adaptive(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,P_syn,Seed_syn)
    case default
        print*, 'invalid synchrotron integral case, check your chosen model!'
        stop
    end select
end subroutine get_syn_selected

subroutine get_syn_simpson(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                   P_syn,Seed_syn)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) ::P_syn(Num_nu),Seed_syn(Num_nu)

real(8),allocatable,dimension (:) :: dN1,ddN
allocate (dN1(Num_gam_e),ddN(Num_gam_e-1))

    if (Num_gam_e > 2) then
        h = log(gam_e(2))-log(gam_e(1))
        do I_gam_e=3,Num_gam_e
            if (abs((log(gam_e(I_gam_e))-log(gam_e(I_gam_e-1)))-h) > 1d-6*max(abs(h),1d-30)) then
                deallocate(dN1,ddN)
                call get_syn(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed,P_syn,Seed_syn)
                return
            end if
        end do
    end if

    factor=(3.62d0/pi)**2
    Temp_syn=dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    Rariv2=R_loc*R_loc
    dN1=dN_gam_e/(gam_e*gam_e)
    ddN=dN1(1:Num_gam_e-1)-dN1(2:Num_gam_e)
    
    h = log(gam_e(2))-log(gam_e(1))
    
    !$ call omp_set_dynamic(.true.)
    !$OMP PARALLEL num_threads(n_threads), private(I_nu,I_gam_e,V_cal,dInteg,Tau,Vc,x,Fx,P_v,simpson_sum, &
    !$OMP& val,gam_e_mean2,ratio_v)
    !$OMP DO SCHEDULE(STATIC)
    do I_nu=1,Num_nu
       V_cal=V_seed(I_nu)
       dInteg=zero
       Tau=zero

       ! \int f(x)dx \prox (h/3)[f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + ... + f(xn)]
       simpson_sum = zero
       do I_gam_e=1,Num_gam_e
           Vc = (4.2d6)*gam_e(I_gam_e)**2*DB
           x = V_cal/Vc
           ratio_v = Vc/V_cal
           Fx = 1.81d0*exp(-x)/sqrt(ratio_v**(2d0/3d0)+factor)
           val = dN_gam_e(I_gam_e) * Fx * gam_e(I_gam_e)
           if (I_gam_e == 1 .or. I_gam_e == Num_gam_e) then
               simpson_sum = simpson_sum + val
           else if (mod(I_gam_e,2) == 0) then
               simpson_sum = simpson_sum + 4d0 * val
           else
               simpson_sum = simpson_sum + two * val
           endif
       end do
       dInteg = (h/3.0d0) * simpson_sum
       ! ====================  [SSA]  ======================
       do I_gam_e=1,Num_gam_e-1
          gam_e_mean2=(gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0
          Vc=(4.2d6)*gam_e_mean2*DB
          x=V_cal/Vc
          ratio_v=Vc/V_cal
          Fx=1.81d0*exp(-x)/sqrt(ratio_v**(2d0/3d0)+factor)
          Tau=Tau+gam_e_mean2*ddN(I_gam_e)*Fx
       end do
       ! ===================================================
       P_v=Temp_syn*DB*dInteg
       Tau=1.046d4*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal)
       if ((Tau-1d-4) < 1d-5) Tau=1d-4
       P_syn(I_nu)=P_v*(one-exp(-Tau))/Tau
       Seed_syn(I_nu)=P_syn(I_nu)/(Rariv2*V_cal)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    temp_para=4d0*pi*Para_c*Para_h
    Seed_syn=Seed_syn/temp_para

deallocate (dN1,ddN)
end subroutine get_syn_simpson

subroutine get_forward_cooling(index_Y,Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,R_loc,R_Gamma_loc, &
                               beta_Gam,dNe,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,Seed_syn,dEl)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: index_Y,Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,R_loc,R_Gamma_loc,beta_Gam,dNe
real(8), intent(inout) :: Gam_e_max
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dEl(Num_gam_e)
real(8) :: Compton(Num_gam_e),dot_gam_e(Num_gam_e),dot_gam_e_SSA(Num_gam_e)

    call get_SSA_numerical(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn,dot_gam_e_SSA)
    cooling_scale=one/(beta_Gam*R_Gamma_loc)
    ssa_scale=cooling_scale/para_c
    f_r=1.35d-19*DB**2*cooling_scale/pi

    select case(index_Y)
    case(0)
        dEl=(f_r-dot_gam_e_SSA*ssa_scale)*gam_e
    case(1)
        call get_IC_numerical(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn,dot_gam_e)
        dEl=(f_r+(dot_gam_e-dot_gam_e_SSA)*ssa_scale)*gam_e

    case(2)
        call get_Y_Nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,Compton)
        Q=4d0*pi*R_loc*R_loc*para_c
        Compton=one+Compton/Q/(4d0*R_Gamma_loc*R_Gamma_loc*dNe*Para_m_p_E)
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA*ssa_scale)*gam_e

    case(3)
        call get_Y_Fan(Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,Compton)
        Compton=one+Compton
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA*ssa_scale)*gam_e

    case default
        print*, 'invalid Compton case, check your chosen model!'
        stop
    end select
end subroutine get_forward_cooling

!****************************************************************************************
!**************** get SSA pile-up effect with Ghisellini & Svensson 1991 ****************
!****************************************************************************************
subroutine get_SSA_numerical(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn, dot_gam_e)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dot_gam_e(Num_gam_e)
integer :: low_idx
real(8) :: cell_low,cell_high

    B_cr=4.4d13
    Temp1=2.5042d-22*B_cr/DB
    Temp2=7.787d-22*B_cr/DB
    Cyclotron_nu=para_e*DB/(two*pi*para_m_e*para_c)
    
    dot_gam_e=zero

    !$OMP PARALLEL num_threads(n_threads), private(I_gam_e, gam, gam2, gam3, V_lowlim, V_uplim, &
    !$OMP& I_nu, ssa_sum, low_idx, sigma_prefactor, cell_low, cell_high)
    !$OMP DO SCHEDULE(STATIC)
    do I_gam_e=1,Num_gam_e
       gam=gam_e(I_gam_e)
       gam2=gam*gam
       gam3=gam2*gam
       V_lowlim=Cyclotron_nu/gam
       V_uplim=1.5d0*gam2*Cyclotron_nu
       ssa_sum=zero
       call first_greater_monotonic(V_seed,Num_nu,V_lowlim,low_idx)
       if (low_idx <= Num_nu) then
          do I_nu=max(1,low_idx-1),Num_nu-1
             if (V_seed(I_nu+1) <= V_lowlim) cycle
             cell_low=max(V_seed(I_nu),V_lowlim)
             cell_high=min(V_seed(I_nu+1),V_uplim)
             if (cell_high > cell_low) then
                sigma_prefactor=Temp1*(3d0*V_lowlim)**(5d0/3d0)
                ssa_sum=ssa_sum+electron_ssa_segment(cell_low,cell_high,Seed_syn(I_nu),Seed_syn(I_nu+1), &
                                                     sigma_prefactor,1,Cyclotron_nu,V_uplim)
             end if

             cell_low=max(V_seed(I_nu),V_uplim)
             cell_high=V_seed(I_nu+1)
             if (cell_high > cell_low) then
                sigma_prefactor=Temp2/gam3
                ssa_sum=ssa_sum+electron_ssa_segment(cell_low,cell_high,Seed_syn(I_nu),Seed_syn(I_nu+1), &
                                                     sigma_prefactor,2,Cyclotron_nu,V_uplim)
             end if
          end do
       end if
       dot_gam_e(I_gam_e)=ssa_sum
    end do
    !$OMP END DO
    !$OMP END PARALLEL

end subroutine get_SSA_numerical


!****************************************************************************************
!******************* get Compton cooling effiency fully numerical ***********************
!****************************************************************************************
subroutine get_IC_numerical(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn, dot_gam_e)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dot_gam_e(Num_gam_e)

real(8),allocatable,dimension (:) :: d_nu,photon_number,gam_e_mean,E_seed,x_seed,V_seed_mid
real(8) :: kn_factor

    allocate (d_nu(Num_nu-1),photon_number(Num_nu-1),gam_e_mean(Num_gam_e-1), &
              E_seed(Num_nu-1),x_seed(Num_nu),V_seed_mid(Num_nu-1))

    para_hEme = Para_h/para_m_energy

    x_seed=dlog(V_seed)
    V_seed_mid=dexp(0.5d0*(x_seed(1:Num_nu-1)+x_seed(2:Num_nu)))
    d_nu=V_seed_mid*(x_seed(2:Num_nu)-x_seed(1:Num_nu-1))
    gam_e_mean=(gam_e(1:Num_gam_e-1)+gam_e(2:Num_gam_e))/two
    E_seed=V_seed_mid*para_hEme
    
    dot_gam_e=zero

    do I_nu=1,Num_nu-1
       photon_number(I_nu)=electron_powerlaw_interp(V_seed(I_nu),V_seed(I_nu+1), &
                                                    Seed_syn(I_nu),Seed_syn(I_nu+1),V_seed_mid(I_nu))
    end do
    
    !$ call omp_set_dynamic(.true.)
    !$OMP PARALLEL num_threads(n_threads), &
    !$OMP& private(i_gam_e, dInteg1, game, game_pow, var, I_nu, dInteg2, &
    !$OMP& V_t, E_t2eV, Nu_s, fssc, Vloc, E_Vloc2eV, uplim, temp, q, kn_factor)
    !$OMP DO SCHEDULE(STATIC)
    do i_gam_e=1,Num_gam_e-1
       dInteg1=zero
       game=gam_e_mean(i_gam_e)
       game_pow=game*game
       var=0.25d0/game_pow
       do I_nu=1,Num_nu-1      !frequency circulation for seed photons
          dInteg2=zero
          V_t=V_seed_mid(I_nu)
          E_t2eV=E_seed(I_nu)
          kn_factor=4d0*game*E_t2eV
          uplim=(4d0*game_pow*E_t2eV)/(one+kn_factor)
          do Nu_s=1,Num_nu-1     !frequency circulation for SSC photons
             fssc=zero
             Vloc=V_seed_mid(Nu_s)
             E_Vloc2eV=E_seed(Nu_s)
             if (Vloc > var*V_t .and. Vloc <= V_t) then
                fssc=Vloc/V_t-var
             else
                if (E_Vloc2eV > uplim) exit
                temp=game-E_Vloc2eV
                q=E_Vloc2eV/(kn_factor*temp)
                fssc=two*q*(log(q)-q)+one+q+ &
                0.5d0*(one-q)*(4d0*game*E_t2eV*q)**2/(1+4d0*game*q*E_t2eV)
             end if
             dInteg2=dInteg2+Vloc*fssc*d_nu(Nu_s)
          end do
          dot_gam_e(i_gam_e)=dot_gam_e(i_gam_e)+photon_number(I_nu)/V_t*d_nu(I_nu)*dInteg2
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    dot_gam_e=dot_gam_e/gam_e/gam_e*para_h*Para_h*Para_SigmaT/para_m_energy
    dot_gam_e(Num_gam_e)=0.99*dot_gam_e(Num_gam_e-1)

deallocate (d_nu,photon_number,gam_e_mean,E_seed,x_seed,V_seed_mid)
end subroutine get_IC_numerical


!****************************************************************************************
!********************* get Compton Y with Nakar & Piran 2009 ****************************
!****************************************************************************************
subroutine get_Y_Nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn, Compton)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu)
real(8), intent(out) :: Compton(Num_gam_e)

real(8),allocatable,dimension (:) :: hat_nu,prefix_syn

    allocate (hat_nu(Num_gam_e),prefix_syn(Num_nu))

    Compton=zero
    var_Compensation=zero
    hat_nu=Para_m_energy/Para_h/gam_e
    prefix_syn(1)=zero
    do I_nu=2,Num_nu
       prefix_syn(I_nu)=prefix_syn(I_nu-1)+ &
                        electron_integrate_powerlaw_segment(V_seed(I_nu-1),V_seed(I_nu), &
                                                            P_syn(I_nu-1),P_syn(I_nu))
    end do

    !$OMP PARALLEL num_threads(n_threads), private(I_Compton, I_nu, temp, var_Compensation, V_loc)
    !$OMP DO SCHEDULE(STATIC)
    do I_Compton=1,Num_gam_e
       if (hat_nu(I_Compton) <= V_seed(1)) cycle
       call first_greater_monotonic_window(V_seed,Num_nu,2,hat_nu(I_Compton),I_nu)
       if (I_nu <= Num_nu) then
          V_loc=min(hat_nu(I_Compton),V_seed(I_nu))
          Compton(I_Compton)=prefix_syn(I_nu-1)+ &
                             electron_integrate_powerlaw_segment(V_seed(I_nu-1),V_loc, &
                                 P_syn(I_nu-1), &
                                 electron_powerlaw_interp(V_seed(I_nu-1),V_seed(I_nu), &
                                                          P_syn(I_nu-1),P_syn(I_nu),V_loc))
       else
          Compton(I_Compton)=prefix_syn(Num_nu)
       end if
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
deallocate (hat_nu,prefix_syn)
end subroutine get_Y_Nakar

!****************************************************************************************
!************************ get Compton Y with Fan et al. 2008 ****************************
!****************************************************************************************
subroutine get_Y_Fan(Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e, Compton)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,gam_e(Num_gam_e)
real(8), intent(out) :: Compton(Num_gam_e)


    eta=(Gam_e_m/Gam_e_c)**(p-two)
    if (eta-one > 0.001) eta=one

        do i_gam_e=1,Num_gam_e
            if (Num_gam_e > i_gam_e) then
!*****************************[The general inverse Compton effect]******************************
                hat_gam=5.4246D6/sqrt(DB*gam_e(i_gam_e+1))
                if (Gam_e_m > Gam_e_c) then
                    if (hat_gam < Gam_e_c) then
                        eta_NK=zero
                    else
                        if (hat_gam < Gam_e_m) then
                            if (p>2) then
                                Step1=(p-1)/(p-2)*Gam_e_m-Gam_e_c
                                eta_NK=(hat_gam-Gam_e_c)/Step1
                            else
                                Step1=Gam_e_m**(p-1)*Gam_e_max**(2-p)-(p-1)*Gam_e_m-(2-p)*Gam_e_c
                                eta_NK=(2-p)*(hat_gam-Gam_e_c)/Step1
                            end if
                        else
                            if (p>2) then
                                Step2=Gam_e_m**(p-1)*hat_gam**(2-p)
                                Step3=(p-1)*Gam_e_m-(p-2)*Gam_e_c
                                eta_NK=1-Step2/Step3
                            else
                                Step2=Gam_e_m**(p-1)*Gam_e_max**(2-p)-(p-1)*Gam_e_m-(2-p)*Gam_e_c
                                Step3=Gam_e_m**(p-1)*(Gam_e_max**(2-p)-hat_gam**(2-p))
                                eta_NK=1-Step2/Step3
                            end if
                        end if
                    end if
                else
                    if (hat_gam < Gam_e_m) then
                        eta_NK=zero
                    else
                        if (hat_gam < Gam_e_c) then
                            if (p>2) then
                                Step4=Gam_e_c**(3-p)/(p-2.0)-Gam_e_m**(3-p)
                                eta_NK=(hat_gam**(3-p)-Gam_e_m**(3-p))/Step4
                            else
                                Step4=(3-p)*Gam_e_c*Gam_e_max**(2-p)-Gam_e_c**(3-p)-(2-p)*Gam_e_m**(3-p)
                                eta_NK=(2-p)*(hat_gam**(3-p)-Gam_e_m**(3-p))/Step4
                            end if
                        else
                            if (p>2) then
                                Step5=(3-p)*Gam_e_c*hat_gam**(2-p)
                                Step6=Gam_e_c**(3.0-p)-(p-2)*Gam_e_m**(3.0-p)
                                eta_NK=1-Step5/Step6
                            else
                                Step5=(3-p)*Gam_e_c*(Gam_e_max**(2-p)-hat_gam**(2-p))
                                Step6=(3-p)*Gam_e_c*Gam_e_max**(2-p)-Gam_e_c**(3-p)-(2-p)*Gam_e_m**(3-p)
                                eta_NK=1-Step5/Step6
                            end if
                        end if
                    end if
                end if
                Compton(i_gam_e)=0.5d0*(-1.0+sqrt(1.0+4.0*eta*eta_NK*Epsilon_e/Epsilon_b))
            else
                Compton(i_gam_e)=0.99*Compton(i_gam_e-1)
            end if
        end do


end subroutine get_Y_Fan

subroutine get_nu_a(R_loc,DB,Num_gam_e,gam_e,dN_gam_e, &
                       V_a)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e)
real(8), intent(out) :: V_a
real(8), parameter :: rel_tol=5d-4

real(8),allocatable,dimension (:) :: dN1,ddN1,x_gam
allocate (dN1(Num_gam_e),ddN1(Num_gam_e),x_gam(Num_gam_e))

    factor=(3.62d0/pi)**2
    Rariv2=R_loc*R_loc
    dN1=dN_gam_e/(gam_e*gam_e)
    x_gam=dlog(gam_e)
    call electron_fill_quadratic_slopes(x_gam,dN1,ddN1,Num_gam_e)
    
    V_a_floor=ten**4d0
    V_a_min=ten**(-20d0)
    V_a_cap=one
    do I_gam_e=1,Num_gam_e-1
       Vc=4.2d6*((gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0)*DB
       if (Vc > V_a_cap) V_a_cap=Vc
    end do
    V_a_cap=max(ten**14d0, min(ten**30d0, ten*V_a_cap))

    V_low=V_a_floor
    call evaluate_tau(V_low,Tau_low)
    if (Tau_low <= one) then
        do I_nu=1,26
           if (V_low <= V_a_min) exit
           V_high=V_low
           Tau_high=Tau_low
           V_low=max(V_a_min,V_low/ten)
           call evaluate_tau(V_low,Tau_low)
           if (Tau_low > one) exit
        end do

        if (Tau_low > one) then
           call refine_nu_a_bracket(V_low,Tau_low,V_high,Tau_high,V_a)
        else
           V_a=V_low
        end if
    else
        V_high=V_low
        Tau_high=Tau_low
       do I_nu=1,26
          V_low=V_high
          Tau_low=Tau_high
          if (V_high >= V_a_cap) exit
          V_high=min(V_a_cap, ten*V_high)
          call evaluate_tau(V_high,Tau_high)
          if (Tau_high <= one) exit
       end do

       if (Tau_high > one) then
          V_a=V_high
          print*, 'nu_a_comoving larger than adaptive upper bound!', V_a_cap
       else
          if (Tau_low <= one .or. Tau_low == Tau_high) then
             V_a=V_high
          else
             call refine_nu_a_bracket(V_low,Tau_low,V_high,Tau_high,V_a)
          end if
       end if
    end if

    deallocate (dN1,ddN1,x_gam)

return

contains

subroutine evaluate_tau(V_cal,Tau)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: V_cal
real(8), intent(out) :: Tau
real(8) :: tau_cell,discard

    Tau=zero
    do I_gam_e=1,Num_gam_e-1
       call electron_syn_cell_adaptive(x_gam(I_gam_e),x_gam(I_gam_e+1), &
                                       dN_gam_e(I_gam_e),dN_gam_e(I_gam_e+1), &
                                       dN1(I_gam_e),dN1(I_gam_e+1), &
                                       ddN1(I_gam_e),ddN1(I_gam_e+1), &
                                       V_cal,DB,factor,rel_tol,discard,tau_cell)
       Tau=Tau+tau_cell
    end do
    Tau=1.046d4*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal)
end subroutine evaluate_tau

subroutine refine_nu_a_bracket(V_low,Tau_low,V_high,Tau_high,V_root)
implicit REAL(8)(A-H,O-Z)
real(8), intent(inout) :: V_low,Tau_low,V_high,Tau_high
real(8), intent(out) :: V_root
real(8) :: log_v_low,log_v_high,log_v_mid,log_tau_low,log_tau_high
real(8) :: V_mid,Tau_mid
integer :: I_iter

    log_v_low=dlog(V_low)
    log_v_high=dlog(V_high)

    do I_iter=1,10
        log_tau_low=dlog(Tau_low)
        log_tau_high=dlog(Tau_high)
        if (log_tau_low == log_tau_high) then
            log_v_mid=0.5d0*(log_v_low+log_v_high)
        else
            log_v_mid=log_v_low-log_tau_low*(log_v_high-log_v_low)/(log_tau_high-log_tau_low)
        end if
        log_v_mid=max(log_v_low,min(log_v_high,log_v_mid))
        V_mid=dexp(log_v_mid)
        call evaluate_tau(V_mid,Tau_mid)

        if (Tau_mid > one) then
            V_low=V_mid
            Tau_low=Tau_mid
            log_v_low=log_v_mid
        else
            V_high=V_mid
            Tau_high=Tau_mid
            log_v_high=log_v_mid
        end if

        if (abs(log_v_high-log_v_low) <= 5d-3) exit
    end do

    log_tau_low=dlog(Tau_low)
    log_tau_high=dlog(Tau_high)
    if (log_tau_low == log_tau_high) then
        V_root=dexp(0.5d0*(log_v_low+log_v_high))
    else
        V_root=dexp(log_v_low-log_tau_low*(log_v_high-log_v_low)/(log_tau_high-log_tau_low))
    end if
end subroutine refine_nu_a_bracket
end subroutine get_nu_a

subroutine get_nu_a_nonuniform(R_loc,DB,Num_gam_e,x_edge_log10,dN_x,V_a)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e
real(8), intent(in) :: R_loc,DB,x_edge_log10(Num_gam_e+1),dN_x(Num_gam_e)
real(8), intent(out) :: V_a
real(8), parameter :: rel_tol=5d-4

real(8), allocatable :: x_edge(:),q_y(:),q_left(:),q_right(:)
allocate(x_edge(Num_gam_e+1),q_y(Num_gam_e),q_left(Num_gam_e),q_right(Num_gam_e))

    factor=(3.62d0/pi)**2
    Rariv2=R_loc*R_loc
    x_edge=dlog(ten)*x_edge_log10
    q_y=dN_x/dlog(ten)
    call electron_ppm_interfaces_nonuniform(Num_gam_e,x_edge,q_y,q_left,q_right)

    V_a_floor=ten**4d0
    V_a_min=ten**(-20d0)
    V_a_cap=one
    do I_gam_e=1,Num_gam_e
       gam_mid=dexp(0.5d0*(x_edge(I_gam_e)+x_edge(I_gam_e+1)))
       Vc=4.2d6*gam_mid*gam_mid*DB
       if (Vc > V_a_cap) V_a_cap=Vc
    end do
    V_a_cap=max(ten**14d0, min(ten**30d0, ten*V_a_cap))

    V_low=V_a_floor
    call evaluate_tau(V_low,Tau_low)
    if (Tau_low <= one) then
        do I_nu=1,26
           if (V_low <= V_a_min) exit
           V_high=V_low
           Tau_high=Tau_low
           V_low=max(V_a_min,V_low/ten)
           call evaluate_tau(V_low,Tau_low)
           if (Tau_low > one) exit
        end do

        if (Tau_low > one) then
           call refine_nu_a_bracket_nonuniform(V_low,Tau_low,V_high,Tau_high,V_a)
        else
           V_a=V_low
        end if
    else
        V_high=V_low
        Tau_high=Tau_low
       do I_nu=1,26
          V_low=V_high
          Tau_low=Tau_high
          if (V_high >= V_a_cap) exit
          V_high=min(V_a_cap, ten*V_high)
          call evaluate_tau(V_high,Tau_high)
          if (Tau_high <= one) exit
       end do

       if (Tau_high > one) then
          V_a=V_high
          print*, 'nu_a_comoving larger than adaptive upper bound!', V_a_cap
       else
          if (Tau_low <= one .or. Tau_low == Tau_high) then
             V_a=V_high
          else
             call refine_nu_a_bracket_nonuniform(V_low,Tau_low,V_high,Tau_high,V_a)
          end if
       end if
    end if

    deallocate(x_edge,q_y,q_left,q_right)

return

contains

subroutine evaluate_tau(V_cal,Tau)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: V_cal
real(8), intent(out) :: Tau
real(8) :: tau_cell

    Tau=zero
    do I_gam_e=1,Num_gam_e
       call electron_tau_ppm_cell_adaptive(x_edge(I_gam_e),x_edge(I_gam_e+1), &
                                           q_y(I_gam_e),q_left(I_gam_e),q_right(I_gam_e), &
                                           V_cal,DB,factor,rel_tol,tau_cell)
       Tau=Tau+tau_cell
    end do
    Tau=1.046d4*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal)
end subroutine evaluate_tau

subroutine electron_ppm_value_derivative_nonuniform(x,cell_lo,cell_hi,qc,q_l,q_r,q_val,dqdx)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x,cell_lo,cell_hi,qc,q_l,q_r
real(8), intent(out) :: q_val,dqdx
real(8) :: dx_cell,xi,q6,bcoef

    dx_cell=cell_hi-cell_lo
    if (dx_cell <= zero) then
        q_val=qc
        dqdx=zero
        return
    end if

    xi=(x-cell_lo)/dx_cell
    xi=max(zero,min(one,xi))
    q6=6d0*qc-3d0*(q_l+q_r)
    bcoef=q_r-q_l+q6
    q_val=q_l+bcoef*xi-q6*xi*xi
    dqdx=(bcoef-two*q6*xi)/dx_cell
end subroutine electron_ppm_value_derivative_nonuniform

real(8) function electron_tau_ppm_integrand_x(x,cell_lo,cell_hi,qc,q_l,q_r,V_cal,DB,factor)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x,cell_lo,cell_hi,qc,q_l,q_r,V_cal,DB,factor
real(8) :: gam,q_val,dqdx,d_dN1_dx

    call electron_ppm_value_derivative_nonuniform(x,cell_lo,cell_hi,qc,q_l,q_r,q_val,dqdx)
    gam=dexp(x)
    d_dN1_dx=(dqdx-3d0*q_val)/(gam*gam*gam)
    electron_tau_ppm_integrand_x=-d_dN1_dx*gam*gam*electron_syn_fx(gam,V_cal,DB,factor)
end function electron_tau_ppm_integrand_x

subroutine electron_tau_ppm_gauss_cell(x0,x1,cell_lo,cell_hi,qc,q_l,q_r,V_cal,DB,factor,t2,t3)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,cell_lo,cell_hi,qc,q_l,q_r,V_cal,DB,factor
real(8), intent(out) :: t2,t3
real(8) :: xm,dx,w2,w3a,x2a,x2b,x3a,x3b

    if (x1 <= x0) then
        t2=zero
        t3=zero
        return
    end if

    xm=0.5d0*(x0+x1)
    dx=0.5d0*(x1-x0)
    w2=one/dsqrt(3d0)
    x2a=xm-dx*w2
    x2b=xm+dx*w2
    t2=dx*(electron_tau_ppm_integrand_x(x2a,cell_lo,cell_hi,qc,q_l,q_r,V_cal,DB,factor)+ &
           electron_tau_ppm_integrand_x(x2b,cell_lo,cell_hi,qc,q_l,q_r,V_cal,DB,factor))

    w3a=dsqrt(3d0/5d0)
    x3a=xm-dx*w3a
    x3b=xm+dx*w3a
    t3=dx*((5d0/9d0)*electron_tau_ppm_integrand_x(x3a,cell_lo,cell_hi,qc,q_l,q_r,V_cal,DB,factor)+ &
           (8d0/9d0)*electron_tau_ppm_integrand_x(xm ,cell_lo,cell_hi,qc,q_l,q_r,V_cal,DB,factor)+ &
           (5d0/9d0)*electron_tau_ppm_integrand_x(x3b,cell_lo,cell_hi,qc,q_l,q_r,V_cal,DB,factor))
end subroutine electron_tau_ppm_gauss_cell

subroutine electron_tau_ppm_cell_adaptive(x0,x1,qc,q_l,q_r,V_cal,DB,factor,rel_tol,tau_int)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: x0,x1,qc,q_l,q_r,V_cal,DB,factor,rel_tol
real(8), intent(out) :: tau_int
real(8) :: t2,t3,xm,ref_t,err_t,t3_l,t3_r

    call electron_tau_ppm_gauss_cell(x0,x1,x0,x1,qc,q_l,q_r,V_cal,DB,factor,t2,t3)
    ref_t=max(abs(t3),1d-30)
    err_t=abs(t3-t2)/ref_t
    if (err_t <= rel_tol) then
        tau_int=t3
    else
        xm=0.5d0*(x0+x1)
        call electron_tau_ppm_gauss_cell(x0,xm,x0,x1,qc,q_l,q_r,V_cal,DB,factor,t2,t3_l)
        call electron_tau_ppm_gauss_cell(xm,x1,x0,x1,qc,q_l,q_r,V_cal,DB,factor,t2,t3_r)
        tau_int=t3_l+t3_r
    end if
end subroutine electron_tau_ppm_cell_adaptive

subroutine refine_nu_a_bracket_nonuniform(V_low,Tau_low,V_high,Tau_high,V_root)
implicit REAL(8)(A-H,O-Z)
real(8), intent(inout) :: V_low,Tau_low,V_high,Tau_high
real(8), intent(out) :: V_root
real(8) :: log_v_low,log_v_high,log_v_mid,log_tau_low,log_tau_high
real(8) :: V_mid,Tau_mid
integer :: I_iter

    log_v_low=dlog(V_low)
    log_v_high=dlog(V_high)
    do I_iter=1,10
        log_tau_low=dlog(Tau_low)
        log_tau_high=dlog(Tau_high)
        if (log_tau_low == log_tau_high) then
            log_v_mid=0.5d0*(log_v_low+log_v_high)
        else
            log_v_mid=log_v_low-log_tau_low*(log_v_high-log_v_low)/(log_tau_high-log_tau_low)
        end if
        log_v_mid=max(log_v_low,min(log_v_high,log_v_mid))
        V_mid=dexp(log_v_mid)
        call evaluate_tau(V_mid,Tau_mid)
        if (Tau_mid > one) then
            V_low=V_mid
            Tau_low=Tau_mid
            log_v_low=log_v_mid
        else
            V_high=V_mid
            Tau_high=Tau_mid
            log_v_high=log_v_mid
        end if
        if (abs(log_v_high-log_v_low) <= 5d-3) exit
    end do

    log_tau_low=dlog(Tau_low)
    log_tau_high=dlog(Tau_high)
    if (log_tau_low == log_tau_high) then
        V_root=dexp(0.5d0*(log_v_low+log_v_high))
    else
        V_root=dexp(log_v_low-log_tau_low*(log_v_high-log_v_low)/(log_tau_high-log_tau_low))
    end if
    V_root=min(max(V_root,V_low),V_high)
end subroutine refine_nu_a_bracket_nonuniform
end subroutine get_nu_a_nonuniform

end module get_Y
