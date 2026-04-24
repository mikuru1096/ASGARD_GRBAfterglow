module electron_reverse_kernel
    use constants
    use dynamics_common, only: dynamics_external_density_base, dynamics_reverse_gamma_extrema
    use electron_transport_common, only: electron_prepare_implicit_coeffs_common, electron_backward_sweep_common
    use electron_radiation_kernel, only: get_syn_selected
    use electron_cooling_kernel, only: get_IC_numerical
    use electron_y_kernel, only: get_Y_Nakar, get_Y_Fan
    implicit none
contains

subroutine electron_reverse_evolve(Delta_0,e_r,b_r,p_r,f_e_r,eta_0,Epsilon_e,Epsilon_b,z,A_star,dNe_ISM,para_m_ej, &
                                   T_cross,R_cross,R_Tobs,R_Gamma,R,B3,V_seed,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger, &
                                   n_threads,gam_e,dN_gam_e)
    implicit none
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer :: I_tobs,I_gam_e,L1,L,i
    real(8), intent(in) :: Delta_0,e_r,b_r,p_r,f_e_r,eta_0,Epsilon_e,Epsilon_b,z,A_star,dNe_ISM,para_m_ej,T_cross,R_cross
    real(8), intent(in) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),B3(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R)
    real(8), parameter :: reverse_gamma_c_coeff=7.7d8, reverse_synch_b_coeff=0.39d0, reverse_adv_coeff=1.35d-19
    real(8) :: factor2,dB,gamma34,Gam_e_max,Gam_e_m,Gam_e_c,dNe,DB_min,Gam_e_max_max,d_x,R_loc,R_Gamma_loc,Delta
    real(8) :: R_n4,beta4,beta2,eta,f_r,dDR,dDD,CFL,Q0,Q,Q1,Qshell,cooling_scale
    real(8), allocatable :: dEl(:),principal(:),x(:),dF1(:),up(:),temp1(:),temp2(:),temp3(:),dN_x(:),para_minus_gam_e_p(:)
    real(8), allocatable :: dB3_serial(:),P_syn(:),Seed_syn(:),cooling_aux(:),Compton(:)

    allocate(dEl(Num_gam_e),principal(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e),up(Num_gam_e-1), &
             temp1(Num_gam_e-1),temp2(Num_gam_e),temp3(Num_gam_e-1),dN_x(Num_gam_e), &
             para_minus_gam_e_p(Num_gam_e),dB3_serial(Num_R),P_syn(Num_nu),Seed_syn(Num_nu), &
             cooling_aux(Num_gam_e),Compton(Num_gam_e))
    dB3_serial=B3

    factor2=(p_r-two)/(p_r-one)*e_r*Para_m_p_div_m_e
    if (p_r < 2.05d0) factor2=0.05d0/1.05d0*e_r*Para_m_p_div_m_e
    dB3_serial(1)=dB3_serial(min(2,Num_R))
    dB=dB3_serial(1); gamma34=1.001d0
    call dynamics_reverse_gamma_extrema(dB,gamma34,factor2,f_e_r,Gam_e_max,Gam_e_m)
    Gam_e_c=reverse_gamma_c_coeff/(one+dsqrt(e_r/b_r))/R_Gamma(1)/dB**2/(R_Tobs(1)/two)

    call dynamics_external_density_base(A_star,dNe_ISM,R(1),dNe)
    DB_min=reverse_synch_b_coeff*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)

    do I_gam_e=1,Num_gam_e
        gam_e(I_gam_e)=3d0*ten**(dlog10(Gam_e_max_max)*(I_gam_e-1)/(Num_gam_e-1))
        dN_gam_e(I_gam_e,1)=zero
        if (Gam_e_m > Gam_e_c) then
            Q1=1d10*Gam_e_c
            if ((Gam_e_c-gam_e(I_gam_e)) <= one) then
                if (Gam_e_m > gam_e(I_gam_e)) then
                    dN_gam_e(I_gam_e,1)=Q1*gam_e(I_gam_e)**(-2)
                else if (Gam_e_max > gam_e(I_gam_e)) then
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e_m**(p_r-one)*gam_e(I_gam_e)**(-(p_r+one))
                end if
            end if
        else
            Q1=1d10*Gam_e_m**(p_r-one)
            if (Gam_e_m <= gam_e(I_gam_e)) then
                if (Gam_e_c > gam_e(I_gam_e)) then
                    dN_gam_e(I_gam_e,1)=Q1*gam_e(I_gam_e)**(-p_r)
                else if (Gam_e_max > gam_e(I_gam_e)) then
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e_c*gam_e(I_gam_e)**(-(p_r+one))
                end if
            end if
        end if
    end do

    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
    d_x=dlog10(gam_e(2)/gam_e(1))
    para_minus_gam_e_p=one/(gam_e-one)**p_r

    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1); R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        Delta=max(Delta_0,R_loc/Eta_0**2)
        R_n4=para_m_ej/(4d0*pi*Para_m_p*R_loc*R_loc*Eta_0*Delta)
        beta4=dsqrt(one-one/eta_0**2); beta2=dsqrt(one-one/R_Gamma_loc**2)
        gamma34=(one-beta2*beta4)*eta_0*R_Gamma_loc
        dB=(dB3_serial(I_tobs)+dB3_serial(I_tobs-1))/two
        call dynamics_reverse_gamma_extrema(dB,gamma34,factor2,f_e_r,Gam_e_max,Gam_e_m)
        Gam_e_c=reverse_gamma_c_coeff*(one+z)/R_Gamma_loc/dB**2/R_Tobs(I_tobs)
        eta=(Gam_e_m/Gam_e_c)**(p_r-two); if (eta > one+0.001d0) eta=one
        f_r=reverse_adv_coeff/beta2/R_Gamma_loc*dB**2/pi
        dDR=0.7d0/(f_r*Gam_e_max+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1); L1=max(100,min(1000,int(dDD/dDR))); dDR=dDD/L1; CFL=dDR/d_x
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)

        call get_syn_selected(index_syn_intger,R(I_tobs-1),dB,Num_gam_e,Num_nu,n_threads, &
                              gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_syn,Seed_syn)
        cooling_aux=zero
        Compton=zero
        select case(index_Y)
        case(0)
            dEl=f_r*gam_e
        case(1)
            cooling_scale=one/(beta2*R_Gamma_loc*Para_c)
            call get_IC_numerical(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn,cooling_aux)
            dEl=(f_r+cooling_aux*cooling_scale)*gam_e
        case(2)
            Qshell=4d0*pi*R(I_tobs-1)*R(I_tobs-1)*Para_c
            call get_Y_Nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn,cooling_aux)
            Compton=one+cooling_aux/Qshell/(4d0*R_Gamma_loc*R_Gamma_loc*R_n4*Para_m_p_E)
            Gam_e_max=Gam_e_max/dsqrt(Compton(Num_gam_e))
            dEl=f_r*Compton*gam_e
        case(3)
            call get_Y_Fan(e_r,b_r,p_r,dB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e,Compton)
            Compton=one+Compton
            Gam_e_max=Gam_e_max/dsqrt(Compton(Num_gam_e))
            dEl=f_r*Compton*gam_e
        case default
            print*, 'invalid Compton case, check your chosen model!'
            stop
        end select
        dEl(Num_gam_e)=dEl(Num_gam_e)+0.1d0

        Q0=4d0*pi*R_n4*(p_r-one)*(Gam_e_m-one)**(p_r-one)*f_e_r
        do L=1,L1
            R_loc=R_loc+dDR
            if (R_cross >= R_loc) then
                Q=Q0*R_loc*R_loc
            else
                Q=zero
            end if
            dF1=zero
            where(gam_e < Gam_e_max .and. gam_e > Gam_e_m) dF1=Q*para_minus_gam_e_p*gam_e*dlog(ten)
            temp3=((dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/two+one/R_loc)/dlog(ten)
            up=-CFL*temp3
            call electron_prepare_implicit_coeffs_common(Num_gam_e,one,up,principal,temp1)
            temp2=(dN_x+dDR*dF1)/principal
            call electron_backward_sweep_common(Num_gam_e,temp1,temp2,x)
            dN_x=x
            if (L == L1) dN_gam_e(:,I_tobs)=dN_x/gam_e/dlog(ten)
        end do
    end do

    deallocate(dEl,principal,x,dF1,up,temp1,temp2,temp3,dN_x,para_minus_gam_e_p,dB3_serial,P_syn,Seed_syn,cooling_aux,Compton)
end subroutine electron_reverse_evolve

end module electron_reverse_kernel
