include 'Constants.f90'
subroutine dynamics_cip(Boundary,Num_gam_e, R_Tobs,R_Gamma,R,gam_e,dN_gam_e)
    !$ use omp_lib
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    PARAMETER(Num_R=900,Num_R1=Num_R-1)
    PARAMETER(Num_RK4=4)
    EXTERNAL F

    real(8),dimension (15) :: Boundary
    real(8),dimension (Num_RK4) :: Y,D,B,C,G,E
    real(8),dimension (Num_R) :: R_Tobs,R_Gamma,R_m,R
    real(8),dimension (Num_R,Num_gam_e) :: dN_gam_e
    real(8),dimension (Num_gam_e) :: gam_e
    real(8),allocatable,dimension (:) :: temp_N_e,f_r_times_gam_e,dEl,para_minus_gam_e_p,dEl1, &
            x,temp_11,dN_x,ddN_x,temp_12,cal_para_minus_gam_e_p,dF1,temp_a,temp_b

    integer,allocatable,dimension (:) :: NP

    allocate (temp_N_e(Num_gam_e),f_r_times_gam_e(Num_gam_e),dEl(Num_gam_e), &
            para_minus_gam_e_p(Num_gam_e),dEl1(Num_gam_e),x(Num_gam_e),temp_11(Num_gam_e), &
            dN_x(Num_gam_e),ddN_x(Num_gam_e),temp_12(Num_gam_e),cal_para_minus_gam_e_p(Num_gam_e), &
            dF1(Num_gam_e),temp_a(Num_gam_e),temp_b(Num_gam_e),NP(Num_gam_e))

    !f2py intent(in) :: Boundary,Num_gam_e
    !f2py intent(out) :: R_Tobs,R_Gamma,R,gam_e,dN_gam_e
    !f2py depend(Num_gam_e) :: gam_e,dN_gam_e

    !***********************[Parameter Initial]**********************
    Eta_0=Boundary(1)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    E_iso=Boundary(14)
    T_log10_duration=Boundary(15)
    !***********************[Parameter Initial]**********************
    !   Y(1)=Gamma,Y(2)=m,Y(3)=thermal energy,Y(4)=R

    Y(1)=Eta_0-0.001
    Y(2)=Boundary(2)
    Y(3)=Boundary(3)
    Y(4)=Boundary(4)
    T00=Y(4)*(one/dsqrt(one-one/Eta_0**2)-one)/Para_c
    M=4
    EPS=1.0D-5
    !**********************[Time bin]**********************
    DM_0=E_iso/((Eta_0-one)*4d0*pi*Para_m_p_E)
    R_dec_ISM=(dNe_ISM*Eta_0/DM_0)**(-one/3d0)
    if (A_star > zero) then
        R_dec_wind=DM_0/(2.0d35*A_star*Eta_0)
        R_dec=min(R_dec_wind,R_dec_ISM)
    else
        R_dec=R_dec_ISM
    end if
    t_dec=R_dec/(two*Eta_0*Eta_0*Para_c)
    Grid_Tobs_bin=min(-5d0,dlog10(t_dec*0.1))
    T_log10=T_log10_duration-Grid_Tobs_bin
    !   log time bin
    do I_tobs=1,Num_R
        T=T00+ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        if (I_tobs < one) then
            H=ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        else
            H=ten**(Grid_Tobs_bin+T_log10*I_tobs/Num_R1)-ten**(Grid_Tobs_bin+T_log10*(I_tobs-one)/Num_R1)
        end if

        call GRKT2(T,H,Y,M,F,EPS,D,B,C,G,E,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Epsilon_b,p,z)
        R_Tobs(I_tobs)=T*(one+z)
        R_Gamma(I_tobs)=Y(1)
        R_m(I_tobs)=Y(2)/Para_m_p
        R(I_tobs)=Y(4)
    end do

    !*****************Part 1: given the boundary consition [Using the analytical approximation]*********************
    if (A_star >= zero) then
        dNe_wind=A_star*3.0d35/R(1)**2
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
    end if
    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe_ISM*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    Gam_e_m=(p-two)/(p-one)*Epsilon_e*1836d0*(R_Gamma(1)-one)+one
    if (p-2.05 < 0.001) then
        Gam_e_m=0.05d0/1.05d0*Epsilon_e*1836d0*(R_Gamma(1)-one)+one
    end if
    Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    do I_gam_e=1,Num_gam_e
        Gam_e(I_gam_e)=3d0*ten**(dlog10(Gam_e_max_max)*(I_gam_e-1)/(Num_gam_e-1))
        if (Gam_e_m > Gam_e_c) then
            Q1=R_m(1)*Gam_e_c
            if ((Gam_e_c-Gam_e(I_gam_e)) > one) then
                dN_gam_e(1,I_gam_e)=zero
            else
                if (Gam_e_m > Gam_e(I_gam_e)) then
                    dN_gam_e(1,I_gam_e)=Q1*Gam_e(I_gam_e)**(-2)
                else
                    if (Gam_e_max > Gam_e(I_gam_e)) then
                        dN_gam_e(1,I_gam_e)=Q1*Gam_e_m**(p-one)*Gam_e(I_gam_e)**(-(p+one))
                    else
                        dN_gam_e(1,I_gam_e)=zero
                    end if
                end if
            end if
        else
            Q1=R_m(1)*Gam_e_m**(p-one)
            if (Gam_e_m > Gam_e(I_gam_e)) then
                dN_gam_e(1,I_gam_e)=zero
            else
                if (Gam_e_c > Gam_e(I_gam_e)) then
                    dN_gam_e(1,I_gam_e)=Q1*Gam_e(I_gam_e)**(-p)
                else
                    if (Gam_e_max > Gam_e(I_gam_e)) then
                        dN_gam_e(1,I_gam_e)=Q1*Gam_e_c*Gam_e(I_gam_e)**(-(p+one))
                    else
                        dN_gam_e(1,I_gam_e)=zero
                    end if
                end if
            end if
        end if
    end do
    dN_x=dN_gam_e(1,:)*gam_e*dlog(ten)
    d_x=dlog10(gam_e(2)/gam_e(1))
    d_x2=d_x*d_x
    d_x3=d_x2*d_x
    ddN_x(2:Num_gam_e)=(dN_x(2:Num_gam_e)-dN_x(1:Num_gam_e-1))/d_x
    ddN_x(1)=ddN_x(2)

    !*******************Part 1 is completed [has been checked and there is no bug]**********************************
    !*******************Part 2: To calculate the electron distribution**********************************************
    para_minus_gam_e_p=one/(gam_e-one)**p
    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        if (A_star >= 0.0) then
            dNe_wind=A_star*3.0d35/R_loc**2
            if (dNe_wind <= dNe_ISM/4d0) then
                dNe=dNe_ISM
            else
                dNe=dNe_wind
            end if
        else
            dNe=dNe_ISM
        end if
        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        Gam_e_m=(p-two)/(p-one)*Epsilon_e*1836d0*(R_Gamma_loc-one)+one
        if (p-2.05 < 0.001) then
            Gam_e_m=0.05d0/1.05d0*Epsilon_e*1836d0*(R_Gamma_loc-one)+one
        end if
        Gam_e_m_p=(p-one)*(Gam_e_m-one)**(p-one)
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)
        eta=(Gam_e_m/Gam_e_c)**(p-two)
        if (eta-one > 0.001) eta=one
        beta_Gam=dsqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        f_r_times_gam_e=f_r*gam_e
        dDR=3d-1*dlog(ten)*d_x/(f_r*Gam_e_max+one/R(I_tobs))
        !***********************[Here we have presented the choice on Delta_r]******************************************
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=Int(dDD/dDR)
        L1=L1+10
        dDR=dDD/L1
        temp_N_e=dN_gam_e(I_tobs-1,:)
        dN_x=temp_N_e*gam_e*dlog(ten)
        !$OMP PARALLEL  NUM_THREADS(8)
        !$OMP DO
        do i_gam_e=1,Num_gam_e
            if (gam_e(i_gam_e) < Gam_e_max .and. gam_e(i_gam_e) > Gam_e_m) then
                cal_para_minus_gam_e_p(i_gam_e)=para_minus_gam_e_p(i_gam_e)
            else
                cal_para_minus_gam_e_p(i_gam_e)=zero
            end if
            if (i_gam_e < Num_gam_e) then
                !*****************************[The general inverse Compton effect]******************************
                hat_gam=6.5d6/sqrt(DB*(gam_e(i_gam_e)+gam_e(i_gam_e+1))/two)
                if (Gam_e_m > Gam_e_c) then
                    if (hat_gam < Gam_e_c) then
                        eta_NK=zero
                    else
                        if (hat_gam < Gam_e_m) then
                            Step1=(p-one)/(p-two)*Gam_e_m-Gam_e_c
                            eta_NK=(hat_gam-Gam_e_c)/Step1
                        else
                            Step2=Gam_e_m**(p-one)*hat_gam**(two-p)
                            Step3=(p-one)*Gam_e_m-(p-two)*Gam_e_c
                            eta_NK=one-Step2/Step3
                        end if
                    end if
                else
                    if (hat_gam < Gam_e_m) then
                        eta_NK=zero
                    else
                        if (hat_gam < Gam_e_c) then
                            Step4=Gam_e_c**(3d0-p)/(p-two)-Gam_e_m**(3d0-p)
                            eta_NK=(hat_gam**(3d0-p)-Gam_e_m**(3d0-p))/Step4
                        else
                            Step5=(3d0-p)*Gam_e_c*hat_gam**(two-p)
                            Step6=Gam_e_c**(3d0-p)-(p-two)*Gam_e_m**(3d0-p)
                            eta_NK=one-Step5/Step6
                        end if
                    end if
                end if
                Compton=one+(-one+dsqrt(one+4d0*eta*eta_NK*Epsilon_e/Epsilon_b))/two
                !                Compton=1d0
                dEl(i_gam_e)=f_r_times_gam_e(i_gam_e)*Compton
            else
                dEl(i_gam_e)=f_r_times_gam_e(i_gam_e)*Compton
            end if
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        do L=1,L1
            R_loc=R_loc+dDR
            dEl1=(dEl+one/R_loc)/log(ten) !
            if (A_star >= zero) then
                dNe_wind=A_star*3.0d35/R_loc**2
                if (dNe_wind <= dNe_ISM/4d0) then
                    dNe=dNe_ISM
                else
                    dNe=dNe_wind
                end if
            else
                dNe=dNe_ISM
            end if
            Q=4d0*pi*R_loc*R_loc*dNe*Gam_e_m_p  !here Q is Q_0*\gamma_m**p
            dF1=Q*cal_para_minus_gam_e_p*gam_e*dlog(ten)

            temp_12=dN_x
            temp_11=ddN_x
            dN_x=temp_12+temp_12*dEl*dDR+dF1*dDR
            dN_x=0.75d0*temp_12+0.25d0*dN_x+0.25d0*(dN_x*dEl*dDR+dF1*dDR)
            dN_x=temp_12/3d0+2d0*dN_x/3d0+2d0*(dN_x*dEl*dDR+dF1*dDR)/3d0

            do i_gam_e=2,Num_gam_e-1
                ddN_x(i_gam_e)=temp_11(i_gam_e)+(dN_x(i_gam_e+1)- &
                        dN_x(i_gam_e-1)-temp_12(i_gam_e+1)+temp_12(i_gam_e-1))/d_x/two+ddN_x(i_gam_e)*dEl(i_gam_e)*dDR
            end do
            ddN_x(1)=ddN_x(2)
            ddN_x(Num_gam_e)=ddN_x(Num_gam_e-1)

            x=dEl1*dDR

            do i_gam_e=1,Num_gam_e-1
                temp_a(i_gam_e)=(ddN_x(i_gam_e)+ddN_x(i_gam_e+1))/d_x2+two*(dN_x(i_gam_e)-dN_x(i_gam_e+1))/d_x3
                temp_b(i_gam_e)=3d0*(dN_x(i_gam_e+1)-dN_x(i_gam_e))/d_x2-(two*ddN_x(i_gam_e)+ddN_x(i_gam_e+1))/d_x
            end do
            temp_a(Num_gam_e)=temp_a(Num_gam_e-1)
            temp_b(Num_gam_e)=temp_b(Num_gam_e-1)

            temp_12=dN_x
            temp_11=ddN_x
            dN_x=x*(x*(temp_a*x+temp_b)+temp_11)+temp_12
            ddN_x=x*(3d0*temp_a*x+two*temp_b)+temp_11

            temp_12=dN_x
            temp_11=ddN_x
            dN_x=temp_12+temp_12*dEl*dDR+dF1*dDR
            dN_x=0.75d0*temp_12+0.25d0*dN_x+0.25d0*(dN_x*dEl*dDR+dF1*dDR)
            dN_x=temp_12/3d0+two*dN_x/3d0+two*(dN_x*dEl*dDR+dF1*dDR)/3d0

            isgn=1
            ee=0.01
            call Nindex (NP, Num_gam_e, dN_x, ee)

            do i=2,Num_gam_e-1
                if (NP(i)==-isgn .and. (x(i+1) - x(i-1) >= 0)) then
                    call Linear(dN_x(i),ddN_x(i),x(i),dEl1(i-1),dN_x(i+isgn),ddN_x(i+isgn),dN_x(i),ddN_x(i),1)
                end if

                if (NP(i+isgn)==isgn .and. (x(i+1) -x(i-1) >= 0)) then
                    call Linear(dN_x(i),ddN_x(i),x(i),dEl1(i-1),dN_x(i+isgn),ddN_x(i+isgn),dN_x(i),ddN_x(i),-1)
                end if
            end do

            do i_gam_e=1,Num_gam_e-1
                if (dN_x(i_gam_e) < 0.0) then
                    dN_x(i_gam_e)=0.0
                end if
            end do

        end do
        temp_N_e=dN_x/gam_e/log(ten)
        !        write (*,*) temp_N_e(20)
        do i_gam_e=1,Num_gam_e
            if (temp_N_e(i_gam_e) < zero) then
                dN_gam_e(I_tobs,i_gam_e)=zero
            else
                dN_gam_e(I_tobs,i_gam_e)=temp_N_e(i_gam_e)
            end if
        end do
    end do

    return
end subroutine

subroutine Nindex(NP, Num_gam_e, dN_x, ee)
    IMPLICIT REAL(8)(A-H,O-Z)
    DIMENSION NP(Num_gam_e),dN_x(Num_gam_e)

    do i=3,Num_gam_e-2
        if (abs(dN_x(i)-dN_x(i-1)/(dN_x(i+1)-dN_x(i))) < ee .and. &
                abs(dN_x(i-1)-dN_x(i-2)/(dN_x(i+1)-dN_x(i))) < ee) then
            NP(i)=-1
        else if (abs(dN_x(i)-dN_x(i-1)/(dN_x(i+2)-dN_x(i+1))) >= 1./ee .and. &
                abs(dN_x(i)-dN_x(i-1)/(dN_x(i+1)-dN_x(i))) >= 1./ee) then
            NP(i)=1
        end if
    end do

    return
end subroutine Nindex

subroutine Linear(f_star, f_hat_star, xi, dr, fu, fu_hat, fd, fd_hat, ip)
    IMPLICIT REAL(8)(A-H,O-Z)

    f1=fu-(dr-xi)*fu_hat
    f2=fd+xi*fd_hat
    if (ip*(fd-fu)*(f2-f1) < 0.0 ) then
        f_star=f2
        f_hat_star=fd_hat
    else
        f_star=f1
        f_hat_star=fu_hat
    end if

    return
end subroutine Linear

!**********************[Dynamic]**********************
SUBROUTINE F(T, Y, M, D, E_e, E_iso, Eta_0, dNe_ISM, A_star, E_b, p, z)
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    DIMENSION Y(M),D(M)

    if (A_star >= 0.0) then
        dNe_wind=A_star*3.0d35/Y(4)**2
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
    end if

    Epe=E_e
    qq=0.0
    t_1=700.0/(1d0+z)
    t_2=1500.0/(1d0+z)
    if (T >= t_1 .and. T<= t_2) then
        A=1.07d50
    else
        A=0.0
    end if
    A=zero
    dB=0.39d0*dsqrt((E_b*dNe)*(Y(1)**2-one))
    gam_c=7.739d8/(dB**2*Y(1)*T)
    gam_m=Epe*1836d0*(p-two)*(Y(1)-one)/(p-one)+one
    if ((gam_c-gam_m) > 0.001) Epe=Epe*(gam_m/gam_c)**(p-two)
    Bgam=dsqrt(one-one/Y(1)**2)

    ! **************************Huang's Dynamic model******************************
    !     DNe0=E_iso/(Eta_0-1.0)/1.5*1.0D3/(4d0*pi)
    !     D00=Bgam*2.997D10/(1.0-Bgam)
    !     D01=Y(1)**2-1.0
    !     D02=DNe0+Epe*Y(2)+2.0*(1.0-Epe)*Y(1)*Y(2)
    !     dM=dNe*Y(4)**2*D00
    !     D(1)=(-D01*dM+A/(4d0*pi)*(1.0+T/T_1)**qq/1.5*1.0D3)/D02
    !     D(2)=dM
    !     D(3)=(Y(1)-1.0)*dM
    !     D(4)=D00
    !! ************************Pe'er's Dynamic model*******************************
    ! *******************Caculate hat_gamma derived by Pe'er **********************
    !             four_v=Bgam*Y(1)
    !             theta=four_v/3.0*(four_v+1.07*four_v**2)
    !    $        /(1.0+four_v+1.07*four_v**2)
    !             cz=theta/(0.24+theta)
    !             hat_gam=(5.0-1.21937*cz+0.18203*cz**2-0.96583*cz**3
    !    $        +2.32513*cz**4-2.39332*cz**5+1.07136*cz**6)/3.0
    !!! *************************new equations by Pe'er*****************************
    !             DNe0=E_iso/(Eta_0-1.0)/1.5*1.0D3/12.5664
    !             D00=Bgam*2.997D10/(1.0-Bgam)
    !             D01=Y(1)**2-1.0
    !             D011=hat_gam*D01-(hat_gam-1.0)*Y(1)*Bgam**2
    !             D02=DNe0+Epe*Y(2)+2.0*(1.0-Epe)*Y(1)*Y(2)
    !             D022=DNe0+Epe*Y(2)+(1.0-Epe)*Y(2)*(2.0*hat_gam*Y(1)
    !    $         -(hat_gam-1.0)*(1.0+Y(1)**(-2)))
    !             dM=DNe*Y(4)**2*D00
    !             D(1)=-(D011/D022)*dM
    !             D(2)=dM
    !             D(3)=(Y(1)-1.0)*dM
    !             D(4)=D00
    !! ************************Zhang's Dynamic model*******************************
    ! *******************Caculate hat_gamma derived by Pe'er **********************
    four_v=Bgam*Y(1)
    theta=four_v/3d0*(four_v+1.07d0*four_v**2)/(one+four_v+1.07d0*four_v**2)
    cz=theta/(0.24d0+theta)
    hat_gam=(5d0-1.21937d0*cz+0.18203d0*cz**2-0.96583d0*cz**3+ &
            2.32513d0*cz**4-2.39332d0*cz**5+1.07136d0*cz**6)/3d0
    ! *************************new equations by Zhang******************************
    DNe0=E_iso/(Eta_0-one)/Para_c**2
    D00=Bgam*Para_c/(one-Bgam)
    D01=Y(1)**2-one
    dM=4d0*pi*dNe*Y(4)**2*D00*Para_m_p
    D011=-dM/D00*Para_c**2*Y(1)*D01*(hat_gam*Y(1)-hat_gam+one)- &
            (hat_gam*D01+one)*Y(1)*(one-hat_gam)*3d0*Y(3)/Y(4)+Y(1)**2*A*(one+T/t_1)**qq/D00
    D022=Y(1)**2*(DNe0+Y(2))*Para_c**2+(hat_gam**2*D01+3d0*hat_gam-two)*Y(3)

    D(1)=D011/D022*D00
    D(2)=dM
    D(3)=((1d0-Epe)*(Y(1)-one)*dM/D00*Para_c**2- &
            (hat_gam-one)*(3d0/Y(4)-one/Y(1)*D(1)/D00)*Y(3))*D00+A*(one+T/t_1)**qq/D00
    D(4)=D00
    return
end subroutine

SUBROUTINE GRKT2(T,H,Y,M,F,EPS,D,B,C,G,E,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0)
    IMPLICIT REAL(8)(A-H,O-Z)
    DIMENSION Y(M),D(M),A(4),B(M),C(M),G(M),E(M)

    HH=H
    N=1
    P=1+EPS
    X=T
    C=Y
    do while (P >= EPS)
        A(1)=HH/2.0
        A(2)=A(1)
        A(3)=HH
        A(4)=HH
        G=Y
        Y=C
        DT=H/N
        T=X
        do J=1,N
            call F(T,Y,M,D,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0)
            E=Y
            B=Y
            do K=1,3
                Y=E+A(K)*D
                B=B+A(K+1)*D/3.0
                TT=T+A(K)
                call F(TT,Y,M,D,Epsilon_e,E_iso,Eta_0,dNe_ISM,A_star,Eb,pp,z0)
            end do
            Y=B+HH*D/6.0
            T=T+DT
        end do
        P=0.0
        do I=1,M
            Q=2.0*abs(Y(I)-G(I))/(Y(I)+G(I))
            if (Q > P) P=Q
        end do
        HH=HH/2.0
        N=N+N
    end do
    T=X

    return
end subroutine GRKT2

