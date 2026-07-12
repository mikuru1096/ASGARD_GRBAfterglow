subroutine dynamics_forward(Boundary,n,Num_R,index_dyn,R_Tobs,R_Gamma,R,R_m)
    use constants
    use dynamics_density_profile, only: density_moment, set_density_profile
    implicit none
    integer, intent(in) :: n,Num_R,index_dyn
    integer, parameter :: r0_slot = 27
    integer :: I_tobs, Num_R1, M
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(out), dimension(Num_R) :: R_Tobs,R_Gamma,R,R_m
    real(8) :: Eta_0,Epsilon_e,Epsilon_b,p,z,dNe_ISM,A_star,E_iso,tdur_log,f_e
    real(8) :: einj1,einj2,E_inj,einj_q,R_tr,f_jump,f_wide,R0,EPS,dm0,rdec,tdec,gam0
    real(8) :: tbin,tspan,T,H,dtgrid,rdecism,rdecwind,mstart
    real(8),dimension(4) :: Y,D,B,C,E

    ! Boundary 只在入口解包一次；后续计算使用具名标量。
    ! Boundary is unpacked once at the entry; later stages use named scalars.
    Eta_0=Boundary(1); Epsilon_e=Boundary(5); Epsilon_b=Boundary(6); p=Boundary(7); z=Boundary(8)
    dNe_ISM=Boundary(11); A_star=Boundary(12); E_iso=Boundary(14); tdur_log=Boundary(15); f_e=Boundary(16)
    einj1=Boundary(17); einj2=Boundary(18); E_inj=Boundary(19); einj_q=Boundary(20)
    R_tr=Boundary(21); f_jump=Boundary(22); f_wide=Boundary(23)
    if (n >= r0_slot) then
        R0=Boundary(r0_slot)
    else
        R0=Boundary(n)
    end if
    call set_density_profile(Boundary,n)
    call density_moment(A_star,dNe_ISM,Boundary(4),R0,R_tr,f_jump,f_wide,mstart)
    if (index_dyn == 3) mstart=4d0*pi*Para_m_p*mstart

    Num_R1=Num_R-1
    gam0=Eta_0-0.001d0
    ! Y(1) 保存正四速度 Gamma*beta，使 RK 状态跨越相对论与非相对论区。
    ! Y(1) stores positive four-velocity Gamma*beta across both dynamical regimes.
    Y=[dsqrt((gam0-1d0)*(gam0+1d0)),mstart,Boundary(3),Boundary(4)]
    M=4; EPS=1.0D-5

    ! 减速时间只决定 observer-time 网格；动力学仍由 RK 状态推进。
    ! The deceleration time only sets the observer-time grid; dynamics still come from RK advancement.
    dm0=E_iso/((Eta_0-1d0)*4d0*pi*Para_m_p_E)
    rdecism=(dNe_ISM*Eta_0/dm0)**(-1d0/3d0)
    if (A_star > 0d0) then
        rdecwind=dm0/(2.0d35*A_star*Eta_0)
        rdec=min(rdecwind,rdecism)
    else
        rdec=rdecism
    end if
    tdec=rdec/(2d0*Eta_0*Eta_0*Para_c)
    tbin=min(-5d0,dlog10(tdec*0.1d0))
    tspan=tdur_log-tbin

    ! 每个输出点从当前时间推进到目标观测时间，并写公开输出数组。
    ! Each output step advances from the current time to the target observer time.
    T=0d0
    do I_tobs=1,Num_R
        dtgrid=1d1**(tbin+tspan*(I_tobs-1d0)/Num_R1)
        H=dtgrid-T
        call forward_rk4()
        R_Tobs(I_tobs)=T*(1d0+z); R_Gamma(I_tobs)=dsqrt(1d0+Y(1)*Y(1)); R(I_tobs)=Y(4)
        if (index_dyn == 3) then
            R_m(I_tobs)=Y(2)
        else
            R_m(I_tobs)=4d0*pi*Para_m_p*Y(2)
        end if
    end do

contains

! 自适应 RK4：正时间用 log(time) 步长，初始零时间段用线性时间步长。
! Adaptive RK4: use log(time) steps after positive time, and linear steps near the initial 0d0 time.
subroutine forward_rk4()
    implicit none
    integer :: N, J, K
    real(8) :: HH, DT, P, X, TT, S, HS, SS, T_phys
    real(8), parameter :: logstep=0.25d0
    logical :: trial_ok
    real(8), dimension(4) :: A, ref, trial

    X = T
    C = Y
    ref = C
    if (X > 0d0) then
        HS = log((X+H)/X)
        N = max(1, ceiling(abs(HS)/logstep))
        HH = HS/N
        P = 1d0+EPS
        do while (P >= EPS)
            A = [0.5d0*HH, 0.5d0*HH, HH, HH]
            Y = C
            S = 0d0
            trial_ok = .true.
            do J = 1, N
                T_phys = X*exp(S)
                call forward_rhs(T_phys, Y, D)
                D = T_phys*D
                E = Y
                B = Y
                do K = 1, 3
                    Y = E+A(K)*D
                    trial_ok = Y(1) > 0d0
                    if (.not. trial_ok) exit
                    B = B+A(K+1)*D/3.0d0
                    SS = S+A(K)
                    T_phys = X*exp(SS)
                    call forward_rhs(T_phys, Y, D)
                    D = T_phys*D
                end do
                if (.not. trial_ok) exit
                Y = B+HH*D/6.0d0
                trial_ok = Y(1) > 0d0
                if (.not. trial_ok) exit
                S = S+HH
            end do
            if (trial_ok) then
                trial = Y
                call forward_error(trial, ref, P, trial_ok)
            else
                P = huge(1d0)
            end if
            if (trial_ok) ref = trial
            HH = 0.5d0*HH
            N = N+N
        end do
        T = X+H
        return
    end if

    HH = H
    N = 1
    P = 1d0+EPS
    do while (P >= EPS)
        A = [0.5d0*HH, 0.5d0*HH, HH, HH]
        Y = C
        T = X
        DT = HH
        trial_ok = .true.
        do J = 1, N
            call forward_rhs(T, Y, D)
            E = Y
            B = Y
            do K = 1, 3
                Y = E+A(K)*D
                trial_ok = Y(1) > 0d0
                if (.not. trial_ok) exit
                B = B+A(K+1)*D/3.0d0
                TT = T+A(K)
                call forward_rhs(TT, Y, D)
            end do
            if (.not. trial_ok) exit
            Y = B+HH*D/6.0d0
            trial_ok = Y(1) > 0d0
            if (.not. trial_ok) exit
            T = T+DT
        end do
        if (trial_ok) then
            trial = Y
            call forward_error(trial, ref, P, trial_ok)
        else
            P = huge(1d0)
        end if
        if (trial_ok) ref = trial
        HH = HH/2.0d0
        N = N+N
    end do
    T = X+H
end subroutine forward_rk4

subroutine forward_error(trial, ref, error, valid)
    implicit none
    integer :: I
    real(8), intent(in), dimension(M) :: trial, ref
    real(8), intent(out) :: error
    logical, intent(out) :: valid
    real(8) :: ratio

    error = 0d0
    valid = .true.
    do I = 1, M
        if (trial(I) <= 0d0 .or. ref(I) <= 0d0) then
            valid = .false.
            error = huge(1d0)
            return
        end if
        ratio = 2d0*abs(trial(I)-ref(I))/(trial(I)+ref(I))
        if (ratio > error) error = ratio
    end do
end subroutine forward_error

! 正向激波 RHS：外介质密度、能量注入、辐射修正和动力学分支都在这里定义。
! Forward-shock RHS: external density, energy injection, radiative correction, and dynamics branch live here.
subroutine forward_rhs(T_rhs, Y_rhs, D_rhs)
    use dynamics_density_profile, only: density_profile
    implicit none
    real(8), intent(in), dimension(M) :: Y_rhs
    real(8), intent(in) :: T_rhs
    real(8), intent(out), dimension(M) :: D_rhs
    real(8), parameter :: bsyn_coeff=0.39d0, gcool_coeff=7.739d8
    real(8) :: dNe,t_1,t_2,A,Epe,dB,gam_c,gam_m,D00,D01,dM,hat_gam,four_v,theta,cz
    real(8) :: gamma,gamone,gamdot,DNe0,D011,D02,D022

    call density_profile(A_star,dNe_ISM,Y_rhs(4),R0,1,R_tr,f_jump,f_wide,dNe)
    t_1=einj1/(1d0+z); t_2=einj2/(1d0+z)
    if (T_rhs >= t_1 .and. T_rhs <= t_2) then
        A=E_inj*(T_rhs/t_1)**einj_q
    else
        A=0d0
    end if

    ! u=Gamma*beta: Gamma-1=u^2/(Gamma+1), beta/(1-beta)=u*(Gamma+u).
    ! 这些恒等式在深度非相对论区避免相消，并保持 RHS 定义域。
    four_v=Y_rhs(1)
    D01=four_v*four_v
    gamma=dsqrt(1d0+D01)
    gamone=D01/(gamma+1d0)
    D00=four_v*(gamma+four_v)*Para_c

    Epe=Epsilon_e
    dB=bsyn_coeff*dsqrt((Epsilon_b*dNe)*D01)
    gam_c=gcool_coeff/(dB**2*gamma*T_rhs)
    gam_m=Epe/f_e*Para_m_p_DIV_m_e*(p-2d0)*gamone/(p-1d0)+1d0
    if ((gam_c-gam_m) > 0.001d0 .and. p >= 2d0) Epe=Epe*(gam_m/gam_c)**(p-2d0)

    if (index_dyn == 3) then
        dM=4d0*pi*dNe*Y_rhs(4)**2*D00*Para_m_p
    else
        dM=dNe*Y_rhs(4)**2*D00
    end if
    if (index_dyn == 2 .or. index_dyn == 3) then
        theta=four_v/3d0*(four_v+1.07d0*four_v**2)/(1d0+four_v+1.07d0*four_v**2)
        cz=theta/(0.24d0+theta)
        hat_gam=(5d0-1.21937d0*cz+0.18203d0*cz**2-0.96583d0*cz**3+2.32513d0*cz**4-2.39332d0*cz**5+1.07136d0*cz**6)/3d0
    end if
    D_rhs(2)=dM; D_rhs(4)=D00

    select case(index_dyn)
    case(1)
        DNe0=E_iso/(Eta_0-1d0)/1.5d0*1.0D3/(4d0*pi)
        D02=DNe0+Epe*Y_rhs(2)+2d0*(1d0-Epe)*gamma*Y_rhs(2)
        gamdot=(-D01*dM+A/(4d0*pi)/1.5d0*1.0D3)/D02
        D_rhs(3)=gamone*dM
    case(2)
        DNe0=E_iso/(Eta_0-1d0)/1.5d0*1.0D3/(4d0*pi)
        D011=D01*(hat_gam-(hat_gam-1d0)/gamma)
        D022=DNe0+Epe*Y_rhs(2)+(1d0-Epe)*Y_rhs(2)*(2d0*hat_gam*gamma- &
             (hat_gam-1d0)*(1d0+gamma**(-2)))
        gamdot=-(D011/D022)*dM
        D_rhs(3)=gamone*dM
    case(3)
        DNe0=E_iso/(Eta_0-1d0)/Para_c**2
        D011=-dM/D00*Para_c**2*gamma*D01*(hat_gam*gamma-hat_gam+1d0)- &
             (hat_gam*D01+1d0)*gamma*(1d0-hat_gam)*3d0*Y_rhs(3)/Y_rhs(4)+gamma**2*A/D00
        D022=gamma**2*(DNe0+Y_rhs(2))*Para_c**2+(hat_gam**2*D01+3d0*hat_gam-2d0)*Y_rhs(3)
        gamdot=D011/D022*D00
        D_rhs(3)=((1d0-Epe)*gamone*dM/D00*Para_c**2- &
                 (hat_gam-1d0)*(3d0/Y_rhs(4)-gamdot/(gamma*D00))*Y_rhs(3))*D00+A
    end select
    D_rhs(1)=gamma/four_v*gamdot
end subroutine forward_rhs
end subroutine dynamics_forward
