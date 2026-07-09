! 反向激波动力学 f2py 入口。
! Reverse-shock dynamics f2py entry.
! first RS 是第 0 个 ejecta 分支；density jump 触发后续分支，同一输出循环内同步推进和记录。
! The first RS is branch 0; density jumps trigger later branches in the same output loop.
subroutine dynamics_reverse(Delta_t,e_r,b_r,p_r,fer,sigma_r,Boundary,n,Num_R, &
                            T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gmcross,b3ordx, &
                            R_Tobs,R_Gamma,R,M2,M3,B3,U3_th,V3_comoving,Gamma34_inst, &
                            branch_m3,branch_u3,branch_v3,branch_b3, &
                            total_m3,total_u3,total_v3,total_b3, &
                            total_p,total_w, &
                            avg_gc,avg_p3,avg_g43,avg_brs, &
                            total_ud,diss_e,inj_e, &
                            branch_gm,branch_gc,branch_g43, &
                            branch_comp,branch_brs,branch_ud, &
                            nu_m,nu_c, &
                            event_on,start_r,end_r, &
                            start_t,end_t)
    use constants
    use dynamics_density_profile, only: set_density_profile, density_profile, density_moment, &
                                        jump_count, jump_radius, &
                                        jump_factor, jump_width
    use reverse_jump_conditions, only: reverse_contact
    use reverse_shock_mhd_jump, only: rs_vegas_ud, rs_vegas_comp, rs_mag_internal
    use reverse_shock_state, only: rhs_phase, wait_phase, precross_phase, &
                                   postcross_phase, rs_db3, rs_tcross, rs_rcross, rs_e3cross, &
                                   rs_gam20, rs_u3cross, rs_v3cross, rs_m3cross, rs_gammcross, &
                                   rs_b3ordered, rs_nstate
    use reverse_rhs_module, only: reverse_dynamics_rhs
    implicit none
    integer, intent(in) :: n,Num_R
    integer, parameter :: ir0 = 27, jmax = 8
    integer :: nstate,it,event_j,njump,Num_R1
    real(8), dimension(n), intent(in) :: Boundary
    real(8), intent(in) :: Delta_t,e_r,b_r,p_r,fer,sigma_r
    real(8), intent(out) :: T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gmcross,b3ordx
    real(8), dimension(Num_R), intent(out) :: R_Tobs,R,M2,M3,B3,R_Gamma,U3_th,V3_comoving,Gamma34_inst
    real(8), dimension(Num_R), intent(out) :: total_m3,total_u3,total_v3,total_b3,total_p,total_w
    real(8), dimension(Num_R), intent(out) :: avg_gc,avg_p3,avg_g43,avg_brs,total_ud,diss_e,inj_e,nu_m,nu_c
    real(8), dimension(jmax,Num_R), intent(out) :: branch_m3,branch_u3,branch_v3,branch_b3,branch_gm
    real(8), dimension(jmax,Num_R), intent(out) :: branch_gc,branch_g43,branch_comp,branch_brs,branch_ud
    logical, dimension(jmax), intent(out) :: event_on
    real(8), dimension(jmax), intent(out) :: start_r,end_r,start_t,end_t
    real(8) :: Eta_0,Epsilon_e,Epsilon_b,p_f,z,nism,astar,eiso,tdur,f_e
    real(8) :: R_tr,f_jump,f_wide,R0
    real(8) :: delta0,mej,v3s,m2init,m3init,dmref,rdec,tbase,tdec,tbin,tspan,dB3
    real(8) :: dtgrid,rdecism,rdecwind,tnow,tout,tevent,dbwait,dbevent
    real(8) :: u2init,u4init,dinit,n4init,g34init,n3init,compinit
    real(8) :: prev_r,prev_g,prev_t,curr_r,curr_g,curr_t
    real(8), dimension(jmax) :: src_prev,src_cur,jump_r,jump_f,jump_w,diss_prev,gm_prev
    real(8), dimension(:), allocatable :: y_prev,y_cur,Y,ywait,yevent
    logical, dimension(jmax) :: event_done
    logical :: ready

    ! Y 状态向量布局：
    ! Y state-vector layout:
    !   1: 接触面 Lorentz 因子，2: 半径，3: 扫过的外介质质量。
    !   1: contact Lorentz factor, 2: radius, 3: swept external mass.
    !   4: first RS 扫过 ejecta 分数，5: first RS U3/(Mej c^2)，6: first RS V3/v3s。
    !   4: first-RS swept ejecta fraction, 5: first-RS U3/(Mej c^2), 6: first-RS V3/v3s.
    !   每个 density jump 后面五段依次存 multiple-RS 分支质量、热能、体积、累计耗散能、
    !   以及 gamma_m 加权的累计耗散能。
    !   For each density jump, the next five blocks store branch mass, thermal energy,
    !   volume, cumulative dissipated energy, and gamma_m-weighted cumulative dissipated energy.
    call unpack()
    call set_density_profile(Boundary,n)
    njump=jump_count
    jump_r=jump_radius
    jump_f=jump_factor
    jump_w=jump_width
    nstate=6+5*njump
    allocate(y_prev(nstate),y_cur(nstate),Y(nstate),ywait(nstate),yevent(nstate))

    Y=0d0
    call init_shell()
    call init_dec()
    call init_cross()
    tbase=Y(2)*(1d0/dsqrt(1d0-1d0/Eta_0**2)-1d0)/Para_c
    tdec=rdec/(2d0*Eta_0*Eta_0*Para_c)
    tbin=min(log10(tbase)-2d0,dlog10(tdec*0.1d0))
    tspan=tdur-tbin
    Num_R1=Num_R-1
    tnow=tbase
    call init_output()
    call init_event()
    call init_scan()

    ! 观测时间网格只走一遍：先推进 first RS，写 first 输出，再用同一对相邻状态扫描
    ! density-jump multiple-RS 事件和分支诊断。
    ! Single pass over the observer grid: advance first RS, store first outputs,
    ! then scan multiple-RS events and diagnostics from the same adjacent state pair.
    do it=1,Num_R
        call advance_out(it)
        call save_prim(it)
        call save_r3(it)
        call save_g34(it)
        call step_scan(it)
        call save_multi(it)
    end do
    call close_events()

    deallocate(y_prev,y_cur,Y,ywait,yevent)

contains

    ! 在初始半径设置 ejecta shell 和 first region-3 状态。
    ! Seed the ejecta shell and first region-3 state at the initial radius.
    ! v3s 只是体积归一化；Y(6) 存无量纲共动体积。
    ! v3s is only a volume normalization; Y(6) stores dimensionless comoving volume.
    subroutine init_shell()
    implicit none

        delta0=Delta_t*para_c
        if (sigma_r <= 0d0) then
            mej=eiso/eta_0/para_c**2
        else
            mej=eiso/(1d0+sigma_r)/eta_0/para_c**2
        end if

        call density_moment(astar,nism,R(1),R0,R_tr,f_jump,f_wide,m2init)
        m2init=4d0*pi*para_m_p*m2init

        m3init=1d1
        R_Gamma(1)=Eta_0-0.001d0
        u2init=dsqrt(R_Gamma(1)*R_Gamma(1)-1d0)
        u4init=dsqrt(Eta_0*Eta_0-1d0)
        dinit=max(delta0,R(1)/Eta_0**2)
        n4init=mej/(4d0*pi*Para_m_p*R(1)*R(1)*Eta_0*dinit)
        g34init=(R_Gamma(1)*R_Gamma(1)+Eta_0*Eta_0-1d0)/(Eta_0*R_Gamma(1)+u2init*u4init)
        compinit=rs_vegas_comp(g34init,sigma_r)
        n3init=compinit*n4init
        v3s=mej/(n3init*Para_m_p)
        Y(1:6)=[R_Gamma(1),R(1),m2init,m3init/mej, &
                rs_mag_internal(g34init,sigma_r)*m3init/mej, &
                m3init/(n3init*Para_m_p)/v3s]
    end subroutine init_shell

    ! 减速尺度只用来设置对数输出时间网格。
    ! The deceleration scale only sets the logarithmic output-time grid.
    subroutine init_dec()
    implicit none

        dmref=eiso/((Eta_0-1d0)*4d0*pi*Para_m_p_E)
        rdecism=(nism*Eta_0/dmref)**(-1d0/3d0)
        if (astar > 0d0) then
            rdecwind=dmref/(2.0d35*astar*Eta_0)
            rdec=min(rdecwind,rdecism)
        else
            rdec=rdecism
        end if
    end subroutine init_dec

    ! 当 Y(4) 到达 1 时，只写一次 crossing 诊断量。
    ! Crossing diagnostics are written exactly once when Y(4) reaches 1.
    subroutine init_cross()
    implicit none

        T_cross=-1d0
        R_cross=0d0
        e3_cross=0d0
        gam20=1d0
        U3_cross=0d0
        V3_cross=0d0
        M3_cross=0d0
        gmcross=0d0
        b3ordx=0d0
    end subroutine init_cross

    ! 输出编号先转成 lab-frame 时间，再推进共享动力学状态。
    ! Convert output index to lab-frame time, then advance the shared dynamical state.
    subroutine advance_out(i_out)
    implicit none
    integer, intent(in) :: i_out

        dtgrid=1d1**(tbin+tspan*(i_out-1d0)/Num_R1)
        tout=tbase+dtgrid
        call advance_to(tout)
    end subroutine advance_out

    ! 写 public ABI 中主 RS 的接触面运动学和扫过质量。
    ! Store public contact kinematics and swept mass for the main RS dynamics.
    subroutine save_prim(i_out)
    implicit none
    integer, intent(in) :: i_out

        R_Tobs(i_out)=tout*(1d0+z)
        R_Gamma(i_out)=Y(1)
        R(i_out)=Y(2)
        M2(i_out)=Y(3)
    end subroutine save_prim

    ! 写 first RS 的 region-3 热状态。multiple-RS 输出绑定 density-jump 分支，
    ! 后面由 save_multi 写。
    ! Store the first-RS region-3 thermal state. Multiple-RS output is branch-based
    ! and is written later by save_multi.
    subroutine save_r3(i_out)
    implicit none
    integer, intent(in) :: i_out

        M3(i_out)=Y(4)*mej
        U3_th(i_out)=Y(5)*mej*para_c**2
        V3_comoving(i_out)=Y(6)*v3s
        B3(i_out)=dB3
    end subroutine save_r3

    ! 写 RS 电子注入使用的瞬时 shock-front gamma34。
    ! Store instantaneous shock-front gamma34 used by RS electron injection.
    subroutine save_g34(i_out)
    implicit none
    integer, intent(in) :: i_out

        Gamma34_inst(i_out)=(Y(1)*Y(1)+Eta_0*Eta_0-1d0)/(Eta_0*Y(1)+dsqrt(Y(1)*Y(1)-1d0)*u4init)
    end subroutine save_g34

    ! crossing 前用 swept ejecta fraction Y(4) 作自变量推进，因此 crossing 事件固定在 Y(4)=1。
    ! Before crossing, integrate with swept ejecta fraction Y(4) as the independent variable,
    ! so the crossing event is pinned at Y(4)=1.
    subroutine advance_m3(ttin)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    integer :: I, J, N, nbracket
    real(8), intent(in) :: ttin
    real(8), parameter :: hmax = 5d-2, rk_eps = 1d-5
    logical :: crossing_first, has_reference
    real(8) :: H_bound, h_event, H_lo, H_hi, H_mid, hest, HH, P, Q
    real(8), dimension(nstate) :: ytry,G,D
    real(8), dimension(rs_nstate) :: rs
    real(8) :: dbtry,ttry,tref

        rs = 0d0
        rs(rs_db3) = dB3
        rhs_phase = precross_phase
        call reverse_dynamics_rhs(rs,tnow,Y,D,nstate,mej,v3s,delta0, &
                                  eta_0,astar,nism,R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e, &
                                  p_f,f_e,e_r,b_r,p_r,fer,sigma_r)
        dB3 = rs(rs_db3)
        rhs_phase = 0
        H_bound = 1d0-Y(4)
        hest = D(4)*(ttin-tnow)
        H_hi = min(hmax, H_bound, hest)
        ttry = tnow
        ytry = Y
        dbtry = dB3
        nbracket = max(2, ceiling(H_hi/hmax))
        HH = H_hi/dble(nbracket)
        do J = 1, nbracket
            call rk_m3(dbtry, ttry, ytry, HH)
        end do
        do while (ttry < ttin .and. H_hi < H_bound)
            H_hi = min(H_bound, 2d0*H_hi)
            ttry = tnow
            ytry = Y
            dbtry = dB3
            nbracket = max(2, ceiling(H_hi/hmax))
            HH = H_hi/dble(nbracket)
            do J = 1, nbracket
                call rk_m3(dbtry, ttry, ytry, HH)
            end do
        end do
        crossing_first = (ttry < ttin)
        if (crossing_first) then
            h_event = H_bound
        else
            H_lo = 0d0
            do I = 1, 30
                H_mid = 0.5d0*(H_lo+H_hi)
                ttry = tnow
                ytry = Y
                dbtry = dB3
                nbracket = max(2, ceiling(H_mid/hmax))
                HH = H_mid/dble(nbracket)
                do J = 1, nbracket
                    call rk_m3(dbtry, ttry, ytry, HH)
                end do
                if (ttry >= ttin) then
                    H_hi = H_mid
                else
                    H_lo = H_mid
                end if
            end do
            h_event = H_hi
        end if

        N = max(1, ceiling(h_event/hmax))
        P = 1d0+rk_eps
        has_reference = .false.
        do while (P >= rk_eps)
            ttry = tnow
            ytry = Y
            dbtry = dB3
            HH = h_event/dble(N)
            do J = 1, N
                call rk_m3(dbtry, ttry, ytry, HH)
            end do
            if (has_reference) then
                P = 0d0
                do I = 1, 3
                    if (.not. ieee_is_finite(ytry(I)) .or. .not. ieee_is_finite(G(I))) then
                        P = huge(1d0)
                        exit
                    end if
                    if (ytry(I)+G(I) <= 0d0) then
                        P = huge(1d0)
                        exit
                    end if
                    Q = 2d0*abs(ytry(I)-G(I))/(ytry(I)+G(I))
                    if (Q > P) P = Q
                end do
                Q = 2d0*abs(ttry-tref)/(ttry+tref)
                if (Q > P) P = Q
            else
                P = 1d0+rk_eps
                has_reference = .true.
            end if
            if (P < rk_eps) then
                dB3 = dbtry
                tnow = ttry
                Y = ytry
            else
                G = ytry
                tref = ttry
                N = N+N
            end if
        end do

        if (crossing_first) then
            ! crossing 后做一次 RHS 计算，记录冻结的 crossing 半径、热状态、gamma_m 和有序磁场。
            ! After crossing, 1 RHS evaluation records the frozen crossing radius,
            ! thermal state, gamma_m, and ordered magnetic field.
            Y(4) = 1d0
            rs = [dB3,T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross, &
                         gmcross,b3ordx]
            rhs_phase = postcross_phase
            call reverse_dynamics_rhs(rs,tnow,Y,D,nstate,mej,v3s,delta0, &
                                      eta_0,astar,nism,R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e, &
                                      p_f,f_e,e_r,b_r,p_r,fer,sigma_r)
            dB3=rs(rs_db3); T_cross=rs(rs_tcross); R_cross=rs(rs_rcross)
            e3_cross=rs(rs_e3cross); gam20=rs(rs_gam20)
            U3_cross=rs(rs_u3cross); V3_cross=rs(rs_v3cross)
            M3_cross=rs(rs_m3cross); gmcross=rs(rs_gammcross)
            b3ordx=rs(rs_b3ordered)
            rhs_phase = 0
        end if

    end subroutine advance_m3

    ! swept-mass fraction 下的一步 RK4。第 4 个状态量由步长直接推进，
    ! 其他状态量使用 reverse_dynamics_rhs 给出的 dY/dM3。
    ! One RK4 step in swept-mass fraction. Component 4 advances by the step itself;
    ! the other components use dY/dM3 from reverse_dynamics_rhs.
    subroutine rk_m3(db_step, t_step, y_step, h_step)
    implicit none
    real(8), dimension(nstate), intent(inout) :: y_step
    real(8), intent(inout) :: db_step,t_step
    real(8), intent(in) :: h_step
    real(8), dimension(nstate) :: Y0,d_step,K1,K2,K3,K4
    real(8), dimension(rs_nstate) :: rs_step
    real(8) :: T0,L1,L2,L3,L4,scale

        Y0 = y_step
        T0 = t_step
        rs_step = 0d0
        rs_step(rs_db3) = db_step
        rhs_phase = precross_phase
        call reverse_dynamics_rhs(rs_step,T0,Y0,d_step,nstate,mej,v3s,delta0, &
                                  eta_0,astar,nism,R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e, &
                                  p_f,f_e,e_r,b_r,p_r,fer,sigma_r)
        db_step = rs_step(rs_db3)
        scale = 1d0/d_step(4)
        K1 = d_step*scale
        K1(4) = 0d0
        L1 = scale

        y_step = Y0+0.5d0*h_step*K1
        y_step(4) = Y0(4)+0.5d0*h_step
        t_step = T0+0.5d0*h_step*L1
        rs_step(rs_db3) = db_step
        call reverse_dynamics_rhs(rs_step,t_step,y_step,d_step,nstate,mej,v3s,delta0, &
                                  eta_0,astar,nism,R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e, &
                                  p_f,f_e,e_r,b_r,p_r,fer,sigma_r)
        db_step = rs_step(rs_db3)
        scale = 1d0/d_step(4)
        K2 = d_step*scale
        K2(4) = 0d0
        L2 = scale

        y_step = Y0+0.5d0*h_step*K2
        y_step(4) = Y0(4)+0.5d0*h_step
        t_step = T0+0.5d0*h_step*L2
        rs_step(rs_db3) = db_step
        call reverse_dynamics_rhs(rs_step,t_step,y_step,d_step,nstate,mej,v3s,delta0, &
                                  eta_0,astar,nism,R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e, &
                                  p_f,f_e,e_r,b_r,p_r,fer,sigma_r)
        db_step = rs_step(rs_db3)
        scale = 1d0/d_step(4)
        K3 = d_step*scale
        K3(4) = 0d0
        L3 = scale

        y_step = Y0+h_step*K3
        y_step(4) = Y0(4)+h_step
        t_step = T0+h_step*L3
        rs_step(rs_db3) = db_step
        call reverse_dynamics_rhs(rs_step,t_step,y_step,d_step,nstate,mej,v3s,delta0, &
                                  eta_0,astar,nism,R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e, &
                                  p_f,f_e,e_r,b_r,p_r,fer,sigma_r)
        db_step = rs_step(rs_db3)
        scale = 1d0/d_step(4)
        K4 = d_step*scale
        K4(4) = 0d0
        L4 = scale

        y_step = Y0+h_step*(K1+2d0*K2+2d0*K3+K4)/6.0d0
        y_step(4) = Y0(4)+h_step
        t_step = T0+h_step*(L1+2d0*L2+2d0*L3+L4)/6.0d0
        rhs_phase = 0
    end subroutine rk_m3

    ! log(time) RK driver，覆盖 waiting、pre-crossing 和 post-crossing 三个阶段。
    ! 如果 pre-crossing 步内碰到 Y(4)=1，就在事件点切步。
    ! Log-time RK driver for waiting, pre-crossing, and post-crossing phases.
    ! A pre-crossing step that reaches Y(4)=1 is split at the event.
    subroutine advance_log(db_step, tcs, rcs, e3_hot, gam20_step, &
                                        u3s, v3sstep, m3sstep, &
                                        gms, b3ords, &
                                        t_step, h_step, y_step, phase)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    integer, intent(in) :: phase
    integer :: I, J, N
    real(8), intent(inout) :: db_step, tcs, rcs, e3_hot, gam20_step
    real(8), intent(inout) :: u3s, v3sstep, m3sstep, gms, b3ords
    real(8), dimension(nstate), intent(inout) :: y_step
    real(8), intent(inout) :: t_step
    real(8), intent(in) :: h_step
    logical :: has_reference
    real(8), parameter :: rk_eps = 1d-5, mtol = 1d-7, ttol = 1d-12
    real(8), dimension(nstate) :: C,G,D,y_trial
    real(8), dimension(rs_nstate) :: rs_try,rs_trial
    real(8) :: HH,P,Q,X,S,HS,s_trial,H_lo,H_hi,H_mid,h_event,h_post,t_phys

        N = 1
        HS = log((t_step+h_step)/t_step)
        HH = HS
        P = 1d0+rk_eps
        X = t_step
        C = y_step
        has_reference = .false.
        do while (P >= rk_eps)
            y_step = C
            S = 0d0
            rs_try=[db_step,tcs,rcs,e3_hot,gam20_step, &
                        u3s,v3sstep,m3sstep,gms,b3ords]
            do J = 1, N
                select case (phase)
                case (wait_phase)
                call rk_log(wait_phase, rs_try, X, S, HH, y_step)
                case (precross_phase)
                    if (rs_try(rs_tcross) >= 0d0 .or. y_step(4) >= 1d0) then
                        call rk_log(postcross_phase, rs_try, X, S, HH, y_step)
                    else
                        y_trial = y_step
                        s_trial = S
                        rs_trial = rs_try
                        call rk_log(precross_phase, rs_trial, X, s_trial, HH, y_trial)
                        if (y_trial(4) < 1d0) then
                            rs_try = rs_trial
                            S = s_trial
                            y_step = y_trial
                        else
                            H_lo = 0d0
                            H_hi = HH
                            do I = 1, 60
                                H_mid = 0.5d0*(H_lo+H_hi)
                                y_trial = y_step
                                s_trial = S
                                rs_trial = rs_try
                                call rk_log(precross_phase, rs_trial, &
                                                                   X, s_trial, H_mid, y_trial)
                                if (y_trial(4) >= 1d0) then
                                    H_hi = H_mid
                                else
                                    H_lo = H_mid
                                end if
                                if (abs(y_trial(4)-1d0) <= mtol) exit
                                if (H_hi-H_lo <= ttol*abs(HH)) exit
                            end do
                            h_event = H_hi
                            call rk_log(precross_phase, rs_try, X, S, h_event, y_step)
                            y_step(4) = 1d0
                            rhs_phase = postcross_phase
                            t_phys = X*exp(S)
                            call reverse_dynamics_rhs(rs_try,t_phys,y_step,D,nstate,mej,v3s,delta0, &
                                                      eta_0,astar,nism,R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e, &
                                                      p_f,f_e,e_r,b_r,p_r,fer,sigma_r)
                            rhs_phase = 0
                            h_post = HH-h_event
                            if (h_post > 0d0) then
                                call rk_log(postcross_phase, rs_try, X, S, h_post, y_step)
                            end if
                        end if
                    end if
                case (postcross_phase)
                    call rk_log(postcross_phase, rs_try, X, S, HH, y_step)
                end select
            end do
            if (has_reference) then
                P = 0d0
                do I = 1, 3
                    if (.not. ieee_is_finite(y_step(I)) .or. .not. ieee_is_finite(G(I))) then
                        P = huge(1d0)
                        exit
                    end if
                    if (y_step(I)+G(I) <= 0d0) then
                        P = huge(1d0)
                        exit
                    end if
                    Q = 2d0*abs(y_step(I)-G(I))/(y_step(I)+G(I))
                    if (Q > P) P = Q
                end do
            else
                P = 1d0+rk_eps
                has_reference = .true.
            end if
            if (P < rk_eps) then
                db_step=rs_try(rs_db3); tcs=rs_try(rs_tcross); rcs=rs_try(rs_rcross)
                e3_hot=rs_try(rs_e3cross); gam20_step=rs_try(rs_gam20)
                u3s=rs_try(rs_u3cross); v3sstep=rs_try(rs_v3cross)
                m3sstep=rs_try(rs_m3cross); gms=rs_try(rs_gammcross)
                b3ords=rs_try(rs_b3ordered)
            else
                G = y_step
                HH = 0.5d0*HH
                N = N+N
            end if
        end do
        t_step = X

    end subroutine advance_log

    ! log(time) 下的 RK4 stage。rhs_phase 选择物理分支，
    ! 状态向量仍使用同一个 Y 布局。
    ! RK4 stage in log(time). rhs_phase selects the physical branch;
    ! the state vector keeps the same Y layout.
    subroutine rk_log(phase, rs_step, x_base, s_step, h_phase, y_phase)
    implicit none
    integer, intent(in) :: phase
    integer :: K
    real(8), dimension(rs_nstate), intent(inout) :: rs_step
    real(8), dimension(nstate), intent(inout) :: y_phase
    real(8), intent(inout) :: s_step
    real(8), intent(in) :: x_base, h_phase
    real(8), dimension(nstate) :: d_phase,b_phase,e_phase
    real(8), dimension(4) :: a_phase
    real(8) :: s_stage,t_phys

        a_phase = [0.5d0*h_phase, 0.5d0*h_phase, h_phase, h_phase]
        rhs_phase = phase
        t_phys = x_base*exp(s_step)
        call reverse_dynamics_rhs(rs_step,t_phys,y_phase,d_phase,nstate,mej,v3s,delta0, &
                                  eta_0,astar,nism,R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e, &
                                  p_f,f_e,e_r,b_r,p_r,fer,sigma_r)
        d_phase = t_phys*d_phase
        e_phase = y_phase
        b_phase = y_phase
        do K = 1, 3
            y_phase = e_phase+a_phase(K)*d_phase
            b_phase = b_phase+a_phase(K+1)*d_phase/3.0d0
            s_stage = s_step+a_phase(K)
            t_phys = x_base*exp(s_stage)
            call reverse_dynamics_rhs(rs_step,t_phys,y_phase,d_phase,nstate,mej,v3s,delta0, &
                                      eta_0,astar,nism,R_tr,f_jump,f_wide,R0,Epsilon_b,Epsilon_e, &
                                      p_f,f_e,e_r,b_r,p_r,fer,sigma_r)
            d_phase = t_phys*d_phase
        end do
        y_phase = b_phase+h_phase*d_phase/6.0d0
        s_step = s_step+h_phase
        rhs_phase = 0
    end subroutine rk_log

    ! 推进到目标 lab-frame 时间。磁化 ejecta 在 pressure 和 fast-mode 条件满足前
    ! 可能先处于 waiting phase。
    ! Advance to a lab-frame target time. Magnetized ejecta may first stay in a
    ! waiting phase before pressure and fast-mode conditions allow an RS.
    subroutine advance_to(ttin)
    implicit none
    real(8), intent(in) :: ttin
    real(8) :: t_local,h_local

        do while (tnow < ttin)
            if (T_cross < 0d0 .and. Y(4) < 1d0) then
                if (pressure_ok(Y)) then
                    call advance_m3(ttin)
                else
                    call wait_trial(ttin,ywait,dbwait)
                    ready=pressure_ok(ywait)
                    if (.not. ready) then
                        Y=ywait
                        dB3=dbwait
                        tnow=ttin
                    else
                        call wait_root(ttin,tevent,yevent,dbevent)
                        Y=yevent
                        dB3=dbevent
                        tnow=tevent
                    end if
                end if
            else
                h_local=ttin-tnow
                t_local=tnow
                call advance_log(dB3,T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross, &
                                             gmcross,b3ordx,t_local,h_local,Y,precross_phase)
                tnow=ttin
            end if
        end do
    end subroutine advance_to

    ! waiting phase 的 trial 推进，用于二分定位；不提交 host Y。
    ! Waiting-phase trial used for bisection; it does not commit the host Y.
    subroutine wait_trial(end_t,y_out,db_out)
    implicit none
    real(8), intent(in) :: end_t
    real(8), dimension(nstate), intent(out) :: y_out
    real(8), intent(out) :: db_out
    real(8) :: t_local,h_local,tctry,rctry,e3try,gam20_trial
    real(8) :: u3try,v3try,m3try,gmtry,b3ordtry

        y_out=Y
        db_out=dB3
        tctry=T_cross; rctry=R_cross; e3try=e3_cross; gam20_trial=gam20
        u3try=U3_cross; v3try=V3_cross; m3try=M3_cross
        gmtry=gmcross; b3ordtry=b3ordx
        t_local=tnow
        h_local=end_t-tnow
        call advance_log(db_out,tctry,rctry,e3try,gam20_trial, &
                                     u3try,v3try,m3try,gmtry, &
                                     b3ordtry,t_local,h_local,y_out,wait_phase)
    end subroutine wait_trial

    ! 定位 waiting phase 到 shock-ready phase 的过渡时间。
    ! Locate the transition from waiting phase to shock-ready phase.
    subroutine wait_root(end_t,event_t,event_y,event_db)
    implicit none
    integer :: iter
    real(8), intent(in) :: end_t
    real(8), dimension(nstate), intent(out) :: event_y
    real(8), intent(out) :: event_t,event_db
    real(8) :: t_lo,t_hi,t_mid

        t_lo=tnow
        t_hi=end_t
        event_y=ywait
        event_db=dbwait
        do iter=1,32
            t_mid=0.5d0*(t_lo+t_hi)
            call wait_trial(t_mid,ywait,dbwait)
            if (pressure_ok(ywait)) then
                t_hi=t_mid
                event_y=ywait
                event_db=dbwait
            else
                t_lo=t_mid
            end if
        end do
        event_t=t_hi
    end subroutine wait_root

    ! 磁化 RS 形成条件：pressure balance，加上 shock frame 中上游速度超过 fast magnetosonic 条件。
    ! Magnetized RS formation condition: pressure balance plus upstream speed
    ! above the fast magnetosonic condition in the shock frame.
    logical function pressure_ok(state)
    implicit none
    real(8), dimension(nstate), intent(in) :: state
    real(8) :: radius,gamma_cd,u_cd,u4,gamma43,shell_d,n1,n4
    real(8) :: u_down,u_up,sigma_cr
    logical :: p_ready,fast_ok

        if (sigma_r <= 0d0) then
            pressure_ok=.true.
            return
        end if
        radius=state(2)
        gamma_cd=state(1)
        u_cd=dsqrt(gamma_cd*gamma_cd-1d0)
        u4=dsqrt(Eta_0*Eta_0-1d0)
        gamma43=(gamma_cd*gamma_cd+Eta_0*Eta_0-1d0)/(Eta_0*gamma_cd+u_cd*u4)
        shell_d=max(delta0,radius/(Eta_0*Eta_0))
        call density_profile(astar,nism,radius,R0,1,R_tr,f_jump,f_wide,n1)
        n4=mej/(4d0*pi*Para_m_p*radius*radius*Eta_0*shell_d)
        sigma_cr=2d0*(4d0*gamma_cd+3d0)*(gamma_cd-1d0)*n1/(3d0*n4)
        p_ready=(sigma_cr >= sigma_r)
        u_down=rs_vegas_ud(gamma43,sigma_r)
        u_up=dsqrt((1d0+u_down*u_down)*(gamma43-1d0)*(gamma43+1d0))+u_down*gamma43
        fast_ok=(u_up*u_up > sigma_r)
        pressure_ok=p_ready .and. fast_ok
    end function pressure_ok


    ! multiple-RS 输出和 first-RS 输出分开初始化。只有 density-jump source 发生耗散后，
    ! 对应分支才会变成非零。
    ! Multiple-RS outputs are initialized separately from first-RS outputs.
    ! A branch becomes non-0d0 only after its density-jump source dissipates.
    subroutine init_output()
    implicit none

        call init_branch()
        call init_total()
        call init_diag()
        call init_history()
    end subroutine init_output

    subroutine init_branch()
    implicit none

        branch_m3=0d0
        branch_u3=0d0
        branch_v3=0d0
        branch_b3=0d0
        branch_gm=0d0
        branch_gc=0d0
        branch_g43=1d0
        branch_comp=1d0
        branch_brs=0d0
        branch_ud=0d0
    end subroutine init_branch

    subroutine init_total()
    implicit none

        total_m3=0d0
        total_u3=0d0
        total_v3=0d0
        total_b3=0d0
        total_p=0d0
        total_w=0d0
    end subroutine init_total

    subroutine init_diag()
    implicit none

        avg_gc=0d0
        avg_p3=0d0
        avg_g43=1d0
        avg_brs=0d0
        total_ud=0d0
        diss_e=0d0
        inj_e=0d0
        nu_m=0d0
        nu_c=0d0
    end subroutine init_diag

    subroutine init_history()
    implicit none

        diss_prev=0d0
        gm_prev=0d0
    end subroutine init_history

    subroutine init_event()
    implicit none

        event_on=.false.
        event_done=.false.
        start_r=0d0
        end_r=0d0
        start_t=0d0
        end_t=0d0
    end subroutine init_event

    ! 第一次输出推进前，初始化事件扫描的左端点。
    ! Initialize the left endpoint for event scanning before the first output advance.
    subroutine init_scan()
    implicit none

        prev_r=Y(2)
        prev_g=Y(1)
        prev_t=tbase*(1d0+z)
        y_prev=Y
        src_prev=-1d0
        do event_j=1,njump
            call event_source(event_j,prev_r,prev_g, &
                                                y_prev,src_prev(event_j))
        end do
    end subroutine init_scan

    ! first-RS 推进后，构造事件扫描的当前右端点。
    ! Build the current right endpoint for event scanning after first-RS advancement.
    subroutine step_scan(i_out)
    implicit none
    integer, intent(in) :: i_out

        curr_r=Y(2)
        curr_g=Y(1)
        curr_t=R_Tobs(i_out)
        y_cur=Y
        src_cur=-1d0
        do event_j=1,njump
            call event_source(event_j,curr_r,curr_g, &
                                                y_cur,src_cur(event_j))
        end do
        call scan_events()
        prev_r=curr_r
        prev_g=curr_g
        prev_t=curr_t
        src_prev=src_cur
        y_prev=y_cur
    end subroutine step_scan

    ! 如果分支到采样网格末端仍未关闭，就在最后一个已知物理状态关闭，
    ! 不外推编造 endpoint。
    ! If a branch stays open through the sampled grid, close it at the last known
    ! physical state rather than inventing an extrapolated endpoint.
    subroutine close_events()
    implicit none

        do event_j=1,njump
            if (event_on(event_j) .and. .not. event_done(event_j)) then
                event_done(event_j)=.true.
                end_r(event_j)=prev_r
                end_t(event_j)=prev_t
            end if
        end do
    end subroutine close_events

    ! 在前后两个 first-RS 状态之间扫描每个 density-jump window。
    ! source > 0 表示 multiple-RS 分支物理上处于 active。
    ! Scan each density-jump window between previous and current first-RS states.
    ! source > 0 means a multiple-RS branch is physically active.
    subroutine scan_events()
    implicit none
    integer :: j,k,nscan
    real(8) :: rl,rh,gl,gh,tl,th,sl,sh
    real(8) :: width,ol,oh,drs
    real(8) :: r_root,t_root
    real(8), dimension(nstate) :: yl,yh
    real(8) :: fs

        do j=1,njump
            width=jump_w(j)*jump_r(j)
            ol=max(prev_r,jump_r(j)-4d0*width)
            oh=min(curr_r,jump_r(j))
            drs=max(0d0,oh-ol)
            nscan=max(1,ceiling(drs/(width/16d0)))
            rl=prev_r; gl=prev_g
            tl=prev_t; sl=src_prev(j)
            yl=y_prev
            do k=1,nscan
                fs=dble(k)/dble(nscan)
                rh=prev_r+(curr_r-prev_r)*fs
                gh=prev_g+(curr_g-prev_g)*fs
                th=prev_t+(curr_t-prev_t)*fs
                yh=y_prev+fs*(y_cur-y_prev)
                call event_source(j,rh,gh,yh,sh)
                if (.not. event_on(j)) then
                    if (sl > 0d0) then
                        event_on(j)=.true.
                        start_r(j)=rl
                        start_t(j)=tl
                    else if (sl <= 0d0 .and. sh > 0d0) then
                        call event_root(j,rl,rh,gl,gh,tl,th,yl,yh,r_root,t_root)
                        event_on(j)=.true.
                        start_r(j)=r_root
                        start_t(j)=t_root
                    end if
                end if
                if (event_on(j) .and. .not. event_done(j)) then
                    if (sl > 0d0 .and. sh <= 0d0) then
                        call event_root(j,rl,rh,gl,gh,tl,th,yl,yh,r_root,t_root)
                        event_done(j)=.true.
                        if (r_root >= jump_r(j)) r_root=nearest(jump_r(j),-1d0)
                        end_r(j)=r_root
                        end_t(j)=t_root
                    end if
                end if
                rl=rh; gl=gh; tl=th; sl=sh; yl=yh
            end do
        end do
    end subroutine scan_events

    ! 对 source 变号点做半径二分；这个点对应 multiple-RS 分支开始或结束。
    ! Bisect in radius for the source sign change that starts or ends a multiple-RS branch.
    subroutine event_root(j,rli,rhi,gli,ghi,tli,t_hi,yli,yhi,r0,t0root)
    implicit none
    integer, intent(in) :: j
    integer :: iter
    real(8), intent(in) :: rli,rhi,gli,ghi,tli,t_hi
    real(8), dimension(nstate), intent(in) :: yli,yhi
    real(8), intent(out) :: r0,t0root
    real(8) :: rl,rh,rm,gl,gh,gm,tl,th,tm,sl,sm,frac
    real(8), dimension(nstate) :: ymid

        rl=rli; rh=rhi
        gl=gli; gh=ghi
        tl=tli; th=t_hi
        call event_source(j,rl,gl,yli,sl)
        do iter=1,80
            rm=0.5d0*(rl+rh)
            frac=(rm-rli)/(rhi-rli)
            gm=gli+frac*(ghi-gli)
            tm=tli+frac*(t_hi-tli)
            ymid=yli+frac*(yhi-yli)
            call event_source(j,rm,gm,ymid,sm)
            if (sl*sm <= 0d0) then
                rh=rm; gh=gm; th=tm
            else
                rl=rm; gl=gm; tl=tm; sl=sm
            end if
        end do
        r0=0.5d0*(rl+rh)
        t0root=0.5d0*(tl+th)
    end subroutine event_root

    ! multiple-RS mechanical source：新 shock 下游热能减去上游分支绝热压缩后的能量。
    ! Multiple-RS mechanical source: newly shocked downstream thermal energy
    ! minus adiabatically compressed upstream branch energy.
    subroutine event_source(j,radius,gb,state,source)
    implicit none
    integer, intent(in) :: j
    integer :: parent_m,parent_u,parent_v
    real(8), intent(in) :: radius,gb
    real(8), dimension(nstate), intent(in) :: state
    real(8), intent(out) :: source
    real(8) :: dens_all,dens_bump,n1,dens_ex,dens_pre,n4,e4,p4,p3,e3_hot,e_ad
    real(8) :: contact_g,shock_g,comp,beta_rs

        call branch_density(radius,j,dens_all,dens_bump)
        if (dens_bump <= 0d0 .or. gb <= 1d0) then
            source=-1d0
            return
        end if
        n1=dens_all
        dens_ex=dens_bump
        dens_pre=n1-dens_ex
        if (dens_pre <= 0d0) error stop 'secondary reverse event source found non-positive pre-bump density'
        n4=4d0*gb*dens_pre
        e4=4d0*gb*(gb-1d0)*dens_pre*Para_m_p*para_c**2
        p4=e4/3d0
        if (j > 1) then
            parent_m=6+j-1
            parent_u=6+njump+j-1
            parent_v=6+2*njump+j-1
            if (parent_ready(j,state(parent_m),state(parent_u),state(parent_v))) then
                n4=state(parent_m)*mej/(Para_m_p*state(parent_v)*v3s)
                e4=state(parent_u)*mej*para_c**2/(state(parent_v)*v3s)
                p4=e4/3d0
            end if
        end if
        call reverse_contact(gb,n1,n4,e4,p4,contact_g,p3,shock_g,comp,beta_rs)
        if (comp <= 0d0) then
            source=-1d0
            return
        end if
        e3_hot=3d0*p3
        e_ad=e4*comp**(4d0/3d0)
        source=e3_hot-e_ad
    end subroutine event_source

    ! 把分支状态变量写入 public arrays，并累计 total diagnostics。
    ! Convert branch state variables into public arrays and total diagnostics.
    subroutine save_multi(i_out)
    implicit none
    integer, intent(in) :: i_out
    integer :: j
    real(8) :: de,dg,wd,wg
    real(8) :: mt,ut,vt,bt

        wd=0d0; wg=0d0
        do j=1,njump
            call save_branch(i_out,j,de,dg)
            call sum_diag(i_out,j,de,dg,wd,wg)
            call sum_total(i_out,j)
        end do

        if (wd > 0d0) then
            avg_gc(i_out)=avg_gc(i_out)/wd
            avg_p3(i_out)=avg_p3(i_out)/wd
            avg_g43(i_out)=1d0+avg_g43(i_out)/wd
            avg_brs(i_out)=avg_brs(i_out)/wd
        end if

        mt=total_m3(i_out)
        ut=total_u3(i_out)
        vt=total_v3(i_out)
        if (vt <= 0d0) return

        bt=dsqrt(8d0*pi*b_r*ut/vt)
        total_b3(i_out)=bt
        total_p(i_out)=ut/(3d0*vt)
        total_w(i_out)=mt*para_c**2/vt+ut/vt+total_p(i_out)
        call freqs(i_out,wd,wg,bt)
    end subroutine save_multi

    ! 写一个 multiple-RS 分支，并返回本输出步新增的耗散能。
    ! Store 1 multiple-RS branch and return the newly dissipated energy.
    subroutine save_branch(i_out,jb,de,dg)
    implicit none
    integer, intent(in) :: i_out,jb
    integer :: im,iu,iv,ie,ig
    real(8), intent(out) :: de,dg
    real(8) :: ecum,gcum

        im=6+jb
        iu=6+njump+jb
        iv=6+2*njump+jb
        ie=6+3*njump+jb
        ig=6+4*njump+jb

        branch_m3(jb,i_out)=Y(im)*mej
        branch_u3(jb,i_out)=Y(iu)*mej*para_c**2
        branch_v3(jb,i_out)=Y(iv)*v3s
        if (branch_v3(jb,i_out) > 0d0) &
            branch_b3(jb,i_out)=dsqrt(8d0*pi*b_r*branch_u3(jb,i_out)/branch_v3(jb,i_out))

        ecum=Y(ie)*mej*para_c**2
        gcum=Y(ig)*mej*para_c**2
        de=ecum-diss_prev(jb)
        dg=gcum-gm_prev(jb)
        if (de < 0d0) error stop 'multiple RS dissipated energy must not decrease'

        diss_e(i_out)=diss_e(i_out)+de
        inj_e(i_out)=inj_e(i_out)+e_r*de
        if (de > 0d0) branch_gm(jb,i_out)=dg/de

        diss_prev(jb)=ecum
        gm_prev(jb)=gcum
    end subroutine save_branch

    ! 总质量、总热能和总体积是各分支直接求和。
    ! Total mass, thermal energy, and volume are direct sums over branches.
    subroutine sum_total(i_out,jb)
    implicit none
    integer, intent(in) :: i_out,jb

        total_m3(i_out)=total_m3(i_out)+branch_m3(jb,i_out)
        total_u3(i_out)=total_u3(i_out)+branch_u3(jb,i_out)
        total_v3(i_out)=total_v3(i_out)+branch_v3(jb,i_out)
    end subroutine sum_total

    ! shock diagnostics 按新增耗散能加权，而不是按总分支质量加权；
    ! 否则旧的 inactive 物质会主导平均值。
    ! Shock diagnostics are weighted by newly dissipated energy, not by total
    ! branch mass; otherwise inactive old material would dominate averages.
    subroutine sum_diag(i_out,jb,de,dg,wd,wg)
    implicit none
    integer, intent(in) :: i_out,jb
    real(8), intent(in) :: de,dg
    real(8), intent(inout) :: wd,wg
    logical :: hasup
    real(8) :: n1,n4,e4,p4,p3,e3_hot,e_ad,comp,beta_rs
    real(8) :: contact_g,shock_g,n3,gmj,b_i,gemax

        if (de <= 0d0) return
        call get_upstream(R(i_out),i_out,jb,hasup,n1,n4,e4,p4)
        if (.not. hasup) return

        call reverse_contact(R_Gamma(i_out),n1,n4,e4,p4,contact_g,p3,shock_g,comp,beta_rs)
        if (comp <= 0d0) return

        e3_hot=3d0*p3
        e_ad=e4*comp**(4d0/3d0)
        n3=comp*n4
        if (p_r <= 2d0) error stop 'secondary RS injection requires p_r > 2'
        b_i=dsqrt(8d0*pi*b_r*e3_hot)
        gmj=1d0+e_r/fer*(p_r-2d0)/(p_r-1d0)*(e3_hot-e_ad)/(n3*Para_m_e*para_c**2)
        gemax=3d0*Para_m_energy/dsqrt(8d0*b_i*Para_e**3)
        if (gmj >= gemax) error stop 'secondary RS electron injection exceeds maximum'

        call sum_weight(i_out,jb,de,dg,contact_g,p3,shock_g, &
                                                       comp,beta_rs,e3_hot,e_ad,wd,wg)
    end subroutine sum_diag

    ! 第 1 个分支的上游是 bump 前的 FS 下游物质。更后面的分支在 parent 有有限状态后，
    ! 可以使用前一个 multiple-RS 分支作为上游。
    ! Branch 1 upstream is the pre-bump FS downstream material. Later branches
    ! can use the previous multiple-RS branch after it has finite state.
    subroutine get_upstream(radius,i_out,jb,hasup,n1,n4,e4,p4)
    implicit none
    integer, intent(in) :: i_out,jb
    logical, intent(out) :: hasup
    real(8), intent(in) :: radius
    real(8), intent(out) :: n1,n4,e4,p4
    integer :: ipar
    real(8) :: dens_all,dens_bump,dens_ex,dens_pre

        hasup=.false.
        call branch_density(radius,jb,dens_all,dens_bump)
        if (dens_bump <= 0d0 .or. R_Gamma(i_out) <= 1d0) return

        n1=dens_all
        dens_ex=dens_bump
        dens_pre=n1-dens_ex
        if (dens_pre <= 0d0) error stop 'secondary branch output found non-positive pre-bump density'
        n4=4d0*R_Gamma(i_out)*dens_pre
        e4=4d0*R_Gamma(i_out)*(R_Gamma(i_out)-1d0)*dens_pre*Para_m_p*para_c**2
        p4=e4/3d0
        if (jb > 1) then
            ipar = jb-1
            if (parent_ready(jb,branch_m3(ipar,i_out)/mej, &
                                                    branch_u3(ipar,i_out)/(mej*para_c**2), &
                                                    branch_v3(ipar,i_out)/v3s)) then
                n4=branch_m3(ipar,i_out)/(Para_m_p*branch_v3(ipar,i_out))
                e4=branch_u3(ipar,i_out)/branch_v3(ipar,i_out)
                p4=e4/3d0
            end if
        end if

        hasup=.true.
    end subroutine get_upstream

    ! 累计当前输出时刻的能量加权分支诊断量。
    ! Accumulate energy-weighted branch diagnostics for this output time.
    subroutine sum_weight(i_out,jb,de,dg,contact_g,p3,shock_g,comp,beta_rs,e3_hot,e_ad,wd,wg)
    implicit none
    integer, intent(in) :: i_out,jb
    real(8), intent(in) :: de,dg,contact_g,p3,shock_g,comp,beta_rs,e3_hot,e_ad
    real(8), intent(inout) :: wd,wg

        wd=wd+de
        wg=wg+dg
        avg_gc(i_out)=avg_gc(i_out)+de*contact_g
        avg_p3(i_out)=avg_p3(i_out)+de*p3
        avg_g43(i_out)=avg_g43(i_out)+de*(shock_g-1d0)
        avg_brs(i_out)=avg_brs(i_out)+de*beta_rs
        total_ud(i_out)=total_ud(i_out)+(e3_hot-e_ad)
        branch_gc(jb,i_out)=contact_g
        branch_g43(jb,i_out)=shock_g
        branch_comp(jb,i_out)=comp
        branch_brs(jb,i_out)=beta_rs
        branch_ud(jb,i_out)=e3_hot-e_ad
    end subroutine sum_weight

    ! 计算合并 multiple-RS 成分的 synchrotron break frequencies。
    ! Compute synchrotron break frequencies for the combined multiple-RS component.
    subroutine freqs(i_out,wd,wg,b)
    implicit none
    integer, intent(in) :: i_out
    real(8), intent(in) :: wd,wg,b
    real(8) :: dop,gcool

        if (wd > 0d0 .and. avg_gc(i_out) > 1d0) then
            dop=avg_gc(i_out)*(1d0-dsqrt(1d0-avg_gc(i_out)**(-2)))*(1d0+z)
            nu_m(i_out)=4.2d6*b*(wg/wd)**2/dop
        end if
        gcool=7.7d8*(1d0+z)/(R_Gamma(i_out)*b*b*R_Tobs(i_out))
        dop=R_Gamma(i_out)*(1d0-dsqrt(1d0-R_Gamma(i_out)**(-2)))*(1d0+z)
        nu_c(i_out)=4.2d6*b*gcool*gcool/dop
    end subroutine freqs

    ! 把当前 density profile 拆成 smooth baseline 加一个选中的 Gaussian enhancement。
    ! 每个 enhancement 最多对应一个分支。
    ! Split the active density profile into smooth baseline plus 1 selected
    ! Gaussian enhancement. Each enhancement owns at most 1 branch.
    subroutine branch_density(radius,j,dens_all,dens_bump)
    implicit none
    integer, intent(in) :: j
    integer :: k
    real(8), intent(in) :: radius
    real(8), intent(out) :: dens_all,dens_bump
    real(8) :: x,width,prof,nall,nbase,enh

        call density_profile(astar,nism,radius,R0,1,R_tr,f_jump,f_wide,nall)
        enh=1d0; dens_bump=0d0
        do k=1,njump
            x=radius-jump_r(k)
            width=jump_w(k)*jump_r(k)
            prof=(jump_f(k)-1d0)*dexp(-(x*x)/(2d0*width*width))
            enh=enh+prof
        end do
        nbase=nall/enh
        do k=1,njump
            x=radius-jump_r(k)
            width=jump_w(k)*jump_r(k)
            prof=(jump_f(k)-1d0)*dexp(-(x*x)/(2d0*width*width))
            if (k == j .and. x >= -4d0*width .and. x < 0d0) dens_bump=nbase*prof
        end do
        dens_all=nall
    end subroutine branch_density

    ! 判断 chained multiple RS 中 parent 分支是否已经可用。
    ! Parent branch eligibility for chained multiple RS.
    logical function parent_ready(j,pm,pe,pv)
    implicit none
    integer, intent(in) :: j
    real(8), intent(in) :: pm,pe,pv
        parent_ready=.false.
        if (j <= 1) return
        if (pm <= 0d0 .or. pe <= 0d0 .or. pv <= 0d0) return
        parent_ready=.true.
    end function parent_ready


    ! Boundary 解包留在本地，保证 public f2py 参数列表和顺序不变，
    ! 同时内部实现可以使用具名标量。
    ! Boundary unpacking stays local so the public f2py argument list and order
    ! remain unchanged while the implementation uses named scalars.
    subroutine unpack()
    implicit none

        Eta_0=Boundary(1)
        R(1)=Boundary(4)
        Epsilon_e=Boundary(5)
        Epsilon_b=Boundary(6)
        p_f=Boundary(7)
        z=Boundary(8)
        nism=Boundary(11)
        astar=Boundary(12)
        eiso=Boundary(14)
        tdur=Boundary(15)
        f_e=Boundary(16)
        R_tr=Boundary(21)
        f_jump=Boundary(22)
        f_wide=Boundary(23)
        if (n >= ir0) then
            R0=Boundary(ir0)
        else
            R0=Boundary(n)
        end if
    end subroutine unpack

end subroutine dynamics_reverse
