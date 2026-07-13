module electron_reverse_kernel
    use constants
    use dynamics_density_profile, only: density_profile
    use electron_coord_common, only: build_fourvel_grid, coord_from_xg, &
                                                 dxg_dcoord, &
                                                 coord_fourvel, &
                                                 fourvel_scale
    use electron_injection_profiles, only: kinetic_edges, &
                                           kinetic_coord, &
                                           log_edges
    use electron_common, only: tail_factor
    use electron_shell_transport, only: resolve_solver, &
                                               shell_coord_step, &
                                               coord_to_dgamma, &
                                               solver_dg, solver_fullhide
    use electron_transport_common, only: build_piece_u, remap_edges, &
                                         dnx_dgamma, trace_piece_u, u_edges, &
                                         flux_split_step, x_from_u
    use electron_dg_transport, only: dg_mesh, dg_char_step, &
                                               dg_advance_step, &
                                               dg_build_mesh, dg_integral, &
                                               dg_limit_positive, &
                                               dg_kinetic_source, &
                                               dg_project_state, dg_project_cells, &
                                               dg_scale_content, dg_tail_fraction
    use electron_radiation_kernel, only: syn_state, nua_fromtau
    use electron_cooling_kernel, only: forward_cooling
    implicit none
    integer, parameter :: reverse_dg_base_substeps = 10
contains

! 反向激波电子演化主驱动：注入 -> 同步/IC 冷却 -> 1D 或 DG 输运推进。
! Reverse-shock electron driver: injection -> synchrotron/IC cooling -> 1D or DG transport.
    subroutine electron_reverse_evolve(Delta_0,e_r,b_r,p_r,f_e_r,eta_0,Epsilon_e,Epsilon_b,z,A_star,dNe_ISM,para_m_ej, &
                                       R_tr,f_jump,f_wide,R0, &
                                       T_cross,R_cross,U3_cross,V3_cross,M3_cross,R_Tobs,R_Gamma,R,B3,M3_shell, &
                                       U3_shell,V3_shell,V_seed, &
                                       Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads,gam_e,dN_gam_e, &
                                       solver_id)
    implicit none
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads
    integer, intent(in), optional :: solver_id
    integer :: I_tobs,I_gam_e,L1,L,active_solver,i_empty,i_edge
    real(8), intent(in) :: Delta_0,e_r,b_r,p_r,f_e_r,eta_0,Epsilon_e,Epsilon_b,z,A_star,dNe_ISM,para_m_ej
    real(8), intent(in) :: R_tr,f_jump,f_wide,R0
    real(8), intent(in) :: T_cross,R_cross,U3_cross,V3_cross,M3_cross
    real(8), intent(in), dimension(Num_R) :: R_Tobs, R_Gamma, R, B3, M3_shell, U3_shell, V3_shell
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(out), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_gam_e,Num_R) :: dN_gam_e
    real(8), parameter :: gc_coeff=7.7d8, bsyn_coeff=0.39d0, adv_coeff=1.35d-19, trans_break=2d0
    real(8), parameter :: tail_thresh=1d-10, tail_power=2d0
    real(8) :: factor2, dB, gamma34, gmax, gm, gc, dNe, DB_min, gmax0, gmin, d_x, rloc, gloc, Delta, rn4
    real(8) :: beta4, beta2, u2, u4, f_r, dDR, dDD, Qshell, coolscale, R_step_mid, vol_lo, vol_hi
    real(8) :: thermloss, adrate, dgscale, coord_scale, dglow, dgmid, dghigh
    real(8) :: injection_rate, inj_hi, inj_width, mass_lo, mass_hi, rhi
    real(8) :: pc_scale,pc_shift,step_scale,step_shift
    real(8), allocatable, dimension(:) :: dEl,x,dF1,temp3,dN_x,x_edge,coord_edge
    real(8), allocatable, dimension(:) :: pc_log,pc_map,pc_work,pc_back,pc_u,pc_a,pc_b
    real(8), allocatable, dimension(:,:) :: pc_affine
    real(8), allocatable, dimension(:) :: dB3_serial,P_syn,Seed_syn,cooling_aux,Compton
    type(dg_mesh) :: dg_mesh,dg_new_mesh
    real(8), allocatable, dimension(:) :: dg_state,dg_work,dg_dEl,dg_source
    logical :: crossed_reverse,pc_ready

    allocate(dEl(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e),temp3(Num_gam_e-1),dN_x(Num_gam_e), &
             x_edge(Num_gam_e+1),coord_edge(Num_gam_e+1), &
             pc_log(Num_gam_e),pc_map(Num_gam_e+1), &
             pc_work(Num_gam_e+1),pc_back(Num_gam_e+1), &
             pc_u(Num_gam_e+1),pc_a(Num_gam_e),pc_b(Num_gam_e), &
             pc_affine(Num_gam_e+1,1), &
             dB3_serial(Num_R),P_syn(Num_nu),Seed_syn(Num_nu), &
             cooling_aux(Num_gam_e),Compton(Num_gam_e))
    dB3_serial=B3

    ! 先决定输运后端；后续 shell 循环保留同一物理源项和冷却定义。
    ! Resolve the transport backend first; later shell loops keep the same source and cooling definitions.
    active_solver=resolve_solver(solver_id)
    crossed_reverse=(T_cross > 0d0 .and. R_cross > 0d0 .and. M3_cross > 0d0 .and. V3_cross > 0d0)
    pc_ready=.false.
    pc_scale=1d0
    pc_shift=0d0
    if (maxval(B3) <= 0d0 .or. maxval(M3_shell) <= M3_shell(1)) then
        do i_empty=1,Num_gam_e
            if (Num_gam_e == 1) then
                gam_e(i_empty)=1d0
            else
                gam_e(i_empty)=dexp(dlog(1d1)*dble(i_empty-1)/dble(Num_gam_e-1))
            end if
        end do
        dN_gam_e=0d0
        deallocate(dEl,x,dF1,temp3,dN_x,x_edge,coord_edge,pc_log,pc_map, &
                   pc_work,pc_back,pc_u,pc_a, &
                   pc_b,pc_affine,dB3_serial,P_syn,Seed_syn,cooling_aux,Compton)
        return
    end if

    ! 初始 gamma 网格覆盖整段 RS 可能达到的最高注入能量。
    ! The initial gamma grid spans the largest injection energy reachable over the RS history.
    factor2=(p_r-2d0)/(p_r-1d0)*e_r*Para_m_p_DIV_m_e
    if (p_r < 2.05d0) factor2=0.05d0/1.05d0*e_r*Para_m_p_DIV_m_e
    beta4=dsqrt(1d0-1d0/eta_0**2); u4=dsqrt(eta_0*eta_0-1d0)
    dB3_serial(1)=dB3_serial(min(2,Num_R))
    dB=dB3_serial(1); gamma34=1.001d0
    gmax=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
    gm=factor2*(gamma34-1d0)/f_e_r+1d0
    dgscale=fourvel_scale
    dglow=1d0
    dgmid=gm
    dghigh=gmax
    gc=gc_coeff/(1d0+dsqrt(e_r/b_r))/R_Gamma(1)/dB**2/(R_Tobs(1)/2d0)

    call density_profile(A_star,dNe_ISM,R(1),R0,1,R_tr,f_jump,f_wide,dNe)
    DB_min=bsyn_coeff*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-1d0)))
    gmax0=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    gmin=1d0

    ! 主循环每个输出壳层先重建局部 shock 状态，再推进电子谱。
    ! The main loop rebuilds local shock state for each output shell before advancing the spectrum.
    do I_tobs=2,Num_R
        gloc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/2d0
        beta2=dsqrt(1d0-1d0/gloc**2)
        u2=dsqrt(gloc*gloc-1d0)
        gamma34=(gloc*gloc+eta_0*eta_0-1d0)/(eta_0*gloc+u2*u4)
        dB=(dB3_serial(I_tobs)+dB3_serial(I_tobs-1))/2d0
        gmax=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        gm=factor2*(gamma34-1d0)/f_e_r+1d0
    end do
    if (gmax0 <= gmin) error stop "electron_reverse_evolve: reverse electron grid maximum must exceed minimum."

    do I_gam_e=1,Num_gam_e
        if (Num_gam_e == 1) then
            gam_e(I_gam_e)=gmin
        else
            gam_e(I_gam_e)=gmin*dexp(dlog(gmax0/gmin)*(I_gam_e-1)/(Num_gam_e-1))
        end if
        dN_gam_e(I_gam_e,1)=0d0
    end do

    coord_scale=dgscale*dgscale-1d0
    call build_fourvel_grid(Num_gam_e,gmin,gmax0,dgscale,gam_e,coord_edge,x_edge)
    dN_x=0d0
    d_x=coord_edge(2)-coord_edge(1)
    rloc=R(1)
    if (active_solver == solver_dg) call init_dg_state()

    do I_tobs=2,Num_R
        rloc=R(I_tobs-1)
        gloc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/2d0
        Delta=max(Delta_0,rloc/Eta_0**2)
        rn4=para_m_ej/(4d0*pi*Para_m_p*rloc*rloc*Eta_0*Delta)
        beta2=dsqrt(1d0-1d0/gloc**2)
        u2=dsqrt(gloc*gloc-1d0)
        gamma34=(gloc*gloc+eta_0*eta_0-1d0)/(eta_0*gloc+u2*u4)
        dB=(dB3_serial(I_tobs)+dB3_serial(I_tobs-1))/2d0
        gmax=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        gm=factor2*(gamma34-1d0)/f_e_r+1d0
        gc=gc_coeff*(1d0+z)/gloc/dB**2/R_Tobs(I_tobs)
        f_r=adv_coeff/beta2/gloc*dB**2/pi
        dDR=0.7d0/(f_r*gmax+1.333d0/(R(I_tobs)+R(I_tobs-1)))
        dDD=R(I_tobs)-R(I_tobs-1)
        if (dB <= 0d0) then
            dN_gam_e(:,I_tobs)=dN_gam_e(:,I_tobs-1)
            cycle
        end if
        thermloss=0d0
        if (crossed_reverse .and. R(I_tobs) > R_cross) then
            if (R(I_tobs-1) < R_cross) then
                vol_lo=V3_cross
                vol_hi=V3_shell(I_tobs)
            else
                vol_lo=V3_shell(I_tobs-1)
                vol_hi=V3_shell(I_tobs)
            end if
            if (vol_lo <= 0d0 .or. vol_hi <= 0d0) &
                error stop "electron_reverse_evolve: post-crossing comoving volume must be positive."
            if (R(I_tobs-1) < R_cross) then
                thermloss=dlog(vol_hi/vol_lo)/(3d0*(R(I_tobs)-R_cross))
            else
                thermloss=dlog(vol_hi/vol_lo)/(3d0*dDD)
            end if
            if (thermloss < 0d0) &
                error stop "electron_reverse_evolve: post-crossing comoving volume must not decrease."
        end if
        injection_rate=0d0
        if ((.not. crossed_reverse) .or. R(I_tobs-1) < R_cross) then
            if (crossed_reverse) then
                inj_hi=min(R(I_tobs),R_cross)
            else
                inj_hi=R(I_tobs)
            end if
            if (inj_hi > R(I_tobs-1)) then
                mass_lo=M3_shell(I_tobs-1)
                mass_hi=M3_shell(I_tobs)
                if (crossed_reverse .and. R(I_tobs) > R_cross) mass_hi=M3_cross
                if (mass_hi < mass_lo) error stop "electron_reverse_evolve: reverse swept mass must not decrease."
                inj_width=inj_hi-R(I_tobs-1)
                injection_rate=f_e_r*(mass_hi-mass_lo)/(Para_m_p*inj_width)
            end if
        end if
        L1=reverse_transport_substeps(dDR,dDD,active_solver)
        dDR=dDD/L1
        call compute_cooling(I_tobs)
        call advance_shell(I_tobs)
    end do

    deallocate(dEl,x,dF1,temp3,dN_x,x_edge,coord_edge,pc_log,pc_map, &
               pc_work,pc_back,pc_u,pc_a, &
               pc_b,pc_affine,dB3_serial,P_syn,Seed_syn,cooling_aux,Compton)
    if (allocated(dg_state)) deallocate(dg_state)
    if (allocated(dg_work)) deallocate(dg_work,dg_dEl,dg_source)

contains

    ! 当前壳层同步辐射和 IC 冷却；index_Y 控制公开冷却模式。
    ! Synchrotron and IC cooling for the current shell; index_Y selects the public cooling mode.
    subroutine compute_cooling(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs
    real(8), dimension(Num_nu,1) :: P_syn_column,Seed_syn_column
    real(8), dimension(Num_gam_e,1) :: cooling_aux_column
    real(8), dimension(Num_nu) :: P_emit_tmp,Tau_syn_tmp

        call syn_state(index_syn_intger,R(I_tobs-1),dB,Num_gam_e,Num_nu,n_threads, &
                                    gam_e,dN_gam_e(:,I_tobs-1),V_seed,P_emit_tmp,P_syn,Seed_syn,Tau_syn_tmp)
        P_syn_column(:,1)=P_syn
        Seed_syn_column(:,1)=Seed_syn
        call forward_cooling(0,index_Y,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0, &
                             Num_gam_e,Num_nu,1,n_threads,gam_e,V_seed, &
                             P_syn_column,Seed_syn_column,Seed_syn_column,cooling_aux_column,dEl)
        cooling_aux=cooling_aux_column(:,1)
        select case(index_Y)
        case(0)
            dEl=f_r*gam_e
        case(1)
            coolscale=1d0/(beta2*gloc*Para_c)
            dEl=(f_r+cooling_aux*coolscale)*gam_e
        case(2)
            Qshell=4d0*pi*R(I_tobs-1)*R(I_tobs-1)*Para_c
            Compton=1d0+cooling_aux/Qshell/(4d0*gloc*gloc*rn4*Para_m_p_E)
            gmax=gmax/dsqrt(Compton(Num_gam_e))
            dEl=f_r*Compton*gam_e
        case default
            error stop 'electron_reverse_evolve: index_Y must be 0, 1, or 2.'
        end select
    end subroutine compute_cooling

    ! 推进一个输出壳层：pre-crossing 有注入，post-crossing 只做冷却和膨胀。
    ! Advance a single output shell: pre-crossing injects particles; post-crossing only cools and expands.
    subroutine advance_shell(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs
    integer :: i_node
    logical :: postonly

        postonly = (crossed_reverse .and. R(I_tobs - 1) >= R_cross)
        if (postonly .and. active_solver == solver_fullhide) then
            call advance_postcross(I_tobs)
            return
        end if
        if (active_solver == solver_dg) then
            call remesh_dg_state()
            if (index_Y /= 0) then
                do i_node=1,dg_mesh%ntot
                    dg_dEl(i_node)=log_interp(Num_gam_e,x_edge,dEl,dg_mesh%x_gamma(i_node))
                end do
            end if
            do L=1,L1
                rhi=rloc+dDR
                if (postonly) then
                    call advance_dgpart(I_tobs,dDR,.false.)
                else if (crossed_reverse .and. rloc < R_cross .and. rhi > R_cross) then
                    call advance_dgpart(I_tobs,R_cross-rloc,.true.)
                    rloc=R_cross
                    call advance_dgpart(I_tobs,rhi-rloc,.false.)
                else
                    call advance_dgpart(I_tobs,dDR,(.not. crossed_reverse) .or. rhi <= R_cross)
                end if
                rloc=rhi
            end do
            call dg_project_cells(dg_mesh,dg_state,Num_gam_e,coord_edge,dF1)
            call coord_to_dgamma(Num_gam_e,coord_edge,coord_scale, &
                                                               gam_e,dF1,dN_gam_e(:,I_tobs))
            where (dN_gam_e(:,I_tobs) <= 0d0) dN_gam_e(:,I_tobs)=0d0
            dN_x=dN_gam_e(:,I_tobs)*gam_e
            return
        end if

        do L=1,L1
            rhi=rloc+dDR
            if (crossed_reverse .and. rloc < R_cross .and. rhi > R_cross) then
                call advance_fullhidepart(I_tobs,R_cross-rloc,.true.)
                rloc=R_cross
                call advance_fullhidepart(I_tobs,rhi-rloc,.false.)
            else
                call advance_fullhidepart(I_tobs,dDR,(.not. crossed_reverse) .or. rhi <= R_cross)
            end if
            rloc=rhi
        end do
        call coord_to_dgamma(Num_gam_e,coord_edge,coord_scale, &
                                                        gam_e,dN_x,dN_gam_e(:,I_tobs))
    end subroutine advance_shell

    ! 在 RS 穿越事件单侧推进 DG 谱，使注入和绝热系数不跨越物理分支。
    ! Advance DG on one side of RS crossing so injection and adiabatic coefficients stay branch-local.
    subroutine advance_dgpart(I_tobs,drseg,preseg)
    implicit none
    integer, intent(in) :: I_tobs
    real(8), intent(in) :: drseg
    logical, intent(in) :: preseg
    real(8) :: srcseg

        R_step_mid=rloc+0.5d0*drseg
        if (index_Y == 0) call prepare_substep_state(I_tobs,R_step_mid)
        if (preseg) then
            srcseg=injection_rate
            adrate=1d0/R_step_mid
        else
            srcseg=0d0
            adrate=thermloss
        end if
        if (srcseg > 0d0) then
            call dg_kinetic_source(dg_mesh,srcseg,p_r,gm,gmax,dg_source)
        else
            dg_source=0d0
        end if
        call dg_scale_content(dg_mesh,srcseg,dg_source)
        if (index_Y == 0) then
            call dg_char_step(dg_mesh,drseg,f_r,adrate,dg_source,dg_state,dg_work)
        else
            call dg_advance_step(dg_mesh,adrate,drseg,dg_dEl,dg_source,dg_state,dg_work)
            call dg_limit_positive(dg_mesh,dg_work)
        end if
        dg_state=dg_work
        if (srcseg > 0d0 .and. dglow <= 1d0) dglow=gm
        if (dglow > 1d0) call advance_dg_front(dglow,drseg)
        if (srcseg > 0d0) dglow=min(dglow,gm)
        dglow=max(1d0,dglow)
        call advance_dg_inject(srcseg,drseg)
        if (dghigh > 1d0) call advance_dg_front(dghigh,drseg)
        if (srcseg > 0d0) dghigh=max(dghigh,gmax)
    end subroutine advance_dgpart

    ! Advance one side of the exact reverse-shock crossing event.
    subroutine advance_fullhidepart(I_tobs,drseg,preseg)
    implicit none
    integer, intent(in) :: I_tobs
    real(8), intent(in) :: drseg
    logical, intent(in) :: preseg

        R_step_mid=rloc+0.5d0*drseg
        if (index_Y == 0) call prepare_substep_state(I_tobs,R_step_mid)
        if (preseg) then
            call kinetic_coord(Num_gam_e,coord_edge,coord_scale, &
                                                       gm,gmax,injection_rate,p_r,dF1)
            adrate=1d0/R_step_mid
        else
            dF1=0d0
            adrate=thermloss
        end if
        call shell_coord_step(Num_gam_e,drseg,coord_edge,coord_scale, &
                                                   dEl,adrate,dF1,dN_x,x)
        dN_x=x
    end subroutine advance_fullhidepart

    ! Crossing 后没有新注入，沿特征线搬运既有电子谱。
    ! After crossing there is no new injection; advect the existing spectrum along characteristics.
    subroutine advance_postcross(I_tobs)
    implicit none
    integer, intent(in) :: I_tobs

        dF1=0d0
        if (.not. pc_ready) then
            dF1=dN_x
            pc_log=dF1*(coord_edge(2:Num_gam_e+1)-coord_edge(1:Num_gam_e))/ &
                         (x_edge(2:Num_gam_e+1)-x_edge(1:Num_gam_e))
            if (index_Y /= 0) then
                call u_edges(Num_gam_e,x_edge,pc_u)
                pc_map=pc_u
            end if
            pc_ready=.true.
        end if
        do L=1,L1
            R_step_mid=rloc+0.5d0*dDR
            if (index_Y == 0) then
                call prepare_substep_state(I_tobs,R_step_mid)
                if (abs(thermloss) <= 1d-30) then
                    step_scale=1d0
                    step_shift=-f_r*dDR
                else
                    step_scale=dexp(-thermloss*dDR)
                    step_shift=(f_r/thermloss)*(step_scale-1d0)
                end if
                pc_shift=pc_scale*step_shift+pc_shift
                pc_scale=pc_scale*step_scale
            else
                call build_piece_u(Num_gam_e,x_edge,gam_e,dEl,thermloss, &
                                                       pc_u,pc_a,pc_b)
                call trace_piece_u(Num_gam_e,1,pc_u,pc_u, &
                                   pc_a,pc_b,[dDR],pc_affine)
                pc_back=pc_affine(:,1)
                do i_edge=1,Num_gam_e+1
                    pc_work(i_edge)=postcross_umap(dexp(-pc_back(i_edge)))
                end do
                pc_map=pc_work
            end if
            rloc=rloc+dDR
        end do
        if (index_Y == 0) then
            do i_edge=1,Num_gam_e+1
                pc_back(i_edge)=x_from_u(pc_scale*dexp(-x_edge(i_edge))+pc_shift)
            end do
        else
            do i_edge=1,Num_gam_e+1
                pc_back(i_edge)=x_from_u(pc_map(i_edge))
            end do
        end if
        call remap_edges(Num_gam_e,x_edge,pc_back,pc_log,dN_x,.true.)
        call dnx_dgamma(Num_gam_e,x_edge,gam_e,dN_x,dN_gam_e(:,I_tobs))
        where (dN_gam_e(:,I_tobs) <= 0d0) dN_gam_e(:,I_tobs)=0d0
    end subroutine advance_postcross

    real(8) function postcross_umap(u_eval) result(u_cross)
    implicit none
    integer :: left,right,mid
    real(8), intent(in) :: u_eval
    real(8) :: frac

        if (u_eval >= pc_u(1)) then
            frac=(u_eval-pc_u(1))/(pc_u(2)-pc_u(1))
            u_cross=pc_map(1)+frac*(pc_map(2)-pc_map(1))
            return
        end if
        if (u_eval <= pc_u(Num_gam_e+1)) then
            frac=(u_eval-pc_u(Num_gam_e))/(pc_u(Num_gam_e+1)-pc_u(Num_gam_e))
            u_cross=pc_map(Num_gam_e)+ &
                    frac*(pc_map(Num_gam_e+1)-pc_map(Num_gam_e))
            return
        end if
        left=1
        right=Num_gam_e
        do while (left < right)
            mid=(left+right)/2
            if (pc_u(mid+1) <= u_eval) then
                right=mid
            else
                left=mid+1
            end if
        end do
        frac=(u_eval-pc_u(left))/(pc_u(left+1)-pc_u(left))
        u_cross=pc_map(left)+frac*(pc_map(left+1)-pc_map(left))
    end function postcross_umap

    subroutine init_dg_state()
    implicit none

        call dg_build_mesh(x_edge(1),dg_active_xmax(),dlog(dg_low_break()), &
                                                    dlog(dg_inject_break()), &
                                                    dlog(dg_upper_break()),dgscale,dg_mesh)
        allocate(dg_state(dg_mesh%ntot))
        dg_state=0d0
    end subroutine init_dg_state

    subroutine remesh_dg_state()
    implicit none

        call dg_build_mesh(x_edge(1),dg_active_xmax(),dlog(dg_low_break()), &
                                                    dlog(dg_inject_break()), &
                                                    dlog(dg_upper_break()),dgscale,dg_new_mesh)
        call ensure_dg_work(dg_new_mesh%ntot)
        call dg_project_state(dg_mesh,dg_state,dg_new_mesh,dg_work)
        if (size(dg_state) /= dg_new_mesh%ntot) then
            deallocate(dg_state)
            allocate(dg_state(dg_new_mesh%ntot))
        end if
        dg_state=dg_work
        dg_mesh=dg_new_mesh
    end subroutine remesh_dg_state

    subroutine ensure_dg_work(ntot)
    implicit none
    integer, intent(in) :: ntot

        if (allocated(dg_work)) then
            if (size(dg_work) /= ntot) deallocate(dg_work,dg_dEl,dg_source)
        end if
        if (.not. allocated(dg_work)) allocate(dg_work(ntot),dg_dEl(ntot),dg_source(ntot))
    end subroutine ensure_dg_work

    subroutine prepare_substep_state(I_tobs,radius_eval)
    implicit none
    integer, intent(in) :: I_tobs
    real(8), intent(in) :: radius_eval

        gloc=shell_linear_value(I_tobs,R_Gamma,radius_eval)
        beta2=dsqrt(1d0-1d0/gloc**2)
        u2=dsqrt(gloc*gloc-1d0)
        gamma34=(gloc*gloc+eta_0*eta_0-1d0)/(eta_0*gloc+u2*u4)
        dB=shell_linear_value(I_tobs,dB3_serial,radius_eval)
        gmax=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        gm=factor2*(gamma34-1d0)/f_e_r+1d0
        f_r=adv_coeff/beta2/gloc*dB**2/pi
        dEl=f_r*gam_e
    end subroutine prepare_substep_state

    real(8) function shell_linear_value(i_shell,values,radius_eval) result(value)
    implicit none
    integer, intent(in) :: i_shell
    real(8), intent(in), dimension(Num_R) :: values
    real(8), intent(in) :: radius_eval
    real(8) :: frac

        frac=(radius_eval-R(i_shell-1))/(R(i_shell)-R(i_shell-1))
        value=values(i_shell-1)+frac*(values(i_shell)-values(i_shell-1))
    end function shell_linear_value

    real(8) function dg_upper_break() result(gamma_break)
    implicit none

        if ((.not. crossed_reverse) .or. rloc < R_cross) then
            gamma_break=gmax
        else
            gamma_break=dghigh
        end if
    end function dg_upper_break

    real(8) function dg_active_xmax() result(x_max)
    implicit none
    real(8) :: gamma_grid_max,tail_fraction

        gamma_grid_max=dexp(x_edge(Num_gam_e+1))
        x_max=dlog(min(gamma_grid_max,tail_factor*dg_upper_break()))
        if (allocated(dg_state)) then
            if (x_max < dg_mesh%x_gamma(dg_mesh%ntot)) then
                call dg_tail_fraction(dg_mesh,dg_state,dexp(x_max), &
                                                        tail_power,tail_fraction)
                if (tail_fraction > tail_thresh) x_max=dg_mesh%x_gamma(dg_mesh%ntot)
            end if
        end if
    end function dg_active_xmax

    real(8) function dg_low_break() result(gamma_break)
    implicit none

        gamma_break=1d0
        if (dglow > 1d0) gamma_break=min(dglow,gm)
        ! Resolve the trans-relativistic kinetic-source singularity when gamma_m is nearly 1.
        if (gm < trans_break .and. gmax > trans_break) &
            gamma_break=trans_break
    end function dg_low_break

    subroutine advance_dg_inject(source_norm,drseg)
    implicit none
    real(8), intent(in) :: source_norm,drseg

        if (source_norm > 0d0) then
            dgmid=gm
        else if (dgmid > 1d0) then
            call advance_dg_front(dgmid,drseg)
        end if
        dgmid=max(1d0,dgmid)
    end subroutine advance_dg_inject

    real(8) function dg_inject_break() result(gamma_break)
    implicit none

        gamma_break=gm
        if (crossed_reverse .and. rloc >= R_cross .and. dgmid > 1d0) gamma_break=dgmid
    end function dg_inject_break

    subroutine advance_dg_front(gamma_front,drseg)
    implicit none
    real(8), intent(inout) :: gamma_front
    real(8), intent(in) :: drseg
    real(8) :: exp_b,loss_front

        if (index_Y == 0) then
            if (abs(adrate*drseg) <= 1d-10) then
                gamma_front=gamma_front/(1d0+f_r*gamma_front*drseg)
            else
                exp_b=dexp(-adrate*drseg)
                gamma_front=gamma_front*exp_b/(1d0+(f_r/adrate)*gamma_front*(1d0-exp_b))
            end if
        else
            loss_front=log_interp(Num_gam_e,x_edge,dEl,dlog(gamma_front))
            gamma_front=gamma_front*dexp(-(loss_front+adrate)*drseg)
        end if
    end subroutine advance_dg_front
end subroutine electron_reverse_evolve

! Secondary RS 源项历史：逐半径壳层调用同步辐射和SSA核，返回给统一EATS投影使用。
subroutine multiple_synch(index_syn_intger,Num_nu,Num_R,Num_gam_e, &
                                                  n_threads,R,R_Gamma,B3,gam_e, &
                                                  dN_gam_e,V_seed,z,L_syn_spec,Seed_syn,Nu_a)
    implicit none
    integer, intent(in) :: index_syn_intger,Num_nu,Num_R,Num_gam_e,n_threads
    integer :: I_tobs
    real(8), intent(in), dimension(Num_R) :: R,R_Gamma,B3
    real(8), intent(in), dimension(Num_gam_e) :: gam_e
    real(8), intent(in), dimension(Num_gam_e,Num_R) :: dN_gam_e
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in) :: z
    real(8), intent(out), dimension(Num_nu,Num_R) :: L_syn_spec,Seed_syn
    real(8), intent(out), dimension(Num_R) :: Nu_a

    L_syn_spec=0d0
    Seed_syn=0d0
    Nu_a=0d0
    if (n_threads > 1) then
        !$OMP PARALLEL DO SCHEDULE(STATIC) num_threads(n_threads) private(I_tobs)
        do I_tobs=1,Num_R
            call write_shell(I_tobs,1)
        end do
        !$OMP END PARALLEL DO
    else
        do I_tobs=1,Num_R
            call write_shell(I_tobs,n_threads)
        end do
    end if

contains

    subroutine write_shell(I_shell,syn_threads)
    implicit none
    integer, intent(in) :: I_shell,syn_threads
    real(8), dimension(Num_nu) :: P_emit_tmp,Tau_syn_tmp
    real(8) :: doppler_den

        if (B3(I_shell) <= 0d0) return
        call syn_state(index_syn_intger,R(I_shell),B3(I_shell),Num_gam_e,Num_nu,syn_threads, &
                                    gam_e,dN_gam_e(:,I_shell),V_seed,P_emit_tmp,L_syn_spec(:,I_shell), &
                                    Seed_syn(:,I_shell),Tau_syn_tmp)
        doppler_den=R_Gamma(I_shell)*(1d0-dsqrt(1d0-R_Gamma(I_shell)**(-2)))*(1d0+z)
        call nua_fromtau(Num_nu,V_seed,Tau_syn_tmp,Nu_a(I_shell))
        Nu_a(I_shell)=Nu_a(I_shell)/doppler_den
    end subroutine write_shell
end subroutine multiple_synch

! 二次激波分支电子演化：每个 branch 固定 DG mesh，壳层 gamma_m 只定义注入源。
! Secondary-shock branch evolution: fix one DG mesh per branch and use shell gamma_m only for injection.
subroutine branch_reaccel(e_r,b_r,p_r,f_e_r,z,R_Tobs,R_Gamma,R, &
                                                           B3_branch,M3_branch,U3_branch,V3_branch, &
                                                           Gam_m_branch,Gamma43_branch,Comp_branch,Parent_branch, &
                                                           V_seed,Num_jump,Num_nu,Num_R,Num_gam_e,index_syn_intger, &
                                                           n_threads,gam_e,dN_total,Branch_L_syn_spec,L_syn_spec, &
                                                           Branch_seed_energy,Branch_reaccel_energy,solver_id)
    implicit none
    integer, intent(in) :: Num_jump,Num_nu,Num_R,Num_gam_e,index_syn_intger,n_threads
    integer, intent(in), dimension(Num_jump) :: Parent_branch
    integer, intent(in), optional :: solver_id
    integer :: I_tobs,I_jump,I_gam_e,L1,L,parent,active_solver,maxNode
    real(8), intent(in) :: e_r,b_r,p_r,f_e_r,z
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_jump,Num_R) :: B3_branch, M3_branch, U3_branch, V3_branch
    real(8), intent(in), dimension(Num_jump,Num_R) :: Gam_m_branch, Gamma43_branch, Comp_branch
    real(8), intent(out), dimension(Num_gam_e) :: gam_e
    real(8), intent(out), dimension(Num_gam_e,Num_R) :: dN_total
    real(8), intent(out), dimension(Num_jump,Num_nu,Num_R) :: Branch_L_syn_spec
    real(8), intent(out), dimension(Num_nu,Num_R) :: L_syn_spec
    real(8), intent(out), dimension(Num_jump,Num_R) :: Branch_seed_energy,Branch_reaccel_energy
    real(8), parameter :: secondary_adv_coeff=1.35d-19
    real(8) :: dB,gmax,gm,gmax0,gmin,d_x,rloc,gloc,beta2
    real(8) :: f_r,dDR,dDD,injection_rate,mass_lo,mass_hi,mass_delta,fresh_mass,adrate
    real(8) :: parent_mass,transfer_mass,transfer_fraction,boost_factor,seed_energy,out_energy
    logical :: dissipative_shell
    real(8), allocatable, dimension(:) :: dEl,x,dF1,temp3,x_edge,coordEdge
    real(8), allocatable, dimension(:,:) :: dN_x,dN_work,dN_branch,spec_branch,branchState
    real(8), allocatable, dimension(:) :: dN_seed,dN_boost,dN_reaccel,nu_a_dummy
    real(8), allocatable, dimension(:,:,:) :: dN_gamma_branch
    real(8), allocatable, dimension(:,:) :: seed_dummy
    real(8), allocatable, dimension(:) :: branch_mass_available,fresh_mass_branch
    type(dg_mesh), allocatable, dimension(:) :: branch_mesh

    active_solver=resolve_solver(solver_id)
    allocate(dEl(Num_gam_e),x(Num_gam_e),dF1(Num_gam_e),temp3(Num_gam_e-1),x_edge(Num_gam_e+1), &
             dN_x(Num_jump,Num_gam_e),dN_work(Num_jump,Num_gam_e),dN_seed(Num_gam_e),dN_boost(Num_gam_e), &
             dN_reaccel(Num_gam_e),dN_gamma_branch(Num_jump,Num_gam_e,Num_R), &
             dN_branch(Num_gam_e,Num_R),spec_branch(Num_nu,Num_R),seed_dummy(Num_nu,Num_R), &
             nu_a_dummy(Num_R),branch_mass_available(Num_jump),fresh_mass_branch(Num_jump))

    call reaccel_grid()
    call log_edges(Num_gam_e,gam_e,x_edge)
    if (active_solver == solver_dg) then
        allocate(branch_mesh(Num_jump))
        call build_meshes(branch_mesh)
        maxNode=0
        do I_jump=1,Num_jump
            maxNode=max(maxNode,branch_mesh(I_jump)%ntot)
        end do
        allocate(branchState(maxNode,Num_jump),coordEdge(Num_gam_e+1))
        branchState=0d0
        do I_gam_e=1,Num_gam_e+1
            coordEdge(I_gam_e)=coord_from_xg(coord_fourvel,fourvel_scale,x_edge(I_gam_e))
        end do
    end if
    d_x=dlog(gam_e(2)/gam_e(1))
    dN_gamma_branch=0d0; dN_total=0d0; dN_x=0d0; branch_mass_available=0d0
    fresh_mass_branch=0d0
    Branch_L_syn_spec=0d0; L_syn_spec=0d0; Branch_seed_energy=0d0; Branch_reaccel_energy=0d0

    do I_tobs=2,Num_R
        if (active_solver /= solver_dg) &
            dN_work=dN_gamma_branch(:,:,I_tobs-1)*spread(gam_e,1,Num_jump)
        call transfer_parent(I_tobs)
        do I_jump=1,Num_jump
            call advance_reaccel(I_tobs,I_jump)
        end do
    end do
    do I_jump=1,Num_jump
        dN_branch=dN_gamma_branch(I_jump,:,:)
        dN_total=dN_total+dN_branch
        if (.not. any(dN_branch > 0d0)) cycle
        call multiple_synch(index_syn_intger,Num_nu,Num_R,Num_gam_e,n_threads,R,R_Gamma, &
                                                    B3_branch(I_jump,:),gam_e,dN_branch,V_seed,z, &
                                                    spec_branch,seed_dummy,nu_a_dummy)
        Branch_L_syn_spec(I_jump,:,:)=spec_branch
        L_syn_spec=L_syn_spec+spec_branch
    end do
    if (active_solver == solver_dg) deallocate(branch_mesh,branchState,coordEdge)
    deallocate(dEl,x,dF1,temp3,x_edge,dN_x,dN_work,dN_seed,dN_boost,dN_reaccel,dN_gamma_branch, &
               dN_branch,spec_branch,seed_dummy,nu_a_dummy,branch_mass_available,fresh_mass_branch)

contains

    subroutine reaccel_grid()
    implicit none

        gmin=1d0
        gmax0=0d0
        do I_tobs=2,Num_R
            do I_jump=1,Num_jump
                dB=(B3_branch(I_jump,I_tobs)+B3_branch(I_jump,I_tobs-1))/2d0
                if (dB > 0d0) then
                    gmax=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
                    gmax0=max(gmax0,gmax)
                end if
            end do
        end do
        if (gmax0 <= gmin) &
            error stop "branch_reaccel: empty electron grid."
        do I_gam_e=1,Num_gam_e
            if (Num_gam_e == 1) then
                gam_e(I_gam_e)=gmin
            else
                gam_e(I_gam_e)=gmin*dexp(dlog(gmax0/gmin)*(I_gam_e-1)/(Num_gam_e-1))
            end if
        end do
    end subroutine reaccel_grid

    ! 一次扫描各分支的活动壳层，以最大 gamma_m 固定 mesh；未耗散分支保持未构造。
    ! Fix active-branch meshes from maximum gamma_m; leave non-dissipative branches unbuilt.
    subroutine build_meshes(meshes)
    implicit none
    type(dg_mesh), intent(out), dimension(Num_jump) :: meshes
    integer :: i_branch,i_shell
    real(8) :: gm_ref,upper_break

        do i_branch=1,Num_jump
            gm_ref=0d0
            do i_shell=1,Num_R
                if (Gam_m_branch(i_branch,i_shell) > 1d0) &
                    gm_ref=max(gm_ref,Gam_m_branch(i_branch,i_shell))
            end do
            if (gm_ref <= 1d0) cycle
            upper_break=min(gmax0,20d0*gm_ref)
            call dg_build_mesh(x_edge(1),x_edge(Num_gam_e+1),dlog(gm_ref), &
                               dlog(upper_break),dlog(dsqrt(upper_break*gmax0)), &
                               fourvel_scale,meshes(i_branch))
        end do
    end subroutine build_meshes

    ! 仅在真实 parent transfer 时把相关 DG 状态投到公共 gamma 网格，转移后守恒投回。
    ! Project only real parent transfers to the common gamma grid and conservatively return both DG states.
    subroutine transfer_parent(i_shell)
    implicit none
    integer, intent(in) :: i_shell

        fresh_mass_branch=0d0
        do I_jump=1,Num_jump
            mass_delta=M3_branch(I_jump,i_shell)-M3_branch(I_jump,i_shell-1)
            if (mass_delta < 0d0) error stop "branch_reaccel: branch mass decreased."
            dissipative_shell=Gam_m_branch(I_jump,i_shell) > 1d0 .and. Gamma43_branch(I_jump,i_shell) > 1d0
            fresh_mass=0d0
            if (dissipative_shell) fresh_mass=mass_delta
            if (Parent_branch(I_jump) > 0 .and. dissipative_shell .and. mass_delta > 0d0) then
                parent=Parent_branch(I_jump)
                parent_mass=branch_mass_available(parent)
                if (parent_mass > 0d0) then
                    transfer_mass=min(mass_delta,parent_mass)
                    transfer_fraction=transfer_mass/parent_mass
                    if (active_solver == solver_dg) then
                        call export_dg(parent,dN_work(parent,:))
                        call export_dg(I_jump,dN_work(I_jump,:))
                    end if
                    dN_seed=dN_work(parent,:)*transfer_fraction
                    dN_work(parent,:)=dN_work(parent,:)-dN_seed
                    boost_factor=Gamma43_branch(I_jump,i_shell)
                    if (Comp_branch(I_jump,i_shell) > 1d0) &
                        boost_factor=boost_factor*Comp_branch(I_jump,i_shell)**(1d0/3d0)
                    call boost_log(boost_factor,dN_seed,dN_boost)
                    call dsa_reaccel(p_r,dN_boost,dN_reaccel)
                    seed_energy=log_energy(dN_seed)
                    out_energy=log_energy(dN_reaccel)
                    Branch_seed_energy(I_jump,i_shell)=Branch_seed_energy(I_jump,i_shell)+seed_energy
                    Branch_reaccel_energy(I_jump,i_shell)=Branch_reaccel_energy(I_jump,i_shell)+out_energy
                    dN_work(I_jump,:)=dN_work(I_jump,:)+dN_reaccel
                    if (active_solver == solver_dg) then
                        call import_dg(parent,dN_work(parent,:))
                        call import_dg(I_jump,dN_work(I_jump,:))
                    end if
                    branch_mass_available(parent)=branch_mass_available(parent)-transfer_mass
                    branch_mass_available(I_jump)=branch_mass_available(I_jump)+transfer_mass
                    fresh_mass=mass_delta-transfer_mass
                end if
            end if
            fresh_mass_branch(I_jump)=fresh_mass
            branch_mass_available(I_jump)=branch_mass_available(I_jump)+fresh_mass
            if (active_solver /= solver_dg) dN_x(I_jump,:)=dN_work(I_jump,:)
        end do
    end subroutine transfer_parent

    ! 把常驻 DG 状态守恒投影到公共 dN/dln(gamma) 网格；该输出不参与后续推进。
    ! Conservatively export persistent DG state to the common dN/dln(gamma) grid without changing it.
    subroutine export_dg(jumpIndex,dNx)
    implicit none
    integer, intent(in) :: jumpIndex
    real(8), intent(out), dimension(Num_gam_e) :: dNx
    real(8), dimension(Num_gam_e) :: dNcoord,dNdgamma
    real(8) :: dgContent,gridContent
    integer :: nodeCount

        nodeCount=branch_mesh(jumpIndex)%ntot
        call dg_integral(branch_mesh(jumpIndex),branchState(1:nodeCount,jumpIndex),dgContent)
        call dg_project_cells(branch_mesh(jumpIndex),branchState(1:nodeCount,jumpIndex), &
                              Num_gam_e,coordEdge,dNcoord)
        call coord_to_dgamma(Num_gam_e,coordEdge,branch_mesh(jumpIndex)%coord_scale, &
                             gam_e,dNcoord,dNdgamma)
        dNx=dNdgamma*gam_e
        gridContent=sum(dNx*(x_edge(2:Num_gam_e+1)-x_edge(1:Num_gam_e)))
        if (dgContent > 0d0) then
            if (gridContent <= 0d0) error stop "branch_reaccel: DG output projection lost positive content."
            dNx=dNx*(dgContent/gridContent)
        end if
    end subroutine export_dg

    ! 只在 parent transfer 后把公共网格状态守恒导回对应 branch mesh。
    ! Import a common-grid state conservatively to its branch mesh only after parent transfer.
    subroutine import_dg(jumpIndex,dNx)
    implicit none
    integer, intent(in) :: jumpIndex
    real(8), intent(in), dimension(Num_gam_e) :: dNx
    real(8) :: gridContent
    integer :: node,nodeCount

        nodeCount=branch_mesh(jumpIndex)%ntot
        do node=1,nodeCount
            branchState(node,jumpIndex)=log_interp(Num_gam_e,x_edge,dNx, &
                                                   branch_mesh(jumpIndex)%x_gamma(node))* &
                                         branch_mesh(jumpIndex)%dxgamma_dcoord(node)
        end do
        gridContent=sum(dNx*(x_edge(2:Num_gam_e+1)-x_edge(1:Num_gam_e)))
        call dg_scale_content(branch_mesh(jumpIndex),gridContent,branchState(1:nodeCount,jumpIndex))
    end subroutine import_dg

    ! 常驻 DG 状态在本壳层直接推进；壳层输出仅由 export_dg 生成。
    ! Advance persistent DG state directly in this shell; export_dg alone creates the shell output.
    subroutine advance_reaccel(i_shell,jump_index)
    implicit none
    integer, intent(in) :: i_shell,jump_index
    real(8), dimension(Num_gam_e-1) :: face_speed
    real(8), allocatable, dimension(:) :: dg_adiabatic,dg_source_norm

        if (M3_branch(jump_index,i_shell) <= 0d0 .and. M3_branch(jump_index,i_shell-1) <= 0d0) then
            if (active_solver == solver_dg .and. branch_mesh(jump_index)%ntot > 0) &
                call export_dg(jump_index,dN_x(jump_index,:))
            dN_gamma_branch(jump_index,:,i_shell)=dN_x(jump_index,:)/gam_e
            return
        end if
        call prepare_branch_shell(i_shell,jump_index)
        call branch_inject(i_shell,jump_index)
        if (active_solver == solver_dg) then
            allocate(dg_adiabatic(L1),dg_source_norm(L1))
            do L=1,L1
                rloc=rloc+dDR
                dg_adiabatic(L)=adrate
                if (injection_rate > 0d0) then
                    dg_source_norm(L)=injection_rate
                else
                    dg_source_norm(L)=0d0
                end if
            end do
            call dg_sequence(branch_mesh(jump_index),L1,dDR,f_r,dg_adiabatic,dg_source_norm, &
                             p_r,gm,gmax,branchState(1:branch_mesh(jump_index)%ntot,jump_index))
            call export_dg(jump_index,dN_x(jump_index,:))
            dN_gamma_branch(jump_index,:,i_shell)=dN_x(jump_index,:)/gam_e
            deallocate(dg_adiabatic,dg_source_norm)
            return
        end if

        do L=1,L1
            rloc=rloc+dDR
            if (injection_rate > 0d0) then
                call kinetic_edges(Num_gam_e,x_edge,gm,gmax, &
                                                                          injection_rate,p_r,dF1)
            else
                dF1=0d0
            end if
            face_speed=((dEl(2:Num_gam_e)+dEl(1:Num_gam_e-1))/2d0+adrate)
            call flux_split_step(Num_gam_e,dDR,d_x,face_speed,dF1,dN_x(jump_index,:),x,.true.)
            dN_x(jump_index,:)=x
        end do
        dN_gamma_branch(jump_index,:,i_shell)=dN_x(jump_index,:)/gam_e
    end subroutine advance_reaccel

    subroutine prepare_branch_shell(i_shell,jump_index)
    implicit none
    integer, intent(in) :: i_shell,jump_index

        rloc=R(i_shell-1)
        gloc=(R_Gamma(i_shell)+R_Gamma(i_shell-1))/2d0
        beta2=dsqrt(1d0-1d0/gloc**2)
        dB=(B3_branch(jump_index,i_shell)+B3_branch(jump_index,i_shell-1))/2d0
        if (dB <= 0d0) error stop "branch_reaccel: active branch requires B3 > 0."
        if (Gam_m_branch(jump_index,i_shell-1) > 1d0) then
            gm=(Gam_m_branch(jump_index,i_shell)+Gam_m_branch(jump_index,i_shell-1))/2d0
        else
            gm=Gam_m_branch(jump_index,i_shell)
        end if
        gmax=3d0*Para_m_energy/dsqrt(8d0*dB*Para_e**3)
        f_r=secondary_adv_coeff/beta2/gloc*dB**2/pi
        dDD=R(i_shell)-R(i_shell-1)
        dDR=0.7d0/(f_r*gmax+1d0/R(i_shell-1))
        L1=reverse_transport_substeps(dDR,dDD,active_solver)
        dDR=dDD/L1
        dEl=f_r*gam_e
        if (V3_branch(jump_index,i_shell) <= 0d0 .or. V3_branch(jump_index,i_shell-1) <= 0d0) then
            adrate=0d0
        else
            adrate=dlog(V3_branch(jump_index,i_shell)/V3_branch(jump_index,i_shell-1))/(3d0*dDD)
        end if
    end subroutine prepare_branch_shell

    subroutine branch_inject(i_shell,jump_index)
    implicit none
    integer, intent(in) :: i_shell,jump_index

        mass_lo=M3_branch(jump_index,i_shell-1)
        mass_hi=M3_branch(jump_index,i_shell)
        if (mass_hi < mass_lo) error stop "branch swept mass decreased."
        injection_rate=f_e_r*fresh_mass_branch(jump_index)/(Para_m_p*dDD)
    end subroutine branch_inject

    subroutine boost_log(boost,dN_in,dN_out)
    implicit none
    real(8), intent(in), dimension(Num_gam_e) :: dN_in
    real(8), intent(in) :: boost
    real(8), intent(out), dimension(Num_gam_e) :: dN_out
    integer :: i_src,i_hi
    real(8) :: x_src,pos,frac

        if (boost <= 1d0) error stop "branch_reaccel: boost must exceed unity."
        dN_out=0d0
        do I_gam_e=1,Num_gam_e
            x_src=dlog(gam_e(I_gam_e)/boost)
            if (x_src < dlog(gam_e(1)) .or. x_src > dlog(gam_e(Num_gam_e))) cycle
            pos=(x_src-dlog(gam_e(1)))/d_x+1d0
            i_src=int(pos)
            if (i_src < 1) cycle
            if (i_src >= Num_gam_e) then
                dN_out(I_gam_e)=dN_in(Num_gam_e)
            else
                i_hi=i_src+1
                frac=pos-dble(i_src)
                dN_out(I_gam_e)=(1d0-frac)*dN_in(i_src)+frac*dN_in(i_hi)
            end if
        end do
    end subroutine boost_log

    subroutine dsa_reaccel(p,dN_seed_log,dN_out_log)
    implicit none
    real(8), intent(in), dimension(Num_gam_e) :: dN_seed_log
    real(8), intent(in) :: p
    real(8), intent(out), dimension(Num_gam_e) :: dN_out_log
    integer :: i
    real(8) :: integral,dN_dgamma

        if (p <= 2d0) error stop "branch_reaccel: DSA requires p > 2."
        integral=0d0
        dN_out_log=0d0
        do i=1,Num_gam_e
            dN_dgamma=dN_seed_log(i)/(gam_e(i))
            integral=integral+dN_dgamma*gam_e(i)**(p-1d0)*(gam_e(i)*d_x)
            dN_out_log(i)=(p-1d0)*gam_e(i)**(-p)*integral*gam_e(i)
        end do
    end subroutine dsa_reaccel

    real(8) function log_energy(dN_log)
    implicit none
    real(8), intent(in), dimension(Num_gam_e) :: dN_log

        log_energy=sum((gam_e-1d0)*Para_m_e*Para_c**2*dN_log)*d_x
    end function log_energy
end subroutine branch_reaccel

integer function reverse_transport_substeps(candidate_dr,shell_dr,solver) result(nstep)
    implicit none
    integer, intent(in) :: solver
    real(8), intent(in) :: candidate_dr,shell_dr

    if (solver == solver_dg) then
        nstep=reverse_dg_base_substeps
    else
        nstep=max(100,min(1000,int(shell_dr/candidate_dr)))
    end if
end function reverse_transport_substeps

! 在预建 mesh 上推进单个壳层；gamma_m/gamma_max 只用于真实源项投影。
! Advance one shell on a prebuilt mesh; gamma_m/gamma_max enter only the physical source projection.
subroutine dg_sequence(mesh,num_step,dR,rad_coeff,adrate_step,srcstep,p,gamma_m,gamma_max,state)
    implicit none
    type(dg_mesh), intent(in) :: mesh
    integer, intent(in) :: num_step
    integer :: step
    real(8), intent(in) :: dR,rad_coeff
    real(8), intent(in), dimension(num_step) :: adrate_step,srcstep
    real(8), intent(in) :: p,gamma_m,gamma_max
    real(8), intent(inout), dimension(mesh%ntot) :: state
    real(8), allocatable, dimension(:) :: advanced,source_nodes,source_template
    logical :: has_source

    allocate(advanced(mesh%ntot),source_nodes(mesh%ntot),source_template(mesh%ntot))
    has_source = maxval(srcstep) > 0d0
    if (has_source) then
        call dg_kinetic_source(mesh,1d0,p,gamma_m,gamma_max,source_template)
        call dg_scale_content(mesh,1d0,source_template)
    endif
    do step=1,num_step
        if (srcstep(step) > 0d0) then
            source_nodes=srcstep(step)*source_template
        else
            source_nodes=0d0
        end if
        call dg_char_step(mesh,dR,rad_coeff,adrate_step(step), &
                                                       source_nodes,state,advanced)
        state=advanced
    end do
    deallocate(advanced,source_nodes,source_template)
end subroutine dg_sequence

real(8) function log_interp(Num_gam_e,x_edge,values,x_eval) result(value)
    implicit none
    integer, intent(in) :: Num_gam_e
    integer :: i_cell,lo,hi,mid
    real(8), intent(in), dimension(Num_gam_e+1) :: x_edge
    real(8), intent(in), dimension(Num_gam_e) :: values
    real(8), intent(in) :: x_eval
    real(8) :: x_left,x_right

    x_left=0.5d0*(x_edge(1)+x_edge(2))
    if (x_eval <= x_left) then
        value=values(1)
        return
    end if
    x_right=0.5d0*(x_edge(Num_gam_e)+x_edge(Num_gam_e+1))
    if (x_eval >= x_right) then
        value=values(Num_gam_e)
        return
    end if
    lo=1
    hi=Num_gam_e-1
    do while (lo < hi)
        mid=(lo+hi)/2
        x_right=0.5d0*(x_edge(mid+1)+x_edge(mid+2))
        if (x_eval <= x_right) then
            hi=mid
        else
            lo=mid+1
        end if
    end do
    i_cell=lo
    x_left=0.5d0*(x_edge(i_cell)+x_edge(i_cell+1))
    x_right=0.5d0*(x_edge(i_cell+1)+x_edge(i_cell+2))
    value=values(i_cell)+(values(i_cell+1)-values(i_cell))*(x_eval-x_left)/(x_right-x_left)
end function log_interp

end module electron_reverse_kernel
