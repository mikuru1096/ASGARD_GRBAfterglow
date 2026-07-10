! 反向激波强子 light path：质子注入、冷却输运、质子同步辐射。
! Reverse-shock hadronic light path: proton injection, transport, and proton synchrotron.
subroutine reverse_hadronic_1d(R_Tobs,R_Gamma,R,shell_energy_inj_erg,B_field_g,V_seed, &
                                  include_proton_synch,Num_nu,Num_R,num_gam_p, &
                                  gam_p,dN_gam_p,P_had_syn,Seed_had_syn)
    use constants
    use hadronic_base
    use hadronic_shell, only: shell_geom
    use hadronic_transport
    use hadronic_rad, only: proton_syn
    implicit none
    integer, intent(in) :: include_proton_synch,Num_nu,Num_R,num_gam_p
    real(8), intent(in), dimension(Num_R) :: R_Tobs,R_Gamma,R,shell_energy_inj_erg
    real(8), intent(in), dimension(Num_R) :: B_field_g
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(out), dimension(num_gam_p) :: gam_p
    real(8), intent(out), dimension(num_gam_p,Num_R) :: dN_gam_p
    real(8), intent(out), dimension(Num_nu,Num_R) :: P_had_syn,Seed_had_syn
    integer :: ir
    real(8) :: gmax,tdyn,dt,dr,gmin,ebudget
    real(8), dimension(num_gam_p) :: dnprev,dnnext,qinj
    real(8), dimension(num_gam_p) :: lossad,losssyn,losstot

    call init_grid()
    dnprev=0d0

    dN_gam_p=0d0
    P_had_syn=0d0
    Seed_had_syn=0d0

    do ir=1,Num_R
        call advance_shell(ir)
        call emit_syn(ir)
        dnprev=dnnext
    end do

contains

    ! 用所有 RS shell 的最大可加速能量建立统一 proton gamma grid。
    ! Build a single proton gamma grid from the maximum acceleration limit over RS shells.
    subroutine init_grid()
        implicit none

        tdyn=dyn_time(R(1),R_Gamma(1))
        gmax=proton_limit(B_field_g(1),tdyn,1d1)
        do ir=2,Num_R
            tdyn=dyn_time(R(ir),R_Gamma(ir))
            gmax=max(gmax,proton_limit(B_field_g(ir),tdyn,1d1))
        end do
        if (gmax <= 1d0+1d-3) error stop "reverse hadronic gamma_p_max must exceed the injection grid minimum."
        call build_grid(num_gam_p,1d0+1d-3,gmax,gam_p)
    end subroutine init_grid

    ! 单个 RS shell：注入 proton source，再用 log-gamma transport 推进。
    ! One RS shell: inject proton source, then advance in log-gamma transport.
    subroutine advance_shell(ir)
        implicit none
        integer, intent(in) :: ir

        call shell_geom(Num_R,R,R_Gamma,ir,dr,dt)
        tdyn=dyn_time(R(ir),R_Gamma(ir))
        if (shell_energy_inj_erg(ir) < 0d0) error stop "reverse hadronic shell injection energy must be non-negative."
        ebudget=shell_energy_inj_erg(ir)
        gmin=max(gam_p(1),R_Gamma(ir))
        call proton_inject(num_gam_p,gam_p,2.2d0,ebudget,gmin,gam_p(num_gam_p),qinj)
        call proton_loss(num_gam_p,gam_p,B_field_g(ir),tdyn,lossad,losssyn,losstot)
        call advance_loggamma(num_gam_p,gam_p,dnprev,qinj,losstot,dt,dnnext)
        dN_gam_p(:,ir)=dnnext
    end subroutine advance_shell

    ! 可选 proton synch 输出，使用当前 shell transport 结果。
    ! Optional proton synchrotron output from the current shell spectrum.
    subroutine emit_syn(ir)
        implicit none
        integer, intent(in) :: ir

        if (include_proton_synch /= 0) then
            call proton_syn(R(ir),B_field_g(ir),num_gam_p,Num_nu,gam_p,dnnext, &
                            V_seed,P_had_syn(:,ir),Seed_had_syn(:,ir))
        end if
    end subroutine emit_syn
end subroutine reverse_hadronic_1d
