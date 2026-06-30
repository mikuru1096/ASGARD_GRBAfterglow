!f2py: skip
module electron_cooling_y_kernel
  use constants
  use electron_radiation_kernel, only: first_greater_monotonic_window, electron_powerlaw_interp, &
                                       electron_integrate_powerlaw_segment, electron_log_gauss2_interval
  private

  public :: electron_cooling_y_nakar, electron_cooling_y_fan

  integer, save :: y_nakar_num_gam_cache=0, y_nakar_num_nu_cache=0
  integer, allocatable, save :: y_nakar_idx_cache(:)
  real(8), allocatable, save :: y_nakar_hat_nu_cache(:), y_nakar_prefix_cache(:), y_nakar_v_cache(:), y_nakar_vloc_cache(:), &
                                y_nakar_vg1(:), y_nakar_vg2(:), y_nakar_wg1(:), y_nakar_wg2(:)

contains
! 确保Nakar Y参数工作数组已分配（缓存hat_nu、频率段Gauss节点和查找区间）。
subroutine ensure_y_nakar_workspace(Num_gam_e,Num_nu,gam_e,V_seed)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu)
logical :: rebuild

    rebuild=.not. allocated(y_nakar_hat_nu_cache)
    if (.not. rebuild) rebuild = (y_nakar_num_gam_cache /= Num_gam_e)
    if (.not. rebuild) rebuild = (y_nakar_num_nu_cache /= Num_nu)
    if (.not. rebuild) rebuild = any(y_nakar_hat_nu_cache /= Para_m_energy/Para_h/gam_e)
    if (.not. rebuild) rebuild = any(y_nakar_v_cache /= V_seed)
    if (.not. rebuild) return

    if (allocated(y_nakar_hat_nu_cache)) deallocate(y_nakar_hat_nu_cache,y_nakar_prefix_cache,y_nakar_v_cache, &
                                                    y_nakar_vloc_cache,y_nakar_idx_cache,y_nakar_vg1,y_nakar_vg2, &
                                                    y_nakar_wg1,y_nakar_wg2)
    allocate(y_nakar_hat_nu_cache(Num_gam_e),y_nakar_prefix_cache(Num_nu),y_nakar_v_cache(Num_nu), &
             y_nakar_vloc_cache(Num_gam_e),y_nakar_idx_cache(Num_gam_e),y_nakar_vg1(Num_nu-1), &
             y_nakar_vg2(Num_nu-1),y_nakar_wg1(Num_nu-1),y_nakar_wg2(Num_nu-1))
    y_nakar_hat_nu_cache=Para_m_energy/Para_h/gam_e
    y_nakar_v_cache=V_seed
    do I_nu=1,Num_nu-1
        call electron_log_gauss2_interval(V_seed(I_nu),V_seed(I_nu+1), &
                                          y_nakar_vg1(I_nu),y_nakar_vg2(I_nu), &
                                          y_nakar_wg1(I_nu),y_nakar_wg2(I_nu))
    end do
    do I_Compton=1,Num_gam_e
        if (y_nakar_hat_nu_cache(I_Compton) <= V_seed(1)) then
            y_nakar_idx_cache(I_Compton)=0
            y_nakar_vloc_cache(I_Compton)=V_seed(1)
        else
            call first_greater_monotonic_window(V_seed,Num_nu,2,y_nakar_hat_nu_cache(I_Compton),I_nu)
            y_nakar_idx_cache(I_Compton)=I_nu
            if (I_nu <= Num_nu) y_nakar_vloc_cache(I_Compton)=min(y_nakar_hat_nu_cache(I_Compton),V_seed(I_nu))
        end if
    end do
    y_nakar_num_gam_cache=Num_gam_e
    y_nakar_num_nu_cache=Num_nu
end subroutine ensure_y_nakar_workspace

! Nakar+2009 Compton Y参数：Y(γ) = ∫_{ν̂(γ)}^{ν_max} P_syn(ν) dν，谱形依赖。
subroutine electron_cooling_y_nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn, Compton)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu)
real(8), intent(out) :: Compton(Num_gam_e)
integer :: I_Compton,I_nu

    call ensure_y_nakar_workspace(Num_gam_e,Num_nu,gam_e,V_seed)

    Compton=zero
    y_nakar_prefix_cache(1)=zero
    do I_nu=2,Num_nu
       y_nakar_prefix_cache(I_nu)=y_nakar_prefix_cache(I_nu-1)+ &
            y_nakar_wg1(I_nu-1)*electron_powerlaw_interp(V_seed(I_nu-1),V_seed(I_nu), &
                                                          P_syn(I_nu-1),P_syn(I_nu),y_nakar_vg1(I_nu-1))+ &
            y_nakar_wg2(I_nu-1)*electron_powerlaw_interp(V_seed(I_nu-1),V_seed(I_nu), &
                                                          P_syn(I_nu-1),P_syn(I_nu),y_nakar_vg2(I_nu-1))
    end do

    do I_Compton=1,Num_gam_e
       ! Inlined from accumulate_y_nakar_point
       Compton(I_Compton)=zero
       I_nu=y_nakar_idx_cache(I_Compton)
       if (I_nu == 0) cycle
       if (I_nu <= Num_nu) then
          Compton(I_Compton)=y_nakar_prefix_cache(I_nu-1)+ &
                      electron_integrate_powerlaw_segment(V_seed(I_nu-1),y_nakar_vloc_cache(I_Compton), &
                          P_syn(I_nu-1), &
                          electron_powerlaw_interp(V_seed(I_nu-1),V_seed(I_nu), &
                                                   P_syn(I_nu-1),P_syn(I_nu),y_nakar_vloc_cache(I_Compton)))
       else
          Compton(I_Compton)=y_nakar_prefix_cache(Num_nu)
       end if
    end do
end subroutine electron_cooling_y_nakar

! Fan+2008 Compton Y参数：解析分段η_NK(γ) × η_rad，含快/慢冷却和谱指数分支。
subroutine electron_cooling_y_fan(Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e, Compton)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,gam_e(Num_gam_e)
real(8), intent(out) :: Compton(Num_gam_e)
integer :: i_gam_e

    eta=(Gam_e_m/Gam_e_c)**(p-two)
    if (eta-one > 0.001) eta=one

    do i_gam_e=1,Num_gam_e-1
        hat_gam=5.4246D6/sqrt(DB*gam_e(i_gam_e+1))
        if (Gam_e_m > Gam_e_c) then
            if (hat_gam < Gam_e_c) then
                eta_NK=zero
            else if (hat_gam < Gam_e_m) then
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
        else if (hat_gam < Gam_e_m) then
            eta_NK=zero
        else if (hat_gam < Gam_e_c) then
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
        Compton(i_gam_e)=0.5d0*(-1.0+sqrt(1.0+4.0*eta*eta_NK*Epsilon_e/Epsilon_b))
    end do
    Compton(Num_gam_e)=0.99*Compton(Num_gam_e-1)
end subroutine electron_cooling_y_fan

end module electron_cooling_y_kernel
