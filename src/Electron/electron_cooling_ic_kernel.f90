!f2py: skip
module electron_cooling_ic_kernel
  use constants
  use radiation_common, only: compute_simpson_weights
  use electron_radiation_kernel, only: electron_powerlaw_interp
  private

  public :: electron_cooling_ic_loss, electron_cooling_ic_loss_emissivity_budget

  integer, save :: ic_num_gam_cache=0, ic_num_nu_cache=0
  logical, save :: ic_grid_cache_ready=.false.
  real(8), allocatable, save :: ic_d_nu_cache(:), ic_gam_e_mean_cache(:), &
                                ic_e_seed_cache(:), ic_x_seed_cache(:), ic_v_seed_mid_cache(:)

contains
! 确保IC网格缓存已计算（种子频率中点值、间距、电子能量中点值等）。
subroutine ensure_ic_grid_cache(Num_gam_e,Num_nu,gam_e,V_seed)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu)

    if (ic_grid_cache_current()) return
    call rebuild_ic_grid_cache()

contains
logical function ic_grid_cache_current()
implicit REAL(8)(A-H,O-Z)

    ic_grid_cache_current=.false.
    if (.not. ic_grid_cache_ready) return
    if (ic_num_gam_cache /= Num_gam_e) return
    if (ic_num_nu_cache /= Num_nu) return
    if (.not. ic_seed_grid_current()) return
    if (.not. ic_gamma_grid_current()) return
    ic_grid_cache_current=.true.
end function ic_grid_cache_current

logical function ic_seed_grid_current()
implicit REAL(8)(A-H,O-Z)
integer :: I_nu

    ic_seed_grid_current=.false.
    if (.not. allocated(ic_x_seed_cache)) return
    do I_nu=1,Num_nu
        if (ic_x_seed_cache(I_nu) /= dlog(V_seed(I_nu))) return
    end do
    ic_seed_grid_current=.true.
end function ic_seed_grid_current

logical function ic_gamma_grid_current()
implicit REAL(8)(A-H,O-Z)
integer :: I_gam_e

    ic_gamma_grid_current=.false.
    if (.not. allocated(ic_gam_e_mean_cache)) return
    do I_gam_e=1,Num_gam_e-1
        if (ic_gam_e_mean_cache(I_gam_e) /= (gam_e(I_gam_e)+gam_e(I_gam_e+1))/two) return
    end do
    ic_gamma_grid_current=.true.
end function ic_gamma_grid_current

subroutine rebuild_ic_grid_cache()
implicit REAL(8)(A-H,O-Z)

    if (allocated(ic_d_nu_cache)) deallocate(ic_d_nu_cache,ic_gam_e_mean_cache,ic_e_seed_cache, &
                                             ic_x_seed_cache,ic_v_seed_mid_cache)
    allocate(ic_d_nu_cache(Num_nu-1),ic_gam_e_mean_cache(Num_gam_e-1),ic_e_seed_cache(Num_nu-1), &
             ic_x_seed_cache(Num_nu),ic_v_seed_mid_cache(Num_nu-1))

    para_hEme=Para_h/para_m_energy
    ic_x_seed_cache=dlog(V_seed)
    ic_v_seed_mid_cache=dexp(0.5d0*(ic_x_seed_cache(1:Num_nu-1)+ic_x_seed_cache(2:Num_nu)))
    ic_d_nu_cache=ic_v_seed_mid_cache*(ic_x_seed_cache(2:Num_nu)-ic_x_seed_cache(1:Num_nu-1))
    ic_gam_e_mean_cache=(gam_e(1:Num_gam_e-1)+gam_e(2:Num_gam_e))/two
    ic_e_seed_cache=ic_v_seed_mid_cache*para_hEme
    ic_num_gam_cache=Num_gam_e
    ic_num_nu_cache=Num_nu
    ic_grid_cache_ready=.true.
end subroutine rebuild_ic_grid_cache
end subroutine ensure_ic_grid_cache

! 数值计算逆康普顿（IC）冷却率：双重积分（种子光子×散射截面），含Jones/Blumenthal核。
subroutine electron_cooling_ic_loss(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn, dot_gam_e)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dot_gam_e(Num_gam_e)
real(8) :: photon_number(Num_nu-1)

    call ensure_ic_grid_cache(Num_gam_e,Num_nu,gam_e,V_seed)
    dot_gam_e=zero

    do I_nu=1,Num_nu-1
       photon_number(I_nu)=electron_powerlaw_interp(V_seed(I_nu),V_seed(I_nu+1), &
                                                    Seed_syn(I_nu),Seed_syn(I_nu+1),ic_v_seed_mid_cache(I_nu))
    end do

    do i_gam_e=1,Num_gam_e-1
       call accumulate_ic_gamma_loss(i_gam_e,dot_gam_e(i_gam_e))
    end do

    dot_gam_e=dot_gam_e/gam_e/gam_e*para_h*Para_h*Para_SigmaT/para_m_energy
    dot_gam_e(Num_gam_e)=0.99*dot_gam_e(Num_gam_e-1)

contains

subroutine accumulate_ic_gamma_loss(i_gam_e,dot_val)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: i_gam_e
real(8), intent(out) :: dot_val
integer :: I_nu,Nu_s
real(8) :: game,game_pow,var,dInteg2,V_t,E_t2eV,Vloc,E_Vloc2eV,uplim,temp,q,fssc,kn_factor

    dot_val=zero
    game=ic_gam_e_mean_cache(i_gam_e)
    game_pow=game*game
    var=0.25d0/game_pow
    do I_nu=1,Num_nu-1
       dInteg2=zero
       V_t=ic_v_seed_mid_cache(I_nu)
       E_t2eV=ic_e_seed_cache(I_nu)
       kn_factor=4d0*game*E_t2eV
       uplim=(4d0*game_pow*E_t2eV)/(one+kn_factor)
       do Nu_s=1,Num_nu-1
          fssc=zero
          Vloc=ic_v_seed_mid_cache(Nu_s)
          E_Vloc2eV=ic_e_seed_cache(Nu_s)
          if (Vloc > var*V_t .and. Vloc <= V_t) then
             fssc=Vloc/V_t-var
          else
             if (E_Vloc2eV > uplim) exit
             temp=game-E_Vloc2eV
             if (temp <= zero) exit
             q=E_Vloc2eV/(kn_factor*temp)
             if (q <= zero) cycle
             if (q >= one) exit
             fssc=two*q*(log(q)-q)+one+q+ &
                  0.5d0*(one-q)*(4d0*game*E_t2eV*q)**2/(1+4d0*game*q*E_t2eV)
          end if
          dInteg2=dInteg2+Vloc*fssc*ic_d_nu_cache(Nu_s)
       end do
       dot_val=dot_val+photon_number(I_nu)/V_t*ic_d_nu_cache(I_nu)*dInteg2
    end do
end subroutine accumulate_ic_gamma_loss
end subroutine electron_cooling_ic_loss

! IC cooling from the same Jones/KN emissivity kernel used by radiation_ssc_spectrum.
subroutine electron_cooling_ic_loss_emissivity_budget(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn,dot_gam_over_gam)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dot_gam_over_gam(Num_gam_e)
real(8), allocatable :: seed_weights(:),obs_weights(:),e_seed(:),inv_v_seed(:)
integer :: i_gam,work_items,thread_count
logical :: use_parallel

    allocate(seed_weights(Num_nu),obs_weights(Num_nu),e_seed(Num_nu),inv_v_seed(Num_nu))
    call compute_simpson_weights(seed_weights,Num_nu)
    call compute_simpson_weights(obs_weights,Num_nu)
    h_nu=dlog(V_seed(2))-dlog(V_seed(1))
    h_nu_third=h_nu/3.0d0
    para_hEme=Para_h/para_m_energy
    e_seed=V_seed*para_hEme
    inv_v_seed=one/V_seed
    temp_norm_ic=0.75d0*Para_c*Para_h*Para_SigmaT/Para_m_energy
    dot_gam_over_gam=zero

    work_items=Num_gam_e*Num_nu*Num_nu
    thread_count=max(1,n_threads)
    use_parallel=(n_threads > 1 .and. work_items >= 8192)
    !$OMP PARALLEL DO if(use_parallel) num_threads(thread_count) schedule(static) &
    !$OMP& private(i_gam)
    do i_gam=1,Num_gam_e
        call accumulate_budget_gamma(i_gam,dot_gam_over_gam(i_gam))
    end do
    !$OMP END PARALLEL DO

    deallocate(seed_weights,obs_weights,e_seed,inv_v_seed)

contains

subroutine accumulate_budget_gamma(i_gam,loss_val)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: i_gam
real(8), intent(out) :: loss_val
integer :: i_obs,i_seed
real(8) :: gam,gam2,seed_sum,power_log

    gam=gam_e(i_gam)
    gam2=gam*gam
    loss_val=zero
    do i_obs=1,Num_nu
        if (gam_e(Num_gam_e) <= e_seed(i_obs)) cycle
        seed_sum=zero
        do i_seed=1,Num_nu
            if (Seed_syn(i_seed) <= zero) cycle
            if (i_seed < i_obs) then
                seed_sum=seed_sum+seed_weights(i_seed)*Seed_syn(i_seed)*low_seed_kernel(gam,i_obs,i_seed)/gam2
            else
                seed_sum=seed_sum+seed_weights(i_seed)*Seed_syn(i_seed) * &
                         max(zero,V_seed(i_obs)*inv_v_seed(i_seed)-0.25d0/gam2)/gam2
            end if
        end do
        power_log=temp_norm_ic*V_seed(i_obs)*h_nu_third*seed_sum
        loss_val=loss_val+obs_weights(i_obs)*power_log
    end do
    loss_val=h_nu_third*loss_val/gam
end subroutine accumulate_budget_gamma

real(8) function low_seed_kernel(gam,i_obs,i_seed)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: gam
integer, intent(in) :: i_obs,i_seed
real(8) :: temp,q,log_q,q_gamma,kn_coeff

    low_seed_kernel=zero
    temp=gam-e_seed(i_obs)
    if (temp <= zero) return
    q=V_seed(i_obs)/(4.0d0*gam*temp*V_seed(i_seed))
    if (q <= zero .or. q >= one) return
    log_q=dlog(q)
    q_gamma=e_seed(i_obs)/temp
    kn_coeff=q_gamma*q_gamma/(two*(one+q_gamma))
    low_seed_kernel=two*q*(log_q-q)+one+q+kn_coeff*(one-q)
end function low_seed_kernel
end subroutine electron_cooling_ic_loss_emissivity_budget

end module electron_cooling_ic_kernel
