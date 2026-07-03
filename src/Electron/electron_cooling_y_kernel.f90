!f2py: skip
module electron_y_kernel
  use constants
  use electron_radiation_kernel, only: greater_window, pl_interp, &
                                       pl_integral, log_gauss2
  private

  public :: electron_y_nakar, electron_y_fan

  integer, save :: nakar_ng=0,nakar_nnu=0
  integer, allocatable, save, dimension(:) :: idx_cache
  real(8), allocatable, save, dimension(:) :: hat_cache,prefix_cache,v_cache,vloc_cache,vg1_cache,vg2_cache &
      & ,wg1_cache,wg2_cache

contains
! 准备 Nakar Y 积分缓存：hat_nu、频段 Gauss 节点和查找区间。
! Prepare Nakar-Y integration cache: hat_nu, frequency-bin Gauss nodes, and lookup windows.
subroutine ensure_y_nakar(ng,nnu,gam,vseed)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: ng,nnu
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in), dimension(nnu) :: vseed
logical :: rebuild

    rebuild=.not. allocated(hat_cache)
    if (.not. rebuild) rebuild = (nakar_ng /= ng)
    if (.not. rebuild) rebuild = (nakar_nnu /= nnu)
    if (.not. rebuild) rebuild = any(hat_cache /= Para_m_energy/Para_h/gam)
    if (.not. rebuild) rebuild = any(v_cache /= vseed)
    if (.not. rebuild) return

    if (allocated(hat_cache)) deallocate(hat_cache,prefix_cache,v_cache,vloc_cache,idx_cache,vg1_cache,vg2_cache, &
                                         wg1_cache,wg2_cache)
    allocate(hat_cache(ng),prefix_cache(nnu),v_cache(nnu),vloc_cache(ng),idx_cache(ng),vg1_cache(nnu-1), &
             vg2_cache(nnu-1),wg1_cache(nnu-1),wg2_cache(nnu-1))
    hat_cache=Para_m_energy/Para_h/gam
    v_cache=vseed
    do inu=1,nnu-1
        call log_gauss2(vseed(inu),vseed(inu+1),vg1_cache(inu),vg2_cache(inu), &
                                          wg1_cache(inu),wg2_cache(inu))
    end do
    do icomp=1,ng
        if (hat_cache(icomp) <= vseed(1)) then
            idx_cache(icomp)=0
            vloc_cache(icomp)=vseed(1)
        else
            call greater_window(vseed,nnu,2,hat_cache(icomp),inu)
            idx_cache(icomp)=inu
            if (inu <= nnu) vloc_cache(icomp)=min(hat_cache(icomp),vseed(inu))
        end if
    end do
    nakar_ng=ng
    nakar_nnu=nnu
end subroutine ensure_y_nakar

! Nakar+2009 Y 参数：从 hat_nu(gamma) 到最高频积分同步光谱。
! Nakar+2009 Y parameter: integrate the synchrotron spectrum above hat_nu(gamma).
subroutine electron_y_nakar(ng,nnu,nthr,gam,vseed,psyn,comp)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: ng,nnu,nthr
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in), dimension(nnu) :: vseed,psyn
real(8), intent(out), dimension(ng) :: comp
integer :: icomp,inu

    call ensure_y_nakar(ng,nnu,gam,vseed)

    comp=0d0
    prefix_cache(1)=0d0
    do inu=2,nnu
       prefix_cache(inu)=prefix_cache(inu-1)+ &
            wg1_cache(inu-1)*pl_interp(vseed(inu-1),vseed(inu),psyn(inu-1),psyn(inu),vg1_cache(inu-1))+ &
            wg2_cache(inu-1)*pl_interp(vseed(inu-1),vseed(inu),psyn(inu-1),psyn(inu),vg2_cache(inu-1))
    end do

    do icomp=1,ng
       comp(icomp)=0d0
       inu=idx_cache(icomp)
       if (inu == 0) cycle
       if (inu <= nnu) then
          comp(icomp)=prefix_cache(inu-1)+pl_integral(vseed(inu-1),vloc_cache(icomp), &
                      psyn(inu-1),pl_interp(vseed(inu-1),vseed(inu),psyn(inu-1),psyn(inu), &
                      vloc_cache(icomp)))
       else
          comp(icomp)=prefix_cache(nnu)
       end if
    end do
end subroutine electron_y_nakar

! Fan+2008 Y 参数：解析分段 eta_KN(gamma) x eta_rad，区分快/慢冷却。
! Fan+2008 Y parameter: analytic eta_KN(gamma) x eta_rad branches for fast/slow cooling.
subroutine electron_y_fan(ee,eb,p,db,gm,gc,gmax,ng,gam,comp)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: ng
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in) :: ee,eb,p,db,gm,gc,gmax
real(8), intent(out), dimension(ng) :: comp
integer :: ig

    eta=(gm/gc)**(p-2d0)
    if (eta-1d0 > 0.001) eta=1d0

    do ig=1,ng-1
        hatgam=5.4246d6/sqrt(db*gam(ig+1))
        if (gm > gc) then
            if (hatgam < gc) then
                etakn=0d0
            else if (hatgam < gm) then
                if (p>2) then
                    Step1=(p-1)/(p-2)*gm-gc
                    etakn=(hatgam-gc)/Step1
                else
                    Step1=gm**(p-1)*gmax**(2-p)-(p-1)*gm-(2-p)*gc
                    etakn=(2-p)*(hatgam-gc)/Step1
                end if
            else
                if (p>2) then
                    Step2=gm**(p-1)*hatgam**(2-p)
                    Step3=(p-1)*gm-(p-2)*gc
                    etakn=1-Step2/Step3
                else
                    Step2=gm**(p-1)*gmax**(2-p)-(p-1)*gm-(2-p)*gc
                    Step3=gm**(p-1)*(gmax**(2-p)-hatgam**(2-p))
                    etakn=1-Step2/Step3
                end if
            end if
        else if (hatgam < gm) then
            etakn=0d0
        else if (hatgam < gc) then
            if (p>2) then
                Step4=gc**(3-p)/(p-2d0)-gm**(3-p)
                etakn=(hatgam**(3-p)-gm**(3-p))/Step4
            else
                Step4=(3-p)*gc*gmax**(2-p)-gc**(3-p)-(2-p)*gm**(3-p)
                etakn=(2-p)*(hatgam**(3-p)-gm**(3-p))/Step4
            end if
        else
            if (p>2) then
                Step5=(3-p)*gc*hatgam**(2-p)
                Step6=gc**(3d0-p)-(p-2)*gm**(3d0-p)
                etakn=1-Step5/Step6
            else
                Step5=(3-p)*gc*(gmax**(2-p)-hatgam**(2-p))
                Step6=(3-p)*gc*gmax**(2-p)-gc**(3-p)-(2-p)*gm**(3-p)
                etakn=1-Step5/Step6
            end if
        end if
        comp(ig)=0.5d0*(-1d0+sqrt(1d0+4d0*eta*etakn*ee/eb))
    end do
    comp(ng)=0.99*comp(ng-1)
end subroutine electron_y_fan

end module electron_y_kernel
