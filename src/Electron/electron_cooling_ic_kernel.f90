!f2py: skip
module electron_ic_kernel
  use constants
  use rad_common, only: compute_simpson_weights
  use electron_radiation_kernel, only: pl_interp
  private

  public :: electron_ic_loss, electron_ic_budget

  integer, save :: ng_cache=0,nnu_cache=0
  logical, save :: grid_ready=.false.
  real(8), allocatable, save, dimension(:) :: dnu_cache,gmean_cache,eseed_cache,xseed_cache,vmid_cache

contains
! 更新 IC 积分网格缓存；缓存只依赖电子能格和种子光子频格。
! Refresh the IC quadrature cache; it depends only on electron and seed-photon grids.
subroutine ensure_ic_grid(ng,nnu,gam,vseed)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: ng,nnu
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in), dimension(nnu) :: vseed

    if (ic_grid_current()) return
    call rebuild_ic_grid()

contains
logical function ic_grid_current()
implicit REAL(8)(A-H,O-Z)

    ic_grid_current=.false.
    if (.not. grid_ready) return
    if (ng_cache /= ng) return
    if (nnu_cache /= nnu) return
    if (.not. ic_seed_current()) return
    if (.not. ic_gamma_current()) return
    ic_grid_current=.true.
end function ic_grid_current

logical function ic_seed_current()
implicit REAL(8)(A-H,O-Z)
integer :: inu

    ic_seed_current=.false.
    if (.not. allocated(xseed_cache)) return
    do inu=1,nnu
        if (xseed_cache(inu) /= dlog(vseed(inu))) return
    end do
    ic_seed_current=.true.
end function ic_seed_current

logical function ic_gamma_current()
implicit REAL(8)(A-H,O-Z)
integer :: ig

    ic_gamma_current=.false.
    if (.not. allocated(gmean_cache)) return
    do ig=1,ng-1
        if (gmean_cache(ig) /= (gam(ig)+gam(ig+1))/2d0) return
    end do
    ic_gamma_current=.true.
end function ic_gamma_current

subroutine rebuild_ic_grid()
implicit REAL(8)(A-H,O-Z)

    if (allocated(dnu_cache)) deallocate(dnu_cache,gmean_cache,eseed_cache,xseed_cache,vmid_cache)
    allocate(dnu_cache(nnu-1),gmean_cache(ng-1),eseed_cache(nnu-1),xseed_cache(nnu),vmid_cache(nnu-1))

    para_hEme=Para_h/para_m_energy
    xseed_cache=dlog(vseed)
    vmid_cache=dexp(0.5d0*(xseed_cache(1:nnu-1)+xseed_cache(2:nnu)))
    dnu_cache=vmid_cache*(xseed_cache(2:nnu)-xseed_cache(1:nnu-1))
    gmean_cache=(gam(1:ng-1)+gam(2:ng))/2d0
    eseed_cache=vmid_cache*para_hEme
    ng_cache=ng
    nnu_cache=nnu
    grid_ready=.true.
end subroutine rebuild_ic_grid
end subroutine ensure_ic_grid

! Jones/Blumenthal IC 冷却率，沿种子光子频率和散射后频率做双重积分。
! Jones/Blumenthal IC cooling rate from a double integral over seed and scattered frequencies.
subroutine electron_ic_loss(ng,nnu,nthr,gam,vseed,seed,loss)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: ng,nnu,nthr
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in), dimension(nnu) :: vseed,seed
real(8), intent(out), dimension(ng) :: loss
real(8), dimension(nnu-1) :: photons

    call ensure_ic_grid(ng,nnu,gam,vseed)
    loss=0d0

    do inu=1,nnu-1
       photons(inu)=pl_interp(vseed(inu),vseed(inu+1),seed(inu),seed(inu+1),vmid_cache(inu))
    end do

    do ig=1,ng-1
       call accumulate_ic_loss(ig,loss(ig))
    end do

    loss=loss/gam/gam*para_h*Para_h*Para_SigmaT/para_m_energy
    loss(ng)=0.99*loss(ng-1)

contains

subroutine accumulate_ic_loss(ig,rate)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: ig
real(8), intent(out) :: rate
integer :: inu,is
real(8) :: gmid,g2,var,integ,vt,et,vloc,ev,uplim,temp,q,fssc,kn

    rate=0d0
    gmid=gmean_cache(ig)
    g2=gmid*gmid
    var=0.25d0/g2
    do inu=1,nnu-1
       integ=0d0
       vt=vmid_cache(inu)
       et=eseed_cache(inu)
       kn=4d0*gmid*et
       uplim=(4d0*g2*et)/(1d0+kn)
       do is=1,nnu-1
          fssc=0d0
          vloc=vmid_cache(is)
          ev=eseed_cache(is)
          if (vloc > var*vt .and. vloc <= vt) then
             fssc=vloc/vt-var
          else
             if (ev > uplim) exit
             temp=gmid-ev
             if (temp <= 0d0) exit
             q=ev/(kn*temp)
             if (q <= 0d0) cycle
             if (q >= 1d0) exit
             fssc=2d0*q*(log(q)-q)+1d0+q+0.5d0*(1d0-q)*(4d0*gmid*et*q)**2/(1d0+4d0*gmid*q*et)
          end if
          integ=integ+vloc*fssc*dnu_cache(is)
       end do
       rate=rate+photons(inu)/vt*dnu_cache(inu)*integ
    end do
end subroutine accumulate_ic_loss
end subroutine electron_ic_loss

! 用 SSC 辐射同源的 Jones/KN 发射率核计算 IC 能量预算。
! Compute the IC energy budget with the same Jones/KN emissivity kernel used by SSC radiation.
subroutine electron_ic_budget(ng,nnu,nthr,gam,vseed,seed,loss)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: ng,nnu,nthr
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in), dimension(nnu) :: vseed,seed
real(8), intent(out), dimension(ng) :: loss
real(8), allocatable, dimension(:) :: wseed,wobs,eseed,inverseed
integer :: ig,work,nt
logical :: doomp

    allocate(wseed(nnu),wobs(nnu),eseed(nnu),inverseed(nnu))
    call compute_simpson_weights(wseed,nnu)
    call compute_simpson_weights(wobs,nnu)
    hstep=dlog(vseed(2))-dlog(vseed(1))
    hthird=hstep/3d0
    para_hEme=Para_h/para_m_energy
    eseed=vseed*para_hEme
    inverseed=1d0/vseed
    cnorm=0.75d0*Para_c*Para_h*Para_SigmaT/Para_m_energy
    loss=0d0

    work=ng*nnu*nnu
    nt=max(1,nthr)
    doomp=(nthr > 1 .and. work >= 8192)
    !$OMP PARALLEL DO if(doomp) num_threads(nt) schedule(static) &
    !$OMP& private(ig)
    do ig=1,ng
        call accumulate_budget(ig,loss(ig))
    end do
    !$OMP END PARALLEL DO

    deallocate(wseed,wobs,eseed,inverseed)

contains

subroutine accumulate_budget(ig,rate)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: ig
real(8), intent(out) :: rate
integer :: iobs,iseed
real(8) :: ge,g2,seedsum,power

    ge=gam(ig)
    g2=ge*ge
    rate=0d0
    do iobs=1,nnu
        if (gam(ng) <= eseed(iobs)) cycle
        seedsum=0d0
        do iseed=1,nnu
            if (seed(iseed) <= 0d0) cycle
            if (iseed < iobs) then
                seedsum=seedsum+wseed(iseed)*seed(iseed)*low_seed_kernel(ge,iobs,iseed)/g2
            else
                seedsum=seedsum+wseed(iseed)*seed(iseed)*max(0d0,vseed(iobs)*inverseed(iseed)-0.25d0/g2)/g2
            end if
        end do
        power=cnorm*vseed(iobs)*hthird*seedsum
        rate=rate+wobs(iobs)*power
    end do
    rate=hthird*rate/ge
end subroutine accumulate_budget

real(8) function low_seed_kernel(gam,i_obs,i_seed)
implicit REAL(8)(A-H,O-Z)
real(8), intent(in) :: gam
integer, intent(in) :: i_obs,i_seed
real(8) :: temp,q,logq,qgam,kncoeff

    low_seed_kernel=0d0
    temp=gam-eseed(i_obs)
    if (temp <= 0d0) return
    q=vseed(i_obs)/(4d0*gam*temp*vseed(i_seed))
    if (q <= 0d0 .or. q >= 1d0) return
    logq=dlog(q)
    qgam=eseed(i_obs)/temp
    kncoeff=qgam*qgam/(2d0*(1d0+qgam))
    low_seed_kernel=2d0*q*(logq-q)+1d0+q+kncoeff*(1d0-q)
end function low_seed_kernel
end subroutine electron_ic_budget

end module electron_ic_kernel
