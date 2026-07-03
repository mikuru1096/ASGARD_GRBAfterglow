!f2py: skip
module electron_ssa_kernel
  use constants
  use electron_radiation_kernel, only: greater_window, pl_interp, &
                                       log_gauss2
  private

  public :: electron_ssa_loss

  integer, parameter :: lowseg=1,highseg=2
  integer, save :: nnu_cache=0
  logical, save :: seed_ready=.false.
  real(8), allocatable, save, dimension(:) :: v_cache,vlow_cache,vhigh_cache,vg1_cache,vg2_cache,wg1_cache,wg2_cache
  !$omp threadprivate(nnu_cache,seed_ready,v_cache,vlow_cache,vhigh_cache)
  !$omp threadprivate(vg1_cache,vg2_cache,wg1_cache,wg2_cache)

contains
! 刷新 SSA 种子频率缓存；线程私有，避免 OpenMP 列计算共享临时状态。
! Refresh the SSA seed-frequency cache; it is thread-private for OpenMP column work.
subroutine ensure_ssa_cache(nnu,vseed)
implicit none
integer, intent(in) :: nnu
real(8), intent(in), dimension(nnu) :: vseed
integer :: inu
logical :: match

    match=.false.
    if (seed_ready) then
        if (nnu_cache == nnu) then
            if (allocated(v_cache)) then
                match=all(v_cache == vseed)
            end if
        end if
    end if

    if (match) return

    if (allocated(v_cache)) deallocate(v_cache,vlow_cache,vhigh_cache,vg1_cache,vg2_cache,wg1_cache,wg2_cache)
    allocate(v_cache(nnu),vlow_cache(nnu-1),vhigh_cache(nnu-1),vg1_cache(nnu-1),vg2_cache(nnu-1), &
             wg1_cache(nnu-1),wg2_cache(nnu-1))
    v_cache=vseed
    vlow_cache=vseed(1:nnu-1)
    vhigh_cache=vseed(2:nnu)
    do inu=1,nnu-1
        call log_gauss2(vlow_cache(inu),vhigh_cache(inu),vg1_cache(inu),vg2_cache(inu), &
                                          wg1_cache(inu),wg2_cache(inu))
    end do
    nnu_cache=nnu
    seed_ready=.true.
end subroutine ensure_ssa_cache

! 把 SSA 种子频率游标移到第一个超过低频阈值的网格点。
! Move the SSA seed cursor to the first grid point above the low-frequency bound.
subroutine advance_ssa_cursor(nnu,vmin,cursor)
implicit none
integer, intent(in) :: nnu
real(8), intent(in) :: vmin
integer, intent(inout) :: cursor

    if (cursor < 1) cursor=1
    if (cursor > nnu+1) cursor=nnu+1

    do while (cursor > 1)
        if (v_cache(cursor-1) <= vmin) exit
        cursor=cursor-1
    end do

    do while (cursor <= nnu)
        if (v_cache(cursor) > vmin) exit
        cursor=cursor+1
    end do
end subroutine advance_ssa_cursor

! 为每个电子能格预计算 SSA 积分段和截面前因子。
! Precompute SSA integration segments and cross-section prefactors for each electron bin.
subroutine build_ssa_geometry(db,ng,nnu,gam,lowpos,lowfirst,lowlast,highfirst,siglow,sighigh,cyclonu)
implicit none
integer, intent(in) :: ng,nnu
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in) :: db
integer, intent(out), dimension(ng) :: lowpos,lowfirst,lowlast,highfirst
real(8), intent(out), dimension(ng) :: siglow,sighigh
real(8), intent(out) :: cyclonu
integer :: ig,cursor,upper
real(8) :: bcr,pref1,pref2,ge,g2,g3,vmin,vmax

    bcr=4.4d13
    pref1=2.5042d-22*bcr/db
    pref2=7.787d-22*bcr/db
    cyclonu=para_e*db/(2d0*pi*para_m_e*para_c)

    cursor=1
    do ig=1,ng
       ge=gam(ig)
       g2=ge*ge
       g3=g2*ge
       vmin=cyclonu/ge
       vmax=1.5d0*g2*cyclonu

       call advance_ssa_cursor(nnu,vmin,cursor)
       lowpos(ig)=cursor

       if (cursor <= nnu) then
          call greater_window(v_cache,nnu,cursor,vmax,upper)
          siglow(ig)=pref1*(3d0*vmin)**(5d0/3d0)
          sighigh(ig)=pref2/g3
          lowfirst(ig)=max(1,cursor-1)
          lowlast(ig)=min(nnu-1,upper-1)
          highfirst(ig)=max(1,upper-1)
       else
          siglow(ig)=0d0
          sighigh(ig)=0d0
          lowfirst(ig)=1
          lowlast(ig)=0
          highfirst(ig)=nnu
       end if
    end do
end subroutine build_ssa_geometry

! SSA 冷却率：同一电子能格上批量处理多个 chi 列的种子光子场。
! SSA cooling rate: process several chi-column seed photon fields on the same electron grid.
subroutine electron_ssa_loss(db,ng,nnu,nchi,nthr,gam,vseed,seed,loss)
!$ use omp_lib
implicit none
integer, intent(in) :: ng,nnu,nchi,nthr
real(8), intent(in), dimension(ng) :: gam
real(8), intent(in), dimension(nnu) :: vseed
real(8), intent(in), dimension(nnu,nchi) :: seed
real(8), intent(in) :: db
real(8), intent(out), dimension(ng,nchi) :: loss
integer, parameter :: work_min=512
integer, dimension(ng) :: lowpos,lowfirst,lowlast,highfirst
integer :: inu,ic,ig,work,nt
logical :: doomp
real(8), dimension(nnu-1) :: vlow,vhigh,vg1,vg2,wg1,wg2
real(8), dimension(ng) :: siglow,sighigh
real(8), dimension(nnu-1,nchi) :: seedg1,seedg2
real(8) :: cyclonu
real(8), dimension(0:nnu-1,nchi) :: lowpref
real(8), dimension(nnu-1,nchi) :: high1,high2

    call ensure_ssa_cache(nnu,vseed)
    vlow=vlow_cache
    vhigh=vhigh_cache
    vg1=vg1_cache
    vg2=vg2_cache
    wg1=wg1_cache
    wg2=wg2_cache
    lowpref(0,:)=0d0
    do ic=1,nchi
        do inu=1,nnu-1
            seedg1(inu,ic)=pl_interp(vlow(inu),vhigh(inu),seed(inu,ic),seed(inu+1,ic),vg1(inu))
            seedg2(inu,ic)=pl_interp(vlow(inu),vhigh(inu),seed(inu,ic),seed(inu+1,ic),vg2(inu))
            lowpref(inu,ic)=lowpref(inu-1,ic)+para_h*para_c*(wg1(inu)*seedg1(inu,ic)*vg1(inu)**(-2d0/3d0)+ &
                              wg2(inu)*seedg2(inu,ic)*vg2(inu)**(-2d0/3d0))
            high1(inu,ic)=wg1(inu)*seedg1(inu,ic)
            high2(inu,ic)=wg2(inu)*seedg2(inu,ic)
        end do
    end do
    call build_ssa_geometry(db,ng,nnu,gam,lowpos,lowfirst,lowlast,highfirst,siglow,sighigh,cyclonu)

    loss=0d0
    work=ng*nchi*nnu
    nt=max(1,nthr)
    doomp=(nthr > 1 .and. work >= work_min)
    !$OMP PARALLEL DO collapse(2) if(doomp) num_threads(nt) schedule(static) &
    !$OMP& private(ic,ig)
    do ic=1,nchi
        do ig=1,ng
            call accumulate_ssa_cell(ig,ic,loss(ig,ic))
        end do
    end do
    !$OMP END PARALLEL DO

contains

subroutine accumulate_ssa_cell(ig,ic,rate)
implicit none
integer, intent(in) :: ig,ic
real(8), intent(out) :: rate
integer :: inu,lfirst,llast,hfirst
real(8) :: ge,g2,vmin,vmax,invmax,highpref,total,vlo,vhi

    rate=0d0
    if (lowpos(ig) > nnu) return

    ge=gam(ig)
    g2=ge*ge
    vmin=cyclonu/ge
    vmax=1.5d0*g2*cyclonu
    invmax=1d0/vmax
    highpref=sighigh(ig)*cyclonu*para_h*para_c
    total=0d0

    lfirst=lowfirst(ig)
    llast=lowlast(ig)
    if (lfirst <= llast) then
        vlo=max(vlow(lfirst),vmin)
        vhi=min(vhigh(lfirst),vmax)
        if (vhi > vlo) then
            if (vlo /= vlow(lfirst) .or. vhi /= vhigh(lfirst)) then
                total=total+clipped_ssa_segment(vlo,vhi,lfirst,ic,siglow(ig),lowseg,cyclonu,vmax)
                lfirst=lfirst+1
            end if
        else
            lfirst=lfirst+1
        end if
    end if
    if (llast >= lfirst) then
        vlo=max(vlow(llast),vmin)
        vhi=min(vhigh(llast),vmax)
        if (vhi > vlo) then
            if (vlo /= vlow(llast) .or. vhi /= vhigh(llast)) then
                total=total+clipped_ssa_segment(vlo,vhi,llast,ic,siglow(ig),lowseg,cyclonu,vmax)
                llast=llast-1
            end if
        else
            llast=llast-1
        end if
    end if
    if (llast >= lfirst) then
        total=total+siglow(ig)*(lowpref(llast,ic)-lowpref(lfirst-1,ic))
    end if

    hfirst=highfirst(ig)
    if (hfirst <= nnu-1) then
        vlo=max(vlow(hfirst),vmax)
        vhi=vhigh(hfirst)
        if (vhi > vlo) then
            if (vlo /= vlow(hfirst)) then
                total=total+clipped_ssa_segment(vlo,vhi,hfirst,ic,sighigh(ig),highseg,cyclonu,vmax)
                hfirst=hfirst+1
            end if
        end if
    end if
    do inu=hfirst,nnu-1
        total=total+highpref*(high1(inu,ic)*dexp(-vg1(inu)*invmax)+high2(inu,ic)*dexp(-vg2(inu)*invmax))
    end do

    rate=total
end subroutine accumulate_ssa_cell

real(8) function clipped_ssa_segment(vlo,vhi,inu,ic,sigpref,mode,cyclonu,vmax)
implicit none
integer, intent(in) :: inu,ic,mode
real(8), intent(in) :: vlo,vhi,sigpref,cyclonu,vmax
integer :: iq
real(8), dimension(2) :: vg,wg
real(8) :: seedval,sigval

    call log_gauss2(vlo,vhi,vg(1),vg(2),wg(1),wg(2))
    clipped_ssa_segment=0d0

    do iq=1,2
        seedval=pl_interp(vlow(inu),vhigh(inu),seed(inu,ic),seed(inu+1,ic),vg(iq))
        if (mode == lowseg) then
            sigval=sigpref*vg(iq)**(-5d0/3d0)
        else
            sigval=sigpref*(cyclonu/vg(iq))*dexp(-vg(iq)/vmax)
        end if
        clipped_ssa_segment=clipped_ssa_segment+wg(iq)*sigval*seedval*para_h*vg(iq)*para_c
    end do
end function clipped_ssa_segment
end subroutine electron_ssa_loss

end module electron_ssa_kernel
