!f2py: skip
module hadronic_pg
    use constants
    use hadronic_base, only: proton_m
    implicit none
    private
    public :: pg_hummer

    real(8), parameter :: mbarn_cm2 = 1.0d-30

    real(8), parameter, dimension(28) :: it1_er= [ 0.34d0, 0.64d0, 0.64d0, 0.64d0, 0.75d0, 0.75d0, 0.75d0, 0.75d0, 0.77d0, 1.03d0, 1.03d0, 1.03d0, 1.04d0, 1.04d0, 1.04d0, 1.04d0, 1.05d0, 1.05d0, 1.05d0, 1.05d0, 1.45d0, 1.45d0, 1.45d0, 1.45d0, 1.56d0, 1.56d0, 1.56d0, 1.56d0 ]
    real(8), parameter, dimension(28) :: it1_bout= [ 1.00d0, 0.67d0, 0.33d0, 0.33d0, 0.52d0, 0.42d0, 0.42d0, 0.06d0, 0.45d0, 0.75d0, 0.14d0, 0.14d0, 0.64d0, 0.22d0, 0.22d0, 0.14d0, 0.14d0, 0.55d0, 0.55d0, 0.31d0, 0.14d0, 0.13d0, 0.13d0, 0.73d0, 0.37d0, 0.39d0, 0.39d0, 0.24d0 ]
    real(8), parameter, dimension(28) :: it1_g= [ 66.69d0, 4.260d0, 4.260d0, 4.260d0, 27.01d0, 27.01d0, 27.01d0, 27.01d0, 6.590d0, 3.190d0, 3.190d0, 3.190d0, 15.72d0, 15.72d0, 15.72d0, 15.72d0, 21.68d0, 21.68d0, 21.68d0, 21.68d0, 3.030d0, 3.030d0, 3.030d0, 3.030d0, 16.46d0, 16.46d0, 16.46d0, 16.46d0 ]
    real(8), parameter, dimension(28) :: it1_chi= [ 0.22d0, 0.29d0, 0.14d0, 0.19d0, 0.31d0, 0.17d0, 0.18d0, 0.22d0, 0.32d0, 0.35d0, 0.23d0, 0.17d0, 0.35d0, 0.24d0, 0.17d0, 0.23d0, 0.35d0, 0.24d0, 0.16d0, 0.23d0, 0.38d0, 0.29d0, 0.15d0, 0.23d0, 0.39d0, 0.30d0, 0.15d0, 0.23d0 ]
    real(8), parameter, dimension(28) :: it1_chipn= [ 0.78d0, 0.71d0, 1.00d0, 0.67d0, 0.69d0, 1.00d0, 0.65d0, 0.56d0, 0.68d0, 0.65d0, 1.00d0, 0.60d0, 0.65d0, 1.00d0, 0.59d0, 0.54d0, 0.65d0, 1.00d0, 0.60d0, 0.54d0, 0.62d0, 1.00d0, 0.56d0, 0.54d0, 0.61d0, 1.00d0, 0.55d0, 0.54d0 ]
    real(8), parameter, dimension(28) :: it1_pi0= [ 0.667d0, 0.333d0, 0.333d0, 0.333d0, 0.333d0, 0.333d0, 0.333d0, 0.667d0, 0.333d0, 0.333d0, 0.333d0, 0.333d0, 0.333d0, 0.333d0, 0.333d0, 0.667d0, 0.667d0, 0.067d0, 0.400d0, 0.333d0, 0.667d0, 0.067d0, 0.400d0, 0.333d0, 0.667d0, 0.067d0, 0.400d0, 0.333d0 ]
    real(8), parameter, dimension(28) :: it1_pip= [ 0.333d0, 0.667d0, 0.167d0, 0.611d0, 0.667d0, 0.167d0, 0.611d0, 1.000d0, 0.667d0, 0.667d0, 0.167d0, 0.611d0, 0.667d0, 0.167d0, 0.611d0, 1.000d0, 0.333d0, 0.533d0, 0.422d0, 1.000d0, 0.333d0, 0.533d0, 0.422d0, 1.000d0, 0.333d0, 0.533d0, 0.422d0, 1.000d0 ]
    real(8), parameter, dimension(28) :: it1_pim= [ 0.000d0, 0.000d0, 0.500d0, 0.056d0, 0.000d0, 0.500d0, 0.056d0, 0.333d0, 0.000d0, 0.000d0, 0.500d0, 0.056d0, 0.000d0, 0.500d0, 0.056d0, 0.333d0, 0.000d0, 0.400d0, 0.178d0, 0.667d0, 0.000d0, 0.400d0, 0.178d0, 0.667d0, 0.000d0, 0.400d0, 0.178d0, 0.667d0 ]
    real(8), parameter, dimension(28) :: it1_pro= [ 0.667d0, 0.333d0, 0.000d0, 0.778d0, 0.333d0, 0.000d0, 0.778d0, 0.333d0, 0.333d0, 0.333d0, 0.000d0, 0.778d0, 0.333d0, 0.000d0, 0.778d0, 0.333d0, 0.667d0, 0.000d0, 0.622d0, 0.667d0, 0.667d0, 0.000d0, 0.622d0, 0.667d0, 0.667d0, 0.000d0, 0.622d0, 0.667d0 ]
    real(8), parameter, dimension(28) :: it1_ntr= [ 0.333d0, 0.667d0, 0.000d0, 0.222d0, 0.667d0, 0.000d0, 0.222d0, 0.667d0, 0.667d0, 0.667d0, 0.000d0, 0.222d0, 0.667d0, 0.000d0, 0.222d0, 0.667d0, 0.333d0, 0.000d0, 0.378d0, 0.333d0, 0.333d0, 0.000d0, 0.378d0, 0.333d0, 0.333d0, 0.000d0, 0.378d0, 0.333d0 ]
    real(8), parameter, dimension(9) :: it2_em= [ 0.17d0, 0.56d0, 10.0d0, 0.4d0, 1.58d0, 10.0d0, 0.4d0, 1.58d0, 10.0d0 ]
    real(8), parameter, dimension(9) :: it2_ex= [ 0.56d0, 10.0d0, 1.0d6, 1.58d0, 10.0d0, 1.0d6, 1.58d0, 10.0d0, 1.0d6 ]
    real(8), parameter, dimension(9) :: it2_chi= [ 0.13d0, 0.05d0, 0.001d0, 0.08d0, 0.02d0, 0.001d0, 0.2d0, 0.2d0, 0.2d0 ]
    real(8), parameter, dimension(9) :: it2_chipn= [ 0.87d0, 0.95d0, 0.999d0, 0.72d0, 0.78d0, 0.799d0, 1.00d0, 1.00d0, 1.00d0 ]
    real(8), parameter, dimension(9) :: it2_pi0= [ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.1667d0, 0.1667d0, 0.1667d0 ]
    real(8), parameter, dimension(9) :: it2_pip= [ 1.0d0, 1.0d0, 1.0d0, 0.25d0, 0.25d0, 0.25d0, 0.75d0, 0.75d0, 0.75d0 ]
    real(8), parameter, dimension(9) :: it2_pim= [ 0.0d0, 0.0d0, 0.0d0, 0.75d0, 0.75d0, 0.75d0, 0.0833d0, 0.0833d0, 0.0833d0 ]
    real(8), parameter, dimension(9) :: it2_pro= [ 0.000d0, 0.000d0, 0.000d0, 0.833d0, 0.833d0, 0.833d0, 0.0d0, 0.0d0, 0.0d0 ]
    real(8), parameter, dimension(9) :: it2_ntr= [ 1.000d0, 1.000d0, 1.000d0, 0.167d0, 0.167d0, 0.167d0, 0.0d0, 0.0d0, 0.0d0 ]
    real(8), parameter, dimension(14) :: it3_sig= [ 60.0d0, 60.0d0, 85.0d0, 85.0d0, 120.0d0, 120.0d0, 120.0d0, 120.0d0, 120.0d0, 120.0d0, 120.0d0, 120.0d0, 120.0d0, 120.0d0 ]
    real(8), parameter, dimension(14) :: it3_em= [ 0.5d0, 0.5d0, 0.9d0, 0.9d0, 1.5d0, 1.5d0, 5.0d0, 5.0d0, 50.0d0, 50.0d0, 500.0d0, 500.0d0, 5000.0d0, 5000.0d0 ]
    real(8), parameter, dimension(14) :: it3_ex= [ 0.9d0, 0.9d0, 1.5d0, 1.5d0, 5.0d0, 5.0d0, 50.0d0, 50.0d0, 500.0d0, 500.0d0, 5000.0d0, 5000.0d0, 5.0d6, 5.0d6 ]
    real(8), parameter, dimension(14) :: it3_chi= [ 0.1d0, 0.4d0, 0.15d0, 0.35d0, 0.15d0, 0.35d0, 0.07d0, 0.35d0, 0.02d0, 0.5d0, 0.007d0, 0.5d0, 0.002d0, 0.6d0 ]
    real(8), parameter, dimension(14) :: it3_chipn= [ 0.73d0, 1.00d0, 0.66d0, 1.00d0, 0.61d0, 1.00d0, 0.51d0, 1.00d0, 0.55d0, 1.00d0, 0.56d0, 1.00d0, 0.56d0, 1.00d0 ]
    real(8), parameter, dimension(14) :: it3_pi0= [ 0.32d0, 0.17d0, 0.42d0, 0.19d0, 0.59d0, 0.16d0, 1.38d0, 0.16d0, 3.01d0, 0.20d0, 5.13d0, 0.27d0, 7.59d0, 0.26d0 ]
    real(8), parameter, dimension(14) :: it3_pip= [ 0.34d0, 0.29d0, 0.31d0, 0.35d0, 0.57d0, 0.21d0, 1.37d0, 0.25d0, 2.86d0, 0.21d0, 4.68d0, 0.29d0, 6.80d0, 0.27d0 ]
    real(8), parameter, dimension(14) :: it3_pim= [ 0.04d0, 0.05d0, 0.07d0, 0.08d0, 0.30d0, 0.13d0, 1.11d0, 0.23d0, 2.64d0, 0.14d0, 4.57d0, 0.12d0, 6.65d0, 0.13d0 ]
    real(8), parameter, dimension(14) :: it3_pro= [ 0.46d0, 0.00d0, 0.49d0, 0.00d0, 0.65d0, 0.00d0, 0.72d0, 0.00d0, 0.71d0, 0.00d0, 0.72d0, 0.00d0, 0.71d0, 0.00d0 ]
    real(8), parameter, dimension(14) :: it3_ntr= [ 0.54d0, 0.00d0, 0.51d0, 0.00d0, 0.35d0, 0.00d0, 0.28d0, 0.00d0, 0.29d0, 0.00d0, 0.28d0, 0.00d0, 0.29d0, 0.00d0 ]

contains

! Hummer 2010 p-gamma 主算子：生成 pi 源、重子再注入、质子/中子损失和光子 sink。
! Hummer 2010 p-gamma operator: pi sources, baryon reinjection, nucleon losses, and photon sink.
subroutine pg_hummer(ngp,nnu,ehad,hden,eph, &
                                           phden,nden,qpi0, &
                                           qpip,qpim, &
                                           qpro,qntr, &
                                           ploss,nloss,phloss)
    integer, intent(in) :: ngp,nnu
    real(8), intent(in), dimension(ngp) :: ehad,hden,nden
    real(8), intent(in), dimension(nnu) :: eph,phden
    real(8), intent(out), dimension(ngp) :: qpi0,qpip &
        & ,qpim,qpro,qntr &
        & ,ploss,nloss
    real(8), intent(out), dimension(nnu) :: phloss
    integer :: i_gam
    real(8) :: dln_ep, dln_eg, gp, src_p, src_n, ep
    real(8), dimension(ngp) :: plog,nlog
    real(8), dimension(nnu) :: phlog
    real(8), dimension(28) :: conv_res
    real(8), dimension(9) :: conv_dir
    real(8), dimension(14) :: conv_mul

    dln_ep = dlog(ehad(2)/ehad(1))
    dln_eg = dlog(eph(2)/eph(1))
    plog = ehad * hden
    nlog = ehad * nden
    phlog = eph * phden

    qpi0 = 0d0
    qpip = 0d0
    qpim = 0d0
    qpro = 0d0
    qntr = 0d0
    ploss = 0d0
    nloss = 0d0
    phloss = 0d0

    do i_gam=1,ngp
        src_p = plog(i_gam)
        src_n = nlog(i_gam)
        if (src_p <= 0d0 .and. src_n <= 0d0) cycle
        ep = ehad(i_gam)
        gp = ep / proton_m
        call rates_res(gp,nnu,eph,phlog,dln_eg,conv_res)
        call rates_dir(gp,nnu,eph,phlog,dln_eg,conv_dir)
        call rates_mul(gp,nnu,eph,phlog,dln_eg,conv_mul)

        call deposit_pions(28,i_gam,dln_ep,src_p,src_n,it1_chi,it1_pi0,it1_pip,it1_pim, &
                                        it1_pi0,it1_pim,it1_pip,conv_res,qpi0, &
                                        qpip,qpim)
        call deposit_pions(9,i_gam,dln_ep,src_p,src_n,it2_chi,it2_pi0,it2_pip,it2_pim, &
                                        it2_pi0,it2_pim,it2_pip,conv_dir,qpi0, &
                                        qpip,qpim)
        call deposit_pions(14,i_gam,dln_ep,src_p,src_n,it3_chi,it3_pi0,it3_pip,it3_pim, &
                                        it3_pi0,it3_pim,it3_pip,conv_mul,qpi0, &
                                        qpip,qpim)

        call deposit_baryons(28,i_gam,dln_ep,src_p,src_n,it1_chipn,it1_pro,it1_ntr,it1_ntr,it1_pro, &
                                         conv_res,qpro,qntr)
        call deposit_baryons(9,i_gam,dln_ep,src_p,src_n,it2_chipn,it2_pro,it2_ntr,it2_ntr,it2_pro, &
                                         conv_dir,qpro,qntr)
        call deposit_baryons(14,i_gam,dln_ep,src_p,src_n,it3_chipn,it3_pro,it3_ntr,it3_ntr,it3_pro, &
                                         conv_mul,qpro,qntr)

        ploss(i_gam) = loss_res(conv_res) + loss_dir(conv_dir) + &
                                  loss_mul(conv_mul)
        nloss(i_gam) = ploss(i_gam)

        call accum_phloss(gp,src_p+src_n)
    end do

    qpi0 = qpi0 / ehad
    qpip = qpip / ehad
    qpim = qpim / ehad
    qpro = qpro / ehad
    qntr = qntr / ehad

contains

    subroutine accum_phloss(gp,srclog)
        real(8), intent(in) :: gp,srclog
        integer :: i_nu,it
        real(8) :: y

        do i_nu=1,nnu
            y = gp * eph(i_nu)
            do it=1,28
                if (it1_pro(it) > 1d-4 .or. it1_ntr(it) > 1d-4) then
                    phloss(i_nu) = phloss(i_nu) + &
                                             srclog * phloss_res(it,y,dln_ep)
                end if
            end do
            do it=1,9
                if (it2_pro(it) > 1d-4 .or. it2_ntr(it) > 1d-4) then
                    phloss(i_nu) = phloss(i_nu) + &
                                             srclog * phloss_dir(it,y,dln_ep)
                end if
            end do
            do it=1,14
                if (it3_pro(it) > 1d-4 .or. it3_ntr(it) > 1d-4) then
                    phloss(i_nu) = phloss(i_nu) + &
                                             srclog * phloss_mul(it,y,dln_ep)
                end if
            end do
        end do
    end subroutine accum_phloss
end subroutine pg_hummer

! pi 族沉积：按 chi 能量偏移把 p/n 诱发的 pi0/pi+/pi- 多重数投到 hadron grid。
! Pion-family deposition: map p/n-induced pi0/pi+/pi- yields onto the hadron grid by chi shift.
subroutine deposit_pions(nfam,i_gam,dln_ep,src_p,src_n,chiarr,mpi0p,mpipp,mpimp, &
                                      mpi0n,mpipn,mpimn,conv,qpi0,qpip,qpim)
    integer, intent(in) :: nfam,i_gam
    real(8), intent(in), dimension(nfam) :: chiarr,mpi0p,mpipp,mpimp
    real(8), intent(in) :: dln_ep,src_p,src_n
    real(8), intent(in), dimension(nfam) :: mpi0n,mpipn,mpimn,conv
    real(8), intent(inout), dimension(:) :: qpi0,qpip,qpim
    integer :: it
    do it=1,nfam
        if (conv(it) <= 0d0) cycle
        call deposit_shift(qpi0,i_gam,chiarr(it),dln_ep, &
                                         conv(it)*(src_p*mpi0p(it)+src_n*mpi0n(it)))
        call deposit_shift(qpip,i_gam,chiarr(it),dln_ep, &
                                         conv(it)*(src_p*mpipp(it)+src_n*mpipn(it)))
        call deposit_shift(qpim,i_gam,chiarr(it),dln_ep, &
                                         conv(it)*(src_p*mpimp(it)+src_n*mpimn(it)))
    end do
end subroutine deposit_pions

! 重子沉积：按 chi 能量偏移把 p/n 出射分量投到 proton/neutron 再注入谱。
! Baryon deposition: map p/n outgoing components to proton/neutron reinjection spectra by chi shift.
subroutine deposit_baryons(nfam,i_gam,dln_ep,src_p,src_n,chiarr,mprop,mntrp,mpron,mntrn, &
                                       conv,qpro,qntr)
    integer, intent(in) :: nfam,i_gam
    real(8), intent(in), dimension(nfam) :: chiarr,mprop,mntrp
    real(8), intent(in) :: dln_ep,src_p,src_n
    real(8), intent(in), dimension(nfam) :: mpron,mntrn,conv
    real(8), intent(inout), dimension(:) :: qpro,qntr
    integer :: it
    do it=1,nfam
        if (conv(it) <= 0d0) cycle
        call deposit_shift(qpro,i_gam,chiarr(it),dln_ep, &
                                         conv(it)*(src_p*mprop(it)+src_n*mpron(it)))
        call deposit_shift(qntr,i_gam,chiarr(it),dln_ep, &
                                         conv(it)*(src_p*mntrp(it)+src_n*mntrn(it)))
    end do
end subroutine deposit_baryons

! 对数能量沉积 primitive：把一个权重点线性分配到相邻两个能量格点。
! Log-energy deposition primitive: split a weighted packet between adjacent grid cells.
subroutine deposit_shift(target,i_gam,chi,dln_ep,weight)
    real(8), intent(inout), dimension(:) :: target
    integer, intent(in) :: i_gam
    real(8), intent(in) :: chi,dln_ep,weight
    integer :: ileft,iright
    real(8) :: pos,frac
    if (weight <= 0d0 .or. chi <= 0d0) return
    pos = dble(i_gam) + dlog(chi)/dln_ep
    ileft = floor(pos)
    frac = pos - dble(ileft)
    iright = ileft + 1
    if (ileft >= 1 .and. ileft <= size(target)) target(ileft) = target(ileft) + (1d0-frac)*weight
    if (iright >= 1 .and. iright <= size(target)) target(iright) = target(iright) + frac*weight
end subroutine deposit_shift

! 共振区 rate：Delta-resonance 通道对目标光子场卷积。
! Resonance rate: convolve Delta-resonance channels with the target photon field.
subroutine rates_res(gp,nnu,eph,phlog,dln_eg,conv)
    integer, intent(in) :: nnu
    real(8), intent(in), dimension(nnu) :: eph,phlog
    real(8), intent(in) :: gp,dln_eg
    real(8), intent(out), dimension(28) :: conv
    integer :: it,i_nu
    real(8) :: y,sumk
    if (gp <= 0d0) error stop "rates_res: gp must be positive."
    do it=1,28
        sumk = 0d0
        do i_nu=1,nnu
            y = gp * eph(i_nu)
            if (eph(i_nu) <= 0d0) error stop "rates_res: photon energy grid must be positive."
            if (2d0*y >= it1_er(it)) sumk = sumk + phlog(i_nu) / (2d0*y*y)
        end do
        conv(it) = Para_c * mbarn_cm2 * dln_eg * (it1_er(it)*it1_bout(it)*it1_g(it)) * sumk
    end do
end subroutine rates_res

! 直接区 rate：direct pion-production 通道对目标光子场卷积。
! Direct rate: convolve direct pion-production channels with the target photon field.
subroutine rates_dir(gp,nnu,eph,phlog,dln_eg,conv)
    integer, intent(in) :: nnu
    real(8), intent(in), dimension(nnu) :: eph,phlog
    real(8), intent(in) :: gp,dln_eg
    real(8), intent(out), dimension(9) :: conv
    integer :: it,i_nu
    real(8) :: y,sumk
    do it=1,9
        sumk = 0d0
        do i_nu=1,nnu
            y = gp * eph(i_nu)
            sumk = sumk + phlog(i_nu) * kernel_dir(it,y)
        end do
        conv(it) = Para_c * mbarn_cm2 * dln_eg * sumk
    end do
end subroutine rates_dir

! 高能多 pi 区 rate：multipion 通道对目标光子场卷积。
! Multipion rate: convolve high-energy multipion channels with the target photon field.
subroutine rates_mul(gp,nnu,eph,phlog,dln_eg,conv)
    integer, intent(in) :: nnu
    real(8), intent(in), dimension(nnu) :: eph,phlog
    real(8), intent(in) :: gp,dln_eg
    real(8), intent(out), dimension(14) :: conv
    integer :: it,i_nu
    real(8) :: y,sumk
    do it=1,14
        sumk = 0d0
        do i_nu=1,nnu
            y = gp * eph(i_nu)
            sumk = sumk + phlog(i_nu) * kernel_mul(it,y)
        end do
        conv(it) = Para_c * mbarn_cm2 * dln_eg * it3_sig(it) * sumk
    end do
end subroutine rates_mul

! 共振区 nucleon loss：只累加有出射 p/n 的有效通道。
! Resonance nucleon loss: sum only channels with outgoing p/n.
real(8) function loss_res(conv)
    real(8), intent(in), dimension(28) :: conv
    integer :: it
    loss_res = 0d0
    do it=1,28
        if (it1_pro(it) > 1d-4 .or. it1_ntr(it) > 1d-4) loss_res = loss_res + conv(it)
    end do
end function loss_res

! 直接区 nucleon loss：只累加有出射 p/n 的有效通道。
! Direct nucleon loss: sum only channels with outgoing p/n.
real(8) function loss_dir(conv)
    real(8), intent(in), dimension(9) :: conv
    integer :: it
    loss_dir = 0d0
    do it=1,9
        if (it2_pro(it) > 1d-4 .or. it2_ntr(it) > 1d-4) loss_dir = loss_dir + conv(it)
    end do
end function loss_dir

! 高能多 pi 区 nucleon loss：只累加有出射 p/n 的有效通道。
! Multipion nucleon loss: sum only channels with outgoing p/n.
real(8) function loss_mul(conv)
    real(8), intent(in), dimension(14) :: conv
    integer :: it
    loss_mul = 0d0
    do it=1,14
        if (it3_pro(it) > 1d-4 .or. it3_ntr(it) > 1d-4) then
            loss_mul = loss_mul + conv(it)
        end if
    end do
end function loss_mul

! 共振区 photon sink kernel。
! Resonance photon-sink kernel.
real(8) function phloss_res(it,y,dln_ep)
    integer, intent(in) :: it
    real(8), intent(in) :: y,dln_ep
    phloss_res = Para_c * mbarn_cm2 * dln_ep * &
                                         (it1_er(it)*it1_bout(it)*it1_g(it)) * kernel_res(it,y)
end function phloss_res

! 直接区 photon sink kernel。
! Direct photon-sink kernel.
real(8) function phloss_dir(it,y,dln_ep)
    integer, intent(in) :: it
    real(8), intent(in) :: y,dln_ep
    phloss_dir = Para_c * mbarn_cm2 * dln_ep * kernel_dir(it,y)
end function phloss_dir

! 高能多 pi 区 photon sink kernel。
! Multipion photon-sink kernel.
real(8) function phloss_mul(it,y,dln_ep)
    integer, intent(in) :: it
    real(8), intent(in) :: y,dln_ep
    phloss_mul = Para_c * mbarn_cm2 * dln_ep * &
                                         it3_sig(it) * kernel_mul(it,y)
end function phloss_mul

! 共振区 cross-section kernel：阶跃形式 sigma proportional to 1/y^2，阈值 epsilon_r/2。
! Resonance cross-section kernel: step form with sigma proportional to 1/y^2 above epsilon_r/2.
real(8) function kernel_res(it,y)
    integer, intent(in) :: it
    real(8), intent(in) :: y
    if (y <= 0d0) error stop "kernel_res: y must be positive."
    if (2d0*y >= it1_er(it)) then
        kernel_res = 1d0 / (2d0 * y * y)
    else
        kernel_res = 0d0
    end if
end function kernel_res

! 高能多 pi cross-section kernel：分段平台加 1/y^2 尾部。
! Multipion cross-section kernel: piecewise plateau plus 1/y^2 tail.
real(8) function kernel_mul(it,y)
    integer, intent(in) :: it
    real(8), intent(in) :: y
    if (y <= 0d0) error stop "kernel_mul: y must be positive."
    if (y < 0.5d0*it3_em(it)) then
        kernel_mul = 0d0
    else if (y < 0.5d0*it3_ex(it)) then
        kernel_mul = 1d0 - it3_em(it)*it3_em(it)/(4d0*y*y)
    else
        kernel_mul = (it3_ex(it)*it3_ex(it)-it3_em(it)*it3_em(it))/(4d0*y*y)
    end if
end function kernel_mul

! 直接区 cross-section kernel：由分段 I(z) 的差分给出。
! Direct cross-section kernel: difference of the piecewise integral I(z).
real(8) function kernel_dir(it,y)
    integer, intent(in) :: it
    real(8), intent(in) :: y
    real(8) :: zz
    if (y <= 0d0) error stop "kernel_dir: y must be positive."
    if (y < 0.5d0*it2_em(it)) then
        kernel_dir = 0d0
    else if (y < 0.5d0*it2_ex(it)) then
        zz = 2d0*y
        kernel_dir = (idir(it,zz)-idir(it,it2_em(it))) / (2d0*y*y)
    else
        kernel_dir = (idir(it,it2_ex(it))-idir(it,it2_em(it))) / (2d0*y*y)
    end if
end function kernel_dir

! 直接区截面积分 I(z)：低 z 用四次多项式，高 z 用解析近似。
! Direct-region integral I(z): quartic fit at low z and analytic approximation at high z.
real(8) function idir(it,z)
    integer, intent(in) :: it
    real(8), intent(in) :: z
    real(8) :: x
    idir = 0d0
    if (z <= 0d0) error stop "idir: z must be positive."
    x = dlog10(0.5d0*z)
    if (it <= 3) then
        if (z >= 0.17d0 .and. z < 0.96d0) then
            idir = 35.95d0 + 84.08d0*x + 110.76d0*x*x + 102.73d0*x*x*x + 40.47d0*x*x*x*x
        else if (z >= 0.96d0) then
            idir = 30.20d0 + 40.55d0*x + 2.031d0*x*x - 0.3879d0*x*x*x + 0.02504d0*x*x*x*x
        end if
    else
        if (z >= 0.4d0) idir = -3.4083d0 + 16.2864d0/z + 40.7160d0*dlog(z)
    end if
end function idir

end module hadronic_pg
