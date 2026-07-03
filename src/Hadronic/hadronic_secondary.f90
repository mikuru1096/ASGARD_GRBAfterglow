!f2py: skip
module hadronic_secondary
    use constants
    use hadronic_base, only: check_grid
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    real(8), parameter :: c_cgs = Para_c
    real(8), parameter :: sigt_cgs = Para_sigmaT
    real(8), parameter :: megev = Para_m_energy*Para_erg2eV*1d-9
    real(8), parameter :: mpigev = 1.396d-1
    real(8), parameter :: mmugev = 1.0566d-1
    real(8), parameter :: bcrit = 4.41d13
    real(8), parameter :: erggev = 6.24d2

    public :: secondary_calc

contains

! 次级 π/μ 辐射主入口。
! Main entry for pion/muon secondary radiation.
subroutine secondary_calc(nhad,ehad,nph,eph,pip,pim,muml,mumr,mupl,mupr,phden, &
                          iphmin,bfield,pisyn,musyn,piic,muic,dlnsyn,kpion,kmuon, &
                          dlnic,depi,jmpi,demu,jmmu)
    integer, intent(in) :: nhad,nph,iphmin
    real(8), intent(in), dimension(nhad) :: ehad,pip,pim,muml,mumr,mupl,mupr
    real(8), intent(in), dimension(nph) :: eph,phden
    real(8), intent(in) :: bfield
    real(8), intent(out), dimension(nph) :: pisyn,musyn,piic,muic
    real(8), intent(out) :: dlnsyn,dlnic
    real(8), intent(out), dimension(nph,nhad) :: kpion,kmuon
    integer, intent(out), dimension(nhad) :: depi,jmpi,demu,jmmu

    call init_syn(nhad,ehad,nph,eph,bfield,dlnsyn,kpion,kmuon)
    call init_ic(nhad,ehad,nph,eph,iphmin,dlnic,depi,jmpi,demu,jmmu)

    call apply_rad(nhad,pip,pim,muml,mumr,mupl,mupr,nph,phden,dlnsyn,kpion,kmuon, &
                   dlnic,depi,jmpi,demu,jmmu,pisyn,musyn,piic,muic)
end subroutine secondary_calc

! 应用已预计算的同步辐射和 IC 映射。
! Apply precomputed synchrotron and IC maps.
subroutine apply_rad(nhad,pip,pim,muml,mumr,mupl,mupr,nph,phden,dlnsyn,kpion,kmuon, &
                     dlnic,depi,jmpi,demu,jmmu,pisyn,musyn,piic,muic)
    integer, intent(in) :: nhad,nph
    integer, intent(in), dimension(nhad) :: depi,jmpi,demu,jmmu
    real(8), intent(in), dimension(nhad) :: pip,pim,muml,mumr,mupl,mupr
    real(8), intent(in), dimension(nph) :: phden
    real(8), intent(in) :: dlnsyn,dlnic
    real(8), intent(in), dimension(nph,nhad) :: kpion,kmuon
    real(8), intent(out), dimension(nph) :: pisyn,musyn,piic,muic

    call pion_syn(nhad,pip,pim,nph,dlnsyn,kpion,pisyn)
    call muon_syn(nhad,muml,mumr,mupl,mupr,nph,dlnsyn,kmuon,musyn)
    call pion_ic(nph,phden,nhad,pip,pim,dlnic,depi,jmpi,piic)
    call muon_ic(nph,phden,nhad,muml,mumr,mupl,mupr,dlnic,demu,jmmu,muic)
end subroutine apply_rad

! 建立同步辐射核矩阵：每列对应一个次级粒子能格。
! Build synchrotron matrices: each column maps a secondary energy bin.
subroutine init_syn(nhad,ehad,nph,eph,bfield,dln,kpion,kmuon)
    integer, intent(in) :: nhad,nph
    real(8), intent(in), dimension(nhad) :: ehad
    real(8), intent(in), dimension(nph) :: eph
    real(8), intent(in) :: bfield
    real(8), intent(out), dimension(nph,nhad) :: kpion,kmuon
    real(8), intent(out) :: dln
    integer :: i,j

    call check_grid(nhad,ehad,"ehad",dln)
    call check_grid(nph,eph,"eph")
    if (bfield <= 0d0) error stop "secondary synchrotron requires bfield > 0."

    do i=1,nph
        do j=1,nhad
            kpion(i,j) = syn_kernel(eph(i),ehad(j),mpigev,bfield)
            kmuon(i,j) = syn_kernel(eph(i),ehad(j),mmugev,bfield)
        end do
    end do
end subroutine init_syn

! 建立 IC 索引映射：保存能量偏移和每个粒子能格的最高光子格。
! Build IC index maps: store energy offsets and the last photon bin per particle bin.
subroutine init_ic(nhad,ehad,nph,eph,iphmin,dln,depi,jmpi,demu,jmmu)
    integer, intent(in) :: nhad,nph,iphmin
    real(8), intent(in), dimension(nhad) :: ehad
    real(8), intent(in), dimension(nph) :: eph
    real(8), intent(out) :: dln
    integer, intent(out), dimension(nhad) :: depi,jmpi,demu,jmmu
    real(8) :: dlnhad,dlnph

    call check_grid(nhad,ehad,"ehad",dlnhad)
    call check_grid(nph,eph,"eph",dlnph)
    if (dabs(dlnhad-dlnph) > dmax1(1d-12,1d-10*dabs(dlnhad))) then
        error stop "secondary IC requires matched logarithmic spacing."
    end if

    dln = dlnhad
    call build_ic(nhad,ehad,dln,nph,iphmin,mpigev,depi,jmpi)
    call build_ic(nhad,ehad,dln,nph,iphmin,mmugev,demu,jmmu)
end subroutine init_ic

! π 介子同步辐射：总 π 密度乘以同步辐射核矩阵。
! Pion synchrotron: total pion density multiplied by the synchrotron matrix.
subroutine pion_syn(nhad,pip,pim,nph,dln,kpion,pisyn)
    integer, intent(in) :: nhad,nph
    real(8), intent(in), dimension(nhad) :: pip,pim
    real(8), intent(in), dimension(nph,nhad) :: kpion
    real(8), intent(in) :: dln
    real(8), intent(out), dimension(nph) :: pisyn
    real(8), dimension(nhad) :: pitot

    call check_density(nhad,pip,"pip")
    call check_density(nhad,pim,"pim")
    pitot = pip + pim
    pisyn = dln*matmul(kpion,pitot)
end subroutine pion_syn

! μ 子同步辐射：四个手征/电荷分量先合并再投影。
! Muon synchrotron: merge four helicity/charge components before projection.
subroutine muon_syn(nhad,muml,mumr,mupl,mupr,nph,dln,kmuon,musyn)
    integer, intent(in) :: nhad,nph
    real(8), intent(in), dimension(nhad) :: muml,mumr,mupl,mupr
    real(8), intent(in), dimension(nph,nhad) :: kmuon
    real(8), intent(in) :: dln
    real(8), intent(out), dimension(nph) :: musyn
    real(8), dimension(nhad) :: mutot

    call check_density(nhad,muml,"muml")
    call check_density(nhad,mumr,"mumr")
    call check_density(nhad,mupl,"mupl")
    call check_density(nhad,mupr,"mupr")
    mutot = muml + mumr + mupl + mupr
    musyn = dln*matmul(kmuon,mutot)
end subroutine muon_syn

! π 介子 IC：目标光子密度按预计算索引卷积。
! Pion IC: convolve target photons through the precomputed index map.
subroutine pion_ic(nph,phden,nhad,pip,pim,dln,depi,jmpi,piic)
    integer, intent(in) :: nph,nhad
    integer, intent(in), dimension(nhad) :: depi,jmpi
    real(8), intent(in), dimension(nph) :: phden
    real(8), intent(in), dimension(nhad) :: pip,pim
    real(8), intent(in) :: dln
    real(8), intent(out), dimension(nph) :: piic
    real(8), dimension(nhad) :: pitot
    real(8) :: cpi

    call check_density(nph,phden,"phden")
    call check_density(nhad,pip,"pip")
    call check_density(nhad,pim,"pim")
    pitot = pip + pim
    cpi = ic_coeff(mpigev)
    call ic_channel(nph,phden,nhad,pitot,depi,jmpi,dln,cpi,piic)
end subroutine pion_ic

! μ 子 IC：四个 μ 分量共享同一 IC 映射。
! Muon IC: the four muon components use the same IC map.
subroutine muon_ic(nph,phden,nhad,muml,mumr,mupl,mupr,dln,demu,jmmu,muic)
    integer, intent(in) :: nph,nhad
    integer, intent(in), dimension(nhad) :: demu,jmmu
    real(8), intent(in), dimension(nph) :: phden
    real(8), intent(in), dimension(nhad) :: muml,mumr,mupl,mupr
    real(8), intent(in) :: dln
    real(8), intent(out), dimension(nph) :: muic
    real(8), dimension(nhad) :: mutot
    real(8) :: cmu

    call check_density(nph,phden,"phden")
    call check_density(nhad,muml,"muml")
    call check_density(nhad,mumr,"mumr")
    call check_density(nhad,mupl,"mupl")
    call check_density(nhad,mupr,"mupr")
    mutot = muml + mumr + mupl + mupr
    cmu = ic_coeff(mmugev)
    call ic_channel(nph,phden,nhad,mutot,demu,jmmu,dln,cmu,muic)
end subroutine muon_ic

! 超相对论同步辐射核：保留 AM3 分段近似。
! Ultrarelativistic synchrotron kernel: keep the AM3 piecewise fit.
real(8) function syn_kernel(eph,epart,mass,bfield)
    real(8), intent(in) :: eph,epart,mass,bfield
    real(8) :: norm,bdim,mratio,xbar,x2,y,poly

    norm = (3d0*dsqrt(3d0)/pi) * sigt_cgs * c_cgs * bfield * &
           bcrit/(8d0*pi) * erggev/mass
    bdim = bfield/bcrit
    mratio = megev/mass
    xbar = eph*mass/(3d0*epart*epart*bdim*mratio*mratio)
    x2 = 2d0*xbar

    if (x2 < 1d-2) then
        syn_kernel = norm*1.80842d0*xbar**(1d0/3d0)*2d0**(-2d0/3d0)
    else if (x2 < 1d0) then
        y = dlog10(x2)
        poly = -0.35775237d0 - 0.83695385d0*y - 1.1449608d0*y*y - 0.68137283d0*y**3 - &
               0.22754737d0*y**4 - 0.031967334d0*y**5
        syn_kernel = norm*(1d1**poly)/2d0
    else if (x2 < 1d1) then
        y = dlog10(x2)
        poly = -0.35842494d0 - 0.79652041d0*y - 1.6113032d0*y*y + 0.26055213d0*y**3 - &
               1.6979017d0*y**4 + 0.032955035d0*y**5
        syn_kernel = norm*(1d1**poly)/2d0
    else if (x2 < 1d2) then
        syn_kernel = norm*(pi/4d0)*dexp(-x2)*(1d0-99d0/(162d0*x2))
    else
        syn_kernel = 0d0
    end if
end function syn_kernel

! 为一个粒子种类建立 IC 偏移和上界索引。
! Build IC offset and upper-index arrays for a particle species.
subroutine build_ic(nhad,ehad,dln,nph,iphmin,mass,de,jmax)
    integer, intent(in) :: nhad,nph,iphmin
    real(8), intent(in), dimension(nhad) :: ehad
    real(8), intent(in) :: dln,mass
    integer, intent(out), dimension(nhad) :: de,jmax
    integer :: i,cand,jtop
    real(8) :: gamma

    do i=1,nhad
        gamma = ehad(i)/mass
        de(i) = int(dlog(gamma*gamma)/dln)
        jtop = int(dlog(0.5d0*mass/gamma)/dln) + iphmin
        cand = jtop + de(i)
        if (cand > nph) then
            jmax(i) = nph
        else
            jmax(i) = cand
        end if
    end do
end subroutine build_ic

! 单个 IC 通道卷积：粒子能格决定读取哪个目标光子格。
! IC channel convolution: each particle bin selects a target photon bin.
subroutine ic_channel(nph,phden,nhad,hden,de,jmax,dln,coeff,rate)
    integer, intent(in) :: nph,nhad
    integer, intent(in), dimension(nhad) :: de,jmax
    real(8), intent(in), dimension(nph) :: phden
    real(8), intent(in), dimension(nhad) :: hden
    real(8), intent(in) :: dln,coeff
    real(8), intent(out), dimension(nph) :: rate
    integer :: i,j,j0,isrc
    real(8) :: z

    rate = 0d0
    !$OMP PARALLEL DO if(nph*nhad >= 8192) schedule(static) private(i,j,j0,isrc,z)
    do j=1,nph
        j0 = j - 1
        z = 0d0
        do i=1,nhad
            if (j0 < de(i) .or. j0 > jmax(i)) cycle
            isrc = j - de(i)
            if (isrc < 1 .or. isrc > nph) then
                error stop "secondary IC kernel maps outside photon grid."
            end if
            z = z + phden(isrc)*hden(i)
        end do
        rate(j) = z*dln*coeff
    end do
    !$OMP END PARALLEL DO
end subroutine ic_channel

! IC 系数：Thomson 截面和质量缩放。
! IC coefficient: Thomson cross section with mass scaling.
real(8) function ic_coeff(mass)
    real(8), intent(in) :: mass
    real(8) :: mratio

    mratio = mass/megev
    ic_coeff = c_cgs*sigt_cgs/(mratio*mratio)
end function ic_coeff

! 检查内部密度数组有限性，避免 NaN/Inf 进入矩阵投影。
! Check internal density arrays before matrix projection.
subroutine check_density(ngrid,values,name)
    integer, intent(in) :: ngrid
    character(*), intent(in) :: name
    real(8), intent(in), dimension(ngrid) :: values
    integer :: i

    do i=1,ngrid
        if (.not. ieee_is_finite(values(i))) then
            error stop trim(name)//" must contain finite densities."
        end if
    end do
end subroutine check_density

end module hadronic_secondary
