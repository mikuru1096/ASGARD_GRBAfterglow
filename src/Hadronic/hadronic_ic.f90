!f2py: skip
module hadronic_ic
    use constants
    use hadronic_base, only: check_grid
    implicit none
    private

    real(8), parameter :: c_cgs = Para_c
    real(8), parameter :: sigt_cgs = Para_sigmaT
    real(8), parameter :: megev = Para_m_energy*Para_erg2eV*1d-9
    real(8), parameter :: mpgev = Para_m_p_E*Para_erg2eV*1d-9
    real(8), parameter :: mpigev = 1.396d-1
    real(8), parameter :: mmugev = 1.0566d-1

    public :: ic_init
    public :: ic_operator

contains

! 初始化强子 IC 核：验证网格，并预计算 p/pi/mu 的映射索引。
! Initialize hadronic IC kernels: validate grids and precompute p/pi/mu index maps.
subroutine ic_init(nhad,ehad,nph,eph, &
                                                  iphmin,dln, &
                                                  dep,jmp,depi,jmpi,demu,jmmu)
    integer, intent(in) :: nhad,nph,iphmin
    real(8), intent(in), dimension(nhad) :: ehad
    real(8), intent(in), dimension(nph) :: eph
    real(8), intent(out) :: dln
    integer, intent(out), dimension(nhad) :: dep,jmp,depi,jmpi
    integer, intent(out), dimension(nhad) :: demu,jmmu
    real(8) :: dlnhad,dlnph

    call check_grid(nhad,ehad,"hadronic IC hadron grid",dlnhad)
    call check_grid(nph,eph,"hadronic IC photon grid",dlnph)
    if (dabs(dlnhad-dlnph) > dmax1(1d-12,1d-10*dabs(dlnhad))) then
        error stop "hadronic IC requires hadron/photon grids with the same logarithmic spacing."
    end if

    dln = dlnhad
    call build_ic(nhad,ehad,dln,nph,iphmin,mpgev,dep,jmp)
    call build_ic(nhad,ehad,dln,nph,iphmin,mpigev,depi,jmpi)
    call build_ic(nhad,ehad,dln,nph,iphmin,mmugev,demu,jmmu)
end subroutine ic_init

! 使用预计算索引直接计算强子 IC；调用方已完成网格初始化。
! Apply hadronic IC with precomputed maps; the caller has already initialized the grids.
subroutine ic_operator(nhad,nph,phden, &
                                                     pden,pipden,pimden, &
                                                     mumlden,mumrden, &
                                                     muplden,muprden, &
                                                     dln,dep,jmp,depi,jmpi, &
                                                     demu,jmmu,epsp,epspi,epsmu, &
                                                     cp,cpi,cmu)
    integer, intent(in) :: nhad,nph
    real(8), intent(in), dimension(nph) :: phden
    real(8), intent(in), dimension(nhad) :: pden,pipden,pimden
    real(8), intent(in), dimension(nhad) :: mumlden,mumrden
    real(8), intent(in), dimension(nhad) :: muplden,muprden
    real(8), intent(in) :: dln
    integer, intent(in), dimension(nhad) :: dep,jmp,depi,jmpi
    integer, intent(in), dimension(nhad) :: demu,jmmu
    real(8), intent(out), dimension(nph) :: epsp,epspi,epsmu
    real(8), intent(out) :: cp,cpi,cmu

    call apply_ic(nhad,nph,phden,pden, &
                                           pipden,pimden,mumlden, &
                                           mumrden,muplden,muprden, &
                                           dln,dep,jmp,depi,jmpi,demu,jmmu, &
                                           epsp,epspi,epsmu,cp, &
                                           cpi,cmu)
end subroutine ic_operator

! 分通道累加强子 IC：质子、带电 pi 和 muon 各自使用相同的卷积结构。
! Accumulate hadronic IC by channel: proton, charged pion, and muon share the same convolution.
subroutine apply_ic(nhad,nph,phden,pden, &
                                             pipden,pimden,mumlden, &
                                             mumrden,muplden,muprden, &
                                             dln,dep,jmp,depi,jmpi,demu,jmmu, &
                                             epsp,epspi,epsmu,cp, &
                                             cpi,cmu)
    integer, intent(in) :: nhad,nph
    real(8), intent(in), dimension(nph) :: phden
    real(8), intent(in), dimension(nhad) :: pden,pipden,pimden
    real(8), intent(in), dimension(nhad) :: mumlden,mumrden
    real(8), intent(in), dimension(nhad) :: muplden,muprden
    real(8), intent(in) :: dln
    integer, intent(in), dimension(nhad) :: dep,jmp,depi,jmpi
    integer, intent(in), dimension(nhad) :: demu,jmmu
    real(8), intent(out), dimension(nph) :: epsp,epspi,epsmu
    real(8), intent(out) :: cp,cpi,cmu
    real(8), dimension(nhad) :: hsum

    call apply_species(mpgev,pden,dep,jmp,cp,epsp)

    ! 合并带电 pi 电荷对。 / Merge the charged-pion pair.
    hsum = pipden + pimden
    call apply_species(mpigev,hsum,depi,jmpi,cpi,epspi)

    ! 合并 muon 电荷与螺旋度态。 / Merge muon charge and helicity states.
    hsum = mumlden + mumrden + muplden + muprden
    call apply_species(mmugev,hsum,demu,jmmu,cmu,epsmu)

contains

    subroutine apply_species(mass,hden,de,jmax,coeff,eps)
        real(8), intent(in), dimension(nhad) :: hden
        real(8), intent(in) :: mass
        integer, intent(in), dimension(nhad) :: de,jmax
        real(8), intent(out), dimension(nph) :: eps
        real(8), intent(out) :: coeff

        coeff = ic_coeff(mass)
        call ic_channel(nph,phden,nhad,hden,de,jmax,dln,coeff,eps)
    end subroutine apply_species
end subroutine apply_ic

! 构建单一粒子种类的 IC 映射：de 是能量偏移，jmax 是可贡献的最大光子索引。
! Build the IC map for a single species: de is the energy shift, jmax is the last contributing photon index.
subroutine build_ic(nhad,ehad,dln,nph,iphmin,mass,de,jmax)
    integer, intent(in) :: nhad,nph,iphmin
    real(8), intent(in), dimension(nhad) :: ehad
    real(8), intent(in) :: dln,mass
    integer, intent(out), dimension(nhad) :: de,jmax
    integer :: i,cand,jtop
    real(8) :: gamma

    do i=1,nhad
        if (ehad(i) <= 0d0) then
            error stop "hadronic IC requires positive hadron energies."
        end if
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

! 计算单个强子种类的 IC 通道：在 photon grid 上做离散卷积。
! Compute a single hadronic IC channel: perform the discrete convolution on the photon grid.
subroutine ic_channel(nph,phden,nhad,hden,de,jmax,dln,coeff,eps)
    integer, intent(in) :: nph,nhad
    integer, intent(in), dimension(nhad) :: de,jmax
    real(8), intent(in), dimension(nph) :: phden
    real(8), intent(in), dimension(nhad) :: hden
    real(8), intent(in) :: dln,coeff
    real(8), intent(out), dimension(nph) :: eps
    integer :: i,j,isrc,j0
    real(8) :: z

    eps = 0d0
    !$OMP PARALLEL DO if(nph*nhad >= 8192) schedule(static) private(i,j,isrc,j0,z)
    do j=1,nph
        j0 = j - 1
        z = 0d0
        do i=1,nhad
            if (j0 < de(i) .or. j0 > jmax(i)) cycle
            isrc = j - de(i)
            if (isrc < 1 .or. isrc > nph) then
                error stop "hadronic IC kernel maps to an out-of-grid photon source index."
            end if
            z = z + phden(isrc)*hden(i)
        end do
        eps(j) = z*dln*coeff
    end do
    !$OMP END PARALLEL DO
end subroutine ic_channel

! IC 前因子：sigma_T * c / mratio^2。
! IC prefactor: sigma_T * c / mratio^2.
real(8) function ic_coeff(mass)
    real(8), intent(in) :: mass
    real(8) :: mratio

    mratio = mass/megev
    ic_coeff = c_cgs*sigt_cgs/(mratio*mratio)
end function ic_coeff

end module hadronic_ic
