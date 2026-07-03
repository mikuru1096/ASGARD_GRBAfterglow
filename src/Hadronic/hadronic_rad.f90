!f2py: skip
module hadronic_rad
    use constants
    use hadronic_base
    use syn_polarization, only: synchrotron_polarized_components
    implicit none
    private

    public :: syn_kernel
    public :: proton_syn
    public :: syn_polarization

contains

! 同步核无量纲频率 x = nu / nu_c。
! Dimensionless synchrotron-kernel frequency x = nu / nu_c.
real(8) function syn_xarg(vnu,epart,mass,bfield)
    implicit none
    real(8), intent(in) :: vnu,epart,mass,bfield
    real(8), parameter :: bcrit=4.41d13
    real(8) :: megev,eph,bdim,mratio,xbar

    megev=Para_m_energy*Para_erg2eV*1d-9
    eph=Para_h_GeV*vnu
    bdim=bfield/bcrit
    mratio=megev/mass
    xbar=eph*mass/(3d0*epart*epart*bdim*mratio*mratio)
    syn_xarg=2d0*xbar
end function syn_xarg

! 任意带电强子同步辐射核。
! Synchrotron kernel for a charged hadron with arbitrary mass.
real(8) function syn_mass(vnu,epart,mass,bfield)
    implicit none
    real(8), intent(in) :: vnu,epart,mass,bfield
    real(8) :: xarg,y,poly

    xarg=syn_xarg(vnu,epart,mass,bfield)

    if (xarg < 1d-2) then
        syn_mass=1.80842d0*(0.5d0*xarg)**(1d0/3d0)*2d0**(-2d0/3d0)
    else if (xarg < 1d0) then
        y=dlog10(xarg)
        poly=-0.35775237d0 &
            -0.83695385d0*y &
            -1.1449608d0*y*y &
            -0.68137283d0*y**3 &
            -0.22754737d0*y**4 &
            -0.031967334d0*y**5
        syn_mass=1d1**poly/2d0
    else if (xarg < 1d1) then
        y=dlog10(xarg)
        poly=-0.35842494d0 &
            -0.79652041d0*y &
            -1.6113032d0*y*y &
            +0.26055213d0*y**3 &
            -1.6979017d0*y**4 &
            +0.032955035d0*y**5
        syn_mass=1d1**poly/2d0
    else if (xarg < 1d2) then
        syn_mass=(pi/4d0)*dexp(-xarg)*(1d0-99d0/(162d0*xarg))
    else
        syn_mass=0d0
    end if
end function syn_mass

! 质子同步辐射核：先把 gamma 转成 proton energy。
! Proton synchrotron kernel; converts gamma to proton energy first.
real(8) function syn_kernel(vnu,gp,bfield)
    implicit none
    real(8), intent(in) :: vnu,gp,bfield
    real(8) :: mpgev,eprot

    mpgev=Para_m_p_E*Para_erg2eV*1d-9
    eprot=gp*mpgev
    syn_kernel=syn_mass(vnu,eprot,mpgev,bfield)
end function syn_kernel

! 质子同步辐射状态：卷积 proton spectrum，输出 power 和 seed photon density。
! Proton synchrotron state: convolve proton spectrum into power and seed photons.
subroutine proton_syn(rad,bfield,ngp,nnu,gp,dn,nuseed,psyn,seedsyn)
    implicit none
    integer, intent(in) :: ngp,nnu
    real(8), intent(in), dimension(ngp) :: gp,dn
    real(8), intent(in), dimension(nnu) :: nuseed
    real(8), intent(in) :: rad,bfield
    real(8), intent(out), dimension(nnu) :: psyn,seedsyn
    integer :: inu,ig
    real(8) :: pref,r2,vcal,integ,fx
    real(8), dimension(ngp) :: dnmid,dlng,gmid
    real(8) :: norm

    pref=dsqrt(3d0)*Para_e*Para_e*Para_e/Para_m_p_E
    r2=rad*rad
    psyn=0d0
    seedsyn=0d0

    ! 构建 proton-bin midpoint quadrature。 / Build proton-bin midpoint quadrature.
    do ig=1,ngp-1
        gmid(ig)=dsqrt(gp(ig)*gp(ig+1))
        dnmid(ig)=(dn(ig)+dn(ig+1))/2d0
        dlng(ig)=dlog(gp(ig+1)/gp(ig))
    end do
    ! 对每个频率卷积 synchrotron kernel。 / Convolve the synchrotron kernel per frequency.
    do inu=1,nnu
        vcal=nuseed(inu)
        integ=0d0
        do ig=1,ngp-1
            fx=syn_kernel(vcal,gmid(ig),bfield)
            integ=integ+dnmid(ig)*fx*gmid(ig)*dlng(ig)
        end do
        psyn(inu)=pref*bfield*integ
        seedsyn(inu)=psyn(inu)/(r2*vcal)
    end do
    ! 转为 seed photon density。 / Convert power to seed photon density.
    norm=4d0*pi*Para_c*Para_h
    seedsyn=seedsyn/norm
end subroutine proton_syn

! 强子同步偏振率：卷积 perpendicular/parallel 两个偏振发射核。
! Hadronic synchrotron polarization from perpendicular/parallel polarized kernels.
subroutine syn_polarization(nhad,ehad,dens,nph,nuph, &
                                              mass,bfield,pidx,pinu)
    implicit none
    integer, intent(in) :: nhad,nph
    real(8), intent(in), dimension(nhad) :: ehad,dens
    real(8), intent(in), dimension(nph) :: nuph
    real(8), intent(in) :: mass,bfield,pidx
    real(8), intent(out), dimension(nph) :: pinu
    integer :: iph,ihad
    real(8) :: dperp,dpara,ptotal,emid,dmid,dlne,xarg,kperp,kpara

    call check_grid(nhad,ehad,"ehad",dlne)

    do iph=1,nph
        call polar_integral(iph)
        ptotal=dperp+dpara
        if (ptotal > 0d0) then
            pinu(iph)=(dperp-dpara)/ptotal
        else
            pinu(iph)=0d0
        end if
    end do

contains

    subroutine polar_integral(iph)
    implicit none
    integer, intent(in) :: iph

        dperp=0d0
        dpara=0d0
        do ihad=1,nhad-1
            dmid=0.5d0*(dens(ihad)+dens(ihad+1))
            if (dmid <= 0d0) cycle
            emid=dsqrt(ehad(ihad)*ehad(ihad+1))
            xarg=syn_xarg(nuph(iph),emid,mass,bfield)
            call synchrotron_polarized_components(xarg,kperp,kpara)
            dperp=dperp+dmid*kperp*emid*dlne
            dpara=dpara+dmid*kpara*emid*dlne
        end do
    end subroutine polar_integral
end subroutine syn_polarization

end module hadronic_rad
