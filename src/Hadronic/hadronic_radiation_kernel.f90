!f2py: skip
module hadronic_radiation_kernel
    use constants
    use hadronic_common
    implicit none
    private

    public :: hadronic_syn_kernel_ultrarel
    public :: hadronic_get_proton_syn_state
    public :: hadronic_syn_polarization_fraction

contains

! 任意带电强子同步辐射核函数：输入光子频率、粒子能量和粒子质量。
real(8) function hadronic_syn_kernel_ultrarel_mass(V_nu,particle_energy_gev,particle_mass_gev,B_field_g)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: V_nu,particle_energy_gev,particle_mass_gev,B_field_g
    real(8), parameter :: B_crit=4.41d13
    real(8) :: mass_e_gev,energy_photon_gev,b_dimless,mass_ratio,xbar,x_arg,y,poly

    mass_e_gev=Para_m_energy*Para_erg2eV*1d-9
    energy_photon_gev=Para_h_GeV*V_nu
    b_dimless=B_field_g/B_crit
    mass_ratio=mass_e_gev/particle_mass_gev
    xbar=energy_photon_gev*particle_mass_gev/(3d0*particle_energy_gev*particle_energy_gev*b_dimless*mass_ratio*mass_ratio)
    x_arg=two*xbar

    if (x_arg < 1d-2) then
        hadronic_syn_kernel_ultrarel_mass=1.80842d0*xbar**(1d0/3d0)*two**(-2d0/3d0)
    else if (x_arg < one) then
        y=dlog10(x_arg)
        poly=-0.35775237d0 &
            -0.83695385d0*y &
            -1.1449608d0*y*y &
            -0.68137283d0*y**3 &
            -0.22754737d0*y**4 &
            -0.031967334d0*y**5
        hadronic_syn_kernel_ultrarel_mass=ten**poly/two
    else if (x_arg < ten) then
        y=dlog10(x_arg)
        poly=-0.35842494d0 &
            -0.79652041d0*y &
            -1.6113032d0*y*y &
            +0.26055213d0*y**3 &
            -1.6979017d0*y**4 &
            +0.032955035d0*y**5
        hadronic_syn_kernel_ultrarel_mass=ten**poly/two
    else if (x_arg < 1d2) then
        hadronic_syn_kernel_ultrarel_mass=(pi/4d0)*dexp(-x_arg)*(one-99d0/(162d0*x_arg))
    else
        hadronic_syn_kernel_ultrarel_mass=zero
    end if
end function hadronic_syn_kernel_ultrarel_mass

! 质子同步辐射超相对论核函数：分段解析/多项式近似计算 F(x)。
real(8) function hadronic_syn_kernel_ultrarel(V_nu,gam_p,B_field_g)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: V_nu,gam_p,B_field_g
    real(8) :: mass_p_gev,energy_proton_gev

    mass_p_gev=Para_m_p_E*Para_erg2eV*1d-9
    energy_proton_gev=gam_p*mass_p_gev
    hadronic_syn_kernel_ultrarel=hadronic_syn_kernel_ultrarel_mass(V_nu,energy_proton_gev,mass_p_gev,B_field_g)
end function hadronic_syn_kernel_ultrarel

! 计算质子同步辐射发射功率和种子光子场，对质子谱卷积同步核函数。
subroutine hadronic_get_proton_syn_state(R_loc,B_field_g,Num_gam_p,Num_nu,gam_p,dN_gam_p,V_seed,P_had_syn,Seed_had_syn)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p,Num_nu
    real(8), intent(in) :: R_loc,B_field_g,gam_p(Num_gam_p),dN_gam_p(Num_gam_p),V_seed(Num_nu)
    real(8), intent(out) :: P_had_syn(Num_nu),Seed_had_syn(Num_nu)
    integer :: I_nu,I_gam
    real(8) :: temp_syn,Rariv2,V_cal,dInteg,Fx
    real(8) :: dN_mid(Num_gam_p),dln_gam(Num_gam_p),temp_para,gam_mid(Num_gam_p)

    temp_syn=dsqrt(3d0)*Para_e*Para_e*Para_e/Para_m_p_E
    if (R_loc <= zero) error stop "hadronic_get_proton_syn_state: radius must be positive."
    if (B_field_g <= zero) error stop "hadronic_get_proton_syn_state: magnetic field must be positive."
    Rariv2=R_loc*R_loc
    P_had_syn=zero
    Seed_had_syn=zero

    do I_gam=1,Num_gam_p-1
        gam_mid(I_gam)=dsqrt(gam_p(I_gam)*gam_p(I_gam+1))
        dN_mid(I_gam)=(dN_gam_p(I_gam)+dN_gam_p(I_gam+1))/two
        dln_gam(I_gam)=dlog(gam_p(I_gam+1)/gam_p(I_gam))
    end do

    do I_nu=1,Num_nu
        V_cal=V_seed(I_nu)
        if (V_cal <= zero) error stop "hadronic_get_proton_syn_state: frequency grid must be positive."
        dInteg=zero
        do I_gam=1,Num_gam_p-1
            Fx=hadronic_syn_kernel_ultrarel(V_cal,gam_mid(I_gam),B_field_g)
            dInteg=dInteg+dN_mid(I_gam)*Fx*gam_mid(I_gam)*dln_gam(I_gam)
        end do
        P_had_syn(I_nu)=temp_syn*B_field_g*dInteg
        Seed_had_syn(I_nu)=P_had_syn(I_nu)/(Rariv2*V_cal)
    end do

    temp_para=4d0*pi*Para_c*Para_h
    Seed_had_syn=Seed_had_syn/temp_para
end subroutine hadronic_get_proton_syn_state

! 输出专用强子同步频率依赖偏振率核：由粒子谱局域斜率和同步发射权重给出Pi_nu。
subroutine hadronic_syn_polarization_fraction(num_had,hadron_energy_gev,density_per_gev,num_ph,photon_frequency_hz, &
                                              particle_mass_gev,magnetic_field_g,p_index,Pi_nu)
    implicit none
    integer, intent(in) :: num_had,num_ph
    real(8), intent(in) :: hadron_energy_gev(num_had),density_per_gev(num_had),photon_frequency_hz(num_ph)
    real(8), intent(in) :: particle_mass_gev,magnetic_field_g,p_index
    real(8), intent(out) :: Pi_nu(num_ph)
    integer :: i_ph,i_had
    real(8) :: dInteg,dPol,e_mid,density_mid,dln_e,Fx,p_eff,pi_eff

    call hadronic_validate_log_grid(num_had,hadron_energy_gev,"hadron_energy_gev",dln_e)
    if (particle_mass_gev <= zero) error stop "hadronic_syn_polarization_fraction requires particle_mass_gev > 0."
    if (magnetic_field_g <= zero) error stop "hadronic_syn_polarization_fraction requires magnetic_field_g > 0."

    do i_ph=1,num_ph
        if (photon_frequency_hz(i_ph) <= zero) error stop "hadronic_syn_polarization_fraction requires positive frequency."
        dInteg=zero
        dPol=zero
        do i_had=1,num_had-1
            density_mid=0.5d0*(density_per_gev(i_had)+density_per_gev(i_had+1))
            if (density_mid <= zero) cycle
            e_mid=dsqrt(hadron_energy_gev(i_had)*hadron_energy_gev(i_had+1))
            Fx=hadronic_syn_kernel_ultrarel_mass(photon_frequency_hz(i_ph),e_mid,particle_mass_gev,magnetic_field_g)
            if (density_per_gev(i_had) > zero .and. density_per_gev(i_had+1) > zero) then
                p_eff=-dlog(density_per_gev(i_had+1)/density_per_gev(i_had))/ &
                       dlog(hadron_energy_gev(i_had+1)/hadron_energy_gev(i_had))
            else
                p_eff=p_index
            end if
            pi_eff=(p_eff+1d0)/(p_eff+7d0/3d0)
            dInteg=dInteg+density_mid*Fx*e_mid*dln_e
            dPol=dPol+density_mid*Fx*pi_eff*e_mid*dln_e
        end do
        if (dInteg > zero) then
            Pi_nu(i_ph)=dPol/dInteg
        else
            Pi_nu(i_ph)=(p_index+1d0)/(p_index+7d0/3d0)
        end if
    end do
end subroutine hadronic_syn_polarization_fraction

end module hadronic_radiation_kernel
